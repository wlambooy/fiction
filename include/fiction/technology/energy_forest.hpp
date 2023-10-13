//
// Created by Willem Lambooy on 11/10/2023.
//

#ifndef FICTION_ENERGY_FOREST_HPP
#define FICTION_ENERGY_FOREST_HPP

#include "fiction/algorithms/path_finding/distance.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/layouts/bounding_box.hpp"
#include "fiction/layouts/cell_level_layout.hpp"
#include "fiction/technology/cell_ports.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/technology/sidb_nm_position.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/hash.hpp"
#include "fiction/utils/layout_utils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <utility>
#include <variant>

namespace fiction
{

template <typename Lyt>
class energy_forest
{
  public:
    energy_forest(const std::pair<double, double>& bbox_min, const std::pair<double, double>& bbox_max,
                  const std::vector<cell<Lyt>>& sidb_order, const sidb_simulation_parameters& physical_parameters) :
            phys_params{physical_parameters},
            mu_bounds_with_error{phys_params.ef_params.error - phys_params.mu_minus,
                                 -phys_params.ef_params.error - phys_params.mu_minus,
                                 phys_params.ef_params.error - phys_params.mu_plus(),
                                 -phys_params.ef_params.error - phys_params.mu_plus()},
            root_collection{make_root_collection(bbox_min, bbox_max, sidb_order, phys_params)}
    {
        initialize_energy_data_structures();
    }

    void update(const uint64_t w, const uint64_t i, int8_t delta) noexcept
    {
        if (delta == 0)
        {
            return;
        }

        workers[w].local_potential_expressions[i]->update(delta);
    }

    bool check_population_stability(const uint64_t w) noexcept
    {
        if (*workers[w].system_energy - ground_state_energy > 0)
        {
            return false;
        }

        for (std::shared_ptr<local_potential_expression> lpe : workers[w].local_potential_expressions)
        {
            if (!lpe->check_stability(mu_bounds_with_error))
            {
                return false;
            }
        }

        return true;
    }

    inline void challenge_ground_state_energy(const double candidate, std::mutex& mutex) noexcept
    {
        const std::lock_guard lock{mutex};

        const double candidate_with_error = candidate + phys_params.ef_params.energy_error;

        if (candidate_with_error < ground_state_energy)
        {
            ground_state_energy = candidate_with_error;
        }
    }

    inline void reset_ground_state_energy() noexcept
    {
        ground_state_energy = std::numeric_limits<double>::max();
    }

    uint64_t add_worker(std::mutex& mutex)
    {
        worker w{};

        w.system_energy = std::make_shared<double>(*workers[0].system_energy);

        w.local_potential_expressions.reserve(root_collection->sidbs.size());

        for (uint64_t i = 0; i < root_collection->sidbs.size(); ++i)
        {
            w.local_potential_expressions.push_back(workers[0].local_potential_expressions[i]->copy(w.system_energy));

            bool duplicate = false;

            for (const auto& lpe : w.unique_local_potential_expressions)
            {
                duplicate |= w.local_potential_expressions.back() == lpe;
            }

            if (!duplicate)
            {
                w.unique_local_potential_expressions.push_back(w.local_potential_expressions.back());
            }
        }

        // update local potential expression update sets
        std::map<std::shared_ptr<local_potential_expression>, std::shared_ptr<local_potential_expression>> map{};

        for (uint64_t i = 0; i < w.unique_local_potential_expressions.size(); ++i)
        {
            map[workers[0].unique_local_potential_expressions[i]] = w.unique_local_potential_expressions[i];
        }

        for (uint64_t i = 0; i < w.unique_local_potential_expressions.size(); ++i)
        {
            std::set<std::pair<std::shared_ptr<local_potential_expression>, uint64_t>> new_update_set{};

            for (auto& [lpe, j] : w.unique_local_potential_expressions[i]->update_set)
            {
                new_update_set.emplace(map.at(lpe), j);
            }

            w.unique_local_potential_expressions[i]->update_set = new_update_set;
        }

        {
            const std::lock_guard lock{mutex};

            workers.push_back(w);

            return workers.size() - 1;
        }
    }

  protected:
    const sidb_simulation_parameters& phys_params;

    const std::array<double, 4> mu_bounds_with_error;

  private:
    struct sidb_collection;

    struct quadtree
    {
        using nm_pos = std::pair<double, double>;

        struct bounding_box
        {
            const std::pair<double, double> min{};
            const std::pair<double, double> max{};
        };

        const bounding_box bb{};

        std::vector<std::shared_ptr<quadtree>> children{};

        std::set<std::pair<nm_pos, uint64_t>> sidbs{};
        std::pair<double, double>             center_of_mass{};

        std::map<port_direction::cardinal, typename sidb_collection::ptr> hemiquadrants{};

        quadtree(const bounding_box& bbox) : bb{bbox} {}

        [[nodiscard]] inline constexpr uint8_t quadrant_of_pos(const nm_pos& pos) const noexcept
        {
            return (pos.second < (bb.min.second + bb.max.second) / 2 ? 0 : 2) +
                   static_cast<uint8_t>(pos.first >= (bb.min.first + bb.max.first) / 2);
        }

        [[nodiscard]] inline constexpr bounding_box quadrant_bbox(const uint8_t quad) const noexcept
        {
            const double x_mid = (bb.min.first + bb.max.first) / 2;
            const double y_mid = (bb.min.second + bb.max.second) / 2;

            switch (quad)
            {
                case 0: return bounding_box{bb.min, {x_mid, y_mid}};
                case 1: return bounding_box{{x_mid, bb.min.second}, {bb.max.first, y_mid}};
                case 2: return bounding_box{{bb.min.first, y_mid}, {x_mid, bb.max.second}};
                default: return bounding_box{{x_mid, y_mid}, bb.max};
            }
        }

        void insert(const std::pair<nm_pos, uint64_t> pos_ix_pair) noexcept
        {
            const auto& [x, y] = pos_ix_pair.first;

            if (sidbs.empty())
            {
                sidbs          = {pos_ix_pair};
                center_of_mass = {x, y};
                return;
            }

            if (children.empty())
            {
                for (uint8_t i = 0; i < 4; ++i)
                {
                    children.push_back(std::make_shared<quadtree>(quadrant_bbox(i)));
                }

                children[quadrant_of_pos(center_of_mass)]->insert(*sidbs.cbegin());
            }

            sidbs.insert(pos_ix_pair);

            children[quadrant_of_pos({x, y})]->insert(pos_ix_pair);

            center_of_mass = {(center_of_mass.first * static_cast<double>(sidbs.size() - 1) + x) /
                                  static_cast<double>(sidbs.size()),
                              (center_of_mass.second * static_cast<double>(sidbs.size() - 1) + y) /
                                  static_cast<double>(sidbs.size())};
        }

        void make_hemiquadrants(std::shared_ptr<typename sidb_collection::fresh_t> fresh) noexcept
        {
            if (children.empty())
            {
                return;
            }

            for (auto& qt : children)
            {
                qt->make_hemiquadrants(fresh);
            }

            for (uint8_t i = 0; i < 2; i++)
            {
                for (uint8_t j = 0; j < 2; j++)
                {
                    hemiquadrants.insert(
                        {static_cast<port_direction::cardinal>(4 * i + 2 * j),
                         std::make_shared<sidb_collection>(fresh, children[static_cast<uint8_t>(i != j) * (i + 1)],
                                                           children[3 - static_cast<uint8_t>(i == j) * (2 - i)])});
                }
            }
        }
    };

    struct sidb_collection
    {
        using ptr = const std::shared_ptr<sidb_collection>;

        struct fresh_t
        {
            uint64_t fresh{};

            uint64_t get()
            {
                return fresh++;
            }
        };

        std::shared_ptr<fresh_t>        fresh{};
        const uint64_t                  id{};
        const double                    width{};
        const std::set<uint64_t>        sidbs{};
        const std::pair<double, double> center_of_mass{};

        const std::optional<
            std::variant<std::shared_ptr<quadtree>, std::pair<std::shared_ptr<quadtree>, std::shared_ptr<quadtree>>>>
            qt{};

        static inline std::set<uint64_t> get_indices(const std::set<std::pair<std::pair<double, double>, uint64_t>> ps)
        {
            std::set<uint64_t> indices{};

            for (const auto& [_, ix] : ps)
            {
                indices.emplace(ix);
            }

            return indices;
        }

        explicit sidb_collection(std::shared_ptr<fresh_t> fresh_val, std::shared_ptr<quadtree> q) :
                fresh{fresh_val},
                id{fresh->get()},
                width{q->bb.max.first - q->bb.min.first},
                sidbs{get_indices(q->sidbs)},
                center_of_mass{q->center_of_mass},
                qt{q}
        {}

        static inline std::set<uint64_t> get_sidbs_union(const std::set<uint64_t>& s1,
                                                         const std::set<uint64_t>& s2) noexcept
        {
            std::set<uint64_t> sidbs_union{};

            for (const auto& s : {s1, s2})
            {
                sidbs_union.insert(s.cbegin(), s.cend());
            }

            return sidbs_union;
        }

        explicit sidb_collection(std::shared_ptr<fresh_t> fresh_val, std::shared_ptr<quadtree> q1,
                                 std::shared_ptr<quadtree> q2) :
                fresh{fresh_val},
                id{fresh->get()},
                width{std::max(
                    std::max(q1->bb.max.first, q2->bb.max.first) - std::min(q1->bb.min.first, q2->bb.min.first),
                    std::max(q1->bb.max.second, q2->bb.max.second) - std::min(q1->bb.min.second, q2->bb.min.second))},
                sidbs{get_sidbs_union(get_indices(q1->sidbs), get_indices(q2->sidbs))},
                center_of_mass{(q1->center_of_mass.first * static_cast<double>(q1->sidbs.size()) +
                                q2->center_of_mass.first * static_cast<double>(q2->sidbs.size())) /
                                   static_cast<double>(sidbs.size()),
                               (q1->center_of_mass.second * static_cast<double>(q1->sidbs.size()) +
                                q2->center_of_mass.second * static_cast<double>(q2->sidbs.size())) /
                                   static_cast<double>(sidbs.size())},
                qt{std::make_pair(q1, q2)}
        {}

        [[nodiscard]] double dist_to(const sidb_collection& other) const noexcept
        {
            return std::hypot(center_of_mass.first - other.center_of_mass.first,
                              center_of_mass.second - other.center_of_mass.second);
        }

        [[nodiscard]] bool is_well_separated_from(const sidb_collection& other, const double& theta) const noexcept
        {
            return std::min(width, other.width) / dist_to(other) < theta;
        }

        [[nodiscard]] double compute_potential_with(const sidb_collection&            other,
                                                    const sidb_simulation_parameters& ps) const noexcept
        {
            const double r = dist_to(other);

            return ps.k() / (r * 1E-9) * std::exp(-r / ps.lambda_tf) * physical_constants::ELEMENTARY_CHARGE;
        }

        std::optional<std::pair<ptr, ptr>> split() const noexcept
        {
            if (sidbs.size() <= 1 || !qt.has_value())
            {
                return std::nullopt;
            }

            if (!std::holds_alternative<std::shared_ptr<quadtree>>(qt.value()))
            {
                const auto& [q1, q2] =
                    std::get<std::pair<std::shared_ptr<quadtree>, std::shared_ptr<quadtree>>>(qt.value());

                if (q1->sidbs.empty())
                {
                    return std::make_shared<sidb_collection>(fresh, q2)->split();
                }

                if (q2->sidbs.empty())
                {
                    return std::make_shared<sidb_collection>(fresh, q1)->split();
                }

                return std::make_pair(std::make_shared<sidb_collection>(fresh, q1),
                                      std::make_shared<sidb_collection>(fresh, q2));
            }

            const quadtree& q = *std::get<std::shared_ptr<quadtree>>(qt.value());

            for (uint8_t i = 0; i < 8; i += 2)
            {
                if (q.hemiquadrants.at(static_cast<port_direction::cardinal>(i))->sidbs.empty())
                {
                    return q.hemiquadrants.at(static_cast<port_direction::cardinal>((i + 4) % 8))->split();
                }
            }

            const double ns_diff =
                std::hypot(q.hemiquadrants.at(port_direction::cardinal::NORTH)->center_of_mass.first -
                               q.hemiquadrants.at(port_direction::cardinal::SOUTH)->center_of_mass.first,
                           q.hemiquadrants.at(port_direction::cardinal::NORTH)->center_of_mass.second -
                               q.hemiquadrants.at(port_direction::cardinal::SOUTH)->center_of_mass.second);

            const double ew_diff =
                std::hypot(q.hemiquadrants.at(port_direction::cardinal::EAST)->center_of_mass.first -
                               q.hemiquadrants.at(port_direction::cardinal::WEST)->center_of_mass.first,
                           q.hemiquadrants.at(port_direction::cardinal::EAST)->center_of_mass.second -
                               q.hemiquadrants.at(port_direction::cardinal::WEST)->center_of_mass.second);

            if (ns_diff >= ew_diff)
            {
                return std::make_pair(q.hemiquadrants.at(port_direction::cardinal::NORTH),
                                      q.hemiquadrants.at(port_direction::cardinal::SOUTH));
            }

            return std::make_pair(q.hemiquadrants.at(port_direction::cardinal::EAST),
                                  q.hemiquadrants.at(port_direction::cardinal::WEST));
        }
    };

    typename sidb_collection::ptr root_collection;

    struct energy_term
    {
        typename sidb_collection::ptr c1{};
        typename sidb_collection::ptr c2{};
        const double                  chargeless_potential{};

        explicit energy_term(typename sidb_collection::ptr col1, typename sidb_collection::ptr col2,
                             const std::optional<double> chargeless_pot = std::nullopt) :
                c1{col1},
                c2{col2},
                chargeless_potential{chargeless_pot.value()}
        {}

        explicit energy_term(typename sidb_collection::ptr col1, typename sidb_collection::ptr col2,
                             const sidb_simulation_parameters& ps) :
                c1{col1},
                c2{col2},
                chargeless_potential{c1->compute_potential_with(*c2, ps)}
        {}

        std::shared_ptr<energy_term> copy(const bool swap) const noexcept
        {
            return swap ? std::make_shared<energy_term>(c2, c1, chargeless_potential) :
                          std::make_shared<energy_term>(c1, c2, chargeless_potential);
        }
    };

    struct local_potential_expression
    {
        std::shared_ptr<double> system_energy{};

        explicit local_potential_expression(std::shared_ptr<double> system_energy_p) :
                system_energy{std::move(system_energy_p)}
        {}

        std::optional<std::shared_ptr<energy_term>> inner_et{};
        int64_t                                     charge{};
        uint64_t                                    num_sidbs{1};
        std::vector<std::shared_ptr<energy_term>>   ets{};
        double                                      outer_pot{};
        double                                      inner_pot{};

        std::set<std::pair<std::shared_ptr<local_potential_expression>, uint64_t>> update_set{};

        explicit local_potential_expression(
            std::shared_ptr<double> system_energy_copy, std::optional<std::shared_ptr<energy_term>> inner_et_copy,
            int64_t charge_copy, uint64_t num_sidbs_copy, std::vector<std::shared_ptr<energy_term>> ets_copy,
            double outer_pot_copy, double inner_pot_copy,
            std::set<std::pair<std::shared_ptr<local_potential_expression>, uint64_t>> update_set_copy) :
                system_energy{std::move(system_energy_copy)},
                inner_et{inner_et_copy},
                charge{charge_copy},
                num_sidbs{num_sidbs_copy},
                ets{ets_copy},
                outer_pot{outer_pot_copy},
                inner_pot{inner_pot_copy},
                update_set{update_set_copy}
        {}

        void set_inner_energy_term(const uint64_t ix) noexcept
        {
            inner_et  = ets[ix];
            num_sidbs = ets[ix]->c1->sidbs.size() + ets[ix]->c2->sidbs.size();

            ets.erase(std::next(ets.begin(), static_cast<int64_t>(ix)));
        }

        std::set<uint64_t> to_greatest_lower_bound(const sidb_simulation_parameters& ps) noexcept
        {
            uint64_t minimum_size = std::numeric_limits<uint64_t>::max();

            for (const auto& et : ets)
            {
                if (et->c1->sidbs.size() < minimum_size)
                {
                    minimum_size = et->c1->sidbs.size();
                }
            }

            // if precisely one partner is found, set this pair as inner ET and return their union
            std::optional<std::pair<bool, uint64_t>> partner_ix{};

            for (uint64_t i = 0; i < ets.size(); ++i)
            {
                if (ets[i]->c1->sidbs.size() == minimum_size)
                {
                    if (partner_ix.has_value())
                    {
                        partner_ix = {true, i};
                        break;
                    }

                    partner_ix = {false, i};
                }
            }

            if (!partner_ix.value().first)
            {
                set_inner_energy_term(partner_ix.value().second);
            }

            for (auto& loc_et : ets)
            {
                while (loc_et->c1->sidbs.size() > num_sidbs)
                {
                    const auto [cl, cr] = loc_et->c1->split().value();

                    loc_et = std::make_shared<energy_term>(
                        cl->sidbs.count(*ets[partner_ix.value().second]->c1->sidbs.cbegin()) != 0 ? cl : cr, loc_et->c2,
                        ps);
                }
            }

            if (partner_ix.value().first)  // partner_ix must have value
            {
                return std::set<uint64_t>{};
            }

            return sidb_collection::get_sidbs_union(inner_et.value()->c1->sidbs, inner_et.value()->c2->sidbs);
        }

        void external_update(const uint64_t index, const int8_t delta) noexcept
        {
            outer_pot += static_cast<double>(delta) * ets[index]->chargeless_potential;
        }

        void update(const int8_t delta) noexcept
        {
            *system_energy += delta * outer_pot;

            // update internal energy
            if (inner_et.has_value())
            {
                inner_pot += delta * inner_et.value()->chargeless_potential;
            }

            charge += delta;

            // update other local potential expressions
            for (auto& [lpe, i] : update_set)
            {
                lpe->external_update(i, delta);
            }
        }

        bool check_stability(const std::array<double, 4>& mu_bounds) noexcept
        {
            const auto   n           = static_cast<double>(charge);
            const auto   norm_charge = n / static_cast<double>(num_sidbs);

            if (norm_charge == 0 && -outer_pot > mu_bounds[1] && -outer_pot < mu_bounds[2])
            {
                return true;
            }

            const double incr_ratio  = (n + 1) / n;

            if (norm_charge == -1 && -inner_pot * incr_ratio - outer_pot < mu_bounds[0])
            {
                return true;
            }

            const double decr_ratio = (n - 1) / n;

            return (norm_charge == 1 && -inner_pot * decr_ratio - outer_pot > mu_bounds[3]) ||
                   (-inner_pot - outer_pot > mu_bounds[1] && -inner_pot - outer_pot < mu_bounds[2] &&
                    ((norm_charge < 0 && -inner_pot * incr_ratio - outer_pot < mu_bounds[0]) ||
                     -inner_pot * decr_ratio - outer_pot > mu_bounds[3]));
        }

        std::shared_ptr<local_potential_expression> copy(std::shared_ptr<double> system_energy_p) const noexcept
        {
            return std::make_shared<local_potential_expression>(system_energy_p, inner_et, charge, num_sidbs, ets,
                                                                outer_pot, inner_pot, update_set);
        }
    };

    struct worker
    {
        std::vector<std::shared_ptr<local_potential_expression>> local_potential_expressions{};

        std::vector<std::shared_ptr<local_potential_expression>> unique_local_potential_expressions{};

        std::shared_ptr<double> system_energy{};
    };

    std::vector<worker> workers{};

    double ground_state_energy{std::numeric_limits<double>::max()};

    struct potential_tree
    {
        struct node_t
        {
            std::shared_ptr<energy_term>                                                et{};
            std::pair<std::shared_ptr<potential_tree>, std::shared_ptr<potential_tree>> subtrees{};
        };

        std::optional<node_t> node{};

        void walk_tree(const uint64_t i, std::shared_ptr<local_potential_expression> lpe,
                       const sidb_simulation_parameters& ps) noexcept
        {
            if (!node.has_value())
            {
                return;
            }

            bool is_leaf = node.value().subtrees.first == nullptr && node.value().subtrees.second == nullptr;

            if (!is_leaf && !node.value().subtrees.second)
            {
                node.value().subtrees.first->walk_tree(i, lpe, ps);
                return;
            }

            if (!is_leaf && !node.value().subtrees.first)
            {
                node.value().subtrees.second->walk_tree(i, lpe, ps);
                return;
            }

            const bool c1_contains_i = node.value().et->c1->sidbs.count(i) != 0;

            if (!c1_contains_i && node.value().et->c2->sidbs.count(i) == 0)
            {
                return;
            }

            std::shared_ptr<sidb_collection> c_t = c1_contains_i ? node.value().et->c1 : node.value().et->c2;
            std::shared_ptr<sidb_collection> c_b = c1_contains_i ? node.value().et->c2 : node.value().et->c1;

            if ((c_t->sidbs.size() == 1 || c_t->is_well_separated_from(*c_b, ps.ef_params.local_theta)) &&
                (c_b->sidbs.size() == 1 || c_b->is_well_separated_from(*c_t, ps.ef_params.global_theta)))
            {
                lpe->ets.push_back(node.value().et->copy(!c1_contains_i));
                return;
            }

            if (is_leaf)
            {
                if (c_t->sidbs.size() == 1)
                {
                    c_t.swap(c_b);
                }

                const auto [cl, cr] = c_t->split().value();  // split must have value

                node.value().subtrees = {std::make_unique<potential_tree>(), std::make_unique<potential_tree>()};

                node.value().subtrees.first->node =
                    potential_tree::node_t{std::make_shared<energy_term>(cl, c_b, ps),
                                           {std::unique_ptr<potential_tree>{}, std::unique_ptr<potential_tree>{}}};
                node.value().subtrees.second->node =
                    potential_tree::node_t{std::make_shared<energy_term>(cr, c_b, ps),
                                           {std::unique_ptr<potential_tree>{}, std::unique_ptr<potential_tree>{}}};
            }

            node.value().subtrees.first->walk_tree(i, lpe, ps);
            node.value().subtrees.second->walk_tree(i, lpe, ps);
        }
    };

    // V ( c1 , c2 )
    potential_tree potential_between_collections(typename sidb_collection::ptr c1,
                                                 typename sidb_collection::ptr c2) const noexcept
    {
        if (c1->sidbs.empty() || c2->sidbs.empty())
        {
            return potential_tree{};
        }

        typename potential_tree::node_t pot_tree_node{std::make_shared<energy_term>(c1, c2, phys_params)};

        if (c1->is_well_separated_from(*c2, phys_params.ef_params.global_theta) ||
            (c1->sidbs.size() == 1 && c2->sidbs.size() == 1))
        {
            return potential_tree{pot_tree_node};
        }

        const std::array<typename sidb_collection::ptr, 2> a = {c1, c2};
        const std::array<std::optional<std::pair<typename sidb_collection::ptr, typename sidb_collection::ptr>>, 2>
            res = {c1->split(), c2->split()};

        for (uint8_t i = 0; i < 2; ++i)
        {
            if (res[i].has_value() && (!res[1 - i].has_value() || !(res[1 - i].value().first->sidbs.empty() ||
                                                                    res[1 - i].value().second->sidbs.empty())))
            {
                continue;
            }

            pot_tree_node.subtrees = {
                std::make_unique<potential_tree>(potential_between_collections(a[i], res[1 - i].value().first)),
                std::make_unique<potential_tree>(potential_between_collections(a[i], res[1 - i].value().second))};

            return potential_tree{pot_tree_node};
        }

        double  max_diff{};
        uint8_t max_diff_ix{};

        for (uint8_t i = 0; i < 2; ++i)
        {
            std::array<double, 2> dists{};

            for (uint8_t j = 0; j < 2; ++j)
            {
                const sidb_collection& ca = *a[i];
                const sidb_collection& cb = *(j == 0 ? res[1 - i].value().first : res[1 - i].value().second);

                *std::next(dists.begin(), j) = ca.dist_to(cb);
            }

            const double diff = dists[0] - dists[1];

            if (diff > max_diff)
            {
                max_diff    = diff;
                max_diff_ix = i;
            }
        }

        pot_tree_node.subtrees = {std::make_unique<potential_tree>(potential_between_collections(
                                      a[max_diff_ix], res[1 - max_diff_ix].value().first)),
                                  std::make_unique<potential_tree>(potential_between_collections(
                                      a[max_diff_ix], res[1 - max_diff_ix].value().second))};

        return potential_tree{pot_tree_node};
    }

    struct energy_tree
    {
        struct node_t
        {
            potential_tree                                                        pot_tree{};
            std::pair<std::unique_ptr<energy_tree>, std::unique_ptr<energy_tree>> subtrees{};
        };

        std::optional<node_t> node{};

        void construct_lpe(const uint64_t i, std::shared_ptr<local_potential_expression> lpe,
                           const sidb_simulation_parameters& ps) noexcept
        {
            if (!node.has_value())
            {
                return;
            }

            node.value().pot_tree.walk_tree(i, lpe, ps);

            if (node.value().pot_tree.node.has_value())
            {
                node.value().pot_tree.node.value().et->c1->sidbs.count(i) != 0 ?
                    node.value().subtrees.first->construct_lpe(i, lpe, ps) :
                    node.value().subtrees.second->construct_lpe(i, lpe, ps);
                return;
            }

            node.value().subtrees.first->construct_lpe(i, lpe, ps);
            node.value().subtrees.second->construct_lpe(i, lpe, ps);
        }
    };

    energy_tree global_energy_tree{};

    // E ( c )
    energy_tree energy_for_collection(typename sidb_collection::ptr c) const noexcept
    {
        const std::optional<std::pair<typename sidb_collection::ptr, typename sidb_collection::ptr>> res = c->split();

        if (!res.has_value())
        {
            return energy_tree{};
        }

        return energy_tree{
            typename energy_tree::node_t{potential_between_collections(res.value().first, res.value().second),
                                         {std::make_unique<energy_tree>(energy_for_collection(res.value().first)),
                                          std::make_unique<energy_tree>(energy_for_collection(res.value().second))}}};
    }

    void initialize_energy_data_structures() noexcept
    {
        // initialise global energy tree
        global_energy_tree = energy_for_collection(root_collection);
        //        global_energy_tree.set_charges(cell_charge);

        worker w{};

        w.system_energy = std::make_shared<double>(0.0);
        //        w.system_energy = std::make_shared<double>(global_energy_tree.compute_energy());
        //        w.system_energy = std::make_shared<double>(system_energy);

        // initialise local potential expressions
        w.local_potential_expressions.reserve(root_collection->sidbs.size());
        w.unique_local_potential_expressions.reserve(root_collection->sidbs.size());

        std::map<uint64_t, std::shared_ptr<local_potential_expression>> reference_store{};

        for (uint64_t i = 0; i < root_collection->sidbs.size(); ++i)
        {
            if (reference_store.count(i) != 0)
            {
                w.local_potential_expressions.push_back(reference_store.at(i));
                continue;
            }

            std::shared_ptr<local_potential_expression> lpe =
                std::make_shared<local_potential_expression>(w.system_energy);

            global_energy_tree.construct_lpe(i, lpe, phys_params);

            const std::set<uint64_t>& meet = lpe->to_greatest_lower_bound(phys_params);

            for (const uint64_t j : meet)
            {
                if (j == i)
                {
                    continue;
                }

                reference_store[j] = lpe;

                if (j > i)
                {
                    continue;
                }

                w.local_potential_expressions[j] = lpe;
            }

            w.local_potential_expressions.push_back(lpe);
            w.unique_local_potential_expressions.push_back(lpe);
        }

        // initialise update sets
        for (uint64_t i = 0; i < root_collection->sidbs.size(); ++i)
        {
            if (reference_store.count(i) != 0)
            {
                continue;
            }

            for (uint64_t j = 0; j < w.local_potential_expressions[i]->ets.size(); ++j)
            {
                for (const uint64_t k : w.local_potential_expressions[i]->ets[j]->c2->sidbs)
                {
                    if (reference_store.count(k) == 0)
                    {
                        w.local_potential_expressions[k]->update_set.emplace(w.local_potential_expressions[i], j);
                    }
                }
            }
        }

        workers.push_back(w);
    }

    static typename sidb_collection::ptr make_root_collection(const std::pair<double, double>&  bbox_min,
                                                              const std::pair<double, double>&  bbox_max,
                                                              const std::vector<cell<Lyt>>&     sidb_order,
                                                              const sidb_simulation_parameters& ps) noexcept
    {
        std::vector<std::pair<std::pair<double, double>, uint64_t>> pairs{};

        for (uint64_t i = 0; i < sidb_order.size(); ++i)
        {
            pairs.emplace_back(sidb_nm_position<Lyt>(ps, sidb_order[i]), i);
        }

        std::shared_ptr<quadtree> qt = std::make_shared<quadtree>(typename quadtree::bounding_box{bbox_min, bbox_max});

        for (const auto& p : pairs)
        {
            qt->insert(p);
        }

        typename sidb_collection::ptr root_collection =
            std::make_shared<sidb_collection>(std::make_shared<typename sidb_collection::fresh_t>(), qt);

        qt->make_hemiquadrants(root_collection->fresh);

        return root_collection;
    }
};

}  // namespace fiction

#endif  // FICTION_ENERGY_FOREST_HPP

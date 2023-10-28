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
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <utility>
#include <variant>

namespace fiction
{

enum class energy_forest_action
{
    ADD_WORKER,
    RESET
};

template <typename Lyt>
class energy_forest_worker;

template <typename Lyt>
class energy_forest
{
  public:
    energy_forest(const std::pair<double, double>& bbox_min, const std::pair<double, double>& bbox_max,
                  const std::vector<cell<Lyt>>& sidb_order, const sidb_simulation_parameters& physical_parameters) :
            phys_params{physical_parameters},
            mu_bounds_with_error{phys_params.ef_params.stability_error - phys_params.mu_minus,
                                 -phys_params.ef_params.stability_error - phys_params.mu_minus,
                                 phys_params.ef_params.stability_error - phys_params.mu_plus(),
                                 -phys_params.ef_params.stability_error - phys_params.mu_plus()},
            root_collection{std::move(make_root_collection(bbox_min, bbox_max, sidb_order, phys_params))},
            global_energy_tree{energy_for_collection(root_collection)}
    {}

    inline void challenge_ground_state_energy(const double candidate) noexcept
    {
        const std::lock_guard lock{mutex};

        const double candidate_with_error = candidate + phys_params.ef_params.energy_error;

        if (candidate_with_error < ground_state_energy)
        {
            ground_state_energy = candidate_with_error;
        }
    }

    inline void update_phys_params(const sidb_simulation_parameters& params) noexcept
    {
        if (phys_params.mu_minus != params.mu_minus ||
            phys_params.ef_params.stability_error != params.ef_params.stability_error)
        {
            mu_bounds_with_error = {params.ef_params.stability_error - params.mu_minus,
                                    -params.ef_params.stability_error - params.mu_minus,
                                    params.ef_params.stability_error - params.mu_plus(),
                                    -params.ef_params.stability_error - params.mu_plus()};
        }

        if (phys_params.ef_params.energy_error != params.ef_params.energy_error)
        {
            const std::lock_guard lock{mutex};

            ground_state_energy += params.ef_params.energy_error - phys_params.ef_params.energy_error;
        }

        phys_params = params;
    }

    inline void reset_ground_state_energy() noexcept
    {
        const std::lock_guard lock{mutex};

        ground_state_energy = std::numeric_limits<double>::max();
    }

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

        void insert(const std::pair<nm_pos, uint64_t>& pos_ix_pair) noexcept
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

        void make_hemiquadrants() noexcept
        {
            if (children.empty())
            {
                return;
            }

            for (auto& qt : children)
            {
                qt->make_hemiquadrants();
            }

            for (uint8_t i = 0; i < 2; i++)
            {
                for (uint8_t j = 0; j < 2; j++)
                {
                    hemiquadrants.insert(
                        {static_cast<port_direction::cardinal>(4 * i + 2 * j),
                         std::make_shared<sidb_collection>(children[static_cast<uint8_t>(i != j) * (i + 1)],
                                                           children[3 - static_cast<uint8_t>(i == j) * (2 - i)])});
                }
            }
        }
    };

    struct sidb_collection
    {
        using ptr = const std::shared_ptr<sidb_collection>;

        const double                    width{};
        const std::set<uint64_t>        sidbs{};
        const std::pair<double, double> center_of_mass{};

        const std::optional<
            std::variant<std::shared_ptr<quadtree>, std::pair<std::shared_ptr<quadtree>, std::shared_ptr<quadtree>>>>
            qt{};

        static inline std::set<uint64_t> get_indices(const std::set<std::pair<std::pair<double, double>, uint64_t>>& ps)
        {
            std::set<uint64_t> indices{};

            for (const auto& [_, ix] : ps)
            {
                indices.emplace(ix);
            }

            return indices;
        }

        explicit sidb_collection(const std::shared_ptr<quadtree> q) :
                width{q->bb.max.first - q->bb.min.first},
                sidbs{get_indices(q->sidbs)},
                center_of_mass{q->center_of_mass},
                qt{std::move(q)}
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

        explicit sidb_collection(const std::shared_ptr<quadtree> q1, const std::shared_ptr<quadtree> q2) :
                width{1 / std::sqrt(2) *
                      std::hypot(std::max(q1->bb.max.first, q2->bb.max.first) -
                                     std::min(q1->bb.min.first, q2->bb.min.first),
                                 std::max(q1->bb.max.second, q2->bb.max.second) -
                                     std::min(q1->bb.min.second, q2->bb.min.second))},
                sidbs{get_sidbs_union(get_indices(q1->sidbs), get_indices(q2->sidbs))},
                center_of_mass{(q1->center_of_mass.first * static_cast<double>(q1->sidbs.size()) +
                                q2->center_of_mass.first * static_cast<double>(q2->sidbs.size())) /
                                   static_cast<double>(sidbs.size()),
                               (q1->center_of_mass.second * static_cast<double>(q1->sidbs.size()) +
                                q2->center_of_mass.second * static_cast<double>(q2->sidbs.size())) /
                                   static_cast<double>(sidbs.size())},
                qt{std::make_pair(std::move(q1), std::move(q2))}
        {}

        explicit sidb_collection(ptr c1, ptr c2) :
                width{std::min(c1->width, c2->width)},
                sidbs{get_sidbs_union(c1->sidbs, c2->sidbs)},
                center_of_mass{(c1->center_of_mass.first * static_cast<double>(c1->sidbs.size()) +
                                c2->center_of_mass.first * static_cast<double>(c2->sidbs.size())) /
                                   static_cast<double>(sidbs.size()),
                               (c1->center_of_mass.second * static_cast<double>(c1->sidbs.size()) +
                                c2->center_of_mass.second * static_cast<double>(c2->sidbs.size())) /
                                   static_cast<double>(sidbs.size())}
        {}

        [[nodiscard]] double dist_to(const sidb_collection& other) const noexcept
        {
            return std::hypot(center_of_mass.first - other.center_of_mass.first,
                              center_of_mass.second - other.center_of_mass.second);
        }

        [[nodiscard]] bool maximum_is_well_separated_from(const sidb_collection& other, const double theta) const noexcept
        {
            return std::max(width, other.width) / dist_to(other) < theta;
        }

        [[nodiscard]] bool minimum_is_well_separated_from(const sidb_collection& other,
                                                          const double           theta) const noexcept
        {
            return std::min(width, other.width) / dist_to(other) < theta;
        }

        [[nodiscard]] bool is_well_separated_from(const sidb_collection& other,
                                                  const double           theta) const noexcept
        {
            return width / dist_to(other) < theta;
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
                    return std::make_shared<sidb_collection>(q2)->split();
                }

                if (q2->sidbs.empty())
                {
                    return std::make_shared<sidb_collection>(q1)->split();
                }

                return std::make_pair(std::make_shared<sidb_collection>(q1), std::make_shared<sidb_collection>(q2));
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

    struct energy_term
    {
        typename sidb_collection::ptr c1{};
        typename sidb_collection::ptr c2{};
        const double                  chargeless_potential{};

        explicit energy_term(typename sidb_collection::ptr col1, typename sidb_collection::ptr col2,
                             const std::optional<double> chargeless_pot = std::nullopt) :
                c1{std::move(col1)},
                c2{std::move(col2)},
                chargeless_potential{chargeless_pot.value()}
        {}

        explicit energy_term(typename sidb_collection::ptr col1, typename sidb_collection::ptr col2,
                             const sidb_simulation_parameters& ps) :
                c1{std::move(col1)},
                c2{std::move(col2)},
                chargeless_potential{c1->compute_potential_with(*c2, ps)}
        {}

        explicit energy_term(typename sidb_collection::ptr col1, const std::shared_ptr<energy_term> et,
                             const sidb_simulation_parameters& ps) :
                c1{std::move(col1)},
                c2{std::make_shared<sidb_collection>(et->c1, et->c2)},
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
            uint64_t num_sidbs_copy, std::vector<std::shared_ptr<energy_term>> ets_copy,
            std::set<std::pair<std::shared_ptr<local_potential_expression>, uint64_t>> update_set_copy) :
                system_energy{std::move(system_energy_copy)},
                inner_et{inner_et_copy},
                num_sidbs{num_sidbs_copy},
                ets{ets_copy},
                update_set{update_set_copy}
        {}

        void set_inner_energy_term(const uint64_t ix, const sidb_simulation_parameters& ps,
                                   const std::vector<cell<Lyt>>& sidb_order) noexcept
        {
//            inner_et  = ets[ix]; //added
            num_sidbs = ets[ix]->c1->sidbs.size() + ets[ix]->c2->sidbs.size();

            std::vector<std::pair<uint64_t, std::pair<double, double>>> sidb_locs{};

            for (const uint64_t i : sidb_collection::get_sidbs_union(ets[ix]->c1->sidbs, ets[ix]->c2->sidbs))
            {
                sidb_locs.push_back({i, sidb_nm_position<Lyt>(ps, sidb_order[i])});
            }

            double cum_pot = 0;

            for (const auto& [db1, pos1] : sidb_locs)
            {
                for (const auto& [db2, pos2] : sidb_locs)
                {
                    if (db1 >= db2)
                    {
                        continue;
                    }

                    const double r = std::hypot(pos1.first - pos2.first, pos1.second - pos2.second);
                    cum_pot +=
                        ps.k() / (r * 1E-9) * std::exp(-r / ps.lambda_tf) * physical_constants::ELEMENTARY_CHARGE;
                }
            }

            inner_et = std::make_shared<energy_term>(ets[ix]->c1, ets[ix]->c2,
                                                     2 * cum_pot / static_cast<double>(num_sidbs * (num_sidbs - 1)));

            ets.erase(std::next(ets.begin(), static_cast<int64_t>(ix)));
        }

        std::set<uint64_t> compute_greatest_lower_bound(const sidb_simulation_parameters& ps,
                                                        const std::vector<cell<Lyt>>&     sidb_order) noexcept
        {
            uint64_t minimum_size = std::numeric_limits<uint64_t>::max();

            for (const auto& et : ets)
            {
                if (et->c1->sidbs.size() < minimum_size)
                {
                    minimum_size = et->c1->sidbs.size();
                }
            }  // design flaw: minimum_size is always 1.

            // if precisely one partner is found, set this pair as the inner energy term and return their union
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
                set_inner_energy_term(partner_ix.value().second, ps, sidb_order);
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
            outer_pot += delta * ets[index]->chargeless_potential;
        }

        void update(const int8_t delta) noexcept
        {
//            *system_energy += delta * outer_pot;  // delta * inner_pot?

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

        bool check_stability(const std::array<double, 4>& mu_bounds) const noexcept
        {
            const auto n           = static_cast<double>(charge);
            const auto norm_charge = n / static_cast<double>(num_sidbs);

            if (norm_charge == 0 && -outer_pot > mu_bounds[1] && -outer_pot < mu_bounds[2])
            {
                return true;
            }

            const double incr_ratio = (n + 1) / n;

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
            return std::make_shared<local_potential_expression>(system_energy_p, inner_et, num_sidbs, ets, update_set);
        }

        void reset_energy() noexcept
        {
            inner_pot = 0.0;
            outer_pot = 0.0;
            charge    = 0;
        }
    };

    struct potential_tree
    {
        struct node_t
        {
            std::shared_ptr<energy_term>                                                et{};
            std::pair<std::shared_ptr<potential_tree>, std::shared_ptr<potential_tree>> subtrees{};
        };

        std::optional<node_t> node{};

        void walk_tree(const uint64_t i, std::shared_ptr<local_potential_expression>& lpe,
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

            std::shared_ptr<sidb_collection> c_t = std::move(c1_contains_i ? node.value().et->c1 : node.value().et->c2);
            std::shared_ptr<sidb_collection> c_b = std::move(c1_contains_i ? node.value().et->c2 : node.value().et->c1);

            const bool accept_t =
                c_t->sidbs.size() == 1 || c_t->is_well_separated_from(*c_b, ps.ef_params.global_theta);

            if (accept_t && (c_b->sidbs.size() == 1 || c_b->is_well_separated_from(*c_t, ps.ef_params.local_theta)))
            {
                lpe->ets.push_back(node.value().et->copy(!c1_contains_i));
                return;
            }

            if (is_leaf)
            {

                if (accept_t)
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

    struct energy_tree
    {
        struct node_t
        {
            potential_tree                                                        pot_tree{};
            std::pair<std::unique_ptr<energy_tree>, std::unique_ptr<energy_tree>> subtrees{};
        };

        std::optional<node_t> node{};

        void construct_lpe(const uint64_t i, std::shared_ptr<local_potential_expression>& lpe,
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

    sidb_simulation_parameters phys_params;
    std::array<double, 4>      mu_bounds_with_error;

    std::mutex mutex{};

    double ground_state_energy{std::numeric_limits<double>::max()};

    typename sidb_collection::ptr root_collection;

    energy_tree global_energy_tree;

    std::shared_ptr<energy_forest_worker<Lyt>> primary_worker{};

    inline void set_primary_worker(const std::shared_ptr<energy_forest_worker<Lyt>> w) noexcept
    {
        primary_worker = std::move(w);
    }

  private:
    // V ( c1 , c2 )
    potential_tree potential_between_collections(typename sidb_collection::ptr c1,
                                                 typename sidb_collection::ptr c2) const noexcept
    {
        if (c1->sidbs.empty() || c2->sidbs.empty())
        {
            return potential_tree{};
        }

        typename potential_tree::node_t pot_tree_node{std::make_shared<energy_term>(c1, c2, phys_params)};

        if (c2->minimum_is_well_separated_from(*c1, phys_params.ef_params.global_theta) ||
            (c1->sidbs.size() == 1 && c2->sidbs.size() == 1))
        {
            return potential_tree{pot_tree_node};
        }

        const std::array<typename sidb_collection::ptr, 2> a = {std::move(c1), std::move(c2)};
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

        qt->make_hemiquadrants();

        return std::make_shared<sidb_collection>(qt);
    }
};

template <typename Lyt>
class energy_forest_worker
{
  public:
    std::shared_ptr<energy_forest<Lyt>> ef{};

    // constructor for secondary worker
    explicit energy_forest_worker(std::shared_ptr<energy_forest<Lyt>> p) : ef{std::move(p)} {}

    // constructor for primary worker
    explicit energy_forest_worker(const Lyt& lyt, const std::vector<cell<Lyt>>& sidb_order,
                                  const std::vector<sidb_charge_state>& cell_charge,
                                  const sidb_simulation_parameters&     phys_params)
    {
        // create energy forest
        const bounding_box_2d bb{lyt};

        std::pair<double, double> min = sidb_nm_position<Lyt>(phys_params, bb.get_min());
        std::pair<double, double> max = sidb_nm_position<Lyt>(phys_params, bb.get_max());

        if (bb.get_x_size() != bb.get_y_size())
        {
            const double padding = (max.first - min.first - max.second + min.second) / 2;

            if (padding > 0)
            {
                min.second -= padding;
                max.second += padding;
            }
            else
            {
                min.first += padding;
                max.first -= padding;
            }
        }

        ef = std::make_shared<energy_forest<Lyt>>(min, max, sidb_order, phys_params);

        // initialise local potential expressions
        local_potential_expressions.reserve(ef->root_collection->sidbs.size());
        unique_local_potential_expressions.reserve(ef->root_collection->sidbs.size());

        if (ef->root_collection->sidbs.size() <= 2)
        {
            for (uint64_t i = 0; i < ef->root_collection->sidbs.size(); ++i)
            {
                std::shared_ptr<local_pot_expr> lpe = std::make_shared<local_pot_expr>(system_energy);
                ef->global_energy_tree.construct_lpe(i, lpe, ef->phys_params);

                local_potential_expressions.push_back(lpe);
                unique_local_potential_expressions.push_back(lpe);
            }

            initialize_charges(cell_charge);

            return;
        }

        std::map<uint64_t, std::shared_ptr<local_pot_expr>> reference_store{};

        for (uint64_t i = 0; i < ef->root_collection->sidbs.size(); ++i)
        {
            if (reference_store.count(i) != 0)
            {
                local_potential_expressions.push_back(reference_store.at(i));
                continue;
            }

            std::shared_ptr<local_pot_expr> lpe = std::make_shared<local_pot_expr>(system_energy);

            ef->global_energy_tree.construct_lpe(i, lpe, ef->phys_params);

            for (const uint64_t j : lpe->compute_greatest_lower_bound(ef->phys_params, sidb_order))
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

                local_potential_expressions[j] = lpe;
            }

            local_potential_expressions.push_back(lpe);
            unique_local_potential_expressions.push_back(lpe);
        }

        apply_greatest_lower_bounds(reference_store);

        initialize_update_sets(reference_store);

        initialize_charges(cell_charge);
    }

    using local_pot_expr = struct energy_forest<Lyt>::local_potential_expression;

    void update(const uint64_t i, int8_t delta) noexcept
    {
        if (delta == 0)
        {
            return;
        }

        local_potential_expressions[i]->update(delta);
    }

    bool check_population_stability() const noexcept
    {
        //        {
        //            const std::lock_guard lock{ef->mutex};
        //
        //            if (*system_energy - ef->ground_state_energy > 0)
        //            {
        //                return false;
        //            }
        //        }

        for (const auto& lpe : unique_local_potential_expressions)
        {
            if (!lpe->check_stability(ef->mu_bounds_with_error))
            {
                return false;
            }
        }

        return true;
    }

    std::shared_ptr<energy_forest_worker> make_new_worker(const std::vector<sidb_charge_state>& cell_charge)
    {
        std::shared_ptr<energy_forest_worker> w = std::make_shared<energy_forest_worker>(ef);

        w->local_potential_expressions.reserve(ef->root_collection->sidbs.size());

        for (uint64_t i = 0; i < ef->primary_worker->unique_local_potential_expressions.size(); ++i)
        {
            w->unique_local_potential_expressions.push_back(
                ef->primary_worker->unique_local_potential_expressions[i]->copy(w->system_energy));
        }

        for (uint64_t i = 0; i < ef->root_collection->sidbs.size(); ++i)
        {
            uint64_t unique_ix = 0;
            for (; ef->primary_worker->local_potential_expressions[i] !=
                   ef->primary_worker->unique_local_potential_expressions[unique_ix];
                 unique_ix++)
            {}

            w->local_potential_expressions.push_back(w->unique_local_potential_expressions[unique_ix]);
        }

        // update local potential expression update sets
        std::map<std::shared_ptr<local_pot_expr>, std::shared_ptr<local_pot_expr>> map{};

        for (uint64_t i = 0; i < w->unique_local_potential_expressions.size(); ++i)
        {
            map[ef->primary_worker->unique_local_potential_expressions[i]] = w->unique_local_potential_expressions[i];
        }

        for (auto& lpe : w->unique_local_potential_expressions)
        {
            std::set<std::pair<std::shared_ptr<local_pot_expr>, uint64_t>> new_update_set{};

            for (const auto& [update_lpe, i] : lpe->update_set)
            {
                new_update_set.emplace(map.at(update_lpe), i);
            }

            lpe->update_set = new_update_set;
        }

        w->initialize_charges(cell_charge);

        return w;
    }

    void reset_energy_forest(const std::vector<sidb_charge_state>& cell_charge) noexcept
    {
        *system_energy = 0.0;

        for (uint64_t i = 0; i < unique_local_potential_expressions.size(); ++i)
        {
            unique_local_potential_expressions[i]->reset_energy();
        }

        initialize_charges(cell_charge);

        ef->reset_ground_state_energy();
    }

  private:
    void apply_greatest_lower_bounds(const std::map<uint64_t, std::shared_ptr<local_pot_expr>>& ref_store)
    {
        if (unique_local_potential_expressions.size() < 2 || ref_store.empty())
        {
            return;
        }

        for (auto& lpe : unique_local_potential_expressions)
        {
            bool keep_going = true;

            while (keep_going)
            {
                for (uint64_t i = 0; i < lpe->ets.size(); ++i)
                {
                    bool found = false;

                    for (const auto& [_, cluster_lpe] : ref_store)
                    {
                        if (lpe->ets[i]->c2->sidbs == cluster_lpe->inner_et.value()->c1->sidbs)
                        {
                            lpe->ets[i] = std::make_shared<typename energy_forest<Lyt>::energy_term>(
                                lpe->ets[i]->c1, cluster_lpe->inner_et.value(), ef->phys_params);
                        }
                        else if (lpe->ets[i]->c2->sidbs == cluster_lpe->inner_et.value()->c2->sidbs)
                        {
                            lpe->ets.erase(std::next(lpe->ets.begin(), static_cast<int64_t>(i)));
                            found = true;
                            break;
                        }
                    }

                    if (found)
                    {
                        keep_going = true;
                        break;
                    }

                    keep_going = false;
                }
            }
        }
    }

    void initialize_update_sets(const std::map<uint64_t, std::shared_ptr<local_pot_expr>>& ref_store)
    {
        for (const auto& lpe : unique_local_potential_expressions)
        {
            for (uint64_t i = 0; i < lpe->ets.size(); ++i)
            {
                for (const uint64_t j : lpe->ets[i]->c2->sidbs)
                {
                    if (ref_store.count(j) == 0)
                    {
                        local_potential_expressions[j]->update_set.emplace(lpe, i);
                    }
                }
            }
        }
    }

    void initialize_charges(const std::vector<sidb_charge_state>& cell_charge) noexcept
    {
        for (uint64_t i = 0; i < ef->root_collection->sidbs.size(); ++i)
        {
            update(i, charge_state_to_sign(cell_charge[i]));
        }
    }

    std::vector<std::shared_ptr<local_pot_expr>> local_potential_expressions{};
    std::vector<std::shared_ptr<local_pot_expr>> unique_local_potential_expressions{};

    std::shared_ptr<double> system_energy{std::make_shared<double>(0.0)};
};

}  // namespace fiction

#endif  // FICTION_ENERGY_FOREST_HPP
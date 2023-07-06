//
// Created by Willem Lambooy on 06/04/2023.
//

#ifndef FICTION_COMPOSIM_HPP
#define FICTION_COMPOSIM_HPP

#include "fiction/algorithms/simulation/sidb/quicksim.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/hash.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fiction
{

struct composim_params
{
    quicksim_params qs_params{};
    uint64_t        trials{20};
    uint64_t        max_charge_layouts{256};
    uint64_t        max_permutation_logarithm{30};
    bool            debug{false};
};

template <typename Lyt>
struct composim_stats
{
    mockturtle::stopwatch<>::duration time_total{0};

    mockturtle::stopwatch<>::duration time_quicksim{0};

    uint64_t quicksim_attempts{0};

    // TODO: use a variable size bit array to avoid collisions
    std::unordered_map<uint64_t, const charge_distribution_surface<Lyt>*> valid_lyts;

    bool failed{false};

    void report(std::ostream& out = std::cout, std::optional<uint64_t> trials = std::nullopt)
    {
        if (failed)
        {
            out << "\n[w] The log of permutations to try is too large, returning only QuickSim results" << std::endl;
        }

        if (trials.has_value())
        {
            out << "\n[i] QuickSim stats:\n";
            out << fmt::format("[i] Total runtime: {:.2f} s\n", mockturtle::to_seconds(time_quicksim));
            out << fmt::format("[i] Total no. attempts: {}\n", quicksim_attempts);
            out << fmt::format("[i] Average runtime per attempt: {:.2f} ms\n",
                               1000 * mockturtle::to_seconds(time_quicksim) / static_cast<double>(quicksim_attempts));
            out << fmt::format("[i] Average runtime per trial: {:.2f} s\n",
                               mockturtle::to_seconds(time_quicksim) / static_cast<double>(trials.value()));
            out << fmt::format("[i] Average runtime per distinct charge layout: {:.2f} s\n",
                               mockturtle::to_seconds(time_quicksim) / static_cast<double>(valid_lyts.size()));
            out << fmt::format("[i] Average no. attempts per trial: {:.2f}\n",
                               static_cast<double>(quicksim_attempts) / static_cast<double>(trials.value()));
            out << fmt::format("[i] Average no. distinct charge layouts per trial: {:.2f}",
                               static_cast<double>(valid_lyts.size()) / static_cast<double>(trials.value()))
                << std::endl;
        }
        else
        {
            out << fmt::format("\n[i] non-QuickSim runtime: {:.2f} s\n",
                               mockturtle::to_seconds(time_total - time_quicksim));
            out << fmt::format("[i] total runtime: {:.2f} s", mockturtle::to_seconds(time_total)) << std::endl;
        }

        if (!valid_lyts.empty())
        {
            std::vector<const charge_distribution_surface<Lyt>*> sort_by_energy{};

            for (const auto& [_, lyt] : valid_lyts)
            {
                sort_by_energy.emplace_back(lyt);
            }

            std::sort(sort_by_energy.begin(), sort_by_energy.end(),
                      [](const auto* cds_p1, const auto* cds_p2)
                      { return cds_p1->get_system_energy() < cds_p2->get_system_energy(); });

            out << fmt::format("\n[i] lowest energy state: {:.4f} meV \n", sort_by_energy[0]->get_system_energy())
                << std::endl;
            for (const auto& lyt : sort_by_energy)
            {
                out << fmt::format("[i] energy {:.6f}", lyt->get_system_energy()) << std::endl;
            }
        }
        else
        {
            std::cout << "no state found" << std::endl;
        }

        std::cout << "_________________________________________________________" << std::endl;
    }
};

namespace detail
{

template <typename Lyt>
class composim_impl
{
  public:
    composim_impl(const Lyt& lyt, const composim_params& p, composim_stats<Lyt>* pst) :
            layout{lyt},
            ps{p},
            stats_pointer{pst}
    {}

    sidb_simulation_result<Lyt> run()
    {
        sidb_simulation_result<Lyt> res{};

        res.algorithm_name = "CompoSim";
        res.additional_simulation_parameters.emplace_back("iteration_steps", ps.qs_params.interation_steps);
        res.additional_simulation_parameters.emplace_back("alpha", ps.qs_params.alpha);
        res.additional_simulation_parameters.emplace_back("global_potential", ps.qs_params.global_potential);
        res.physical_parameters = ps.qs_params.phys_params;

        // handle empty layouts
        if (layout.num_cells() == 0)
        {
            return res;
        }

        // run CompoSim and track time
        {
            mockturtle::stopwatch stop{st.time_total};
            collect_ground_states();
        }

        if (ps.debug)
        {
            st.report();
        }

        // write simulation runtime and physically valid charge layouts
        if (stats_pointer)
        {
            *stats_pointer = st;
        }

        res.simulation_runtime = st.time_total;

        for (const auto& [_, cds_p] : st.valid_lyts)
        {
            res.charge_distributions.emplace_back(charge_distribution_surface<Lyt>{*cds_p});
        }

        return res;
    }

  private:
    const Lyt& layout;

    const composim_params& ps;

    composim_stats<Lyt>* stats_pointer;

    composim_stats<Lyt> st{};

    using db = std::pair<bool, uint64_t>;

    using midpoint = std::pair<double, double>;

    static inline constexpr midpoint compute_midpoint(const midpoint& x, const midpoint& y) noexcept
    {
        return std::make_pair((x.first + y.first) / 2, (x.second + y.second) / 2);
    }

    using group = std::pair<std::set<db>, midpoint>;

    struct group_hash
    {
        size_t operator()(const group& grp) const
        {
            std::set<uint64_t> indices;
            transform(grp.first.cbegin(), grp.first.cend(), std::inserter(indices, indices.end()),
                      [](const db& p) { return p.second; });
            return std::hash<std::set<uint64_t>>()(indices);
        }
    };

    using groupset = std::unordered_set<group, group_hash>;

    // TODO: log to variable buffer
    void collect_ground_states()
    {
        // track the cumulative QuickSim runtime
        {
            mockturtle::stopwatch stop{st.time_quicksim};

            if (ps.debug)
            {
                std::cout << std::endl;
            }

            // collect different charge layouts produced by QuickSim
            for (const auto i : mockturtle::range(ps.trials))
            {
                auto quicksim_lyts = quicksim<Lyt>(layout, ps.qs_params).charge_distributions;

                uint64_t attempts = 1;

                while (quicksim_lyts.empty())
                {
                    attempts++;

                    const auto& result = quicksim<Lyt>(layout, ps.qs_params);

                    quicksim_lyts.insert(quicksim_lyts.end(), result.charge_distributions.cbegin(),
                                         result.charge_distributions.cend());

                    if (i == 0 && attempts > 5000)
                    {
                        std::cout << "More than 5000 QuickSim attempts failed, try running with different parameters"
                                  << std::endl;
                    }
                }

                for (const auto& lyt : quicksim_lyts)
                {
                    st.valid_lyts[lyt.get_charge_index().first] = new charge_distribution_surface<Lyt>{lyt};
                }

                st.quicksim_attempts += attempts;

                if (ps.debug)
                {
                    std::cout << fmt::format("QuickSim trial {} finished in {} attempt{}", i + 1, attempts,
                                             attempts == 1 ? "" : "s")
                              << std::endl;
                }
            }
        }

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] FOUND {} {}", st.valid_lyts.size(),
                                     st.valid_lyts.size() == 1 ? "QUICKSIM CHARGE LAYOUT" :
                                                                 "DISTINCT QUICKSIM CHARGE LAYOUTS")
                      << std::endl;

            st.report(std::cout, ps.trials);
        }

        // only execute the CompoSim algorithm if multiple distinct charge layouts have been obtained
        if (st.valid_lyts.size() == 1)
        {
            return;
        }

        // use first charge layout for accessing constant properties of the layout
        const auto& fst_cds = *st.valid_lyts.cbegin()->second;

        // initialise vectors that hold indices of DBs that have the same charge state assignment in all layouts
        std::vector<db> const_dbs{};

        // iterate over the DBs in the layout
        for (uint64_t i = 0; i < layout.num_cells(); ++i)
        {
            // get charge state assigned to this DB in the first charge layout
            const auto cs = fst_cds.get_charge_state_by_index(i);

            // check if all other layouts have the same charge state assigned to this DB
            if (std::all_of(++st.valid_lyts.cbegin(), st.valid_lyts.cend(),
                            [&i, &cs](const auto& p) { return p.second->get_charge_state_by_index(i) == cs; }))
            {
                const_dbs.emplace_back(cs == sidb_charge_state::NEGATIVE, i);
                // (assuming QuickSim does not produce charge layouts with positive charges)
            }
        }

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] FOUND {} CONSTANT DB{}", const_dbs.size(), const_dbs.size() == 1 ? "" : "S")
                      << std::endl;
        }

        // All other DBs can form groups with each other; groups are initially considered to be pairs of DBs,
        // where distinct pairs of a DB and its closest non-constant neighbor are added to the set of groups
        // if their charge states are dissimilar in all charge layouts. We annotate the group elements with a flag
        // that is true iff the DB has a negative charge state in the first charge layout.

        // initialise the set of groups, where groups are (lexicographically ordered) sets of (flag, DB index) pairs
        // along with their midpoint coordinate
        groupset groups{};

        // track considered closest neighbors such that we don't add duplicate groups
        std::vector<uint64_t> considered{};

        // iterate over the non-constant and non-considered DBs in the layout
        for (uint64_t i = 0; i < layout.num_cells() - 1; ++i)
        {
            // skip if DB i is labelled constant, or if we already considered this DB as the closest neighbor of DB i' <
            // i
            if (is_const_db(const_dbs, i) || std::find(considered.cbegin(), considered.cend(), i) != considered.cend())
            {
                continue;
            }

            // find the closest non-constant neighboring DB to DB i
            std::pair<double, int64_t> closest_neighbor = {std::numeric_limits<double>::max(), -1};
            for (uint64_t j = i + 1; j < layout.num_cells(); ++j)
            {
                // skip if DB j is labelled constant
                if (is_const_db(const_dbs, j))
                {
                    continue;
                }

                // update closest neighbor if one is found with a distance less than the closest neighbor thus far
                const double dist = fst_cds.get_nm_distance_by_indices(i, j);
                if (dist < closest_neighbor.first)
                {
                    closest_neighbor = {dist, j};
                }
            }

            // move on if no closest non-constant neighbor was found
            if (closest_neighbor.second < 0)
            {
                continue;
            }

            // no need for signedness anymore
            const auto j = static_cast<uint64_t>(closest_neighbor.second);

            // check if they make a group and add them to our set of groups if true
            if (std::all_of(st.valid_lyts.cbegin(), st.valid_lyts.cend(),
                            [&i, &j](const auto& p) {
                                return p.second->get_charge_state_by_index(i) != p.second->get_charge_state_by_index(j);
                            }))
            {
                // obtain flag
                const auto flag = fst_cds.get_charge_state_by_index(i) == sidb_charge_state::NEGATIVE;

                // insert into the set of groups
                groups.emplace(std::make_pair(std::set<db>({std::make_pair(flag, i), std::make_pair(!flag, j)}),
                                              compute_midpoint(fst_cds.get_all_sidb_locations_in_nm()[i],
                                                               fst_cds.get_all_sidb_locations_in_nm()[j])));
            }

            // add DB j to the vector of considered DBs
            considered.emplace_back(j);
        }

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] FOUND {} PAIR{}", groups.size(), groups.size() == 1 ? "" : "S")
                      << std::endl;
        }

        // Now we can create groups of groups, where groups X and Y can be grouped iff for arbitrary elements x in X and
        // y in Y, x has either the same charge state as y in all charge layouts, or not the same charge state in all
        // charge layouts. The annotations let us know which grouped DBs should be assigned the same charge state later.

        reduce_groups(groups);

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] REDUCED TO {} GROUP{}", groups.size(), groups.size() == 1 ? "" : "S")
                      << std::endl;
        }

        // Finally, we try all combinations of group assignments and add the physically valid ones to our store. For
        // each combination, we set the constant DBs and try all possible configurations for the other non-grouped DBs.

        try_all_permutations(groups, const_dbs);
    }

    // TODO: documentation

    static inline constexpr bool is_const_db(const std::vector<db>& const_dbs, const uint64_t& i)
    {
        return std::binary_search(const_dbs.cbegin(), const_dbs.cend(), std::make_pair(false, i),
                                  [](const db& db_i, const db& db_j) { return db_i.second < db_j.second; });
    }

    group merge_groups(const group& g1, const group& g2) const
    {
        std::set<db> merged;
        std::set_union(g1.first.cbegin(), g1.first.cend(), g2.first.cbegin(), g2.first.cend(),
                       std::inserter(merged, merged.begin()));
        return std::make_pair(merged, compute_midpoint(g1.second, g2.second));
    }

    void reduce_groups(groupset& groups) const
    {
        auto groups_size = groups.size() + 1;

        while (groups.size() != groups_size)
        {
            groupset merged_groups;

            std::vector<std::pair<typename groupset::const_iterator, bool>> considered;

            for (auto it1 = groups.cbegin(); it1 != groups.cend(); ++it1)
            {
                const auto considered_it = std::find_if(
                    considered.cbegin(), considered.cend(),
                    [&it1](const std::pair<typename groupset::const_iterator, bool>& p) { return p.first == it1; });
                if (considered_it != considered.cend())
                {
                    if (!considered_it->second)
                    {
                        merged_groups.emplace(*it1);
                    }
                    continue;
                }

                std::pair<double, typename groupset::const_iterator> closest_group{std::numeric_limits<double>::max(),
                                                                                   groups.cend()};
                for (auto it2 = std::next(it1); it2 != groups.cend(); ++it2)
                {
                    const double dist =
                        std::hypot(it1->second.first - it2->second.first, it1->second.second - it2->second.second);
                    if (dist < closest_group.first)
                    {
                        closest_group = std::make_pair(dist, it2);
                    }
                }

                if (closest_group.second == groups.cend())
                {
                    merged_groups.emplace(*it1);
                    continue;
                }

                const auto& pair_i = *it1->first.cbegin();
                const auto& pair_j = *closest_group.second->first.cbegin();
                if (pair_i.second == pair_j.second ||
                    std::all_of(st.valid_lyts.cbegin(), st.valid_lyts.cend(),
                                [&i = pair_i.second, &j = pair_j.second](const auto& p) {
                                    return p.second->get_charge_state_by_index(i) ==
                                           p.second->get_charge_state_by_index(j);
                                }))
                {
                    merged_groups.emplace(merge_groups(*it1, *closest_group.second));
                    considered.emplace_back(closest_group.second, true);
                }
                else
                {
                    merged_groups.emplace(*it1);
                    considered.emplace_back(closest_group.second, false);
                }
            }

            groups_size = groups.size();

            groups = merged_groups;
        }
    }

    static std::set<db> get_grouped_dbs(const groupset& groups)
    {
        std::set<db> indices;
        for (const auto& g : groups)
        {
            std::set_union(indices.begin(), indices.end(), g.first.cbegin(), g.first.cend(),
                           std::inserter(indices, indices.begin()));
        }
        return indices;
    }

    void try_all_permutations(const groupset& groups, const std::vector<db>& const_dbs)
    {
        // get indices of grouped DBs
        const std::set<db>& grouped_dbs = get_grouped_dbs(groups);

        // initialise ordered set of "other" (non-grouped and non-constant) DB indices
        std::set<uint64_t> other;
        auto               k = 0;
        std::generate_n(std::inserter(other, other.begin()), layout.num_cells(), [&k]() { return k++; });

        // remove grouped DB indices from this set
        for (const auto& [_, ix] : grouped_dbs)
        {
            other.erase(ix);
        }

        // remove constant DB indices fromm this set
        for (const auto& [_, ix] : const_dbs)
        {
            other.erase(ix);
        }

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] FOUND {} OTHER DB{}", other.size(), other.size() == 1 ? "" : "S")
                      << std::endl;
        }

        // obtain the 2-log of the number of permutations to try
        const uint64_t n = groups.size() + other.size();

        // check if the number of permutations to try is larger than the given maximum logarithm
        if (n > ps.max_permutation_logarithm)
        {
            st.failed = true;
            return;
        }

        // constrain number of threads when there are more threads available than number of permutations to try
        // also run with one thread if the number of threads is initially set to zero
        const auto     num_bitstrings = 1ul << n;
        const uint64_t num_threads    = std::min(std::max(ps.qs_params.number_threads, 1ul), num_bitstrings);

        // define the bit string ranges per thread
        std::vector<std::pair<uint64_t, uint64_t>> ranges;
        const uint64_t                             chunk_size = std::max(num_bitstrings / num_threads, 1ul);
        uint64_t                                   start      = 0;
        uint64_t                                   end        = chunk_size - 1;

        for (uint64_t i = 0; i < num_threads; ++i)
        {
            ranges.emplace_back(start, end);
            start = end + 1;
            end   = i == num_threads - 2 ? num_bitstrings - 1 : start + chunk_size - 1;
        }

        std::vector<std::thread> threads{};
        threads.reserve(num_threads);
        std::mutex mutex{};

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] TRYING {} BIT STRINGS WITH {} THREAD{} (log = {})", num_bitstrings,
                                     num_threads, num_threads == 1 ? "" : "S", n)
                      << std::endl;
        }

        st.valid_lyts.erase(++st.valid_lyts.cbegin(), st.valid_lyts.cend());
        uint64_t valid_lyt_count = 1;

        std::map<double, uint64_t, std::greater<>> energy_map{
            {st.valid_lyts.cbegin()->second->get_system_energy(),
                            st.valid_lyts.cbegin()->second->get_charge_index().first}};

        const auto max_charge_lyts = std::max(ps.max_charge_layouts, 1ul);

        for (const auto& range : ranges)
        {
            threads.emplace_back(
                [this, &mutex, &range, &n, &groups, &const_dbs, &other, &valid_lyt_count, &energy_map, &max_charge_lyts]
                {
                    // create new charge distribution surface with all charges states set to NEGATIVE
                    charge_distribution_surface<Lyt> cds{*st.valid_lyts.cbegin()->second};

                    // change charge state for the constant neutral DBs
                    for (const auto& [negative, ix] : const_dbs)
                    {
                        if (!negative)
                        {
                            cds.assign_charge_state_by_cell_index(ix, sidb_charge_state::NEUTRAL, false);
                        }
                    }

                    // iterate over all bit strings of length n
                    for (uint64_t i = range.first; i <= range.second; ++i)
                    {
                        // iterate over the bits in the bitstring
                        for (uint64_t j = 0; j < n; j++)
                        {
                            // bit == 1 <-> assign negative charge state
                            const auto negative = static_cast<bool>(i & (1ul << j));

                            // check if this bit changed w.r.t. the previous bit string
                            if (i > range.first && (negative == static_cast<bool>((i - 1ul) & (1ul << j))))
                            {
                                continue;
                            }

                            // assign the charge state to a DB index in the set of other DB indices if the current
                            // bit number is greater than the number of groups, otherwise assign to a group of DBs
                            if (j >= groups.size())
                            {
                                cds.assign_charge_state_by_cell_index(
                                    *std::next(other.cbegin(), static_cast<int64_t>(j - groups.size())),
                                    negative ? sidb_charge_state::NEGATIVE : sidb_charge_state::NEUTRAL, false);
                            }
                            else
                            {
                                for (const auto& [flag, ix] :
                                     std::next(groups.cbegin(), static_cast<int64_t>(j))->first)
                                {
                                    cds.assign_charge_state_by_cell_index(
                                        ix, flag == negative ? sidb_charge_state::NEGATIVE : sidb_charge_state::NEUTRAL,
                                        false);
                                }
                            }
                        }

                        // compute system energy
                        cds.charge_distribution_to_index();
                        cds.update_after_charge_change();

                        // store this charge layout if it is physically valid
                        if (cds.is_physically_valid())
                        {
                            const std::lock_guard lock{mutex};

                            valid_lyt_count++;

                            const auto charge_index = cds.get_charge_index().first;
                            const auto energy       = cds.get_system_energy();

                            if (st.valid_lyts.size() >= max_charge_lyts)
                            {
                                const auto& [highest_energy, highest_energy_charge_index] = *energy_map.cbegin();

                                if (highest_energy > energy)
                                {
                                    st.valid_lyts.erase(highest_energy_charge_index);
                                    energy_map.erase(energy_map.cbegin());

                                    st.valid_lyts[charge_index] = new charge_distribution_surface<Lyt>{cds};
                                    energy_map.emplace(energy, charge_index);
                                }
                            }
                            else
                            {
                                st.valid_lyts[cds.get_charge_index().first] = new charge_distribution_surface<Lyt>{cds};
                                energy_map.emplace(energy, charge_index);
                            }
                        }
                    }
                });
        }

        for (auto& thread : threads)
        {
            thread.join();
        }

        if (ps.debug)
        {
            std::cout << fmt::format("\n[i] FOUND {} VALID LAYOUT{}", valid_lyt_count, valid_lyt_count == 1 ? "" : "S")
                      << std::endl;
        }
    }
};

}  // namespace detail

template <typename Lyt>
sidb_simulation_result<Lyt> composim(const Lyt& lyt, const composim_params& ps = composim_params{},
                                     composim_stats<Lyt>* pst = nullptr)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt must be an SiDB layout");

    return detail::composim_impl<Lyt>{lyt, ps, pst}.run();
}

}  // namespace fiction

#endif  // FICTION_COMPOSIM_HPP

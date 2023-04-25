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

// TODO: remove boost dependency
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/global_fun.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>

namespace fiction
{

namespace detail
{

// TODO: use a variable size bit array to avoid collisions
template <typename Lyt>
uint64_t get_charge_index(const charge_distribution_surface<Lyt>& cds) noexcept
{
    return cds.get_charge_index().first;
}

template <typename Lyt>
double get_system_energy(const charge_distribution_surface<Lyt>& cds) noexcept
{
    return cds.get_system_energy();
}

template <typename Lyt>
using charge_layout_set = boost::multi_index::multi_index_container<
    charge_distribution_surface<Lyt>,
    boost::multi_index::indexed_by<
        boost::multi_index::hashed_unique<
            boost::multi_index::global_fun<
                const charge_distribution_surface<Lyt>&,
                uint64_t,
                &get_charge_index<Lyt>
            >
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::global_fun<
                const charge_distribution_surface<Lyt>&,
                double,
                &get_system_energy<Lyt>
            >
        >
    >
>;

}  // namespace detail

struct composim_params
{
    quicksim_params qs_params{};
    uint64_t        trials{20};
    uint64_t        max_permutation_logarithm{30};
    uint64_t        max_charge_layouts{256};
};

template <typename Lyt>
struct composim_stats
{
    mockturtle::stopwatch<>::duration time_total{0};

    mockturtle::stopwatch<>::duration time_quicksim{0};

    uint64_t quicksim_attempts{0};

    detail::charge_layout_set<Lyt> valid_lyts{};

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
            out << fmt::format("\n[i] lowest energy state: {:.4f} meV \n",
                               valid_lyts.template get<1>().cbegin()->get_system_energy()) << std::endl;
            for (const auto& lyt : valid_lyts.template get<1>())
            {
                out << fmt::format("[i] energy {:.6f}", lyt.get_system_energy()) << std::endl;
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

    void run()
    {
        // handle empty layouts
        if (layout.num_cells() == 0)
        {
            return;
        }

        // run CompoSim and track time
        {
            mockturtle::stopwatch stop{st.time_total};
            collect_ground_states();
        }

        // write simulation runtime and physically valid charge layouts
        if (stats_pointer)
        {
            *stats_pointer = st;
        }
    }

  private:
    const Lyt& layout;

    const composim_params& ps;

    composim_stats<Lyt>* stats_pointer;

    composim_stats<Lyt> st{};

    using db = std::pair<bool, uint64_t>;

    using midpoint = std::pair<double, double>;

    static midpoint compute_midpoint(const midpoint& x, const midpoint& y) noexcept
    {
        return std::make_pair((x.first + y.first) / 2 , (x.second + y.second) / 2);
    }

    using group = std::pair<std::set<db>, midpoint>;

    struct group_hash
    {
        size_t operator() (const group& grp) const
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

            std::cout << std::endl;

            // collect different charge layouts produced by QuickSim
            for (auto i : mockturtle::range(ps.trials))
            {
                uint64_t attempts = 0;

                quicksim_stats<Lyt> quicksim_stats{};

                while (quicksim_stats.valid_lyts.empty())
                {
                    attempts++;

                    quicksim<Lyt>(layout, ps.qs_params, &quicksim_stats);

                    if (i == 0 && attempts > 5000)
                    {
                        std::cout << "More than 5000 QuickSim attempts failed, try running with different parameters"
                                  << std::endl;
                    }
                }

                for (const auto& lyt : quicksim_stats.valid_lyts)
                {
                    lyt.charge_distribution_to_index();
                    st.valid_lyts.insert(lyt);
                }

                st.quicksim_attempts += attempts;
                std::cout << fmt::format("QuickSim trial {} finished in {} attempt{}",
                                         i + 1, attempts, attempts == 1 ? "" : "s") << std::endl;
            }
        }

        std::cout << fmt::format("\n[i] FOUND {} {}", st.valid_lyts.size(),
                                 st.valid_lyts.size() == 1 ? "QUICKSIM CHARGE LAYOUT" :
                                                             "DISTINCT QUICKSIM CHARGE LAYOUTS") << std::endl;
        st.report(std::cout, ps.trials);

        // only execute the CompoSim algorithm if multiple distinct charge layouts have been obtained
        if (st.valid_lyts.size() == 1)
        {
            return;
        }

        // use first charge layout for accessing constant properties of the layout
        auto& fst_cds = *st.valid_lyts.cbegin();

        // initialise vectors that hold indices of DBs that have the same charge state assignment in all layouts
        std::vector<db> const_dbs;

        // iterate over the DBs in the layout
        for (uint64_t i = 0; i < layout.num_cells(); ++i)
        {
            // get charge state assigned to this DB in the first charge layout
            auto cs = fst_cds.get_charge_state_by_index(i);

            // check if all other layouts have the same charge state assigned to this DB
            if (std::all_of(++st.valid_lyts.cbegin(), st.valid_lyts.cend(),
                            [&i, &cs](const charge_distribution_surface<Lyt>& cds)
                            { return cds.get_charge_state_by_index(i) == cs; }))
            {
                const_dbs.push_back(std::make_pair(cs == sidb_charge_state::NEGATIVE, i));
                // (assuming QuickSim does not produce charge layouts with positive charges)
            }
        }

        std::cout << fmt::format("\n[i] FOUND {} CONSTANT DB{}", const_dbs.size(), const_dbs.size() == 1 ? "" : "S")
                  << std::endl;

        // All other DBs can form groups with each other; groups are initially considered to be pairs of DBs,
        // where distinct pairs of a DB and its closest non-constant neighbor are added to the set of groups
        // if their charge states are dissimilar in all charge layouts. We annotate the group elements with a flag
        // that is true iff the DB has a negative charge state in the first charge layout.

        // initialise the set of groups, where groups are (lexicographically ordered) sets of (flag, DB index) pairs
        // along with their midpoint coordinate
        groupset groups;

        // track considered closest neighbors such that we don't add duplicate groups
        std::vector<uint64_t> considered;

        // iterate over the non-constant and non-considered DBs in the layout
        for (uint64_t i = 0; i < layout.num_cells() - 1; ++i)
        {
            // skip if DB i is labelled constant, or if we already considered this DB as the closest neighbor of DB i' < i
            if (is_const_db(const_dbs, i) ||
                std::find(considered.cbegin(), considered.cend(), i) != considered.cend())
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
                double dist = fst_cds.get_nm_distance_by_indices(i, j);
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
            auto j = static_cast<uint64_t>(closest_neighbor.second);

            // check if they make a group and add them to our set of groups if true
            if (std::all_of(st.valid_lyts.cbegin(), st.valid_lyts.cend(),
                            [&i, &j](const auto& cds)
                            { return cds.get_charge_state_by_index(i) != cds.get_charge_state_by_index(j); }))
            {
                // obtain flag
                auto flag = fst_cds.get_charge_state_by_index(i) == sidb_charge_state::NEGATIVE;

                // insert into the set of groups
                groups.insert(std::make_pair(std::set<db>({std::make_pair(flag, i),
                                                           std::make_pair(!flag, j)}),
                                             compute_midpoint(fst_cds.get_all_sidb_locations_in_nm()[i],
                                                              fst_cds.get_all_sidb_locations_in_nm()[j])));
            }

            // add DB j to the vector of considered DBs
            considered.push_back(j);
        }

        std::cout << fmt::format("\n[i] FOUND {} PAIR{}", groups.size(), groups.size() == 1 ? "" : "S") << std::endl;

        // Now we can create groups of groups, where groups X and Y can be grouped iff for arbitrary elements x in X and
        // y in Y, x has either the same charge state as y in all charge layouts, or not the same charge state in all
        // charge layouts. The annotations let us know which grouped DBs should be assigned the same charge state later.

        reduce_groups(groups);

        std::cout << fmt::format("\n[i] REDUCED TO {} GROUP{}", groups.size(), groups.size() == 1 ? "" : "S")
                  << std::endl;

        // Finally, we try all combinations of group assignments and add the physically valid ones to our store. For
        // each combination, we set the constant DBs and try all possible configurations for the other non-grouped DBs.

        try_all_permutations(groups, const_dbs);
    }

    // TODO: documentation

    bool is_const_db(const std::vector<db>& const_dbs, const uint64_t& i)
    {
        return std::binary_search(const_dbs.cbegin(), const_dbs.cend(), std::make_pair(false, i),
                                  [](const db& db_i, const db& db_j) { return db_i.second < db_j.second; });
    }

    group merge_groups (const group& g1, const group& g2)
    {
        std::set<db> merged;
        std::set_union(g1.first.cbegin(), g1.first.cend(), g2.first.cbegin(), g2.first.cend(),
                       std::inserter(merged, merged.begin()));
        return std::make_pair(merged, compute_midpoint(g1.second, g2.second));
    }

    void reduce_groups(groupset& groups)
    {
        auto groups_size = groups.size() + 1;

        while (groups.size() != groups_size)
        {
            groupset merged_groups;

            std::vector<std::pair<typename groupset::const_iterator, bool>> considered;

            for (auto it1 = groups.cbegin(); it1 != groups.cend(); ++it1)
            {
                auto considered_it = std::find_if(considered.cbegin(), considered.cend(),
                                                  [&it1](const std::pair<typename groupset::const_iterator, bool>& p)
                                                  { return p.first == it1; });
                if (considered_it != considered.cend())
                {
                    if (!considered_it->second)
                    {
                        merged_groups.insert(*it1);
                    }
                    continue;
                }

                std::pair<double, typename groupset::const_iterator> closest_group{std::numeric_limits<double>::max(),
                                                                                   groups.cend()};
                for (auto it2 = std::next(it1); it2 != groups.cend(); ++it2)
                {
                    double dist = std::hypot(it1->second.first - it2->second.first,
                                             it1->second.second - it2->second.second);
                    if (dist < closest_group.first)
                    {
                        closest_group = std::make_pair(dist, it2);
                    }
                }

                if (closest_group.second == groups.cend())
                {
                    merged_groups.insert(*it1);
                    continue;
                }

                const auto& pair_i = *it1->first.cbegin();
                const auto& pair_j = *closest_group.second->first.cbegin();
                if (pair_i.second == pair_j.second ||
                    std::all_of(st.valid_lyts.cbegin(), st.valid_lyts.cend(),
                                [&i = pair_i.second, &j = pair_j.second](const auto& cds)
                                { return cds.get_charge_state_by_index(i) == cds.get_charge_state_by_index(j); }))
                {
                    merged_groups.insert(merge_groups(*it1, *closest_group.second));
                    considered.push_back(std::make_pair(closest_group.second, true));
                }
                else
                {
                    merged_groups.insert(*it1);
                    considered.push_back(std::make_pair(closest_group.second, false));
                }

            }

            groups_size = groups.size();

            groups = merged_groups;
        }
    }

    std::set<db> get_grouped_dbs(const groupset& groups)
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
        std::set<db> grouped_dbs = get_grouped_dbs(groups);

        // initialise ordered set of "other" (non-grouped and non-constant) DB indices
        std::set<uint64_t> other;
        auto k = 0;
        std::generate_n(std::inserter(other, other.begin()), layout.num_cells(), [&k](){ return k++; });

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

        std::cout << fmt::format("\n[i] FOUND {} OTHER DB{}", other.size(), other.size() == 1 ? "" : "S") << std::endl;

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
        const auto num_bitstrings = static_cast<uint64_t>(std::pow(2, n));
        const uint64_t num_threads = std::min(std::max(ps.qs_params.number_threads, 1ul), num_bitstrings);

        // define the bit string ranges per thread
        std::vector<std::pair<uint64_t, uint64_t>> ranges;
        const uint64_t chunk_size = std::max(num_bitstrings / num_threads, 1ul);
        uint64_t start = 0;
        uint64_t end = chunk_size - 1;

        for (uint64_t i = 0; i < num_threads; ++i)
        {
            ranges.emplace_back(start, end);
            start = end + 1;
            end = i == num_threads - 2 ? num_bitstrings - 1 : start + chunk_size - 1;
        }

        std::vector<std::thread> threads{};
        threads.reserve(num_threads);
        std::mutex mutex{};

        std::cout << fmt::format("\n[i] TRYING {} BIT STRINGS WITH {} THREAD{} (log = {})",
                                 num_bitstrings, num_threads, num_threads == 1 ? "" : "S", n) << std::endl;

        uint64_t valid_lyt_count = 0;

        for (auto& range : ranges)
        {
            threads.emplace_back(
                [&mutex, &range, &n, &groups, &const_dbs, &other, &valid_lyt_count,
                 &valid_lyts = st.valid_lyts, &max_size = ps.max_charge_layouts]
                {
                    // create new charge distribution surface
                    auto cds = charge_distribution_surface<Lyt>{*valid_lyts.begin()};

                    // reset all charge states to NEUTRAL
                    cds.set_all_charge_states(sidb_charge_state::NEUTRAL);

                    // change charge state for the constant negative DBs
                    for (const auto& [negative, ix] : const_dbs)
                    {
                        if (negative)
                        {
                            cds.assign_charge_by_cell_index(ix, sidb_charge_state::NEGATIVE);
                        }
                    }

                    // iterate over all bit strings of length n
                    for (uint64_t i = range.first; i <= range.second; ++i)
                    {
                        // iterate over the bits in the bitstring
                        for (uint64_t j = 0; j < n; j++)
                        {
                            // bit == 1 <-> assign negative charge state
                            const auto negative = static_cast<bool>(i & (1 << j));

                            // check if this bit changed w.r.t. the previous bit string
                            if (i > range.first && (negative == static_cast<bool>((i - 1) & (1 << j))))
                            {
                                continue;
                            }

                            // assign the charge state to a DB index in the set of other DB indices if the current
                            // bit number is greater than the number of groups, otherwise assign to a group of DBs
                            if (j >= groups.size())
                            {
                                cds.assign_charge_by_cell_index(*std::next(other.cbegin(),
                                                                           static_cast<int64_t>(j - groups.size())),
                                                                negative ? sidb_charge_state::NEGATIVE :
                                                                           sidb_charge_state::NEUTRAL);
                            }
                            else
                            {
                                for (const auto& [flag, ix] : std::next(groups.cbegin(), static_cast<int64_t>(j))->first)
                                {
                                    cds.assign_charge_by_cell_index(ix, flag == negative ? sidb_charge_state::NEGATIVE :
                                                                                           sidb_charge_state::NEUTRAL);
                                }
                            }
                        }

                        // compute system energy
                        cds.update_after_charge_change();

                        // store this charge layout if it is physically valid
                        if (cds.is_physically_valid())
                        {
                            const std::lock_guard lock{mutex};

                            valid_lyt_count++;

                            if (valid_lyts.size() >= max_size)
                            {
                                auto& highest_energy_cds = *valid_lyts.template get<1>().crbegin();
                                if (highest_energy_cds.get_system_energy() > cds.get_system_energy())
                                {
                                    valid_lyts.erase(highest_energy_cds.get_charge_index().first);
                                    valid_lyts.insert(charge_distribution_surface<Lyt>{cds});
                                }
                            }
                            else
                            {
                                valid_lyts.insert(charge_distribution_surface<Lyt>{cds});
                            }
                        }
                    }
                });
        }

        for (auto& thread : threads)
        {
            thread.join();
        }
        
        std::cout << fmt::format("\n[i] FOUND {} VALID LAYOUT{}", valid_lyt_count, valid_lyt_count == 1 ? "" : "S")
                  << std::endl;
    }
};

}  // namespace detail

template <typename Lyt>
void composim(const Lyt& lyt, const composim_params& ps = composim_params{}, composim_stats<Lyt>* pst = nullptr)
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt must be an SiDB layout");

    detail::composim_impl<Lyt>{lyt, ps, pst}.run();
}

}  // namespace fiction

#endif  // FICTION_COMPOSIM_HPP

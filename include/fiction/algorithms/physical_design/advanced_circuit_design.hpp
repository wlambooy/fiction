//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_circuit_operational.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/gate_design_utils.hpp"
#include "fiction/utils/math_utils.hpp"

#include <phmap.h>

#include <mutex>
#include <thread>

#if (PROGRESS_BARS)
#include <mockturtle/utils/progress_bar.hpp>
#endif

#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fiction
{

/**
 * This struct stores the parameters to design an SiDB circuit on a defective surface.
 *
 * @tparam CellLyt SiDB cell-level layout type.
 */
template <typename CellLyt>
struct advanced_circuit_design_params
{
    enum sub_circuit_creation_mode
    {
        CONNECTED_GATES,
        ADJACENT_GATES,
        ALL_GATES
    };
    enum quantization_mode
    {
        VISUALIZATION_ONLY,
        SOFTEN_PRUNING
    };
    /**
     * Parameters for the *exact* placement and routing algorithm.
     */
    exact_physical_design_params exact_design_parameters = {};
    /**
     * This struct holds parameters to design SiDB gates.
     */
    std::pair<cell<CellLyt>, cell<CellLyt>> canvas{};
    uint64_t                                maximum_number_of_canvas_sidbs = 4;
    /** This variable specifies the radius in nanometers around the center of the hexagon where atomic defects are
     * incorporated into the gate design. (unit: nm)
     */

    sidb_simulation_parameters sim_params{};

    std::optional<CellLyt> defect_surface{};
    double                 influence_radius_charged_defects = 15;

    uint64_t num_trials = 100000000000000;

    double            quantization_step = 7.5;
    quantization_mode quantize_mode     = quantization_mode::SOFTEN_PRUNING;

    double selectivity           = 0.5;
    double selectivity_tolerance = 0.5;

    double success_rate_ceiling = 0.95;

    uint64_t maximum_repeated_discrimination_attempts = 5;
    uint64_t maximum_discrimination_attempts          = 15;

    sub_circuit_creation_mode sub_circuit_mode = sub_circuit_creation_mode::CONNECTED_GATES;

    uint64_t available_threads = std::thread::hardware_concurrency();

    bool print_found_circuits = false;
};

/**
 * Statistics for the on-the-fly defect-aware circuit design.
 */
template <typename GateLyt>
struct advanced_circuit_design_stats
{
    /**
     * The total runtime of the operational domain computation.
     */
    mockturtle::stopwatch<>::duration time_total{0};
    /**
     * The `stats` of the *exact* algorithm.
     */
    exact_physical_design_stats exact_stats{};
};

namespace detail
{

template <typename CellLyt, typename GateLyt>
class advanced_circuit_design_impl
{
    using sidb_count_map_t = std::unordered_map<mockturtle::node<GateLyt>, std::pair<uint8_t, uint8_t>>;

  public:
    advanced_circuit_design_impl(const GateLyt& gate_lyt, const sidb_count_map_t& sidb_count_map,
                                 const advanced_circuit_design_params<CellLyt>& design_params,
                                 advanced_circuit_design_stats<GateLyt>&        st) :
            gate_layout{gate_lyt},
            sidb_counts{sidb_count_map},
            params{design_params},
            stats{st}
    {
        operational_params.simulation_parameters = params.sim_params;

        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;
        // operational_params.input_bdl_iterator_params =
        //     params.design_gate_params.operational_params.input_bdl_iterator_params;
    }

    [[nodiscard]] std::vector<CellLyt> design_sidb_layouts(const uint64_t num_gates_to_design)
    {
        // initialize
        collect_initial_gate_designs();

        circuit->tighten_gate_influence_bounds_until_fixpoint();

        while (circuit_design_level < num_gates_to_design - 1)
        {
            // prune by assessing gate design combinations for increasingly large sets of connected gates
            if (!prune_gate_designs(
                    circuit_design_level,
                    std::string_view{
                        std::to_string(circuit_design_level + 1) +
                        (circuit_design_level == num_gates_to_design - 1 ?
                             " | GLOBAL" :
                             (circuit_design_level == 0 ?
                                  " GATE" :
                                  (params.sub_circuit_mode == advanced_circuit_design_params<
                                                                  CellLyt>::sub_circuit_creation_mode::CONNECTED_GATES ?
                                       " CONNECTED GATES" :
                                   params.sub_circuit_mode == advanced_circuit_design_params<
                                                                  CellLyt>::sub_circuit_creation_mode::ADJACENT_GATES ?
                                       " ADJACENT GATES" :
                                       " GATES")))}))
            {
                circuit_design_level = 0;

                continue;
            }

            ++circuit_design_level;
        }

        const std::vector<CellLyt>& result = exhaustively_enumerate_gate_design_combinations();

        std::cout << result.size() << " OPERATIONAL CIRCUITS FOUND" << std::endl;

        return result;
    }

    [[nodiscard]] std::vector<CellLyt> design_circuit()
    {
        const mockturtle::stopwatch stop{stats.time_total};

        circuit.emplace(gate_layout, params.sim_params, sidb_counts);

        uint64_t number_of_gates_to_design = 0;

        gate_layout.foreach_node(
            [&](const mockturtle::node<GateLyt>& n)
            {
                if (!gate_layout.is_constant(n))
                {
                    number_of_gates_to_design++;
                }
            });

        return design_sidb_layouts(number_of_gates_to_design);
    }

  private:
    /**
     * Gate-level layout.
     */
    const GateLyt gate_layout;
    /**
     *
     */
    const sidb_count_map_t sidb_counts;
    /**
     * Parameters for the on-the-fly circuit design.
     */
    const advanced_circuit_design_params<CellLyt> params{};
    /**
     * Statistics for the on-the-fly circuit design.
     */
    advanced_circuit_design_stats<GateLyt>&           stats;
    std::optional<sidb_bdl_circuit<CellLyt, GateLyt>> circuit{};

    is_circuit_operational_params operational_params{};

    template <typename T>
    using foreach_node = std::unordered_map<mockturtle::node<GateLyt>, T>;

    uint64_t circuit_design_level = 0;

    using gate_design_t = typename sidb_bdl_circuit<CellLyt, GateLyt>::canvas_combination;

    using gate_designs_per_node = foreach_node<std::vector<gate_design_t>>;

    struct SetHash
    {
        std::size_t operator()(const std::pair<mockturtle::node<GateLyt>, gate_design_t>& p) const
        {
            std::size_t seed = 0;

            hash_combine(seed, static_cast<uint64_t>(p.first));

            for (const uint64_t v : p.second)
            {
                hash_combine(seed, v);
            }

            return seed;
        }
    };

    struct SetSetHash
    {
        std::size_t
        operator()(const phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash>& s) const
        {
            std::size_t seed = 0;
            for (const auto& [e1, e2] : s)
            {
                hash_combine(seed, e1);
                hash_combine(seed, e2);
            }
            return seed;
        }
    };

    phmap::flat_hash_set<phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash>, SetSetHash>
        excluded_combinations{};

    std::mutex excluded_combinations_mutex;

    void collect_initial_gate_designs()
    {
        gate_designs_per_node gate_designs{};

        try
        {
            gate_layout.foreach_node(
                [&, this](const auto& n, [[maybe_unused]] auto i)
                {
                    if (gate_layout.is_constant(n))
                    {
                        return;
                    }

                    std::cout << gate_layout.get_tile(n) << std::endl;

                    std::vector<gate_design_t> all_combinations;

                    // // Reserve an estimated capacity if possible (optional optimization) todo
                    // all_combinations.reserve(estimate_total_combinations(max_sidbs, total_positions));

                    for (std::size_t num_sidbs = sidb_counts.at(n).first; num_sidbs <= sidb_counts.at(n).second;
                         ++num_sidbs)
                    {
                        auto combinations = determine_all_combinations_of_distributing_k_entities_on_n_positions(
                            num_sidbs, circuit->canvasses.at(n).positions.size());
                        std::move(combinations.begin(), combinations.end(), std::back_inserter(all_combinations));
                    }

                    std::cout << "tile: " << gate_layout.get_tile(n)
                              << "\t|\tNUM COMBINATIONS: " << all_combinations.size() << std::endl;

                    std::shuffle(all_combinations.begin(), all_combinations.end(),
                                 std::mt19937(std::random_device()()));

                    circuit->gate_designs[n] = std::move(all_combinations);
                });
        }

        catch (const gate_design_exception<tt, GateLyt>& _)
        {
            throw unsuccessful_gate_design_error("Gate design was unsuccessful");
        }
    }

    void collect_nodes(const tile<GateLyt>&                             t,
                       phmap::flat_hash_set<mockturtle::node<GateLyt>>& collected_nodes) const noexcept
    {
        const mockturtle::node<GateLyt>& n = gate_layout.get_node(t);

        uint8_t i = 3;

        // std::cout << "\n\n\nstarting for " << t << std::endl;

        while (true)
        {
            const auto search_in_radius = [&]
            {
                int8_t x = 0;
                int8_t y = 0;

                bool    b = true;
                uint8_t j = 0;

                while (j++ < i)
                {
                    ++(b ? x : y);

                    b = !b;
                }

                std::vector<mockturtle::node<GateLyt>> surrounding_nodes{};

                for (int8_t xi = t.x < x ? -t.x : -x; xi <= x; ++xi)
                {
                    for (int8_t yi = t.y < y ? -t.y : -y; yi <= y; ++yi)
                    {
                        if (const mockturtle::node<GateLyt>& other_n = gate_layout.get_node({t.x + xi, t.y + yi});
                            other_n != 0 && other_n != n)
                        {
                            surrounding_nodes.push_back(other_n);
                        }
                    }
                }

                return surrounding_nodes;
            };

            if (const std::vector<mockturtle::node<GateLyt>>& surrounding_nodes = search_in_radius();
                surrounding_nodes.size() > 0)
            {
                for (const mockturtle::node<GateLyt>& other_n : surrounding_nodes)
                {
                    // std::cout << "collected: " << gate_layout.get_tile(other_n) << std::endl;
                    collected_nodes.emplace(other_n);
                }

                return;
            }

            i += 1;
        }
    }

    struct VectorHash
    {
        template <typename T>
        std::size_t operator()(const std::vector<T>& v) const
        {
            std::size_t seed = v.size();
            for (const auto& i : v)
            {
                seed ^= std::hash<T>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    using sub_circuit_map =
        std::unordered_map<std::vector<mockturtle::node<GateLyt>>, sidb_bdl_sub_circuit<CellLyt, GateLyt>, VectorHash>;

    using sub_circuit_ptr_map =
        std::unordered_map<std::vector<mockturtle::node<GateLyt>>, sidb_bdl_sub_circuit<CellLyt, GateLyt>*, VectorHash>;

    void build_sub_circuits(const phmap::flat_hash_set<mockturtle::node<GateLyt>>& root_collection,
                            std::vector<mockturtle::node<GateLyt>>&                path,
                            phmap::flat_hash_set<mockturtle::node<GateLyt>>& visited, const size_t depth,
                            const size_t max_depth, sub_circuit_map& result) const
    {
        if (depth >= max_depth)
        {
            // // todo: does this work for disjoint regions?
            // for (const auto& node : visited)
            // {
            //     if (gate_layout.is_wire(node))
            //     {
            //         continue;
            //     }
            //
            //     // make sure to collect the whole logic clump  todo: does this work for logic clumps with a y spacer?
            //     // std::cout << "collecting other logic clump components" << std::endl;
            //     // collect_nodes(gate_layout.get_tile(node), visited);
            //
            //     // make sure we have an observable output
            //     bool observable_output_found = false;
            //
            //     for (const auto& other_node : visited)
            //     {
            //         if (gate_layout.is_wire(other_node) &&
            //             gate_layout.get_tile(other_node).y > gate_layout.get_tile(node).y)
            //         {
            //             observable_output_found = true;
            //
            //             break;
            //         }
            //     }
            //
            //     if (!observable_output_found)
            //     {
            //         // std::cout << "additionally collecting observable output" << std::endl;
            //         collect_nodes(gate_layout.get_tile(*std::max_element(
            //                           visited.cbegin(), visited.cend(), [&](const auto& n1, const auto& n2)
            //                           { return gate_layout.get_tile(n1).y < gate_layout.get_tile(n2).y; })),
            //                       visited);
            //     }
            //
            //     break;
            // }

            // check duplicate path

            for (const auto& [other_path, _] : result)
            {
                if (std::all_of(other_path.cbegin(), other_path.cend(),
                                [&](const auto& n) { return visited.count(n) > 0; }))
                {
                    return;
                }
            }

            // path.clear();
            //
            // for (const auto& n : visited)
            // {
            //     path.emplace_back(n);
            // }

            if (path.size() <= max_depth + 1 &&
                std::any_of(path.cbegin(), path.cend(), [&](const auto& n) { return gate_layout.is_wire(n); }))
            {
                result.insert({path, sidb_bdl_sub_circuit<CellLyt, GateLyt>{std::cref(*circuit), path}});
            }

            return;
        }

        const tile<GateLyt> current_tile = gate_layout.get_tile(path.back());

        phmap::flat_hash_set<mockturtle::node<GateLyt>> collection = root_collection;
        collect_nodes(current_tile, collection);

        for (const auto& next_node : collection)
        {
            if (visited.count(next_node) > 0)
            {
                continue;
            }

            visited.insert(next_node);
            path.push_back(next_node);

            build_sub_circuits(root_collection, path, visited, depth + 1, max_depth, result);

            path.pop_back();
            visited.erase(next_node);
        }
    }
    void build_sub_circuits_from_root(const mockturtle::node<GateLyt>& root, const size_t max_depth,
                                      sub_circuit_map& result) const
    {
        std::vector<mockturtle::node<GateLyt>>          path    = {root};
        phmap::flat_hash_set<mockturtle::node<GateLyt>> visited = {root};

        // if (max_depth == 0)
        // {
        //     result.insert({path, sidb_bdl_sub_circuit<CellLyt, GateLyt>{std::cref(*circuit), path}});
        //     return;
        // }

        phmap::flat_hash_set<mockturtle::node<GateLyt>> root_collection;
        collect_nodes(gate_layout.get_tile(root), root_collection);

        build_sub_circuits(root_collection, path, visited, 0, max_depth, result);

        std::cout << "RESULT: (root = " << gate_layout.get_tile(root) << ')' << std::endl;
        for (const auto& [v, _] : result)
        {
            for (const auto nn : v)
            {
                std::cout << gate_layout.get_tile(nn) << '\t';
            }
            std::cout << std::endl;
        }
        std::cout << "END" << std::endl;
    }

    [[nodiscard]] uint64_t collect_indices_to_trial(const mockturtle::node<GateLyt>&              n,
                                                    const std::vector<mockturtle::node<GateLyt>>& node_vec,
                                                    foreach_node<std::vector<uint64_t>>& sampled_indices) const noexcept
    {
        const size_t                       num_nodes = node_vec.size() - 1;  // skip n
        std::vector<std::vector<uint64_t>> trialable_indices_per_node;
        std::vector<uint64_t>              bases;  // base for each position (number of valid indices per node)

        uint64_t num_possible_trials = 1;
        for (const mockturtle::node<GateLyt>& node : node_vec)
        {
            if (node == n)
            {
                continue;
            }

            std::vector<uint64_t> trialable_indices;

            for (uint64_t i = 0; i < circuit->gate_designs.at(node).size(); ++i)
            {
                trialable_indices.emplace_back(i);
            }

            assert(!trialable_indices.empty() && "No valid indices to trial for node.");

            bases.push_back(trialable_indices.size());
            num_possible_trials *= trialable_indices.size();
            trialable_indices_per_node.emplace_back(std::move(trialable_indices));
        }

        // Step 1: Sample unique indices in range [0, num_possible_trials)
        phmap::flat_hash_set<uint64_t>          sampled_flat_indices;
        std::mt19937                            rng(std::random_device{}());
        std::uniform_int_distribution<uint64_t> dist(0, num_possible_trials - 1);

        while (sampled_flat_indices.size() < num_possible_trials)
        {
            sampled_flat_indices.insert(dist(rng));
        }

        // Step 2: For each flat index, compute the mixed-base digit vector
        std::vector<std::vector<uint64_t>> combinations;
        combinations.reserve(num_possible_trials);

        for (uint64_t flat_index : sampled_flat_indices)
        {
            std::vector<uint64_t> combination;
            for (size_t i = 0; i < num_nodes; ++i)
            {
                const uint64_t base  = bases[i];
                const uint64_t digit = flat_index % base;
                flat_index /= base;
                combination.push_back(trialable_indices_per_node[i][digit]);
            }
            combinations.emplace_back(std::move(combination));
        }

        // Step 3: Populate sampled_indices per node
        for (uint64_t node_idx = 0, other_node_idx = 0; node_idx < node_vec.size(); ++node_idx)
        {
            if (node_vec.at(node_idx) == n)
            {
                continue;
            }

            std::vector<uint64_t> node_trials;
            node_trials.reserve(num_possible_trials);

            for (const auto& combo : combinations)
            {
                node_trials.emplace_back(combo[other_node_idx]);
            }

            ++other_node_idx;

            sampled_indices[node_vec[node_idx]] = std::move(node_trials);
        }

        return num_possible_trials;
    }

    [[nodiscard]] bool powerset_check(
        const phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash>& input) const noexcept
    {
        const uint64_t combination_size = input.size();

        std::vector<decltype(input.cbegin())> iterators;
        iterators.reserve(combination_size);

        for (auto it = input.cbegin(); it != input.cend(); ++it)
        {
            iterators.push_back(it);
        }

        // Iterate over all non-empty subsets using bitmask
        for (uint64_t mask = 1; mask < 1ULL << combination_size; ++mask)
        {
            if ((mask & (mask - 1)) == 0)
            {
                continue;  // combination must be size > 1
            }

            phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash> subset;

            for (size_t i = 0; i < combination_size; ++i)
            {
                if (mask & (1ULL << i))
                {
                    subset.insert(*iterators[i]);
                }
            }

            if (excluded_combinations.find(subset) != excluded_combinations.cend())
            {
                return false;
            }
        }

        return true;
    }

    sidb_technology::cell_type cell_type_of_canvas(const mockturtle::node<GateLyt>& node) const noexcept
    {
        if (gate_layout.is_pi(node))
        {
            return sidb_technology::cell_type::INPUT;
        }

        if (gate_layout.is_po(node))
        {
            if (sidb_counts.at(node) == std::make_pair<uint8_t, uint8_t>(2, 2))
            {
                return sidb_technology::cell_type::OUTPUT;
            }

            if (sidb_counts.at(node) == std::make_pair<uint8_t, uint8_t>(1, 1))
            {
                return sidb_technology::cell_type::OUTPUT_PERTURBER;
            }

            assert(false);
        }

        if (gate_layout.is_buf(node))
        {
            return sidb_technology::cell_type::NORMAL;
        }

        return sidb_technology::cell_type::LOGIC;
    }

    [[nodiscard]] bool perform_trial(const sidb_bdl_sub_circuit<CellLyt, GateLyt>& sub_circuit,
                                     const mockturtle::node<GateLyt>& n, const gate_design_t& gate_design,
                                     const std::vector<mockturtle::node<GateLyt>>& node_vec,
                                     const foreach_node<std::vector<uint64_t>>&    indices_to_trial,
                                     const uint64_t trial_number, const CellLyt& cell_lyt,
                                     uint64_t& total_number_of_simulator_calls) noexcept
    {
        phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash> current_combination{};

        if (node_vec.size() > 1)
        {
            current_combination.emplace(n, gate_design);

            for (const mockturtle::node<GateLyt>& node : node_vec)
            {
                if (node != n)
                {
                    current_combination.emplace(
                        node, circuit->gate_designs.at(node).at(indices_to_trial.at(node).at(trial_number)));
                }
            }

            const std::lock_guard lock{excluded_combinations_mutex};

            if (!powerset_check(current_combination))
            {
                return false;
            }
        }

        CellLyt cell_lyt_clone = cell_lyt.clone();

        // assign random gate design to other gates
        for (const mockturtle::node<GateLyt>& node : node_vec)
        {
            if (node == n)
            {
                continue;
            }

            assert(trial_number < indices_to_trial.at(node).size() &&
                   "Trial number exceeds the number of available samples");

            // select a random gate implementation
            for (const uint64_t gate_design_cell_index :
                 circuit->gate_designs.at(node).at(indices_to_trial.at(node).at(trial_number)))
            {
                cell_lyt_clone.assign_cell_type(circuit->canvasses.at(node).positions.at(gate_design_cell_index),
                                                cell_type_of_canvas(node));
            }
        }

        // if (op_assessment.simulator_invocations > 0)
        // {
        ++total_number_of_simulator_calls;
        // }

        if (is_circuit_operational<CellLyt, GateLyt, local_external_potential_type::BOUNDED>(
                sidb_cell_level_bdl_circuit<CellLyt, GateLyt>{cell_lyt_clone, sub_circuit}, operational_params) ==
            operational_status::OPERATIONAL)
        {
            return true;
        }

        if (node_vec.size() > 1)
        {
            const std::lock_guard lock{excluded_combinations_mutex};

            excluded_combinations.insert(std::move(current_combination));
        }

        return false;
    }

    struct gate_fitness_assessment
    {
        uint64_t gate_index;
        double   fitness;
        bool     selected{false};

        gate_fitness_assessment() = default;

        gate_fitness_assessment(const uint64_t gate_index_, const double fitness_) :
                gate_index(gate_index_),
                fitness(fitness_)
        {}
    };

    void make_pruning_assessments(const sub_circuit_ptr_map& sub_circuits, const mockturtle::node<GateLyt>& n,
                                  const tile<GateLyt>& t, std::vector<bool>& gate_design_assessments) noexcept
    {
        const uint64_t num_gate_designs = circuit->gate_designs.at(n).size();
        const uint64_t num_threads      = params.available_threads;
        // std::min(params.available_threads, num_gate_designs ? num_gate_designs : uint64_t{1});

        // 2) Precompute, for each (gate_design j, sub s), the sampled_indices and number of trials,
        //    and prepare a CellLyt configured for that gate_design & sub-circuit (cloned skeleton + assigned logic
        //    cells).
        struct PairPrep
        {
            // keep the node_vec (copy) and sampled_indices for that pair
            std::vector<mockturtle::node<GateLyt>> node_vec;
            foreach_node<std::vector<uint64_t>>    sampled_indices;
            uint64_t                               actual_num_trials{0};
            CellLyt cell_lyt{};  // cell layout adjusted for that gate_design on this sub-circuit
            const sidb_bdl_sub_circuit<CellLyt, GateLyt>* sub_ptr{nullptr};
        };

        // We'll store in a flat vector indexed by: idx = j * num_subs + s_idx
        const uint64_t num_subs  = sub_circuits.size();
        const size_t   num_pairs = static_cast<size_t>(num_gate_designs) * num_subs;

        std::vector<PairPrep> prep(num_pairs);

        std::atomic<uint64_t> next_gate{0};
        std::atomic<uint64_t> total_number_of_trials_all{0};

        std::vector<std::thread> prep_threads;
        prep_threads.reserve(num_threads);

        for (uint64_t thread_id = 0; thread_id < num_threads; ++thread_id)
        {
            prep_threads.emplace_back(
                [&]
                {
                    uint64_t local_trials = 0;

                    while (true)
                    {
                        const uint64_t j = next_gate.fetch_add(1, std::memory_order_relaxed);
                        if (j >= num_gate_designs)
                            break;

                        const gate_design_t& gate_design = circuit->gate_designs.at(n).at(j);

                        uint32_t s_idx = 0;
                        for (const auto& [node_vec, sub_circuit_ptr] : sub_circuits)
                        {
                            const size_t idx = j * num_subs + s_idx;

                            PairPrep pp;
                            pp.node_vec = node_vec;
                            pp.sub_ptr  = sub_circuit_ptr;

                            for (const uint64_t gate_design_cell_index : gate_design)
                            {
                                pp.cell_lyt.assign_cell_type(
                                    circuit->canvasses.at(n).positions.at(gate_design_cell_index),
                                    cell_type_of_canvas(n));
                            }

                            pp.actual_num_trials = collect_indices_to_trial(n, pp.node_vec, pp.sampled_indices);

                            local_trials += pp.actual_num_trials;

                            prep[idx] = std::move(pp);

                            ++s_idx;
                        }
            }

            total_number_of_trials_all.fetch_add(local_trials, std::memory_order_relaxed);
        });
}

for (auto& th : prep_threads)
{
    if (th.joinable())
        th.join();
}

        // 3) Build a job list where each job is one trial.
        struct TrialJob
        {
            uint32_t gate_design_index;  // j
            uint32_t sub_index;          // s
            uint64_t trial_number;       // 0..actual_num_trials-1 (for that pair)
        };

        std::vector<TrialJob> jobs;
        jobs.reserve(static_cast<size_t>(total_number_of_trials_all));

        // Also create a compact array that stores actual_num_trials per pair for easy lookup later
        // and also track which pairs had zero trials (they automatically make gate design fail).
        std::vector<uint64_t> actual_num_trials_flat;
        actual_num_trials_flat.reserve(prep.size());
        for (uint64_t j = 0; j < num_gate_designs; ++j)
        {
            for (uint32_t s = 0; s < static_cast<uint32_t>(sub_circuits.size()); ++s)
            {
                size_t   idx    = j * sub_circuits.size() + s;
                uint64_t actual = prep[idx].actual_num_trials;
                actual_num_trials_flat.push_back(actual);

                for (uint64_t tr = 0; tr < actual; ++tr)
                {
                    jobs.push_back({static_cast<uint32_t>(j), s, tr});
                }
            }
        }

        // 4) Array of per-(j,s) atomic flags marking whether that pair has found an operational combination.
        //    Flat storage: flag_idx = j * num_subs + s
        // Use dynamically allocated atomic array (std::atomic is not copyable, but array new works).
        std::unique_ptr<std::atomic<uint8_t>[]> pair_found_operational(new std::atomic<uint8_t>[num_pairs]);
        for (size_t i = 0; i < num_pairs; ++i)
        {
            pair_found_operational[i].store(0, std::memory_order_relaxed);
        }

        // 5) Set up worker threads consuming jobs (one job == one trial).
        std::vector<std::thread> threads;
        threads.reserve(num_threads);

        std::atomic<uint64_t> job_index{0};
        std::atomic<uint64_t> processed_jobs{0};
        std::atomic<uint64_t> global_simulator_calls{0};

#if (PROGRESS_BARS)
        mockturtle::progress_bar bar{static_cast<uint32_t>(std::min<uint64_t>(jobs.size(), UINT32_MAX)),
                                     "[i] Finding gate designs to prune for tile " + fmt::format("({},{})", t.x, t.y) +
                                         ": |{0}|"};
#endif

        // std::mutex cout_mutex;

        for (uint64_t thread_id = 0; thread_id < num_threads; ++thread_id)
        {
            threads.emplace_back(
                [&, thread_id]
                {
                    uint64_t local_simulator_calls = 0;

                    while (true)
                    {
                        const uint64_t my_job = job_index.fetch_add(1, std::memory_order_relaxed);
                        if (my_job >= jobs.size())
                            break;

                        const TrialJob& job = jobs[my_job];
                        const size_t    pair_idx =
                            static_cast<size_t>(job.gate_design_index) * sub_circuits.size() + job.sub_index;

                        // If another thread already found an operational combo for this (j,s) pair, skip.
                        if (pair_found_operational[pair_idx].load(std::memory_order_acquire))
                        {
                            processed_jobs.fetch_add(1, std::memory_order_relaxed);
#if (PROGRESS_BARS)
                            if (thread_id == 0)
                                bar(processed_jobs.load());
#endif
                            continue;
                        }

                        // Grab necessary precomputed items for this pair:
                        const PairPrep& pp = prep[pair_idx];

                        // If pair has zero trials skip — but we wouldn't have jobs for zero-trial pairs by
                        // construction. Call perform_trial with the prepared arguments.
                        const gate_design_t& gate_design = circuit->gate_designs.at(n).at(job.gate_design_index);

                        const bool operational =
                            perform_trial(*pp.sub_ptr, n, gate_design, pp.node_vec, pp.sampled_indices,
                                          job.trial_number, pp.cell_lyt, local_simulator_calls);

                        if (operational)
                        {
                            // mark that this (j,s) has an operational combination
                            pair_found_operational[pair_idx].store(1, std::memory_order_release);
                        }

                        processed_jobs.fetch_add(1, std::memory_order_relaxed);
#if (PROGRESS_BARS)
                        if (thread_id == 0)
                            bar(processed_jobs.load());
#endif
                    }

                    // aggregate this thread's simulator calls to global counter
                    global_simulator_calls.fetch_add(local_simulator_calls, std::memory_order_relaxed);

                    // const std::lock_guard lock{cout_mutex};
                    // std::cout << fmt::format("thread {} finished :: processed jobs: {} | local sim calls: {}",
                    // thread_id,
                    //                          job_index.load(), local_simulator_calls)
                    //           << std::endl;
                });
        }

        // join
        for (auto& th : threads)
        {
            if (th.joinable())
                th.join();
        }

        // 6) Reduce results to gate_design_assessments[j]
        for (uint64_t j = 0; j < num_gate_designs; ++j)
        {
            bool exists_sub_circuit_with_no_operational_combination = false;

            for (uint64_t s = 0; s < sub_circuits.size(); ++s)
            {
                // If for this pair nobody found an operational trial, then prune
                if (pair_found_operational[j * sub_circuits.size() + s].load(std::memory_order_acquire) == 0)
                {
                    exists_sub_circuit_with_no_operational_combination = true;
                    break;
                }
            }

            gate_design_assessments[j] = !exists_sub_circuit_with_no_operational_combination;
        }

        std::cout << fmt::format("\n{} combinations | {} simulations | cache hits ~ {:.4f}%",
                                 total_number_of_trials_all.load(), global_simulator_calls.load(),
                                 total_number_of_trials_all.load() == 0 ?
                                     0.0 :
                                     (100.0 * (1.0 - static_cast<double>(global_simulator_calls.load()) /
                                                         static_cast<double>(total_number_of_trials_all.load()))))
                  << std::endl;
    }
    /**
     * todo
     */
    bool prune_gate_designs(const uint64_t level, const std::string_view& level_str)
    {
        std::cout << "\n\nPRUNING GATE DESIGNS\tLEVEL: " << level_str << std::endl;

        sub_circuit_map                   all_sub_circuits{};
        foreach_node<sub_circuit_ptr_map> sub_circuits{};

        foreach_node<std::vector<bool>> gate_design_assessments{};

        foreach_node<CellLyt> sub_circuit_skeletons{};

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (gate_layout.is_constant(n))
                {
                    return;
                }

                build_sub_circuits_from_root(n, level, all_sub_circuits);
            });

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (gate_layout.is_constant(n))
                {
                    return;
                }

                sub_circuits[n];

                for (auto&& [node_vec, sub_circuit] : all_sub_circuits)
                {
                    if (std::find(node_vec.cbegin(), node_vec.cend(), n) != node_vec.end())
                    {
                        sub_circuits[n].insert({node_vec, &sub_circuit});
                    }
                }
            });

        std::map<mockturtle::node<GateLyt>, uint64_t> num_gate_designs_pruned{};

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (gate_layout.is_constant(n))
                {
                    return;
                }

                gate_design_assessments[n].resize(circuit->gate_designs.at(n).size());

                const tile<GateLyt>& t = gate_layout.get_tile(n);

                const uint64_t number_of_layouts = sub_circuits.at(n).size();

                std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                std::cout << fmt::format("Checking gate design combinations for {} sub-circuit{} for {} gate designs",
                                         number_of_layouts, number_of_layouts > 1 ? "s" : "",
                                         circuit->gate_designs.at(n).size())
                          << std::endl;

                make_pruning_assessments(sub_circuits.at(n), n, t, gate_design_assessments[n]);

                std::vector<gate_design_t> selected_gate_implementations{};

                for (uint64_t ix = 0; ix < gate_design_assessments.at(n).size(); ++ix)
                {
                    if (gate_design_assessments.at(n).at(ix))
                    {
                        selected_gate_implementations.push_back(std::move(circuit->gate_designs[n][ix]));
                    }
                    else
                    {
                        ++num_gate_designs_pruned[n];
                    }
                }

                circuit->gate_designs[n] = std::move(selected_gate_implementations);

                std::cout << "PRUNED " << num_gate_designs_pruned[n] << " out of "
                          << gate_design_assessments.at(n).size() << std::endl;

                if (circuit->gate_designs.at(n).empty())
                {
                    throw std::runtime_error{fmt::format("All gate designs pruned for tile {}", t)};
                }

                if (num_gate_designs_pruned.at(n) == 0)
                {
                    return;
                }

                if (circuit->gate_designs.at(n).size() == 1)
                {
                    CellLyt lyt{};
                    for (const uint64_t gate_design_cell_index : circuit->gate_designs.at(n).front())
                    {
                        lyt.assign_cell_type(circuit->canvasses.at(n).positions.at(gate_design_cell_index),
                                                     cell_type_of_canvas(n));
                    }
                    print_layout(lyt);
                }

                circuit->tighten_gate_influence_bounds_until_fixpoint();
            });

        std::stringstream ss{};

        uint64_t pruned_total    = 0;
        uint64_t remaining_total = 0;

        ss << "\n\n**********************************\n";
        ss << std::left << std::setw(10) << "TILE"
           << " | " << std::right << std::setw(8) << "#PRUNED"
           << " | " << std::right << std::setw(10) << "#REMAINING"
           << "\n";
        ss << "----------------------------------\n";

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (gate_layout.is_constant(n))
                {
                    return;
                }

                const auto num_pruned = num_gate_designs_pruned.at(n);

                const auto& tile      = gate_layout.get_tile(n);
                const auto  remaining = circuit->gate_designs.at(n).size();

                ss << std::left << std::setw(10) << tile << " | " << std::right << std::setw(8) << num_pruned << " | "
                   << std::right << std::setw(10) << remaining << "\n";

                pruned_total += num_pruned;
                remaining_total += remaining;
            });
        ss << "----------------------+----------- +\n";
        ss << std::right << std::setw(21) << pruned_total << " | " << std::right << std::setw(10) << remaining_total
           << "\n";
        ss << "**********************************\n";

        if (pruned_total == 0)
        {
            return true;
        }

        std::cout << ss.str();

        return false;

        return pruned_total == 0;
    }
    /**
     * todo
     */
    std::vector<CellLyt> exhaustively_enumerate_gate_design_combinations() const noexcept
    {
        std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;

        std::vector<CellLyt> result{};

        // ---- Precompute gate list and radix sizes ----
        std::vector<std::pair<mockturtle::node<GateLyt>, const std::vector<gate_design_t>*>> gates;
        std::vector<uint64_t>                                                                radix;

        gates.reserve(circuit->gate_designs.size());
        radix.reserve(circuit->gate_designs.size());

        uint64_t total_num_combinations = 1;

        for (const auto& [n, designs] : circuit->gate_designs)
        {
            gates.emplace_back(n, &designs);
            radix.push_back(designs.size());
            total_num_combinations *= designs.size();
        }

        std::cout << fmt::format("{} combinations to explore", total_num_combinations) << std::endl;

        const uint64_t max_threads = std::min<uint64_t>(params.available_threads, total_num_combinations);

        std::atomic<uint64_t> processed_jobs{0};
        std::atomic<uint64_t> next_job{0};

#if (PROGRESS_BARS)
        mockturtle::progress_bar bar{static_cast<uint32_t>(std::min<uint64_t>(total_num_combinations, UINT32_MAX)),
                                     "[i] Finding operational circuits: |{0}|"};
#endif

        std::mutex               result_mutex;
        std::vector<std::thread> workers;

        uint64_t next_printed_result = 0;

        const auto tcm = std::make_unique<thread_count_manager>(params.available_threads - max_threads);

        auto worker = [&](const uint64_t thread_id)
        {
            std::vector<uint64_t> indices(radix.size());

            while (true)
            {
                const uint64_t job = next_job.fetch_add(1, std::memory_order_relaxed);
                if (job >= total_num_combinations)
                    break;

                uint64_t tmp = job;

                // decode mixed radix index
                for (uint64_t i = 0; i < radix.size(); ++i)
                {
                    indices[i] = tmp % radix[i];
                    tmp /= radix[i];
                }

                phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash> current_combination;

                for (uint64_t i = 0; i < gates.size(); ++i)
                {
                    const auto& [n, designs_ptr] = gates[i];
                    current_combination.emplace(n, designs_ptr->at(indices[i]));
                }

                if (!powerset_check(current_combination))
                {
                    processed_jobs.fetch_add(1, std::memory_order_relaxed);
                    continue;
                }

                CellLyt operational_circuit_candidate{};

                for (uint64_t i = 0; i < gates.size(); ++i)
                {
                    const auto& [n, designs_ptr] = gates[i];
                    const auto& design           = designs_ptr->at(indices[i]);

                    for (const uint64_t gate_design_cell_index : design)
                    {
                        operational_circuit_candidate.assign_cell_type(
                            circuit->canvasses.at(n).positions.at(gate_design_cell_index), cell_type_of_canvas(n));
                    }
                }

                if (is_circuit_operational(
                        sidb_cell_level_bdl_circuit<CellLyt, GateLyt>{operational_circuit_candidate,
                                                                      sidb_bdl_sub_circuit<CellLyt, GateLyt>{*circuit}},
                        operational_params, tcm) == operational_status::OPERATIONAL)
                {
                    std::lock_guard lock{result_mutex};
                    std::cout << " tryghfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff "
                                 "\n\n\n\n\n\n\\n";
                    print_layout(operational_circuit_candidate);
                    result.push_back(operational_circuit_candidate);
                }

                processed_jobs.fetch_add(1, std::memory_order_relaxed);

#if (PROGRESS_BARS)
                if (thread_id == 0)
                {
                    bar(processed_jobs.load());

                    if (params.print_found_circuits)
                    {
                        uint64_t size;

                        {
                            std::lock_guard lock{result_mutex};
                            size = result.size();
                        }

                        for (uint64_t ix = next_printed_result; ix < size; ++ix)
                        {
                            print_layout(result.at(ix));
                        }

                        next_printed_result += size - next_printed_result;
                }
            }
#endif
        }

        tcm->return_threads(1);
    };

    for (uint64_t i = 0; i < max_threads; ++i)
    {
        workers.emplace_back(worker, i);
    }

    for (auto& w : workers)
    {
        w.join();
    }

    std::cout << fmt::format("{} combinations | {} simulations | cache hits ~ {:.4f}%",
                             total_num_combinations,
                             processed_jobs.load(),
                             1.0 - static_cast<double>(processed_jobs.load()) /
                                       static_cast<double>(total_num_combinations))
              << std::endl;

    return result;
    }
};

}  // namespace detail

/**
 *
 *
 * @tparam Ntk The type of the input network.
 * @tparam CellLyt SiDB cell-level layout type.
 * @tparam GateLyt Gate-level layout type.
 * todo
 * @param ntk The input network to be mapped onto the defective surface.
 * @param lattice_tiling The lattice tiling used for the circuit design.
 * @param params The parameters used for designing the circuit, encapsulated in an
 * `advanced_circuit_design_params` object.
 * @param stats Pointer to a structure for collecting statistics. If nullptr, statistics are not collected.
 * @return A `sidb_defect_surface<CellLyt>` representing the designed circuit on the defective surface.
 */
template <typename CellLyt, typename GateLyt>
[[nodiscard]] std::vector<CellLyt> advanced_circuit_design(
    const GateLyt&                                                                                      gate_lyt,
    const std::unordered_map<mockturtle::node<cart_odd_row_gate_clk_lyt>, std::pair<uint8_t, uint8_t>>& sidb_count_map,
    const advanced_circuit_design_params<CellLyt>& params = {}, advanced_circuit_design_stats<GateLyt>* stats = nullptr)
{
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    // static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt is not a hexagonal");
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");

    advanced_circuit_design_stats<GateLyt> st{};

    const auto perform_design_task_and_write_stats = [&](detail::advanced_circuit_design_impl<CellLyt, GateLyt>&& p)
    {
        const auto result = p.design_circuit();

        if (stats)
        {
            *stats = st;
        }

        return result;
    };

    return perform_design_task_and_write_stats(
        detail::advanced_circuit_design_impl<CellLyt, GateLyt>{gate_lyt, sidb_count_map, params, st});
}

}  // namespace fiction

#endif  // advanced_CIRCUIT_DESIGN_HPP

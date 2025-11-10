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
    /**
     * The gate-level layout after P&R.
     */
    std::optional<GateLyt> gate_layout{};
};

namespace detail
{

template <typename Ntk, typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
class advanced_circuit_design_impl
{
  public:
    advanced_circuit_design_impl(const Ntk& ntk, const advanced_circuit_design_params<CellLyt>& design_params,
                                 advanced_circuit_design_stats<GateLyt>& st) :
            network{ntk},
            params{design_params},
            stats{st}
    {
        operational_params.simulation_parameters = params.sim_params;
        // operational_params.input_bdl_iterator_params =
        //     params.design_gate_params.operational_params.input_bdl_iterator_params;
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;
    }

    advanced_circuit_design_impl(const GateLyt& gate_lyt, const advanced_circuit_design_params<CellLyt>& design_params,
                                 advanced_circuit_design_stats<GateLyt>& st) :
            gate_layout{gate_lyt},
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

        if (network)
        {
            stats.gate_layout = exact<GateLyt>(*network, params.exact_design_parameters, &stats.exact_stats);

            if (!stats.gate_layout.has_value() || bounding_box_2d<GateLyt>{*stats.gate_layout}.get_x_size() > 5 ||
                bounding_box_2d<GateLyt>{*stats.gate_layout}.get_y_size() > 5 ||
                bounding_box_2d<GateLyt>{*stats.gate_layout}.get_x_size() < 1 ||
                bounding_box_2d<GateLyt>{*stats.gate_layout}.get_y_size() < 1)
            {
                // P&R was unsuccessful
                std::cout << "UNSUCCESS" << std::endl;
                return {};
            }
        }
        else if (gate_layout)
        {
            stats.gate_layout = gate_layout;
        }
        else
        {
            std::cout << "NO NETWORK OR GATE LAYOUT GIVEN" << std::endl;
            return {};
        }

        circuit.emplace(*stats.gate_layout, params.sim_params, params.canvas);

        uint64_t number_of_gates_to_design = 0;

        stats.gate_layout->foreach_node(
            [&](const mockturtle::node<GateLyt>& n)
            {
                if (!skip_physical_design_for_node(*stats.gate_layout, n))
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
    std::optional<GateLyt> gate_layout;
    /**
     * Network.
     */
    std::optional<Ntk> network;
    /**
     * Parameters for the on-the-fly circuit design.
     */
    const advanced_circuit_design_params<CellLyt> params{};
    /**
     * Statistics for the on-the-fly circuit design.
     */
    advanced_circuit_design_stats<GateLyt>&                                stats;
    std::optional<sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>> circuit{};

    is_circuit_operational_params operational_params{};

    template <typename T>
    using foreach_node = std::unordered_map<mockturtle::node<GateLyt>, T>;

    using gate_designs_per_node = foreach_node<std::vector<typename SkeletonGateLibrary::fcn_gate>>;

    uint64_t circuit_design_level = 0;

    using gate_design_t = typename sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::canvas_combination;

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
            stats.gate_layout->foreach_node(
                [&, this](const auto& n, [[maybe_unused]] auto i)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = stats.gate_layout->get_tile(n);

                    const std::pair<cell<CellLyt>, cell<CellLyt>> actual_canvas = params.canvas;
                    // t == tile<GateLyt>{1, 2} ? params.canvas :
                    //                            std::make_pair(cell<CellLyt>{11, 10}, cell<CellLyt>{13, 13});

                    const std::vector<cell<CellLyt>>& all_sidbs_in_canvas = all_coordinates_in_spanned_area(
                        circuit->relative_to_absolute_canvas_position(circuit->gate_layout, actual_canvas.first, t),
                        circuit->relative_to_absolute_canvas_position(circuit->gate_layout, actual_canvas.second, t));

                    std::vector<gate_design_t> all_combinations;

                    // // Reserve an estimated capacity if possible (optional optimization) todo
                    // all_combinations.reserve(estimate_total_combinations(max_sidbs, total_positions));

                    for (std::size_t num_sidbs = 0;
                         // num_sidbs <= (t == tile<GateLyt>{1, 2} ? params.maximum_number_of_canvas_sidbs : 2);
                         num_sidbs <= params.maximum_number_of_canvas_sidbs; ++num_sidbs)
                    {
                        auto combinations = determine_all_combinations_of_distributing_k_entities_on_n_positions(
                            num_sidbs, all_sidbs_in_canvas.size());
                        std::move(combinations.begin(), combinations.end(), std::back_inserter(all_combinations));
                    }

                    std::cout << "NUM COMBINATIONS: " << all_combinations.size() << std::endl;

                    std::shuffle(all_combinations.begin(), all_combinations.end(),
                                 std::mt19937(std::random_device()()));

                    circuit->all_canvas_positions[n] = std::move(all_sidbs_in_canvas);
                    circuit->gate_designs[n]         = std::move(all_combinations);
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
        switch (params.sub_circuit_mode)
        {
            case advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::CONNECTED_GATES:
                for (const auto& in_t : stats.gate_layout->incoming_data_flow(t))
                {
                    if (!stats.gate_layout->is_pi(stats.gate_layout->get_node(in_t)))
                    {
                        collected_nodes.emplace(stats.gate_layout->get_node(stats.gate_layout->below(in_t)));
                    }
                }

                for (const auto& out_t : stats.gate_layout->outgoing_data_flow(t))
                {
                    if (!stats.gate_layout->is_po(stats.gate_layout->get_node(out_t)))
                    {
                        collected_nodes.emplace(stats.gate_layout->get_node(stats.gate_layout->below(out_t)));
                    }
                }

                if (stats.gate_layout->is_wire_tile(t))
                {
                    if (const auto above_t = stats.gate_layout->above(t);
                        t != above_t && stats.gate_layout->is_wire_tile(above_t))
                    {
                        collect_nodes(above_t, collected_nodes);
                    }
                }

                break;
            case advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::ADJACENT_GATES:
                stats.gate_layout->foreach_adjacent_coordinate(
                    t,
                    [&](const auto& c)
                    {
                        if (const auto& n = stats.gate_layout->get_node(c);
                            !skip_physical_design_for_node(*stats.gate_layout, n))
                        {
                            collected_nodes.emplace(n);
                        }
                    });
                break;
            default:  // advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::ALL_GATES:
                stats.gate_layout->foreach_node(
                    [&](const auto& n)
                    {
                        if (skip_physical_design_for_node(*stats.gate_layout, n))
                        {
                            return;
                        }

                        if (const tile<GateLyt>& other_t = stats.gate_layout->get_tile(n); other_t != t)
                        {
                            collected_nodes.emplace(n);
                        }
                    });
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

    using gate_lyt_window_map =
        std::unordered_map<std::vector<mockturtle::node<GateLyt>>,
                           sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>, VectorHash>;

    void build_subcircuits(const phmap::flat_hash_set<mockturtle::node<GateLyt>>& root_collection,
                           std::vector<mockturtle::node<GateLyt>>&                path,
                           phmap::flat_hash_set<mockturtle::node<GateLyt>>& visited, const size_t depth,
                           const size_t max_depth, gate_lyt_window_map& result) const
    {
        if (depth == max_depth)
        {
            // Get tiles corresponding to path nodes
            std::vector<tile<GateLyt>> tiles{};
            for (const auto& node : path)
            {
                tiles.push_back(stats.gate_layout->get_tile(node));
            }

            for (const auto& [other_path, _] : result)
            {
                if (std::all_of(other_path.cbegin(), other_path.cend(),
                                [&](const auto& n) { return visited.count(n) > 0; }))
                {
                    return;
                }
            }

            result.insert(
                {path, sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{std::cref(*circuit), tiles}});

            return;
        }

        const tile<GateLyt> current_tile = stats.gate_layout->get_tile(path.back());

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

            build_subcircuits(root_collection, path, visited, depth + 1, max_depth, result);

            path.pop_back();
            visited.erase(next_node);
        }
    }
    void build_subcircuits_from_root(const mockturtle::node<GateLyt> root, const size_t max_depth,
                                     gate_lyt_window_map& result) const
    {
        const tile<GateLyt> current_tile = stats.gate_layout->get_tile(root);

        phmap::flat_hash_set<mockturtle::node<GateLyt>> root_collection;
        collect_nodes(current_tile, root_collection);

        std::vector<mockturtle::node<GateLyt>>          path    = {root};
        phmap::flat_hash_set<mockturtle::node<GateLyt>> visited = {root};

        build_subcircuits(root_collection, path, visited, 0, max_depth, result);
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
    };

    [[nodiscard]] bool perform_trial(const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& sub_circuit,
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
                cell_lyt_clone.assign_cell_type(circuit->all_canvas_positions.at(node).at(gate_design_cell_index),
                                                sidb_technology::cell_type::LOGIC);
            }
        }

        // if (op_assessment.simulator_invocations > 0)
        // {
        ++total_number_of_simulator_calls;
        // }

        if (is_circuit_operational<CellLyt, GateLyt, local_external_potential_type::BOUNDED, SkeletonGateLibrary>(
                sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{cell_lyt_clone, sub_circuit},
                operational_params) == operational_status::OPERATIONAL)
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

    void make_pruning_assessments(const gate_lyt_window_map& gate_lyt_windows, const mockturtle::node<GateLyt>& n,
                                  const tile<GateLyt>& t, std::vector<bool>& gate_design_assessments) noexcept
    {
        const uint64_t num_gate_designs = circuit->gate_designs.at(n).size();
        const uint64_t num_threads      = params.available_threads;
        // std::min(params.available_threads, num_gate_designs ? num_gate_designs : uint64_t{1});

        // 1) Collect only the sub-circuits that include node 'n' into a vector (ordered, indexed).
        struct SubEntry
        {
            std::vector<mockturtle::node<GateLyt>>                             node_vec;
            const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>* sub_ptr;
        };
        std::vector<SubEntry> sub_entries;
        sub_entries.reserve(gate_lyt_windows.size());

        for (const auto& [node_vec, sub_circuit] : gate_lyt_windows)
        {
            if (std::find(node_vec.cbegin(), node_vec.cend(), n) != node_vec.cend())
            {
                sub_entries.push_back({node_vec, &sub_circuit});
            }
        }

        const size_t num_subs = sub_entries.size();

        // 2) Precompute, for each (gate_design j, sub s), the sampled_indices and number of trials,
        //    and prepare a CellLyt configured for that gate_design & sub-circuit (cloned skeleton + assigned logic
        //    cells).
        struct PairPrep
        {
            // keep the node_vec (copy) and sampled_indices for that pair
            std::vector<mockturtle::node<GateLyt>> node_vec;
            foreach_node<std::vector<uint64_t>>    sampled_indices;
            uint64_t                               actual_num_trials{0};
            CellLyt cell_lyt;  // cell layout adjusted for that gate_design on this sub-circuit
            const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>* sub_ptr{nullptr};
        };

        // We'll store in a flat vector indexed by: idx = j * num_subs + s_idx
        std::vector<PairPrep> prep;
        prep.reserve(static_cast<size_t>(num_gate_designs) * num_subs);

        uint64_t total_number_of_trials_all = 0;

        for (uint64_t j = 0; j < num_gate_designs; ++j)
        {
            const gate_design_t& gate_design = circuit->gate_designs.at(n).at(j);

            for (size_t s = 0; s < num_subs; ++s)
            {
                PairPrep pp;
                pp.node_vec = sub_entries[s].node_vec;
                pp.sub_ptr  = sub_entries[s].sub_ptr;

                // clone the skeleton for the sub-circuit and assign cells for this gate design
                pp.cell_lyt = pp.sub_ptr->skeleton.clone();

                for (const uint64_t gate_design_cell_index : gate_design)
                {
                    // original code used circuit->all_canvas_positions.at(n).at(gate_design_cell_index)
                    pp.cell_lyt.assign_cell_type(circuit->all_canvas_positions.at(n).at(gate_design_cell_index),
                                                 sidb_technology::cell_type::LOGIC);
                }

                // determine how many trials (and the sampled indices) for this (j, s)
                pp.actual_num_trials = collect_indices_to_trial(n, pp.node_vec, pp.sampled_indices);

                total_number_of_trials_all += pp.actual_num_trials;

                prep.push_back(std::move(pp));
            }
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
            for (uint32_t s = 0; s < static_cast<uint32_t>(num_subs); ++s)
            {
                size_t   idx    = j * num_subs + s;
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
        const auto num_pairs = static_cast<size_t>(num_gate_designs) * num_subs;
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

                        const TrialJob& job   = jobs[my_job];
                        const size_t pair_idx = static_cast<size_t>(job.gate_design_index) * num_subs + job.sub_index;

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

            for (uint64_t s = 0; s < num_subs; ++s)
            {
                // If for this pair nobody found an operational trial, then prune
                if (pair_found_operational[j * num_subs + s].load(std::memory_order_acquire) == 0)
                {
                    exists_sub_circuit_with_no_operational_combination = true;
                    break;
                }
            }

            gate_design_assessments[j] = !exists_sub_circuit_with_no_operational_combination;
        }

        std::cout << fmt::format("\n{} combinations | {} simulations | cache hits ~ {:.4f}%",
                                 total_number_of_trials_all, global_simulator_calls.load(),
                                 total_number_of_trials_all == 0 ?
                                     0.0 :
                                     (100.0 * (1.0 - static_cast<double>(global_simulator_calls.load()) /
                                                         static_cast<double>(total_number_of_trials_all))))
                  << std::endl;
    }
    /**
     * todo
     */
    bool prune_gate_designs(const uint64_t level, const std::string_view& level_str)
    {
        std::cout << "\n\nPRUNING GATE DESIGNS\tLEVEL: " << level_str << std::endl;

        gate_lyt_window_map gate_lyt_windows{};

        foreach_node<std::vector<bool>> gate_design_assessments{};

        foreach_node<CellLyt> sub_circuit_skeletons{};

        stats.gate_layout->foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(*stats.gate_layout, n))
                {
                    return;
                }

                build_subcircuits_from_root(n, level, gate_lyt_windows);
            });

        stats.gate_layout->foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(*stats.gate_layout, n))
                {
                    return;
                }

                gate_design_assessments[n].resize(circuit->gate_designs.at(n).size());
            });

        std::map<mockturtle::node<GateLyt>, uint64_t> num_gate_designs_pruned{};

        stats.gate_layout->foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(*stats.gate_layout, n))
                {
                    return;
                }

                gate_design_assessments[n].resize(circuit->gate_designs.at(n).size());

                const tile<GateLyt>& t = stats.gate_layout->get_tile(n);

                uint64_t number_of_layouts = 0;

                for (const auto& [node_vec, _] : gate_lyt_windows)
                {
                    if (std::find(node_vec.cbegin(), node_vec.cend(), n) == node_vec.cend())
                    {
                        continue;
                    }

                    number_of_layouts++;
                }

                std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                std::cout << fmt::format("Checking gate design combinations for {} sub-circuit{} for {} gate designs",
                                         number_of_layouts, number_of_layouts > 1 ? "s" : "",
                                         circuit->gate_designs.at(n).size())
                          << std::endl;

                make_pruning_assessments(gate_lyt_windows, n, t, gate_design_assessments[n]);

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

                if (circuit->gate_designs[n].empty())
                {
                    throw std::runtime_error{fmt::format("All gate designs pruned for tile {}", t)};
                }

                if (num_gate_designs_pruned.at(n) == 0)
                {
                    return;
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

        circuit->gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(circuit->gate_layout, n))
                {
                    return;
                }

                const auto num_pruned = num_gate_designs_pruned.at(n);

                const auto& tile      = stats.gate_layout->get_tile(n);
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

        std::vector<CellLyt>               result{};
        std::vector<std::vector<uint64_t>> all_combinations;

        // === 1. Generate all combinations first ===
        std::vector<uint64_t> indices(circuit->gate_designs.size(), 0);  // todo reserve?
        bool                  stop = false;
        while (!stop)
        {
            phmap::flat_hash_set<std::pair<mockturtle::node<GateLyt>, gate_design_t>, SetHash> current_combination{};

            for (uint64_t i = 0; i < indices.size(); ++i)
            {
                const auto& [n, designs] = *std::next(circuit->gate_designs.cbegin(), static_cast<int64_t>(i));
                current_combination.emplace(n, designs.at(indices.at(i)));
            }

            if (powerset_check(current_combination))
            {
                all_combinations.push_back(indices);
            }

            // Increment indices like an odometer f
            for (uint64_t i = 0; i < indices.size(); ++i)
            {
                if (++indices[i] < std::next(circuit->gate_designs.cbegin(), static_cast<int64_t>(i))->second.size())
                {
                    break;  // No carry needed
                }
                indices[i] = 0;  // Reset this index and carry over to the next
                if (i == indices.size() - 1)
                {
                    stop = true;  // Stop when the last index overflows
                }
            }
        }

        uint64_t total_num_combinations = 1;

        for (const auto [_, designs] : circuit->gate_designs)
        {
            total_num_combinations *= designs.size();
        }

        std::cout << fmt::format("{} combinations | {} simulations | cache hits ~ {:.4f}%", total_num_combinations,
                                 all_combinations.size(),
                                 1.0 - static_cast<double>(all_combinations.size()) /
                                           static_cast<double>(total_num_combinations))
                  << std::endl;

        const auto            total_jobs  = all_combinations.size();
        const uint64_t        max_threads = std::min<uint64_t>(params.available_threads, total_jobs);
        std::atomic<uint64_t> processed_jobs{0};

        // === 2. Progress bar setup ===
#if (PROGRESS_BARS)
        mockturtle::progress_bar bar{static_cast<uint32_t>(std::min<uint64_t>(total_jobs, UINT32_MAX)),
                                     "[i] Finding operational circuits: |{0}|"};
#endif

        // === 3. Thread management ===
        std::mutex               result_mutex;
        std::vector<std::thread> workers;
        std::atomic<uint64_t>    next_job{0};

        uint64_t next_printed_result = 0;

        // Thread count manager
        const auto tcm = std::make_unique<thread_count_manager>(params.available_threads - max_threads);

        auto worker = [&](const uint64_t thread_id)
        {
            while (true)
            {
                const uint64_t job_index = next_job.fetch_add(1);
                if (job_index >= total_jobs)
                    break;

                const auto& indices = all_combinations[job_index];

                CellLyt operational_circuit_candidate{circuit->skeleton.clone()};
                for (uint64_t i = 0; i < circuit->gate_designs.size(); i++)
                {
                    const auto& [n, op_gate_designs_for_gate] =
                        *std::next(circuit->gate_designs.cbegin(), static_cast<int64_t>(i));

                    for (const uint64_t gate_design_cell_index : op_gate_designs_for_gate.at(indices.at(i)))
                    {
                        operational_circuit_candidate.assign_cell_type(
                            circuit->all_canvas_positions.at(n).at(gate_design_cell_index),
                            sidb_technology::cell_type::LOGIC);
                    }
                }

                if (is_circuit_operational(
                        sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                            operational_circuit_candidate,
                            sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{*circuit}},
                        operational_params, tcm) == operational_status::OPERATIONAL)
                {
                    std::lock_guard lock{result_mutex};
                    result.push_back(operational_circuit_candidate);
                }

                processed_jobs.fetch_add(1);

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

            // Return one thread's worth of work to the pool
            tcm->return_threads(1);
        };

        // === 4. Spawn threads ===
        for (uint64_t i = 0; i < max_threads; ++i)
        {
            workers.emplace_back(worker, i);
        }

        for (auto& w : workers)
        {
            w.join();
        }

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
template <typename Ntk, typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
[[nodiscard]] std::vector<CellLyt> advanced_circuit_design(const std::variant<Ntk, GateLyt>& ntk_or_gate_lyt,
                                                           const advanced_circuit_design_params<CellLyt>& params = {},
                                                           advanced_circuit_design_stats<GateLyt>* stats = nullptr)
{
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt is not a hexagonal");
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");

    advanced_circuit_design_stats<GateLyt> st{};

    const auto perform_design_task_and_write_stats =
        [&](detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt, SkeletonGateLibrary>&& p)
    {
        const auto result = p.design_circuit();

        if (stats)
        {
            *stats = st;
        }

        return result;
    };

    if (std::holds_alternative<Ntk>(ntk_or_gate_lyt))
    {
        return perform_design_task_and_write_stats(
            detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt, SkeletonGateLibrary>{
                std::get<Ntk>(ntk_or_gate_lyt), params, st});
    }

    return perform_design_task_and_write_stats(
        detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt, SkeletonGateLibrary>{
            std::get<GateLyt>(ntk_or_gate_lyt), params, st});
}

}  // namespace fiction

#endif  // advanced_CIRCUIT_DESIGN_HPP

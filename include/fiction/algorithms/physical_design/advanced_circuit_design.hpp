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
        // operational_params.input_bdl_iterator_params =
        //     params.design_gate_params.operational_params.input_bdl_iterator_params;
    }

    [[nodiscard]] std::optional<CellLyt> design_sidb_layout(const uint64_t num_gates_to_design)
    {
        // initialize
        collect_initial_gate_designs();

        circuit->tighten_gate_influence_bounds_until_fixpoint();

        //
        // operational_params.termination_cond =
        //     is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        while (circuit_design_level < num_gates_to_design - 1)
        {
            // prune by assessing gate design combinations for increasingly large sets of connected gates
            prune_gate_designs(
                circuit_design_level,
                std::string_view{
                    std::to_string(circuit_design_level + 1) +
                    (circuit_design_level == 0 ?
                         " GATE" :
                         (params.sub_circuit_mode ==
                                  advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::CONNECTED_GATES ?
                              " CONNECTED GATES" :
                          params.sub_circuit_mode ==
                                  advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::ADJACENT_GATES ?
                              " ADJACENT GATES" :
                              " GATES"))});

            ++circuit_design_level;
        }

        // prune at the global level (all gates are considered together)
        if (const std::optional<CellLyt>& maybe_lyt = prune_gate_designs(circuit_design_level, "GLOBAL");
            maybe_lyt.has_value())
        {
            return *maybe_lyt;
        }

        std::cout << "exhaustive" << std::endl;

        return std::nullopt;  // exhaustively_enumerate_gate_design_combinations();
    }

    [[nodiscard]] std::optional<CellLyt> design_circuit()
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
                return std::nullopt;
            }
        }
        else if (gate_layout)
        {
            stats.gate_layout = gate_layout;
        }
        else
        {
            std::cout << "NO NETWORK OR GATE LAYOUT GIVEN" << std::endl;
            return std::nullopt;
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

        if (const std::optional<CellLyt>& maybe_sidb_lyt = design_sidb_layout(number_of_gates_to_design);
            maybe_sidb_lyt.has_value())
        {
            std::cout << "\n\nSUCCESS! GENERATED OPERATIONAL CIRCUIT:" << std::endl;
            print_layout(maybe_sidb_lyt.value());

            /// todo not necessary,,,, right?
            // // // add defects to the circuit.
            // // params.sidb_on_the_fly_gate_library_parameters.defect_surface.foreach_sidb_defect(
            // //     [&sidbs_and_defects](const auto& defect)
            // //     { sidbs_and_defects.assign_sidb_defect(defect.first, defect.second); });

            return maybe_sidb_lyt.value();
        }

        std::cout << "\n\nFAILURE: NO OPERATIONAL CIRCUIT COULD BE GENERATED" << std::endl;

        return std::nullopt;
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

                    const std::vector<cell<CellLyt>>& all_sidbs_in_canvas = all_coordinates_in_spanned_area(
                        circuit->relative_to_absolute_canvas_position(circuit->gate_layout, params.canvas.first, t),
                        circuit->relative_to_absolute_canvas_position(circuit->gate_layout, params.canvas.second, t));

                    std::vector<typename sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::canvas_combination>
                        all_combinations;

                    // // Reserve an estimated capacity if possible (optional optimization) todo
                    // all_combinations.reserve(estimate_total_combinations(max_sidbs, total_positions));

                    for (std::size_t num_sidbs = 0; num_sidbs <= params.maximum_number_of_canvas_sidbs; ++num_sidbs)
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

    void collect_nodes(const tile<GateLyt>&                           t,
                       std::unordered_set<mockturtle::node<GateLyt>>& collected_nodes) const noexcept
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

    void build_subcircuits(const std::unordered_set<mockturtle::node<GateLyt>>& root_collection,
                           std::vector<mockturtle::node<GateLyt>>&              path,
                           std::unordered_set<mockturtle::node<GateLyt>>& visited, const size_t depth,
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

        std::unordered_set<mockturtle::node<GateLyt>> collection = root_collection;
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

        std::unordered_set<mockturtle::node<GateLyt>> root_collection;
        collect_nodes(current_tile, root_collection);

        std::vector<mockturtle::node<GateLyt>>        path    = {root};
        std::unordered_set<mockturtle::node<GateLyt>> visited = {root};

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
        std::unordered_set<uint64_t>            sampled_flat_indices;
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

    [[nodiscard]] bool perform_trial(const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& sub_circuit,
                                     const mockturtle::node<GateLyt>&                                   n,
                                     const std::vector<mockturtle::node<GateLyt>>&                      node_vec,
                                     foreach_node<std::vector<uint64_t>>&& indices_to_trial,
                                     const uint64_t trial_number, double& logic_match_average_over_inputs,
                                     CellLyt&                                     cell_lyt_clone,
                                     const std::unique_ptr<thread_count_manager>& tcm) const noexcept
    {
        // std::cout << "PERFORMING TRIAL | number = " << trial_number << " | ";
        std::flush(std::cout);

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

        const sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary> c{cell_lyt_clone, sub_circuit};

        // sub-circuit logic match assessment
        const circuit_operational_assessment<CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
            is_circuit_operational<CellLyt, GateLyt, local_external_potential_type::BOUNDED, SkeletonGateLibrary>(
                c, operational_params, tcm);

        logic_match_average_over_inputs +=
            op_assessment.status == operational_status::OPERATIONAL ? 1.0 : 0.0;  // todo not needed

        // std::cout << (op_assessment.status == operational_status::OPERATIONAL ? "" : "NON-") << "OPERATIONAL"
        //           << std::endl;

        return op_assessment.status == operational_status::OPERATIONAL;

        // assert(op_assessment.assessment_per_input.has_value() && "ALL_INPUT_COMBINATIONS_ENUMERATED is not set.");
        //
        // uint64_t successful_input_combinations  = 0;
        // double   logic_match_sum_for_all_inputs = 0.0;
        //
        // // sum logic match over the different input combinations
        // for (const typename circuit_operational_assessment<CellLyt, local_external_potential_type::BOUNDED>::
        //          operational_assessment_for_input& op_assessment_for_input : *op_assessment.assessment_per_input)
        // {
        //     if (op_assessment_for_input.status == operational_status::OPERATIONAL)
        //     {
        //         successful_input_combinations++;
        //     }
        //
        //     logic_match_sum_for_all_inputs += op_assessment_for_input.logic_match;
        // }
        //
        // logic_match_average_over_inputs +=
        //     logic_match_sum_for_all_inputs / static_cast<double>(op_assessment.assessment_per_input->size());
        //
        // return successful_input_combinations;
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

    void make_trial_based_fitness_assessments(const gate_lyt_window_map&       gate_lyt_windows,
                                              const mockturtle::node<GateLyt>& n, const tile<GateLyt>& t,
                                              const uint64_t min_bound, const uint64_t max_bound,
                                              std::vector<gate_fitness_assessment>& gate_fitness_assessments,
                                              const bool global_pruning, const uint64_t num_input_combinations,
                                              std::mutex& lyt_mutex, std::optional<CellLyt>& maybe_lyt) const noexcept
    {
        const uint64_t num_gate_designs = max_bound - min_bound;

        const uint64_t num_threads = std::min(params.available_threads, num_gate_designs);

        const uint64_t chunk_size = (num_gate_designs + num_threads - 1) / num_threads;  // Ceiling division

#if (PROGRESS_BARS)
        mockturtle::progress_bar bar{static_cast<uint32_t>(std::min(chunk_size, num_gate_designs)),
                                     "[i] Determining successful trial ratio for tile " +
                                         fmt::format("({},{})", t.x, t.y) + ": |{0}|"};
#endif

        std::vector<std::thread> threads{};
        threads.reserve(num_threads);

        std::unique_ptr<thread_count_manager> tcm =
            std::make_unique<thread_count_manager>(params.available_threads - num_threads);

        for (uint64_t i = 0; i < num_threads; ++i)
        {
            threads.emplace_back(
                [&gate_fitness_assessments, &maybe_lyt, &lyt_mutex, &tcm, &n, i,
#if (PROGRESS_BARS)
                 &bar,
#endif
                 min_bound, max_bound, chunk_size, this, &gate_lyt_windows, global_pruning]
                {
                    const uint64_t start_index = min_bound + i * chunk_size;
                    const uint64_t end_index   = std::min(start_index + chunk_size, max_bound);

                    for (uint64_t j = start_index; j < end_index; ++j)
                    {
                        uint64_t total_number_of_trials = 0;

                        double successful_trials = 0;

                        foreach_node<std::vector<uint64_t>> sampled_indices{};

                        // for each sub-circuit
                        for (const auto& [node_vec, gate_lyt_window] : gate_lyt_windows)
                        {
                            if (std::find(node_vec.cbegin(), node_vec.cend(), n) == node_vec.cend())
                            {
                                continue;
                            }

                            CellLyt cell_lyt{gate_lyt_window.skeleton.clone()};

                            // select the j-th gate implementation for n
                            for (const uint64_t gate_design_cell_index : circuit->gate_designs.at(n).at(j))
                            {
                                cell_lyt.assign_cell_type(
                                    circuit->all_canvas_positions.at(n).at(gate_design_cell_index),
                                    sidb_technology::cell_type::LOGIC);
                            }

                            const uint64_t actual_num_trials =
                                collect_indices_to_trial(n, node_vec, sampled_indices);

                            total_number_of_trials += actual_num_trials;

                            uint64_t current_trial = 0;

                            while (current_trial < actual_num_trials)
                            {
                                if (global_pruning)
                                {
                                    std::lock_guard lock{lyt_mutex};

                                    if (maybe_lyt.has_value())
                                    {
                                        return;  // quit if an operational circuit was already found
                                    }
                                }

                                CellLyt cell_lyt_clone = cell_lyt.clone();

                                const bool operational =
                                    perform_trial(gate_lyt_window, n, node_vec, std::move(sampled_indices),
                                                  current_trial, successful_trials, cell_lyt_clone, tcm);

                                current_trial++;

                                if (!operational)
                                {
                                    continue;
                                }

                                if (global_pruning)
                                {
                                    // all input combinations are operational: operational circuit found

                                    std::lock_guard lock{lyt_mutex};

                                    if (!maybe_lyt.has_value())
                                    {
                                        maybe_lyt.emplace(cell_lyt_clone);
                                    }

                                    return;
                                }

                                // return right away -- cannot prune

                                break;
                            }

                            if (successful_trials == 0)
                            {
                                break;
                            }
                        }

#if (PROGRESS_BARS)
                        if (i == 0)
                        {
                            bar(j - min_bound);
                        }
#endif

                        const double successful_trial_ratio =
                            successful_trials / static_cast<double>(total_number_of_trials);

                        gate_fitness_assessments[j].gate_index = j;
                        gate_fitness_assessments[j].fitness    = successful_trial_ratio;
                        gate_fitness_assessments[j].selected   = false;
                    }

                    tcm->return_threads(1);
                });
        }

        for (auto& thread : threads)
        {
            if (thread.joinable())
            {
                thread.join();
            }
        }
    }
    /**
     * todo
     */
    std::optional<CellLyt> prune_gate_designs(const uint64_t level, const std::string_view& level_str)
    {
        std::cout << "\n\nSTARTING TO PRUNE GATE DESIGNS\tLEVEL: " << level_str << std::endl;

        const bool global_pruning = level_str == "GLOBAL";

        std::optional<CellLyt> maybe_lyt{};

        const uint64_t num_input_combinations = 1 << stats.gate_layout->num_pis();

        gate_lyt_window_map gate_lyt_windows{};

        foreach_node<std::vector<gate_fitness_assessment>> gate_fitness_assessments{};

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

        std::mutex lyt_mutex{};

        bool fixpoint = false;

        while (!fixpoint)
        {
            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

                    gate_fitness_assessments[n].resize(circuit->gate_designs.at(n).size());
                });

            std::unordered_map<mockturtle::node<GateLyt>, uint64_t> num_gate_designs_pruned{};

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

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
                    std::cout << fmt::format("Trialing for {} sub-circuit{} for {} gate implementations",
                                             number_of_layouts, number_of_layouts > 1 ? "s" : "",
                                             circuit->gate_designs.at(n).size())
                              << std::endl;

                    uint64_t min_bound = 0;                                   // inclusive
                    uint64_t max_bound = circuit->gate_designs.at(n).size();  // exclusive
                                                                              //
                                                                              // uint64_t repeated_attempt_number = 1;
                                                                              // uint64_t attempt_number          = 1;
                                                                              //
                                                                              // while (true)
                                                                              // {
                    make_trial_based_fitness_assessments(gate_lyt_windows, n, t, min_bound, max_bound,
                                                         gate_fitness_assessments[n], global_pruning,
                                                         num_input_combinations, lyt_mutex, maybe_lyt);

                    if (global_pruning)
                    {
                        if (maybe_lyt.has_value())
                        {
                            return;  // quit if an operational circuit has been found
                        }
                    }

                    for (uint64_t gate_index = 0; gate_index < gate_fitness_assessments.at(n).size(); ++gate_index)
                    {
                        if (gate_fitness_assessments.at(n).at(gate_index).fitness != 0.0)
                        {
                            gate_fitness_assessments[n][gate_index].selected = true;
                        }
                        else
                        {
                            ++num_gate_designs_pruned[n];
                        }
                    }

                    std::cout << "PRUNED " << num_gate_designs_pruned[n] << " out of "
                              << gate_fitness_assessments.at(n).size() << std::endl;

                    //
                    //     if (discriminate_fitness_assessments(
                    //             global_pruning, selectivity, quantization_step, success_rate_ceiling,
                    //             circuit->gate_designs[n], gate_fitness_assessments[n], min_bound, max_bound,
                    //             repeated_attempt_number, attempt_number, completed_assessment[n]))
                    //     {
                    //         break;
                    //     }
                    // }
                });

            std::cout << std::endl;

            if (global_pruning)
            {
                if (maybe_lyt.has_value())
                {
                    return maybe_lyt.value();  // return the operational circuit
                }
            }

            uint64_t pruned_total = 0;
            std::cout << "\n=================\n" << "  TILE  | #PRUNED" << std::endl;
            for (const auto& [n, num_pruned] : num_gate_designs_pruned)
            {
                std::cout << stats.gate_layout->get_tile(n) << " |    " << num_pruned << std::endl;

                pruned_total += num_pruned;
            }
            std::cout << "----------------- +\n" << "             " << pruned_total << std::endl;

            if (pruned_total == 0)
            {
                fixpoint = true;

                continue;
            }

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n) || gate_fitness_assessments.at(n).empty())
                    {
                        return;
                    }

                    std::vector<typename sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::canvas_combination>
                        selected_gate_implementations{};

                    for (const gate_fitness_assessment& fitness_assessment : gate_fitness_assessments.at(n))
                    {
                        if (fitness_assessment.selected)
                        {
                            selected_gate_implementations.push_back(
                                std::move(circuit->gate_designs[n][fitness_assessment.gate_index]));
                        }
                    }

                    circuit->gate_designs[n] = std::move(selected_gate_implementations);
                });

            circuit->tighten_gate_influence_bounds_until_fixpoint();
        }

        return std::nullopt;
    }
    // /**
    //  * todo
    //  */
    // std::optional<CellLyt> exhaustively_enumerate_gate_design_combinations() const noexcept
    // {
    //     std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;
    //
    //     // operational_params.print = true;  // after comment: reenable function constness
    //
    //     std::vector<uint64_t> indices(circuit->gate_designs.size(), 0);
    //
    //     while (true)
    //     {
    //         CellLyt operational_circuit_candidate{};
    //         for (uint64_t i = 0; i < circuit->gate_designs.size(); i++)
    //         {
    //             const auto& [n, op_gate_designs_for_gate] =
    //                 *std::next(circuit->gate_designs.cbegin(), static_cast<int64_t>(i));
    //             // select a random gate implementation for the tile that connects as input to n
    //             assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(operational_circuit_candidate,
    //                                                                op_gate_designs_for_gate.at(indices.at(i)),
    //                                                                *stats.gate_layout,
    //                                                                stats.gate_layout->get_tile(n));
    //         }
    //
    //         std::cout << "trying combination: ";
    //         for (uint64_t i = 0; i < circuit->gate_designs.size(); i++)
    //         {
    //             std::cout << indices.at(i) << " ";
    //         }
    //         std::cout << std::endl;
    //
    //         if (is_circuit_operational(
    //                 sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
    //                     operational_circuit_candidate,
    //                     sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{*circuit}},
    //                 operational_params)
    //                 .status == operational_status::OPERATIONAL)
    //         {
    //             std::cout << "\n\nFINAL GENERATED CIRCUIT:" << std::endl;
    //             print_layout(operational_circuit_candidate);
    //
    //             return operational_circuit_candidate;
    //         }
    //
    //         // Increment indices like an odometer
    //         for (uint64_t i = 0; i < indices.size(); ++i)
    //         {
    //             if (++indices[i] < std::next(circuit->gate_designs.cbegin(), static_cast<int64_t>(i))->second.size())
    //             {
    //                 break;  // No carry needed
    //             }
    //
    //             indices[i] = 0;  // Reset this index and carry over to the next
    //
    //             if (i == indices.size() - 1)
    //             {
    //                 return std::nullopt;  // Stop when the last index overflows
    //             }
    //         }
    //     }
    //
    //     return std::nullopt;
    // }
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
[[nodiscard]] std::optional<CellLyt> advanced_circuit_design(const std::variant<Ntk, GateLyt>& ntk_or_gate_lyt,
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

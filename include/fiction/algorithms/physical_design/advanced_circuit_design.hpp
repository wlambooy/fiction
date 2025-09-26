//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_circuit_operational.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/technology/sidb_on_the_fly_gate_library.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/gate_design_utils.hpp"

#include <fiction/io/write_svg_layout.hpp>

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
    /**
     * Parameters for the *exact* placement and routing algorithm.
     */
    exact_physical_design_params exact_design_parameters = {};
    /**
     * This struct holds parameters to design SiDB gates.
     */
    design_sidb_gates_params<CellLyt> design_gate_params{};
    /** This variable specifies the radius in nanometers around the center of the hexagon where atomic defects are
     * incorporated into the gate design. (unit: nm)
     */

    std::optional<CellLyt> defect_surface{};
    double                 influence_radius_charged_defects = 15;

    uint64_t num_trials                      = 500;
    double   quantization_factor             = 0.075;
    double   selectivity                     = 0.5;
    double   success_rate_ceiling            = 0.95;
    uint64_t maximum_discrimination_attempts = 5;

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
    advanced_circuit_design_impl(const Ntk& ntk, advanced_circuit_design_params<CellLyt>& design_params,
                                 const GateLyt& tiling, advanced_circuit_design_stats<GateLyt>& st) :
            lattice_tiling{tiling},
            network{ntk},
            params{design_params},
            stats{st}
    {
        operational_params.simulation_parameters = params.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params =
            params.design_gate_params.operational_params.input_bdl_iterator_params;
    }

    [[nodiscard]] std::optional<CellLyt> design_sidb_layout(const uint64_t num_gates_to_design)
    {
        circuit.emplace(*stats.gate_layout,
                        params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params);

        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;

        // initialize
        collect_initial_gate_designs();

        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        while (++circuit_design_level < num_gates_to_design - 1)
        {
            // prune by assessing gate design combinations for increasingly large sets of connected gates
            prune_gate_designs(circuit_design_level,
                               std::string_view{std::to_string(circuit_design_level) + " GATE CONNECTION" +
                                                (circuit_design_level > 1 ? "S" : "")});
        }

        // prune at the global level (all gates are considered together)
        if (const std::optional<CellLyt>& maybe_lyt = prune_gate_designs(circuit_design_level, "GLOBAL");
            maybe_lyt.has_value())
        {
            return *maybe_lyt;
        }

        return exhaustively_enumerate_gate_design_combinations();
    }

    [[nodiscard]] std::optional<CellLyt> design_circuit()
    {
        const mockturtle::stopwatch stop{stats.time_total};

        stats.gate_layout = exact<GateLyt>(network, params.exact_design_parameters, &stats.exact_stats);

        if (!stats.gate_layout.has_value() || bounding_box_2d<GateLyt>{*stats.gate_layout}.get_x_size() > 5 ||
            bounding_box_2d<GateLyt>{*stats.gate_layout}.get_y_size() > 5 ||
            bounding_box_2d<GateLyt>{*stats.gate_layout}.get_x_size() < 1 ||
            bounding_box_2d<GateLyt>{*stats.gate_layout}.get_y_size() < 1)
        {
            // P&R was unsuccessful
            std::cout << "UNSUCCESS" << std::endl;
            return std::nullopt;
        }

        print_skeleton_gate_layout(*stats.gate_layout);

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
    GateLyt lattice_tiling;
    /**
     * Network.
     */
    Ntk network;
    /**
     * Parameters for the on-the-fly circuit design.
     */
    advanced_circuit_design_params<CellLyt> params{};
    /**
     * Statistics for the on-the-fly circuit design.
     */
    advanced_circuit_design_stats<GateLyt>&                                stats;
    std::optional<sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>> circuit{};

    is_circuit_operational_params operational_params{};

    template <typename T>
    using foreach_node = std::unordered_map<mockturtle::node<GateLyt>, T>;

    using gate_designs_per_node = foreach_node<std::vector<typename SkeletonGateLibrary::fcn_gate>>;

    gate_designs_per_node gate_designs{};

    uint64_t circuit_design_level = 0;

    void print_skeleton_gate_layout(const GateLyt& gate_layout) const noexcept
    {
        CellLyt lyt = apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout);

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_layout, n))
                {
                    return;
                }

                const auto& t = gate_layout.get_tile(n);

                for (const cell<CellLyt>& relative_c : all_coordinates_in_spanned_area(
                         params.design_gate_params.canvas.first, params.design_gate_params.canvas.second))
                {
                    const cell<CellLyt> absolute_c =
                        relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                           SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                            gate_layout, t, relative_c);
                    lyt.assign_cell_type(absolute_c, sidb_technology::cell_type::LOGIC);
                }
            });

        std::cout << "Skeleton looks like:" << std::endl;
        print_layout(lyt);
        std::cout << std::endl;
    }

    void collect_initial_gate_designs()
    {
        const sidb_on_the_fly_gate_library_params<CellLyt> on_the_fly_params{params.design_gate_params,
                                                                             params.influence_radius_charged_defects};

        try
        {
            stats.gate_layout->foreach_node(
                [&, this](const auto& n, [[maybe_unused]] auto i)
                {
                    if (!skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        gate_designs[n] = sidb_on_the_fly_gate_library<SkeletonGateLibrary::gate_x_size(),
                                                                       SkeletonGateLibrary::gate_y_size()>::
                            template set_up_gates<GateLyt, CellLyt, sidb_on_the_fly_gate_library_params<CellLyt>,
                                                  local_external_potential_type::BOUNDED, SkeletonGateLibrary>(
                                *stats.gate_layout, stats.gate_layout->get_tile(n), on_the_fly_params,
                                params.defect_surface, std::make_optional(*circuit), operational_params);
                    }
                });
        }

        catch (const gate_design_exception<tt, GateLyt>& _)
        {
            throw unsuccessful_gate_design_error("Gate design was unsuccessful");
        }
    }

    static void
    collect_connecting_nodes(const GateLyt& gate_lyt, const tile<GateLyt>& t,
                             std::vector<mockturtle::node<GateLyt>>& connecting_to_t,
                             const std::optional<tile<GateLyt>>&     maybe_do_not_collect = std::nullopt) noexcept
    {
        for (const auto& in_t : gate_lyt.incoming_data_flow(t))
        {
            if (!gate_lyt.is_pi(gate_lyt.get_node(in_t)) &&
                (!maybe_do_not_collect.has_value() || gate_lyt.below(in_t) != *maybe_do_not_collect))
            {
                connecting_to_t.emplace_back(gate_lyt.get_node(gate_lyt.below(in_t)));
            }
        }

        for (const auto& out_t : gate_lyt.outgoing_data_flow(t))
        {
            if (!gate_lyt.is_po(gate_lyt.get_node(out_t)) &&
                (!maybe_do_not_collect.has_value() || gate_lyt.below(out_t) != *maybe_do_not_collect))
            {
                connecting_to_t.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
            }
        }

        if (gate_lyt.is_wire_tile(t))
        {
            if (const auto above_t = gate_lyt.above(t); t != above_t && gate_lyt.is_wire_tile(above_t))
            {
                collect_connecting_nodes(gate_lyt, above_t, connecting_to_t, maybe_do_not_collect);
            }
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

    void build_subcircuits(const std::vector<mockturtle::node<GateLyt>>   root_connections,
                           std::vector<mockturtle::node<GateLyt>>&        path,
                           std::unordered_set<mockturtle::node<GateLyt>>& visited, const size_t depth,
                           const size_t max_depth, gate_lyt_window_map& result) const
    {
        // todo: check out symmetry

        if (depth == max_depth)
        {
            // Get tiles corresponding to path nodes
            std::vector<tile<GateLyt>>                    tiles{};
            std::unordered_set<mockturtle::node<GateLyt>> nodes{};
            for (const auto& node : path)
            {
                tiles.push_back(stats.gate_layout->get_tile(node));
                nodes.insert(node);
            }

            // filter duplicates (i.e.: RAB, RBA)
            for (const auto& [other_path, _] : result)
            {
                if (path.front() == other_path.front() &&
                    std::all_of(std::next(other_path.cbegin(), 1), other_path.cend(),
                                [&](const auto& n) { return nodes.count(n) > 0; }))
                {
                    return;
                }
            }

            result.insert(
                {path, sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                           std::cref(*circuit), tiles,
                           params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params}});

            return;
        }

        const tile<GateLyt> current_tile = stats.gate_layout->get_tile(path.back());

        std::vector<mockturtle::node<GateLyt>> connections = root_connections;
        collect_connecting_nodes(*stats.gate_layout, current_tile, connections);

        for (const auto& next_node : connections)
        {
            if (visited.count(next_node) > 0)
            {
                continue;
            }

            visited.insert(next_node);
            path.push_back(next_node);

            build_subcircuits(root_connections, path, visited, depth + 1, max_depth, result);

            path.pop_back();
            visited.erase(next_node);
        }
    }
    void build_subcircuits_from_root(const mockturtle::node<GateLyt> root, const size_t max_depth,
                                     gate_lyt_window_map& result) const
    {
        const tile<GateLyt> current_tile = stats.gate_layout->get_tile(root);

        std::vector<mockturtle::node<GateLyt>> root_connections;
        collect_connecting_nodes(*stats.gate_layout, current_tile, root_connections);

        std::vector<mockturtle::node<GateLyt>>        path    = {root};
        std::unordered_set<mockturtle::node<GateLyt>> visited = {root};

        build_subcircuits(root_connections, path, visited, 0, max_depth, result);
    }

    void collect_indices_to_trial(const std::vector<mockturtle::node<GateLyt>>& node_vec, uint64_t& num_trials,
                                  foreach_node<std::vector<uint64_t>>& sampled_indices) const noexcept
    {
        const size_t                       num_nodes = node_vec.size() - 1;  // skip the first node
        std::vector<std::vector<uint64_t>> trialable_indices_per_node;
        std::vector<uint64_t>              bases;  // base for each position (number of valid indices per node)

        uint64_t num_possible_trials = 1;
        for (auto node_vec_it = std::next(node_vec.cbegin(), 1); node_vec_it != node_vec.cend(); ++node_vec_it)
        {
            std::vector<uint64_t> trialable_indices;

            for (uint64_t i = 0; i < gate_designs.at(*node_vec_it).size(); ++i)
            {
                trialable_indices.emplace_back(i);
            }

            assert(!trialable_indices.empty() && "No valid indices to trial for node.");

            bases.push_back(trialable_indices.size());
            num_possible_trials *= trialable_indices.size();
            trialable_indices_per_node.emplace_back(std::move(trialable_indices));
        }

        const uint64_t actual_num_trials = std::min(num_trials, num_possible_trials);

        // Step 1: Sample unique indices in range [0, num_possible_trials)
        std::unordered_set<uint64_t>            sampled_flat_indices;
        std::mt19937                            rng(std::random_device{}());
        std::uniform_int_distribution<uint64_t> dist(0, num_possible_trials - 1);

        while (sampled_flat_indices.size() < actual_num_trials)
        {
            sampled_flat_indices.insert(dist(rng));
        }

        // Step 2: For each flat index, compute the mixed-base digit vector
        std::vector<std::vector<uint64_t>> combinations;
        combinations.reserve(actual_num_trials);

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
        size_t node_idx = 1;  // skip first node
        for (auto node_vec_it = std::next(node_vec.cbegin(), 1); node_vec_it != node_vec.cend();
             ++node_vec_it, ++node_idx)
        {
            std::vector<uint64_t> node_trials;
            node_trials.reserve(actual_num_trials);

            for (const auto& combo : combinations)
            {
                node_trials.emplace_back(combo[node_idx - 1]);  // offset by one
            }

            sampled_indices[*node_vec_it] = std::move(node_trials);
        }

        num_trials = actual_num_trials;  // update trial count to actual number sampled
    }

    [[nodiscard]] uint64_t perform_trial(const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& sub_circuit,
                                         const std::vector<mockturtle::node<GateLyt>>&                      node_vec,
                                         foreach_node<std::vector<uint64_t>>&& indices_to_trial,
                                         const uint64_t trial_number, double& logic_match_average_over_inputs,
                                         CellLyt& cell_lyt_clone) const noexcept
    {
        // assign random gate design to other gates
        for (typename std::vector<mockturtle::node<GateLyt>>::const_iterator node_vec_it =
                 std::next(node_vec.cbegin(), 1);
             node_vec_it != node_vec.cend(); ++node_vec_it)
        {
            assert(trial_number < indices_to_trial.at(*node_vec_it).size() &&
                   "Trial number exceeds the number of available samples");

            // select a random gate implementation
            assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(
                cell_lyt_clone, gate_designs.at(*node_vec_it).at(indices_to_trial.at(*node_vec_it).at(trial_number)),
                *stats.gate_layout, stats.gate_layout->get_tile(*node_vec_it));
        }

        const sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary> c{cell_lyt_clone, sub_circuit};

        // sub-circuit logic match assessment
        const circuit_operational_assessment<CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
            is_circuit_operational<CellLyt, GateLyt, local_external_potential_type::BOUNDED, SkeletonGateLibrary>(
                c, operational_params);

        assert(op_assessment.assessment_per_input.has_value() && "ALL_INPUT_COMBINATIONS_ENUMERATED is not set.");

        uint64_t successful_input_combinations  = 0;
        double   logic_match_sum_for_all_inputs = 0.0;

        // sum logic match over the different input combinations
        for (const typename circuit_operational_assessment<CellLyt, local_external_potential_type::BOUNDED>::
                 operational_assessment_for_input& op_assessment_for_input : *op_assessment.assessment_per_input)
        {
            if (op_assessment_for_input.status == operational_status::OPERATIONAL)
            {
                successful_input_combinations++;
            }

            logic_match_sum_for_all_inputs += op_assessment_for_input.logic_match;
        }

        logic_match_average_over_inputs +=
            logic_match_sum_for_all_inputs / static_cast<double>(op_assessment.assessment_per_input->size());

        return successful_input_combinations;
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
                                              const uint64_t num_trials, const uint64_t min_bound,
                                              const uint64_t                        max_bound,
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

        for (uint64_t i = 0; i < num_threads; ++i)
        {
            threads.emplace_back(
                [&gate_fitness_assessments, &maybe_lyt, &lyt_mutex, num_input_combinations, &t, &n, i,
#if (PROGRESS_BARS)
                 &bar,
#endif
                 min_bound, max_bound, chunk_size, this, &gate_lyt_windows, global_pruning, num_trials]
                {
                    const uint64_t start_index = min_bound + i * chunk_size;
                    const uint64_t end_index   = std::min(start_index + chunk_size, max_bound);

                    for (uint64_t j = start_index; j < end_index; ++j)
                    {
                        CellLyt cell_lyt{};

                        // select the j-th gate implementation for n
                        assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(cell_lyt, gate_designs.at(n).at(j),
                                                                           *stats.gate_layout, t);

                        uint64_t total_number_of_trials = 0;

                        double successful_trials = 0;

                        foreach_node<std::vector<uint64_t>> sampled_indices{};

                        // for each sub-circuit
                        for (const auto& [node_vec, gate_lyt_window] : gate_lyt_windows)
                        {
                            if (n != node_vec.front())
                            {
                                continue;
                            }

                            uint64_t actual_num_trials = num_trials;

                            collect_indices_to_trial(node_vec, actual_num_trials, sampled_indices);

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

                                const uint64_t successful_input_combinations =
                                    perform_trial(gate_lyt_window, node_vec, std::move(sampled_indices), current_trial,
                                                  successful_trials, cell_lyt_clone);

                                if (global_pruning && successful_input_combinations == num_input_combinations)
                                {
                                    // all input combinations are operational: operational circuit found

                                    std::lock_guard lock{lyt_mutex};

                                    if (!maybe_lyt.has_value())
                                    {
                                        maybe_lyt.emplace(cell_lyt_clone);
                                    }

                                    return;
                                }

                                current_trial++;
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

    static void apply_quantization(std::vector<gate_fitness_assessment>& sorted_fitness_assessments,
                                   const double quantization_factor, const double success_rate_ceiling) noexcept
    {
        // apply ceiling
        for (auto& fitness_assessment : sorted_fitness_assessments)
        {
            fitness_assessment.fitness = std::min(success_rate_ceiling, fitness_assessment.fitness);
        }

        if (sorted_fitness_assessments.empty() || quantization_factor <= 0.0)
        {
            return;
        }

        const double min_fitness = sorted_fitness_assessments.front().fitness;
        const double max_fitness = sorted_fitness_assessments.back().fitness;

        // If all values are identical, nothing to do.
        if (std::abs(min_fitness - max_fitness) < std::numeric_limits<double>::epsilon())
        {
            return;
        }

        const double scale = std::floor(1.0 / quantization_factor);

        for (auto& fitness_assessment : sorted_fitness_assessments)
        {
            fitness_assessment.fitness =
                std::min(success_rate_ceiling, std::round(fitness_assessment.fitness * scale) / scale);
        }
    }

    static void print_success_rate_distribution(const std::vector<gate_fitness_assessment>& fitness_assessments,
                                                const uint64_t first_passing_gate_ix) noexcept
    {
        const double selectivity_threshold = fitness_assessments.at(first_passing_gate_ix).fitness;

        std::cout << "\nDetermined success threshold: " << std::fixed << std::setprecision(1)
                  << selectivity_threshold * 100
                  << fmt::format("% — pruning {} out of {} gate designs — reduced pool by {:.1f}%\n\n",
                                 first_passing_gate_ix, fitness_assessments.size(),
                                 static_cast<double>(first_passing_gate_ix) /
                                     static_cast<double>(fitness_assessments.size()) * 100);

        // Map: success_rate -> count of gates with this rate
        std::map<double, uint64_t> count_by_success_rate;

        // Count how many gates per success rate
        for (const gate_fitness_assessment& fitness_assessment : fitness_assessments)
        {
            ++count_by_success_rate[fitness_assessment.fitness];
        }

        // Total number of gates
        const uint64_t total_gates = fitness_assessments.size();

        assert(total_gates > 0 && "The distribution is empty.");

        std::cout << "Success Rate | Status | Gate Count | Graph\n";
        std::cout
            << "-------------|--------|------------|--------------------------------------------------------------\n";

        for (const auto& [rate, count] : count_by_success_rate)
        {
            constexpr size_t max_bar_length = 60;

            const bool        is_kept = rate >= selectivity_threshold;
            const std::string status  = is_kept ? "KEPT" : "PRUNED";

            // Bar length proportional to total count
            const auto bar_len = 1 + static_cast<size_t>(static_cast<double>(count) / static_cast<double>(total_gates) *
                                                         (max_bar_length - 1));

            std::string bar(bar_len, is_kept ? 'o' : 'X');

            std::cout << std::right << std::setw(11) << std::fixed << std::setprecision(1) << rate * 100 << "%"
                      << " | " << std::setw(6) << status << " | " << std::setw(10) << count << " | " << std::left
                      << std::setw(max_bar_length) << bar << '\n';
        }

        std::cout << std::endl;
    }

    [[nodiscard]] bool discriminate_fitness_assessments(
        const double selectivity, const double quantization_factor, const double success_rate_ceiling,
        std::vector<typename SkeletonGateLibrary::fcn_gate>& remaining_gate_designs,
        std::vector<gate_fitness_assessment>& gate_fitness_assessments, uint64_t& min_bound, uint64_t& max_bound,
        uint64_t& attempt_number, bool& big_fixpoint) const noexcept
    {
        std::sort(gate_fitness_assessments.begin(), gate_fitness_assessments.end(),
                  [](const auto& lhs, const auto& rhs) { return lhs.fitness < rhs.fitness; });

        auto threshold_ix =
            static_cast<uint64_t>(std::round(static_cast<double>(gate_fitness_assessments.size()) * selectivity));

        const auto threshold_val_is_above_success_rate_ceiling = [&]
        {
            return gate_fitness_assessments.at(threshold_ix).fitness >
                   success_rate_ceiling - std::numeric_limits<double>::epsilon();
        };

        while (threshold_ix > 0 && threshold_val_is_above_success_rate_ceiling())
        {
            --threshold_ix;
        }

        const auto lb_ix = static_cast<uint64_t>(std::distance(
            gate_fitness_assessments.cbegin(),
            std::lower_bound(gate_fitness_assessments.cbegin(), gate_fitness_assessments.cend(),
                             gate_fitness_assessments.at(threshold_ix).fitness,

                             [](const gate_fitness_assessment& fitness_assessment, const double& val)
                             { return fitness_assessment.fitness < val - std::numeric_limits<double>::epsilon(); })));
        const auto ub_ix = static_cast<uint64_t>(std::distance(
            gate_fitness_assessments.cbegin(),
            std::upper_bound(gate_fitness_assessments.cbegin(), gate_fitness_assessments.cend(),
                             gate_fitness_assessments.at(threshold_ix).fitness,
                             [](const double val, const gate_fitness_assessment& fitness_assessment)
                             { return val + std::numeric_limits<double>::epsilon() < fitness_assessment.fitness; })));

        if (ub_ix - lb_ix == 1 || threshold_val_is_above_success_rate_ceiling() ||
            attempt_number == params.maximum_discrimination_attempts)
        {
            apply_quantization(gate_fitness_assessments, quantization_factor, success_rate_ceiling);

            const uint64_t first_passing_ix = ub_ix - threshold_ix >= threshold_ix - lb_ix ? lb_ix : ub_ix;

            print_success_rate_distribution(gate_fitness_assessments, first_passing_ix);

            for (uint64_t gate_index = first_passing_ix; gate_index < gate_fitness_assessments.size(); ++gate_index)
            {
                gate_fitness_assessments.at(gate_index).selected = true;
            }

            if (first_passing_ix != 0)
            {
                big_fixpoint = false;
            }

            return true;
        }

        const auto reorder_in_place = [](std::vector<typename SkeletonGateLibrary::fcn_gate>& data,
                                         const std::vector<gate_fitness_assessment>&          order)
        {
            std::vector visited(data.size(), false);

            for (size_t i = 0; i < data.size(); ++i)
            {
                if (visited[i] || order[i].gate_index == i)
                {
                    continue;  // already in place or visited
                }

                size_t j = i;

                typename SkeletonGateLibrary::fcn_gate temp = std::move(data[i]);

                // Follow the cycle
                while (!visited[j])
                {
                    visited[j] = true;

                    size_t next = order[j].gate_index;

                    data[j] = std::move(next == i ? temp : data[next]);

                    j = next;
                }
            }
        };

        reorder_in_place(remaining_gate_designs, gate_fitness_assessments);

        for (uint64_t i = 0; i < gate_fitness_assessments.size(); ++i)
        {
            gate_fitness_assessments[i].gate_index = i;
        }

        if (min_bound != lb_ix || max_bound != ub_ix)
        {
            attempt_number = 0;
        }
        else
        {
            ++attempt_number;
        }

        min_bound = lb_ix;
        max_bound = ub_ix;

        return false;
    }

    /**
     * todo
     */
    std::optional<CellLyt> prune_gate_designs(const uint64_t level, const std::string_view& level_str) noexcept
    {
        std::cout << "\n\nSTARTING TO PRUNE GATE DESIGNS\tLEVEL: " << level_str << std::endl;

        bool global_pruning = false;

        if (level_str == "GLOBAL")
        {
            global_pruning = true;
        }

        std::optional<CellLyt> maybe_lyt{};

        const uint64_t num_input_combinations = 1 << stats.gate_layout->num_pis();

        // todo
        const uint64_t num_trials           = std::max(uint64_t{1}, params.num_trials / level);
        const double   quantization_factor  = params.quantization_factor * std::pow(0.8, level - 1);
        const double   selectivity          = params.selectivity * std::pow(1.2, level - 1);
        const double   success_rate_ceiling = params.success_rate_ceiling + static_cast<double>(level - 1) * 0.01;

        std::cout << "\nNumber of trials: " << num_trials << std::endl;
        std::cout << "Quantization factor: " << fmt::format("{:.3f}", quantization_factor) << std::endl;
        std::cout << "Selectivity: " << fmt::format("{:.2f}", selectivity) << std::endl;
        std::cout << "Success_rate_ceiling: " << fmt::format("{:.2f}", success_rate_ceiling) << std::endl;

        gate_lyt_window_map gate_lyt_windows{};

        foreach_node<std::vector<gate_fitness_assessment>> gate_fitness_assessments{};

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

        bool big_fixpoint = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (!skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        gate_fitness_assessments[n].resize(gate_designs.at(n).size());
                    }
                });

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n) || gate_designs.at(n).size() == 1)
                    {
                        return;
                    }

                    const tile<GateLyt>& t = stats.gate_layout->get_tile(n);

                    uint64_t number_of_layouts = 0;

                    for (const auto& [node_vec, _] : gate_lyt_windows)
                    {
                        assert(node_vec.size() > 1 && "The connected nodes vector cannot be singleton");

                        if (n != node_vec.front())
                        {
                            continue;
                        }

                        number_of_layouts++;
                    }

                    std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                    std::cout << fmt::format("Performing {} trials for {} sub-circuit{} for {} gate implementations",
                                             num_trials, number_of_layouts, number_of_layouts > 1 ? "s" : "",
                                             gate_designs.at(n).size())
                              << std::endl;

                    uint64_t min_bound = 0;                          // inclusive
                    uint64_t max_bound = gate_designs.at(n).size();  // exclusive

                    uint64_t attempt_number = 0;

                    while (true)
                    {
                        make_trial_based_fitness_assessments(gate_lyt_windows, n, t, num_trials, min_bound, max_bound,
                                                             gate_fitness_assessments[n], global_pruning,
                                                             num_input_combinations, lyt_mutex, maybe_lyt);

                        if (global_pruning)
                        {
                            if (maybe_lyt.has_value())
                            {
                                return;  // quit if an operational circuit has been found
                            }
                        }

                        if (discriminate_fitness_assessments(selectivity, quantization_factor, success_rate_ceiling,
                                                             gate_designs[n], gate_fitness_assessments[n], min_bound,
                                                             max_bound, attempt_number, big_fixpoint))
                        {
                            break;
                        }
                    }
                });

            std::cout << std::endl;

            if (global_pruning)
            {
                if (maybe_lyt.has_value())
                {
                    return maybe_lyt.value();  // return the operational circuit
                }
            }

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n) || gate_designs.at(n).size() == 1)
                    {
                        return;
                    }

                    std::vector<typename SkeletonGateLibrary::fcn_gate> selected_gate_implementations{};

                    for (const gate_fitness_assessment& fitness_assessment : gate_fitness_assessments.at(n))
                    {
                        if (fitness_assessment.selected)
                        {
                            selected_gate_implementations.push_back(
                                std::move(gate_designs[n][fitness_assessment.gate_index]));
                        }
                    }

                    gate_designs[n] = std::move(selected_gate_implementations);
                });
        }

        return std::nullopt;
    }
    /**
     * todo
     */
    std::optional<CellLyt> exhaustively_enumerate_gate_design_combinations() noexcept
    {
        std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;

        operational_params.print = true;  // after comment: reenable function constness

        std::vector<uint64_t> indices(gate_designs.size(), 0);

        while (true)
        {
            CellLyt operational_circuit_candidate{};
            for (uint64_t i = 0; i < gate_designs.size(); i++)
            {
                const auto& [n, op_gate_designs_for_gate] = *std::next(gate_designs.cbegin(), static_cast<int64_t>(i));
                // select a random gate implementation for the tile that connects as input to n
                assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(operational_circuit_candidate,
                                                                   op_gate_designs_for_gate.at(indices.at(i)),
                                                                   *stats.gate_layout, stats.gate_layout->get_tile(n));
            }

            std::cout << "trying combination: ";
            for (uint64_t i = 0; i < gate_designs.size(); i++)
            {
                std::cout << indices.at(i) << " ";
            }
            std::cout << std::endl;

            if (is_circuit_operational(
                    sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                        operational_circuit_candidate,
                        sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{*circuit}},
                    operational_params)
                    .status == operational_status::OPERATIONAL)
            {
                std::cout << "\n\nFINAL GENERATED CIRCUIT:" << std::endl;
                print_layout(operational_circuit_candidate);

                return operational_circuit_candidate;
            }

            // Increment indices like an odometer
            for (uint64_t i = 0; i < indices.size(); ++i)
            {
                if (++indices[i] < std::next(gate_designs.cbegin(), static_cast<int64_t>(i))->second.size())
                {
                    break;  // No carry needed
                }

                indices[i] = 0;  // Reset this index and carry over to the next

                if (i == indices.size() - 1)
                {
                    return std::nullopt;  // Stop when the last index overflows
                }
            }
        }

        return std::nullopt;
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
[[nodiscard]] std::optional<CellLyt> advanced_circuit_design(const Ntk& ntk, const GateLyt& lattice_tiling,
                                                             advanced_circuit_design_params<CellLyt>& params = {},
                                                             advanced_circuit_design_stats<GateLyt>*  stats  = nullptr)
{
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt is not a hexagonal");
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");

    advanced_circuit_design_stats<GateLyt> st{};

    detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt, SkeletonGateLibrary> p{ntk, params, lattice_tiling, st};

    const auto result = p.design_circuit();

    if (stats)
    {
        *stats = st;
    }

    return result;
}

}  // namespace fiction

#endif  // advanced_CIRCUIT_DESIGN_HPP

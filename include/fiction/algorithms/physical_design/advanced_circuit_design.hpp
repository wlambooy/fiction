//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_circuit_operational.hpp"
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
    design_sidb_gates_params<CellLyt> design_gate_params{};
    /** This variable specifies the radius in nanometers around the center of the hexagon where atomic defects are
     * incorporated into the gate design. (unit: nm)
     */

    std::optional<CellLyt> defect_surface{};
    double                 influence_radius_charged_defects = 15;

    uint64_t num_trials = 500;

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
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;

        // initialize
        collect_initial_gate_designs();

        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        while (++circuit_design_level < num_gates_to_design - 1)
        {
            // prune by assessing gate design combinations for increasingly large sets of connected gates
            prune_gate_designs(
                circuit_design_level,
                std::string_view{
                    std::to_string(circuit_design_level + 1) +
                    (params.sub_circuit_mode ==
                             advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::CONNECTED_GATES ?
                         " CONNECTED GATES" :
                     params.sub_circuit_mode ==
                             advanced_circuit_design_params<CellLyt>::sub_circuit_creation_mode::ADJACENT_GATES ?
                         " ADJACENT GATES" :
                         "GATES")});
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

        circuit.emplace(*stats.gate_layout, params.design_gate_params.operational_params.simulation_parameters,
                        params.design_gate_params.canvas,
                        params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params);

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
    const advanced_circuit_design_params<CellLyt> params{};
    /**
     * Statistics for the on-the-fly circuit design.
     */
    advanced_circuit_design_stats<GateLyt>&                                stats;
    std::optional<sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>> circuit{};

    is_circuit_operational_params operational_params{};

    template <typename T>
    using foreach_node = std::unordered_map<mockturtle::node<GateLyt>, T>;

    using gate_design_t         = std::vector<uint64_t>;
    using gate_designs_per_node = foreach_node<std::vector<gate_design_t>>;

    gate_designs_per_node gate_designs{};

    uint64_t circuit_design_level = 0;

    [[nodiscard]] std::vector<gate_design_t>
    design_gates_randomly(const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& sub_circuit,
                          const mockturtle::node<GateLyt>&                                   n) const noexcept
    {
        const std::vector<cell<CellLyt>>& all_canvas_positions = circuit->all_canvas_positions.at(n);
        const uint64_t                    number_of_layouts    = binomial_coefficient(
            all_canvas_positions.size(), params.design_gate_params.number_of_canvas_sidbs);  // todo sum 0 to num

        std::vector<gate_design_t> randomly_designed_gate_layouts = {};

        std::vector<std::thread> threads{};
        threads.reserve(params.design_gate_params.available_threads);

        std::mutex mutex_to_protect_designed_gate_layouts{};  // used to control access to shared resources

        // const auto check_if_gate_design_is_already_present = [&](const Lyt& gate_design)
        // {
        //     for (const Lyt& stored_gate_design : randomly_designed_gate_layouts)
        //     {
        //         if (std::all_of(all_sidbs_in_canvas.cbegin(), all_sidbs_in_canvas.cend(), [&](const cell<Lyt>& sidb)
        //                         { return gate_design.get_cell_type(sidb) == stored_gate_design.get_cell_type(sidb);
        //                         }))
        //         {
        //             return true;
        //         }
        //     }
        //
        //     return false;
        // };

        std::atomic<uint64_t> num_solutions_found = 0;

        const uint64_t max_number_of_solutions =
            std::min(params.design_gate_params.maximum_number_of_solutions, number_of_layouts);

#if (PROGRESS_BARS)
        // initialize a progress bar
        mockturtle::progress_bar bar{
            static_cast<uint32_t>(max_number_of_solutions),
            fmt::format("[i] looking for {} random operational gate designs: ", max_number_of_solutions) +
                "|{0}|                            "};
#endif

        for (uint64_t z = 0u; z < params.design_gate_params.available_threads; z++)
        {
            threads.emplace_back(
                [this, &num_solutions_found, &max_number_of_solutions,
                 &sub_circuit,  // &check_if_gate_design_is_already_present,
                 &mutex_to_protect_designed_gate_layouts, &all_canvas_positions, &randomly_designed_gate_layouts
#if (PROGRESS_BARS)
                 ,
                 &bar
#endif
            ]
                {
                    while (num_solutions_found < max_number_of_solutions)
                    {

                        static std::random_device              rd;
                        static std::mt19937                    gen(rd());
                        static std::uniform_int_distribution<> dist(
                            1, static_cast<int32_t>(params.design_gate_params.number_of_canvas_sidbs));

                        const auto number_of_sidbs_of_final_layout = static_cast<uint64_t>(dist(gen));

                        static std::mt19937_64 generator{std::random_device{}()};

                        // instantiate distribution
                        std::uniform_int_distribution<std::size_t> distributions{0, all_canvas_positions.size() - 1};

                        // container for the random samples
                        phmap::btree_set<uint64_t> canvas_positions{};

                        for (std::size_t i = 0; i < number_of_sidbs_of_final_layout; ++i)
                        {
                            canvas_positions.insert(distributions(generator));
                        }

                        CellLyt random_lyt = sub_circuit.skeleton.clone();
                        for (const uint64_t gate_design_cell_index : canvas_positions)
                        {
                            random_lyt.assign_cell_type(all_canvas_positions.at(gate_design_cell_index),
                                                        sidb_technology::cell_type::LOGIC);
                        }

                        if (is_circuit_operational<CellLyt, GateLyt, local_external_potential_type::BOUNDED,
                                                   SkeletonGateLibrary>(
                                sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{random_lyt,
                                                                                                   sub_circuit},
                                operational_params)
                                .status != operational_status::OPERATIONAL)
                        {
                            continue;
                        }

                        const std::lock_guard lock{mutex_to_protect_designed_gate_layouts};

                        randomly_designed_gate_layouts.emplace_back(canvas_positions.cbegin(), canvas_positions.cend());

                        ++num_solutions_found;

#if (PROGRESS_BARS)
                        if (num_solutions_found < params.design_gate_params.maximum_number_of_solutions)
                        {
                            // update the progress bar
                            bar(num_solutions_found.load());
                        }
#endif
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

        return randomly_designed_gate_layouts;
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
                    if (skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

                    std::cout << "starting gate design for tile " << stats.gate_layout->get_tile(n)
                              << "\t|\tnode function:";
                    print_binary(stats.gate_layout->node_function(n));
                    std::cout << std::endl;

                    gate_designs[n] = design_gates_randomly(
                        sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                            std::cref(*circuit),
                            {stats.gate_layout->get_tile(n)},
                            params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params},
                        n);

                    std::cout << "number of gate layouts found: " << gate_designs.at(n).size() << std::endl;
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
                {path, sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                           std::cref(*circuit), tiles,
                           params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params}});

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

    void collect_indices_to_trial(const mockturtle::node<GateLyt>&              n,
                                  const std::vector<mockturtle::node<GateLyt>>& node_vec, uint64_t& num_trials,
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

            for (uint64_t i = 0; i < gate_designs.at(node).size(); ++i)
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
        for (uint64_t node_idx = 0, other_node_idx = 0; node_idx < node_vec.size(); ++node_idx)
        {
            if (node_vec.at(node_idx) == n)
            {
                continue;
            }

            std::vector<uint64_t> node_trials;
            node_trials.reserve(actual_num_trials);

            for (const auto& combo : combinations)
            {
                node_trials.emplace_back(combo[other_node_idx]);
            }

            ++other_node_idx;

            sampled_indices[node_vec[node_idx]] = std::move(node_trials);
        }

        num_trials = actual_num_trials;  // update trial count to actual number sampled
    }

    [[nodiscard]] uint64_t perform_trial(const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& sub_circuit,
                                         const mockturtle::node<GateLyt>&                                   n,
                                         const std::vector<mockturtle::node<GateLyt>>&                      node_vec,
                                         const foreach_node<std::vector<uint64_t>>& indices_to_trial,
                                         const uint64_t trial_number, double& logic_match_average_over_inputs,
                                         CellLyt&                                     cell_lyt_clone,
                                         const std::unique_ptr<thread_count_manager>& tcm) const noexcept
    {
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
                 gate_designs.at(node).at(indices_to_trial.at(node).at(trial_number)))
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
        // number of gate designs we need to evaluate (indices in [min_bound, max_bound))
        const uint64_t num_gate_designs = max_bound - min_bound;

        // ---------------------------
        // 1) Collect only the sub-circuits that include node 'n'
        // ---------------------------
        struct SubEntry
        {
            std::vector<mockturtle::node<GateLyt>>                             node_vec;
            const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>* sub_ptr;
        };
        std::vector<SubEntry> sub_entries;
        sub_entries.reserve(gate_lyt_windows.size());

        for (const auto& [node_vec, gate_lyt_window] : gate_lyt_windows)
        {
            if (std::find(node_vec.cbegin(), node_vec.cend(), n) != node_vec.cend())
            {
                // store pointer to the gate_lyt_window (sub-circuit)
                sub_entries.push_back({node_vec, &gate_lyt_window});
            }
        }

        const size_t num_subs = sub_entries.size();
        if (num_subs == 0)
        {
            // no sub-circuits include node n: every gate design has zero trials -> fitness 0
            for (uint64_t j = 0; j < num_gate_designs; ++j)
            {
                const uint64_t idx                       = min_bound + j;
                gate_fitness_assessments[idx].gate_index = idx;
                gate_fitness_assessments[idx].fitness    = 0.0;
                gate_fitness_assessments[idx].selected   = false;
            }
            return;
        }

        // ---------------------------
        // 2) Precompute per (gate_design j, sub s) the sampled indices, actual_num_trials and a prepared CellLyt
        //    We'll store them in a flat vector indexed by (j_idx * num_subs + s)
        //    where j_idx runs 0..num_gate_designs-1 corresponding to absolute index (min_bound + j_idx)
        // ---------------------------
        struct PairPrep
        {
            std::vector<mockturtle::node<GateLyt>> node_vec;
            foreach_node<std::vector<uint64_t>>    sampled_indices;
            uint64_t                               actual_num_trials{0};
            CellLyt cell_lyt;  // skeleton + assigned logic cells for this gate design & sub-circuit
            const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>* sub_ptr{nullptr};
        };

        std::vector<PairPrep> prep;
        prep.reserve(static_cast<size_t>(num_gate_designs) * num_subs);

        // total number of trial jobs across all pairs
        uint64_t total_number_of_trials_all = 0;

        for (uint64_t j_idx = 0; j_idx < num_gate_designs; ++j_idx)
        {
            const uint64_t j = min_bound + j_idx;  // absolute index in circuit->gate_designs

            const gate_design_t& gate_design = gate_designs.at(n).at(j);

            for (size_t s = 0; s < num_subs; ++s)
            {
                PairPrep pp;
                pp.node_vec = sub_entries[s].node_vec;
                pp.sub_ptr  = sub_entries[s].sub_ptr;

                // clone skeleton for this sub-circuit and assign cells for this gate design
                pp.cell_lyt = pp.sub_ptr->skeleton.clone();

                for (const uint64_t gate_design_cell_index : gate_design)
                {
                    // original code used stats.gate_layout and t to assign cells; keep same behavior:
                    // assign_gate previously used assign_cell_type by position; here we mimic assignment by
                    // using the stored positions in circuit->all_canvas_positions as in forked code.
                    pp.cell_lyt.assign_cell_type(circuit->all_canvas_positions.at(n).at(gate_design_cell_index),
                                                 sidb_technology::cell_type::LOGIC);
                }

                // determine how many trials (and sampled indices) for this (j, s)
                // We must call the same helper as original: collect_indices_to_trial(n, node_vec, actual_num_trials,
                // sampled_indices) but that helper in original used a reference to actual_num_trials and the
                // sampled_indices. To preserve behavior we call it similarly.
                uint64_t actual_num_trials = num_trials;
                collect_indices_to_trial(n, pp.node_vec, actual_num_trials, pp.sampled_indices);
                pp.actual_num_trials = actual_num_trials;

                total_number_of_trials_all += pp.actual_num_trials;
                prep.push_back(std::move(pp));
            }
        }

        // ---------------------------
        // 3) Build a job list where each job is one trial
        // ---------------------------
        struct TrialJob
        {
            uint32_t gate_design_jidx;  // j index relative to min_bound: 0..num_gate_designs-1
            uint32_t sub_index;         // s index 0..num_subs-1
            uint64_t trial_number;      // 0..actual_num_trials-1 for that pair
        };

        std::vector<TrialJob> jobs;
        jobs.reserve(static_cast<size_t>(total_number_of_trials_all));

        // also keep total trials per gate design for final reduction
        std::vector<uint64_t> total_trials_per_j;
        total_trials_per_j.assign(static_cast<size_t>(num_gate_designs), 0);

        for (uint64_t j_idx = 0; j_idx < num_gate_designs; ++j_idx)
        {
            for (uint32_t s = 0; s < static_cast<uint32_t>(num_subs); ++s)
            {
                const size_t   idx    = static_cast<size_t>(j_idx) * num_subs + s;
                const uint64_t actual = prep[idx].actual_num_trials;
                total_trials_per_j[j_idx] += actual;

                for (uint64_t tr = 0; tr < actual; ++tr)
                {
                    jobs.push_back({static_cast<uint32_t>(j_idx), s, tr});
                }
            }
        }

        // number of worker threads to spawn (keep original policy)
        const uint64_t num_threads = std::min(params.available_threads, jobs.size());

#if (PROGRESS_BARS)
        // progress bar will be driven by processed_jobs counter below
        mockturtle::progress_bar bar{static_cast<uint32_t>(std::min(jobs.size(), static_cast<uint64_t>(UINT32_MAX))),
                                     "[i] Determining successful trial ratio for tile " +
                                         fmt::format("({},{})", t.x, t.y) + ": |{0}|"};
#endif

        // ---------------------------
        // 5) Per-j accumulators for successful trials (count of successful trials)
        //    We'll use atomic<uint64_t> counts because perform_trial will be called with a local double
        //    and we will add the integer count (0 or 1) to this atomic. This preserves the final
        //    successful_trial_ratio = successful_trials / total_number_of_trials behavior.
        // ---------------------------
        std::vector<double> successful_trial_counts;
        successful_trial_counts.resize(num_gate_designs);
        std::vector<std::mutex> successful_trial_mutexes(num_gate_designs);

        // job index and processed jobs counters
        std::atomic<uint64_t> job_index{0};
        std::atomic<uint64_t> processed_jobs{0};

        // global flag indicating we've found an operational circuit (for global pruning)
        std::atomic<uint8_t> global_found_operational{0};

        // A thread_count_manager shared as before
        std::unique_ptr<thread_count_manager> tcm =
            // circuit_design_level > 2 ?
            // std::make_unique<thread_count_manager>(params.available_threads - num_threads) :
            nullptr;

        // launch workers
        std::vector<std::thread> threads;
        threads.reserve(num_threads);

        for (uint64_t thread_id = 0; thread_id < num_threads; ++thread_id)
        {
            threads.emplace_back(
                [&, thread_id]
                {
                    // A small local value to collect per-thread simulator calls if needed by perform_trial
                    // implementation. (we do not rely on it for correctness here)
                    while (true)
                    {
                        const uint64_t my_job = job_index.fetch_add(1, std::memory_order_relaxed);
                        if (my_job >= jobs.size())
                            break;

                        // if global pruning is enabled and some other thread already found an operational circuit ->
                        // exit early
                        if (global_pruning && global_found_operational.load(std::memory_order_acquire))
                        {
                            // We still increment processed_jobs for progress bookkeeping
                            processed_jobs.fetch_add(1, std::memory_order_relaxed);
#if (PROGRESS_BARS)
                            if (thread_id == 0)
                                bar(processed_jobs.load());
#endif
                            break;
                        }

                        const TrialJob& job      = jobs[my_job];
                        const size_t    pair_idx = static_cast<size_t>(job.gate_design_jidx) * num_subs + job.sub_index;

                        // Grab the precomputed pair
                        const PairPrep& pp = prep[pair_idx];

                        // Prepare cell layout clone for this trial (like original code)
                        CellLyt cell_lyt_clone = pp.cell_lyt.clone();

                        // local accumulator passed to perform_trial (preserve the original semantics):
                        double local_successful_trials = 0.0;

                        // Call perform_trial with the same semantics as the original function.
                        // Note: perform_trial in the original added to 'successful_trials' and returned
                        // the number of successful input combinations for the trial. We keep the same call.
                        const uint64_t successful_input_combinations =
                            perform_trial(*pp.sub_ptr, n, pp.node_vec, pp.sampled_indices, job.trial_number,
                                          local_successful_trials, cell_lyt_clone, tcm);

                        // If this trial yields an entirely operational input-combination (all input combinations work),
                        // handle global pruning: set maybe_lyt under mutex and notify other threads to stop.
                        if (global_pruning && successful_input_combinations == num_input_combinations)
                        {
                            std::lock_guard lock{lyt_mutex};
                            if (!maybe_lyt.has_value())
                            {
                                maybe_lyt.emplace(cell_lyt_clone);
                            }
                            // mark global found so other threads stop pulling jobs
                            global_found_operational.store(1, std::memory_order_release);

                            // We increment processed_jobs and break out (other threads will observe global flag).
                            processed_jobs.fetch_add(1, std::memory_order_relaxed);
#if (PROGRESS_BARS)
                            if (thread_id == 0)
                                bar(processed_jobs.load());
#endif
                            break;
                        }

                        {
                            const std::lock_guard lock{successful_trial_mutexes[job.gate_design_jidx]};

                            successful_trial_counts[job.gate_design_jidx] += local_successful_trials;
                        }

                        // If for this pair some thread already discovered operational combination concurrently, mark
                        // accordingly (if perform_trial determined pair operational in some other way, it would have
                        // resulted in successful_input_combinations == num_input_combinations and handled above). if
                        // current thread finds pair operational by other criteria (not expected), we could mark
                        // pair_found_operational. For safety: if local_successes_as_uint > 0 AND actual_num_trials for
                        // this pair equals 1 and local_successes_as_uint equals actual num, we could mark pair found —
                        // but we avoid changing pair_found_operational here except in the explicit global_pruning case.

                        processed_jobs.fetch_add(1, std::memory_order_relaxed);
#if (PROGRESS_BARS)
                        if (thread_id == 0)
                            bar(processed_jobs.load());
#endif

                        // indicate this worker will return one thread slot (mimic original behavior)
                        if (tcm)
                        {
                            tcm->return_threads(1);
                        }
                    }  // while jobs
                });
        }

        // join workers
        for (auto& th : threads)
        {
            if (th.joinable())
                th.join();
        }

        // ---------------------------
        // 6) Reduce results to gate_fitness_assessments (map back j_idx -> absolute j)
        // ---------------------------
        for (uint64_t j_idx = 0; j_idx < num_gate_designs; ++j_idx)
        {
            const uint64_t abs_j        = min_bound + j_idx;
            const uint64_t total_trials = total_trials_per_j[j_idx];

            gate_fitness_assessments[abs_j].gate_index = abs_j;
            gate_fitness_assessments[abs_j].fitness =
                successful_trial_counts[j_idx] / static_cast<double>(total_trials);
            gate_fitness_assessments[abs_j].selected = false;
        }
    }

    static void apply_quantization(std::vector<gate_fitness_assessment>& sorted_fitness_assessments,
                                   const double quantization_step, const double success_rate_ceiling) noexcept
    {
        // apply ceiling
        for (auto& fitness_assessment : sorted_fitness_assessments)
        {
            fitness_assessment.fitness = std::min(success_rate_ceiling, fitness_assessment.fitness);
        }

        if (sorted_fitness_assessments.empty() || quantization_step <= 0.0)
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

        const double scale = std::floor(100.0 / quantization_step);

        for (auto& fitness_assessment : sorted_fitness_assessments)
        {
            fitness_assessment.fitness =
                std::min(success_rate_ceiling, std::round(fitness_assessment.fitness * scale) / scale);
        }
    }

    struct DoubleComparator
    {
        bool operator()(const double lhs, const double rhs) const
        {
            return lhs + std::numeric_limits<double>::epsilon() < rhs;
        }
    };

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
        std::map<double, uint64_t, DoubleComparator> count_by_success_rate;

        // Count how many gates per success rate
        for (const gate_fitness_assessment& fitness_assessment : fitness_assessments)
        {
            ++count_by_success_rate[fitness_assessment.fitness];
        }

        // Total number of gates
        const uint64_t total_gates = fitness_assessments.size();

        assert(total_gates > 0 && "The distribution is empty.");

        std::cout << "Success Rate | Status | Gate Count | Graph\n";
        std::cout << "-------------|--------|------------|-------------------------------------"
                     "-------------------------\n";

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
        bool global_pruning, const double selectivity, const double quantization_step,
        const double success_rate_ceiling, std::vector<gate_design_t>& remaining_gate_designs,
        std::vector<gate_fitness_assessment>& gate_fitness_assessments, uint64_t& min_bound, uint64_t& max_bound,
        uint64_t& repeated_attempt_number, uint64_t& attempt_number, bool& completed_assessment) const noexcept
    {
        std::sort(gate_fitness_assessments.begin(), gate_fitness_assessments.end(),
                  [](const auto& lhs, const auto& rhs) { return lhs.fitness < rhs.fitness; });

        auto threshold_ix =
            static_cast<uint64_t>(std::round(static_cast<double>(gate_fitness_assessments.size()) * selectivity));

        const auto threshold_val_is_above_success_rate_ceiling = [&](const uint8_t offset)
        {
            return gate_fitness_assessments.at(threshold_ix - offset).fitness >
                   success_rate_ceiling - std::numeric_limits<double>::epsilon();
        };

        while (threshold_ix > 0 && threshold_val_is_above_success_rate_ceiling(1))
        {
            --threshold_ix;
        }

        if (params.quantize_mode == advanced_circuit_design_params<CellLyt>::quantization_mode::SOFTEN_PRUNING)
        {
            apply_quantization(gate_fitness_assessments, quantization_step, success_rate_ceiling);
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

        if (ub_ix - lb_ix == 1 || threshold_val_is_above_success_rate_ceiling(0) ||
            repeated_attempt_number >= params.maximum_repeated_discrimination_attempts ||
            attempt_number >= params.maximum_discrimination_attempts)
        {
            const uint64_t first_passing_ix = ub_ix == gate_fitness_assessments.size() ||
                                                      threshold_val_is_above_success_rate_ceiling(0) ||
                                                      (lb_ix != 0 && ub_ix - threshold_ix >= threshold_ix - lb_ix) ?
                                                  lb_ix :
                                                  ub_ix;

            if (params.quantize_mode == advanced_circuit_design_params<CellLyt>::quantization_mode::VISUALIZATION_ONLY)
            {
                apply_quantization(gate_fitness_assessments, quantization_step, success_rate_ceiling);
            }

            print_success_rate_distribution(gate_fitness_assessments, first_passing_ix);

            for (uint64_t gate_index = first_passing_ix; gate_index < gate_fitness_assessments.size(); ++gate_index)
            {
                gate_fitness_assessments.at(gate_index).selected = true;
            }

            if ((global_pruning && first_passing_ix == 0) ||
                (!global_pruning &&
                 static_cast<double>(first_passing_ix) / static_cast<double>(gate_fitness_assessments.size()) <
                     (1.0 - params.selectivity_tolerance) * selectivity))
            {
                std::cout << "ASSESSMENT COMPLETED\n" << std::endl;

                completed_assessment = true;
            }

            return true;
        }

        const auto reorder_in_place =
            [](std::vector<gate_design_t>& data, const std::vector<gate_fitness_assessment>& order)
        {
            std::vector visited(data.size(), false);

            for (size_t i = 0; i < data.size(); ++i)
            {
                if (visited[i] || order[i].gate_index == i)
                {
                    continue;  // already in place or visited
                }

                size_t j = i;

                gate_design_t temp = std::move(data[i]);

                // Follow the cycle
                while (!visited[j])
                {
                    visited[j] = true;

                    const size_t next = order[j].gate_index;

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
            repeated_attempt_number = 1;
        }
        else
        {
            ++repeated_attempt_number;
        }

        ++attempt_number;

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

        const bool global_pruning = level_str == "GLOBAL";

        std::optional<CellLyt> maybe_lyt{};

        const uint64_t num_input_combinations = 1 << stats.gate_layout->num_pis();

        // todo
        const uint64_t num_trials           = std::max(uint64_t{1}, params.num_trials / level);
        const double   quantization_step    = params.quantization_step * std::pow(0.8, level - 1);
        const double   selectivity          = params.selectivity * std::pow(1.2, level - 1);
        const double   success_rate_ceiling = params.success_rate_ceiling + static_cast<double>(level - 1) * 0.01;

        std::cout << "\nNumber of trials:     " << num_trials << std::endl;
        std::cout << "Quantization step:    " << fmt::format("{:.2f}", quantization_step) << std::endl;
        std::cout << "Selectivity:          " << fmt::format("{:.2f}", selectivity) << std::endl;
        std::cout << "Success rate ceiling: " << fmt::format("{:.2f}", success_rate_ceiling) << std::endl;

        gate_lyt_window_map gate_lyt_windows{};

        foreach_node<std::vector<gate_fitness_assessment>> gate_fitness_assessments{};

        foreach_node<bool> completed_assessment{};

        stats.gate_layout->foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(*stats.gate_layout, n))
                {
                    return;
                }

                completed_assessment[n] = gate_designs.at(n).size() == 1;

                build_subcircuits_from_root(n, level, gate_lyt_windows);
            });

        std::mutex lyt_mutex{};

        while (!std::all_of(completed_assessment.cbegin(), completed_assessment.cend(),
                            [](const auto& kv) { return kv.second; }))
        {
            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

                    if (completed_assessment.at(n))
                    {
                        gate_fitness_assessments[n].clear();
                    }
                    else
                    {
                        gate_fitness_assessments[n].resize(gate_designs.at(n).size());
                    }
                });

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(*stats.gate_layout, n) || completed_assessment.at(n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = stats.gate_layout->get_tile(n);

                    uint64_t number_of_layouts = 0;

                    for (const auto& [node_vec, _] : gate_lyt_windows)
                    {
                        assert(node_vec.size() > 1 && "The connected nodes vector cannot be singleton");

                        if (std::find(node_vec.cbegin(), node_vec.cend(), n) == node_vec.cend())
                        {
                            continue;
                        }

                        number_of_layouts++;
                    }

                    std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                    std::cout << fmt::format("Performing {} trials for {} sub-circuit{} for {} "
                                             "gate implementations",
                                             num_trials, number_of_layouts, number_of_layouts > 1 ? "s" : "",
                                             gate_designs.at(n).size())
                              << std::endl;

                    uint64_t min_bound = 0;                          // inclusive
                    uint64_t max_bound = gate_designs.at(n).size();  // exclusive

                    uint64_t repeated_attempt_number = 1;
                    uint64_t attempt_number          = 1;

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

                        if (discriminate_fitness_assessments(
                                global_pruning, selectivity, quantization_step, success_rate_ceiling, gate_designs[n],
                                gate_fitness_assessments[n], min_bound, max_bound, repeated_attempt_number,
                                attempt_number, completed_assessment[n]))
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
                    if (skip_physical_design_for_node(*stats.gate_layout, n) || gate_fitness_assessments.at(n).empty())
                    {
                        return;
                    }

                    std::vector<gate_design_t> selected_gate_implementations{};

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
    std::optional<CellLyt> exhaustively_enumerate_gate_design_combinations() const noexcept
    {
        std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;

        // operational_params.print = true;  // after comment: reenable function constness

        std::vector<uint64_t> indices(gate_designs.size(), 0);

        while (true)
        {
            CellLyt operational_circuit_candidate{circuit->skeleton};
            for (uint64_t i = 0; i < gate_designs.size(); i++)
            {
                const auto& [n, op_gate_designs_for_gate] = *std::next(gate_designs.cbegin(), static_cast<int64_t>(i));
                // select a random gate implementation for the tile that connects as input to n
                for (const uint64_t gate_design_cell_index : op_gate_designs_for_gate.at(indices.at(i)))
                {
                    operational_circuit_candidate.assign_cell_type(
                        circuit->all_canvas_positions.at(n).at(gate_design_cell_index),
                        sidb_technology::cell_type::LOGIC);
                }
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
                                                             const advanced_circuit_design_params<CellLyt>& params = {},
                                                             advanced_circuit_design_stats<GateLyt>* stats = nullptr)
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

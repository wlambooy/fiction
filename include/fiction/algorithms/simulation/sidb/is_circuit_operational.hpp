//
// Created by Jan Drewniok on 11.09.23.
//

#ifndef FICTION_IS_CIRCUIT_OPERATIONAL_HPP
#define FICTION_IS_CIRCUIT_OPERATIONAL_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/simulation/sidb/clustercomplete.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_result.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_bdl_circuit.hpp"
#include "fiction/technology/sidb_bdl_skeletons.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"

#include <kitty/bit_operations.hpp>
// #include <kitty/print.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiction
{

/**
 * Parameters for the `is_operational` algorithm. TODO
 */
struct is_circuit_operational_params
{
    /**
     * The termination condition for assessment of the operational status of the given layout.
     */
    enum class termination_condition : uint8_t
    {
        /**
         * The assessment for the given layout terminates either when it is found to be operational for all input
         * combinations, or an input combination is found for which the layout is not operational.
         */
        ON_FIRST_NON_OPERATIONAL,
        /**
         * The operational status is assessed for all input combinations.
         */
        ALL_INPUT_COMBINATIONS_ASSESSED
    };
    /**
     * Selector for the different ways to handle obtained simulation results.
     */
    enum class simulation_results_mode : uint8_t
    {
        /**
         * The simulation results for each input pattern are returned for operational gates.
         */
        KEEP_SIMULATION_RESULTS,
        /**
         * The simulation results are discarded after the operational status was assessed.
         */
        DISCARD_SIMULATION_RESULTS
    };
    /**
     * The simulation parameters for the physical simulation of the ground state.
     */
    sidb_simulation_parameters simulation_parameters{};
    /**
     * Parameters for the BDL input iterator.
     */
    bdl_input_iterator_params input_bdl_iterator_params{};
    /**
     * Condition to decide when to terminate the assessment of the operational status of the given layout.
     */
    termination_condition termination_cond = termination_condition::ON_FIRST_NON_OPERATIONAL;
    /**
     * Simulation results that are used to certify the status `OPERATIONAL` are not kept by default.
     */
    simulation_results_mode simulation_results_retention = simulation_results_mode::DISCARD_SIMULATION_RESULTS;
    bool                    print                        = false;
};

/**
 * This struct is used to collect results from the operational status assessment.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * todo
 */
template <typename Lyt, local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED>
struct circuit_operational_assessment
{
    /**
     * Typedef for the simulation results type, i.e., all physically valid charge distributions obtained during the
     * simulation.
     */
    using simulation_results_t = std::vector<charge_distribution_surface<Lyt, ExtPotType>>;
    /**
     * This struct collects the information for a specific input combination that was obtained during the assessment.
     */
    struct operational_assessment_for_input
    {
        /**
         * Standard constructor that only sets the operational status.
         *
         * @param op_status The operational status to set.
         */
        explicit operational_assessment_for_input(const operational_status op_status) noexcept : status{op_status} {}
        /**
         * The assessed operational status of the given layout under one input combination.
         */
        operational_status status;
        /**
         * todo
         */
        double logic_match = 1.0;
        /**
         * The charge distributions obtained for the input combination that was tested.
         */
        std::optional<simulation_results_t> simulation_results{};
    };
    /**
     * Standard constructor that only sets the operational status.
     *
     * @param op_status The operational status to set.
     */
    explicit circuit_operational_assessment(const operational_status op_status) noexcept : status{op_status} {}
    /**
     * The assessed operational status of the given layout. The status `OPERATIONAL` is given if and only the layout is
     * operational under all input combinations.
     */
    operational_status status;
    /**
     * When the termination condition is set to `ALL_INPUT_COMBINATIONS_ASSESSED`, the operational status for each
     * respective input combination is stored here, sorted by their binary representation. When the simulation retention
     * is set to `KEEP_SIMULATION_RESULTS`, this optional structure is also populated.
     */
    std::optional<std::vector<operational_assessment_for_input>> assessment_per_input{};
    /**
     * The number of input combinations tested.
     */
    std::size_t simulator_invocations{0};
    /**
     * Extracts the simulation results contained in this operational assessment through moves.
     *
     * @return A vector containing the simulation results for each respective input that was assessed.
     */
    std::vector<simulation_results_t> extract_simulation_results_per_input() const noexcept
    {
        assert(assessment_per_input.has_value() && "Assessment results per input are not present.");
        assert(!assessment_per_input.value().empty() && "No input combinations were assessed.");
        assert(assessment_per_input.value().front().simulation_results.has_value() &&
               "Simulation results were not retained during assessment.");

        std::vector<simulation_results_t> simulation_results_per_input{};
        simulation_results_per_input.reserve(assessment_per_input.value().size());

        for (operational_assessment_for_input assessment : assessment_per_input.value())
        {
            simulation_results_per_input.push_back(std::move(assessment.simulation_results.value()));
        }

        return simulation_results_per_input;
    }
};

namespace detail
{

struct thread_count_manager
{
    std::mutex mutex{};
    uint64_t count;

    explicit thread_count_manager(const uint64_t thread_count) : count{thread_count} {}

    uint64_t reserve_threads() noexcept
    {
        const std::lock_guard lock{mutex};

        const uint64_t reserved_threads = count;

        count = 0;

        return reserved_threads;
    }

    void return_threads(const uint64_t num) noexcept
    {
        if (num > 0)
        {
            const std::lock_guard lock{mutex};

            count += num;
        }
    }
};

/**
 * Implementation of the `is_operational` algorithm for a given SiDB layout.
 *
 * This class provides an implementation of the `is_operational` algorithm for
 * a specified SiDB layout and parameters. It checks whether the SiDB layout is operational
 * by simulating its behavior for different input combinations and comparing the results
 * to expected outputs from a truth table.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @tparam TT Type of the truth table.
 * @tparam todo.
 */
template <typename Lyt, typename GateLyt,
          local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED,
          typename SkeletonGateLibrary             = sidb_bdl_skeleton_1>
class is_circuit_operational_impl
{
  public:
    using operational_assessment_for_input =
        typename circuit_operational_assessment<Lyt, ExtPotType>::operational_assessment_for_input;
    /**
     * Constructor to initialize the algorithm with a layout and parameters.
     *
     * @param lyt The SiDB cell-level layout to be checked.
     * @param spec Expected Boolean function of the layout given as a multi-output truth table.
     * @param params Parameters for the `is_operational` algorithm.
     */
    is_circuit_operational_impl(
        const sidb_cell_level_bdl_circuit<Lyt, GateLyt, SkeletonGateLibrary>& implemented_bdl_circuit,
        const is_circuit_operational_params& params, const std::unique_ptr<thread_count_manager>& tcm) :
            implemented_circuit{implemented_bdl_circuit},
            parameters{params},
            thread_counter{tcm}
    {}
    /**
     * Run the `is_operational` algorithm.
     *
     * This function executes the operational status checking algorithm for the given SiDB layout
     * and parameters provided during initialization.
     *
     * @return Pair with the first element indicating the operational status (either `OPERATIONAL` or `NON_OPERATIONAL`)
     * and the second element indicating the reason if it is non-operational.
     */
    [[nodiscard]] circuit_operational_assessment<Lyt, ExtPotType> run() noexcept
    {
        circuit_operational_assessment<Lyt, ExtPotType> operational_assessment_results{operational_status::OPERATIONAL};

        // when `termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED` is set, the results of the operational status
        // assessment are also stored for each input separately
        std::vector<operational_assessment_for_input> assessment_results_per_input{};

        if (parameters.termination_cond ==
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED)
        {
            assessment_results_per_input.reserve(1 << implemented_circuit.circuit.num_inputs);
        }

        // when `simulation_results_mode::KEEP_SIMULATION_RESULTS` is set, the simulation results must be collected for
        // each input combination
        std::vector<typename circuit_operational_assessment<Lyt, ExtPotType>::simulation_results_t> sim_res_per_input{};
        if (parameters.simulation_results_retention ==
            is_circuit_operational_params::simulation_results_mode::KEEP_SIMULATION_RESULTS)
        {
            sim_res_per_input.reserve(1 << implemented_circuit.circuit.num_inputs);
        }

        bdl_input_iterator<Lyt> bii{implemented_circuit.cell_layout, parameters.input_bdl_iterator_params};
        bii = 0;

        // number of different input combinations
        for (auto i = 0u; i < 1 << implemented_circuit.circuit.num_inputs; ++i, ++bii)
        {
            operational_assessment_for_input assessment_results_for_this_input_combination{
                operational_status::OPERATIONAL};

            ++operational_assessment_results.simulator_invocations;

            // performs physical simulation of a given SiDB layout at a given input combination
            auto simulation_results = physical_simulation_of_layout(bii);

            if (!simulation_results.has_value())
            {
                continue;
            }

            // if no physically valid charge distributions were found, the layout is non-operational
            if (simulation_results->charge_distributions.empty())
            {
                operational_assessment_results.status = operational_status::NON_OPERATIONAL;

                if (parameters.termination_cond ==
                    is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL)
                {
                    return operational_assessment_results;
                }

                // all input combinations are being assessed

                assessment_results_for_this_input_combination.status = operational_status::NON_OPERATIONAL;

                assessment_results_for_this_input_combination.simulation_results.emplace();

                assessment_results_per_input.push_back(std::move(assessment_results_for_this_input_combination));

                continue;
            }

            auto   status      = operational_status::NON_OPERATIONAL;
            double logic_match = 0.0;

            if constexpr (ExtPotType == local_external_potential_type::BOUNDED)
            {
                std::vector<uint64_t> ground_state_in_environment_count(simulation_results->charge_distributions.size(),
                                                                        0);

                uint64_t total_environment_count = 0;

                if (simulation_results->charge_distributions.size() > 1)
                {
                    for (auto j = 0u; j < 1 << implemented_circuit.circuit.super_circuit.gate_layout.num_pis(); ++j)
                    {
                        if (!implemented_circuit.circuit.get_simulated_bdl_wires_for_input_indices(i, j)
                                 .get()
                                 .has_value())
                        {
                            continue;
                        }

                        total_environment_count++;

                        bool at_least_one_physically_valid_ground_state_for_super_circuit_input_combination = false;

                        std::vector<std::pair<double, uint64_t>> sort_energies{};
                        sort_energies.reserve(simulation_results->charge_distributions.size());

                        for (uint64_t cds_ix = 0; cds_ix < simulation_results->charge_distributions.size(); ++cds_ix)
                        {
                            sort_energies.emplace_back(
                                implemented_circuit.circuit
                                    .get_energy_of_expected_charge_distribution_with_sub_circuit_charge_distribution(
                                        i, j, simulation_results->charge_distributions.at(cds_ix),
                                        at_least_one_physically_valid_ground_state_for_super_circuit_input_combination),
                                cds_ix);
                        }

                        if (!at_least_one_physically_valid_ground_state_for_super_circuit_input_combination)
                        {
                            total_environment_count = 0;

                            break;
                        }

                        std::sort(sort_energies.begin(), sort_energies.end(),
                                  [&](const std::pair<double, uint64_t>& a, const std::pair<double, uint64_t>& b)
                                  { return a.first < b.first; });

                        ground_state_in_environment_count[sort_energies.front().second]++;
                    }
                }
                else
                {
                    ground_state_in_environment_count.push_back(1);
                    total_environment_count = 1;
                }

                if (total_environment_count != 0)
                {
                    for (uint64_t cds_ix = 0; cds_ix < simulation_results->charge_distributions.size(); ++cds_ix)
                    {
                        if (ground_state_in_environment_count.at(cds_ix) == 0)
                        {
                            continue;
                        }

                        const operational_assessment_for_input& op_assessment =
                            assess_logic_match_of_charge_distribution(
                                simulation_results->charge_distributions.at(cds_ix), i);

                        logic_match += static_cast<double>(ground_state_in_environment_count.at(cds_ix)) /
                                       static_cast<double>(total_environment_count) * op_assessment.logic_match;

                        if (op_assessment.status == operational_status::OPERATIONAL)
                        {
                            status = operational_status::OPERATIONAL;
                        }
                    }
                }
            }
            else
            {
                const auto ground_states = simulation_results->groundstates();

                status = operational_status::OPERATIONAL;

                double logic_match_sum = 0.0;

                for (const auto& gs : ground_states)
                {
                    const operational_assessment_for_input& op_assessment =
                        assess_logic_match_of_charge_distribution(gs, i);

                    logic_match_sum += op_assessment.logic_match;

                    if (op_assessment.status == operational_status::NON_OPERATIONAL)
                    {
                        status = operational_status::NON_OPERATIONAL;

                        // break;
                    }
                }

                logic_match = logic_match_sum / static_cast<double>(ground_states.size());
            }

            // std::cout << "final logic match: " << logic_match << std::endl;

            if (status == operational_status::NON_OPERATIONAL)
            {
                // the input combination is not operational

                operational_assessment_results.status = operational_status::NON_OPERATIONAL;

                if (parameters.termination_cond ==
                    is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL)
                {
                    return operational_assessment_results;
                }
            }

            if (parameters.termination_cond ==
                is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED)
            {
                assessment_results_for_this_input_combination.status      = status;
                assessment_results_for_this_input_combination.logic_match = logic_match;
            }

            if (parameters.print)
            {
                std::cout << "STATUS: " << (status == operational_status::NON_OPERATIONAL ? "NON-" : "")
                          << "OPERATIONAL" << std::endl;
                std::cout << fmt::format("logic match: {:.3f}", logic_match) << std::endl;
            }

            // store the assessment results for this input combination when the termination condition is set to
            // `termination_condition::ALL_INPUT_COMBINATION_ASSESSED` or the simulation result retention is set to
            // `simulation_results_mode::KEEP_SIMULATION_RESULTS`
            if (parameters.termination_cond ==
                    is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED ||
                parameters.simulation_results_retention ==
                    is_circuit_operational_params::simulation_results_mode::KEEP_SIMULATION_RESULTS)
            {
                // save simulation results when the simulation result retention is set to
                // `simulation_results_mode::KEEP_SIMULATION_RESULTS`
                if (parameters.simulation_results_retention ==
                    is_circuit_operational_params::simulation_results_mode::KEEP_SIMULATION_RESULTS)
                {
                    assessment_results_for_this_input_combination.simulation_results =
                        std::move(simulation_results->charge_distributions);
                }

                assessment_results_per_input.push_back(std::move(assessment_results_for_this_input_combination));
            }
        }

        // store the assessment results for all input combinations when the termination condition is set to
        // `termination_condition::ALL_INPUT_COMBINATION_ASSESSED` or the simulation result retention is set to
        // `simulation_results_mode::KEEP_SIMULATION_RESULTS`
        if (parameters.termination_cond ==
                is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED ||
            parameters.simulation_results_retention ==
                is_circuit_operational_params::simulation_results_mode::KEEP_SIMULATION_RESULTS)
        {
            operational_assessment_results.assessment_per_input = std::move(assessment_results_per_input);
        }

        // note: when all input combinations are assessed per termination condition, the assessment can yield the layout
        // is non-operational, yet we do not give a reason
        return operational_assessment_results;
    }
    /**
     * todo

     * @param given_cds The charge distribution surface to be checked for operation.
     * @param input_pattern Input pattern represented by the position of perturbers.
     * @return Pair with the first element indicating the operational status (either `OPERATIONAL` or `NON_OPERATIONAL`)
     * and the second element indicating the reason if it is non-operational.
     */
    [[nodiscard]] operational_assessment_for_input
    assess_logic_match_of_charge_distribution(const charge_distribution_surface<Lyt, ExtPotType>& given_cds,
                                              const uint64_t input_pattern) noexcept
    {
        const bool print_it = parameters.print;  // std::round(static_cast<double>(std::rand()) / (RAND_MAX + 1.0) *
                                                 // 10000) == 5000;
        if (print_it)
        {
            // std::cout << std::endl;
            // print_layout(given_cds);
            std::cout << "input_pattern: " << input_pattern << std::endl;
        }
        operational_assessment_for_input op_assessment{operational_status::NON_OPERATIONAL};

        uint64_t successful_bdl_pairs_count = 0;

        const auto count_logic_matching_bdl_pairs_in_wire_range =
            [&](const typename std::vector<bdl_pair<cell<Lyt>>>::const_iterator begin,
                const typename std::vector<bdl_pair<cell<Lyt>>>::const_iterator end, const bool signal) noexcept
        {
            for (auto it = begin; it != end; ++it)

            {
                assert(it->type != sidb_technology::cell_type::INPUT && "input BDL pairs are handled separately");

                successful_bdl_pairs_count += static_cast<uint64_t>((signal && encodes_bit_one(given_cds, *it)) ||
                                                                    (!signal && encodes_bit_zero(given_cds, *it)));

                if (print_it)
                {
                    std::cout << "bdl pair with " << it->upper.x << ',' << it->upper.y << " and " << it->lower.x << ','
                              << it->lower.y << ": "
                              << static_cast<uint64_t>((signal && encodes_bit_one(given_cds, *it)) ||
                                                       (!signal && encodes_bit_zero(given_cds, *it)))
                              << std::endl;
                }
            }
        };

        std::unordered_map<tile<GateLyt>, std::unordered_map<tile<GateLyt>, bool>> expected_signal_at_gate_connection{};

        uint64_t current_input_number = 0;

        for (uint64_t wire_ix = 0; wire_ix < implemented_circuit.circuit.bdl_wires.size(); ++wire_ix)
        {
            const bdl_wire<Lyt>& wire = implemented_circuit.circuit.bdl_wires.at(wire_ix);

            assert((wire.port.dir == port_direction::SOUTH || wire.port.dir == port_direction::EAST ||
                    wire.port.dir == port_direction::NONE) &&
                   "Wrong port direction; only row clocking is supported");

            const auto& [upper_tile, lower_tile] = implemented_circuit.circuit.gate_connections.at(wire_ix);

            typename std::vector<bdl_pair<cell<Lyt>>>::const_iterator successful_bdl_pairs_counting_start_it =
                wire.pairs.cbegin();

            if (implemented_circuit.circuit.gate_layout.is_pi_tile(upper_tile))
            {
                assert((expected_signal_at_gate_connection.count(lower_tile) == 0 ||
                        expected_signal_at_gate_connection.at(lower_tile).count(upper_tile) == 0) &&
                       "PI is visited twice");

                const bool current_bit_set =
                    (input_pattern &
                     (uint64_t{1ull} << (implemented_circuit.circuit.num_inputs - 1 - current_input_number++))) != 0ull;

                expected_signal_at_gate_connection[lower_tile].insert({upper_tile, current_bit_set});

                const bdl_pair<cell<Lyt>>& input_pair = wire.pairs.front();

                assert(input_pair.type == sidb_technology::cell_type::INPUT &&
                       "BDL wire connecting to a PI does not start with an input BDL pair");

                // successful_bdl_pairs_count += static_cast<uint64_t>(
                //     (current_bit_set && given_cds.get_charge_state(input_pair.lower) == sidb_charge_state::NEGATIVE)
                //     ||
                //     (!current_bit_set && given_cds.get_charge_state(input_pair.upper) ==
                //     sidb_charge_state::NEGATIVE));

                if (print_it)
                {
                    std::cout << "input: " << current_input_number - 1 << ": "
                              << static_cast<uint64_t>(
                                     (current_bit_set &&
                                      given_cds.get_charge_state(input_pair.lower) == sidb_charge_state::NEGATIVE) ||
                                     (!current_bit_set &&
                                      given_cds.get_charge_state(input_pair.upper) == sidb_charge_state::NEGATIVE))
                              << std::endl;
                }

                successful_bdl_pairs_counting_start_it = std::next(wire.pairs.cbegin(), 1);
            }

            assert(expected_signal_at_gate_connection.count(lower_tile) != 0 &&
                   expected_signal_at_gate_connection.at(lower_tile).count(upper_tile) != 0 &&
                   "Tile is visited before the incoming tile that connects it");

            const bool expected_signal_for_wire = expected_signal_at_gate_connection.at(lower_tile).at(upper_tile);

            count_logic_matching_bdl_pairs_in_wire_range(successful_bdl_pairs_counting_start_it, wire.pairs.cend(),
                                                         expected_signal_for_wire);

            const uint32_t num_inputs = implemented_circuit.circuit.gate_layout
                                            .node_function(implemented_circuit.circuit.gate_layout.get_node(lower_tile))
                                            .num_vars();

            assert(expected_signal_at_gate_connection.at(lower_tile).size() <= num_inputs &&
                   "Number of tiles visited connecting to the current tile exceeds the number of inputs to the node "
                   "function");

            if (implemented_circuit.circuit.gate_layout.is_po_tile(lower_tile) ||
                expected_signal_at_gate_connection.at(lower_tile).size() < num_inputs)
            {
                continue;
            }

            const std::vector<tile<GateLyt>>& outgoing_tiles =
                implemented_circuit.circuit.gate_layout.outgoing_data_flow(lower_tile);

            assert(!outgoing_tiles.empty() && "Non-PO tile does not have outgoing data flow");

            if constexpr (has_is_fanout_v<GateLyt>)
            {
                if (implemented_circuit.circuit.gate_layout.is_fanout(
                        implemented_circuit.circuit.gate_layout.get_node(lower_tile)))
                {
                    for (const tile<GateLyt>& lower_lower_t : outgoing_tiles)
                    {
                        expected_signal_at_gate_connection[lower_lower_t].insert(
                            {lower_tile, expected_signal_for_wire});
                    }

                    continue;
                }
            }

            assert(outgoing_tiles.size() == 1 && "Tile with single-output gate has more than one outgoing tile.");

            const tile<GateLyt>& first_input_tile = expected_signal_at_gate_connection.at(lower_tile).cbegin()->first;

            assert((num_inputs == 1 ||
                    (expected_signal_at_gate_connection.at(lower_tile).size() == 2 &&
                     first_input_tile.y ==
                         std::next(expected_signal_at_gate_connection.at(lower_tile).cbegin(), 1)->first.y &&
                     first_input_tile.x !=
                         std::next(expected_signal_at_gate_connection.at(lower_tile).cbegin(), 1)->first.x)) &&
                   "Only row clocking is supported; inputs to a tile must be on the same y with differing x");

            auto tt_inp = static_cast<uint8_t>(expected_signal_at_gate_connection.at(lower_tile).at(first_input_tile));

            if (num_inputs == 2)
            {
                const tile<GateLyt>& second_input_tile =
                    std::next(expected_signal_at_gate_connection.at(lower_tile).cbegin(), 1)->first;
                const auto tt_inp_second =
                    static_cast<uint8_t>(expected_signal_at_gate_connection.at(lower_tile).at(second_input_tile));

                // tt_inp <- 2 * L_in + R_in
                if (first_input_tile.x < second_input_tile.x)
                {
                    tt_inp = 2 * tt_inp + tt_inp_second;
                }
                else
                {
                    tt_inp += 2 * tt_inp_second;
                }
            }

            expected_signal_at_gate_connection[outgoing_tiles.front()].insert(
                {lower_tile, kitty::get_bit(implemented_circuit.circuit.gate_layout.node_function(
                                                implemented_circuit.circuit.gate_layout.get_node(lower_tile)),
                                            tt_inp)});
        }

        op_assessment.logic_match =
            static_cast<double>(successful_bdl_pairs_count) /
            static_cast<double>(implemented_circuit.circuit.num_bdl_pairs - implemented_circuit.circuit.num_inputs);

        if (successful_bdl_pairs_count ==
            implemented_circuit.circuit.num_bdl_pairs - implemented_circuit.circuit.num_inputs)
        {
            op_assessment.status = operational_status::OPERATIONAL;
        }

        if (print_it)
        {
            std::cout << "successful bdl pairs:" << successful_bdl_pairs_count << std::endl;
            std::cout << "total number of bdl pairs: " << implemented_circuit.circuit.num_bdl_pairs << std::endl;
            std::cout << fmt::format("logic match: {:.3f}\n", op_assessment.logic_match) << std::endl;
        }
        return op_assessment;
    }

  private:
    const sidb_cell_level_bdl_circuit<Lyt, GateLyt, SkeletonGateLibrary>& implemented_circuit{};
    /**
     * Parameters for the `is_operational` algorithm.
     */
    const is_circuit_operational_params&  parameters;
    const std::unique_ptr<thread_count_manager>& thread_counter;

    /**
     * This function conducts physical simulation of the given SiDB layout.
     * The simulation results are stored in the `sim_result` variable.
     *
     * @param bdl_iterator BDL input iterator representing the SiDB layout with a given input
     * combination.
     * @return Simulation results.
     */
    [[nodiscard]] std::optional<sidb_simulation_result<Lyt, ExtPotType>>
    physical_simulation_of_layout(const bdl_input_iterator<Lyt>& bdl_iterator) noexcept
    {
#if (FICTION_ALGLIB_ENABLED)
        // perform ClusterComplete exact simulation
        if constexpr (ExtPotType == local_external_potential_type::BOUNDED)
        {
            clustercomplete_params<cell<Lyt>, ExtPotType> cc_params{parameters.simulation_parameters};
            cc_params.available_threads = 1 + (thread_counter ? thread_counter->reserve_threads() : 0);

            Lyt cell_lyt{};

            std::unordered_map<cell<Lyt>, std::array<double, 2>> skeleton_influence_bounds_map{};

            (*bdl_iterator)
                .foreach_cell(
                    [&](const auto& c)
                    {
                        if (const auto ct = (*bdl_iterator).get_cell_type(c);
                            ct != sidb_technology::cell_type::OUTPUT_PERTURBER ||
                            implemented_circuit.circuit.super_circuit.skeleton.get_cell_type(c) ==
                                sidb_technology::cell_type::OUTPUT_PERTURBER)
                        {
                            cell_lyt.assign_cell_type(c, ct);

                            skeleton_influence_bounds_map.insert(
                                {c, std::array<double, 2>{std::numeric_limits<double>::infinity(),
                                                          -std::numeric_limits<double>::infinity()}});
                        }
                    });

            for (auto i = 0u; i < 1 << implemented_circuit.circuit.super_circuit.num_inputs; ++i)
            {
                const std::optional<charge_distribution_surface<Lyt>>& maybe_expected_charge_distribution_for_input =
                    implemented_circuit.circuit.get_simulated_bdl_wires_for_input_indices(
                        bdl_iterator.get_current_input_index(), i);

                if (!maybe_expected_charge_distribution_for_input.has_value())
                {
                    continue;
                }

                const charge_distribution_surface<Lyt>& expected_charge_distribution_for_input =
                    *maybe_expected_charge_distribution_for_input;

                for (auto& [c, bounds] : skeleton_influence_bounds_map)
                {
                    assert(expected_charge_distribution_for_input.get_local_internal_potential(c).has_value() &&
                           "c is not part of the layout");

                    const double pot_val_for_input =
                        *expected_charge_distribution_for_input.get_local_internal_potential(c);

                    bounds[0] = std::min(bounds[0], pot_val_for_input);
                    bounds[1] = std::max(bounds[1], pot_val_for_input);
                }
            }

            if (std::isinf(skeleton_influence_bounds_map.cbegin()->second.front()))
            {
                // input combination cannot exist in super circuit

                return std::nullopt;
            }

            for (auto& [c, bounds] : skeleton_influence_bounds_map)
            {
                cc_params.local_external_potential[c] = bounds;
            }

            const auto& sim_res = clustercomplete<Lyt, ExtPotType>(cell_lyt, cc_params);

            if (thread_counter)
            {
                thread_counter->return_threads(cc_params.available_threads - 1);
            }

            return sim_res;
        }
        else
        {
            if (parameters.print)
            {
                std::cout << "\nstarting exact simulation task (#SiDBs: " << (*bdl_iterator).num_cells() << ")"
                          << std::endl;
            }

            clustercomplete_params<cell<Lyt>> cc_params{parameters.simulation_parameters};
            const auto&                       res = clustercomplete(*bdl_iterator, cc_params);

            if (parameters.print)
            {
                std::cout << "exact simulation terminated in " << res.simulation_runtime.count() << " seconds\n\n";
                if (res.charge_distributions.empty())
                {
                    std::cout << "NO CHARGE DISTRIBUTIONS FOUNDS" << std::endl;
                    return res;
                }
                print_layout(res.groundstates().front());
                std::cout << std::endl;
            }

            return res;
        }

#else   // FICTION_ALGLIB_ENABLED
        assert(false && "ALGLIB must be enabled if ClusterComplete is to be used");
        return {}
#endif  // FICTION_ALGLIB_ENABLED
    }
    /**
     * This function returns `true` if `0` is encoded in the charge state of the given BDL pair. `false` otherwise.
     * Assumes row clocking.
     *
     * @param ground_state The ground state charge distribution surface.
     * @param bdl BDL pair to be evaluated.
     * @return `true` if `0` is encoded, `false` otherwise.
     */
    [[nodiscard]] bool encodes_bit_zero(const charge_distribution_surface<Lyt, ExtPotType>& ground_state,
                                        const bdl_pair<cell<Lyt>>&                          bdl) const noexcept
    {
        return static_cast<bool>((ground_state.get_charge_state(bdl.upper) == sidb_charge_state::NEGATIVE) &&
                                 (ground_state.get_charge_state(bdl.lower) == sidb_charge_state::NEUTRAL));
    }

    /**
     * This function returns `true` if `1` is encoded in the charge state of the given BDL pair. `false` otherwise.
     * Assumes row clocking.
     *
     * @param ground_state The ground state charge distribution surface.
     * @param bdl BDL pair to be evaluated.
     * @return `true` if `1` is encoded, `false` otherwise.
     */
    [[nodiscard]] bool encodes_bit_one(const charge_distribution_surface<Lyt, ExtPotType>& ground_state,
                                       const bdl_pair<cell<Lyt>>&                          bdl) const noexcept
    {
        return static_cast<bool>((ground_state.get_charge_state(bdl.upper) == sidb_charge_state::NEUTRAL) &&
                                 (ground_state.get_charge_state(bdl.lower) == sidb_charge_state::NEGATIVE));
    }
};

}  // namespace detail

/**
 * Determine the operational status of an SiDB layout.
 *
 * This function checks the operational status of a given SiDB layout using the `is_operational` algorithm. It
 * determines whether the SiDB layout is operational and returns the correct result for all \f$2^n\f$ input
 * combinations.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @tparam TT Type of the truth table.
 * @param lyt The SiDB cell-level layout to be checked.
 * @param spec Expected Boolean function of the layout given as a multi-output truth table.
 * @param params Parameters for the `is_operational` algorithm. todo
 * @return A datatype containing the operational status of the gate-level layout (either `OPERATIONAL` or
 * `NON_OPERATIONAL`) along with auxiliary statistics.
 */
template <typename Lyt, typename GateLyt,
          local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED,
          typename SkeletonGateLibrary>
[[nodiscard]] circuit_operational_assessment<Lyt, ExtPotType>
is_circuit_operational(const sidb_cell_level_bdl_circuit<Lyt, GateLyt, SkeletonGateLibrary>& implemented_circuit,
                       const is_circuit_operational_params& params = {}, const std::unique_ptr<detail::thread_count_manager>& tcm = nullptr) noexcept
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");

    assert(implemented_circuit.cell_layout.num_pis() > 0 && "lyt needs input cells");
    assert(implemented_circuit.cell_layout.num_pos() > 0 && "lyt needs output cells");

    assert(implemented_circuit.circuit.gate_layout.num_pis() * 2 == implemented_circuit.cell_layout.num_pis() &&
           "Each PI in the gate lyt needs to be implemented by a BDL pair");
    assert(implemented_circuit.circuit.gate_layout.num_pos() * 2 == implemented_circuit.cell_layout.num_pos() &&
           "Each PO in the gate lyt needs to be implemented by a BDL pair");

    detail::is_circuit_operational_impl<Lyt, GateLyt, ExtPotType, SkeletonGateLibrary> p{implemented_circuit, params,
                                                                                         tcm};

    const auto& assessment_result = p.run();

    if (params.print)
    {
        std::cout << "\n\nOVERALL STATUS: "
                  << (assessment_result.status == operational_status::NON_OPERATIONAL ? "NON-" : "") << "OPERATIONAL"
                  << std::endl;

        std::cout << std::endl;
    }

    return assessment_result;
}

}  // namespace fiction

#endif  // FICTION_IS_CIRCUIT_OPERATIONAL_HPP

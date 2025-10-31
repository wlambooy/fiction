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
#include <memory>
#include <optional>
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
        const is_circuit_operational_params&                                  params) :
            implemented_circuit{implemented_bdl_circuit},
            parameters{params}
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

        // perform quickcell pruning
        Lyt only_canvasses{};
        implemented_circuit.cell_layout.foreach_cell(
            [&](const cell<Lyt>& c)
            {
                if (implemented_circuit.cell_layout.get_cell_type(c) == sidb_technology::cell_type::LOGIC)
                {
                    only_canvasses.assign_cell_type(c, sidb_technology::cell_type::LOGIC);
                }
            });

        std::set<uint64_t> consistent_super_circuit_input_indices_set{};
        for (auto i = 0u; i < 1 << implemented_circuit.circuit.num_inputs; ++i)
        {
            if (implemented_circuit.circuit.input_index_possible_in_super_circuit(i))
            {
                const std::vector<uint64_t>& consistent_super_circuit_input_indices_for_input =
                    implemented_circuit.circuit.get_consistent_super_circuit_input_indices(i);

                consistent_super_circuit_input_indices_set.insert(
                    consistent_super_circuit_input_indices_for_input.cbegin(),
                    consistent_super_circuit_input_indices_for_input.cend());
            }
        }

        const std::vector<uint64_t> consistent_super_circuit_input_indices{
            consistent_super_circuit_input_indices_set.cbegin(), consistent_super_circuit_input_indices_set.cend()};

        for (uint64_t super_circuit_input_index_ix = 0;
             super_circuit_input_index_ix < consistent_super_circuit_input_indices.size();
             ++super_circuit_input_index_ix)
        {
            clustercomplete_params<cell<Lyt>, local_external_potential_type::BOUNDED> cc_params{
                parameters.simulation_parameters};
            cc_params.available_threads = 1;

            const auto collect_gate_influence_bounds = [&](const cell<Lyt>& c)
            {
                const mockturtle::node<GateLyt>& n = implemented_circuit.circuit.super_circuit.gate_layout.get_node(
                    implemented_circuit.circuit.super_circuit.skeleton_with_canvasses
                        .template get_cell_tile<tile<GateLyt>>(c));

                std::array<double, 2> gate_design_influence_bound_sum = {0, 0};

                implemented_circuit.circuit.super_circuit.gate_layout.foreach_node(
                    [&](const auto& super_circuit_n)
                    {
                        if (skip_physical_design_for_node(implemented_circuit.circuit.super_circuit.gate_layout,
                                                          super_circuit_n) ||
                            std::find(
                                implemented_circuit.circuit.tiles.cbegin(), implemented_circuit.circuit.tiles.cend(),
                                implemented_circuit.circuit.super_circuit.gate_layout.get_tile(super_circuit_n)) !=
                                implemented_circuit.circuit.tiles.cend())
                        {
                            return;
                        }

                        const std::array<double, 2>& influence_bounds_from_super_circuit_n =
                            implemented_circuit.circuit.super_circuit.get_gate_design_influence_bounds(
                                consistent_super_circuit_input_indices.at(super_circuit_input_index_ix), n, c,
                                super_circuit_n);

                        gate_design_influence_bound_sum[0] += influence_bounds_from_super_circuit_n[0];
                        gate_design_influence_bound_sum[1] += influence_bounds_from_super_circuit_n[1];
                    });

                return gate_design_influence_bound_sum;
            };

            only_canvasses.foreach_cell(
                [&](const auto& c)
                {
                    cc_params.local_external_potential[c] = collect_gate_influence_bounds(c);

                    const double skeleton_influence =
                        *implemented_circuit.circuit.super_circuit
                             .get_simulated_bdl_wires_for_input_index(
                                 consistent_super_circuit_input_indices.at(super_circuit_input_index_ix))
                             .get()
                             .get_local_internal_potential(c);

                    cc_params.local_external_potential[c][0] += skeleton_influence;
                    cc_params.local_external_potential[c][1] += skeleton_influence;
                });

            const auto sim_res =
                clustercomplete<Lyt, local_external_potential_type::BOUNDED>(only_canvasses, cc_params);

            if (sim_res.charge_distributions.empty())
            {
                // first pruning: physical infeasibility of canvas layouts

                operational_assessment_results.status = operational_status::NON_OPERATIONAL;

                return operational_assessment_results;
            }

            // second pruning: physical infeasibility of skeleton
            // NOTE: instead of the implementation below it might be better to only check the sub-circuit skeleton
            bdl_input_iterator<Lyt> bii{implemented_circuit.circuit.super_circuit.skeleton,
                                        parameters.input_bdl_iterator_params};
            bii = consistent_super_circuit_input_indices.at(super_circuit_input_index_ix);

            charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED> simulated_bdl_wires{
                implemented_circuit.circuit.super_circuit
                    .get_simulated_bdl_wires_for_input_index(
                        consistent_super_circuit_input_indices.at(super_circuit_input_index_ix))
                    .get()};  // todo: check clone behavior

            simulated_bdl_wires.template update_local_internal_potential<true>();  // todo do this earlier

            typename charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED>::
                local_external_potential_map_t& bounded_influence_from_canvasses =
                    simulated_bdl_wires.get_local_external_potentials_reference();

            // first collect bounded influence from other canvasses
            (*bii).foreach_cell([&](const auto& c)
                                { bounded_influence_from_canvasses[c] = collect_gate_influence_bounds(c); });

            bool at_least_one_physically_valid = false;

            // then collect influence for each simulated charge distribution of the sub-circuit canvasses
            for (uint64_t cds_ix = 0; cds_ix < sim_res.charge_distributions.size(); ++cds_ix)
            {
                const charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED>& cds =
                    sim_res.charge_distributions.at(cds_ix);

                std::unordered_map<cell<Lyt>, double> sub_circuit_canvas_influences{};
                sub_circuit_canvas_influences.reserve((*bii).num_cells());

                (*bii).foreach_cell(
                    [&](const auto& c)
                    {
                        double sub_circuit_canvas_influence_sum = 0.0;

                        cds.foreach_cell(
                            [&](const cell<Lyt>& canvas_c)
                            {
                                sub_circuit_canvas_influence_sum +=
                                    simulated_bdl_wires.get_chargeless_potential_between_sidbs(c, canvas_c) *
                                    charge_state_to_sign(cds.get_charge_state(canvas_c));
                            });

                        bounded_influence_from_canvasses[c][0] += sub_circuit_canvas_influence_sum;
                        bounded_influence_from_canvasses[c][1] += sub_circuit_canvas_influence_sum;

                        sub_circuit_canvas_influences[c] = sub_circuit_canvas_influence_sum;
                    });

                simulated_bdl_wires.update_local_external_potential();
                simulated_bdl_wires.determine_effective_charge_transition_thresholds();

                simulated_bdl_wires.template validity_check<true>();

                if (simulated_bdl_wires.is_physically_valid())
                {
                    // second pruning: the skeleton is not physically valid under any of the simulated sub-circuit
                    // canvas charge distributions

                    at_least_one_physically_valid = true;

                    break;
                }

                if (cds_ix < sim_res.charge_distributions.size() - 1)
                {
                    (*bii).foreach_cell(
                        [&](const auto& c)
                        {
                            bounded_influence_from_canvasses[c][0] -= sub_circuit_canvas_influences[c];
                            bounded_influence_from_canvasses[c][1] -= sub_circuit_canvas_influences[c];
                        });
                }
            }

            if (!at_least_one_physically_valid)
            {
                operational_assessment_results.status = operational_status::NON_OPERATIONAL;

                return operational_assessment_results;
            }
        }

        bdl_input_iterator<Lyt> bii{implemented_circuit.cell_layout, parameters.input_bdl_iterator_params};
        bii = 0;

        // number of different input combinations
        for (auto i = 0u; i < 1 << implemented_circuit.circuit.num_inputs; ++i, ++bii)
        {
            operational_assessment_for_input assessment_results_for_this_input_combination{
                operational_status::OPERATIONAL};

            // performs physical simulation of a given SiDB layout at a given input combination
            auto maybe_results = physical_simulation_of_layout(*bii, bii.get_current_input_index());

            if (!maybe_results.has_value())
            {
                continue;
            }

            ++operational_assessment_results.simulator_invocations;

            auto& [simulation_results, skeleton_influence_bounds] = *maybe_results;

            // if no physically valid charge distributions were found, the layout is non-operational
            if (simulation_results.charge_distributions.empty())
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

            const operational_assessment_for_input& op_assessment_for_input =
                determine_status_and_logic_match(i, std::move(skeleton_influence_bounds), simulation_results);

            if (op_assessment_for_input.status == operational_status::NON_OPERATIONAL)
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
                assessment_results_for_this_input_combination.status      = op_assessment_for_input.status;
                assessment_results_for_this_input_combination.logic_match = op_assessment_for_input.logic_match;
            }

            if (parameters.print)
            {
                std::cout << "STATUS: "
                          << (op_assessment_for_input.status == operational_status::NON_OPERATIONAL ? "NON-" : "")
                          << "OPERATIONAL" << std::endl;
                std::cout << fmt::format("logic match: {:.3f}", op_assessment_for_input.logic_match) << std::endl;
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
                        std::move(simulation_results.charge_distributions);
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

  private:
    const sidb_cell_level_bdl_circuit<Lyt, GateLyt, SkeletonGateLibrary>& implemented_circuit{};
    /**
     * Parameters for the `is_operational` algorithm.
     */
    const is_circuit_operational_params& parameters;

    using pair_t = std::pair<sidb_simulation_result<Lyt, ExtPotType>,
                             typename charge_distribution_surface<
                                 Lyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t>;
    /**
     * This function conducts physical simulation of the given SiDB layout.
     * The simulation results are stored in the `sim_result` variable.
     *
     * @param todo
     * @return Simulation results.
     */
    template <bool consider_internal_skeleton = false>
    [[nodiscard]] std::optional<pair_t> physical_simulation_of_layout(const Lyt&     lyt,
                                                                      const uint64_t input_index) noexcept
    {
#if (FICTION_ALGLIB_ENABLED)
        // perform ClusterComplete exact simulation
        if constexpr (ExtPotType == local_external_potential_type::BOUNDED)
        {
            if (!implemented_circuit.circuit.input_index_possible_in_super_circuit(input_index))
            {
                // input combination cannot exist in super circuit

                return std::nullopt;
            }

            clustercomplete_params<cell<Lyt>, ExtPotType> cc_params{parameters.simulation_parameters};

            const typename charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED>::
                local_external_potential_map_t& skeleton_influence_bounds =
                    implemented_circuit.circuit.template collect_influence_bounds<consider_internal_skeleton>(
                        lyt, input_index, cc_params.local_external_potential);

            Lyt cell_lyt_without_internal_output_perturber{};

            lyt.foreach_cell(
                [&](const auto& c)
                {
                    if (const auto& ct = implemented_circuit.circuit.is_not_internal_output_perturber(lyt, c);
                        ct.has_value())
                    {
                        cell_lyt_without_internal_output_perturber.assign_cell_type(c, *ct);
                    }
                });

            cc_params.available_threads = 1;

            const auto& sim_res =
                clustercomplete<Lyt, ExtPotType>(cell_lyt_without_internal_output_perturber, cc_params);

            return std::make_pair(std::move(sim_res), std::move(skeleton_influence_bounds));
        }
        else
        {
            if (parameters.print)
            {
                std::cout << "\nstarting exact simulation task (#SiDBs: " << lyt.num_cells() << ")" << std::endl;
            }

            clustercomplete_params<cell<Lyt>> cc_params{parameters.simulation_parameters};
            const auto&                       res = clustercomplete(lyt, cc_params);

            if (parameters.print)
            {
                std::cout << "exact simulation terminated in " << res.simulation_runtime.count() << " seconds\n\n";
                if (res.charge_distributions.empty())
                {
                    std::cout << "NO CHARGE DISTRIBUTIONS FOUNDS" << std::endl;
                    return std::make_optional<pair_t>(
                        {std::move(res),
                         typename charge_distribution_surface<
                             Lyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t{}});
                }
                print_layout(res.groundstates().front());
                std::cout << std::endl;
            }

            return std::make_optional<pair_t>(
                {std::move(res), typename charge_distribution_surface<
                                     Lyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t{}});
        }

#else   // FICTION_ALGLIB_ENABLED
        assert(false && "ALGLIB must be enabled if ClusterComplete is to be used");
        return {}
#endif  // FICTION_ALGLIB_ENABLED
    }

    [[nodiscard]] std::optional<std::vector<std::array<double, 2>>> get_energy_bounds_per_super_circuit_input(
        const uint64_t sub_circuit_input_index,
        const typename charge_distribution_surface<
            Lyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t& skeleton_influence_bounds,
        charge_distribution_surface<Lyt, ExtPotType>&                                     cds) noexcept
    {
        std::vector<std::array<double, 2>> energy_bounds_per_super_circuit_input{};
        energy_bounds_per_super_circuit_input.reserve(
            implemented_circuit.circuit.get_consistent_super_circuit_input_indices(sub_circuit_input_index)
                .get()
                .size());

        typename charge_distribution_surface<
            Lyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t& influence_bounds =
            cds.get_local_external_potentials_reference();
        cds.foreach_cell(
            [&](const auto& c)
            {
                if (implemented_circuit.circuit.is_not_internal_output_perturber(cds, c))
                {
                    influence_bounds[c][0] -= skeleton_influence_bounds.at(c)[0];
                    influence_bounds[c][1] -= skeleton_influence_bounds.at(c)[1];
                }
            });

        for (const uint64_t i :
             implemented_circuit.circuit.get_consistent_super_circuit_input_indices(sub_circuit_input_index).get())
        {
            cds.foreach_cell(
                [&](const auto& c)
                {
                    if (implemented_circuit.circuit.is_not_internal_output_perturber(cds, c))
                    {
                        const double skeleton_influence =
                            implemented_circuit.circuit.get_skeleton_influence(sub_circuit_input_index, i, c);

                        influence_bounds[c][0] += skeleton_influence;
                        influence_bounds[c][1] += skeleton_influence;
                    }
                });

            cds.recompute_electrostatic_potential_energy();

            cds.determine_effective_charge_transition_thresholds();
            cds.validity_check();

            if (!cds.is_physically_valid())
            {
                return std::nullopt;
            }

            energy_bounds_per_super_circuit_input.push_back(cds.get_electrostatic_potential_energy());

            if (i == implemented_circuit.circuit.get_consistent_super_circuit_input_indices(sub_circuit_input_index)
                         .get()
                         .back())
            {
                break;
            }

            cds.foreach_cell(
                [&](const auto& c)
                {
                    if (implemented_circuit.circuit.is_not_internal_output_perturber(cds, c))
                    {
                        const double skeleton_influence =
                            implemented_circuit.circuit.get_skeleton_influence(sub_circuit_input_index, i, c);

                        influence_bounds[c][0] -= skeleton_influence;
                        influence_bounds[c][1] -= skeleton_influence;
                    }
                });
        }

        return std::make_optional(std::move(energy_bounds_per_super_circuit_input));
    }

    [[nodiscard]] std::optional<std::vector<std::vector<uint64_t>>>
    get_possible_ground_state_indices_per_super_circuit_input(
        const uint64_t                                                   sub_circuit_input_index,
        const std::vector<charge_distribution_surface<Lyt, ExtPotType>>& simulated_charge_distributions,
        std::vector<std::optional<std::vector<std::array<double, 2>>>>&& energy_bounds) noexcept
    {
        std::vector<std::vector<uint64_t>> possible_ground_state_indices_per_super_circuit_input{};
        possible_ground_state_indices_per_super_circuit_input.reserve(
            implemented_circuit.circuit.get_consistent_super_circuit_input_indices(sub_circuit_input_index)
                .get()
                .size());

        for (uint64_t super_circuit_input_index_ix = 0;
             super_circuit_input_index_ix <
             implemented_circuit.circuit.get_consistent_super_circuit_input_indices(sub_circuit_input_index)
                 .get()
                 .size();
             ++super_circuit_input_index_ix)
        {
            std::vector<uint64_t> sorted_indices{};
            sorted_indices.reserve(simulated_charge_distributions.size());

            for (uint64_t cds_ix = 0; cds_ix < simulated_charge_distributions.size(); ++cds_ix)
            {
                if (energy_bounds.at(cds_ix))
                {
                    sorted_indices.push_back(cds_ix);
                }
            }

            if (sorted_indices.empty())
            {
                return std::nullopt;
            }

            // Sort by energy[0] ascending, then energy[1] descending
            std::sort(sorted_indices.begin(), sorted_indices.end(),
                      [&](const auto& a, const auto& b)
                      {
                          const std::array<double, 2>& e1 = energy_bounds.at(a)->at(super_circuit_input_index_ix);
                          const std::array<double, 2>& e2 = energy_bounds.at(b)->at(super_circuit_input_index_ix);

                          if (std::abs(e1[0] - e2[0]) < constants::ERROR_MARGIN)
                          {
                              return e1[1] > e2[1];
                          }

                          return e1[0] < e2[0];
                      });

            std::vector possible_ground_state_indices{sorted_indices.front()};
            possible_ground_state_indices.reserve(sorted_indices.size());

            double max_energy = energy_bounds.at(sorted_indices.front())->at(super_circuit_input_index_ix)[1];

            for (uint64_t ix = 1; ix < sorted_indices.size(); ++ix)
            {
                const std::array<double, 2>& bounded_energy =
                    energy_bounds.at(sorted_indices.at(ix))->at(super_circuit_input_index_ix);

                if (bounded_energy[0] > max_energy - constants::ERROR_MARGIN)
                {
                    break;  // Done: all future entries start after max_energy
                }

                possible_ground_state_indices.push_back(sorted_indices.at(ix));

                max_energy = std::max(max_energy, bounded_energy[1]);
            }

            possible_ground_state_indices_per_super_circuit_input.push_back(std::move(possible_ground_state_indices));
        }

        return possible_ground_state_indices_per_super_circuit_input;
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
        const charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED>& simulated_bdl_wires =
            implemented_circuit.circuit.super_circuit
                .get_simulated_bdl_wires_for_input_index(
                    implemented_circuit.circuit.get_consistent_super_circuit_input_indices(input_pattern).get().front())
                .get();

        for (const bdl_wire<Lyt>& wire : implemented_circuit.circuit.bdl_wires)
        {
            for (const bdl_pair<cell<Lyt>>& pair : wire.pairs)
            {
                const sidb_charge_state upper_cs = given_cds.get_charge_state(pair.upper);
                const sidb_charge_state lower_cs = given_cds.get_charge_state(pair.lower);

                if (upper_cs != sidb_charge_state::NONE)
                {
                    if (upper_cs != simulated_bdl_wires.get_charge_state(pair.upper))
                    {
                        return operational_assessment_for_input{operational_status::NON_OPERATIONAL};
                    }
                }

                if (lower_cs != sidb_charge_state::NONE)
                {
                    if (lower_cs != simulated_bdl_wires.get_charge_state(pair.lower))
                    {
                        return operational_assessment_for_input{operational_status::NON_OPERATIONAL};
                    }
                }
            }
        }

        return operational_assessment_for_input{operational_status::OPERATIONAL};
    }

    [[nodiscard]] operational_assessment_for_input determine_status_and_logic_match(
        const uint64_t i,
        typename charge_distribution_surface<
            Lyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t&& skeleton_influence_bounds,
        sidb_simulation_result<Lyt, ExtPotType>&                                           simulation_results) noexcept
    {
        operational_assessment_for_input op_ass{operational_status::OPERATIONAL};

        if constexpr (ExtPotType == local_external_potential_type::BOUNDED)
        {
            std::vector<std::optional<std::vector<std::array<double, 2>>>> energy_bounds{};
            energy_bounds.reserve(simulation_results.charge_distributions.size());

            for (charge_distribution_surface<Lyt, ExtPotType>& cds : simulation_results.charge_distributions)
            {
                energy_bounds.push_back(get_energy_bounds_per_super_circuit_input(i, skeleton_influence_bounds, cds));
            }

            const std::optional<std::vector<std::vector<uint64_t>>>& maybe_possible_ground_state_indices =
                get_possible_ground_state_indices_per_super_circuit_input(i, simulation_results.charge_distributions,
                                                                          std::move(energy_bounds));

            if (!maybe_possible_ground_state_indices.has_value())
            {
                op_ass.status = operational_status::NON_OPERATIONAL;

                return op_ass;
            }

            for (const std::vector<uint64_t>& possible_ground_state_indices : *maybe_possible_ground_state_indices)
            {
                bool at_least_one_operational = false;

                for (const uint64_t cds_ix : possible_ground_state_indices)
                {
                    if (assess_logic_match_of_charge_distribution(simulation_results.charge_distributions.at(cds_ix), i)
                            .status == operational_status::OPERATIONAL)
                    {
                        at_least_one_operational = true;

                        break;
                    }
                }

                if (!at_least_one_operational)
                {
                    op_ass.status = operational_status::NON_OPERATIONAL;

                    break;
                }
            }
        }
        else
        {
            const auto ground_states = simulation_results.groundstates();

            double logic_match_sum = 0.0;

            for (const auto& gs : ground_states)
            {
                const operational_assessment_for_input& op_assessment =
                    assess_logic_match_of_charge_distribution(gs, i);

                logic_match_sum += op_assessment.logic_match;

                if (op_assessment.status == operational_status::NON_OPERATIONAL)
                {
                    op_ass.status = operational_status::NON_OPERATIONAL;

                    // break;
                }
            }

            op_ass.logic_match = logic_match_sum / static_cast<double>(ground_states.size());
        }

        return op_ass;
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
                       const is_circuit_operational_params&                                  params = {}) noexcept
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

    detail::is_circuit_operational_impl<Lyt, GateLyt, ExtPotType, SkeletonGateLibrary> p{implemented_circuit, params};

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

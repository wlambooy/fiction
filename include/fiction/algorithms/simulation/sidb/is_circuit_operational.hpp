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

#include <fiction/io/print_layout.hpp>

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

// #define  PRINT_DEBUG

namespace detail
{

struct thread_count_manager
{
    std::mutex mutex{};
    uint64_t   count;

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
          local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED>
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
    is_circuit_operational_impl(const sidb_cell_level_bdl_circuit<Lyt, GateLyt>& implemented_bdl_circuit,
                                const is_circuit_operational_params&             params,
                                const std::unique_ptr<thread_count_manager>&     tcm) :
            implemented_circuit{implemented_bdl_circuit},
            parameters{params},
            thread_counter{tcm},
            inp_pairs{detect_bdl_pairs<Lyt>(implemented_circuit.cell_layout, sidb_technology::cell_type::INPUT,
                                            detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance})}

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
    [[nodiscard]] operational_status run() noexcept
    {
        /* // perform quickcell pruning
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
        for (auto i = 0u; i < 1 << implemented_circuit.circuit.gate_layout.num_pis(); ++i)
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

        for (const uint64_t super_circuit_input_index : consistent_super_circuit_input_indices)
        {
            clustercomplete_params<cell<Lyt>, local_external_potential_type::BOUNDED> cc_params{
                parameters.simulation_parameters};
            cc_params.available_threads = 1;

            const auto collect_gate_influence_bounds = [&](const cell<Lyt>& c)
            {
                const mockturtle::node<GateLyt>& n = implemented_circuit.circuit.super_circuit.gate_layout.get_node(
                    implemented_circuit.circuit.super_circuit.canvasses_lyt
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
                                super_circuit_input_index, n, c, super_circuit_n);

                        gate_design_influence_bound_sum[0] += influence_bounds_from_super_circuit_n[0];
                        gate_design_influence_bound_sum[1] += influence_bounds_from_super_circuit_n[1];
                    });

                return gate_design_influence_bound_sum;
            };

            only_canvasses.foreach_cell(
                [&](const auto& c)
                {
                    cc_params.local_external_potential[c] = collect_gate_influence_bounds(c);

                    // const double skeleton_influence =
                    //     *implemented_circuit.circuit.super_circuit
                    //          .get_simulated_bdl_wires_for_input_index(super_circuit_input_index)
                    //          .get()
                    //          .get_local_internal_potential(c);
                    //
                    // cc_params.local_external_potential[c][0] += skeleton_influence;
                    // cc_params.local_external_potential[c][1] += skeleton_influence;
                });

            const auto sim_res =
                clustercomplete<Lyt, local_external_potential_type::BOUNDED>(only_canvasses, cc_params);

            if (sim_res.charge_distributions.empty())
            {
                // first pruning: physical infeasibility of canvas layouts

                return operational_status::NON_OPERATIONAL;
            }

            // second pruning: physical infeasibility of skeleton
            // NOTE: instead of the implementation below it might be better to only check the sub-circuit skeleton
            bdl_input_iterator<Lyt> bii{implemented_circuit.circuit.super_circuit.canvasses_lyt,
                                        parameters.input_bdl_iterator_params};
            bii = super_circuit_input_index;

            // charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED> simulated_bdl_wires{
            //     implemented_circuit.circuit.super_circuit
            //         .get_simulated_bdl_wires_for_input_index(super_circuit_input_index)
            //         .get()};  // todo: check clone behavior
            //
            typename charge_distribution_surface<Lyt, local_external_potential_type::BOUNDED>::
                local_external_potential_map_t& bounded_influence_from_canvasses{};
                    // simulated_bdl_wires.get_local_external_potentials_reference();

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
                return operational_status::NON_OPERATIONAL;
            }
        }*/

        std::vector<bool> flip_chart{};
        flip_chart.reserve(inp_pairs.size());

        for (const auto& pair : inp_pairs)
        {
            bool flip = false;

            if (pair.upper.y == pair.lower.y)
            {
                const tile<GateLyt>& t =
                    implemented_circuit.circuit.super_circuit.canvasses_lyt.template get_cell_tile<tile<GateLyt>>(
                        pair.upper);

                for (int8_t x = -2; x < 0; ++x)  // todo make parametric?
                {
                    if (implemented_circuit.circuit.super_circuit.gate_layout.get_node(
                            tile<GateLyt>{t.x + x, t.y + 1}) != 0)
                    {
                        flip = true;
                        break;
                    }
                }
            }

            flip_chart.push_back(flip);
        }

        Lyt lyt{implemented_circuit.cell_layout.clone()};

        // number of different input combinations
        for (auto iix = 0u; iix < 1 << inp_pairs.size(); ++iix)
        {
#ifdef PRINT_DEBUG
            std::cout << "input index: " << iix << std::endl;
#endif
            if (!implemented_circuit.circuit.input_index_possible_in_super_circuit(iix))
            {
                // input combination cannot exist in super circuit
#ifdef PRINT_DEBUG
                std::cout << "SKIP" << std::endl;
#endif
                continue;
            }

            for (uint64_t input_number = inp_pairs.size() - 1; input_number < inp_pairs.size(); --input_number)
            {
                if (flip_chart.at(input_number) ^
                    ((iix & (uint64_t{1ull} << (inp_pairs.size() - 1 - input_number))) != 0ull))
                {
                    lyt.assign_cell_type(inp_pairs.at(input_number).lower, technology<Lyt>::cell_type::INPUT);
                    lyt.assign_cell_type(inp_pairs.at(input_number).upper, technology<Lyt>::cell_type::EMPTY);
                }
                else
                {
                    lyt.assign_cell_type(inp_pairs.at(input_number).lower, technology<Lyt>::cell_type::EMPTY);
                    lyt.assign_cell_type(inp_pairs.at(input_number).upper, technology<Lyt>::cell_type::INPUT);
                }
            }

            for (const uint64_t super_circuit_input_index :
                 implemented_circuit.circuit.get_consistent_super_circuit_input_indices(iix).get())
            {
                // performs physical simulation of a given SiDB layout at a given input combination

#ifdef PRINT_DEBUG
                std::cout << "SUPER CIRCUIT INPUT INDEX: " << super_circuit_input_index << std::endl;
                print_layout(lyt);
                std::cout << "to simulate ^^^^^^^^^" << std::endl;
#endif
                if (const auto& sim_res = physical_simulation_of_layout(lyt, iix, super_circuit_input_index);
                    determine_status_and_logic_match(iix, super_circuit_input_index, sim_res) ==
                    operational_status::NON_OPERATIONAL)
                {
                    return operational_status::NON_OPERATIONAL;
                }
            }
        }

        return operational_status::OPERATIONAL;
    }

  private:
    const sidb_cell_level_bdl_circuit<Lyt, GateLyt>& implemented_circuit{};
    /**
     * Parameters for the `is_operational` algorithm.
     */
    const is_circuit_operational_params& parameters;

    const std::unique_ptr<thread_count_manager>& thread_counter;

    const std::vector<bdl_pair<cell<Lyt>>> inp_pairs{};

    using pair_t = std::pair<sidb_simulation_result<Lyt, ExtPotType>,
                             typename charge_distribution_surface<Lyt>::local_external_potential_map_t>;
    /**
     * This function conducts physical simulation of the given SiDB layout.
     * The simulation results are stored in the `sim_result` variable.
     *
     * @param todo
     * @return Simulation results.
     */
    template <bool consider_internal_skeleton = false>
    [[nodiscard]] sidb_simulation_result<Lyt, ExtPotType>
    physical_simulation_of_layout(const Lyt& lyt, const uint64_t input_index,
                                  const uint64_t super_circuit_input_index) noexcept
    {
#if (FICTION_ALGLIB_ENABLED)
        // perform ClusterComplete exact simulation
        if constexpr (ExtPotType == local_external_potential_type::BOUNDED)
        {
            clustercomplete_params<cell<Lyt>, ExtPotType> cc_params{parameters.simulation_parameters};

            implemented_circuit.circuit.collect_influence_bounds(lyt, input_index, super_circuit_input_index,
                                                                 cc_params.local_external_potential);

            // Lyt cell_lyt_without_internal_perturbers{};
            //
            // lyt.foreach_cell(
            //     [&](const auto& c)
            //     {
            //         if (const auto& ct = implemented_circuit.circuit.is_not_internal_perturber(lyt, c);
            //         ct.has_value())
            //         {
            //             cell_lyt_without_internal_perturbers.assign_cell_type(c, *ct);
            //         }
            //     });

            cc_params.available_threads = 1;

            if (parameters.print)
            {
                for (const auto& [c, b] : cc_params.local_external_potential)
                {
                    std::cout << c.x << ' ' << c.y << " : " << b[0] << ' ' << b[1] << std::endl;
                }
                print_layout(lyt);
            }

            const auto& sim_res = clustercomplete<Lyt, ExtPotType>(lyt, cc_params);

            if (parameters.print)
            {
                std::cout << "RES START" << std::endl;
                for (const auto& c : sim_res.charge_distributions)
                {
                    print_layout(c);
                }
                std::cout << "RES END" << std::endl;
            }

            return sim_res;
        }
        else
        {
            if (parameters.print)
            {
                std::cout << "\nstarting exact simulation task (#SiDBs: " << lyt.num_cells() << ")" << std::endl;
            }

            clustercomplete_params<cell<Lyt>> cc_params{parameters.simulation_parameters};

            cc_params.available_threads = 1 + (thread_counter ? thread_counter->reserve_threads() : 0);

            const auto& res = clustercomplete(lyt, cc_params);

            if (thread_counter)
            {
                thread_counter->return_threads(cc_params.available_threads - 1);
            }

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

    [[nodiscard]] std::vector<uint64_t> get_possible_ground_state_indices(
        const std::vector<charge_distribution_surface<Lyt, ExtPotType>>& simulated_charge_distributions) noexcept
    {
        std::vector<uint64_t> possible_ground_state_indices{};
        possible_ground_state_indices.reserve(simulated_charge_distributions.size());

        for (uint64_t ix = 0; ix < simulated_charge_distributions.size(); ++ix)
        {
            possible_ground_state_indices.push_back(ix);
        }

        return possible_ground_state_indices;

        // todo mistake in exact reasoning! this considers energy bounds on the subsystem, but should reason on energy
        // todo   bounds of the whole system.... should disabled for now.

        //         std::vector<std::pair<uint64_t, std::array<double, 2>>> sorted_indices{};
        //         sorted_indices.reserve(simulated_charge_distributions.size());
        //
        //         for (uint64_t cds_ix = 0; cds_ix < simulated_charge_distributions.size(); ++cds_ix)
        //         {
        //             sorted_indices.emplace_back(cds_ix,
        //                                         simulated_charge_distributions.at(cds_ix).get_electrostatic_potential_energy());
        //         }
        //
        //         // Sort by energy[0] ascending, then energy[1] descending
        //         std::sort(sorted_indices.begin(), sorted_indices.end(),
        //                   [&](const auto& a, const auto& b)
        //                   {
        //                       const std::array<double, 2>& e1 = a.second;
        //                       const std::array<double, 2>& e2 = b.second;
        //
        //                       if (std::abs(e1[0] - e2[0]) < constants::ERROR_MARGIN)
        //                       {
        //                           return e1[1] > e2[1];
        //                       }
        //
        //                       return e1[0] < e2[0];
        //                   });
        //
        //         std::vector possible_ground_state_indices{sorted_indices.front().first};
        //
        //         // todo ...
        //
        //         uint64_t i = 1;
        //         while (i < std::min(decltype(sorted_indices.size()){4}, sorted_indices.size()))
        //         {
        // //             bool found_pos = false;
        // //             bool found_neg = false;
        // //             for (const auto& db :
        // simulated_charge_distributions.at(sorted_indices.at(i).first).get_sidb_order())
        // //             {
        // //                 switch (
        // simulated_charge_distributions.at(sorted_indices.at(i).first).get_charge_state(db))
        // //                 {
        // //                     case sidb_charge_state::POSITIVE: found_pos = true; break;
        // //                     case sidb_charge_state::NEGATIVE: found_neg = true; break;
        // //                     default: break;
        // //                 }
        // //                 if (found_neg || found_pos)
        // // break;
        // //             }
        // //             if (found_pos || !found_neg)
        // //             {
        //                 possible_ground_state_indices.emplace_back(sorted_indices.at(i++).first);
        //             // }
        //             // else
        //             // {
        //             //     break;
        //             // }
        //         }
        //
        //         // todo very not exact!!!
        //         // possible_ground_state_indices.reserve(sorted_indices.size());
        //         //
        //         // double max_energy = sorted_indices.front().second[1];
        //         //
        //         // for (uint64_t ix = 1; ix < sorted_indices.size(); ++ix)
        //         // {
        //         //     const std::array<double, 2>& bounded_energy = sorted_indices.at(ix).second;
        //         //
        //         //     if (bounded_energy[0] > max_energy - constants::ERROR_MARGIN)
        //         //     {
        //         //         break;  // Done: all future entries start after max_energy
        //         //     }
        //         //
        //         //     possible_ground_state_indices.push_back(sorted_indices.at(ix).first);
        //         //
        //         //     max_energy = std::max(max_energy, bounded_energy[1]);
        //         // }
        //
        //         return possible_ground_state_indices;
    }
    /**
    * todo

    * @param given_cds The charge distribution surface to be checked for operation.
    * @param input_pattern Input pattern represented by the position of perturbers.
    * @return Pair with the first element indicating the operational status (either `OPERATIONAL` or `NON_OPERATIONAL`)
    * and the second element indicating the reason if it is non-operational.
    */
    [[nodiscard]] operational_status
    assess_logic_match_of_charge_distribution(const charge_distribution_surface<Lyt, ExtPotType>& given_cds,
                                              const uint64_t                                      input_pattern,
                                              const uint64_t super_circuit_input_pattern) noexcept
    {
        const auto is_input_bit_set = [&](const uint64_t input_number)
        {
            return (super_circuit_input_pattern & (uint64_t{1ull} << (2 - 1 - input_number))) != 0ull;
        };  // todo hardcoded 2 input here

        std::vector<bdl_pair<cell<Lyt>>> input_bdl_pairs = inp_pairs;

#ifdef PRINT_DEBUG
        std::cout << "layout:" << std::endl;
        print_layout(given_cds);
        std::cout << "inputs: " << is_input_bit_set(0) << '\t' << is_input_bit_set(1) << std::endl;
        std::cout << "input pattern: " << input_pattern << std::endl;
        std::cout << "super circuit input pattern: " << super_circuit_input_pattern << std::endl;
        std::cout << "input pairs: " << input_bdl_pairs.size() << std::endl;
#endif

        // todo: this is very barebones
        for (auto&& p : detect_bdl_pairs(implemented_circuit.cell_layout, sidb_technology::cell_type::NORMAL,
                                         detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance}))
        {
            const tile<GateLyt>& t =
                implemented_circuit.circuit.super_circuit.canvasses_lyt.template get_cell_tile<tile<GateLyt>>(p.upper);

            if (t.y == 0)
            {
                continue;
            }

            for (int8_t x = -3; x < 3; ++x)  // todo make parametric?
            {
                if (implemented_circuit.circuit.super_circuit.gate_layout.is_pi_tile(tile<GateLyt>{t.x + x, t.y - 1}))
                {
                    input_bdl_pairs.push_back(std::move(p));

                    break;
                }
            }
        }

#ifdef PRINT_DEBUG
        std::cout << "augmented input pairs: " << input_bdl_pairs.size() << std::endl;

        std::cout << "bounding min: " << bounding_box_2d<Lyt>{given_cds}.get_min().x << ", "
                  << bounding_box_2d<Lyt>{given_cds}.get_min().y << std::endl;
        std::cout << "bounding max: " << bounding_box_2d<Lyt>{given_cds}.get_max().x << ", "
                  << bounding_box_2d<Lyt>{given_cds}.get_max().y << std::endl;
#endif

        for (uint8_t i = 0; i < input_bdl_pairs.size(); ++i)
        {
            bool south_east = true;

            const tile<GateLyt>& t =
                implemented_circuit.circuit.super_circuit.canvasses_lyt.template get_cell_tile<tile<GateLyt>>(
                    input_bdl_pairs.at(i).upper);

            for (int8_t x = -2; x < 0; ++x)  // todo make parametric?
            {
                if (implemented_circuit.circuit.super_circuit.gate_layout.get_node(tile<GateLyt>{t.x + x, t.y + 1}) !=
                    0)
                {
                    south_east = false;
                    break;
                }
            }

            const bool flip = !south_east && input_bdl_pairs.at(i).upper.y == input_bdl_pairs.at(i).lower.y;

            const bool bit_set = is_input_bit_set(south_east ? 0 : 1);  // todo make work with > 2 inputs

            if (bit_set ^ flip)
            {
                if (!(given_cds.get_charge_state(input_bdl_pairs.at(i).upper) == sidb_charge_state::NONE ||
                      given_cds.get_charge_state(input_bdl_pairs.at(i).upper) == sidb_charge_state::NEUTRAL) ||
                    given_cds.get_charge_state(input_bdl_pairs.at(i).lower) != sidb_charge_state::NEGATIVE)
                {
#ifdef PRINT_DEBUG
                    std::cout << "a" << std::endl;
                    if (given_cds.get_charge_state(input_bdl_pairs.at(i).lower) != sidb_charge_state::NEGATIVE)
                    {
                        std::cout << "NEG failure at index " << uint64_t{i} << ", lower @ "
                                  << input_bdl_pairs.at(i).lower.x << ", " << input_bdl_pairs.at(i).lower.y
                                  << std::endl;
                    }
                    else
                    {
                        std::cout << "NEUT failure at index " << uint64_t{i} << ", upper @ "
                                  << input_bdl_pairs.at(i).upper.x << ", " << input_bdl_pairs.at(i).upper.y
                                  << std::endl;
                    }
#endif
                    return operational_status::NON_OPERATIONAL;
                }
            }
            else
            {
                if (given_cds.get_charge_state(input_bdl_pairs.at(i).upper) != sidb_charge_state::NEGATIVE ||
                    !(given_cds.get_charge_state(input_bdl_pairs.at(i).lower) == sidb_charge_state::NEUTRAL ||
                      given_cds.get_charge_state(input_bdl_pairs.at(i).lower) == sidb_charge_state::NONE))
                {
#ifdef PRINT_DEBUG
                    std::cout << "b" << std::endl;
                    if (given_cds.get_charge_state(input_bdl_pairs.at(i).upper) != sidb_charge_state::NEGATIVE)
                    {
                        std::cout << "NEG failure at index " << uint64_t{i} << ", upper @ "
                                  << input_bdl_pairs.at(i).upper.x << ", " << input_bdl_pairs.at(i).upper.y
                                  << std::endl;
                    }
                    else
                    {
                        std::cout << "NEUT failure at index " << uint64_t{i} << ", lower @ "
                                  << input_bdl_pairs.at(i).lower.x << ", " << input_bdl_pairs.at(i).lower.y
                                  << std::endl;
                    }
#endif
                    return operational_status::NON_OPERATIONAL;
                }
            }
        }

        const std::vector<bdl_pair<cell<Lyt>>>& bdl_pairs =
            detect_bdl_pairs(implemented_circuit.cell_layout, sidb_technology::cell_type::OUTPUT,
                             detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance});

        const bool bit_set = kitty::get_bit(implemented_circuit.circuit.function, input_pattern);

#ifdef PRINT_DEBUG
        std::cout << "BDL pairs: " << bdl_pairs.size() << std::endl;
        std::cout << "function: ";
        kitty::print_binary(implemented_circuit.circuit.function);
        std::cout << "\tbit set: " << bit_set << std::endl;
#endif

        for (uint8_t i = 0; i < bdl_pairs.size(); ++i)
        {
            bool flip = false;

            if (bdl_pairs.at(i).upper.y == bdl_pairs.at(i).lower.y)
            {
                const tile<GateLyt>& t =
                    implemented_circuit.circuit.super_circuit.canvasses_lyt.template get_cell_tile<tile<GateLyt>>(
                        bdl_pairs.at(i).upper);

                for (int8_t x = -2; x < 0; ++x)  // todo make parametric?
                {
                    if (implemented_circuit.circuit.super_circuit.gate_layout.get_node(
                            tile<GateLyt>{t.x + x, t.y + 1}) != 0)
                    {
                        flip = true;
                        break;
                    }
                }
            }

            if (bit_set ^ flip)
            {
                if (!(given_cds.get_charge_state(bdl_pairs.at(i).upper) == sidb_charge_state::NONE ||
                      given_cds.get_charge_state(bdl_pairs.at(i).upper) == sidb_charge_state::NEUTRAL) ||
                    given_cds.get_charge_state(bdl_pairs.at(i).lower) != sidb_charge_state::NEGATIVE)
                {
#ifdef PRINT_DEBUG
                    std::cout << "c" << std::endl;
#endif
                    return operational_status::NON_OPERATIONAL;
                }
            }
            else
            {
                if (given_cds.get_charge_state(bdl_pairs.at(i).upper) != sidb_charge_state::NEGATIVE ||
                    !(given_cds.get_charge_state(bdl_pairs.at(i).lower) == sidb_charge_state::NEUTRAL ||
                      given_cds.get_charge_state(bdl_pairs.at(i).lower) == sidb_charge_state::NONE))
                {
#ifdef PRINT_DEBUG
                    std::cout << "d" << std::endl;
#endif
                    return operational_status::NON_OPERATIONAL;
                }
            }
        }

        if (const cell<Lyt> op{bounding_box_2d<Lyt>{implemented_circuit.cell_layout}.get_max()};
            implemented_circuit.cell_layout.get_cell_type(op) == sidb_technology::cell_type::OUTPUT_PERTURBER &&
            given_cds.get_charge_state(op) != sidb_charge_state::NEGATIVE)
        {
#ifdef PRINT_DEBUG
            std::cout << "e" << std::endl;
#endif
            return operational_status::NON_OPERATIONAL;
        }

        return operational_status::OPERATIONAL;

        /*
        const auto is_bit_set = [&](const uint64_t input_number)
        {
            return (input_pattern &
                    (uint64_t{1ull} << (implemented_circuit.circuit.function.num_vars() - 1 - input_number))) != 0ull;
        };

        const std::vector<bdl_pair<cell<Lyt>>>& input_bdl_pairs =
            detect_bdl_pairs(implemented_circuit.circuit.canvasses_lyt, sidb_technology::cell_type::INPUT);

        uint8_t bit_position = 0;

        for (uint8_t i = 0; i < input_bdl_pairs.size(); ++i)
        {
            if (is_bit_set(i))
            {
                if (given_cds.get_charge_state(input_bdl_pairs.at(i).lower) != sidb_charge_state::NEGATIVE)
                {
                    return operational_status::NON_OPERATIONAL;
                }
            }
            else
            {
                if (given_cds.get_charge_state(input_bdl_pairs.at(i).upper) != sidb_charge_state::NEGATIVE)
                {
                    return operational_status::NON_OPERATIONAL;
                }
            }

            if (i == 0)
            {
                bit_position += 2 * static_cast<uint8_t>(is_bit_set(i));
            }
            else if (i == 1)
            {
                bit_position += static_cast<uint8_t>(is_bit_set(i));
            }
            else
            {
                assert(false && "todo");
            }
        }

        const std::vector<bdl_pair<cell<Lyt>>>& output_bdl_pairs =
            detect_bdl_pairs(implemented_circuit.circuit.canvasses_lyt, sidb_technology::cell_type::OUTPUT);

        for (uint8_t i = 0; i < output_bdl_pairs.size(); ++i)
        {
            const mockturtle::node<GateLyt>& po_n = implemented_circuit.circuit.super_circuit.gate_layout.get_node(
                implemented_circuit.circuit.canvasses_lyt.template get_cell_tile<tile<GateLyt>>(
                    output_bdl_pairs.at(i).upper));

            std::vector<mockturtle::node<GateLyt>> logic_n{};
            implemented_circuit.circuit.super_circuit.gate_layout.foreach_fanin(po_n,
                                                                  [&logic_n](const auto& n) { logic_n.push_back(n); });
            assert(logic_n.size() == 1 && "precisely one node must connect to the PO");

            if (kitty::get_bit(implemented_circuit.circuit.super_circuit.gate_layout.node_function(logic_n.front()),
        bit_position))
            {
                if (given_cds.get_charge_state(output_bdl_pairs.at(i).upper) != sidb_charge_state::NEUTRAL ||
                    given_cds.get_charge_state(output_bdl_pairs.at(i).lower) != sidb_charge_state::NEGATIVE)
                {
                    return operational_status::NON_OPERATIONAL;
                }
            }
            else
            {
                if (given_cds.get_charge_state(output_bdl_pairs.at(i).upper) != sidb_charge_state::NEGATIVE ||
                    given_cds.get_charge_state(output_bdl_pairs.at(i).lower) != sidb_charge_state::NEUTRAL)
                {
                    return operational_status::NON_OPERATIONAL;
                }
            }
        }

        return operational_status::OPERATIONAL;*/
    }

    [[nodiscard]] operational_status
    determine_status_and_logic_match(const uint64_t i, const uint64_t s_i,
                                     const sidb_simulation_result<Lyt, ExtPotType>& simulation_results) noexcept
    {
        if constexpr (ExtPotType == local_external_potential_type::BOUNDED)
        {
            for (const uint64_t cds_ix : get_possible_ground_state_indices(simulation_results.charge_distributions))
            {
                if (assess_logic_match_of_charge_distribution(simulation_results.charge_distributions.at(cds_ix), i,
                                                              s_i) == operational_status::OPERATIONAL)
                {
#ifdef PRINT_DEBUG
                    std::cout << "verdict: operational\n" << std::endl;
#endif
                    return operational_status::OPERATIONAL;
                }
            }

#ifdef PRINT_DEBUG
            std::cout << "verdict: NON-operational\n" << std::endl;
#endif
            return operational_status::NON_OPERATIONAL;
        }
        else
        {
            const auto ground_states = simulation_results.groundstates();

            for (const auto& gs : ground_states)
            {
                if (assess_logic_match_of_charge_distribution(gs, i, s_i) == operational_status::NON_OPERATIONAL)
                {
                    return operational_status::NON_OPERATIONAL;

                    break;
                }
            }

            return operational_status::OPERATIONAL;
        }
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
          local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED>
[[nodiscard]] operational_status
is_circuit_operational(const sidb_cell_level_bdl_circuit<Lyt, GateLyt>&     implemented_circuit,
                       const is_circuit_operational_params&                 params = {},
                       const std::unique_ptr<detail::thread_count_manager>& tcm    = nullptr) noexcept
{
    static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");

    // assert(implemented_circuit.cell_layout.num_pis() > 0 && "lyt needs input cells");
    // assert(implemented_circuit.cell_layout.num_pos() > 0 && "lyt needs output cells");

    // assert(implemented_circuit.circuit.gate_layout.num_pis() * 2 == implemented_circuit.cell_layout.num_pis() &&
    //        "Each PI in the gate lyt needs to be implemented by a BDL pair");
    // assert(implemented_circuit.circuit.gate_layout.num_pos() * 2 == implemented_circuit.cell_layout.num_pos() &&
    //        "Each PO in the gate lyt needs to be implemented by a BDL pair");

    // todo
    if (implemented_circuit.cell_layout.num_cells() == 0 ||
        detect_bdl_pairs(implemented_circuit.cell_layout, std::nullopt,
                         detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance})
            .empty())
    {
        return operational_status::OPERATIONAL;
    }

#ifdef PRINT_DEBUG
    std::cout << std::endl;
    print_layout(implemented_circuit.circuit.canvasses_lyt);
    print_layout(implemented_circuit.cell_layout);
    std::cout << "tiles: ";
    for (const auto& n : implemented_circuit.circuit.nodes)
    {
        std::cout << implemented_circuit.circuit.super_circuit.gate_layout.get_tile(n) << '\t';
    }
    std::cout << std::endl;
    implemented_circuit.cell_layout.foreach_cell(
        [&](const auto& c)
        {
            std::string type{};
            switch (implemented_circuit.cell_layout.get_cell_type(c))
            {
                case sidb_technology::cell_type::INPUT: type = "INPUT"; break;
                case sidb_technology::cell_type::OUTPUT: type = "OUTPUT"; break;
                case sidb_technology::cell_type::OUTPUT_PERTURBER: type = "OUTPUT_PERTURBER"; break;
                case sidb_technology::cell_type::NORMAL: type = "NORMAL"; break;
                case sidb_technology::cell_type::LOGIC: type = "LOGIC"; break;
                default: type = "UNKNOWN";
            };
            std::cout << type << " @ " << c.x << ", " << c.y << std::endl;
        });
    std::cout << "num consistent super circuit inputs:" << std::endl;

    for (uint64_t i = 0; i < implemented_circuit.circuit.consistent_super_circuit_input_indices_per_input.size(); ++i)
    {
        const std::vector<uint64_t>& s_ix =
            implemented_circuit.circuit.consistent_super_circuit_input_indices_per_input.at(i);

        std::cout << "sub circuit input index: " << i << '\t';
        for (const uint64_t sixx : s_ix)
        {
            std::cout << sixx << ", ";
        }
        std::cout << std::endl;
    }
#endif

    detail::is_circuit_operational_impl<Lyt, GateLyt, ExtPotType> p{implemented_circuit, params, tcm};

    const auto& status = p.run();

    if (params.print)
    {
        std::cout << "\n\nOVERALL STATUS: " << (status == operational_status::NON_OPERATIONAL ? "NON-" : "")
                  << "OPERATIONAL" << std::endl;

        std::cout << std::endl;
    }

    return status;
}

}  // namespace fiction

#endif  // FICTION_IS_CIRCUIT_OPERATIONAL_HPP

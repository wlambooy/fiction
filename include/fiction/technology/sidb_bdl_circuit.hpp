//
// Created by Willem Lambooy on 06/06/2025.
//

#ifndef SIDB_BDL_CIRCUIT_HPP
#define SIDB_BDL_CIRCUIT_HPP

#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/technology/sidb_bdl_skeletons.hpp"
#include "fiction/technology/sidb_defect_surface.hpp"
#include "fiction/traits.hpp"

#include <array>
#include <set>
#include <utility>
#include <vector>

#include <signal.h>

namespace fiction
{

/**
 * This struct stores the parameters to design an SiDB circuit on a defective surface.
 *
 * @tparam CellLyt SiDB cell-level layout type.
 */
template <typename CellLyt>
struct sidb_bdl_circuit_params
{
    /**
     * This struct holds parameters to design SiDB gates.
     */
    design_sidb_gates_params<CellLyt> design_gate_params{};
    /** This variable specifies the radius in nanometers around the center of the hexagon where atomic defects are
     * incorporated into the gate design. (unit: nm)
     */

    std::optional<sidb_defect_surface<CellLyt>> defect_surface{};
    double                                      influence_radius_charged_defects = 15;

    uint64_t num_trials          = 500;
    double   quantization_factor = 0.075;
    double   selectivity         = 0.5;

    uint64_t num_trials_for_double_scope          = 100;
    double   quantization_factor_for_double_scope = 0.005;
    double   selectivity_for_double_scope         = 0.6;

    uint64_t num_trials_for_global_scope          = 20;
    double   quantization_factor_for_global_scope = 0.025;
    double   selectivity_for_global_scope         = 0.8;

    double excited_state_alpha = 1.0;

    uint64_t available_threads = std::thread::hardware_concurrency();
};

template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
class sidb_bdl_circuit
{
  public:
    explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const sidb_bdl_circuit_params<CellLyt>& bdl_circuit_parameters,
                              const std::optional<tile<GateLyt>>& this_tile = std::nullopt) noexcept :
            gate_layout{gate_lyt.clone()},
            bdl_circuit_params{bdl_circuit_parameters},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt)},
            bdl_wires{detect_bdl_wires(
                skeleton,
                bdl_circuit_params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params)},
            num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_circuit_params.design_gate_params.operational_params
                                                          .input_bdl_iterator_params.bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            gate_connections{get_gate_connections(bdl_wires, skeleton, gate_lyt)},
            gate_tile{this_tile},
            num_gates_to_design{get_number_of_gates_to_design(gate_lyt)}
    {
        operational_params.simulation_parameters =
            bdl_circuit_params.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params =
            bdl_circuit_params.design_gate_params.operational_params.input_bdl_iterator_params;
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;
        // is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;
        operational_params.excited_state_alpha = bdl_circuit_params.excited_state_alpha;
    }

    // explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const tile<GateLyt>& t, const tile<GateLyt>& connecting_t,
    //                           const sidb_bdl_circuit_params<CellLyt>& bdl_wire_params) noexcept :
    //         gate_layout{create_gate_lyt_window_for_gate_connection(gate_lyt, t, connecting_t)},
    //         skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout)},
    //         bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
    //         num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
    //         input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
    //                                                   bdl_wire_params.bdl_pairs_params)},
    //         num_inputs{input_bdl_pairs.size()},
    //         gate_connections{get_gate_connections(bdl_wires, skeleton, gate_layout)}
    // {}
    //
    // explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const tile<GateLyt>& t, const tile<GateLyt>& connecting_t,
    //                           const tile<GateLyt>&           connecting_to_connecting_t,
    //                           const detect_bdl_wires_params& bdl_wire_params) noexcept :
    //         gate_layout{
    //             create_gate_lyt_window_for_two_gate_connections(gate_lyt, t, connecting_t,
    //             connecting_to_connecting_t)},
    //         skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout)},
    //         bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
    //         num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
    //         input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
    //                                                   bdl_wire_params.bdl_pairs_params)},
    //         num_inputs{input_bdl_pairs.size()},
    //         gate_connections{get_gate_connections(bdl_wires, skeleton, gate_layout)}
    // {}

    [[nodiscard]] std::optional<CellLyt> design_circuit()
    {
        // initialize
        collect_initial_gate_designs();

        while (++circuit_design_level < num_gates_to_design)
        {
            // prune by assessing gate design combinations for increasingly large sets of connected gates
            prune_gate_designs();
        }

        // prune at the global level (all gates are considered together)
        if (const std::optional<CellLyt>& maybe_lyt = prune_gate_designs(); maybe_lyt.has_value())
        {
            return maybe_lyt.value();
        }

        return exhaustively_enumerate_gate_design_combinations();
    }

    /**
     * SiDB gate-level layout.
     */
    const GateLyt gate_layout;

    const sidb_bdl_circuit_params<CellLyt> bdl_circuit_params;

    is_circuit_operational_params<CellLyt> operational_params{};

    const CellLyt skeleton;

    const std::vector<bdl_wire<CellLyt>>                       bdl_wires{};
    const uint64_t                                             num_bdl_pairs{};
    const std::vector<bdl_pair<cell<CellLyt>>>                 input_bdl_pairs{};
    const uint64_t                                             num_inputs{};
    const std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};

    const std::optional<tile<GateLyt>> gate_tile{};

    const uint64_t num_gates_to_design{};

    using gate_designs_per_node =
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename SkeletonGateLibrary::fcn_gate>>;

    gate_designs_per_node gate_designs{};

    uint64_t circuit_design_level = 0;

  private:
    void collect_initial_gate_designs()
    {
        gate_layout->foreach_node(
            [&, this](const auto& n, [[maybe_unused]] auto i)
            {
                if (!skip_physical_design_for_node(*gate_layout, n))
                {
                    gate_designs[n] =  // design_gates_for_node(n);
                        SkeletonGateLibrary::template set_up_gates<GateLyt, CellLyt,
                                                                   sidb_on_the_fly_gate_library_params<CellLyt>,

                                                                   local_external_potential_type::BOUNDED>(
                            gate_layout, gate_layout.get_tile(n),
                            {bdl_circuit_params.design_gate_params,
                             bdl_circuit_params.influence_radius_charged_defects},
                            bdl_circuit_params.defect_surface, std::make_optional(this), operational_params);
                }
            });
    }

    // std::vector<typename SkeletonGateLibrary::fcn_gate> design_gates_for_node(const mockturtle::node<GateLyt>& n)
    // {
    //     const auto t = gate_layout.get_tile(n);
    //     const auto f = gate_layout.node_function(n);
    //     const auto p = SkeletonGateLibrary::determine_port_routing(gate_layout, t);
    //
    //     auto center_cell = relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
    //                                                           SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
    //         gate_layout, t,
    //         cell<CellLyt>{SkeletonGateLibrary::gate_x_size() / 2, SkeletonGateLibrary::gate_y_size() / 2});
    //     auto absolute_cell =
    //         relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
    //         SkeletonGateLibrary::gate_y_size(),
    //                                            GateLyt, CellLyt>(gate_layout, t, cell<CellLyt>{0, 0});
    //
    //     const auto cell_list = sidb_bdl_skeleton_1{}.set_up_gate(gate_layout, t);
    //     if (cell_list == SkeletonGateLibrary::EMPTY_GATE)
    //     {
    //         return {SkeletonGateLibrary::EMPTY_GATE};
    //     }
    //
    //     const auto skeleton = cell_list_to_cell_level_layout<CellLyt>(cell_list);
    //
    //     const auto design_gates = [&]()
    //     {
    //         std::cout << "starting gate design for tile " << tile << "\t|\tnode function:";
    //         for (const tt& tt : spec)
    //         {
    //             std::cout << '\t';
    //             kitty::print_binary(tt);
    //         }
    //         std::cout << std::endl;
    //
    //         const auto found_gate_layouts =
    //             design_sidb_gates<LytSkeleton, TT, ExtPotType, GateLyt, SkeletonGateLibrary>(
    //                 skeleton, spec, parameters.design_gate_params, nullptr, std::make_optional(std::move(circuit)),
    //                 super_circuit, op_params);
    //     }
    //
    //     try
    //     {
    //         if constexpr (fiction::has_is_fanout_v<GateLyt>)
    //         {
    //             if (gate_layout.is_fanout(n))
    //             {
    //                 if (gate_layout.fanout_size(n) == 2)
    //                 {
    //                     if constexpr (is_sidb_defect_surface_v<CellLyt>)
    //                     {
    //                         if (bdl_circuit_params.defect_surface.has_value())
    //                         {
    //                             const auto skeleton_with_defects = add_defect_to_skeleton(
    //                                 bdl_circuit_params.defect_surface.value(), skeleton,
    //                                 bdl_circuit_params.influence_radius_charged_defects, center_cell, absolute_cell);
    //
    //                             return design_gates<CellLyt, tt, CellLyt, GateLyt,
    //                                                 local_external_potential_type::BOUNDED, SkeletonGateLibrary>(
    //                                 skeleton_with_defects, create_fan_out_tt(), params, p, t,
    //                                 make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                                     lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params));
    //                         }
    //                     }
    //                     return design_gates<CellLyt, tt, CellLyt, GateLyt, local_external_potential_type::BOUNDED,
    //                                         SkeletonGateLibrary>(
    //                         skeleton, create_fan_out_tt(), params, p, t,
    //                         make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                             lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params));
    //                 }
    //             }
    //         }
    //         if constexpr (fiction::has_is_buf_v<GateLyt>)
    //         {
    //             if (gate_layout.is_buf(n))
    //             {
    //                 if (gate_layout.is_ground_layer(t))
    //                 {
    //                     // crossing case
    //                     if (const auto at = gate_layout.above(t); (t != at) && gate_layout.is_wire_tile(at))
    //                     {
    //                         // two possible options: actual crossover and (parallel) hourglass wire
    //                         const auto pa = SkeletonGateLibrary::determine_port_routing(gate_layout, at);
    //
    //                         const auto spec = TWO_IN_TWO_OUT_MAP.at({p, pa});
    //
    //                         auto complex_gate_param               = params;
    //                         complex_gate_param.design_gate_params = params.design_gate_params_complex_gates;
    //
    //                         complex_gate_param.design_gate_params.operational_params.cc_map =
    //                             params.design_gate_params.operational_params.cc_map;
    //
    //                         if constexpr (is_sidb_defect_surface_v<CellLyt>)
    //                         {
    //                             if (defect_surface.has_value())
    //                             {
    //                                 const auto skeleton_with_defects = add_defect_to_skeleton(
    //                                     defect_surface.value(), skeleton, params.influence_radius_charged_defects,
    //                                     center_cell, absolute_cell);
    //
    //                                 return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType,
    //                                 SkeletonGateLibrary>(
    //                                     skeleton_with_defects, spec, complex_gate_param, p, t,
    //                                     make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                                         lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
    //                                     super_circuit, op_params);
    //                             }
    //                         }
    //
    //                         return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
    //                             skeleton, spec, complex_gate_param, p, t,
    //                             make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                                 lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
    //                             super_circuit, op_params);
    //                     }
    //
    //                     if constexpr (is_sidb_defect_surface_v<CellLyt>)
    //                     {
    //                         if (defect_surface.has_value())
    //                         {
    //                             const auto skeleton_with_defects = add_defect_to_skeleton(
    //                                 defect_surface.value(), skeleton, params.influence_radius_charged_defects,
    //                                 center_cell, absolute_cell);
    //                             return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
    //                                 skeleton_with_defects, std::vector<tt>{f}, params, p, t,
    //                                 make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                                     lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
    //                                 super_circuit, op_params);
    //                         }
    //                     }
    //
    //                     return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
    //                         skeleton, std::vector<tt>{f}, params, p, t,
    //                         make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                             lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
    //                         super_circuit, op_params);
    //                 }
    //                 return {SkeletonGateLibrary::EMPTY_GATE};
    //             }
    //         }
    //
    //         if constexpr (is_sidb_defect_surface_v<CellLyt>)
    //         {
    //             if (bdl_circuit_params.defect_surface.has_value())
    //             {
    //                 const auto skeleton_with_defects = add_defect_to_skeleton(
    //                     bdl_circuit_params.defect_surface.value(), skeleton,
    //                     bdl_circuit_params.influence_radius_charged_defects, center_cell, absolute_cell);
    //
    //                 return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
    //                     skeleton_with_defects, std::vector<tt>{f}, params, p, t,
    //                     make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                         lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
    //                     super_circuit, op_params);
    //             }
    //         }
    //
    //         return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
    //             skeleton, std::vector<tt>{f}, params, p, t,
    //             make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
    //                 lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
    //             super_circuit, op_params);
    //     }
    //
    //     catch (const std::out_of_range&)
    //     {
    //         throw unsupported_gate_orientation_exception(t, p);
    //     }
    //
    //     throw unsupported_gate_type_exception(t);
    // }

    void prune_gate_designs() noexcept {}

    std::optional<CellLyt> exhaustively_enumerate_gate_design_combinations() const noexcept
    {
        std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;

        is_circuit_operational_params operational_params{};
        operational_params.simulation_parameters =
            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params = params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                                           .operational_params.input_bdl_iterator_params;
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        // operational_params.print = true;

        std::vector<uint64_t> indices(operational_gate_designs.size(), 0);

        while (true)
        {
            CellLyt operational_circuit_candidate{};
            for (uint64_t i = 0; i < operational_gate_designs.size(); i++)
            {
                const auto& [n, op_gate_designs_for_gate] =
                    *std::next(operational_gate_designs.cbegin(), static_cast<int64_t>(i));
                // select a random gate implementation for the tile that connects as input to n
                assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(operational_circuit_candidate,
                                                                   op_gate_designs_for_gate.at(indices.at(i)), gate_lyt,
                                                                   gate_lyt.get_tile(n));
            }

            std::cout << "trying combination: ";
            for (uint64_t i = 0; i < operational_gate_designs.size(); i++)
            {
                std::cout << indices.at(i) << " ";
            }
            std::cout << std::endl;

            if (is_circuit_operational(
                    sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{operational_circuit_candidate,
                                                                                       *circuit},
                    operational_params)
                    .status == operational_status::OPERATIONAL)
            {
                lyt.emplace(std::move(operational_circuit_candidate));

                std::cout << "\n\nFINAL GENERATED CIRCUIT:" << std::endl;
                print_layout(lyt.value());

                return true;
            }

            // Increment indices like an odometer
            for (uint64_t i = 0; i < indices.size(); ++i)
            {
                if (++indices[i] < std::next(operational_gate_designs.cbegin(), static_cast<int64_t>(i))->second.size())
                {
                    break;  // No carry needed
                }

                indices[i] = 0;  // Reset this index and carry over to the next

                if (i == indices.size() - 1)
                {
                    return false;  // Stop when the last index overflows
                }
            }
        }

        return false;
    }

    static void augment_gate_lyt_window(const GateLyt& gate_lyt, GateLyt& gate_lyt_window,
                                        const tile<GateLyt>&           current_t,
                                        const std::set<tile<GateLyt>>& connecting_inputs  = {},
                                        const std::set<tile<GateLyt>>& connecting_outputs = {}) noexcept
    {
        std::vector<mockturtle::signal<GateLyt>> inputs_to_current_t{};

        for (auto in_t : gate_lyt.incoming_data_flow(current_t))
        {
            if (connecting_inputs.count(gate_lyt.below(in_t)) == 0)
            {
                if (!gate_lyt_window.is_empty_tile(in_t))
                {
                    in_t.z = 1 - in_t.z;
                }

                gate_lyt_window.create_pi("", in_t);
            }

            inputs_to_current_t.push_back(static_cast<mockturtle::signal<GateLyt>>(in_t));
        }

        assert(gate_lyt_window.is_empty_tile(current_t) && "tile on which node is to be created is already populated");
        gate_lyt_window.create_node(inputs_to_current_t, gate_lyt.node_function(gate_lyt.get_node(current_t)),
                                    current_t);

        for (auto out_t : gate_lyt.outgoing_data_flow(current_t))
        {
            if (connecting_outputs.count(gate_lyt.below(out_t)) == 0)
            {
                if (!gate_lyt_window.is_empty_tile(out_t))
                {
                    out_t.z = 1 - out_t.z;
                }

                gate_lyt_window.create_po(static_cast<mockturtle::signal<GateLyt>>(current_t), "", out_t);
            }
        }

        if (gate_lyt.is_wire_tile(current_t))
        {
            if (const auto above_current_t = gate_lyt.above(current_t);
                current_t != above_current_t && gate_lyt.is_wire_tile(above_current_t))
            {
                // upper_t hosts a crossing or double wire

                augment_gate_lyt_window(gate_lyt, gate_lyt_window, above_current_t, connecting_inputs,
                                        connecting_outputs);
            }
        }
    }

    [[nodiscard]] static GateLyt
    create_gate_lyt_window_for_connection_sequence(const GateLyt& gate_lyt, std::vector<tile<GateLyt>> tiles) noexcept
    {
        GateLyt gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};

        // Sort tiles in ascending y, then x, then z for deterministic ordering
        std::sort(tiles.begin(), tiles.end());

        const auto same_clock_zone = [](const tile<GateLyt>& a, const tile<GateLyt>& b)
        {
            return a.y == b.y;  // Adjust if clock zone definition differs from y==y
        };

        // Process each tile in sorted order
        for (size_t i = 0; i < tiles.size(); ++i)
        {
            std::set<tile<GateLyt>> inputs, outputs;

            // --- Determine inputs ---
            for (size_t j = 0; j < i; ++j)
            {
                // Rule: in "shared clock zone" cases, skip direct non-adjacent wiring
                if (same_clock_zone(tiles[j], tiles[i]))
                    continue;

                // If immediate predecessor is in same zone, also take inputs from its inputs
                // (this is the "both feed lower" branch in your original)
                if (i > 0 && same_clock_zone(tiles[i - 1], tiles[i]))
                {
                    if (j == i - 1)  // direct predecessor in same zone → don't connect
                        continue;
                }

                inputs.insert(tiles[j]);
            }

            // --- Determine outputs ---
            for (size_t j = i + 1; j < tiles.size(); ++j)
            {
                if (same_clock_zone(tiles[i], tiles[j]))
                    continue;

                // If immediate successor is in same zone, also give outputs to both
                // (this is the "upper gives output to both" branch)
                if (j == i + 1 && j + 1 < tiles.size() && same_clock_zone(tiles[j], tiles[j + 1]))
                {
                    outputs.insert(tiles[j]);
                    outputs.insert(tiles[j + 1]);
                    break;  // handled both at once
                }

                outputs.insert(tiles[j]);
                break;  // only connect to first in next zone
            }

            augment_gate_lyt_window(gate_lyt, gate_lyt_window, tiles[i], inputs, outputs);
        }

        return gate_lyt_window;
    }

    [[nodiscard]] static GateLyt create_gate_lyt_window_for_gate_connection(const GateLyt&       gate_lyt,
                                                                            const tile<GateLyt>& t,
                                                                            const tile<GateLyt>& connecting_t) noexcept
    {
        GateLyt              gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};
        const tile<GateLyt>& upper_t = t.y < connecting_t.y ? t : connecting_t;
        const tile<GateLyt>& lower_t = t.y < connecting_t.y ? connecting_t : t;

        augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {lower_t});
        augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {upper_t}, {});

        return gate_lyt_window;
    }

    [[nodiscard]] static GateLyt
    create_gate_lyt_window_for_two_gate_connections(const GateLyt& gate_lyt, const tile<GateLyt>& t,
                                                    const tile<GateLyt>& connecting_t,
                                                    const tile<GateLyt>& connecting_to_connecting_t) noexcept
    {
        GateLyt gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};

        std::array<tile<GateLyt>, 3> tiles_sorted{{t, connecting_t, connecting_to_connecting_t}};
        std::sort(tiles_sorted.begin(), tiles_sorted.end());

        const tile<GateLyt>& upper_t  = tiles_sorted[0];
        const tile<GateLyt>& middle_t = tiles_sorted[1];
        const tile<GateLyt>& lower_t  = tiles_sorted[2];

        assert(upper_t.y <= middle_t.y && middle_t.y <= lower_t.y && "tiles are not sorted by row number");

        if (middle_t.y == lower_t.y || middle_t.y == upper_t.y)
        {
            // a clock zone is shared

            assert(upper_t.y != lower_t.y && "all tiles are in the same clockzone");

            if (middle_t.y < lower_t.y)
            {
                // case: lower_t follows the others and gets input from both
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {lower_t});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, middle_t, {}, {lower_t});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {upper_t, middle_t}, {});
            }
            else
            {
                // case: upper_t precedes the others and gives output to both
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {middle_t, lower_t});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, middle_t, {upper_t}, {});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {upper_t}, {});
            }

            return gate_lyt_window;
        }

        augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {middle_t});
        augment_gate_lyt_window(gate_lyt, gate_lyt_window, middle_t, {upper_t}, {lower_t});
        augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {middle_t}, {});

        return gate_lyt_window;
    }
    /**
     *
     */
    [[nodiscard]] static std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>>
    get_gate_connections(const std::vector<bdl_wire<CellLyt>>& bdl_wires, const CellLyt& lyt,
                         const GateLyt& gate_lyt) noexcept
    {
        std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};
        gate_connections.reserve(bdl_wires.size());

        for (const bdl_wire<CellLyt>& wire : bdl_wires)
        {
            assert(!wire.pairs.empty() && "BDL wire is empty");

            const bool right_to_left = wire.pairs.front().upper.x > wire.pairs.front().lower.x;
            const bool is_input_wire = wire.pairs.front().type == sidb_technology::cell_type::INPUT;

            std::set<tile<GateLyt>> tiles_in_wire = {};

            for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
            {
                // std::cout << "\ncell: " << pair.upper.x << " " << pair.upper.y << " at tile "
                //           << lyt.get_cell_tile(pair.upper).x << ',' << lyt.get_cell_tile(pair.upper).y << ','
                //           << lyt.get_cell_tile(pair.upper).z << std::endl;
                // std::cout << "cell: " << pair.lower.x << " " << pair.lower.y << " at tile "
                //           << lyt.get_cell_tile(pair.lower).x << ',' << lyt.get_cell_tile(pair.lower).y << ','
                //           << lyt.get_cell_tile(pair.lower).z << std::endl;
                assert(lyt.get_cell_tile(pair.upper) == lyt.get_cell_tile(pair.lower));

                tiles_in_wire.insert(
                    [&](const auto& this_t)
                    {
                        const tile<GateLyt> t = {this_t.x, this_t.y};
                        auto                z = 0;

                        // std::cout << "current tile: " << t.x << ", " << t.y << std::endl;

                        if (const auto at = gate_lyt.above(t); at != t && gate_lyt.is_wire_tile(at))
                        {
                            assert(gate_lyt.is_wire_tile(t) && "stacked tiles are not both wires");

                            const bool is_input_wire_to_this_tile = is_input_wire || !tiles_in_wire.empty();

                            // stacked tiles are both populated with wire tiles  ---  CROSSING or DOUBLE WIRE

                            // todo: does NOT work for crossing / double wire in connection .. ?

                            if (gate_lyt.below(gate_lyt.north_west(at)) ==
                                gate_lyt.below(gate_lyt.incoming_data_flow(at).front()))
                            {
                                // upper tile connects to north west

                                if (gate_lyt.below(gate_lyt.south_east(at)) ==
                                    gate_lyt.below(gate_lyt.outgoing_data_flow(at).front()))
                                {
                                    // upper tile connects to south east  ---  CROSSING CASE (upper tile goes L->R)
                                    if (!right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                                else
                                {
                                    // lower tile connects to south east  ---  DOUBLE WIRE CASE
                                    if (is_input_wire_to_this_tile != right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                            }
                            else
                            {
                                // lower tile connects to north west

                                if (gate_lyt.below(gate_lyt.south_east(at)) ==
                                    gate_lyt.below(gate_lyt.outgoing_data_flow(at).front()))
                                {
                                    // upper tile connects to south east  ---  DOUBLE WIRE CASE
                                    if (is_input_wire_to_this_tile == right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                                else
                                {
                                    // lower tile connects to south east  ---  CROSSING CASE (upper tile goes R->L)
                                    if (right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                            }
                        }
                        else if (!is_input_wire && !tiles_in_wire.empty())
                        {
                            const std::vector<tile<GateLyt>>& outgoing_tiles =
                                gate_lyt.outgoing_data_flow(*tiles_in_wire.cbegin());

                            if (outgoing_tiles.size() == 1)
                            {
                                z = outgoing_tiles.front().z;
                            }
                            else
                            {
                                assert(outgoing_tiles.size() == 2 && "fan out must be either 1 or 2");

                                if (right_to_left)
                                {
                                    z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.front().z :
                                                                                             outgoing_tiles.back().z;
                                }
                                else
                                {
                                    z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.back().z :
                                                                                             outgoing_tiles.front().z;
                                }
                            }
                        }

                        // std::cout << "z: " << z << std::endl;

                        return tile<GateLyt>{t.x, t.y, z};
                    }(lyt.get_cell_tile(pair.upper)));
            }

            assert((tiles_in_wire.size() == 2 ||
                    (tiles_in_wire.size() == 1 &&
                     (is_input_wire || wire.pairs.back().type == sidb_technology::cell_type::OUTPUT))) &&
                   "Gate connection is malformed");

            const auto get_tile_pair = [&]
            {
                // std::cout << "\ntile in wire: " << tiles_in_wire.cbegin()->x << ',' << tiles_in_wire.cbegin()->y <<
                // ','
                // << tiles_in_wire.cbegin()->z << std::endl;
                if (tiles_in_wire.size() == 2)
                {
                    // std::cout << "tile TO wire: " << tiles_in_wire.crbegin()->x << ',' << tiles_in_wire.crbegin()->y
                    // << ',' << tiles_in_wire.crbegin()->z << std::endl;
                    const tile<GateLyt>& first_tile  = *tiles_in_wire.cbegin();
                    const tile<GateLyt>& second_tile = *tiles_in_wire.crbegin();

                    if (first_tile.y < second_tile.y || (first_tile.y == second_tile.y && first_tile.x < second_tile.x))
                    {
                        return std::make_pair(first_tile, second_tile);
                    }

                    return std::make_pair(second_tile, first_tile);
                }

                assert(tiles_in_wire.size() == 1 &&
                       (is_input_wire || wire.pairs.back().type == sidb_technology::cell_type::OUTPUT) &&
                       "Gate connection is malformed");

                const tile<GateLyt>& tile_in_wire = *tiles_in_wire.cbegin();

                // todo: does NOT work for crossing / double wire in connection

                if (wire.pairs.front().type == sidb_technology::cell_type::INPUT)
                {

                    const tile<GateLyt>& upper_tile = right_to_left ?
                                                          gate_lyt.below(gate_lyt.north_east(tile_in_wire)) :
                                                          gate_lyt.below(gate_lyt.north_west(tile_in_wire));

                    auto z = 0;

                    const std::vector<tile<GateLyt>>& incoming_tiles = gate_lyt.incoming_data_flow(tile_in_wire);

                    if (incoming_tiles.size() == 1)
                    {
                        z = incoming_tiles.front().z;
                    }
                    else if (incoming_tiles.size() == 2)
                    {
                        // assert(incoming_tiles.size() == 2 && "fan in must be either 1 or 2");

                        if (right_to_left)
                        {
                            z = incoming_tiles.front().x < incoming_tiles.back().x ? incoming_tiles.back().z :
                                                                                     incoming_tiles.front().z;
                        }
                        else
                        {
                            z = incoming_tiles.front().x < incoming_tiles.back().x ? incoming_tiles.front().z :
                                                                                     incoming_tiles.back().z;
                        }
                    }

                    // std::cout << "tile INP wire: " << upper_tile.x << ',' << upper_tile.y << ',' << z << std::endl;
                    return std::make_pair(tile<GateLyt>{upper_tile.x, upper_tile.y, z}, tile_in_wire);
                }

                const tile<GateLyt>& lower_tile = right_to_left ? gate_lyt.below(gate_lyt.south_west(tile_in_wire)) :
                                                                  gate_lyt.below(gate_lyt.south_east(tile_in_wire));

                auto z = 0;

                const std::vector<tile<GateLyt>>& outgoing_tiles = gate_lyt.outgoing_data_flow(tile_in_wire);

                if (outgoing_tiles.size() == 1)
                {
                    z = outgoing_tiles.front().z;
                }
                else if (outgoing_tiles.size() == 2)
                {
                    // assert(outgoing_tiles.size() == 2 && "fan out must be either 1 or 2");

                    if (right_to_left)
                    {
                        z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.front().z :
                                                                                 outgoing_tiles.back().z;
                    }
                    else
                    {
                        z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.back().z :
                                                                                 outgoing_tiles.front().z;
                    }
                }
                // std::cout << "tile OUT wire: " << lower_tile.x << ',' << lower_tile.y << ',' << z << std::endl;
                return std::make_pair(tile_in_wire, tile<GateLyt>{lower_tile.x, lower_tile.y, z});
            };

            std::pair<tile<GateLyt>, tile<GateLyt>> tile_pair_at_gate_connection = get_tile_pair();
            // std::cout << fmt::format("tile pair at gate connection: {},{},{}   {},{},{}",
            //                          tile_pair_at_gate_connection.first.x, tile_pair_at_gate_connection.first.y,
            //                          tile_pair_at_gate_connection.first.z, tile_pair_at_gate_connection.second.x,
            //                          tile_pair_at_gate_connection.second.y, tile_pair_at_gate_connection.second.z)
            //           << std::endl;
            // assert(tile_pair_at_gate_connection.first.y == tile_pair_at_gate_connection.second.y - 1 &&
            //        "tiles are not represent a row clocked gate connection"); todo

            gate_connections.push_back(std::move(tile_pair_at_gate_connection));
        }

        return gate_connections;
    }

    [[nodiscard]] static uint64_t get_number_of_bdl_pairs(const std::vector<bdl_wire<CellLyt>>& bdl_wires) noexcept
    {
        uint64_t total = 0;

        for (const bdl_wire<CellLyt>& wire : bdl_wires)
        {
            total += wire.pairs.size();
        }

        return total;
    }

    [[nodiscard]] static uint64_t get_number_of_gates_to_design(const GateLyt& gate_lyt) noexcept
    {
        uint64_t total = 0;

        gate_lyt.foreach_node(
            [&](const mockturtle::node<GateLyt>& n)
            {
                if (!skip_physical_design_for_node(gate_lyt, n))
                {
                    total++;
                }
            });

        return total;
    }
};
template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
struct sidb_cell_level_bdl_circuit
{
    /**
     * SiDB cell-level layout.
     */
    const CellLyt&                                          cell_layout;
    sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary> circuit{};

    // todo construct using std::unordered_map<tile<GateLyt>, canvas_positions>
    explicit sidb_cell_level_bdl_circuit(
        const CellLyt& lyt, const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& bdl_circuit) noexcept :
            cell_layout{lyt},
            circuit{bdl_circuit}
    {}
};

}  // namespace fiction

#endif  // SIDB_BDL_CIRCUIT_HPP

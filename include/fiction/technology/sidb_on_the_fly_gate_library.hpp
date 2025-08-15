//
// Created by Jan Drewniok 20.09.23.
//

#ifndef FICTION_SIDB_GATE_LIBRARY_HPP
#define FICTION_SIDB_GATE_LIBRARY_HPP

#include "fiction/algorithms/physical_design/compare_designed_sidb_gates.hpp"
#include "fiction/algorithms/physical_design/design_sidb_gates.hpp"
#include "fiction/algorithms/simulation/sidb/compare_by_ground_state_isolation.hpp"
#include "fiction/algorithms/simulation/sidb/is_circuit_operational.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/layouts/bounding_box.hpp"
#include "fiction/technology/cell_ports.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/fcn_gate_library.hpp"
#include "fiction/technology/is_sidb_gate_design_impossible.hpp"
#include "fiction/technology/sidb_bdl_circuit.hpp"
#include "fiction/technology/sidb_bdl_skeletons.hpp"
#include "fiction/technology/sidb_nm_distance.hpp"
#include "fiction/technology/sidb_on_the_fly_gate_library.hpp"
#include "fiction/technology/sidb_skeleton_gate_library.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/layout_utils.hpp"
#include "fiction/utils/truth_table_utils.hpp"

#include <phmap.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fiction
{
/**
 * This struct encapsulates parameters for the parameterized SiDB gate library.
 *
 * @tparam Lyt Cell-level layout type.
 */
template <typename Lyt>
struct sidb_on_the_fly_gate_library_params
{
    /**
     * This struct holds parameters to design SiDB gates.
     */
    design_sidb_gates_params<Lyt> design_gate_params{};
    /** This variable specifies the radius in nanometers around the center of the hexagon where atomic defects are
     * incorporated into the gate design. (unit: nm)
     */
    double influence_radius_charged_defects = 15;
};

/**
 * A parameterized gate library for SiDB technology. It allows the design of SiDB gates tailored to given atomic
 * defects, thus enabling the design of SiDB circuits in the presence of atomic defects. The skeleton (i.e., the
 * pre-defined input and output wires) is hexagonal.
 *
 * @tparam GateSizeX Width of a hexagon.
 * @tparam GateSizeY Height of a hexagon.
 */
template <uint16_t GateSizeX, uint16_t GateSizeY>
class sidb_on_the_fly_gate_library : public fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>
{
  public:
    explicit sidb_on_the_fly_gate_library() = delete;
    /**
     * Overrides the corresponding function in fcn_gate_library. Given a tile `t`, this function takes all necessary
     * information from the stored grid into account to design the correct fcn_gate representation for that tile. In
     * case there is no possible SiDB design, the blacklist is updated and an error fcn gate is returned.
     *
     * @tparam GateLyt Pointy-top hexagonal gate-level layout type.
     * @tparam CellLyt SiDB cell-level layout type.
     * @tparam Params Type of the parameter used for the gate library.
     * todo
     * @param lyt Layout that hosts tile `t`.
     * @param t Tile to be realized as a Bestagon gate.
     * @param params Parameter to design SiDB gates.
     * @param defect_surface Optional atomic defect surface in case atomic defects are present.
     * @return Bestagon gate representation of `t` including mirroring.
     */
    template <typename GateLyt, typename CellLyt, typename Params,
              local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED,
              typename SkeletonGateLibrary             = sidb_skeleton_bestagon_mini_library>
    static std::vector<typename sidb_on_the_fly_gate_library::fcn_gate> set_up_gates(
        const GateLyt& lyt, const tile<GateLyt>& t, Params& params,
        const std::optional<CellLyt>&                                                 defect_surface = std::nullopt,
        const std::optional<sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>>& super_circuit  = std::nullopt,
        const std::optional<is_circuit_operational_params>&                           op_params      = std::nullopt)
    {
        static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt must be a gate-level layout");
        static_assert(has_cube_coord_v<CellLyt>, "CellLyt must be based on cube coordinates");
        static_assert(is_cell_level_layout_v<CellLyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<CellLyt>, "Lyt is not an SiDB layout");

        const auto n = lyt.get_node(t);
        const auto f = lyt.node_function(n);
        const auto p = determine_port_routing(lyt, t);

        // center cell of the Bestagon tile. IMPORTANT: There is no center for the specified Bestagon library. The
        // middle is at 22.66666 (34*2/3). However, this is not an integer and does not specify a cell. Cell close to it
        // is chosen.
        auto center_cell =
            relative_to_absolute_cell_position<typename sidb_on_the_fly_gate_library::gate_x_size(),
                                               typename sidb_on_the_fly_gate_library::gate_y_size(), GateLyt, CellLyt>(
                lyt, t,
                cell<CellLyt>{typename sidb_on_the_fly_gate_library::gate_x_size() / 2,
                              typename sidb_on_the_fly_gate_library::gate_y_size() / 2});
        // center cell of the current tile
        auto absolute_cell =
            relative_to_absolute_cell_position<typename sidb_on_the_fly_gate_library::gate_x_size(),
                                               typename sidb_on_the_fly_gate_library::gate_y_size(), GateLyt, CellLyt>(
                lyt, t, cell<CellLyt>{0, 0});

        const auto cell_list = sidb_bdl_skeleton_1{}.set_up_gate(lyt, t);
        if (cell_list == typename sidb_on_the_fly_gate_library::EMPTY_GATE)
        {
            return {typename sidb_on_the_fly_gate_library::EMPTY_GATE};
        }

        const auto skeleton = cell_list_to_cell_level_layout<CellLyt>(cell_list);

        try
        {
            if constexpr (fiction::has_is_fanout_v<GateLyt>)
            {
                if (lyt.is_fanout(n))
                {
                    if (lyt.fanout_size(n) == 2)
                    {
                        if constexpr (is_sidb_defect_surface_v<CellLyt>)
                        {
                            if (defect_surface.has_value())
                            {
                                const auto skeleton_with_defects = add_defect_to_skeleton(
                                    defect_surface.value(), skeleton, params.influence_radius_charged_defects,
                                    center_cell, absolute_cell);

                                return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                                    skeleton_with_defects, create_fan_out_tt(), params, p, t,
                                    make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                                        lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                                    super_circuit, op_params);
                            }
                        }
                        return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                            skeleton, create_fan_out_tt(), params, p, t,
                            make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                                lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                            super_circuit, op_params);
                    }
                }
            }
            if constexpr (fiction::has_is_buf_v<GateLyt>)
            {
                if (lyt.is_buf(n))
                {
                    if (lyt.is_ground_layer(t))
                    {
                        // crossing case
                        if (const auto at = lyt.above(t); (t != at) && lyt.is_wire_tile(at))
                        {
                            // two possible options: actual crossover and (parallel) hourglass wire
                            const auto pa = determine_port_routing(lyt, at);

                            const auto spec = TWO_IN_TWO_OUT_MAP.at({p, pa});

                            auto complex_gate_param               = params;
                            complex_gate_param.design_gate_params = params.design_gate_params_complex_gates;

                            complex_gate_param.design_gate_params.operational_params.cc_map =
                                params.design_gate_params.operational_params.cc_map;

                            if constexpr (is_sidb_defect_surface_v<CellLyt>)
                            {
                                if (defect_surface.has_value())
                                {
                                    const auto skeleton_with_defects = add_defect_to_skeleton(
                                        defect_surface.value(), skeleton, params.influence_radius_charged_defects,
                                        center_cell, absolute_cell);

                                    return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                                        skeleton_with_defects, spec, complex_gate_param, p, t,
                                        make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                                            lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                                        super_circuit, op_params);
                                }
                            }

                            return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                                skeleton, spec, complex_gate_param, p, t,
                                make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                                    lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                                super_circuit, op_params);
                        }

                        if constexpr (is_sidb_defect_surface_v<CellLyt>)
                        {
                            if (defect_surface.has_value())
                            {
                                const auto skeleton_with_defects = add_defect_to_skeleton(
                                    defect_surface.value(), skeleton, params.influence_radius_charged_defects,
                                    center_cell, absolute_cell);
                                return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                                    skeleton_with_defects, std::vector<tt>{f}, params, p, t,
                                    make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                                        lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                                    super_circuit, op_params);
                            }
                        }

                        return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                            skeleton, std::vector<tt>{f}, params, p, t,
                            make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                                lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                            super_circuit, op_params);
                    }
                    return {fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::EMPTY_GATE};
                }
            }

            if constexpr (is_sidb_defect_surface_v<CellLyt>)
            {
                if (defect_surface.has_value())
                {
                    const auto skeleton_with_defects =
                        add_defect_to_skeleton(defect_surface.value(), skeleton,
                                               params.influence_radius_charged_defects, center_cell, absolute_cell);

                    return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                        skeleton_with_defects, std::vector<tt>{f}, params, p, t,
                        make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                            lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                        super_circuit, op_params);
                }
            }

            return design_gates<CellLyt, tt, CellLyt, GateLyt, ExtPotType, SkeletonGateLibrary>(
                skeleton, std::vector<tt>{f}, params, p, t,
                make_bdl_circuit_for_tile<CellLyt, GateLyt, SkeletonGateLibrary>(
                    lyt, t, op_params->input_bdl_iterator_params.bdl_wire_params),
                super_circuit, op_params);
        }

        catch (const std::out_of_range&)
        {
            throw unsupported_gate_orientation_exception(t, p);
        }

        throw unsupported_gate_type_exception(t);
    }

  private:
    template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
    [[nodiscard]] static sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>
    make_bdl_circuit_for_tile(const GateLyt& gate_lyt, const tile<GateLyt>& t,
                              const detect_bdl_wires_params& bdl_wire_params) noexcept
    {
        GateLyt gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};

        const mockturtle::node<GateLyt>& n = gate_lyt.get_node(t);

        std::vector<mockturtle::signal<GateLyt>> inputs_to_t{};

        for (const auto& in_t : gate_lyt.incoming_data_flow(t))
        {
            assert(gate_lyt_window.is_empty_tile(in_t) && "tile on which PI is to be created is already populated");
            // std::cout << "created pi at " << in_t.x << " " << in_t.y << " " << in_t.z << std::endl;
            gate_lyt_window.create_pi("", in_t);
            // std::cout << "num pis" << gate_lyt_window.num_pis() << std::endl;

            inputs_to_t.push_back(static_cast<mockturtle::signal<GateLyt>>(in_t));
        }

        // std::cout << "created node at " << t.x << " " << t.y << " " << t.z << std::endl;
        assert(gate_lyt_window.is_empty_tile(t) && "tile on which node is to be created is already populated");
        if (gate_lyt.is_pi(gate_lyt.get_node(t)))
        {
            assert(inputs_to_t.size() == 0 && "PI tile has inputs");
            gate_lyt_window.create_pi("", t);
        }
        else if (gate_lyt.is_po(gate_lyt.get_node(t)))
        {
            assert(inputs_to_t.size() != 1 && "PO tile has not precisely one input");
            gate_lyt_window.create_po(static_cast<mockturtle::signal<GateLyt>>(inputs_to_t.front()), "", t);
        }
        else
        {
            gate_lyt_window.create_node(inputs_to_t, gate_lyt.node_function(n), t);
        }

        for (const auto& out_t : gate_lyt.outgoing_data_flow(t))
        {
            assert(gate_lyt_window.is_empty_tile(out_t) && "tile on which PO is to be created is already populated");

            // std::cout << "created po at " << out_t.x << " " << out_t.y << " " << out_t.z << std::endl;
            gate_lyt_window.create_po(static_cast<mockturtle::signal<GateLyt>>(t), "", out_t);
        }

        if (gate_lyt.is_buf(n))
        {
            if (const auto above_t = gate_lyt.above(t); t != above_t && gate_lyt.is_wire_tile(above_t))
            {
                // upper_t hosts a crossing or double wire

                std::vector<mockturtle::signal<GateLyt>> inputs_to_above_t{};

                for (const auto& in_t : gate_lyt.incoming_data_flow(above_t))
                {
                    assert(gate_lyt_window.is_empty_tile(in_t) &&
                           "tile on which PI is to be created is already populated");
                    // std::cout << "created pi at " << in_t.x << " " << in_t.y << " " << in_t.z << std::endl;

                    gate_lyt_window.create_pi("", in_t);
                    // std::cout << "num pis" << gate_lyt_window.num_pis() << std::endl;

                    inputs_to_above_t.push_back(static_cast<mockturtle::signal<GateLyt>>(in_t));
                }

                assert(gate_lyt_window.is_empty_tile(above_t) &&
                       "tile on which node is to be created is already populated");
                // std::cout << "created node at " << above_t.x << " " << above_t.y << " " << above_t.z << std::endl;
                gate_lyt_window.create_node(inputs_to_above_t, gate_lyt.node_function(gate_lyt.get_node(above_t)),
                                            above_t);

                for (const auto& out_t : gate_lyt.outgoing_data_flow(above_t))
                {
                    assert(gate_lyt_window.is_empty_tile(out_t) &&
                           "tile on which PO is to be created is already populated");
                    // std::cout << "created po at " << out_t.x << " " << out_t.y << " " << out_t.z << std::endl;
                    gate_lyt_window.create_po(static_cast<mockturtle::signal<GateLyt>>(above_t), "", out_t);
                }
            }
        }

        return sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{gate_lyt_window, bdl_wire_params,
                                                                       std::make_optional(t)};
    }
    /**
     * This function designs an SiDB gate for a given Boolean function at a given tile and a given rotation. If
     atomic
     * defects exist, they are incorporated into the design process.
     *
     * An exception is thrown in case there is no possible gate design.
     *
     * @tparam LytSkeleton The cell-level layout of the skeleton.
     * @tparam TT Truth table type.
     * @tparam CellLyt The cell-level layout.
     * @tparam GateLyt The gate-level layout.
     * todo
     * @param skeleton Skeleton with atomic defects if available.
     * @param spec Expected Boolean function of the layout given as a multi-output truth table.
     * @param parameters Parameters for the SiDB gate design process.
     * @param p The list of ports and their directions.
     * @param tile The specific tile on which the gate should be designed.
     * @return An `fcn_gate` object.
     */
    template <typename LytSkeleton, typename TT, typename CellLyt, typename GateLyt,
              local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED,
              typename SkeletonGateLibrary             = sidb_on_the_fly_gate_library>
    [[nodiscard]] static std::vector<typename sidb_on_the_fly_gate_library::fcn_gate>
    design_gates(const LytSkeleton& skeleton, const std::vector<TT>& spec,
                 const sidb_on_the_fly_gate_library_params<CellLyt>& parameters, const port_list<port_direction>& p,
                 const tile<GateLyt>& tile, sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>&& circuit,
                 const std::optional<sidb_bdl_circuit<LytSkeleton, GateLyt, SkeletonGateLibrary>>& super_circuit,
                 const std::optional<is_circuit_operational_params>&                               op_params)
    {
        static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
        static_assert(has_cube_coord_v<CellLyt>, "CellLyt is not based on cube coordinates");

        const auto create_fcn_gates = [&](const auto& found_gate_layouts)
        {
            std::cout << "number of gate layouts found: " << found_gate_layouts.size() << std::endl;

            if (found_gate_layouts.empty())
            {
                throw gate_design_exception<tt, GateLyt>(tile, create_id_tt(), p);
            }

            std::vector<typename sidb_on_the_fly_gate_library::fcn_gate> gates{};
            gates.reserve(found_gate_layouts.size());

            designed_sidb_gates<LytSkeleton, ExtPotType> sorted_gates{};
            sorted_gates.gate_layouts.reserve(found_gate_layouts.size());

            for (const auto& gate : found_gate_layouts)
            {
                gates.emplace_back(cell_list_to_gate<char>(
                    cell_level_layout_to_list(std::move(gate), circuit.gate_layout, tile, false)));
                sorted_gates.gate_layouts.emplace_back(gate);
            }

            // order_designed_sidb_gates({std::make_shared<compare_by_minimum_ground_state_isolation<LytSkeleton>>(),
            //                            std::make_shared<compare_by_average_ground_state_isolation<LytSkeleton>>()},
            //                           sorted_gates);
            //
            // print_layout(sorted_gates.gate_layouts.front());

            return gates;
        };

        const auto params = is_sidb_gate_design_impossible_params{
            parameters.design_gate_params.operational_params.simulation_parameters};

        if constexpr (is_sidb_defect_surface_v<LytSkeleton>)
        {
            if (is_sidb_gate_design_impossible(skeleton, spec, params))
            {
                throw gate_design_exception<tt, GateLyt>(tile, spec.front(), p);
            }
        }

        std::cout << "starting gate design for tile " << tile << "\t|\tnode function:";
        for (const tt& tt : spec)
        {
            std::cout << '\t';
            kitty::print_binary(tt);
        }
        std::cout << std::endl;

        const auto found_gate_layouts = design_sidb_gates<LytSkeleton, TT, ExtPotType, GateLyt, SkeletonGateLibrary>(
            skeleton, spec, parameters.design_gate_params, nullptr, std::make_optional(std::move(circuit)),
            super_circuit, op_params);

        return create_fcn_gates(found_gate_layouts);
    }
    /**
     * This function takes a defect surface and a skeleton and adds defects from the surrounding area
     * to the skeleton. The defects within a specified distance from the center cell are taken into account.
     * The resulting skeleton with added defects is returned.
     *
     * @tparam CellLyt SiDB defect surface type.
     * @tparam Params Type of Parameters.
     * @param skeleton The skeleton to which defects will be added.
     * @param center_cell The coordinates of the center cell.
     * @param absolute_cell The coordinates of the skeleton's absolute cell.
     * @param parameters Parameters for defect handling.
     * @return The updated skeleton with added defects from the surrounding area.
     */
    template <typename CellLyt>
    [[nodiscard]] static CellLyt
    add_defect_to_skeleton(const CellLyt& defect_surface, const CellLyt& skeleton, const double influence_distance,
                           const cell<CellLyt>& center_cell, const cell<CellLyt>& absolute_cell)
    {
        static_assert(is_sidb_defect_surface_v<CellLyt>, "CellLyt is not a defect surface");
        static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
        static_assert(has_cube_coord_v<CellLyt>, "CellLyt is not based on cube coordinates");

        auto skeleton_with_defect = skeleton;

        defect_surface.foreach_sidb_defect(
            [&skeleton_with_defect, &center_cell, &absolute_cell, &influence_distance](const auto& cd)
            {
                // all defects (charged) in a distance of influence_radius_charged_defects from the center are taken
                // into account.
                if (sidb_nm_distance(CellLyt{}, center_cell, cd.first) < influence_distance)
                {
                    const auto relative_defect_position = cd.first - absolute_cell;
                    skeleton_with_defect.assign_sidb_defect(relative_defect_position, cd.second);
                }
            });

        const auto bb = bounding_box_2d(skeleton_with_defect);
        skeleton_with_defect.resize(bb.get_max());

        return skeleton_with_defect;
    }
    /**
     * Determines the port directions of a given tile.
     *
     * @tparam Lyt Pointy-top hexagonal gate-level layout type.
     * @param lyt Given tile `t` for which the port directions are determined.
     * @return port directions of the given tile are returned as `port_list`.
     */
    template <typename Lyt>
    [[nodiscard]] static port_list<port_direction> determine_port_routing(const Lyt& lyt, const tile<Lyt>& t) noexcept
    {
        port_list<port_direction> p{};

        // determine incoming connector ports
        if (lyt.has_north_eastern_incoming_signal(t))
        {
            p.inp.emplace(port_direction::cardinal::NORTH_EAST);
        }
        if (lyt.has_north_western_incoming_signal(t))
        {
            p.inp.emplace(port_direction::cardinal::NORTH_WEST);
        }

        // determine outgoing connector ports
        if (lyt.has_south_eastern_outgoing_signal(t))
        {
            p.out.emplace(port_direction::cardinal::SOUTH_EAST);
        }
        if (lyt.has_south_western_outgoing_signal(t))
        {
            p.out.emplace(port_direction::cardinal::SOUTH_WEST);
        }

        // gates without connector ports

        // 1-input functions
        if (const auto n = lyt.get_node(t); lyt.is_pi(n) || lyt.is_po(n) || lyt.is_buf(n) || lyt.is_inv(n))
        {
            if (lyt.has_no_incoming_signal(t))
            {
                p.inp.emplace(port_direction::cardinal::NORTH_WEST);
            }
            if (lyt.has_no_outgoing_signal(t))
            {
                p.out.emplace(port_direction::cardinal::SOUTH_EAST);
            }
        }
        else  // 2-input functions
        {
            if (lyt.has_no_incoming_signal(t))
            {
                p.inp.emplace(port_direction::cardinal::NORTH_WEST);
                p.inp.emplace(port_direction::cardinal::NORTH_EAST);
            }
            if (lyt.has_no_outgoing_signal(t))
            {
                p.out.emplace(port_direction::cardinal::SOUTH_EAST);
            }
        }

        return p;
    }

    using double_port_gate_function_map =
        phmap::flat_hash_map<std::pair<port_list<port_direction>, port_list<port_direction>>, std::vector<tt>>;

    static inline const double_port_gate_function_map TWO_IN_TWO_OUT_MAP = {
        {{{{port_direction(port_direction::cardinal::NORTH_WEST)},
           {port_direction(port_direction::cardinal::SOUTH_WEST)}},
          {{port_direction(port_direction::cardinal::NORTH_EAST)},
           {port_direction(port_direction::cardinal::SOUTH_EAST)}}},
         create_double_wire_tt()},
        {{{{port_direction(port_direction::cardinal::NORTH_EAST)},
           {port_direction(port_direction::cardinal::SOUTH_EAST)}},
          {{port_direction(port_direction::cardinal::NORTH_WEST)},
           {port_direction(port_direction::cardinal::SOUTH_WEST)}}},
         create_double_wire_tt()},
        {{{{port_direction(port_direction::cardinal::NORTH_WEST)},
           {port_direction(port_direction::cardinal::SOUTH_EAST)}},
          {{port_direction(port_direction::cardinal::NORTH_EAST)},
           {port_direction(port_direction::cardinal::SOUTH_WEST)}}},
         create_crossing_wire_tt()},
        {{{{port_direction(port_direction::cardinal::NORTH_EAST)},
           {port_direction(port_direction::cardinal::SOUTH_WEST)}},
          {{port_direction(port_direction::cardinal::NORTH_WEST)},
           {port_direction(port_direction::cardinal::SOUTH_EAST)}}},
         create_crossing_wire_tt()},
    };
};

}  // namespace fiction

#endif  // FICTION_SIDB_ON_THE_FLY_GATE_LIBRARY_HPP

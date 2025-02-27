//
// Created by Willem Lambooy on 18/02/2025.
//

#ifndef SKELETON_INFLUENCE_BOUNDS_HPP
#define SKELETON_INFLUENCE_BOUNDS_HPP

#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/io/print_layout.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/layout_utils.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiction
{

template <typename CellType>
struct skeleton_influence_bounds_params
{
    /**
     * The simulation parameters for the physical simulation of the ground state.
     */
    sidb_simulation_parameters simulation_parameters{};
    /**
     * Canvas spanned by the northwest and southeast cell.
     */
    std::pair<CellType, CellType> canvas{};
    /**
     * Parameters to detect BDL wires.
     */
    detect_bdl_wires_params bdl_wire_params{};
};

namespace detail
{

template <typename CellLyt, typename SkeletonGateLibrary, typename GateLyt>
class skeleton_influence_bounds_impl
{
  public:
    skeleton_influence_bounds_impl(const GateLyt& gate_layout, const tile<GateLyt>& t,
                                   const skeleton_influence_bounds_params<cell<CellLyt>>& ps) noexcept :

            params{ps},
            gate_lyt{gate_layout},
            current_tile{t},
            bdl_wires_of_designed_gates_current_tile_with_other_tile{
                obtain_bdl_wires_for_all_tile_pairs(gate_layout, t)}
    // designed_gate{apply_io_labels_to_bdl_pairs(gate_layout, t, ps)},
    // input_wires_of_designed_gate{detect_bdl_wires<CellLyt>(designed_gate, {}, bdl_wire_selection::INPUT)},
    // output_wires_of_designed_gate{detect_bdl_wires<CellLyt>(designed_gate, {}, bdl_wire_selection::OUTPUT)},
    {}

    [[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>> run() noexcept
    {
        const CellLyt designed_gate =
            apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile});

        // print_layout(designed_gate);

        designed_gate.foreach_cell(
            [&](const cell<CellLyt>& absolute_c)
            {
                const cell<CellLyt> absolute_offset =
                    relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                       SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                        gate_lyt, current_tile, {0, 0});

                // std::cout << absolute_c.x - absolute_offset.x << "," << absolute_c.y - absolute_offset.y << " @ ";

                collect_influence_bounds_for_sidb({absolute_c.x - absolute_offset.x, absolute_c.y - absolute_offset.y},
                                                  absolute_c);
                // if (sidb_is_in_io_wires(relative_c))
                // {
                // }
                // else
                // {
                //     collect_influence_bounds_for_sidb(relative_c);
                // }
            });

        for (const cell<CellLyt>& relative_c :
             all_coordinates_in_spanned_area(params.canvas.first, params.canvas.second))
        {
            const cell<CellLyt> absolute_c =
                relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                   SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                    gate_lyt, current_tile, relative_c);
            collect_influence_bounds_for_sidb(relative_c, absolute_c);
        }

        return skeleton_influence_bounds_map;
    }

  private:
    // static [[nodiscard]] CellLyt
    // make_designed_gate_with_io_marked_bdl_wires(const GateLyt& gate_layout, const tile<GateLyt>& t,
    //                                             const skeleton_influence_bounds_params<cell<CellLyt>>& ps) noexcept
    // {
    //     CellLyt designed_gate = apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout, {t}, ps);
    //
    //     designed_gate.foreach_cell(
    //         [&](cell<CellLyt>& sidb)
    //         {
    //             if (sidb.x >= ps.canvas.first.x && sidb.x <= ps.canvas.second.x && sidb.y >= ps.canvas.first.y &&
    //                 sidb.y <= ps.canvas.second.y)
    //             {
    //                 return;
    //             }
    //
    //             // this code assumes north to south row clocking with hexagonal tiling
    //             if (sidb.y < ps.canvas.first.y)
    //             {
    //                 sidb.assign_cell_type(sidb_technology::cell_type::INPUT);
    //             }
    //
    //             sidb.assign_cell_type(sidb_technology::cell_type::OUTPUT);
    //         });
    //
    //     return designed_gate;
    // }
    //
    // struct sidb_identification
    // {
    //     enum class sidb_role : uint8_t
    //     {
    //         CANVAS,
    //         INPUT_WIRE,
    //         OUTPUT_WIRE
    //     };
    //
    //     cell<CellLyt> relative_position;
    //     cell<CellLyt> absolute_position;
    //     sidb_role     role;
    //
    //     sidb_identification(const skeleton_influence_bounds_impl* parent_class, const GateLyt& gate_lyt,
    //                         const tile<GateLyt>& tile, const cell<CellLyt>& sidb) noexcept :
    //             relative_position{sidb},
    //             absolute_position{
    //                 relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
    //                                                    SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
    //                     gate_lyt, tile, relative_position)},
    //             role{determine_sidb_role(parent_class, sidb)}
    //     {}
    //
    //     static [[nodiscard]] sidb_role determine_sidb_role(const skeleton_influence_bounds_impl* parent_class)
    //     noexcept
    //     {
    //         for (const bdl_wire<CellLyt>& wire : parent_class->input_wires_of_designed_gate)
    //         {
    //             for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
    //             {
    //                 if (pair.lower == relative_position || pair.upper == relative_position)
    //                 {
    //                     return sidb_role::INPUT_WIRE;
    //                 }
    //             }
    //         }
    //
    //         for (const bdl_wire<CellLyt>& wire : parent_class->output_wires_of_designed_gate)
    //         {
    //             for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
    //             {
    //                 if (pair.lower == relative_position || pair.upper == relative_position)
    //                 {
    //                     return sidb_role::OUTPUT_WIRE;
    //                 }
    //             }
    //         }
    //
    //         return sidb_role::CANVAS;
    //     }
    // };

    // void collect_influence_bounds_for_sidb(const sidb_identification& sidb_id) noexcept
    void collect_influence_bounds_for_sidb(const cell<CellLyt>& relative_sidb,
                                           const cell<CellLyt>& absolute_sidb) noexcept
    {
        // std::cout << absolute_sidb.x << " " << absolute_sidb.y << std::endl;

        std::array<double, 2> bounds{{0, 0}};

        uint64_t node_counter = 0;

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (gate_lyt.is_constant(n))
                {
                    return;
                }

                const tile<GateLyt>& t = gate_lyt.get_tile(n);

                if (t == current_tile)
                {
                    return;
                }

                for (bdl_wire<CellLyt>& wire :
                     detect_bdl_wires(apply_gate_library<CellLyt, SkeletonGateLibrary>(gate_lyt, std::set{t})))
                {
                    std::array<double, 2> potential_from_wire{{0, 0}};

                    for (bdl_pair<cell<CellLyt>>& pair : wire.pairs)
                    {
                        CellLyt cell_lyt_with_c_and_bound_pair{};

                        cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb,
                                                                        sidb_technology::cell_type::NORMAL);

                        make_pair_sidbs_respect_io_connections(
                            absolute_sidb, bdl_wires_of_designed_gates_current_tile_with_other_tile.at(node_counter),
                            pair);

                        cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);
                        cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);

                        const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                                       params.simulation_parameters};

                        potential_from_wire[0] +=
                            cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.lower);
                        potential_from_wire[1] +=
                            cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.upper);
                    }

                    std::sort(potential_from_wire.begin(), potential_from_wire.end());

                    bounds[0] += potential_from_wire[0];
                    bounds[1] += potential_from_wire[1];
                }

                node_counter++;
            });

        skeleton_influence_bounds_map.insert({relative_sidb, std::move(bounds)});
    }

    void make_pair_sidbs_respect_io_connections(const cell<CellLyt>&                  absolute_sidb,
                                                const std::vector<bdl_wire<CellLyt>>& bdl_wires_in_gates_together,
                                                bdl_pair<cell<CellLyt>>&              pair) const noexcept
    {
        // if (sidb_id.role == sidb_identification::sidb_role::CANVAS)
        // {
        //     return;
        // }
        //
        // // the code below assumes north to south row clocking with hexagonal tiling
        //
        // if (sidb_id.role == sidb_identification::sidb_role::INPUT_WIRE)
        // {
        //     if (gate_lyt.north_west(tile) != t && gate_lyt.north_east(tile) != t)
        //     {
        //         return;
        //     }
        //
        //     output_wires_of_other_designed_gate{detect_bdl_wires<CellLyt>(designed_gate, {},
        //     bdl_wire_selection::OUTPUT)},
        //
        // }
        //
        // if (gate_lyt.south_west(tile) != t && gate_lyt.south_east(tile) != t)
        // {
        //     return;
        // }
        //
        // input_wires_of_other_designed_gate{detect_bdl_wires<CellLyt>(designed_gate, {}, bdl_wire_selection::INPUT)},
        //
        //
        // if (t.north_west)
        // {
        //     lyt.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);
        //     lyt.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);
        // }

        // check if the biggest wire is unique, otherwise the gates are not in connection
        if (bdl_wires_in_gates_together.front().pairs.size() == bdl_wires_in_gates_together.at(1).pairs.size() ||
            !is_upper_in_wire(pair.upper, bdl_wires_in_gates_together.front()))
        {
            return;
        }

        if (!is_upper_in_wire(absolute_sidb, bdl_wires_in_gates_together.front()))
        {
            if (!is_lower_in_wire(absolute_sidb, bdl_wires_in_gates_together.front()))
            {
                return;
            }

            pair.lower = pair.upper;
        }

        if (!is_lower_in_wire(absolute_sidb, bdl_wires_in_gates_together.front()))
        {
            return;
        }

        pair.upper = pair.lower;
    }
    const skeleton_influence_bounds_params<cell<CellLyt>>& params;
    GateLyt                                                gate_lyt{};
    const tile<GateLyt>&                                   current_tile{};
    const std::vector<std::vector<bdl_wire<CellLyt>>>      bdl_wires_of_designed_gates_current_tile_with_other_tile{};
    // const std::vector<bdl_wire<CellLyt>>                     input_wires_of_designed_gate{};
    // const std::vector<bdl_wire<CellLyt>>                     output_wires_of_designed_gate{};
    std::unordered_map<cell<CellLyt>, std::array<double, 2>> skeleton_influence_bounds_map{};

    [[nodiscard]] static std::vector<std::vector<bdl_wire<CellLyt>>>
    obtain_bdl_wires_for_all_tile_pairs(const GateLyt& gate_lyt, const tile<GateLyt>& current_tile) noexcept
    {
        std::vector<std::vector<bdl_wire<CellLyt>>> bdl_wires_of_designed_gates_current_tile_with_other_tile{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (gate_lyt.is_constant(n))
                {
                    return;
                }

                const tile<GateLyt>& t = gate_lyt.get_tile(n);

                if (t == current_tile)
                {
                    return;
                }

                CellLyt designed_gates_together =
                    apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile, t});

                // print_layout(designed_gates_together);
                // std::cout << "\n";

                std::vector<bdl_wire<CellLyt>> bdl_wires_of_designed_gates_together =
                    detect_bdl_wires(designed_gates_together);

                // sort by wire length, first will be biggest
                std::sort(bdl_wires_of_designed_gates_together.begin(), bdl_wires_of_designed_gates_together.end(),
                          [](const bdl_wire<CellLyt>& lhs, const bdl_wire<CellLyt>& rhs)
                          { return lhs.pairs.size() < rhs.pairs.size(); });

                bdl_wires_of_designed_gates_current_tile_with_other_tile.push_back(
                    std::move(bdl_wires_of_designed_gates_together));
            });

        return bdl_wires_of_designed_gates_current_tile_with_other_tile;
    }

    [[nodiscard]] static bool is_upper_in_wire(const cell<CellLyt>& sidb, const bdl_wire<CellLyt>& wire) noexcept
    {

        for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
        {
            if (pair.upper == sidb)
            {
                return true;
            }
        }

        return false;
    }

    [[nodiscard]] static bool is_lower_in_wire(const cell<CellLyt>& sidb, const bdl_wire<CellLyt>& wire) noexcept
    {

        for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
        {
            if (pair.lower == sidb)
            {
                return true;
            }
        }

        return false;
    }
};

}  // namespace detail

template <typename CellLyt, typename SkeletonGateLibrary, typename GateLyt>
[[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>>
skeleton_influence_bounds(const GateLyt& gate_layout, const tile<GateLyt>& tile,
                          const skeleton_influence_bounds_params<cell<CellLyt>>& params = {}) noexcept
{
    return detail::skeleton_influence_bounds_impl<CellLyt, SkeletonGateLibrary, GateLyt>{gate_layout, tile, params}
        .run();
}

}  // namespace fiction

#endif  // SKELETON_INFLUENCE_BOUNDS_HPP

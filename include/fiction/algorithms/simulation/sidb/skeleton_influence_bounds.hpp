//
// Created by Willem Lambooy on 18/02/2025.
//

#ifndef SKELETON_INFLUENCE_BOUNDS_HPP
#define SKELETON_INFLUENCE_BOUNDS_HPP

#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/layout_utils.hpp"

#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>

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
    std::pair<CellType, CellType> canvas_complex_gates{};
    /**
     * Parameters to detect BDL wires.
     */
    detect_bdl_wires_params bdl_wire_params{};
    bool                    absolute_positions = false;
};

namespace detail
{

template <typename CellLyt, typename SkeletonGateLibrary, typename GateLyt>
class skeleton_influence_bounds_impl
{
  public:
    skeleton_influence_bounds_impl(const GateLyt& gate_layout, const std::vector<tile<GateLyt>>& ts,
                                   const skeleton_influence_bounds_params<cell<CellLyt>>& ps) noexcept :

            params{ps},
            gate_lyt{gate_layout},
            tiles_of_interest{ts},
            cell_lyt_with_all_skeletons{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt)},
            cell_lyt_of_tiles_of_interest{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(
                gate_lyt, std::set(tiles_of_interest.cbegin(), tiles_of_interest.cend()))},
            bdl_wires_of_skeleton_with_other_tile_for_all_tiles_of_interest{
                obtain_bdl_wires_for_all_tile_pairs(gate_lyt, tiles_of_interest)}
    {}

    [[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>> run() noexcept
    {
        uint64_t tile_counter = 0;
        for (const tile<GateLyt>& current_tile : tiles_of_interest)
        {
            bool first_one_passed = false;
            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (first_one_passed)
                    {
                        return;
                    }

                    if (!skip_physical_design_for_node(gate_lyt, n) && tiles_of_interest.size() == 1)
                    {
                        if (tiles_of_interest.size() == 1 && tiles_of_interest.front() == gate_lyt.get_tile(n))
                        {
                            std::cout << "Skeleton looks like:" << std::endl;
                            CellLyt lyt = apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt);
                            gate_lyt.foreach_node(
                                [&](const auto& nn)
                                {
                                    if (skip_physical_design_for_node(gate_lyt, nn))
                                    {
                                        return;
                                    }
                                    const auto& t = gate_lyt.get_tile(nn);
                                    const auto  canvas =
                                        is_complex_gate<GateLyt>(gate_lyt, nn) ?
                                             params.canvas_complex_gates :
                                             //      make_gate_design_params_for_complex_gates<design_sidb_gates_params<CellLyt>,
                                            //                                                CellLyt>()
                                            //         .canvas :
                                            params.canvas;

                                    for (const cell<CellLyt>& relative_c :
                                         all_coordinates_in_spanned_area(canvas.first, canvas.second))
                                    {
                                        const cell<CellLyt> absolute_c =
                                            relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                                               SkeletonGateLibrary::gate_y_size(),
                                                                               GateLyt, CellLyt>(gate_lyt, t,
                                                                                                 relative_c);
                                        lyt.assign_cell_type(absolute_c, sidb_technology::cell_type::LOGIC);
                                    }
                                });

                            print_layout(lyt);
                            std::cout << std::endl;
                        }

                        first_one_passed = true;
                    }
                });

            const CellLyt designed_gate =
                apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile});

            const charge_distribution_surface<CellLyt> cds{designed_gate};

            designed_gate.foreach_cell(
                [&](const cell<CellLyt>& absolute_c)
                {
                    const cell<CellLyt> absolute_offset =
                        relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                           SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                            gate_lyt, current_tile, {0, 0});

                    collect_influence_bounds_for_sidb(
                        {absolute_c.x - absolute_offset.x, absolute_c.y - absolute_offset.y}, absolute_c, tile_counter);
                });

            const auto canvas = is_complex_gate<GateLyt>(gate_lyt, gate_lyt.get_node(current_tile)) ?
                                    params.canvas_complex_gates :
                                    params.canvas;

            for (const cell<CellLyt>& relative_c : all_coordinates_in_spanned_area(canvas.first, canvas.second))
            {
                const cell<CellLyt> absolute_c =
                    relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                       SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                        gate_lyt, current_tile, relative_c);
                collect_influence_bounds_for_sidb(relative_c, absolute_c, tile_counter);
            }

            tile_counter++;
        }

        return skeleton_influence_bounds_map;
    }

  private:
    const skeleton_influence_bounds_params<cell<CellLyt>>& params;
    GateLyt                                                gate_lyt{};
    const std::vector<tile<GateLyt>>&                      tiles_of_interest{};
    const CellLyt                                          cell_lyt_with_all_skeletons{};
    const CellLyt                                          cell_lyt_of_tiles_of_interest{};
    const std::vector<std::vector<std::vector<bdl_wire<CellLyt>>>>
        bdl_wires_of_skeleton_with_other_tile_for_all_tiles_of_interest{};
    std::unordered_map<cell<CellLyt>, std::array<double, 2>> skeleton_influence_bounds_map{};

    void collect_influence_bounds_for_sidb(const cell<CellLyt>& relative_sidb, const cell<CellLyt>& absolute_sidb,
                                           const uint64_t tile_of_interest_counter) noexcept
    {
        std::array<double, 2> bounds{{0, 0}};

        uint64_t node_counter = 0;

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const tile<GateLyt>& t = gate_lyt.get_tile(n);

                if (std::find(tiles_of_interest.cbegin(), tiles_of_interest.cend(), t) != tiles_of_interest.cend())
                {
                    return;
                }

                const CellLyt cell_lyt_of_other_tile =
                    apply_gate_library<CellLyt, SkeletonGateLibrary>(gate_lyt, std::set{t});

                for (bdl_wire<CellLyt>& wire : detect_bdl_wires(cell_lyt_of_other_tile, params.bdl_wire_params))
                {
                    std::array<double, 2> potential_from_wire{{0, 0}};

                    std::vector<cell<CellLyt>> w1{}, w2{};

                    for (bdl_pair<cell<CellLyt>>& pair : wire.pairs)
                    {
                        if (pair.type == sidb_technology::cell_type::INPUT ||
                            is_contained_in_tile_of_interest(pair.lower))
                        {
                            continue;
                        }

                        CellLyt cell_lyt_with_c_and_bound_pair{};

                        cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb,
                                                                        sidb_technology::cell_type::NORMAL);

                        // todo ENABLE ... OR NOT ???
                        // make_pair_sidbs_respect_io_connections(
                        //     absolute_sidb,
                        //     bdl_wires_of_skeleton_with_other_tile_for_all_tiles_of_interest.at(tile_of_interest_counter)
                        //         .at(node_counter),
                        //     pair);

                        cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);
                        cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);

                        const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                                       params.simulation_parameters};

                        potential_from_wire[0] +=
                            -cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.lower);

                        // the upper SiDB in a pair can be an output perturber, we need to check for it
                        if (!is_output_perturber_in_tile_of_interest(pair.upper))
                        {
                            potential_from_wire[1] +=
                                -cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.upper);
                        }
                    }

                    std::sort(potential_from_wire.begin(), potential_from_wire.end());

                    bounds[0] += potential_from_wire[0];
                    bounds[1] += potential_from_wire[1];
                }

                node_counter++;
            });

        // include influence from the output perturbers if they are part of a different tile
        cell_lyt_with_all_skeletons.foreach_cell(
            [&](const cell<CellLyt>& c)
            {
                if (cell_lyt_with_all_skeletons.get_cell_type(c) == sidb_technology::OUTPUT_PERTURBER &&
                    !is_contained_in_tile_of_interest(c))
                {
                    CellLyt cell_lyt_with_c_and_bound_pair{};

                    cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb, sidb_technology::cell_type::NORMAL);

                    cell_lyt_with_c_and_bound_pair.assign_cell_type(c, sidb_technology::cell_type::NORMAL);

                    const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                                   params.simulation_parameters};

                    bounds[0] += -cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, c);
                    bounds[1] += -cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, c);
                }
            });

        // include influence from the input perturbers to the whole layout
        for (bdl_wire<CellLyt>& wire :
             detect_bdl_wires(cell_lyt_with_all_skeletons, params.bdl_wire_params, bdl_wire_selection::INPUT))
        {
            std::array<double, 2> potential_from_wire{{0, 0}};

            for (bdl_pair<cell<CellLyt>>& pair : wire.pairs)
            {
                if (cell_lyt_with_all_skeletons.get_cell_type(pair.upper) != sidb_technology::cell_type::INPUT ||
                    is_contained_in_tile_of_interest(pair.lower))
                {
                    continue;
                }

                CellLyt cell_lyt_with_c_and_bound_pair{};

                cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb, sidb_technology::cell_type::NORMAL);

                cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);
                cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);

                const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                               params.simulation_parameters};

                potential_from_wire[0] +=
                    -cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.lower);
                potential_from_wire[1] +=
                    -cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.upper);
            }

            std::sort(potential_from_wire.begin(), potential_from_wire.end());

            bounds[0] += potential_from_wire[0];
            bounds[1] += potential_from_wire[1];
        }

        skeleton_influence_bounds_map.insert(
            {params.absolute_positions ? absolute_sidb : relative_sidb, std::move(bounds)});
    }

    bool is_contained_in_tile_of_interest(const cell<CellLyt>& absolute_sidb) noexcept
    {
        bool found = false;

        cell_lyt_of_tiles_of_interest.foreach_cell(
            [&](const cell<CellLyt>& sidb_in_tile_of_interest)
            {
                if (!found && absolute_sidb == sidb_in_tile_of_interest &&
                    cell_lyt_of_tiles_of_interest.get_cell_type(sidb_in_tile_of_interest) !=
                        sidb_technology::cell_type::EMPTY)
                {
                    found = true;
                }
            });

        return found;
    }

    bool is_output_perturber_in_tile_of_interest(const cell<CellLyt>& absolute_sidb) noexcept
    {
        bool found = false;

        cell_lyt_of_tiles_of_interest.foreach_cell(
            [&](const cell<CellLyt>& sidb_in_tile_of_interest)
            {
                if (!found && absolute_sidb == sidb_in_tile_of_interest &&
                    cell_lyt_of_tiles_of_interest.get_cell_type(sidb_in_tile_of_interest) ==
                        sidb_technology::cell_type::OUTPUT_PERTURBER)
                {
                    found = true;
                }
            });

        return found;
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

    void make_pair_sidbs_respect_io_connections(const cell<CellLyt>&                  absolute_sidb,
                                                const std::vector<bdl_wire<CellLyt>>& bdl_wires_in_gates_together,
                                                bdl_pair<cell<CellLyt>>&              pair) const noexcept
    {
        // check if the biggest wire is unique, otherwise the gates are not in connection
        if (bdl_wires_in_gates_together.front().pairs.size() == bdl_wires_in_gates_together.at(1).pairs.size() ||
            !is_upper_in_wire(pair.upper, bdl_wires_in_gates_together.front()))
        {
            return;
        }

        if (is_upper_in_wire(absolute_sidb, bdl_wires_in_gates_together.front()))
        {
            pair.lower = pair.upper;

            return;
        }

        if (!is_lower_in_wire(absolute_sidb, bdl_wires_in_gates_together.front()))
        {
            return;
        }

        pair.upper = pair.lower;
    }

    [[nodiscard]] static std::vector<std::vector<std::vector<bdl_wire<CellLyt>>>>
    obtain_bdl_wires_for_all_tile_pairs(const GateLyt&                    gate_lyt,
                                        const std::vector<tile<GateLyt>>& tiles_of_interest) noexcept
    {
        std::vector<std::vector<std::vector<bdl_wire<CellLyt>>>>
            bdl_wires_of_skeleton_with_other_tile_for_all_tiles_of_interest{};

        for (const tile<GateLyt>& current_tile : tiles_of_interest)
        {
            std::vector<std::vector<bdl_wire<CellLyt>>> bdl_wires_of_skeleton_current_tile_with_other_tile{};

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (gate_lyt.is_constant(n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = gate_lyt.get_tile(n);

                    if (std::find(tiles_of_interest.cbegin(), tiles_of_interest.cend(), t) != tiles_of_interest.cend())
                    {
                        return;
                    }

                    CellLyt designed_gates_together =
                        apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile, t});

                    std::vector<bdl_wire<CellLyt>> bdl_wires_of_designed_gates_together =
                        detect_bdl_wires(designed_gates_together);

                    // sort by wire length, first will be biggest
                    std::sort(bdl_wires_of_designed_gates_together.begin(), bdl_wires_of_designed_gates_together.end(),
                              [](const bdl_wire<CellLyt>& lhs, const bdl_wire<CellLyt>& rhs)
                              { return lhs.pairs.size() > rhs.pairs.size(); });

                    bdl_wires_of_skeleton_current_tile_with_other_tile.push_back(
                        std::move(bdl_wires_of_designed_gates_together));
                });

            bdl_wires_of_skeleton_with_other_tile_for_all_tiles_of_interest.push_back(
                std::move(bdl_wires_of_skeleton_current_tile_with_other_tile));
        }

        return bdl_wires_of_skeleton_with_other_tile_for_all_tiles_of_interest;
    }
};

}  // namespace detail

template <typename CellLyt, typename SkeletonGateLibrary, typename GateLyt>
[[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>>
skeleton_influence_bounds(const GateLyt& gate_layout, const std::vector<tile<GateLyt>>& tile,
                          const skeleton_influence_bounds_params<cell<CellLyt>>& params = {}) noexcept
{
    return detail::skeleton_influence_bounds_impl<CellLyt, SkeletonGateLibrary, GateLyt>{gate_layout, tile, params}
        .run();
}

}  // namespace fiction

#endif  // SKELETON_INFLUENCE_BOUNDS_HPP

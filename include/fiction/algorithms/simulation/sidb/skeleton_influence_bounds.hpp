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
                gate_lyt, std::set(tiles_of_interest.cbegin(), tiles_of_interest.cend()))}
    {
    }

    [[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>> run() noexcept
    {
        for (const tile<GateLyt>& current_tile : tiles_of_interest)
        {
            if (current_tile.y == 0 && tiles_of_interest.size() == 1)
            {
                std::cout << "Skeleton looks like:" << std::endl;
                print_layout(apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt));
                std::cout << std::endl;
            }

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
                        {absolute_c.x - absolute_offset.x, absolute_c.y - absolute_offset.y}, absolute_c);
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
        }

        return skeleton_influence_bounds_map;
    }

  private:
    void collect_influence_bounds_for_sidb(const cell<CellLyt>& relative_sidb,
                                           const cell<CellLyt>& absolute_sidb) noexcept
    {
        std::array<double, 2> bounds{{0, 0}};

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

    const skeleton_influence_bounds_params<cell<CellLyt>>&   params;
    GateLyt                                                  gate_lyt{};
    const std::vector<tile<GateLyt>>&                        tiles_of_interest{};
    const CellLyt                                            cell_lyt_with_all_skeletons{};
    const CellLyt                                            cell_lyt_of_tiles_of_interest{};
    std::unordered_map<cell<CellLyt>, std::array<double, 2>> skeleton_influence_bounds_map{};
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

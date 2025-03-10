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
    {}

    [[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>> run() noexcept
    {
        if (current_tile.y == 0)
        {
            std::cout << "Skeleton looks like:" << std::endl;
            print_layout(apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt));
            std::cout << std::endl;
        }

        const CellLyt designed_gate =
            apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile});

        designed_gate.foreach_cell(
            [&](const cell<CellLyt>& absolute_c)
            {
                const cell<CellLyt> absolute_offset =
                    relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                       SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                        gate_lyt, current_tile, {0, 0});

                collect_influence_bounds_for_sidb({absolute_c.x - absolute_offset.x, absolute_c.y - absolute_offset.y},
                                                  absolute_c);
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
    void collect_influence_bounds_for_sidb(const cell<CellLyt>& relative_sidb,
                                           const cell<CellLyt>& absolute_sidb) noexcept
    {
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
    const skeleton_influence_bounds_params<cell<CellLyt>>&   params;
    GateLyt                                                  gate_lyt{};
    const tile<GateLyt>&                                     current_tile{};
    const std::vector<std::vector<bdl_wire<CellLyt>>>        bdl_wires_of_designed_gates_current_tile_with_other_tile{};
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

                std::vector<bdl_wire<CellLyt>> bdl_wires_of_designed_gates_together =
                    detect_bdl_wires(designed_gates_together);

                // sort by wire length, first will be biggest
                std::sort(bdl_wires_of_designed_gates_together.begin(), bdl_wires_of_designed_gates_together.end(),
                          [](const bdl_wire<CellLyt>& lhs, const bdl_wire<CellLyt>& rhs)
                          { return lhs.pairs.size() > rhs.pairs.size(); });

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

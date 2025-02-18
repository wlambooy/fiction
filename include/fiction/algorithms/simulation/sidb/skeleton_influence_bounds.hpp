//
// Created by Willem Lambooy on 18/02/2025.
//

#ifndef SKELETON_INFLUENCE_BOUNDS_HPP
#define SKELETON_INFLUENCE_BOUNDS_HPP

#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/layout_utils.hpp"

#include <algorithm>
#include <array>
#include <unordered_map>

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
            gate_lyt_copy{gate_layout},
            tile{t},
            params{ps}
    {
        gate_lyt_copy.clear_tile(t);
    }

    [[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>> run() const noexcept
    {
        std::unordered_map<cell<CellLyt>, std::array<double, 2>> skeleton_influence_bounds_map{};

        for (const cell<CellLyt>& relative_c :
             all_coordinates_in_spanned_area(params.canvas.first, params.canvas.second))
        {
            const cell<CellLyt> absolute_c{
                relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                   SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                    gate_lyt_copy, tile, relative_c)};

            std::array<double, 2> bounds{{0, 0}};

            for (const bdl_wire<CellLyt>& wire :
                 detect_bdl_wires(apply_gate_library<CellLyt, SkeletonGateLibrary>(gate_lyt_copy)))
            {
                std::array<double, 2> potential_from_wire{{0, 0}};

                for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
                {
                    CellLyt cell_lyt_with_c_and_pair{};
                    cell_lyt_with_c_and_pair.assign_cell_type(absolute_c, sidb_technology::cell_type::NORMAL);
                    cell_lyt_with_c_and_pair.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);
                    cell_lyt_with_c_and_pair.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);

                    const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_pair,
                                                                                   params.simulation_parameters};

                    potential_from_wire[0] +=
                        cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_c, pair.lower);
                    potential_from_wire[1] +=
                        cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_c, pair.upper);
                }

                std::sort(potential_from_wire.begin(), potential_from_wire.end());

                bounds[0] += potential_from_wire[0];
                bounds[1] += potential_from_wire[1];
            }

            skeleton_influence_bounds_map.insert({relative_c, std::move(bounds)});
        }

        return skeleton_influence_bounds_map;
    }

  private:
    GateLyt                                                gate_lyt_copy{};
    const tile<GateLyt>&                                   tile{};
    const skeleton_influence_bounds_params<cell<CellLyt>>& params;
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

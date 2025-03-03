//
// Created by Willem Lambooy on 28/02/2025.
//

#include "fiction/traits.hpp"
#include "mockturtle/traits.hpp"

#ifndef GATE_DESIGN_UTILS_HPP
#define GATE_DESIGN_UTILS_HPP

namespace fiction
{

/**
 * This function assigns a given FCN gate implementation to the total cell layout.
 *
 * @tparam CellLyt Type of the returned cell-level layout.
 * @tparam GateLibrary Type of the gate library to apply.
 * @tparam GateLyt Type of the gate-level layout to apply the library to.
 * @param c Top-left cell of the tile where the gate is placed.
 * @param g Gate implementation.
 * @param n Corresponding node in the gate-level layout.
 */
template <typename CellLyt, typename GateLibrary, typename GateLyt>
void assign_gate(CellLyt& cell_lyt, const cell<CellLyt>& c, const typename GateLibrary::fcn_gate& g,
                 const GateLyt& gate_lyt, const mockturtle::node<GateLyt>& n)
{
    const auto start_x = c.x;
    const auto start_y = c.y;
    const auto layer   = c.z;

    for (auto y = 0ul; y < g.size(); ++y)
    {
        for (auto x = 0ul; x < g[y].size(); ++x)
        {
            const cell<CellLyt> pos{start_x + x, start_y + y, layer};
            const auto          type{g[y][x]};

            if (!technology<CellLyt>::is_empty_cell(type))
            {
                cell_lyt.assign_cell_type(pos, type);
            }

            // set IO names
            if (technology<CellLyt>::is_input_cell(type) || technology<CellLyt>::is_output_cell(type))
            {
                cell_lyt.assign_cell_name(pos, gate_lyt.get_name(n));
            }
        }
    }
}

}  // namespace fiction

#endif  // GATE_DESIGN_UTILS_HPP

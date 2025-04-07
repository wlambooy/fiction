//
// Created by Willem Lambooy on 28/02/2025.
//

#ifndef GATE_DESIGN_UTILS_HPP
#define GATE_DESIGN_UTILS_HPP

#include "fiction/traits.hpp"
#include "mockturtle/traits.hpp"

#include <cstdint>

namespace fiction
{
/**
 * This function assigns a given FCN gate implementation to the total cell layout.
 *
 * @tparam CellLyt Type of the returned cell-level layout.
 * @tparam GateLibrary Type of the gate library to apply.
 * @tparam GateLyt Type of the gate-level layout to apply the library to.
 * todo
 * @param c Top-left cell of the tile where the gate is placed.
 * @param g Gate implementation.
 * @param n Corresponding node in the gate-level layout.
 */
template <typename CellLyt, typename GateLibrary, typename GateLyt>
static void assign_gate(CellLyt& cell_lyt, const cell<CellLyt>& c, const typename GateLibrary::fcn_gate& g,
                        const GateLyt& gate_lyt, const mockturtle::node<GateLyt>& n)
{
    if (gate_lyt.is_pi(n) || gate_lyt.is_po(n))
    {
        return;
    }

    const auto start_x = c.x;
    const auto start_y = c.y;
    const auto layer   = c.z;

    for (decltype(c.y) y = 0ul; static_cast<uint64_t>(y) < g.size(); ++y)
    {
        for (decltype(c.x) x = 0ul; static_cast<uint64_t>(x) < g[static_cast<uint64_t>(y)].size(); ++x)
        {
            const cell<CellLyt> pos{start_x + x, start_y + y, layer};
            const auto          type{g[static_cast<uint64_t>(y)][static_cast<uint64_t>(x)]};

            if (technology<CellLyt>::is_empty_cell(type))
            {
                continue;
            }

            // overwrites always make a NORMAL type cell
            if (!technology<CellLyt>::is_empty_cell(cell_lyt.get_cell_type(pos)))
            {
                cell_lyt.assign_cell_type(pos, sidb_technology::cell_type::NORMAL);

                continue;
            }

            cell_lyt.assign_cell_type(pos, type);

            // set IO names
            if (technology<CellLyt>::is_input_cell(type) || technology<CellLyt>::is_output_cell(type))
            {
                cell_lyt.assign_cell_name(pos, gate_lyt.get_name(n));
            }
        }
    }
}

template <typename GateLyt>
static bool skip_physical_design_for_node(const GateLyt& gate_lyt, const mockturtle::node<GateLyt>& n) noexcept
{
    return gate_lyt.is_constant(n) || gate_lyt.is_pi(n) || gate_lyt.is_po(n) || gate_lyt.is_constant(n) ||
           (gate_lyt.is_buf(n) && !gate_lyt.is_ground_layer(gate_lyt.get_tile(n)));
}

template <typename GateLyt>
static bool is_complex_gate(const GateLyt& gate_lyt, const mockturtle::node<GateLyt>& n) noexcept
{
    const auto t  = gate_lyt.get_tile(n);
    const auto at = gate_lyt.above(t);

    return gate_lyt.is_buf(n) && t != at && gate_lyt.is_wire_tile(at);
}

template <typename Params, typename CellLyt>
static Params make_gate_design_params_for_complex_gates(const uint64_t canvas_sidb_complex_gates = 0) noexcept
{
    Params design_gate_params{};
    design_gate_params.number_of_sidbs = canvas_sidb_complex_gates;
    design_gate_params.design_mode     = design_sidb_gates_params<CellLyt>::design_sidb_gates_mode::RANDOM;
    design_gate_params.canvas          = {{10, 9}, {24, 19}};
    // design_gate_params.operational_params.op_condition_kinks =
    // is_operational_params<cell<CellLyt>>::operational_condition_kinks::TOLERATE_KINKS;
    design_gate_params.operational_params.op_condition_kinks =
        is_operational_params<cell<CellLyt>>::operational_condition_kinks::REJECT_KINKS;
    design_gate_params.operational_params.strategy_to_analyze_operational_status =
        is_operational_params<cell<CellLyt>>::operational_analysis_strategy::FILTER_THEN_SIMULATION;

    return design_gate_params;
}

}  // namespace fiction

#endif  // GATE_DESIGN_UTILS_HPP

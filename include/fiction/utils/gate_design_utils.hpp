//
// Created by Willem Lambooy on 28/02/2025.
//

#ifndef GATE_DESIGN_UTILS_HPP
#define GATE_DESIGN_UTILS_HPP

#include "fiction/traits.hpp"
#include "fiction/utils/layout_utils.hpp"
#include "mockturtle/traits.hpp"

#include <fiction/layouts/bounding_box.hpp>

#include <cstdint>

namespace fiction
{

/**
 * This exception is thrown when an error occurs during the design of an SiDB gate.
 * It provides information about the tile, truth table, and port list associated with the error.
 *
 * @tparam TT The type representing the truth table.
 * @tparam GateLyt The type representing the gate-level layout.
 */
template <typename TT, typename GateLyt>
class gate_design_exception : public std::exception
{
  public:
    /**
     * Constructor for the gate_design_exception class.
     *
     * @param ti The tile associated with the error.
     * @param spec The truth table associated with the error.
     * @param portlist The port list associated with the error.
     */
    explicit gate_design_exception(const tile<GateLyt>& ti, const TT& spec,
                                   const port_list<port_direction>& portlist) noexcept :
            error_tile{ti},
            truth_table{spec},
            p{portlist}
    {}
    /**
     * Get the tile associated with the exception.
     */
    [[nodiscard]] tile<GateLyt> which_tile() const noexcept
    {
        return error_tile;
    }
    /**
     * Get the truth table associated with the exception.
     */
    [[nodiscard]] TT which_truth_table() const noexcept
    {
        return truth_table;
    }
    /**
     * Get the port list associated with the exception.
     */
    [[nodiscard]] port_list<port_direction> which_port_list() const noexcept
    {
        return p;
    }

  private:
    /**
     * The tile associated with the error.
     */
    const tile<GateLyt> error_tile{};
    /**
     * The truth table associated with the error.
     */
    const TT truth_table{};
    /**
     * The port list associated with the error.
     */
    const port_list<port_direction> p;
};
/**
 * Exception thrown if the gate design was unsuccessful. Depending on the given gate design parameters and the defect
 * density, the gate design may fail.
 */
class unsuccessful_gate_design_error : public std::runtime_error
{
  public:
    /**
     * This explicit constructor initializes the base `std::runtime_error` class
     * with the provided error message, ensuring that the exception contains
     * detailed information about the reason for the gate design failure.
     *
     * @param msg A descriptive message explaining why the gate design failed.
     */
    explicit unsuccessful_gate_design_error(const std::string_view& msg) noexcept : std::runtime_error(msg.data()) {}
};
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
static void assign_gate(CellLyt& cell_lyt, const typename GateLibrary::fcn_gate& g, const GateLyt& gate_lyt,
                        const tile<GateLyt>& t, const std::optional<typename technology<CellLyt>::cell_type>& only_cell_type = std::nullopt)
{
    const mockturtle::node<GateLyt>& n = gate_lyt.get_node(t);

    // physical design is skipped for input and output ports
    if ((gate_lyt.is_pi(n)/* && bounding_box_2d<GateLyt>{gate_lyt}.get_y_size() >= 4*/) ||
        (gate_lyt.is_po(n)/* && bounding_box_2d<GateLyt>{gate_lyt}.get_y_size() >= 5*/))
    {
        return;
    }

    // retrieve the top-leftmost cell in tile t
    const auto c =
        relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt, CellLyt>(
            gate_lyt, t, cell<CellLyt>{0, 0});

    const auto start_x = c.x;
    const auto start_y = c.y;
    const auto layer   = c.z;

    for (decltype(c.y) y = 0ul; static_cast<uint64_t>(y) < g.size(); ++y)
    {
        for (decltype(c.x) x = 0ul; static_cast<uint64_t>(x) < g[static_cast<uint64_t>(y)].size(); ++x)
        {
            const cell<CellLyt> pos{start_x + x, start_y + y, layer};
            const auto          type{g[static_cast<uint64_t>(y)][static_cast<uint64_t>(x)]};

            if (technology<CellLyt>::is_empty_cell(type) || (only_cell_type.has_value() && *only_cell_type != type))
            {
                continue;
            }

            // overwrites always make a NORMAL type cell
            if (!technology<CellLyt>::is_empty_cell(cell_lyt.get_cell_type(pos)))
            {
                // the cell tile is only overwritten for input cells and output perturber cells
                if (technology<CellLyt>::is_input_cell(cell_lyt.get_cell_type(pos)) ||
                    technology<CellLyt>::is_output_perturber_cell(cell_lyt.get_cell_type(pos)))
                {
                    cell_lyt.assign_cell_tile(pos, typename CellLyt::clock_zone{t.x, t.y, t.z});
                }

                cell_lyt.assign_cell_type(pos, sidb_technology::cell_type::NORMAL);

                continue;
            }

            cell_lyt.assign_cell_tile(pos, typename CellLyt::clock_zone{t.x, t.y, t.z});

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
    if (gate_lyt.is_pi(n))
    {
        return true;//bounding_box_2d<GateLyt>{gate_lyt}.get_y_size() >= 4;
    }

    if (gate_lyt.is_po(n))
    {
        return true;//bounding_box_2d<GateLyt>{gate_lyt}.get_y_size() >= 5;
    }

    return gate_lyt.is_constant(n) || (gate_lyt.is_buf(n) && !gate_lyt.is_ground_layer(gate_lyt.get_tile(n)));
}

template <typename GateLyt>
static bool is_complex_gate(const GateLyt& gate_lyt, const mockturtle::node<GateLyt>& n) noexcept
{
    const auto t  = gate_lyt.get_tile(n);
    const auto at = gate_lyt.above(t);

    return gate_lyt.is_buf(n) && t != at && gate_lyt.is_wire_tile(at);
}

// template <typename Params, typename CellLyt>
// static Params make_gate_design_params_for_complex_gates(const uint64_t canvas_sidb_complex_gates = 0) noexcept
// {
//     Params design_gate_params{};
//     design_gate_params.number_of_sidbs = canvas_sidb_complex_gates;
//     design_gate_params.design_mode     = design_sidb_gates_params<CellLyt>::design_sidb_gates_mode::RANDOM;
//     design_gate_params.canvas          = {{10, 9}, {24, 19}};
//     // design_gate_params.operational_params.op_condition_kinks =
//     // is_operational_params<cell<CellLyt>>::operational_condition_kinks::TOLERATE_KINKS;
//     design_gate_params.operational_params.op_condition_kinks =
//         is_operational_params<cell<CellLyt>>::operational_condition_kinks::REJECT_KINKS;
//     design_gate_params.operational_params.strategy_to_analyze_operational_status =
//         is_operational_params<cell<CellLyt>>::operational_analysis_strategy::FILTER_THEN_SIMULATION;
//
//     return design_gate_params;
// }

}  // namespace fiction

#endif  // GATE_DESIGN_UTILS_HPP

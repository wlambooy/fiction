//
// Created by marcel on 28.06.21.
//

#ifndef FICTION_APPLY_GATE_LIBRARY_HPP
#define FICTION_APPLY_GATE_LIBRARY_HPP

#include "fiction/traits.hpp"
#include "fiction/utils/gate_design_utils.hpp"
#include "fiction/utils/layout_utils.hpp"
#include "fiction/utils/name_utils.hpp"

#include <optional>

#if (PROGRESS_BARS)
#include <mockturtle/utils/progress_bar.hpp>

#include <cstdint>
#endif
#include <mockturtle/traits.hpp>

#include <algorithm>
#include <functional>
#include <set>
#include <type_traits>

// data types cannot properly be converted to bit field types
#pragma GCC diagnostic push
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wuseless-cast"
#endif
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace fiction
{
class sidb_skeleton_bestagon_mini_library;

namespace detail
{

template <typename CellLyt, typename GateLibrary, typename GateLyt>
class apply_gate_library_impl
{
public:
    explicit apply_gate_library_impl(const GateLyt& lyt) :
            gate_lyt{lyt},
            cell_lyt{determine_aspect_ratio_for_cell_level_layout(gate_lyt)}
    {
        cell_lyt.set_tile_size_x(GateLibrary::gate_x_size());
        cell_lyt.set_tile_size_y(GateLibrary::gate_y_size());

        // if GateLyt and CellLyt are based on the same coordinate type, copy the clocking scheme over
        if constexpr (std::is_same_v<coordinate<CellLyt>, coordinate<GateLyt>>)
        {
            cell_lyt.replace_clocking_scheme(gate_lyt.get_clocking_scheme());
        }
        // otherwise, try to find a matching clocking scheme (this will discard overwritten clock numbers)
        else
        {
            if (const auto clk_scheme = get_clocking_scheme<CellLyt>(gate_lyt.get_clocking_scheme().name);
                clk_scheme.has_value())
            {
                cell_lyt.replace_clocking_scheme(clk_scheme.value());
            }
        }
    }

    /**
     * Run the cell layout generation process.
     *
     * This function performs the cell layout generation process based on the gate library and the gate-level layout
     * information provided by `GateLibrary` and `gate_lyt`. It iterates through the nodes in the gate-level layout and
     * maps gates to cell implementations based on their corresponding positions and types. Optionally, it performs
     * post-layout optimization and sets the layout name if certain conditions are met.
     *
     * @param defect_lyt Optional defect surface.
     * @param whitelist A whitelist for the tiles to which gates are assigned. When `std::nullopt` (default), all tiles
     * are considered that are not blacklisted.
     * @param blacklist A blacklist for the tiles to which gates are not assigned. When `std::nullopt` (default), all
     * tiles are considered (unless a whitelist is given). The blacklist has priority over the whitelist.
     * @return A `CellLyt` object representing the generated cell layout.
     */
    [[nodiscard]] CellLyt run_static_gate_library(const std::optional<CellLyt>& defect_surface = std::nullopt, const std::optional<std::set<tile<GateLyt>>>& whitelist,
                                                  const std::optional<std::set<tile<GateLyt>>>& blacklist)
    {
#if (PROGRESS_BARS)
        // initialize a progress bar
        mockturtle::progress_bar bar{static_cast<uint32_t>(gate_lyt.size()), "[i] applying gate library: |{0}|"};
#endif

        gate_lyt.foreach_node(
            [&, this](const auto& n, [[maybe_unused]] auto i)
            {
                if (!gate_lyt.is_constant(n))
                {
                    if (const auto t = gate_lyt.get_tile(n);
                        (!blacklist.has_value() || blacklist.value().count(t) == 0) &&
                        (!whitelist.has_value() || whitelist.value().count(t) != 0))
                    {

                        // retrieve the top-leftmost cell in tile t
                        const auto c =
                            relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(),
                                                               GateLyt, CellLyt>(gate_lyt, t, cell<CellLyt>{0, 0});

                        if constexpr (std::is_same_v<GateLibrary, sidb_skeleton_bestagon_mini_library>)
                        {
                            assign_gate<CellLyt, GateLibrary, GateLyt>(
                                cell_lyt, c, GateLibrary{}.set_up_gate(gate_lyt, t), gate_lyt, n);
                        }
                        else
                        {
                            assign_gate<CellLyt, GateLibrary, GateLyt>(
                                cell_lyt, c, GateLibrary::set_up_gate(gate_lyt, t), gate_lyt, n);
                        }
                    }
                }
#if (PROGRESS_BARS)
                // update progress
                bar(i);
#endif
            });

        // perform post-layout optimization if necessary
        if constexpr (has_post_layout_optimization_v<GateLibrary, CellLyt>)
        {
            GateLibrary::post_layout_optimization(cell_lyt);
        }
        // if available, recover layout name
        if constexpr (has_get_layout_name_v<GateLyt> && has_set_layout_name_v<CellLyt>)
        {
            cell_lyt.set_layout_name(gate_lyt.get_layout_name());
        }

        if constexpr (is_sidb_defect_surface_v<CellLyt>)
        {
            if (defect_surface.has_value())
            {
                // due to issue with windows-2019 Visual Studio 16 2019 and v142. It doesn't compile without using
                // "copy_lyt". When using "cell_lyt.assign_sidb_defect(...)" inside the lambda function, it results in
                // the error: "error C2059: syntax error: '.'".
                auto copy_lyt = cell_lyt.clone();
                // copy the original defects over to the circuit since they are gone when converting the gate-level
                // layout to the cell-level layout.
                defect_surface.value().foreach_sidb_defect([this, &copy_lyt](const auto& def)
                                                           { copy_lyt.assign_sidb_defect(def.first, def.second); });
                return copy_lyt;
            }
        }

        return cell_lyt;
    }
    /**
     * Run the cell layout generation process.
     *
     * This function performs the cell layout generation process based on the SiDB on-the-fly gate library and the
     * gate-level layout information provided by `GateLibrary` and `gate_lyt`. It iterates through the nodes in the
     * gate-level layout and maps gates to cell implementations based on their corresponding positions and types.
     * Optionally, it performs post-layout optimization and sets the layout name if certain conditions are met.
     *
     * @tparam Params Type of the Parameters used for the SiDB on-the-fly gate library.
     * @param params Parameters used for the SiDB on-the-fly gate library.
     * @param defect_surface Optional defect surface.
     * @param whitelist A whitelist for the tiles to which gates are assigned. When `std::nullopt` (default), all tiles
     * are considered that are not blacklisted.
     * @param blacklist A blacklist for the tiles to which gates are not assigned. When `std::nullopt` (default), all
     * tiles are considered (unless a whitelist is given). The blacklist has priority over the whitelist.
     * @return A `CellLyt` object representing the generated cell layout.
     */
    template <typename Params>
    [[nodiscard]] auto run_parameterized_gate_library(const Params&                                       params,
                                                      const std::optional<CellLyt>& defect_surface = std::nullopt,
                                                         const std::optional<std::set<tile<GateLyt>>>& whitelist,
                                                         const std::optional<std::set<tile<GateLyt>>>& blacklist)
    {
#if (PROGRESS_BARS)
        // initialize a progress bar
        mockturtle::progress_bar bar{static_cast<uint32_t>(gate_lyt.size()), "[i] applying gate library: |{0}|"};
#endif
        // perform post-layout optimization if necessary
        if constexpr (has_post_layout_optimization_v<GateLibrary, CellLyt>)
        {
            GateLibrary::post_layout_optimization(gate_lyt);
        }

        gate_lyt.foreach_node(
            [&, this](const auto& n, [[maybe_unused]] auto i)
            {
                if (!gate_lyt.is_constant(n))
                {

                    if (const auto t = gate_lyt.get_tile(n);
                        (!blacklist.has_value() || blacklist.value().count(t) == 0) &&
                        (!whitelist.has_value() || whitelist.value().count(t) != 0))
                    {
                        // retrieve the top-leftmost cell in tile t
                        const auto c =
                            relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(),
                                                               GateLyt, CellLyt>(gate_lyt, t, cell<CellLyt>{0, 0});

                        const auto gate =
                            GateLibrary::template set_up_gate<GateLyt, CellLyt, Params>(gate_lyt, t, params, defect_surface);

                        assign_gate<CellLyt, GateLibrary, GateLyt>(cell_lyt, c, gate, gate_lyt, n);
                    }
                }
#if (PROGRESS_BARS)
                // update progress
                bar(i);
#endif
            });

        // if available, recover layout name
        cell_lyt.set_layout_name(get_name(gate_lyt));

        if constexpr (is_sidb_defect_surface_v<CellLyt>)
        {
            if (defect_surface.has_value())
            {
                // due to issue with windows-2019 Visual Studio 16 2019 and v142. It doesn't compile without using
                // "copy_lyt". When using "cell_lyt.assign_sidb_defect(...)" inside the lambda function, it results in
                // the error: "error C2059: syntax error: '.'".
                auto copy_lyt = cell_lyt.clone();
                // copy the original defects over to the circuit since they are gone when converting the gate-level
                // layout to the cell-level layout.
                defect_surface.value().foreach_sidb_defect([this, &copy_lyt](const auto& def)
                                                           { copy_lyt.assign_sidb_defect(def.first, def.second); });
                return copy_lyt;
            }
        }

        return cell_lyt;
    }

private:
    /**
     * Gate-level layout.
     */
    GateLyt gate_lyt;
    /**
     * Cell-level layout.
     */
    CellLyt cell_lyt;

    /**
     * Computes the (inclusively) bounding coordinate for a cell-level layout that is derived from the dimensions of the
     * given gate-level layout, while respecting tiling geometry in which even and odd rows/columns do not line up.
     *
     * @param gate_lyt Gate-level layout of which the dimensions are read.
     * @return Aspect ratio for a cell-level layout that corresponds to the dimensions of the given gate-level layout.
     */
    static aspect_ratio<CellLyt> determine_aspect_ratio_for_cell_level_layout(const GateLyt& gate_lyt) noexcept
    {
        const std::function<cell<CellLyt>(GateLyt, tile<GateLyt>, cell<CellLyt>)> rel_to_abs_cell_pos =
            relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                               CellLyt>;

        const cell<CellLyt> max_rel_coord = {GateLibrary::gate_x_size() - 1, GateLibrary::gate_y_size() - 1};

        const cell<CellLyt> first_odd_tile = {gate_lyt.x() != 0 ? 1 : 0, gate_lyt.y() != 0 ? 1 : 0};

        const auto max_coord_even_x = rel_to_abs_cell_pos(gate_lyt, {0, gate_lyt.y()}, max_rel_coord);
        const auto max_coord_odd_x  = rel_to_abs_cell_pos(gate_lyt, {first_odd_tile.x, gate_lyt.y()}, max_rel_coord);
        const auto max_coord_even_y = rel_to_abs_cell_pos(gate_lyt, {gate_lyt.x(), 0}, max_rel_coord);
        const auto max_coord_odd_y  = rel_to_abs_cell_pos(gate_lyt, {gate_lyt.x(), first_odd_tile.y}, max_rel_coord);

        return {std::max(max_coord_even_y.x, max_coord_odd_y.x), std::max(max_coord_even_x.y, max_coord_odd_x.y)};
    }
};

}  // namespace detail

/**
 * Applies a gate library to a given gate-level layout and, thereby, creates and returns a cell-level layout. The gate
 * library type should provide all functions specified in fcn_gate_library. It is, thus, easiest to extend
 * fcn_gate_library to implement a new gate library. Examples are `qca_one_library`, `inml_topolinano_library`, and
 * `sidb_bestagon_library`.
 *
 * May pass through, and thereby throw, an `unsupported_gate_type_exception` or an
 * `unsupported_gate_orientation_exception`.
 *
 * @tparam CellLyt Type of the returned cell-level layout.
 * @tparam GateLibrary Type of the gate library to apply.
 * @tparam GateLyt Type of the gate-level layout to apply the library to.
 * @param lyt The gate-level layout.
 * @param whitelist A whitelist for the tiles to which gates are assigned. When `std::nullopt` (default), all tiles are
 * considered that are not blacklisted.
 * @param blacklist A blacklist for the tiles to which gates are not assigned. When `std::nullopt` (default), all tiles
 * are considered (unless a whitelist is given). The blacklist has priority over the whitelist.
 * @return A cell-level layout that implements `lyt`'s gate types with building blocks defined in `GateLibrary`.
 */
template <typename CellLyt, typename GateLibrary, typename GateLyt>
[[nodiscard]] CellLyt apply_gate_library(const GateLyt&                                lyt,
                                         const std::optional<std::set<tile<GateLyt>>>& whitelist = {},
                                         const std::optional<std::set<tile<GateLyt>>>& blacklist = {})
{
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(!has_siqad_coord_v<CellLyt>, "CellLyt cannot have SiQAD coordinates");
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(mockturtle::has_is_constant_v<GateLyt>, "GateLyt does not implement the is_constant function");
    static_assert(mockturtle::has_foreach_node_v<GateLyt>, "GateLyt does not implement the foreach_node function");

    static_assert(std::is_same_v<technology<CellLyt>, technology<GateLibrary>>,
                  "CellLyt and GateLibrary must implement the same technology");

    detail::apply_gate_library_impl<CellLyt, GateLibrary, GateLyt> p{lyt};

    return p.run_static_gate_library({}, whitelist, blacklist);
}

/**
 * Applies a gate library to a given gate-level layout and maps the SiDB and defect locations onto a defect surface. The
 * gate library type should provide all functions specified in fcn_gate_library. It is, thus, easiest to extend
 * fcn_gate_library to implement a new gate library. Examples are `qca_one_library`, `inml_topolinano_library`, and
 * `sidb_bestagon_library`.
 *
 * May pass through, and thereby throw, an `unsupported_gate_type_exception` or an
 * `unsupported_gate_orientation_exception`.
 *
 * @tparam DefectLyt Type of the returned cell-level layout.
 * @tparam GateLibrary Type of the gate library to apply.
 * @tparam GateLyt Type of the gate-level layout to apply the library to.
 * @param lyt The gate-level layout.
 * @return A cell-level layout that implements `lyt`'s gate types with building blocks defined in `GateLibrary`.
 */
template <typename DefectLyt, typename GateLibrary, typename GateLyt>
[[nodiscard]] DefectLyt apply_gate_library_to_defective_surface(const GateLyt& lyt, const DefectLyt& defect_surface)
{
    static_assert(is_cell_level_layout_v<DefectLyt>, "DefectLyt is not a cell-level layout");
    static_assert(is_sidb_defect_surface_v<DefectLyt>, "DefectLyt is not an SiDB defect surface");
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(mockturtle::has_is_constant_v<GateLyt>, "GateLyt does not implement the is_constant function");
    static_assert(mockturtle::has_foreach_node_v<GateLyt>, "GateLyt does not implement the foreach_node function");

    static_assert(std::is_same_v<technology<DefectLyt>, technology<GateLibrary>>,
                  "DefectLyt and GateLibrary must implement the same technology");

    detail::apply_gate_library_impl<DefectLyt, GateLibrary, GateLyt> p{lyt};

    return p.run_static_gate_library(defect_surface, {}, {});
}
/**
 * Applies a parameterized gate library to a given
 * gate-level layout and, thereby, creates and returns a cell-level layout.
 *
 * May pass through, and thereby throw, an `unsupported_gate_type_exception`, an
 * `unsupported_gate_orientation_exception` and any further custom exceptions of the gate libraries.
 *
 * @tparam CellLyt Type of the returned cell-level layout.
 * @tparam GateLibrary Type of the gate library to apply.
 * @tparam GateLyt Type of the gate-level layout to apply the library to.
 * @tparam Params Type of the parameter used for SiDB on-the-fly gate library.
 * @param lyt The gate-level layout.
 * @param params Parameter for the gate library.
 * @param whitelist A whitelist for the tiles to which gates are assigned. When `std::nullopt` (default), all tiles are
 * considered that are not blacklisted.
 * @param blacklist A blacklist for the tiles to which gates are not assigned. When `std::nullopt` (default), all tiles
 * are considered (unless a whitelist is given). The blacklist has priority over the whitelist.
 * @return A cell-level layout that implements `lyt`'s gate types with building blocks defined in `GateLibrary`.
 */
template <typename CellLyt, typename GateLibrary, typename GateLyt, typename Params>
[[nodiscard]] CellLyt apply_parameterized_gate_library(const GateLyt& lyt, Params& params,
                                                       const std::optional<std::set<tile<GateLyt>>>& whitelist = {},
                                                       const std::optional<std::set<tile<GateLyt>>>& blacklist = {})
{
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(has_cube_coord_v<CellLyt>, "CellLyt must be based on cube coordinates");
    static_assert(mockturtle::has_is_constant_v<GateLyt>, "GateLyt does not implement the is_constant function");
    static_assert(mockturtle::has_foreach_node_v<GateLyt>, "GateLyt does not implement the foreach_node function");

    static_assert(std::is_same_v<technology<CellLyt>, technology<GateLibrary>>,
                  "CellLyt and GateLibrary must implement the same technology");

    detail::apply_gate_library_impl<CellLyt, GateLibrary, GateLyt> p{lyt};

    return p.template run_parameterized_gate_library<Params>(params, whitelist, blacklist);
}


/**
 * Applies a defect-aware parameterized gate library to a given
 * gate-level layout and, thereby, creates and returns a cell-level layout.
 *
 * May pass through, and thereby throw, an `unsupported_gate_type_exception`, an
 * `unsupported_gate_orientation_exception` and any further custom exceptions of the gate libraries.
 *
 * @tparam DefectLyt Type of the returned cell-level layout.
 * @tparam GateLibrary Type of the gate library to apply.
 * @tparam GateLyt Type of the gate-level layout to apply the library to.
 * @tparam Params Type of the parameter used for SiDB on-the-fly gate library.
 * @param lyt The gate-level layout.
 * @param params Parameter for the gate library.
 * @param defect_surface Defect surface.
 * @return A cell-level layout that implements `lyt`'s gate types with building blocks defined in `GateLibrary`.
 */
template <typename DefectLyt, typename GateLibrary, typename GateLyt, typename Params>
[[nodiscard]] DefectLyt apply_parameterized_gate_library_to_defective_surface(const GateLyt& lyt, const Params& params,
                                                                              const DefectLyt& defect_surface)
{
    static_assert(is_cell_level_layout_v<DefectLyt>, "DefectLyt is not a cell-level layout");
    static_assert(is_sidb_defect_surface_v<DefectLyt>, "DefectLyt is not an SiDB defect surface");
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(has_cube_coord_v<DefectLyt>, "DefectLyt must be based on cube coordinates");
    static_assert(mockturtle::has_is_constant_v<GateLyt>, "GateLyt does not implement the is_constant function");
    static_assert(mockturtle::has_foreach_node_v<GateLyt>, "GateLyt does not implement the foreach_node function");

    static_assert(std::is_same_v<technology<DefectLyt>, technology<GateLibrary>>,
                  "DefectLyt and GateLibrary must implement the same technology");

    detail::apply_gate_library_impl<DefectLyt, GateLibrary, GateLyt> p{lyt};

    // Run the gate library with the parameters
    const DefectLyt result = p.template run_parameterized_gate_library<Params>(params, defect_surface, {}, {});

    return result;
}

}  // namespace fiction

#pragma GCC diagnostic pop

#endif  // FICTION_APPLY_GATE_LIBRARY_HPP

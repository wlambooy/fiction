//
// Created by Willem Lambooy on 06/06/2025.
//

#ifndef SIDB_BDL_CIRCUIT_HPP
#define SIDB_BDL_CIRCUIT_HPP

#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/traits.hpp"
#include "kitty/print.hpp"

#include <array>
#include <set>
#include <utility>
#include <vector>

#include <signal.h>

namespace fiction
{

template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
class sidb_bdl_circuit
{
  public:
    explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const detect_bdl_wires_params& bdl_wire_params,
                              const std::optional<tile<GateLyt>>& this_tile = std::nullopt) noexcept :
            gate_layout{gate_lyt.clone()},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            gate_connections{get_gate_connections(bdl_wires, skeleton, gate_lyt)},
            gate_tile{this_tile}
    {}
    /**
     * SiDB gate-level layout.
     */
    const GateLyt gate_layout;

    const CellLyt skeleton;

    const std::vector<bdl_wire<CellLyt>>                       bdl_wires{};
    const uint64_t                                             num_bdl_pairs{};
    const std::vector<bdl_pair<cell<CellLyt>>>                 input_bdl_pairs{};
    const uint64_t                                             num_inputs{};
    const std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};

    const std::optional<tile<GateLyt>> gate_tile{};

    static uint64_t get_number_of_bdl_pairs(const std::vector<bdl_wire<CellLyt>>& bdl_wires) noexcept
    {
        uint64_t total = 0;

        for (const bdl_wire<CellLyt>& wire : bdl_wires)
        {
            total += wire.pairs.size();
        }

        return total;
    }

  private:
    /**
     *
     */
    static std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>>
    get_gate_connections(const std::vector<bdl_wire<CellLyt>>& bdl_wires, const CellLyt& lyt,
                         const GateLyt& gate_lyt) noexcept
    {
        std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};
        gate_connections.reserve(bdl_wires.size());

        for (const bdl_wire<CellLyt>& wire : bdl_wires)
        {
            assert(!wire.pairs.empty() && "BDL wire is empty");

            const bool right_to_left = wire.pairs.front().upper.x > wire.pairs.front().lower.x;
            const bool is_input_wire = wire.pairs.front().type == sidb_technology::cell_type::INPUT;

            std::set<tile<GateLyt>> tiles_in_wire = {};

            for (const bdl_pair<cell<CellLyt>>& pair : wire.pairs)
            {
                // std::cout << "\ncell: " << pair.upper.x << " " << pair.upper.y << " at tile "
                //           << lyt.get_cell_tile(pair.upper).x << ',' << lyt.get_cell_tile(pair.upper).y << ','
                //           << lyt.get_cell_tile(pair.upper).z << std::endl;
                // std::cout << "cell: " << pair.lower.x << " " << pair.lower.y << " at tile "
                //           << lyt.get_cell_tile(pair.lower).x << ',' << lyt.get_cell_tile(pair.lower).y << ','
                //           << lyt.get_cell_tile(pair.lower).z << std::endl;
                assert(lyt.get_cell_tile(pair.upper) == lyt.get_cell_tile(pair.lower));

                tiles_in_wire.insert(
                    [&](const auto& this_t)
                    {
                        const tile<GateLyt> t = {this_t.x, this_t.y};
                        auto                z = 0;

                        // std::cout << "current tile: " << t.x << ", " << t.y << std::endl;

                        if (const auto at = gate_lyt.above(t); at != t && gate_lyt.is_wire_tile(at))
                        {
                            assert(gate_lyt.is_wire_tile(t) && "stacked tiles are not both wires");

                            const bool is_input_wire_to_this_tile = is_input_wire || !tiles_in_wire.empty();

                            // stacked tiles are both populated with wire tiles  ---  CROSSING or DOUBLE WIRE

                            // todo: does NOT work for crossing / double wire in connection .. ?

                            if (gate_lyt.below(gate_lyt.north_west(at)) ==
                                gate_lyt.below(gate_lyt.incoming_data_flow(at).front()))
                            {
                                // upper tile connects to north west

                                if (gate_lyt.below(gate_lyt.south_east(at)) ==
                                    gate_lyt.below(gate_lyt.outgoing_data_flow(at).front()))
                                {
                                    // upper tile connects to south east  ---  CROSSING CASE (upper tile goes L->R)
                                    if (!right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                                else
                                {
                                    // lower tile connects to south east  ---  DOUBLE WIRE CASE
                                    if (is_input_wire_to_this_tile != right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                            }
                            else
                            {
                                // lower tile connects to north west

                                if (gate_lyt.below(gate_lyt.south_east(at)) ==
                                    gate_lyt.below(gate_lyt.outgoing_data_flow(at).front()))
                                {
                                    // upper tile connects to south east  ---  DOUBLE WIRE CASE
                                    if (is_input_wire_to_this_tile == right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                                else
                                {
                                    // lower tile connects to south east  ---  CROSSING CASE (upper tile goes R->L)
                                    if (right_to_left)
                                    {
                                        z = 1;
                                    }
                                }
                            }
                        }
                        else if (!is_input_wire && !tiles_in_wire.empty())
                        {
                            const std::vector<tile<GateLyt>>& outgoing_tiles =
                                gate_lyt.outgoing_data_flow(*tiles_in_wire.cbegin());

                            if (outgoing_tiles.size() == 1)
                            {
                                z = outgoing_tiles.front().z;
                            }
                            else
                            {
                                assert(outgoing_tiles.size() == 2 && "fan out must be either 1 or 2");

                                if (right_to_left)
                                {
                                    z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.front().z :
                                                                                             outgoing_tiles.back().z;
                                }
                                else
                                {
                                    z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.back().z :
                                                                                             outgoing_tiles.front().z;
                                }
                            }
                        }

                        // std::cout << "z: " << z << std::endl;

                        return tile<GateLyt>{t.x, t.y, z};
                    }(lyt.get_cell_tile(pair.upper)));
            }

            assert((tiles_in_wire.size() == 2 ||
                    (tiles_in_wire.size() == 1 &&
                     (is_input_wire || wire.pairs.back().type == sidb_technology::cell_type::OUTPUT))) &&
                   "Gate connection is malformed");

            const auto get_tile_pair = [&]
            {
                // std::cout << "\ntile in wire: " << tiles_in_wire.cbegin()->x << ',' << tiles_in_wire.cbegin()->y <<
                // ','
                // << tiles_in_wire.cbegin()->z << std::endl;
                if (tiles_in_wire.size() == 2)
                {
                    // std::cout << "tile TO wire: " << tiles_in_wire.crbegin()->x << ',' << tiles_in_wire.crbegin()->y
                    // << ',' << tiles_in_wire.crbegin()->z << std::endl;
                    const tile<GateLyt>& first_tile  = *tiles_in_wire.cbegin();
                    const tile<GateLyt>& second_tile = *tiles_in_wire.crbegin();

                    if (first_tile.y < second_tile.y || (first_tile.y == second_tile.y && first_tile.x < second_tile.x))
                    {
                        return std::make_pair(first_tile, second_tile);
                    }

                    return std::make_pair(second_tile, first_tile);
                }

                assert(tiles_in_wire.size() == 1 &&
                       (is_input_wire || wire.pairs.back().type == sidb_technology::cell_type::OUTPUT) &&
                       "Gate connection is malformed");

                const tile<GateLyt>& tile_in_wire = *tiles_in_wire.cbegin();

                // todo: does NOT work for crossing / double wire in connection

                if (wire.pairs.front().type == sidb_technology::cell_type::INPUT)
                {

                    const tile<GateLyt>& upper_tile = right_to_left ?
                                                          gate_lyt.below(gate_lyt.north_east(tile_in_wire)) :
                                                          gate_lyt.below(gate_lyt.north_west(tile_in_wire));

                    auto z = 0;

                    const std::vector<tile<GateLyt>>& incoming_tiles = gate_lyt.incoming_data_flow(tile_in_wire);

                    if (incoming_tiles.size() == 1)
                    {
                        z = incoming_tiles.front().z;
                    }
                    else if (incoming_tiles.size() == 2)
                    {
                        // assert(incoming_tiles.size() == 2 && "fan in must be either 1 or 2");

                        if (right_to_left)
                        {
                            z = incoming_tiles.front().x < incoming_tiles.back().x ? incoming_tiles.back().z :
                                                                                     incoming_tiles.front().z;
                        }
                        else
                        {
                            z = incoming_tiles.front().x < incoming_tiles.back().x ? incoming_tiles.front().z :
                                                                                     incoming_tiles.back().z;
                        }
                    }

                    // std::cout << "tile INP wire: " << upper_tile.x << ',' << upper_tile.y << ',' << z << std::endl;
                    return std::make_pair(tile<GateLyt>{upper_tile.x, upper_tile.y, z}, tile_in_wire);
                }

                const tile<GateLyt>& lower_tile = right_to_left ? gate_lyt.below(gate_lyt.south_west(tile_in_wire)) :
                                                                  gate_lyt.below(gate_lyt.south_east(tile_in_wire));

                auto z = 0;

                const std::vector<tile<GateLyt>>& outgoing_tiles = gate_lyt.outgoing_data_flow(tile_in_wire);

                if (outgoing_tiles.size() == 1)
                {
                    z = outgoing_tiles.front().z;
                }
                else if (outgoing_tiles.size() == 2)
                {
                    // assert(outgoing_tiles.size() == 2 && "fan out must be either 1 or 2");

                    if (right_to_left)
                    {
                        z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.front().z :
                                                                                 outgoing_tiles.back().z;
                    }
                    else
                    {
                        z = outgoing_tiles.front().x < outgoing_tiles.back().x ? outgoing_tiles.back().z :
                                                                                 outgoing_tiles.front().z;
                    }
                }
                // std::cout << "tile OUT wire: " << lower_tile.x << ',' << lower_tile.y << ',' << z << std::endl;
                return std::make_pair(tile_in_wire, tile<GateLyt>{lower_tile.x, lower_tile.y, z});
            };

            std::pair<tile<GateLyt>, tile<GateLyt>> tile_pair_at_gate_connection = get_tile_pair();
            // std::cout << fmt::format("tile pair at gate connection: {},{},{}   {},{},{}",
            //                          tile_pair_at_gate_connection.first.x, tile_pair_at_gate_connection.first.y,
            //                          tile_pair_at_gate_connection.first.z, tile_pair_at_gate_connection.second.x,
            //                          tile_pair_at_gate_connection.second.y, tile_pair_at_gate_connection.second.z)
            //           << std::endl;
            // assert(tile_pair_at_gate_connection.first.y == tile_pair_at_gate_connection.second.y - 1 &&
            //        "tiles are not represent a row clocked gate connection"); todo

            gate_connections.push_back(std::move(tile_pair_at_gate_connection));
        }

        return gate_connections;
    }
};

template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
class sidb_bdl_sub_circuit
{
  public:
    explicit sidb_bdl_sub_circuit(
        const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& bdl_super_circuit) noexcept :
            super_circuit{bdl_super_circuit},
            gate_layout{super_circuit.gate_layout},
            skeleton{super_circuit.skeleton},
            bdl_wires{super_circuit.bdl_wires},
            num_bdl_pairs{super_circuit.num_bdl_pairs},
            input_bdl_pairs{super_circuit.input_bdl_pairs},
            num_inputs{super_circuit.num_inputs},
            gate_connections{filter_gate_connections_from_super_circuit(super_circuit, bdl_wires)}
    {}

    explicit sidb_bdl_sub_circuit(const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& bdl_super_circuit,
                                  const std::vector<tile<GateLyt>>&                              tiles,
                                  const detect_bdl_wires_params&                                 bdl_wire_params,
                                  const std::optional<tile<GateLyt>>& this_tile = std::nullopt) noexcept :
            super_circuit{bdl_super_circuit},
            gate_layout{create_gate_lyt_window_for_connection_sequence(super_circuit.gate_layout, tiles)},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            gate_connections{filter_gate_connections_from_super_circuit(super_circuit, bdl_wires)},
            gate_tile{this_tile}
    {}

    sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary> super_circuit{};

    /**
     * SiDB gate-level layout.
     */
    const GateLyt gate_layout;

    const CellLyt skeleton;

    const std::vector<bdl_wire<CellLyt>>                       bdl_wires{};
    const uint64_t                                             num_bdl_pairs{};
    const std::vector<bdl_pair<cell<CellLyt>>>                 input_bdl_pairs{};
    const uint64_t                                             num_inputs{};
    const std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};

    const std::optional<tile<GateLyt>> gate_tile{};

  private:
    static void augment_gate_lyt_window(const GateLyt& gate_lyt, GateLyt& gate_lyt_window,
                                        const tile<GateLyt>&           current_t,
                                        const std::set<tile<GateLyt>>& connecting_inputs  = {},
                                        const std::set<tile<GateLyt>>& connecting_outputs = {}) noexcept
    {
        std::vector<mockturtle::signal<GateLyt>> inputs_to_current_t{};

        for (auto in_t : gate_lyt.incoming_data_flow(current_t))
        {
            if (connecting_inputs.count(gate_lyt.below(in_t)) == 0)
            {
                if (!gate_lyt_window.is_empty_tile(in_t))
                {
                    in_t.z = 1 - in_t.z;
                }

                gate_lyt_window.create_pi("", in_t);
            }

            inputs_to_current_t.push_back(static_cast<mockturtle::signal<GateLyt>>(in_t));
        }

        assert(gate_lyt_window.is_empty_tile(current_t) && "tile on which node is to be created is already populated");
        gate_lyt_window.create_node(inputs_to_current_t, gate_lyt.node_function(gate_lyt.get_node(current_t)),
                                    current_t);

        for (auto out_t : gate_lyt.outgoing_data_flow(current_t))
        {
            if (connecting_outputs.count(gate_lyt.below(out_t)) == 0)
            {
                if (!gate_lyt_window.is_empty_tile(out_t))
                {
                    out_t.z = 1 - out_t.z;
                }

                gate_lyt_window.create_po(static_cast<mockturtle::signal<GateLyt>>(current_t), "", out_t);
            }
        }

        if (gate_lyt.is_wire_tile(current_t))
        {
            if (const auto above_current_t = gate_lyt.above(current_t);
                current_t != above_current_t && gate_lyt.is_wire_tile(above_current_t))
            {
                // upper_t hosts a crossing or double wire

                augment_gate_lyt_window(gate_lyt, gate_lyt_window, above_current_t, connecting_inputs,
                                        connecting_outputs);
            }
        }
    }

    [[nodiscard]] static GateLyt
    create_gate_lyt_window_for_connection_sequence(const GateLyt& gate_lyt, std::vector<tile<GateLyt>> tiles) noexcept
    {
        GateLyt gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};

        // Sort tiles in ascending y, then x, then z for deterministic ordering
        std::sort(tiles.begin(), tiles.end());

        const auto same_clock_zone = [](const tile<GateLyt>& a, const tile<GateLyt>& b)
        {
            return a.y == b.y;  // Adjust if clock zone definition differs from y==y
        };

        // Process each tile in sorted order
        for (size_t i = 0; i < tiles.size(); ++i)
        {
            std::set<tile<GateLyt>> inputs, outputs;

            // --- Determine inputs ---
            for (size_t j = 0; j < i; ++j)
            {
                if (tiles[i].y - tiles[j].y > 1)
                {
                    // tiles are not directly connected
                    continue;
                }

                if (same_clock_zone(tiles[j], tiles[i]))
                {
                    // tiles are in the same row and thus not connected
                    continue;
                }

                inputs.insert(tiles[j]);
            }

            // --- Determine outputs ---
            for (size_t j = i + 1; j < tiles.size(); ++j)
            {
                if (tiles[j].y - tiles[i].y > 1)
                {
                    // tiles are not directly connected
                    continue;
                }

                if (same_clock_zone(tiles[j], tiles[i]))
                {
                    // tiles are in the same row and thus not connected
                    continue;
                }

                outputs.insert(tiles[j]);
            }

            augment_gate_lyt_window(gate_lyt, gate_lyt_window, tiles[i], inputs, outputs);
        }

        return gate_lyt_window;
    }

    static std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> filter_gate_connections_from_super_circuit(
        const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& super_circuit,
        const std::vector<bdl_wire<CellLyt>>&                          sub_circuit_bdl_wires) noexcept
    {
        assert(!sub_circuit_bdl_wires.empty() && "There are no BDL wires in the sub-circuit.");

        std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> sub_circuit_gate_connections{};
        sub_circuit_gate_connections.reserve(sub_circuit_bdl_wires.size());

        for (const bdl_wire<CellLyt>& sub_circuit_bdl_wire : sub_circuit_bdl_wires)
        {
            uint64_t super_circuit_bdl_wires_ix = 0;

            for (; super_circuit_bdl_wires_ix < super_circuit.bdl_wires.size(); ++super_circuit_bdl_wires_ix)
            {
                if (sub_circuit_bdl_wire.pairs.front().upper ==
                        super_circuit.bdl_wires.at(super_circuit_bdl_wires_ix).pairs.front().upper ||
                    sub_circuit_bdl_wire.pairs.back().upper ==
                        super_circuit.bdl_wires.at(super_circuit_bdl_wires_ix).pairs.back().upper)
                {
                    sub_circuit_gate_connections.push_back(
                        super_circuit.gate_connections.at(super_circuit_bdl_wires_ix));

                    break;
                }
            }

            assert(super_circuit_bdl_wires_ix < super_circuit.bdl_wires.size() &&
                   "Sub-circuit BDL wire could not be matched to one of the super circuit.");
        }

        return sub_circuit_gate_connections;
    }
};

template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
struct sidb_cell_level_bdl_circuit
{
    /**
     * SiDB cell-level layout.
     */
    const CellLyt&                                                     cell_layout;
    const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& circuit{};

    // todo construct using std::unordered_map<tile<GateLyt>, canvas_positions> ?
    explicit sidb_cell_level_bdl_circuit(
        const CellLyt& lyt, const sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& bdl_circuit) noexcept :
            cell_layout{lyt},
            circuit{bdl_circuit}
    {}
};

}  // namespace fiction

#endif  // SIDB_BDL_CIRCUIT_HPP

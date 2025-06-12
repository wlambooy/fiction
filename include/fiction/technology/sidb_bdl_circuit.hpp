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

    explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const tile<GateLyt>& t, const tile<GateLyt>& connecting_t,
                              const detect_bdl_wires_params& bdl_wire_params) noexcept :
            gate_layout{create_gate_lyt_window_for_gate_connection(gate_lyt, t, connecting_t)},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            gate_connections{get_gate_connections(bdl_wires, skeleton, gate_layout)}
    {}

    explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const tile<GateLyt>& t, const tile<GateLyt>& connecting_t,
                              const tile<GateLyt>&           connecting_to_connecting_t,
                              const detect_bdl_wires_params& bdl_wire_params) noexcept :
            gate_layout{
                create_gate_lyt_window_for_two_gate_connections(gate_lyt, t, connecting_t, connecting_to_connecting_t)},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            gate_connections{get_gate_connections(bdl_wires, skeleton, gate_layout)}
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

                // std::cout << "created pi at " << in_t.x << " " << in_t.y << " " << in_t.z << std::endl;
            }

            inputs_to_current_t.push_back(static_cast<mockturtle::signal<GateLyt>>(in_t));
        }

        // std::cout << "created node at " << current_t.x << " " << current_t.y << " " << current_t.z << " with inputs: ";
        // for (const auto& s : inputs_to_current_t)
        // {
        //     const auto& t = static_cast<tile<GateLyt>>(s);
        //     std::cout << t.x << " " << t.y << " " << t.z << " \t ";
        // }
        // std::cout << std::endl;
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

                // std::cout << "created po at " << out_t.x << " " << out_t.y << " " << out_t.z << " with input "
                //           << current_t.x << " " << current_t.y << " " << current_t.z << std::endl;

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
    };

    [[nodiscard]] static GateLyt create_gate_lyt_window_for_gate_connection(const GateLyt&       gate_lyt,
                                                                            const tile<GateLyt>& t,
                                                                            const tile<GateLyt>& connecting_t) noexcept
    {
        GateLyt gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};
        // std::cout << "start num pis" << gate_lyt_window.num_pis() << std::endl;

        const tile<GateLyt>& upper_t = t.y < connecting_t.y ? t : connecting_t;
        const tile<GateLyt>& lower_t = t.y < connecting_t.y ? connecting_t : t;

        augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {lower_t});
        augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {upper_t}, {});

        return gate_lyt_window;
    }

    [[nodiscard]] static GateLyt
    create_gate_lyt_window_for_two_gate_connections(const GateLyt& gate_lyt, const tile<GateLyt>& t,
                                                    const tile<GateLyt>& connecting_t,
                                                    const tile<GateLyt>& connecting_to_connecting_t) noexcept
    {
        GateLyt gate_lyt_window{{gate_lyt.x(), gate_lyt.y(), gate_lyt.z()}, row_clocking<GateLyt>()};
        // std::cout << "start num pis" << gate_lyt_window.num_pis() << std::endl;

        std::array<tile<GateLyt>, 3> tiles_sorted{{t, connecting_t, connecting_to_connecting_t}};
        std::sort(tiles_sorted.begin(), tiles_sorted.end());

        const tile<GateLyt>& upper_t  = tiles_sorted[0];
        const tile<GateLyt>& middle_t = tiles_sorted[1];
        const tile<GateLyt>& lower_t  = tiles_sorted[2];

        assert(upper_t.y <= middle_t.y && middle_t.y <= lower_t.y && "tiles are not sorted by row number");

        if (middle_t.y == lower_t.y || middle_t.y == upper_t.y)
        {
            // a clock zone is shared

            assert(upper_t.y != lower_t.y && "all tiles are in the same clockzone");

            if (middle_t.y < lower_t.y)
            {
                // case: lower_t follows the others and gets input from both
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {lower_t});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, middle_t, {}, {lower_t});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {upper_t, middle_t}, {});
            }
            else
            {
                // case: upper_t precedes the others and gives output to both
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {middle_t, lower_t});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, middle_t, {upper_t}, {});
                augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {upper_t}, {});
            }

            return gate_lyt_window;
        }

        augment_gate_lyt_window(gate_lyt, gate_lyt_window, upper_t, {}, {middle_t});
        augment_gate_lyt_window(gate_lyt, gate_lyt_window, middle_t, {upper_t}, {lower_t});
        augment_gate_lyt_window(gate_lyt, gate_lyt_window, lower_t, {middle_t}, {});

        return gate_lyt_window;
    }
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
                    else
                    {
                        assert(incoming_tiles.size() == 2 && "fan in must be either 1 or 2");

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
                // std::cout << "tile OUT wire: " << lower_tile.x << ',' << lower_tile.y << ',' << z << std::endl;
                return std::make_pair(tile_in_wire, tile<GateLyt>{lower_tile.x, lower_tile.y, z});
            };

            std::pair<tile<GateLyt>, tile<GateLyt>> tile_pair_at_gate_connection = get_tile_pair();
            // std::cout << fmt::format("tile pair at gate connection: {},{},{}   {},{},{}",
            //                          tile_pair_at_gate_connection.first.x, tile_pair_at_gate_connection.first.y,
            //                          tile_pair_at_gate_connection.first.z, tile_pair_at_gate_connection.second.x,
            //                          tile_pair_at_gate_connection.second.y, tile_pair_at_gate_connection.second.z)
            //           << std::endl;
            assert(tile_pair_at_gate_connection.first.y == tile_pair_at_gate_connection.second.y - 1 &&
                   "tiles are not represent a row clocked gate connection");

            gate_connections.push_back(std::move(tile_pair_at_gate_connection));
        }

        return gate_connections;
    }

    static uint64_t get_number_of_bdl_pairs(const std::vector<bdl_wire<CellLyt>>& bdl_wires) noexcept
    {
        uint64_t total = 0;

        for (const bdl_wire<CellLyt>& wire : bdl_wires)
        {
            total += wire.pairs.size();
        }

        return total;
    }
};
template <typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
struct sidb_cell_level_bdl_circuit
{
    /**
     * SiDB cell-level layout.
     */
    const CellLyt&                                          cell_layout;
    sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary> circuit{};

    // todo construct using std::unordered_map<tile<GateLyt>, canvas_positions>
    explicit sidb_cell_level_bdl_circuit(
        const CellLyt& lyt, const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& bdl_circuit) noexcept :
            cell_layout{lyt},
            circuit{bdl_circuit}
    {}
};

}  // namespace fiction

#endif  // SIDB_BDL_CIRCUIT_HPP

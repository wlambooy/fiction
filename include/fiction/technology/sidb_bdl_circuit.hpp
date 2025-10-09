//
// Created by Willem Lambooy on 06/06/2025.
//

#ifndef SIDB_BDL_CIRCUIT_HPP
#define SIDB_BDL_CIRCUIT_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/traits.hpp"

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
    explicit sidb_bdl_circuit(const GateLyt& gate_lyt, const sidb_simulation_parameters& simulation_parameters,
                              const std::pair<cell<CellLyt>, cell<CellLyt>>& rel_canvas,
                              const detect_bdl_wires_params& bdl_wire_params, const bool print_skeleton = true) noexcept
            :
            gate_layout{gate_lyt.clone()},
            sim_params{simulation_parameters},
            canvas{rel_canvas},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            output_perturbers{get_output_perturbers(skeleton, gate_lyt)},
            gate_connections{get_gate_connections(bdl_wires, skeleton, gate_lyt)}
    {
        make_skeleton_with_canvasses(gate_layout, skeleton, canvas, print_skeleton);
    }
    /**
     * SiDB gate-level layout.
     */
    const GateLyt gate_layout;

    const sidb_simulation_parameters sim_params;

    const std::pair<cell<CellLyt>, cell<CellLyt>> canvas;

    const CellLyt skeleton;

    const std::vector<bdl_wire<CellLyt>>                       bdl_wires{};
    const uint64_t                                             num_bdl_pairs{};
    const std::vector<bdl_pair<cell<CellLyt>>>                 input_bdl_pairs{};
    const uint64_t                                             num_inputs{};
    const std::vector<cell<CellLyt>>                           output_perturbers{};
    const std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};

    [[nodiscard]] static cell<CellLyt> relative_to_absolute_canvas_position(const GateLyt&       gate_lyt,
                                                                            const cell<CellLyt>& rel_pos,
                                                                            const tile<GateLyt>& t) noexcept
    {
        cell<CellLyt> absolute_c =
            relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(), SkeletonGateLibrary::gate_y_size(),
                                               GateLyt, CellLyt>(gate_lyt, t, rel_pos);

        const bool nw = gate_lyt.has_north_western_incoming_signal(t);
        const bool ne = gate_lyt.has_north_eastern_incoming_signal(t);

        if (nw && !ne)
        {
            absolute_c.x -= SkeletonGateLibrary::gate_x_size() / 6;
        }

        if (!nw && ne)
        {
            absolute_c.x += SkeletonGateLibrary::gate_x_size() / 6;
        }

        return absolute_c;
    }

    static CellLyt make_skeleton_with_canvasses(
        const GateLyt& gate_lyt, const CellLyt& skeleton, const std::pair<cell<CellLyt>, cell<CellLyt>>& canvas,
        const bool                                       print          = true,
        const std::optional<std::vector<tile<GateLyt>>>& tile_whitelist = std::nullopt) noexcept
    {
        CellLyt lyt = skeleton.clone();

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const auto& t = gate_lyt.get_tile(n);

                if (tile_whitelist.has_value() &&
                    std::find(tile_whitelist->cbegin(), tile_whitelist->cend(), t) == tile_whitelist->cend())
                {
                    return;
                }

                for (const cell<CellLyt>& relative_c : all_coordinates_in_spanned_area(canvas.first, canvas.second))
                {
                    lyt.assign_cell_type(relative_to_absolute_canvas_position(gate_lyt, relative_c, t),
                                         sidb_technology::cell_type::LOGIC);
                }
            });

        if (print)
        {
            std::cout << "Skeleton looks like:" << std::endl;
            print_layout(lyt);
            std::cout << std::endl;
        }

        return lyt;
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

    /**
     *
     */
    [[nodiscard]] static std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>>
    get_gate_connections(const std::vector<bdl_wire<CellLyt>>& bdl_wires, const CellLyt& lyt,
                         const GateLyt& gate_lyt) noexcept
    {
        std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};
        gate_connections.reserve(bdl_wires.size());

        for (const bdl_wire<CellLyt>& wire : bdl_wires)
        {
            assert(!wire.pairs.empty() && "BDL wire is empty");

            const bool right_to_left = wire.pairs.front().upper.x > wire.pairs.back().lower.x;
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

  private:
    [[nodiscard]] static std::vector<cell<CellLyt>> get_output_perturbers(const CellLyt& skeleton,
                                                                          const GateLyt& gate_lyt)
    {
        std::vector<cell<CellLyt>> output_perturbers{};
        output_perturbers.reserve(gate_lyt.num_pos());

        skeleton.foreach_cell(
            [&](const cell<CellLyt>& c)
            {
                if (skeleton.get_cell_type(c) == sidb_technology::cell_type::OUTPUT_PERTURBER)
                {
                    output_perturbers.push_back(c);
                }
            });

        return output_perturbers;
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
            gate_connections{super_circuit.gate_connections}
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
            gate_connections{sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::get_gate_connections(bdl_wires, skeleton, gate_layout)},// filter_gate_connections_from_super_circuit(super_circuit, bdl_wires)},
            gate_tile{this_tile},
            super_circuit_simulated_bdl_wires_per_input{simulate_bdl_wires_of_super_circuit(
                super_circuit, tiles, input_bdl_pairs, num_inputs,
                sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::make_skeleton_with_canvasses(
                    gate_layout, skeleton, super_circuit.canvas, false))}
    {}

    const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& super_circuit{};

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

    [[nodiscard]] std::reference_wrapper<const std::optional<charge_distribution_surface<CellLyt>>>
    get_simulated_bdl_wires_for_input_indices(const uint64_t sub_circuit_input_index,
                                              const uint64_t super_circuit_input_index) const noexcept
    {
        return std::cref(
            super_circuit_simulated_bdl_wires_per_input.at(sub_circuit_input_index).at(super_circuit_input_index));
    }

    [[nodiscard]] double get_energy_of_expected_charge_distribution_with_sub_circuit_charge_distribution(
        const uint64_t sub_circuit_input_index, const uint64_t super_circuit_input_index,
        charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>& sub_circuit_charge_distribution,
        bool& at_least_one_physically_valid_ground_state_for_super_circuit_input_combination) const noexcept
    {
        const charge_distribution_surface<CellLyt>& simulated_bdl_wires =
            *get_simulated_bdl_wires_for_input_indices(sub_circuit_input_index, super_circuit_input_index).get();

        typename charge_distribution_surface<
            CellLyt, local_external_potential_type::BOUNDED>::local_external_potential_map_t& pot_from_super_circuit = sub_circuit_charge_distribution.get_local_external_potentials_reference();

        for (const cell<CellLyt>& c : sub_circuit_charge_distribution.get_sidb_order())
        {
            assert(simulated_bdl_wires.get_local_internal_potential(c).has_value() && "c is not part of the layout");

            pot_from_super_circuit[c] = {simulated_bdl_wires.get_local_internal_potential(c).value(),
                                         simulated_bdl_wires.get_local_internal_potential(c).value()};
        }

        sub_circuit_charge_distribution.update_local_external_potential();
        sub_circuit_charge_distribution.determine_effective_charge_transition_thresholds();
        sub_circuit_charge_distribution.recompute_electrostatic_potential_energy();
        sub_circuit_charge_distribution.validity_check();
        at_least_one_physically_valid_ground_state_for_super_circuit_input_combination |=
            sub_circuit_charge_distribution.is_physically_valid();

        return sub_circuit_charge_distribution.get_electrostatic_potential_energy()[0];

        // todo: defects ... ?
    }

  private:
    const std::vector<std::vector<std::optional<charge_distribution_surface<CellLyt>>>>
        super_circuit_simulated_bdl_wires_per_input{};

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

    [[nodiscard]] static std::vector<std::vector<std::optional<charge_distribution_surface<CellLyt>>>>
    simulate_bdl_wires_of_super_circuit(const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& super_circuit,
                                        const std::vector<tile<GateLyt>>&           sub_circuit_tiles,
                                        const std::vector<bdl_pair<cell<CellLyt>>>& input_bdl_pairs,
                                        const uint64_t                              num_sub_circuit_inputs,
                                        const CellLyt& sub_circuit_skeleton_with_canvasses) noexcept
    {
        std::vector<std::vector<std::optional<charge_distribution_surface<CellLyt>>>>
            super_circuit_simulated_bdl_wires_per_sub_circuit_input{};
        super_circuit_simulated_bdl_wires_per_sub_circuit_input.reserve(1 << num_sub_circuit_inputs);

        for (uint64_t sub_circuit_input_index = 0; sub_circuit_input_index < 1 << num_sub_circuit_inputs;
             ++sub_circuit_input_index)
        {
            std::vector<std::optional<charge_distribution_surface<CellLyt>>> super_circuit_simulated_bdl_wires{};

            for (uint64_t super_circuit_input_index = 0; super_circuit_input_index < 1 << super_circuit.num_inputs;
                 ++super_circuit_input_index)
            {
                super_circuit_simulated_bdl_wires.emplace_back();

                bool mismatch = false;

                charge_distribution_surface<CellLyt> current_cds = charge_distribution_surface<CellLyt>{
                    sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::make_skeleton_with_canvasses(
                        super_circuit.gate_layout, super_circuit.skeleton, super_circuit.canvas, false,
                        sub_circuit_tiles),
                    super_circuit.sim_params, sidb_charge_state::NEUTRAL, cds_configuration::CHARGE_LOCATION_ONLY};

                for (const cell<CellLyt>& p : super_circuit.output_perturbers)
                {
                    current_cds.assign_charge_state(p, sidb_charge_state::NEGATIVE,
                                                    charge_index_mode::KEEP_CHARGE_INDEX);
                }

                const auto assign_logic_state_to_bdl_pairs =
                    [&](const typename std::vector<bdl_pair<cell<CellLyt>>>::const_iterator begin,
                        const typename std::vector<bdl_pair<cell<CellLyt>>>::const_iterator end,
                        const bool                                                          signal) noexcept
                {
                    for (auto it = begin; it != end; ++it)

                    {
                        assert(it->type != sidb_technology::cell_type::INPUT &&
                               "input BDL pairs are handled separately");

                        current_cds.assign_charge_state(signal ? it->lower : it->upper, sidb_charge_state::NEGATIVE,
                                                        charge_index_mode::KEEP_CHARGE_INDEX);
                    }
                };

                std::unordered_map<tile<GateLyt>, std::unordered_map<tile<GateLyt>, bool>>
                    expected_signal_at_gate_connection{};

                uint64_t current_input_number             = 0;
                uint64_t current_sub_circuit_input_number = 0;

                const auto is_bit_set =
                    [&](const uint64_t input_index, uint64_t& input_number, const uint64_t number_of_inputs)
                { return (input_index & (uint64_t{1ull} << (number_of_inputs - 1 - input_number++))) != 0ull; };

                for (uint64_t wire_ix = 0; wire_ix < super_circuit.bdl_wires.size(); ++wire_ix)
                {
                    const bdl_wire<CellLyt>& wire = super_circuit.bdl_wires.at(wire_ix);

                    assert((wire.port.dir == port_direction::SOUTH || wire.port.dir == port_direction::EAST ||
                            wire.port.dir == port_direction::NONE) &&
                           "Wrong port direction; only row clocking is supported");

                    const auto& [upper_tile, lower_tile] = super_circuit.gate_connections.at(wire_ix);

                    typename std::vector<bdl_pair<cell<CellLyt>>>::const_iterator
                        assign_logic_state_to_bdl_pairs_start_it = wire.pairs.cbegin();

                    // if (circuit->gate_layout.is_pi_tile(upper_tile)) todo
                    if (wire.pairs.front().type == sidb_technology::cell_type::INPUT)
                    {
                        // assert((expected_signal_at_gate_connection.count(lower_tile) == 0 ||
                        //         expected_signal_at_gate_connection.at(lower_tile).count(upper_tile) == 0) &&
                        //        "PI is visited twice"); todo

                        const bool current_bit_set =
                            is_bit_set(super_circuit_input_index, current_input_number, super_circuit.num_inputs);

                        expected_signal_at_gate_connection[lower_tile].insert({upper_tile, current_bit_set});

                        const bdl_pair<cell<CellLyt>>& input_pair = wire.pairs.front();

                        assert(input_pair.type == sidb_technology::cell_type::INPUT &&
                               "BDL wire connecting to a PI does not start with an input BDL pair");

                        current_cds.assign_charge_state(current_bit_set ? input_pair.lower : input_pair.upper,
                                                        sidb_charge_state::NEGATIVE,
                                                        charge_index_mode::KEEP_CHARGE_INDEX);

                        assign_logic_state_to_bdl_pairs_start_it = std::next(wire.pairs.cbegin(), 1);
                    }

                    // assert(expected_signal_at_gate_connection.count(lower_tile) != 0 &&
                    //        expected_signal_at_gate_connection.at(lower_tile).count(upper_tile) != 0 &&
                    //        "Tile is visited before the incoming tile that connects it"); todo

                    const bool expected_signal_for_wire =
                        expected_signal_at_gate_connection.at(lower_tile).at(upper_tile);

                    if (current_sub_circuit_input_number < num_sub_circuit_inputs &&
                        wire.pairs.front().upper == input_bdl_pairs.at(current_sub_circuit_input_number).upper)
                    {
                        const bool bitset = is_bit_set(sub_circuit_input_index, current_sub_circuit_input_number,
                                                       num_sub_circuit_inputs);
                        if (expected_signal_for_wire != bitset)
                        {
                            // input mismatches with circuit

                            mismatch = true;

                            break;
                        }
                    }

                    assign_logic_state_to_bdl_pairs(assign_logic_state_to_bdl_pairs_start_it, wire.pairs.cend(),
                                                    expected_signal_for_wire);

                    const uint32_t num_inputs_to_lower_tile =
                        super_circuit.gate_layout.node_function(super_circuit.gate_layout.get_node(lower_tile))
                            .num_vars();

                    assert(expected_signal_at_gate_connection.at(lower_tile).size() <= num_inputs_to_lower_tile &&
                           "Number of tiles visited connecting to the current tile exceeds the number of inputs to the "
                           "node function");

                    if (super_circuit.gate_layout.is_po_tile(lower_tile) ||
                        expected_signal_at_gate_connection.at(lower_tile).size() < num_inputs_to_lower_tile)
                    {
                        continue;
                    }

                    const std::vector<tile<GateLyt>>& outgoing_tiles =
                        super_circuit.gate_layout.outgoing_data_flow(lower_tile);

                    assert(!outgoing_tiles.empty() && "Non-PO tile does not have outgoing data flow");

                    if constexpr (has_is_fanout_v<GateLyt>)
                    {
                        if (super_circuit.gate_layout.is_fanout(super_circuit.gate_layout.get_node(lower_tile)))
                        {
                            for (const tile<GateLyt>& lower_lower_t : outgoing_tiles)
                            {
                                expected_signal_at_gate_connection[lower_lower_t].insert(
                                    {lower_tile, expected_signal_for_wire});
                            }

                            continue;
                        }
                    }

                    assert(outgoing_tiles.size() == 1 &&
                           "Tile with single-output gate has more than one outgoing tile.");

                    const tile<GateLyt>& first_input_tile =
                        expected_signal_at_gate_connection.at(lower_tile).cbegin()->first;

                    assert((num_inputs_to_lower_tile == 1 ||
                            (expected_signal_at_gate_connection.at(lower_tile).size() == 2 &&
                             first_input_tile.y ==
                                 std::next(expected_signal_at_gate_connection.at(lower_tile).cbegin(), 1)->first.y &&
                             first_input_tile.x !=
                                 std::next(expected_signal_at_gate_connection.at(lower_tile).cbegin(), 1)->first.x)) &&
                           "Only row clocking is supported; inputs to a tile must be on the same y with differing x");

                    auto tt_inp =
                        static_cast<uint8_t>(expected_signal_at_gate_connection.at(lower_tile).at(first_input_tile));

                    if (num_inputs_to_lower_tile == 2)
                    {
                        const tile<GateLyt>& second_input_tile =
                            std::next(expected_signal_at_gate_connection.at(lower_tile).cbegin(), 1)->first;
                        const auto tt_inp_second = static_cast<uint8_t>(
                            expected_signal_at_gate_connection.at(lower_tile).at(second_input_tile));

                        // tt_inp <- 2 * L_in + R_in
                        if (first_input_tile.x < second_input_tile.x)
                        {
                            tt_inp = 2 * tt_inp + tt_inp_second;
                        }
                        else
                        {
                            tt_inp += 2 * tt_inp_second;
                        }
                    }

                    expected_signal_at_gate_connection[outgoing_tiles.front()].insert(
                        {lower_tile, kitty::get_bit(super_circuit.gate_layout.node_function(
                                                        super_circuit.gate_layout.get_node(lower_tile)),
                                                    tt_inp)});
                }

                if (mismatch)
                {
                    continue;
                }

                // the charge index is not updated

                sub_circuit_skeleton_with_canvasses.foreach_cell(
                    [&](const auto& c)
                    {
                        if (const auto ct = sub_circuit_skeleton_with_canvasses.get_cell_type(c);
                            ct != sidb_technology::cell_type::OUTPUT_PERTURBER ||
                            super_circuit.skeleton.get_cell_type(c) == sidb_technology::cell_type::OUTPUT_PERTURBER)
                        {
                            current_cds.assign_charge_state(c, sidb_charge_state::NEUTRAL,
                                                            charge_index_mode::KEEP_CHARGE_INDEX);
                        }
                    });

                current_cds.initialize_matrices_for_electrostatic_calculation();
                current_cds.update_local_internal_potential();

                // todo: defect influence
                // if constexpr (is_sidb_defect_surface_v<CellLyt>)
                // {
                //     CellLyt::foreach_sidb_defect([&current_cds](const auto cd)
                //                                  { current_cds.add_sidb_defect_to_potential_landscape(cd.first,
                //                                  cd.second); });
                // }

                super_circuit_simulated_bdl_wires[super_circuit_input_index].emplace(std::move(current_cds));
            }

            super_circuit_simulated_bdl_wires_per_sub_circuit_input.push_back(
                std::move(super_circuit_simulated_bdl_wires));
        }

        return super_circuit_simulated_bdl_wires_per_sub_circuit_input;
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

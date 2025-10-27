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
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/traits.hpp"

#include <fiction/io/write_sqd_layout.hpp>

#include <phmap.h>

#include <cmath>
#include <limits>
#include <set>
#include <unordered_map>
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
                              const detect_bdl_wires_params&                 bdl_wire_params = {},
                              const bool                                     print_skeleton  = true) noexcept :
            gate_layout{gate_lyt.clone()},
            sim_params{simulation_parameters},
            canvas{rel_canvas},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt)},
            skeleton_with_canvasses{make_skeleton_with_canvasses(gate_layout, skeleton, canvas, print_skeleton)},
            output_perturbers{get_output_perturbers(skeleton, gate_lyt)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            gate_connections{get_gate_connections(bdl_wires, skeleton, gate_lyt)},
            simulated_bdl_wires_per_input{simulate_bdl_wires_for_each_input(gate_layout, sim_params,
                                                                            skeleton_with_canvasses, output_perturbers,
                                                                            bdl_wires, num_inputs, gate_connections)},
            all_cells_at_node_per_input{make_all_cells_at_node_map_per_input(gate_lyt, skeleton, canvas)}
    {
        write_sqd_layout(skeleton_with_canvasses, "test.sqd");
        set_initial_gate_design_influence_bounds();
    }
    /**
     * SiDB gate-level layout.
     */
    const GateLyt gate_layout;

    const sidb_simulation_parameters sim_params;

    const std::pair<cell<CellLyt>, cell<CellLyt>> canvas;

    const CellLyt skeleton;
    const CellLyt skeleton_with_canvasses;

    const std::vector<cell<CellLyt>>                           output_perturbers{};
    const std::vector<bdl_wire<CellLyt>>                       bdl_wires{};
    const uint64_t                                             num_bdl_pairs{};
    const std::vector<bdl_pair<cell<CellLyt>>>                 input_bdl_pairs{};
    const uint64_t                                             num_inputs{};
    const std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>> gate_connections{};

    /**
     * A canvas combination is a combination of canvas positions as a vector of canvas position indices.
     */
    using canvas_combination = std::vector<std::size_t>;

    std::unordered_map<mockturtle::node<GateLyt>, std::vector<cell<CellLyt>>>      all_canvas_positions{};
    std::unordered_map<mockturtle::node<GateLyt>, std::vector<canvas_combination>> gate_designs{};

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
            absolute_c.x -= 2;  // SkeletonGateLibrary::gate_x_size() / 8;
        }

        if (!nw && ne)
        {
            absolute_c.x += 2;  // SkeletonGateLibrary::gate_x_size() / 8;
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

                for (const cell<CellLyt>& c :
                     all_coordinates_in_spanned_area(relative_to_absolute_canvas_position(gate_lyt, canvas.first, t),
                                                     relative_to_absolute_canvas_position(gate_lyt, canvas.second, t)))
                {
                    lyt.assign_cell_type(c, sidb_technology::cell_type::LOGIC);
                    lyt.assign_cell_tile(c, {gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y, 0});
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
                assert(lyt.template get_cell_tile<tile<GateLyt>>(pair.upper) ==
                       lyt.template get_cell_tile<tile<GateLyt>>(pair.lower));

                tiles_in_wire.insert(
                    [&](const auto& this_t)
                    {
                        const tile<GateLyt> t = {this_t.x, this_t.y};
                        auto                z = 0;

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

                        return tile<GateLyt>{t.x, t.y, z};
                    }(lyt.template get_cell_tile<tile<GateLyt>>(pair.upper)));
            }

            assert((tiles_in_wire.size() == 2 ||
                    (tiles_in_wire.size() == 1 &&
                     (is_input_wire || wire.pairs.back().type == sidb_technology::cell_type::OUTPUT))) &&
                   "Gate connection is malformed");

            const auto get_tile_pair = [&]
            {
                if (tiles_in_wire.size() == 2)
                {
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

                    if (const std::vector<tile<GateLyt>>& incoming_tiles = gate_lyt.incoming_data_flow(tile_in_wire);
                        incoming_tiles.size() == 1)
                    {
                        z = incoming_tiles.front().z;
                    }
                    else if (incoming_tiles.size() == 2)
                    {
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

                return std::make_pair(tile_in_wire, tile<GateLyt>{lower_tile.x, lower_tile.y, z});
            };

            std::pair<tile<GateLyt>, tile<GateLyt>> tile_pair_at_gate_connection = get_tile_pair();

            gate_connections.push_back(std::move(tile_pair_at_gate_connection));
        }

        return gate_connections;
    }

    [[nodiscard]] std::reference_wrapper<
        const charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>
    get_simulated_bdl_wires_for_input_index(const uint64_t input_index) const noexcept
    {
        return std::cref(simulated_bdl_wires_per_input.at(input_index));
    }

    [[nodiscard]] std::reference_wrapper<const std::array<double, 2>>
    get_gate_design_influence_bounds(const uint64_t input_index, const mockturtle::node<GateLyt>& n,
                                     const cell<CellLyt>& c, const mockturtle::node<GateLyt>& other_n) const noexcept
    {
        return std::cref(gate_design_influence_bounds_per_input.at(input_index).at(n).at(c).at(other_n));
    }

    void tighten_gate_influence_bounds_until_fixpoint(
        const uint64_t available_threads = std::thread::hardware_concurrency()) noexcept
    {
        std::cout << '\n';

        std::unordered_map<mockturtle::node<GateLyt>, uint64_t> gate_design_counts{};
        gate_design_counts.reserve(gate_designs.size());

        for (const auto& [n, designs] : gate_designs)
        {
            gate_design_counts[n] = designs.size();
        }

        std::vector<uint64_t> input_indices{};
        input_indices.reserve(1 << num_inputs);
        for (uint64_t i = 0; i < 1 << num_inputs; input_indices.push_back(i++))
        {}

        // tighten bounds until fixed point

        bool big_fixpoint = false;

        while (!big_fixpoint)
        {
            std::cout << "\nSTARTING FIXED POINT ITERATION " << std::endl;

            big_fixpoint = true;

            std::shuffle(input_indices.begin(), input_indices.end(), std::mt19937(std::random_device()()));

            for (const uint64_t input_index : input_indices)
            {
                std::cout << "NARROWING INFLUENCE BOUNDS FOR INPUT INDEX " << input_index << std::endl;

                bool fixpoint = false;

                while (!fixpoint)
                {
                    fixpoint = true;

                    gate_layout.foreach_node(
                        [&](const auto& n)
                        {
                            if (skip_physical_design_for_node(gate_layout, n))
                            {
                                return;
                            }

                            std::unordered_map<mockturtle::node<GateLyt>, std::set<uint64_t>> gate_nums_to_prune{};
                            std::mutex mutex_to_protect_gate_nums_to_prune{};

                            const uint64_t num_cells = all_cells_at_node_per_input.at(input_index).at(n).size();

                            const uint64_t num_threads = std::min(available_threads, num_cells);

                            const uint64_t chunk_size =
                                (num_cells + num_threads - 1) / num_threads;  // Ceiling division

                            std::vector<std::thread> threads{};
                            threads.reserve(num_threads);

                            for (uint64_t i = 0; i < num_threads; ++i)
                            {
                                threads.emplace_back(
                                    [&, i]
                                    {
                                        const uint64_t cell_start_index = i * chunk_size;
                                        const uint64_t cell_end_index =
                                            std::min(cell_start_index + chunk_size, num_cells);

                                        tighten_gate_influence_bounds(input_index, n, cell_start_index, cell_end_index,
                                                                      gate_nums_to_prune,
                                                                      mutex_to_protect_gate_nums_to_prune, fixpoint);
                                    });
                            }

                            for (auto& thread : threads)
                            {
                                if (thread.joinable())
                                {
                                    thread.join();
                                }
                            }

                            gate_layout.foreach_node(
                                [&](const auto& other_n)
                                {
                                    if (skip_physical_design_for_node(gate_layout, other_n))
                                    {
                                        return;
                                    }

                                    if (gate_nums_to_prune.count(other_n) == 0)
                                    {
                                        return;
                                    }

                                    const tile<GateLyt>& other_t = gate_layout.get_tile(other_n);

                                    if (gate_designs.at(other_n).size() == gate_nums_to_prune.at(other_n).size())
                                    {
                                        throw std::runtime_error{
                                            fmt::format("All gate designs pruned for tile {}", other_t)};
                                    }

                                    for (auto it = gate_nums_to_prune.at(other_n).rbegin();
                                         it != gate_nums_to_prune.at(other_n).rend(); ++it)
                                    {
                                        std::swap(gate_designs[other_n][*it], gate_designs[other_n].back());

                                        gate_designs[other_n].pop_back();
                                    }

                                    std::cout << "pruned " << gate_nums_to_prune.at(other_n).size() << " from tile "
                                              << other_t << " (remaining: " << gate_designs.at(other_n).size() << ")"
                                              << std::endl;

                                    fixpoint     = false;
                                    big_fixpoint = false;
                                });
                        });
                }
            }
        }

        uint64_t pruned_total = 0;

        std::cout << "\n=================\n" << "  TILE  | #PRUNED" << std::endl;

        for (const auto& [n, designs] : gate_designs)
        {
            const uint64_t pruned = gate_design_counts.at(n) - designs.size();

            std::cout << gate_layout.get_tile(n) << " |    " << pruned << std::endl;

            pruned_total += pruned;
        }

        std::cout << "----------------- +\n" << "             " << pruned_total << std::endl;
    }

  private:
    const std::vector<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>
        simulated_bdl_wires_per_input;

    const std::vector<std::unordered_map<mockturtle::node<GateLyt>, std::vector<cell<CellLyt>>>>
        all_cells_at_node_per_input;

    // forall super circuit input index, for all nodes, for all cells at node, for all other nodes, bounds on influence
    // from possible gate designs of other node (excl. skeleton)
    std::vector<std::unordered_map<
        mockturtle::node<GateLyt>,
        phmap::flat_hash_map<cell<CellLyt>, std::unordered_map<mockturtle::node<GateLyt>, std::array<double, 2>>>>>
        gate_design_influence_bounds_per_input;

    [[nodiscard]] static std::vector<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>
    simulate_bdl_wires_for_each_input(
        const GateLyt& gate_layout, const sidb_simulation_parameters& sim_params,
        const CellLyt& skeleton_with_canvasses, const std::vector<cell<CellLyt>>& output_perturbers,
        const std::vector<bdl_wire<CellLyt>>& bdl_wires, const uint64_t num_inputs,
        const std::vector<std::pair<tile<GateLyt>, tile<GateLyt>>>& gate_connections) noexcept
    {
        std::vector<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>
            simulated_bdl_wires_per_input{};
        simulated_bdl_wires_per_input.reserve(1 << num_inputs);

        bdl_input_iterator<CellLyt> bii{skeleton_with_canvasses};

        for (uint64_t input_index = 0; input_index < 1 << num_inputs; ++input_index, ++bii)
        {
            simulated_bdl_wires_per_input.emplace_back((*bii).clone(), sim_params, sidb_charge_state::NONE,
                                                       cds_configuration::CHARGE_LOCATION_ONLY);

            charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>& cds =
                simulated_bdl_wires_per_input.back();
            cds.initialize_matrices_for_electrostatic_calculation();
            cds.get_local_external_potentials_reference().reserve(cds.num_cells());

            for (const cell<CellLyt>& p : output_perturbers)
            {
                cds.assign_charge_state(p, sidb_charge_state::NEGATIVE, charge_index_mode::KEEP_CHARGE_INDEX);
            }

            const auto assign_logic_state_to_bdl_pairs =
                [&](const typename std::vector<bdl_pair<cell<CellLyt>>>::const_iterator begin,
                    const typename std::vector<bdl_pair<cell<CellLyt>>>::const_iterator end, const bool signal) noexcept
            {
                for (auto it = begin; it != end; ++it)

                {
                    assert(it->type != sidb_technology::cell_type::INPUT && "input BDL pairs are handled separately");

                    cds.assign_charge_state(signal ? it->upper : it->lower, sidb_charge_state::NEUTRAL,
                                            charge_index_mode::KEEP_CHARGE_INDEX);
                    cds.assign_charge_state(signal ? it->lower : it->upper, sidb_charge_state::NEGATIVE,
                                            charge_index_mode::KEEP_CHARGE_INDEX);
                }
            };

            std::unordered_map<tile<GateLyt>, std::unordered_map<tile<GateLyt>, bool>>
                expected_signal_at_gate_connection{};

            uint64_t current_input_number = 0;

            const auto is_bit_set = [&](uint64_t& input_number, const uint64_t number_of_inputs)
            { return (input_index & (uint64_t{1ull} << (number_of_inputs - 1 - input_number++))) != 0ull; };

            for (uint64_t wire_ix = 0; wire_ix < bdl_wires.size(); ++wire_ix)
            {
                const bdl_wire<CellLyt>& wire = bdl_wires.at(wire_ix);

                assert((wire.port.dir == port_direction::SOUTH || wire.port.dir == port_direction::EAST ||
                        wire.port.dir == port_direction::NONE) &&
                       "Wrong port direction; only row clocking is supported");

                const auto& [upper_tile, lower_tile] = gate_connections.at(wire_ix);

                typename std::vector<bdl_pair<cell<CellLyt>>>::const_iterator assign_logic_state_to_bdl_pairs_start_it =
                    wire.pairs.cbegin();

                if (wire.pairs.front().type == sidb_technology::cell_type::INPUT)
                {
                    const bool current_bit_set = is_bit_set(current_input_number, num_inputs);

                    expected_signal_at_gate_connection[lower_tile].insert({upper_tile, current_bit_set});

                    const bdl_pair<cell<CellLyt>>& input_pair = wire.pairs.front();

                    assert(input_pair.type == sidb_technology::cell_type::INPUT &&
                           "BDL wire connecting to a PI does not start with an input BDL pair");

                    cds.assign_charge_state(current_bit_set ? input_pair.upper : input_pair.lower,
                                            sidb_charge_state::NEUTRAL, charge_index_mode::KEEP_CHARGE_INDEX);
                    cds.assign_charge_state(current_bit_set ? input_pair.lower : input_pair.upper,
                                            sidb_charge_state::NEGATIVE, charge_index_mode::KEEP_CHARGE_INDEX);

                    assign_logic_state_to_bdl_pairs_start_it = std::next(wire.pairs.cbegin(), 1);
                }

                const bool expected_signal_for_wire = expected_signal_at_gate_connection.at(lower_tile).at(upper_tile);

                assign_logic_state_to_bdl_pairs(assign_logic_state_to_bdl_pairs_start_it, wire.pairs.cend(),
                                                expected_signal_for_wire);

                const uint32_t num_inputs_to_lower_tile =
                    gate_layout.node_function(gate_layout.get_node(lower_tile)).num_vars();

                assert(expected_signal_at_gate_connection.at(lower_tile).size() <= num_inputs_to_lower_tile &&
                       "Number of tiles visited connecting to the current tile exceeds the number of inputs to the "
                       "node function");

                if (gate_layout.is_po_tile(lower_tile) ||
                    expected_signal_at_gate_connection.at(lower_tile).size() < num_inputs_to_lower_tile)
                {
                    continue;
                }

                const std::vector<tile<GateLyt>>& outgoing_tiles = gate_layout.outgoing_data_flow(lower_tile);

                assert(!outgoing_tiles.empty() && "Non-PO tile does not have outgoing data flow");

                if constexpr (has_is_fanout_v<GateLyt>)
                {
                    if (gate_layout.is_fanout(gate_layout.get_node(lower_tile)))
                    {
                        for (const tile<GateLyt>& lower_lower_t : outgoing_tiles)
                        {
                            expected_signal_at_gate_connection[lower_lower_t].insert(
                                {lower_tile, expected_signal_for_wire});
                        }

                        continue;
                    }
                }

                assert(outgoing_tiles.size() == 1 && "Tile with single-output gate has more than one outgoing tile.");

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
                    const auto tt_inp_second =
                        static_cast<uint8_t>(expected_signal_at_gate_connection.at(lower_tile).at(second_input_tile));

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
                    {lower_tile, kitty::get_bit(gate_layout.node_function(gate_layout.get_node(lower_tile)), tt_inp)});
            }
        }

        return simulated_bdl_wires_per_input;
    }

    [[nodiscard]] static std::vector<std::unordered_map<mockturtle::node<GateLyt>, std::vector<cell<CellLyt>>>>
    make_all_cells_at_node_map_per_input(const GateLyt& gate_lyt, const CellLyt& skeleton,
                                         const std::pair<cell<CellLyt>, cell<CellLyt>>& canvas) noexcept
    {
        std::vector<std::unordered_map<mockturtle::node<GateLyt>, std::vector<cell<CellLyt>>>>
            all_cells_at_node_map_per_input{};
        all_cells_at_node_map_per_input.reserve(1 << gate_lyt.num_pis());

        bdl_input_iterator<CellLyt> bii{skeleton};

        for (auto i = 0u; i < 1 << gate_lyt.num_pis(); ++i, ++bii)
        {
            std::unordered_map<mockturtle::node<GateLyt>, std::vector<cell<CellLyt>>> all_cells_at_node_map{};

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(gate_lyt, n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = gate_lyt.get_tile(n);

                    std::vector<cell<CellLyt>> cells = all_coordinates_in_spanned_area(
                        relative_to_absolute_canvas_position(gate_lyt, canvas.first, t),
                        relative_to_absolute_canvas_position(gate_lyt, canvas.second, t));

                    (*bii).foreach_cell(
                        [&](const cell<CellLyt>& c)
                        {
                            if ((*bii).template get_cell_tile<tile<GateLyt>>(c) == t)
                            {
                                cells.push_back(c);
                            }
                        });

                    all_cells_at_node_map[n] = std::move(cells);
                });

            all_cells_at_node_map_per_input.push_back(std::move(all_cells_at_node_map));
        }

        return all_cells_at_node_map_per_input;
    }

    void
    tighten_gate_influence_bounds(const uint64_t input_index, const mockturtle::node<GateLyt>& n,
                                  const uint64_t cell_start_index, const uint64_t cell_end_index,
                                  std::unordered_map<mockturtle::node<GateLyt>, std::set<uint64_t>>& gate_nums_to_prune,
                                  std::mutex& mutex_to_protect_gate_nums_to_prune, bool& fixpoint) noexcept
    {
        const auto update_bounds = [&](std::array<double, 2>& current_bounds, const std::array<double, 2>& new_bounds)
        {
            if (std::isinf(current_bounds[0]))
            {
                fixpoint = false;

                current_bounds = new_bounds;

                return;
            }

            if (current_bounds[0] - new_bounds[0] < -std::numeric_limits<double>::epsilon())
            {
                fixpoint = false;

                current_bounds[0] = new_bounds[0];
            }

            if (current_bounds[1] - new_bounds[1] > std::numeric_limits<double>::epsilon())
            {
                fixpoint = false;

                current_bounds[1] = new_bounds[1];
            }
        };

        charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED> simulated_bdl_wires{
            simulated_bdl_wires_per_input.at(input_index)};  // todo test whether .clone() matters here

        gate_layout.foreach_node(
            [&](const auto& other_n)
            {
                if (skip_physical_design_for_node(gate_layout, other_n))
                {
                    return;
                }

                typename charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>::
                    local_external_potential_map_t& bounded_influence_from_other_canvasses =
                        simulated_bdl_wires.get_local_external_potentials_reference();

                const auto collect_bounds =
                    [&](const cell<CellLyt>& c, const mockturtle::node<GateLyt>& containing_node)
                {
                    std::array<double, 2> bounds{0, 0};

                    gate_layout.foreach_node(
                        [&](const auto& other_other_n)
                        {
                            if (skip_physical_design_for_node(gate_layout, other_other_n) || other_n == other_other_n)
                            {
                                return;
                            }

                            const std::array<double, 2>& influence_bounds_from_other_n =
                                get_gate_design_influence_bounds(input_index, containing_node, c, other_other_n);

                            bounds[0] += influence_bounds_from_other_n[0];
                            bounds[1] += influence_bounds_from_other_n[1];
                        });

                    bounded_influence_from_other_canvasses[c] = std::move(bounds);
                };

                simulated_bdl_wires.foreach_cell(
                    [&](const auto& c)
                    {
                        collect_bounds(
                            c, gate_layout.get_node(simulated_bdl_wires.template get_cell_tile<tile<GateLyt>>(c)));
                    });

                simulated_bdl_wires.update_local_external_potential();
                simulated_bdl_wires.determine_effective_charge_transition_thresholds();

                for (uint64_t j = cell_start_index; j < cell_end_index; ++j)
                {
                    const cell<CellLyt>& c = all_cells_at_node_per_input.at(input_index).at(n).at(j);

                    std::array<double, 2> bounds{std::numeric_limits<double>::infinity(),
                                                 -std::numeric_limits<double>::infinity()};

                    uint64_t gate_num = 0;
                    for (const canvas_combination& gate : gate_designs.at(other_n))
                    {
                        CellLyt canvas_of_other_n{};

                        for (const uint64_t gate_design_cell_index : gate)
                        {
                            if (const cell<CellLyt>& canvas_c =
                                    all_canvas_positions.at(other_n).at(gate_design_cell_index);
                                canvas_c != c)
                            {
                                canvas_of_other_n.assign_cell_type(canvas_c, sidb_technology::cell_type::LOGIC);
                            }
                        }

                        charge_distribution_surface<CellLyt> canvas_cds{canvas_of_other_n, sim_params,
                                                                        sidb_charge_state::NEGATIVE,
                                                                        cds_configuration::CHARGE_LOCATION_ONLY};
                        canvas_cds.initialize_matrices_for_electrostatic_calculation();

                        canvas_cds.add_sidb_defect_to_potential_landscape(
                            c, sidb_defect{sidb_defect_type::DB, 0, sim_params.epsilon_r, sim_params.lambda_tf});

                        const auto max_index = canvas_cds.get_max_charge_index();

                        bool at_least_one_charge_index_valid = false;

                        for (uint64_t charge_index = 0; charge_index <= max_index; charge_index++)
                        {
                            canvas_cds.assign_charge_index(charge_index,
                                                           charge_distribution_mode::UPDATE_CHARGE_DISTRIBUTION);

                            canvas_cds.foreach_cell(
                                [&](const cell<CellLyt>& canvas_c)
                                {
                                    simulated_bdl_wires.assign_charge_state(canvas_c,
                                                                            canvas_cds.get_charge_state(canvas_c),
                                                                            charge_index_mode::KEEP_CHARGE_INDEX);
                                });

                            simulated_bdl_wires.template update_after_charge_change<true>(
                                dependent_cell_mode::FIXED, energy_calculation::KEEP_OLD_ENERGY_VALUE);

                            if (!simulated_bdl_wires.is_physically_valid())
                            {
                                continue;
                            }

                            at_least_one_charge_index_valid = true;

                            canvas_cds.update_local_defect_potential();

                            const double pot_at_c = *canvas_cds.get_local_defect_potential(c);

                            bounds[0] = std::min(bounds[0], pot_at_c);
                            bounds[1] = std::max(bounds[1], pot_at_c);
                        }

                        canvas_cds.foreach_cell(
                            [&](const cell<CellLyt>& canvas_c)
                            {
                                simulated_bdl_wires.assign_charge_state(canvas_c, sidb_charge_state::NONE,
                                                                        charge_index_mode::KEEP_CHARGE_INDEX);
                            });

                        if (!at_least_one_charge_index_valid)
                        {
                            const std::lock_guard guard{mutex_to_protect_gate_nums_to_prune};

                            gate_nums_to_prune[other_n].emplace(gate_num);
                        }

                        ++gate_num;
                    }

                    update_bounds(gate_design_influence_bounds_per_input[input_index][n][c][other_n], bounds);
                }
            });
    }

    void set_initial_gate_design_influence_bounds() noexcept
    {
        for (uint64_t input_index = 0; input_index < 1 << num_inputs; ++input_index)
        {
            gate_design_influence_bounds_per_input.emplace_back();

            // initialize all bounds to (-inf,inf)

            gate_layout.foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(gate_layout, n))
                    {
                        return;
                    }

                    std::unordered_map<mockturtle::node<GateLyt>, std::array<double, 2>> bounds_per_node{};

                    gate_layout.foreach_node(
                        [&](const auto& other_n)
                        {
                            if (skip_physical_design_for_node(gate_layout, other_n))
                            {
                                return;
                            }

                            bounds_per_node[other_n] = std::array<double, 2>{-std::numeric_limits<double>::infinity(),
                                                                             std::numeric_limits<double>::infinity()};
                        });

                    for (const cell<CellLyt>& c : all_cells_at_node_per_input.at(input_index).at(n))
                    {
                        gate_design_influence_bounds_per_input[input_index][n][c] = bounds_per_node;
                    }
                });
        }
    }

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
            gate_connections{super_circuit.gate_connections},
            tiles{get_all_tiles(gate_layout)}
    {}

    explicit sidb_bdl_sub_circuit(const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& bdl_super_circuit,
                                  const std::vector<tile<GateLyt>>&                              sub_circuit_tiles,
                                  const detect_bdl_wires_params& bdl_wire_params = {}) noexcept :
            super_circuit{bdl_super_circuit},
            gate_layout{create_gate_lyt_window_for_tiles(super_circuit.gate_layout, sub_circuit_tiles)},
            skeleton{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout)},
            bdl_wires{detect_bdl_wires(skeleton, bdl_wire_params)},
            num_bdl_pairs{sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>::get_number_of_bdl_pairs(bdl_wires)},
            input_bdl_pairs{detect_bdl_pairs<CellLyt>(skeleton, sidb_technology::cell_type::INPUT,
                                                      bdl_wire_params.bdl_pairs_params)},
            num_inputs{input_bdl_pairs.size()},
            tiles{sub_circuit_tiles},
            consistent_super_circuit_input_indices_per_input{
                collect_consistent_super_circuit_input_indices(super_circuit, input_bdl_pairs)},
            super_circuit_simulated_bdl_wires_per_input{simulate_bdl_wires_of_super_circuit(
                super_circuit, input_bdl_pairs, consistent_super_circuit_input_indices_per_input,
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

    const std::vector<tile<GateLyt>>   tiles{};

    [[nodiscard]] std::optional<sidb_technology::cell_type>
    is_not_internal_output_perturber(const CellLyt& lyt, const cell<CellLyt>& c) const noexcept
    {
        if (const auto ct = lyt.get_cell_type(c);
            ct != sidb_technology::cell_type::OUTPUT_PERTURBER ||
            super_circuit.skeleton.get_cell_type(c) == sidb_technology::cell_type::OUTPUT_PERTURBER)
        {
            return ct;
        }

        return std::nullopt;
    }

    [[nodiscard]] bool input_index_possible_in_super_circuit(const uint64_t input_index) const noexcept
    {
        return !consistent_super_circuit_input_indices_per_input.at(input_index).empty();
    }
    [[nodiscard]] std::reference_wrapper<const std::vector<uint64_t>>
    get_consistent_super_circuit_input_indices(const uint64_t input_index) const noexcept
    {
        return std::cref(consistent_super_circuit_input_indices_per_input.at(input_index));
    }
    [[nodiscard]] double get_skeleton_influence(const uint64_t       sub_circuit_input_index,
                                                const uint64_t       super_circuit_input_index,
                                                const cell<CellLyt>& c) const noexcept
    {
        return *super_circuit_simulated_bdl_wires_per_input.at(sub_circuit_input_index)
                    .at(super_circuit_input_index)
                    ->get_local_internal_potential(c);
    }

    typename charge_distribution_surface<CellLyt,
                                         local_external_potential_type::BOUNDED>::local_external_potential_map_t
    collect_influence_bounds(const CellLyt& lyt, const uint64_t input_index,
                             typename charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>::
                                 local_external_potential_map_t& influence_bounds) const noexcept
    {
        const std::vector<uint64_t>& consistent_super_circuit_input_indices =
            consistent_super_circuit_input_indices_per_input.at(input_index);

        influence_bounds.reserve(lyt.num_cells());

        typename charge_distribution_surface<CellLyt,
                                             local_external_potential_type::BOUNDED>::local_external_potential_map_t
            skeleton_influence_bounds{};
        skeleton_influence_bounds.reserve(lyt.num_cells());

        lyt.foreach_cell(
            [&](const auto& c)
            {
                if (!is_not_internal_output_perturber(lyt, c))
                {
                    return;
                }

                const mockturtle::node<GateLyt>& n = super_circuit.gate_layout.get_node(
                    super_circuit.skeleton_with_canvasses.template get_cell_tile<tile<GateLyt>>(c));

                std::array<double, 2> bounds{std::numeric_limits<double>::infinity(),
                                             -std::numeric_limits<double>::infinity()};

                std::array<double, 2> skeleton_bounds{std::numeric_limits<double>::infinity(),
                                                      -std::numeric_limits<double>::infinity()};

                for (uint64_t super_circuit_input_index_ix = 0;
                     super_circuit_input_index_ix < consistent_super_circuit_input_indices.size();
                     ++super_circuit_input_index_ix)
                {
                    const charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>&
                        simulated_bdl_wires =
                            *super_circuit_simulated_bdl_wires_per_input.at(input_index)
                                 .at(consistent_super_circuit_input_indices.at(super_circuit_input_index_ix));

                    assert(simulated_bdl_wires.get_local_internal_potential(c).has_value() &&
                           "c is not part of the layout");

                    const double skeleton_influence = *simulated_bdl_wires.get_local_internal_potential(c);

                    std::array<double, 2> gate_design_influence_bound_sum = {0, 0};

                    super_circuit.gate_layout.foreach_node(
                        [&](const auto& super_circuit_n)
                        {
                            if (skip_physical_design_for_node(super_circuit.gate_layout, super_circuit_n) ||
                                std::find(tiles.cbegin(), tiles.cend(),
                                          super_circuit.gate_layout.get_tile(super_circuit_n)) != tiles.cend())
                            {
                                return;
                            }

                            const std::array<double, 2>& influence_bounds_from_super_circuit_n =
                                super_circuit.get_gate_design_influence_bounds(
                                    consistent_super_circuit_input_indices.at(super_circuit_input_index_ix), n, c,
                                    super_circuit_n);

                            gate_design_influence_bound_sum[0] += influence_bounds_from_super_circuit_n[0];
                            gate_design_influence_bound_sum[1] += influence_bounds_from_super_circuit_n[1];
                        });

                    bounds[0] = std::min(bounds[0], skeleton_influence + gate_design_influence_bound_sum[0]);
                    bounds[1] = std::max(bounds[1], skeleton_influence + gate_design_influence_bound_sum[1]);

                    skeleton_bounds[0] = std::min(skeleton_bounds[0], skeleton_influence);
                    skeleton_bounds[1] = std::max(skeleton_bounds[1], skeleton_influence);
                }

                influence_bounds[c]          = std::move(bounds);
                skeleton_influence_bounds[c] = std::move(skeleton_bounds);
            });

        return skeleton_influence_bounds;
    }

  private:
    const std::vector<std::vector<uint64_t>> consistent_super_circuit_input_indices_per_input{};
    const std::vector<
        std::vector<std::optional<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>>>
        super_circuit_simulated_bdl_wires_per_input{};

    [[nodiscard]] static std::vector<tile<GateLyt>> get_all_tiles(const GateLyt& gate_lyt) noexcept
    {
        std::vector<tile<GateLyt>> tiles{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (!skip_physical_design_for_node(gate_lyt, n))
                {
                    tiles.push_back(gate_lyt.get_tile(n));
                }
            });

        return tiles;
    }

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
                while (!gate_lyt_window.is_empty_tile(in_t))
                {
                    if (in_t.z == gate_lyt_window.z())
                    {
                        gate_lyt_window.resize(
                            aspect_ratio<GateLyt>{gate_lyt_window.x(), gate_lyt_window.y(), gate_lyt_window.z() + 1});
                    }

                    ++in_t.z;
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
                while (!gate_lyt_window.is_empty_tile(out_t))
                {
                    if (out_t.z == gate_lyt_window.z())
                    {
                        gate_lyt_window.resize(
                            aspect_ratio<GateLyt>{gate_lyt_window.x(), gate_lyt_window.y(), gate_lyt_window.z() + 1});
                    }

                    ++out_t.z;
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

    [[nodiscard]] static GateLyt create_gate_lyt_window_for_tiles(const GateLyt&             gate_lyt,
                                                                  std::vector<tile<GateLyt>> tiles) noexcept
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

    [[nodiscard]] static std::vector<std::vector<uint64_t>> collect_consistent_super_circuit_input_indices(
        const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& super_circuit,
        const std::vector<bdl_pair<cell<CellLyt>>>&                    input_bdl_pairs) noexcept
    {
        const auto is_bit_set = [&](const uint64_t input_index, const uint64_t input_number)
        { return (input_index & (uint64_t{1ull} << (input_bdl_pairs.size() - 1 - input_number))) != 0ull; };

        std::vector<std::vector<uint64_t>> consistent_super_circuit_input_indices_per_sub_circuit_index{};
        consistent_super_circuit_input_indices_per_sub_circuit_index.reserve(1 << input_bdl_pairs.size());

        for (uint64_t sub_circuit_iix = 0; sub_circuit_iix < 1 << input_bdl_pairs.size(); ++sub_circuit_iix)
        {
            std::vector<uint64_t> consistent_super_circuit_input_indices{};
            consistent_super_circuit_input_indices.reserve(1 << super_circuit.num_inputs);

            for (uint64_t super_circuit_iix = 0; super_circuit_iix < 1 << super_circuit.num_inputs; ++super_circuit_iix)
            {
                const charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>&
                    simulated_bdl_wires = super_circuit.get_simulated_bdl_wires_for_input_index(super_circuit_iix);

                bool pass = true;

                for (uint64_t sub_circuit_input_number = 0; sub_circuit_input_number < input_bdl_pairs.size();
                     ++sub_circuit_input_number)
                {
                    const bool require_1 = is_bit_set(sub_circuit_iix, sub_circuit_input_number);
                    const bool is_1 =
                        simulated_bdl_wires.get_charge_state(input_bdl_pairs.at(sub_circuit_input_number).lower) ==
                        sidb_charge_state::NEGATIVE;

                    if (!((require_1 && is_1) || (!require_1 && !is_1)))
                    {
                        pass = false;
                        break;
                    }
                }

                if (pass)
                {
                    consistent_super_circuit_input_indices.push_back(super_circuit_iix);
                }
            }

            consistent_super_circuit_input_indices_per_sub_circuit_index.push_back(
                std::move(consistent_super_circuit_input_indices));
        }

        return consistent_super_circuit_input_indices_per_sub_circuit_index;
    }

    [[nodiscard]] static std::vector<
        std::vector<std::optional<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>>>
    simulate_bdl_wires_of_super_circuit(
        const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>& super_circuit,
        const std::vector<bdl_pair<cell<CellLyt>>>&                    input_bdl_pairs,
        const std::vector<std::vector<uint64_t>>&                      consistent_super_circuit_input_indices_per_input,
        const CellLyt&                                                 sub_circuit_skeleton_with_canvasses) noexcept
    {
        std::vector<
            std::vector<std::optional<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>>>
            super_circuit_simulated_bdl_wires_per_sub_circuit_input{};
        super_circuit_simulated_bdl_wires_per_sub_circuit_input.reserve(1 << input_bdl_pairs.size());

        for (uint64_t sub_circuit_input_index = 0; sub_circuit_input_index < 1 << input_bdl_pairs.size();
             ++sub_circuit_input_index)
        {
            std::vector<std::optional<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>>
                super_circuit_simulated_bdl_wires{};
            super_circuit_simulated_bdl_wires.resize(1 << super_circuit.num_inputs);

            for (const uint64_t consistent_super_circuit_iix :
                 consistent_super_circuit_input_indices_per_input.at(sub_circuit_input_index))
            {
                charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED> simulated_bdl_wires_copy{
                    super_circuit.get_simulated_bdl_wires_for_input_index(consistent_super_circuit_iix)};

                // neutralise sub-circuit and store charge information in external potential

                simulated_bdl_wires_copy.initialize_matrices_for_electrostatic_calculation();

                // typename charge_distribution_surface<CellLyt>::local_external_potential_map_t
                //     pot_from_sub_circuit_skeleton_to_super_circuit{};
                // pot_from_sub_circuit_skeleton_to_super_circuit.reserve(cds_with_simulated_bdl_wires.num_cells() -
                //                                                        sub_circuit_skeleton_with_canvasses.num_cells());

                sub_circuit_skeleton_with_canvasses.foreach_cell(
                    [&](const auto& c)
                    {
                        if ((sub_circuit_skeleton_with_canvasses.get_cell_type(c) ==
                                 sidb_technology::cell_type::OUTPUT_PERTURBER &&
                             super_circuit.skeleton.get_cell_type(c) != sidb_technology::cell_type::OUTPUT_PERTURBER) ||
                            simulated_bdl_wires_copy.get_charge_state(c) == sidb_charge_state::NEUTRAL)
                        {
                            return;
                        }

                        // for (const cell<CellLyt>& other_c : cds_with_simulated_bdl_wires->get_sidb_order())
                        // {
                        //     if (other_c == c)
                        //     {
                        //         return;
                        //     }
                        //
                        //     pot_from_sub_circuit_skeleton_to_super_circuit[other_c] +=
                        //         cds_with_simulated_bdl_wires->get_potential_between_sidbs(other_c, c);
                        // }

                        simulated_bdl_wires_copy.assign_charge_state(c, sidb_charge_state::NEUTRAL,
                                                                     charge_index_mode::KEEP_CHARGE_INDEX);
                    });

                // cds_with_simulated_bdl_wires->assign_local_external_potential_map(
                //     pot_from_sub_circuit_skeleton_to_super_circuit);
                simulated_bdl_wires_copy.template update_local_internal_potential<true>();

                // todo: defect influence
                // if constexpr (is_sidb_defect_surface_v<CellLyt>)
                // {
                //     CellLyt::foreach_sidb_defect([&current_cds](const auto cd)
                //                                  { current_cds.add_sidb_defect_to_potential_landscape(cd.first,
                //                                  cd.second); });
                // }

                super_circuit_simulated_bdl_wires[consistent_super_circuit_iix].emplace(
                    std::move(simulated_bdl_wires_copy));
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

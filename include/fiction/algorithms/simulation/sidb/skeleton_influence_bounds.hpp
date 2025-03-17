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
    bool                    absolute_positions = false;
};

namespace detail
{

template <typename CellLyt, typename SkeletonGateLibrary, typename GateLyt>
class skeleton_influence_bounds_impl
{
  public:
    skeleton_influence_bounds_impl(const GateLyt& gate_layout, const std::vector<tile<GateLyt>>& ts,
                                   const skeleton_influence_bounds_params<cell<CellLyt>>& ps) noexcept :

            params{ps},
            gate_lyt{gate_layout},
            tiles_of_interest{ts},
            cell_lyt_with_all_skeletons{apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt)}
    // bdl_wires_of_designed_gates_current_tile_with_other_tile{
    //     obtain_bdl_wires_for_all_tile_pairs(gate_layout, ts.front())}
    {}

    [[nodiscard]] std::unordered_map<cell<CellLyt>, std::array<double, 2>> run() noexcept
    {
        for (const tile<GateLyt>& current_tile : tiles_of_interest)
        {
            if (current_tile.y == 0 && tiles_of_interest.size() == 1)
            {
                std::cout << "Skeleton looks like:" << std::endl;
                print_layout(apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt));
                std::cout << std::endl;
            }

            const CellLyt designed_gate =
                apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile});

            std::cout << "SKELETON OF GATE OF INTEREST:" << std::endl;
            print_layout(designed_gate);

            const charge_distribution_surface<CellLyt> cds{designed_gate};

            first_sidb  = cds.get_sidb_order().front();
            second_sidb = cds.get_sidb_order().front().y == cds.get_sidb_order().at(1).y ? cds.get_sidb_order().at(2) :
                                                                                           cds.get_sidb_order().at(1);

            designed_gate.foreach_cell(
                [&](const cell<CellLyt>& absolute_c)
                {
                    const cell<CellLyt> absolute_offset =
                        relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                           SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                            gate_lyt, current_tile, {0, 0});

                    collect_influence_bounds_for_sidb(
                        current_tile, {absolute_c.x - absolute_offset.x, absolute_c.y - absolute_offset.y}, absolute_c);
                });

            third_sidb = relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                            SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                gate_lyt, current_tile,
                all_coordinates_in_spanned_area(params.canvas.first, params.canvas.second).front());
            fourth_sidb = relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                             SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                gate_lyt, current_tile,
                all_coordinates_in_spanned_area(params.canvas.first, params.canvas.second).back());

            for (const cell<CellLyt>& relative_c :
                 all_coordinates_in_spanned_area(params.canvas.first, params.canvas.second))
            {
                const cell<CellLyt> absolute_c =
                    relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                       SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                        gate_lyt, current_tile, relative_c);
                collect_influence_bounds_for_sidb(current_tile, relative_c, absolute_c);
            }

            for (const auto& [cell_maps, sidb] :
                 std::array<std::pair<std::pair<std::unordered_map<cell<CellLyt>, sidb_charge_state>,
                                                std::unordered_map<cell<CellLyt>, sidb_charge_state>>,
                                      cell<CellLyt>>,
                            4>{{{{cell_map1_lb, cell_map1_ub}, first_sidb},
                                {{cell_map2_lb, cell_map2_ub}, second_sidb},
                                {{cell_map3_lb, cell_map3_ub}, third_sidb},
                                {{cell_map4_lb, cell_map4_ub}, fourth_sidb}}})
            {
                CellLyt lyt{};
                designed_gate.foreach_cell([&](const cell<CellLyt>& absolute_c)
                                           { lyt.assign_cell_type(absolute_c, sidb_technology::cell_type::NORMAL); });

                for (const cell<CellLyt>& relative_c :
                     all_coordinates_in_spanned_area(params.canvas.first, params.canvas.second))
                {
                    const cell<CellLyt> absolute_c =
                        relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                           SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                            gate_lyt, current_tile, relative_c);

                    lyt.assign_cell_type(absolute_c, sidb_technology::cell_type::NORMAL);
                }

                for (const auto& [c, _] : cell_maps.first)
                {
                    lyt.assign_cell_type(c, sidb_technology::cell_type::NORMAL);
                }

                charge_distribution_surface<CellLyt> cds_lb{lyt}, cds_ub{lyt};

                cds_lb.assign_all_charge_states(sidb_charge_state::NEUTRAL);
                cds_ub.assign_all_charge_states(sidb_charge_state::NEUTRAL);

                for (const auto& cc : {first_sidb, second_sidb, third_sidb, fourth_sidb})
                {
                    cds_lb.assign_charge_state(cc, cc == sidb ? sidb_charge_state::POSITIVE : sidb_charge_state::NONE);
                    cds_ub.assign_charge_state(cc, cc == sidb ? sidb_charge_state::POSITIVE : sidb_charge_state::NONE);
                }

                for (const auto& [c, cs] : cell_maps.first)
                {
                    cds_lb.assign_charge_state(c, cs);
                }

                for (const auto& [c, cs] : cell_maps.second)
                {
                    cds_ub.assign_charge_state(c, cs);
                }

                std::cout << "\n\nLOWER BOUND (" << sidb.x << "," << sidb.y << ")" << std::endl;
                print_layout(cds_lb);
                cds_lb.update_after_charge_change();
                std::cout << "LOC_POT = " << cds_lb.get_local_potential(sidb).value() << std::endl;
                // std::cout << "LOC_POT (stored) = " << skeleton_influence_bounds_map.at(sidb).front()
                std::cout << "LOC_POT (stored) = " << skeleton_influence_bounds_map.at(to_rel_pos.at(sidb)).front()
                          << std::endl;

                std::cout << "\n\nUPPER BOUND (" << sidb.x << "," << sidb.y << ")" << std::endl;
                print_layout(cds_ub);
                cds_ub.update_after_charge_change();
                std::cout << "LOC_POT = " << cds_ub.get_local_potential(sidb).value() << std::endl;
                // std::cout << "LOC_POT (stored) = " << skeleton_influence_bounds_map.at(sidb).back()
                std::cout << "LOC_POT (stored) = " << skeleton_influence_bounds_map.at(to_rel_pos.at(sidb)).back()
                          << std::endl;

                std::cout << std::endl;
            }
        }

        return skeleton_influence_bounds_map;
    }

  private:
    void collect_influence_bounds_for_sidb(const tile<GateLyt>& tile_of_interest, const cell<CellLyt>& relative_sidb,
                                           const cell<CellLyt>& absolute_sidb) noexcept
    {
        if (absolute_sidb == first_sidb || absolute_sidb == second_sidb || absolute_sidb == third_sidb ||
            absolute_sidb == fourth_sidb)
        {
            std::cout << absolute_sidb.x << ',' << absolute_sidb.y << '=' << relative_sidb.x << ',' << relative_sidb.y
                      << std::endl;
            to_rel_pos.insert({absolute_sidb, params.absolute_positions ? absolute_sidb : relative_sidb});
        }

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

                if (std::find(tiles_of_interest.cbegin(), tiles_of_interest.cend(), t) != tiles_of_interest.cend())
                {
                    return;
                }

                const CellLyt cell_lyt_of_other_tile =
                    apply_gate_library<CellLyt, SkeletonGateLibrary>(gate_lyt, std::set{t});

                for (bdl_wire<CellLyt>& wire : detect_bdl_wires(cell_lyt_of_other_tile, params.bdl_wire_params))
                {
                    std::array<double, 2> potential_from_wire{{0, 0}};

                    std::vector<cell<CellLyt>> w1{}, w2{};  //, w11{}, w22{};

                    for (bdl_pair<cell<CellLyt>>& pair : wire.pairs)
                    {
                        // todo replace with pair cell type?
                        if (pair.type == sidb_technology::cell_type::INPUT ||
                            is_contained_in_tile_of_interest(tile_of_interest, pair.lower))
                        {
                            continue;
                        }

                        CellLyt cell_lyt_with_c_and_bound_pair{};

                        cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb,
                                                                        sidb_technology::cell_type::NORMAL);

                        cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);
                        cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);

                        const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                                       params.simulation_parameters};

                        if (absolute_sidb == first_sidb || absolute_sidb == second_sidb ||
                            absolute_sidb == third_sidb || absolute_sidb == fourth_sidb)
                        {
                            w1.push_back(pair.lower);
                            w2.push_back(pair.upper);
                        }

                        // if (tiles_of_interest.size() == 1)
                        // {
                        //     make_pair_sidbs_respect_io_connections(
                        //         absolute_sidb,
                        //         bdl_wires_of_designed_gates_current_tile_with_other_tile.at(node_counter), pair);
                        // }

                        // if (absolute_sidb == first_sidb || absolute_sidb == second_sidb ||
                        //     absolute_sidb == third_sidb || absolute_sidb == fourth_sidb)
                        // {
                        //     w11.push_back(pair.lower);
                        //     w22.push_back(pair.upper);
                        // }

                        potential_from_wire[0] +=
                            cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.lower);

                        if (!is_output_perturber_in_tile_of_interest(tile_of_interest, pair.upper))
                        {
                            potential_from_wire[1] +=
                                cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.upper);
                        }
                    }

                    for (const auto& [cell_maps, sidb] : std::array{
                             std::make_pair(std::make_pair(std::ref(cell_map1_lb), std::ref(cell_map1_ub)), first_sidb),
                             std::make_pair(std::make_pair(std::ref(cell_map2_lb), std::ref(cell_map2_ub)),
                                            second_sidb),
                             std::make_pair(std::make_pair(std::ref(cell_map3_lb), std::ref(cell_map3_ub)), third_sidb),
                             std::make_pair(std::make_pair(std::ref(cell_map4_lb), std::ref(cell_map4_ub)),
                                            fourth_sidb)})
                    {
                        if (absolute_sidb != sidb)
                        {
                            continue;
                        }

                        // bool found  = false;
                        // bool w1_win = false;
                        // for (uint8_t k = 0; k < w11.size(); ++k)
                        // {
                        //     if (!found && w11[k] == w22[k])
                        //     {
                        //         found = true;
                        //
                        //         if (w1[k] == w11[k])
                        //         {
                        //             std::cout << "lowers win" << std::endl;
                        //             w1_win = true;
                        //         }
                        //         else
                        //         {
                        //             std::cout << "uppers win" << std::endl;
                        //         }
                        //     }
                        // }
                        //
                        // if (found)
                        // {
                        //     std::cout << "CODE DOESN'T WORK" << std::endl;
                        //     if (w1_win)
                        //     {
                        //         for (uint8_t k = 0; k < w2.size(); ++k)
                        //         {
                        //             cell_maps.first.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                        //             cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                        //         }
                        //         for (uint8_t k = 0; k < w1.size(); ++k)
                        //         {
                        //             cell_maps.first.insert({w1.at(k), sidb_charge_state::NEGATIVE});
                        //             cell_maps.second.insert({w1.at(k), sidb_charge_state::NEGATIVE});
                        //         }
                        //     }
                        //     else
                        //     {
                        //         for (uint8_t k = 0; k < w2.size(); ++k)
                        //         {
                        //             cell_maps.first.insert({w2.at(k), sidb_charge_state::NEGATIVE});
                        //             cell_maps.second.insert({w2.at(k), sidb_charge_state::NEGATIVE});
                        //         }
                        //         for (uint8_t k = 0; k < w1.size(); ++k)
                        //         {
                        //             cell_maps.first.insert({w1.at(k), sidb_charge_state::NEUTRAL});
                        //             cell_maps.second.insert({w1.at(k), sidb_charge_state::NEUTRAL});
                        //         }
                        //     }
                        // }

                        if (potential_from_wire[0] > potential_from_wire[1])
                        {
                            // 0 orientation is stronger
                            for (uint8_t k = 0; k < w2.size(); ++k)
                            {
                                if (!is_output_perturber_in_tile_of_interest(tile_of_interest, w2.at(k)))
                                {
                                    cell_maps.first.insert({w2.at(k), sidb_charge_state::NEGATIVE});
                                }
                                else
                                {
                                    cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                                }
                                cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                            }
                            for (uint8_t k = 0; k < w1.size(); ++k)
                            {
                                cell_maps.first.insert({w1.at(k), sidb_charge_state::NEUTRAL});
                                cell_maps.second.insert({w1.at(k), sidb_charge_state::NEGATIVE});
                            }
                        }
                        else
                        {
                            // 1 orientation is stronger, LB gets lowers <- DB-, uppers <- DB0
                            //                            UB gets lowers <- DB0, uppers <- DB-
                            for (uint8_t k = 0; k < w1.size(); ++k)
                            {
                                cell_maps.first.insert({w1.at(k), sidb_charge_state::NEGATIVE});
                                cell_maps.second.insert({w1.at(k), sidb_charge_state::NEUTRAL});
                            }
                            for (uint8_t k = 0; k < w2.size(); ++k)
                            {
                                cell_maps.first.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                                if (!is_output_perturber_in_tile_of_interest(tile_of_interest, w2.at(k)))
                                {
                                    cell_maps.second.insert({w2.at(k), sidb_charge_state::NEGATIVE});
                                }
                                else
                                {
                                    cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                                }
                            }
                        }
                    }

                    std::sort(potential_from_wire.begin(), potential_from_wire.end());

                    bounds[0] += potential_from_wire[0];
                    bounds[1] += potential_from_wire[1];
                }

                node_counter++;
            });

        // include influence from the output perturbers if they are part of a different tile
        cell_lyt_with_all_skeletons.foreach_cell(
            [&](const cell<CellLyt>& c)
            {
                if (is_contained_in_tile_of_interest(tile_of_interest, c))
                {
                    return;
                }

                if (cell_lyt_with_all_skeletons.get_cell_type(c) == sidb_technology::OUTPUT_PERTURBER)
                {
                    CellLyt cell_lyt_with_c_and_bound_pair{};

                    cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb, sidb_technology::cell_type::NORMAL);

                    cell_lyt_with_c_and_bound_pair.assign_cell_type(c, sidb_technology::cell_type::NORMAL);

                    const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                                   params.simulation_parameters};

                    bounds[0] += cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, c);
                    bounds[1] += cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, c);

                    for (const auto& [cell_maps, sidb] : std::array{
                             std::make_pair(std::make_pair(std::ref(cell_map1_lb), std::ref(cell_map1_ub)), first_sidb),
                             std::make_pair(std::make_pair(std::ref(cell_map2_lb), std::ref(cell_map2_ub)),
                                            second_sidb),
                             std::make_pair(std::make_pair(std::ref(cell_map3_lb), std::ref(cell_map3_ub)), third_sidb),
                             std::make_pair(std::make_pair(std::ref(cell_map4_lb), std::ref(cell_map4_ub)),
                                            fourth_sidb)})
                    {
                        if (absolute_sidb != sidb)
                        {
                            continue;
                        }

                        cell_maps.first.insert({c, sidb_charge_state::NEGATIVE});
                        cell_maps.second.insert({c, sidb_charge_state::NEGATIVE});
                    }
                }
            });

        // include influence from the input perturbers to the whole layout
        for (bdl_wire<CellLyt>& wire :
             detect_bdl_wires(cell_lyt_with_all_skeletons, params.bdl_wire_params, bdl_wire_selection::INPUT))
        {
            std::array<double, 2> potential_from_wire{{0, 0}};

            std::vector<cell<CellLyt>> w1{}, w2{};  //, w11{}, w22{};

            for (bdl_pair<cell<CellLyt>>& pair : wire.pairs)
            {
                if (is_contained_in_tile_of_interest(tile_of_interest, pair.lower))
                {
                    continue;
                }

                CellLyt cell_lyt_with_c_and_bound_pair{};

                cell_lyt_with_c_and_bound_pair.assign_cell_type(absolute_sidb, sidb_technology::cell_type::NORMAL);

                cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.upper, sidb_technology::cell_type::NORMAL);
                cell_lyt_with_c_and_bound_pair.assign_cell_type(pair.lower, sidb_technology::cell_type::NORMAL);

                const charge_distribution_surface<CellLyt> cds_with_c_and_pair{cell_lyt_with_c_and_bound_pair,
                                                                               params.simulation_parameters};

                if (absolute_sidb == first_sidb || absolute_sidb == second_sidb || absolute_sidb == third_sidb ||
                    absolute_sidb == fourth_sidb)
                {
                    w1.push_back(pair.lower);
                    w2.push_back(pair.upper);
                }

                potential_from_wire[0] +=
                    cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.lower);
                potential_from_wire[1] +=
                    cds_with_c_and_pair.get_chargeless_potential_between_sidbs(absolute_sidb, pair.upper);
            }

            for (const auto& [cell_maps, sidb] : std::array{
                     std::make_pair(std::make_pair(std::ref(cell_map1_lb), std::ref(cell_map1_ub)), first_sidb),
                     std::make_pair(std::make_pair(std::ref(cell_map2_lb), std::ref(cell_map2_ub)), second_sidb),
                     std::make_pair(std::make_pair(std::ref(cell_map3_lb), std::ref(cell_map3_ub)), third_sidb),
                     std::make_pair(std::make_pair(std::ref(cell_map4_lb), std::ref(cell_map4_ub)), fourth_sidb)})
            {
                if (absolute_sidb != sidb)
                {
                    continue;
                }

                if (potential_from_wire[0] > potential_from_wire[1])
                {
                    // 0 orientation is stronger
                    for (uint8_t k = 0; k < w2.size(); ++k)
                    {
                        if (!is_output_perturber_in_tile_of_interest(tile_of_interest, w2.at(k)))
                        {
                            cell_maps.first.insert({w2.at(k), sidb_charge_state::NEGATIVE});
                        }
                        else
                        {
                            cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                        }
                        cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                    }
                    for (uint8_t k = 0; k < w1.size(); ++k)
                    {
                        cell_maps.first.insert({w1.at(k), sidb_charge_state::NEUTRAL});
                        cell_maps.second.insert({w1.at(k), sidb_charge_state::NEGATIVE});
                    }
                }
                else
                {
                    // 1 orientation is stronger, LB gets lowers <- DB-, uppers <- DB0
                    //                            UB gets lowers <- DB0, uppers <- DB-
                    for (uint8_t k = 0; k < w1.size(); ++k)
                    {
                        cell_maps.first.insert({w1.at(k), sidb_charge_state::NEGATIVE});
                        cell_maps.second.insert({w1.at(k), sidb_charge_state::NEUTRAL});
                    }
                    for (uint8_t k = 0; k < w2.size(); ++k)
                    {
                        cell_maps.first.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                        if (!is_output_perturber_in_tile_of_interest(tile_of_interest, w2.at(k)))
                        {
                            cell_maps.second.insert({w2.at(k), sidb_charge_state::NEGATIVE});
                        }
                        else
                        {
                            cell_maps.second.insert({w2.at(k), sidb_charge_state::NEUTRAL});
                        }
                    }
                }
            }

            std::sort(potential_from_wire.begin(), potential_from_wire.end());

            bounds[0] += potential_from_wire[0];
            bounds[1] += potential_from_wire[1];
        }

        if (skeleton_influence_bounds_map.count(params.absolute_positions ? absolute_sidb : relative_sidb) != 0)
        {
            std::cout << "HELP IT'S OVERWRITING" << std::endl;
        }

        std::cout << fmt::format("storing ({}, {}) at {},{}", bounds.front(), bounds.back(),
                                 (params.absolute_positions ? absolute_sidb : relative_sidb).x,
                                 (params.absolute_positions ? absolute_sidb : relative_sidb).y)
                  << std::endl;

        skeleton_influence_bounds_map.insert(
            {params.absolute_positions ? absolute_sidb : relative_sidb, std::move(bounds)});  // TODO: negate???
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

    bool is_contained_in_tile_of_interest(const tile<GateLyt>& tile_of_interest,
                                          const cell<CellLyt>& absolute_sidb) noexcept
    {
        const CellLyt cell_lyt_of_tile_of_interest =
            apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{tile_of_interest});

        bool found = false;

        cell_lyt_of_tile_of_interest.foreach_cell(
            [&](const cell<CellLyt>& sidb_in_tile_of_interest)
            {
                if (!found && absolute_sidb == sidb_in_tile_of_interest &&
                    cell_lyt_of_tile_of_interest.get_cell_type(sidb_in_tile_of_interest) !=
                        sidb_technology::cell_type::EMPTY)
                {
                    found = true;
                }
            });

        return found;
    }

    bool is_output_perturber_in_tile_of_interest(const tile<GateLyt>& tile_of_interest,
                                                 const cell<CellLyt>& absolute_sidb) noexcept
    {
        const CellLyt cell_lyt_of_tile_of_interest =
            apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{tile_of_interest});

        bool found = false;

        cell_lyt_of_tile_of_interest.foreach_cell(
            [&](const cell<CellLyt>& sidb_in_tile_of_interest)
            {
                if (!found && absolute_sidb == sidb_in_tile_of_interest &&
                    cell_lyt_of_tile_of_interest.get_cell_type(sidb_in_tile_of_interest) ==
                        sidb_technology::cell_type::OUTPUT_PERTURBER)
                {
                    found = true;
                }
            });

        return found;
    }

    const skeleton_influence_bounds_params<cell<CellLyt>>&   params;
    GateLyt                                                  gate_lyt{};
    const std::vector<tile<GateLyt>>&                        tiles_of_interest{};
    // const std::vector<std::vector<bdl_wire<CellLyt>>> bdl_wires_of_designed_gates_current_tile_with_other_tile{};
    const CellLyt                                            cell_lyt_with_all_skeletons{};
    std::unordered_map<cell<CellLyt>, std::array<double, 2>> skeleton_influence_bounds_map{};

    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map1_lb{};
    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map1_ub{};

    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map2_lb{};
    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map2_ub{};

    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map3_lb{};
    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map3_ub{};

    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map4_lb{};
    std::unordered_map<cell<CellLyt>, sidb_charge_state> cell_map4_ub{};
    cell<CellLyt>                                        first_sidb{};
    cell<CellLyt>                                        second_sidb{};

    cell<CellLyt>                                    third_sidb{};
    cell<CellLyt>                                    fourth_sidb{};
    std::unordered_map<cell<CellLyt>, cell<CellLyt>> to_rel_pos{};

    // [[nodiscard]] static std::vector<std::vector<bdl_wire<CellLyt>>>
    // obtain_bdl_wires_for_all_tile_pairs(const GateLyt& gate_lyt, const tile<GateLyt>& current_tile) noexcept
    // {
    //     std::vector<std::vector<bdl_wire<CellLyt>>> bdl_wires_of_designed_gates_current_tile_with_other_tile{};
    //
    //     gate_lyt.foreach_node(
    //         [&](const auto& n)
    //         {
    //             if (gate_lyt.is_constant(n))
    //             {
    //                 return;
    //             }
    //
    //             const tile<GateLyt>& t = gate_lyt.get_tile(n);
    //
    //             if (t == current_tile)
    //             {
    //                 return;
    //             }
    //
    //             CellLyt designed_gates_together =
    //                 apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_lyt, std::set{current_tile, t});
    //
    //             std::vector<bdl_wire<CellLyt>> bdl_wires_of_designed_gates_together =
    //                 detect_bdl_wires(designed_gates_together, params.bdl_wire_params);
    //
    //             // sort by wire length, first will be biggest
    //             std::sort(bdl_wires_of_designed_gates_together.begin(), bdl_wires_of_designed_gates_together.end(),
    //                       [](const bdl_wire<CellLyt>& lhs, const bdl_wire<CellLyt>& rhs)
    //                       { return lhs.pairs.size() > rhs.pairs.size(); });
    //
    //             bdl_wires_of_designed_gates_current_tile_with_other_tile.push_back(
    //                 std::move(bdl_wires_of_designed_gates_together));
    //         });
    //
    //     return bdl_wires_of_designed_gates_current_tile_with_other_tile;
    // }

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
skeleton_influence_bounds(const GateLyt& gate_layout, const std::vector<tile<GateLyt>>& tile,
                          const skeleton_influence_bounds_params<cell<CellLyt>>& params = {}) noexcept
{
    return detail::skeleton_influence_bounds_impl<CellLyt, SkeletonGateLibrary, GateLyt>{gate_layout, tile, params}
        .run();
}

}  // namespace fiction

#endif  // SKELETON_INFLUENCE_BOUNDS_HPP

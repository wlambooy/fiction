//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/technology/sidb_defect_surface.hpp"
#include "fiction/technology/sidb_on_the_fly_gate_library.hpp"
#include "fiction/technology/sidb_skeleton_bestagon_library.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/gate_design_utils.hpp"

#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/views/window_view.hpp>

#include <cstdint>
#include <optional>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiction
{

/**
 * This struct stores the parameters to design an SiDB circuit on a defective surface.
 *
 * @tparam CellLyt SiDB cell-level layout type.
 */
template <typename CellLyt>
struct advanced_circuit_design_params
{
    /**
     * Parameters for the SiDB on-the-fly gate library.
     */
    sidb_on_the_fly_gate_library_params<CellLyt> sidb_on_the_fly_gate_library_parameters = {};
    /**
     * Parameters for the *exact* placement and routing algorithm.
     */
    exact_physical_design_params exact_design_parameters  = {};
    uint64_t                     num_trials               = 500;
    double                       minimum_passing_fraction = 0.3;
};

/**
 * Statistics for the on-the-fly defect-aware circuit design.
 */
template <typename GateLyt>
struct advanced_circuit_design_stats
{
    /**
     * The total runtime of the operational domain computation.
     */
    mockturtle::stopwatch<>::duration time_total{0};
    /**
     * The `stats` of the *exact* algorithm.
     */
    exact_physical_design_stats exact_stats{};
    /**
     * The gate-level layout after P&R.
     */
    std::optional<GateLyt> gate_layout{};
};

namespace detail
{

template <typename Ntk, typename CellLyt, typename GateLyt, typename GateLibrary, typename SkeletonGateLibrary>
class advanced_circuit_design_impl
{
  public:
    advanced_circuit_design_impl(const Ntk& ntk, advanced_circuit_design_params<CellLyt>& design_params,
                                 const GateLyt& tiling, advanced_circuit_design_stats<GateLyt>& st) :
            lattice_tiling{tiling},
            network{ntk},
            params{design_params},
            stats{st}
    {}

    [[nodiscard]] sidb_defect_surface<CellLyt> design_circuit_on_defective_surface()
    {
        const mockturtle::stopwatch stop{stats.time_total};

        std::optional<GateLyt> gate_lyt = std::nullopt;

        CellLyt lyt{};

        std::unordered_map<mockturtle::node<Ntk>, std::vector<typename GateLibrary::fcn_gate>>
            operational_gate_designs{};

        // generating the blacklist based on neutral defects. The long-range electrostatic influence of charged defects
        // is not considered as gates are designed on-the-fly.
        auto black_list = sidb_surface_analysis<SkeletonGateLibrary>(
            lattice_tiling, params.sidb_on_the_fly_gate_library_parameters.defect_surface, std::make_pair(0, 0));

        while (!gate_lyt.has_value())
        {
            // P&R with *exact* and the pre-determined blacklist
            gate_lyt =
                exact_with_blacklist<GateLyt>(network, black_list, params.exact_design_parameters, &stats.exact_stats);

            if (!gate_lyt.has_value())
            {
                // P&R was unsuccessful
                break;
            }

            operational_gate_designs.clear();

            try
            {
                gate_lyt->foreach_node(
                    [&, this](const auto& n, [[maybe_unused]] auto i)
                    {
                        if (!gate_lyt->is_constant(n))
                        {
                            const auto t = gate_lyt->get_tile(n);

                            operational_gate_designs[n] =
                                GateLibrary::template set_up_gates<GateLyt, CellLyt,
                                                                   sidb_on_the_fly_gate_library_params<CellLyt>>(
                                    *gate_lyt, t, params.sidb_on_the_fly_gate_library_parameters);

                            // std::cout << "OPERATIONAL GATE IMPLEMENTATIONS FOR TILE " << gate_lyt->get_tile(n)
                            //           << std::endl;
                            //
                            // for (const auto& gate : operational_gate_designs.at(n))
                            // {
                            //     CellLyt lytt{};
                            //
                            //     assign_gate<CellLyt, GateLibrary, GateLyt>(
                            //         lytt, cell<CellLyt>{0, 0},
                            //         // relative_to_absolute_cell_position<GateLibrary::gate_x_size(),
                            //         // GateLibrary::gate_y_size(),
                            //         //                                    GateLyt, CellLyt>(*gate_lyt,
                            //         //                                    gate_lyt->get_tile(n),
                            //         //                                                      cell<CellLyt>{0, 0}),
                            //         gate, *gate_lyt, n);
                            //
                            //     print_layout(lytt);
                            //     std::cout << std::endl << std::endl;
                            // }
                        }
                    });
            }

            // on-the-fly gate design was unsuccessful at a certain tile. Hence, this tile-gate pair is added to the
            // blacklist and the process is rerun.
            catch (const gate_design_exception<tt, GateLyt>& e)
            {
                black_list[e.which_tile()][e.which_truth_table()].push_back(e.which_port_list());
            }
        }

        if (prune_gate_designs(*gate_lyt, operational_gate_designs))
        {
            // assign gates
            gate_lyt->foreach_node(
                [&](const auto& n)
                {
                    if (!gate_lyt->is_constant(n))
                    {
                        // select a random gate implementation for the tile that connects as input to n
                        assign_gate<CellLyt, GateLibrary, GateLyt>(
                            lyt,
                            relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(),
                                                               GateLyt, CellLyt>(*gate_lyt, gate_lyt->get_tile(n),
                                                                                 cell<CellLyt>{0, 0}),
                            operational_gate_designs.at(n).front(), *gate_lyt, n);
                    }
                });

            std::cout << "FINAL GENERATED CIRCUIT:" << std::endl;
            print_layout(lyt);
        }
        else
        {
            std::cout << "FAILURE: NO OPERATIONAL CIRCUIT COULD BE GENERATED" << std::endl;
        }

        stats.gate_layout = std::optional{gate_lyt};

        sidb_defect_surface<CellLyt> sidbs_and_defects{lyt};

        // add defects to the circuit.
        params.sidb_on_the_fly_gate_library_parameters.defect_surface.foreach_sidb_defect(
            [&sidbs_and_defects](const auto& defect)
            { sidbs_and_defects.assign_sidb_defect(defect.first, defect.second); });

        return sidbs_and_defects;
    }

  private:
    /**
     * Gate-level layout.
     */
    GateLyt lattice_tiling;
    /**
     * Network.
     */
    Ntk network;
    /**
     * Parameters for the on-the-fly circuit design.
     */
    advanced_circuit_design_params<CellLyt> params{};
    /**
     * Statistics for the on-the-fly circuit design.
     */
    advanced_circuit_design_stats<GateLyt>& stats;
    /**
     *
     */
    bool prune_gate_designs(const GateLyt& gate_lyt,
                            std::unordered_map<mockturtle::node<Ntk>, std::vector<typename GateLibrary::fcn_gate>>&
                                operational_gate_designs) const noexcept
    {

        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS" << std::endl;

        std::unordered_map<mockturtle::node<Ntk>, cell<CellLyt>>                      top_left_corner_abs{};
        std::unordered_map<mockturtle::node<Ntk>, std::vector<mockturtle::node<Ntk>>> inputs_to_node{};
        std::unordered_map<mockturtle::node<Ntk>, uint64_t>                           number_of_variables_to_node{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (gate_lyt.is_constant(n))
                {
                    return;
                }

                // std::cout << "CURRENT NODE: " << n << std::endl;
                // std::cout << "CURRENT TILE: " << gate_lyt.get_tile(n) << std::endl;

                top_left_corner_abs.insert(
                    {n,
                     relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                                        CellLyt>(gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0})});

                std::vector<mockturtle::node<Ntk>> inputs_to_n{};

                const auto& incoming_tiles = gate_lyt.is_pi(n) ? gate_lyt.outgoing_data_flow(gate_lyt.get_tile(n)) :
                                                                 gate_lyt.incoming_data_flow(gate_lyt.get_tile(n));

                for (const auto& t : incoming_tiles)
                {
                    inputs_to_n.emplace_back(gate_lyt.get_node(t));
                    // std::cout << " OTHER NODE: " << inputs_to_n.back() << std::endl;
                    // std::cout << " OTHER TILE: " << t << std::endl;
                }
                //
                // gate_lyt.foreach_fanin(n,
                //                        [&](const mockturtle::node<Ntk>& other_n)
                //                        {
                //                            if (!gate_lyt.is_dead(other_n))
                //                            {
                //                                std::cout << " OTHER NODE: " << other_n << std::endl;
                //                                 std::cout << " OTHER TILE: " << gate_lyt.get_tile(other_n) <<
                //                                 std::endl;
                //                                inputs_to_n.emplace_back(other_n);
                //                            }
                //                            // else
                //                            // {
                //                            //     std::cout << " other NODE: " << other_n << std::endl;
                //                            //      std::cout << " other TILE: " << gate_lyt.get_tile(other_n) <<
                //                            std::endl;
                //                            // }
                //                        });

                // if (gate_lyt.is_pi(n))
                // {
                //     number_of_variables_to_node.insert({n, 1});
                // }
                // else
                // {
                //     number_of_variables_to_node.insert({n, inputs_to_n.size()});
                // }

                inputs_to_node.insert({n, std::move(inputs_to_n)});
            });

        bool big_fixpoint = false;

        bool exit_by_failure = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            std::cout << "Starting main fixpoint iteration" << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || gate_lyt.is_constant(n))
                    {
                        return;
                    }

                    std::cout << fmt::format("Starting pruning for tile {}\n", gate_lyt.get_tile(n));

                    // std::vector<mockturtle::node<Ntk>> tiles_selected_for_operational_assessment{n.children};
                    // std::vector<mockturtle::signal<Ntk>> output_signals{};
                    // gate_lyt.foreach_fanout(n,
                    //                         [&](const auto& output_n) {
                    //                             output_signals.emplace_back(
                    //                                 static_cast<mockturtle::signal<Ntk>>(gate_lyt.get_tile(output_n)));
                    //                         });

                    // std::vector<mockturtle::node<Ntk>> inputs_to_n{};

                    // const std::vector<tile<GateLyt>>& outgoing_tiles =
                    //     gate_lyt.outgoing_data_flow(gate_lyt.get_tile(n));
                    // std::vector<mockturtle::signal<Ntk>> output_signals{};
                    //
                    // for (const auto& t : outgoing_tiles)
                    // {
                    //     output_signals.emplace_back(static_cast<mockturtle::signal<Ntk>>(t));
                    //     // std::cout << " OTHER NODE: " << inputs_to_n.back() << std::endl;
                    //     // std::cout << " OTHER TILE: " << t << std::endl;
                    // }

                    // uint32_t total_number_of_inputs = 0;
                    //
                    // for (const auto& other_n : inputs_to_node.at(n))
                    // {
                    //     total_number_of_inputs += number_of_variables_to_node.at(other_n);
                    // }

                    //     std::accumulate(
                    //                             inputs_to_node.at(n).cbegin(), inputs_to_node.at(n).cend(),
                    //                             uint32_t{0},
                    //                             [&](const uint32_t               sum,
                    //                                                            const mockturtle::node<Ntk>&
                    //                                                            sub_child)
                    //                             {
                    //                                 std::cout << " OTHER NODE: " << sub_child << std::endl;
                    //     std::cout << " OTHER TILE: " << gate_lyt.get_tile(sub_child) << std::endl;
                    //                                 for (const auto& [k,v] : number_of_variables_to_node)
                    //                                 {
                    //                                 std::cout << " key: " << k << "\tvalue: " << v << std::endl;
                    //                                 }
                    //                             std::cout << "now comes: " <<
                    //                             number_of_variables_to_node.at(sub_child) << std::endl;
                    // return sum + static_cast<uint32_t>(number_of_variables_to_node.at(sub_child));
                    //                         })}

                    // gate_lyt.foreach_fanin(n,
                    //                        [&](const auto& before_n)
                    //                        {
                    //                            if (!gate_lyt->is_constant(before_n))  // necessary?
                    //                            {
                    //                                tiles_selected_for_operational_assessment.emplace_back(before_n);
                    //                            }
                    //                        });

                    uint64_t original_size = operational_gate_designs.at(n).size();

                    uint64_t current_tile_gate_implementation_index = 0;

                    while (current_tile_gate_implementation_index < operational_gate_designs.at(n).size())
                    {
                        // uint64_t current_iteration_gate = 0;

                        // std::vector<mockturtle::node<Ntk>> other_n_order{};
                        // other_n_order.reserve(tiles_selected_for_operational_assessment.size());
                        //
                        // std::vector<uint64_t> gate_design_indices{};
                        // gate_design_indices.reserve(tiles_selected_for_operational_assessment.size());
                        //
                        // for (const tile<GateLyt>& other_n : tiles_selected_for_operational_assessment)
                        // {
                        //     gate_design_indices.push_back(0);
                        //     other_n_order.push_back(other_n);
                        // }

                        CellLyt cell_lyt{};

                        // select the first gate implementation for n
                        assign_gate<CellLyt, GateLibrary, GateLyt>(
                            cell_lyt, top_left_corner_abs.at(n),
                            *std::next(operational_gate_designs.at(n).cbegin(),
                                       static_cast<int64_t>(current_tile_gate_implementation_index)),
                            gate_lyt, n);

                        uint64_t current_trial = 0, successful_trials = 0;

                        // do `num_trials` random trials of looking for gate implementations of the other gates
                        // connecting as inputs to the current; when less than `minimum_passing_fraction` * `num_trials`
                        // number of trials are found to be operational, the current chosen gate for n gets discarded,
                        // and we move to the next
                        while (current_trial < params.num_trials)
                        {
                            current_trial++;

                            CellLyt cell_lyt_clone = cell_lyt.clone();

                            for (const auto& other_n : inputs_to_node.at(n))
                            {
                                std::random_device rd;         // a seed source for the random number engine
                                std::mt19937       gen(rd());  // mersenne_twister_engine seeded with rd()
                                std::uniform_int_distribution<uint64_t> distrib{
                                    0, operational_gate_designs.at(other_n).size() - 1};

                                // select a random gate implementation for the tile that connects as input to n
                                assign_gate<CellLyt, GateLibrary, GateLyt>(
                                    cell_lyt_clone, top_left_corner_abs.at(other_n),
                                    operational_gate_designs.at(other_n).at(distrib(gen)), gate_lyt, n);
                            }

                            // std::cout << "IMPLEMENTATION COMBINATION TO ASSESS:" << std::endl;
                            // print_layout(cell_lyt_clone);

                            // const auto& [reduced_ntk, num_pis] =
                            //     make_network_around_nodes(n, inputs_to_node.at(n));

                            //
                            // std::vector<mockturtle::node<Ntk>> copy_vec{};
                            // std::copy(inputs_to_node.at(n).cbegin(), inputs_to_node.at(n).cend(),
                            //           std::back_inserter(copy_vec));
                            // // copy_vec.push_back(n);
                            // // for (const mockturtle::signal<Ntk>& s : output_signals)
                            // // {
                            // //     copy_vec.emplace_back(gate_lyt.get_node(static_cast<tile<GateLyt>>(s)));
                            // // }
                            //
                            // if (params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                            //         .cc_map.has_value())
                            // {
                            //     std::cout << "YES IT HAS A VALUE | size = "
                            //               << params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                            //                      .operational_params.cc_map.value()
                            //                      .size()
                            //               << std::endl;
                            // }
                            //
                            // std::cout << "TOTAL NUMBER OF INPUTS: " << total_number_of_inputs << std::endl;

                            // mockturtle::window_view<Ntk> window{
                            //                           reduced_ntk,
                            //                           reduced_ntk_pis,
                            //                           {n},
                            //                           {n}};
                            //
                            // mockturtle::default_simulator<kitty::dynamic_truth_table> sim{total_number_of_inputs};
                            //
                            // Ntk ntk;
                            // for (const auto& other_n : inputs_to_node.at(n))
                            // {
                            //     std::vector<mockturtle::node<Ntk>> pis{};
                            //
                            //
                            //
                            //     gate_lyt.foreach_fanin(other_n, [&](const auto& input_to_other_n) { if
                            //     (!gate_lyt.is_constant(input_to_other_n)) {pis.emplace_back(ntk.create_pi());}});
                            // }
                            //
                            //
                            // switch_node_function(ntk, )

                            if (is_operational(
                                    cell_lyt_clone, create_fan_out_tt(),
                                    // mockturtle::simulate<kitty::dynamic_truth_table>(reduced_ntk, num_pis),
                                    params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                        .operational_params)
                                    .status == operational_status::OPERATIONAL)
                            {
                                // std::cout << (gate_lyt.is_fanout(n) ? "FANOUT" : "BUFFER") << " IS OPERATIONAL!!!!!!"
                                // << std::endl;
                                successful_trials++;
                            }
                            else
                            {
                                // std::cout << (gate_lyt.is_fanout(n) ? "FANOUT" : "BUFFER") << " IS NON-OPERATIONAL"
                                // << std::endl;
                            }
                        }

                        if (const double passing_fraction =
                                static_cast<double>(successful_trials) / static_cast<double>(params.num_trials);
                            passing_fraction < params.minimum_passing_fraction)
                        {
                            operational_gate_designs[n].erase(
                                std::next(operational_gate_designs.at(n).begin(),
                                          static_cast<int64_t>(current_tile_gate_implementation_index)));

                            std::cout
                                << fmt::format(
                                       "pruned gate design for ({},{}) since {} out of {} trials were successful "
                                       "(percentage = {:.1f} < {:.1f})  |  current gate index = {}, remaining = {}",
                                       gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y, successful_trials,
                                       params.num_trials, passing_fraction * 100, params.minimum_passing_fraction * 100,
                                       current_tile_gate_implementation_index,
                                       operational_gate_designs.at(n).size() - current_tile_gate_implementation_index)
                                << std::endl;
                        }
                        else
                        {
                            current_tile_gate_implementation_index++;

                            std::cout
                                << fmt::format(
                                       "*KEPT* gate design for ({},{}) since {} out of {} trials were successful "
                                       "(percentage = {:.1f} ≥ {:.1f})  |  current gate index = {}, remaining = {}",
                                       gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y, successful_trials,
                                       params.num_trials, passing_fraction * 100, params.minimum_passing_fraction * 100,
                                       current_tile_gate_implementation_index,
                                       operational_gate_designs.at(n).size() - current_tile_gate_implementation_index)
                                << std::endl;
                        }

                        // while (fixpoint && current_iteration_gate < other_n_order.size())
                        // {
                        //     if (gate_design_indices.at(current_iteration_gate) ==
                        //     operational_gate_designs.at(other_n_order.at(current_iteration_gate)))
                        //     {
                        //
                        //     }
                        //
                        //     CellLyt cell_lyt{};
                        //
                        //     assign_gate<CellLyt, sidb_on_the_fly_gate_library, GateLyt>(
                        //         cell_lyt, c, operational_gate_designs.at(t).at(counter++), gate_lyt, n);
                        //
                        //     // if (gate_design_indices.at(current_iteration_gate) ==
                        //     other_n_order[current_iteration_gate]())
                        //     // if ()
                        //     // {
                        //     //     fixpoint = false;
                        //     //     continue;
                        //     // }
                        //
                        //
                        // }
                    }

                    if (operational_gate_designs.at(n).empty())
                    {
                        std::cout << "ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (operational_gate_designs.at(n).size() < original_size)
                    {
                        big_fixpoint = false;
                        std::cout << "SET FIXPOINT TO FALSE" << std::endl;
                    }

                    // assign_gate<CellLyt, sidb_on_the_fly_gate_library, GateLyt>(cell_lyt, c,
                    // operational_gate_designs.at(t).at(counter++), gate_lyt, n);
                });
        }

        return !exit_by_failure;
    }

    // [[nodiscard]] std::vector<kitty::dynamic_truth_table>
    // determine_truth_table_of_extended_window(const mockturtle::node<Ntk>&              n,
    //                                          const std::vector<mockturtle::node<Ntk>>& inputs) noexcept
    // {
    //     const mockturtle::window_view<Ntk> singleton_window{network, inputs, n.children, {n}};
    //
    //     return mockturtle::simulate<kitty::dynamic_truth_table, Ntk>(singleton_window);
    // }

    // std::pair<Ntk, uint32_t> make_network_around_nodes(const GateLyt& gate_lyt, const mockturtle::node<Ntk>& node,
    // const std::vector<mockturtle::node<Ntk>>& fan_in_nodes) const noexcept
    // {
    //     Ntk reduced_ntk;
    //
    //     for (const auto& fi_n : fan_in_nodes)
    //     {
    //         for (const tile<GateLyt>& t : gate_lyt.incoming_data_flow(gate_lyt.get_tile(fi_n)))
    //         {
    //             const mockturtle::node<Ntk>& fi_fi_n = gate_lyt.get_node(t);
    //
    //             if (gate_lyt.is_pi(fi_fi_n))
    //             {
    //                 reduced_ntk.create_pi();
    //             }
    //         }
    //         // gate_lyt.foreach_fanin(fi_n, [&](const auto& fi_fi_n)
    //         // {
    //         //
    //         // });
    //     }
    // }
};

}  // namespace detail

/**
 *
 *
 * @tparam Ntk The type of the input network.
 * @tparam CellLyt SiDB cell-level layout type.
 * @tparam GateLyt Gate-level layout type.
 * todo
 * @param ntk The input network to be mapped onto the defective surface.
 * @param lattice_tiling The lattice tiling used for the circuit design.
 * @param params The parameters used for designing the circuit, encapsulated in an
 * `advanced_circuit_design_params` object.
 * @param stats Pointer to a structure for collecting statistics. If nullptr, statistics are not collected.
 * @return A `sidb_defect_surface<CellLyt>` representing the designed circuit on the defective surface.
 */
template <typename Ntk, typename CellLyt, typename GateLyt, typename GateLibrary, typename SkeletonGateLibrary>
[[nodiscard]] sidb_defect_surface<CellLyt>
advanced_circuit_design(const Ntk& ntk, const GateLyt& lattice_tiling,
                        advanced_circuit_design_params<CellLyt>& params = {},
                        advanced_circuit_design_stats<GateLyt>*  stats  = nullptr)
{
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt is not a hexagonal");
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");

    advanced_circuit_design_stats<GateLyt> st{};

    detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt, GateLibrary, SkeletonGateLibrary> p{ntk, params,
                                                                                                    lattice_tiling, st};

    const auto result = p.design_circuit_on_defective_surface();

    if (stats)
    {
        *stats = st;
    }

    return result;
}

}  // namespace fiction

#endif  // advanced_CIRCUIT_DESIGN_HPP

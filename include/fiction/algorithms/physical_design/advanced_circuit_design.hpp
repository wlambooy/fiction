//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/skeleton_influence_bounds.hpp"
#include "fiction/technology/sidb_defect_surface.hpp"
#include "fiction/technology/sidb_on_the_fly_gate_library.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/gate_design_utils.hpp"

#if (PROGRESS_BARS)
#include <mockturtle/utils/progress_bar.hpp>
#endif

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
    exact_physical_design_params exact_design_parameters = {};
    uint64_t                     num_trials              = 500;
    double                       selectivity             = 0.5;
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
        std::unordered_map<mockturtle::node<Ntk>, std::vector<mockturtle::node<Ntk>>> gate_connections_to_simulate{};
        std::unordered_map<mockturtle::node<Ntk>, is_operational_params<cell<CellLyt>>>
                                                                         operational_params_for_joint_simulation{};
        std::unordered_map<mockturtle::node<Ntk>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (gate_lyt.is_constant(n))
                {
                    return;
                }

                top_left_corner_abs.insert(
                    {n,
                     relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                                        CellLyt>(gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0})});

                std::vector<mockturtle::node<Ntk>> inputs_to_n{};
                std::vector<tile<GateLyt>>         tiles_to_simulate_together{{gate_lyt.get_tile(n)}};

                const auto& incoming_tiles = gate_lyt.is_pi(n) ? gate_lyt.outgoing_data_flow(gate_lyt.get_tile(n)) :
                                                                 gate_lyt.incoming_data_flow(gate_lyt.get_tile(n));

                for (const auto& t : incoming_tiles)
                {
                    inputs_to_n.emplace_back(gate_lyt.get_node(t));
                    tiles_to_simulate_together.emplace_back(t);
                }

                is_operational_params operational_params =
                    params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params;

                operational_params.cc_map = skeleton_influence_bounds<CellLyt, SkeletonGateLibrary, GateLyt>(
                    gate_lyt, tiles_to_simulate_together,
                    skeleton_influence_bounds_params<cell<CellLyt>>{
                        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                            .simulation_parameters,
                        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.canvas,
                        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                            .input_bdl_iterator_params.bdl_wire_params,
                        true});

                gate_connections_to_simulate.insert({n, std::move(inputs_to_n)});
                operational_params_for_joint_simulation.insert({n, std::move(operational_params)});
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

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{
                        static_cast<uint32_t>(operational_gate_designs.at(n).size()),
                        "[i] Determining successful trial ratio for tile " +
                            fmt::format("({},{})", gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y) + ": |{0}|"};
#endif

                    for (uint64_t current_tile_gate_implementation_index = 0;
                         current_tile_gate_implementation_index < operational_gate_designs.at(n).size();
                         ++current_tile_gate_implementation_index)
                    {
                        CellLyt cell_lyt{};

                        // select the first gate implementation for n
                        assign_gate<CellLyt, GateLibrary, GateLyt>(
                            cell_lyt, top_left_corner_abs.at(n),
                            *std::next(operational_gate_designs.at(n).cbegin(),
                                       static_cast<int64_t>(current_tile_gate_implementation_index)),
                            gate_lyt, n);

                        uint64_t current_trial = 0, successful_trials = 0;

                        while (current_trial < params.num_trials)
                        {
                            current_trial++;

                            CellLyt cell_lyt_clone = cell_lyt.clone();

                            for (const auto& other_n : gate_connections_to_simulate.at(n))
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

                            if (is_operational(
                                    cell_lyt_clone, create_fan_out_tt(),
                                    // mockturtle::simulate<kitty::dynamic_truth_table>(reduced_ntk, num_pis),
                                    operational_params_for_joint_simulation.at(n))
                                    .status == operational_status::OPERATIONAL)
                            {
                                successful_trials++;
                            }
                        }

                        const double successful_trial_ratio =
                            static_cast<double>(successful_trials) / static_cast<double>(params.num_trials);

                        successful_trial_ratio_per_gate_implementation.emplace_back(
                            successful_trial_ratio, current_tile_gate_implementation_index);
#if (PROGRESS_BARS)
                        // update progress
                        bar(current_tile_gate_implementation_index);
#endif
                    }

                    std::sort(successful_trial_ratio_per_gate_implementation.begin(),
                              successful_trial_ratio_per_gate_implementation.end(),
                              [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

                    const auto first_passing_gate_ix = static_cast<uint64_t>(
                        params.selectivity *
                        static_cast<double>(successful_trial_ratio_per_gate_implementation.size()));

                    std::cout << fmt::format(
                                     "\nDetermined passing ratio: {:.1f}%",
                                     successful_trial_ratio_per_gate_implementation.at(first_passing_gate_ix).first *
                                         100)
                              << std::endl;

                    std::cout << fmt::format("Pruning {} out of {} gate designs", first_passing_gate_ix,
                                             successful_trial_ratio_per_gate_implementation.size())
                              << std::endl;

                    selected_gate_implementation_indices[n].clear();

                    for (uint64_t current_tile_gate_implementation_index = 0;
                         current_tile_gate_implementation_index < operational_gate_designs.at(n).size();
                         ++current_tile_gate_implementation_index)
                    {
                        if (current_tile_gate_implementation_index < first_passing_gate_ix)
                        {
                            std::cout
                                << fmt::format(
                                       "pruned gate design for ({},{}) since {:.1f}% of the {} trials were successful",
                                       gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y,
                                       100 * successful_trial_ratio_per_gate_implementation
                                                 .at(current_tile_gate_implementation_index)
                                                 .first,
                                       params.num_trials)
                                << std::endl;
                        }
                        else
                        {
                            std::cout
                                << fmt::format(
                                       "*KEPT* gate design for ({},{}) since {:.1f}% of the {} trials were successful",
                                       gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y,
                                       100 * successful_trial_ratio_per_gate_implementation
                                                 .at(current_tile_gate_implementation_index)
                                                 .first,
                                       params.num_trials)
                                << std::endl;
                            selected_gate_implementation_indices[n].emplace_back(
                                successful_trial_ratio_per_gate_implementation
                                    .at(current_tile_gate_implementation_index)
                                    .second);
                        }
                    }

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() < operational_gate_designs.at(n).size())
                    {
                        big_fixpoint = false;
                    }
                });

            std::cout << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || gate_lyt.is_constant(n))
                    {
                        return;
                    }

                    std::vector<typename GateLibrary::fcn_gate> selected_gate_implementations{};

                    for (const uint64_t selected_gate_implementation_index : selected_gate_implementation_indices.at(n))
                    {
                        selected_gate_implementations.push_back(
                            std::move(operational_gate_designs[n][selected_gate_implementation_index]));
                    }

                    operational_gate_designs[n] = std::move(selected_gate_implementations);
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

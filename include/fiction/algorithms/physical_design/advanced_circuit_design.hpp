//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_circuit_operational.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/technology/sidb_defect_surface.hpp"
#include "fiction/technology/sidb_on_the_fly_mini_gate_library.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/utils/gate_design_utils.hpp"

#include <fiction/io/write_svg_layout.hpp>

#include <mutex>
#include <thread>

#if (PROGRESS_BARS)
#include <mockturtle/utils/progress_bar.hpp>
#endif

#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
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
    exact_physical_design_params            exact_design_parameters = {};
    std::vector<kitty::dynamic_truth_table> spec{};
    uint64_t                                num_trials                           = 500;
    double                                  quantization_factor                  = 0.075;
    double                                  selectivity                          = 0.5;
    uint64_t                                num_trials_for_double_scope          = 100;
    double                                  quantization_factor_for_double_scope = 0.005;
    double                                  selectivity_for_double_scope         = 0.6;
    uint64_t                                num_trials_for_global_scope          = 20;
    double                                  quantization_factor_for_global_scope = 0.025;
    double                                  selectivity_for_global_scope         = 0.8;
    double                                  excited_state_alpha                  = 1.0;
    uint64_t                                available_threads                    = std::thread::hardware_concurrency();
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

    [[nodiscard]] std::optional<sidb_defect_surface<CellLyt>> design_circuit_on_defective_surface()
    {
        const mockturtle::stopwatch stop{stats.time_total};

        std::optional<GateLyt> gate_lyt = std::nullopt;

        CellLyt lyt{};

        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename GateLibrary::fcn_gate>>
            operational_gate_designs{};

        // // generating the blacklist based on neutral defects. The long-range electrostatic influence of charged
        // defects
        // // is not considered as gates are designed on-the-fly.
        // auto black_list = sidb_surface_analysis<SkeletonGateLibrary>(
        //     lattice_tiling, params.sidb_on_the_fly_gate_library_parameters.defect_surface, std::make_pair(0, 0));

        while (!gate_lyt.has_value())
        {
            // P&R with *exact* and the pre-determined blacklist
            gate_lyt =
                // exact_with_blacklist<GateLyt>(network, black_list, params.exact_design_parameters,
                // &stats.exact_stats);
                exact_with_blacklist<GateLyt>(network, {}, params.exact_design_parameters, &stats.exact_stats);

            if (!gate_lyt.has_value())
            {
                // P&R was unsuccessful
                std::cout << "UNSUCCESS" << std::endl;
                break;
            }

            circuit.emplace(*gate_lyt, params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                           .operational_params.input_bdl_iterator_params.bdl_wire_params);

            is_circuit_operational_params operational_params{};
            operational_params.simulation_parameters = params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                                           .operational_params.simulation_parameters;
            operational_params.input_bdl_iterator_params =
                params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                    .input_bdl_iterator_params;
            operational_params.termination_cond =
                is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;
            // is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;
            operational_params.excited_state_alpha = params.excited_state_alpha;

            operational_gate_designs.clear();

            try
            {
                gate_lyt->foreach_node(
                    [&, this](const auto& n, [[maybe_unused]] auto i)
                    {
                        if (!skip_physical_design_for_node(*gate_lyt, n))
                        {
                            const auto t = gate_lyt->get_tile(n);

                            operational_gate_designs[n] =
                                GateLibrary::template set_up_gates<GateLyt, CellLyt,
                                                                   sidb_on_the_fly_gate_library_params<CellLyt>,
                                                                   local_external_potential_type::BOUNDED>(
                                    *gate_lyt, t, params.sidb_on_the_fly_gate_library_parameters, std::nullopt,
                                    std::make_optional(*circuit), operational_params);
                        }
                    });
            }

            catch (const gate_design_exception<tt, GateLyt>& e)
            {
                throw unsuccessful_gate_design_error("Gate design was unsuccessful");
            }

            // // on-the-fly gate design was unsuccessful at a certain tile. Hence, this tile-gate pair is added to the
            // // blacklist and the process is rerun.
            // catch (const gate_design_exception<tt, GateLyt>& e)
            // {
            //     black_list[e.which_tile()][e.which_truth_table()].push_back(e.which_port_list());
            // }
        }

        std::optional<sidb_defect_surface<CellLyt>> sidbs_and_defects{};

        if (!prune_gate_designs_by_gate_connections(*gate_lyt, operational_gate_designs) ||
            !prune_gate_designs_by_two_gate_connections(*gate_lyt, operational_gate_designs) ||
            !prune_gate_designs_at_global_level(*gate_lyt, operational_gate_designs, lyt) ||
            !look_for_operational_circuit_exhaustively(*gate_lyt, operational_gate_designs, lyt))
        {
            std::cout << "\n\nFAILURE: NO OPERATIONAL CIRCUIT COULD BE GENERATED" << std::endl;
        }
        else
        {
            std::cout << "\n\nSUCCESS! GENERATED OPERATIONAL CIRCUIT:" << std::endl;
            print_layout(lyt);

            sidbs_and_defects.emplace(lyt);
        }

        stats.gate_layout = std::optional{gate_lyt};

        // // add defects to the circuit.
        // params.sidb_on_the_fly_gate_library_parameters.defect_surface.foreach_sidb_defect(
        //     [&sidbs_and_defects](const auto& defect)
        //     { sidbs_and_defects.assign_sidb_defect(defect.first, defect.second); });

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
    advanced_circuit_design_stats<GateLyt>&                                stats;
    std::optional<sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>> circuit{};

    static void
    collect_connecting_nodes(const GateLyt& gate_lyt, const tile<GateLyt>& t,
                             std::vector<mockturtle::node<GateLyt>>& connecting_to_t,
                             const std::optional<tile<GateLyt>>&     maybe_do_not_collect = std::nullopt) noexcept
    {
        for (const auto& in_t : gate_lyt.incoming_data_flow(t))
        {
            if (!gate_lyt.is_pi(gate_lyt.get_node(in_t)) &&
                (!maybe_do_not_collect.has_value() || gate_lyt.below(in_t) != *maybe_do_not_collect))
            {
                connecting_to_t.emplace_back(gate_lyt.get_node(gate_lyt.below(in_t)));
                std::cout << "inc tile: " << gate_lyt.below(in_t).x << ',' << gate_lyt.below(in_t).y << std::endl;
            }
        }

        for (const auto& out_t : gate_lyt.outgoing_data_flow(t))
        {
            if (!gate_lyt.is_po(gate_lyt.get_node(out_t)) &&
                (!maybe_do_not_collect.has_value() || gate_lyt.below(out_t) != *maybe_do_not_collect))
            {
                connecting_to_t.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
                std::cout << "out tile: " << gate_lyt.below(out_t).x << ',' << gate_lyt.below(out_t).y << std::endl;
            }
        }

        if (gate_lyt.is_wire_tile(t))
        {
            if (const auto above_t = gate_lyt.above(t); t != above_t && gate_lyt.is_wire_tile(above_t))
            {
                collect_connecting_nodes(gate_lyt, above_t, connecting_to_t, maybe_do_not_collect);
            }
        }
    }

    static void apply_quantization(std::vector<std::pair<double, uint64_t>>& sorted_success_ratios,
                                   const double                              quantization_factor) noexcept
    {
        if (sorted_success_ratios.empty() || quantization_factor <= 0.0)
        {
            return;
        }

        const double min_sr = sorted_success_ratios.front().first;
        const double max_sr = sorted_success_ratios.back().first;

        // If all values are identical, nothing to do.
        if (std::abs(min_sr - max_sr) < std::numeric_limits<double>::epsilon())
        {
            return;
        }

        const double scale = std::floor(1.0 / quantization_factor);

        for (auto& [success_rate, _] : sorted_success_ratios)
        {
            success_rate = std::round(success_rate * scale) / scale;
        }
    }

    [[nodiscard]] static uint64_t determine_first_passing_gate_ix(
        const std::vector<std::pair<double, uint64_t>>& successful_trial_ratio_per_gate_implementation,
        const double                                    selectivity) noexcept
    {
        auto first_passing_gate_ix = static_cast<int64_t>(
            selectivity * static_cast<double>(successful_trial_ratio_per_gate_implementation.size()));

        const bool same_as_first =
            std::abs(
                successful_trial_ratio_per_gate_implementation.front().first -
                successful_trial_ratio_per_gate_implementation.at(static_cast<uint64_t>(first_passing_gate_ix)).first) <
            std::numeric_limits<double>::epsilon();

        const bool same_as_last =
            std::abs(
                successful_trial_ratio_per_gate_implementation.back().first -
                successful_trial_ratio_per_gate_implementation.at(static_cast<uint64_t>(first_passing_gate_ix)).first) <
            std::numeric_limits<double>::epsilon();

        if (same_as_first && same_as_last)
        {
            return 0;
        }

        // todo update comment
        // if the success ratio of the first passing gate index is the same lowest and there also a higher success
        // ratio dir = +1, otherwise, dir = -1

        // for the case of dir = -1
        // walk down (i.e., decrementing success ratio) to the first gate implementation with a different success
        // ratio opposite for dir = +1

        std::array<int64_t, 2> first_passing_gate_up_or_down{{first_passing_gate_ix, first_passing_gate_ix}};

        for (const int64_t dir : same_as_first || same_as_last ? std::vector<int64_t>{same_as_last ? -1 : +1} :
                                                                 std::vector<int64_t>{{-1, +1}})
        {
            const uint64_t index = static_cast<uint64_t>(dir + 1) / 2;

            for (; std::abs(successful_trial_ratio_per_gate_implementation
                                .at(static_cast<uint64_t>(first_passing_gate_up_or_down[index] + dir))
                                .first -
                            successful_trial_ratio_per_gate_implementation
                                .at(static_cast<uint64_t>(first_passing_gate_up_or_down[index]))
                                .first) < std::numeric_limits<double>::epsilon();
                 first_passing_gate_up_or_down[index] += dir)
            {}

            if (dir == +1)
            {
                ++first_passing_gate_up_or_down[index];
            }
        }

        if (same_as_last)
        {
            return static_cast<uint64_t>(first_passing_gate_up_or_down[0]);
        }

        if (same_as_first)
        {
            return static_cast<uint64_t>(first_passing_gate_up_or_down[1]);
        }

        return static_cast<uint64_t>(
            first_passing_gate_up_or_down[first_passing_gate_up_or_down[1] - first_passing_gate_ix >=
                                                  first_passing_gate_ix - first_passing_gate_up_or_down[0] ?
                                              0 :
                                              1]);
    }

    static void print_success_rate_distribution(const std::vector<std::pair<double, uint64_t>>& success_ratios,
                                                const uint64_t first_passing_gate_ix) noexcept
    {
        const double selectivity_threshold = success_ratios.at(first_passing_gate_ix).first;

        std::cout << "\nDetermined success threshold: " << std::fixed << std::setprecision(1)
                  << selectivity_threshold * 100
                  << fmt::format("% — pruning {} out of {} gate designs — reduced pool by {:.1f}%\n\n",
                                 first_passing_gate_ix, success_ratios.size(),
                                 static_cast<double>(first_passing_gate_ix) /
                                     static_cast<double>(success_ratios.size()) * 100);

        // Map: success_rate -> count of gates with this rate
        std::map<double, uint64_t> count_by_success_rate;

        // Count how many gates per success rate
        for (const auto& [rate, _] : success_ratios)
        {
            count_by_success_rate[rate]++;
        }

        // Total number of gates
        const uint64_t total_gates = success_ratios.size();

        assert(total_gates > 0 && "The distribution is empty.");

        std::cout << "Success Rate | Status | Gate Count | Graph\n";
        std::cout
            << "-------------|--------|------------|--------------------------------------------------------------\n";

        for (const auto& [rate, count] : count_by_success_rate)
        {
            constexpr size_t max_bar_length = 60;

            const bool        is_kept = rate >= selectivity_threshold;
            const std::string status  = is_kept ? "KEPT" : "PRUNED";

            // Bar length proportional to total count
            const auto bar_len = 1 + static_cast<size_t>(static_cast<double>(count) / static_cast<double>(total_gates) *
                                                         (max_bar_length - 1));

            std::string bar(bar_len, is_kept ? 'o' : 'X');

            std::cout << std::right << std::setw(11) << std::fixed << std::setprecision(1) << rate * 100 << "%"
                      << " | " << std::setw(6) << status << " | " << std::setw(10) << count << " | " << std::left
                      << std::setw(max_bar_length) << bar << '\n';
        }

        std::cout << std::endl;
    }

    [[nodiscard]] std::vector<uint64_t> select_gate_implementations_by_successful_trial_ratio(
        std::vector<std::pair<double, uint64_t>>&& successful_trial_ratio_per_gate_implementation,
        const double quantization_factor, const double selectivity) const noexcept
    {
        std::sort(successful_trial_ratio_per_gate_implementation.begin(),
                  successful_trial_ratio_per_gate_implementation.end(),
                  [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

        apply_quantization(successful_trial_ratio_per_gate_implementation, quantization_factor);

        const uint64_t first_passing_gate_ix =
            determine_first_passing_gate_ix(successful_trial_ratio_per_gate_implementation, selectivity);

        print_success_rate_distribution(successful_trial_ratio_per_gate_implementation, first_passing_gate_ix);

        std::vector<uint64_t> selected_gate_implementation_indices{};
        selected_gate_implementation_indices.reserve(successful_trial_ratio_per_gate_implementation.size() -
                                                     first_passing_gate_ix);

        for (uint64_t current_tile_gate_implementation_index = first_passing_gate_ix;
             current_tile_gate_implementation_index < successful_trial_ratio_per_gate_implementation.size();
             ++current_tile_gate_implementation_index)
        {
            selected_gate_implementation_indices.push_back(
                successful_trial_ratio_per_gate_implementation.at(current_tile_gate_implementation_index).second);
        }

        return selected_gate_implementation_indices;
    }
    /**
     *
     */
    bool prune_gate_designs_by_gate_connections(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename GateLibrary::fcn_gate>>&
            operational_gate_designs) const noexcept
    {
        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS BY GATE CONNECTIONS" << std::endl;

        std::unordered_map<
            mockturtle::node<GateLyt>,
            std::unordered_map<mockturtle::node<GateLyt>, sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>>>
            gate_lyt_window_for_joint_simulation{};  // todo: optimise through symmetry
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        is_circuit_operational_params operational_params{};
        operational_params.simulation_parameters =
            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params = params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                                           .operational_params.input_bdl_iterator_params;
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;
        operational_params.excited_state_alpha = params.excited_state_alpha;

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const tile<GateLyt> t = gate_lyt.get_tile(n);

                // if (t.x != 1 || t.y != 1)
                //     return;

                std::cout << "\nNOW DOING TILE " << t.x << ',' << t.y << ',' << t.z << std::endl;
                std::cout << "tt: ";
                kitty::print_binary(gate_lyt.node_function(n));
                if (const auto& at = gate_lyt.above(t); t != at && gate_lyt.is_wire_tile(at))
                {

                    std::cout << "\t\tabove tt:";
                    kitty::print_binary(gate_lyt.node_function(gate_lyt.get_node(at)));
                }
                std::cout << std::endl;

                std::vector<mockturtle::node<GateLyt>> connecting_to_n{};

                collect_connecting_nodes(gate_lyt, t, connecting_to_n);

                for (const auto& connecting_n : connecting_to_n)
                {
                    const tile<GateLyt>& connecting_t = gate_lyt.get_tile(connecting_n);

                    // std::cout << "connecting " << gate_lyt.get_tile(connecting_n).x << ',' <<
                    // gate_lyt.get_tile(connecting_n).y << ',' << gate_lyt.get_tile(connecting_n).z << std::endl;

                    assert(gate_lyt_window_for_joint_simulation.count(n) == 0 ||
                           gate_lyt_window_for_joint_simulation.at(n).count(connecting_n) == 0);

                    gate_lyt_window_for_joint_simulation[n].insert(
                        {connecting_n, sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                                           gate_lyt, t, connecting_t,
                                           params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                               .operational_params.input_bdl_iterator_params.bdl_wire_params}});
                }
            });

        bool big_fixpoint = false;

        bool exit_by_failure = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || skip_physical_design_for_node(gate_lyt, n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = gate_lyt.get_tile(n);

                    // if (t.y == 1)
                    //     return;

                    // if (t.x != 1 || t.y != 1)
                    //     return;

                    std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                    std::cout << fmt::format(
                                     "Performing {} trials for {} gate connection{} for {} gate implementations",
                                     params.num_trials, gate_lyt_window_for_joint_simulation.at(n).size(),
                                     gate_lyt_window_for_joint_simulation.at(n).size() > 1 ? "s" : "",
                                     operational_gate_designs.at(n).size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};
                    std::mutex                               mutex_to_protect_successful_trial_ratios;

                    // const uint64_t num_threads = 1;
                    const uint64_t num_threads =
                        std::min(params.available_threads, operational_gate_designs.at(n).size());
                    //  const uint64_t num_threads = std::min(uint64_t{20},
                    //  operational_gate_designs.at(n).size());

                    const uint64_t chunk_size =
                        (operational_gate_designs.at(n).size() + num_threads - 1) /
                        num_threads;  // Ceiling division

                    std::vector<std::thread> threads{};
                    threads.reserve(num_threads);

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{static_cast<uint32_t>(std::min(
                                                     chunk_size, operational_gate_designs.at(n).size())),
                                                 "[i] Determining successful trial ratio for tile " +
                                                     fmt::format("({},{})", t.x, t.y) + ": |{0}|"};
#endif

                    for (uint64_t i = 0; i < num_threads; ++i)
                    {
                        threads.emplace_back(
                            [&successful_trial_ratio_per_gate_implementation, &t, i, chunk_size,
                             &operational_gate_designs, &n, &gate_lyt, this, &operational_params,
#if (PROGRESS_BARS)

                             &bar,
#endif
                             &gate_lyt_window_for_joint_simulation, &mutex_to_protect_successful_trial_ratios]
                            {
                                const uint64_t start_index = i * chunk_size;
                                const uint64_t end_index   = std::min(
                                    start_index + chunk_size, operational_gate_designs.at(n).size());

                                for (uint64_t j = start_index; j < end_index; ++j)
                                {
                                    // std::cout << "gate design ix: " << j << std::endl;

                                    CellLyt cell_lyt{};

                                    // select the first gate implementation for n
                                    assign_gate<CellLyt, GateLibrary, GateLyt>(
                                        cell_lyt,
                                        *std::next(operational_gate_designs.at(n).cbegin(),
                                                   static_cast<int64_t>(j)),
                                        gate_lyt, t);

                                    double successful_trials = 0;

                                    // uint64_t gate_connection_ix = 0;
                                    for (const auto& [connecting_n, gate_lyt_window] :
                                         gate_lyt_window_for_joint_simulation.at(n))
                                    {
                                        // std::cout << "connecting " << gate_lyt.get_tile(connecting_n).x << ','
                                        //           << gate_lyt.get_tile(connecting_n).y << ','
                                        //           << gate_lyt.get_tile(connecting_n).z << std::endl;
                                        // GateLyt gate_lyt_window_for_gate_connection{};

                                        // std::cout << "gate connection ix: " << gate_connection_ix++ << std::endl;

                                        uint64_t current_trial = 0;

                                        // double   avg_time  = 0;
                                        // uint64_t sample_ix = 0;

                                        while (current_trial < params.num_trials)
                                        {
                                            current_trial++;
                                            // std::cout << "trial number: " << current_trial << std::endl;

                                            CellLyt cell_lyt_clone = cell_lyt.clone();

                                            std::random_device rd;         // a seed source for the random number engine
                                            std::mt19937       gen(rd());  // mersenne_twister_engine seeded with rd()
                                            std::uniform_int_distribution<uint64_t> distrib{
                                                0, operational_gate_designs.at(connecting_n).size() - 1};

                                            // select a random gate implementation for the tile that connects as
                                            // input to n
                                            assign_gate<CellLyt, GateLibrary, GateLyt>(
                                                cell_lyt_clone,
                                                operational_gate_designs.at(connecting_n)
                                                    .at(distrib(gen)),
                                                gate_lyt, gate_lyt.get_tile(connecting_n));

                                            // mockturtle::stopwatch<>::duration time_counter{};
                                            // {
                                            //     mockturtle::stopwatch stop{time_counter};

                                            // std::cout << "TT here:";
                                            // kitty::print_binary(gate_lyt_window_for_joint_simulation.at(n)
                                            //             .at(connecting_n)
                                            //             .first.node_function(gate_lyt_window_for_joint_simulation.at(n)
                                            //             .at(connecting_n)
                                            //             .first.get_node({1,2,0})));
                                            // std::cout << std::endl;
                                            //
                                            // std::cout << "is po here: " <<
                                            // gate_lyt_window_for_joint_simulation.at(n)
                                            //             .at(connecting_n)
                                            //             .first.is_po_tile({1,2,0}) << std::endl;

                                            // std::cout
                                            //     << "gate level inputs: "
                                            //     << gate_lyt_window_for_joint_simulation.at(n)
                                            //            .at(connecting_n)
                                            //            .first.gate_layout.num_pis()
                                            //     << std::endl;
                                            // std::cout << "cell level inputs: " << cell_lyt_clone.num_pis() <<
                                            // std::endl;

                                            const circuit_operational_assessment<
                                                CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
                                                is_circuit_operational<CellLyt, GateLyt,
                                                                       local_external_potential_type::BOUNDED,
                                                                       SkeletonGateLibrary>(
                                                    sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                                                        cell_lyt_clone, gate_lyt_window},
                                                    operational_params,
                                                    std::make_optional<std::reference_wrapper<
                                                        const sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>>>(
                                                        std::cref(*circuit)));

                                            assert(op_assessment.assessment_per_input.has_value() &&
                                                   "ALL_COMBINATIONS_ENUMERATED is not set.");

                                            // if (op_assessment.status == operational_status::OPERATIONAL)
                                            // {
                                            double logic_match_sum_for_all_inputs = 0.0;

                                            for (const typename circuit_operational_assessment<
                                                     CellLyt, local_external_potential_type::BOUNDED>::
                                                     operational_assessment_for_input& op_assessment_for_input :
                                                 *op_assessment.assessment_per_input)
                                            {
                                                logic_match_sum_for_all_inputs +=
                                                    // op_assessment_for_input.status ==
                                                    //         operational_status::OPERATIONAL ?
                                                    //     1.0 :
                                                    //     0.0;
                                                    op_assessment_for_input.logic_match;
                                            }

                                            successful_trials +=
                                                logic_match_sum_for_all_inputs /
                                                static_cast<double>(op_assessment.assessment_per_input->size());
                                            // }
                                            // }

                                            // if (i == 0)
                                            // {
                                            //     sample_ix++;
                                            //
                                            //     if (avg_time == 0)
                                            //     {
                                            //         avg_time = time_counter.count();
                                            //     }
                                            //     else
                                            //     {
                                            //         avg_time =
                                            //             ((sample_ix - 1) * avg_time + time_counter.count()) /
                                            //             sample_ix;
                                            //     }
                                            //
                                            //     std::cout << fmt::format("time taken: {:.3f} s | avg time: {:.3f}
                                            //     s",
                                            //                              time_counter.count() / 1e9, avg_time /
                                            //                              1e9)
                                            //               << std::endl;
                                            // }
                                        }
                                    }

#if (PROGRESS_BARS)
                                    if (i == 0)
                                    {
                                        bar(j);
                                    }
#endif

                                    const double successful_trial_ratio =
                                        successful_trials /
                                        static_cast<double>(params.num_trials *
                                                            gate_lyt_window_for_joint_simulation.at(n).size());

                                    const std::lock_guard lock{mutex_to_protect_successful_trial_ratios};

                                    successful_trial_ratio_per_gate_implementation.emplace_back(successful_trial_ratio,
                                                                                                j);
                                }
                            });
                    }

                    for (auto& thread : threads)
                    {
                        if (thread.joinable())
                        {
                            thread.join();
                        }
                    }

                    selected_gate_implementation_indices[n] = select_gate_implementations_by_successful_trial_ratio(
                        std::move(successful_trial_ratio_per_gate_implementation), params.quantization_factor,
                        params.selectivity);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() <
                             operational_gate_designs.at(n).size())
                    {
                        big_fixpoint = false;
                    }
                });

            std::cout << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || skip_physical_design_for_node(gate_lyt, n))
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

    bool prune_gate_designs_by_two_gate_connections(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename GateLibrary::fcn_gate>>&
            operational_gate_designs) const noexcept
    {
        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS BY TWO GATE CONNECTIONS" << std::endl;

        std::unordered_map<
            mockturtle::node<GateLyt>,
            std::unordered_map<
                mockturtle::node<GateLyt>,
                std::unordered_map<mockturtle::node<GateLyt>, sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>>>>
            gate_lyt_window_for_joint_simulation{};  // todo: optimise through symmetry
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        is_circuit_operational_params operational_params{};
        operational_params.simulation_parameters =
            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params = params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                                           .operational_params.input_bdl_iterator_params;
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;
        operational_params.excited_state_alpha = params.excited_state_alpha;

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const tile<GateLyt> t = gate_lyt.get_tile(n);

                // if (t.x != 0 || t.y != 2)
                //     return;

                std::cout << "\nNOW DOING TILE " << t.x << ',' << t.y << ',' << t.z << std::endl;
                std::cout << "tt: ";
                kitty::print_binary(gate_lyt.node_function(n));
                if (const auto& at = gate_lyt.above(t); t != at && gate_lyt.is_wire_tile(at))
                {

                    std::cout << "\t\tabove tt:";
                    kitty::print_binary(gate_lyt.node_function(gate_lyt.get_node(at)));
                }
                std::cout << std::endl;

                std::vector<mockturtle::node<GateLyt>> connecting_to_n{};

                collect_connecting_nodes(gate_lyt, t, connecting_to_n);

                std::unordered_map<mockturtle::node<GateLyt>, std::vector<mockturtle::node<GateLyt>>>
                    second_gate_connections_to_connecting_n{};

                for (const auto& connecting_n : connecting_to_n)
                {
                    const tile<GateLyt>& connecting_t = gate_lyt.get_tile(connecting_n);

                    std::vector<mockturtle::node<GateLyt>> second_gate_connections{};

                    collect_connecting_nodes(gate_lyt, connecting_t, second_gate_connections, std::make_optional(t));

                    second_gate_connections_to_connecting_n[connecting_n] = std::move(second_gate_connections);
                }

                for (const auto& connecting_n : connecting_to_n)
                {
                    const tile<GateLyt>& connecting_t = gate_lyt.get_tile(connecting_n);

                    // std::cout << "connecting " << connecting_t.x << ',' << connecting_t.y << ',' << connecting_t.z
                    //           << std::endl;

                    for (const auto& connecting_to_connecting_n :
                         second_gate_connections_to_connecting_n.at(connecting_n))
                    {
                        const tile<GateLyt>& connecting_to_connecting_t = gate_lyt.get_tile(connecting_to_connecting_n);

                        // std::cout << "connecting further " << connecting_to_connecting_t.x << ','
                        //           << connecting_to_connecting_t.y << ',' << connecting_to_connecting_t.z <<
                        //           std::endl;

                        gate_lyt_window_for_joint_simulation[n][connecting_n].insert(
                            {connecting_to_connecting_n,
                             sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                                 gate_lyt, t, connecting_t, connecting_to_connecting_t,
                                 params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                                     .input_bdl_iterator_params.bdl_wire_params}});
                    }

                    // assert(gate_lyt_window_for_joint_simulation.count(n) == 0 ||
                    //        gate_lyt_window_for_joint_simulation.at(n).count(connecting_n) == 0);

                    for (const auto& other_connecting_n : connecting_to_n)
                    {
                        const tile<GateLyt>& other_connecting_t = gate_lyt.get_tile(other_connecting_n);

                        if (connecting_t <= other_connecting_t)
                        {
                            continue;
                        }

                        std::cout << "connecting with " << other_connecting_t.x << ',' << other_connecting_t.y << ','
                                  << other_connecting_t.z << std::endl;

                        gate_lyt_window_for_joint_simulation[n][connecting_n].insert(
                            {other_connecting_n,
                             sidb_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                                 gate_lyt, t, connecting_t, other_connecting_t,
                                 params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                                     .input_bdl_iterator_params.bdl_wire_params}});
                    }
                }
            });

        bool big_fixpoint = false;

        bool exit_by_failure = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || skip_physical_design_for_node(gate_lyt, n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = gate_lyt.get_tile(n);

                    // if (t.x != 2 || t.y != 3)
                    //     return;

                    uint64_t number_of_layouts = 0;

                    for (const auto& [_, map] : gate_lyt_window_for_joint_simulation.at(n))
                    {
                        number_of_layouts += map.size();
                    }

                    std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                    std::cout << fmt::format(
                                     "Performing {} trials for {} gate connection{} for {} gate implementations",
                                     params.num_trials_for_double_scope, number_of_layouts,
                                     number_of_layouts > 1 ? "s" : "",
                                     operational_gate_designs.at(n).size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};
                    std::mutex                               mutex_to_protect_successful_trial_ratios;

                    // const uint64_t num_threads = 1;
                    const uint64_t num_threads =
                        std::min(params.available_threads, operational_gate_designs.at(n).size());
                    //  const uint64_t num_threads = std::min(uint64_t{20},
                    //  operational_gate_designs.at(n).size());

                    const uint64_t chunk_size =
                        (operational_gate_designs.at(n).size() + num_threads - 1) /
                        num_threads;  // Ceiling division

                    std::vector<std::thread> threads{};
                    threads.reserve(num_threads);

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{static_cast<uint32_t>(std::min(
                                                     chunk_size, operational_gate_designs.at(n).size())),
                                                 "[i] Determining successful trial ratio for tile " +
                                                     fmt::format("({},{})", t.x, t.y) + ": |{0}|"};
#endif

                    for (uint64_t i = 0; i < num_threads; ++i)
                    {
                        threads.emplace_back(
                            [&successful_trial_ratio_per_gate_implementation, &t, i, chunk_size,
                             &operational_gate_designs, &n, &gate_lyt, this, &operational_params,
#if (PROGRESS_BARS)
                             &bar,
#endif
                             &gate_lyt_window_for_joint_simulation, number_of_layouts,
                             &mutex_to_protect_successful_trial_ratios]
                            {
                                const uint64_t start_index = i * chunk_size;
                                const uint64_t end_index   = std::min(
                                    start_index + chunk_size, operational_gate_designs.at(n).size());

                                for (uint64_t j = start_index; j < end_index; ++j)
                                {
                                    // std::cout << "gate design ix: " << j << std::endl;

                                    CellLyt cell_lyt{};

                                    // select the first gate implementation for n
                                    assign_gate<CellLyt, GateLibrary, GateLyt>(
                                        cell_lyt,
                                        *std::next(operational_gate_designs.at(n).cbegin(),
                                                   static_cast<int64_t>(j)),
                                        gate_lyt, t);

                                    double successful_trials = 0;

                                    // uint64_t gate_connection_ix = 0;
                                    for (const auto& [connecting_n, connecting_to_connection_map] :
                                         gate_lyt_window_for_joint_simulation.at(n))
                                    {
                                        // std::cout << "connecting " << gate_lyt.get_tile(connecting_n).x << ','
                                        //           << gate_lyt.get_tile(connecting_n).y << ','
                                        //           << gate_lyt.get_tile(connecting_n).z << std::endl;

                                        for (const auto& [connecting_to_connecting_n, gate_lyt_window] :
                                             connecting_to_connection_map)
                                        {
                                            // std::cout << "connecting further "
                                            //           << gate_lyt.get_tile(connecting_to_connecting_n).x << ','
                                            //           << gate_lyt.get_tile(connecting_to_connecting_n).y << ','
                                            //           << gate_lyt.get_tile(connecting_to_connecting_n).z <<
                                            //           std::endl;

                                            // GateLyt gate_lyt_window_for_gate_connection{};

                                            // std::cout << "gate connection ix: " << gate_connection_ix++ << std::endl;

                                            uint64_t current_trial = 0;

                                            // double   avg_time  = 0;
                                            // uint64_t sample_ix = 0;

                                            while (current_trial < params.num_trials_for_double_scope)
                                            {
                                                current_trial++;
                                                // std::cout << "trial number: " << current_trial << std::endl;

                                                CellLyt cell_lyt_clone = cell_lyt.clone();

                                                std::random_device rd;   // a seed source for the random number engine
                                                std::mt19937 gen(rd());  // mersenne_twister_engine seeded with rd()
                                                std::uniform_int_distribution<uint64_t> distrib{
                                                    0, operational_gate_designs.at(connecting_n).size() -
                                                           1};

                                                // select a random gate implementation for the tile that connects as
                                                // input to n
                                                assign_gate<CellLyt, GateLibrary, GateLyt>(
                                                    cell_lyt_clone,
                                                    operational_gate_designs.at(connecting_n)
                                                        .at(distrib(gen)),
                                                    gate_lyt, gate_lyt.get_tile(connecting_n));

                                                std::random_device rd2;    // a seed source for the random number engine
                                                std::mt19937 gen2(rd2());  // mersenne_twister_engine seeded with rd()
                                                std::uniform_int_distribution<uint64_t> distrib2{
                                                    0, operational_gate_designs.at(connecting_to_connecting_n)
                                                               .size() -
                                                           1};

                                                // select a random gate implementation for the tile that connects as
                                                // input to n
                                                assign_gate<CellLyt, GateLibrary, GateLyt>(
                                                    cell_lyt_clone,
                                                    operational_gate_designs.at(connecting_to_connecting_n)
                                                        .at(distrib2(gen2)),
                                                    gate_lyt, gate_lyt.get_tile(connecting_to_connecting_n));

                                                // mockturtle::stopwatch<>::duration time_counter{};
                                                // {
                                                //     mockturtle::stopwatch stop{time_counter};

                                                // std::cout << "TT here:";
                                                // kitty::print_binary(gate_lyt_window_for_joint_simulation.at(n)
                                                //             .at(connecting_n)
                                                //             .first.node_function(gate_lyt_window_for_joint_simulation.at(n)
                                                //             .at(connecting_n)
                                                //             .first.get_node({1,2,0})));
                                                // std::cout << std::endl;
                                                //
                                                // std::cout << "is po here: " <<
                                                // gate_lyt_window_for_joint_simulation.at(n)
                                                //             .at(connecting_n)
                                                //             .first.is_po_tile({1,2,0}) << std::endl;

                                                // std::cout
                                                //     << "gate level inputs: " << gate_lyt_window.gate_layout.num_pis()
                                                //     << std::endl;
                                                // std::cout << "cell level inputs: " << cell_lyt_clone.num_pis()
                                                //           << std::endl;

                                                const circuit_operational_assessment<
                                                    CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
                                                    is_circuit_operational<CellLyt, GateLyt,
                                                                           local_external_potential_type::BOUNDED,
                                                                           SkeletonGateLibrary>(
                                                        sidb_cell_level_bdl_circuit<CellLyt, GateLyt,
                                                                                    SkeletonGateLibrary>{
                                                            cell_lyt_clone, gate_lyt_window},
                                                        operational_params,
                                                        std::make_optional<
                                                            std::reference_wrapper<const sidb_bdl_circuit<
                                                                CellLyt, GateLyt, SkeletonGateLibrary>>>(
                                                            std::cref(*circuit)));

                                                assert(op_assessment.assessment_per_input.has_value() &&
                                                       "ALL_COMBINATIONS_ENUMERATED is not set.");

                                                // if (op_assessment.status == operational_status::OPERATIONAL)
                                                // {
                                                double logic_match_sum_for_all_inputs = 0.0;

                                                for (const typename circuit_operational_assessment<
                                                         CellLyt, local_external_potential_type::BOUNDED>::
                                                         operational_assessment_for_input& op_assessment_for_input :
                                                     *op_assessment.assessment_per_input)
                                                {
                                                    logic_match_sum_for_all_inputs +=
                                                        // op_assessment_for_input.status ==
                                                        //         operational_status::OPERATIONAL ?
                                                        //     1.0 :
                                                        //     0.0;
                                                        op_assessment_for_input.logic_match;
                                                }

                                                successful_trials +=
                                                    logic_match_sum_for_all_inputs /
                                                    static_cast<double>(op_assessment.assessment_per_input->size());
                                                // }
                                                // }

                                                // if (i == 0)
                                                // {
                                                //     sample_ix++;
                                                //
                                                //     if (avg_time == 0)
                                                //     {
                                                //         avg_time = time_counter.count();
                                                //     }
                                                //     else
                                                //     {
                                                //         avg_time =
                                                //             ((sample_ix - 1) * avg_time + time_counter.count()) /
                                                //             sample_ix;
                                                //     }
                                                //
                                                //     std::cout << fmt::format("time taken: {:.3f} s | avg time: {:.3f}
                                                //     s",
                                                //                              time_counter.count() / 1e9, avg_time /
                                                //                              1e9)
                                                //               << std::endl;
                                                // }
                                            }
                                        }
                                    }

#if (PROGRESS_BARS)
                                    if (i == 0)
                                    {
                                        bar(j);
                                    }
#endif

                                    const double successful_trial_ratio =
                                        successful_trials /
                                        static_cast<double>(params.num_trials_for_double_scope * number_of_layouts);

                                    const std::lock_guard lock{mutex_to_protect_successful_trial_ratios};

                                    successful_trial_ratio_per_gate_implementation.emplace_back(successful_trial_ratio,
                                                                                                j);
                                }
                            });
                    }

                    for (auto& thread : threads)
                    {
                        if (thread.joinable())
                        {
                            thread.join();
                        }
                    }

                    selected_gate_implementation_indices[n] = select_gate_implementations_by_successful_trial_ratio(
                        std::move(successful_trial_ratio_per_gate_implementation),
                        params.quantization_factor_for_double_scope, params.selectivity_for_double_scope);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() <
                             operational_gate_designs.at(n).size())
                    {
                        big_fixpoint = false;
                    }
                });

            std::cout << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || skip_physical_design_for_node(gate_lyt, n))
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

    bool prune_gate_designs_at_global_level(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename GateLibrary::fcn_gate>>&
                 operational_gate_designs,
        CellLyt& lyt) noexcept
    {
        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS WITH GLOBAL SCOPE" << std::endl;

        bool operational_circuit_found = false;

        is_circuit_operational_params operational_params{};
        operational_params.simulation_parameters =
            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params = params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                                                           .operational_params.input_bdl_iterator_params;
        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        operational_params.print = true;

        const uint64_t num_input_combinations = 1 << gate_lyt.num_pis();

        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        bool big_fixpoint = false;

        bool exit_by_failure = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (operational_circuit_found || exit_by_failure || skip_physical_design_for_node(gate_lyt, n))
                    {
                        return;
                    }

                    if (operational_gate_designs.at(n).size() == 1)
                    {
                        selected_gate_implementation_indices[n] = {0};

                        return;
                    }

                    const tile<GateLyt>& t = gate_lyt.get_tile(n);

                    std::cout << fmt::format("\nStarting pruning for tile {}\n", t);
                    std::cout << fmt::format("Performing {} trials for {} gate implementations",
                                             params.num_trials_for_global_scope,
                                             operational_gate_designs.at(n).size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{
                        static_cast<uint32_t>(operational_gate_designs.at(n).size()),
                        "[i] Determining successful trial ratio for tile " + fmt::format("({},{})", t.x, t.y) +
                            ": |{0}|\t"};
#endif

                    for (uint64_t j = 0;
                         !operational_circuit_found && j < operational_gate_designs.at(n).size(); ++j)
                    {
                        CellLyt cell_lyt{};

                        // select the first gate implementation for n
                        assign_gate<CellLyt, GateLibrary, GateLyt>(
                            cell_lyt,
                            *std::next(operational_gate_designs.at(n).cbegin(), static_cast<int64_t>(j)),
                            gate_lyt, t);

                        double successful_trials = 0;

                        uint64_t current_trial = 0;

                        while (!operational_circuit_found && current_trial < params.num_trials_for_global_scope)
                        {
                            current_trial++;

                            CellLyt cell_lyt_clone = cell_lyt.clone();

                            gate_lyt.foreach_node(
                                [&](const auto& other_n)
                                {
                                    if (skip_physical_design_for_node(gate_lyt, other_n) || other_n == n)
                                    {
                                        return;
                                    }

                                    std::random_device rd;         // a seed source for the random number engine
                                    std::mt19937       gen(rd());  // mersenne_twister_engine seeded with rd()
                                    std::uniform_int_distribution<uint64_t> distrib{
                                        0, operational_gate_designs.at(other_n).size() - 1};

                                    // select a random gate implementation for the tile that connects as
                                    // input to n
                                    assign_gate<CellLyt, GateLibrary, GateLyt>(
                                        cell_lyt_clone,
                                        operational_gate_designs.at(other_n).at(distrib(gen)), gate_lyt,
                                        gate_lyt.get_tile(other_n));
                                });

                            const circuit_operational_assessment<CellLyt>& op_assessment_for_all_inputs =
                                is_circuit_operational(
                                    sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{cell_lyt_clone,
                                                                                                       *circuit},
                                    operational_params);

                            assert(op_assessment_for_all_inputs.assessment_per_input.has_value() &&
                                   "Operational assessment termination condition was not set to "
                                   "ALL_INPUT_COMBINATIONS_ENUMERATED");

                            uint64_t successful_input_combinations = 0;
                            double   logic_match_sum               = 0.0;

                            for (const typename circuit_operational_assessment<
                                     CellLyt>::operational_assessment_for_input& op_assessment :
                                 op_assessment_for_all_inputs.assessment_per_input.value())
                            {
                                if (op_assessment.status == operational_status::OPERATIONAL)
                                {
                                    successful_input_combinations++;
                                }

                                logic_match_sum += op_assessment.logic_match;
                            }

                            if (successful_input_combinations == num_input_combinations)
                            {
                                // all input combinations are operational---operational circuit found

                                lyt = cell_lyt_clone;

                                operational_circuit_found = true;
                            }
                            else
                            {
                                successful_trials += logic_match_sum / static_cast<double>(num_input_combinations);
                            }
                        }

                        if (operational_circuit_found)
                        {
                            break;
                        }

#if (PROGRESS_BARS)
                        bar(j);
#endif

                        const double successful_trial_ratio =
                            successful_trials /
                            static_cast<double>(params.num_trials_for_global_scope * num_input_combinations);

                        successful_trial_ratio_per_gate_implementation.emplace_back(successful_trial_ratio, j);
                    }

                    if (operational_circuit_found)
                    {
                        return;
                    }

                    selected_gate_implementation_indices[n] = select_gate_implementations_by_successful_trial_ratio(
                        std::move(successful_trial_ratio_per_gate_implementation),
                        params.quantization_factor_for_global_scope, params.selectivity_for_global_scope);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() <
                             operational_gate_designs.at(n).size())
                    {
                        big_fixpoint = false;
                    }
                });

            if (operational_circuit_found || big_fixpoint)
            {
                return true;
            }

            if (exit_by_failure)
            {
                return false;
            }

            std::cout << std::endl;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (skip_physical_design_for_node(gate_lyt, n))
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
    /**
     * todo
     */
    bool look_for_operational_circuit_exhaustively(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename GateLibrary::fcn_gate>>&
                 operational_gate_designs,
        CellLyt& lyt) const noexcept
    {
        std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;

        std::vector<uint64_t> indices(operational_gate_designs.size(), 0);

        while (true)
        {
            CellLyt operational_circuit_candidate{};
            for (uint64_t i = 0; i < operational_gate_designs.size(); i++)
            {
                const auto& [n, op_gate_designs_for_gate] =
                    *std::next(operational_gate_designs.cbegin(), static_cast<int64_t>(i));
                // select a random gate implementation for the tile that connects as input to n
                assign_gate<CellLyt, GateLibrary, GateLyt>(operational_circuit_candidate,
                                                           op_gate_designs_for_gate.at(indices.at(i)),
                                                           gate_lyt, gate_lyt.get_tile(n));
            }

            std::cout << "trying combination: ";
            for (uint64_t i = 0; i < operational_gate_designs.size(); i++)
            {
                std::cout << indices.at(i) << " ";
            }
            std::cout << std::endl;

            if (is_operational(
                    operational_circuit_candidate, params.spec,
                    [](auto ps)
                    {
                        ps.termination_cond = decltype(ps)::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;
                        ps.print            = true;
                        return ps;
                    }(params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params),
                    std::make_optional(gate_lyt))
                    .status == operational_status::OPERATIONAL)
            {
                lyt = operational_circuit_candidate;

                std::cout << "\n\nFINAL GENERATED CIRCUIT:" << std::endl;
                print_layout(lyt);

                return true;
            }

            // Increment indices like an odometer
            for (uint64_t i = 0; i < indices.size(); ++i)
            {
                if (++indices[i] <
                    std::next(operational_gate_designs.cbegin(), static_cast<int64_t>(i))->second.size())
                {
                    break;  // No carry needed
                }

                indices[i] = 0;  // Reset this index and carry over to the next

                if (i == indices.size() - 1)
                {
                    return false;  // Stop when the last index overflows
                }
            }
        }

        return false;
    }
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
[[nodiscard]] std::optional<sidb_defect_surface<CellLyt>>
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

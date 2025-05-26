//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/algorithms/simulation/sidb/skeleton_influence_bounds.hpp"
#include "fiction/technology/sidb_defect_surface.hpp"
#include "fiction/technology/sidb_on_the_fly_mini_gate_library.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
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
    uint64_t                                num_trials                   = 500;
    double                                  selectivity                  = 0.5;
    uint64_t                                num_trials_for_global_scope  = 20;
    double                                  selectivity_for_global_scope = 0.8;
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

        std::unordered_map<mockturtle::node<GateLyt>, sidb_on_the_fly_mini_gate_library::designed_fcn_gates>
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
                                    *gate_lyt, t, params.sidb_on_the_fly_gate_library_parameters);
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

        if (!prune_gate_designs_by_gate_connections(*gate_lyt, operational_gate_designs) ||
            !prune_gate_designs_at_global_level(*gate_lyt, operational_gate_designs, lyt) ||
            !look_for_operational_circuit_exhaustively(*gate_lyt, operational_gate_designs, lyt))
        {
            std::cout << "\n\nFAILURE: NO OPERATIONAL CIRCUIT COULD BE GENERATED" << std::endl;
        }
        else
        {
            std::cout << "\n\nSUCCESS! GENERATED OPERATIONAL CIRCUIT:" << std::endl;
            print_layout(lyt);
        }

        stats.gate_layout = std::optional{gate_lyt};

        sidb_defect_surface<CellLyt> sidbs_and_defects{lyt};

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
    advanced_circuit_design_stats<GateLyt>& stats;

    using tt = kitty::dynamic_truth_table;

    [[nodiscard]] static std::vector<tt> determine_truth_table_of_gate_connection(
        const GateLyt& gate_lyt, const mockturtle::node<GateLyt>& n, const std::vector<tt>& spec_n,
        const mockturtle::node<GateLyt>& connecting_n, const std::vector<tt>& spec_connecting_n) noexcept
    {
        const tile<GateLyt> t            = {gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y, 0};
        const tile<GateLyt> connecting_t = {gate_lyt.get_tile(connecting_n).x, gate_lyt.get_tile(connecting_n).y, 0};

        bool connecting_n_is_incoming = false;

        std::vector<tile<GateLyt>> incoming_data_flow       = gate_lyt.incoming_data_flow(t);
        std::vector<tile<GateLyt>> incoming_data_flow_above = gate_lyt.incoming_data_flow(gate_lyt.above(t));

        incoming_data_flow.insert(incoming_data_flow.end(), std::make_move_iterator(incoming_data_flow_above.begin()),
                                  std::make_move_iterator(incoming_data_flow_above.end()));

        for (const auto& incoming_t : incoming_data_flow)
        {
            if (connecting_t.x == incoming_t.x && connecting_t.y == incoming_t.y)
            {
                connecting_n_is_incoming = true;
                break;
            }
        }

        const std::unordered_map<mockturtle::node<GateLyt>, uint8_t> num_fan_in{
            {n, incoming_data_flow.size()},
            {connecting_n, gate_lyt.incoming_data_flow(connecting_t).size() +
                               gate_lyt.incoming_data_flow(gate_lyt.above(connecting_t)).size()}};

        const uint8_t num_fan_in_combined = num_fan_in.at(n) + num_fan_in.at(connecting_n) - 1;

        std::vector<tt> vars{};
        for (uint8_t i = 0; i < num_fan_in_combined; ++i)
        {
            tt var{num_fan_in_combined};
            kitty::create_nth_var(var, i);
            vars.push_back(std::move(var));
        }

        const mockturtle::node<GateLyt> upper_n = connecting_n_is_incoming ? connecting_n : n;
        const mockturtle::node<GateLyt> lower_n = connecting_n_is_incoming ? n : connecting_n;

        const tile<GateLyt> upper_t = gate_lyt.get_tile(upper_n);
        const tile<GateLyt> lower_t = gate_lyt.get_tile(lower_n);

        const std::vector<tt> upper_spec = connecting_n_is_incoming ? spec_connecting_n : spec_n;
        const std::vector<tt> lower_spec = connecting_n_is_incoming ? spec_n : spec_connecting_n;

        // right to left means NE to SW, left to right means NW to SE
        const bool right_to_left = upper_t == gate_lyt.north_east(lower_t);

        // the first inputs go to the upper tile
        std::vector<tt> input_to_upper_n{};
        for (uint8_t i = 0; i < num_fan_in.at(upper_n); ++i)
        {
            input_to_upper_n.push_back(std::move(vars[i]));
        }

        // if the upper tile has a two-output function, we collect this direct output---this is always the first output
        std::optional<tt> upper_n_direct_output{};

        if (gate_lyt.outgoing_data_flow(upper_t).size() + gate_lyt.outgoing_data_flow(gate_lyt.above(upper_t)).size() ==
            2)
        {
            // for R->L, the direct output is determined by the second truth table in the spec; for L->R, the first
            upper_n_direct_output.emplace(kitty::compose_truth_table<tt, tt>(
                right_to_left ? upper_spec.back() : upper_spec.front(), input_to_upper_n));
        }

        std::vector<tt> lower_n_outputs{};

        if (upper_n_direct_output.has_value())
        {
            lower_n_outputs.push_back(std::move(upper_n_direct_output.value()));
        }

        // compute the output of the upper tile going to the lower one; the remaining truth table in the spec is used
        const tt mid_tt = kitty::compose_truth_table<tt, tt>(right_to_left ? upper_spec.front() : upper_spec.back(),
                                                             input_to_upper_n);

        // the remaining input variables go to the lower tile
        std::vector<tt> input_to_lower_n{};

        for (uint8_t i = num_fan_in.at(upper_n); i < num_fan_in_combined; ++i)
        {
            input_to_lower_n.push_back(std::move(vars[i]));
        }

        // the output of the upper tile is the first input in case of L->R, otherwise (R->L), it is the second
        input_to_lower_n.emplace(right_to_left ? input_to_lower_n.cend() : input_to_lower_n.cbegin(), mid_tt);

        // use the spec of the lower tile to determine the lower tile outputs
        for (const tt& lower_tt : lower_spec)
        {
            lower_n_outputs.emplace_back(kitty::compose_truth_table<tt, tt>(lower_tt, input_to_lower_n));
        }

        return lower_n_outputs;
    }

    static void print_success_rate_histogram(const std::vector<std::pair<double, uint64_t>>& success_ratios,
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
        std::cout << "-------------|--------|------------|--------------------------------------------------------------\n";

        for (const auto& [rate, count] : count_by_success_rate)
        {
            constexpr size_t max_bar_length = 60;

            const bool        is_kept = rate >= selectivity_threshold;
            const std::string status  = is_kept ? "KEPT" : "PRUNED";

            // Bar length proportional to total count
            const auto bar_len = 1 + static_cast<size_t>(static_cast<double>(count) / static_cast<double>(total_gates) *
                                                         (max_bar_length - 1));

            std::string bar(bar_len, is_kept ? 'O' : 'X');

            std::cout << std::right << std::setw(11) << std::fixed << std::setprecision(1) << rate * 100 << "%"
                      << " | " << std::setw(6) << status << " | " << std::setw(10) << count << " | " << std::left
                      << std::setw(max_bar_length) << bar << '\n';
        }

        std::cout << std::endl;
    }

    [[nodiscard]] uint64_t determine_first_passing_gate_ix(
        const std::vector<std::pair<double, uint64_t>>& successful_trial_ratio_per_gate_implementation,
        const double                                    selectivity) const noexcept
    {
        uint64_t first_passing_gate_ix = static_cast<uint64_t>(
            selectivity * static_cast<double>(successful_trial_ratio_per_gate_implementation.size()));

        // walk down (i.e., decrementing success ratio) to the first gate implementation with a different success ratio
        for (; first_passing_gate_ix > 0 &&
               std::abs(successful_trial_ratio_per_gate_implementation.at(first_passing_gate_ix - 1).first -
                        successful_trial_ratio_per_gate_implementation.at(first_passing_gate_ix).first) <
                   std::numeric_limits<double>::epsilon();
             --first_passing_gate_ix)
        {}

        print_success_rate_histogram(successful_trial_ratio_per_gate_implementation, first_passing_gate_ix);

        return first_passing_gate_ix;
    }

    [[nodiscard]] std::vector<uint64_t> select_gate_implementations_by_successful_trial_ratio(
        std::vector<std::pair<double, uint64_t>>&& successful_trial_ratio_per_gate_implementation,
        const double                               selectivity) const noexcept
    {
        std::sort(successful_trial_ratio_per_gate_implementation.begin(),
                  successful_trial_ratio_per_gate_implementation.end(),
                  [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

        const uint64_t first_passing_gate_ix =
            determine_first_passing_gate_ix(successful_trial_ratio_per_gate_implementation, selectivity);

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
        std::unordered_map<mockturtle::node<GateLyt>, sidb_on_the_fly_mini_gate_library::designed_fcn_gates>&
            operational_gate_designs) const noexcept
    {

        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS BY GATE CONNECTIONS" << std::endl;

        std::unordered_map<mockturtle::node<GateLyt>, cell<CellLyt>> top_left_corner_abs{};
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<mockturtle::node<GateLyt>>>
            gate_connections_to_simulate{};
        std::unordered_map<mockturtle::node<GateLyt>,
                           std::unordered_map<mockturtle::node<GateLyt>,
                                              std::pair<std::vector<tt>, is_operational_params<cell<CellLyt>>>>>
            spec_and_operational_params_for_joint_simulation{};
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const tile<GateLyt> t = gate_lyt.get_tile(n);

                top_left_corner_abs.insert(
                    {n,
                     relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                                        CellLyt>(gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0})});

                std::vector<mockturtle::node<GateLyt>> connecting_to_n{};

                for (const auto& in_t : gate_lyt.incoming_data_flow(t))
                {
                    if (!gate_lyt.is_pi(gate_lyt.get_node(in_t)))
                    {
                        connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(in_t)));
                    }
                }

                for (const auto& out_t : gate_lyt.outgoing_data_flow(t))
                {
                    if (!gate_lyt.is_po(gate_lyt.get_node(out_t)))
                    {
                        connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
                    }
                }

                if (gate_lyt.is_buf(n))
                {
                    for (const auto& in_t : gate_lyt.incoming_data_flow(gate_lyt.above(t)))
                    {
                        if (!gate_lyt.is_pi(gate_lyt.get_node(in_t)))
                        {
                            connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(in_t)));
                        }
                    }

                    for (const auto& out_t : gate_lyt.outgoing_data_flow(gate_lyt.above(t)))
                    {
                        if (!gate_lyt.is_po(gate_lyt.get_node(out_t)))
                        {
                            connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
                        }
                    }
                }

                for (const auto& connecting_n : connecting_to_n)
                {
                    is_operational_params operational_params =
                        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params;

                    operational_params.termination_cond =
                        is_operational_params<cell<CellLyt>>::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

                    operational_params.cc_map = skeleton_influence_bounds<CellLyt, SkeletonGateLibrary, GateLyt>(
                        gate_lyt, {t, gate_lyt.get_tile(connecting_n)},
                        skeleton_influence_bounds_params<cell<CellLyt>>{
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                                .simulation_parameters,
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.canvas,
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params_complex_gates.canvas,
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                                .input_bdl_iterator_params.bdl_wire_params,
                            true});

                    spec_and_operational_params_for_joint_simulation[n].insert(
                        {connecting_n,
                         {determine_truth_table_of_gate_connection(gate_lyt, n, operational_gate_designs.at(n).function,
                                                                   connecting_n,
                                                                   operational_gate_designs.at(connecting_n).function),
                          std::move(operational_params)}});
                }

                gate_connections_to_simulate.insert({n, std::move(connecting_to_n)});
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

                    std::cout << fmt::format("\nStarting pruning for tile {}", gate_lyt.get_tile(n)) << std::endl;
                    std::cout << fmt::format(
                                     "Performing {} trials for {} gate connection{} for {} gate implementations",
                                     params.num_trials, gate_connections_to_simulate.at(n).size(),
                                     gate_connections_to_simulate.at(n).size() > 1 ? "s" : "",
                                     operational_gate_designs.at(n).designed_gates.size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};
                    std::mutex                               mutex_to_protect_successful_trial_ratios;

                    // const uint64_t num_threads = 1;
                    const uint64_t num_threads = std::min(static_cast<uint64_t>(std::thread::hardware_concurrency()),
                                                          operational_gate_designs.at(n).designed_gates.size());
                    //  const uint64_t num_threads = std::min(uint64_t{20},
                    //  operational_gate_designs.at(n).designed_gates.size());

                    const uint64_t chunk_size =
                        (operational_gate_designs.at(n).designed_gates.size() + num_threads - 1) /
                        num_threads;  // Ceiling division

                    std::vector<std::thread> threads{};
                    threads.reserve(num_threads);

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{
                        static_cast<uint32_t>(
                            std::min(chunk_size, operational_gate_designs.at(n).designed_gates.size())),
                        "[i] Determining successful trial ratio for tile " +
                            fmt::format("({},{})", gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y) + ": |{0}|"};
#endif

                    for (uint64_t i = 0; i < num_threads; ++i)
                    {
                        threads.emplace_back(
                            [&successful_trial_ratio_per_gate_implementation, i, chunk_size, &operational_gate_designs,
                             &n, &top_left_corner_abs, &gate_lyt, this, &gate_connections_to_simulate,
#if (PROGRESS_BARS)

                             &bar,
#endif
                             &spec_and_operational_params_for_joint_simulation,
                             &mutex_to_protect_successful_trial_ratios]
                            {
                                const uint64_t start_index = i * chunk_size;
                                const uint64_t end_index   = std::min(
                                    start_index + chunk_size, operational_gate_designs.at(n).designed_gates.size());

                                for (uint64_t j = start_index; j < end_index; ++j)
                                {
                                    // std::cout << "gate design ix: " << j << std::endl;

                                    CellLyt cell_lyt{};

                                    // select the first gate implementation for n
                                    assign_gate<CellLyt, GateLibrary, GateLyt>(
                                        cell_lyt, top_left_corner_abs.at(n),
                                        *std::next(operational_gate_designs.at(n).designed_gates.cbegin(),
                                                   static_cast<int64_t>(j)),
                                        gate_lyt, n);

                                    double successful_trials = 0;

                                    // uint64_t gate_connection_ix = 0;
                                    for (const auto& connecting_n : gate_connections_to_simulate.at(n))
                                    {
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
                                                0, operational_gate_designs.at(connecting_n).designed_gates.size() - 1};

                                            // select a random gate implementation for the tile that connects as
                                            // input to n
                                            assign_gate<CellLyt, GateLibrary, GateLyt>(
                                                cell_lyt_clone, top_left_corner_abs.at(connecting_n),
                                                operational_gate_designs.at(connecting_n)
                                                    .designed_gates.at(distrib(gen)),
                                                gate_lyt, connecting_n);

                                            // mockturtle::stopwatch<>::duration time_counter{};
                                            // {
                                            //     mockturtle::stopwatch stop{time_counter};

                                            const operational_assessment<
                                                CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
                                                is_operational<CellLyt, tt, local_external_potential_type::BOUNDED>(
                                                    cell_lyt_clone,
                                                    spec_and_operational_params_for_joint_simulation.at(n)
                                                        .at(connecting_n)
                                                        .first,
                                                    spec_and_operational_params_for_joint_simulation.at(n)
                                                        .at(connecting_n)
                                                        .second);

                                            assert(op_assessment.assessment_per_input.has_value() &&
                                                   "ALL_COMBINATIONS_ENUMERATED is not set.");

                                            if (op_assessment.status == operational_status::OPERATIONAL)
                                            {
                                                double valid_groundstate_probabilty_sum_for_all_inputs = 0.0;

                                                for (const typename operational_assessment<
                                                         CellLyt, local_external_potential_type::BOUNDED>::
                                                         operational_assessment_for_input& op_assessment_for_input :
                                                     *op_assessment.assessment_per_input)
                                                {
                                                    valid_groundstate_probabilty_sum_for_all_inputs +=
                                                        op_assessment_for_input.valid_ground_state_probability;
                                                }

                                                successful_trials +=
                                                    valid_groundstate_probabilty_sum_for_all_inputs /
                                                    static_cast<double>(op_assessment.assessment_per_input->size());
                                            }
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
                                            //     std::cout << fmt::format("time taken: {:.3f} s | avg time: {:.3f} s",
                                            //                              time_counter.count() / 1e9, avg_time / 1e9)
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
                                                            gate_connections_to_simulate.at(n).size());

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
                        std::move(successful_trial_ratio_per_gate_implementation), params.selectivity);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() <
                             operational_gate_designs.at(n).designed_gates.size())
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
                            std::move(operational_gate_designs[n].designed_gates[selected_gate_implementation_index]));
                    }

                    operational_gate_designs[n].designed_gates = std::move(selected_gate_implementations);
                });
        }

        return !exit_by_failure;
    }

    bool prune_gate_designs_by_two_gate_connections(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, sidb_on_the_fly_mini_gate_library::designed_fcn_gates>&
            operational_gate_designs) const noexcept
    {

        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS BY TWO GATE CONNECTIONS" << std::endl;

        std::unordered_map<mockturtle::node<GateLyt>, cell<CellLyt>> top_left_corner_abs{};
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<mockturtle::node<GateLyt>>>
            gate_connections_to_simulate{};
        std::unordered_map<mockturtle::node<GateLyt>,
                           std::unordered_map<mockturtle::node<GateLyt>,
                                              std::pair<std::vector<tt>, is_operational_params<cell<CellLyt>>>>>
            spec_and_operational_params_for_joint_simulation{};
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const tile<GateLyt> t = gate_lyt.get_tile(n);

                top_left_corner_abs.insert(
                    {n,
                     relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                                        CellLyt>(gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0})});

                std::vector<mockturtle::node<GateLyt>> connecting_to_n{};

                for (const auto& in_t : gate_lyt.incoming_data_flow(t))
                {
                    if (!gate_lyt.is_pi(gate_lyt.get_node(in_t)))
                    {
                        connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(in_t)));
                    }
                }

                for (const auto& out_t : gate_lyt.outgoing_data_flow(t))
                {
                    if (!gate_lyt.is_po(gate_lyt.get_node(out_t)))
                    {
                        connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
                    }
                }

                if (gate_lyt.is_buf(n))
                {
                    for (const auto& in_t : gate_lyt.incoming_data_flow(gate_lyt.above(t)))
                    {
                        if (!gate_lyt.is_pi(gate_lyt.get_node(in_t)))
                        {
                            connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(in_t)));
                        }
                    }

                    for (const auto& out_t : gate_lyt.outgoing_data_flow(gate_lyt.above(t)))
                    {
                        if (!gate_lyt.is_po(gate_lyt.get_node(out_t)))
                        {
                            connecting_to_n.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
                        }
                    }
                }

                for (const auto& connecting_n : connecting_to_n)
                {
                    is_operational_params operational_params =
                        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params;

                    operational_params.termination_cond =
                        is_operational_params<cell<CellLyt>>::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

                    operational_params.cc_map = skeleton_influence_bounds<CellLyt, SkeletonGateLibrary, GateLyt>(
                        gate_lyt, {t, gate_lyt.get_tile(connecting_n)},
                        skeleton_influence_bounds_params<cell<CellLyt>>{
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                                .simulation_parameters,
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.canvas,
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params_complex_gates.canvas,
                            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params
                                .input_bdl_iterator_params.bdl_wire_params,
                            true});

                    spec_and_operational_params_for_joint_simulation[n].insert(
                        {connecting_n,
                         {determine_truth_table_of_gate_connection(gate_lyt, n, operational_gate_designs.at(n).function,
                                                                   connecting_n,
                                                                   operational_gate_designs.at(connecting_n).function),
                          std::move(operational_params)}});
                }

                gate_connections_to_simulate.insert({n, std::move(connecting_to_n)});
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

                    std::cout << fmt::format("\nStarting pruning for tile {}", gate_lyt.get_tile(n)) << std::endl;
                    std::cout << fmt::format(
                                     "Performing {} trials for {} gate connection{} for {} gate implementations",
                                     params.num_trials, gate_connections_to_simulate.at(n).size(),
                                     gate_connections_to_simulate.at(n).size() > 1 ? "s" : "",
                                     operational_gate_designs.at(n).designed_gates.size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};
                    std::mutex                               mutex_to_protect_successful_trial_ratios;

                    // const uint64_t num_threads = 1;
                    const uint64_t num_threads = std::min(static_cast<uint64_t>(std::thread::hardware_concurrency()),
                                                          operational_gate_designs.at(n).designed_gates.size());
                    //  const uint64_t num_threads = std::min(uint64_t{20},
                    //  operational_gate_designs.at(n).designed_gates.size());

                    const uint64_t chunk_size =
                        (operational_gate_designs.at(n).designed_gates.size() + num_threads - 1) /
                        num_threads;  // Ceiling division

                    std::vector<std::thread> threads{};
                    threads.reserve(num_threads);

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{
                        static_cast<uint32_t>(
                            std::min(chunk_size, operational_gate_designs.at(n).designed_gates.size())),
                        "[i] Determining successful trial ratio for tile " +
                            fmt::format("({},{})", gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y) + ": |{0}|\t"};
#endif

                    for (uint64_t i = 0; i < num_threads; ++i)
                    {
                        threads.emplace_back(
                            [&successful_trial_ratio_per_gate_implementation, i, chunk_size, &operational_gate_designs,
                             &n, &top_left_corner_abs, &gate_lyt, this, &gate_connections_to_simulate,
#if (PROGRESS_BARS)

                             &bar,
#endif
                             &spec_and_operational_params_for_joint_simulation,
                             &mutex_to_protect_successful_trial_ratios]
                            {
                                const uint64_t start_index = i * chunk_size;
                                const uint64_t end_index   = std::min(
                                    start_index + chunk_size, operational_gate_designs.at(n).designed_gates.size());

                                for (uint64_t j = start_index; j < end_index; ++j)
                                {
                                    // std::cout << "gate design ix: " << j << std::endl;

                                    CellLyt cell_lyt{};

                                    // select the first gate implementation for n
                                    assign_gate<CellLyt, GateLibrary, GateLyt>(
                                        cell_lyt, top_left_corner_abs.at(n),
                                        *std::next(operational_gate_designs.at(n).designed_gates.cbegin(),
                                                   static_cast<int64_t>(j)),
                                        gate_lyt, n);

                                    double successful_trials = 0;

                                    // uint64_t gate_connection_ix = 0;
                                    for (const auto& connecting_n : gate_connections_to_simulate.at(n))
                                    {
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
                                                0, operational_gate_designs.at(connecting_n).designed_gates.size() - 1};

                                            // select a random gate implementation for the tile that connects as
                                            // input to n
                                            assign_gate<CellLyt, GateLibrary, GateLyt>(
                                                cell_lyt_clone, top_left_corner_abs.at(connecting_n),
                                                operational_gate_designs.at(connecting_n)
                                                    .designed_gates.at(distrib(gen)),
                                                gate_lyt, connecting_n);

                                            // mockturtle::stopwatch<>::duration time_counter{};
                                            // {
                                            //     mockturtle::stopwatch stop{time_counter};

                                            const operational_assessment<
                                                CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
                                                is_operational<CellLyt, tt, local_external_potential_type::BOUNDED>(
                                                    cell_lyt_clone,
                                                    spec_and_operational_params_for_joint_simulation.at(n)
                                                        .at(connecting_n)
                                                        .first,
                                                    spec_and_operational_params_for_joint_simulation.at(n)
                                                        .at(connecting_n)
                                                        .second);

                                            assert(op_assessment.assessment_per_input.has_value() &&
                                                   "ALL_COMBINATIONS_ENUMERATED is not set.");

                                            if (op_assessment.status == operational_status::OPERATIONAL)
                                            {
                                                double valid_groundstate_probabilty_sum_for_all_inputs = 0.0;

                                                for (const typename operational_assessment<
                                                         CellLyt, local_external_potential_type::BOUNDED>::
                                                         operational_assessment_for_input& op_assessment_for_input :
                                                     *op_assessment.assessment_per_input)
                                                {
                                                    valid_groundstate_probabilty_sum_for_all_inputs +=
                                                        op_assessment_for_input.valid_ground_state_probability;
                                                }

                                                successful_trials +=
                                                    valid_groundstate_probabilty_sum_for_all_inputs /
                                                    static_cast<double>(op_assessment.assessment_per_input->size());
                                            }
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
                                            //     std::cout << fmt::format("time taken: {:.3f} s | avg time: {:.3f} s",
                                            //                              time_counter.count() / 1e9, avg_time / 1e9)
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
                                                            gate_connections_to_simulate.at(n).size());

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
                        std::move(successful_trial_ratio_per_gate_implementation), params.selectivity);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() <
                             operational_gate_designs.at(n).designed_gates.size())
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
                            std::move(operational_gate_designs[n].designed_gates[selected_gate_implementation_index]));
                    }

                    operational_gate_designs[n].designed_gates = std::move(selected_gate_implementations);
                });
        }

        return !exit_by_failure;
    }

    bool prune_gate_designs_at_global_level(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, sidb_on_the_fly_mini_gate_library::designed_fcn_gates>&
                 operational_gate_designs,
        CellLyt& lyt) noexcept
    {
        std::cout << "\nSTARTING TO PRUNE GATE DESIGNS WITH GLOBAL SCOPE" << std::endl;

        bool operational_circuit_found = false;

        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.termination_cond =
            decltype(params.sidb_on_the_fly_gate_library_parameters.design_gate_params
                         .operational_params)::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        // params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.print = true;

        const uint64_t num_input_combinations = 1 << gate_lyt.num_pis();

        std::unordered_map<mockturtle::node<GateLyt>, cell<CellLyt>>         top_left_corner_abs{};
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_lyt, n))
                {
                    return;
                }

                const tile<GateLyt> t = gate_lyt.get_tile(n);

                top_left_corner_abs.insert(
                    {n,
                     relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                                        CellLyt>(gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0})});
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
                    if (operational_circuit_found || exit_by_failure || skip_physical_design_for_node(gate_lyt, n) ||
                        operational_gate_designs.at(n).designed_gates.size() == 1)
                    {
                        return;
                    }

                    std::cout << fmt::format("\nStarting pruning for tile {}\n", gate_lyt.get_tile(n));
                    std::cout << fmt::format("Performing {} trials for {} gate implementations",
                                             params.num_trials_for_global_scope,
                                             operational_gate_designs.at(n).designed_gates.size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{
                        static_cast<uint32_t>(operational_gate_designs.at(n).designed_gates.size()),
                        "[i] Determining successful trial ratio for tile " +
                            fmt::format("({},{})", gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y) + ": |{0}|\t"};
#endif

                    for (uint64_t j = 0;
                         !operational_circuit_found && j < operational_gate_designs.at(n).designed_gates.size(); ++j)
                    {
                        CellLyt cell_lyt{};

                        // select the first gate implementation for n
                        assign_gate<CellLyt, GateLibrary, GateLyt>(
                            cell_lyt, top_left_corner_abs.at(n),
                            *std::next(operational_gate_designs.at(n).designed_gates.cbegin(), static_cast<int64_t>(j)),
                            gate_lyt, n);

                        uint64_t successful_trials = 0;

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
                                        0, operational_gate_designs.at(other_n).designed_gates.size() - 1};

                                    // select a random gate implementation for the tile that connects as
                                    // input to n
                                    assign_gate<CellLyt, GateLibrary, GateLyt>(
                                        cell_lyt_clone, top_left_corner_abs.at(other_n),
                                        operational_gate_designs.at(other_n).designed_gates.at(distrib(gen)), gate_lyt,
                                        other_n);
                                });

                            const operational_assessment<CellLyt>& op_assessment_for_all_inputs = is_operational(
                                cell_lyt_clone, params.spec,
                                params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params);

                            assert(op_assessment_for_all_inputs.assessment_per_input.has_value() &&
                                   "Operational assessment termination condition was not set to "
                                   "ALL_INPUT_COMBINATIONS_ENUMERATED");

                            uint64_t successful_input_combinations = 0;

                            for (const typename operational_assessment<CellLyt>::operational_assessment_for_input&
                                     op_assessment : op_assessment_for_all_inputs.assessment_per_input.value())
                            {
                                if (op_assessment.status == operational_status::OPERATIONAL)
                                {
                                    successful_input_combinations++;
                                }
                            }

                            if (successful_input_combinations == num_input_combinations)
                            {
                                // all input combinations are operational---operational circuit found

                                lyt = cell_lyt_clone;

                                operational_circuit_found = true;
                            }
                            else
                            {
                                successful_trials += successful_input_combinations;
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
                            static_cast<double>(successful_trials) /
                            static_cast<double>(params.num_trials_for_global_scope * num_input_combinations);

                        successful_trial_ratio_per_gate_implementation.emplace_back(successful_trial_ratio, j);
                    }

                    if (operational_circuit_found)
                    {
                        return;
                    }

                    selected_gate_implementation_indices[n] = select_gate_implementations_by_successful_trial_ratio(
                        std::move(successful_trial_ratio_per_gate_implementation), params.selectivity_for_global_scope);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE " << gate_lyt.get_tile(n)
                                  << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() <
                             operational_gate_designs.at(n).designed_gates.size())
                    {
                        big_fixpoint = false;
                    }
                });

            if (operational_circuit_found)
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
                            std::move(operational_gate_designs[n].designed_gates[selected_gate_implementation_index]));
                    }

                    operational_gate_designs[n].designed_gates = std::move(selected_gate_implementations);
                });
        }

        return false;
    }
    /**
     * todo
     */
    bool look_for_operational_circuit_exhaustively(
        const GateLyt& gate_lyt,
        std::unordered_map<mockturtle::node<GateLyt>, sidb_on_the_fly_mini_gate_library::designed_fcn_gates>&
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
                assign_gate<CellLyt, GateLibrary, GateLyt>(
                    operational_circuit_candidate,
                    relative_to_absolute_cell_position<GateLibrary::gate_x_size(), GateLibrary::gate_y_size(), GateLyt,
                                                       CellLyt>(gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0}),
                    op_gate_designs_for_gate.designed_gates.at(indices.at(i)), gate_lyt, n);
            }

            std::cout << "trying combination: ";
            for (uint64_t i = 0; i < operational_gate_designs.size(); i++)
            {
                std::cout << indices.at(i) << " ";
            }
            std::cout << std::endl;

            if (is_operational(operational_circuit_candidate, params.spec,
                               params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params)
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
                    std::next(operational_gate_designs.cbegin(), static_cast<int64_t>(i))->second.designed_gates.size())
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

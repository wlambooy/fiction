//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/algorithms/simulation/sidb/is_circuit_operational.hpp"
#include "fiction/algorithms/simulation/sidb/is_operational.hpp"
#include "fiction/technology/sidb_on_the_fly_gate_library.hpp"
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
#include <mutex>
#include <optional>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
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
     * Parameters for the *exact* placement and routing algorithm.
     */
    exact_physical_design_params exact_design_parameters = {};
    /**
     * This struct holds parameters to design SiDB gates.
     */
    design_sidb_gates_params<CellLyt> design_gate_params{};
    /** This variable specifies the radius in nanometers around the center of the hexagon where atomic defects are
     * incorporated into the gate design. (unit: nm)
     */

    std::optional<CellLyt> defect_surface{};
    double                 influence_radius_charged_defects = 15;

    uint64_t num_trials           = 500;
    double   quantization_factor  = 0.075;
    double   selectivity          = 0.5;
    double   success_rate_ceiling = 0.95;

    uint64_t available_threads = std::thread::hardware_concurrency();
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

template <typename Ntk, typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
class advanced_circuit_design_impl
{
  public:
    advanced_circuit_design_impl(const Ntk& ntk, advanced_circuit_design_params<CellLyt>& design_params,
                                 const GateLyt& tiling, advanced_circuit_design_stats<GateLyt>& st) :
            lattice_tiling{tiling},
            network{ntk},
            params{design_params},
            stats{st}
    {
        operational_params.simulation_parameters = params.design_gate_params.operational_params.simulation_parameters;
        operational_params.input_bdl_iterator_params =
            params.design_gate_params.operational_params.input_bdl_iterator_params;
    }

    [[nodiscard]] std::optional<CellLyt> design_sidb_layout(const uint64_t num_gates_to_design)
    {
        circuit.emplace(*stats.gate_layout,
                        params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params);

        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ON_FIRST_NON_OPERATIONAL;

        // initialize
        collect_initial_gate_designs();

        operational_params.termination_cond =
            is_circuit_operational_params::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED;

        const auto quit_with_error = []
        {
            std::cout << "\nGate design was unsuccessful." << std::endl;
            return std::nullopt;
        };

        while (++circuit_design_level < num_gates_to_design - 1)
        {
            // prune by assessing gate design combinations for increasingly large sets of connected gates
            if (std::get<bool>(prune_gate_designs(
                    circuit_design_level, std::string_view{std::to_string(circuit_design_level) + " GATE CONNECTION" +
                                                           (circuit_design_level > 1 ? "S" : "")})))
            {
                return quit_with_error();
            }
        }

        // prune at the global level (all gates are considered together)
        if (const std::variant<bool, CellLyt>& maybe_lyt = prune_gate_designs(circuit_design_level, "GLOBAL");
            std::holds_alternative<CellLyt>(maybe_lyt))
        {
            return std::get<CellLyt>(maybe_lyt);
        }
        else if (std::get<bool>(maybe_lyt))
        {
            return quit_with_error();
        }

        return exhaustively_enumerate_gate_design_combinations();
    }

    [[nodiscard]] std::optional<CellLyt> design_circuit()
    {
        const mockturtle::stopwatch stop{stats.time_total};

        stats.gate_layout = exact<GateLyt>(network, params.exact_design_parameters, &stats.exact_stats);

        if (!stats.gate_layout.has_value() || bounding_box_2d<GateLyt>{*stats.gate_layout}.get_x_size() > 5 ||
            bounding_box_2d<GateLyt>{*stats.gate_layout}.get_y_size() > 5 ||
            bounding_box_2d<GateLyt>{*stats.gate_layout}.get_x_size() < 1 ||
            bounding_box_2d<GateLyt>{*stats.gate_layout}.get_y_size() < 1)
        {
            // P&R was unsuccessful
            std::cout << "UNSUCCESS" << std::endl;
            return std::nullopt;
        }

        print_skeleton_gate_layout(*stats.gate_layout);

        uint64_t number_of_gates_to_design = 0;

        stats.gate_layout->foreach_node(
            [&](const mockturtle::node<GateLyt>& n)
            {
                if (!skip_physical_design_for_node(*stats.gate_layout, n))
                {
                    number_of_gates_to_design++;
                }
            });

        if (const std::optional<CellLyt>& maybe_sidb_lyt = design_sidb_layout(number_of_gates_to_design);
            maybe_sidb_lyt.has_value())
        {
            std::cout << "\n\nSUCCESS! GENERATED OPERATIONAL CIRCUIT:" << std::endl;
            print_layout(maybe_sidb_lyt.value());

            /// todo not necessary,,,, right?
            // // // add defects to the circuit.
            // // params.sidb_on_the_fly_gate_library_parameters.defect_surface.foreach_sidb_defect(
            // //     [&sidbs_and_defects](const auto& defect)
            // //     { sidbs_and_defects.assign_sidb_defect(defect.first, defect.second); });

            return maybe_sidb_lyt.value();
        }

        std::cout << "\n\nFAILURE: NO OPERATIONAL CIRCUIT COULD BE GENERATED" << std::endl;

        return std::nullopt;
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

    is_circuit_operational_params operational_params{};

    using gate_designs_per_node =
        std::unordered_map<mockturtle::node<GateLyt>, std::vector<typename SkeletonGateLibrary::fcn_gate>>;

    gate_designs_per_node gate_designs{};

    uint64_t circuit_design_level = 0;

    void print_skeleton_gate_layout(const GateLyt& gate_layout) const noexcept
    {
        CellLyt lyt = apply_gate_library<CellLyt, SkeletonGateLibrary, GateLyt>(gate_layout);

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(gate_layout, n))
                {
                    return;
                }

                const auto& t = gate_layout.get_tile(n);

                for (const cell<CellLyt>& relative_c : all_coordinates_in_spanned_area(
                         params.design_gate_params.canvas.first, params.design_gate_params.canvas.second))
                {
                    const cell<CellLyt> absolute_c =
                        relative_to_absolute_cell_position<SkeletonGateLibrary::gate_x_size(),
                                                           SkeletonGateLibrary::gate_y_size(), GateLyt, CellLyt>(
                            gate_layout, t, relative_c);
                    lyt.assign_cell_type(absolute_c, sidb_technology::cell_type::LOGIC);
                }
            });

        std::cout << "Skeleton looks like:" << std::endl;
        print_layout(lyt);
        std::cout << std::endl;
    }

    void collect_initial_gate_designs()
    {
        const sidb_on_the_fly_gate_library_params<CellLyt> on_the_fly_params{params.design_gate_params,
                                                                             params.influence_radius_charged_defects};

        try
        {
            stats.gate_layout->foreach_node(
                [&, this](const auto& n, [[maybe_unused]] auto i)
                {
                    if (!skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        gate_designs[n] = sidb_on_the_fly_gate_library<SkeletonGateLibrary::gate_x_size(),
                                                                       SkeletonGateLibrary::gate_y_size()>::
                            template set_up_gates<GateLyt, CellLyt, sidb_on_the_fly_gate_library_params<CellLyt>,
                                                  local_external_potential_type::BOUNDED, SkeletonGateLibrary>(
                                *stats.gate_layout, stats.gate_layout->get_tile(n), on_the_fly_params,
                                params.defect_surface, std::make_optional(*circuit), operational_params);
                    }
                });
        }

        catch (const gate_design_exception<tt, GateLyt>& e)
        {
            throw unsuccessful_gate_design_error("Gate design was unsuccessful");
        }
    }

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
            }
        }

        for (const auto& out_t : gate_lyt.outgoing_data_flow(t))
        {
            if (!gate_lyt.is_po(gate_lyt.get_node(out_t)) &&
                (!maybe_do_not_collect.has_value() || gate_lyt.below(out_t) != *maybe_do_not_collect))
            {
                connecting_to_t.emplace_back(gate_lyt.get_node(gate_lyt.below(out_t)));
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
                                   const double quantization_factor, const double success_rate_ceiling) noexcept
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
            success_rate = std::min(success_rate_ceiling, std::round(success_rate * scale) / scale);
        }
    }

    [[nodiscard]] static uint64_t determine_first_passing_gate_ix(
        const std::vector<std::pair<double, uint64_t>>& successful_trial_ratio_per_gate_implementation,
        const double                                    selectivity) noexcept
    {
        const auto first_passing_gate_ix = static_cast<int64_t>(
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
        const double quantization_factor, const double selectivity, const double success_rate_ceiling) const noexcept
    {
        std::sort(successful_trial_ratio_per_gate_implementation.begin(),
                  successful_trial_ratio_per_gate_implementation.end(),
                  [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

        apply_quantization(successful_trial_ratio_per_gate_implementation, quantization_factor, success_rate_ceiling);

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

    struct VectorHash
    {
        template <typename T>
        std::size_t operator()(const std::vector<T>& v) const
        {
            std::size_t seed = v.size();
            for (const auto& i : v)
            {
                seed ^= std::hash<T>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
    using gate_lyt_window_map =
        std::unordered_map<std::vector<mockturtle::node<GateLyt>>,
                           sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>, VectorHash>;

    void build_subcircuits(const std::vector<mockturtle::node<GateLyt>>   root_connections,
                           std::vector<mockturtle::node<GateLyt>>&        path,
                           std::unordered_set<mockturtle::node<GateLyt>>& visited, const size_t depth,
                           const size_t max_depth, gate_lyt_window_map& result) const
    {
        // todo: check out symmetry

        if (depth == max_depth)
        {
            // Get tiles corresponding to path nodes
            std::vector<tile<GateLyt>>                    tiles{};
            std::unordered_set<mockturtle::node<GateLyt>> nodes{};
            for (const auto& node : path)
            {
                tiles.push_back(stats.gate_layout->get_tile(node));
                nodes.insert(node);
            }

            // filter duplicates (i.e.: RAB, RBA)
            for (const auto& [other_path, _] : result)
            {
                if (path.front() == other_path.front() &&
                    std::all_of(std::next(other_path.cbegin(), 1), other_path.cend(),
                                [&](const auto& n) { return nodes.count(n) > 0; }))
                {
                    return;
                }
            }

            result.insert(
                {path, sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                           std::cref(*circuit), tiles,
                           params.design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params}});

            return;
        }

        const tile<GateLyt> current_tile = stats.gate_layout->get_tile(path.back());

        std::vector<mockturtle::node<GateLyt>> connections = root_connections;
        collect_connecting_nodes(*stats.gate_layout, current_tile, connections);

        for (const auto& next_node : connections)
        {
            if (visited.count(next_node) > 0)
            {
                continue;
            }

            visited.insert(next_node);
            path.push_back(next_node);

            build_subcircuits(root_connections, path, visited, depth + 1, max_depth, result);

            path.pop_back();
            visited.erase(next_node);
        }
    }
    void build_subcircuits_from_root(const mockturtle::node<GateLyt> root, const size_t max_depth,
                                     gate_lyt_window_map& result) const
    {
        const tile<GateLyt> current_tile = stats.gate_layout->get_tile(root);

        std::vector<mockturtle::node<GateLyt>> root_connections;
        collect_connecting_nodes(*stats.gate_layout, current_tile, root_connections);

        std::vector<mockturtle::node<GateLyt>>        path    = {root};
        std::unordered_set<mockturtle::node<GateLyt>> visited = {root};

        build_subcircuits(root_connections, path, visited, 0, max_depth, result);
    }
    /**
     * todo
     */
    std::variant<bool, CellLyt> prune_gate_designs(const uint64_t level, const std::string_view& level_str) noexcept
    {
        std::cout << "\n\nSTARTING TO PRUNE GATE DESIGNS\tLEVEL: " << level_str << std::endl;

        bool global_pruning = false;

        if (level_str == "GLOBAL")
        {
            global_pruning = true;
        }

        std::optional<CellLyt> maybe_lyt{};

        const uint64_t num_input_combinations = 1 << stats.gate_layout->num_pis();

        // todo
        const uint64_t num_trials           = params.num_trials / level;
        const double   quantization_factor  = params.quantization_factor * std::pow(1.2, level - 1);
        const double   selectivity          = params.selectivity / std::pow(3.0, level - 1);
        const double   success_rate_ceiling = params.success_rate_ceiling + static_cast<double>(level - 1) * 0.01;

        std::cout << "\nNumber of trials: " << num_trials << std::endl;
        std::cout << "Quantization factor: " << fmt::format("{:.3f}", quantization_factor) << std::endl;
        std::cout << "Selectivity: " << fmt::format("{:.2f}", selectivity) << std::endl;
        std::cout << "Success_rate_ceiling: " << fmt::format("{:.2f}", success_rate_ceiling) << std::endl;

        gate_lyt_window_map gate_lyt_windows{};

        stats.gate_layout->foreach_node(
            [&](const auto& n)
            {
                if (skip_physical_design_for_node(*stats.gate_layout, n))
                {
                    return;
                }

                build_subcircuits_from_root(n, level, gate_lyt_windows);
            });

        std::unordered_map<mockturtle::node<GateLyt>, std::vector<uint64_t>> selected_gate_implementation_indices{};

        bool big_fixpoint = false;

        bool exit_by_failure = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            std::cout << "\nStarting main fixpoint iteration\n" << std::endl;

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

                    const tile<GateLyt>& t = stats.gate_layout->get_tile(n);

                    uint64_t number_of_layouts = 0;

                    for (const auto& [node_vec, gate_lyt_window] : gate_lyt_windows)
                    {
                        assert(node_vec.size() > 1 && "The connected nodes vector cannot be singleton");

                        if (n == node_vec.front())
                        {
                            number_of_layouts++;
                        }
                    }

                    std::cout << fmt::format("\nStarting pruning for tile {}", t) << std::endl;
                    std::cout << fmt::format("Performing {} trials for {} sub-circuit{} for {} gate implementations",
                                             num_trials, number_of_layouts, number_of_layouts > 1 ? "s" : "",
                                             gate_designs.at(n).size())
                              << std::endl;

                    std::vector<std::pair<double, uint64_t>> successful_trial_ratio_per_gate_implementation{};
                    std::mutex                               mutex_to_protect_successful_trial_ratios;

                    const uint64_t num_threads = std::min(params.available_threads, gate_designs.at(n).size());
                    const uint64_t chunk_size =
                        (gate_designs.at(n).size() + num_threads - 1) / num_threads;  // Ceiling division

                    std::vector<std::thread> threads{};
                    threads.reserve(num_threads);

#if (PROGRESS_BARS)
                    mockturtle::progress_bar bar{static_cast<uint32_t>(std::min(chunk_size, gate_designs.at(n).size())),
                                                 "[i] Determining successful trial ratio for tile " +
                                                     fmt::format("({},{})", t.x, t.y) + ": |{0}|"};
#endif

                    std::mutex lyt_mutex{};

                    for (uint64_t i = 0; i < num_threads; ++i)
                    {
                        threads.emplace_back(
                            [&successful_trial_ratio_per_gate_implementation, &t, i, chunk_size, &n, this,
#if (PROGRESS_BARS)
                             &bar,
#endif
                             &gate_lyt_windows, number_of_layouts, &mutex_to_protect_successful_trial_ratios,
                             num_input_combinations, global_pruning, num_trials, &maybe_lyt, &lyt_mutex]
                            {
                                const uint64_t start_index = i * chunk_size;
                                const uint64_t end_index =
                                    std::min(start_index + chunk_size, gate_designs.at(n).size());

                                for (uint64_t j = start_index; j < end_index; ++j)
                                {
                                    CellLyt cell_lyt{};

                                    // select the first gate implementation for n
                                    assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(
                                        cell_lyt, *std::next(gate_designs.at(n).cbegin(), static_cast<int64_t>(j)),
                                        *stats.gate_layout, t);

                                    double successful_trials = 0;

                                    for (const auto& [node_vec, gate_lyt_window] : gate_lyt_windows)
                                    {
                                        if (n != node_vec.front())
                                        {
                                            continue;
                                        }

                                        uint64_t current_trial = 0;

                                        while (current_trial < num_trials)
                                        {
                                            if (global_pruning)
                                            {
                                                std::lock_guard lock{lyt_mutex};

                                                if (maybe_lyt.has_value())
                                                {
                                                    return;  // quit if an operational circuit has already been found
                                                }
                                            }

                                            current_trial++;

                                            CellLyt cell_lyt_clone = cell_lyt.clone();

                                            // assign random gate design to other gates
                                            for (typename std::vector<mockturtle::node<GateLyt>>::const_iterator
                                                     node_vec_it = std::next(node_vec.cbegin(), 1);
                                                 node_vec_it != node_vec.cend(); ++node_vec_it)
                                            {
                                                std::random_device rd;   // a seed source for the random number engine
                                                std::mt19937 gen(rd());  // mersenne_twister_engine seeded with rd()
                                                std::uniform_int_distribution<uint64_t> distrib{
                                                    0, gate_designs.at(*node_vec_it).size() - 1};

                                                // select a random gate implementation
                                                assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(
                                                    cell_lyt_clone, gate_designs.at(*node_vec_it).at(distrib(gen)),
                                                    *stats.gate_layout, stats.gate_layout->get_tile(*node_vec_it));
                                            }

                                            // sub-circuit logic match assessment
                                            const circuit_operational_assessment<
                                                CellLyt, local_external_potential_type::BOUNDED>& op_assessment =
                                                is_circuit_operational<CellLyt, GateLyt,
                                                                       local_external_potential_type::BOUNDED,
                                                                       SkeletonGateLibrary>(
                                                    sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                                                        cell_lyt_clone, gate_lyt_window},
                                                    operational_params);

                                            assert(op_assessment.assessment_per_input.has_value() &&
                                                   "ALL_INPUT_COMBINATIONS_ENUMERATED is not set.");

                                            uint64_t successful_input_combinations  = 0;
                                            double   logic_match_sum_for_all_inputs = 0.0;

                                            // sum logic match over the different input combinations
                                            for (const typename circuit_operational_assessment<
                                                     CellLyt, local_external_potential_type::BOUNDED>::
                                                     operational_assessment_for_input& op_assessment_for_input :
                                                 *op_assessment.assessment_per_input)
                                            {
                                                if (global_pruning &&
                                                    op_assessment_for_input.status == operational_status::OPERATIONAL)
                                                {
                                                    successful_input_combinations++;
                                                }

                                                logic_match_sum_for_all_inputs += op_assessment_for_input.logic_match;
                                            }

                                            if (global_pruning &&
                                                successful_input_combinations == num_input_combinations)
                                            {
                                                // all input combinations are operational---operational circuit found

                                                std::lock_guard lock{lyt_mutex};

                                                if (!maybe_lyt.has_value())
                                                {
                                                    maybe_lyt.emplace(cell_lyt_clone);
                                                }

                                                return;
                                            }

                                            successful_trials +=
                                                logic_match_sum_for_all_inputs /
                                                static_cast<double>(op_assessment.assessment_per_input->size());
                                        }
                                    }

#if (PROGRESS_BARS)
                                    if (i == 0)
                                    {
                                        bar(j);
                                    }
#endif

                                    const double successful_trial_ratio =
                                        successful_trials / static_cast<double>(num_trials * number_of_layouts);

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

                    if (global_pruning)
                    {
                        if (maybe_lyt.has_value())
                        {
                            return;  // quit if an operational circuit has been found
                        }
                    }

                    selected_gate_implementation_indices[n] = select_gate_implementations_by_successful_trial_ratio(
                        std::move(successful_trial_ratio_per_gate_implementation), quantization_factor, selectivity,
                        success_rate_ceiling);

                    if (selected_gate_implementation_indices.at(n).empty())
                    {
                        std::cout << "ERROR: ALL GATE IMPLEMENTATIONS ARE PRUNED FOR TILE "
                                  << stats.gate_layout->get_tile(n) << std::endl;
                        exit_by_failure = true;
                    }
                    else if (selected_gate_implementation_indices.at(n).size() < gate_designs.at(n).size())
                    {
                        big_fixpoint = false;
                    }
                });

            std::cout << std::endl;

            if (global_pruning)
            {
                if (maybe_lyt.has_value())
                {
                    return maybe_lyt.value();  // return the operational circuit
                }
            }

            stats.gate_layout->foreach_node(
                [&](const auto& n)
                {
                    if (exit_by_failure || skip_physical_design_for_node(*stats.gate_layout, n))
                    {
                        return;
                    }

                    std::vector<typename SkeletonGateLibrary::fcn_gate> selected_gate_implementations{};

                    for (const uint64_t selected_gate_implementation_index : selected_gate_implementation_indices.at(n))
                    {
                        selected_gate_implementations.push_back(
                            std::move(gate_designs[n][selected_gate_implementation_index]));
                    }

                    gate_designs[n] = std::move(selected_gate_implementations);
                });
        }

        return exit_by_failure;
    }
    /**
     * todo
     */
    std::optional<CellLyt> exhaustively_enumerate_gate_design_combinations() const noexcept
    {
        std::cout << "\n\nLOOKING FOR OPERATIONAL CIRCUIT EXHAUSTIVELY" << std::endl;

        // operational_params.print = true;

        std::vector<uint64_t> indices(gate_designs.size(), 0);

        while (true)
        {
            CellLyt operational_circuit_candidate{};
            for (uint64_t i = 0; i < gate_designs.size(); i++)
            {
                const auto& [n, op_gate_designs_for_gate] = *std::next(gate_designs.cbegin(), static_cast<int64_t>(i));
                // select a random gate implementation for the tile that connects as input to n
                assign_gate<CellLyt, SkeletonGateLibrary, GateLyt>(operational_circuit_candidate,
                                                                   op_gate_designs_for_gate.at(indices.at(i)),
                                                                   *stats.gate_layout, stats.gate_layout->get_tile(n));
            }

            std::cout << "trying combination: ";
            for (uint64_t i = 0; i < gate_designs.size(); i++)
            {
                std::cout << indices.at(i) << " ";
            }
            std::cout << std::endl;

            if (is_circuit_operational(
                    sidb_cell_level_bdl_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{
                        operational_circuit_candidate,
                        sidb_bdl_sub_circuit<CellLyt, GateLyt, SkeletonGateLibrary>{*circuit}},
                    operational_params)
                    .status == operational_status::OPERATIONAL)
            {
                std::cout << "\n\nFINAL GENERATED CIRCUIT:" << std::endl;
                print_layout(operational_circuit_candidate);

                return operational_circuit_candidate;
            }

            // Increment indices like an odometer
            for (uint64_t i = 0; i < indices.size(); ++i)
            {
                if (++indices[i] < std::next(gate_designs.cbegin(), static_cast<int64_t>(i))->second.size())
                {
                    break;  // No carry needed
                }

                indices[i] = 0;  // Reset this index and carry over to the next

                if (i == indices.size() - 1)
                {
                    return std::nullopt;  // Stop when the last index overflows
                }
            }
        }

        return std::nullopt;
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
template <typename Ntk, typename CellLyt, typename GateLyt, typename SkeletonGateLibrary>
[[nodiscard]] std::optional<CellLyt> advanced_circuit_design(const Ntk& ntk, const GateLyt& lattice_tiling,
                                                             advanced_circuit_design_params<CellLyt>& params = {},
                                                             advanced_circuit_design_stats<GateLyt>*  stats  = nullptr)
{
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt is not a hexagonal");
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");

    advanced_circuit_design_stats<GateLyt> st{};

    detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt, SkeletonGateLibrary> p{ntk, params, lattice_tiling, st};

    const auto result = p.design_circuit();

    if (stats)
    {
        *stats = st;
    }

    return result;
}

}  // namespace fiction

#endif  // advanced_CIRCUIT_DESIGN_HPP

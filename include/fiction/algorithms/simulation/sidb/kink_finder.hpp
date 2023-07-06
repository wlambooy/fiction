//
// Created by Willem Lambooy on 15/05/2023.
//

#ifndef FICTION_KINK_FINDER_HPP
#define FICTION_KINK_FINDER_HPP

#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/composim.hpp"
#include "fiction/algorithms/simulation/sidb/quickexact.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_result.hpp"
#include "fiction/io/print_layout.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/technology/cell_ports.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/fcn_gate_library.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/hash.hpp"
#include "fiction/utils/layout_utils.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/hash.hpp>
#include <mockturtle/utils/stopwatch.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiction
{

class sidb_bestagon_library_wrapper
{
  public:
    using library = sidb_bestagon_library;

    explicit sidb_bestagon_library_wrapper() = delete;

    // index of input 0 is before input 1
    static inline const std::array<std::pair<port_direction::cardinal, coordinate<sidb_cell_clk_lyt>>, 4>
        INPUT_PERTURBERS{std::make_pair(port_direction::cardinal::NORTH_WEST, coordinate<sidb_cell_clk_lyt>{10, 2}),
                         std::make_pair(port_direction::cardinal::NORTH_WEST, coordinate<sidb_cell_clk_lyt>{12, 4}),
                         std::make_pair(port_direction::cardinal::NORTH_EAST, coordinate<sidb_cell_clk_lyt>{48, 2}),
                         std::make_pair(port_direction::cardinal::NORTH_EAST, coordinate<sidb_cell_clk_lyt>{46, 4})};

    static inline const std::array<std::pair<port_direction::cardinal, coordinate<sidb_cell_clk_lyt>>, 2>
        OUTPUT_PERTURBERS{std::make_pair(port_direction::cardinal::SOUTH_EAST, coordinate<sidb_cell_clk_lyt>{46, 40}),
                          std::make_pair(port_direction::cardinal::SOUTH_WEST, coordinate<sidb_cell_clk_lyt>{12, 40})};

    static inline const std::array<std::pair<port_direction::cardinal, coordinate<sidb_cell_clk_lyt>>, 2> OUTPUTS{
        std::make_pair(port_direction::cardinal::SOUTH_EAST, coordinate<sidb_cell_clk_lyt>{42, 38}),
        std::make_pair(port_direction::cardinal::SOUTH_WEST, coordinate<sidb_cell_clk_lyt>{16, 38})};
};

/**
 * This struct stores the parameters for the `expected ground state simulation` algorithm.
 */
struct expected_gs_params
{
    /**
     * All parameters for physical SiDB simulations.
     */
    sidb_simulation_parameters physical_parameters{};

    int64_t exact_sim_max_db{27};

    std::vector<double> global_potentials{-0.05};

    std::optional<std::ostream*> kinked_layout_ostream{&std::cout};

    uint64_t number_threads{std::thread::hardware_concurrency()};
};

namespace kinkfinder
{

struct key
{
    kitty::dynamic_truth_table                                      tt{};
    std::pair<port_list<port_direction>, port_list<port_direction>> pl{};
    std::map<port_direction::cardinal, bool>                        in{};
    double                                                          gp{};

    std::size_t operator()(const key& key) const noexcept
    {
        std::size_t h = 0;

        for (const auto& p : key.in)
        {
            hash_combine(h, p);
        }

        hash_combine(h, key.pl);
        hash_combine(h, key.gp);
        kitty::hash_combine(h, kitty::hash<kitty::dynamic_truth_table>()(key.tt));

        return h;
    }

    bool operator==(const key& other) const noexcept
    {
        return tt == other.tt && pl == other.pl && in == other.in && gp == other.gp;
    }
};

struct charge_map_container
{
    std::map<offset::ucoord_t, sidb_charge_state> charge_map{};
    std::map<port_direction::cardinal, bool>      output{};
};

using store = std::unordered_map<key, charge_map_container, key>;

store INIT = store{};

}  // namespace kinkfinder

namespace detail
{

constexpr const char* BIG_TEXT_EXPECTED =
    "  ______  ____   _    _  _   _  _____     _____   _    _ __     __ _____  _____  _____    "
    "        _       _   __     __   _____  _   _ __      __       _       _____  _____     "
    "______ __   __ _____   ______  _____  _______  ______  _____      _____  _____    ____   "
    "_    _  _   _  _____           _____  _______      _______  ______ \n"
    " |  ____|/ __ \\ | |  | || \\ | ||  __ \\   |  __ \\ | |  | |\\ \\   / // ____||_   _|/ "
    "____|    /\\    | |     | |  \\ \\   / /  |_   _|| \\ | |\\ \\    / //\\    | |     |_   "
    "_||  __ \\   |  ____|\\ \\ / /|  __ \\ |  ____|/ ____||__   __||  ____||  __ \\    / "
    "____||  __ \\  / __ \\ | |  | || \\ | ||  __ \\         / ____||__   __| /\\ |__   __||  "
    "____|\n"
    " | |__  | |  | || |  | ||  \\| || |  | |  | |__) || |__| | \\ \\_/ /| (___    | | | |     "
    "   /  \\   | |     | |   \\ \\_/ /     | |  |  \\| | \\ \\  / //  \\   | |       | |  | | "
    " | |  | |__    \\ V / | |__) || |__  | |        | |   | |__   | |  | |  | |  __ | |__) || "
    "|  | || |  | ||  \\| || |  | | ______| (___     | |   /  \\   | |   | |__   \n"
    " |  __| | |  | || |  | || . ` || |  | |  |  ___/ |  __  |  \\   /  \\___ \\   | | | |     "
    "  / /\\ \\  | |     | |    \\   /      | |  | . ` |  \\ \\/ // /\\ \\  | |       | |  | | "
    " | |  |  __|    > <  |  ___/ |  __| | |        | |   |  __|  | |  | |  | | |_ ||  _  / | "
    "|  | || |  | || . ` || |  | ||______|\\___ \\    | |  / /\\ \\  | |   |  __|  \n"
    " | |    | |__| || |__| || |\\  || |__| |  | |     | |  | |   | |   ____) | _| |_| |____  "
    "/ ____ \\ | |____ | |____ | |      _| |_ | |\\  |   \\  // ____ \\ | |____  _| |_ | |__| "
    "|  | |____  / . \\ | |     | |____| |____    | |   | |____ | |__| |  | |__| || | \\ \\ | "
    "|__| || |__| || |\\  || |__| |        ____) |   | | / ____ \\ | |   | |____ \n"
    " |_|     \\____/  \\____/ |_| \\_||_____/   |_|     |_|  |_|   |_|  |_____/ "
    "|_____|\\_____|/_/    \\_\\|______||______||_|     |_____||_| \\_|    \\//_/    "
    "\\_\\|______||_____||_____/   |______|/_/ \\_\\|_|     |______|\\_____|   |_|   "
    "|______||_____/    \\_____||_|  \\_\\ \\____/  \\____/ |_| \\_||_____/        |_____/    "
    "|_|/_/    \\_\\|_|   |______|\n";

constexpr const char* BIG_TEXT_KINKED =
    "  ______  ____   _    _  _   _  _____     _____   _    _ __     __ _____  _____  _____    "
    "        _       _   __     __  __      __       _       _____  _____     _  __ _____  _   "
    "_  _  __ ______  _____      _____  _____    ____   _    _  _   _  _____           _____  "
    "_______      _______  ______ \n"
    " |  ____|/ __ \\ | |  | || \\ | ||  __ \\   |  __ \\ | |  | |\\ \\   / // ____||_   _|/ "
    "____|    /\\    | |     | |  \\ \\   / /  \\ \\    / //\\    | |     |_   _||  __ \\   | "
    "|/ /|_   _|| \\ | || |/ /|  ____||  __ \\    / ____||  __ \\  / __ \\ | |  | || \\ | ||  "
    "__ \\         / ____||__   __| /\\ |__   __||  ____|\n"
    " | |__  | |  | || |  | ||  \\| || |  | |  | |__) || |__| | \\ \\_/ /| (___    | | | |     "
    "   /  \\   | |     | |   \\ \\_/ /    \\ \\  / //  \\   | |       | |  | |  | |  | ' /   "
    "| |  |  \\| || ' / | |__   | |  | |  | |  __ | |__) || |  | || |  | ||  \\| || |  | | "
    "______| (___     | |   /  \\   | |   | |__   \n"
    " |  __| | |  | || |  | || . ` || |  | |  |  ___/ |  __  |  \\   /  \\___ \\   | | | |     "
    "  / /\\ \\  | |     | |    \\   /      \\ \\/ // /\\ \\  | |       | |  | |  | |  |  <    "
    "| |  | . ` ||  <  |  __|  | |  | |  | | |_ ||  _  / | |  | || |  | || . ` || |  | "
    "||______|\\___ \\    | |  / /\\ \\  | |   |  __|  \n"
    " | |    | |__| || |__| || |\\  || |__| |  | |     | |  | |   | |   ____) | _| |_| |____  "
    "/ ____ \\ | |____ | |____ | |        \\  // ____ \\ | |____  _| |_ | |__| |  | . \\  _| "
    "|_ | |\\  || . \\ | |____ | |__| |  | |__| || | \\ \\ | |__| || |__| || |\\  || |__| |    "
    "    ____) |   | | / ____ \\ | |   | |____ \n"
    " |_|     \\____/  \\____/ |_| \\_||_____/   |_|     |_|  |_|   |_|  |_____/ "
    "|_____|\\_____|/_/    \\_\\|______||______||_|         \\//_/    "
    "\\_\\|______||_____||_____/   |_|\\_\\|_____||_| \\_||_|\\_\\|______||_____/    "
    "\\_____||_|  \\_\\ \\____/  \\____/ |_| \\_||_____/        |_____/    |_|/_/    \\_\\|_|  "
    " |______|\n";

template <typename GateLibrary, typename GateLyt>
class kink_finder_impl
{
  public:
    kink_finder_impl(const GateLyt& gate_level_layout, const expected_gs_params& parameter,
                     kinkfinder::store& simulation_store, std::mutex& mutex) :
            gate_lyt{gate_level_layout},
            params{parameter},
            ss{simulation_store},
            mutex{mutex}
    {}

    sidb_simulation_result<sidb_cell_clk_lyt_siqad> run()
    {
        sidb_simulation_result<sqd_lyt> simulation_result{};
        simulation_result.algorithm_name      = "kink_finder";
        simulation_result.physical_parameters = params.physical_parameters;

        mockturtle::stopwatch<>::duration time_counter{};
        {
            const mockturtle::stopwatch stop{time_counter};

            const auto     num_bitstrings = 1ul << gate_lyt.num_pis();
            const uint64_t num_threads    = std::min(std::max(params.number_threads, 1ul), num_bitstrings);

            max_subsim_threads = static_cast<uint64_t>(
                std::max(1l, static_cast<int64_t>(num_threads) - static_cast<int64_t>(num_bitstrings)));

            // define the bit string ranges per thread
            std::vector<std::pair<uint64_t, uint64_t>> ranges{};
            const uint64_t                             chunk_size = std::max(num_bitstrings / num_threads, 1ul);
            uint64_t                                   start      = 0;
            uint64_t                                   end        = chunk_size - 1;

            for (uint64_t i = 0; i < num_threads; ++i)
            {
                ranges.emplace_back(start, end);
                start = end + 1;
                end   = i == num_threads - 2 ? num_bitstrings - 1 : start + chunk_size - 1;
            }

            std::vector<std::thread> threads{};
            threads.reserve(num_threads);

            for (const auto& range : ranges)
            {
                threads.emplace_back(
                    [this, &range, &simulation_result]
                    {
                        for (uint64_t pi_input = range.first; pi_input <= range.second; ++pi_input)
                        {
                            if (this->found_solution())
                            {
                                return;
                            }

                            std::map<tile<GateLyt>, double> empty_global_potential_map{};

                            const std::map<ucoord_t, charge_map_t>& expected_charge_maps =
                                this->get_charge_maps(pi_input, empty_global_potential_map);

                            const auto& expected_charge_lyt = this->apply_charge_maps(pi_input, expected_charge_maps);

                            if (!expected_charge_lyt.is_physically_valid())
                            {
                                const std::lock_guard lock{mutex};

                                if (this->found_solution())
                                {
                                    return;
                                }

                                simulation_result.charge_distributions.emplace_back(expected_charge_lyt);

                                expected_gs_fail = true;

                                return;
                            }

                            const uint64_t num_non_const_nodes{gate_lyt.num_gates() + gate_lyt.num_wires()};

                            for (const auto g_pot : params.global_potentials)
                            {

                                for (uint64_t kinked_node_index = 0; kinked_node_index < num_non_const_nodes;
                                     ++kinked_node_index)
                                {
                                    if (this->found_solution())
                                    {
                                        return;
                                    }

                                    uint64_t count{0};

                                    std::map<tile<GateLyt>, double> global_potential_map{};

                                    gate_lyt.foreach_node(
                                        [this, &kinked_node_index, &count, &global_potential_map, &g_pot](const auto& n)
                                        {
                                            if (!gate_lyt.is_constant(n))
                                            {
                                                if (count == kinked_node_index)
                                                {
                                                    global_potential_map[gate_lyt.get_tile(n)] = g_pot;
                                                }

                                                count++;
                                            }
                                        });

                                    const auto& kinked_charge_maps =
                                        this->get_charge_maps(pi_input, global_potential_map);

                                    if (this->found_solution())
                                    {
                                        return;
                                    }

                                    const auto& charge_lyt = this->apply_charge_maps(pi_input, kinked_charge_maps);

                                    if (charge_lyt.is_physically_valid() &&
                                        charge_lyt.get_system_energy() < expected_charge_lyt.get_system_energy())
                                    {
                                        const std::lock_guard lock{mutex};

                                        if (this->found_solution())
                                        {
                                            return;
                                        }

                                        simulation_result.charge_distributions.emplace_back(charge_lyt);

                                        kinked_gs_energy.emplace(expected_charge_lyt.get_system_energy());

                                        return;
                                    }
                                }
                            }
                        }
                    });
            }

            for (auto& thread : threads)
            {
                thread.join();
            }
        }

        if (found_solution())
        {
            std::cout << fmt::format("FOUND PHYSICALLY {} GROUND STATE",
                                     expected_gs_fail ? "INVALID EXPECTED" : "VALID KINKED")
                      << std::endl;

            if (params.kinked_layout_ostream.has_value())
            {
                *params.kinked_layout_ostream.value()
                    << (expected_gs_fail ? BIG_TEXT_EXPECTED : BIG_TEXT_KINKED) << std::endl;

                print_layout(simulation_result.charge_distributions.front(), *params.kinked_layout_ostream.value());
            }

            std::cout << fmt::format("\n\nPHYSICALLY VALID: {}\n{} {} eV\n",
                                     simulation_result.charge_distributions.front().is_physically_valid(),
                                     expected_gs_fail ? "ENERGY:" : "KINKED GROUND STATE ENERGY:  ",
                                     simulation_result.charge_distributions.front().get_system_energy());

            if (kinked_gs_energy.has_value())
            {
                std::cout << fmt::format("EXPECTED GROUND STATE ENERGY: {} eV\n", kinked_gs_energy.value());
            }
        }
        else
        {
            std::cout << fmt::format("PASSED\n");
        }

        std::cout << fmt::format(
                         "\nTOTAL RUNTIME: {}ms\nINTERNAL RUNTIME: {}ms\n",
                         std::chrono::duration_cast<std::chrono::milliseconds>(time_counter).count(),
                         std::chrono::duration_cast<std::chrono::milliseconds>(time_counter).count() -
                             static_cast<int64_t>(
                                 static_cast<double>(
                                     std::chrono::duration_cast<std::chrono::milliseconds>(subsim_time).count()) /
                                 static_cast<double>(params.number_threads -
                                                     (max_subsim_threads == 1 ? 0 : max_subsim_threads))))
                  << std::endl;

        simulation_result.simulation_runtime = time_counter;

        return simulation_result;
    }

  private:
    using sidb_lyt = sidb_cell_clk_lyt;
    using sqd_lyt  = sidb_cell_clk_lyt_siqad;

    using ucoord_t = coordinate<sidb_lyt>;

    using charge_map_t = std::map<ucoord_t, sidb_charge_state>;
    using io_t         = std::map<port_direction::cardinal, bool>;

    std::mutex&           mutex;
    bool                  expected_gs_fail{false};
    std::optional<double> kinked_gs_energy{};

    const GateLyt& gate_lyt;
    /**
     * Parameter used for the simulation.
     */
    const expected_gs_params& params;

    uint64_t max_subsim_threads{1};
    /**
     * Total sub-simulator runtime.
     */
    std::chrono::duration<double> subsim_time{};

    using input_queue_t = std::vector<std::pair<tile<GateLyt>, kinkfinder::key>>;
    using new_input_t   = std::vector<std::pair<tile<GateLyt>, std::pair<port_direction::cardinal, bool>>>;

    kinkfinder::store& ss;

    [[nodiscard]] inline constexpr bool found_solution() const noexcept
    {
        return expected_gs_fail || kinked_gs_energy.has_value();
    }

    inline constexpr void propagate(const tile<GateLyt>& t, const port_direction pd, const bool& value,
                                    new_input_t& new_i) const
    {
        if (pd.po)
        {
            return;
        }

        switch (pd.dir)
        {
            case port_direction::cardinal::NORTH:
                new_i.emplace_back(gate_lyt.north(t), std::make_pair(port_direction::cardinal::SOUTH, value));
                break;
            case port_direction::cardinal::NORTH_EAST:
                new_i.emplace_back(gate_lyt.north_east(t), std::make_pair(port_direction::cardinal::SOUTH_WEST, value));
                break;
            case port_direction::cardinal::EAST:
                new_i.emplace_back(gate_lyt.east(t), std::make_pair(port_direction::cardinal::WEST, value));
                break;
            case port_direction::cardinal::SOUTH_EAST:
                new_i.emplace_back(gate_lyt.south_east(t), std::make_pair(port_direction::cardinal::NORTH_WEST, value));
                break;
            case port_direction::cardinal::SOUTH:
                new_i.emplace_back(gate_lyt.south(t), std::make_pair(port_direction::cardinal::NORTH, value));
                break;
            case port_direction::cardinal::SOUTH_WEST:
                new_i.emplace_back(gate_lyt.south_west(t), std::make_pair(port_direction::cardinal::NORTH_EAST, value));
                break;
            case port_direction::cardinal::WEST:
                new_i.emplace_back(gate_lyt.west(t), std::make_pair(port_direction::cardinal::EAST, value));
                break;
            case port_direction::cardinal::NORTH_WEST:
                new_i.emplace_back(gate_lyt.north_west(t), std::make_pair(port_direction::cardinal::SOUTH_EAST, value));
                break;
        }
    }

    inline bool check_key_exists(const kinkfinder::key& key, const tile<GateLyt>& t, const ucoord_t& offset,
                                 std::map<ucoord_t, charge_map_t>& cms, new_input_t& new_in) const
    {
        if (ss.count(key))
        {
            for (const auto& out_pd : key.pl.first.out)
            {
                propagate(t, out_pd, ss.at(key).output.at(static_cast<port_direction::cardinal>(out_pd.dir)), new_in);
            }

            cms[offset] = ss.at(key).charge_map;

            return true;
        }

        return false;
    }

    sidb_simulation_result<sqd_lyt> simulate_gate(const tile<GateLyt>& t, const kinkfinder::key& key) const
    {
        sidb_lyt lyt{{GateLibrary::library::gate_x_size(), GateLibrary::library::gate_y_size()}};

        // physical design
        const auto& g = GateLibrary::library::set_up_gate(gate_lyt, t);

        for (auto y = 0ul; y < g.size(); ++y)
        {
            for (auto x = 0ul; x < g[y].size(); ++x)
            {
                const auto type{g[y][x]};

                if (!technology<sidb_lyt>::is_empty_cell(type))
                {
                    lyt.assign_cell_type(cell<sidb_lyt>{x, y}, type);
                }
            }
        }

        // input perturbators
        for (const auto& [dir, b] : key.in)
        {
            for (auto ix = static_cast<uint64_t>(b); ix < GateLibrary::INPUT_PERTURBERS.size(); ix += 2)
            {
                if (GateLibrary::INPUT_PERTURBERS[ix].first == dir)
                {
                    lyt.assign_cell_type(GateLibrary::INPUT_PERTURBERS[ix].second, sidb_lyt::cell_type::NORMAL);
                    break;
                }
            }
        }

        // output perturbators
        for (const auto& out_pd : key.pl.first.out)
        {
            for (const auto& [dir, coord] : GateLibrary::OUTPUT_PERTURBERS)
            {
                if (dir == out_pd.dir)
                {
                    lyt.assign_cell_type(coord, sidb_lyt::cell_type::NORMAL);
                    break;
                }
            }
        }

        sqd_lyt conv_lyt = convert_to_siqad_coordinates(lyt);

        // simulate
        return params.exact_sim_max_db < 0 || conv_lyt.num_cells() <= params.exact_sim_max_db ?
                   quickexact<sqd_lyt>(conv_lyt, quickexact_params<sqd_lyt>{params.physical_parameters,
                                                                            automatic_base_number_detection::ON,
                                                                            {},
                                                                            key.gp}) :
                   composim<sqd_lyt>(conv_lyt, composim_params{quicksim_params{params.physical_parameters, 7, 0.8,
                                                                               key.gp, max_subsim_threads},
                                                               10, 1});
    }

    static std::map<double, std::unordered_map<uint64_t, const charge_distribution_surface<sqd_lyt>*>>
    get_charge_occurrences(const std::vector<charge_distribution_surface<sqd_lyt>>& cds_vec)
    {
        std::map<double, std::unordered_map<uint64_t, const charge_distribution_surface<sqd_lyt>*>>
            charge_occurrences{};

        for (const auto& lyt : cds_vec)
        {
            charge_occurrences[lyt.get_system_energy()][lyt.get_charge_index().first] = &lyt;
        }

        return charge_occurrences;
    }

    void handle_tile(const tile<GateLyt>& t, const kinkfinder::key& key, std::map<ucoord_t, charge_map_t>& charge_maps,
                     new_input_t& new_inputs)
    {
        const auto& n = gate_lyt.get_node(t);

        if (gate_lyt.is_constant(n) || found_solution())
        {
            return;
        }

        if (static_cast<ucoord_t>(t).z == 1)
        {
            return;
        }

        const auto offset =
            relative_to_absolute_cell_position<GateLibrary::library::gate_x_size(), GateLibrary::library::gate_y_size(),
                                               GateLyt, sidb_lyt>(gate_lyt, t, cell<sidb_lyt>{0, 0});

        if (check_key_exists(key, t, offset, charge_maps, new_inputs))
        {
            return;
        }

        const auto& single_gate_simulation = simulate_gate(t, key);

        {
            const std::lock_guard lock{mutex};

            subsim_time += single_gate_simulation.simulation_runtime;
        }

        if (single_gate_simulation.charge_distributions.empty())
        {
            const std::lock_guard lock{mutex};

            std::cout << "ERROR: no ground state found" << std::endl;

            return;
        }

        if (found_solution() || check_key_exists(key, t, offset, charge_maps, new_inputs))
        {
            return;
        }

        const auto& charge_occurrences = get_charge_occurrences(single_gate_simulation.charge_distributions);

        if (charge_occurrences.cbegin()->second.size() != 1)
        {
            const std::lock_guard lock{mutex};

            std::cout << "WARNING: degenerate ground states:\n" << std::endl;

            for (const auto& [_, lyt_p] : charge_occurrences.cbegin()->second)
            {
                print_layout(*lyt_p);
            }

            std::cout << std::endl;
        }

        // get ground state
        const auto* gs = charge_occurrences.cbegin()->second.cbegin()->second;

        if (found_solution() || check_key_exists(key, t, offset, charge_maps, new_inputs))
        {
            return;
        }

        //        if (params.exact_sim_max_db >= 0 && gs->num_cells() > params.exact_sim_max_db)
        //        {
        //            const std::lock_guard lock{mutex};
        //
        //            std::cout << "\n\nINPUT:\n";
        //            for (const auto& [dir, b] : key.in)
        //            {
        //                std::cout << fmt::format("{} -> {}\n", port_direction{dir}, b);
        //            }
        //
        //            std::cout
        //                << fmt::format(
        //                       "\nNUM DBS: {}\nPORTS: {}\nGROUND PORTS: {}\nGLOBAL POTENTIAL: {}\nSUBSIM RUNTIME:
        //                       {}\nCHARGE " "LAYOUTS: {}\nENERGY:{}\n", gs->num_cells(), key.pl.first, key.pl.second,
        //                       key.gp, single_gate_simulation.simulation_runtime.count(),
        //                       single_gate_simulation.charge_distributions.size(), gs->get_system_energy())
        //                << std::endl;
        //            print_layout(*gs);
        //        }

        kinkfinder::charge_map_container cmc{};

        // propagate outputs
        for (const auto& out_pd : key.pl.first.out)
        {
            for (const auto& [dir, coord] : GateLibrary::OUTPUTS)
            {
                if (dir == out_pd.dir)
                {
                    cmc.output[dir] = gs->get_charge_state(siqad::to_siqad_coord(coord)) == sidb_charge_state::NEGATIVE;

                    propagate(t, out_pd, cmc.output[dir], new_inputs);

                    break;
                }
            }
        }

        if (found_solution() || check_key_exists(key, t, offset, charge_maps, new_inputs))
        {
            return;
        }

        // store to charge map
        gs->foreach_cell(
            [&cm = cmc.charge_map, &cds = *gs](const auto& c)
            {
                if (!cds.is_empty_cell(c))
                {
                    cm[siqad::to_fiction_coord<ucoord_t>(c)] = cds.get_charge_state(c);
                }
            });

        {
            const std::lock_guard lock{mutex};

            if (check_key_exists(key, t, offset, charge_maps, new_inputs))
            {
                return;
            }

            ss[key] = cmc;
        }

        charge_maps[offset] = cmc.charge_map;
    }

    void merge_into_input_queue(input_queue_t& input_queue, const new_input_t& new_inputs,
                                std::map<tile<GateLyt>, double>& global_potential_map) const
    {
        for (const auto& p : new_inputs)
        {
            auto it = std::find_if(input_queue.begin(), input_queue.end(),
                                   [&t = p.first](const auto& p) { return p.first == t; });
            if (it == input_queue.end())
            {
                input_queue.emplace_back(p.first, kinkfinder::key{gate_lyt.node_function(gate_lyt.get_node(p.first)),
                                                                  determine_port_routing(gate_lyt, p.first),
                                                                  io_t{p.second}, global_potential_map[p.first]});
            }
            else
            {
                it->second.in[p.second.first] = p.second.second;
            }
        }
    }

    std::map<ucoord_t, charge_map_t> get_charge_maps(uint64_t                         pi_input,
                                                     std::map<tile<GateLyt>, double>& global_potential_map)
    {
        input_queue_t                    input_queue{};
        std::map<ucoord_t, charge_map_t> charge_maps{};

        gate_lyt.foreach_pi(
            [this, &input_queue, &charge_maps, &global_potential_map, &pi_input](const typename GateLyt::node& n,
                                                                                 const uint64_t                i)
            {
                const auto& t  = gate_lyt.get_tile(n);
                const auto& pl = determine_port_routing(gate_lyt, t);

                new_input_t new_inputs{};

                this->handle_tile(t,
                                  kinkfinder::key{gate_lyt.node_function(n), pl,
                                                  io_t{{pl.first.inp.empty() ? GateLibrary::INPUT_PERTURBERS[0].first :
                                                                               static_cast<port_direction::cardinal>(
                                                                                   pl.first.inp.cbegin()->dir),
                                                        static_cast<bool>((1ul << i) & pi_input)}},
                                                  global_potential_map[t]},
                                  charge_maps, new_inputs);

                merge_into_input_queue(input_queue, new_inputs, global_potential_map);
            });

        while (!input_queue.empty())
        {
            std::optional<std::pair<uint64_t, new_input_t>> input_queue_mod{};

            for (auto i = 0; i < input_queue.size(); ++i)
            {
                if (input_queue[i].second.pl.first.inp.size() == input_queue[i].second.in.size())
                {
                    new_input_t new_inputs{};

                    handle_tile(input_queue[i].first, input_queue[i].second, charge_maps, new_inputs);

                    if (found_solution())
                    {
                        return charge_maps;
                    }

                    input_queue_mod = std::make_pair(i, new_inputs);

                    break;
                }
            }

            if (!input_queue_mod.has_value())
            {
                const std::lock_guard lock{mutex};

                std::cout << "ERROR: a set of inputs could not be saturated" << std::endl;

                return charge_maps;
            }

            input_queue.erase(std::next(input_queue.cbegin(), input_queue_mod.value().first));

            merge_into_input_queue(input_queue, input_queue_mod.value().second, global_potential_map);
        }

        return charge_maps;
    }

    [[nodiscard]] charge_distribution_surface<sqd_lyt>
    apply_charge_maps(uint64_t pi_input, const std::map<ucoord_t, charge_map_t>& charge_maps) const
    {
        charge_distribution_surface<sqd_lyt> charge_lyt{
            convert_to_siqad_coordinates(get_cell_level_layout(gate_lyt, pi_input)), params.physical_parameters,
            sidb_charge_state::NEUTRAL};

        // apply charge maps
        for (const auto& [offset, cm] : charge_maps)
        {
            for (const auto& [c, cs] : cm)
            {
                charge_lyt.assign_charge_state(siqad::to_siqad_coord(offset) + siqad::to_siqad_coord(c), cs, false);
            }
        }

        charge_lyt.update_after_charge_change();

        return charge_distribution_surface<sqd_lyt>{charge_lyt};
    }

    template <typename Lyt>
    [[nodiscard]] static std::pair<port_list<port_direction>, port_list<port_direction>>
    determine_port_routing(const Lyt& lyt, const tile<Lyt>& t) noexcept
    {
        port_list<port_direction> p_combined{};
        port_list<port_direction> p_ground{};

        const auto& n = lyt.get_node(t);

        const bool is_pi = lyt.is_pi(n);
        const bool is_po = lyt.is_po(n);

        for (const auto& ti : lyt.z() == 0 || lyt.is_constant(lyt.get_node(lyt.above(t))) ?
                                  std::vector{t} :
                                  std::vector{t, lyt.above(t)})
        {
            for (auto* p : ti == t ? std::vector{&p_combined, &p_ground} : std::vector{&p_combined})
            {
                // determine incoming connector ports
                if (lyt.has_north_eastern_incoming_signal(ti))
                {
                    p->inp.emplace(port_direction::cardinal::NORTH_EAST, is_pi, is_po);
                }
                if (lyt.has_north_western_incoming_signal(ti))
                {
                    p->inp.emplace(port_direction::cardinal::NORTH_WEST, is_pi, is_po);
                }

                // determine outgoing connector ports
                if (lyt.has_south_eastern_outgoing_signal(ti))
                {
                    p->out.emplace(port_direction::cardinal::SOUTH_EAST, is_pi, is_po);
                }
                if (lyt.has_south_western_outgoing_signal(ti))
                {
                    p->out.emplace(port_direction::cardinal::SOUTH_WEST, is_pi, is_po);
                }

                // default PI and PO port routing
                if (is_pi || is_po)
                {
                    if (lyt.has_no_incoming_signal(t))
                    {
                        p->inp.emplace(port_direction::cardinal::NORTH_WEST, is_pi, is_po);
                    }
                    if (lyt.has_no_outgoing_signal(t))
                    {
                        p->out.emplace(port_direction::cardinal::SOUTH_EAST, is_pi, is_po);
                    }
                }
            }
        }

        return std::make_pair(p_combined, p_ground);
    }

    [[nodiscard]] static sidb_lyt get_cell_level_layout(const GateLyt& gate_lyt, uint64_t pi_input)
    {
        sidb_lyt cell_lyt = apply_gate_library<sidb_lyt, typename GateLibrary::library>(gate_lyt);

        std::vector<std::pair<tile<GateLyt>, coordinate<sidb_lyt>>> perturbers{};

        gate_lyt.foreach_pi(
            [&gate_lyt, &pi_input, &perturbers](const auto& n, const auto i)
            {
                const tile<GateLyt>&            t       = gate_lyt.get_tile(n);
                const std::set<port_direction>& inp_pds = determine_port_routing(gate_lyt, t).first.inp;

                const auto value = static_cast<uint64_t>(static_cast<bool>((1ul << i) & pi_input));

                if (inp_pds.empty())
                {
                    perturbers.emplace_back(t, GateLibrary::INPUT_PERTURBERS[value].second);
                }
                else
                {
                    for (auto ix = value; ix < GateLibrary::INPUT_PERTURBERS.size(); ix += 2)
                    {
                        if (GateLibrary::INPUT_PERTURBERS[ix].first == inp_pds.cbegin()->dir)
                        {
                            perturbers.emplace_back(t, GateLibrary::INPUT_PERTURBERS[ix].second);
                            break;
                        }
                    }
                }
            });

        gate_lyt.foreach_po(
            [&gate_lyt, &perturbers](const auto& n)
            {
                const tile<GateLyt>& t =
                    static_cast<tile<GateLyt>>(n);  // TODO: notify marcel; gate_lyt.get_tile(n) does not give the same
                const std::set<port_direction>& out_pds = determine_port_routing(gate_lyt, t).first.out;

                if (out_pds.empty())
                {
                    perturbers.emplace_back(t, GateLibrary::OUTPUT_PERTURBERS.front().second);
                }
                else
                {
                    for (const auto& [dir, coord] : GateLibrary::OUTPUT_PERTURBERS)
                    {
                        if (dir == out_pds.cbegin()->dir)
                        {
                            perturbers.emplace_back(t, coord);
                            break;
                        }
                    }
                }
            });

        for (auto& [tile, coord] : perturbers)
        {
            cell_lyt.assign_cell_type(
                relative_to_absolute_cell_position<GateLibrary::library::gate_x_size(),
                                                   GateLibrary::library::gate_y_size(), GateLyt, sidb_lyt>(gate_lyt,
                                                                                                           tile, coord),
                sidb_lyt::cell_type::NORMAL);
        }

        return cell_lyt;
    }
};

}  // namespace detail

template <typename GateLibrary, typename GateLyt>
sidb_simulation_result<sidb_cell_clk_lyt_siqad>
kink_finder(const GateLyt& gate_level_layout, const expected_gs_params& params = {},
            typename kinkfinder::store& ss    = kinkfinder::INIT,
            std::unique_ptr<std::mutex> mutex = std::make_unique<std::mutex>())
{
    return detail::kink_finder_impl<GateLibrary, GateLyt>{gate_level_layout, params, ss, *mutex}.run();
}

}  // namespace fiction

#endif  // FICTION_KINK_FINDER_HPP

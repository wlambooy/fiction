//
// Created by Willem Lambooy on 15/05/2023.
//

#ifndef FICTION_EXPECTED_GROUND_STATE_SIMULATION_HPP
#define FICTION_EXPECTED_GROUND_STATE_SIMULATION_HPP

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

    std::optional<std::ostream*> kinked_layout_ostream{&std::cout};

    uint64_t number_threads{std::thread::hardware_concurrency()};
};

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
class expected_ground_state_simulation_impl
{
  public:
    expected_ground_state_simulation_impl(const GateLyt& gate_level_layout, const expected_gs_params& parameter) :
            gate_lyt{gate_level_layout},
            params{parameter}
    {}

    sidb_simulation_result<sidb_cell_clk_lyt_siqad> run()
    {
        sidb_simulation_result<sqd_lyt> simulation_result{};
        simulation_result.algorithm_name      = "expected_ground_state_simulation";
        simulation_result.physical_parameters = params.physical_parameters;

        mockturtle::stopwatch<>::duration time_counter{};
        {
            const mockturtle::stopwatch stop{time_counter};

            const auto     num_bitstrings = 1ul << gate_lyt.num_pis();
            const uint64_t num_threads    = std::min(std::max(params.number_threads, 1ul), num_bitstrings);

            max_subsim_threads = std::max(1ul, num_threads - num_bitstrings);

            // define the bit string ranges per thread
            std::vector<std::pair<uint64_t, uint64_t>> ranges;
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

                            std::map<typename GateLyt::node, double> empty_global_potential_map{};

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

                            for (uint64_t kinked_node_index = 0; kinked_node_index < num_non_const_nodes;
                                 ++kinked_node_index)
                            {
                                if (this->found_solution())
                                {
                                    return;
                                }

                                uint64_t count{0};

                                std::map<typename GateLyt::node, double> global_potential_map{};

                                gate_lyt.foreach_node(
                                    [this, &kinked_node_index, &count, &global_potential_map](const auto& n)
                                    {
                                        if (!gate_lyt.is_constant(n))
                                        {
                                            if (count == kinked_node_index)
                                            {
                                                global_potential_map[n] = -0.05;
                                            }

                                            count++;
                                        }
                                    });

                                const auto& kinked_charge_maps = this->get_charge_maps(pi_input, global_potential_map);

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
                             static_cast<long>(
                                 static_cast<double>(
                                     std::chrono::duration_cast<std::chrono::milliseconds>(subsim_time).count()) /
                                 static_cast<double>(params.number_threads - (max_subsim_threads == 1 ? 0 : max_subsim_threads))))
                  << std::endl;

        simulation_result.simulation_runtime = time_counter;

        return simulation_result;
    }

  private:
    using sidb_lyt = sidb_cell_clk_lyt;
    using sqd_lyt  = sidb_cell_clk_lyt_siqad;

    using ucoord_t     = coordinate<sidb_lyt>;
    using charge_map_t = std::map<ucoord_t, sidb_charge_state>;
    using io_t         = std::map<port_direction::cardinal, bool>;

    std::mutex            mutex{};
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

    struct simulation_store_key
    {
        kitty::dynamic_truth_table tt{};
        port_list<port_direction>  pl{};
        io_t                       in{};
        double                     gp{};

        std::size_t operator()(const simulation_store_key& key) const noexcept
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

        bool operator==(const simulation_store_key& other) const noexcept
        {
            return tt == other.tt && pl == other.pl && in == other.in;
        }
    };

    struct charge_map_container
    {
        charge_map_t charge_map{};
        io_t         output{};
    };

    using simulation_store = std::unordered_map<simulation_store_key, charge_map_container, simulation_store_key>;

    simulation_store ss{};

    inline constexpr bool found_solution()
    {
        return expected_gs_fail || kinked_gs_energy.has_value();
    }

    static inline constexpr void propagate(const GateLyt& gate_lyt, const tile<GateLyt>& t, const port_direction& pd,
                                           const bool& value, std::map<tile<GateLyt>, io_t>& input_map)
    {
        if (pd.po)
        {
            return;
        }

        switch (pd.dir)
        {
            case port_direction::cardinal::NORTH:
                input_map[gate_lyt.north(t)][port_direction::cardinal::SOUTH] = value;
                break;
            case port_direction::cardinal::NORTH_EAST:
                input_map[gate_lyt.north_east(t)][port_direction::cardinal::SOUTH_WEST] = value;
                break;
            case port_direction::cardinal::EAST:
                input_map[gate_lyt.east(t)][port_direction::cardinal::WEST] = value;
                break;
            case port_direction::cardinal::SOUTH_EAST:
                input_map[gate_lyt.south_east(t)][port_direction::cardinal::NORTH_WEST] = value;
                break;
            case port_direction::cardinal::SOUTH:
                input_map[gate_lyt.south(t)][port_direction::cardinal::NORTH] = value;
                break;
            case port_direction::cardinal::SOUTH_WEST:
                input_map[gate_lyt.south_west(t)][port_direction::cardinal::NORTH_EAST] = value;
                break;
            case port_direction::cardinal::WEST:
                input_map[gate_lyt.west(t)][port_direction::cardinal::EAST] = value;
                break;
            case port_direction::cardinal::NORTH_WEST:
                input_map[gate_lyt.north_west(t)][port_direction::cardinal::SOUTH_EAST] = value;
                break;
        }
    }

    bool check_key_exists(const simulation_store_key& key, const tile<GateLyt>& t, const port_list<port_direction>& pl,
                          const ucoord_t offset, std::map<ucoord_t, charge_map_t>& cms,
                          std::map<tile<GateLyt>, io_t>& im)
    {
        if (ss.count(key))
        {
            for (const auto& out_pd : pl.out)
            {
                propagate(gate_lyt, t, out_pd, ss[key].output[static_cast<port_direction::cardinal>(out_pd.dir)], im);
            }

            cms[offset] = ss[key].charge_map;

            return true;
        }

        return false;
    }

    void handle_node(const typename GateLyt::node& n, std::map<tile<GateLyt>, io_t>& input_map,
                     std::map<ucoord_t, charge_map_t>& charge_maps, const double& global_potential,
                     const std::optional<bool>& pi_input = std::nullopt)
    {
        if (gate_lyt.is_constant(n))
        {
            return;
        }

        const auto& t = gate_lyt.get_tile(n);

        if (static_cast<ucoord_t>(t).z == 1)
        {
            return;
        }

        const auto offset =
            relative_to_absolute_cell_position<GateLibrary::library::gate_x_size(), GateLibrary::library::gate_y_size(),
                                               GateLyt, sidb_lyt>(gate_lyt, t, cell<sidb_lyt>{0, 0});

        const auto& pl = determine_port_routing(gate_lyt, t);

        const auto& inp = pi_input.has_value() ?
                              io_t{{pl.inp.empty() ? GateLibrary::INPUT_PERTURBERS[0].first :
                                                     static_cast<port_direction::cardinal>(pl.inp.cbegin()->dir),
                                    pi_input.value()}} :
                              input_map[t];

        const auto key = simulation_store_key{gate_lyt.node_function(n), pl, inp, global_potential};

        if (check_key_exists(key, t, pl, offset, charge_maps, input_map))
        {
            return;
        }

        sidb_lyt single_gate_cell_lyt{{GateLibrary::library::gate_x_size(), GateLibrary::library::gate_y_size()}};

        // physical design
        const auto& g = GateLibrary::library::set_up_gate(gate_lyt, t);

        for (auto y = 0ul; y < g.size(); ++y)
        {
            for (auto x = 0ul; x < g[y].size(); ++x)
            {
                const auto type{g[y][x]};

                if (!technology<sidb_lyt>::is_empty_cell(type))
                {
                    single_gate_cell_lyt.assign_cell_type(cell<sidb_lyt>{x, y}, type);
                }
            }
        }

        // input perturbators
        for (const auto& [dir, b] : inp)
        {
            for (auto ix = static_cast<uint64_t>(b); ix < GateLibrary::INPUT_PERTURBERS.size(); ix += 2)
            {
                if (GateLibrary::INPUT_PERTURBERS[ix].first == dir)
                {
                    single_gate_cell_lyt.assign_cell_type(GateLibrary::INPUT_PERTURBERS[ix].second,
                                                          sidb_lyt::cell_type::NORMAL);
                    break;
                }
            }
        }

        // output perturbators
        for (const auto& out_pd : pl.out)
        {
            for (const auto& [dir, coord] : GateLibrary::OUTPUT_PERTURBERS)
            {
                if (dir == out_pd.dir)
                {
                    single_gate_cell_lyt.assign_cell_type(coord, sidb_lyt::cell_type::NORMAL);
                    break;
                }
            }
        }

        // simulate
        sqd_lyt conv_lyt = convert_to_siqad_coordinates(single_gate_cell_lyt);

        const auto& single_gate_simulation =
            params.exact_sim_max_db < 0 || conv_lyt.num_cells() <= params.exact_sim_max_db ?
                quickexact<sqd_lyt>(conv_lyt, quickexact_params<sqd_lyt>{params.physical_parameters,
                                                                         automatic_base_number_detection::ON,
                                                                         {},
                                                                         global_potential}) :
                composim<sqd_lyt>(conv_lyt, composim_params{quicksim_params{params.physical_parameters, 10, 0.8,
                                                                            global_potential, max_subsim_threads},
                                                            20, 1});

        {
            const std::lock_guard lock{mutex};

            subsim_time += single_gate_simulation.simulation_runtime;
        }

        if (check_key_exists(key, t, pl, offset, charge_maps, input_map))
        {
            return;
        }
// global potential -> check how many should be kept      |     apply charge maps checks for phys validity

        // get ground state
        const auto gs = std::min_element(
            single_gate_simulation.charge_distributions.cbegin(), single_gate_simulation.charge_distributions.cend(),
            [](const auto& cd1, const auto& cd2) { return cd1.get_system_energy() < cd2.get_system_energy(); });

        if (gs == single_gate_simulation.charge_distributions.cend())
        {
            std::cout << "ERROR: no ground state found" << std::endl;
            return;
        }

        if (check_key_exists(key, t, pl, offset, charge_maps, input_map))
        {
            return;
        }

//        if (params.exact_sim_max_db >= 0 && conv_lyt.num_cells() > params.exact_sim_max_db)
//        {
//            const std::lock_guard lock{mutex};
//            std::cout << fmt::format(
//                             "\n\nPORTS: {}\nGLOBAL POTENTIAL: {}\nCS RUNTIME: {}\nCHARGE LAYOUTS: {}\nENERGY:{}\n", pl,
//                             global_potential, single_gate_simulation.simulation_runtime.count(),
//                             single_gate_simulation.charge_distributions.size(), gs->get_system_energy())
//                      << std::endl;
//            print_layout(*gs);
//        }

        charge_map_container cmc{};

        // propagate outputs
        for (const auto& out_pd : pl.out)
        {
            for (const auto& [dir, coord] : GateLibrary::OUTPUTS)
            {
                if (dir == out_pd.dir)
                {
                    cmc.output[dir] = gs->get_charge_state(siqad::to_siqad_coord(coord)) == sidb_charge_state::NEGATIVE;

                    propagate(gate_lyt, t, out_pd, cmc.output[dir], input_map);

                    break;
                }
            }
        }

        if (check_key_exists(key, t, pl, offset, charge_maps, input_map))
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

            if (ss.count(key))
            {
                for (const auto& out_pd : pl.out)
                {
                    propagate(gate_lyt, t, out_pd, ss[key].output[static_cast<port_direction::cardinal>(out_pd.dir)],
                              input_map);
                }

                charge_maps[offset] = ss[key].charge_map;

                return;
            }

            ss[key] = cmc;
        }

        charge_maps[offset] = cmc.charge_map;
    }

    std::map<ucoord_t, charge_map_t> get_charge_maps(uint64_t                                  pi_input,
                                                     std::map<typename GateLyt::node, double>& global_potential_map)
    {
        std::map<tile<GateLyt>, io_t>    input_map{};
        std::map<ucoord_t, charge_map_t> charge_maps{};

        gate_lyt.foreach_pi(
            [this, &input_map, &charge_maps, &global_potential_map, &pi_input](const typename GateLyt::node& n,
                                                                               const uint64_t                i)
            {
                this->handle_node(n, input_map, charge_maps, global_potential_map[n],
                                  std::optional{static_cast<bool>((1ul << i) & pi_input)});
            });

        gate_lyt.foreach_node(
            [this, &input_map, &charge_maps, &global_potential_map](const typename GateLyt::node& n)
            {
                if (!gate_lyt.is_pi(n))
                {
                    this->handle_node(n, input_map, charge_maps, global_potential_map[n]);
                }
            });

        return charge_maps;
    }

    charge_distribution_surface<sqd_lyt> apply_charge_maps(uint64_t                                pi_input,
                                                           const std::map<ucoord_t, charge_map_t>& charge_maps)
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
    [[nodiscard]] static port_list<port_direction> determine_port_routing(const Lyt& lyt, const tile<Lyt>& t) noexcept
    {
        port_list<port_direction> p{};

        for (const auto& ti : lyt.z() == 0 || lyt.is_constant(lyt.get_node(lyt.above(t))) ?
                                  std::vector{t} :
                                  std::vector{t, lyt.above(t)})
        {
            // determine incoming connector ports
            if (lyt.has_north_eastern_incoming_signal(ti))
            {
                p.inp.emplace(port_direction::cardinal::NORTH_EAST);
            }
            if (lyt.has_north_western_incoming_signal(ti))
            {
                p.inp.emplace(port_direction::cardinal::NORTH_WEST);
            }

            // determine outgoing connector ports
            if (lyt.has_south_eastern_outgoing_signal(ti))
            {
                p.out.emplace(port_direction::cardinal::SOUTH_EAST);
            }
            if (lyt.has_south_western_outgoing_signal(ti))
            {
                p.out.emplace(port_direction::cardinal::SOUTH_WEST);
            }

            // gates without connector ports

            // 1-input functions
            if (const auto n = lyt.get_node(ti); lyt.is_pi(n) || lyt.is_po(n) || lyt.is_buf(n) || lyt.is_inv(n))
            {
                if (lyt.has_no_incoming_signal(ti))
                {
                    p.inp.emplace(port_direction::cardinal::NORTH_WEST);
                }
                if (lyt.has_no_outgoing_signal(ti))
                {
                    p.out.emplace(port_direction::cardinal::SOUTH_EAST);
                }
            }
            else  // 2-input functions
            {
                if (lyt.has_no_incoming_signal(ti))
                {
                    p.inp.emplace(port_direction::cardinal::NORTH_WEST);
                    p.inp.emplace(port_direction::cardinal::NORTH_EAST);
                }
                if (lyt.has_no_outgoing_signal(ti))
                {
                    p.out.emplace(port_direction::cardinal::SOUTH_EAST);
                }
            }
        }

        return p;
    }

    static sidb_lyt get_cell_level_layout(const GateLyt& gate_lyt, uint64_t pi_input)
    {
        sidb_lyt cell_lyt = apply_gate_library<sidb_lyt, typename GateLibrary::library>(gate_lyt);

        std::vector<std::pair<tile<GateLyt>, coordinate<sidb_lyt>>> perturbers{};

        gate_lyt.foreach_pi(
            [&gate_lyt, &pi_input, &perturbers](const auto& n, const auto i)
            {
                const tile<GateLyt>&            t       = gate_lyt.get_tile(n);
                const std::set<port_direction>& inp_pds = determine_port_routing(gate_lyt, t).inp;

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
                const std::set<port_direction>& out_pds = determine_port_routing(gate_lyt, t).out;

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
sidb_simulation_result<sidb_cell_clk_lyt_siqad> expected_ground_state_simulation(const GateLyt& gate_level_layout,
                                                                                 const expected_gs_params& params = {})
{
    detail::expected_ground_state_simulation_impl<GateLibrary, GateLyt> p{gate_level_layout, params};
    return p.run();
}

}  // namespace fiction

#endif  // FICTION_EXPECTED_GROUND_STATE_SIMULATION_HPP

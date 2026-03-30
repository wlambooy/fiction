//
// Created by Willem Lambooy on 27/01/2025.
//

#include <fiction/algorithms/physical_design/advanced_circuit_design.hpp>
#if (FICTION_Z3_SOLVER)

#include "fiction_experiments.hpp"

#include <fiction/algorithms/network_transformation/technology_mapping.hpp>
// #include <fiction/algorithms/physical_design/advanced_circuit_design.hpp>
// #include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/io/read_sidb_surface_defects.hpp>
#include <fiction/io/write_sqd_layout.hpp>
#include <fiction/layouts/bounding_box.hpp>
#include <fiction/technology/area.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/sidb_bdl_skeletons.hpp>
#include <fiction/technology/sidb_defect_surface.hpp>
#include <fiction/technology/sidb_defects.hpp>
// #include <fiction/technology/sidb_on_the_fly_gate_library.hpp>
#include <fiction/traits.hpp>
#include <fiction/types.hpp>

#include <fmt/format.h>
#include <lorina/lorina.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>

#include <CLI/CLI.hpp>

// This script conducts defect-aware placement and routing with defect-aware on-the-fly SiDB gate design. Thereby, SiDB
// circuits can be designed in the presence of atomic defects.

// This algorithm was proposed in \"On-the-fly Defect-Aware Design of Circuits based on Silicon Dangling Bond Logic\" by
// J. Drewniok, M. Walter, S. S. H. Ng, K. Walus, and R. Wille in IEEE NANO 2024
// (https://ieeexplore.ieee.org/abstract/document/10628962).

// #define USE_MINI

namespace fs = std::filesystem;

template <typename lyt_t, typename skeleton>
std::pair<advanced_circuit_design_params<lyt_t>, std::string>
parse_params(const int argc, char** argv,
             // design_sidb_gates_params<lyt_t>& design_gate_params,
             const std::optional<lyt_t>& surface_lattice)
{
    advanced_circuit_design_params<lyt_t> params{};

    CLI::App app{"SiDB circuit design parameters"};

    // ---- Temporary signed variables for parsing ----
    int64_t canvas_sidbs_signed = static_cast<int64_t>(params.maximum_number_of_canvas_sidbs);
    // int64_t max_gate_designs_signed            = static_cast<int64_t>(params.maximum_number_of_solutions);
    // int64_t design_gate_threads_signed         = static_cast<int64_t>(design_gate_params.available_threads);
    int64_t num_trials_signed                  = static_cast<int64_t>(params.num_trials);
    int64_t max_discrimination_attempts_signed = static_cast<int64_t>(params.maximum_discrimination_attempts);
    int64_t max_repeated_discrimination_attempts_signed =
        static_cast<int64_t>(params.maximum_repeated_discrimination_attempts);
    int64_t threads_signed = static_cast<int64_t>(params.available_threads);

    std::vector canvas_x{static_cast<int64_t>(params.canvas.first.x), static_cast<int64_t>(params.canvas.second.x)};
    std::vector canvas_y{static_cast<int64_t>(params.canvas.first.y), static_cast<int64_t>(params.canvas.second.y)};

    app.add_option("--canvas_x", canvas_x, "Canvas X range (xmin,xmax)")->delimiter(',');
    app.add_option("--canvas_y", canvas_y, "Canvas Y range (ymin,ymax)")->delimiter(',');

    // Gate design parameters
    app.add_option("--canvas_sidbs", canvas_sidbs_signed, "Number of canvas SiDBs");
    // app.add_option("--max_gate_designs", max_gate_designs_signed, "Maximum number of initial gate designs");
    // app.add_option("--design_gate_threads", design_gate_threads_signed, "Available gate design threads (0 = all)");

    // Circuit design parameters
    app.add_option("--num_trials", num_trials_signed, "Number of trials");

    std::string sub_circuit_mode = "c";
    app.add_option("--sub_circuit_mode", sub_circuit_mode,
                   "Sub-circuit mode: c (connected), n (neighboring/adjacent), a (all)");

    app.add_option("--quantization_step", params.quantization_step, "Quantization step");

    std::string quantization_mode = "v";
    app.add_option("--quantization_mode", quantization_mode, "Quantization mode: v (visual), p (pruning)");

    app.add_option("--selectivity", params.selectivity, "Selectivity");
    app.add_option("--selectivity_tolerance", params.selectivity_tolerance, "Selectivity tolerance");
    app.add_option("--success_rate_ceiling", params.success_rate_ceiling, "Success rate ceiling");

    app.add_option("--max_discrimination_attempts", max_discrimination_attempts_signed,
                   "Maximum number of discrimination attempts");
    app.add_option("--max_repeated_discrimination_attempts", max_repeated_discrimination_attempts_signed,
                   "Maximum number of repeated discrimination attempts");

    std::string gate_design_mode = "r";
    app.add_option("--gate_design_mode", gate_design_mode, "Gate design mode: r (random), e (exhaustive)");

    app.add_option("--threads", threads_signed, "Available threads (0 = all)");

    std::string func = "FREEANDEX";
    app.add_option("--func", func, "Function name");

    app.add_flag("--print_solutions", params.print_found_circuits, "Print found solutions during exhaustive search");

    // Parse CLI args
    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        std::exit(app.exit(e));
    }

    // ---- Convert to unsigned after validation ----
    auto to_uint64_checked = [](const int64_t value, const std::string& name, const bool allow_zero = false) -> uint64_t
    {
        if (value < 0)
        {
            throw std::invalid_argument(name + " cannot be negative");
        }

        if (!allow_zero && value == 0)
        {
            throw std::invalid_argument(name + " must be > 0");
        }

        return static_cast<uint64_t>(value);
    };

    params.maximum_number_of_canvas_sidbs = to_uint64_checked(canvas_sidbs_signed, "canvas_sidbs");
    // design_gate_params.maximum_number_of_solutions = to_uint64_checked(max_gate_designs_signed, "max_gate_designs");
    // design_gate_params.available_threads = to_uint64_checked(design_gate_threads_signed, "design_gate_threads",
    // true);
    params.num_trials = to_uint64_checked(num_trials_signed, "num_trials");
    params.maximum_discrimination_attempts =
        to_uint64_checked(max_discrimination_attempts_signed, "max_discrimination_attempts");
    params.maximum_repeated_discrimination_attempts =
        to_uint64_checked(max_repeated_discrimination_attempts_signed, "max_repeated_discrimination_attempts");
    params.available_threads = to_uint64_checked(threads_signed, "threads", true);

    // Normalize mode strings
    auto string_to_lower = [](std::string& s)
    { std::transform(s.begin(), s.end(), s.begin(), [](const unsigned char c) { return std::tolower(c); }); };

    string_to_lower(sub_circuit_mode);
    string_to_lower(quantization_mode);
    string_to_lower(gate_design_mode);

    // Map sub_circuit_mode
    if (sub_circuit_mode == "connected" || sub_circuit_mode == "c")
    {
        params.sub_circuit_mode = decltype(params)::CONNECTED_GATES;
    }
    else if (sub_circuit_mode == "neighboring" || sub_circuit_mode == "adjacent" || sub_circuit_mode == "n")
    {
        params.sub_circuit_mode = decltype(params)::ADJACENT_GATES;
    }
    else if (sub_circuit_mode == "all" || sub_circuit_mode == "a")
    {
        params.sub_circuit_mode = decltype(params)::ALL_GATES;
    }
    else
    {
        throw std::invalid_argument("Invalid sub_mode: must be c, n, or a");
    }

    // Quantization mode
    if (quantization_mode == "p" || quantization_mode == "pruning")
    {
        params.quantize_mode = advanced_circuit_design_params<lyt_t>::quantization_mode::SOFTEN_PRUNING;
    }

    // // Gate design mode
    // if (gate_design_mode == "e" || gate_design_mode == "exhaustive")
    // {
    //     design_gate_params.design_mode = design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::EXHAUSTIVE;
    // }

    // Argument validation
    auto check_positive = [](auto value, const std::string& name)
    {
        if (value <= 0)
        {
            throw std::invalid_argument(name + " must be > 0");
        }
    };

    check_positive(params.maximum_number_of_canvas_sidbs, "canvas_sidbs");
    // check_positive(design_gate_params.maximum_number_of_solutions, "max_gate_designs");
    // if (design_gate_params.available_threads < 0)
    // {
    //     throw std::invalid_argument("design_gate_threads must be >= 0");
    // }
    check_positive(params.num_trials, "num_trials");

    if (params.quantization_step < 0.0 || params.quantization_step > 10.0)
    {
        throw std::invalid_argument("quantization_step must be in [0, 10]");
    }
    if (params.selectivity <= 0.0 || params.selectivity >= 1.0)
    {
        throw std::invalid_argument("selectivity must be in (0, 1)");
    }
    if (params.selectivity_tolerance < 0.0 || params.selectivity_tolerance > 1.0)
    {
        throw std::invalid_argument("selectivity_tolerance must be in [0, 1]");
    }
    if (params.success_rate_ceiling <= 0.5 || params.success_rate_ceiling > 1.0)
    {
        throw std::invalid_argument("success_rate_ceiling must be in (0.5, 1]");
    }

    check_positive(params.maximum_discrimination_attempts, "max_discrimination_attempts");
    check_positive(params.maximum_repeated_discrimination_attempts, "max_repeated_discrimination_attempts");
    check_positive(params.available_threads, "threads");

    auto normalize_threads = [](uint64_t& threads)
    {
        if (threads == 0)
        {
            threads = std::thread::hardware_concurrency();
        }
    };

    // Threads handling (0 -> all)
    normalize_threads(params.available_threads);
    // normalize_threads(design_gate_params.available_threads);

    // --- Canvas validation ---
    if (canvas_x.size() != 2 || canvas_y.size() != 2)
    {
        throw std::invalid_argument("Canvas ranges must have exactly 2 values each (xmin,xmax and ymin,ymax)");
    }

    const int64_t x_min = canvas_x[0];
    const int64_t x_max = canvas_x[1];
    const int64_t y_min = canvas_y[0];
    const int64_t y_max = canvas_y[1];

    if (x_min < 0 || x_max < x_min || x_max >= static_cast<int64_t>(skeleton::gate_x_size()))
    {
        throw std::invalid_argument("Invalid canvas X range: 0 <= xmin <= xmax < skeleton::gate_x_size() = " +
                                    std::to_string(skeleton::gate_x_size()));
    }
    if (y_min < 0 || y_max < y_min || y_max >= static_cast<int64_t>(skeleton::gate_y_size()))
    {
        throw std::invalid_argument("Invalid canvas Y range: 0 <= ymin <= ymax < skeleton::gate_y_size() = " +
                                    std::to_string(skeleton::gate_y_size()));
    }

    params.canvas = {{static_cast<uint64_t>(x_min), static_cast<uint64_t>(y_min)},
                     {static_cast<uint64_t>(x_max), static_cast<uint64_t>(y_max)}};

    // Fixed/default parameters
    params.exact_design_parameters.scheme        = "ROW4";
    params.exact_design_parameters.crossings     = true;
    params.exact_design_parameters.border_io     = false;
    params.exact_design_parameters.desynchronize = true;
    // params.exact_design_parameters.upper_bound_x = 4;          // 5 x 5 tiles
    // params.exact_design_parameters.upper_bound_y = 4;          // 5 x 5 tiles
    params.exact_design_parameters.upper_bound_x = 11;         // 5 x 5 tiles
    params.exact_design_parameters.upper_bound_y = 30;         // 5 x 5 tiles
    params.exact_design_parameters.timeout       = 3'600'000;  // 1h in ms

    params.defect_surface = surface_lattice;
    // params.design_gate_params = design_gate_params;

    return {params, func};
}

int main(int argc, char* argv[])  // NOLINT
{
    for (int i = 0; i < argc; ++i)
    {
        std::cout << argv[i] << ' ';
    }
    std::cout << '\n' << std::endl;

    using gate_lyt = hex_even_row_gate_clk_lyt;
    using cell_lyt = sidb_cell_clk_lyt_cube;
    using lyt_t    = sidb_defect_surface<cell_lyt>;

    // using skeleton = sidb_bdl_skeleton_original_bestagon;
    // using skeleton = sidb_bdl_skeleton_1;
    using skeleton = sidb_bdl_skeleton_hexamini;
    //
    // /// DESIGN GATE PARAMS
    //
    // design_sidb_gates_params<lyt_t> design_gate_params{};
    //
    // design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{3, -0.32};
    // // design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{3, -0.26};
    // // design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params.threshold_bdl_interdistance =
    // 3;
    // //
    // design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params.bdl_pairs_params.minimum_distance
    // // = 0.5;
    // //
    // design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params.bdl_pairs_params.maximum_distance
    // // = 1.7;
    // // =; sidb_simulation_parameters{2, -0.32};
    // design_gate_params.operational_params.op_condition_positive_charges =
    //     is_operational_params<cell<lyt_t>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
    // design_gate_params.operational_params.op_condition_kinks =
    //     is_operational_params<cell<lyt_t>>::operational_condition_kinks::REJECT_KINKS;
    // design_gate_params.design_mode = design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::RANDOM;
    //
    // // design_gate_params.post_design_process = {
    // //     std::make_unique<compare_by_minimum_ground_state_isolation<lyt_t>>(),
    // //     std::make_unique<compare_by_average_ground_state_isolation<lyt_t>>()};
    //
    // // needs to be changed if a different skeleton is used.
    // // design_gate_params.canvas = {{10, 8}, {31, 19}}; // new_mini_bestagon_august.sqd
    // // design_gate_params.canvas = {{23, 12}, {37, 25}};  // original_bestagon.sqd
    // design_gate_params.canvas = {{8, 11}, {18, 20}};  // new_mini_bestagon.sqd
    // // design_gate_params.canvas = {{5, 8}, {21, 21}};  // new_mini_bestagon.sqd
    //
    // design_gate_params.number_of_canvas_sidbs        = 4;
    // design_gate_params.operational_params.sim_engine = sidb_simulation_engine::CLUSTERCOMPLETE;
    // design_gate_params.termination_cond =
    // design_sidb_gates_params<lyt_t>::termination_condition::OBTAINED_N_SOLUTIONS;
    // design_gate_params.maximum_number_of_solutions = 200;

    // // save atomic defects which their respective physical parameters as experimentally determined by T. R. Huff, T.
    // // Dienel, M. Rashidi, R. Achal, L. Livadaru, J. Croshaw, and R. A. Wolkow, "Electrostatic landscape of a
    // // Hydrogen-terminated Silicon Surface Probed by a Moveable Quantum Dot."
    // const auto stray_db   = sidb_defect{sidb_defect_type::DB, -1, 4.1, 1.8};
    // const auto si_vacancy = sidb_defect{sidb_defect_type::SI_VACANCY, -1, 10.6, 5.9};

    // const auto lyt = read_sqd_layout<sidb_100_cell_clk_lyt>("handmadeXNOR.sqd", "file");
    // const auto sim_res = clustercomplete(lyt);
    // print_layout(sim_res.groundstates().front());
    // auto gs = sim_res.groundstates().front();
    // for (const auto& cell : gs.get_sidb_order())
    // {
    //     std::cout << cell.x << ' ' << cell.y << " : " << *gs.get_local_internal_potential(cell) << std::endl;
    //     if ((cell.y >= 32 && cell.y <= 50 || cell.x > 48) && gs.get_charge_state(cell) ==
    //     sidb_charge_state::NEGATIVE)
    //     {
    //         gs.assign_charge_state(cell, sidb_charge_state::NEUTRAL);
    //     }
    // }
    // print_layout(gs);
    // gs.update_after_charge_change();
    // for (const auto& cell : gs.get_sidb_order())
    // {
    //     std::cout << cell.x << ' ' << cell.y << " : " << *gs.get_local_internal_potential(cell) << std::endl;
    // }

    static const std::string layouts_folder =
        fmt::format("{}/physical_design_with_on_the_fly_gate_design/layouts", EXPERIMENTS_PATH);

    // // read-in the initial defects. Physical parameters of the defects are not stored yet.
    // auto surface_lattice_initial = read_sidb_surface_defects<cell_lyt>(
    //     "../../experiments/physical_design_with_on_the_fly_gate_design/1_percent_with_charged_surface.txt");

    // create an empty surface.
    sidb_defect_surface<cell_lyt> surface_lattice{};

    // // add physical parameters of the defects to the surface_lattice.
    // surface_lattice_initial.foreach_sidb_defect(
    //     [&surface_lattice, &stray_db, &si_vacancy](const auto& cd)
    //     {
    //         if (cd.second.type == sidb_defect_type::DB)
    //         {
    //             surface_lattice.assign_sidb_defect(cd.first, stray_db);
    //         }
    //         else if (cd.second.type == sidb_defect_type::SI_VACANCY)
    //         {
    //             surface_lattice.assign_sidb_defect(cd.first, si_vacancy);
    //         }
    //         else
    //         {
    //             surface_lattice.assign_sidb_defect(cd.first, cd.second);
    //         }
    //     });

    // determine bounding-box of the surface to set the aspect ratio of the surface lattice.
    const auto bb_defect_surface = bounding_box_2d{surface_lattice};
    surface_lattice.resize(bb_defect_surface.get_max());

    experiments::experiment<std::string, double, uint64_t, bool> sidb_circuits_with_defects{
        "sidb_circuits_with_defects", "benchmark", "runtime", "number of aspect ratios", "success"};

    const fs::path b_dir = fmt::format("{}", EXPERIMENTS_PATH);

    const auto& [params, func] = parse_params<lyt_t, skeleton>(argc, argv, surface_lattice);

    advanced_circuit_design_stats<cart_odd_row_gate_clk_lyt> st{};

    cart_odd_row_gate_clk_lyt g{{7, 5}};
    const auto                i1 = g.create_pi("1", {0, 0});
    const auto                i2 = g.create_pi("2", {7, 0});
    const auto                b1 = g.create_buf(i1, {1, 1});
    const auto                b2 = g.create_buf(i2, {5, 1});
    const auto                a1 = g.create_and(i1, i2, {3, 2});
    // const auto                a2 = g.create_and(i1, i2, {4, 2});
    // const auto                a3 = g.create_and(i1, i2, {3, 3});
    const auto o = g.create_po(a1, "o", {4, 4});
    const auto p = g.create_po(o, "p", {5, 5});

    std::unordered_map<mockturtle::node<cart_odd_row_gate_clk_lyt>, std::pair<uint8_t, uint8_t>> sidb_count_map{};
    sidb_count_map[g.get_node({0, 0})] = {2, 2};
    sidb_count_map[g.get_node({7, 0})] = {2, 2};
    sidb_count_map[g.get_node({1, 1})] = {2, 2};
    sidb_count_map[g.get_node({5, 1})] = {2, 2};
    sidb_count_map[g.get_node({3, 2})] = {0, 1};
    sidb_count_map[g.get_node({4, 2})] = {0, 1};
    sidb_count_map[g.get_node({3, 3})] = {0, 1};
    sidb_count_map[g.get_node({4, 4})] = {2, 2};
    sidb_count_map[g.get_node({5, 5})] = {1, 1};

    g.foreach_node(
        [&](const auto& n, [[maybe_unused]] auto i)
        {
            if (g.is_constant(n))
            {
                return;
            }

            const tile<cart_odd_row_gate_clk_lyt>& t = g.get_tile(n);

            const cell<lyt_t> nw_position = relative_to_absolute_cell_position<3, 4, cart_odd_row_gate_clk_lyt, lyt_t>(
                g, t, {g.is_in_odd_row(t) ? 1 : 0, 0});
            const cell<lyt_t> se_position =
                relative_to_absolute_cell_position<3, 4, cart_odd_row_gate_clk_lyt, lyt_t>(g, t, {2, 3});

            const std::vector<cell<lyt_t>>& all_sidbs_in_canvas =
                all_coordinates_in_spanned_area(nw_position, se_position);

            std::cout << "tile: " << t << std::endl;
            for (const auto& c : all_sidbs_in_canvas)
            {
                std::cout << c.x << " " << c.y << std::endl;
            }
            std::cout << std::endl;
        });

    const std::vector<lyt_t>& lyts =
        advanced_circuit_design<lyt_t, cart_odd_row_gate_clk_lyt>(g, sidb_count_map, params, &st);

    // write a SiQAD simulation file
    // write_sqd_layout(*lyt, layout_path.c_str());

    create_directory(b_dir / "exact_benchmarks_layout" / func);

    uint64_t ix = 0;
    for (const auto& lyt : lyts)
    {
        write_sqd_layout(lyt, (b_dir / "exact_benchmarks_layout" / func / (std::to_string(ix++) + ".sqd")).c_str());
    }

    // write runtime to file
    const auto    runtime_path = b_dir / "exact_benchmarks_runtime" /*/ num_in_dir */ / (func + ".txt");
    std::ofstream os{runtime_path, std::ofstream::out};
    if (!os.is_open())
    {
        throw std::ofstream::failure("could not open file");
    }
    const auto runtime_string = fmt::format("{:.2f}", mockturtle::to_seconds(st.time_total));
    os.write(runtime_string.c_str(), static_cast<uint32_t>(runtime_string.size()));

    return EXIT_SUCCESS;

    // sidb_circuits_with_defects(func, mockturtle::to_seconds(st.time_total), st.exact_stats.num_aspect_ratios,
    //                            !lyts.empty());
    // sidb_circuits_with_defects.save();
    // sidb_circuits_with_defects.table();
    //
    // return EXIT_SUCCESS;
}

#else  // FICTION_Z3_SOLVER

#include <cstdlib>
#include <iostream>

int main()  // NOLINT
{
    std::cerr << "[e] Z3 solver is not available, please install Z3 and recompile the code" << std::endl;

    return EXIT_FAILURE;
}

#endif  // FICTION_Z3_SOLVER

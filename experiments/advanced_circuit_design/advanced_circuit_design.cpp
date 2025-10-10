//
// Created by Willem Lambooy on 27/01/2025.
//

#if (FICTION_Z3_SOLVER)

#include "../alice/lib/cli11/CLI11.hpp"
#include "fiction_experiments.hpp"

#include <fiction/algorithms/network_transformation/technology_mapping.hpp>
#include <fiction/algorithms/physical_design/advanced_circuit_design.hpp>
#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
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
#include <fiction/technology/sidb_on_the_fly_gate_library.hpp>
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

// This script conducts defect-aware placement and routing with defect-aware on-the-fly SiDB gate design. Thereby, SiDB
// circuits can be designed in the presence of atomic defects.

// This algorithm was proposed in \"On-the-fly Defect-Aware Design of Circuits based on Silicon Dangling Bond Logic\" by
// J. Drewniok, M. Walter, S. S. H. Ng, K. Walus, and R. Wille in IEEE NANO 2024
// (https://ieeexplore.ieee.org/abstract/document/10628962).

// #define USE_MINI

namespace fs = std::filesystem;

template <typename lyt_t>
advanced_circuit_design_params<lyt_t> parse_params(const int argc, char** argv,
                                                   design_sidb_gates_params<lyt_t>& design_gate_params,
                                                   const std::optional<lyt_t>&      surface_lattice)
{
    advanced_circuit_design_params<lyt_t> params{};

    CLI::App app{"SiDB circuit design parameters"};

    // ---- Temporary signed variables for parsing ----
    int64_t canvas_sidbs_signed;
    int64_t max_gate_designs_signed;
    int64_t design_gate_threads_signed;
    int64_t num_trials_signed;
    int64_t max_discrimination_attempts_signed;
    int64_t max_repeated_discrimination_attempts_signed;
    int64_t threads_signed;

    // Gate design parameters
    app.add_option("--canvas_sidbs", canvas_sidbs_signed, "Number of canvas SiDBs");
    app.add_option("--max_gate_designs", max_gate_designs_signed, "Maximum number of initial gate designs");
    app.add_option("--design_gate_threads", design_gate_threads_signed, "Available gate design threads (0 = all)");

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

    design_gate_params.number_of_canvas_sidbs      = to_uint64_checked(canvas_sidbs_signed, "canvas_sidbs");
    design_gate_params.maximum_number_of_solutions = to_uint64_checked(max_gate_designs_signed, "max_gate_designs");
    design_gate_params.available_threads = to_uint64_checked(design_gate_threads_signed, "design_gate_threads", true);
    params.num_trials                    = to_uint64_checked(num_trials_signed, "num_trials");
    params.maximum_discrimination_attempts =
        to_uint64_checked(max_discrimination_attempts_signed, "max_discrimination_attempts");
    params.maximum_repeated_discrimination_attempts =
        to_uint64_checked(max_repeated_discrimination_attempts_signed, "max_repeated_discrimination_attempts");
    params.available_threads = to_uint64_checked(threads_signed, "threads", true);

    // Normalize mode strings
    std::transform(sub_circuit_mode.begin(), sub_circuit_mode.end(), sub_circuit_mode.begin(),
                   [](const unsigned char c) { return std::tolower(c); });
    std::transform(quantization_mode.begin(), quantization_mode.end(), quantization_mode.begin(),
                   [](const unsigned char c) { return std::tolower(c); });
    std::transform(gate_design_mode.begin(), gate_design_mode.end(), gate_design_mode.begin(),
                   [](const unsigned char c) { return std::tolower(c); });

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

    // Gate design mode
    if (gate_design_mode == "e" || gate_design_mode == "exhaustive")
    {
        design_gate_params.design_mode = design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::EXHAUSTIVE;
    }

    // Argument validation
    auto check_positive = [](auto value, const std::string& name)
    {
        if (value <= 0)
        {
            throw std::invalid_argument(name + " must be > 0");
        }
    };

    check_positive(design_gate_params.number_of_canvas_sidbs, "canvas_sidbs");
    check_positive(design_gate_params.maximum_number_of_solutions, "max_gate_designs");
    if (design_gate_params.available_threads < 0)
    {
        throw std::invalid_argument("design_gate_threads must be >= 0");
    }
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
    normalize_threads(design_gate_params.available_threads);

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

    params.defect_surface     = surface_lattice;
    params.design_gate_params = design_gate_params;

    return params;
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

    using skeleton = sidb_bdl_skeleton_1;

    /// DESIGN GATE PARAMS

    design_sidb_gates_params<lyt_t> design_gate_params{};

    design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{3, -0.32};
    // design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{3, -0.26};
    // design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params.threshold_bdl_interdistance = 3;
    // design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params.bdl_pairs_params.minimum_distance
    // = 0.5;
    // design_gate_params.operational_params.input_bdl_iterator_params.bdl_wire_params.bdl_pairs_params.maximum_distance
    // = 1.7;
    // =; sidb_simulation_parameters{2, -0.32};
    design_gate_params.operational_params.op_condition_positive_charges =
        is_operational_params<cell<lyt_t>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
    design_gate_params.operational_params.op_condition_kinks =
        is_operational_params<cell<lyt_t>>::operational_condition_kinks::REJECT_KINKS;
    design_gate_params.design_mode = design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::RANDOM;

    // design_gate_params.post_design_process = {
    //     std::make_unique<compare_by_minimum_ground_state_isolation<lyt_t>>(),
    //     std::make_unique<compare_by_average_ground_state_isolation<lyt_t>>()};

    // needs to be changed if a different skeleton is used.
    // design_gate_params.canvas = {{10, 8}, {31, 19}}; // new_mini_bestagon_august.sqd
    // design_gate_params.canvas = {{23, 12}, {37, 25}};  // original_bestagon.sqd
    design_gate_params.canvas = {{8, 11}, {18, 20}};  // new_mini_bestagon.sqd
    // design_gate_params.canvas = {{5, 8}, {21, 21}};  // new_mini_bestagon.sqd

    design_gate_params.number_of_canvas_sidbs        = 4;
    design_gate_params.operational_params.sim_engine = sidb_simulation_engine::CLUSTERCOMPLETE;
    design_gate_params.termination_cond = design_sidb_gates_params<lyt_t>::termination_condition::OBTAINED_N_SOLUTIONS;
    design_gate_params.maximum_number_of_solutions = 200;

    // // save atomic defects which their respective physical parameters as experimentally determined by T. R. Huff, T.
    // // Dienel, M. Rashidi, R. Achal, L. Livadaru, J. Croshaw, and R. A. Wolkow, "Electrostatic landscape of a
    // // Hydrogen-terminated Silicon Surface Probed by a Moveable Quantum Dot."
    // const auto stray_db   = sidb_defect{sidb_defect_type::DB, -1, 4.1, 1.8};
    // const auto si_vacancy = sidb_defect{sidb_defect_type::SI_VACANCY, -1, 10.6, 5.9};

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

    const auto lattice_tiling = gate_lyt{{11, 30}};  // todo

    experiments::experiment<std::string, double, uint64_t, bool> sidb_circuits_with_defects{
        "sidb_circuits_with_defects", "benchmark", "runtime", "number of aspect ratios", "success"};

    const std::string base_dir = fmt::format("{}benchmarks_verilog", EXPERIMENTS_PATH);

    for (int n = 2; n <= 4; ++n)
    {
        const fs::path input_dir = base_dir + "/" + std::to_string(n) + "_in";

        for (const auto& entry : fs::directory_iterator(input_dir))
        {
            const auto& verilog_path = entry.path();  // e.g., <something>/benchmarks_verilog/3_in/foo.v
            const auto  num_in_dir   = verilog_path.parent_path().filename();                   // "3_in"
            const auto  name         = verilog_path.stem();                                     // "foo"
            const auto  b_dir        = verilog_path.parent_path().parent_path().parent_path();  // <something>
            const auto  layout_path  = b_dir / "benchmarks_layout" / num_in_dir / (name.string() + ".sqd");
            if (std::ifstream is{layout_path.c_str()}; is.is_open())
            {
                continue;
            }

            const auto benchmark = entry.path().string();
            fmt::print("[attempts] processing {}\n", benchmark);

            mockturtle::xag_network xag{};
            const auto              result = lorina::read_verilog(benchmark, mockturtle::verilog_reader(xag));
            assert(result == lorina::return_code::success);

            // compute depth
            const mockturtle::depth_view depth_xag{xag};

            const technology_mapping_params tech_map_params = all_2_input_functions();

            // parameters for cut rewriting
            mockturtle::cut_rewriting_params cut_params{};
            cut_params.cut_enumeration_ps.cut_size = 4;

            const mockturtle::xag_npn_resynthesis<
                mockturtle::xag_network,                    // the input network type
                mockturtle::xag_network,                    // the database network type
                mockturtle::xag_npn_db_kind::xag_complete>  // the kind of database to use

                resynthesis_function{};

            // rewrite network cuts using the given re-synthesis function
            const auto cut_xag = mockturtle::cut_rewriting(xag, resynthesis_function, cut_params);

            // perform technology mapping
            const auto mapped_network = technology_mapping(cut_xag, tech_map_params);

            const advanced_circuit_design_params<lyt_t> params =
                parse_params<lyt_t>(argc, argv, design_gate_params, surface_lattice);

            advanced_circuit_design_stats<gate_lyt> st{};

            const std::optional<lyt_t>& lyt =
                advanced_circuit_design<decltype(mapped_network), lyt_t, gate_lyt, skeleton>(
                    mapped_network, lattice_tiling, params, &st);

            if (!lyt.has_value())
            {
                sidb_circuits_with_defects(benchmark, 0, st.exact_stats.num_aspect_ratios, false);
                sidb_circuits_with_defects.save();
                sidb_circuits_with_defects.table();
                continue;
            }

            // check equivalence
            const auto miter = mockturtle::miter<mockturtle::klut_network>(mapped_network, st.gate_layout.value());
            // const auto eq    = mockturtle::equivalence_checking(*miter);
            // assert(eq.has_value());

            // determine bounding box and exclude atomic defects
            const auto bb = bounding_box_2d<cell_lyt>(static_cast<cell_lyt>(*lyt));

            // write a SiQAD simulation file
            write_sqd_layout(*lyt, layout_path.c_str());

            // write runtime to file
            const auto    runtime_path = b_dir / "benchmarks_runtime" / num_in_dir / (name.string() + ".txt");
            std::ofstream os{runtime_path, std::ofstream::out};
            if (!os.is_open())
            {
                throw std::ofstream::failure("could not open file");
            }
            const auto runtime_string = fmt::format("{:.2f}", mockturtle::to_seconds(st.time_total));
            os.write(runtime_string.c_str(), static_cast<uint32_t>(runtime_string.size()));

            // compute area
            area_stats                   area_stats{};
            area_params<sidb_technology> area_ps{};
            area(bb, area_ps, &area_stats);

            sidb_circuits_with_defects(benchmark, mockturtle::to_seconds(st.time_total),
                                       st.exact_stats.num_aspect_ratios, lyt.has_value());
            sidb_circuits_with_defects.save();
            sidb_circuits_with_defects.table();
        }
    }

    return EXIT_SUCCESS;
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

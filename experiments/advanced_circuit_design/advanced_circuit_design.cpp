//
// Created by Willem Lambooy on 27/01/2025.
//

#include <fiction/algorithms/simulation/sidb/compare_by_ground_state_isolation.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/technology/sidb_on_the_fly_gate_library.hpp>
#include <fiction/technology/sidb_on_the_fly_mini_gate_library.hpp>
#include <fiction/technology/sidb_skeleton_bestagon_library.hpp>
#include <fiction/technology/sidb_skeleton_bestagon_mini_library.hpp>
#if (FICTION_Z3_SOLVER)

#include "fiction_experiments.hpp"

#include <fiction/algorithms/network_transformation/technology_mapping.hpp>
#include <fiction/algorithms/physical_design/advanced_circuit_design.hpp>
#include <fiction/algorithms/physical_design/design_sidb_gates.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/io/read_sidb_surface_defects.hpp>
#include <fiction/layouts/bounding_box.hpp>
#include <fiction/technology/area.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/sidb_defect_surface.hpp>
#include <fiction/technology/sidb_defects.hpp>
#include <fiction/technology/sidb_on_the_fly_mini_gate_library.hpp>
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
#include <string>

// This script conducts defect-aware placement and routing with defect-aware on-the-fly SiDB gate design. Thereby, SiDB
// circuits can be designed in the presence of atomic defects.

// This algorithm was proposed in \"On-the-fly Defect-Aware Design of Circuits based on Silicon Dangling Bond Logic\" by
// J. Drewniok, M. Walter, S. S. H. Ng, K. Walus, and R. Wille in IEEE NANO 2024
// (https://ieeexplore.ieee.org/abstract/document/10628962).

#define USE_MINI

int main(int argc, char* argv[])  // NOLINT
{
    using gate_lyt = fiction::hex_even_row_gate_clk_lyt;
    using cell_lyt = fiction::sidb_cell_clk_lyt_cube;

#ifdef USE_MINI
    using gate_lib          = fiction::sidb_on_the_fly_mini_gate_library;
    using skeleton_gate_lib = fiction::sidb_skeleton_bestagon_mini_library;
#else
    using gate_lib          = fiction::sidb_on_the_fly_gate_library;
    using skeleton_gate_lib = fiction::sidb_skeleton_bestagon_library;
#endif

    fiction::sidb_defect_surface<cell_lyt> lyt{};

    fiction::design_sidb_gates_params<fiction::sidb_defect_surface<cell_lyt>> design_gate_params{};

    design_gate_params.operational_params.simulation_parameters = fiction::sidb_simulation_parameters{2, -0.32};
    design_gate_params.operational_params.op_condition_positive_charges =
        fiction::is_operational_params<fiction::cell<fiction::sidb_defect_surface<cell_lyt>>>::
            operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
    design_gate_params.operational_params.op_condition_kinks = fiction::is_operational_params<
        fiction::cell<fiction::sidb_defect_surface<cell_lyt>>>::operational_condition_kinks::REJECT_KINKS;

    design_gate_params.post_design_process = {
        std::make_unique<fiction::compare_by_minimum_ground_state_isolation<fiction::sidb_defect_surface<cell_lyt>>>(),
        std::make_unique<fiction::compare_by_average_ground_state_isolation<fiction::sidb_defect_surface<cell_lyt>>>()};

    // needs to be changed if a different skeleton is used.
    if constexpr (std::is_same_v<gate_lib, fiction::sidb_on_the_fly_mini_gate_library>)
    {
        design_gate_params.canvas = {{14, 10}, {21, 19}};  // smaller canvas
    }
    else
    {
        design_gate_params.canvas = {{24, 17}, {34, 28}};  // normal canvas
    }

    design_gate_params.number_of_sidbs               = 4;
    design_gate_params.operational_params.sim_engine = fiction::sidb_simulation_engine::CLUSTERCOMPLETE;
    design_gate_params.termination_cond              = fiction::design_sidb_gates_params<
                     fiction::sidb_defect_surface<cell_lyt>>::termination_condition::OBTAINED_N_SOLUTIONS;
    // design_gate_params.max_num_solutions = 200;
    design_gate_params.design_mode = fiction::design_sidb_gates_params<
        fiction::sidb_defect_surface<cell_lyt>>::design_sidb_gates_mode::EXHAUSTIVE_GATE_DESIGNER;

    // // save atomic defects which their respective physical parameters as experimentally determined by T. R. Huff, T.
    // // Dienel, M. Rashidi, R. Achal, L. Livadaru, J. Croshaw, and R. A. Wolkow, "Electrostatic landscape of a
    // // Hydrogen-terminated Silicon Surface Probed by a Moveable Quantum Dot."
    // const auto stray_db   = fiction::sidb_defect{fiction::sidb_defect_type::DB, -1, 4.1, 1.8};
    // const auto si_vacancy = fiction::sidb_defect{fiction::sidb_defect_type::SI_VACANCY, -1, 10.6, 5.9};

    static const std::string layouts_folder =
        fmt::format("{}/physical_design_with_on_the_fly_gate_design/layouts", EXPERIMENTS_PATH);

    // // read-in the initial defects. Physical parameters of the defects are not stored yet.
    // auto surface_lattice_initial = fiction::read_sidb_surface_defects<cell_lyt>(
    //     "../../experiments/physical_design_with_on_the_fly_gate_design/1_percent_with_charged_surface.txt");

    // create an empty surface.
    fiction::sidb_defect_surface<cell_lyt> surface_lattice{};

    // // add physical parameters of the defects to the surface_lattice.
    // surface_lattice_initial.foreach_sidb_defect(
    //     [&surface_lattice, &stray_db, &si_vacancy](const auto& cd)
    //     {
    //         if (cd.second.type == fiction::sidb_defect_type::DB)
    //         {
    //             surface_lattice.assign_sidb_defect(cd.first, stray_db);
    //         }
    //         else if (cd.second.type == fiction::sidb_defect_type::SI_VACANCY)
    //         {
    //             surface_lattice.assign_sidb_defect(cd.first, si_vacancy);
    //         }
    //         else
    //         {
    //             surface_lattice.assign_sidb_defect(cd.first, cd.second);
    //         }
    //     });

    // determine bounding-box of the surface to set the aspect ratio of the surface lattice.
    const auto bb_defect_surface = fiction::bounding_box_2d{surface_lattice};
    surface_lattice.resize(bb_defect_surface.get_max());

    const auto lattice_tiling = gate_lyt{{11, 30}};

    experiments::experiment<std::string, double, uint64_t, bool> sidb_circuits_with_defects{
        "sidb_circuits_with_defects", "benchmark", "runtime", "number of aspect ratios", "equivalent"};

    // constexpr const uint64_t bench_select = fiction_experiments::fo;
    // constexpr const uint64_t bench_select = fiction_experiments::and3;
    // constexpr const uint64_t bench_select = fiction_experiments::xor3;
    constexpr const uint64_t bench_select = fiction_experiments::xor5_maj;
    // fiction_experiments::all & ~fiction_experiments::parity & ~fiction_experiments::two_bit_add_maj &
    // ~fiction_experiments::b1_r2 & ~fiction_experiments::clpl & ~fiction_experiments::iscas85 &
    // ~fiction_experiments::epfl & ~fiction_experiments::half_adder & ~fiction_experiments::full_adder &
    // ~fiction_experiments::one_bit_add_aoig & ~fiction_experiments::one_bit_add_maj & ~fiction_experiments::cm82a_5;

    const auto create_and3_tt = []
    {
        static constexpr const char* and3_tt_string = "10000000";  // Output 1
        kitty::dynamic_truth_table   and3_tb{3};
        kitty::create_from_binary_string(and3_tb, and3_tt_string);
        return std::vector{and3_tb};
    };

    const auto create_xor3_tt = []
    {
        static constexpr const char* xor3_tt_string = "10010110";  // Output 1
        kitty::dynamic_truth_table   xor3_tb{3};
        kitty::create_from_binary_string(xor3_tb, xor3_tt_string);
        return std::vector{xor3_tb};
    };

    const auto create_xor5_tt = []
    {
        static constexpr const char* xor5_tt_string = "96699669";  // Output 1
        kitty::dynamic_truth_table   xor5_tb{5};
        kitty::create_from_hex_string(xor5_tb, xor5_tt_string);
        return std::vector{xor5_tb};
    };

    std::map<std::string, std::vector<kitty::dynamic_truth_table>> tt_map{{"TEST/fo", fiction::create_fan_out_tt()},
                                                                          {"TEST/and3", create_and3_tt()},
                                                                          {"TEST/xor3", create_xor3_tt()},
                                                                          {"fontes18/xor5_r1", create_xor5_tt()},
                                                                          {"fontes18/xor5Maj", create_xor5_tt()}};

    for (const auto& benchmark : fiction_experiments::all_benchmarks(bench_select))
    {
        fmt::print("[attempts] processing {}\n", benchmark);
        mockturtle::xag_network xag{};

        [[maybe_unused]] const auto read_verilog_result =
            lorina::read_verilog(fiction_experiments::benchmark_path(benchmark), mockturtle::verilog_reader(xag));
        assert(read_verilog_result == lorina::return_code::success);

        // compute depth
        const mockturtle::depth_view depth_xag{xag};

        const fiction::technology_mapping_params tech_map_params = fiction::all_2_input_functions();

        // parameters for cut rewriting
        mockturtle::cut_rewriting_params cut_params{};
        cut_params.cut_enumeration_ps.cut_size = 4;

        const mockturtle::xag_npn_resynthesis<mockturtle::xag_network,                    // the input network type
                                              mockturtle::xag_network,                    // the database network type
                                              mockturtle::xag_npn_db_kind::xag_complete>  // the kind of database to use

            resynthesis_function{};

        // rewrite network cuts using the given re-synthesis function
        const auto cut_xag = mockturtle::cut_rewriting(xag, resynthesis_function, cut_params);

        // perform technology mapping
        const auto mapped_network = fiction::technology_mapping(cut_xag, tech_map_params);

        fiction::advanced_circuit_design_params<cell_lyt> params{};

        design_gate_params.max_num_solutions = 200;
        params.num_trials                    = 30;
        params.selectivity                   = 0.93;

        params.spec = tt_map.at(benchmark);

        if (argc != 1)
        {
            design_gate_params.max_num_solutions = std::stoi(argv[1]);
            params.num_trials                    = std::stoi(argv[2]);
            params.selectivity                   = std::stod(argv[3]);
            params.num_trials_for_global_scope   = std::stoi(argv[4]);
            params.selectivity_for_global_scope  = std::stod(argv[5]);
        }

        params.exact_design_parameters.scheme        = "ROW4";
        params.exact_design_parameters.crossings     = true;
        params.exact_design_parameters.border_io     = false;
        params.exact_design_parameters.desynchronize = true;
        params.exact_design_parameters.upper_bound_x = 11;         // 12 x 31 tiles
        params.exact_design_parameters.upper_bound_y = 30;         // 12 x 31 tiles
        params.exact_design_parameters.timeout       = 3'600'000;  // 1h in ms

        params.sidb_on_the_fly_gate_library_parameters.defect_surface                = surface_lattice;
        params.sidb_on_the_fly_gate_library_parameters.use_skeleton_influence_bounds = true;
        params.sidb_on_the_fly_gate_library_parameters.design_gate_params            = design_gate_params;
        params.sidb_on_the_fly_gate_library_parameters.canvas_sidb_complex_gates = // 5;//
            params.sidb_on_the_fly_gate_library_parameters.design_gate_params.number_of_sidbs;

        fiction::advanced_circuit_design_stats<gate_lyt> st{};

        lyt =
            fiction::advanced_circuit_design<decltype(mapped_network), cell_lyt, gate_lyt, gate_lib, skeleton_gate_lib>(
                mapped_network, lattice_tiling, params, &st);

        params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.print = true;

        std::cout << "\nassessing operational status for each input combination of the generated circuit..."
                  << std::endl;

        if (is_operational(lyt, params.spec,
                           params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params)
                .status == fiction::operational_status::OPERATIONAL)
        {
            std::cout << "\n\nCIRCUIT OPERATION VERIFICATION COMPLETED: PASS" << std::endl;
        }
        else
        {
            std::cout << "\n\nCIRCUIT OPERATION VERIFICATION COMPLETED: FAILED" << std::endl;
        }

        // check equivalence
        const auto miter = mockturtle::miter<mockturtle::klut_network>(mapped_network, st.gate_layout.value());
        const auto eq    = mockturtle::equivalence_checking(*miter);
        assert(eq.has_value());

        // determine bounding box and exclude atomic defects
        const auto bb = fiction::bounding_box_2d<cell_lyt>(static_cast<cell_lyt>(lyt));

        // compute area
        fiction::area_stats                            area_stats{};
        fiction::area_params<fiction::sidb_technology> area_ps{};
        fiction::area(bb, area_ps, &area_stats);

        sidb_circuits_with_defects(benchmark, mockturtle::to_seconds(st.time_total), st.exact_stats.num_aspect_ratios,
                                   *eq);
        sidb_circuits_with_defects.save();
        sidb_circuits_with_defects.table();

        // write a SiQAD simulation file
        // fiction::write_sqd_layout(result, fmt::format("{}/{}.sqd", layouts_folder, benchmark));
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

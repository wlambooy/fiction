//
// Created by Willem Lambooy on 27/01/2025.
//

#if (FICTION_Z3_SOLVER)

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
#include <string>

// This script conducts defect-aware placement and routing with defect-aware on-the-fly SiDB gate design. Thereby, SiDB
// circuits can be designed in the presence of atomic defects.

// This algorithm was proposed in \"On-the-fly Defect-Aware Design of Circuits based on Silicon Dangling Bond Logic\" by
// J. Drewniok, M. Walter, S. S. H. Ng, K. Walus, and R. Wille in IEEE NANO 2024
// (https://ieeexplore.ieee.org/abstract/document/10628962).

// #define USE_MINI

namespace fs = std::filesystem;

int main(int argc, char* argv[])  // NOLINT
{
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
    design_gate_params.canvas = {{5, 12}, {21, 19}};  // new_mini_bestagon.sqd
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

            // write_

            advanced_circuit_design_params<lyt_t> params{};

            // design_gate_params.max_num_solutions = 200;
            // params.num_trials                   = 100;
            // params.selectivity                  = 0.93;

            if (argc == 3)
            {
                design_gate_params.number_of_canvas_sidbs = std::stoull(argv[1]);
            }
            else if (argc == 7)
            {
                design_gate_params.number_of_canvas_sidbs      = std::stoull(argv[1]);
                design_gate_params.maximum_number_of_solutions = std::stoull(argv[2]);
                params.num_trials                              = std::stoull(argv[3]);
                params.quantization_factor                     = std::stod(argv[4]);
                params.selectivity                             = std::stod(argv[5]);
                params.success_rate_ceiling                    = std::stod(argv[6]);
            }
            else if (argc == 8)
            {
                design_gate_params.number_of_canvas_sidbs      = std::stoull(argv[1]);
                design_gate_params.maximum_number_of_solutions = std::stoull(argv[2]);
                params.num_trials                              = std::stoull(argv[3]);
                params.quantization_factor                     = std::stod(argv[4]);
                params.selectivity                             = std::stod(argv[5]);
                params.success_rate_ceiling                    = std::stod(argv[6]);

                if (strncmp(argv[7], "e", 1) == 0)
                {
                    design_gate_params.design_mode =
                        design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::EXHAUSTIVE;
                }
                else
                {
                    params.available_threads = std::stoull(argv[7]);
                }
            }
            else if (argc == 9)
            {
                design_gate_params.number_of_canvas_sidbs      = std::stoull(argv[1]);
                design_gate_params.maximum_number_of_solutions = std::stoull(argv[2]);
                params.num_trials                              = std::stoull(argv[3]);
                params.quantization_factor                     = std::stod(argv[4]);
                params.selectivity                             = std::stod(argv[5]);
                params.success_rate_ceiling                    = std::stod(argv[6]);
                if (strncmp(argv[7], "e", 1) == 0)
                {
                    design_gate_params.design_mode =
                        design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::EXHAUSTIVE;
                    params.available_threads = std::stoull(argv[8]);
                }
                else
                {
                    params.available_threads             = std::stoull(argv[7]);
                    design_gate_params.available_threads = std::stoull(argv[8]);
                }
            }
            else if (argc == 10)
            {
                design_gate_params.number_of_canvas_sidbs      = std::stoull(argv[1]);
                design_gate_params.maximum_number_of_solutions = std::stoull(argv[2]);
                params.num_trials                              = std::stoull(argv[3]);
                params.quantization_factor                     = std::stod(argv[4]);
                params.selectivity                             = std::stod(argv[5]);
                params.success_rate_ceiling                    = std::stod(argv[6]);
                params.available_threads                       = std::stoull(argv[8]);
                design_gate_params.available_threads           = std::stoull(argv[9]);

                design_gate_params.design_mode = design_sidb_gates_params<lyt_t>::design_sidb_gates_mode::EXHAUSTIVE;
            }
            else if (argc != 1)  // todo
            {
                design_gate_params.maximum_number_of_solutions = std::stoull(argv[1]);
                params.num_trials                              = std::stoull(argv[2]);
                params.selectivity                             = std::stod(argv[3]);
            }

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

            // params.sidb_on_the_fly_gate_library_parameters.design_gate_params_complex_gates =
            // 6;  //
            // params.sidb_on_the_fly_gate_library_parameters.design_gate_params.number_of_sidbs;

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

            // params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params.print = true;
            //
            // std::cout << "\nassessing operational status for each input combination of the generated circuit..."
            //           << std::endl;
            //
            // if (is_operational(*lyt, params.spec,
            //                    params.sidb_on_the_fly_gate_library_parameters.design_gate_params.operational_params)
            //         .status != operational_status::OPERATIONAL)
            // {
            //     std::cout << "\n\nCIRCUIT OPERATION VERIFICATION COMPLETED: FAILED" << std::endl;
            //
            //     return EXIT_FAILURE;
            // }
            //
            // std::cout << "\n\nCIRCUIT OPERATION VERIFICATION COMPLETED: PASS" << std::endl;

            // check equivalence
            const auto miter = mockturtle::miter<mockturtle::klut_network>(mapped_network, st.gate_layout.value());
            // const auto eq    = mockturtle::equivalence_checking(*miter);
            // assert(eq.has_value());/

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

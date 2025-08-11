//
// Created by Willem Lambooy on 23/06/2025.
//

#include <fiction/algorithms/physical_design/advanced_circuit_design.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/layouts/hexagonal_layout.hpp>
#include <fiction/technology/sidb_skeleton_bestagon_mini_library.hpp>
#include <fiction/types.hpp>

#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/networks/klut.hpp>

#include <algorithm>
#include <iostream>
#include <optional>
#include <random>

using namespace fiction;

using cell_lyt          = sidb_100_cell_clk_lyt_cube;
using gate_lyt          = hex_even_row_gate_clk_lyt;
using gate_lib          = sidb_on_the_fly_mini_gate_library;
using skeleton_gate_lib = sidb_skeleton_bestagon_mini_library;

void print_skeleton_gate_layout(const gate_lyt&                                 gate_layout,
                                const advanced_circuit_design_params<cell_lyt>& params) noexcept
{
    cell_lyt lyt = apply_gate_library<cell_lyt, skeleton_gate_lib, gate_lyt>(gate_layout);

    gate_layout.foreach_node(
        [&](const auto& nn)
        {
            if (skip_physical_design_for_node(gate_layout, nn))
            {
                return;
            }
            const auto& t = gate_layout.get_tile(nn);
            const auto  canvas =
                is_complex_gate<gate_lyt>(gate_layout, nn) ?
                     params.sidb_on_the_fly_gate_library_parameters.design_gate_params_complex_gates.canvas :
                     params.sidb_on_the_fly_gate_library_parameters.design_gate_params.canvas;

            for (const cell<cell_lyt>& relative_c : all_coordinates_in_spanned_area(canvas.first, canvas.second))
            {
                const cell<cell_lyt> absolute_c =
                    relative_to_absolute_cell_position<skeleton_gate_lib::gate_x_size(),
                                                       skeleton_gate_lib::gate_y_size(), gate_lyt, cell_lyt>(
                        gate_layout, t, relative_c);
                lyt.assign_cell_type(absolute_c, sidb_technology::cell_type::LOGIC);
            }
        });

    std::cout << "Skeleton looks like:" << std::endl;
    print_layout(lyt);
    std::cout << std::endl;
}

[[nodiscard]] kitty::dynamic_truth_table generate_random_truth_table(const uint32_t num_vars) noexcept
{
    const auto num_bits = static_cast<size_t>(1) << num_vars;

    // Use random_device to seed the Mersenne Twister engine
    static std::random_device              rd;
    static std::mt19937                    gen(rd());
    static std::uniform_int_distribution<> dist(0, 1);

    std::string tt_string;
    tt_string.reserve(num_bits);

    while (tt_string.empty() ||
           std::all_of(tt_string.cbegin(), tt_string.cend(), [&](const char c) { return c == tt_string.front(); }))
    {
        tt_string.clear();
        for (size_t i = 0; i < num_bits; ++i)
        {
            tt_string += dist(gen) ? '1' : '0';
        }
    }

    kitty::dynamic_truth_table tt{num_vars};
    std::cout << tt_string << std::endl;
    create_from_binary_string(tt, tt_string);
    return tt;
}

[[nodiscard]] gate_lyt generate_gate_layout_type_1() noexcept
{
    gate_lyt gate_layout{{2, 4}, row_clocking<gate_lyt>()};

    const auto pi1 = gate_layout.create_pi("A", {0, 0});
    const auto pi2 = gate_layout.create_pi("B", {2, 0});

    const auto i1o1_1 = gate_layout.create_node({pi1}, generate_random_truth_table(1), {1, 1});
    const auto i1o1_2 = gate_layout.create_node({pi2}, generate_random_truth_table(1), {2, 1});

    const auto i2o1 = gate_layout.create_node({i1o1_1, i1o1_2}, generate_random_truth_table(2), {1, 2});

    const auto i1o1_3 = gate_layout.create_node({i2o1}, generate_random_truth_table(1), {1, 3});

    gate_layout.create_po(i1o1_3, "O", {0, 4});

    return gate_layout;
}

[[nodiscard]] gate_lyt generate_gate_layout_type_2() noexcept
{
    gate_lyt gate_layout{{2, 4}, row_clocking<gate_lyt>()};

    const auto pi1 = gate_layout.create_pi("A", {2, 0});
    const auto pi2 = gate_layout.create_pi("B", {0, 1});

    const auto i1o1_1 = gate_layout.create_node({pi1}, generate_random_truth_table(1), {2, 1});
    const auto i1o1_2 = gate_layout.create_node({pi2}, generate_random_truth_table(1), {0, 2});
    const auto i1o1_3 = gate_layout.create_node({i1o1_1}, generate_random_truth_table(1), {1, 2});

    const auto i2o1 = gate_layout.create_node({i1o1_2, i1o1_3}, generate_random_truth_table(2), {1, 3});

    gate_layout.create_po(i2o1, "O", {1, 4});

    return gate_layout;
}

[[nodiscard]] gate_lyt generate_gate_layout_type_3() noexcept
{
    gate_lyt gate_layout{{2, 4, 1}, row_clocking<gate_lyt>()};

    const auto pi1 = gate_layout.create_pi("A", {1, 0});
    const auto pi2 = gate_layout.create_pi("B", {1, 0, 1});

    const auto i1o1_1 = gate_layout.create_node({pi1}, generate_random_truth_table(1), {1, 1});
    const auto i1o1_2 = gate_layout.create_node({pi2}, generate_random_truth_table(1), {2, 1});

    const auto i2o1 = gate_layout.create_node({i1o1_1, i1o1_2}, generate_random_truth_table(2), {1, 2});

    const auto i1o1_3 = gate_layout.create_node({i2o1}, generate_random_truth_table(1), {2, 3});

    gate_layout.create_po(i1o1_3, "O", {1, 4});

    return gate_layout;
}

[[nodiscard]] gate_lyt generate_gate_layout_type_4() noexcept
{
    gate_lyt gate_layout{{3, 4}, row_clocking<gate_lyt>()};

    const auto pi1 = gate_layout.create_pi("A", {2, 0});
    const auto pi2 = gate_layout.create_pi("B", {3, 1});

    const auto i1o1_1 = gate_layout.create_node({pi1}, generate_random_truth_table(1), {2, 1});
    const auto i1o1_2 = gate_layout.create_node({i1o1_1}, generate_random_truth_table(1), {1, 2});
    const auto i1o1_3 = gate_layout.create_node({pi2}, generate_random_truth_table(1), {2, 2});

    const auto i2o1 = gate_layout.create_node({i1o1_2, i1o1_3}, generate_random_truth_table(2), {2, 3});

    gate_layout.create_po(i2o1, "O", {2, 4});

    return gate_layout;
}


[[nodiscard]] std::optional<gate_lyt> generate_supertile_gate_layout_for_tt() noexcept
{

}

int main(int argc, char* argv[])
{
    /// DESIGN GATE PARAMS

    design_sidb_gates_params<cell_lyt> design_gate_params{};

    design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{3, -0.32};
    design_gate_params.operational_params.op_condition_positive_charges =
        is_operational_params<cell<cell_lyt>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
    design_gate_params.operational_params.op_condition_kinks =
        is_operational_params<cell<cell_lyt>>::operational_condition_kinks::REJECT_KINKS;
    design_gate_params.design_mode = design_sidb_gates_params<cell_lyt>::design_sidb_gates_mode::RANDOM;

    // design_gate_params.canvas = {{15, 10}, {26, 19}};  // smaller canvas
    design_gate_params.canvas = {{15, 8}, {29, 17}};  // smaller canvas

    design_gate_params.number_of_canvas_sidbs        = 4;
    design_gate_params.operational_params.sim_engine = sidb_simulation_engine::CLUSTERCOMPLETE;
    design_gate_params.termination_cond =
        design_sidb_gates_params<cell_lyt>::termination_condition::OBTAINED_N_SOLUTIONS;
    design_gate_params.maximum_number_of_solutions = 200;

    /// COMPLEX DESIGN GATE PARAMS

    design_sidb_gates_params<cell_lyt> design_gate_params_complex_gates{};

    design_gate_params_complex_gates.operational_params.simulation_parameters = sidb_simulation_parameters{3, -0.32};
    design_gate_params_complex_gates.operational_params.op_condition_positive_charges =
        is_operational_params<cell<cell_lyt>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
    design_gate_params_complex_gates.operational_params.op_condition_kinks =
        is_operational_params<cell<cell_lyt>>::operational_condition_kinks::REJECT_KINKS;
    design_gate_params_complex_gates.design_mode = design_sidb_gates_params<cell_lyt>::design_sidb_gates_mode::RANDOM;

    // design_gate_params_complex_gates.canvas = {{13, 10}, {28, 19}};
    design_gate_params_complex_gates.canvas = {{13, 8}, {31, 17}};

    design_gate_params_complex_gates.number_of_canvas_sidbs        = 4;
    design_gate_params_complex_gates.operational_params.sim_engine = sidb_simulation_engine::CLUSTERCOMPLETE;
    design_gate_params_complex_gates.termination_cond =
        design_sidb_gates_params<cell_lyt>::termination_condition::OBTAINED_N_SOLUTIONS;
    design_gate_params_complex_gates.maximum_number_of_solutions = 200;

    /// CIRCUIT DESIGN PARAMS

    advanced_circuit_design_params<cell_lyt> params{};

    params.num_trials                   = 100;
    params.selectivity                  = 0.93;
    params.num_trials_for_global_scope  = 4;
    params.selectivity_for_global_scope = 0.9;

    if (argc == 2)
    {
        design_gate_params_complex_gates.number_of_canvas_sidbs = std::stoull(argv[1]);
    }
    else if (argc == 3)
    {
        design_gate_params.number_of_canvas_sidbs               = std::stoull(argv[1]);
        design_gate_params_complex_gates.number_of_canvas_sidbs = std::stoull(argv[2]);
    }
    else if (argc == 15)
    {
        design_gate_params.number_of_canvas_sidbs                    = std::stoull(argv[1]);
        design_gate_params_complex_gates.number_of_canvas_sidbs      = std::stoull(argv[2]);
        design_gate_params.maximum_number_of_solutions               = std::stoull(argv[3]);
        design_gate_params_complex_gates.maximum_number_of_solutions = std::stoull(argv[4]);
        params.num_trials                                            = std::stoull(argv[5]);
        params.quantization_factor                                   = std::stod(argv[6]);
        params.selectivity                                           = std::stod(argv[7]);
        params.num_trials_for_double_scope                           = std::stoull(argv[8]);
        params.quantization_factor_for_double_scope                  = std::stod(argv[9]);
        params.selectivity_for_double_scope                          = std::stod(argv[10]);
        params.num_trials_for_global_scope                           = std::stoull(argv[11]);
        params.quantization_factor_for_global_scope                  = std::stod(argv[12]);
        params.selectivity_for_global_scope                          = std::stod(argv[13]);
        params.excited_state_alpha                                   = std::stod(argv[14]);
    }
    else if (argc == 16)
    {
        design_gate_params.number_of_canvas_sidbs                    = std::stoull(argv[1]);
        design_gate_params_complex_gates.number_of_canvas_sidbs      = std::stoull(argv[2]);
        design_gate_params.maximum_number_of_solutions               = std::stoull(argv[3]);
        design_gate_params_complex_gates.maximum_number_of_solutions = std::stoull(argv[4]);
        params.num_trials                                            = std::stoull(argv[5]);
        params.quantization_factor                                   = std::stod(argv[6]);
        params.selectivity                                           = std::stod(argv[7]);
        params.num_trials_for_double_scope                           = std::stoull(argv[8]);
        params.quantization_factor_for_double_scope                  = std::stod(argv[9]);
        params.selectivity_for_double_scope                          = std::stod(argv[10]);
        params.num_trials_for_global_scope                           = std::stoull(argv[11]);
        params.quantization_factor_for_global_scope                  = std::stod(argv[12]);
        params.selectivity_for_global_scope                          = std::stod(argv[13]);
        params.excited_state_alpha                                   = std::stod(argv[14]);
        params.available_threads                                     = std::stoull(argv[15]);
    }
    else if (argc != 1)
    {
        design_gate_params.maximum_number_of_solutions               = std::stoull(argv[1]);
        design_gate_params_complex_gates.maximum_number_of_solutions = std::stoull(argv[2]);
        params.num_trials                                            = std::stoull(argv[3]);
        params.selectivity                                           = std::stod(argv[4]);
        params.num_trials_for_global_scope                           = std::stoull(argv[5]);
        params.selectivity_for_global_scope                          = std::stod(argv[6]);
    }

    params.sidb_on_the_fly_gate_library_parameters.design_gate_params               = design_gate_params;
    params.sidb_on_the_fly_gate_library_parameters.design_gate_params_complex_gates = design_gate_params_complex_gates;
    advanced_circuit_design_stats<gate_lyt> st{};

    const uint64_t number = 500;

    for (uint64_t i = 0; i < number; ++i)
    {
        // Use random_device to seed the Mersenne Twister engine
        static std::random_device              rd;
        static std::mt19937                    gen(rd());
        static std::uniform_int_distribution<> dist(0, 3);

        const gate_lyt& gate_layout =
            (*std::vector<gate_lyt (*)()>{generate_gate_layout_type_1, generate_gate_layout_type_2,
                                          generate_gate_layout_type_3, generate_gate_layout_type_4}
                  .at(dist(gen)))();

        print_skeleton_gate_layout(gate_layout, params);

        std::optional<cell_lyt> lyt{};

        mockturtle::stopwatch<>::duration time_counter{};
        {
            const mockturtle::stopwatch stop{time_counter};

            lyt = fiction::advanced_circuit_design<mockturtle::klut_network, cell_lyt, gate_lyt, gate_lib,
                                                   skeleton_gate_lib>(std::nullopt, gate_layout, params, &st);
        }

        if (!lyt.has_value())
        {
            std::cout << "RESULT: "
                         "FAILURE\tXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                         "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX\nFailed runtime: "
                      << mockturtle::to_seconds(time_counter) << " s" << std::endl;
            continue;
        }

        std::cout << "RESULT: "
                     "SUCCESS\tOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\nOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"
                     "OOOOOOOOOOOOOOOOOOOOOOOOOOOOO\nRuntime: "
                  << mockturtle::to_seconds(time_counter) << " s" << std::endl;
    }

    return EXIT_SUCCESS;
}
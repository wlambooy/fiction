//
// Created by wille on 17/10/2023.
//

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/physical_design/apply_gate_library.hpp>
#include <fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/technology/sidb_bestagon_library.hpp>
#include <fiction/utils/layout_utils.hpp>

using namespace fiction;

TEST_CASE("Benchmark ExGS-EF simulations of Bestagon gates", "[benchmark]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{0, 0}};

    layout.create_or({}, {}, {0, 0});

    sidb_simulation_parameters params{2, -0.32};

    //    params.ef_params = {1.3831063607828464, 1.1793595882585919, 0.14006305601323665, -0.30951723965575634};

    //    params.ef_params = {0.6, 1.15, 0.03, 0};

    params.ef_params = {2.4741532309861984, 0.7261618183793207, 0.041491891025909895, 0};

    BENCHMARK("ExGS-EF simulation of a Bestagon OR-gate (base 2)")
    {
        charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> or_gate{
            convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)), params,
            sidb_charge_state::NEGATIVE};

        return exhaustive_ground_state_simulation(or_gate, params);
    };

    layout.clear_tile({0, 0});
    layout.create_pi("x1", {0, 0});

    params.base = 2;

    BENCHMARK("ExGS-EF simulation of a Bestagon wire (base 2)")
    {
        charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> wire{
            convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)), params,
            sidb_charge_state::NEGATIVE};

        return exhaustive_ground_state_simulation(wire, params);
    };

    params.base = 3;

    BENCHMARK("ExGS-EF simulation of a Bestagon wire (base 3)")
    {
        charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> wire{
            convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)), params,
            sidb_charge_state::NEGATIVE};

        return exhaustive_ground_state_simulation(wire, params);
    };
}

TEST_CASE("Reference benchmark ExGS simulations of Bestagon gates", "[benchmark]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{0, 0}};

    layout.create_or({}, {}, {0, 0});

    sidb_simulation_parameters params{2, -0.32};

    BENCHMARK("ExGS simulation of a Bestagon OR-gate (base 2)")
    {
        charge_distribution_surface<sidb_cell_clk_lyt_siqad, false> or_gate{
            convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)), params,
            sidb_charge_state::NEGATIVE};

        return exhaustive_ground_state_simulation(or_gate, params);
    };

    layout.clear_tile({0, 0});
    layout.create_pi("x1", {0, 0});

    params.base = 2;

    BENCHMARK("ExGS simulation of a Bestagon wire (base 2)")
    {
        charge_distribution_surface<sidb_cell_clk_lyt_siqad, false> wire{
            convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)), params,
            sidb_charge_state::NEGATIVE};

        return exhaustive_ground_state_simulation(wire, params);
    };

    params.base = 3;

    BENCHMARK("ExGS simulation of a Bestagon wire (base 3)")
    {
        charge_distribution_surface<sidb_cell_clk_lyt_siqad, false> wire{
            convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)), params,
            sidb_charge_state::NEGATIVE};

        return exhaustive_ground_state_simulation(wire, params);
    };
}



//ef_params = {1.0812582466265876, 1.4991189374313565, 0.11185634057430065, -0.36101092529278067};
//[I 2023-10-18 00:35:47,804] Trial 2874 finished with value: 0.09209537506103516 and parameters: {'glob_theta': 1.0812582466265876, 'loc_theta': 1.4991189374313565, 'stable_err': 0.11185634057430065, 'energy_err': -0.36101092529278067}. Best is trial 2874 with value: 0.09209537506103516.


//[I 2023-10-18 00:57:56,141] Trial 1374 finished with value: 0.09522557258605957 and parameters: {'glob_theta': 0.8639610673297877, 'loc_theta': 1.4952889834625438, 'stable_err': 0.0051128341500607495, 'energy_err': -0.4979215555936839}. Best is trial 1374 with value: 0.09522557258605957.

//[I 2023-10-18 11:19:50,988] Trial 32020 finished with value: 0.09145760536193848 and parameters: {'glob_theta': 2.4741532309861984, 'loc_theta': 0.7261618183793207, 'stable_err': 0.041491891025909895, 'energy_err': -0.38317175145482807}. Best is trial 32020 with value: 0.09145760536193848.
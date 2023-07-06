//
// Created by Willem Lambooy on 15/05/2023.
//

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/kink_finder.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/gate_level_layout.hpp>
#include <fiction/utils/math_utils.hpp>

using namespace fiction;

using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

TEST_CASE("Single-threaded OR of two ANDs (Bestagon)", "[kink_finder]")
{
    hex_gate_lyt gates{{3, 4}, row_clocking<hex_gate_lyt>()};

    const auto x1 = gates.create_pi("x1", {0, 0});
    const auto x2 = gates.create_pi("x2", {1, 0});
    const auto x3 = gates.create_pi("x3", {2, 0});
    const auto x4 = gates.create_pi("x4", {3, 0});

    const auto a1 = gates.create_and(x1, x2, {0, 1});
    const auto a2 = gates.create_and(x3, x4, {2, 1});

    const auto b1 = gates.create_buf(a1, {1, 2});
    const auto b2 = gates.create_buf(a2, {2, 2});

    const auto o1 = gates.create_or(b1, b2, {1, 3});

    gates.create_po(o1, "f1", {1, 4});

    const expected_gs_params params{sidb_simulation_parameters{2, -0.32}, -1, {-0.05}, std::nullopt, 1};

    const auto simulation_results = kink_finder<sidb_bestagon_library_wrapper>(gates, params);

    CHECK(simulation_results.charge_distributions.size() == 1);
    CHECK(round_to_n_decimal_places(simulation_results.charge_distributions.front().get_system_energy(), 6) ==
          7.349980);
    CHECK(simulation_results.charge_distributions.front().is_physically_valid());
    CHECK(simulation_results.additional_simulation_parameters.empty());
    CHECK(simulation_results.algorithm_name == "kink_finder");
    CHECK(simulation_results.simulation_runtime.count() > 0);
}

TEST_CASE("OR of two ANDs (Bestagon)", "[kink_finder]")
{
    hex_gate_lyt gates{{3, 4}, row_clocking<hex_gate_lyt>()};

    const auto x1 = gates.create_pi("x1", {0, 0});
    const auto x2 = gates.create_pi("x2", {1, 0});
    const auto x3 = gates.create_pi("x3", {2, 0});
    const auto x4 = gates.create_pi("x4", {3, 0});

    const auto a1 = gates.create_and(x1, x2, {0, 1});
    const auto a2 = gates.create_and(x3, x4, {2, 1});

    const auto b1 = gates.create_buf(a1, {1, 2});
    const auto b2 = gates.create_buf(a2, {2, 2});

    const auto o1 = gates.create_or(b1, b2, {1, 3});

    gates.create_po(o1, "f1", {1, 4});

    const expected_gs_params params{sidb_simulation_parameters{2, -0.32}, -1, {-0.05}, std::nullopt};

    const auto simulation_results = kink_finder<sidb_bestagon_library_wrapper>(gates, params);

    REQUIRE(simulation_results.charge_distributions.size() == 1);
    CHECK(simulation_results.charge_distributions.front().is_physically_valid());
    CHECK(simulation_results.simulation_runtime.count() > 0);
}

TEST_CASE("Gates with two POs (Bestagon)", "[kink_finder]")
{

    const expected_gs_params params{sidb_simulation_parameters{2, -0.32}, 27};

    SECTION("Hourglass")
    {
        hex_gate_lyt gates{{1, 4, 1}, row_clocking<hex_gate_lyt>()};

        const auto x2 = gates.create_pi("x2", {0, 0});

        const auto n1 = gates.create_not(x2, {0, 1});
        const auto b1 = gates.create_buf(n1, {1, 2});

        const auto x1 = gates.create_pi("x1", {0, 2});

        const auto c1 = gates.create_buf(x1, {0, 3, 0});
        const auto c2 = gates.create_buf(b1, {0, 3, 1});

        gates.create_po(c1, "f1", {0, 4});
        gates.create_po(c2, "f2", {1, 4});

        const auto simulation_results = kink_finder<sidb_bestagon_library_wrapper>(gates, params);

        CHECK(simulation_results.charge_distributions.empty());
        CHECK(simulation_results.simulation_runtime.count() > 0);
    }

    SECTION("Crossover")
    {
        hex_gate_lyt gates{{1, 3, 1}, row_clocking<hex_gate_lyt>()};

        const auto x1 = gates.create_pi("x1", {0, 0});
        const auto x2 = gates.create_pi("x2", {1, 0});

        const auto c1 = gates.create_buf(x1, {0, 1, 0});
        const auto c2 = gates.create_buf(x2, {0, 1, 1});

        const auto n1 = gates.create_not(c2, {0, 2});

        gates.create_po(n1, "f1", {0, 3});
        gates.create_po(c1, "f2", {1, 2});

        const auto simulation_results = kink_finder<sidb_bestagon_library_wrapper>(gates, params);

        REQUIRE(simulation_results.charge_distributions.size() == 1);
        CHECK(simulation_results.charge_distributions.front().is_physically_valid());
        CHECK(simulation_results.simulation_runtime.count() > 0);
    }
}

TEST_CASE("Double wire and crossover in sequence (Bestagon)", "[kink_finder]")
{
    hex_gate_lyt gates{{1, 4, 1}, row_clocking<hex_gate_lyt>()};

    const auto x1 = gates.create_pi("x1", {0, 0});
    const auto x2 = gates.create_pi("x1", {1, 0});

    const auto a1 = gates.create_buf(x1, {0, 1, 0});
    const auto a2 = gates.create_buf(x2, {0, 1, 1});
    const auto b1 = gates.create_buf(a1, {1, 2});
    const auto b2 = gates.create_buf(a2, {0, 2});
    const auto c1 = gates.create_buf(b2, {0, 3, 0});
    const auto c2 = gates.create_buf(b1, {0, 3, 1});

    gates.create_po(c1, "f1", {0, 4});
    gates.create_po(c2, "f2", {1, 4});

    const auto simulation_results =
        kink_finder<sidb_bestagon_library_wrapper>(gates, expected_gs_params{sidb_simulation_parameters{2, -0.32}});

    REQUIRE(simulation_results.charge_distributions.size() == 1);
    CHECK(simulation_results.charge_distributions.front().is_physically_valid());
    CHECK(simulation_results.simulation_runtime.count() > 0);
}

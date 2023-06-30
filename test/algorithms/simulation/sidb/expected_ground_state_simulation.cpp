//
// Created by Willem Lambooy on 15/05/2023.
//

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/expected_ground_state_simulation.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/layouts/gate_level_layout.hpp>
#include <fiction/utils/math_utils.hpp>

using namespace fiction;

using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

TEST_CASE("Single-threaded OR of two ANDs (Bestagon)", "[expect_gs]")
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

    const expected_gs_params params{sidb_simulation_parameters{2, -0.32}, -1, std::nullopt, 1};

    const auto simulation_results = expected_ground_state_simulation<sidb_bestagon_library_wrapper>(gates, params);

    CHECK(simulation_results.charge_distributions.size() == 1);
    CHECK(round_to_n_decimal_places(simulation_results.charge_distributions.front().get_system_energy(), 6) ==
          7.349980);
    CHECK(simulation_results.charge_distributions.front().is_physically_valid());
    CHECK(simulation_results.additional_simulation_parameters.empty());
    CHECK(simulation_results.algorithm_name == "expected_ground_state_simulation");
    CHECK(simulation_results.simulation_runtime.count() > 0);
}

TEST_CASE("OR of two ANDs (Bestagon)", "[expect_gs]")
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

    const expected_gs_params params{sidb_simulation_parameters{2, -0.32}, -1, std::nullopt};

    const auto simulation_results = expected_ground_state_simulation<sidb_bestagon_library_wrapper>(gates, params);

    REQUIRE(simulation_results.charge_distributions.size() == 1);
    CHECK(simulation_results.charge_distributions.front().is_physically_valid());
    CHECK(simulation_results.simulation_runtime.count() > 0);
}

TEST_CASE("Gates with two POs (Bestagon)", "[expect_gs]")
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

        const auto simulation_results = expected_ground_state_simulation<sidb_bestagon_library_wrapper>(gates, params);

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

        const auto simulation_results = expected_ground_state_simulation<sidb_bestagon_library_wrapper>(gates, params);

        REQUIRE(simulation_results.charge_distributions.size() == 1);
        CHECK(simulation_results.charge_distributions.front().is_physically_valid());
        CHECK(simulation_results.simulation_runtime.count() > 0);
    }
}

TEST_CASE("Many double wires (Bestagon)", "[expect_gs]")
{
    hex_gate_lyt gates{{1, 20, 1}, row_clocking<hex_gate_lyt>()};

    const auto x1 = gates.create_pi("x1", {0, 0});
    const auto x2 = gates.create_pi("x1", {1, 0});

    const auto a1 = gates.create_buf(x1, {0, 1, 0});
    const auto a2 = gates.create_buf(x2, {0, 1, 1});
    const auto b1 = gates.create_buf(a1, {0, 2});
    const auto b2 = gates.create_buf(a2, {1, 2});
    const auto c1 = gates.create_buf(b1, {0, 3, 0});
    const auto c2 = gates.create_buf(b2, {0, 3, 1});
    const auto d1 = gates.create_buf(c1, {0, 4});
    const auto d2 = gates.create_buf(c2, {1, 4});
    const auto e1 = gates.create_buf(d1, {0, 5, 0});
    const auto e2 = gates.create_buf(d2, {0, 5, 1});
    const auto f1 = gates.create_buf(e1, {0, 6});
    const auto f2 = gates.create_buf(e2, {1, 6});
    const auto g1 = gates.create_buf(f1, {0, 7, 0});
    const auto g2 = gates.create_buf(f2, {0, 7, 1});
    const auto h1 = gates.create_buf(g1, {0, 8});
    const auto h2 = gates.create_buf(g2, {1, 8});
    const auto i1 = gates.create_buf(h1, {0, 9, 0});
    const auto i2 = gates.create_buf(h2, {0, 9, 1});
    const auto j1 = gates.create_buf(i1, {0, 10});
    const auto j2 = gates.create_buf(i2, {1, 10});
    const auto k1 = gates.create_buf(j1, {0, 11, 0});
    const auto k2 = gates.create_buf(j2, {0, 11, 1});
    const auto l1 = gates.create_buf(k1, {0, 12});
    const auto l2 = gates.create_buf(k2, {1, 12});
    const auto m1 = gates.create_buf(l1, {0, 13, 0});
    const auto m2 = gates.create_buf(l2, {0, 13, 1});
    const auto n1 = gates.create_buf(m1, {0, 14});
    const auto n2 = gates.create_buf(m2, {1, 14});
    const auto o1 = gates.create_buf(n1, {0, 15, 0});
    const auto o2 = gates.create_buf(n2, {0, 15, 1});
    const auto p1 = gates.create_buf(o1, {0, 16});
    const auto p2 = gates.create_buf(o2, {1, 16});
    const auto q1 = gates.create_buf(p1, {0, 17, 0});
    const auto q2 = gates.create_buf(p2, {0, 17, 1});
    const auto r1 = gates.create_buf(q1, {0, 18});
    const auto r2 = gates.create_buf(q2, {1, 18});
    const auto s1 = gates.create_buf(r1, {0, 19, 0});
    const auto s2 = gates.create_buf(r2, {0, 19, 1});

    gates.create_po(s1, "f1", {0, 20});
    gates.create_po(s2, "f2", {1, 20});

    const auto simulation_results = expected_ground_state_simulation<sidb_bestagon_library_wrapper>(gates, expected_gs_params{sidb_simulation_parameters{2, -0.32}});

    CHECK(simulation_results.charge_distributions.empty());
    CHECK(simulation_results.simulation_runtime.count() > 0);
}

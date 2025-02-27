//
// Created by Willem Lambooy on 18/02/2025.
//

#include <catch2/catch_test_macros.hpp>

#include <fiction/algorithms/simulation/sidb/skeleton_influence_bounds.hpp>
#include <fiction/layouts/gate_level_layout.hpp>
#include <fiction/technology/sidb_on_the_fly_mini_gate_library.hpp>
#include <fiction/technology/sidb_skeleton_bestagon_library.hpp>
#include <fiction/types.hpp>

#include <array>
#include <cstdint>
#include <unordered_map>

using namespace fiction;

using GateLyt = hex_even_row_gate_clk_lyt;
using CellLyt = sidb_cell_clk_lyt_cube;

TEST_CASE("Skeleton influence bounds of Bestagon gates in connection without clocking", "[skeleton-influence-bounds]")
{
    gate_level_layout<GateLyt> gate_lyt{{1, 1}};
    gate_lyt.create_and({}, {}, {0, 0});
    gate_lyt.create_and({}, {}, {1, 1});

    constexpr skeleton_influence_bounds_params<cell<CellLyt>> params{sidb_simulation_parameters{},
                                                                     {{24, 17}, {34, 28}}};

    const std::array<std::unordered_map<cell<CellLyt>, std::array<double, 2>>, 2>& res = {
        skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(gate_lyt, {0, 0}, params),
        skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(gate_lyt, {1, 1}, params)};

    SECTION("Influence of one gate")
    {
        for (const uint8_t i : std::array<uint8_t, 2>{{0, 1}})
        {
            CHECK(res.at(i).size() == (34 - 24 + 1) * (28 - 17 + 1) + 12);

            for (const auto& [c, bounds] : res.at(i))
            {
                CHECK(bounds[0] > 0);
                CHECK(bounds[0] < bounds[1]);
            }
        }
    }

    SECTION("Influence from two gates")
    {
        gate_level_layout<GateLyt> three_gate_lyt{{1, 2}};
        three_gate_lyt.create_and({}, {}, {0, 0});
        three_gate_lyt.create_and({}, {}, {1, 1});
        three_gate_lyt.create_and({}, {}, {1, 2});

        for (const uint8_t i : std::array<uint8_t, 2>{{0, 1}})
        {
            const auto& res_three_gates =
                skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(three_gate_lyt, {i, i}, params);

            CHECK(res_three_gates.size() == (34 - 24 + 1) * (28 - 17 + 1) + 12);

            for (const auto& [c, bounds] : res_three_gates)
            {
                CHECK(bounds[0] > 0);
                CHECK(bounds[0] < bounds[1]);

                REQUIRE(res.at(i).count(c) != 0);
                CHECK(bounds[0] > res.at(i).at(c)[0]);
                CHECK(bounds[1] > res.at(i).at(c)[1]);
            }
        }
    }
}

TEST_CASE("Skeleton influence bounds of Bestagon gates in connection with clocking", "[skeleton-influence-bounds]")
{
    GateLyt gate_lyt_straight_wire{{1, 2}, row_clocking<GateLyt>()};

    constexpr skeleton_influence_bounds_params<cell<CellLyt>> params{sidb_simulation_parameters{},
                                                                     {{24, 17}, {34, 28}}};

    constexpr std::array<tile<GateLyt>, 3> tile_order{{{0, 0}, {1, 1}, {0, 2}}};

    const auto s0 = gate_lyt_straight_wire.create_pi("A", tile_order.at(0));
    const auto s1 = gate_lyt_straight_wire.create_buf(s0, tile_order.at(1));
    gate_lyt_straight_wire.create_po(s1, "b1", tile_order.at(2));

    const std::array<std::unordered_map<cell<CellLyt>, std::array<double, 2>>, 3>& res{
        skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(gate_lyt_straight_wire, tile_order.at(0),
                                                                           params),
        skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(gate_lyt_straight_wire, tile_order.at(1),
                                                                           params),
        skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(gate_lyt_straight_wire, tile_order.at(2),
                                                                           params)};

    SECTION("Influence from two gates")
    {
        for (const auto& r : res)
        {
            for (const auto& [c, bounds] : r)
            {
                std::cout << c.x << "," << c.y << " " << bounds[0] << " " << bounds[1] << std::endl;
                CHECK(bounds[0] > 0);
                CHECK(bounds[0] < bounds[1]);
            }
        }
    }

    SECTION("Influence from three gates")
    {
        gate_level_layout<hex_even_row_gate_clk_lyt> gate_lyt_fo{{1, 2}, row_clocking<hex_even_row_gate_clk_lyt>()};
        const auto                                   s2 = gate_lyt_fo.create_pi("A", {0, 0});
        const auto                                   s3 = gate_lyt_fo.create_buf(s2, {1, 1});
        gate_lyt_fo.create_po(s3, "fo1", {0, 2});
        gate_lyt_fo.create_po(s3, "fo2", {1, 2});

        for (const uint8_t i : std::array<uint8_t, 2>{{0, 2}})
        {
            const auto& res_fo = skeleton_influence_bounds<CellLyt, sidb_skeleton_bestagon_library>(
                gate_lyt_fo, tile_order.at(i), params);

            for (const auto& [c, bounds] : res_fo)
            {
                CHECK(bounds[0] > 0);
                CHECK(bounds[0] < bounds[1]);

                REQUIRE(res.at(i).count(c) != 0);
                CHECK(bounds[0] > res.at(i).at(c)[0]);
                CHECK(bounds[1] > res.at(i).at(c)[1]);
            }
        }
    }
}

// TEST_CASE("Verify skeleton influence bounds in mini Bestagon circuit analytically", "[skeleton-influence-bounds]")
// {
//     GateLyt gate_lyt{{1, 2}, row_clocking<GateLyt>()};
//
//     const auto s  = gate_lyt.create_pi("A", {0, 0});
//     const auto s1 = gate_lyt.create_buf(s, {1, 1});
//     gate_lyt.create_po(s1, "fo1", {0, 2});
//     gate_lyt.create_po(s1, "fo2", {1, 2});
//
//     design_sidb_gates_params<sidb_defect_surface<CellLyt>> design_gate_params{};
//     design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{2, -0.32};
//     design_gate_params.operational_params.op_condition_positive_charges =
//         is_operational_params<cell<CellLyt>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
//     design_gate_params.operational_params.op_condition_kinks =
//         is_operational_params<cell<CellLyt>>::operational_condition_kinks::REJECT_KINKS;
//     design_gate_params.canvas                        = {{14, 10}, {21, 19}};  // smaller canvas
//     design_gate_params.number_of_sidbs               = 4;
//     design_gate_params.operational_params.sim_engine = sidb_simulation_engine::CLUSTERCOMPLETE;
//     design_gate_params.design_mode =
//         design_sidb_gates_params<sidb_defect_surface<CellLyt>>::design_sidb_gates_mode::QUICKCELL;
//
//     sidb_on_the_fly_gate_library_params<CellLyt> params{};
//     params.design_gate_params        = design_gate_params;
//     params.canvas_sidb_complex_gates = 4;
//
//     // for ()/
// }
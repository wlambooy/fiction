//
// Created by Willem Lambooy on 18/02/2025.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers_container_properties.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

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

        // verify local potentials for two skeleton SiDBs
        CHECK_THAT(res.at(0).at({16, 6}).at(0) - 0.000700401, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(0).at({16, 6}).at(1) - 0.000879205, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        CHECK_THAT(res.at(0).at({18, 8}).at(0) - 0.00091755, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(0).at({18, 8}).at(1) - 0.00115617, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        // verify local potentials for two canvas SiDBs
        CHECK_THAT(res.at(0).at({24, 17}).at(0) - 0.00252539, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(0).at({24, 17}).at(1) - 0.00325522, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        CHECK_THAT(res.at(0).at({34, 28}).at(0) - 0.0129194, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(0).at({34, 28}).at(1) - 0.017796, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        // verify local potentials for two skeleton SiDBs
        CHECK_THAT(res.at(1).at({16, 6}).at(0) - 0.0733532, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(1).at({16, 6}).at(1) - 0.0743375, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        CHECK_THAT(res.at(1).at({18, 8}).at(0) - 0.0722706, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(1).at({18, 8}).at(1) - 0.0729778, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        // verify local potentials for two canvas SiDBs
        CHECK_THAT(res.at(1).at({24, 17}).at(0) - 0.0113122, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(1).at({24, 17}).at(1) - 0.015815, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

        CHECK_THAT(res.at(1).at({34, 28}).at(0) - 0.00217837, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
        CHECK_THAT(res.at(1).at({34, 28}).at(1) - 0.00286426, Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
    }

    SECTION("Influence from two gates")
    {
        gate_level_layout<GateLyt> three_gate_lyt{{1, 2}};
        three_gate_lyt.create_and({}, {}, {0, 0});
        three_gate_lyt.create_and({}, {}, {1, 1});
        three_gate_lyt.create_and({}, {}, {1, 2});

        for (const uint8_t i : std::array<uint8_t, 3>{{0, 1}})
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

            // do extra checks for the sandwiched gate
            if (i == 1)
            {
                // verify local potentials for two skeleton SiDBs
                CHECK_THAT(res_three_gates.at({16, 6}).at(0) - 0.0740536,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
                CHECK_THAT(res_three_gates.at({16, 6}).at(1) - 0.0752167,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

                CHECK_THAT(res_three_gates.at({18, 8}).at(0) - 0.0731881,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
                CHECK_THAT(res_three_gates.at({18, 8}).at(1) - 0.074134,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

                // verify local potentials for two canvas SiDBs
                CHECK_THAT(res_three_gates.at({24, 17}).at(0) - 0.0138376,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
                CHECK_THAT(res_three_gates.at({24, 17}).at(1) - 0.0190702,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));

                CHECK_THAT(res_three_gates.at({34, 28}).at(0) - 0.0150978,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
                CHECK_THAT(res_three_gates.at({34, 28}).at(1) - 0.0206602,
                           Catch::Matchers::WithinAbs(0.0, constants::ERROR_MARGIN));
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
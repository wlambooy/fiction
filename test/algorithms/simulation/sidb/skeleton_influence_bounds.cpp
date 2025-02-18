//
// Created by Willem Lambooy on 18/02/2025.
//

#include <catch2/catch_test_macros.hpp>

#include <fiction/algorithms/simulation/sidb/skeleton_influence_bounds.hpp>
#include <fiction/layouts/gate_level_layout.hpp>
#include <fiction/technology/sidb_skeleton_bestagon_library.hpp>
#include <fiction/types.hpp>

using namespace fiction;

TEST_CASE("Bestagon skeleton influence bounds", "[skeleton-influence-bounds]")
{
    gate_level_layout<hex_even_row_gate_clk_lyt> gate_lyt{{1, 1}};
    gate_lyt.create_and({}, {}, {0, 0});
    gate_lyt.create_and({}, {}, {1, 1});

    const auto& res = skeleton_influence_bounds<sidb_100_cell_clk_lyt, sidb_skeleton_bestagon_library>(
        gate_lyt, {0, 0},
        skeleton_influence_bounds_params<cell<sidb_100_cell_clk_lyt>>{sidb_simulation_parameters{},
                                                                      {{24, 17}, {34, 28}}});

    SECTION("Influence of one gate")
    {
        CHECK(res.size() == (34 - 24 + 1) * (28 - 17 + 1));

        for (const auto& [c, bounds] : res)
        {
            CHECK(bounds[0] >= 0);
            CHECK(bounds[0] <= bounds[1]);
        }
    }

    SECTION("Influence onto a sandwiched gate")
    {
        gate_level_layout<hex_even_row_gate_clk_lyt> gate_lyt_sandwich{{2, 2}};
        gate_lyt_sandwich.create_and({}, {}, {0, 0});
        gate_lyt_sandwich.create_and({}, {}, {1, 1});
        gate_lyt_sandwich.create_and({}, {}, {1, 2});

        const auto& res_sandwich = skeleton_influence_bounds<sidb_100_cell_clk_lyt, sidb_skeleton_bestagon_library>(
            gate_lyt_sandwich, {1, 1},
            skeleton_influence_bounds_params<cell<sidb_100_cell_clk_lyt>>{sidb_simulation_parameters{},
                                                                          {{24, 17}, {34, 28}}});

        CHECK(res_sandwich.size() == (34 - 24 + 1) * (28 - 17 + 1));

        for (const auto& [c, bounds] : res_sandwich)
        {
            CHECK(bounds[0] >= 0);
            CHECK(bounds[0] <= bounds[1]);

            REQUIRE(res.count(c) != 0);
            CHECK(bounds[0] > res.at(c)[0]);
            CHECK(bounds[1] > res.at(c)[1]);
        }
    }
}

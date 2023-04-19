//
// Created by Willem Lambooy on 06/04/2023.
//

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/fuzzsim.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>

using namespace fiction;

TEMPLATE_TEST_CASE("Empty layout FuzzSim simulation", "[fuzzsim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    fuzzsim_stats<TestType> fuzzsimstats{};
    const fuzzsim_params fuzzsim_params{quicksim_params{sidb_simulation_parameters{2, -0.30}}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.30);

    fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

    CHECK(fuzzsimstats.valid_lyts.empty());
}

TEMPLATE_TEST_CASE("Single SiDB FuzzSim simulation", "[fuzzsim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);

    fuzzsim_stats<TestType> fuzzsimstats{};
    const fuzzsim_params fuzzsim_params{quicksim_params{sidb_simulation_parameters{2, -0.30}}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.30);

    fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

    CHECK(fuzzsimstats.valid_lyts.size() == 1);
}

template <typename Lyt>
void check_for_absence_of_positive_charges(const fuzzsim_stats<Lyt>& stats) noexcept
{
    REQUIRE(!stats.valid_lyts.empty());

    for (const auto& lyt : stats.valid_lyts)
    {
        CHECK(!lyt.charge_exists(sidb_charge_state::POSITIVE));
    }
}

template <typename Lyt>
void check_for_runtime_measurement(const fuzzsim_stats<Lyt>& stats) noexcept
{
    CHECK(stats.time_total.count() > 0);
}

TEMPLATE_TEST_CASE("FuzzSim simulation of several SiDBs with varying thread counts", "[fuzzsim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({3, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({4, 3, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({6, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({7, 3, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({6, 10, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({7, 10, 0}, TestType::cell_type::NORMAL);

    fuzzsim_stats<TestType>          fuzzsimstats{};
    const sidb_simulation_parameters params{2, -0.30};

    fuzzsim_params fuzzsim_params{quicksim_params{params}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.30);

    SECTION("Default settings")
    {
        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
    }
    SECTION("0 threads")
    {
        fuzzsim_params.qs_params.number_threads = 0;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
    }
    SECTION("1 thread")
    {
        fuzzsim_params.qs_params.number_threads = 1;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
    }
    SECTION("2 threads")
    {
        fuzzsim_params.qs_params.number_threads = 2;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
    }
    SECTION("100 threads")
    {
        fuzzsim_params.qs_params.number_threads = 100;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
    }
}


TEMPLATE_TEST_CASE("FuzzSim simulation of an SiDB layout comprising of 10 SiDBs with varying thread counts",
                   "[fuzzsim]", (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};
    
    lyt.assign_cell_type({-13, -1, 1}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({-9, -1, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({-7, -1, 1}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({-3, -1, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({-1, -1, 1}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({3, -1, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({5, -1, 1}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({9, -1, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({11, -1, 1}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({15, -1, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({17, -1, 1}, TestType::cell_type::NORMAL);

    fuzzsim_stats<TestType>          fuzzsimstats{};
    const sidb_simulation_parameters params{2, -0.32};

    fuzzsim_params fuzzsim_params{quicksim_params{params}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.32);

    const auto check_charge_configuration = [](const fuzzsim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin();

        CHECK(charge_lyt_first.get_charge_state({-13, -1, 1}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({-9, -1, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({-7, -1, 1}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({-3, -1, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({-1, -1, 1}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({3, -1, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({5, -1, 1}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({9, -1, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({11, -1, 1}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({15, -1, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({17, -1, 1}) == sidb_charge_state::NEGATIVE);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.47982940640, fiction::physical_constants::POP_STABILITY_ERR));
    };

    SECTION("Default settings")
    {
        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("0 threads")
    {
        fuzzsim_params.qs_params.number_threads = 0;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("1 thread")
    {
        fuzzsim_params.qs_params.number_threads = 1;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("2 threads")
    {
        fuzzsim_params.qs_params.number_threads = 2;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("100 threads")
    {
        fuzzsim_params.qs_params.number_threads = 100;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
}

TEMPLATE_TEST_CASE("FuzzSim simulation of a Y-shape SiDB arrangement with varying thread counts", "[fuzzsim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({-11, -2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({-10, -1, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({-4, -1, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({-3, -2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({-7, 0, 1}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({-7, 1, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({-7, 3, 0}, TestType::cell_type::NORMAL);

    fuzzsim_stats<TestType>          fuzzsimstats{};
    const sidb_simulation_parameters params{2, -0.32};

    fuzzsim_params fuzzsim_params{quicksim_params{params}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.32);

    const auto check_charge_configuration = [](const fuzzsim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin();

        CHECK(charge_lyt_first.get_charge_state({-11, -2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({-10, -1, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({-3, -2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({-4, -1, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({-7, 0, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({-7, 1, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({-7, 3, 0}) == sidb_charge_state::NEGATIVE);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.3191504062951, fiction::physical_constants::POP_STABILITY_ERR));
    };

    SECTION("Default settings")
    {
        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("0 threads")
    {
        fuzzsim_params.qs_params.number_threads = 0;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("1 thread")
    {
        fuzzsim_params.qs_params.number_threads = 1;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("2 threads")
    {
        fuzzsim_params.qs_params.number_threads = 2;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("100 threads")
    {
        fuzzsim_params.qs_params.number_threads = 100;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
}

TEMPLATE_TEST_CASE("FuzzSim simulation of a Y-shape SiDB OR gate with input 01 and varying thread counts",
                   "[fuzzsim]", (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({6, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({12, 3, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({14, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({10, 5, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({10, 6, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({10, 8, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({16, 1, 0}, TestType::cell_type::NORMAL);

    fuzzsim_stats<TestType>          fuzzsimstats{};
    const sidb_simulation_parameters params{2, -0.28};

    fuzzsim_params fuzzsim_params{quicksim_params{params}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.28);

    const auto check_charge_configuration = [](const fuzzsim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin();

        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.46621669, fiction::physical_constants::POP_STABILITY_ERR));
    };

    SECTION("Default settings")
    {
        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("0 threads")
    {
        fuzzsim_params.qs_params.number_threads = 0;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("1 thread")
    {
        fuzzsim_params.qs_params.number_threads = 1;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("2 threads")
    {
        fuzzsim_params.qs_params.number_threads = 2;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("100 threads")
    {
        fuzzsim_params.qs_params.number_threads = 100;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
}

TEMPLATE_TEST_CASE("FuzzSim simulation of an SiDB BDL pair with varying thread counts", "[fuzzsim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({6, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({8, 2, 0}, TestType::cell_type::NORMAL);

    fuzzsim_stats<TestType>          fuzzsimstats{};
    const sidb_simulation_parameters params{2, -0.25};

    fuzzsim_params fuzzsim_params{quicksim_params{params}};

    REQUIRE(fuzzsim_params.qs_params.phys_params.mu == -0.25);

    const auto check_charge_configuration = [](const fuzzsim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin();

        CHECK((((charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE) &&
                (charge_lyt_first.get_charge_state({8, 2, 0}) == sidb_charge_state::NEUTRAL)) ||
               ((charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEUTRAL) &&
                (charge_lyt_first.get_charge_state({8, 2, 0}) == sidb_charge_state::NEGATIVE))));
    };

    SECTION("Default settings")
    {
        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("0 threads")
    {
        fuzzsim_params.qs_params.number_threads = 0;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("1 thread")
    {
        fuzzsim_params.qs_params.number_threads = 1;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("2 threads")
    {
        fuzzsim_params.qs_params.number_threads = 2;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
    SECTION("100 threads")
    {
        fuzzsim_params.qs_params.number_threads = 100;

        fuzzsim<TestType>(lyt, fuzzsim_params, &fuzzsimstats);

        check_for_absence_of_positive_charges(fuzzsimstats);
        check_for_runtime_measurement(fuzzsimstats);
        check_charge_configuration(fuzzsimstats);
    }
}
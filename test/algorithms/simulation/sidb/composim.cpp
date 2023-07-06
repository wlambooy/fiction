//
// Created by Willem Lambooy on 06/04/2023.
//

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/composim.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>

using namespace fiction;

TEMPLATE_TEST_CASE("Empty layout CompoSim simulation", "[composim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    composim_stats<TestType> composimstats{};
    const composim_params    composim_params{quicksim_params{sidb_simulation_parameters{2, -0.30}}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.30);

    composim<TestType>(lyt, composim_params, &composimstats);

    CHECK(composimstats.valid_lyts.empty());
}

TEMPLATE_TEST_CASE("Single SiDB CompoSim simulation", "[composim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);

    composim_stats<TestType> composimstats{};
    const composim_params    composim_params{quicksim_params{sidb_simulation_parameters{2, -0.30}}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.30);

    const auto res = composim<TestType>(lyt, composim_params, &composimstats);

    CHECK(composim<TestType>(lyt, composim_params).charge_distributions.size() == 1);

    CHECK(composimstats.time_total.count() > 0);
    CHECK(res.simulation_runtime.count() > 0);
}

template <typename Lyt>
void check_for_absence_of_positive_charges(const composim_stats<Lyt>& stats) noexcept
{
    REQUIRE(!stats.valid_lyts.empty());

    for (const auto& lyt : stats.valid_lyts)
    {
        CHECK(!lyt.second->charge_exists(sidb_charge_state::POSITIVE));
    }
}

template <typename Lyt>
void check_for_runtime_measurement(const composim_stats<Lyt>& stats) noexcept
{
    CHECK(stats.time_total.count() > 0);
}

TEMPLATE_TEST_CASE("CompoSim simulation of several SiDBs with varying thread counts", "[composim]",
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

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.30};

    composim_params composim_params{quicksim_params{params}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.30);

    SECTION("Default settings")
    {
        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
}

TEMPLATE_TEST_CASE("CompoSim simulation of an SiDB layout comprising of 10 SiDBs with varying thread counts",
                   "[composim]", (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
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

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.32};

    composim_params composim_params{quicksim_params{params}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.32);

    const auto check_charge_configuration = [](const composim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin()->second;

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
        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
}

TEMPLATE_TEST_CASE("CompoSim simulation of a Y-shape SiDB arrangement with varying thread counts", "[composim]",
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

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.32};

    composim_params composim_params{quicksim_params{params}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.32);

    const auto check_charge_configuration = [](const composim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin()->second;

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
        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
}

TEMPLATE_TEST_CASE("CompoSim simulation of a Y-shape SiDB OR gate with input 01 and varying thread counts",
                   "[composim]", (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
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

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.28};

    composim_params composim_params{quicksim_params{params}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.28);

    const auto check_charge_configuration = [](const composim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin()->second;

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
        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
}

TEMPLATE_TEST_CASE("CompoSim simulation of an SiDB BDL pair with varying thread counts", "[composim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({6, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({8, 2, 0}, TestType::cell_type::NORMAL);

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.25};

    composim_params composim_params{quicksim_params{params}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.25);

    const auto check_charge_configuration = [](const composim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin()->second;

        CHECK((((charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE) &&
                (charge_lyt_first.get_charge_state({8, 2, 0}) == sidb_charge_state::NEUTRAL)) ||
               ((charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEUTRAL) &&
                (charge_lyt_first.get_charge_state({8, 2, 0}) == sidb_charge_state::NEGATIVE))));
    };

    SECTION("Default settings")
    {
        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
}

TEMPLATE_TEST_CASE("CompoSim simulation of an SiDB BDL pair with global potential with varying thread counts",
                   "[composim]", (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({6, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({8, 2, 0}, TestType::cell_type::NORMAL);

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.25};

    composim_params composim_params{quicksim_params{params, 80, 0.7, -0.05}};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.25);

    const auto check_charge_configuration = [](const std::vector<charge_distribution_surface<TestType>>& lyts) noexcept
    {
        REQUIRE(!lyts.empty());

        const auto& charge_lyt_first = *lyts.cbegin();

        CHECK((((charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE) &&
                (charge_lyt_first.get_charge_state({8, 2, 0}) == sidb_charge_state::NEUTRAL)) ||
               ((charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEUTRAL) &&
                (charge_lyt_first.get_charge_state({8, 2, 0}) == sidb_charge_state::NEGATIVE))));
    };

    SECTION("Default settings")
    {
        auto sim_res = composim<TestType>(lyt, composim_params, &composimstats);

        check_charge_configuration(sim_res.charge_distributions);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
    }
}

TEMPLATE_TEST_CASE("CompoSim simulation of a SiDB CX gate with input 01 and varying thread counts", "[composim]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{36, 19}};

    lyt.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({6, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({24, 5, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({26, 4, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({30, 3, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({32, 2, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({36, 1, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({14, 15, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({12, 16, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({8, 17, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({6, 18, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({2, 19, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({24, 15, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({30, 17, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({32, 18, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({36, 19, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({17, 8, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({24, 9, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({15, 10, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({15, 11, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({20, 12, 1}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({14, 13, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({17, 13, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({22, 13, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({20, 13, 1}, TestType::cell_type::NORMAL);

    composim_stats<TestType>         composimstats{};
    const sidb_simulation_parameters params{2, -0.32};

    composim_params composim_params{quicksim_params{params, 100, 0.8, 0}, 20, 1};

    REQUIRE(composim_params.qs_params.phys_params.mu == -0.32);
    REQUIRE(composim_params.max_charge_layouts == 1);

    const auto check_charge_configuration = [](const composim_stats<TestType>& stats) noexcept
    {
        REQUIRE(!stats.valid_lyts.empty());
        CHECK(stats.valid_lyts.size() == 1);

        const auto& charge_lyt_first = *stats.valid_lyts.cbegin()->second;

        CHECK(charge_lyt_first.get_charge_state({0, 0, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({12, 4, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({14, 5, 0}) == sidb_charge_state::NEUTRAL);

        CHECK(charge_lyt_first.get_charge_state({24, 5, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({26, 4, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({30, 3, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({32, 2, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({36, 1, 0}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({14, 15, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({12, 16, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({8, 17, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({6, 18, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({2, 19, 0}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({24, 15, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({26, 16, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({30, 17, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({32, 18, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({36, 19, 0}) == sidb_charge_state::NEGATIVE);

        CHECK(charge_lyt_first.get_charge_state({17, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({24, 9, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({15, 10, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({15, 11, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({20, 12, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({14, 13, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({17, 13, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({22, 13, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({20, 13, 1}) == sidb_charge_state::NEUTRAL);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(1.2818425861, fiction::physical_constants::POP_STABILITY_ERR));
    };

    SECTION("Default settings")
    {
        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("0 threads")
    {
        composim_params.qs_params.number_threads = 0;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("1 thread")
    {
        composim_params.qs_params.number_threads = 1;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("2 threads")
    {
        composim_params.qs_params.number_threads = 2;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
    SECTION("100 threads")
    {
        composim_params.qs_params.number_threads = 100;

        composim<TestType>(lyt, composim_params, &composimstats);

        check_for_absence_of_positive_charges(composimstats);
        check_for_runtime_measurement(composimstats);
        check_charge_configuration(composimstats);
    }
}

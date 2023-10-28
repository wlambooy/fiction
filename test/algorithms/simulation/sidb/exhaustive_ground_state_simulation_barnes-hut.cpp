//
// Created by Jan Drewniok on 18.12.22.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp>
#include <fiction/algorithms/simulation/sidb/quicksim.hpp>
#include <fiction/layouts/cartesian_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/layouts/clocked_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/physical_constants.hpp>
#include <fiction/utils/layout_utils.hpp>
#include <fiction/algorithms/physical_design/apply_gate_library.hpp>
#include <fiction/technology/sidb_bestagon_library.hpp>
#include <fiction/io/read_sqd_layout.hpp>

#include <cstdint>
#include <algorithm>
#include <utility>
#include <iostream>


using namespace fiction;

TEMPLATE_TEST_CASE("Empty layout ExGS simulation", "[exhaustive-ground-state-simulation]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    const sidb_simulation_parameters params{2, -0.32};

    const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

    CHECK(simulation_results.charge_distributions.empty());
    CHECK(simulation_results.additional_simulation_parameters.empty());
    CHECK(simulation_results.algorithm_name == "ExGS");
    CHECK(simulation_results.additional_simulation_parameters.empty());
}

TEMPLATE_TEST_CASE("Single SiDB ExGS simulation", "[exhaustive-ground-state-simulation]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};
    lyt.assign_cell_type({1, 3, 0}, TestType::cell_type::NORMAL);

    const sidb_simulation_parameters params{2, -0.32};

    const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

    REQUIRE(simulation_results.charge_distributions.size() == 1);
    CHECK(simulation_results.charge_distributions.front().get_charge_state_by_index(0) == sidb_charge_state::NEGATIVE);
}

TEMPLATE_TEST_CASE("ExGS simulation of a one BDL pair with one perturber", "[exhaustive-ground-state-simulation]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({4, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({6, 0, 0}, TestType::cell_type::NORMAL);

    const sidb_simulation_parameters params{2, -0.32};

    const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);
    CHECK(simulation_results.charge_distributions.size() == 1);
}

TEMPLATE_TEST_CASE("ExGS simulation of a two-pair BDL wire with one perturber", "[exhaustive-ground-state-simulation]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
{
    TestType lyt{{20, 10}};

    lyt.assign_cell_type({0, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({5, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({7, 0, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({11, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({13, 0, 0}, TestType::cell_type::NORMAL);

    lyt.assign_cell_type({17, 0, 0}, TestType::cell_type::NORMAL);
    lyt.assign_cell_type({19, 0, 0}, TestType::cell_type::NORMAL);

    const sidb_simulation_parameters params{2, -0.32};

    const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

    const auto size_before = simulation_results.charge_distributions.size();

    const auto simulation_results_after = exhaustive_ground_state_simulation<TestType>(lyt, params);
    auto       size_after               = simulation_results_after.charge_distributions.size();

    CHECK(size_before == 1);
    CHECK(size_after == 1);

    REQUIRE(!simulation_results_after.charge_distributions.empty());

    const auto& charge_lyt_first = simulation_results_after.charge_distributions.front();

    CHECK(charge_lyt_first.get_charge_state({0, 0, 0}) == sidb_charge_state::NEGATIVE);
    CHECK(charge_lyt_first.get_charge_state({5, 0, 0}) == sidb_charge_state::NEUTRAL);
    CHECK(charge_lyt_first.get_charge_state({7, 0, 0}) == sidb_charge_state::NEGATIVE);
    CHECK(charge_lyt_first.get_charge_state({11, 0, 0}) == sidb_charge_state::NEUTRAL);
    CHECK(charge_lyt_first.get_charge_state({13, 0, 0}) == sidb_charge_state::NEGATIVE);
    CHECK(charge_lyt_first.get_charge_state({17, 0, 0}) == sidb_charge_state::NEUTRAL);
    CHECK(charge_lyt_first.get_charge_state({19, 0, 0}) == sidb_charge_state::NEGATIVE);

    CHECK_THAT(charge_lyt_first.get_system_energy(),
               Catch::Matchers::WithinAbs(0.2460493219, physical_constants::POP_STABILITY_ERR));
}

TEMPLATE_TEST_CASE("ExGS simulation of a Y-shape SiDB arrangement", "[exhaustive-ground-state-simulation]",
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

    const sidb_simulation_parameters params{2, -0.32};

    const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

    REQUIRE(!simulation_results.charge_distributions.empty());

    const auto& charge_lyt_first = simulation_results.charge_distributions.front();

    CHECK(charge_lyt_first.get_charge_state({-11, -2, 0}) == sidb_charge_state::NEGATIVE);
    CHECK(charge_lyt_first.get_charge_state({-10, -1, 0}) == sidb_charge_state::NEUTRAL);
    CHECK(charge_lyt_first.get_charge_state({-3, -2, 0}) == sidb_charge_state::NEGATIVE);
    CHECK(charge_lyt_first.get_charge_state({-4, -1, 0}) == sidb_charge_state::NEUTRAL);
    CHECK(charge_lyt_first.get_charge_state({-7, 0, 1}) == sidb_charge_state::NEGATIVE);
    CHECK(charge_lyt_first.get_charge_state({-7, 1, 1}) == sidb_charge_state::NEUTRAL);
    CHECK(charge_lyt_first.get_charge_state({-7, 3, 0}) == sidb_charge_state::NEGATIVE);

    CHECK_THAT(charge_lyt_first.get_system_energy(),
               Catch::Matchers::WithinAbs(0.3191788254, physical_constants::POP_STABILITY_ERR));
}

TEMPLATE_TEST_CASE("ExGS simulation of a Y-shape SiDB OR gate with input 01", "[exhaustive-ground-state-simulation]",
                   (cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>))
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

    sidb_simulation_parameters params{2, -0.28};

    SECTION("Standard Physical Parameters")
    {
        const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

        REQUIRE(!simulation_results.charge_distributions.empty());
        const auto& charge_lyt_first = simulation_results.charge_distributions.front();

        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEUTRAL);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.4662582096, physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Increased mu_minus")
    {
        params.mu_minus = -0.1;

        const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

        REQUIRE(!simulation_results.charge_distributions.empty());
        const auto& charge_lyt_first = simulation_results.charge_distributions.front();

        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEUTRAL);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.061037632, physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Decreased mu_minus")
    {
        params.mu_minus = -0.7;

        const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

        REQUIRE(!simulation_results.charge_distributions.empty());
        const auto& charge_lyt_first = simulation_results.charge_distributions.front();

        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEGATIVE);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(2.069954113, physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Decreased lambda_tf")
    {
        params.lambda_tf = 1;

        const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

        REQUIRE(!simulation_results.charge_distributions.empty());
        const auto& charge_lyt_first = simulation_results.charge_distributions.front();

        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEGATIVE);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.5432404075, physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Increased lambda_tf")
    {
        params.lambda_tf = 10;

        const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

        REQUIRE(!simulation_results.charge_distributions.empty());
        const auto& charge_lyt_first = simulation_results.charge_distributions.front();

        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEUTRAL);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.2930574885, physical_constants::POP_STABILITY_ERR));
    }

    SECTION("Increased epsilon_r")
    {
        params.epsilon_r = 10;

        const auto simulation_results = exhaustive_ground_state_simulation<TestType>(lyt, params);

        REQUIRE(!simulation_results.charge_distributions.empty());
        const auto& charge_lyt_first = simulation_results.charge_distributions.front();

        CHECK(charge_lyt_first.get_charge_state({6, 2, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({12, 3, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 8, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 6, 1}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({16, 1, 0}) == sidb_charge_state::NEGATIVE);
        CHECK(charge_lyt_first.get_charge_state({10, 5, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({14, 2, 0}) == sidb_charge_state::NEUTRAL);
        CHECK(charge_lyt_first.get_charge_state({8, 3, 0}) == sidb_charge_state::NEGATIVE);

        CHECK_THAT(charge_lyt_first.get_system_energy(),
                   Catch::Matchers::WithinAbs(0.505173434, physical_constants::POP_STABILITY_ERR));
    }
}

template <typename Lyt, bool b>
std::optional<std::pair<double, uint64_t>> sort_by_energy(std::vector<charge_distribution_surface<Lyt, b>> cls)
{
    const auto comp = [](charge_distribution_surface<Lyt, b>& cl1, charge_distribution_surface<Lyt, b>& cl2) { return cl1.get_system_energy() < cl2.get_system_energy(); };

    std::sort(cls.begin(), cls.end(), comp);

    uint64_t max{0};
    for (const auto& cl : cls)
    {
        max = cl.get_num_validity_checks() > max ? cl.get_num_validity_checks() : max;
//        std::cout << "CHARGE INDEX: " << cl.get_charge_index_and_base().first << "\tENERGY: " << cl.get_system_energy() << std::endl;
    }

    std::cout << "\nEXACT VALIDITY CHECKS: " << max << "\n" << std::endl;

    return cls.empty() ? std::nullopt : std::make_optional(std::make_pair(cls.front().get_system_energy(), cls.front().get_charge_index_and_base().first));
}

TEST_CASE("ExGS simulation of a Bestagon wire (base 2)", "[ExGS]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{0, 0}};

    layout.create_pi("x1", {0, 0});

    sidb_simulation_parameters params{2, -0.32};

    charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> cl{
        convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)),
        params, sidb_charge_state::NEGATIVE};

    auto res = exhaustive_ground_state_simulation(cl, params);

    std::optional<std::pair<double, uint64_t>> gs_info = sort_by_energy(res.charge_distributions);

    REQUIRE(gs_info.has_value());

    CHECK_THAT(gs_info.value().first, Catch::Matchers::WithinAbs(0.5252063592, fiction::physical_constants::POP_STABILITY_ERR));
    CHECK(gs_info.value().second == 5322);
}

TEST_CASE("ExGS simulation of a Bestagon wire (base 3)", "[ExGS]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{0, 0}};

    layout.create_pi("x1", {0, 0});

    sidb_simulation_parameters params{3, -0.32};

    const charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> cl{
        convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)),
        params, sidb_charge_state::NEGATIVE};

    auto res = exhaustive_ground_state_simulation(cl, params);

    std::optional<std::pair<double, uint64_t>> gs_info = sort_by_energy(res.charge_distributions);

    REQUIRE(gs_info.has_value());

    CHECK_THAT(gs_info.value().first, Catch::Matchers::WithinAbs(0.5252063592, fiction::physical_constants::POP_STABILITY_ERR));
    CHECK(gs_info.value().second == 593436);
}

TEST_CASE("ExGS simulation of two Bestagon wires", "[ExGS]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{2, 0}};

    layout.create_pi("x1", {0, 0});
    layout.create_pi("x2", {2, 0});

    sidb_simulation_parameters params{2, -0.32};

    charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> cl{
        convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)),
        params, sidb_charge_state::NEGATIVE};

    auto res = exhaustive_ground_state_simulation(cl, params);

    std::optional<std::pair<double, uint64_t>> gs_info = sort_by_energy(res.charge_distributions);

    REQUIRE(gs_info.has_value());

    CHECK_THAT(gs_info.value().first, Catch::Matchers::WithinAbs(1.0504611338, fiction::physical_constants::POP_STABILITY_ERR));
    CHECK(gs_info.value().second == 53539020);
}

TEST_CASE("ExGS simulation of a Bestagon OR-gate (base 2)", "[ExGS]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{0, 0}};

    layout.create_or({}, {}, {0, 0});

    sidb_simulation_parameters params{2, -0.32};

    charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> cl{
        convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)),
        params, sidb_charge_state::NEGATIVE};

    auto res = exhaustive_ground_state_simulation(cl, params);

    std::optional<std::pair<double, uint64_t>> gs_info = sort_by_energy(res.charge_distributions);

    REQUIRE(gs_info.has_value());

    CHECK_THAT(gs_info.value().first, Catch::Matchers::WithinAbs(0.6021180782, fiction::physical_constants::POP_STABILITY_ERR));
    CHECK(gs_info.value().second == 52458);
}

TEST_CASE("ExGS simulation of a Bestagon OR-gate (base 3)", "[ExGS]")
{
    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;

    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{0, 0}};

    layout.create_or({}, {}, {0, 0});

    sidb_simulation_parameters params{3, -0.32};

    const charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> cl{
        convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)),
        params, sidb_charge_state::NEGATIVE};

    auto res = exhaustive_ground_state_simulation(cl, params);

    std::optional<std::pair<double, uint64_t>> gs_info = sort_by_energy(res.charge_distributions);

    REQUIRE(gs_info.has_value());

    CHECK_THAT(gs_info.value().first, Catch::Matchers::WithinAbs(0.6021180782, fiction::physical_constants::POP_STABILITY_ERR));
    CHECK(gs_info.value().second == 19371261);
}


TEST_CASE("testje", "[ExGS]")
{
    using sidb_layout = cell_level_layout<sidb_technology, clocked_layout<cartesian_layout<siqad::coord_t>>>;

    std::ifstream infile{"../../cmake-build-release/test/4-in-NAND.sqd", std::ifstream::in};
    REQUIRE(infile.is_open());

    sidb_layout layout = read_sqd_layout<sidb_layout>(infile);

    sidb_simulation_parameters params{2, -0.32};

    charge_distribution_surface<sidb_cell_clk_lyt_siqad> cl{layout, params};

//    auto res = quicksim(layout, quicksim_params{params, 10000});

//    using hex_gate_lyt = hex_odd_row_gate_clk_lyt;
//
//    hex_gate_lyt layout{aspect_ratio<hex_gate_lyt>{2, 0}};
//
//    layout.create_pi("x1", {0, 0});
//    layout.create_pi("x2", {2, 0});
//
//    sidb_simulation_parameters params{2, -0.32};
//
//    charge_distribution_surface<sidb_cell_clk_lyt_siqad, true> cl{
//        convert_to_siqad_coordinates(apply_gate_library<sidb_cell_clk_lyt, sidb_bestagon_library>(layout)),
//        params, sidb_charge_state::NEGATIVE};

    cl.assign_charge_index(17753454441793673608ul);
//    cl.assign_two_part_charge_index(15ul, 18410715276690587644ul);
    cl.update_after_charge_change();
    REQUIRE(cl.is_physically_valid());
//    cl.charge_distribution_to_index_general();
//    auto ix = cl.get_charge_index_and_base().first;
//    REQUIRE(ix == 295111876382333861884);
//    cl.assign_charge_index(ix);
//    cl.update_after_charge_change();
//    CHECK(cl.is_physically_valid());

}
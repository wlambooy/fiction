//
// Created by Jan Drewniok on 05.03.25.
//

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_result.hpp>
#include <fiction/types.hpp>

using namespace fiction;

TEST_CASE("Determine the groundstate from simulation results", "[sidb-simulation-result]")
{
    using lattice = sidb_cell_clk_lyt;

    SECTION("Three distinct charge distributions")
    {
        lattice lyt{};

        lyt.assign_cell_type({5, 4}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({5, 5}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({5, 6}, lattice::cell_type::NORMAL);

        charge_distribution_surface cds1{lyt};

        charge_distribution_surface cds2{lyt};
        cds2.assign_all_charge_states(sidb_charge_state::NEUTRAL);
        cds2.update_after_charge_change();

        charge_distribution_surface cds3{lyt};
        cds2.assign_all_charge_states(sidb_charge_state::POSITIVE);
        cds3.update_after_charge_change();

        CHECK_THAT(cds2.get_electrostatic_potential_energy(), Catch::Matchers::WithinAbs(0.0, 0.00001));
        CHECK(cds2.get_electrostatic_potential_energy() < cds3.get_electrostatic_potential_energy());
        CHECK(cds2.get_electrostatic_potential_energy() < cds1.get_electrostatic_potential_energy());

        cds1.assign_charge_index(0, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);
        cds2.assign_charge_index(1, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);
        cds3.assign_charge_index(2, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);

        sidb_simulation_result<lattice> results{};
        results.charge_distributions = {cds1, cds2, cds3};
        results.algorithm_name       = "test";

        const auto ground_state = results.groundstates();
        CHECK(ground_state.size() == 1);
    }
    SECTION("Several charge distributions with degeneracy")
    {
        lattice lyt{};

        lyt.assign_cell_type({5, 4}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({6, 4}, lattice::cell_type::NORMAL);

        charge_distribution_surface cds1{lyt};
        cds1.assign_charge_state({5, 4}, sidb_charge_state::NEUTRAL);
        cds1.assign_charge_state({6, 4}, sidb_charge_state::NEGATIVE);
        cds1.update_after_charge_change();

        charge_distribution_surface cds2{lyt};
        cds2.assign_charge_state({5, 4}, sidb_charge_state::NEGATIVE);
        cds2.assign_charge_state({6, 4}, sidb_charge_state::NEUTRAL);
        cds2.update_after_charge_change();

        charge_distribution_surface cds3{lyt};
        cds2.assign_all_charge_states(sidb_charge_state::POSITIVE);
        cds3.update_after_charge_change();

        // copy cds2 to check for degeneracy.
        charge_distribution_surface cds4{cds2};

        cds1.assign_charge_index(0, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);
        cds2.assign_charge_index(1, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);
        cds3.assign_charge_index(2, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);
        cds4.assign_charge_index(3, charge_distribution_mode::KEEP_CHARGE_DISTRIBUTION);

        CHECK_THAT(cds2.get_electrostatic_potential_energy() - cds1.get_electrostatic_potential_energy(),
                   Catch::Matchers::WithinAbs(0.0, 0.00001));

        sidb_simulation_result<lattice> results{};
        results.charge_distributions = {cds1, cds2, cds3, cds4};
        results.algorithm_name       = "test";

        const auto ground_states = results.groundstates();
        REQUIRE(ground_states.size() == 2);
        CHECK_THAT(ground_states[0].get_electrostatic_potential_energy() -
                       ground_states[1].get_electrostatic_potential_energy(),
                   Catch::Matchers::WithinAbs(0.0, 0.00001));
    }
}

TEST_CASE("Determine the groundstate from simulation results for Si-111 lattice orientation",
          "[sidb-simulation-result]")
{
    using lattice = sidb_111_cell_clk_lyt;

    SECTION("Three charge distributions with a degenerated ground state")
    {
        lattice lyt{};

        lyt.assign_cell_type({0, 0}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({2, 0}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({4, 0}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({0, 3}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({2, 3}, lattice::cell_type::NORMAL);
        lyt.assign_cell_type({4, 3}, lattice::cell_type::NORMAL);

        const sidb_simulation_parameters params{2, -0.30};
        const auto                       results = exhaustive_ground_state_simulation(lyt, params);

        const auto ground_state = results.groundstates();
        REQUIRE(ground_state.size() == 2);
        CHECK_THAT(ground_state.front().get_electrostatic_potential_energy(),
                   Catch::Matchers::WithinAbs(0.29683, 0.00001));
        CHECK_THAT(ground_state.front().get_electrostatic_potential_energy() -
                       ground_state.back().get_electrostatic_potential_energy(),
                   Catch::Matchers::WithinAbs(0.0, 0.00001));
    }
}

TEMPLATE_TEST_CASE("Determine the groundstate of a two BDL pair wire with input 1 applied", "[sidb-simulation-result]",
                   sidb_100_cell_clk_lyt, sidb_cell_clk_lyt)
{
    TestType lyt{};

    lyt.assign_cell_type({2, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({6, 0, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 0, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 0, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({14, 0, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({18, 0, 0}, sidb_technology::cell_type::NORMAL);

    const sidb_simulation_parameters params{2, -0.32};

    const auto results = exhaustive_ground_state_simulation(lyt, params);

    const auto ground_state = results.groundstates();
    REQUIRE(ground_state.size() == 2);
}

TEMPLATE_TEST_CASE("get_ordered_weights_under_bounded_energy() with bounded local potentials", "[sidb-simulation-result]", sidb_100_cell_clk_lyt, sidb_cell_clk_lyt)
{
    using lyt = TestType;
    using cds_t = charge_distribution_surface<lyt, local_external_potential_type::BOUNDED>;
    using result_t = sidb_simulation_result<lyt, local_external_potential_type::BOUNDED>;

    lyt layout{};

    cell<lyt> c0 = {0, 0}, c1 = {2,0}, c2 = {4,0};

    layout.assign_cell_type(c0, lyt::cell_type::NORMAL);
    layout.assign_cell_type(c1, lyt::cell_type::NORMAL);
    layout.assign_cell_type(c2, lyt::cell_type::NORMAL);

    const std::unordered_map<cell<lyt>, std::array<double, 2>> potential_1{
        {c0, {-0.1, 0.1}}, {c1, {-0.1, 0.1}}, {c2, {-0.1, 0.1}}};
    const std::unordered_map<cell<lyt>, std::array<double, 2>> potential_2{
        {c0, {-0.2, 0.2}}, {c1, {-0.1, 0.15}}, {c2, {-0.1, 0.1}}};
    const std::unordered_map<cell<lyt>, std::array<double, 2>> potential_3{
        {c0, {-0.3, 0.3}}, {c1, {-0.2, 0.2}}, {c2, {-0.1, 0.1}}};

    // Charge distribution 1: lowest uncertainty
    cds_t cds1{layout};
    cds1.assign_local_external_potential(potential_1);
    cds1.assign_charge_state(c0, sidb_charge_state::NEGATIVE);
    cds1.assign_charge_state(c1, sidb_charge_state::NEUTRAL);
    cds1.assign_charge_state(c2, sidb_charge_state::NEUTRAL);
    cds1.update_after_charge_change();

    // Charge distribution 2: medium uncertainty
    cds_t cds2{layout};
    cds2.assign_local_external_potential(potential_2);
    cds2.assign_charge_state(c0, sidb_charge_state::NEGATIVE);
    cds2.assign_charge_state(c1, sidb_charge_state::NEGATIVE);
    cds2.assign_charge_state(c2, sidb_charge_state::NEUTRAL);
    cds2.update_after_charge_change();

    // Charge distribution 3: highest uncertainty
    cds_t cds3{layout};
    cds3.assign_local_external_potential(potential_3);
    cds3.assign_charge_state(c0, sidb_charge_state::NEGATIVE);
    cds3.assign_charge_state(c1, sidb_charge_state::NEGATIVE);
    cds3.assign_charge_state(c2, sidb_charge_state::NEGATIVE);
    cds3.update_after_charge_change();

    result_t result{};
    result.charge_distributions = {cds3, cds1, cds2};  // shuffled on purpose

    result.reduce_to_groundstates_under_bounded_energy();
    const auto weights = result.get_ordered_weights_under_bounded_energy();

    REQUIRE(weights.size() == result.charge_distributions.size());

    // Check weight range and normalization
    for (const auto& w : weights)
    {
        CHECK(w >= 0.0);
        CHECK(w <= 1.0);
    }

    const double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    CHECK_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-6));

    // Since cds1 has the smallest energy and least uncertainty, it should get the highest weight
    CHECK(weights.front() > weights[1]);
    CHECK(weights[1] > weights.back());
}

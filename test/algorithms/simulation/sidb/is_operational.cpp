//
// Created by Jan Drewniok on 11.09.23.
//

#include <catch2/catch_template_test_macros.hpp>

#include "utils/blueprints/layout_blueprints.hpp"

#include <fiction/algorithms/iter/bdl_input_iterator.hpp>
#include <fiction/algorithms/physical_design/apply_gate_library.hpp>
#include <fiction/algorithms/simulation/sidb/compare_by_ground_state_isolation.hpp>
#include <fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp>
#include <fiction/algorithms/simulation/sidb/is_operational.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_engine.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/io/print_layout.hpp>
#include <fiction/layouts/cell_level_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/sidb_defects.hpp>
#include <fiction/technology/sidb_on_the_fly_mini_gate_library.hpp>
#include <fiction/traits.hpp>
#include <fiction/types.hpp>
#include <fiction/utils/truth_table_utils.hpp>

#include <cstdint>
#include <optional>
#include <set>
#include <vector>

using namespace fiction;

TEST_CASE("SiQAD OR gate", "[is-operational]")
{
    const auto or_gate = blueprints::siqad_or_gate<sidb_cell_clk_lyt_siqad>();

    const sidb_100_cell_clk_lyt_siqad lat{or_gate};

    auto op_params = is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
        sidb_simulation_parameters{2, -0.32},
        sidb_simulation_engine::QUICKEXACT,
        bdl_input_iterator_params{detect_bdl_wires_params{1.5},
                                  bdl_input_iterator_params::input_bdl_configuration::PERTURBER_ABSENCE_ENCODED},
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::TOLERATE_KINKS,
        {},
        {},
        is_operational_params<
            cell<sidb_100_cell_clk_lyt_siqad>>::termination_condition::ALL_INPUT_COMBINATIONS_ASSESSED,
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::simulation_results_mode::KEEP_SIMULATION_RESULTS};

    SECTION(
        "determine if layout is operational, tolerate kinks, assess all input combinations and keep simulation results")
    {
        const auto assessment_results = is_operational(lat, std::vector<tt>{create_or_tt()}, op_params);
        CHECK(assessment_results.status == operational_status::OPERATIONAL);
        REQUIRE(assessment_results.assessment_per_input.has_value());
        REQUIRE(!assessment_results.assessment_per_input.value().empty());
        CHECK(assessment_results.assessment_per_input.value().front().simulation_results.has_value());
    }

    // from now on, we will discard simulation results
    op_params.simulation_results_retention =
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::simulation_results_mode::DISCARD_SIMULATION_RESULTS;

    SECTION("determine if layout is operational, tolerate kinks, assess all input combinations and discard simulation "
            "results")
    {
        const auto assessment_results = is_operational(lat, std::vector<tt>{create_or_tt()}, op_params);
        CHECK(assessment_results.status == operational_status::OPERATIONAL);
        REQUIRE(assessment_results.assessment_per_input.has_value());
        REQUIRE(assessment_results.assessment_per_input.value().size() == 4);
        for (uint64_t i = 0; i < 4; ++i)
        {
            CHECK(assessment_results.assessment_per_input.value().at(i).status == operational_status::OPERATIONAL);
            CHECK(!assessment_results.assessment_per_input.value().at(i).simulation_results.has_value());
        }

        const auto assessment_results_and = is_operational(lat, std::vector<tt>{create_and_tt()}, op_params);
        CHECK(assessment_results_and.status == operational_status::NON_OPERATIONAL);
        REQUIRE(assessment_results_and.assessment_per_input.has_value());
        REQUIRE(assessment_results_and.assessment_per_input.value().size() == 4);
        CHECK(assessment_results_and.assessment_per_input.value().at(0).status == operational_status::OPERATIONAL);
        CHECK(assessment_results_and.assessment_per_input.value().at(1).status == operational_status::NON_OPERATIONAL);
        CHECK(assessment_results_and.assessment_per_input.value().at(2).status == operational_status::NON_OPERATIONAL);
        CHECK(assessment_results_and.assessment_per_input.value().at(3).status == operational_status::OPERATIONAL);
    }

    SECTION("determine if layout is operational under non-realistic physical parameters, assess all input combinations")
    {
        const auto check_for_non_operationality = [&lat, &op_params]
        {
            const auto assessment_results = is_operational(lat, std::vector<tt>{create_or_tt()}, op_params);
            CHECK(assessment_results.status == operational_status::NON_OPERATIONAL);
            REQUIRE(assessment_results.assessment_per_input.has_value());
            REQUIRE(assessment_results.assessment_per_input.value().size() == 4);
            for (uint64_t i = 0; i < 4; ++i)
            {
                CHECK(assessment_results.assessment_per_input.value().at(i).status ==
                      operational_status::NON_OPERATIONAL);
            }
        };

        op_params.simulation_parameters.epsilon_r = 1.0e-3;
        check_for_non_operationality();

        op_params.op_condition_positive_charges = is_operational_params<
            cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
        check_for_non_operationality();

        op_params.op_condition_positive_charges = is_operational_params<
            cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_positive_charges::REJECT_POSITIVE_CHARGES;
        op_params.simulation_parameters.epsilon_r = 5.6;
    }

    // from now on, we will terminate when the first non-operational input combination is found
    op_params.termination_cond =
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::termination_condition::ON_FIRST_NON_OPERATIONAL;

    SECTION("determine if layout is operational, tolerate kinks and terminate on first non-operational assessment")
    {
        const auto assessment_results = is_operational(lat, std::vector<tt>{create_or_tt()}, op_params);
        CHECK(assessment_results.status == operational_status::OPERATIONAL);
        CHECK(!assessment_results.assessment_per_input.has_value());
    }

    // from now on, we will reject kinks
    op_params.op_condition_kinks =
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS;

    SECTION("determine if layout is operational, reject kinks")
    {
        CHECK(is_operational(lat, std::vector<tt>{create_and_tt()}, op_params).status ==
              operational_status::NON_OPERATIONAL);
    }

    SECTION("determine if kinks induce layout to become non-operational")
    {
        const auto kink_induced_non_operational =
            is_kink_induced_non_operational(lat, std::vector<tt>{create_or_tt()}, op_params);
        CHECK(kink_induced_non_operational);
    }

    const auto input_wires  = detect_bdl_wires(lat, detect_bdl_wires_params{1.5}, bdl_wire_selection::INPUT);
    const auto output_wires = detect_bdl_wires(lat, detect_bdl_wires_params{1.5}, bdl_wire_selection::OUTPUT);

    REQUIRE(input_wires.size() == 2);

    CHECK(input_wires[0].pairs.size() == 2);
    CHECK(input_wires[1].pairs.size() == 2);

    CHECK(output_wires.size() == 1);

    SECTION("use pre-determined I/O pins")
    {
        CHECK(is_operational(lat, std::vector<tt>{create_and_tt()}, op_params, input_wires, output_wires).status ==
              operational_status::NON_OPERATIONAL);
    }

    SECTION("determine if kinks induce layout to become non-operational")
    {
        CHECK(is_kink_induced_non_operational(lat, std::vector<tt>{create_or_tt()}, op_params, input_wires,
                                              output_wires));
    }

    SECTION("determine input patterns for which kinks induce layout to become non-operational")
    {
        const auto kink_induced_non_operational_input_pattern =
            kink_induced_non_operational_input_patterns(lat, std::vector<tt>{create_or_tt()}, op_params);

        CHECK(kink_induced_non_operational_input_pattern.size() == 1);

        op_params.op_condition_kinks =
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::TOLERATE_KINKS;
        CHECK(is_operational(lat, std::vector<tt>{create_or_tt()}, op_params).status ==
              operational_status::OPERATIONAL);
    }
}

TEST_CASE("SiQAD NAND gate", "[is-operational]")
{
    const auto nand_gate = blueprints::siqad_nand_gate<sidb_cell_clk_lyt_siqad>();

    const sidb_100_cell_clk_lyt_siqad lat{nand_gate};

    auto op_params = is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
        sidb_simulation_parameters{2, -0.28},
        sidb_simulation_engine::QUICKEXACT,
        bdl_input_iterator_params{detect_bdl_wires_params{1.5},
                                  bdl_input_iterator_params::input_bdl_configuration::PERTURBER_ABSENCE_ENCODED},
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS,
        {},
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_analysis_strategy::FILTER_THEN_SIMULATION,
        {},
        {}};

    SECTION("Pruning and simulation")
    {
        CHECK(is_operational(lat, std::vector<tt>{create_nand_tt()}, op_params).status ==
              operational_status::OPERATIONAL);
    }
    SECTION("only pruning")
    {
        op_params.strategy_to_analyze_operational_status =
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_analysis_strategy::FILTER_ONLY;
        CHECK(is_operational(lat, std::vector<tt>{create_nand_tt()}, op_params).status ==
              operational_status::OPERATIONAL);
    }

    const auto input_wires  = detect_bdl_wires(lat, detect_bdl_wires_params{2.0}, bdl_wire_selection::INPUT);
    const auto output_wires = detect_bdl_wires(lat, detect_bdl_wires_params{2.0}, bdl_wire_selection::OUTPUT);

    sidb_100_cell_clk_lyt_siqad canvas_lyt{};
    canvas_lyt.assign_cell_type({10, 4, 1}, sidb_technology::cell_type::NORMAL);
    canvas_lyt.assign_cell_type({10, 5, 1}, sidb_technology::cell_type::NORMAL);

    SECTION("use pre-determined I/O pins")
    {
        CHECK(is_operational(lat, std::vector<tt>{create_nand_tt()}, op_params, input_wires, output_wires,
                             std::optional{canvas_lyt})
                  .status == operational_status::OPERATIONAL);
    }
}

TEST_CASE("SiQAD's AND gate with input BDL pairs of different size", "[is-operational]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({0, 0, 1}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 1}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({20, 0, 1}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({19, 1, 1}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({4, 2, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({6, 3, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({14, 3, 1}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({16, 2, 1}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({10, 6, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({10, 7, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({10, 9, 1}, sidb_technology::cell_type::NORMAL);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    CHECK(is_operational(lat, std::vector<tt>{create_and_tt()},
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.28}})
              .status == operational_status::OPERATIONAL);
    CHECK(is_operational(lat, std::vector<tt>{create_and_tt()},
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.1}})
              .status == operational_status::NON_OPERATIONAL);
}

TEST_CASE("Bestagon FO2 gate", "[is-operational]")
{
    const auto lyt = blueprints::bestagon_fo2<sidb_cell_clk_lyt_siqad>();

    SECTION("using QuickExact")
    {
        CHECK(is_operational(lyt, std::vector<tt>{create_fan_out_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lyt, std::vector<tt>{create_fan_out_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::NON_OPERATIONAL);
    }

    SECTION("using QuickSim")
    {
        CHECK(is_operational(lyt, std::vector<tt>{create_fan_out_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKSIM})
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lyt, std::vector<tt>{create_fan_out_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKSIM})
                  .status == operational_status::NON_OPERATIONAL);
    }

#if (FICTION_ALGLIB_ENABLED)

    SECTION("using ClusterComplete")
    {
        CHECK(is_operational(lyt, std::vector<tt>{create_fan_out_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{3, -0.32}, sidb_simulation_engine::CLUSTERCOMPLETE})
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lyt, std::vector<tt>{create_fan_out_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{3, -0.30}, sidb_simulation_engine::CLUSTERCOMPLETE})
                  .status == operational_status::NON_OPERATIONAL);
    }

#endif  // FICTION_ALGLIB_ENABLED
}

TEST_CASE("Bestagon CROSSING gate", "[is-operational]")
{
    const auto lyt = blueprints::bestagon_crossing<sidb_cell_clk_lyt_siqad>();

    CHECK(lyt.num_cells() == 29);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    CHECK(is_operational(lat, create_crossing_wire_tt(),
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.32},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::OPERATIONAL);
    CHECK(is_operational(lat, create_crossing_wire_tt(),
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.30},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::NON_OPERATIONAL);
}

TEST_CASE("Bestagon AND gate", "[is-operational]")
{
    auto lyt = blueprints::bestagon_and<sidb_defect_cell_clk_lyt_siqad>();

    const sidb_simulation_parameters params{2, -0.32};

    SECTION("Without defects")
    {
        CHECK(lyt.num_cells() == 23);

        CHECK(is_operational(lyt, std::vector<tt>{create_and_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lyt, std::vector<tt>{create_and_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::NON_OPERATIONAL);
    }
    SECTION("With defects")
    {
        lyt.assign_sidb_defect({3, 16, 1},
                               sidb_defect{sidb_defect_type::UNKNOWN, -1, params.epsilon_r, params.lambda_tf});
        CHECK(is_operational(
                  lyt, std::vector<tt>{create_and_tt()},
                  is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{params, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::OPERATIONAL);

        // move defect one to the right
        lyt.move_sidb_defect({3, 16, 1}, {4, 16, 1});
        CHECK(is_operational(
                  lyt, std::vector<tt>{create_and_tt()},
                  is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{params, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::OPERATIONAL);

        // move defect one to the right
        lyt.move_sidb_defect({4, 16, 1}, {5, 16, 1});
        CHECK(is_operational(
                  lyt, std::vector<tt>{create_and_tt()},
                  is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{params, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::NON_OPERATIONAL);
    }
    SECTION("Check operation for different values of mu")
    {
        CHECK(is_operational(lyt, std::vector<tt>{create_and_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lyt, std::vector<tt>{create_and_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::NON_OPERATIONAL);
    }
    SECTION("Count the number of non-operational input combinations, accepting kinks")
    {
        const auto op_inputs =
            operational_input_patterns(lyt, std::vector<tt>{create_and_tt()},
                                       is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                           sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT});
        CHECK(op_inputs.size() == 1);
        CHECK(op_inputs == std::set<uint64_t>{3});
    }
}

TEST_CASE("Not working diagonal Wire", "[is-operational]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{};

    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 2, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 5, 0}, sidb_technology::cell_type::NORMAL);

    // canvas SiDB
    lyt.assign_cell_type({14, 6, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({24, 15, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({26, 16, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({30, 17, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({32, 18, 0}, sidb_technology::cell_type::OUTPUT);

    lyt.assign_cell_type({36, 19, 0}, sidb_technology::cell_type::NORMAL);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    CHECK(is_operational(lat, std::vector<tt>{create_id_tt()},
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.32},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::NON_OPERATIONAL);
}

TEMPLATE_TEST_CASE("AND gate on the H-Si(111)-1x1 surface", "[is-operational]", sidb_111_cell_clk_lyt_siqad,
                   cds_sidb_111_cell_clk_lyt_siqad)
{
    const auto lyt = blueprints::and_gate_111<TestType>();

    SECTION("check operation for different values of mu")
    {
        const auto op_inputs =
            operational_input_patterns(lyt, std::vector<tt>{create_and_tt()},
                                       is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                           sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT});
        CHECK(op_inputs.size() == 4);
        CHECK(op_inputs == std::set<uint64_t>{0, 1, 2, 3});
    }
    SECTION("count the number of non-operational input combinations")
    {
        const auto op_inputs =
            operational_input_patterns(lyt, std::vector<tt>{create_and_tt()},
                                       is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                           sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT});
        CHECK(op_inputs.size() == 2);
        CHECK(op_inputs == std::set<uint64_t>{0, 3});
    }

    SECTION("verify the operational status of the AND gate, which is mirrored on the x-axis. Note that the input BDL "
            "pairs are located at the bottom, while the output BDL pairs are at the top.")
    {
        const auto lyt_mirrored_x = blueprints::and_gate_111_mirrored_on_the_x_axis<TestType>();
        const auto op_inputs =
            operational_input_patterns(lyt_mirrored_x, std::vector<tt>{create_and_tt()},
                                       is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                           sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT});
        CHECK(op_inputs.size() == 4);
        CHECK(op_inputs == std::set<uint64_t>{0, 1, 2, 3});
    }
}

TEST_CASE(
    "AND gate with Bestagon structure and kink state on right input wire for input 01 and left input wire for input 10",
    "[is-operational]")
{
    const auto lyt = blueprints::and_gate_with_kink_states<sidb_cell_clk_lyt_siqad>();

    SECTION("allow kink states")
    {
        CHECK(is_operational(
                  lyt, std::vector<tt>{create_and_tt()},
                  is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.32}})
                  .status == operational_status::OPERATIONAL);
    }
    SECTION("reject kink states")
    {
        CHECK(is_operational(lyt, std::vector<tt>{create_and_tt()},
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT,
                                 bdl_input_iterator_params{},
                                 is_operational_params<
                                     cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS})
                  .status == operational_status::NON_OPERATIONAL);
    }
    SECTION("check if is_kink_induced_non_operational returns true")
    {
        // check if the function works correctly even if the parameter is wrong (kinks are accepted).
        CHECK(is_kink_induced_non_operational(
            lyt, std::vector<tt>{create_and_tt()},
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT, bdl_input_iterator_params{},
                is_operational_params<
                    cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::TOLERATE_KINKS}));
    }

    SECTION("check input patterns for which kinks induce the layout to become non-operational")
    {
        CHECK(kink_induced_non_operational_input_patterns(
                  lyt, std::vector<tt>{create_and_tt()},
                  is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                      sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT,
                      bdl_input_iterator_params{},
                      is_operational_params<
                          cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::TOLERATE_KINKS}) ==
              std::set<uint64_t>{1, 2});
    }
}

TEST_CASE("BDL wire", "[is-operational]")
{
    using layout = sidb_cell_clk_lyt_siqad;

    layout lyt{{24, 0}, "BDL wire"};

    lyt.assign_cell_type({0, 0, 0}, sidb_technology::cell_type::INPUT);
    lyt.assign_cell_type({3, 0, 0}, sidb_technology::cell_type::INPUT);

    lyt.assign_cell_type({6, 0, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({8, 0, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({12, 0, 0}, sidb_technology::cell_type::NORMAL);
    lyt.assign_cell_type({14, 0, 0}, sidb_technology::cell_type::NORMAL);

    lyt.assign_cell_type({18, 0, 0}, sidb_technology::cell_type::OUTPUT);
    lyt.assign_cell_type({20, 0, 0}, sidb_technology::cell_type::OUTPUT);

    // output perturber
    lyt.assign_cell_type({24, 0, 0}, sidb_technology::cell_type::NORMAL);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    sidb_simulation_parameters sim_params{};

    sim_params.base = 2;

    const is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>> params{sim_params};

    CHECK(is_operational(lyt, std::vector<tt>{create_id_tt()}, params).status == operational_status::OPERATIONAL);
}

TEST_CASE("Special wire that cannot be pruned, but is non-operational when kinks are rejected", "[is-operational]")
{
    sidb_cell_clk_lyt_siqad lyt{};

    // input wires
    lyt.assign_cell_type({0, 0, 0}, sidb_cell_clk_lyt_siqad::cell_type::INPUT);
    lyt.assign_cell_type({2, 1, 0}, sidb_cell_clk_lyt_siqad::cell_type::INPUT);

    lyt.assign_cell_type({6, 2, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);
    lyt.assign_cell_type({8, 3, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);

    lyt.assign_cell_type({14, 5, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);
    lyt.assign_cell_type({12, 4, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);

    // canvas SiDBs
    lyt.assign_cell_type({11, 7, 0}, sidb_cell_clk_lyt_siqad::cell_type::LOGIC);
    lyt.assign_cell_type({13, 13, 0}, sidb_cell_clk_lyt_siqad::cell_type::LOGIC);

    // output wires
    lyt.assign_cell_type({14, 15, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);
    lyt.assign_cell_type({12, 16, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);

    lyt.assign_cell_type({8, 17, 0}, sidb_cell_clk_lyt_siqad::cell_type::OUTPUT);
    lyt.assign_cell_type({6, 18, 0}, sidb_cell_clk_lyt_siqad::cell_type::OUTPUT);

    lyt.assign_cell_type({2, 19, 0}, sidb_cell_clk_lyt_siqad::cell_type::NORMAL);

    sidb_simulation_parameters sim_params{};

    sim_params.base = 2;

    is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>> params{sim_params};

    SECTION("Rejecting Kinks")
    {
        params.op_condition_kinks =
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS;
        params.strategy_to_analyze_operational_status = is_operational_params<
            cell<sidb_100_cell_clk_lyt_siqad>>::operational_analysis_strategy::FILTER_THEN_SIMULATION;

        CHECK(is_operational(lyt, std::vector<tt>{create_id_tt()}, params).status ==
              operational_status::NON_OPERATIONAL);
    }

    SECTION("Only conducting pruning and tolerating kinks")
    {
        params.op_condition_kinks =
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::TOLERATE_KINKS;
        params.strategy_to_analyze_operational_status =
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_analysis_strategy::FILTER_ONLY;

        CHECK(is_operational(lyt, std::vector<tt>{create_id_tt()}, params).status ==
              operational_status::NON_OPERATIONAL);
    }
}

// to save runtime in the CI, this test is only run in RELEASE mode
#ifdef NDEBUG
TEST_CASE("flipped CX bestagon gate", "[is-operational]")
{
    const auto lyt = blueprints::crossing_bestagon_shape_input_down_output_up<sidb_cell_clk_lyt_siqad>();

    CHECK(is_operational(
              lyt, create_crossing_wire_tt(),
              is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                  sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT, bdl_input_iterator_params{},
                  is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS})
              .status == operational_status::OPERATIONAL);

    const auto kink_induced_non_operational_input_pattern = kink_induced_non_operational_input_patterns(
        lyt, create_crossing_wire_tt(),
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
            sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT, bdl_input_iterator_params{},
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS});

    CHECK(kink_induced_non_operational_input_pattern.empty());

    const auto kink_induced_non_operational = is_kink_induced_non_operational(
        lyt, create_crossing_wire_tt(),
        is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
            sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT, bdl_input_iterator_params{},
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_condition_kinks::REJECT_KINKS});

    CHECK(!kink_induced_non_operational);
}

TEST_CASE("is operational check for Bestagon CX gate", "[is-operational], [quality]")
{
    const auto lyt = blueprints::bestagon_crossing<sidb_cell_clk_lyt_siqad>();

    CHECK(lyt.num_cells() == 29);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    SECTION("without predetermined wires")
    {
        CHECK(is_operational(lat, create_crossing_wire_tt(),
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lat, create_crossing_wire_tt(),
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT})
                  .status == operational_status::NON_OPERATIONAL);
    }

    SECTION("using predetermined wires")
    {
        const auto input_bdl_wires  = detect_bdl_wires(lat, detect_bdl_wires_params{}, bdl_wire_selection::INPUT);
        const auto output_bdl_wires = detect_bdl_wires(lat, detect_bdl_wires_params{}, bdl_wire_selection::OUTPUT);

        CHECK(is_operational(lat, create_crossing_wire_tt(),
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.32}, sidb_simulation_engine::QUICKEXACT},
                             input_bdl_wires, output_bdl_wires)
                  .status == operational_status::OPERATIONAL);
        CHECK(is_operational(lat, create_crossing_wire_tt(),
                             is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{
                                 sidb_simulation_parameters{2, -0.30}, sidb_simulation_engine::QUICKEXACT},
                             input_bdl_wires, output_bdl_wires)
                  .status == operational_status::NON_OPERATIONAL);
        CHECK(!is_kink_induced_non_operational(
            lat, create_crossing_wire_tt(),
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.30},
                                                                     sidb_simulation_engine::QUICKEXACT},
            input_bdl_wires, output_bdl_wires));
    }

    SECTION("using predetermined wires and only applying pruning without simulation")
    {
        const auto input_bdl_wires  = detect_bdl_wires(lat, detect_bdl_wires_params{}, bdl_wire_selection::INPUT);
        const auto output_bdl_wires = detect_bdl_wires(lat, detect_bdl_wires_params{}, bdl_wire_selection::OUTPUT);

        auto op_params = is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.32}};
        op_params.strategy_to_analyze_operational_status =
            is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>::operational_analysis_strategy::FILTER_ONLY;

        CHECK(is_operational(lat, create_crossing_wire_tt(), op_params, input_bdl_wires, output_bdl_wires).status ==
              operational_status::OPERATIONAL);
    }
}

TEST_CASE("is operational check for Bestagon double wire", "[is-operational], [quality]")
{
    const auto lyt = blueprints::bestagon_double_wire<sidb_cell_clk_lyt_siqad>();

    CHECK(lyt.num_cells() == 30);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    CHECK(is_operational(lat, create_double_wire_tt(),
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.32},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::OPERATIONAL);
    CHECK(is_operational(lat, create_double_wire_tt(),
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.30},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::NON_OPERATIONAL);
}

TEST_CASE("is operational check for Bestagon half adder", "[is-operational], [quality]")
{
    const auto lyt = blueprints::bestagon_ha<sidb_cell_clk_lyt_siqad>();

    CHECK(lyt.num_cells() == 26);

    const sidb_100_cell_clk_lyt_siqad lat{lyt};

    CHECK(is_operational(lat, create_half_adder_tt(),
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.32},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::OPERATIONAL);
    CHECK(is_operational(lat, create_half_adder_tt(),
                         is_operational_params<cell<sidb_100_cell_clk_lyt_siqad>>{sidb_simulation_parameters{2, -0.25},
                                                                                  sidb_simulation_engine::QUICKEXACT})
              .status == operational_status::NON_OPERATIONAL);
}
#endif

TEMPLATE_TEST_CASE("is operational check of a designed multi-tile fan-out", "[is-operational], [quality]",
                   // sidb_on_the_fly_gate_library, sidb_on_the_fly_mini_gate_library)
                   sidb_on_the_fly_mini_gate_library)
{
    using GateLyt  = hex_even_row_gate_clk_lyt;
    using CellLyt  = sidb_cell_clk_lyt_cube;
    using gate_lib = TestType;

    GateLyt gate_lyt{{1, 2}, row_clocking<GateLyt>()};

    const auto s  = gate_lyt.create_pi("A", {0, 0});
    const auto s1 = gate_lyt.create_buf(s, {1, 1});
    gate_lyt.create_po(s1, "fo1", {0, 2});
    gate_lyt.create_po(s1, "fo2", {1, 2});

    design_sidb_gates_params<sidb_defect_surface<CellLyt>> design_gate_params{};
    design_gate_params.operational_params.simulation_parameters = sidb_simulation_parameters{2, -0.32};
    design_gate_params.operational_params.op_condition_positive_charges =
        is_operational_params<cell<CellLyt>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES;
    design_gate_params.operational_params.op_condition_kinks =
        is_operational_params<cell<CellLyt>>::operational_condition_kinks::REJECT_KINKS;
    if constexpr (std::is_same_v<gate_lib, sidb_on_the_fly_mini_gate_library>)
    {
        design_gate_params.canvas = {{14, 10}, {21, 19}};  // smaller canvas
    }
    else
    {
        design_gate_params.canvas = {{24, 17}, {34, 28}};  // normal canvas
    }
    design_gate_params.number_of_sidbs               = 4;
    design_gate_params.operational_params.sim_engine = sidb_simulation_engine::CLUSTERCOMPLETE;
    design_gate_params.design_mode =
        design_sidb_gates_params<sidb_defect_surface<CellLyt>>::design_sidb_gates_mode::QUICKCELL;

    sidb_on_the_fly_gate_library_params<CellLyt> params{};
    params.design_gate_params            = design_gate_params;
    params.canvas_sidb_complex_gates     = 4;
    params.use_skeleton_influence_bounds = true;

    auto make_layout = [&]
    {
        try
        {
            CellLyt cell_lyt =
                apply_parameterized_gate_library<CellLyt, gate_lib, GateLyt,
                                                 sidb_on_the_fly_gate_library_params<CellLyt>>(gate_lyt, params);

            // cell_lyt.foreach_coordinate(
            //     [&](auto cell)
            //     {
            //         if (cell_lyt.get_cell_type(cell) == sidb_100_cell_clk_lyt_siqad::cell_type::OUTPUT)
            //         {
            //             std::cout << cell.x << " " << cell.y << std::endl;
            //         }
            //     });
            //
            // std::cout << std::endl;
            // print_sidb_layout(std::cout, cell_lyt, true);

            if constexpr (std::is_same_v<gate_lib, sidb_on_the_fly_mini_gate_library>)
            {

                // input
                cell_lyt.assign_cell_type({22, 1}, CellLyt::technology::cell_type::INPUT);
                cell_lyt.assign_cell_type({24, 3}, CellLyt::technology::cell_type::INPUT);

                // left output
                cell_lyt.assign_cell_type({40, 61}, CellLyt::technology::cell_type::OUTPUT);
                cell_lyt.assign_cell_type({42, 63}, CellLyt::technology::cell_type::OUTPUT);
                // cell_lyt.assign_cell_type({52, 69}, CellLyt::technology::cell_type::NORMAL);  // perturber

                // right output
                cell_lyt.assign_cell_type({76, 61}, CellLyt::technology::cell_type::OUTPUT);
                cell_lyt.assign_cell_type({78, 63}, CellLyt::technology::cell_type::OUTPUT);
                // cell_lyt.assign_cell_type({88, 69}, CellLyt::technology::cell_type::NORMAL);  // perturber
            }
            else
            {
                // input
                cell_lyt.assign_cell_type({40, 2}, CellLyt::technology::cell_type::INPUT);
                cell_lyt.assign_cell_type({42, 4}, CellLyt::technology::cell_type::INPUT);

                // left output
                cell_lyt.assign_cell_type({70, 104}, CellLyt::technology::cell_type::OUTPUT);
                cell_lyt.assign_cell_type({72, 106}, CellLyt::technology::cell_type::OUTPUT);
                // cell_lyt.assign_cell_type({76, 108}, CellLyt::technology::cell_type::NORMAL);  // perturber

                // right output
                cell_lyt.assign_cell_type({130, 104}, CellLyt::technology::cell_type::OUTPUT);
                cell_lyt.assign_cell_type({132, 106}, CellLyt::technology::cell_type::OUTPUT);
                // cell_lyt.assign_cell_type({136, 108}, CellLyt::technology::cell_type::NORMAL);  // perturber
            }

            print_sidb_layout(std::cout, cell_lyt, true);

            return cell_lyt;
        }
        catch (const gate_design_exception<tt, GateLyt>& e)
        {
            std::cout << "Gate design exception at tile: " << e.which_tile() << std::endl;
            return CellLyt{};
        }
    };

    SECTION("With sorting of designed gates")
    {
        design_gate_params.post_design_process = {
            std::make_unique<compare_by_minimum_ground_state_isolation<sidb_defect_surface<CellLyt>>>(),
            std::make_unique<compare_by_average_ground_state_isolation<sidb_defect_surface<CellLyt>>>()};

        CellLyt cell_lyt = make_layout();

        REQUIRE(!cell_lyt.is_empty());

        CHECK(is_operational(cell_lyt, create_fan_out_tt(),
                             is_operational_params<cell<CellLyt>>{
                                 sidb_simulation_parameters{2, -0.32},
                                 sidb_simulation_engine::CLUSTERCOMPLETE,
                                 {},
                                 is_operational_params<cell<CellLyt>>::operational_condition_kinks::REJECT_KINKS,
                                 is_operational_params<
                                     cell<CellLyt>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES})
                  .status == operational_status::OPERATIONAL);
    }
    SECTION("Without sorting of designed gates")
    {
        design_gate_params.post_design_process.clear();

        CellLyt cell_lyt = make_layout();

        REQUIRE(!cell_lyt.is_empty());

        CHECK(is_operational(cell_lyt, create_fan_out_tt(),
                             is_operational_params<cell<CellLyt>>{
                                 sidb_simulation_parameters{2, -0.32},
                                 sidb_simulation_engine::CLUSTERCOMPLETE,
                                 {},
                                 is_operational_params<cell<CellLyt>>::operational_condition_kinks::REJECT_KINKS,
                                 is_operational_params<
                                     cell<CellLyt>>::operational_condition_positive_charges::TOLERATE_POSITIVE_CHARGES})
                  .status == operational_status::OPERATIONAL);
    }
}

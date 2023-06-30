//
// Created by marcel on 28.03.22.
//

#include <catch2/catch_test_macros.hpp>

#include <fiction/technology/sidb_defects.hpp>

using namespace fiction;

TEST_CASE("Charged and uncharged defect types", "[sidb-defects]")
{
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::NONE}));
    CHECK(is_charged_defect(sidb_defect{sidb_defect_type::DB}));
    CHECK(is_charged_defect(sidb_defect{sidb_defect_type::SI_VACANCY}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::SINGLE_DIHYDRIDE}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::DIHYDRIDE_PAIR}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::ONE_BY_ONE}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::THREE_BY_ONE}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::SILOXANE}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::RAISED_SI}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::MISSING_DIMER}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::ETCH_PIT}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::STEP_EDGE}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::GUNK}));
    CHECK(!is_charged_defect(sidb_defect{sidb_defect_type::UNKNOWN}));

    CHECK(!is_neutral_defect(sidb_defect{sidb_defect_type::NONE}));
    CHECK(!is_neutral_defect(sidb_defect{sidb_defect_type::DB}));
    CHECK(!is_neutral_defect(sidb_defect{sidb_defect_type::SI_VACANCY}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::SINGLE_DIHYDRIDE}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::DIHYDRIDE_PAIR}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::ONE_BY_ONE}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::THREE_BY_ONE}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::SILOXANE}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::RAISED_SI}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::MISSING_DIMER}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::ETCH_PIT}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::STEP_EDGE}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::GUNK}));
    CHECK(is_neutral_defect(sidb_defect{sidb_defect_type::UNKNOWN}));
}

TEST_CASE("Defect extent", "[sidb-defects]")
{
    static constexpr std::pair<uint16_t, uint16_t> no_spacing{0, 0};
    static constexpr std::pair<uint16_t, uint16_t> neutral_spacing{SIDB_NEUTRAL_DEFECT_HORIZONTAL_SPACING,
                                                                   SIDB_NEUTRAL_DEFECT_VERTICAL_SPACING};
    static constexpr std::pair<uint16_t, uint16_t> charged_spacing{SIDB_CHARGED_DEFECT_HORIZONTAL_SPACING,
                                                                   SIDB_CHARGED_DEFECT_VERTICAL_SPACING};

    CHECK(defect_extent(sidb_defect{sidb_defect_type::NONE}) == no_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::DB}) == charged_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::SI_VACANCY}) == charged_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::SINGLE_DIHYDRIDE}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::DIHYDRIDE_PAIR}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::ONE_BY_ONE}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::THREE_BY_ONE}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::SILOXANE}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::RAISED_SI}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::MISSING_DIMER}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::ETCH_PIT}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::STEP_EDGE}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::GUNK}) == neutral_spacing);
    CHECK(defect_extent(sidb_defect{sidb_defect_type::UNKNOWN}) == neutral_spacing);
}

TEST_CASE("Compare Defect", "[sidb-defects]")
{
    SECTION("Different types")
    {
        const sidb_defect defect_one{sidb_defect_type::GUNK};
        const sidb_defect defect_two{sidb_defect_type::UNKNOWN};
        CHECK(defect_one != defect_two);
    }

    SECTION("Different charge")
    {
        const sidb_defect defect_one{sidb_defect_type::UNKNOWN, -5};
        const sidb_defect defect_two{sidb_defect_type::UNKNOWN, -1};
        CHECK(defect_one != defect_two);
    }

    SECTION("Different epsilon_r")
    {
        const sidb_defect defect_one{sidb_defect_type::UNKNOWN, -1, 2};
        const sidb_defect defect_two{sidb_defect_type::UNKNOWN, -1, 5};
        CHECK(defect_one != defect_two);
    }

    SECTION("Different lambda_tf")
    {
        const sidb_defect defect_one{sidb_defect_type::UNKNOWN, -1, 2, 4};
        const sidb_defect defect_two{sidb_defect_type::UNKNOWN, -1, 2, 5};
        CHECK(defect_one != defect_two);
    }

    SECTION("Completely different")
    {
        const sidb_defect defect_one{sidb_defect_type::UNKNOWN, -1, 2, 4};
        const sidb_defect defect_two{sidb_defect_type::DB, 5, 5, 0.3};
        CHECK(defect_one != defect_two);
    }

    SECTION("Identical Defects")
    {
        const sidb_defect defect_one{sidb_defect_type::UNKNOWN, -1, 2, 4};
        const sidb_defect defect_two{sidb_defect_type::UNKNOWN, -1, 2, 4};
        CHECK(defect_one == defect_two);
    }
}

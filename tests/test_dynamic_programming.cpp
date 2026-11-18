#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "project/DynamicProgramming.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

TEST_CASE("Matrix chain product computes classic optimal cost") {
    std::vector<int> dimensions = {30, 35, 15, 5, 10, 20, 25};

    MatrixChainResult result = matrixChainOrder(dimensions);

    REQUIRE(result.minimumCost == 15125);
    REQUIRE(result.parenthesization == "((A1 x (A2 x A3)) x ((A4 x A5) x A6))");
}

TEST_CASE("Matrix chain product handles two matrices") {
    std::vector<int> dimensions = {10, 20, 30};

    MatrixChainResult result = matrixChainOrder(dimensions);

    REQUIRE(result.minimumCost == 6000);
    REQUIRE(result.parenthesization == "(A1 x A2)");
}

TEST_CASE("Matrix chain product rejects invalid dimensions") {
    REQUIRE_THROWS_AS(
        matrixChainOrder({}),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        matrixChainOrder({10}),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        matrixChainOrder({10, 0, 20}),
        std::invalid_argument);
}

TEST_CASE("Sequence alignment rewards identical sequences") {
    AlignmentResult result =
        globalSequenceAlignment("ACGT", "ACGT");

    REQUIRE(result.score == 4);
    REQUIRE(result.alignedFirst == "ACGT");
    REQUIRE(result.alignedSecond == "ACGT");
}

TEST_CASE("Sequence alignment supports gaps") {
    const auto defaults = globalSequenceAlignment("GATTACA", "GCATGCU");
    REQUIRE(defaults.score == -1);
    REQUIRE(defaults.alignedFirst == "GATTACA");
    REQUIRE(defaults.alignedSecond == "GCATGCU");

    const auto result = globalSequenceAlignment("GATTACA", "GCATGCU", 1, -1, -1);

    REQUIRE(result.score == 0);
    REQUIRE(result.alignedFirst.size() == result.alignedSecond.size());
    REQUIRE((result.alignedFirst.find('-') != std::string::npos ||
             result.alignedSecond.find('-') != std::string::npos));
}

TEST_CASE("Sequence alignment preserves original characters") {
    AlignmentResult result =
        globalSequenceAlignment("ACG", "AG");

    REQUIRE(result.alignedFirst.size() == result.alignedSecond.size());

    std::string firstWithoutGaps;
    std::string secondWithoutGaps;

    for (char c : result.alignedFirst) {
        if (c != '-') {
            firstWithoutGaps += c;
        }
    }

    for (char c : result.alignedSecond) {
        if (c != '-') {
            secondWithoutGaps += c;
        }
    }

    REQUIRE(firstWithoutGaps == "ACG");
    REQUIRE(secondWithoutGaps == "AG");
}

TEST_CASE("Sequence alignment supports custom scoring") {
    AlignmentResult result =
        globalSequenceAlignment("AA", "AA", 2, -1, -2);

    REQUIRE(result.score == 4);
    REQUIRE(result.alignedFirst == "AA");
    REQUIRE(result.alignedSecond == "AA");
}

TEST_CASE("Matrix chain stores splits and handles a single matrix") {
    const auto single = matrixChainOrder({10, 20});
    REQUIRE(single.minimumCost == 0);
    REQUIRE(single.parenthesization == "A1");

    const auto chain = matrixChainOrder({5, 10, 3, 12});
    REQUIRE(chain.minimumCost == 330);
    REQUIRE(chain.split[1][3] == 2);
    REQUIRE(buildParenthesization(chain.split, 1, 3) == "((A1 x A2) x A3)");
}

TEST_CASE("Matrix chain skips an overflowing split when another split fits") {
    const int large = std::numeric_limits<int>::max();
    const auto chain = matrixChainOrder({large, large, 1, large});
    REQUIRE(chain.minimumCost == 2LL * large * large);
    REQUIRE(chain.split[1][3] == 2);
    REQUIRE_THROWS_AS(matrixChainOrder({large, large, large}), std::overflow_error);
}

TEST_CASE("Sequence alignment handles either empty input") {
    const auto firstEmpty = globalSequenceAlignment("", "AC", 2, -3, -4);
    REQUIRE(firstEmpty.score == -8);
    REQUIRE(firstEmpty.alignedFirst == "--");
    REQUIRE(firstEmpty.alignedSecond == "AC");

    const auto secondEmpty = globalSequenceAlignment("AC", "");
    REQUIRE(secondEmpty.score == -4);
    REQUIRE(secondEmpty.alignedFirst == "AC");
    REQUIRE(secondEmpty.alignedSecond == "--");

    REQUIRE(globalSequenceAlignment("", "").score == 0);
}

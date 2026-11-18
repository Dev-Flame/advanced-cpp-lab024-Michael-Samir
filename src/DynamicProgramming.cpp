#include "project/DynamicProgramming.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

std::string buildParenthesization(
    const std::vector<std::vector<int>>& split,
    int i,
    int j) {

    if (i < 1 || j < i) {
        throw std::invalid_argument("Invalid matrix chain range");
    }
    if (i == j) {
        return "A" + std::to_string(i);
    }
    if (static_cast<std::size_t>(i) >= split.size() ||
        static_cast<std::size_t>(j) >= split[i].size()) {
        throw std::invalid_argument("Missing matrix chain split");
    }
    const int k = split[i][j];
    if (k < i || k >= j) {
        throw std::invalid_argument("Invalid matrix chain split");
    }
    return "(" + buildParenthesization(split, i, k) + " x " +
           buildParenthesization(split, k + 1, j) + ")";
}

MatrixChainResult matrixChainOrder(
    const std::vector<int>& dimensions) {

    if (dimensions.size() < 2 ||
        std::any_of(dimensions.begin(), dimensions.end(), [](int d) { return d <= 0; })) {
        throw std::invalid_argument("At least two positive dimensions are required");
    }

    const std::size_t n = dimensions.size() - 1;
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("Too many matrices for the split table");
    }
    const long long limit = std::numeric_limits<long long>::max();
    std::vector<std::vector<long long>> cost(n + 1, std::vector<long long>(n + 1, 0));
    MatrixChainResult result{0, std::vector<std::vector<int>>(n + 1,
        std::vector<int>(n + 1, 0)), ""};

    // Shorter subchains are ready before any split that needs them.
    for (std::size_t length = 2; length <= n; ++length) {
        for (std::size_t i = 1; i + length - 1 <= n; ++i) {
            const std::size_t j = i + length - 1;
            long long best = -1;  // No representable cost found yet.
            for (std::size_t k = i; k < j; ++k) {
                const long long left = cost[i][k];
                const long long right = cost[k + 1][j];
                const long long outer = 1LL * dimensions[i - 1] * dimensions[j];
                if (left < 0 || right < 0 || left > limit - right ||
                    outer > (limit - left - right) / dimensions[k]) {
                    continue;
                }
                const long long candidate = left + right + outer * dimensions[k];
                if (best < 0 || candidate < best) {
                    best = candidate;
                    result.split[i][j] = static_cast<int>(k);
                }
            }
            cost[i][j] = best;
        }
    }
    result.minimumCost = cost[1][n];
    if (result.minimumCost < 0) {
        throw std::overflow_error("Matrix chain cost exceeds long long");
    }
    result.parenthesization = buildParenthesization(result.split, 1, static_cast<int>(n));
    return result;
}

AlignmentResult globalSequenceAlignment(
    const std::string& first,
    const std::string& second,
    int matchScore,
    int mismatchScore,
    int gapScore) {

    const std::size_t m = first.size();
    const std::size_t n = second.size();
    std::vector<std::vector<long long>> scores(m + 1, std::vector<long long>(n + 1, 0));
    for (std::size_t i = 1; i <= m; ++i) {
        scores[i][0] = scores[i - 1][0] + gapScore;
    }
    for (std::size_t j = 1; j <= n; ++j) {
        scores[0][j] = scores[0][j - 1] + gapScore;
    }
    for (std::size_t i = 1; i <= m; ++i) {
        for (std::size_t j = 1; j <= n; ++j) {
            const int pairScore = first[i - 1] == second[j - 1] ? matchScore : mismatchScore;
            scores[i][j] = std::max({scores[i - 1][j - 1] + pairScore,
                                     scores[i - 1][j] + gapScore,
                                     scores[i][j - 1] + gapScore});
        }
    }

    const long long score = scores[m][n];
    if (score < std::numeric_limits<int>::min() || score > std::numeric_limits<int>::max()) {
        throw std::overflow_error("Alignment score exceeds int");
    }
    AlignmentResult result{static_cast<int>(score), "", ""};
    std::size_t i = m;
    std::size_t j = n;
    // Trace backward, preferring diagonal transitions when scores tie.
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 &&
            scores[i][j] == scores[i - 1][j - 1] +
                (first[i - 1] == second[j - 1] ? matchScore : mismatchScore)) {
            result.alignedFirst += first[--i];
            result.alignedSecond += second[--j];
        } else if (i > 0 && scores[i][j] == scores[i - 1][j] + gapScore) {
            result.alignedFirst += first[--i];
            result.alignedSecond += '-';
        } else {
            result.alignedFirst += '-';
            result.alignedSecond += second[--j];
        }
    }
    std::reverse(result.alignedFirst.begin(), result.alignedFirst.end());
    std::reverse(result.alignedSecond.begin(), result.alignedSecond.end());
    return result;
}

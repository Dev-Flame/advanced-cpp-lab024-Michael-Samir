#include "project/DynamicProgramming.hpp"

#include <iostream>

int main() {
    const auto chain = matrixChainOrder({30, 35, 15, 5, 10, 20, 25});
    std::cout << "Matrix order: " << chain.parenthesization
              << "\nMinimum cost: " << chain.minimumCost << '\n';

    const auto alignment = globalSequenceAlignment("ACGT", "AGT");
    std::cout << "\nAlignment score: " << alignment.score << '\n'
              << alignment.alignedFirst << '\n' << alignment.alignedSecond << '\n';
}

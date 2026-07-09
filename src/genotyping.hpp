#pragma once

#include "htslib/sam.h"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "conversion.hpp"
#include "panman.hpp"
#include <iostream>

namespace genotyping {

void stringSplit(const std::string& str, char delimiter, std::vector<std::string>& out);

struct mutationMatrices {
    std::vector<std::vector<double>> submat;
    std::unordered_map<int64_t, double> insmat;
    std::unordered_map<int64_t, double> delmat;
    double maxInsLogProb;
    double maxDelLogProb;
    bool filled;

    mutationMatrices()
        : submat(4, std::vector<double>(4, 0.0)), maxInsLogProb(100.0), maxDelLogProb(100.0), filled(false) {}
};

void fillMutationMatricesFromFile(mutationMatrices& mutMat, std::ifstream& inf);

void buildMutationMatricesHelper(mutationMatrices& mutMat,
                                 panmanUtils::Tree* tree,
                                 panmanUtils::Node* node,
                                 std::vector<int64_t>& parentBaseCounts,
                                 std::vector<int64_t>& totalBaseCounts,
                                 std::vector<std::vector<int64_t>>& subCount,
                                 std::unordered_map<int64_t, int64_t>& insCount,
                                 std::unordered_map<int64_t, int64_t>& delCount);

std::string applyMutationSpectrum(const std::string& line, const std::vector<std::vector<double>>& scaled_submat);

}  // namespace genotyping

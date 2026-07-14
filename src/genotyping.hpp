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

std::string applyMutationSpectrum(const std::string& line, const std::vector<std::vector<double>>& scaled_submat,
                                  int minDepth, double minQual);

// Consensus gate for haploid genotyping: the called allele (calledIdx into AD; 0 = ref)
// must be the majority (strictly more high-quality reads than all others combined) and
// have >= minDepth high-quality reads, else the position is a no-call. Returns true (no
// filtering) only when AD is missing/uninformative.
bool passesConsensusGate(int calledIdx, const std::vector<int>& ad, int minDepth);
// Same gate applied to a raw VCF sample field ("GT:PL:AD").
bool passesConsensusGate(const std::string& sampleField, int minDepth);

}  // namespace genotyping

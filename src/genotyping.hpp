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

// Ensure correct usage of std::getline
  void stringSplit(const std::string& str, char delimiter, std::vector<std::string>& out);

  // Ensure correct initialization of vectors and maps
  struct mutationMatrices {
    std::vector<std::vector<double>> submat;
    std::unordered_map<int64_t, double> insmat;
    std::unordered_map<int64_t, double> delmat;
    double maxInsLogProb;
    double maxDelLogProb;
    bool filled;

    mutationMatrices() : submat(4, std::vector<double>(4, 0.0)), maxInsLogProb(100.0), maxDelLogProb(100.0), filled(false) {}
  };

  /**
   * Fill mutation matrices from a file
   * 
   * @param mutMat Mutation matrices to fill
   * @param inf Input file stream to read from
   */
  void fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf);

  /**
   * Build mutation matrices helper function
   * 
   * @param mutMat Mutation matrices to fill
   * @param tree Tree to process
   * @param node Current node being processed
   * @param parentBaseCounts Parent base counts
   * @param totalBaseCounts Total base counts
   * @param subCount Substitution counts
   * @param insCount Insertion counts
   * @param delCount Deletion counts
   */
  void buildMutationMatricesHelper(
      mutationMatrices &mutMat,
      panmanUtils::Tree *tree,
      panmanUtils::Node *node,
      std::vector<int64_t> &parentBaseCounts,
      std::vector<int64_t> &totalBaseCounts,
      std::vector<std::vector<int64_t>> &subCount,
      std::unordered_map<int64_t, int64_t> &insCount,
      std::unordered_map<int64_t, int64_t> &delCount);

std::vector<std::vector<double>>
scaleMutationSpectrum(const mutationMatrices &mutMat, double mutationRate);
std::string
applyMutationSpectrum(const std::string &line,
                      const std::vector<std::vector<double>> &scaled_submat);

} // namespace genotyping
#pragma once

/**
 * @file seed_helpers.hpp
 * @brief Seed extraction and FASTA/read helpers shared by tests.
 */

#include "metrics_oracle.hpp"  // indexUtils::SeedCountMap

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ts {

// Extract syncmer seeds as the read/genome paths do: rollingSyncmers(seq, k, s, open=false, t=0),
// counting only flagged syncmers. countDuplicates=false collapses to presence (count 1).
indexUtils::SeedCountMap extractSeeds(std::string_view seq, int k, int s, bool countDuplicates = true);

// Read a FASTA file and return the concatenated sequence (ignores headers/whitespace).
std::string readFasta(const std::string& path);

// Deterministically sample fixed-length reads from a genome (skips reads with >25% N).
std::vector<std::string> generateReads(const std::string& genome, int readLen, int numReads, uint64_t seed);

}  // namespace ts

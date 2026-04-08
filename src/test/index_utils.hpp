#pragma once

/**
 * @file index_utils.hpp
 * @brief Shared utilities for metrics computation and node lookup in tests
 */

#include "panmap_utils.hpp"
#include <absl/container/flat_hash_map.h>
#include <cstdint>
#include <string>
#include <cmath>
#include <algorithm>

namespace indexUtils {

// ============================================================================
// Seed Set Types
// ============================================================================

using SeedCountMap = absl::flat_hash_map<uint64_t, int64_t>;

// ============================================================================
// Ground Truth Metrics (computed directly from seed sets)
// ============================================================================

/**
 * @brief Metrics computed directly from read and genome seed sets
 * 
 * This is the "ground truth" computation - no delta traversal.
 * Used for validation and testing.
 */
struct GroundTruthMetrics {
    // Read stats
    size_t readUniqueSeedCount = 0;
    int64_t totalReadSeedFrequency = 0;
    double readMagnitude = 0.0;
    
    // Genome stats
    size_t genomeUniqueSeedCount = 0;
    int64_t genomeTotalSeedFrequency = 0;
    double genomeMagnitudeSquared = 0.0;
    
    // Intersection stats
    size_t presenceIntersectionCount = 0;
    int64_t jaccardNumerator = 0;
    int64_t weightedJaccardNumerator = 0;
    double cosineNumerator = 0.0;
    
    // Final scores
    double jaccardScore = 0.0;
    double weightedJaccardScore = 0.0;
    double cosineScore = 0.0;
    double presenceJaccardScore = 0.0;
    
    /**
     * @brief Compute all metrics from read and genome seed sets
     * 
     * @param readSeeds Hash map of seed -> count from reads
     * @param genomeSeeds Hash map of seed -> count from genome
     * @return GroundTruthMetrics with all computed values
     */
    static GroundTruthMetrics compute(
        const SeedCountMap& readSeeds,
        const SeedCountMap& genomeSeeds);
};

// ============================================================================
// Node Lookup
// ============================================================================

/**
 * @brief Find node index by string ID in LiteTree
 * 
 * @param liteTree The LiteTree to search
 * @param nodeId The node identifier to find
 * @return Node index if found, -1 otherwise
 */
int32_t findNodeIndex(
    const panmapUtils::LiteTree& liteTree,
    const std::string& nodeId);

} // namespace indexUtils

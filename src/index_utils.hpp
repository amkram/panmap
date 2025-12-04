#pragma once

/**
 * @file index_utils.hpp
 * @brief Shared utilities for index traversal, seed reconstruction, and metrics
 * 
 * This module provides the canonical implementations of core operations used by
 * both the panmap workflow (indexing.cpp, placement.cpp) and correctness tests.
 * 
 * KEY DESIGN PRINCIPLE: Same code for production and testing
 * =========================================================
 * 
 * These functions are the source of truth for:
 * - Delta-based seed reconstruction (validated by test_index_dfs_state_matches_genome)
 * - Ground truth metrics computation (validated by test_placement_metrics_at_truth_node)
 * - Parent/child consistency verification (validated by test_seed_change_parent_child_consistency)
 * 
 * Production Usage (placement.cpp):
 *   - Uses computeChildMetrics() for incremental BFS traversal
 *   - Uses findNodeIndex() for resolving node IDs
 *   - Uses buildPathToNode() for path construction
 * 
 * Test Usage (correctness_tests.cpp):
 *   - Uses GroundTruthMetrics::compute() to verify delta-based results
 *   - Uses extractSeedsFromGenome() for direct genome seed extraction  
 *   - Uses compareSeedMaps() to verify reconstructed seeds match truth
 *   - Uses verifyParentChildConsistency() to validate index integrity
 * 
 * The overlap ensures that bugs found in testing are bugs in production.
 */

#include "panmap_utils.hpp"
#include "seeding.hpp"
#include "index_lite.capnp.h"
#include <absl/container/flat_hash_map.h>
#include <capnp/message.h>
#include <cstdint>
#include <string>
#include <vector>
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
     * This is the canonical reference implementation. All delta-based
     * computations should produce identical results.
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
// Path Building
// ============================================================================

/**
 * @brief Build path from root to target node
 * 
 * Traverses parent pointers from target to root, then reverses.
 * 
 * @param liteTree The LiteTree containing node structure
 * @param targetNodeIdx Index of the target node
 * @return Vector of node pointers from root to target (inclusive)
 */
std::vector<panmapUtils::LiteNode*> buildPathToNode(
    const panmapUtils::LiteTree& liteTree,
    int32_t targetNodeIdx);

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

// ============================================================================
// Seed Reconstruction via Delta Traversal
// ============================================================================

/**
 * @brief Reconstruct seed counts at a target node via delta traversal
 * 
 * Applies seed change deltas from root to target node, reconstructing
 * the exact seed state at that node. This is the core operation that
 * must match direct seed extraction from the genome.
 * 
 * VALIDATED BY: test_index_dfs_state_matches_genome
 * 
 * @param pathToTarget Path from root to target (from buildPathToNode)
 * @param seedChangeHashes Index array of seed hashes
 * @param seedChangeChildCounts Index array of child counts
 * @param nodeChangeOffsets Per-node offsets into change arrays
 * @return SeedCountMap with reconstructed seed counts at target
 */
SeedCountMap reconstructSeedsViaDeltas(
    const std::vector<panmapUtils::LiteNode*>& pathToTarget,
    const ::capnp::List<uint64_t>::Reader& seedChangeHashes,
    const ::capnp::List<int64_t>::Reader& seedChangeChildCounts,
    const ::capnp::List<uint32_t>::Reader& nodeChangeOffsets);

/**
 * @brief Overload accepting raw spans (for use during BFS traversal)
 */
SeedCountMap reconstructSeedsViaDeltas(
    const std::vector<panmapUtils::LiteNode*>& pathToTarget,
    const std::vector<std::tuple<uint64_t, int64_t, int64_t>>& allSeedChanges,
    const std::vector<uint32_t>& nodeChangeOffsets);

// ============================================================================
// Direct Seed Extraction
// ============================================================================

/**
 * @brief Extract seeds directly from a genome sequence
 * 
 * This is the ground truth extraction - syncmers directly from sequence.
 * 
 * @param genome The genome sequence
 * @param k Syncmer k parameter
 * @param s Syncmer s parameter
 * @param open Whether to use open syncmers
 * @param t Syncmer t parameter
 * @return SeedCountMap with hash -> count for each seed
 */
SeedCountMap extractSeedsFromGenome(
    const std::string& genome,
    int k, int s, bool open = false, int t = 0);

/**
 * @brief Extract seeds from multiple read sequences
 * 
 * Combines seeds from all reads, summing counts.
 * 
 * @param reads Vector of read sequences
 * @param k Syncmer k parameter
 * @param s Syncmer s parameter  
 * @param open Whether to use open syncmers
 * @param t Syncmer t parameter
 * @return SeedCountMap with combined seed counts
 */
SeedCountMap extractSeedsFromReads(
    const std::vector<std::string>& reads,
    int k, int s, bool open = false, int t = 0);

// ============================================================================
// Metric Computation from Seed Counts
// ============================================================================

/**
 * @brief Compute genome-only statistics from seed counts
 * 
 * @param seedCounts The seed count map
 * @param[out] uniqueCount Number of unique seeds
 * @param[out] totalFrequency Sum of all counts
 * @param[out] magnitudeSquared Sum of count^2
 */
void computeGenomeStats(
    const SeedCountMap& seedCounts,
    size_t& uniqueCount,
    int64_t& totalFrequency,
    double& magnitudeSquared);

/**
 * @brief Compute read statistics including magnitude
 * 
 * @param seedCounts The seed count map
 * @param[out] uniqueCount Number of unique seeds
 * @param[out] totalFrequency Sum of all counts
 * @param[out] magnitude Square root of sum of count^2
 */
void computeReadStats(
    const SeedCountMap& seedCounts,
    size_t& uniqueCount,
    int64_t& totalFrequency,
    double& magnitude);

// ============================================================================
// Validation Utilities
// ============================================================================

/**
 * @brief Compare two seed count maps for equality
 * 
 * @param a First seed count map
 * @param b Second seed count map
 * @param[out] missing Seeds in b but not in a
 * @param[out] extra Seeds in a but not in b
 * @param[out] countMismatch Seeds in both but with different counts
 * @return true if maps are identical
 */
bool compareSeedMaps(
    const SeedCountMap& a,
    const SeedCountMap& b,
    int& missing,
    int& extra,
    int& countMismatch);

/**
 * @brief Verify parent/child count consistency along a path
 * 
 * Checks that parentCount in seed changes matches the accumulated
 * state from prior nodes. This validates index integrity.
 * 
 * VALIDATED BY: test_seed_change_parent_child_consistency
 * 
 * @param pathToTarget Path from root to target
 * @param seedChangeHashes Index array of seed hashes
 * @param seedChangeParentCounts Index array of parent counts
 * @param seedChangeChildCounts Index array of child counts
 * @param nodeChangeOffsets Per-node offsets into change arrays
 * @return Number of inconsistencies found (0 = valid)
 */
int verifyParentChildConsistency(
    const std::vector<panmapUtils::LiteNode*>& pathToTarget,
    const ::capnp::List<uint64_t>::Reader& seedChangeHashes,
    const ::capnp::List<int64_t>::Reader& seedChangeParentCounts,
    const ::capnp::List<int64_t>::Reader& seedChangeChildCounts,
    const ::capnp::List<uint32_t>::Reader& nodeChangeOffsets);

} // namespace indexUtils

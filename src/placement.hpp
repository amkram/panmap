#pragma once

#include "panman.hpp"
#include "panmap_utils.hpp"
#include "seeding.hpp"

#include <absl/container/flat_hash_map.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <span>
#include <string>
#include <vector>

// OPTIMIZATION: Identity hash for seed hashes (they're already well-distributed)
// This avoids double-hashing and improves cache performance
struct IdentityHash {
    std::size_t operator()(uint64_t key) const noexcept {
        return static_cast<std::size_t>(key);
    }
};

namespace placement {

// Forward declarations
struct PlacementResult;
struct PlacementGlobalState;

// Parameters for lite placement traversal
struct TraversalParams {
  int k = 32;               // k-mer size
  int s = 8;                // syncmer parameter s
  int t = 0;                // t-syncmer parameter
  int l = 1;                // k-minimizer window size (1 = use raw syncmers)
  bool open = false;        // Whether to use open syncmers
  bool verify_scores = false; // Whether to recompute scores from scratch for verification
};

// Track global state during placement
struct PlacementGlobalState {
    // Seed hash frequency map from reads (computed once, used for all nodes)
    // OPTIMIZATION: Use IdentityHash since seed hashes are already well-distributed
    absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;
    
    int64_t totalReadSeedFrequency = 0;  // Sum of all read seed frequencies
    size_t readUniqueSeedCount = 0;       // Number of unique seeds in reads
    double readMagnitude = 0.0;           // Precomputed read magnitude for cosine similarity
    int kmerSize = 31;                    // Set from params

    // MGSR index data (deserialized once upfront)
    std::vector<panmapUtils::LiteNode*> liteNodes;  // For node ID lookups
    
    // Root node pointer for traversal
    panmapUtils::LiteNode* root = nullptr;
    
    // Optional: Full tree for verification mode (nullptr if not using verification)
    panmanUtils::Tree* fullTree = nullptr;
};

// Delta-based metrics computation (efficient forward-only traversal)
struct NodeMetrics {
    // Read-genome interaction metrics (computed incrementally)
    int64_t jaccardNumerator = 0;          // Presence-based intersection count (same as presenceIntersectionCount)
    int64_t weightedJaccardNumerator = 0;  // Σmin(readCount, genomeCount)
    double cosineNumerator = 0.0;          // Σ(readCount × genomeCount) for seeds in both
    size_t presenceIntersectionCount = 0;  // Count of seeds in both read and genome
    double rawMatchScore = 0.0;            // Σ(readCount / genomeCount) for seeds in intersection
    
    // PRE-INDEXED genome-only metrics (loaded from index, constant per node)
    double genomeMagnitudeSquared = 0.0;   // Σ(genomeCount²) for cosine denominator
    size_t genomeUniqueSeedCount = 0;      // Total unique seeds in genome
    int64_t genomeTotalSeedFrequency = 0;  // Σ(genomeCount) for weighted Jaccard denominator

    // Compute final scores
    double getJaccardScore(size_t readUniqueSeedCount) const {
        // PRESENCE-based Jaccard: |intersection| / |union|
        size_t genomeOnlyCount = (genomeUniqueSeedCount > presenceIntersectionCount) ? 
            (genomeUniqueSeedCount - presenceIntersectionCount) : 0;
        size_t totalUnion = readUniqueSeedCount + genomeOnlyCount;
        return (totalUnion > 0) ?
            static_cast<double>(jaccardNumerator) / totalUnion : 0.0;
    }

    double getWeightedJaccardScore(size_t totalReadSeedFrequency) const {
        // Weighted Jaccard: Σmin(r,g) / Σmax(r,g)
        int64_t denominator = totalReadSeedFrequency + genomeTotalSeedFrequency - weightedJaccardNumerator;
        return (denominator > 0) ?
            static_cast<double>(weightedJaccardNumerator) / denominator : 0.0;
    }

    double getCosineScore(double readMagnitude) const {
        // Standard cosine similarity: dot(read, genome) / (||read|| * ||genome||)
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (readMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = cosineNumerator / (readMagnitude * genomeMagnitude);
        // Clamp to [0, 1] - values slightly outside due to floating point errors
        return std::clamp(score, 0.0, 1.0);
    }

    double getPresenceJaccardScore(size_t readUniqueSeedCount) const {
        size_t genomeOnlyCount = (genomeUniqueSeedCount > presenceIntersectionCount) ? 
            (genomeUniqueSeedCount - presenceIntersectionCount) : 0;
        size_t totalUnion = readUniqueSeedCount + genomeOnlyCount;
        return (totalUnion > 0) ?
            static_cast<double>(presenceIntersectionCount) / totalUnion : 0.0;
    }

    // Compute child metrics from parent metrics + seed changes
    static void computeChildMetrics(
        NodeMetrics& childMetrics,
        std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges,
        const PlacementGlobalState& state);
};

// Consolidated score tracking and results for lite placement
struct PlacementResult {
  // Hits-based results
  int64_t maxHitsInAnyGenome = 0;
  uint32_t maxHitsNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedMaxHitsNodeIndices;

  // Raw Seed Match Score
  double bestRawSeedMatchScore = 0.0;
  uint32_t bestRawSeedMatchNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedRawSeedMatchNodeIndices;

  // Jaccard-based results
  double bestJaccardScore = 0.0;
  uint32_t bestJaccardNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedJaccardNodeIndices;

  // Weighted Jaccard results
  double bestWeightedJaccardScore = 0.0;
  uint32_t bestWeightedJaccardNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedWeightedJaccardNodeIndices;

  // Jaccard Index (Presence/Absence)
  double bestJaccardPresenceScore = 0.0;
  uint32_t bestJaccardPresenceNodeIndex = UINT32_MAX;

  // Cosine-based results
  double bestCosineScore = 0.0;
  uint32_t bestCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCosineNodeIndices;

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  
  // Update methods (take nodeIndex, resolve to string at end)
  void updateHitsScore(uint32_t nodeIndex, int64_t hits, panmapUtils::LiteNode* node = nullptr);
  void updateRawSeedMatchScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateJaccardScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateWeightedJaccardScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateJaccardPresenceScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  
  // Lazy resolution: convert winning node indices to string IDs (called ONCE at end)
  void resolveNodeIds(panmapUtils::LiteTree* liteTree);
  
  // Resolved string IDs (only populated after resolveNodeIds() is called)
  std::string bestJaccardNodeId;
  std::string bestCosineNodeId;
  std::string bestWeightedJaccardNodeId;
  std::string bestJaccardPresenceNodeId;
  std::string bestRawSeedMatchNodeId;
  std::string maxHitsNodeId;
  
  // ==========================================================================
  // Alignment data (minimal storage - paths only, data extracted on demand)
  // ==========================================================================
  
  // FASTQ file paths (for re-reading during alignment)
  std::string reads1Path;
  std::string reads2Path;
  
  // Seed hash frequency map from all reads (needed for reference seed filtering)
  absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;
  
  // Index parameters (needed for alignment)
  int k = 31;
  int s = 8;
  int t = 0;
  bool open = false;
};

// LiteTree-based placement (the main entry point for lite placement)
void placeLite(PlacementResult &result, 
               panmapUtils::LiteTree *liteTree,
               ::capnp::MessageReader &liteIndex, 
               const std::string &reads1,
               const std::string &reads2,
               std::string &outputPath,
               bool verify_scores = false,
               panmanUtils::Tree *fullTree = nullptr);

} // namespace placement
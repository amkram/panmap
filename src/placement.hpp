#pragma once

#include "panman.hpp"
#include "panmap_utils.hpp"
#include "seeding.hpp"

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
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
  bool store_diagnostics = false; // Store metric components on nodes for diagnostic output
  double seedMaskFraction = 0.001; // Mask top N fraction of seeds by frequency (0 = disabled, 0.001 = top 0.1%)
  int minSeedQuality = 0;   // Min avg Phred quality for seed region (0 = disabled)
  bool dedupReads = false;  // Deduplicate reads before placement (count each unique sequence once)
  bool pairFilter = true;   // Filter: only count seeds from read-pairs that have at least one matching seed at current node
  int trimStart = 0;        // Trim N bases from start of each read before seeding (for primer removal)
  int trimEnd = 0;          // Trim N bases from end of each read before seeding (for primer removal)
  int minReadSupport = 1;   // Minimum read count for a seed to be included (1 = all seeds, 2 = filter single-read seeds)
  bool hpc = false;          // Homopolymer-compressed seeds (from index)
  
  // ========================================
  // Alignment-based refinement parameters
  // After k-mer scoring, optionally refine top candidates by full alignment
  // ========================================
  bool refineEnabled = false;       // Enable alignment-based refinement
  double refineTopPct = 0.01;       // Top X% of nodes to consider (default 1%)
  int refineMaxTopN = 150;          // Max nodes to align against (caps the %)
  int refineNeighborRadius = 2;     // Expand to neighbors within N branches
  int refineMaxNeighborN = 150;     // Max additional nodes from neighbor expansion
};

// Track global state during placement
struct PlacementGlobalState {
    // Seed hash frequency map from reads (computed once, used for all nodes)
    // OPTIMIZATION: Use IdentityHash since seed hashes are already well-distributed
    absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;
    
    // Log-scaled read counts: log(1 + readCount) for each seed
    // Pre-computed once to avoid repeated log() calls during traversal
    absl::flat_hash_map<size_t, double, IdentityHash> logReadCounts;
    
    int64_t totalReadSeedFrequency = 0;  // Sum of all read seed frequencies
    size_t readUniqueSeedCount = 0;       // Number of unique seeds in reads
    double logReadMagnitude = 0.0;        // Precomputed log-scaled read magnitude: sqrt(Σlog(1+r)²)
    
    // Hash to sequence map for debugging output (optional, populated when verbose)
    absl::flat_hash_map<size_t, std::string, IdentityHash> hashToSequence;
    
    // Precomputed inverse genome counts: 1/g_i for each seed in reads
    // g_i = number of genomes in the pangenome containing seed i (from root)
    absl::flat_hash_map<size_t, double, IdentityHash> seedInverseGenomeCounts;
    double weightedContainmentDenominator = 0.0;  // Σ(1/g_i) for all seeds in reads
    double logContainmentDenominator = 0.0;        // Σ log(1+r_i) for all seeds in reads

    int kmerSize = 31;                    // Set from params
    
    // MGSR index data (deserialized once upfront)
    std::vector<panmapUtils::LiteNode*> liteNodes;  // For node ID lookups
    
    // Root node pointer for traversal
    panmapUtils::LiteNode* root = nullptr;
    
    // Optional: Full tree for verification mode (nullptr if not using verification)
    panmanUtils::Tree* fullTree = nullptr;
    
    // Leave-one-out validation: node index to skip during scoring (UINT32_MAX = none)
    uint32_t skipNodeIndex = UINT32_MAX;
};

// Delta-based metrics computation (efficient forward-only traversal)
struct NodeMetrics {
    // Read-genome interaction metrics (computed incrementally)
    double logRawNumerator = 0.0;          // Σ(log(1+readCount) / genomeCount) for LogRAW metric
    double logCosineNumerator = 0.0;       // Σ(log(1+readCount) × genomeCount) for LogCosine metric
    size_t presenceIntersectionCount = 0;  // Count of seeds present in both read and genome
    double weightedContainmentNumerator = 0.0; // Σ(1/nodeGenomeCount_i) for seeds in reads ∩ genome
    double logContainmentNumerator = 0.0;       // Σ log(1+r_i) for seeds in reads ∩ genome

    // PRE-INDEXED genome-only metrics (loaded from index, updated incrementally)
    double genomeMagnitudeSquared = 0.0;   // Σ(genomeCount²) for cosine denominator
    size_t genomeUniqueSeedCount = 0;      // Total unique seeds in genome

    // LogRAW Score: Σ(log(1+readCount) / genomeCount) for matching seeds
    double getLogRawScore(double logReadMagnitude) const {
        if (logReadMagnitude <= 0.0) return 0.0;
        return logRawNumerator / logReadMagnitude;
    }

    // LogCosine Score: Σ(log(1+readCount) × log(1+genomeCount)) / (logReadMag × logGenomeMag)
    double getLogCosineScore(double logReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (logReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = logCosineNumerator / (logReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }

    // Containment Index: |reads ∩ genome| / |reads|
    double getContainmentScore(size_t readUniqueSeedCount) const {
        return (readUniqueSeedCount > 0) ?
            static_cast<double>(presenceIntersectionCount) / readUniqueSeedCount : 0.0;
    }

    // Weighted Containment: Σ(1/nodeGenomeCount_i for R∩G) / Σ(1/rootGenomeCount_i for R)
    // Upweights rare (discriminative) seeds, downweights common ones
    double getWeightedContainmentScore(double weightedContainmentDenominator) const {
        return (weightedContainmentDenominator > 0.0) ?
            weightedContainmentNumerator / weightedContainmentDenominator : 0.0;
    }

    // Log Containment: Σ(log(1+r_i) for R∩G) / Σ(log(1+r_i) for R)
    // Weights each seed by log read abundance instead of binary presence
    double getLogContainmentScore(double logContainmentDenominator) const {
        return (logContainmentDenominator > 0.0) ?
            logContainmentNumerator / logContainmentDenominator : 0.0;
    }

    // Compute child metrics from parent metrics + seed changes
    static void computeChildMetrics(
        NodeMetrics& childMetrics,
        std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges,
        const PlacementGlobalState& state);
};

// Consolidated score tracking and results for lite placement
struct PlacementResult {
  // LogRAW Score results (log-scaled, coverage-robust)
  double bestLogRawScore = 0.0;
  uint32_t bestLogRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogRawNodeIndices;

  // LogCosine Score results (log-scaled cosine similarity)
  double bestLogCosineScore = 0.0;
  uint32_t bestLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogCosineNodeIndices;

  // Containment Index results
  double bestContainmentScore = 0.0;
  uint32_t bestContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedContainmentNodeIndices;

  // Weighted Containment Index results
  double bestWeightedContainmentScore = 0.0;
  uint32_t bestWeightedContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedWeightedContainmentNodeIndices;

  // Log Containment Index results
  double bestLogContainmentScore = 0.0;
  uint32_t bestLogContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogContainmentNodeIndices;

  // Per-metric alignment-based refinement results
  // Each seed metric's top candidates are refined independently by full alignment
  struct RefinedResult {
    double score = 0.0;
    uint32_t nodeIndex = UINT32_MAX;
    std::string nodeId;
  };
  RefinedResult refinedLogRaw;
  RefinedResult refinedLogCosine;
  RefinedResult refinedContainment;
  RefinedResult refinedWeightedContainment;
  RefinedResult refinedLogContainment;
  bool refinementWasRun = false;          // True if refinement was executed

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  
  // Update methods (take nodeIndex, resolve to string at end)
  void updateLogRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateWeightedContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  
  // Lazy resolution: convert winning node indices to string IDs (called ONCE at end)
  void resolveNodeIds(panmapUtils::LiteTree* liteTree);
  
  // Resolved string IDs (only populated after resolveNodeIds() is called)
  std::string bestLogRawNodeId;
  std::string bestLogCosineNodeId;
  std::string bestContainmentNodeId;
  std::string bestWeightedContainmentNodeId;
  std::string bestLogContainmentNodeId;
  
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
  
  // Read statistics (for diagnostic output)
  size_t readUniqueSeedCount = 0;
  int64_t totalReadSeedFrequency = 0;
  double readMagnitude = 0.0;
};

// LiteTree-based placement (the main entry point for lite placement)
// If skipNodeId is non-empty, that node will be excluded from placement results
// (for leave-one-out validation). Its parent ID is returned via parentOfSkippedNode.
void placeLite(PlacementResult &result, 
               panmapUtils::LiteTree *liteTree,
               ::capnp::MessageReader &liteIndex, 
               const std::string &reads1,
               const std::string &reads2,
               std::string &outputPath,
               bool verify_scores = false,
               panmanUtils::Tree *fullTree = nullptr,
               const std::string &skipNodeId = "",
               std::string *parentOfSkippedNode = nullptr,
               bool store_diagnostics = false,
               double seedMaskFraction = 0.001,
               int minSeedQuality = 0,
               bool dedupReads = false,
               bool pairFilter = true,
               int trimStart = 0,
               int trimEnd = 0,
               int minReadSupport = 1,
               bool refineEnabled = false,
               double refineTopPct = 0.01,
               int refineMaxTopN = 150,
               int refineNeighborRadius = 2,
               int refineMaxNeighborN = 150);

} // namespace placement

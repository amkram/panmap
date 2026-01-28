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
};

// Track global state during placement
// SIMPLIFIED: Only fields needed for LogRaw and LogCosine metrics
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
    int kmerSize = 31;                    // Set from params
    
    // ========================================
    // Paired-end fragment tracking for concordance-based scoring
    // Fragments where BOTH mates have matching seeds are weighted higher
    // Distance between R1/R2 matches is checked against expected fragment size
    // ========================================
    bool isPairedEnd = false;  // True if reads2 was provided
    size_t numFragments = 0;   // Number of read-pairs (or single reads if unpaired)
    uint32_t expectedFragmentSize = 400;   // Expected fragment size (insert size + read lengths)
    uint32_t fragmentSizeTolerance = 150;  // Allowed deviation from expected size
    
    // For each seed hash, list of fragment IDs that contain this seed
    // Used to find which fragments have matches at a given node
    absl::flat_hash_map<uint64_t, std::vector<uint32_t>, IdentityHash> seedToFragments;
    
    // For each fragment, seeds from R1 and R2 separately (for concordance checking)
    // A fragment is "concordant" at a node if both R1 and R2 have ≥1 matching seed
    std::vector<absl::flat_hash_set<uint64_t, IdentityHash>> fragmentR1Seeds;  // Seeds from read 1
    std::vector<absl::flat_hash_set<uint64_t, IdentityHash>> fragmentR2Seeds;  // Seeds from read 2
    
    // For each fragment, store read lengths for distance calculation
    std::vector<std::pair<uint16_t, uint16_t>> fragmentReadLengths;  // (R1 length, R2 length)

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
// SIMPLIFIED: Only LogRaw and LogCosine metrics are computed
struct NodeMetrics {
    // Read-genome interaction metrics (computed incrementally)
    double logRawNumerator = 0.0;          // Σ(log(1+readCount) / genomeCount) for LogRAW metric
    double logCosineNumerator = 0.0;       // Σ(log(1+readCount) × genomeCount) for LogCosine metric
    
    // PRE-INDEXED genome-only metrics (loaded from index, constant per node)
    double genomeMagnitudeSquared = 0.0;   // Σ(genomeCount²) for cosine denominator
    size_t genomeUniqueSeedCount = 0;      // Total unique seeds in genome
    int64_t genomeTotalSeedFrequency = 0;  // Σ(genomeCount) for weighted Jaccard denominator (kept for potential future use)

    // LogRAW Score: Σ(log(1+readCount) / genomeCount) for matching seeds
    // - Log-scales read counts to reduce impact of variable coverage
    // - Divides by genomeCount to penalize repeated seeds
    // Normalized by log-scaled read magnitude for scale-invariance
    double getLogRawScore(double logReadMagnitude) const {
        return (logReadMagnitude > 0.0) ?
            logRawNumerator / logReadMagnitude : 0.0;
    }

    // LogCosine Score: Σ(log(1+readCount) × genomeCount) / (logReadMag × genomeMag)
    // Log-scaled cosine similarity - reduces impact of high-coverage seeds
    double getLogCosineScore(double logReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (logReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = logCosineNumerator / (logReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }

    // Compute child metrics from parent metrics + seed changes
    static void computeChildMetrics(
        NodeMetrics& childMetrics,
        std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges,
        const PlacementGlobalState& state);
};

// Consolidated score tracking and results for lite placement
// SIMPLIFIED: Only LogRaw and LogCosine metrics are tracked
struct PlacementResult {
  // LogRAW Score results (log-scaled, coverage-robust)
  double bestLogRawScore = 0.0;
  uint32_t bestLogRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogRawNodeIndices;

  // LogCosine Score results (log-scaled cosine similarity)
  double bestLogCosineScore = 0.0;
  uint32_t bestLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogCosineNodeIndices;

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  
  // Update methods (take nodeIndex, resolve to string at end)
  void updateLogRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  
  // Lazy resolution: convert winning node indices to string IDs (called ONCE at end)
  void resolveNodeIds(panmapUtils::LiteTree* liteTree);
  
  // Resolved string IDs (only populated after resolveNodeIds() is called)
  std::string bestLogRawNodeId;
  std::string bestLogCosineNodeId;
  
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
               uint32_t expectedFragmentSize = 400,
               uint32_t fragmentSizeTolerance = 150);

} // namespace placement
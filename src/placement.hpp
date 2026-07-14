#pragma once

#include "panman.hpp"
#include "panmap_utils.hpp"
#include "seeding.hpp"

#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <span>
#include <string>
#include <vector>

// Identity hash for seed hashes (already well-distributed): avoids double-hashing.
struct IdentityHash {
    constexpr std::size_t operator()(uint64_t key) const noexcept { return static_cast<std::size_t>(key); }
};

namespace placement {

struct PlacementResult;
struct PlacementGlobalState;

struct TraversalParams {
    int k = 19;                       // k-mer size
    int s = 8;
    int t = 0;
    int l = 3;                        // 1 = raw syncmers, not k-minimizers
    bool open = false;
    bool verify_scores = false;
    bool store_diagnostics = false;
    double seedMaskFraction = 0.001;  // 0 = disabled; else mask top fraction of seeds by frequency
    int minSeedQuality = 0;           // min avg Phred over seed span; 0 = disabled
    bool dedupReads = false;
    bool pairFilter =
        true;           // only count seeds from read-pairs with >=1 matching seed at current node
    int trimStart = 0;  // bases trimmed from read start (primer removal)
    int trimEnd = 0;    // bases trimmed from read end (primer removal)
    int minReadSupport =
        -1;            // -1 = auto (2 if est. coverage > 3x, else 1); 1 = keep all; 2 = drop singletons
    bool hpc = false;

    // Alignment-based refinement: after k-mer scoring, optionally refine top candidates by full alignment
    bool refineEnabled = false;
    double refineTopPct = 0.01;    // Top X% of nodes to consider (default 1%)
    int refineMaxTopN = 150;       // Max nodes to align against (caps the %)
    int refineNeighborRadius = 2;  // Expand to neighbors within N branches
    int refineMaxNeighborN = 150;  // Max additional nodes from neighbor expansion
    bool forceLeaf = false;        // Restrict scoring to leaf nodes only
};

struct PlacementGlobalState {
    // Seed hash frequency map from reads (computed once, used for all nodes)
    absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;

    // Log-scaled read counts: log(1 + readCount) for each seed
    // Pre-computed once to avoid repeated log() calls during traversal
    absl::flat_hash_map<size_t, double, IdentityHash> logReadCounts;

    int64_t totalReadSeedFrequency = 0;
    size_t readUniqueSeedCount = 0;
    double logReadMagnitude = 0.0;       // Precomputed log-scaled read magnitude: sqrt(Σlog(1+r)²)

    // Hash to sequence map for debugging output (optional, populated when verbose)
    absl::flat_hash_map<size_t, std::string, IdentityHash> hashToSequence;

    // Precomputed inverse genome counts: 1/g_i for each seed in reads
    // g_i = number of genomes in the pangenome containing seed i (from root)
    absl::flat_hash_map<size_t, double, IdentityHash> seedInverseGenomeCounts;
    double weightedContainmentDenominator = 0.0;  // Σ(1/g_i) for all seeds in reads
    double logContainmentDenominator = 0.0;       // Σ log(1+r_i) for all seeds in reads

    int kmerSize = 31;

    // MGSR index data (deserialized once upfront)
    std::vector<panmapUtils::LiteNode*> liteNodes;  // For node ID lookups

    panmapUtils::LiteNode* root = nullptr;

    // Full tree for verification mode (nullptr if unused)
    panmanUtils::Tree* fullTree = nullptr;

    // Lite tree that owns the zero-copy seed-change SoA pointers (set before traversal)
    const panmapUtils::LiteTree* liteTree = nullptr;

    // Leave-one-out validation: node index to skip during scoring (UINT32_MAX = none)
    uint32_t skipNodeIndex = UINT32_MAX;

    // Restrict scoring to leaf nodes only
    bool forceLeaf = false;
};

// Effective min read-count support for a seed. configuredMinSupport < 0 means "auto":
// 2 when est. coverage (mean read count over seeds seen in >=2 reads) exceeds 3x, else 1.
int64_t resolveMinReadSupport(const absl::flat_hash_map<size_t, int64_t, IdentityHash>& seedFreqInReads,
                              int configuredMinSupport);

// Fill state's read-derived scoring fields (logReadCounts, totalReadSeedFrequency,
// readUniqueSeedCount, logReadMagnitude, logContainmentDenominator) from seedFreqInReads,
// keeping only seeds with read count >= minSupport. Returns count of dropped seeds.
size_t computeReadSeedMagnitudes(PlacementGlobalState& state, int64_t minSupport);

// Delta-based metrics (forward-only traversal)
struct NodeMetrics {
    double logRawNumerator = 0.0;               // Σ(log(1+readCount) / genomeCount) for LogRAW metric
    double logCosineNumerator = 0.0;            // Σ(log(1+readCount) × genomeCount) for LogCosine metric
    size_t presenceIntersectionCount = 0;       // Count of seeds present in both read and genome
    double weightedContainmentNumerator = 0.0;  // Σ(1/nodeGenomeCount_i) for seeds in reads ∩ genome
    double logContainmentNumerator = 0.0;       // Σ log(1+r_i) for seeds in reads ∩ genome

    // Pre-indexed genome-only metrics (loaded from index, updated incrementally)
    double genomeMagnitudeSquared = 0.0;  // Σ(genomeCount²) for cosine denominator
    size_t genomeUniqueSeedCount = 0;

    // LogRAW Score: Σ(log(1+readCount) / genomeCount) for matching seeds
    constexpr double getLogRawScore(double logReadMagnitude) const {
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
    constexpr double getContainmentScore(size_t readUniqueSeedCount) const {
        return (readUniqueSeedCount > 0) ? static_cast<double>(presenceIntersectionCount) / readUniqueSeedCount : 0.0;
    }

    // Weighted Containment: Σ(1/nodeGenomeCount_i for R∩G) / Σ(1/rootGenomeCount_i for R)
    // Upweights rare (discriminative) seeds, downweights common ones
    constexpr double getWeightedContainmentScore(double weightedContainmentDenominator) const {
        return (weightedContainmentDenominator > 0.0) ? weightedContainmentNumerator / weightedContainmentDenominator
                                                      : 0.0;
    }

    // Log Containment: Σ(log(1+r_i) for R∩G) / Σ(log(1+r_i) for R)
    // Weights each seed by log read abundance instead of binary presence
    constexpr double getLogContainmentScore(double logContainmentDenominator) const {
        return (logContainmentDenominator > 0.0) ? logContainmentNumerator / logContainmentDenominator : 0.0;
    }

    static void computeChildMetrics(NodeMetrics& childMetrics,
                                    uint64_t seedChangeOffset,
                                    uint32_t seedChangeSize,
                                    const PlacementGlobalState& state);
};

struct PlacementResult {
    double bestLogRawScore = 0.0;
    uint32_t bestLogRawNodeIndex = UINT32_MAX;
    std::vector<uint32_t> tiedLogRawNodeIndices;

    double bestLogCosineScore = 0.0;
    uint32_t bestLogCosineNodeIndex = UINT32_MAX;
    std::vector<uint32_t> tiedLogCosineNodeIndices;

    double bestContainmentScore = 0.0;
    uint32_t bestContainmentNodeIndex = UINT32_MAX;
    std::vector<uint32_t> tiedContainmentNodeIndices;

    double bestWeightedContainmentScore = 0.0;
    uint32_t bestWeightedContainmentNodeIndex = UINT32_MAX;
    std::vector<uint32_t> tiedWeightedContainmentNodeIndices;

    double bestLogContainmentScore = 0.0;
    uint32_t bestLogContainmentNodeIndex = UINT32_MAX;
    std::vector<uint32_t> tiedLogContainmentNodeIndices;

    // Per-node seed scores by DFS index (populated only when refinement enabled). Kept here
    // not on the shared LiteNode so concurrent batch placements don't race. Metric order:
    // {logRaw, logCosine, containment, weightedContainment, logContainment}.
    std::vector<std::array<float, 5>> nodeScores;

    // Per-metric alignment-based refinement: each metric's top candidates refined independently by full alignment
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
    bool refinementWasRun = false;

    int64_t totalReadsProcessed = 0;

    // Take nodeIndex, resolve to string at end
    void updateLogRawScore(uint32_t nodeIndex, double score);
    void updateLogCosineScore(uint32_t nodeIndex, double score);
    void updateContainmentScore(uint32_t nodeIndex, double score);
    void updateWeightedContainmentScore(uint32_t nodeIndex, double score);
    void updateLogContainmentScore(uint32_t nodeIndex, double score);

    // Lazy resolution: convert winning node indices to string IDs (called once at end)
    void resolveNodeIds(panmapUtils::LiteTree* liteTree);

    // Resolved string IDs (only populated after resolveNodeIds() is called)
    std::string bestLogRawNodeId;
    std::string bestLogCosineNodeId;
    std::string bestContainmentNodeId;
    std::string bestWeightedContainmentNodeId;
    std::string bestLogContainmentNodeId;

    // Alignment data: paths only, data extracted on demand.

    // FASTQ file paths (for re-reading during alignment)
    std::string reads1Path;
    std::string reads2Path;

    // Seed hash frequency map from all reads (needed for reference seed filtering)
    absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;

    // Index parameters (needed for alignment)
    int k = 19;
    int s = 8;
    int t = 0;
    bool open = false;

    // Read statistics (for diagnostic output)
    size_t readUniqueSeedCount = 0;
    int64_t totalReadSeedFrequency = 0;
    double readMagnitude = 0.0;
};

void placeLite(PlacementResult& result,
               panmapUtils::LiteTree* liteTree,
               ::capnp::MessageReader& liteIndex,
               const std::string& reads1,
               const std::string& reads2,
               std::string& outputPath,
               const TraversalParams& params = {},
               panmanUtils::Tree* fullTree = nullptr);

}  // namespace placement

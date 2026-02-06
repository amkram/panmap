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
// SIMPLIFIED: Only fields needed for LogRaw and LogCosine metrics
struct PlacementGlobalState {
    // Seed hash frequency map from reads (computed once, used for all nodes)
    // OPTIMIZATION: Use IdentityHash since seed hashes are already well-distributed
    absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;
    
    // Log-scaled read counts: log(1 + readCount) for each seed
    // Pre-computed once to avoid repeated log() calls during traversal
    absl::flat_hash_map<size_t, double, IdentityHash> logReadCounts;
    
    // Double-log-scaled read counts: log(1 + log(1 + readCount)) for each seed
    // Even more dampened than single-log for extreme coverage variation
    absl::flat_hash_map<size_t, double, IdentityHash> logLogReadCounts;
    
    int64_t totalReadSeedFrequency = 0;  // Sum of all read seed frequencies
    size_t readUniqueSeedCount = 0;       // Number of unique seeds in reads
    double logReadMagnitude = 0.0;        // Precomputed log-scaled read magnitude: sqrt(Σlog(1+r)²)
    double logLogReadMagnitude = 0.0;     // Precomputed double-log-scaled magnitude: sqrt(Σlog(1+log(1+r))²)
    double readMagnitude = 0.0;           // Precomputed raw read magnitude: sqrt(Σr²)
    double cappedReadMagnitude = 0.0;     // Precomputed capped magnitude: sqrt(Σmin(r,100)²)
    double cappedLogReadMagnitude = 0.0;  // Precomputed capped log magnitude: sqrt(Σlog(1+min(r,100))²)
    
    // Capped read counts for capped metrics
    absl::flat_hash_map<size_t, int64_t, IdentityHash> cappedReadCounts;
    absl::flat_hash_map<size_t, double, IdentityHash> cappedLogReadCounts;
    
    // Hash to sequence map for debugging output (optional, populated when verbose)
    absl::flat_hash_map<size_t, std::string, IdentityHash> hashToSequence;
    
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
struct NodeMetrics {
    // ========================================
    // OLD METRICS (from working version 6cfb6cf)
    // ========================================
    double rawNumerator = 0.0;             // Σ(readCount / genomeCount) for RAW metric
    double cosineNumerator = 0.0;          // Σ(readCount × genomeCount) for Cosine metric
    int64_t weightedJaccardNumerator = 0;  // Σmin(readCount, genomeCount)
    size_t presenceIntersectionCount = 0;  // Count of seeds in both read and genome
    double capCosineNumerator = 0.0;       // Σ(min(readCount,100) × genomeCount) for capped cosine
    double capLogCosineNumerator = 0.0;    // Σ(log(1+min(readCount,100)) × genomeCount)
    
    // Read-genome interaction metrics (computed incrementally)
    double logRawNumerator = 0.0;          // Σ(log(1+readCount) / genomeCount) for LogRAW metric
    double logCosineNumerator = 0.0;       // Σ(log(1+readCount) × genomeCount) for LogCosine metric (OLD FORMULA)
    double logLogRawNumerator = 0.0;       // Σ(log(1+log(1+readCount)) / genomeCount) for LogLogRAW metric
    double logLogCosineNumerator = 0.0;    // Σ(log(1+log(1+readCount)) × genomeCount) for LogLogCosine metric
    
    // COVERAGE TRACKING: Track read seed mass that is covered by this genome
    // Used to penalize genomes that don't cover the read seeds (fixes log_cosine catastrophic failures)
    double matchedLogReadMagnitudeSquared = 0.0;     // Σ(log²(1+readCount)) for seeds PRESENT in genome
    double matchedLogLogReadMagnitudeSquared = 0.0;  // Σ(log²(1+log(1+readCount))) for seeds PRESENT in genome
    double matchedReadMagnitudeSquared = 0.0;        // Σ(readCount²) for seeds PRESENT in genome
    
    // IDF-WEIGHTED COSINE: Log-dampen genome counts to reduce repeat amplification
    // Instead of Σ(log(R) × G), use Σ(log(R) × log(1+G)) - proper TF-IDF style
    double logCosineIdfNumerator = 0.0;              // Σ(log(1+readCount) × log(1+genomeCount))
    double logGenomeMagnitudeSquared = 0.0;          // Σ(log(1+genomeCount)²) for IDF denominator
    
    // PRESENCE SCORE: Low-coverage robust metric (ignores read frequency, just binary presence)
    // Ideal for sparse sampling where most k-mers appear only 1-2 times in reads
    // Formula: Σ(1/genomeCount) for each UNIQUE read k-mer present in genome
    double presenceNumerator = 0.0;                  // Σ(1/G_i) for matched unique seeds
    size_t matchedUniqueReadSeeds = 0;               // Count of unique read k-mers found in genome
    
    // CONCORDANCE TRACKING: Track which read seeds are present in genome
    // Used to compute fragment-level concordance (R1/R2 balance)
    absl::flat_hash_set<uint64_t, IdentityHash> genomeSeedSet;  // Seeds from reads that are in genome
    
    // PRE-INDEXED genome-only metrics (loaded from index, constant per node)
    double genomeMagnitudeSquared = 0.0;   // Σ(genomeCount²) for cosine denominator
    size_t genomeUniqueSeedCount = 0;      // Total unique seeds in genome
    int64_t genomeTotalSeedFrequency = 0;  // Σ(genomeCount) for weighted Jaccard denominator

    // ========================================
    // OLD METRIC GETTERS (verified working formulas from 6cfb6cf)
    // ========================================
    
    // RAW Score: Σ(readCount / genomeCount) / readMagnitude
    double getRawScore(double readMagnitude) const {
        if (readMagnitude <= 0.0) return 0.0;
        return rawNumerator / readMagnitude;
    }
    
    // Cosine Score: Σ(readCount × genomeCount) / (readMag × genomeMag)
    // Standard cosine similarity - the baseline metric
    double getCosineScore(double readMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (readMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = cosineNumerator / (readMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }
    
    // OLD LogCosine Score: Σ(log(1+readCount) × genomeCount) / (logReadMag × genomeMag)
    // This is the ORIGINAL formula that works - log on reads only, raw genome counts
    // NO coverage penalty, NO log on genome side
    double getLogCosineOldScore(double logReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (logReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = logCosineNumerator / (logReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }
    
    // Weighted Containment: Σ(readFreq for matching seeds) / Σ(readFreq for all seeds)
    double getWeightedContainmentScore(int64_t totalReadSeedFrequency) const {
        return (totalReadSeedFrequency > 0) ?
            cosineNumerator / static_cast<double>(totalReadSeedFrequency) : 0.0;
    }
    
    // Cap Cosine: Σ(min(readCount,100) × genomeCount) / (cappedMag × genomeMag)
    double getCapCosineScore(double cappedReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (cappedReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = capCosineNumerator / (cappedReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }
    
    // Cap Log Cosine: Σ(log(1+min(readCount,100)) × genomeCount) / (cappedLogMag × genomeMag)
    double getCapLogCosineScore(double cappedLogReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (cappedLogReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = capLogCosineNumerator / (cappedLogReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }
    
    // ========================================
    // NEW METRIC GETTERS (current version)
    // ========================================

    // LogRAW Score: Σ(log(1+readCount) / genomeCount) for matching seeds
    // - Log-scales read counts to reduce impact of variable coverage
    // - Divides by genomeCount to penalize repeated seeds
    // Normalized by log-scaled read magnitude for scale-invariance
    double getLogRawScore(double logReadMagnitude) const {
        if (logReadMagnitude <= 0.0) return 0.0;
        return logRawNumerator / logReadMagnitude;
    }

    // LogCosine Score: IDF-weighted cosine similarity
    // Formula: Σ(log(1+readCount) × log(1+genomeCount)) / (logReadMag × logGenomeMag)
    // Log-dampens BOTH read counts AND genome counts - true TF-IDF style
    // This prevents repeat amplification: 100 copies contribute log(101)≈4.6, not 100
    double getLogCosineScore(double logReadMagnitude) const {
        double logGenomeMagnitude = std::sqrt(logGenomeMagnitudeSquared);
        if (logReadMagnitude <= 0.0 || logGenomeMagnitude <= 0.0) return 0.0;
        double rawScore = logCosineIdfNumerator / (logReadMagnitude * logGenomeMagnitude);
        
        // Coverage penalty: what fraction of read seed mass is covered by genome?
        double logReadMagnitudeSquared = logReadMagnitude * logReadMagnitude;
        double coverageFraction = matchedLogReadMagnitudeSquared / logReadMagnitudeSquared;
        double coveragePenalty = std::sqrt(std::clamp(coverageFraction, 0.0, 1.0));
        
        return std::clamp(rawScore * coveragePenalty, 0.0, 1.0);
    }

    // LogLogRAW Score: Σ(log(1+log(1+readCount)) / genomeCount) for matching seeds
    // Double-log dampens coverage effects even more than single-log
    double getLogLogRawScore(double logLogReadMagnitude) const {
        if (logLogReadMagnitude <= 0.0) return 0.0;
        return logLogRawNumerator / logLogReadMagnitude;
    }

    // Containment Score: fraction of read seed mass that is present in genome
    // Simple and robust metric: high = genome covers most read seeds
    // Formula: sqrt(matchedMag²) / totalMag = matchedMag / totalMag
    double getContainmentScore(double logReadMagnitude) const {
        if (logReadMagnitude <= 0.0) return 0.0;
        double matchedMag = std::sqrt(std::max(0.0, matchedLogReadMagnitudeSquared));
        return std::clamp(matchedMag / logReadMagnitude, 0.0, 1.0);
    }

    // Presence Score: LOW-COVERAGE ROBUST metric
    // Designed for sparse sampling where most k-mers appear only 1-2 times in reads
    // Key insight: With low coverage, READ FREQUENCY (R) is unreliable (mostly 1)
    // So we ignore R entirely and just count BINARY presence, weighted by uniqueness
    // 
    // Formula: Σ(1/G_i) / numUniqueReadKmers
    //   - 1/G: unique markers (G=1) get full point, repeats (G=100) get 0.01
    //   - Dividing by total unique read seeds gives a rate [0,1]
    //
    // Why this works for low coverage:
    //   - Doesn't penalize low read counts (R is ignored)
    //   - Still strongly prefers unique markers over repeats
    //   - Normalizes by what you COULD have matched, not what the genome HAS
    double getPresenceScore(size_t totalUniqueReadSeeds) const {
        if (totalUniqueReadSeeds == 0) return 0.0;
        // Scale: if every matched seed was unique (G=1), score = matchedCount/total
        // With repeats, score < matchedCount/total due to 1/G weighting
        return presenceNumerator / static_cast<double>(totalUniqueReadSeeds);
    }

    // LogLogCosine Score: Σ(log(1+log(1+readCount)) × genomeCount) / (mag × genomeMag)
    // Double-log dampens coverage effects even more than single-log
    // IMPROVED: Add coverage penalty to prevent catastrophic failures
    double getLogLogCosineScore(double logLogReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (logLogReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double rawScore = logLogCosineNumerator / (logLogReadMagnitude * genomeMagnitude);
        
        // Coverage penalty using logLog magnitude
        double logLogMagSquared = logLogReadMagnitude * logLogReadMagnitude;
        double coverageFraction = matchedLogLogReadMagnitudeSquared / logLogMagSquared;
        double coveragePenalty = std::sqrt(std::clamp(coverageFraction, 0.0, 1.0));
        
        return std::clamp(rawScore * coveragePenalty, 0.0, 1.0);
    }

    // Concordance Score: LogRaw weighted by fragment pair concordance
    // For each fragment, penalizes when R1/R2 seed match ratio is far from 1
    // Soft penalty: sqrt(ratio) where ratio = min/max, so ratio=0 gives weight=0.5 not 0
    // NOTE: Full concordance scoring is O(fragments) per node - too slow for large datasets
    // Returning logRaw as approximation until we implement incremental concordance tracking
    double getConcordanceScore(const PlacementGlobalState& state, double logReadMagnitude) const {
        // Concordance scoring disabled - too expensive (O(fragments) per node)
        // Just return logRaw score as approximation
        return logRawNumerator / std::max(logReadMagnitude, 1e-9);
    }

    // Compute child metrics from parent metrics + seed changes
    static void computeChildMetrics(
        NodeMetrics& childMetrics,
        std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges,
        const PlacementGlobalState& state);
};

// Consolidated score tracking and results for lite placement
struct PlacementResult {
  // ========================================
  // OLD METRICS (from working version 6cfb6cf)
  // ========================================
  
  // Raw Score results
  double bestRawScore = 0.0;
  uint32_t bestRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedRawNodeIndices;

  // Cosine Score results (standard cosine similarity)
  double bestCosineScore = 0.0;
  uint32_t bestCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCosineNodeIndices;

  // OLD LogCosine Score (the working formula - log on reads only)
  double bestLogCosineOldScore = 0.0;
  uint32_t bestLogCosineOldNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogCosineOldNodeIndices;

  // Weighted Containment Score results
  double bestWeightedContainmentScore = 0.0;
  uint32_t bestWeightedContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedWeightedContainmentNodeIndices;

  // Cap Cosine Score results
  double bestCapCosineScore = 0.0;
  uint32_t bestCapCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCapCosineNodeIndices;

  // Cap Log Cosine Score results
  double bestCapLogCosineScore = 0.0;
  uint32_t bestCapLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCapLogCosineNodeIndices;

  // ========================================
  // NEW METRICS (current version)
  // ========================================

  // LogRAW Score results (log-scaled, coverage-robust)
  double bestLogRawScore = 0.0;
  uint32_t bestLogRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogRawNodeIndices;

  // LogCosine Score results (NEW formula with log on both sides + coverage penalty)
  double bestLogCosineScore = 0.0;
  uint32_t bestLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogCosineNodeIndices;

  // LogLogRAW Score results (double-log-scaled)
  double bestLogLogRawScore = 0.0;
  uint32_t bestLogLogRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogLogRawNodeIndices;

  // LogLogCosine Score results (double-log-scaled cosine)
  double bestLogLogCosineScore = 0.0;
  uint32_t bestLogLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogLogCosineNodeIndices;

  // Containment Score results (fraction of read seed mass in genome)
  double bestContainmentScore = 0.0;
  uint32_t bestContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedContainmentNodeIndices;

  // Concordance Score results (LogRaw weighted by fragment pair balance)
  double bestConcordanceScore = 0.0;
  uint32_t bestConcordanceNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedConcordanceNodeIndices;

  // Presence Score results (low-coverage robust: binary presence weighted by 1/G)
  double bestPresenceScore = 0.0;
  uint32_t bestPresenceNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedPresenceNodeIndices;

  // Alignment-based refinement results
  // After k-mer scoring, top candidates are re-scored by full alignment
  double bestRefinedScore = 0.0;          // Best alignment score from refinement
  uint32_t bestRefinedNodeIndex = UINT32_MAX;
  bool refinementWasRun = false;          // True if refinement was executed

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  
  // Update methods for OLD metrics
  void updateRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogCosineOldScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateWeightedContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCapCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCapLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  
  // Update methods for NEW metrics (take nodeIndex, resolve to string at end)
  void updateLogRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogLogRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateConcordanceScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updatePresenceScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  
  // Lazy resolution: convert winning node indices to string IDs (called ONCE at end)
  void resolveNodeIds(panmapUtils::LiteTree* liteTree);
  
  // Resolved string IDs for OLD metrics
  std::string bestRawNodeId;
  std::string bestCosineNodeId;
  std::string bestLogCosineOldNodeId;
  std::string bestWeightedContainmentNodeId;
  std::string bestCapCosineNodeId;
  std::string bestCapLogCosineNodeId;
  
  // Resolved string IDs for NEW metrics (only populated after resolveNodeIds() is called)
  std::string bestLogRawNodeId;
  std::string bestLogCosineNodeId;
  std::string bestLogLogRawNodeId;
  std::string bestLogLogCosineNodeId;
  std::string bestContainmentNodeId;
  std::string bestConcordanceNodeId;
  std::string bestPresenceNodeId;
  std::string bestRefinedNodeId;
  
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
               uint32_t expectedFragmentSize = 400,
               uint32_t fragmentSizeTolerance = 150,
               bool refineEnabled = false,
               double refineTopPct = 0.01,
               int refineMaxTopN = 150,
               int refineNeighborRadius = 2,
               int refineMaxNeighborN = 150);

} // namespace placement
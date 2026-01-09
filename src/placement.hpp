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
struct PlacementGlobalState {
    // Seed hash frequency map from reads (computed once, used for all nodes)
    // OPTIMIZATION: Use IdentityHash since seed hashes are already well-distributed
    absl::flat_hash_map<size_t, int64_t, IdentityHash> seedFreqInReads;
    
    // Log-scaled read counts: log(1 + readCount) for each seed
    // Pre-computed once to avoid repeated log() calls during traversal
    absl::flat_hash_map<size_t, double, IdentityHash> logReadCounts;
    
    // Coverage-weighted read counts: weight based on expected coverage
    // Penalizes singletons (likely errors) and over-represented seeds (amplicon bias)
    absl::flat_hash_map<size_t, double, IdentityHash> coverageWeightedCounts;
    
    // Capped read counts: cap at 100 to limit influence of high-frequency seeds
    // Seeds with freq > 100 are treated as having freq = 100
    absl::flat_hash_map<size_t, int64_t, IdentityHash> cappedReadCounts;
    
    // Log-capped read counts: log(min(readCount, 100)) for log-scaled capped metrics
    absl::flat_hash_map<size_t, double, IdentityHash> cappedLogReadCounts;
    
    // Sigmoid-scaled read counts: adaptive scaling centered on median seed frequency
    // Uses sigmoid: scaledCount = 2*median / (1 + exp(-k * (count/median - 1)))
    // This compresses high-frequency seeds (amplicon bias) while preserving sensitivity at median
    absl::flat_hash_map<size_t, double, IdentityHash> sigmoidReadCounts;
    
    int64_t totalReadSeedFrequency = 0;  // Sum of all read seed frequencies
    size_t readUniqueSeedCount = 0;       // Number of unique seeds in reads
    double readMagnitude = 0.0;           // Precomputed read magnitude for cosine similarity
    double logReadMagnitude = 0.0;        // Precomputed log-scaled read magnitude: sqrt(Σlog(1+r)²)
    double covWeightedMagnitude = 0.0;    // Precomputed coverage-weighted magnitude: sqrt(Σw²)
    double cappedReadMagnitude = 0.0;     // Precomputed capped magnitude: sqrt(Σmin(r,100)²)
    double cappedLogReadMagnitude = 0.0;  // Precomputed capped log magnitude: sqrt(Σlog(1+min(r,100))²)
    double sigmoidReadMagnitude = 0.0;   // Precomputed sigmoid-scaled magnitude: sqrt(Σsig(r)²)
    double medianReadSeedFrequency = 0.0; // Median seed frequency in reads (midpoint for sigmoid)
    double expectedSeedCoverage = 0.0;    // Expected coverage: totalReadSeedFreq / medianGenomeSeedCount
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
    
    // ========================================
    // Reference genome statistics (computed dynamically from index)
    // Used for normalizing metrics to account for variable genome quality
    // ========================================
    double medianGenomeUniqueSeedCount = 0.0;  // Median seed count across all genomes
    double maxGenomeUniqueSeedCount = 0.0;     // Max seed count (represents "complete" genome)
    double meanGenomeUniqueSeedCount = 0.0;    // Mean seed count
    double stdGenomeUniqueSeedCount = 0.0;     // Std dev of seed counts
    double incompletenessThreshold = 0.0;      // Threshold for penalizing incomplete genomes (mean - 2*std)
    
    // ========================================
    // IDF (Inverse Document Frequency) weights for seeds
    // Seeds appearing in fewer genomes are more discriminative
    // ========================================
    uint32_t totalLeafGenomes = 0;                                    // Total leaf genomes (N for IDF)
    absl::flat_hash_map<uint64_t, double, IdentityHash> seedIdfWeights;  // seed hash -> IDF weight
};

// Delta-based metrics computation (efficient forward-only traversal)
struct NodeMetrics {
    // Read-genome interaction metrics (computed incrementally)
    int64_t jaccardNumerator = 0;          // Presence-based intersection count (same as presenceIntersectionCount)
    int64_t weightedJaccardNumerator = 0;  // Σmin(readCount, genomeCount)
    double cosineNumerator = 0.0;          // Σ(readCount × genomeCount) for seeds in both
    size_t presenceIntersectionCount = 0;  // Count of seeds in both read and genome
    double rawMatchScore = 0.0;            // Σ(readCount / genomeCount) for seeds in intersection
    double logRawNumerator = 0.0;          // Σ(log(1+readCount) / genomeCount) for LogRAW metric
    double logCosineNumerator = 0.0;       // Σ(log(1+readCount) × genomeCount) for LogCosine metric
    double idfCosineNumerator = 0.0;       // Σ(readCount × genomeCount × IDF(seed)) for IDF-weighted cosine
    double covCosineNumerator = 0.0;       // Σ(covWeight × genomeCount) for coverage-weighted cosine
    double capCosineNumerator = 0.0;       // Σ(min(readCount,100) × genomeCount) for capped cosine
    double capLogCosineNumerator = 0.0;    // Σ(log(1+min(readCount,100)) × genomeCount) for capped log cosine
    double sigCosineNumerator = 0.0;       // Σ(sigmoidScaled(readCount) × genomeCount) for sigmoid cosine
    
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

    // Containment Index: |reads ∩ genome| / |reads|
    // Measures what fraction of read seeds are present in the genome
    // Unlike Jaccard, extra genome sequence doesn't hurt the score
    double getContainmentScore(size_t readUniqueSeedCount) const {
        return (readUniqueSeedCount > 0) ?
            static_cast<double>(presenceIntersectionCount) / readUniqueSeedCount : 0.0;
    }

    // Weighted Containment: Σ(readFreq for matching seeds) / Σ(readFreq for all seeds)
    // Uses cosineNumerator which equals Σ(readFreq * genomeFreq) ≈ Σ(readFreq) when genomeFreq=1
    double getWeightedContainmentScore(int64_t totalReadSeedFrequency) const {
        // cosineNumerator = Σ(readFreq * genomeFreq) for matching seeds
        // For most seeds genomeFreq=1, so this approximates Σ(readFreq for matching seeds)
        return (totalReadSeedFrequency > 0) ?
            cosineNumerator / static_cast<double>(totalReadSeedFrequency) : 0.0;
    }

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

    // AdjCosine: N-adjusted Cosine that penalizes genomes with significantly low seed counts
    // Only penalizes genomes with seed counts below (mean - 2*stddev) threshold
    // penalty = 1.0 if seedCount >= threshold, else seedCount / threshold
    // This catches truly incomplete genomes (with N's) without penalizing normal variation
    double getAdjCosineScore(double readMagnitude, double threshold) const {
        double baseCosine = getCosineScore(readMagnitude);
        if (threshold <= 0.0) return baseCosine;
        
        // Only penalize if significantly below threshold
        if (genomeUniqueSeedCount >= threshold) {
            return baseCosine;  // No penalty for normal genomes
        }
        
        // Penalize incomplete genomes proportionally
        double penalty = static_cast<double>(genomeUniqueSeedCount) / threshold;
        return baseCosine * penalty;
    }

    // AdjRaw: N-adjusted RAW score that penalizes significantly incomplete genomes
    double getAdjRawScore(double logReadMagnitude, double threshold) const {
        double baseLogRaw = getLogRawScore(logReadMagnitude);
        if (threshold <= 0.0) return baseLogRaw;
        
        // Only penalize if significantly below threshold
        if (genomeUniqueSeedCount >= threshold) {
            return baseLogRaw;  // No penalty for normal genomes
        }
        
        double penalty = static_cast<double>(genomeUniqueSeedCount) / threshold;
        return baseLogRaw * penalty;
    }

    // IDF-weighted Cosine: Upweights rare seeds that appear in fewer genomes
    // IDF(seed) = log(N / df) where N = total genomes, df = genomes containing seed
    // Score = Σ(readCount × genomeCount × IDF) / (||read|| × ||genome||)
    // Uses same normalization as regular cosine for comparability
    double getIdfCosineScore(double readMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (readMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = idfCosineNumerator / (readMagnitude * genomeMagnitude);
        // Don't clamp - IDF-weighted scores can exceed 1.0 for rare seeds
        return std::max(0.0, score);
    }

    // CovCosine: Coverage-normalized cosine similarity
    // Weights each seed based on how its observed count compares to expected coverage:
    // - Singletons (n=1) are down-weighted when expected coverage > 1 (likely errors)
    // - Over-represented seeds (n >> E) are down-weighted (amplicon bias)
    // - Seeds near expected coverage get full weight
    // When expected coverage <= 1, all seeds get weight 1.0 (can't distinguish signal)
    double getCovCosineScore(double covWeightedMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (covWeightedMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = covCosineNumerator / (covWeightedMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }

    // CapCosine: Capped cosine similarity
    // Caps read counts at 100 to limit influence of high-frequency seeds
    // Score = Σ(min(readCount,100) × genomeCount) / (cappedMag × genomeMag)
    double getCapCosineScore(double cappedReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (cappedReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = capCosineNumerator / (cappedReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }

    // CapLogCosine: Capped log-scaled cosine similarity
    // Uses log(1 + min(readCount, 100)) to cap contribution at log(101)
    // Score = Σ(log(1+min(r,100)) × g) / (cappedLogMag × genomeMag)
    double getCapLogCosineScore(double cappedLogReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (cappedLogReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = capLogCosineNumerator / (cappedLogReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
    }

    // SigCosine: Sigmoid-scaled cosine similarity
    // Uses adaptive sigmoid centered on median seed frequency to handle amplicon bias
    // Score = Σ(sigmoid(r) × g) / (sigmoidMag × genomeMag)
    double getSigCosineScore(double sigmoidReadMagnitude) const {
        double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
        if (sigmoidReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
        double score = sigCosineNumerator / (sigmoidReadMagnitude * genomeMagnitude);
        return std::clamp(score, 0.0, 1.0);
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

  // Containment Index results
  double bestContainmentScore = 0.0;
  uint32_t bestContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedContainmentNodeIndices;

  // Weighted Containment Index results  
  double bestWeightedContainmentScore = 0.0;
  uint32_t bestWeightedContainmentNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedWeightedContainmentNodeIndices;

  // LogRAW Score results (log-scaled, coverage-robust)
  double bestLogRawScore = 0.0;
  uint32_t bestLogRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogRawNodeIndices;

  // LogCosine Score results (log-scaled cosine similarity)
  double bestLogCosineScore = 0.0;
  uint32_t bestLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedLogCosineNodeIndices;

  // AdjCosine Score results (N-adjusted cosine - penalizes incomplete genomes)
  double bestAdjCosineScore = 0.0;
  uint32_t bestAdjCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedAdjCosineNodeIndices;

  // AdjRaw Score results (N-adjusted RAW - penalizes incomplete genomes)
  double bestAdjRawScore = 0.0;
  uint32_t bestAdjRawNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedAdjRawNodeIndices;

  // IDF-Cosine Score results (IDF-weighted cosine - upweights rare seeds)
  double bestIdfCosineScore = 0.0;
  uint32_t bestIdfCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedIdfCosineNodeIndices;

  // CovCosine Score results (coverage-normalized cosine - penalizes singletons & over-rep seeds)
  double bestCovCosineScore = 0.0;
  uint32_t bestCovCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCovCosineNodeIndices;

  // CapCosine Score results (capped cosine - caps read freqs at 100)
  double bestCapCosineScore = 0.0;
  uint32_t bestCapCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCapCosineNodeIndices;

  // CapLogCosine Score results (capped log cosine - caps read freqs at 100 before log)
  double bestCapLogCosineScore = 0.0;
  uint32_t bestCapLogCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedCapLogCosineNodeIndices;

  // SigCosine Score results (sigmoid-scaled cosine - adaptive to median seed freq)
  double bestSigCosineScore = 0.0;
  uint32_t bestSigCosineNodeIndex = UINT32_MAX;
  std::vector<uint32_t> tiedSigCosineNodeIndices;

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  
  // Update methods (take nodeIndex, resolve to string at end)
  void updateHitsScore(uint32_t nodeIndex, int64_t hits, panmapUtils::LiteNode* node = nullptr);
  void updateRawSeedMatchScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateJaccardScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateWeightedJaccardScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateJaccardPresenceScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateWeightedContainmentScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateAdjCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateAdjRawScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateIdfCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCovCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCapCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateCapLogCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  void updateSigCosineScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node = nullptr);
  
  // Lazy resolution: convert winning node indices to string IDs (called ONCE at end)
  void resolveNodeIds(panmapUtils::LiteTree* liteTree);
  
  // Resolved string IDs (only populated after resolveNodeIds() is called)
  std::string bestJaccardNodeId;
  std::string bestCosineNodeId;
  std::string bestWeightedJaccardNodeId;
  std::string bestJaccardPresenceNodeId;
  std::string bestRawSeedMatchNodeId;
  std::string bestContainmentNodeId;
  std::string bestWeightedContainmentNodeId;
  std::string bestLogRawNodeId;
  std::string bestLogCosineNodeId;
  std::string bestAdjCosineNodeId;
  std::string bestAdjRawNodeId;
  std::string bestIdfCosineNodeId;
  std::string bestCovCosineNodeId;
  std::string bestCapCosineNodeId;
  std::string bestCapLogCosineNodeId;
  std::string bestSigCosineNodeId;
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
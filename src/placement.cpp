#include "panman.hpp"
#include "placement.hpp"
#include "seeding.hpp"
#include "index_lite.capnp.h"
#include "logging.hpp"
#include "mm_align.h"


#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/enumerable_thread_specific.h>
#include <zlib.h>
#include "kseq.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// KSEQ_INIT already defined in seeding.hpp - don't redefine

namespace {

// Compute the canonical hash of a homopolymer k-mer (all same base)
// Returns the min of forward and reverse complement hash
inline size_t computeHomopolymerHash(char base, int k) {
    size_t baseVal = seeding::chash(base);
    size_t compVal;
    switch (base) {
        case 'A': case 'a': compVal = seeding::chash('T'); break;
        case 'C': case 'c': compVal = seeding::chash('G'); break;
        case 'G': case 'g': compVal = seeding::chash('C'); break;
        case 'T': case 't': compVal = seeding::chash('A'); break;
        default: return 0;
    }
    
    // Forward hash: XOR of rol(val, k-1), rol(val, k-2), ..., rol(val, 0)
    size_t fHash = 0;
    for (int i = 0; i < k; i++) {
        fHash ^= seeding::rol(baseVal, k - i - 1);
    }
    
    // Reverse complement hash: same formula but with complement base
    size_t rHash = 0;
    for (int i = 0; i < k; i++) {
        rHash ^= seeding::rol(compVal, k - i - 1);
    }
    
    // Return canonical (minimum) hash
    return std::min(fHash, rHash);
}

// Get all 4 canonical homopolymer hashes for a given k
inline std::array<size_t, 4> getHomopolymerHashes(int k) {
    return {
        computeHomopolymerHash('A', k),
        computeHomopolymerHash('C', k),
        computeHomopolymerHash('G', k),
        computeHomopolymerHash('T', k)
    };
}

// Compute average Phred quality score over a k-mer span
// Quality string uses ASCII encoding: Q = char - 33 (Phred+33)
inline double avgPhredQuality(const std::string& qual, int64_t startPos, int k) {
    if (qual.empty() || startPos < 0 || startPos + k > static_cast<int64_t>(qual.size())) {
        return 0.0;
    }
    int64_t sum = 0;
    for (int i = 0; i < k; ++i) {
        sum += static_cast<int>(qual[startPos + i]) - 33;
    }
    return static_cast<double>(sum) / k;
}

void extractReadSequences(const std::string& readPath1, const std::string& readPath2, std::vector<std::string>& readSequences) {
  FILE *fp;
  kseq_t *seq;
  fp = fopen(readPath1.c_str(), "r");
  if(!fp){
    output::error("File {} not found", readPath1);
    exit(1);
  }
  seq = kseq_init(fileno(fp));
  int line;
  while ((line = kseq_read(seq)) >= 0) {
    readSequences.push_back(seq->seq.s);
  }
  kseq_destroy(seq);
  fclose(fp);

  if (readPath2.size() > 0) {
    fp = fopen(readPath2.c_str(), "r");
    if(!fp){
      output::error("File {} not found", readPath2);
      exit(1);
    }
    seq = kseq_init(fileno(fp));

    line = 0;
    int forwardReads = readSequences.size();
    while ((line = kseq_read(seq)) >= 0) {
      readSequences.push_back(seq->seq.s);
    }
    kseq_destroy(seq);
    fclose(fp);

    if (readSequences.size() != forwardReads*2){
      output::error("File {} does not contain the same number of reads as {}", readPath2, readPath1);
      exit(1);
    }
    
    //Shuffle reads together, so that pairs are next to each other
    seeding::perfect_shuffle(readSequences);
  }
}

// Extract full FASTQ data including sequences, qualities, and names
void extractFullFastqData(const std::string& readPath1, const std::string& readPath2,
                          std::vector<std::string>& readSequences,
                          std::vector<std::string>& readQuals,
                          std::vector<std::string>& readNames) {
  FILE *fp;
  kseq_t *seq;
  fp = fopen(readPath1.c_str(), "r");
  if(!fp){
    output::error("File {} not found", readPath1);
    exit(1);
  }
  seq = kseq_init(fileno(fp));
  int line;
  while ((line = kseq_read(seq)) >= 0) {
    readSequences.push_back(seq->seq.s);
    readNames.push_back(seq->name.s);
    // Quality may be empty for FASTA files
    readQuals.push_back(seq->qual.l > 0 ? seq->qual.s : std::string(seq->seq.l, 'I'));
  }
  kseq_destroy(seq);
  fclose(fp);

  if (readPath2.size() > 0) {
    fp = fopen(readPath2.c_str(), "r");
    if(!fp){
      output::error("File {} not found", readPath2);
      exit(1);
    }
    seq = kseq_init(fileno(fp));

    line = 0;
    int forwardReads = readSequences.size();
    while ((line = kseq_read(seq)) >= 0) {
      readSequences.push_back(seq->seq.s);
      readNames.push_back(seq->name.s);
      readQuals.push_back(seq->qual.l > 0 ? seq->qual.s : std::string(seq->seq.l, 'I'));
    }
    kseq_destroy(seq);
    fclose(fp);

    if (readSequences.size() != static_cast<size_t>(forwardReads*2)){
      output::error("File {} does not contain the same number of reads as {}", readPath2, readPath1);
      exit(1);
    }
    
    // Shuffle reads together, so that pairs are next to each other
    seeding::perfect_shuffle(readSequences);
    seeding::perfect_shuffle(readQuals);
    seeding::perfect_shuffle(readNames);
  }
}

// Build paired-end fragment tracking data
// NOTE: Full fragment R1/R2 seed tracking is disabled (was O(fragments) per node for concordance)
// Now just tracks basic paired-end info for potential future use
void buildPairedEndFragmentData(
    placement::PlacementGlobalState& state,
    const std::vector<std::string>& allReadSequences,
    const placement::TraversalParams& params) {
    
    // Check if paired-end (must have even number of reads, >= 2)
    if (allReadSequences.size() < 2 || allReadSequences.size() % 2 != 0) {
        state.isPairedEnd = false;
        state.numFragments = allReadSequences.size();
        logging::info("Single-end mode: {} reads", allReadSequences.size());
        return;
    }
    
    state.isPairedEnd = true;
    state.numFragments = allReadSequences.size() / 2;
    
    // Skip expensive fragment seed building - concordance scoring is disabled
    logging::info("Paired-end mode: {} fragments (fragment seed tracking disabled)", 
                  state.numFragments);
}

} // anonymous namespace

using namespace logging;


void placement::NodeMetrics::computeChildMetrics(
    placement::NodeMetrics& childMetrics,
    std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges,
    const placement::PlacementGlobalState& state) {

    const size_t numChanges = seedChanges.size();
    
    // Early exit for empty changes
    if (numChanges == 0) [[unlikely]] return;
    constexpr size_t PREFETCH_DISTANCE = 24; // Increased to 24 for deeper pipeline (L3 latency ~40 cycles)
    constexpr size_t BATCH_SIZE = 64; // Process in cache-line aligned batches
    
    // OPTIMIZATION: Pre-fetch both seedChanges AND hash table buckets aggressively
    const size_t initialPrefetchCount = std::min(numChanges, PREFETCH_DISTANCE);
    for (size_t i = 0; i < initialPrefetchCount; ++i) {
        const auto& [seedHash, _, __] = seedChanges[i];
        __builtin_prefetch(&seedChanges[i], 0, 0);
        state.seedFreqInReads.prefetch(seedHash);
    }
    
    // Process in batches for better instruction pipelining and cache utilization
    for (size_t batch_start = 0; batch_start < numChanges; batch_start += BATCH_SIZE) {
        const size_t batch_end = std::min(batch_start + BATCH_SIZE, numChanges);
        
        for (size_t i = batch_start; i < batch_end; ++i) {
            // Prefetch far ahead - both data and hash table
            if (i + PREFETCH_DISTANCE < numChanges) [[likely]] {
                const auto& [nextHash, _, __] = seedChanges[i + PREFETCH_DISTANCE];
                __builtin_prefetch(&seedChanges[i + PREFETCH_DISTANCE], 0, 0);
                state.seedFreqInReads.prefetch(nextHash);
            }
        
            const auto& [seedHash, parentCount, childCount] = seedChanges[i];
            
            // ========================================
            // GENOME-ONLY METRICS (computed for ALL seed changes)
            // These don't require read lookup - just parent/child counts
            // ========================================
            
            // Genome total seed frequency delta
            const int64_t freqDelta = childCount - parentCount;
            childMetrics.genomeTotalSeedFrequency += freqDelta;
            
            // Genome magnitude squared delta: child² - parent²
            const double magDelta = static_cast<double>(childCount * childCount) 
                                  - static_cast<double>(parentCount * parentCount);
            childMetrics.genomeMagnitudeSquared += magDelta;
            
            // Log-dampened genome magnitude: Σ(log(1+G)²) for IDF-weighted cosine denominator
            const double oldLogG = (parentCount > 0) ? std::log1p(static_cast<double>(parentCount)) : 0.0;
            const double newLogG = (childCount > 0) ? std::log1p(static_cast<double>(childCount)) : 0.0;
            childMetrics.logGenomeMagnitudeSquared += (newLogG * newLogG) - (oldLogG * oldLogG);
            
            // Genome unique seed count delta (branchless)
            const int64_t wasPresent = (parentCount > 0) ? 1 : 0;
            const int64_t isPresent = (childCount > 0) ? 1 : 0;
            childMetrics.genomeUniqueSeedCount += (isPresent - wasPresent);
            
            // ========================================
            // READ-GENOME INTERACTION METRICS (only for seeds in reads)
            // SIMPLIFIED: Only compute LogRaw and LogCosine
            // ========================================
            
            // Fast path: skip unchanged seeds (23% of all seeds per profiler)
            if (freqDelta == 0) [[likely]] continue;
            
            // Single hash lookup: find() returns iterator, check if valid
            auto readIt = state.seedFreqInReads.find(seedHash);
            if (readIt == state.seedFreqInReads.end()) [[likely]] continue;
            
            // Hot path: seed changed AND present in reads (29% of all seeds)
            // Get read count and pre-computed values
            const int64_t readCount = readIt->second;
            const double readCountD = static_cast<double>(readCount);
            
            auto logIt = state.logReadCounts.find(seedHash);
            if (logIt != state.logReadCounts.end()) {
                const double logReadCount = logIt->second;
                
                // ========================================
                // OLD METRICS (from working version 6cfb6cf)
                // ========================================
                
                // Raw numerator: Σ(readCount / genomeCount)
                const double oldRawContrib = (parentCount > 0) ? (readCountD / parentCount) : 0.0;
                const double newRawContrib = (childCount > 0) ? (readCountD / childCount) : 0.0;
                childMetrics.rawNumerator += (newRawContrib - oldRawContrib);
                
                // Cosine numerator: Σ(readCount × genomeCount)
                childMetrics.cosineNumerator += readCountD * freqDelta;
                
                // Weighted Jaccard: Σmin(readCount, genomeCount)
                const int64_t oldMin = std::min(readCount, parentCount);
                const int64_t newMin = std::min(readCount, childCount);
                childMetrics.weightedJaccardNumerator += (newMin - oldMin);
                
                // Capped metrics: cap read count at 100
                const int64_t cappedReadCount = std::min(readCount, int64_t(100));
                const double cappedReadCountD = static_cast<double>(cappedReadCount);
                const double cappedLogReadCount = std::log1p(cappedReadCountD);
                childMetrics.capCosineNumerator += cappedReadCountD * freqDelta;
                childMetrics.capLogCosineNumerator += cappedLogReadCount * freqDelta;
                
                // Presence intersection count (binary: seed in both read and genome)
                const bool wasInGenome = (parentCount > 0);
                const bool isInGenome = (childCount > 0);
                if (!wasInGenome && isInGenome) {
                    childMetrics.presenceIntersectionCount += 1;
                    // Track matched read magnitude for old containment
                    childMetrics.matchedReadMagnitudeSquared += readCountD * readCountD;
                } else if (wasInGenome && !isInGenome) {
                    childMetrics.presenceIntersectionCount -= 1;
                    childMetrics.matchedReadMagnitudeSquared -= readCountD * readCountD;
                }
                
                // ========================================
                // NEW METRICS (current version)
                // ========================================
                
                // Update logRaw numerator: Σ(log(1+readCount) / genomeCount) for seeds in intersection
                const double oldLogContrib = (parentCount > 0) ? (logReadCount / parentCount) : 0.0;
                const double newLogContrib = (childCount > 0) ? (logReadCount / childCount) : 0.0;
                childMetrics.logRawNumerator += (newLogContrib - oldLogContrib);
                
                // Update logCosine numerator: Σ(log(1+readCount) × genomeCount) delta (OLD FORMULA!)
                childMetrics.logCosineNumerator += logReadCount * freqDelta;
                
                // IDF-WEIGHTED COSINE: Use log(1+genomeCount) instead of raw genomeCount
                // This prevents repeat amplification: 100 copies → log(101) ≈ 4.6, not 100
                const double oldLogGenome = (parentCount > 0) ? std::log1p(static_cast<double>(parentCount)) : 0.0;
                const double newLogGenome = (childCount > 0) ? std::log1p(static_cast<double>(childCount)) : 0.0;
                childMetrics.logCosineIdfNumerator += logReadCount * (newLogGenome - oldLogGenome);
                
                // Update logLogRaw and logLogCosine numerators
                auto logLogIt = state.logLogReadCounts.find(seedHash);
                if (logLogIt != state.logLogReadCounts.end()) {
                    const double logLogReadCount = logLogIt->second;
                    
                    // LogLogRaw: Σ(log(1+log(1+readCount)) / genomeCount)
                    const double oldLogLogContrib = (parentCount > 0) ? (logLogReadCount / parentCount) : 0.0;
                    const double newLogLogContrib = (childCount > 0) ? (logLogReadCount / childCount) : 0.0;
                    childMetrics.logLogRawNumerator += (newLogLogContrib - oldLogLogContrib);
                    
                    // LogLogCosine: Σ(log(1+log(1+readCount)) × genomeCount)
                    childMetrics.logLogCosineNumerator += logLogReadCount * freqDelta;
                    
                    // Track logLog matched magnitude for coverage penalty
                    const double logLogReadCountSquared = logLogReadCount * logLogReadCount;
                    if (!wasInGenome && isInGenome) {
                        childMetrics.matchedLogLogReadMagnitudeSquared += logLogReadCountSquared;
                    } else if (wasInGenome && !isInGenome) {
                        childMetrics.matchedLogLogReadMagnitudeSquared -= logLogReadCountSquared;
                    }
                }
                
                // COVERAGE TRACKING: Track matched read seed magnitude squared
                // This seed was/is in reads. Track if it enters/leaves the genome.
                const double logReadCountSquared = logReadCount * logReadCount;
                
                if (!wasInGenome && isInGenome) {
                    // Seed newly appears in genome - add to matched magnitude
                    childMetrics.matchedLogReadMagnitudeSquared += logReadCountSquared;
                    // PRESENCE SCORE: Add 1/G for this newly matched seed
                    childMetrics.presenceNumerator += 1.0 / static_cast<double>(childCount);
                    childMetrics.matchedUniqueReadSeeds += 1;
                } else if (wasInGenome && !isInGenome) {
                    // Seed disappears from genome - remove from matched magnitude
                    childMetrics.matchedLogReadMagnitudeSquared -= logReadCountSquared;
                    // PRESENCE SCORE: Remove contribution (was 1/parentCount)
                    childMetrics.presenceNumerator -= 1.0 / static_cast<double>(parentCount);
                    childMetrics.matchedUniqueReadSeeds -= 1;
                } else if (wasInGenome && isInGenome && parentCount != childCount) {
                    // Seed still present but count changed - update 1/G contribution
                    childMetrics.presenceNumerator -= 1.0 / static_cast<double>(parentCount);
                    childMetrics.presenceNumerator += 1.0 / static_cast<double>(childCount);
                }
                // If seed stays in genome (count changes but still > 0), magnitude unchanged
            }
        }
    }
}


namespace placement {

using ::panmanUtils::Node;
using ::panmanUtils::Tree;

// Forward declarations
class PlacementResult;

// Macro to generate score update functions - reduces repetitive code
// Each function: stores score on node, updates best score/index, tracks ties
#define DEFINE_UPDATE_SCORE_FUNC(FuncName, nodeScoreField, bestScore, bestNodeIndex, tiedIndices) \
void PlacementResult::FuncName(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node) { \
    if (node) { \
        node->nodeScoreField = static_cast<float>(score); \
    } \
    double tolerance = std::max(bestScore * 0.0001, 1e-9); \
    if (score > bestScore + tolerance) { \
        bestScore = score; \
        bestNodeIndex = nodeIndex; \
        tiedIndices.clear(); \
        tiedIndices.push_back(nodeIndex); \
    } else if (score >= bestScore - tolerance && score > 0) { \
        if (tiedIndices.empty() || tiedIndices.back() != bestNodeIndex) { \
            tiedIndices.push_back(bestNodeIndex); \
        } \
        if (nodeIndex != bestNodeIndex) { \
            tiedIndices.push_back(nodeIndex); \
        } \
    } \
}

// Generate score update functions for OLD metrics
DEFINE_UPDATE_SCORE_FUNC(updateRawScore, rawScore, bestRawScore, bestRawNodeIndex, tiedRawNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCosineScore, cosineScore, bestCosineScore, bestCosineNodeIndex, tiedCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogCosineOldScore, logCosineOldScore, bestLogCosineOldScore, bestLogCosineOldNodeIndex, tiedLogCosineOldNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateWeightedContainmentScore, weightedContainmentScore, bestWeightedContainmentScore, bestWeightedContainmentNodeIndex, tiedWeightedContainmentNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCapCosineScore, capCosineScore, bestCapCosineScore, bestCapCosineNodeIndex, tiedCapCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCapLogCosineScore, capLogCosineScore, bestCapLogCosineScore, bestCapLogCosineNodeIndex, tiedCapLogCosineNodeIndices)

// Generate score update functions for NEW metrics
DEFINE_UPDATE_SCORE_FUNC(updateLogRawScore, logRawScore, bestLogRawScore, bestLogRawNodeIndex, tiedLogRawNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogCosineScore, logCosineScore, bestLogCosineScore, bestLogCosineNodeIndex, tiedLogCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogLogRawScore, logLogRawScore, bestLogLogRawScore, bestLogLogRawNodeIndex, tiedLogLogRawNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogLogCosineScore, logLogCosineScore, bestLogLogCosineScore, bestLogLogCosineNodeIndex, tiedLogLogCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateContainmentScore, containmentScore, bestContainmentScore, bestContainmentNodeIndex, tiedContainmentNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateConcordanceScore, concordanceScore, bestConcordanceScore, bestConcordanceNodeIndex, tiedConcordanceNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updatePresenceScore, presenceScore, bestPresenceScore, bestPresenceNodeIndex, tiedPresenceNodeIndices)

#undef DEFINE_UPDATE_SCORE_FUNC

void PlacementResult::resolveNodeIds(panmapUtils::LiteTree* liteTree) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (!liteTree) return;
    
    // Resolve OLD metric node indices to string IDs
    if (bestRawNodeIndex != UINT32_MAX) {
        bestRawNodeId = liteTree->resolveNodeId(bestRawNodeIndex);
    }
    if (bestCosineNodeIndex != UINT32_MAX) {
        bestCosineNodeId = liteTree->resolveNodeId(bestCosineNodeIndex);
    }
    if (bestLogCosineOldNodeIndex != UINT32_MAX) {
        bestLogCosineOldNodeId = liteTree->resolveNodeId(bestLogCosineOldNodeIndex);
    }
    if (bestWeightedContainmentNodeIndex != UINT32_MAX) {
        bestWeightedContainmentNodeId = liteTree->resolveNodeId(bestWeightedContainmentNodeIndex);
    }
    if (bestCapCosineNodeIndex != UINT32_MAX) {
        bestCapCosineNodeId = liteTree->resolveNodeId(bestCapCosineNodeIndex);
    }
    if (bestCapLogCosineNodeIndex != UINT32_MAX) {
        bestCapLogCosineNodeId = liteTree->resolveNodeId(bestCapLogCosineNodeIndex);
    }
    
    // Resolve NEW metric node indices to string IDs
    if (bestLogRawNodeIndex != UINT32_MAX) {
        bestLogRawNodeId = liteTree->resolveNodeId(bestLogRawNodeIndex);
    }
    if (bestLogCosineNodeIndex != UINT32_MAX) {
        bestLogCosineNodeId = liteTree->resolveNodeId(bestLogCosineNodeIndex);
    }
    if (bestLogLogRawNodeIndex != UINT32_MAX) {
        bestLogLogRawNodeId = liteTree->resolveNodeId(bestLogLogRawNodeIndex);
    }
    if (bestLogLogCosineNodeIndex != UINT32_MAX) {
        bestLogLogCosineNodeId = liteTree->resolveNodeId(bestLogLogCosineNodeIndex);
    }
    if (bestContainmentNodeIndex != UINT32_MAX) {
        bestContainmentNodeId = liteTree->resolveNodeId(bestContainmentNodeIndex);
    }
    if (bestConcordanceNodeIndex != UINT32_MAX) {
        bestConcordanceNodeId = liteTree->resolveNodeId(bestConcordanceNodeIndex);
    }
    if (bestPresenceNodeIndex != UINT32_MAX) {
        bestPresenceNodeId = liteTree->resolveNodeId(bestPresenceNodeIndex);
    }
    if (bestRefinedNodeIndex != UINT32_MAX) {
        bestRefinedNodeId = liteTree->resolveNodeId(bestRefinedNodeIndex);
    }
}


// =============================================================================
// ALIGNMENT-BASED REFINEMENT
// After k-mer scoring, refine top candidates by full minimap2 alignment
// =============================================================================

// Get all nodes within phylogenetic distance `radius` of a given node
// Uses BFS traversal following parent/child edges
std::vector<panmapUtils::LiteNode*> getNodesWithinRadius(
    panmapUtils::LiteNode* startNode, 
    int radius, 
    int maxNodes) {
    
    if (!startNode || radius <= 0 || maxNodes <= 0) {
        return {};
    }
    
    std::vector<panmapUtils::LiteNode*> result;
    absl::flat_hash_set<panmapUtils::LiteNode*> visited;
    
    // BFS with distance tracking: (node, distance)
    std::queue<std::pair<panmapUtils::LiteNode*, int>> bfsQueue;
    bfsQueue.push({startNode, 0});
    visited.insert(startNode);
    
    while (!bfsQueue.empty() && static_cast<int>(result.size()) < maxNodes) {
        auto [node, dist] = bfsQueue.front();
        bfsQueue.pop();
        
        // Don't include start node in results (it's already in top candidates)
        if (node != startNode) {
            result.push_back(node);
        }
        
        // Stop expanding if we've reached max radius
        if (dist >= radius) continue;
        
        // Add parent (if exists and not visited)
        if (node->parent && !visited.count(node->parent)) {
            visited.insert(node->parent);
            bfsQueue.push({node->parent, dist + 1});
        }
        
        // Add children
        for (auto* child : node->children) {
            if (child && !visited.count(child)) {
                visited.insert(child);
                bfsQueue.push({child, dist + 1});
            }
        }
    }
    
    return result;
}

// Get node sequence from full tree
// Wrapper to access getStringFromReference on the full tree
std::string getNodeSequenceForRefinement(panmanUtils::Tree* T, const std::string& nodeId) {
    if (!T) return "";
    return T->getStringFromReference(nodeId, false, true);
}

// Score a candidate node by aligning reads to its genome sequence
// Returns sum of primary alignment scores
int64_t scoreNodeByAlignment(
    panmanUtils::Tree* fullTree,
    const std::string& nodeId,
    const std::vector<std::string>& readSequences,
    int kmerSize) {
    
    if (!fullTree || nodeId.empty() || readSequences.empty()) {
        return 0;
    }
    
    // Get genome sequence for this node
    std::string genomeSeq = getNodeSequenceForRefinement(fullTree, nodeId);
    if (genomeSeq.empty()) {
        return 0;
    }
    
    // Prepare read data for alignment
    int n_reads = static_cast<int>(readSequences.size());
    std::vector<const char*> readPtrs(n_reads);
    std::vector<int> readLens(n_reads);
    
    for (int i = 0; i < n_reads; i++) {
        readPtrs[i] = readSequences[i].c_str();
        readLens[i] = static_cast<int>(readSequences[i].size());
    }
    
    // Call minimap2 scoring function
    int64_t score = score_reads_vs_reference(
        genomeSeq.c_str(),
        n_reads,
        readPtrs.data(),
        readLens.data(),
        kmerSize
    );
    
    return score;
}

// Main refinement function: align reads to top candidates and pick best
// Called after BFS traversal to refine placement by full alignment
void refineTopCandidates(
    panmapUtils::LiteTree* liteTree,
    panmanUtils::Tree* fullTree,
    const std::vector<std::string>& readSequences,
    PlacementResult& result,
    const TraversalParams& params) {
    
    if (!fullTree || !liteTree || readSequences.empty()) {
        logging::warn("Refinement skipped: missing tree data or no reads");
        return;
    }
    
    auto time_refine_start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Collect top candidates from ALL metrics (not just LogLogRaw)
    // This ensures we don't miss the correct node if different metrics rank it differently
    absl::flat_hash_set<uint32_t> candidateSet;  // Avoid duplicates
    
    auto addTopFromMetric = [&](auto scoreGetter, const char* metricName) {
        std::vector<std::pair<double, uint32_t>> scoredNodes;
        for (size_t i = 0; i < liteTree->dfsIndexToNode.size(); i++) {
            panmapUtils::LiteNode* node = liteTree->dfsIndexToNode[i];
            if (node) {
                double score = scoreGetter(node);
                if (score > 0) {
                    scoredNodes.push_back({score, static_cast<uint32_t>(i)});
                }
            }
        }
        if (scoredNodes.empty()) return;
        
        std::sort(scoredNodes.begin(), scoredNodes.end(), std::greater<>());
        
        size_t numTop = std::min(
            static_cast<size_t>(scoredNodes.size() * params.refineTopPct),
            static_cast<size_t>(params.refineMaxTopN)
        );
        numTop = std::max(numTop, size_t(1));
        
        size_t added = 0;
        for (size_t i = 0; i < numTop && i < scoredNodes.size(); i++) {
            if (candidateSet.insert(scoredNodes[i].second).second) {
                added++;
            }
        }
        logging::debug("Refinement: {} added {} unique candidates from top {}", metricName, added, numTop);
    };
    
    // Add top candidates from each metric
    addTopFromMetric([](auto* n) { return n->logRawScore; }, "LogRaw");
    addTopFromMetric([](auto* n) { return n->logCosineScore; }, "LogCosine");
    addTopFromMetric([](auto* n) { return n->logLogRawScore; }, "LogLogRaw");
    addTopFromMetric([](auto* n) { return n->presenceScore; }, "Presence");
    addTopFromMetric([](auto* n) { return n->containmentScore; }, "Containment");
    addTopFromMetric([](auto* n) { return n->cosineScore; }, "Cosine");
    
    if (candidateSet.empty()) {
        logging::warn("Refinement skipped: no nodes with positive scores");
        return;
    }
    
    logging::info("Refinement: {} unique candidates from all metrics ({}%)", 
                 candidateSet.size(), params.refineTopPct * 100);
    
    // Step 2: Align reads to each candidate and collect neighbors
    std::vector<std::pair<int64_t, uint32_t>> alignmentScores;  // (score, nodeIndex)
    absl::flat_hash_set<uint32_t> scoredSet;  // Track what we've scored
    
    for (uint32_t nodeIdx : candidateSet) {
        if (scoredSet.count(nodeIdx)) continue;
        scoredSet.insert(nodeIdx);
        
        panmapUtils::LiteNode* node = liteTree->dfsIndexToNode[nodeIdx];
        if (!node) continue;
        
        // Score this candidate
        int64_t score = scoreNodeByAlignment(fullTree, node->identifier, readSequences, params.k);
        alignmentScores.push_back({score, nodeIdx});
        
        // Step 3: Get neighbors within radius and score them too
        auto neighbors = getNodesWithinRadius(node, params.refineNeighborRadius, params.refineMaxNeighborN);
        
        for (auto* neighbor : neighbors) {
            if (neighbor && !scoredSet.count(neighbor->nodeIndex)) {
                scoredSet.insert(neighbor->nodeIndex);
                
                int64_t neighborScore = scoreNodeByAlignment(
                    fullTree, neighbor->identifier, readSequences, params.k);
                alignmentScores.push_back({neighborScore, neighbor->nodeIndex});
            }
        }
    }
    
    // Step 4: Find best by alignment score, with seed score as tiebreaker
    if (alignmentScores.empty()) {
        logging::warn("Refinement produced no alignment scores");
        return;
    }
    
    // Find best alignment score, breaking ties by seed-based score (logRawScore)
    auto bestIt = std::max_element(alignmentScores.begin(), alignmentScores.end(),
        [&liteTree](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first < b.first;  // Higher alignment score wins
            // Tie-break by seed score (higher is better)
            auto* nodeA = liteTree->dfsIndexToNode[a.second];
            auto* nodeB = liteTree->dfsIndexToNode[b.second];
            float scoreA = nodeA ? nodeA->logRawScore : 0.0f;
            float scoreB = nodeB ? nodeB->logRawScore : 0.0f;
            return scoreA < scoreB;  // Higher seed score wins
        });
    result.bestRefinedScore = static_cast<double>(bestIt->first);
    result.bestRefinedNodeIndex = bestIt->second;
    result.refinementWasRun = true;
    
    auto time_refine_end = std::chrono::high_resolution_clock::now();
    auto duration_refine = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_refine_end - time_refine_start);
    
    logging::info("Refinement complete: {} candidates scored, best alignment score = {} ({}ms)", 
                 alignmentScores.size(), result.bestRefinedScore, duration_refine.count());
}


// Placement helper using BFS traversal with delta-based metrics
void placeLiteHelperBFS(
    std::vector<panmapUtils::LiteNode*>& nodes,
    std::vector<placement::NodeMetrics>& parentMetrics,
    std::vector<absl::flat_hash_map<uint64_t, int64_t>>& parentSeedCounts,
    const placement::PlacementGlobalState& state,
    placement::PlacementResult& result,
    const placement::TraversalParams& params,
    std::atomic<size_t>& nodesProcessed) {

    if (nodes.empty()) {
        return;
    }

    std::vector<panmapUtils::LiteNode*> current_nodes = std::move(nodes);
    std::vector<placement::NodeMetrics> current_metrics = std::move(parentMetrics);
    std::vector<absl::flat_hash_map<uint64_t, int64_t>> current_seeds = std::move(parentSeedCounts);
    
    struct alignas(64) ThreadLocalData {  // Cache-line aligned to prevent false sharing
        std::vector<panmapUtils::LiteNode*> next_nodes;
        std::vector<placement::NodeMetrics> next_metrics;
        std::vector<absl::flat_hash_map<uint64_t, int64_t>> next_seeds;
        
        // Thread-local best scores (lock-free during traversal!)
        placement::PlacementResult local_result;
        
        // Pre-allocate buffers to reduce reallocation
        ThreadLocalData() {
            next_nodes.reserve(1024);
            next_metrics.reserve(1024);
            next_seeds.reserve(1024);
        }
        
        void clear() {
            next_nodes.clear();
            next_metrics.clear();
            next_seeds.clear();
        }
    };

    // Use enumerable_thread_specific to manage thread-local storage automatically
    tbb::enumerable_thread_specific<ThreadLocalData> thread_data;
    
    // Get thread count for optimal grain size
    const size_t numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

    while (!current_nodes.empty()) {
        size_t num_nodes = current_nodes.size();
        
        // Adaptive grain size based on level size and thread count
        // Larger chunks = less overhead, but need enough chunks for load balancing
        size_t grain_size = std::max(size_t(1), num_nodes / (numThreads * 8));
        
        // Pre-allocate empty seed counts for verification disabled case
        static const absl::flat_hash_map<uint64_t, int64_t> empty_seed_counts;

        tbb::parallel_for(tbb::blocked_range<size_t>(0, num_nodes, grain_size), 
            [&](const tbb::blocked_range<size_t>& r) {
            
            auto& tls = thread_data.local();
            size_t batch_counter = 0;

            for (size_t i = r.begin(); i != r.end(); ++i) {
                panmapUtils::LiteNode* node = current_nodes[i];
                const placement::NodeMetrics& p_metrics = current_metrics[i];
                const absl::flat_hash_map<uint64_t, int64_t>& p_seed_counts = params.verify_scores ? current_seeds[i] : empty_seed_counts;

                if (!node) continue;
                
                // Aggressive software prefetch for next 4 nodes
                constexpr size_t PREFETCH_DISTANCE = 4;
                for (size_t j = 1; j <= PREFETCH_DISTANCE && i + j < r.end(); ++j) {
                    __builtin_prefetch(current_nodes[i + j], 0, 1);         // Prefetch node pointer
                    __builtin_prefetch(&current_metrics[i + j], 0, 1);      // Prefetch metrics
                    if (current_nodes[i + j]) {
                        __builtin_prefetch(&current_nodes[i + j]->seedChanges, 0, 1);  // Prefetch seed changes
                    }
                }
                
                batch_counter++;
                uint32_t nodeIndex = node->nodeIndex;  // Use index instead of deserializing string!

                placement::NodeMetrics nodeMetrics = p_metrics;
                placement::NodeMetrics::computeChildMetrics(nodeMetrics, node->seedChanges, state);
                
                absl::flat_hash_map<uint64_t, int64_t> incrementalGenomeCounts;
                if (params.verify_scores) {
                    incrementalGenomeCounts = p_seed_counts;
                    for (const auto& [seedHash, parentCount, childCount] : node->seedChanges) {
                        int64_t delta = childCount - parentCount;
                        if (delta != 0) {
                            incrementalGenomeCounts[seedHash] += delta;
                            if (incrementalGenomeCounts[seedHash] <= 0) {
                                incrementalGenomeCounts.erase(seedHash);
                            }
                        }
                    }
                }
                
                // Skip scoring for the excluded leaf node (leave-one-out validation)
                if (nodeIndex != state.skipNodeIndex) {
                    // ========================================
                    // Compute OLD metrics (from working version 6cfb6cf)
                    // ========================================
                    
                    // RAW Score
                    double rawScore = nodeMetrics.getRawScore(state.readMagnitude);
                    
                    // Cosine Score (standard cosine similarity)
                    double cosineScore = nodeMetrics.getCosineScore(state.readMagnitude);
                    
                    // OLD LogCosine Score (the working formula - log on reads only)
                    double logCosineOldScore = nodeMetrics.getLogCosineOldScore(state.logReadMagnitude);
                    
                    // Weighted Containment Score
                    double weightedContainmentScore = nodeMetrics.getWeightedContainmentScore(state.totalReadSeedFrequency);
                    
                    // Cap Cosine Score
                    double capCosineScore = nodeMetrics.getCapCosineScore(state.cappedReadMagnitude);
                    
                    // Cap Log Cosine Score
                    double capLogCosineScore = nodeMetrics.getCapLogCosineScore(state.cappedLogReadMagnitude);
                    
                    // ========================================
                    // Compute NEW metrics (current version)
                    // ========================================
                    
                    // LogRAW Score: log-scaled, coverage-robust metric
                    double logRawScore = nodeMetrics.getLogRawScore(state.logReadMagnitude);
                    
                    // NEW LogCosine Score: IDF-weighted with coverage penalty
                    double logCosineScore = nodeMetrics.getLogCosineScore(state.logReadMagnitude);
                    
                    // LogLogRAW Score: double-log dampened coverage
                    double logLogRawScore = nodeMetrics.getLogLogRawScore(state.logLogReadMagnitude);
                    
                    // LogLogCosine Score: double-log cosine similarity
                    double logLogCosineScore = nodeMetrics.getLogLogCosineScore(state.logLogReadMagnitude);
                    
                    // Containment Score: fraction of read seeds covered
                    double containmentScore = nodeMetrics.getContainmentScore(state.logReadMagnitude);
                    
                    // Concordance Score: R1/R2 balance weighted logRaw
                    double concordanceScore = nodeMetrics.getConcordanceScore(state, state.logReadMagnitude);
                    
                    // Presence Score: low-coverage robust (binary presence weighted by 1/G)
                    double presenceScore = nodeMetrics.getPresenceScore(state.readUniqueSeedCount);
                    
                    // Update thread-local results for OLD metrics
                    tls.local_result.updateRawScore(nodeIndex, rawScore, node);
                    tls.local_result.updateCosineScore(nodeIndex, cosineScore, node);
                    tls.local_result.updateLogCosineOldScore(nodeIndex, logCosineOldScore, node);
                    tls.local_result.updateWeightedContainmentScore(nodeIndex, weightedContainmentScore, node);
                    tls.local_result.updateCapCosineScore(nodeIndex, capCosineScore, node);
                    tls.local_result.updateCapLogCosineScore(nodeIndex, capLogCosineScore, node);
                    
                    // Update thread-local results for NEW metrics (NO LOCKS - parallel performance!)
                    tls.local_result.updateLogRawScore(nodeIndex, logRawScore, node);
                    tls.local_result.updateLogCosineScore(nodeIndex, logCosineScore, node);
                    tls.local_result.updateLogLogRawScore(nodeIndex, logLogRawScore, node);
                    tls.local_result.updateLogLogCosineScore(nodeIndex, logLogCosineScore, node);
                    tls.local_result.updateContainmentScore(nodeIndex, containmentScore, node);
                    tls.local_result.updateConcordanceScore(nodeIndex, concordanceScore, node);
                    tls.local_result.updatePresenceScore(nodeIndex, presenceScore, node);
                    
                    // Store diagnostic data on node if enabled
                    if (params.store_diagnostics) {
                        node->genomeMagnitudeSquared = nodeMetrics.genomeMagnitudeSquared;
                        node->genomeUniqueSeedCount = nodeMetrics.genomeUniqueSeedCount;
                        node->genomeTotalSeedFrequency = nodeMetrics.genomeTotalSeedFrequency;
                        node->presenceIntersectionCount = nodeMetrics.presenceIntersectionCount;
                        node->cosineNumerator = nodeMetrics.cosineNumerator;
                        node->jaccardNumerator = nodeMetrics.presenceIntersectionCount; // Same as intersection
                        node->weightedJaccardNumerator = nodeMetrics.weightedJaccardNumerator;
                    }
                }

                if (!node->children.empty()) {
                    // Batch append to thread-local storage
                    tls.next_nodes.insert(tls.next_nodes.end(), node->children.begin(), node->children.end());
                    tls.next_metrics.insert(tls.next_metrics.end(), node->children.size(), nodeMetrics);
                    if (params.verify_scores) {
                        tls.next_seeds.insert(tls.next_seeds.end(), node->children.size(), incrementalGenomeCounts);
                    }
                }
            }
            nodesProcessed.fetch_add(batch_counter, std::memory_order_relaxed);
        });
        
        size_t currentCount = nodesProcessed.load(std::memory_order_relaxed);
        static std::atomic<size_t> lastReported{0};
        size_t last = lastReported.load(std::memory_order_relaxed);
        if (currentCount / 500 > last / 500) {
            lastReported.store(currentCount, std::memory_order_relaxed);
            logging::info("Processed {} nodes", currentCount);
        }
        
        // Parallel Merge Step - Optimized for minimal synchronization
        std::vector<ThreadLocalData*> active_tls;
        active_tls.reserve(thread_data.size());
        for (auto& tls : thread_data) {
            if (!tls.next_nodes.empty()) active_tls.push_back(&tls);
        }
        
        if (active_tls.empty()) break;
        
        // Calculate offsets for parallel merge
        std::vector<size_t> offsets(active_tls.size() + 1, 0);
        for (size_t i = 0; i < active_tls.size(); ++i) {
            offsets[i+1] = offsets[i] + active_tls[i]->next_nodes.size();
        }
        size_t total_next = offsets.back();
        
        // Resize destination vectors once
        current_nodes.clear(); current_nodes.resize(total_next);
        current_metrics.clear(); current_metrics.resize(total_next);
        if (params.verify_scores) {
            current_seeds.clear(); current_seeds.resize(total_next);
        }
        
        // Parallel copy/move from thread-local buffers to unified arrays
        tbb::parallel_for(tbb::blocked_range<size_t>(0, active_tls.size()), 
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto* src = active_tls[i];
                    size_t offset = offsets[i];
                    
                    std::copy(src->next_nodes.begin(), src->next_nodes.end(), current_nodes.begin() + offset);
                    std::copy(src->next_metrics.begin(), src->next_metrics.end(), current_metrics.begin() + offset);
                    if (params.verify_scores) {
                        std::move(src->next_seeds.begin(), src->next_seeds.end(), current_seeds.begin() + offset);
                    }
                    src->clear();
                }
            });
    }
    
    // ========================================
    // MERGE PHASE: Combine thread-local results into global result
    // ========================================
    logging::info("Merging {} thread-local results...", thread_data.size());
    auto merge_start = std::chrono::high_resolution_clock::now();
    
    for (auto& tls : thread_data) {
        // ========================================
        // Merge OLD metrics
        // ========================================
        
        // Merge RAW scores
        if (tls.local_result.bestRawNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestRawScore * 0.0001, 1e-9);
            if (tls.local_result.bestRawScore > result.bestRawScore + tolerance) {
                result.bestRawScore = tls.local_result.bestRawScore;
                result.bestRawNodeIndex = tls.local_result.bestRawNodeIndex;
                result.tiedRawNodeIndices = tls.local_result.tiedRawNodeIndices;
            } else if (tls.local_result.bestRawScore >= result.bestRawScore - tolerance) {
                result.tiedRawNodeIndices.insert(result.tiedRawNodeIndices.end(),
                    tls.local_result.tiedRawNodeIndices.begin(), tls.local_result.tiedRawNodeIndices.end());
            }
        }
        
        // Merge Cosine scores
        if (tls.local_result.bestCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestCosineScore > result.bestCosineScore + tolerance) {
                result.bestCosineScore = tls.local_result.bestCosineScore;
                result.bestCosineNodeIndex = tls.local_result.bestCosineNodeIndex;
                result.tiedCosineNodeIndices = tls.local_result.tiedCosineNodeIndices;
            } else if (tls.local_result.bestCosineScore >= result.bestCosineScore - tolerance) {
                result.tiedCosineNodeIndices.insert(result.tiedCosineNodeIndices.end(),
                    tls.local_result.tiedCosineNodeIndices.begin(), tls.local_result.tiedCosineNodeIndices.end());
            }
        }
        
        // Merge LogCosineOLD scores (the working formula)
        if (tls.local_result.bestLogCosineOldNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogCosineOldScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogCosineOldScore > result.bestLogCosineOldScore + tolerance) {
                result.bestLogCosineOldScore = tls.local_result.bestLogCosineOldScore;
                result.bestLogCosineOldNodeIndex = tls.local_result.bestLogCosineOldNodeIndex;
                result.tiedLogCosineOldNodeIndices = tls.local_result.tiedLogCosineOldNodeIndices;
            } else if (tls.local_result.bestLogCosineOldScore >= result.bestLogCosineOldScore - tolerance) {
                result.tiedLogCosineOldNodeIndices.insert(result.tiedLogCosineOldNodeIndices.end(),
                    tls.local_result.tiedLogCosineOldNodeIndices.begin(), tls.local_result.tiedLogCosineOldNodeIndices.end());
            }
        }
        
        // Merge WeightedContainment scores
        if (tls.local_result.bestWeightedContainmentNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestWeightedContainmentScore * 0.0001, 1e-9);
            if (tls.local_result.bestWeightedContainmentScore > result.bestWeightedContainmentScore + tolerance) {
                result.bestWeightedContainmentScore = tls.local_result.bestWeightedContainmentScore;
                result.bestWeightedContainmentNodeIndex = tls.local_result.bestWeightedContainmentNodeIndex;
                result.tiedWeightedContainmentNodeIndices = tls.local_result.tiedWeightedContainmentNodeIndices;
            } else if (tls.local_result.bestWeightedContainmentScore >= result.bestWeightedContainmentScore - tolerance) {
                result.tiedWeightedContainmentNodeIndices.insert(result.tiedWeightedContainmentNodeIndices.end(),
                    tls.local_result.tiedWeightedContainmentNodeIndices.begin(), tls.local_result.tiedWeightedContainmentNodeIndices.end());
            }
        }
        
        // Merge CapCosine scores
        if (tls.local_result.bestCapCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestCapCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestCapCosineScore > result.bestCapCosineScore + tolerance) {
                result.bestCapCosineScore = tls.local_result.bestCapCosineScore;
                result.bestCapCosineNodeIndex = tls.local_result.bestCapCosineNodeIndex;
                result.tiedCapCosineNodeIndices = tls.local_result.tiedCapCosineNodeIndices;
            } else if (tls.local_result.bestCapCosineScore >= result.bestCapCosineScore - tolerance) {
                result.tiedCapCosineNodeIndices.insert(result.tiedCapCosineNodeIndices.end(),
                    tls.local_result.tiedCapCosineNodeIndices.begin(), tls.local_result.tiedCapCosineNodeIndices.end());
            }
        }
        
        // Merge CapLogCosine scores
        if (tls.local_result.bestCapLogCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestCapLogCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestCapLogCosineScore > result.bestCapLogCosineScore + tolerance) {
                result.bestCapLogCosineScore = tls.local_result.bestCapLogCosineScore;
                result.bestCapLogCosineNodeIndex = tls.local_result.bestCapLogCosineNodeIndex;
                result.tiedCapLogCosineNodeIndices = tls.local_result.tiedCapLogCosineNodeIndices;
            } else if (tls.local_result.bestCapLogCosineScore >= result.bestCapLogCosineScore - tolerance) {
                result.tiedCapLogCosineNodeIndices.insert(result.tiedCapLogCosineNodeIndices.end(),
                    tls.local_result.tiedCapLogCosineNodeIndices.begin(), tls.local_result.tiedCapLogCosineNodeIndices.end());
            }
        }
        
        // ========================================
        // Merge NEW metrics
        // ========================================
        
        // Merge LogRAW scores
        if (tls.local_result.bestLogRawNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogRawScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogRawScore > result.bestLogRawScore + tolerance) {
                result.bestLogRawScore = tls.local_result.bestLogRawScore;
                result.bestLogRawNodeIndex = tls.local_result.bestLogRawNodeIndex;
                result.tiedLogRawNodeIndices = tls.local_result.tiedLogRawNodeIndices;
            } else if (tls.local_result.bestLogRawScore >= result.bestLogRawScore - tolerance) {
                result.tiedLogRawNodeIndices.insert(
                    result.tiedLogRawNodeIndices.end(),
                    tls.local_result.tiedLogRawNodeIndices.begin(),
                    tls.local_result.tiedLogRawNodeIndices.end()
                );
            }
        }
        
        // Merge LogCosine scores
        if (tls.local_result.bestLogCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogCosineScore > result.bestLogCosineScore + tolerance) {
                result.bestLogCosineScore = tls.local_result.bestLogCosineScore;
                result.bestLogCosineNodeIndex = tls.local_result.bestLogCosineNodeIndex;
                result.tiedLogCosineNodeIndices = tls.local_result.tiedLogCosineNodeIndices;
            } else if (tls.local_result.bestLogCosineScore >= result.bestLogCosineScore - tolerance) {
                result.tiedLogCosineNodeIndices.insert(
                    result.tiedLogCosineNodeIndices.end(),
                    tls.local_result.tiedLogCosineNodeIndices.begin(),
                    tls.local_result.tiedLogCosineNodeIndices.end()
                );
            }
        }
        
        // Merge LogLogRAW scores
        if (tls.local_result.bestLogLogRawNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogLogRawScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogLogRawScore > result.bestLogLogRawScore + tolerance) {
                result.bestLogLogRawScore = tls.local_result.bestLogLogRawScore;
                result.bestLogLogRawNodeIndex = tls.local_result.bestLogLogRawNodeIndex;
                result.tiedLogLogRawNodeIndices = tls.local_result.tiedLogLogRawNodeIndices;
            } else if (tls.local_result.bestLogLogRawScore >= result.bestLogLogRawScore - tolerance) {
                result.tiedLogLogRawNodeIndices.insert(
                    result.tiedLogLogRawNodeIndices.end(),
                    tls.local_result.tiedLogLogRawNodeIndices.begin(),
                    tls.local_result.tiedLogLogRawNodeIndices.end()
                );
            }
        }
        
        // Merge LogLogCosine scores
        if (tls.local_result.bestLogLogCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogLogCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogLogCosineScore > result.bestLogLogCosineScore + tolerance) {
                result.bestLogLogCosineScore = tls.local_result.bestLogLogCosineScore;
                result.bestLogLogCosineNodeIndex = tls.local_result.bestLogLogCosineNodeIndex;
                result.tiedLogLogCosineNodeIndices = tls.local_result.tiedLogLogCosineNodeIndices;
            } else if (tls.local_result.bestLogLogCosineScore >= result.bestLogLogCosineScore - tolerance) {
                result.tiedLogLogCosineNodeIndices.insert(
                    result.tiedLogLogCosineNodeIndices.end(),
                    tls.local_result.tiedLogLogCosineNodeIndices.begin(),
                    tls.local_result.tiedLogLogCosineNodeIndices.end()
                );
            }
        }
        
        // Merge Containment scores
        if (tls.local_result.bestContainmentNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestContainmentScore * 0.0001, 1e-9);
            if (tls.local_result.bestContainmentScore > result.bestContainmentScore + tolerance) {
                result.bestContainmentScore = tls.local_result.bestContainmentScore;
                result.bestContainmentNodeIndex = tls.local_result.bestContainmentNodeIndex;
                result.tiedContainmentNodeIndices = tls.local_result.tiedContainmentNodeIndices;
            } else if (tls.local_result.bestContainmentScore >= result.bestContainmentScore - tolerance) {
                result.tiedContainmentNodeIndices.insert(
                    result.tiedContainmentNodeIndices.end(),
                    tls.local_result.tiedContainmentNodeIndices.begin(),
                    tls.local_result.tiedContainmentNodeIndices.end()
                );
            }
        }
        
        // Merge Concordance scores
        if (tls.local_result.bestConcordanceNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestConcordanceScore * 0.0001, 1e-9);
            if (tls.local_result.bestConcordanceScore > result.bestConcordanceScore + tolerance) {
                result.bestConcordanceScore = tls.local_result.bestConcordanceScore;
                result.bestConcordanceNodeIndex = tls.local_result.bestConcordanceNodeIndex;
                result.tiedConcordanceNodeIndices = tls.local_result.tiedConcordanceNodeIndices;
            } else if (tls.local_result.bestConcordanceScore >= result.bestConcordanceScore - tolerance) {
                result.tiedConcordanceNodeIndices.insert(
                    result.tiedConcordanceNodeIndices.end(),
                    tls.local_result.tiedConcordanceNodeIndices.begin(),
                    tls.local_result.tiedConcordanceNodeIndices.end()
                );
            }
        }
        
        // Merge Presence scores
        if (tls.local_result.bestPresenceNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestPresenceScore * 0.0001, 1e-9);
            if (tls.local_result.bestPresenceScore > result.bestPresenceScore + tolerance) {
                result.bestPresenceScore = tls.local_result.bestPresenceScore;
                result.bestPresenceNodeIndex = tls.local_result.bestPresenceNodeIndex;
                result.tiedPresenceNodeIndices = tls.local_result.tiedPresenceNodeIndices;
            } else if (tls.local_result.bestPresenceScore >= result.bestPresenceScore - tolerance) {
                result.tiedPresenceNodeIndices.insert(
                    result.tiedPresenceNodeIndices.end(),
                    tls.local_result.tiedPresenceNodeIndices.begin(),
                    tls.local_result.tiedPresenceNodeIndices.end()
                );
            }
        }
    }
    
    auto merge_end = std::chrono::high_resolution_clock::now();
    auto merge_ms = std::chrono::duration_cast<std::chrono::milliseconds>(merge_end - merge_start).count();
    logging::info("Thread-local merge completed in {}ms", merge_ms);
}

void placeLite(PlacementResult &result, 
                       panmapUtils::LiteTree *liteTree,
                       ::capnp::MessageReader &liteIndex, 
                       const std::string &reads1,
                       const std::string &reads2,
                       std::string &outputPath,
                       bool verify_scores,
                       panmanUtils::Tree *fullTree,
                       const std::string &skipNodeId,
                       std::string *parentOfSkippedNode,
                       bool store_diagnostics,
                       double seedMaskFraction,
                       int minSeedQuality,
                       bool dedupReads,
                       bool pairFilter,
                       int trimStart,
                       int trimEnd,
                       int minReadSupport,
                       uint32_t expectedFragmentSize,
                       uint32_t fragmentSizeTolerance,
                       bool refineEnabled,
                       double refineTopPct,
                       int refineMaxTopN,
                       int refineNeighborRadius,
                       int refineMaxNeighborN) {
    auto placement_total_start = std::chrono::high_resolution_clock::now();
    logging::info("Starting lite-index placement");
    
    // Handle leave-one-out validation: find skip node and its parent
    uint32_t skipNodeIndex = UINT32_MAX;
    if (!skipNodeId.empty()) {
        auto it = liteTree->allLiteNodes.find(skipNodeId);
        if (it == liteTree->allLiteNodes.end()) {
            logging::err("Skip node '{}' not found in tree", skipNodeId);
            throw std::runtime_error("Skip node not found: " + skipNodeId);
        }
        panmapUtils::LiteNode* skipNode = it->second;
        
        // Validate it's a leaf node
        if (!skipNode->children.empty()) {
            logging::err("Skip node '{}' has {} children - only leaf nodes supported for leave-one-out", 
                        skipNodeId, skipNode->children.size());
            throw std::runtime_error("Skip node must be a leaf: " + skipNodeId);
        }
        
        skipNodeIndex = skipNode->nodeIndex;
        logging::info("Leave-one-out mode: skipping leaf node '{}' (index {})", skipNodeId, skipNodeIndex);
        
        // Return parent node ID if requested
        if (parentOfSkippedNode != nullptr) {
            if (skipNode->parent != nullptr) {
                *parentOfSkippedNode = skipNode->parent->identifier;
                logging::info("Expected placement (parent of skipped node): '{}'", *parentOfSkippedNode);
            } else {
                *parentOfSkippedNode = "";
                logging::warn("Skipped node has no parent (is root?)");
            }
        }
    }
    
    // Store verify_scores flag in a place accessible to score calculation
    if (verify_scores) {
        logging::info("VERIFICATION MODE: Will recompute all scores from scratch at each node");
        if (fullTree == nullptr) {
            logging::err("VERIFICATION MODE requires fullTree pointer but got nullptr!");
            throw std::runtime_error("Verification mode requires full Tree to be loaded");
        }
    }
    
    // State initialization moved to after parameter setup
    
    logging::info("Loading pre-computed hash deltas from index for {} nodes...", liteTree->allLiteNodes.size());
    
    auto time_hash_delta_start = std::chrono::high_resolution_clock::now();
    size_t totalHashDeltas = 0;
    
    auto indexRoot = liteIndex.getRoot<LiteIndex>();
    auto genomeMagSquaredReader = indexRoot.getGenomeMagnitudeSquared();
    auto genomeUniqueSeedReader = indexRoot.getGenomeUniqueSeedCount();
    auto genomeTotalFreqReader = indexRoot.getGenomeTotalSeedFrequency();
    
    // OPTIMIZATION: Process nodes in DFS order with direct vector indexing for O(1) access
    // This eliminates hash map lookup overhead in the hot path
    const size_t numNodes = liteTree->allLiteNodes.size();
    
    // Parallel Pass 1: Count seed changes per node to calculate offsets
    std::vector<size_t> nodeOffsets(numNodes + 1, 0);
    
    std::atomic<size_t> nodesWithChanges(0);

    {

        
        if (!indexRoot.hasSeedChangeHashes() || !indexRoot.hasSeedChangeParentCounts() || 
            !indexRoot.hasSeedChangeChildCounts() || !indexRoot.hasNodeChangeOffsets()) {
             throw std::runtime_error("Index missing required V3 fields (seedChangeHashes, etc). V2 is no longer supported.");
        }

        logging::info("Using VERSION 3 struct-of-arrays format from index (optimal parallel access)");
        auto indexOffsets = indexRoot.getNodeChangeOffsets();
        if (indexOffsets.size() >= numNodes + 1) {
            for (size_t i = 0; i <= numNodes; ++i) {
                nodeOffsets[i] = indexOffsets[i];
            }
        } else {
             throw std::runtime_error(fmt::format("Struct-of-arrays format offsets size mismatch: {} vs {}", indexOffsets.size(), numNodes + 1));
        }
        totalHashDeltas = nodeOffsets[numNodes];
        
        // Allocate backing storage once
        liteTree->allSeedChanges.resize(totalHashDeltas);
        
        // Parallel Pass 2: Populate data and assign spans
        
        // VERSION 3: Struct-of-arrays format
        auto hashesReader = indexRoot.getSeedChangeHashes();
        auto parentCountsReader = indexRoot.getSeedChangeParentCounts();
        auto childCountsReader = indexRoot.getSeedChangeChildCounts();
        
        logging::info("Struct-of-arrays format: {} hashes, {} parent counts, {} child counts",
                        hashesReader.size(), parentCountsReader.size(), childCountsReader.size());
        
        // Phase 2a
        const size_t GRAIN_SIZE = 262144; // 256K items per chunk
        tbb::parallel_for(tbb::blocked_range<size_t>(0, totalHashDeltas, GRAIN_SIZE), 
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    liteTree->allSeedChanges[i] = std::make_tuple(
                        hashesReader[i],
                        parentCountsReader[i],
                        childCountsReader[i]
                    );
                }
            });
        
        // Phase 2b: Assign spans to nodes
        tbb::parallel_for(tbb::blocked_range<size_t>(0, numNodes), 
            [&](const tbb::blocked_range<size_t>& range) {
                size_t localNodesWithChanges = 0;
                for (size_t dfsIndex = range.begin(); dfsIndex < range.end(); ++dfsIndex) {
                    auto* liteNode = liteTree->dfsIndexToNode[dfsIndex];
                    size_t offset = nodeOffsets[dfsIndex];
                    size_t size = nodeOffsets[dfsIndex+1] - offset;
                    
                    if (size > 0) {
                        liteNode->seedChanges = std::span<std::tuple<uint64_t, int64_t, int64_t>>(
                            liteTree->allSeedChanges.data() + offset, size);
                        localNodesWithChanges++;
                    } else {
                        liteNode->seedChanges = std::span<std::tuple<uint64_t, int64_t, int64_t>>();
                    }
                }
                nodesWithChanges.fetch_add(localNodesWithChanges, std::memory_order_relaxed);
            });
    }
    
    auto time_hash_delta_end = std::chrono::high_resolution_clock::now();
    auto duration_hash_delta = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_hash_delta_end - time_hash_delta_start);
    logging::info("Hash delta loading: {}ms", duration_hash_delta.count());
    logging::info("Loaded {} hash deltas from index", totalHashDeltas);
    
    // Set traversal parameters
    TraversalParams params;
    params.k = indexRoot.getK();
    params.s = indexRoot.getS();
    params.t = indexRoot.getT();
    params.l = indexRoot.getL();
    params.open = indexRoot.getOpen();
    params.verify_scores = verify_scores;
    params.store_diagnostics = store_diagnostics;
    params.dedupReads = dedupReads;
    params.pairFilter = pairFilter;
    params.trimStart = trimStart;
    params.trimEnd = trimEnd;
    params.minReadSupport = minReadSupport;
    
    // Set refinement parameters
    params.refineEnabled = refineEnabled;
    params.refineTopPct = refineTopPct;
    params.refineMaxTopN = refineMaxTopN;
    params.refineNeighborRadius = refineNeighborRadius;
    params.refineMaxNeighborN = refineMaxNeighborN;
    
    if (dedupReads) {
        logging::info("Read deduplication enabled (--dedup): counting each unique sequence once");
    }
    if (trimStart > 0 || trimEnd > 0) {
        logging::info("Read trimming enabled: trimStart={}, trimEnd={} (for primer removal)", trimStart, trimEnd);
    }
    if (pairFilter && !reads2.empty()) {
        logging::info("Paired-end concordance scoring enabled: fragments with both mates matching weighted higher");
        logging::info("  Expected fragment size: {}bp ± {}bp", expectedFragmentSize, fragmentSizeTolerance);
    } else if (pairFilter && reads2.empty()) {
        logging::info("Paired-end scoring requested but no R2 reads provided - disabled");
    }
    
    std::vector<std::string> allReadSequences;
    std::vector<std::string> allReadQualities;  // Quality strings (only populated if minSeedQuality > 0)
    
    // Initialize placement global state
    PlacementGlobalState state;
    state.kmerSize = indexRoot.getK();
    state.fullTree = fullTree;  // Store full tree pointer for verification mode
    state.liteNodes = liteTree->dfsIndexToNode;
    state.skipNodeIndex = skipNodeIndex;  // Set skip node for leave-one-out validation
    state.expectedFragmentSize = expectedFragmentSize;
    state.fragmentSizeTolerance = fragmentSizeTolerance;
    
    auto time_read_processing_start = std::chrono::high_resolution_clock::now();
    {
        if (!reads1.empty()) {
            if (minSeedQuality > 0) {
                // Need quality strings for per-seed quality filtering
                std::vector<std::string> readNames;  // Unused but required by function
                extractFullFastqData(reads1, reads2, allReadSequences, allReadQualities, readNames);
                logging::info("Loaded {} reads with quality strings for Q{} filtering", 
                            allReadSequences.size(), minSeedQuality);
            } else {
                extractReadSequences(reads1, reads2, allReadSequences);
            }
        }
        
        // Build paired-end fragment tracking if R2 was provided and pairFilter is enabled
        if (pairFilter && !reads2.empty() && !allReadSequences.empty()) {
            buildPairedEndFragmentData(state, allReadSequences, params);
        }
        
        uint32_t l = indexRoot.getL();
        logging::info("Index parameters: k={}, s={}, t={}, l={}", 
                    params.k, params.s, params.t, l);
        
        if (l == 0) {
            logging::info("l=0: Using raw syncmers instead of k-minimizers");
            
            // Quality-filtered seed extraction path (no read deduplication since quality varies)
            if (minSeedQuality > 0 && !allReadQualities.empty()) {
                auto time_seed_extract_start = std::chrono::high_resolution_clock::now();
                
                size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
                std::vector<absl::flat_hash_map<size_t, int64_t>> threadLocalMaps(num_cpus);
                std::vector<std::atomic<size_t>> localSeedsFiltered(num_cpus);
                
                const size_t avg_seeds_per_read = 50;
                const size_t reads_per_thread = (allReadSequences.size() + num_cpus - 1) / num_cpus;
                const size_t estimated_seeds = avg_seeds_per_read * reads_per_thread;
                const size_t reserve_capacity = (estimated_seeds * 3) / 2;
                
                for (auto& localMap : threadLocalMaps) {
                    localMap.reserve(reserve_capacity);
                }
                for (auto& count : localSeedsFiltered) {
                    count.store(0, std::memory_order_relaxed);
                }
                
                tbb::parallel_for(tbb::blocked_range<size_t>(0, allReadSequences.size()), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        size_t threadId = tbb::this_task_arena::current_thread_index();
                        auto& localMap = threadLocalMaps[threadId];
                        size_t filtered = 0;
                        size_t trimFiltered = 0;
                        
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            const std::string& seq = allReadSequences[i];
                            const std::string& qual = allReadQualities[i];
                            
                            // Calculate trim boundaries for this read
                            const int seqLen = static_cast<int>(seq.size());
                            const int trimStartBp = params.trimStart;
                            const int trimEndBp = params.trimEnd;
                            const int validStart = trimStartBp;
                            const int validEnd = seqLen - trimEndBp - params.k;  // Last valid startPos for k-mer
                            
                            const auto& syncmers = seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                            for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                                if (!isSyncmer) [[unlikely]] continue;
                                
                                // Trim filter: skip seeds that start in trimmed regions
                                if (static_cast<int>(startPos) < validStart || static_cast<int>(startPos) > validEnd) {
                                    ++trimFiltered;
                                    continue;  // Skip seed in primer region
                                }
                                
                                // Quality filter: check average quality over k-mer span
                                double avgQual = avgPhredQuality(qual, startPos, params.k);
                                if (avgQual < static_cast<double>(minSeedQuality)) {
                                    ++filtered;
                                    continue;  // Skip low-quality seed
                                }
                                
                                localMap[kmerHash] += 1;
                            }
                        }
                        localSeedsFiltered[threadId].fetch_add(filtered + trimFiltered, std::memory_order_relaxed);
                    });
                
                size_t totalFiltered = 0;
                for (const auto& count : localSeedsFiltered) {
                    totalFiltered += count.load(std::memory_order_relaxed);
                }
                
                // Merge thread-local maps
                size_t exact_total_seeds = 0;
                for (const auto& localMap : threadLocalMaps) {
                    exact_total_seeds += localMap.size();
                }
                state.seedFreqInReads.reserve((exact_total_seeds * 23) / 20);
                for (const auto& localMap : threadLocalMaps) {
                    for (const auto& [hash, count] : localMap) {
                        state.seedFreqInReads[hash] += count;
                    }
                }
                
                auto time_seed_extract_end = std::chrono::high_resolution_clock::now();
                auto duration_seed_extract = std::chrono::duration_cast<std::chrono::milliseconds>(
                    time_seed_extract_end - time_seed_extract_start);
                logging::info("Seed extraction (Q{} filtered): {}ms", 
                            minSeedQuality, duration_seed_extract.count());
                logging::info("Extracted {} unique syncmers from {} reads (filtered {} low-quality seeds)", 
                            state.seedFreqInReads.size(), allReadSequences.size(), totalFiltered);
                
            } else {
                // Original path: deduplicate reads first for efficiency
                std::vector<std::pair<std::string, size_t>> sortedReads;
            auto time_dedup_start = std::chrono::high_resolution_clock::now();
            {
                sortedReads.reserve(allReadSequences.size());
                for (size_t i = 0; i < allReadSequences.size(); i++) {
                    sortedReads.emplace_back(allReadSequences[i], i);
                }
                tbb::parallel_sort(sortedReads.begin(), sortedReads.end(), 
                    [](const auto& a, const auto& b) { return a.first < b.first; });
            }
            auto time_dedup_end = std::chrono::high_resolution_clock::now();
            auto duration_dedup = std::chrono::duration_cast<std::chrono::milliseconds>(
                time_dedup_end - time_dedup_start);
            logging::info("Read deduplication: {}ms", duration_dedup.count());
            
            // Estimate unique count (typically 50-90% unique for paired-end reads)
            std::vector<std::pair<std::string, std::vector<size_t>>> dupReadsIndex;
            
            if (!sortedReads.empty()) {
                // Parallel identification of unique sequence boundaries
                std::vector<uint8_t> is_new_seq(sortedReads.size());
                is_new_seq[0] = 1;
                
                tbb::parallel_for(tbb::blocked_range<size_t>(1, sortedReads.size()), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            is_new_seq[i] = (sortedReads[i].first != sortedReads[i-1].first) ? 1 : 0;
                        }
                    });
                
                std::vector<size_t> unique_starts;
                unique_starts.reserve(sortedReads.size() / 2);
                for (size_t i = 0; i < is_new_seq.size(); ++i) {
                    if (is_new_seq[i]) unique_starts.push_back(i);
                }
                unique_starts.push_back(sortedReads.size());
                
                dupReadsIndex.resize(unique_starts.size() - 1);
                
                tbb::parallel_for(tbb::blocked_range<size_t>(0, unique_starts.size() - 1), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            size_t start = unique_starts[i];
                            size_t end = unique_starts[i+1];
                            
                            std::vector<size_t> duplicates;
                            duplicates.reserve(end - start);
                            for (size_t j = start; j < end; ++j) {
                                duplicates.push_back(sortedReads[j].second);
                            }
                            dupReadsIndex[i] = std::make_pair(sortedReads[start].first, std::move(duplicates));
                        }
                    });
            }
            
            logging::info("Deduplication: {} reads → {} unique sequences", 
                        allReadSequences.size(), dupReadsIndex.size());
            
            auto time_seed_extract_start = std::chrono::high_resolution_clock::now();
            {
                size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
                std::vector<absl::flat_hash_map<size_t, int64_t>> threadLocalMaps(num_cpus);
                
                // Pre-allocate capacity based on empirical analysis to eliminate malloc overhead
                // Profiler shows 14% time in malloc (7.2% consolidate + 6.9% perturb)
                // Formula: estimate seeds/read * reads/thread * load_factor_headroom
                const size_t avg_seeds_per_read = 50; // Conservative estimate for syncmers
                const size_t reads_per_thread = (dupReadsIndex.size() + num_cpus - 1) / num_cpus;
                const size_t estimated_seeds = avg_seeds_per_read * reads_per_thread;
                // Reserve 1.5x for hash table load factor (typically 0.875 max)
                const size_t reserve_capacity = (estimated_seeds * 3) / 2;
                
                for (auto& localMap : threadLocalMaps) {
                    localMap.reserve(reserve_capacity);
                }
                
                tbb::parallel_for(tbb::blocked_range<size_t>(0, dupReadsIndex.size()), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        size_t threadId = tbb::this_task_arena::current_thread_index();
                        auto& localMap = threadLocalMaps[threadId];
                        
                        // Process reads in this thread's range
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            const auto& [seq, duplicates] = dupReadsIndex[i];
                            // If dedupReads is true, count each unique sequence once
                            // Otherwise, use the actual multiplicity (number of duplicate reads)
                            const int64_t multiplicity = params.dedupReads ? 1 : static_cast<int64_t>(duplicates.size());
                            
                            // Calculate trim boundaries for this read
                            const int seqLen = static_cast<int>(seq.size());
                            const int trimStartBp = params.trimStart;
                            const int trimEndBp = params.trimEnd;
                            const int validStart = trimStartBp;
                            const int validEnd = seqLen - trimEndBp - params.k;
                            
                            const auto& syncmers = seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                            for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                                if (!isSyncmer) [[unlikely]] continue;
                                
                                // Trim filter: skip seeds that start in trimmed regions
                                if (static_cast<int>(startPos) < validStart || static_cast<int>(startPos) > validEnd) {
                                    continue;  // Skip seed in primer region
                                }
                                
                                localMap[kmerHash] += multiplicity;
                            }
                        }
                    });
                
                // Calculate exact size needed for final merge to avoid any resizing
                size_t exact_total_seeds = 0;
                for (const auto& localMap : threadLocalMaps) {
                    exact_total_seeds += localMap.size();
                }
                
                // Reserve with headroom for hash table load factor
                // absl::flat_hash_map maintains load factor ~0.875, so reserve 1.15x
                state.seedFreqInReads.reserve((exact_total_seeds * 23) / 20);
                
                // Parallel merge: each thread merges its own map into the final result
                // Use atomic operations or partitioned hash ranges
                // For simplicity, merge serially since it's fast (< 100ms typically)
                // The extraction phase is the real bottleneck, not the merge
                for (const auto& localMap : threadLocalMaps) {
                    for (const auto& [hash, count] : localMap) {
                        state.seedFreqInReads[hash] += count;
                    }
                }
            }
            auto time_seed_extract_end = std::chrono::high_resolution_clock::now();
            auto duration_seed_extract = std::chrono::duration_cast<std::chrono::milliseconds>(
                time_seed_extract_end - time_seed_extract_start);
            logging::info("Seed extraction: {}ms", duration_seed_extract.count());
            
            logging::info("Extracted {} unique syncmers from {} reads ({} unique patterns)", 
                        state.seedFreqInReads.size(), allReadSequences.size(), dupReadsIndex.size());
            }  // end of else block for non-quality-filtered path
            
        } else {
            // l > 0: Use k-minimizers (consecutive syncmers combined)
            
            // Quality-filtered k-minimizer extraction (no read deduplication)
            if (minSeedQuality > 0 && !allReadQualities.empty()) {
                auto time_kminimizer_start = std::chrono::high_resolution_clock::now();
                
                size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
                std::vector<absl::flat_hash_map<size_t, int64_t>> threadLocalMaps(num_cpus);
                std::vector<std::atomic<size_t>> localSeedsFiltered(num_cpus);
                
                const size_t avg_kminmers_per_read = 15;
                const size_t reads_per_thread = (allReadSequences.size() + num_cpus - 1) / num_cpus;
                const size_t estimated_kminmers = avg_kminmers_per_read * reads_per_thread;
                const size_t reserve_capacity = (estimated_kminmers * 3) / 2;
                
                for (auto& localMap : threadLocalMaps) {
                    localMap.reserve(reserve_capacity);
                }
                for (auto& count : localSeedsFiltered) {
                    count.store(0, std::memory_order_relaxed);
                }
                
                tbb::parallel_for(tbb::blocked_range<size_t>(0, allReadSequences.size()), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        size_t threadId = tbb::this_task_arena::current_thread_index();
                        auto& localMap = threadLocalMaps[threadId];
                        size_t filtered = 0;
                        
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            const std::string& seq = allReadSequences[i];
                            const std::string& qual = allReadQualities[i];
                            
                            const auto& syncmers = seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                            
                            if (syncmers.size() < static_cast<size_t>(params.l)) continue;
                            
                            // For l=1, just use syncmer hashes with quality filtering
                            if (params.l == 1) {
                                for (const auto& [seedHash, isReverse, isSyncmer, startPos] : syncmers) {
                                    // Quality filter
                                    double avgQual = avgPhredQuality(qual, startPos, params.k);
                                    if (avgQual < static_cast<double>(minSeedQuality)) {
                                        ++filtered;
                                        continue;
                                    }
                                    localMap[seedHash] += 1;
                                }
                                continue;
                            }
                            
                            // For l > 1: build k-minimizers from quality-filtered syncmers
                            // Pre-compute which syncmers pass quality filter
                            std::vector<bool> syncmerPassesQuality(syncmers.size());
                            for (size_t j = 0; j < syncmers.size(); ++j) {
                                int64_t startPos = std::get<3>(syncmers[j]);
                                double avgQual = avgPhredQuality(qual, startPos, params.k);
                                syncmerPassesQuality[j] = (avgQual >= static_cast<double>(minSeedQuality));
                            }
                            
                            // Build k-minimizers, skipping any window with low-quality syncmers
                            size_t forwardRolledHash = 0;
                            size_t reverseRolledHash = 0;
                            
                            // First k-min-mer
                            bool firstWindowValid = true;
                            for (size_t j = 0; j < static_cast<size_t>(params.l); ++j) {
                                if (!syncmerPassesQuality[j]) {
                                    firstWindowValid = false;
                                    break;
                                }
                                forwardRolledHash = seeding::rol(forwardRolledHash, params.k) ^ std::get<0>(syncmers[j]);
                                reverseRolledHash = seeding::rol(reverseRolledHash, params.k) ^ std::get<0>(syncmers[params.l-j-1]);
                            }
                            
                            if (firstWindowValid && forwardRolledHash != reverseRolledHash) {
                                size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
                                localMap[minHash] += 1;
                            } else if (!firstWindowValid) {
                                ++filtered;
                            }
                            
                            // Rest of k-min-mers using rolling hash
                            for (size_t j = 1; j < syncmers.size() - params.l + 1; ++j) {
                                // Check if new syncmer at end of window passes quality
                                if (!syncmerPassesQuality[j + params.l - 1]) {
                                    ++filtered;
                                    // Reset hash (need to rebuild for next valid window)
                                    forwardRolledHash = 0;
                                    reverseRolledHash = 0;
                                    continue;
                                }
                                
                                // Check if all syncmers in current window pass quality
                                bool windowValid = true;
                                for (size_t w = j; w < j + static_cast<size_t>(params.l); ++w) {
                                    if (!syncmerPassesQuality[w]) {
                                        windowValid = false;
                                        break;
                                    }
                                }
                                
                                if (!windowValid) {
                                    ++filtered;
                                    continue;
                                }
                                
                                // Rebuild hash for this window (simpler than incremental for quality filtering)
                                forwardRolledHash = 0;
                                reverseRolledHash = 0;
                                for (size_t w = 0; w < static_cast<size_t>(params.l); ++w) {
                                    forwardRolledHash = seeding::rol(forwardRolledHash, params.k) ^ std::get<0>(syncmers[j+w]);
                                    reverseRolledHash = seeding::rol(reverseRolledHash, params.k) ^ std::get<0>(syncmers[j+params.l-w-1]);
                                }
                                
                                if (forwardRolledHash != reverseRolledHash) {
                                    size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
                                    localMap[minHash] += 1;
                                }
                            }
                        }
                        localSeedsFiltered[threadId].fetch_add(filtered, std::memory_order_relaxed);
                    });
                
                size_t totalFiltered = 0;
                for (const auto& count : localSeedsFiltered) {
                    totalFiltered += count.load(std::memory_order_relaxed);
                }
                
                // Merge thread-local maps
                size_t exact_total_seeds = 0;
                for (const auto& localMap : threadLocalMaps) {
                    exact_total_seeds += localMap.size();
                }
                state.seedFreqInReads.reserve((exact_total_seeds * 23) / 20);
                for (const auto& localMap : threadLocalMaps) {
                    for (const auto& [hash, count] : localMap) {
                        state.seedFreqInReads[hash] += count;
                    }
                }
                
                auto time_kminimizer_end = std::chrono::high_resolution_clock::now();
                auto duration_kminimizer = std::chrono::duration_cast<std::chrono::milliseconds>(
                    time_kminimizer_end - time_kminimizer_start);
                logging::info("K-minimizer extraction (Q{} filtered): {}ms", 
                            minSeedQuality, duration_kminimizer.count());
                logging::info("Extracted {} unique k-minimizers from {} reads (filtered {} low-quality seeds)", 
                            state.seedFreqInReads.size(), allReadSequences.size(), totalFiltered);
                
            } else {
                // Original path: deduplicate reads first for efficiency
                auto time_kminimizer_start = std::chrono::high_resolution_clock::now();
            
            // Deduplicate reads first (same as l=0 case)
            std::vector<std::pair<std::string, size_t>> sortedReads;
            sortedReads.reserve(allReadSequences.size());
            for (size_t i = 0; i < allReadSequences.size(); i++) {
                sortedReads.emplace_back(allReadSequences[i], i);
            }
            tbb::parallel_sort(sortedReads.begin(), sortedReads.end(), 
                [](const auto& a, const auto& b) { return a.first < b.first; });

            std::vector<std::pair<std::string, std::vector<size_t>>> dupReadsIndex;
            if (!sortedReads.empty()) {
                std::vector<uint8_t> is_new_seq(sortedReads.size());
                is_new_seq[0] = 1;
                tbb::parallel_for(tbb::blocked_range<size_t>(1, sortedReads.size()), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            is_new_seq[i] = (sortedReads[i].first != sortedReads[i-1].first) ? 1 : 0;
                        }
                    });
                
                std::vector<size_t> unique_starts;
                unique_starts.reserve(sortedReads.size() / 2);
                for (size_t i = 0; i < is_new_seq.size(); ++i) {
                    if (is_new_seq[i]) unique_starts.push_back(i);
                }
                unique_starts.push_back(sortedReads.size());
                
                dupReadsIndex.resize(unique_starts.size() - 1);
                tbb::parallel_for(tbb::blocked_range<size_t>(0, unique_starts.size() - 1), 
                    [&](const tbb::blocked_range<size_t>& range) {
                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            size_t start = unique_starts[i];
                            size_t end = unique_starts[i+1];
                            std::vector<size_t> duplicates;
                            duplicates.reserve(end - start);
                            for (size_t j = start; j < end; ++j) {
                                duplicates.push_back(sortedReads[j].second);
                            }
                            dupReadsIndex[i] = std::make_pair(sortedReads[start].first, std::move(duplicates));
                        }
                    });
            }

            logging::info("Deduplication: {} reads → {} unique sequences", 
                        allReadSequences.size(), dupReadsIndex.size());

            // Extract k-min-mers in parallel
            size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
            std::vector<absl::flat_hash_map<size_t, int64_t>> threadLocalMaps(num_cpus);
            
            // Estimate capacity
            const size_t avg_kminmers_per_read = 15;
            const size_t reads_per_thread = (dupReadsIndex.size() + num_cpus - 1) / num_cpus;
            const size_t estimated_kminmers = avg_kminmers_per_read * reads_per_thread;
            const size_t reserve_capacity = (estimated_kminmers * 3) / 2;
            
            for (auto& localMap : threadLocalMaps) {
                localMap.reserve(reserve_capacity);
            }

            tbb::parallel_for(tbb::blocked_range<size_t>(0, dupReadsIndex.size()), 
                [&](const tbb::blocked_range<size_t>& range) {
                    size_t threadId = tbb::this_task_arena::current_thread_index();
                    auto& localMap = threadLocalMaps[threadId];
                    
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto& [seq, duplicates] = dupReadsIndex[i];
                        // If dedupReads is true, count each unique sequence once
                        const int64_t multiplicity = params.dedupReads ? 1 : static_cast<int64_t>(duplicates.size());
                        
                        const auto& syncmers = seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                        
                        if (syncmers.size() < params.l) continue;

                        // Special case for l=1: just use syncmer hashes directly
                        if (params.l == 1) {
                            for (const auto& syncmer : syncmers) {
                                size_t seedHash = std::get<0>(syncmer);
                                localMap[seedHash] += multiplicity;
                            }
                            continue;
                        }

                        size_t forwardRolledHash = 0;
                        size_t reverseRolledHash = 0;
                        
                        // First k-min-mer
                        for (size_t j = 0; j < params.l; ++j) {
                            forwardRolledHash = seeding::rol(forwardRolledHash, params.k) ^ std::get<0>(syncmers[j]);
                            reverseRolledHash = seeding::rol(reverseRolledHash, params.k) ^ std::get<0>(syncmers[params.l-j-1]);
                        }

                        if (forwardRolledHash != reverseRolledHash) {
                            size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
                            localMap[minHash] += multiplicity;
                        }

                        // Rest of k-min-mers
                        for (size_t j = 1; j < syncmers.size() - params.l + 1; ++j) {
                            const size_t& prevSyncmerHash = std::get<0>(syncmers[j-1]);
                            const size_t& nextSyncmerHash = std::get<0>(syncmers[j+params.l-1]);
                            
                            forwardRolledHash = seeding::rol(forwardRolledHash, params.k) ^ seeding::rol(prevSyncmerHash, params.k * params.l) ^ nextSyncmerHash;
                            reverseRolledHash = seeding::ror(reverseRolledHash, params.k) ^ seeding::ror(prevSyncmerHash, params.k) ^ seeding::rol(nextSyncmerHash, params.k * (params.l-1));

                            if (forwardRolledHash != reverseRolledHash) {
                                size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
                                localMap[minHash] += multiplicity;
                            }
                        }
                    }
                });

            // Merge results
            size_t exact_total_seeds = 0;
            for (const auto& localMap : threadLocalMaps) {
                exact_total_seeds += localMap.size();
            }
            state.seedFreqInReads.reserve((exact_total_seeds * 23) / 20);
            
            for (const auto& localMap : threadLocalMaps) {
                for (const auto& [hash, count] : localMap) {
                    state.seedFreqInReads[hash] += count;
                }
            }
            
            auto time_kminimizer_end = std::chrono::high_resolution_clock::now();
            auto duration_kminimizer = std::chrono::duration_cast<std::chrono::milliseconds>(
                time_kminimizer_end - time_kminimizer_start);
            logging::info("K-minimizer extraction: {}ms", duration_kminimizer.count());
            logging::info("Extracted {} unique k-minimizers from {} reads", 
                        state.seedFreqInReads.size(), allReadSequences.size());
            }  // end of else block for non-quality-filtered path (l > 0)
        }
    }
    auto time_read_processing_end = std::chrono::high_resolution_clock::now();
    auto duration_read_processing = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_read_processing_end - time_read_processing_start);
    logging::info("Total read processing: {}ms", duration_read_processing.count());
    
    // ==========================================================================
    // Seed frequency analysis and masking
    // Identifies and optionally removes over-represented seeds (primer/adapter artifacts)
    // ==========================================================================
    if (!allReadSequences.empty()) {
        const size_t totalReads = allReadSequences.size();
        size_t uniqueSeeds = state.seedFreqInReads.size();
        
        // Build hash-to-sequence map ONLY in debug mode (very slow for large datasets)
        // Skip for normal runs - the seed_freq.tsv will just show hashes without sequences
        // This is O(reads * seeds) and can take minutes for 1x coverage data
        
        // Remove homopolymer seeds (all same base - uninformative)
        auto homoHashes = getHomopolymerHashes(params.k);
        size_t homoRemoved = 0;
        for (size_t homoHash : homoHashes) {
            auto it = state.seedFreqInReads.find(homoHash);
            if (it != state.seedFreqInReads.end()) {
                logging::info("Removing homopolymer seed: hash={:016x}, count={}", homoHash, it->second);
                state.seedFreqInReads.erase(it);
                homoRemoved++;
            }
        }
        if (homoRemoved > 0) {
            logging::info("Removed {} homopolymer seeds", homoRemoved);
            uniqueSeeds = state.seedFreqInReads.size();
        }
        
        // Count seeds above various thresholds (no sorting needed)
        size_t above50pct = 0, above20pct = 0, above10pct = 0, above5pct = 0, above1pct = 0;
        double maxFreq = 0.0;
        
        for (const auto& [hash, count] : state.seedFreqInReads) {
            double frac = static_cast<double>(count) / totalReads;
            if (frac > maxFreq) maxFreq = frac;
            if (frac > 0.50) above50pct++;
            if (frac > 0.20) above20pct++;
            if (frac > 0.10) above10pct++;
            if (frac > 0.05) above5pct++;
            if (frac > 0.01) above1pct++;
        }
        
        // Log seed frequency analysis
        logging::info("=== Seed Frequency Analysis ===");
        logging::info("  Total reads: {}, Unique seeds: {}", totalReads, uniqueSeeds);
        logging::info("  Max seed frequency: {:.2f}%", maxFreq * 100.0);
        logging::info("  Seeds >50%% of reads: {}", above50pct);
        logging::info("  Seeds >20%% of reads: {}", above20pct);
        logging::info("  Seeds >10%% of reads: {}", above10pct);
        logging::info("  Seeds >5%% of reads: {}", above5pct);
        logging::info("  Seeds >1%% of reads: {}", above1pct);
        logging::info("===============================");
        
        // Only sort if we need to mask seeds or output diagnostics (slow for large datasets)
        std::vector<std::pair<size_t, int64_t>> sortedSeeds;
        bool needSorting = (seedMaskFraction > 0.0) || params.store_diagnostics;
        
        if (needSorting) {
            sortedSeeds.reserve(uniqueSeeds);
            for (const auto& [hash, count] : state.seedFreqInReads) {
                sortedSeeds.emplace_back(hash, count);
            }
            std::sort(sortedSeeds.begin(), sortedSeeds.end(), 
                      [](const auto& a, const auto& b) { return a.second > b.second; });
            
            // Show top 5 most frequent seeds if any are above 5%
            if (above5pct > 0) {
                logging::info("  Top {} high-frequency seeds:", std::min(size_t(5), above5pct));
                for (size_t i = 0; i < std::min(size_t(5), sortedSeeds.size()); i++) {
                    double frac = static_cast<double>(sortedSeeds[i].second) / totalReads;
                    if (frac < 0.01) break;  // Stop if below 1%
                    logging::info("    {:.2f}%  (hash: {:016x})", frac * 100.0, sortedSeeds[i].first);
                }
            }
        }
        
        // Apply percentile-based masking: remove top N% of seeds by frequency
        // This removes potential primer/adapter artifacts without needing to calibrate a threshold
        if (seedMaskFraction > 0.0 && !sortedSeeds.empty()) {
            const size_t numToMask = static_cast<size_t>(seedMaskFraction * uniqueSeeds);
            
            if (numToMask > 0) {
                int64_t maskedFrequency = 0;
                double minMaskedFrac = 0.0, maxMaskedFrac = 0.0;
                
                for (size_t i = 0; i < numToMask && i < sortedSeeds.size(); i++) {
                    state.seedFreqInReads.erase(sortedSeeds[i].first);
                    maskedFrequency += sortedSeeds[i].second;
                    double frac = static_cast<double>(sortedSeeds[i].second) / totalReads;
                    if (i == 0) maxMaskedFrac = frac;
                    minMaskedFrac = frac;
                }
                
                logging::info("  MASKING: Removed top {} seeds ({:.2f}%% of unique seeds)", 
                            numToMask, seedMaskFraction * 100.0);
                logging::info("  Masked seed freq range: {:.2f}%% - {:.2f}%% of reads", 
                            minMaskedFrac * 100.0, maxMaskedFrac * 100.0);
                logging::info("  Total masked frequency: {}", maskedFrequency);
            } else {
                logging::info("  Seed masking: fraction too small to mask any seeds");
            }
        }
        logging::info("===============================");
        
        // Dump seed frequencies to TSV file for analysis (only in debug mode - slow for large datasets)
        // Format: hash, sequence, count, fraction, masked (0/1)
        // If hashToSequence is populated, include the k-mer sequence
        if (!outputPath.empty() && params.store_diagnostics) {
            std::string seedFreqPath = outputPath + ".seed_freq.tsv";
            std::ofstream seedFreqFile(seedFreqPath);
            if (seedFreqFile.is_open()) {
                bool haveSequences = !state.hashToSequence.empty();
                if (haveSequences) {
                    seedFreqFile << "hash\tsequence\tcount\tfraction\tmasked\n";
                } else {
                    seedFreqFile << "hash\tcount\tfraction\tmasked\n";
                }
                size_t numToMask = static_cast<size_t>(seedMaskFraction * uniqueSeeds);
                for (size_t i = 0; i < sortedSeeds.size(); i++) {
                    double frac = static_cast<double>(sortedSeeds[i].second) / totalReads;
                    int masked = (i < numToMask) ? 1 : 0;
                    seedFreqFile << std::hex << sortedSeeds[i].first << std::dec << "\t";
                    if (haveSequences) {
                        auto it = state.hashToSequence.find(sortedSeeds[i].first);
                        if (it != state.hashToSequence.end()) {
                            seedFreqFile << it->second << "\t";
                        } else {
                            seedFreqFile << ".\t";
                        }
                    }
                    seedFreqFile << sortedSeeds[i].second << "\t"
                                 << std::fixed << std::setprecision(6) << frac << "\t"
                                 << masked << "\n";
                }
                seedFreqFile.close();
                logging::info("Seed frequencies written to: {}", seedFreqPath);
            }
        }
    }
    
    auto time_magnitude_start = std::chrono::high_resolution_clock::now();
    {
        // Compute log-scaled and double-log-scaled magnitudes for all metrics
        state.totalReadSeedFrequency = 0;
        
        // Pre-size maps for efficiency
        state.logReadCounts.reserve(state.seedFreqInReads.size());
        state.logLogReadCounts.reserve(state.seedFreqInReads.size());
        
        // Compute log-scaled values (needed for LogRaw and LogCosine)
        double logMagSquared = 0.0;
        double logLogMagSquared = 0.0;
        double rawMagSquared = 0.0;
        double cappedMagSquared = 0.0;
        double cappedLogMagSquared = 0.0;
        size_t filteredSeedCount = 0;
        size_t lowSupportSeeds = 0;
        
        // IMPROVED: Apply minimum read support filter for noise reduction
        // Seeds seen in only 1 read are more likely to be errors or random matches
        const int64_t minSupport = params.minReadSupport;
        
        for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
            state.totalReadSeedFrequency += readCount;
            
            // Filter seeds below minimum support threshold
            if (readCount < minSupport) {
                lowSupportSeeds++;
                continue;  // Skip this seed for scoring purposes
            }
            
            // Pre-compute raw magnitude: Σ(readCount²)
            const double readCountD = static_cast<double>(readCount);
            rawMagSquared += readCountD * readCountD;
            
            // Pre-compute log(1 + readCount) for each seed - avoids repeated log() in hot path
            double logCount = std::log1p(readCountD);
            state.logReadCounts[seedHash] = logCount;
            logMagSquared += logCount * logCount;
            
            // Pre-compute log(1 + log(1 + readCount)) for double-log metrics
            double logLogCount = std::log1p(logCount);
            state.logLogReadCounts[seedHash] = logLogCount;
            logLogMagSquared += logLogCount * logLogCount;
            
            // Pre-compute capped values: cap read count at 100
            const int64_t cappedCount = std::min(readCount, int64_t(100));
            const double cappedCountD = static_cast<double>(cappedCount);
            state.cappedReadCounts[seedHash] = cappedCount;
            cappedMagSquared += cappedCountD * cappedCountD;
            
            const double cappedLogCount = std::log1p(cappedCountD);
            state.cappedLogReadCounts[seedHash] = cappedLogCount;
            cappedLogMagSquared += cappedLogCount * cappedLogCount;
            
            filteredSeedCount++;
        }
        
        state.readUniqueSeedCount = filteredSeedCount;
        state.readMagnitude = std::sqrt(rawMagSquared);
        state.logReadMagnitude = std::sqrt(logMagSquared);
        state.logLogReadMagnitude = std::sqrt(logLogMagSquared);
        state.cappedReadMagnitude = std::sqrt(cappedMagSquared);
        state.cappedLogReadMagnitude = std::sqrt(cappedLogMagSquared);
        
        if (lowSupportSeeds > 0) {
            logging::info("Filtered {} seeds with read count < {} (kept {} seeds)", 
                         lowSupportSeeds, minSupport, filteredSeedCount);
        }
        logging::info("Precomputed magnitudes: log={:.6f}, loglog={:.6f}, total read seed frequency: {}, unique read seeds: {}", 
                    state.logReadMagnitude, state.logLogReadMagnitude, state.totalReadSeedFrequency, state.readUniqueSeedCount);
    }
    auto time_magnitude_end = std::chrono::high_resolution_clock::now();
    auto duration_magnitude = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_magnitude_end - time_magnitude_start);
    logging::info("Read magnitude precomputation: {}ms", duration_magnitude.count());
    
    state.root = liteTree->root;

    placement::NodeMetrics rootMetrics;
    // CRITICAL: Root starts with empty genome (all genome metrics = 0).
    // The root's seed changes will apply deltas to compute the actual root state.

    logging::info("Starting BFS placement traversal with {} tree nodes",
                  liteTree->allLiteNodes.size());

    auto time_traversal_start = std::chrono::high_resolution_clock::now();
    if (liteTree && liteTree->root) {
        std::vector<panmapUtils::LiteNode*> root_level_nodes;
        std::vector<placement::NodeMetrics> root_level_metrics;
        std::vector<absl::flat_hash_map<uint64_t, int64_t>> root_level_seed_counts;

        // Start with the root node itself (not just its children)
        root_level_nodes.push_back(liteTree->root);
        root_level_metrics.push_back(rootMetrics);
        if (params.verify_scores) {
            root_level_seed_counts.emplace_back();
        }
        
        std::atomic<size_t> nodesProcessed(0);
        placeLiteHelperBFS(root_level_nodes, root_level_metrics, root_level_seed_counts, state, result, params, nodesProcessed);

        // Lazy ID resolution: Only deserialize winning node IDs (eliminates ~1800ms overhead!)
        result.resolveNodeIds(liteTree);
        
        // =======================================================================
        // ALIGNMENT-BASED REFINEMENT (optional)
        // After k-mer scoring, refine top candidates by full minimap2 alignment
        // =======================================================================
        if (params.refineEnabled && state.fullTree != nullptr) {
            TraversalParams refineParams = params;  // Copy params for refinement
            refineTopCandidates(liteTree, state.fullTree, allReadSequences, result, refineParams);
            
            // Re-resolve node IDs to include refinement result
            result.resolveNodeIds(liteTree);
        } else if (params.refineEnabled && state.fullTree == nullptr) {
            logging::warn("Refinement enabled but fullTree is null - skipping alignment refinement");
        }
        
    } else {
        logging::err("LiteTree or root is null, cannot perform placement");
        return;
    }
    auto time_traversal_end = std::chrono::high_resolution_clock::now();
    auto duration_traversal = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_traversal_end - time_traversal_start);
    logging::info("Tree traversal: {}ms", duration_traversal.count());
    
    result.totalReadsProcessed = allReadSequences.size();
    
    // ==========================================================================
    // Store minimal alignment data (paths + seed frequencies only)
    // Full read data will be re-extracted during alignment to save memory
    // ==========================================================================
    {
        result.reads1Path = reads1;
        result.reads2Path = reads2;
        result.seedFreqInReads = std::move(state.seedFreqInReads);
        result.k = params.k;
        result.s = params.s;
        result.t = params.t;
        result.open = params.open;
        
        // Store read statistics for diagnostic output
        result.readUniqueSeedCount = state.readUniqueSeedCount;
        result.totalReadSeedFrequency = state.totalReadSeedFrequency;
        result.readMagnitude = state.logReadMagnitude;
        
        logging::info("Stored alignment parameters: k={}, s={}, t={}, open={}", 
                     result.k, result.s, result.t, result.open);
    }
    
    auto time_write_start = std::chrono::high_resolution_clock::now();
    std::string placementsFilePath = outputPath;
    std::ofstream placementsFile(placementsFilePath);
    if (placementsFile.is_open()) {
        placementsFile << "metric\tscore\tnodes\n";
        
        // ========================================
        // OLD METRICS (from working version 6cfb6cf)
        // ========================================
        
        // Raw score
        placementsFile << "raw\t" << std::fixed << std::setprecision(6) << result.bestRawScore << "\t";
        if (!result.tiedRawNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedRawNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedRawNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestRawNodeId;
        }
        placementsFile << "\n";
        
        // Cosine score (standard)
        placementsFile << "cosine\t" << std::fixed << std::setprecision(6) << result.bestCosineScore << "\t";
        if (!result.tiedCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestCosineNodeId;
        }
        placementsFile << "\n";
        
        // Log Cosine OLD (the working formula - log on reads only)
        placementsFile << "log_cosine_v1\t" << std::fixed << std::setprecision(6) << result.bestLogCosineOldScore << "\t";
        if (!result.tiedLogCosineOldNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogCosineOldNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogCosineOldNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogCosineOldNodeId;
        }
        placementsFile << "\n";
        
        // Weighted Containment score
        placementsFile << "weighted_containment\t" << std::fixed << std::setprecision(6) << result.bestWeightedContainmentScore << "\t";
        if (!result.tiedWeightedContainmentNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedWeightedContainmentNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedWeightedContainmentNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestWeightedContainmentNodeId;
        }
        placementsFile << "\n";
        
        // Cap Cosine score
        placementsFile << "cap_cosine\t" << std::fixed << std::setprecision(6) << result.bestCapCosineScore << "\t";
        if (!result.tiedCapCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedCapCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedCapCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestCapCosineNodeId;
        }
        placementsFile << "\n";
        
        // Cap Log Cosine score
        placementsFile << "cap_log_cosine\t" << std::fixed << std::setprecision(6) << result.bestCapLogCosineScore << "\t";
        if (!result.tiedCapLogCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedCapLogCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedCapLogCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestCapLogCosineNodeId;
        }
        placementsFile << "\n";
        
        // ========================================
        // NEW METRICS (current version)
        // ========================================
        
        // Log Raw score
        placementsFile << "log_raw\t" << std::fixed << std::setprecision(6) << result.bestLogRawScore << "\t";
        if (!result.tiedLogRawNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogRawNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogRawNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogRawNodeId;
        }
        placementsFile << "\n";
        
        // Log Cosine NEW (IDF-weighted with coverage penalty)
        placementsFile << "log_cosine_v2\t" << std::fixed << std::setprecision(6) << result.bestLogCosineScore << "\t";
        if (!result.tiedLogCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogCosineNodeId;
        }
        placementsFile << "\n";
        
        // LogLog Raw score
        placementsFile << "loglog_raw\t" << std::fixed << std::setprecision(6) << result.bestLogLogRawScore << "\t";
        if (!result.tiedLogLogRawNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogLogRawNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogLogRawNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogLogRawNodeId;
        }
        placementsFile << "\n";
        
        // LogLog Cosine score
        placementsFile << "loglog_cosine\t" << std::fixed << std::setprecision(6) << result.bestLogLogCosineScore << "\t";
        if (!result.tiedLogLogCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogLogCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogLogCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogLogCosineNodeId;
        }
        placementsFile << "\n";
        
        // Containment score
        placementsFile << "containment\t" << std::fixed << std::setprecision(6) << result.bestContainmentScore << "\t";
        if (!result.tiedContainmentNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedContainmentNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedContainmentNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestContainmentNodeId;
        }
        placementsFile << "\n";
        
        // Concordance score
        placementsFile << "concordance\t" << std::fixed << std::setprecision(6) << result.bestConcordanceScore << "\t";
        if (!result.tiedConcordanceNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedConcordanceNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedConcordanceNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestConcordanceNodeId;
        }
        placementsFile << "\n";
        
        // Presence score (low-coverage robust)
        placementsFile << "presence\t" << std::fixed << std::setprecision(6) << result.bestPresenceScore << "\t";
        if (!result.tiedPresenceNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedPresenceNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedPresenceNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestPresenceNodeId;
        }
        placementsFile << "\n";
        
        // Refinement score (if refinement was run)
        if (result.refinementWasRun) {
            placementsFile << "refined\t" << std::fixed << std::setprecision(0) << result.bestRefinedScore << "\t";
            placementsFile << result.bestRefinedNodeId;
            placementsFile << "\n";
        }
        
        placementsFile.close();;
        logging::info("Wrote placement results to {}", placementsFilePath);
    } else {
        logging::err("Failed to open placements file: {}", placementsFilePath);
    }
    auto time_write_end = std::chrono::high_resolution_clock::now();
    auto duration_write = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_write_end - time_write_start);
    logging::info("Result file writing: {}ms", duration_write.count());
    
    auto placement_total_end = std::chrono::high_resolution_clock::now();
    auto duration_placement_total = std::chrono::duration_cast<std::chrono::milliseconds>(
        placement_total_end - placement_total_start);
    
    logging::info("\n=== PLACEMENT INTERNAL TIMING BREAKDOWN ===");
    logging::info("HashDeltaLoading: {}ms", duration_hash_delta.count());
    logging::info("ReadProcessing: {}ms", duration_read_processing.count());
    logging::info("ReadMagnitude: {}ms", duration_magnitude.count());
    logging::info("TreeTraversal: {}ms", duration_traversal.count());
    logging::info("FileWriting: {}ms", duration_write.count());
    logging::info("PlacementTotal: {}ms", duration_placement_total.count());
    logging::info("==========================================\n");
    
    logging::info("Best LogRaw score: {:.6f} (node: {})", 
                 result.bestLogRawScore, result.bestLogRawNodeId);
    if (result.tiedLogRawNodeIndices.size() > 1) {
        logging::info("  {} nodes tied for best LogRaw score", result.tiedLogRawNodeIndices.size());
    }
    
    logging::info("Best LogCosine score: {:.6f} (node: {})", 
                 result.bestLogCosineScore, result.bestLogCosineNodeId);
    if (result.tiedLogCosineNodeIndices.size() > 1) {
        logging::info("  {} nodes tied for best LogCosine score", result.tiedLogCosineNodeIndices.size());
    }
    
    if (result.refinementWasRun) {
        logging::info("Best Refined score: {:.0f} (node: {})", 
                     result.bestRefinedScore, result.bestRefinedNodeId);
    }
}

} // namespace placement


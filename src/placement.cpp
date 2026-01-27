#include "panman.hpp"
#include "placement.hpp"
#include "seeding.hpp"
#include "index_lite.capnp.h"
#include "logging.hpp"


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
    std::cerr << "Error: File " << readPath1 << " not found" << std::endl;
    exit(0);
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
      std::cerr << "Error: File " << readPath2 << " not found" << std::endl;
      exit(0);
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
      std::cerr << "Error: File " << readPath2 << " does not contain the same number of reads as " << readPath1 << std::endl;
      exit(0);
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
    std::cerr << "Error: File " << readPath1 << " not found" << std::endl;
    exit(0);
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
      std::cerr << "Error: File " << readPath2 << " not found" << std::endl;
      exit(0);
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
      std::cerr << "Error: File " << readPath2 << " does not contain the same number of reads as " << readPath1 << std::endl;
      exit(0);
    }
    
    // Shuffle reads together, so that pairs are next to each other
    seeding::perfect_shuffle(readSequences);
    seeding::perfect_shuffle(readQuals);
    seeding::perfect_shuffle(readNames);
  }
}

// Build paired-end fragment tracking data
// After reads are shuffled, R1 is at even indices, R2 at odd indices
// Fragment i corresponds to reads 2i (R1) and 2i+1 (R2)
void buildPairedEndFragmentData(
    placement::PlacementGlobalState& state,
    const std::vector<std::string>& allReadSequences,
    const placement::TraversalParams& params) {
    
    // Check if paired-end (must have even number of reads, >= 2)
    if (allReadSequences.size() < 2 || allReadSequences.size() % 2 != 0) {
        state.isPairedEnd = false;
        state.numFragments = allReadSequences.size();
        logging::info("Single-end mode: {} reads (paired-end fragment tracking disabled)", 
                      allReadSequences.size());
        return;
    }
    
    state.isPairedEnd = true;
    state.numFragments = allReadSequences.size() / 2;
    
    logging::info("Paired-end mode: {} fragments (building R1/R2 seed sets for concordance scoring)", 
                  state.numFragments);
    
    // Allocate fragment seed sets
    state.fragmentR1Seeds.resize(state.numFragments);
    state.fragmentR2Seeds.resize(state.numFragments);
    state.fragmentReadLengths.resize(state.numFragments);
    
    auto time_fragment_build_start = std::chrono::high_resolution_clock::now();
    
    // Process fragments in parallel
    size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    
    // Thread-local data for seedToFragments
    std::vector<std::vector<std::pair<uint64_t, uint32_t>>> threadLocalSeedFragments(num_cpus);
    
    tbb::parallel_for(tbb::blocked_range<size_t>(0, state.numFragments), 
        [&](const tbb::blocked_range<size_t>& range) {
            size_t threadId = tbb::this_task_arena::current_thread_index();
            auto& localSeedFragments = threadLocalSeedFragments[threadId];
            
            for (size_t fragId = range.begin(); fragId < range.end(); ++fragId) {
                size_t r1Idx = fragId * 2;
                size_t r2Idx = fragId * 2 + 1;
                
                const std::string& r1Seq = allReadSequences[r1Idx];
                const std::string& r2Seq = allReadSequences[r2Idx];
                
                // Store read lengths for fragment size estimation
                state.fragmentReadLengths[fragId] = {
                    static_cast<uint16_t>(r1Seq.size()),
                    static_cast<uint16_t>(r2Seq.size())
                };
                
                // Calculate trim boundaries
                const int trimStartBp = params.trimStart;
                const int trimEndBp = params.trimEnd;
                
                // Extract R1 seeds
                auto& r1Seeds = state.fragmentR1Seeds[fragId];
                {
                    const int seqLen = static_cast<int>(r1Seq.size());
                    const int validStart = trimStartBp;
                    const int validEnd = seqLen - trimEndBp - params.k;
                    
                    const auto& syncmers = seeding::rollingSyncmers(r1Seq, params.k, params.s, params.open, params.t, false);
                    for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                        if (!isSyncmer) continue;
                        if (static_cast<int>(startPos) < validStart || static_cast<int>(startPos) > validEnd) continue;
                        r1Seeds.insert(kmerHash);
                        localSeedFragments.emplace_back(kmerHash, static_cast<uint32_t>(fragId));
                    }
                }
                
                // Extract R2 seeds
                auto& r2Seeds = state.fragmentR2Seeds[fragId];
                {
                    const int seqLen = static_cast<int>(r2Seq.size());
                    const int validStart = trimStartBp;
                    const int validEnd = seqLen - trimEndBp - params.k;
                    
                    const auto& syncmers = seeding::rollingSyncmers(r2Seq, params.k, params.s, params.open, params.t, false);
                    for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                        if (!isSyncmer) continue;
                        if (static_cast<int>(startPos) < validStart || static_cast<int>(startPos) > validEnd) continue;
                        r2Seeds.insert(kmerHash);
                        localSeedFragments.emplace_back(kmerHash, static_cast<uint32_t>(fragId));
                    }
                }
            }
        });
    
    // Merge thread-local seedToFragments into global map
    size_t totalPairs = 0;
    for (const auto& local : threadLocalSeedFragments) {
        totalPairs += local.size();
    }
    state.seedToFragments.reserve(totalPairs / 2);  // Estimate unique seeds
    
    for (const auto& local : threadLocalSeedFragments) {
        for (const auto& [hash, fragId] : local) {
            state.seedToFragments[hash].push_back(fragId);
        }
    }
    
    auto time_fragment_build_end = std::chrono::high_resolution_clock::now();
    auto duration_fragment_build = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_fragment_build_end - time_fragment_build_start);
    
    // Compute statistics
    size_t totalR1Seeds = 0, totalR2Seeds = 0, overlappingFragments = 0;
    for (size_t i = 0; i < state.numFragments; ++i) {
        totalR1Seeds += state.fragmentR1Seeds[i].size();
        totalR2Seeds += state.fragmentR2Seeds[i].size();
        
        // Check for overlap (shared seeds between R1 and R2)
        for (const auto& seed : state.fragmentR1Seeds[i]) {
            if (state.fragmentR2Seeds[i].contains(seed)) {
                ++overlappingFragments;
                break;  // Just count once per fragment
            }
        }
    }
    
    logging::info("[PLACEMENT TIMING]     Paired-end fragment build: {}ms", duration_fragment_build.count());
    logging::info("  Fragment stats: {} R1 seeds, {} R2 seeds, {}/{} fragments with R1/R2 overlap ({:.1f}%)",
                  totalR1Seeds, totalR2Seeds, 
                  overlappingFragments, state.numFragments,
                  100.0 * overlappingFragments / state.numFragments);
    logging::info("  Seed-to-fragment map: {} unique seeds tracked", state.seedToFragments.size());
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
            
            // Genome unique seed count delta (branchless)
            const int64_t wasPresent = (parentCount > 0) ? 1 : 0;
            const int64_t isPresent = (childCount > 0) ? 1 : 0;
            childMetrics.genomeUniqueSeedCount += (isPresent - wasPresent);
            
            // ========================================
            // READ-GENOME INTERACTION METRICS (only for seeds in reads)
            // ========================================
            
            // Fast path: skip unchanged seeds (23% of all seeds per profiler)
            if (freqDelta == 0) [[likely]] continue;
            
            // Single hash lookup: find() returns iterator, check if valid
            auto readIt = state.seedFreqInReads.find(seedHash);
            if (readIt == state.seedFreqInReads.end()) [[likely]] continue;
            
            // Hot path: seed changed AND present in reads (29% of all seeds)
            const int64_t readCount = readIt->second;
            
            // If parentCount == 0: became present, add 1
            // If childCount == 0: became absent, subtract 1
            // Otherwise: no change
            const int64_t becamePresent = (parentCount == 0) & (childCount != 0);
            const int64_t becameAbsent = (childCount == 0) & (parentCount != 0);
            const int64_t presenceDelta = becamePresent - becameAbsent;
            
            childMetrics.presenceIntersectionCount += presenceDelta;
            childMetrics.jaccardNumerator += presenceDelta;
            
            // Update raw match score: Σ(readCount / genomeCount) for seeds in intersection
            // Delta = readCount/childCount - readCount/parentCount (treating 1/0 as 0)
            const double oldContrib = (parentCount > 0) ? (static_cast<double>(readCount) / parentCount) : 0.0;
            const double newContrib = (childCount > 0) ? (static_cast<double>(readCount) / childCount) : 0.0;
            childMetrics.rawMatchScore += (newContrib - oldContrib);

            // Update logRaw numerator: Σ(log(1+readCount) / genomeCount) for seeds in intersection
            // Log-scales read counts to reduce impact of variable coverage
            // Look up pre-computed log value from state
            auto logIt = state.logReadCounts.find(seedHash);
            if (logIt != state.logReadCounts.end()) {
                const double logReadCount = logIt->second;
                const double oldLogContrib = (parentCount > 0) ? (logReadCount / parentCount) : 0.0;
                const double newLogContrib = (childCount > 0) ? (logReadCount / childCount) : 0.0;
                childMetrics.logRawNumerator += (newLogContrib - oldLogContrib);
                
                // Update logCosine numerator: Σ(log(1+readCount) × genomeCount) delta
                childMetrics.logCosineNumerator += logReadCount * freqDelta;
            }

            // Update weighted Jaccard numerator: delta of min(read, genome)
            // Using branchless min (compiles to CMOV on x86-64)
            const int64_t oldMin = std::min(readCount, parentCount);
            const int64_t newMin = std::min(readCount, childCount);
            childMetrics.weightedJaccardNumerator += (newMin - oldMin);

            // Update cosine numerator: Σ(r_i * g_i) delta
            childMetrics.cosineNumerator += static_cast<double>(readCount) * freqDelta;

            // Update IDF-weighted cosine numerator: Σ(r_i * g_i * IDF(seed)) delta
            // IDF weight upweights rare seeds that appear in fewer genomes
            auto idfIt = state.seedIdfWeights.find(seedHash);
            if (idfIt != state.seedIdfWeights.end()) {
                const double idfWeight = idfIt->second;
                childMetrics.idfCosineNumerator += static_cast<double>(readCount) * freqDelta * idfWeight;
            } else {
                // Seed not in IDF map (appears in all genomes?), use weight 1.0
                childMetrics.idfCosineNumerator += static_cast<double>(readCount) * freqDelta;
            }

            // Update coverage-weighted cosine numerator: Σ(covWeight_i * g_i) delta
            // Uses pre-computed coverage-weighted counts that penalize singletons and over-rep seeds
            auto covIt = state.coverageWeightedCounts.find(seedHash);
            if (covIt != state.coverageWeightedCounts.end()) {
                childMetrics.covCosineNumerator += covIt->second * freqDelta;
            } else {
                // Fallback: use raw read count (shouldn't happen if maps are consistent)
                childMetrics.covCosineNumerator += static_cast<double>(readCount) * freqDelta;
            }

            // Update capped cosine numerator: Σ(min(readCount,100) × genomeCount) delta
            // Caps read frequency at 100 to limit influence of high-frequency seeds
            auto capIt = state.cappedReadCounts.find(seedHash);
            if (capIt != state.cappedReadCounts.end()) {
                childMetrics.capCosineNumerator += static_cast<double>(capIt->second) * freqDelta;
            } else {
                // Fallback: use capped raw read count
                childMetrics.capCosineNumerator += static_cast<double>(std::min(readCount, int64_t(100))) * freqDelta;
            }

            // Update capped log cosine numerator: Σ(log(1+min(readCount,100)) × genomeCount) delta
            auto capLogIt = state.cappedLogReadCounts.find(seedHash);
            if (capLogIt != state.cappedLogReadCounts.end()) {
                childMetrics.capLogCosineNumerator += capLogIt->second * freqDelta;
            } else {
                // Fallback: compute capped log value
                childMetrics.capLogCosineNumerator += std::log1p(static_cast<double>(std::min(readCount, int64_t(100)))) * freqDelta;
            }

            // Update sigmoid cosine numerator: Σ(sigmoidScaled(readCount) × genomeCount) delta
            // Uses pre-computed sigmoid-scaled counts centered on median seed frequency
            auto sigIt = state.sigmoidReadCounts.find(seedHash);
            if (sigIt != state.sigmoidReadCounts.end()) {
                childMetrics.sigCosineNumerator += sigIt->second * freqDelta;
            } else {
                // Fallback: use raw read count (shouldn't happen if maps are consistent)
                childMetrics.sigCosineNumerator += static_cast<double>(readCount) * freqDelta;
            }
        }
    }
    
}


namespace placement {

using ::panmanUtils::Node;
using ::panmanUtils::Tree;

// Forward declarations
class PlacementResult;

// Macro to generate score update functions - reduces ~370 lines of repetitive code to ~30 lines
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

// Generate all 13 score update functions using the macro
DEFINE_UPDATE_SCORE_FUNC(updateRawSeedMatchScore, rawSeedMatchScore, bestRawSeedMatchScore, bestRawSeedMatchNodeIndex, tiedRawSeedMatchNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCosineScore, cosineScore, bestCosineScore, bestCosineNodeIndex, tiedCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateContainmentScore, containmentScore, bestContainmentScore, bestContainmentNodeIndex, tiedContainmentNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateWeightedContainmentScore, weightedContainmentScore, bestWeightedContainmentScore, bestWeightedContainmentNodeIndex, tiedWeightedContainmentNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogRawScore, logRawScore, bestLogRawScore, bestLogRawNodeIndex, tiedLogRawNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogCosineScore, logCosineScore, bestLogCosineScore, bestLogCosineNodeIndex, tiedLogCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateAdjCosineScore, adjCosineScore, bestAdjCosineScore, bestAdjCosineNodeIndex, tiedAdjCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateAdjRawScore, adjRawScore, bestAdjRawScore, bestAdjRawNodeIndex, tiedAdjRawNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateIdfCosineScore, idfCosineScore, bestIdfCosineScore, bestIdfCosineNodeIndex, tiedIdfCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCovCosineScore, covCosineScore, bestCovCosineScore, bestCovCosineNodeIndex, tiedCovCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCapCosineScore, capCosineScore, bestCapCosineScore, bestCapCosineNodeIndex, tiedCapCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateCapLogCosineScore, capLogCosineScore, bestCapLogCosineScore, bestCapLogCosineNodeIndex, tiedCapLogCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateSigCosineScore, sigCosineScore, bestSigCosineScore, bestSigCosineNodeIndex, tiedSigCosineNodeIndices)

#undef DEFINE_UPDATE_SCORE_FUNC

void PlacementResult::resolveNodeIds(panmapUtils::LiteTree* liteTree) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (!liteTree) return;
    
    // Resolve best node indices to string IDs
    if (bestJaccardNodeIndex != UINT32_MAX) {
        bestJaccardNodeId = liteTree->resolveNodeId(bestJaccardNodeIndex);
    }
    if (bestWeightedJaccardNodeIndex != UINT32_MAX) {
        bestWeightedJaccardNodeId = liteTree->resolveNodeId(bestWeightedJaccardNodeIndex);
    }
    if (bestCosineNodeIndex != UINT32_MAX) {
        bestCosineNodeId = liteTree->resolveNodeId(bestCosineNodeIndex);
    }
    if (bestJaccardPresenceNodeIndex != UINT32_MAX) {
        bestJaccardPresenceNodeId = liteTree->resolveNodeId(bestJaccardPresenceNodeIndex);
    }
    if (bestRawSeedMatchNodeIndex != UINT32_MAX) {
        bestRawSeedMatchNodeId = liteTree->resolveNodeId(bestRawSeedMatchNodeIndex);
    }
    if (maxHitsNodeIndex != UINT32_MAX) {
        maxHitsNodeId = liteTree->resolveNodeId(maxHitsNodeIndex);
    }
    if (bestContainmentNodeIndex != UINT32_MAX) {
        bestContainmentNodeId = liteTree->resolveNodeId(bestContainmentNodeIndex);
    }
    if (bestWeightedContainmentNodeIndex != UINT32_MAX) {
        bestWeightedContainmentNodeId = liteTree->resolveNodeId(bestWeightedContainmentNodeIndex);
    }
    if (bestLogRawNodeIndex != UINT32_MAX) {
        bestLogRawNodeId = liteTree->resolveNodeId(bestLogRawNodeIndex);
    }
    if (bestLogCosineNodeIndex != UINT32_MAX) {
        bestLogCosineNodeId = liteTree->resolveNodeId(bestLogCosineNodeIndex);
    }
    if (bestAdjCosineNodeIndex != UINT32_MAX) {
        bestAdjCosineNodeId = liteTree->resolveNodeId(bestAdjCosineNodeIndex);
    }
    if (bestAdjRawNodeIndex != UINT32_MAX) {
        bestAdjRawNodeId = liteTree->resolveNodeId(bestAdjRawNodeIndex);
    }
    if (bestIdfCosineNodeIndex != UINT32_MAX) {
        bestIdfCosineNodeId = liteTree->resolveNodeId(bestIdfCosineNodeIndex);
    }
    if (bestCovCosineNodeIndex != UINT32_MAX) {
        bestCovCosineNodeId = liteTree->resolveNodeId(bestCovCosineNodeIndex);
    }
    if (bestCapCosineNodeIndex != UINT32_MAX) {
        bestCapCosineNodeId = liteTree->resolveNodeId(bestCapCosineNodeIndex);
    }
    if (bestCapLogCosineNodeIndex != UINT32_MAX) {
        bestCapLogCosineNodeId = liteTree->resolveNodeId(bestCapLogCosineNodeIndex);
    }
    if (bestSigCosineNodeIndex != UINT32_MAX) {
        bestSigCosineNodeId = liteTree->resolveNodeId(bestSigCosineNodeIndex);
    }
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
                
                // NOTE: genomeMagnitudeSquared, genomeUniqueSeedCount, genomeTotalSeedFrequency
                // are computed incrementally in computeChildMetrics() - NO per-node storage needed!

                // Skip scoring for the excluded leaf node (leave-one-out validation)
                if (nodeIndex != state.skipNodeIndex) {
                    // Only compute RAW and Cosine - they correctly weight by read frequency
                    // Jaccard/WJ/Presence are disabled: they degenerate to counting when genome freq=1
                    double cosineScore = nodeMetrics.getCosineScore(state.readMagnitude);
                    
                    // Containment Index: |reads ∩ genome| / |reads|
                    // This measures what fraction of read seeds are covered by the genome
                    // Unlike Cosine, it doesn't reward larger genomes with more seeds
                    double containmentScore = nodeMetrics.getContainmentScore(state.readUniqueSeedCount);
                    double weightedContainmentScore = nodeMetrics.getWeightedContainmentScore(state.totalReadSeedFrequency);
                    
                    // LogRAW Score: log-scaled, coverage-robust metric
                    double logRawScore = nodeMetrics.getLogRawScore(state.logReadMagnitude);
                    
                    // LogCosine Score: log-scaled cosine similarity
                    double logCosineScore = nodeMetrics.getLogCosineScore(state.logReadMagnitude);
                    
                    // AdjCosine Score: N-adjusted cosine (penalizes significantly incomplete genomes)
                    // Only penalizes genomes below (mean - 2*std) threshold
                    double adjCosineScore = nodeMetrics.getAdjCosineScore(state.readMagnitude, state.incompletenessThreshold);
                    
                    // AdjRaw Score: N-adjusted RAW (penalizes significantly incomplete genomes)
                    double adjRawScore = nodeMetrics.getAdjRawScore(state.logReadMagnitude, state.incompletenessThreshold);
                    
                    // Update thread-local results (NO LOCKS - parallel performance!)
                    tls.local_result.updateCosineScore(nodeIndex, cosineScore, node);
                    // Raw score: Σ(readCount / genomeCount) - rewards unique hits
                    tls.local_result.updateRawSeedMatchScore(nodeIndex, nodeMetrics.rawMatchScore, node);
                    // Containment metrics
                    tls.local_result.updateContainmentScore(nodeIndex, containmentScore, node);
                    tls.local_result.updateWeightedContainmentScore(nodeIndex, weightedContainmentScore, node);
                    // LogRAW: log-scaled, coverage-robust
                    tls.local_result.updateLogRawScore(nodeIndex, logRawScore, node);
                    // LogCosine: log-scaled cosine
                    tls.local_result.updateLogCosineScore(nodeIndex, logCosineScore, node);
                    // AdjCosine: N-adjusted cosine
                    tls.local_result.updateAdjCosineScore(nodeIndex, adjCosineScore, node);
                    // AdjRaw: N-adjusted RAW
                    tls.local_result.updateAdjRawScore(nodeIndex, adjRawScore, node);
                    // IdfCosine: IDF-weighted cosine (upweights rare seeds)
                    double idfCosineScore = nodeMetrics.getIdfCosineScore(state.readMagnitude);
                    tls.local_result.updateIdfCosineScore(nodeIndex, idfCosineScore, node);
                    
                    // CovCosine: Coverage-weighted cosine (penalizes singletons & over-rep seeds)
                    double covCosineScore = nodeMetrics.getCovCosineScore(state.covWeightedMagnitude);
                    tls.local_result.updateCovCosineScore(nodeIndex, covCosineScore, node);
                    
                    // CapCosine: Capped cosine (caps read freqs at 100)
                    double capCosineScore = nodeMetrics.getCapCosineScore(state.cappedReadMagnitude);
                    tls.local_result.updateCapCosineScore(nodeIndex, capCosineScore, node);
                    
                    // CapLogCosine: Capped log cosine (caps read freqs at 100 before log)
                    double capLogCosineScore = nodeMetrics.getCapLogCosineScore(state.cappedLogReadMagnitude);
                    tls.local_result.updateCapLogCosineScore(nodeIndex, capLogCosineScore, node);
                    
                    // SigCosine: Sigmoid-scaled cosine (adaptive to median seed frequency)
                    double sigCosineScore = nodeMetrics.getSigCosineScore(state.sigmoidReadMagnitude);
                    tls.local_result.updateSigCosineScore(nodeIndex, sigCosineScore, node);
                    
                    // Store diagnostic data on node if enabled
                    if (params.store_diagnostics) {
                        node->jaccardNumerator = nodeMetrics.jaccardNumerator;
                        node->weightedJaccardNumerator = nodeMetrics.weightedJaccardNumerator;
                        node->cosineNumerator = nodeMetrics.cosineNumerator;
                        node->presenceIntersectionCount = nodeMetrics.presenceIntersectionCount;
                        node->genomeMagnitudeSquared = nodeMetrics.genomeMagnitudeSquared;
                        node->genomeUniqueSeedCount = nodeMetrics.genomeUniqueSeedCount;
                        node->genomeTotalSeedFrequency = nodeMetrics.genomeTotalSeedFrequency;
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
    // This is the KEY optimization: all threads computed scores independently (lock-free),
    // now we do a single-pass merge with minimal synchronization
    logging::info("Merging {} thread-local results...", thread_data.size());
    auto merge_start = std::chrono::high_resolution_clock::now();
    
    for (auto& tls : thread_data) {
        // Merge Jaccard scores
        if (tls.local_result.bestJaccardNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestJaccardScore * 0.0001, 1e-9);
            if (tls.local_result.bestJaccardScore > result.bestJaccardScore + tolerance) {
                result.bestJaccardScore = tls.local_result.bestJaccardScore;
                result.bestJaccardNodeIndex = tls.local_result.bestJaccardNodeIndex;
                result.tiedJaccardNodeIndices = tls.local_result.tiedJaccardNodeIndices;
            } else if (tls.local_result.bestJaccardScore >= result.bestJaccardScore - tolerance) {
                // Merge tied nodes
                result.tiedJaccardNodeIndices.insert(
                    result.tiedJaccardNodeIndices.end(),
                    tls.local_result.tiedJaccardNodeIndices.begin(),
                    tls.local_result.tiedJaccardNodeIndices.end()
                );
            }
        }
        
        // Merge Weighted Jaccard scores
        if (tls.local_result.bestWeightedJaccardNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestWeightedJaccardScore * 0.0001, 1e-9);
            if (tls.local_result.bestWeightedJaccardScore > result.bestWeightedJaccardScore + tolerance) {
                result.bestWeightedJaccardScore = tls.local_result.bestWeightedJaccardScore;
                result.bestWeightedJaccardNodeIndex = tls.local_result.bestWeightedJaccardNodeIndex;
                result.tiedWeightedJaccardNodeIndices = tls.local_result.tiedWeightedJaccardNodeIndices;
            } else if (tls.local_result.bestWeightedJaccardScore >= result.bestWeightedJaccardScore - tolerance) {
                result.tiedWeightedJaccardNodeIndices.insert(
                    result.tiedWeightedJaccardNodeIndices.end(),
                    tls.local_result.tiedWeightedJaccardNodeIndices.begin(),
                    tls.local_result.tiedWeightedJaccardNodeIndices.end()
                );
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
                result.tiedCosineNodeIndices.insert(
                    result.tiedCosineNodeIndices.end(),
                    tls.local_result.tiedCosineNodeIndices.begin(),
                    tls.local_result.tiedCosineNodeIndices.end()
                );
            }
        }
        
        // Merge Presence Jaccard scores
        if (tls.local_result.bestJaccardPresenceNodeIndex != UINT32_MAX) {
            if (tls.local_result.bestJaccardPresenceScore > result.bestJaccardPresenceScore) {
                result.bestJaccardPresenceScore = tls.local_result.bestJaccardPresenceScore;
                result.bestJaccardPresenceNodeIndex = tls.local_result.bestJaccardPresenceNodeIndex;
            }
        }
        
        // Merge Raw Seed Match scores
        if (tls.local_result.bestRawSeedMatchNodeIndex != UINT32_MAX) {
            if (tls.local_result.bestRawSeedMatchScore > result.bestRawSeedMatchScore) {
                result.bestRawSeedMatchScore = tls.local_result.bestRawSeedMatchScore;
                result.bestRawSeedMatchNodeIndex = tls.local_result.bestRawSeedMatchNodeIndex;
                result.tiedRawSeedMatchNodeIndices = tls.local_result.tiedRawSeedMatchNodeIndices;
            } else if (tls.local_result.bestRawSeedMatchScore == result.bestRawSeedMatchScore) {
                result.tiedRawSeedMatchNodeIndices.insert(
                    result.tiedRawSeedMatchNodeIndices.end(),
                    tls.local_result.tiedRawSeedMatchNodeIndices.begin(),
                    tls.local_result.tiedRawSeedMatchNodeIndices.end()
                );
            }
        }
        
        // Merge Hits scores
        if (tls.local_result.maxHitsNodeIndex != UINT32_MAX) {
            if (tls.local_result.maxHitsInAnyGenome > result.maxHitsInAnyGenome) {
                result.maxHitsInAnyGenome = tls.local_result.maxHitsInAnyGenome;
                result.maxHitsNodeIndex = tls.local_result.maxHitsNodeIndex;
                result.tiedMaxHitsNodeIndices = tls.local_result.tiedMaxHitsNodeIndices;
            } else if (tls.local_result.maxHitsInAnyGenome == result.maxHitsInAnyGenome) {
                result.tiedMaxHitsNodeIndices.insert(
                    result.tiedMaxHitsNodeIndices.end(),
                    tls.local_result.tiedMaxHitsNodeIndices.begin(),
                    tls.local_result.tiedMaxHitsNodeIndices.end()
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
        
        // Merge Weighted Containment scores
        if (tls.local_result.bestWeightedContainmentNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestWeightedContainmentScore * 0.0001, 1e-9);
            if (tls.local_result.bestWeightedContainmentScore > result.bestWeightedContainmentScore + tolerance) {
                result.bestWeightedContainmentScore = tls.local_result.bestWeightedContainmentScore;
                result.bestWeightedContainmentNodeIndex = tls.local_result.bestWeightedContainmentNodeIndex;
                result.tiedWeightedContainmentNodeIndices = tls.local_result.tiedWeightedContainmentNodeIndices;
            } else if (tls.local_result.bestWeightedContainmentScore >= result.bestWeightedContainmentScore - tolerance) {
                result.tiedWeightedContainmentNodeIndices.insert(
                    result.tiedWeightedContainmentNodeIndices.end(),
                    tls.local_result.tiedWeightedContainmentNodeIndices.begin(),
                    tls.local_result.tiedWeightedContainmentNodeIndices.end()
                );
            }
        }
        
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
        
        // Merge AdjCosine scores (N-adjusted cosine)
        if (tls.local_result.bestAdjCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestAdjCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestAdjCosineScore > result.bestAdjCosineScore + tolerance) {
                result.bestAdjCosineScore = tls.local_result.bestAdjCosineScore;
                result.bestAdjCosineNodeIndex = tls.local_result.bestAdjCosineNodeIndex;
                result.tiedAdjCosineNodeIndices = tls.local_result.tiedAdjCosineNodeIndices;
            } else if (tls.local_result.bestAdjCosineScore >= result.bestAdjCosineScore - tolerance) {
                result.tiedAdjCosineNodeIndices.insert(
                    result.tiedAdjCosineNodeIndices.end(),
                    tls.local_result.tiedAdjCosineNodeIndices.begin(),
                    tls.local_result.tiedAdjCosineNodeIndices.end()
                );
            }
        }
        
        // Merge AdjRaw scores (N-adjusted RAW)
        if (tls.local_result.bestAdjRawNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestAdjRawScore * 0.0001, 1e-9);
            if (tls.local_result.bestAdjRawScore > result.bestAdjRawScore + tolerance) {
                result.bestAdjRawScore = tls.local_result.bestAdjRawScore;
                result.bestAdjRawNodeIndex = tls.local_result.bestAdjRawNodeIndex;
                result.tiedAdjRawNodeIndices = tls.local_result.tiedAdjRawNodeIndices;
            } else if (tls.local_result.bestAdjRawScore >= result.bestAdjRawScore - tolerance) {
                result.tiedAdjRawNodeIndices.insert(
                    result.tiedAdjRawNodeIndices.end(),
                    tls.local_result.tiedAdjRawNodeIndices.begin(),
                    tls.local_result.tiedAdjRawNodeIndices.end()
                );
            }
        }
        
        // Merge IdfCosine scores (IDF-weighted cosine)
        if (tls.local_result.bestIdfCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestIdfCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestIdfCosineScore > result.bestIdfCosineScore + tolerance) {
                result.bestIdfCosineScore = tls.local_result.bestIdfCosineScore;
                result.bestIdfCosineNodeIndex = tls.local_result.bestIdfCosineNodeIndex;
                result.tiedIdfCosineNodeIndices = tls.local_result.tiedIdfCosineNodeIndices;
            } else if (tls.local_result.bestIdfCosineScore >= result.bestIdfCosineScore - tolerance) {
                result.tiedIdfCosineNodeIndices.insert(
                    result.tiedIdfCosineNodeIndices.end(),
                    tls.local_result.tiedIdfCosineNodeIndices.begin(),
                    tls.local_result.tiedIdfCosineNodeIndices.end()
                );
            }
        }
        
        // Merge CovCosine scores (coverage-weighted cosine)
        if (tls.local_result.bestCovCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestCovCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestCovCosineScore > result.bestCovCosineScore + tolerance) {
                result.bestCovCosineScore = tls.local_result.bestCovCosineScore;
                result.bestCovCosineNodeIndex = tls.local_result.bestCovCosineNodeIndex;
                result.tiedCovCosineNodeIndices = tls.local_result.tiedCovCosineNodeIndices;
            } else if (tls.local_result.bestCovCosineScore >= result.bestCovCosineScore - tolerance) {
                result.tiedCovCosineNodeIndices.insert(
                    result.tiedCovCosineNodeIndices.end(),
                    tls.local_result.tiedCovCosineNodeIndices.begin(),
                    tls.local_result.tiedCovCosineNodeIndices.end()
                );
            }
        }
        
        // Merge CapCosine scores (capped cosine - caps read freqs at 100)
        if (tls.local_result.bestCapCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestCapCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestCapCosineScore > result.bestCapCosineScore + tolerance) {
                result.bestCapCosineScore = tls.local_result.bestCapCosineScore;
                result.bestCapCosineNodeIndex = tls.local_result.bestCapCosineNodeIndex;
                result.tiedCapCosineNodeIndices = tls.local_result.tiedCapCosineNodeIndices;
            } else if (tls.local_result.bestCapCosineScore >= result.bestCapCosineScore - tolerance) {
                result.tiedCapCosineNodeIndices.insert(
                    result.tiedCapCosineNodeIndices.end(),
                    tls.local_result.tiedCapCosineNodeIndices.begin(),
                    tls.local_result.tiedCapCosineNodeIndices.end()
                );
            }
        }
        
        // Merge CapLogCosine scores (capped log cosine - caps read freqs at 100 before log)
        if (tls.local_result.bestCapLogCosineNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestCapLogCosineScore * 0.0001, 1e-9);
            if (tls.local_result.bestCapLogCosineScore > result.bestCapLogCosineScore + tolerance) {
                result.bestCapLogCosineScore = tls.local_result.bestCapLogCosineScore;
                result.bestCapLogCosineNodeIndex = tls.local_result.bestCapLogCosineNodeIndex;
                result.tiedCapLogCosineNodeIndices = tls.local_result.tiedCapLogCosineNodeIndices;
            } else if (tls.local_result.bestCapLogCosineScore >= result.bestCapLogCosineScore - tolerance) {
                result.tiedCapLogCosineNodeIndices.insert(
                    result.tiedCapLogCosineNodeIndices.end(),
                    tls.local_result.tiedCapLogCosineNodeIndices.begin(),
                    tls.local_result.tiedCapLogCosineNodeIndices.end()
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
                       uint32_t expectedFragmentSize,
                       uint32_t fragmentSizeTolerance) {
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
        
        // VERSION 3: Struct-of-arrays format - THREE separate primitive arrays
        // This is OPTIMAL for parallel access! Each thread reads contiguous memory.
        auto hashesReader = indexRoot.getSeedChangeHashes();
        auto parentCountsReader = indexRoot.getSeedChangeParentCounts();
        auto childCountsReader = indexRoot.getSeedChangeChildCounts();
        
        logging::info("Struct-of-arrays format: {} hashes, {} parent counts, {} child counts",
                        hashesReader.size(), parentCountsReader.size(), childCountsReader.size());
        
        // Phase 2a: Parallel read of THREE primitive arrays (optimal cache-friendly access!)
        // Each array is contiguous in memory, so threads can read sequential chunks without contention.
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
        
        // Phase 2b: Assign spans to nodes (embarrassingly parallel, no I/O)
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
    logging::info("[PLACEMENT TIMING]   Hash delta loading: {}ms", duration_hash_delta.count());
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
    
    // Load IDF weights for seed discrimination
    // Check if IDF data is present by checking list sizes (totalLeafGenomes is scalar, always present)
    if (indexRoot.hasIdfSeedHashes() && indexRoot.hasIdfGenomeCounts() && indexRoot.getTotalLeafGenomes() > 0) {
        state.totalLeafGenomes = indexRoot.getTotalLeafGenomes();
        auto idfHashes = indexRoot.getIdfSeedHashes();
        auto idfCounts = indexRoot.getIdfGenomeCounts();
        
        if (idfHashes.size() == idfCounts.size() && state.totalLeafGenomes > 0) {
            state.seedIdfWeights.reserve(idfHashes.size());
            double logN = std::log(static_cast<double>(state.totalLeafGenomes));
            
            for (size_t i = 0; i < idfHashes.size(); ++i) {
                // IDF = log(N / df) where df = number of genomes containing this seed
                // Add 1 to df to avoid log(inf) for seeds in all genomes
                double idf = logN - std::log(static_cast<double>(idfCounts[i]));
                state.seedIdfWeights[idfHashes[i]] = std::max(0.0, idf);  // Clamp to non-negative
            }
            
            logging::info("Loaded IDF weights for {} seeds from {} leaf genomes", 
                         state.seedIdfWeights.size(), state.totalLeafGenomes);
        }
    }

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
                logging::info("[PLACEMENT TIMING]     Seed extraction (Q{} filtered): {}ms", 
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
            logging::info("[PLACEMENT TIMING]     Read deduplication: {}ms", duration_dedup.count());
            
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
            logging::info("[PLACEMENT TIMING]     Seed extraction: {}ms", duration_seed_extract.count());
            
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
                logging::info("[PLACEMENT TIMING]     K-minimizer extraction (Q{} filtered): {}ms", 
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
            logging::info("[PLACEMENT TIMING]     K-minimizer extraction: {}ms", duration_kminimizer.count());
            logging::info("Extracted {} unique k-minimizers from {} reads", 
                        state.seedFreqInReads.size(), allReadSequences.size());
            }  // end of else block for non-quality-filtered path (l > 0)
        }
    }
    auto time_read_processing_end = std::chrono::high_resolution_clock::now();
    auto duration_read_processing = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_read_processing_end - time_read_processing_start);
    logging::info("[PLACEMENT TIMING]   Total read processing: {}ms", duration_read_processing.count());
    
    // ==========================================================================
    // Seed frequency analysis and masking
    // Identifies and optionally removes over-represented seeds (primer/adapter artifacts)
    // ==========================================================================
    if (!allReadSequences.empty()) {
        const size_t totalReads = allReadSequences.size();
        size_t uniqueSeeds = state.seedFreqInReads.size();
        
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
        
        // Compute frequency distribution statistics
        std::vector<std::pair<size_t, int64_t>> sortedSeeds;
        sortedSeeds.reserve(uniqueSeeds);
        for (const auto& [hash, count] : state.seedFreqInReads) {
            sortedSeeds.emplace_back(hash, count);
        }
        std::sort(sortedSeeds.begin(), sortedSeeds.end(), 
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Count seeds above various thresholds
        size_t above50pct = 0, above20pct = 0, above10pct = 0, above5pct = 0, above1pct = 0;
        double maxFreq = 0.0;
        
        for (const auto& [hash, count] : sortedSeeds) {
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
        
        // Show top 5 most frequent seeds if any are above 5%
        if (above5pct > 0 && !sortedSeeds.empty()) {
            logging::info("  Top {} high-frequency seeds:", std::min(size_t(5), above5pct));
            for (size_t i = 0; i < std::min(size_t(5), sortedSeeds.size()); i++) {
                double frac = static_cast<double>(sortedSeeds[i].second) / totalReads;
                if (frac < 0.01) break;  // Stop if below 1%
                logging::info("    {:.2f}%  (hash: {:016x})", frac * 100.0, sortedSeeds[i].first);
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
        
        // Dump seed frequencies to TSV file for analysis
        // Format: hash, count, fraction, masked (0/1)
        if (!outputPath.empty()) {
            std::string seedFreqPath = outputPath + ".seed_freq.tsv";
            std::ofstream seedFreqFile(seedFreqPath);
            if (seedFreqFile.is_open()) {
                seedFreqFile << "hash\tcount\tfraction\tmasked\n";
                size_t numToMask = static_cast<size_t>(seedMaskFraction * uniqueSeeds);
                for (size_t i = 0; i < sortedSeeds.size(); i++) {
                    double frac = static_cast<double>(sortedSeeds[i].second) / totalReads;
                    int masked = (i < numToMask) ? 1 : 0;
                    seedFreqFile << std::hex << sortedSeeds[i].first << std::dec << "\t"
                                 << sortedSeeds[i].second << "\t"
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

        state.totalReadSeedFrequency = 0;
        state.readUniqueSeedCount = state.seedFreqInReads.size();
        
        // Pre-size maps for efficiency
        state.logReadCounts.reserve(state.seedFreqInReads.size());
        state.cappedReadCounts.reserve(state.seedFreqInReads.size());
        state.cappedLogReadCounts.reserve(state.seedFreqInReads.size());
        
        // Compute read magnitude, total frequency, and log-scaled values
        // Note: seedFreqInReads map serves as both frequency storage and O(1) membership test
        double logMagSquared = 0.0;
        double cappedMagSquared = 0.0;
        double cappedLogMagSquared = 0.0;
        constexpr int64_t CAP_VALUE = 100;  // Cap read frequencies at 100
        
        for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
            state.readMagnitude += static_cast<double>(readCount * readCount);
            state.totalReadSeedFrequency += readCount;
            
            // Pre-compute log(1 + readCount) for each seed - avoids repeated log() in hot path
            double logCount = std::log1p(static_cast<double>(readCount));
            state.logReadCounts[seedHash] = logCount;
            logMagSquared += logCount * logCount;
            
            // Capped values: cap at 100 to limit influence of high-frequency seeds
            int64_t cappedCount = std::min(readCount, CAP_VALUE);
            state.cappedReadCounts[seedHash] = cappedCount;
            cappedMagSquared += static_cast<double>(cappedCount * cappedCount);
            
            // Capped log: log(1 + min(readCount, 100))
            double cappedLogCount = std::log1p(static_cast<double>(cappedCount));
            state.cappedLogReadCounts[seedHash] = cappedLogCount;
            cappedLogMagSquared += cappedLogCount * cappedLogCount;
        }
        state.readMagnitude = std::sqrt(state.readMagnitude);
        state.logReadMagnitude = std::sqrt(logMagSquared);
        state.cappedReadMagnitude = std::sqrt(cappedMagSquared);
        state.cappedLogReadMagnitude = std::sqrt(cappedLogMagSquared);
        logging::info("Precomputed read magnitude: {:.6f}, log-scaled magnitude: {:.6f}, total read seed frequency: {}, unique read seeds: {}", 
                    state.readMagnitude, state.logReadMagnitude, state.totalReadSeedFrequency, state.readUniqueSeedCount);
        logging::info("Capped (max=100) magnitude: {:.6f}, capped log magnitude: {:.6f}",
                    state.cappedReadMagnitude, state.cappedLogReadMagnitude);
        
        // ============================================================================
        // Compute sigmoid-scaled read counts with median seed frequency as midpoint
        // Sigmoid adaptively compresses high-frequency seeds (amplicon bias) while
        // preserving sensitivity around the median coverage level
        // Formula: scaledCount = 2*median / (1 + exp(-k * (count/median - 1)))
        // ============================================================================
        {
            // First, compute median seed frequency in reads
            std::vector<int64_t> freqs;
            freqs.reserve(state.seedFreqInReads.size());
            for (const auto& [hash, count] : state.seedFreqInReads) {
                freqs.push_back(count);
            }
            std::sort(freqs.begin(), freqs.end());
            
            if (!freqs.empty()) {
                size_t mid = freqs.size() / 2;
                if (freqs.size() % 2 == 0) {
                    state.medianReadSeedFrequency = (freqs[mid - 1] + freqs[mid]) / 2.0;
                } else {
                    state.medianReadSeedFrequency = static_cast<double>(freqs[mid]);
                }
            } else {
                state.medianReadSeedFrequency = 1.0;  // Fallback
            }
            
            const double median = state.medianReadSeedFrequency;
            const double k = 2.0;  // Steepness parameter (controls transition sharpness)
            
            // Pre-size sigmoid map
            state.sigmoidReadCounts.reserve(state.seedFreqInReads.size());
            double sigMagSquared = 0.0;
            
            for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
                // Sigmoid scaling centered on median:
                // scaledCount = 2*median / (1 + exp(-k * (count/median - 1)))
                // When count = median: scaledCount = median (midpoint of sigmoid)
                // When count << median: scaledCount → 0 (but approaches linearly near 0)
                // When count >> median: scaledCount → 2*median (asymptotic cap)
                double ratio = static_cast<double>(readCount) / median;
                double scaledCount = 2.0 * median / (1.0 + std::exp(-k * (ratio - 1.0)));
                
                state.sigmoidReadCounts[seedHash] = scaledCount;
                sigMagSquared += scaledCount * scaledCount;
            }
            
            state.sigmoidReadMagnitude = std::sqrt(sigMagSquared);
            logging::info("Sigmoid scaling: median read seed freq={:.1f}, steepness k={:.1f}, sigmoid magnitude={:.6f}",
                         median, k, state.sigmoidReadMagnitude);
        }
    }
    auto time_magnitude_end = std::chrono::high_resolution_clock::now();
    auto duration_magnitude = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_magnitude_end - time_magnitude_start);
    logging::info("[PLACEMENT TIMING]   Read magnitude precomputation: {}ms", duration_magnitude.count());
    
    // Compute reference genome statistics (median, max seed counts for N-adjustment)
    {
        std::vector<uint32_t> seedCounts;
        seedCounts.reserve(genomeUniqueSeedReader.size());
        for (size_t i = 0; i < genomeUniqueSeedReader.size(); ++i) {
            uint32_t count = genomeUniqueSeedReader[i];
            if (count > 0) {  // Skip root node if it has 0 seeds
                seedCounts.push_back(count);
            }
        }
        
        if (!seedCounts.empty()) {
            // Compute median
            std::sort(seedCounts.begin(), seedCounts.end());
            size_t mid = seedCounts.size() / 2;
            if (seedCounts.size() % 2 == 0) {
                state.medianGenomeUniqueSeedCount = (seedCounts[mid - 1] + seedCounts[mid]) / 2.0;
            } else {
                state.medianGenomeUniqueSeedCount = seedCounts[mid];
            }
            state.maxGenomeUniqueSeedCount = seedCounts.back();
            
            // Compute mean and stddev
            double sum = 0.0;
            for (uint32_t c : seedCounts) sum += c;
            state.meanGenomeUniqueSeedCount = sum / seedCounts.size();
            
            double sumSquares = 0.0;
            for (uint32_t c : seedCounts) {
                double diff = c - state.meanGenomeUniqueSeedCount;
                sumSquares += diff * diff;
            }
            state.stdGenomeUniqueSeedCount = std::sqrt(sumSquares / seedCounts.size());
            
            // Compute incompleteness threshold: mean - 2*std
            // Only genomes below this threshold are penalized (likely have N's or missing data)
            state.incompletenessThreshold = state.meanGenomeUniqueSeedCount - 2.0 * state.stdGenomeUniqueSeedCount;
            if (state.incompletenessThreshold < 0) state.incompletenessThreshold = 0;
            
            logging::info("Genome seed count stats: median={:.0f}, max={:.0f}, mean={:.1f}, std={:.1f}, threshold={:.1f}",
                        state.medianGenomeUniqueSeedCount, state.maxGenomeUniqueSeedCount,
                        state.meanGenomeUniqueSeedCount, state.stdGenomeUniqueSeedCount,
                        state.incompletenessThreshold);
        }
    }
    
    // ============================================================================
    // Compute coverage-weighted seed counts
    // Weight function penalizes singletons (likely errors) and over-represented seeds (amplicon bias)
    // Only applies when expected coverage > 1 (low coverage: all seeds get weight 1.0)
    // ============================================================================
    {
        // Expected coverage = total seed frequency / expected unique seeds per genome
        // Using median genome seed count as the expectation for unique seeds in the sample
        if (state.medianGenomeUniqueSeedCount > 0) {
            state.expectedSeedCoverage = static_cast<double>(state.totalReadSeedFrequency) / 
                                          state.medianGenomeUniqueSeedCount;
        } else {
            state.expectedSeedCoverage = 1.0;  // Fallback: no penalty
        }
        
        const double E = state.expectedSeedCoverage;
        logging::info("Expected seed coverage: {:.2f}x (total seeds: {}, expected unique: {:.0f})",
                     E, state.totalReadSeedFrequency, state.medianGenomeUniqueSeedCount);
        
        // Pre-size coverageWeightedCounts map
        state.coverageWeightedCounts.reserve(state.seedFreqInReads.size());
        
        double covMagSquared = 0.0;
        for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
            double weight = 1.0;  // Default: full weight
            
            if (E > 1.0) {
                const double n = static_cast<double>(readCount);
                
                // Singleton penalty: down-weight seeds appearing only once
                // Weight = 1/log2(E) for singletons, which decreases as coverage increases
                // At E=10: singletons get ~0.30; at E=100: ~0.15
                double singletonPenalty = 1.0;
                if (readCount == 1) {
                    singletonPenalty = std::min(1.0, 1.0 / std::log2(E));
                }
                
                // Over-representation penalty: down-weight seeds appearing >> expected
                // Uses hyperbolic decay: weight = E/n for n > E, with floor at 0.25
                double overrepPenalty = 1.0;
                if (n > E) {
                    overrepPenalty = std::max(0.25, E / n);
                }
                
                weight = singletonPenalty * overrepPenalty;
            }
            // else: E <= 1.0 means low coverage, keep weight = 1.0 (can't distinguish signal)
            
            // Store weighted count (weight * readCount for cosine numerator contribution)
            double weightedCount = weight * static_cast<double>(readCount);
            state.coverageWeightedCounts[seedHash] = weightedCount;
            covMagSquared += weightedCount * weightedCount;
        }
        
        state.covWeightedMagnitude = std::sqrt(covMagSquared);
        
        // Log weight distribution summary
        if (E > 1.0) {
            size_t singletons = 0, overrep = 0;
            double minWeight = 1.0, maxWeight = 0.0;
            for (const auto& [hash, wc] : state.coverageWeightedCounts) {
                auto it = state.seedFreqInReads.find(hash);
                if (it != state.seedFreqInReads.end()) {
                    double w = wc / static_cast<double>(it->second);
                    if (it->second == 1) singletons++;
                    if (static_cast<double>(it->second) > E) overrep++;
                    minWeight = std::min(minWeight, w);
                    maxWeight = std::max(maxWeight, w);
                }
            }
            logging::info("Coverage weights: {} singletons, {} over-represented (>E), weight range [{:.3f}, {:.3f}]",
                         singletons, overrep, minWeight, maxWeight);
            logging::info("Coverage-weighted magnitude: {:.6f} (vs raw: {:.6f})",
                         state.covWeightedMagnitude, state.readMagnitude);
        } else {
            logging::info("Low coverage (E={:.2f}): all seeds get weight 1.0", E);
        }
    }
    
    state.root = liteTree->root;

    placement::NodeMetrics rootMetrics;
    // CRITICAL: Root starts with empty genome (all genome metrics = 0) but with
    // read-based metrics pre-initialized. This represents the state BEFORE the root node
    // processes its mutations (i.e., a hypothetical empty ancestor genome).
    // The root's seed changes will then apply deltas to compute the actual root state.
    //
    // Pre-initialize read-based components:
    // - weightedJaccardDenominator starts with sum of read frequencies (max(read, 0) = read for all seeds)
    // All genome-based components start at 0 and will be updated via computeChildMetrics.
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        // Note: weightedJaccardDenominator is now computed from pre-indexed genomeTotalSeedFrequency
    }

    uint32_t currentDfsIndex = 0;

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
        
    } else {
        logging::err("LiteTree or root is null, cannot perform placement");
        return;
    }
    auto time_traversal_end = std::chrono::high_resolution_clock::now();
    auto duration_traversal = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_traversal_end - time_traversal_start);
    logging::info("[PLACEMENT TIMING]   Tree traversal: {}ms", duration_traversal.count());
    
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
        result.readMagnitude = state.readMagnitude;
        
        logging::info("Stored alignment parameters: k={}, s={}, t={}, open={}", 
                     result.k, result.s, result.t, result.open);
    }
    
    auto time_write_start = std::chrono::high_resolution_clock::now();
    std::string placementsFilePath = outputPath;
    std::ofstream placementsFile(placementsFilePath);
    if (placementsFile.is_open()) {
        placementsFile << "metric\tscore\thits\tnodes\n";
        
        // Raw score - join all tied nodes with comma
        placementsFile << "raw\t" << result.bestRawSeedMatchScore 
                      << "\t" << result.maxHitsInAnyGenome << "\t";
        if (!result.tiedRawSeedMatchNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedRawSeedMatchNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedRawSeedMatchNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestRawSeedMatchNodeId;
        }
        placementsFile << "\n";
        
        // Jaccard score - join all tied nodes with comma
        placementsFile << "jaccard\t" << std::fixed << std::setprecision(6) << result.bestJaccardScore << "\t\t";
        if (!result.tiedJaccardNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedJaccardNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedJaccardNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestJaccardNodeId;
        }
        placementsFile << "\n";
        
        // Cosine score - join all tied nodes with comma
        placementsFile << "cosine\t" << std::fixed << std::setprecision(6) << result.bestCosineScore << "\t\t";
        if (!result.tiedCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestCosineNodeId;
        }
        placementsFile << "\n";
        
        // Weighted Jaccard score - join all tied nodes with comma
        placementsFile << "weighted_jaccard\t" << std::fixed << std::setprecision(6) << result.bestWeightedJaccardScore << "\t\t";
        if (!result.tiedWeightedJaccardNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedWeightedJaccardNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedWeightedJaccardNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestWeightedJaccardNodeId;
        }
        placementsFile << "\n";
        
        // Log Raw score
        placementsFile << "log_raw\t" << std::fixed << std::setprecision(6) << result.bestLogRawScore << "\t\t";
        if (!result.tiedLogRawNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogRawNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogRawNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogRawNodeId;
        }
        placementsFile << "\n";
        
        // Log Cosine score
        placementsFile << "log_cosine\t" << std::fixed << std::setprecision(6) << result.bestLogCosineScore << "\t\t";
        if (!result.tiedLogCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedLogCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedLogCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestLogCosineNodeId;
        }
        placementsFile << "\n";
        
        // IDF Cosine score
        placementsFile << "idf_cosine\t" << std::fixed << std::setprecision(6) << result.bestIdfCosineScore << "\t\t";
        if (!result.tiedIdfCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedIdfCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedIdfCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestIdfCosineNodeId;
        }
        placementsFile << "\n";
        
        // Coverage Cosine score
        placementsFile << "cov_cosine\t" << std::fixed << std::setprecision(6) << result.bestCovCosineScore << "\t\t";
        if (!result.tiedCovCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedCovCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedCovCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestCovCosineNodeId;
        }
        placementsFile << "\n";
        
        // Cap Cosine score
        placementsFile << "cap_cosine\t" << std::fixed << std::setprecision(6) << result.bestCapCosineScore << "\t\t";
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
        placementsFile << "cap_log_cosine\t" << std::fixed << std::setprecision(6) << result.bestCapLogCosineScore << "\t\t";
        if (!result.tiedCapLogCosineNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedCapLogCosineNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedCapLogCosineNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestCapLogCosineNodeId;
        }
        placementsFile << "\n";
        
        // Containment score
        placementsFile << "containment\t" << std::fixed << std::setprecision(6) << result.bestContainmentScore << "\t\t";
        if (!result.tiedContainmentNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedContainmentNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedContainmentNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestContainmentNodeId;
        }
        placementsFile << "\n";
        
        // Weighted Containment score
        placementsFile << "weighted_containment\t" << std::fixed << std::setprecision(6) << result.bestWeightedContainmentScore << "\t\t";
        if (!result.tiedWeightedContainmentNodeIndices.empty()) {
            for (size_t i = 0; i < result.tiedWeightedContainmentNodeIndices.size(); ++i) {
                if (i > 0) placementsFile << ",";
                placementsFile << liteTree->resolveNodeId(result.tiedWeightedContainmentNodeIndices[i]);
            }
        } else {
            placementsFile << result.bestWeightedContainmentNodeId;
        }
        placementsFile << "\n";
        
        placementsFile.close();
        logging::info("Wrote placement results to {}", placementsFilePath);
    } else {
        logging::err("Failed to open placements file: {}", placementsFilePath);
    }
    auto time_write_end = std::chrono::high_resolution_clock::now();
    auto duration_write = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_write_end - time_write_start);
    logging::info("[PLACEMENT TIMING]   Result file writing: {}ms", duration_write.count());
    
    auto placement_total_end = std::chrono::high_resolution_clock::now();
    auto duration_placement_total = std::chrono::duration_cast<std::chrono::milliseconds>(
        placement_total_end - placement_total_start);
    
    logging::info("Placement completed. Best Jaccard score: {} (node: {})", 
                 result.bestJaccardScore, result.bestJaccardNodeId);
    if (result.tiedJaccardNodeIndices.size() > 1) {
        double jaccard_tol = std::max(result.bestJaccardScore * 0.0001, 1e-9);
        logging::info("  {} nodes tied for best Jaccard score (within relative tolerance {:.2e})", 
                     result.tiedJaccardNodeIndices.size(), jaccard_tol);
    }
    
    logging::info("\n=== PLACEMENT INTERNAL TIMING BREAKDOWN ===");
    logging::info("HashDeltaLoading: {}ms", duration_hash_delta.count());
    logging::info("ReadProcessing: {}ms", duration_read_processing.count());
    logging::info("ReadMagnitude: {}ms", duration_magnitude.count());
    logging::info("TreeTraversal: {}ms", duration_traversal.count());
    logging::info("FileWriting: {}ms", duration_write.count());
    logging::info("PlacementTotal: {}ms", duration_placement_total.count());
    logging::info("==========================================\n");
    
    logging::info("Best Cosine score: {} (node: {})", 
                 result.bestCosineScore, result.bestCosineNodeId);
    if (result.tiedCosineNodeIndices.size() > 1) {
        double cosine_tol = std::max(result.bestCosineScore * 0.0001, 1e-9);
        logging::info("  {} nodes tied for best Cosine score (within relative tolerance {:.2e})", 
                     result.tiedCosineNodeIndices.size(), cosine_tol);
    }
    
    logging::info("Best Weighted Jaccard: {} (node: {})", 
                 result.bestWeightedJaccardScore, result.bestWeightedJaccardNodeId);
    if (result.tiedWeightedJaccardNodeIndices.size() > 1) {
        double weighted_tol = std::max(result.bestWeightedJaccardScore * 0.0001, 1e-9);
        logging::info("  {} nodes tied for best Weighted Jaccard (within relative tolerance {:.2e})", 
                     result.tiedWeightedJaccardNodeIndices.size(), weighted_tol);
    }
    
    logging::info("Best Raw matches: {} (node: {})", 
                 result.bestRawSeedMatchScore, result.bestRawSeedMatchNodeId);
    if (result.tiedRawSeedMatchNodeIndices.size() > 1) {
        logging::info("  {} nodes tied for best Raw score", result.tiedRawSeedMatchNodeIndices.size());
    }
    
    logging::info("Best LogRAW score: {:.6f} (node: {})", 
                 result.bestLogRawScore, result.bestLogRawNodeId);
    if (result.tiedLogRawNodeIndices.size() > 1) {
        logging::info("  {} nodes tied for best LogRAW score", result.tiedLogRawNodeIndices.size());
    }
}

} // namespace placement


#include "panman.hpp"
#include "placement.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"
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
constexpr size_t computeHomopolymerHash(char base, int k) {
    size_t baseVal = seeding::chash(base);
    size_t compVal;
    switch (base) {
        case 'A':
        case 'a': compVal = seeding::chash('T'); break;
        case 'C':
        case 'c': compVal = seeding::chash('G'); break;
        case 'G':
        case 'g': compVal = seeding::chash('C'); break;
        case 'T':
        case 't': compVal = seeding::chash('A'); break;
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
constexpr std::array<size_t, 4> getHomopolymerHashes(int k) {
    return {computeHomopolymerHash('A', k),
            computeHomopolymerHash('C', k),
            computeHomopolymerHash('G', k),
            computeHomopolymerHash('T', k)};
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

void extractReadSequences(const std::string& readPath1,
                          const std::string& readPath2,
                          std::vector<std::string>& readSequences) {
    mgsr::FastqFile fq1(readPath1);
    kseq_t* seq = kseq_init(fileno(fq1.fp));
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
    }
    kseq_destroy(seq);

    if (readPath2.size() > 0) {
        mgsr::FastqFile fq2(readPath2);
        seq = kseq_init(fileno(fq2.fp));

        line = 0;
        int forwardReads = readSequences.size();
        while ((line = kseq_read(seq)) >= 0) {
            readSequences.push_back(seq->seq.s);
        }
        kseq_destroy(seq);

        if (readSequences.size() != static_cast<size_t>(forwardReads) * 2) {
            output::error("File {} does not contain the same number of reads as {}", readPath2, readPath1);
            exit(1);
        }

        // Shuffle reads together, so that pairs are next to each other
        seeding::perfect_shuffle(readSequences);
    }
}

// Extract full FASTQ data including sequences, qualities, and names
void extractFullFastqData(const std::string& readPath1,
                          const std::string& readPath2,
                          std::vector<std::string>& readSequences,
                          std::vector<std::string>& readQuals,
                          std::vector<std::string>& readNames) {
    mgsr::FastqFile fq1(readPath1);
    kseq_t* seq = kseq_init(fileno(fq1.fp));
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
        readNames.push_back(seq->name.s);
        // Quality may be empty for FASTA files
        readQuals.push_back(seq->qual.l > 0 ? seq->qual.s : std::string(seq->seq.l, 'I'));
    }
    kseq_destroy(seq);

    if (readPath2.size() > 0) {
        mgsr::FastqFile fq2(readPath2);
        seq = kseq_init(fileno(fq2.fp));

        line = 0;
        int forwardReads = readSequences.size();
        while ((line = kseq_read(seq)) >= 0) {
            readSequences.push_back(seq->seq.s);
            readNames.push_back(seq->name.s);
            readQuals.push_back(seq->qual.l > 0 ? seq->qual.s : std::string(seq->seq.l, 'I'));
        }
        kseq_destroy(seq);

        if (readSequences.size() != static_cast<size_t>(forwardReads) * 2) {
            output::error("File {} does not contain the same number of reads as {}", readPath2, readPath1);
            exit(1);
        }

        // Shuffle reads together, so that pairs are next to each other
        seeding::perfect_shuffle(readSequences);
        seeding::perfect_shuffle(readQuals);
        seeding::perfect_shuffle(readNames);
    }
}

}  // anonymous namespace

void placement::NodeMetrics::computeChildMetrics(placement::NodeMetrics& childMetrics,
                                                 std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges,
                                                 const placement::PlacementGlobalState& state) {
    const size_t numChanges = seedChanges.size();

    // Early exit for empty changes
    if (numChanges == 0) [[unlikely]]
        return;
    constexpr size_t PREFETCH_DISTANCE = 24;  // Increased to 24 for deeper pipeline (L3 latency ~40 cycles)
    constexpr size_t BATCH_SIZE = 64;         // Process in cache-line aligned batches

    // OPTIMIZATION: Pre-fetch both seedChanges AND hash table buckets aggressively
    const size_t initialPrefetchCount = std::min(numChanges, PREFETCH_DISTANCE);
    for (size_t i = 0; i < initialPrefetchCount; ++i) {
        const auto& [seedHash, _, __] = seedChanges[i];
        __builtin_prefetch(&seedChanges[i], 0, 0);
        state.seedFreqInReads.prefetch(seedHash);
        state.seedInverseGenomeCounts.prefetch(seedHash);
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
                state.seedInverseGenomeCounts.prefetch(nextHash);
            }

            const auto& [seedHash, parentCount, childCount] = seedChanges[i];

            const int64_t freqDelta = childCount - parentCount;

            // Pre-compute log-scaled genome counts for cosine metric
            const double logChild = (childCount > 0) ? std::log1p(static_cast<double>(childCount)) : 0.0;
            const double logParent = (parentCount > 0) ? std::log1p(static_cast<double>(parentCount)) : 0.0;

            // ========================================
            // GENOME-ONLY METRICS (computed for ALL seed changes, even if not in reads)
            // ========================================

            // Genome magnitude squared delta: log(1+child)² - log(1+parent)²
            // Use log-scale to match the log-scaled read vector for proper cosine similarity
            const double magDelta = logChild * logChild - logParent * logParent;
            childMetrics.genomeMagnitudeSquared += magDelta;

            // Genome unique seed count delta (branchless)
            const int64_t wasPresent = (parentCount > 0) ? 1 : 0;
            const int64_t isPresent = (childCount > 0) ? 1 : 0;
            childMetrics.genomeUniqueSeedCount += (isPresent - wasPresent);

            // Fast path: skip unchanged seeds
            if (freqDelta == 0) [[likely]]
                continue;

            auto logIt = state.logReadCounts.find(seedHash);
            if (logIt == state.logReadCounts.end()) [[likely]]
                continue;

            const double logReadCount = logIt->second;

            // Presence intersection count delta
            const int64_t becamePresent = (parentCount == 0) & (childCount != 0);
            const int64_t becameAbsent = (childCount == 0) & (parentCount != 0);
            const int64_t presenceDelta = becamePresent - becameAbsent;
            childMetrics.presenceIntersectionCount += presenceDelta;

            // LogRaw numerator: Σ(log(1+readCount) / genomeCount)
            const double oldLogContrib = (parentCount > 0) ? (logReadCount / parentCount) : 0.0;
            const double newLogContrib = (childCount > 0) ? (logReadCount / childCount) : 0.0;
            childMetrics.logRawNumerator += (newLogContrib - oldLogContrib);

            // LogCosine numerator: Σ(log(1+readCount) × log(1+genomeCount))
            // Both vectors in log-space for proper cosine similarity
            childMetrics.logCosineNumerator += logReadCount * (logChild - logParent);

            // Weighted containment numerator: Σ(1/nodeGenomeCount_i) for seeds in reads ∩ genome
            // Unlike regular containment which uses binary presence, this weights
            // each seed by 1/(node's genome count), rewarding genomes where matching
            // seeds are rare/unique. Delta: 1/childCount - 1/parentCount for seeds in reads.
            {
                const double oldWcContrib = (parentCount > 0) ? (1.0 / parentCount) : 0.0;
                const double newWcContrib = (childCount > 0) ? (1.0 / childCount) : 0.0;
                childMetrics.weightedContainmentNumerator += (newWcContrib - oldWcContrib);
            }

            // Log containment numerator: Σ(log(1+r_i)) for seeds in reads ∩ genome
            // Like containment but weights each seed by its log read abundance
            childMetrics.logContainmentNumerator += presenceDelta * logReadCount;
        }
    }
}

namespace placement {

using ::panmanUtils::Node;
using ::panmanUtils::Tree;

class PlacementResult;

// Macro to generate score update functions - reduces repetitive code
// Each function: stores score on node, updates best score/index, tracks ties
#define DEFINE_UPDATE_SCORE_FUNC(FuncName, nodeScoreField, bestScore, bestNodeIndex, tiedIndices)   \
    void PlacementResult::FuncName(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node) { \
        if (node) {                                                                                 \
            node->nodeScoreField = static_cast<float>(score);                                       \
        }                                                                                           \
        double tolerance = std::max(bestScore * 0.0001, 1e-9);                                      \
        if (score > bestScore + tolerance) {                                                        \
            bestScore = score;                                                                      \
            bestNodeIndex = nodeIndex;                                                              \
            tiedIndices.clear();                                                                    \
            tiedIndices.push_back(nodeIndex);                                                       \
        } else if (score >= bestScore - tolerance && score > 0) {                                   \
            if (tiedIndices.empty() || tiedIndices.back() != bestNodeIndex) {                       \
                tiedIndices.push_back(bestNodeIndex);                                               \
            }                                                                                       \
            if (nodeIndex != bestNodeIndex) {                                                       \
                tiedIndices.push_back(nodeIndex);                                                   \
            }                                                                                       \
        }                                                                                           \
    }

// Generate score update functions
DEFINE_UPDATE_SCORE_FUNC(updateLogRawScore, logRawScore, bestLogRawScore, bestLogRawNodeIndex, tiedLogRawNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(
    updateLogCosineScore, logCosineScore, bestLogCosineScore, bestLogCosineNodeIndex, tiedLogCosineNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateContainmentScore,
                         containmentScore,
                         bestContainmentScore,
                         bestContainmentNodeIndex,
                         tiedContainmentNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateWeightedContainmentScore,
                         weightedContainmentScore,
                         bestWeightedContainmentScore,
                         bestWeightedContainmentNodeIndex,
                         tiedWeightedContainmentNodeIndices)
DEFINE_UPDATE_SCORE_FUNC(updateLogContainmentScore,
                         logContainmentScore,
                         bestLogContainmentScore,
                         bestLogContainmentNodeIndex,
                         tiedLogContainmentNodeIndices)

#undef DEFINE_UPDATE_SCORE_FUNC

// Helper to deduplicate a tied-indices vector and pick the lowest index as best
static void finalizeTiedIndices(std::vector<uint32_t>& tied, uint32_t& bestIndex) {
    if (tied.empty()) return;
    // Sort and deduplicate
    std::sort(tied.begin(), tied.end());
    tied.erase(std::unique(tied.begin(), tied.end()), tied.end());
    // Deterministic tie-breaking: lowest node index wins (DFS order is stable)
    bestIndex = tied.front();
}

void PlacementResult::resolveNodeIds(panmapUtils::LiteTree* liteTree) {
    if (!liteTree) return;

    // Deterministic tie-breaking before resolution
    finalizeTiedIndices(tiedLogRawNodeIndices, bestLogRawNodeIndex);
    finalizeTiedIndices(tiedLogCosineNodeIndices, bestLogCosineNodeIndex);
    finalizeTiedIndices(tiedContainmentNodeIndices, bestContainmentNodeIndex);
    finalizeTiedIndices(tiedWeightedContainmentNodeIndices, bestWeightedContainmentNodeIndex);
    finalizeTiedIndices(tiedLogContainmentNodeIndices, bestLogContainmentNodeIndex);

    if (bestLogRawNodeIndex != UINT32_MAX) {
        bestLogRawNodeId = liteTree->resolveNodeId(bestLogRawNodeIndex);
    }
    if (bestLogCosineNodeIndex != UINT32_MAX) {
        bestLogCosineNodeId = liteTree->resolveNodeId(bestLogCosineNodeIndex);
    }
    if (bestContainmentNodeIndex != UINT32_MAX) {
        bestContainmentNodeId = liteTree->resolveNodeId(bestContainmentNodeIndex);
    }
    if (bestWeightedContainmentNodeIndex != UINT32_MAX) {
        bestWeightedContainmentNodeId = liteTree->resolveNodeId(bestWeightedContainmentNodeIndex);
    }
    if (bestLogContainmentNodeIndex != UINT32_MAX) {
        bestLogContainmentNodeId = liteTree->resolveNodeId(bestLogContainmentNodeIndex);
    }
    // Resolve per-metric refined node IDs
    auto resolveRefined = [&](RefinedResult& r) {
        if (r.nodeIndex != UINT32_MAX) {
            r.nodeId = liteTree->resolveNodeId(r.nodeIndex);
        }
    };
    resolveRefined(refinedLogRaw);
    resolveRefined(refinedLogCosine);
    resolveRefined(refinedContainment);
    resolveRefined(refinedWeightedContainment);
    resolveRefined(refinedLogContainment);
}

// After k-mer scoring, refine top candidates by full minimap2 alignment

// Get all nodes within phylogenetic distance `radius` of a given node
// Uses BFS traversal following parent/child edges
std::vector<panmapUtils::LiteNode*> getNodesWithinRadius(panmapUtils::LiteNode* startNode, int radius, int maxNodes) {
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
// Returns negative total edit distance (higher = fewer errors = better)
// When pairedEnd=true, reads are interleaved R1_0,R2_0,R1_1,R2_1,...
int64_t scoreNodeByAlignment(panmanUtils::Tree* fullTree,
                             const std::string& nodeId,
                             const std::vector<std::string>& readSequences,
                             int kmerSize,
                             bool pairedEnd) {
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

    // Call minimap2 scoring function with paired-end awareness
    int64_t score =
        score_reads_vs_reference(genomeSeq.c_str(), n_reads, readPtrs.data(), readLens.data(), kmerSize, pairedEnd);

    return score;
}

// Main refinement function: per-metric candidate selection, shared alignment scoring
// Each seed metric independently nominates its top candidates + neighbors.
// All unique candidates are aligned once, then each metric picks its best.
void refineTopCandidates(panmapUtils::LiteTree* liteTree,
                         panmanUtils::Tree* fullTree,
                         const std::vector<std::string>& readSequences,
                         PlacementResult& result,
                         const TraversalParams& params,
                         bool pairedEnd) {
    if (!fullTree || !liteTree || readSequences.empty()) {
        logging::warn("Refinement skipped: missing tree data or no reads");
        return;
    }

    auto time_refine_start = std::chrono::high_resolution_clock::now();

    // Step 1: For each metric, collect its own top candidate set
    // Returns the set of candidate node indices for one metric
    auto getTopForMetric = [&](auto scoreGetter, const char* metricName) -> absl::flat_hash_set<uint32_t> {
        absl::flat_hash_set<uint32_t> metricCandidates;
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
        if (scoredNodes.empty()) return metricCandidates;

        std::sort(scoredNodes.begin(), scoredNodes.end(), std::greater<>());

        size_t numTop = std::min(static_cast<size_t>(scoredNodes.size() * params.refineTopPct),
                                 static_cast<size_t>(params.refineMaxTopN));
        numTop = std::max(numTop, size_t(1));

        for (size_t i = 0; i < numTop && i < scoredNodes.size(); i++) {
            metricCandidates.insert(scoredNodes[i].second);
        }
        logging::debug("Refinement: {} selected {} candidates", metricName, metricCandidates.size());
        return metricCandidates;
    };

    // Collect per-metric candidate sets
    auto logRawCands = getTopForMetric([](auto* n) { return n->logRawScore; }, "LogRaw");
    auto logCosineCands = getTopForMetric([](auto* n) { return n->logCosineScore; }, "LogCosine");
    auto containCands = getTopForMetric([](auto* n) { return n->containmentScore; }, "Containment");
    auto wContainCands = getTopForMetric([](auto* n) { return n->weightedContainmentScore; }, "WeightedContainment");
    auto logContainCands = getTopForMetric([](auto* n) { return n->logContainmentScore; }, "LogContainment");

    // Always include the unrefined best node for each metric in its candidate set.
    // At high coverage, seed scores saturate and the correct node may fall outside
    // the top-N% cutoff. Including the seed-phase winner ensures refinement always
    // considers the best unrefined placement.
    if (result.bestLogRawNodeIndex != UINT32_MAX) logRawCands.insert(result.bestLogRawNodeIndex);
    if (result.bestLogCosineNodeIndex != UINT32_MAX) logCosineCands.insert(result.bestLogCosineNodeIndex);
    if (result.bestContainmentNodeIndex != UINT32_MAX) containCands.insert(result.bestContainmentNodeIndex);
    if (result.bestWeightedContainmentNodeIndex != UINT32_MAX)
        wContainCands.insert(result.bestWeightedContainmentNodeIndex);
    if (result.bestLogContainmentNodeIndex != UINT32_MAX) logContainCands.insert(result.bestLogContainmentNodeIndex);

    // Step 2: Expand each metric's candidates with neighbors, build union for alignment
    absl::flat_hash_set<uint32_t> allCandidates;  // Union for alignment scoring

    // For each metric, expand candidates with neighbors and track membership
    struct MetricCandidateSet {
        absl::flat_hash_set<uint32_t> expanded;  // candidates + their neighbors
    };

    auto expandWithNeighbors = [&](const absl::flat_hash_set<uint32_t>& baseCands) -> MetricCandidateSet {
        MetricCandidateSet mcs;
        for (uint32_t nodeIdx : baseCands) {
            mcs.expanded.insert(nodeIdx);
            allCandidates.insert(nodeIdx);

            panmapUtils::LiteNode* node = liteTree->dfsIndexToNode[nodeIdx];
            if (!node) continue;
            auto neighbors = getNodesWithinRadius(node, params.refineNeighborRadius, params.refineMaxNeighborN);
            for (auto* neighbor : neighbors) {
                if (neighbor) {
                    mcs.expanded.insert(neighbor->nodeIndex);
                    allCandidates.insert(neighbor->nodeIndex);
                }
            }
        }
        return mcs;
    };

    auto logRawExpanded = expandWithNeighbors(logRawCands);
    auto logCosineExpanded = expandWithNeighbors(logCosineCands);
    auto containExpanded = expandWithNeighbors(containCands);
    auto wContainExpanded = expandWithNeighbors(wContainCands);
    auto logContainExpanded = expandWithNeighbors(logContainCands);

    if (allCandidates.empty()) {
        logging::warn("Refinement skipped: no nodes with positive scores");
        return;
    }

    logging::info("Refinement: {} unique candidates total from all metrics ({}%)",
                  allCandidates.size(),
                  params.refineTopPct * 100);

    // Step 3: Align reads to ALL unique candidates once (shared work)
    absl::flat_hash_map<uint32_t, int64_t> alignmentScoreMap;  // nodeIndex -> score

    // Log seed winners
    logging::info("Refinement seed winners: LogRaw={} LogCosine={} Contain={} WContain={} LogContain={}",
                  result.bestLogRawNodeIndex,
                  result.bestLogCosineNodeIndex,
                  result.bestContainmentNodeIndex,
                  result.bestWeightedContainmentNodeIndex,
                  result.bestLogContainmentNodeIndex);

    for (uint32_t nodeIdx : allCandidates) {
        panmapUtils::LiteNode* node = liteTree->dfsIndexToNode[nodeIdx];
        if (!node) continue;

        int64_t score = scoreNodeByAlignment(fullTree, node->identifier, readSequences, params.k, pairedEnd);
        alignmentScoreMap[nodeIdx] = score;
    }

    if (alignmentScoreMap.empty()) {
        logging::warn("Refinement produced no alignment scores");
        return;
    }

    // Step 4: For each metric, find the best candidate from its own expanded set
    auto findBestForMetric = [&](const MetricCandidateSet& mcs, auto seedScoreGetter) -> std::pair<int64_t, uint32_t> {
        int64_t bestScore = INT64_MIN;
        uint32_t bestIdx = UINT32_MAX;
        for (uint32_t nodeIdx : mcs.expanded) {
            auto it = alignmentScoreMap.find(nodeIdx);
            if (it == alignmentScoreMap.end()) continue;
            int64_t score = it->second;
            if (score > bestScore || (score == bestScore && bestIdx != UINT32_MAX)) {
                // Tie-break by seed score
                if (score == bestScore) {
                    auto* nodeA = liteTree->dfsIndexToNode[nodeIdx];
                    auto* nodeB = liteTree->dfsIndexToNode[bestIdx];
                    float seedA = nodeA ? seedScoreGetter(nodeA) : 0.0f;
                    float seedB = nodeB ? seedScoreGetter(nodeB) : 0.0f;
                    if (seedA <= seedB) continue;  // Current best still wins
                }
                bestScore = score;
                bestIdx = nodeIdx;
            }
        }
        return {bestScore, bestIdx};
    };

    auto [lrScore, lrIdx] = findBestForMetric(logRawExpanded, [](auto* n) { return n->logRawScore; });
    auto [lcScore, lcIdx] = findBestForMetric(logCosineExpanded, [](auto* n) { return n->logCosineScore; });
    auto [cScore, cIdx] = findBestForMetric(containExpanded, [](auto* n) { return n->containmentScore; });
    auto [wcScore, wcIdx] = findBestForMetric(wContainExpanded, [](auto* n) { return n->weightedContainmentScore; });
    auto [lgcScore, lgcIdx] = findBestForMetric(logContainExpanded, [](auto* n) { return n->logContainmentScore; });

    // Store results
    if (lrIdx != UINT32_MAX) {
        result.refinedLogRaw = {static_cast<double>(lrScore), lrIdx, ""};
    }
    if (lcIdx != UINT32_MAX) {
        result.refinedLogCosine = {static_cast<double>(lcScore), lcIdx, ""};
    }
    if (cIdx != UINT32_MAX) {
        result.refinedContainment = {static_cast<double>(cScore), cIdx, ""};
    }
    if (wcIdx != UINT32_MAX) {
        result.refinedWeightedContainment = {static_cast<double>(wcScore), wcIdx, ""};
    }
    if (lgcIdx != UINT32_MAX) {
        result.refinedLogContainment = {static_cast<double>(lgcScore), lgcIdx, ""};
    }
    result.refinementWasRun = true;

    auto time_refine_end = std::chrono::high_resolution_clock::now();
    auto duration_refine = std::chrono::duration_cast<std::chrono::milliseconds>(time_refine_end - time_refine_start);

    logging::info(
        "Refinement complete: {} candidates scored in {}ms", alignmentScoreMap.size(), duration_refine.count());
    logging::info("  refined_log_raw: {} | refined_containment: {} | refined_log_cosine: {}", lrScore, cScore, lcScore);
}

// Placement helper using BFS traversal with delta-based metrics
void placeLiteHelperBFS(std::vector<panmapUtils::LiteNode*>& nodes,
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

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, num_nodes, grain_size), [&](const tbb::blocked_range<size_t>& r) {
                auto& tls = thread_data.local();
                size_t batch_counter = 0;

                for (size_t i = r.begin(); i != r.end(); ++i) {
                    panmapUtils::LiteNode* node = current_nodes[i];
                    const placement::NodeMetrics& p_metrics = current_metrics[i];
                    const absl::flat_hash_map<uint64_t, int64_t>& p_seed_counts =
                        params.verify_scores ? current_seeds[i] : empty_seed_counts;

                    if (!node) continue;

                    // Aggressive software prefetch for next 4 nodes
                    constexpr size_t PREFETCH_DISTANCE = 4;
                    for (size_t j = 1; j <= PREFETCH_DISTANCE && i + j < r.end(); ++j) {
                        __builtin_prefetch(current_nodes[i + j], 0, 1);     // Prefetch node pointer
                        __builtin_prefetch(&current_metrics[i + j], 0, 1);  // Prefetch metrics
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

                    // Skip scoring for excluded nodes (leave-one-out) and internal nodes (force-leaf)
                    bool isLeaf = node->children.empty();
                    if (nodeIndex != state.skipNodeIndex && (!state.forceLeaf || isLeaf)) {
                        // Compute scores
                        double logRawScore = nodeMetrics.getLogRawScore(state.logReadMagnitude);
                        double logCosineScore = nodeMetrics.getLogCosineScore(state.logReadMagnitude);
                        double containmentScore = nodeMetrics.getContainmentScore(state.readUniqueSeedCount);
                        double weightedContainmentScore =
                            nodeMetrics.getWeightedContainmentScore(state.weightedContainmentDenominator);
                        double logContainmentScore =
                            nodeMetrics.getLogContainmentScore(state.logContainmentDenominator);

                        // Update thread-local results (NO LOCKS - parallel performance!)
                        tls.local_result.updateLogRawScore(nodeIndex, logRawScore, node);
                        tls.local_result.updateLogCosineScore(nodeIndex, logCosineScore, node);
                        tls.local_result.updateContainmentScore(nodeIndex, containmentScore, node);
                        tls.local_result.updateWeightedContainmentScore(nodeIndex, weightedContainmentScore, node);
                        tls.local_result.updateLogContainmentScore(nodeIndex, logContainmentScore, node);
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
            logging::debug("Processed {} nodes", currentCount);
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
            offsets[i + 1] = offsets[i] + active_tls[i]->next_nodes.size();
        }
        size_t total_next = offsets.back();

        current_nodes.clear();
        current_nodes.resize(total_next);
        current_metrics.clear();
        current_metrics.resize(total_next);
        if (params.verify_scores) {
            current_seeds.clear();
            current_seeds.resize(total_next);
        }

        // Parallel copy/move from thread-local buffers to unified arrays
        tbb::parallel_for(tbb::blocked_range<size_t>(0, active_tls.size()), [&](const tbb::blocked_range<size_t>& r) {
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
    logging::debug("Merging {} thread-local results...", thread_data.size());
    auto merge_start = std::chrono::high_resolution_clock::now();

    for (auto& tls : thread_data) {
        // Merge LogRAW scores
        if (tls.local_result.bestLogRawNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogRawScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogRawScore > result.bestLogRawScore + tolerance) {
                result.bestLogRawScore = tls.local_result.bestLogRawScore;
                result.bestLogRawNodeIndex = tls.local_result.bestLogRawNodeIndex;
                result.tiedLogRawNodeIndices = tls.local_result.tiedLogRawNodeIndices;
            } else if (tls.local_result.bestLogRawScore >= result.bestLogRawScore - tolerance) {
                result.tiedLogRawNodeIndices.insert(result.tiedLogRawNodeIndices.end(),
                                                    tls.local_result.tiedLogRawNodeIndices.begin(),
                                                    tls.local_result.tiedLogRawNodeIndices.end());
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
                result.tiedLogCosineNodeIndices.insert(result.tiedLogCosineNodeIndices.end(),
                                                       tls.local_result.tiedLogCosineNodeIndices.begin(),
                                                       tls.local_result.tiedLogCosineNodeIndices.end());
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
                result.tiedContainmentNodeIndices.insert(result.tiedContainmentNodeIndices.end(),
                                                         tls.local_result.tiedContainmentNodeIndices.begin(),
                                                         tls.local_result.tiedContainmentNodeIndices.end());
            }
        }
        // Merge Weighted Containment scores
        if (tls.local_result.bestWeightedContainmentNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestWeightedContainmentScore * 0.0001, 1e-9);
            if (tls.local_result.bestWeightedContainmentScore > result.bestWeightedContainmentScore + tolerance) {
                result.bestWeightedContainmentScore = tls.local_result.bestWeightedContainmentScore;
                result.bestWeightedContainmentNodeIndex = tls.local_result.bestWeightedContainmentNodeIndex;
                result.tiedWeightedContainmentNodeIndices = tls.local_result.tiedWeightedContainmentNodeIndices;
            } else if (tls.local_result.bestWeightedContainmentScore >=
                       result.bestWeightedContainmentScore - tolerance) {
                result.tiedWeightedContainmentNodeIndices.insert(
                    result.tiedWeightedContainmentNodeIndices.end(),
                    tls.local_result.tiedWeightedContainmentNodeIndices.begin(),
                    tls.local_result.tiedWeightedContainmentNodeIndices.end());
            }
        }
        // Merge Log Containment scores
        if (tls.local_result.bestLogContainmentNodeIndex != UINT32_MAX) {
            double tolerance = std::max(result.bestLogContainmentScore * 0.0001, 1e-9);
            if (tls.local_result.bestLogContainmentScore > result.bestLogContainmentScore + tolerance) {
                result.bestLogContainmentScore = tls.local_result.bestLogContainmentScore;
                result.bestLogContainmentNodeIndex = tls.local_result.bestLogContainmentNodeIndex;
                result.tiedLogContainmentNodeIndices = tls.local_result.tiedLogContainmentNodeIndices;
            } else if (tls.local_result.bestLogContainmentScore >= result.bestLogContainmentScore - tolerance) {
                result.tiedLogContainmentNodeIndices.insert(result.tiedLogContainmentNodeIndices.end(),
                                                            tls.local_result.tiedLogContainmentNodeIndices.begin(),
                                                            tls.local_result.tiedLogContainmentNodeIndices.end());
            }
        }
    }

    auto merge_end = std::chrono::high_resolution_clock::now();
    auto merge_ms = std::chrono::duration_cast<std::chrono::milliseconds>(merge_end - merge_start).count();
    logging::debug("Thread-local merge completed in {}ms", merge_ms);
}

void placeLite(PlacementResult& result,
               panmapUtils::LiteTree* liteTree,
               ::capnp::MessageReader& liteIndex,
               const std::string& reads1,
               const std::string& reads2,
               std::string& outputPath,
               const TraversalParams& callerParams,
               panmanUtils::Tree* fullTree) {
    auto placement_total_start = std::chrono::high_resolution_clock::now();
    logging::debug("Starting lite-index placement");

    uint32_t skipNodeIndex = UINT32_MAX;

    if (callerParams.verify_scores) {
        logging::info("VERIFICATION MODE: Will recompute all scores from scratch at each node");
        if (fullTree == nullptr) {
            logging::err("VERIFICATION MODE requires fullTree pointer but got nullptr!");
            throw std::runtime_error("Verification mode requires full Tree to be loaded");
        }
    }

    // State initialization moved to after parameter setup

    auto time_hash_delta_start = std::chrono::high_resolution_clock::now();
    size_t totalHashDeltas = 0;

    auto indexRoot = liteIndex.getRoot<LiteIndex>();

    if (!liteTree->seedChangesLoaded) {
        logging::debug("Loading pre-computed hash deltas from index for {} nodes...", liteTree->allLiteNodes.size());

        // OPTIMIZATION: Process nodes in DFS order with direct vector indexing for O(1) access
        // This eliminates hash map lookup overhead in the hot path
        const size_t numNodes = liteTree->allLiteNodes.size();

        // Parallel Pass 1: Count seed changes per node to calculate offsets
        std::vector<size_t> nodeOffsets(numNodes + 1, 0);

        std::atomic<size_t> nodesWithChanges(0);

        {
            if (!indexRoot.hasSeedChangeHashes() || !indexRoot.hasSeedChangeParentCounts() ||
                !indexRoot.hasSeedChangeChildCounts() || !indexRoot.hasNodeChangeOffsets()) {
                throw std::runtime_error(
                    "Index missing required V3 fields (seedChangeHashes, etc). V2 is no longer supported.");
            }

            logging::debug("Using VERSION 3 struct-of-arrays format from index (optimal parallel access)");
            auto indexOffsets = indexRoot.getNodeChangeOffsets();
            if (indexOffsets.size() >= numNodes + 1) {
                for (size_t i = 0; i <= numNodes; ++i) {
                    nodeOffsets[i] = indexOffsets[i];
                }
            } else {
                throw std::runtime_error(fmt::format(
                    "Struct-of-arrays format offsets size mismatch: {} vs {}", indexOffsets.size(), numNodes + 1));
            }
            totalHashDeltas = nodeOffsets[numNodes];

            // Allocate backing storage once
            liteTree->allSeedChanges.resize(totalHashDeltas);

            // Parallel Pass 2: Populate data and assign spans

            // VERSION 3: Struct-of-arrays format (with optional overflow arrays for large indices)
            auto hashesReader = indexRoot.getSeedChangeHashes();
            auto parentCountsReader = indexRoot.getSeedChangeParentCounts();
            auto childCountsReader = indexRoot.getSeedChangeChildCounts();

            size_t numSegs = hashesReader.size();
            logging::debug("Struct-of-arrays format: {} hashes in {} segments", totalHashDeltas, numSegs);

            // Phase 2a
            constexpr size_t GRAIN_SIZE = 262144;  // 256K items per chunk
            constexpr size_t CAPNP_SPLIT = 500'000'000;
            for (uint32_t seg = 0; seg < numSegs; seg++) {
                size_t segStart = static_cast<size_t>(seg) * CAPNP_SPLIT;
                size_t segEnd = std::min(segStart + CAPNP_SPLIT, totalHashDeltas);
                auto hashes = hashesReader[seg];
                auto parents = parentCountsReader[seg];
                auto children = childCountsReader[seg];
                tbb::parallel_for(tbb::blocked_range<size_t>(segStart, segEnd, GRAIN_SIZE),
                                  [&](const tbb::blocked_range<size_t>& range) {
                                      for (size_t i = range.begin(); i < range.end(); ++i) {
                                          uint32_t segOffset = i - segStart;
                                          liteTree->allSeedChanges[i] = std::make_tuple(
                                              hashes[segOffset], parents[segOffset], children[segOffset]);
                                      }
                                  });
            }

            // Phase 2b: Assign spans to nodes
            tbb::parallel_for(tbb::blocked_range<size_t>(0, numNodes), [&](const tbb::blocked_range<size_t>& range) {
                size_t localNodesWithChanges = 0;
                for (size_t dfsIndex = range.begin(); dfsIndex < range.end(); ++dfsIndex) {
                    auto* liteNode = liteTree->dfsIndexToNode[dfsIndex];
                    size_t offset = nodeOffsets[dfsIndex];
                    size_t size = nodeOffsets[dfsIndex + 1] - offset;

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

        logging::debug("Loaded {} hash deltas from index", totalHashDeltas);
        liteTree->seedChangesLoaded = true;
    } else {
        logging::info("Seed changes already loaded, skipping index deserialization");
    }

    auto time_hash_delta_end = std::chrono::high_resolution_clock::now();

    // Set traversal parameters
    TraversalParams params = callerParams;
    // Override k/s/t/l/open/hpc from the index (authoritative source)
    params.k = indexRoot.getK();
    params.s = indexRoot.getS();
    params.t = indexRoot.getT();
    params.l = indexRoot.getL();
    params.open = indexRoot.getOpen();
    params.hpc = indexRoot.getHpc();

    if (params.dedupReads) {
        logging::info("Read deduplication enabled (--dedup): counting each unique sequence once");
    }
    if (params.trimStart > 0 || params.trimEnd > 0) {
        logging::info("Read trimming enabled: params.trimStart={}, params.trimEnd={} (for primer removal)",
                      params.trimStart,
                      params.trimEnd);
    }
    if (params.pairFilter && reads2.empty()) {
        logging::info("Paired-end scoring requested but no R2 reads provided - disabled");
    }

    std::vector<std::string> allReadSequences;
    std::vector<std::string> allReadQualities;  // Quality strings (only populated if params.minSeedQuality > 0)

    // Initialize placement global state
    PlacementGlobalState state;
    state.kmerSize = indexRoot.getK();
    state.fullTree = fullTree;  // Store full tree pointer for verification mode
    state.liteNodes = liteTree->dfsIndexToNode;
    state.skipNodeIndex = skipNodeIndex;  // Set skip node for leave-one-out validation
    state.forceLeaf = params.forceLeaf;   // Restrict scoring to leaf nodes only
    if (params.forceLeaf) {
        logging::info("Force-leaf mode: only leaf nodes will be considered for placement");
    }

    auto time_read_processing_start = std::chrono::high_resolution_clock::now();
    {
        if (!reads1.empty()) {
            if (params.minSeedQuality > 0) {
                // Need quality strings for per-seed quality filtering
                std::vector<std::string> readNames;  // Unused but required by function
                extractFullFastqData(reads1, reads2, allReadSequences, allReadQualities, readNames);
                logging::info("Loaded {} reads with quality strings for Q{} filtering",
                              allReadSequences.size(),
                              params.minSeedQuality);
            } else {
                extractReadSequences(reads1, reads2, allReadSequences);
            }
        }

        // Apply homopolymer compression to reads if the index was built with HPC
        if (params.hpc && !allReadSequences.empty()) {
            logging::info("HPC mode: compressing read sequences");
            tbb::parallel_for(tbb::blocked_range<size_t>(0, allReadSequences.size()),
                              [&](const tbb::blocked_range<size_t>& range) {
                                  for (size_t i = range.begin(); i < range.end(); ++i) {
                                      allReadSequences[i] = seeding::hpcCompress(allReadSequences[i]);
                                  }
                              });
        }

        uint32_t l = indexRoot.getL();
        logging::debug("Index parameters: k={}, s={}, t={}, l={}{}",
                      params.k,
                      params.s,
                      params.t,
                      l,
                      params.hpc ? ", hpc=on" : "");

        if (l == 0) {
            logging::info("l=0: Using raw syncmers instead of k-minimizers");

            // Quality-filtered seed extraction path (no read deduplication since quality varies)
            if (params.minSeedQuality > 0 && !allReadQualities.empty()) {
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

                tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, allReadSequences.size()),
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

                            const auto& syncmers =
                                seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                            for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                                if (!isSyncmer) [[unlikely]]
                                    continue;

                                // Trim filter: skip seeds that start in trimmed regions
                                if (static_cast<int>(startPos) < validStart || static_cast<int>(startPos) > validEnd) {
                                    ++trimFiltered;
                                    continue;  // Skip seed in primer region
                                }

                                // Quality filter: check average quality over k-mer span
                                double avgQual = avgPhredQuality(qual, startPos, params.k);
                                if (avgQual < static_cast<double>(params.minSeedQuality)) {
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
                logging::info(
                    "Seed extraction (Q{} filtered): {}ms", params.minSeedQuality, duration_seed_extract.count());
                logging::info("Extracted {} unique syncmers from {} reads (filtered {} low-quality seeds)",
                              state.seedFreqInReads.size(),
                              allReadSequences.size(),
                              totalFiltered);

            } else {
                // Original path: deduplicate reads first for efficiency
                std::vector<std::pair<std::string, size_t>> sortedReads;
                auto time_dedup_start = std::chrono::high_resolution_clock::now();
                {
                    sortedReads.reserve(allReadSequences.size());
                    for (size_t i = 0; i < allReadSequences.size(); i++) {
                        sortedReads.emplace_back(allReadSequences[i], i);
                    }
                    tbb::parallel_sort(sortedReads.begin(), sortedReads.end(), [](const auto& a, const auto& b) {
                        return a.first < b.first;
                    });
                }
                auto time_dedup_end = std::chrono::high_resolution_clock::now();
                auto duration_dedup =
                    std::chrono::duration_cast<std::chrono::milliseconds>(time_dedup_end - time_dedup_start);
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
                                              is_new_seq[i] =
                                                  (sortedReads[i].first != sortedReads[i - 1].first) ? 1 : 0;
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
                                              size_t end = unique_starts[i + 1];

                                              std::vector<size_t> duplicates;
                                              duplicates.reserve(end - start);
                                              for (size_t j = start; j < end; ++j) {
                                                  duplicates.push_back(sortedReads[j].second);
                                              }
                                              dupReadsIndex[i] =
                                                  std::make_pair(sortedReads[start].first, std::move(duplicates));
                                          }
                                      });
                }

                logging::debug(
                    "Deduplication: {} reads → {} unique sequences", allReadSequences.size(), dupReadsIndex.size());

                auto time_seed_extract_start = std::chrono::high_resolution_clock::now();
                {
                    size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
                    std::vector<absl::flat_hash_map<size_t, int64_t>> threadLocalMaps(num_cpus);

                    // Pre-allocate capacity based on empirical analysis to eliminate malloc overhead
                    // Profiler shows 14% time in malloc (7.2% consolidate + 6.9% perturb)
                    // Formula: estimate seeds/read * reads/thread * load_factor_headroom
                    const size_t avg_seeds_per_read = 50;  // Conservative estimate for syncmers
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
                                              const int64_t multiplicity =
                                                  params.dedupReads ? 1 : static_cast<int64_t>(duplicates.size());

                                              // Calculate trim boundaries for this read
                                              const int seqLen = static_cast<int>(seq.size());
                                              const int trimStartBp = params.trimStart;
                                              const int trimEndBp = params.trimEnd;
                                              const int validStart = trimStartBp;
                                              const int validEnd = seqLen - trimEndBp - params.k;

                                              const auto& syncmers = seeding::rollingSyncmers(
                                                  seq, params.k, params.s, params.open, params.t, false);
                                              for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                                                  if (!isSyncmer) [[unlikely]]
                                                      continue;

                                                  // Trim filter: skip seeds that start in trimmed regions
                                                  if (static_cast<int>(startPos) < validStart ||
                                                      static_cast<int>(startPos) > validEnd) {
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
                              state.seedFreqInReads.size(),
                              allReadSequences.size(),
                              dupReadsIndex.size());
            }  // end of else block for non-quality-filtered path

        } else {
            // l > 0: Use k-minimizers (consecutive syncmers combined)

            // Quality-filtered k-minimizer extraction (no read deduplication)
            if (params.minSeedQuality > 0 && !allReadQualities.empty()) {
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

                tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, allReadSequences.size()),
                    [&](const tbb::blocked_range<size_t>& range) {
                        size_t threadId = tbb::this_task_arena::current_thread_index();
                        auto& localMap = threadLocalMaps[threadId];
                        size_t filtered = 0;

                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            const std::string& seq = allReadSequences[i];
                            const std::string& qual = allReadQualities[i];

                            const auto& syncmers =
                                seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);

                            if (syncmers.size() < static_cast<size_t>(params.l)) continue;

                            // For l=1, just use syncmer hashes with quality filtering
                            if (params.l == 1) {
                                for (const auto& [seedHash, isReverse, isSyncmer, startPos] : syncmers) {
                                    // Quality filter
                                    double avgQual = avgPhredQuality(qual, startPos, params.k);
                                    if (avgQual < static_cast<double>(params.minSeedQuality)) {
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
                                syncmerPassesQuality[j] = (avgQual >= static_cast<double>(params.minSeedQuality));
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
                                forwardRolledHash =
                                    seeding::rol(forwardRolledHash, params.k) ^ std::get<0>(syncmers[j]);
                                reverseRolledHash =
                                    seeding::rol(reverseRolledHash, params.k) ^ std::get<0>(syncmers[params.l - j - 1]);
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
                                    forwardRolledHash =
                                        seeding::rol(forwardRolledHash, params.k) ^ std::get<0>(syncmers[j + w]);
                                    reverseRolledHash = seeding::rol(reverseRolledHash, params.k) ^
                                                        std::get<0>(syncmers[j + params.l - w - 1]);
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
                auto duration_kminimizer =
                    std::chrono::duration_cast<std::chrono::milliseconds>(time_kminimizer_end - time_kminimizer_start);
                logging::debug(
                    "K-minimizer extraction (Q{} filtered): {}ms", params.minSeedQuality, duration_kminimizer.count());
                logging::debug("Extracted {} unique k-minimizers from {} reads (filtered {} low-quality seeds)",
                              state.seedFreqInReads.size(),
                              allReadSequences.size(),
                              totalFiltered);

            } else {
                // Original path: deduplicate reads first for efficiency
                auto time_kminimizer_start = std::chrono::high_resolution_clock::now();

                // Deduplicate reads first (same as l=0 case)
                std::vector<std::pair<std::string, size_t>> sortedReads;
                sortedReads.reserve(allReadSequences.size());
                for (size_t i = 0; i < allReadSequences.size(); i++) {
                    sortedReads.emplace_back(allReadSequences[i], i);
                }
                tbb::parallel_sort(sortedReads.begin(), sortedReads.end(), [](const auto& a, const auto& b) {
                    return a.first < b.first;
                });

                std::vector<std::pair<std::string, std::vector<size_t>>> dupReadsIndex;
                if (!sortedReads.empty()) {
                    std::vector<uint8_t> is_new_seq(sortedReads.size());
                    is_new_seq[0] = 1;
                    tbb::parallel_for(tbb::blocked_range<size_t>(1, sortedReads.size()),
                                      [&](const tbb::blocked_range<size_t>& range) {
                                          for (size_t i = range.begin(); i < range.end(); ++i) {
                                              is_new_seq[i] =
                                                  (sortedReads[i].first != sortedReads[i - 1].first) ? 1 : 0;
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
                                              size_t end = unique_starts[i + 1];
                                              std::vector<size_t> duplicates;
                                              duplicates.reserve(end - start);
                                              for (size_t j = start; j < end; ++j) {
                                                  duplicates.push_back(sortedReads[j].second);
                                              }
                                              dupReadsIndex[i] =
                                                  std::make_pair(sortedReads[start].first, std::move(duplicates));
                                          }
                                      });
                }

                logging::debug(
                    "Deduplication: {} reads → {} unique sequences", allReadSequences.size(), dupReadsIndex.size());

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

                tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, dupReadsIndex.size()), [&](const tbb::blocked_range<size_t>& range) {
                        size_t threadId = tbb::this_task_arena::current_thread_index();
                        auto& localMap = threadLocalMaps[threadId];

                        for (size_t i = range.begin(); i < range.end(); ++i) {
                            const auto& [seq, duplicates] = dupReadsIndex[i];
                            // If dedupReads is true, count each unique sequence once
                            const int64_t multiplicity =
                                params.dedupReads ? 1 : static_cast<int64_t>(duplicates.size());

                            const auto& syncmers =
                                seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);

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
                                forwardRolledHash =
                                    seeding::rol(forwardRolledHash, params.k) ^ std::get<0>(syncmers[j]);
                                reverseRolledHash =
                                    seeding::rol(reverseRolledHash, params.k) ^ std::get<0>(syncmers[params.l - j - 1]);
                            }

                            if (forwardRolledHash != reverseRolledHash) {
                                size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
                                localMap[minHash] += multiplicity;
                            }

                            // Rest of k-min-mers
                            for (size_t j = 1; j < syncmers.size() - params.l + 1; ++j) {
                                const size_t& prevSyncmerHash = std::get<0>(syncmers[j - 1]);
                                const size_t& nextSyncmerHash = std::get<0>(syncmers[j + params.l - 1]);

                                forwardRolledHash = seeding::rol(forwardRolledHash, params.k) ^
                                                    seeding::rol(prevSyncmerHash, params.k * params.l) ^
                                                    nextSyncmerHash;
                                reverseRolledHash = seeding::ror(reverseRolledHash, params.k) ^
                                                    seeding::ror(prevSyncmerHash, params.k) ^
                                                    seeding::rol(nextSyncmerHash, params.k * (params.l - 1));

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
                auto duration_kminimizer =
                    std::chrono::duration_cast<std::chrono::milliseconds>(time_kminimizer_end - time_kminimizer_start);
                logging::debug("K-minimizer extraction: {}ms", duration_kminimizer.count());
                logging::debug("Extracted {} unique k-minimizers from {} reads",
                              state.seedFreqInReads.size(),
                              allReadSequences.size());
            }  // end of else block for non-quality-filtered path (l > 0)
        }
    }
    auto time_read_processing_end = std::chrono::high_resolution_clock::now();
    auto duration_read_processing =
        std::chrono::duration_cast<std::chrono::milliseconds>(time_read_processing_end - time_read_processing_start);
    logging::debug("Total read processing: {}ms", duration_read_processing.count());

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
        logging::debug("=== Seed Frequency Analysis ===");
        logging::debug("  Total reads: {}, Unique seeds: {}", totalReads, uniqueSeeds);
        logging::debug("  Max seed frequency: {:.2f}%", maxFreq * 100.0);
        logging::debug("  Seeds >50%% of reads: {}", above50pct);
        logging::debug("  Seeds >20%% of reads: {}", above20pct);
        logging::debug("  Seeds >10%% of reads: {}", above10pct);
        logging::debug("  Seeds >5%% of reads: {}", above5pct);
        logging::debug("  Seeds >1%% of reads: {}", above1pct);
        logging::debug("===============================");

        // Only sort if we need to mask seeds or output diagnostics (slow for large datasets)
        std::vector<std::pair<size_t, int64_t>> sortedSeeds;
        bool needSorting = (params.seedMaskFraction > 0.0) || params.store_diagnostics;

        if (needSorting) {
            sortedSeeds.reserve(uniqueSeeds);
            for (const auto& [hash, count] : state.seedFreqInReads) {
                sortedSeeds.emplace_back(hash, count);
            }
            std::sort(sortedSeeds.begin(), sortedSeeds.end(), [](const auto& a, const auto& b) {
                return a.second > b.second;
            });

            // Show top 5 most frequent seeds if any are above 5%
            if (above5pct > 0) {
                logging::debug("  Top {} high-frequency seeds:", std::min(size_t(5), above5pct));
                for (size_t i = 0; i < std::min(size_t(5), sortedSeeds.size()); i++) {
                    double frac = static_cast<double>(sortedSeeds[i].second) / totalReads;
                    if (frac < 0.01) break;  // Stop if below 1%
                    logging::debug("    {:.2f}%  (hash: {:016x})", frac * 100.0, sortedSeeds[i].first);
                }
            }
        }

        // Apply percentile-based masking: remove top N% of seeds by frequency
        // This removes potential primer/adapter artifacts without needing to calibrate a threshold
        if (params.seedMaskFraction > 0.0 && !sortedSeeds.empty()) {
            const size_t numToMask = static_cast<size_t>(params.seedMaskFraction * uniqueSeeds);

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
                              numToMask,
                              params.seedMaskFraction * 100.0);
                logging::info("  Masked seed freq range: {:.2f}%% - {:.2f}%% of reads",
                              minMaskedFrac * 100.0,
                              maxMaskedFrac * 100.0);
                logging::info("  Total masked frequency: {}", maskedFrequency);
            } else {
                logging::info("  Seed masking: fraction too small to mask any seeds");
            }
        }
        logging::debug("===============================");

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
                size_t numToMask = static_cast<size_t>(params.seedMaskFraction * uniqueSeeds);
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
                    seedFreqFile << sortedSeeds[i].second << "\t" << std::fixed << std::setprecision(6) << frac << "\t"
                                 << masked << "\n";
                }
                seedFreqFile.close();
                logging::info("Seed frequencies written to: {}", seedFreqPath);
            }
        }
    }

    auto time_magnitude_start = std::chrono::high_resolution_clock::now();
    {
        // Compute log-scaled magnitudes for scoring
        state.totalReadSeedFrequency = 0;

        state.logReadCounts.reserve(state.seedFreqInReads.size());

        double logMagSquared = 0.0;
        double logCountSum = 0.0;  // Σ log(1+r_i) for log containment denominator
        size_t filteredSeedCount = 0;
        size_t lowSupportSeeds = 0;

        const int64_t minSupport = params.minReadSupport;

        for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
            state.totalReadSeedFrequency += readCount;

            if (readCount < minSupport) {
                lowSupportSeeds++;
                continue;
            }

            // Pre-compute log(1 + readCount) for each seed
            const double readCountD = static_cast<double>(readCount);
            double logCount = std::log1p(readCountD);
            state.logReadCounts[seedHash] = logCount;
            logMagSquared += logCount * logCount;
            logCountSum += logCount;

            filteredSeedCount++;
        }

        state.readUniqueSeedCount = filteredSeedCount;
        state.logReadMagnitude = std::sqrt(logMagSquared);
        state.logContainmentDenominator = logCountSum;

        if (lowSupportSeeds > 0) {
            logging::info("Filtered {} seeds with read count < {} (kept {} seeds)",
                          lowSupportSeeds,
                          minSupport,
                          filteredSeedCount);
        }
        logging::debug("Precomputed magnitude: log={:.6f}, total read seed frequency: {}, unique read seeds: {}",
                      state.logReadMagnitude,
                      state.totalReadSeedFrequency,
                      state.readUniqueSeedCount);
    }
    auto time_magnitude_end = std::chrono::high_resolution_clock::now();
    auto duration_magnitude =
        std::chrono::duration_cast<std::chrono::milliseconds>(time_magnitude_end - time_magnitude_start);
    logging::debug("Read magnitude precomputation: {}ms", duration_magnitude.count());

    state.root = liteTree->root;

    // Precompute inverse genome counts from root's seed changes
    // g_i = number of genomes containing seed i (root's childCount for that seed)
    // Used by weighted containment metric to upweight rare/discriminative seeds
    if (liteTree->root) {
        for (const auto& [seedHash, parentCount, childCount] : liteTree->root->seedChanges) {
            if (childCount > 0 && state.logReadCounts.contains(seedHash)) {
                double invCount = 1.0 / static_cast<double>(childCount);
                state.seedInverseGenomeCounts[seedHash] = invCount;
                state.weightedContainmentDenominator += invCount;
            }
        }
        logging::debug("Weighted containment: {} seeds with inverse genome counts, denominator={:.6f}",
                      state.seedInverseGenomeCounts.size(),
                      state.weightedContainmentDenominator);
    }

    placement::NodeMetrics rootMetrics;
    // CRITICAL: Root starts with empty genome (all genome metrics = 0).
    // The root's seed changes will apply deltas to compute the actual root state.

    logging::debug("Starting BFS placement traversal with {} tree nodes", liteTree->allLiteNodes.size());

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
        placeLiteHelperBFS(
            root_level_nodes, root_level_metrics, root_level_seed_counts, state, result, params, nodesProcessed);

        // Lazy ID resolution: Only deserialize winning node IDs (eliminates ~1800ms overhead!)
        result.resolveNodeIds(liteTree);

        // =======================================================================
        // ALIGNMENT-BASED REFINEMENT (optional)
        // After k-mer scoring, refine top candidates by full minimap2 alignment
        // =======================================================================
        if (params.refineEnabled && state.fullTree != nullptr) {
            TraversalParams refineParams = params;  // Copy params for refinement
            bool pairedEnd = !reads2.empty();
            refineTopCandidates(liteTree, state.fullTree, allReadSequences, result, refineParams, pairedEnd);

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
    auto duration_traversal =
        std::chrono::duration_cast<std::chrono::milliseconds>(time_traversal_end - time_traversal_start);
    logging::debug("Tree traversal: {}ms", duration_traversal.count());

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

        logging::debug(
            "Stored alignment parameters: k={}, s={}, t={}, open={}", result.k, result.s, result.t, result.open);
    }

    auto time_write_start = std::chrono::high_resolution_clock::now();
    std::string placementsFilePath = outputPath;
    std::ofstream placementsFile(placementsFilePath);
    if (placementsFile.is_open()) {
        placementsFile << "metric\tscore\tnodes\n";

        auto writeMetric =
            [&](const char* name, double score, const std::vector<uint32_t>& tied, const std::string& bestId) {
                placementsFile << name << "\t" << std::fixed << std::setprecision(6) << score << "\t";
                if (!tied.empty()) {
                    for (size_t i = 0; i < tied.size(); ++i) {
                        if (i > 0) placementsFile << ",";
                        placementsFile << liteTree->resolveNodeId(tied[i]);
                    }
                } else {
                    placementsFile << bestId;
                }
                placementsFile << "\n";
            };
        writeMetric("log_raw", result.bestLogRawScore, result.tiedLogRawNodeIndices, result.bestLogRawNodeId);
        writeMetric(
            "log_cosine", result.bestLogCosineScore, result.tiedLogCosineNodeIndices, result.bestLogCosineNodeId);
        writeMetric("containment",
                    result.bestContainmentScore,
                    result.tiedContainmentNodeIndices,
                    result.bestContainmentNodeId);
        writeMetric("weighted_containment",
                    result.bestWeightedContainmentScore,
                    result.tiedWeightedContainmentNodeIndices,
                    result.bestWeightedContainmentNodeId);
        writeMetric("log_containment",
                    result.bestLogContainmentScore,
                    result.tiedLogContainmentNodeIndices,
                    result.bestLogContainmentNodeId);

        // Per-metric refinement scores (if refinement was run)
        if (result.refinementWasRun) {
            const std::pair<const char*, const PlacementResult::RefinedResult*> refined[] = {
                {"refined_log_raw", &result.refinedLogRaw},
                {"refined_log_cosine", &result.refinedLogCosine},
                {"refined_containment", &result.refinedContainment},
                {"refined_weighted_containment", &result.refinedWeightedContainment},
                {"refined_log_containment", &result.refinedLogContainment}};
            for (const auto& [name, r] : refined) {
                if (r->nodeIndex != UINT32_MAX)
                    placementsFile << name << "\t" << std::fixed << std::setprecision(0) << r->score << "\t"
                                   << r->nodeId << "\n";
            }
        }

        placementsFile.close();
        logging::debug("Wrote placement results to {}", placementsFilePath);
    } else {
        logging::err("Failed to open placements file: {}", placementsFilePath);
    }
    auto time_write_end = std::chrono::high_resolution_clock::now();
    auto duration_write = std::chrono::duration_cast<std::chrono::milliseconds>(time_write_end - time_write_start);
    logging::debug("Result file writing: {}ms", duration_write.count());

    auto placement_total_end = std::chrono::high_resolution_clock::now();
    auto duration_placement_total =
        std::chrono::duration_cast<std::chrono::milliseconds>(placement_total_end - placement_total_start);

    logging::debug("Placement complete in {}ms", duration_placement_total.count());
    logging::debug("Best LogRaw score: {:.6f} (node: {})", result.bestLogRawScore, result.bestLogRawNodeId);
    if (result.tiedLogRawNodeIndices.size() > 1) {
        logging::debug("  {} nodes tied for best LogRaw score", result.tiedLogRawNodeIndices.size());
    }

    if (result.refinementWasRun) {
        logging::info("Per-metric refined placements:");
        const std::pair<const char*, const PlacementResult::RefinedResult*> refined[] = {
            {"refined_log_raw", &result.refinedLogRaw},
            {"refined_log_cosine", &result.refinedLogCosine},
            {"refined_containment", &result.refinedContainment},
            {"refined_weighted_containment", &result.refinedWeightedContainment},
            {"refined_log_containment", &result.refinedLogContainment}};
        for (const auto& [name, r] : refined)
            if (r->nodeIndex != UINT32_MAX) logging::info("  {}: score={:.0f} node={}", name, r->score, r->nodeId);
    }
}

}  // namespace placement

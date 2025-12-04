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

// Perfect shuffle for interleaving paired-end reads
template <typename T>
void perfect_shuffle(std::vector<T>& v) {
  const size_t n = v.size();
  if (n < 2) return;
  
  std::vector<T> canvas(n);
  for (size_t i = 0; i < n / 2; i++) {
    canvas[i*2] = std::move(v[i]);
    canvas[i*2+1] = std::move(v[i + n/2]);
  }
  v = std::move(canvas);
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
    perfect_shuffle(readSequences);
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
    perfect_shuffle(readSequences);
    perfect_shuffle(readQuals);
    perfect_shuffle(readNames);
  }
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

            // Update weighted Jaccard numerator: delta of min(read, genome)
            // Using branchless min (compiles to CMOV on x86-64)
            const int64_t oldMin = std::min(readCount, parentCount);
            const int64_t newMin = std::min(readCount, childCount);
            childMetrics.weightedJaccardNumerator += (newMin - oldMin);

            // Update cosine numerator: Σ(r_i * g_i) delta
            childMetrics.cosineNumerator += static_cast<double>(readCount) * freqDelta;
        }
    }
    
}


namespace placement {

using ::panmanUtils::Node;
using ::panmanUtils::Tree;

// Forward declarations
class PlacementResult;



/**
 * @brief Update Jaccard similarity score for a node
 * 
 * @param nodeId The node ID being evaluated
 * @param jaccardScore The Jaccard similarity score
 * @param node Pointer to the LiteNode (optional, for storing score on node)
 */
void PlacementResult::updateJaccardScore(uint32_t nodeIndex, double jaccardScore, panmapUtils::LiteNode* node) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (node) {
        node->jaccardScore = static_cast<float>(jaccardScore);
    }
    
    // Update thread-local best (no locks needed!)
    double tolerance = std::max(bestJaccardScore * 0.0001, 1e-9);
    if (jaccardScore > bestJaccardScore + tolerance) {
        bestJaccardScore = jaccardScore;
        bestJaccardNodeIndex = nodeIndex;
        tiedJaccardNodeIndices.clear();
        tiedJaccardNodeIndices.push_back(nodeIndex);
    } else if (jaccardScore >= bestJaccardScore - tolerance && jaccardScore > 0) {
        if (tiedJaccardNodeIndices.empty() || tiedJaccardNodeIndices.back() != bestJaccardNodeIndex) {
            tiedJaccardNodeIndices.push_back(bestJaccardNodeIndex);
        }
        if (nodeIndex != bestJaccardNodeIndex) {
            tiedJaccardNodeIndices.push_back(nodeIndex);
        }
    }
}

/**
 * Update the weighted Jaccard similarity score tracking
 * @param nodeId The node ID being evaluated  
 * @param weightedJaccardScore The weighted Jaccard similarity score
 * @param node Pointer to the LiteNode (optional, for storing score on node)
 */
void PlacementResult::updateWeightedJaccardScore(uint32_t nodeIndex, double weightedJaccardScore, panmapUtils::LiteNode* node) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (node) {
        node->weightedJaccardScore = static_cast<float>(weightedJaccardScore);
    }
    
    // Update thread-local best
    double tolerance = std::max(bestWeightedJaccardScore * 0.0001, 1e-9);
    if (weightedJaccardScore > bestWeightedJaccardScore + tolerance) {
        bestWeightedJaccardScore = weightedJaccardScore;
        bestWeightedJaccardNodeIndex = nodeIndex;
        tiedWeightedJaccardNodeIndices.clear();
        tiedWeightedJaccardNodeIndices.push_back(nodeIndex);
    } else if (weightedJaccardScore >= bestWeightedJaccardScore - tolerance && weightedJaccardScore > 0) {
        if (tiedWeightedJaccardNodeIndices.empty() || tiedWeightedJaccardNodeIndices.back() != bestWeightedJaccardNodeIndex) {
            tiedWeightedJaccardNodeIndices.push_back(bestWeightedJaccardNodeIndex);
        }
        if (nodeIndex != bestWeightedJaccardNodeIndex) {
            tiedWeightedJaccardNodeIndices.push_back(nodeIndex);
        }
    }
}

/**
 * @brief Update Raw Seed Match score for a node
 * 
 * @param nodeId The node ID being evaluated
 * @param score The raw seed match score (sum of read frequencies for matched seeds)
 * @param node Pointer to the LiteNode (optional, for storing score on node)
 */
void PlacementResult::updateRawSeedMatchScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (node) {
        node->rawSeedMatchScore = static_cast<float>(score);
    }
    
    // Update thread-local best with tolerance
    double tolerance = std::max(bestRawSeedMatchScore * 0.0001, 1e-9);
    if (score > bestRawSeedMatchScore + tolerance) {
        bestRawSeedMatchScore = score;
        bestRawSeedMatchNodeIndex = nodeIndex;
        tiedRawSeedMatchNodeIndices.clear();
        tiedRawSeedMatchNodeIndices.push_back(nodeIndex);
    } else if (score >= bestRawSeedMatchScore - tolerance && score > 0) {
        if (tiedRawSeedMatchNodeIndices.empty() || tiedRawSeedMatchNodeIndices.back() != bestRawSeedMatchNodeIndex) {
            tiedRawSeedMatchNodeIndices.push_back(bestRawSeedMatchNodeIndex);
        }
        if (nodeIndex != bestRawSeedMatchNodeIndex) {
            tiedRawSeedMatchNodeIndices.push_back(nodeIndex);
        }
    }
}

/**
 * @brief Update hits score for a node
 * 
 * @param nodeId The node ID being evaluated
 * @param hits The number of hits for this node
 * @param node Pointer to the LiteNode (optional, for storing score on node)
 */
void PlacementResult::updateHitsScore(uint32_t nodeIndex, int64_t hits, panmapUtils::LiteNode* node) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (node) {
        node->hitsScore = static_cast<float>(hits);
    }
    
    // Update thread-local best
    if (hits > maxHitsInAnyGenome) {
        maxHitsInAnyGenome = hits;
        maxHitsNodeIndex = nodeIndex;
        tiedMaxHitsNodeIndices.clear();
        tiedMaxHitsNodeIndices.push_back(nodeIndex);
    } else if (hits == maxHitsInAnyGenome && hits > 0) {
        if (tiedMaxHitsNodeIndices.empty() || tiedMaxHitsNodeIndices.back() != maxHitsNodeIndex) {
            tiedMaxHitsNodeIndices.push_back(maxHitsNodeIndex);
        }
        if (nodeIndex != maxHitsNodeIndex) {
            tiedMaxHitsNodeIndices.push_back(nodeIndex);
        }
    }
}

/**
 * @brief Update Jaccard (Presence/Absence) score for a node
 * 
 * @param nodeId The node ID being evaluated
 * @param score The Jaccard score based on presence/absence of seeds
 * @param node Pointer to the LiteNode (optional, for storing score on node)
 */
void PlacementResult::updateJaccardPresenceScore(uint32_t nodeIndex, double score, panmapUtils::LiteNode* node) {
    // Store score on node first
    if (node) {
        node->jaccardPresenceScore = static_cast<float>(score);
    }
    
    // Early exit for non-competitive scores
    if (score <= bestJaccardPresenceScore) [[likely]] return;
    
    // Only update if score is better (simple best tracking, no mutex needed)
    if (score > bestJaccardPresenceScore) {
        bestJaccardPresenceScore = score;
        bestJaccardPresenceNodeIndex = nodeIndex;
    }
}

/**
 * @brief Update cosine similarity score for a node
 * 
 * @param nodeId The node ID being evaluated
 * @param cosineScore The cosine similarity score
 * @param node Pointer to the LiteNode (optional, for storing score on node)
 */
void PlacementResult::updateCosineScore(uint32_t nodeIndex, double cosineScore, panmapUtils::LiteNode* node) {
    // LOCK-FREE: Store score on node + update thread-local best
    if (node) {
        node->cosineScore = static_cast<float>(cosineScore);
    }
    
    // Update thread-local best
    double tolerance = std::max(bestCosineScore * 0.0001, 1e-9);
    if (cosineScore > bestCosineScore + tolerance) {
        bestCosineScore = cosineScore;
        bestCosineNodeIndex = nodeIndex;
        tiedCosineNodeIndices.clear();
        tiedCosineNodeIndices.push_back(nodeIndex);
    } else if (cosineScore >= bestCosineScore - tolerance && cosineScore > 0) {
        if (tiedCosineNodeIndices.empty() || tiedCosineNodeIndices.back() != bestCosineNodeIndex) {
            tiedCosineNodeIndices.push_back(bestCosineNodeIndex);
        }
        if (nodeIndex != bestCosineNodeIndex) {
            tiedCosineNodeIndices.push_back(nodeIndex);
        }
    }
}

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

                double jaccardScore = nodeMetrics.getJaccardScore(state.readUniqueSeedCount);
                double weightedJaccardScore = nodeMetrics.getWeightedJaccardScore(state.totalReadSeedFrequency);
                double cosineScore = nodeMetrics.getCosineScore(state.readMagnitude);
                double presenceJaccardScore = nodeMetrics.getPresenceJaccardScore(state.readUniqueSeedCount);
                
                // Update thread-local results (NO LOCKS - parallel performance!)
                tls.local_result.updateJaccardScore(nodeIndex, jaccardScore, node);
                tls.local_result.updateWeightedJaccardScore(nodeIndex, weightedJaccardScore, node);
                tls.local_result.updateCosineScore(nodeIndex, cosineScore, node);
                tls.local_result.updateJaccardPresenceScore(nodeIndex, presenceJaccardScore, node);
                // Raw score: Σ(readCount / genomeCount) - rewards unique hits
                tls.local_result.updateRawSeedMatchScore(nodeIndex, nodeMetrics.rawMatchScore, node);
                // Hits score uses presence-based intersection count (unique seeds matched)
                tls.local_result.updateHitsScore(nodeIndex, nodeMetrics.presenceIntersectionCount, node);

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
                       panmanUtils::Tree *fullTree) {
    auto placement_total_start = std::chrono::high_resolution_clock::now();
    logging::info("Starting lite-index placement");
    
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
    
    std::vector<std::string> allReadSequences;
    
    // Initialize placement global state
    PlacementGlobalState state;
    state.kmerSize = indexRoot.getK();
    state.fullTree = fullTree;  // Store full tree pointer for verification mode
    state.liteNodes = liteTree->dfsIndexToNode;

    auto time_read_processing_start = std::chrono::high_resolution_clock::now();
    {
        if (!reads1.empty()) {
            extractReadSequences(reads1, reads2, allReadSequences);
        }
        
        uint32_t l = indexRoot.getL();
        logging::info("Index parameters: k={}, s={}, t={}, l={}", 
                    params.k, params.s, params.t, l);
        
        if (l == 0) {
            logging::info("l=0: Using raw syncmers instead of k-minimizers");
            
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
                            const int64_t multiplicity = static_cast<int64_t>(duplicates.size());
                            
                            const auto& syncmers = seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                            for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                                if (!isSyncmer) [[unlikely]] continue;
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
            
        } else {
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
                        const int64_t multiplicity = static_cast<int64_t>(duplicates.size());
                        
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
        }
    }
    auto time_read_processing_end = std::chrono::high_resolution_clock::now();
    auto duration_read_processing = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_read_processing_end - time_read_processing_start);
    logging::info("[PLACEMENT TIMING]   Total read processing: {}ms", duration_read_processing.count());
    
    auto time_magnitude_start = std::chrono::high_resolution_clock::now();
    {

        state.totalReadSeedFrequency = 0;
        state.readUniqueSeedCount = state.seedFreqInReads.size();
        
        // Compute read magnitude and total frequency
        // Note: seedFreqInReads map serves as both frequency storage and O(1) membership test
        for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
            state.readMagnitude += static_cast<double>(readCount * readCount);
            state.totalReadSeedFrequency += readCount;
        }
        state.readMagnitude = std::sqrt(state.readMagnitude);
        logging::info("Precomputed read magnitude: {:.6f}, total read seed frequency: {}, unique read seeds: {}", 
                    state.readMagnitude, state.totalReadSeedFrequency, state.readUniqueSeedCount);
    }
    auto time_magnitude_end = std::chrono::high_resolution_clock::now();
    auto duration_magnitude = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_magnitude_end - time_magnitude_start);
    logging::info("[PLACEMENT TIMING]   Read magnitude precomputation: {}ms", duration_magnitude.count());
    
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
}

} // namespace placement


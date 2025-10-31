#include "panman.hpp"
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <iomanip>  // For std::fixed and std::setprecision
#include <algorithm> // Added for std::min and std::transform
#include <sstream>  // For std::ostringstream
#include "placement.hpp"
#include "mgsr.hpp"
#include "seeding.hpp"
#include "caller_logging.hpp"
#include "index.capnp.h"
#include "mgsr_index.capnp.h"
#include <capnp/common.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <capnp/serialize.h>
#include "logging.hpp"
#include "state.hpp"
#include "seeding.hpp"
#include "gap_map.hpp"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_group.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdio>
#include <fcntl.h> // For O_RDONLY, O_WRONLY, etc.
#include <algorithm> // For std::min
#include <limits>  // Add this include for std::numeric_limits
#include <stack>   // For stack-based traversal

#include "progress_state.hpp"
#include "indexing.hpp"
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>

using namespace coordinates;
using namespace logging;
using namespace state;

// Forward declarations for MGSR-style functions  
void updateSimilarityMetricsAdd(size_t hash,
                               placement::MgsrGenomeState& genomeState,
                               const placement::PlacementGlobalState& state,
                               const std::string& nodeId = "");

void updateSimilarityMetricsDel(size_t hash,
                               placement::MgsrGenomeState& genomeState,
                               const placement::PlacementGlobalState& state,
                               const std::string& nodeId = "");

void calculateMgsrStyleScores(const std::string& nodeId,
                             const placement::MgsrGenomeState& genomeState,
                             const placement::PlacementGlobalState& state,
                             placement::PlacementResult& result,
                             const placement::TraversalParams& params);

void updateMetricsAfterBacktrack(placement::MgsrGenomeState& genomeState,
                                const placement::PlacementGlobalState& state);

void placeLiteHelper(panmapUtils::LiteNode* node,
                    placement::MgsrGenomeState& genomeState,
                    uint32_t& dfsIndex,
                    const placement::PlacementGlobalState& state,
                    placement::PlacementResult& result,
                    const placement::TraversalParams& params);

std::tuple<double, double, size_t> calculateScoresFromGenomeState(
    const std::string& nodeId,
    placement::PlacementGlobalState& state,
    placement::PlacementResult& result,
    const placement::TraversalParams& params,
    const GenomeState& genomeState);

// Implementation of GenomeState methods
void GenomeState::applyNodeChanges(const std::vector<std::pair<size_t, int64_t>>& seedDeltas) {
    for (const auto& [seedHash, countDelta] : seedDeltas) {
        if (countDelta > 0) {
            // Add seeds
            auto [it, inserted] = seedCounts.try_emplace(seedHash, 0);
            it->second += countDelta;
            uniqueSeeds.insert(seedHash);
        } else if (countDelta < 0) {
            // Remove seeds
            auto it = seedCounts.find(seedHash);
            if (it != seedCounts.end()) {
                it->second += countDelta; // countDelta is negative
                if (it->second <= 0) {
                    seedCounts.erase(it);
                    uniqueSeeds.erase(seedHash);
                }
            }
        }
    }
    metricsValid = false; // Invalidate cached metrics
}

void GenomeState::computeMetrics(const placement::PlacementGlobalState& state) const {
    if (metricsValid) return;
    
    // Reset metrics
    weightedJaccardNumerator = 0;
    jaccardNumerator = 0; 
    jaccardDenominator = 0;
    cosineNumerator = 0.0;
    cosineDenominator = 0.0;
    
    // Compute metrics by iterating through read seeds
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        auto genomeIt = seedCounts.find(seedHash);
        int64_t genomeCount = (genomeIt != seedCounts.end()) ? genomeIt->second : 0;
        
        if (genomeCount > 0) {
            // Seed is in both read and genome (intersection)
            jaccardNumerator += readCount; // For regular Jaccard: sum read frequencies in intersection
            weightedJaccardNumerator += std::min(readCount, genomeCount); // min for weighted Jaccard
            cosineNumerator += static_cast<double>(readCount) * genomeCount; // dot product
        }
        
        // For denominators (union logic)
        jaccardDenominator += readCount; // All read frequencies (will add genome-only later if needed)
    }
    
    // For cosine denominator: compute genome magnitude
    double genomeMagnitudeSquared = 0.0;
    for (const auto& [seedHash, genomeCount] : seedCounts) {
        genomeMagnitudeSquared += static_cast<double>(genomeCount) * genomeCount;
    }
    cosineDenominator = state.readMagnitude * std::sqrt(genomeMagnitudeSquared);
    
    metricsValid = true;
}

void applyNodeChangesToGenomeState(GenomeState& genomeState, 
                                  const auto& nodeChanges, 
                                  const placement::PlacementGlobalState& state) {
    // Process seed deletions first - these are POSITIONS, not seed indices
    auto deletions = nodeChanges.getSeedDeletions();
    for (uint32_t position : deletions) {
        // Look up which seed is at this position
        auto posIt = genomeState.positionToSeedIndex.find(position);
        if (posIt != genomeState.positionToSeedIndex.end()) {
            uint32_t seedIdx = posIt->second;
            if (seedIdx < state.seedInfo.size()) {
                auto seedReader = state.seedInfo[seedIdx];
                size_t seedHash = seedReader.getHash();
                
                // Remove from hash-based tracking
                auto it = genomeState.seedCounts.find(seedHash);
                if (it != genomeState.seedCounts.end()) {
                    it->second--;
                    if (it->second <= 0) {
                        genomeState.seedCounts.erase(it);
                        genomeState.uniqueSeeds.erase(seedHash);
                    }
                }
            }
            // Remove from position map
            genomeState.positionToSeedIndex.erase(posIt);
        }
    }
    
    // Process seed additions (called seedInsubIndices in Cap'n Proto schema)
    // These may be true additions OR substitutions (if position already had a seed)
    auto additions = nodeChanges.getSeedInsubIndices();
    for (uint32_t seedIdx : additions) {
        if (seedIdx < state.seedInfo.size()) {
            auto seedReader = state.seedInfo[seedIdx];
            size_t seedHash = seedReader.getHash();
            uint32_t position = seedReader.getStartPos();
            
            // Check if this position already has a seed (shouldn't after deletions, but be safe)
            auto posIt = genomeState.positionToSeedIndex.find(position);
            if (posIt != genomeState.positionToSeedIndex.end()) {
                // Substitution case: remove old seed first
                uint32_t oldSeedIdx = posIt->second;
                if (oldSeedIdx < state.seedInfo.size()) {
                    auto oldSeedReader = state.seedInfo[oldSeedIdx];
                    size_t oldSeedHash = oldSeedReader.getHash();
                    
                    auto oldIt = genomeState.seedCounts.find(oldSeedHash);
                    if (oldIt != genomeState.seedCounts.end()) {
                        oldIt->second--;
                        if (oldIt->second <= 0) {
                            genomeState.seedCounts.erase(oldIt);
                            genomeState.uniqueSeeds.erase(oldSeedHash);
                        }
                    }
                }
            }
            
            // Add new seed
            genomeState.positionToSeedIndex[position] = seedIdx;
            auto [it, inserted] = genomeState.seedCounts.try_emplace(seedHash, 0);
            it->second++;
            genomeState.uniqueSeeds.insert(seedHash);
        }
    }
    
    genomeState.metricsValid = false; // Invalidate cached metrics
}

std::tuple<double, double, size_t> calculateScoresFromGenomeState(
    const std::string& nodeId,
    placement::PlacementGlobalState& state,
    placement::PlacementResult& result,
    const placement::TraversalParams& params,
    const GenomeState& genomeState) {
    
    // Compute similarity metrics from genome state
    genomeState.computeMetrics(state);
    
    // Calculate final similarity scores
    double jaccardScore = 0.0;
    if (genomeState.jaccardDenominator > 0) {
        jaccardScore = static_cast<double>(genomeState.jaccardNumerator) / genomeState.jaccardDenominator;
    }
    
    double cosineScore = 0.0;
    if (genomeState.cosineDenominator > 0.0) {
        cosineScore = genomeState.cosineNumerator / genomeState.cosineDenominator;
    }
    
    double weightedJaccardScore = 0.0;
    // Calculate weighted Jaccard denominator (max logic)
    int64_t weightedJaccardDenominator = 0;
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        auto genomeIt = genomeState.seedCounts.find(seedHash);
        int64_t genomeCount = (genomeIt != genomeState.seedCounts.end()) ? genomeIt->second : 0;
        weightedJaccardDenominator += std::max(readCount, genomeCount);
    }
    if (weightedJaccardDenominator > 0) {
        weightedJaccardScore = static_cast<double>(genomeState.weightedJaccardNumerator) / weightedJaccardDenominator;
    }
    
    // Presence/absence Jaccard
    size_t intersectionCount = 0;
    for (const auto& [seedHash, _] : state.seedFreqInReads) {
        if (genomeState.uniqueSeeds.find(seedHash) != genomeState.uniqueSeeds.end()) {
            intersectionCount++;
        }
    }
    double presenceAbsenceScore = 0.0;
    size_t unionCount = state.seedFreqInReads.size(); // All read seeds are in union
    if (unionCount > 0) {
        presenceAbsenceScore = static_cast<double>(intersectionCount) / unionCount;
    }
    
    // Update best scores in result
    // Debug logging for ad-hoc scoring path: print components for a few nodes to detect anomalies
    static int dbgCount2 = 0;
    if (dbgCount2 < 10) {
        logging::info("DEBUG(calcScores): Node {}: jaccard={:.9f} (num={}, denom={}), weightedJaccard={:.9f} (num={}, denom={}), cosine={:.9f} (num={:.6f}, denom={:.6f}), presence={:.9f}",
                      nodeId,
                      jaccardScore, genomeState.jaccardNumerator, genomeState.jaccardDenominator,
                      weightedJaccardScore, genomeState.weightedJaccardNumerator, weightedJaccardDenominator,
                      cosineScore, genomeState.cosineNumerator, genomeState.cosineDenominator,
                      presenceAbsenceScore);

        if (genomeState.cosineNumerator > 0.0 && !(cosineScore > 0.0)) {
            logging::warn("Possible cosine anomaly (calcScores) for node {}: numerator={} denom={} => cosine={}",
                         nodeId, genomeState.cosineNumerator, genomeState.cosineDenominator, cosineScore);
        }
        dbgCount2++;
    }

    result.updateJaccardScore(nodeId, jaccardScore);
    result.updateJaccardPresenceScore(nodeId, presenceAbsenceScore);
    result.updateCosineScore(nodeId, cosineScore);
    result.updateWeightedJaccardScore(nodeId, weightedJaccardScore);
    result.updateRawSeedMatchScore(nodeId, intersectionCount);
    result.updateHitsScore(nodeId, intersectionCount);
    
    return {jaccardScore, cosineScore, genomeState.uniqueSeeds.size()};
}

// MGSR-style seed addition (matches mgsr::mgsrPlacer::addSeedAtPosition)
void addSeedToGenomeState(uint64_t seedIndex, 
                         placement::MgsrGenomeState& genomeState,
                         const placement::PlacementGlobalState& state,
                         std::vector<std::pair<uint64_t, uint8_t>>& backtrack,
                         const std::string& nodeId) {
    
    auto seedReader = state.seedInfo[seedIndex];
    uint64_t pos = seedReader.getStartPos();
    size_t hash = seedReader.getHash();
    
    // Check if this seed hash already exists in genome (at different position)
    auto existingCountIt = genomeState.currentSeedCounts.find(hash);
    bool seedHashAlreadyInGenome = (existingCountIt != genomeState.currentSeedCounts.end() && existingCountIt->second > 0);
    int64_t existingCount = seedHashAlreadyInGenome ? existingCountIt->second : 0;
    
    auto posIt = genomeState.positionMap.find(pos);
    if (posIt != genomeState.positionMap.end()) {
        // Position exists - substitution
        uint64_t oldSeedIndex = posIt->second;
        auto oldSeedReader = state.seedInfo[oldSeedIndex];
        size_t oldHash = oldSeedReader.getHash();
        
        if (seedHashAlreadyInGenome && hash != oldHash) {
            logging::critical("SUBSTITUTION at pos {}: Replacing oldHash={} (seedIdx={}) with newHash={} (seedIdx={}) which ALREADY has count={}!",
                             pos, oldHash, oldSeedIndex, hash, seedIndex, existingCount);
            logging::critical("  This means newHash appears at {} existing position(s) and now also at pos {}", existingCount, pos);
        }
        
        backtrack.emplace_back(oldSeedIndex, 3);  // SUB
        
        // CRITICAL: Remove old seed's metrics first
        updateSimilarityMetricsDel(oldHash, genomeState, state, nodeId);
        
        // Remove old seed from hash map
        auto oldHashIt = genomeState.hashToPositionMap.find(oldHash);
        if (oldHashIt != genomeState.hashToPositionMap.end()) {
            auto& positions = oldHashIt->second;
            positions.erase(std::remove(positions.begin(), positions.end(), posIt), positions.end());
            if (positions.empty()) {
                genomeState.hashToPositionMap.erase(oldHashIt);
            }
        }
        
        // Update position with new seed
        posIt->second = seedIndex;
        genomeState.hashToPositionMap[hash].push_back(posIt);
        
    } else {
        // New position - addition
        posIt = genomeState.positionMap.emplace(pos, seedIndex).first;
        backtrack.emplace_back(seedIndex, 1);  // ADD
        genomeState.hashToPositionMap[hash].push_back(posIt);
        
        // If this hash already exists in the genome at a different position,
        // we're adding a duplicate occurrence. Log this for debugging.
        if (seedHashAlreadyInGenome) {
            logging::critical("ADDITION at pos {}: Adding newHash={} (seedIdx={}) which ALREADY has count={}!",
                             pos, hash, seedIndex, existingCount);
            logging::critical("  This means newHash appears at {} existing position(s) and now also at pos {}", existingCount, pos);
        }
    }
    
    // Update similarity metrics for new seed
    // Note: If the seed hash already existed in genome, this will increment the count,
    // which is correct for handling duplicate seeds (same k-mer at multiple positions)
    updateSimilarityMetricsAdd(hash, genomeState, state, nodeId);
}

// MGSR-style seed deletion (matches mgsr::mgsrPlacer::delSeedAtPosition)
void delSeedFromGenomeState(uint64_t pos,
                           placement::MgsrGenomeState& genomeState,
                           const placement::PlacementGlobalState& state,
                           std::vector<std::pair<uint64_t, uint8_t>>& backtrack,
                           const std::string& nodeId) {
    
    auto posIt = genomeState.positionMap.find(pos);
    if (posIt == genomeState.positionMap.end()) {
        // This can happen legitimately if the position was never seeded or already deleted
        // logging::debug("Attempted to delete seed at position {} but no seed found", pos);
        return;
    }
    
    uint64_t seedIndex = posIt->second;
    auto seedReader = state.seedInfo[seedIndex];
    size_t hash = seedReader.getHash();
    
    backtrack.emplace_back(seedIndex, 2);  // DEL
    
    // Remove from hash map
    auto hashIt = genomeState.hashToPositionMap.find(hash);
    if (hashIt != genomeState.hashToPositionMap.end()) {
        auto& positions = hashIt->second;
        positions.erase(std::remove(positions.begin(), positions.end(), posIt), positions.end());
        if (positions.empty()) {
            genomeState.hashToPositionMap.erase(hashIt);
        }
    }
    
    // Remove from position map
    genomeState.positionMap.erase(posIt);
    
    // Update similarity metrics
    updateSimilarityMetricsDel(hash, genomeState, state, nodeId);
}

// Backtrack functions (for restoring parent state)
void addSeedToGenomeStateBacktrack(uint64_t seedIndex, 
                                  placement::MgsrGenomeState& genomeState,
                                  const placement::PlacementGlobalState& state) {
    auto seedReader = state.seedInfo[seedIndex];
    uint64_t pos = seedReader.getStartPos();
    size_t hash = seedReader.getHash();
    
    auto posIt = genomeState.positionMap.emplace(pos, seedIndex).first;
    genomeState.hashToPositionMap[hash].push_back(posIt);
    
    updateSimilarityMetricsAdd(hash, genomeState, state);
}

void delSeedFromGenomeStateBacktrack(uint64_t pos,
                                    placement::MgsrGenomeState& genomeState,
                                    const placement::PlacementGlobalState& state) {
    auto posIt = genomeState.positionMap.find(pos);
    if (posIt == genomeState.positionMap.end()) return;
    
    uint64_t seedIndex = posIt->second;
    auto seedReader = state.seedInfo[seedIndex];
    size_t hash = seedReader.getHash();
    
    // Remove from hash map
    auto hashIt = genomeState.hashToPositionMap.find(hash);
    if (hashIt != genomeState.hashToPositionMap.end()) {
        auto& positions = hashIt->second;
        positions.erase(std::remove(positions.begin(), positions.end(), posIt), positions.end());
        if (positions.empty()) {
            genomeState.hashToPositionMap.erase(hashIt);
        }
    }
    
    genomeState.positionMap.erase(posIt);
    updateSimilarityMetricsDel(hash, genomeState, state);
}

// Update similarity metrics when adding a seed  
void updateSimilarityMetricsAdd(size_t hash,
                               placement::MgsrGenomeState& genomeState,
                               const placement::PlacementGlobalState& state,
                               const std::string& nodeId) {
    // Update seed counts
    auto [it, inserted] = genomeState.currentSeedCounts.try_emplace(hash, 0);
    int64_t oldGenomeCount = it->second;
    it->second++;
    int64_t newGenomeCount = it->second;
    
    // Debug logging for specific nodes
    double oldWJaccNum = genomeState.currentWeightedJaccardNumerator;
    double oldCosNum = genomeState.currentCosineNumerator;
    
    if ((nodeId == "node_15" || nodeId == "node_2") && oldGenomeCount > 0) {
        logging::debug("    updateSimilarityMetricsAdd[{}]: hash={}, oldGenomeCount={} -> {}", 
                         nodeId, hash, oldGenomeCount, newGenomeCount);
    }
    
    bool wasInGenome = (oldGenomeCount > 0);
    bool nowInGenome = (newGenomeCount > 0);
    
    if (!wasInGenome && nowInGenome) {
        genomeState.currentUniqueSeeds.insert(hash);
    }
    
    // Check if seed is in reads
    auto readIt = state.seedFreqInReads.find(hash);
    bool inReads = (readIt != state.seedFreqInReads.end());
    int64_t readCount = inReads ? readIt->second : 0;
    
    if (inReads) {
        // Seed is in both reads and genome
        
        // Jaccard numerator: if becoming present, add read frequency
        if (!wasInGenome && nowInGenome) {
            genomeState.currentJaccardNumerator += readCount;
            genomeState.currentPresenceIntersectionCount++;
        }
        
        // Jaccard denominator: update with genome count delta (read part is constant)
        // Denominator = total_read_freq + genome_only_freq
        // When seed is in both, adding genome doesn't change denominator (already counted in reads)
        // So no change to denominator for read-matching seeds
        
        // Weighted Jaccard numerator: update min(read, genome)
        int64_t oldMin = std::min(readCount, oldGenomeCount);
        int64_t newMin = std::min(readCount, newGenomeCount);
        int64_t wJaccNumDelta = (newMin - oldMin);
        genomeState.currentWeightedJaccardNumerator += wJaccNumDelta;
        
        // Weighted Jaccard denominator: update sum of max(read, genome)
        int64_t oldMax = (oldGenomeCount == 0) ? readCount : std::max(readCount, oldGenomeCount);
        int64_t newMax = std::max(readCount, newGenomeCount);
        genomeState.currentWeightedJaccardDenominator += (newMax - oldMax);
        
        // Cosine: update dot product
        double oldDot = static_cast<double>(readCount * oldGenomeCount);
        double newDot = static_cast<double>(readCount * newGenomeCount);
        double cosNumDelta = (newDot - oldDot);
        genomeState.currentCosineNumerator += cosNumDelta;
        
        // Debug logging for node_15
        if (nodeId == "node_15") {
            logging::debug("      [ADD-READ] readCount={}, oldMin={}, newMin={}, wJaccDelta={} ({} -> {}), oldDot={}, newDot={}, cosDelta={} ({} -> {})",
                         readCount, oldMin, newMin, wJaccNumDelta, 
                         oldWJaccNum, genomeState.currentWeightedJaccardNumerator,
                         oldDot, newDot, cosNumDelta,
                         oldCosNum, genomeState.currentCosineNumerator);
        }
        
        // Presence union: if first time in genome, no change (already in reads)
        // Union stays same when seed transitions from read-only to both
    } else {
        // Genome-only seed (not in reads)
        
        // Jaccard denominator: genome-only seeds contribute their frequency delta
        genomeState.currentJaccardDenominator += (newGenomeCount - oldGenomeCount);
        
        // Weighted Jaccard denominator: genome-only seeds contribute their count
        if (!wasInGenome) {
            genomeState.currentWeightedJaccardDenominator += newGenomeCount;
        } else {
            genomeState.currentWeightedJaccardDenominator += (newGenomeCount - oldGenomeCount);
        }
        
        // Presence union: if first time in genome, increment union
        if (!wasInGenome && nowInGenome) {
            genomeState.currentPresenceUnionCount++;
        }
    }
    
    // Update genome magnitude for cosine: maintain squared magnitude to avoid rounding errors
    double oldMagContribution = static_cast<double>(oldGenomeCount * oldGenomeCount);
    double newMagContribution = static_cast<double>(newGenomeCount * newGenomeCount);
    double magnitudeSquaredDelta = newMagContribution - oldMagContribution;
    
    // Update squared magnitude directly (no sqrt/square roundtrip)
    genomeState.currentGenomeMagnitudeSquared += magnitudeSquaredDelta;
    genomeState.currentGenomeMagnitude = std::sqrt(genomeState.currentGenomeMagnitudeSquared);
}

// Update similarity metrics when deleting a seed
void updateSimilarityMetricsDel(size_t hash,
                               placement::MgsrGenomeState& genomeState,
                               const placement::PlacementGlobalState& state,
                               const std::string& nodeId) {
    
    auto countIt = genomeState.currentSeedCounts.find(hash);
    if (countIt == genomeState.currentSeedCounts.end()) return;
    
    int64_t oldGenomeCount = countIt->second;
    countIt->second--;
    int64_t newGenomeCount = countIt->second;
    
    // Debug logging for specific nodes
    double oldWJaccNum = genomeState.currentWeightedJaccardNumerator;
    double oldCosNum = genomeState.currentCosineNumerator;
    
    if ((nodeId == "node_15" || nodeId == "node_2") && newGenomeCount >= 0) {
        logging::debug("    updateSimilarityMetricsDel[{}]: hash={}, oldGenomeCount={} -> {}", 
                         nodeId, hash, oldGenomeCount, newGenomeCount);
    }
    
    bool wasInGenome = (oldGenomeCount > 0);
    bool nowInGenome = (newGenomeCount > 0);
    
    if (newGenomeCount <= 0) {
        genomeState.currentSeedCounts.erase(countIt);
        genomeState.currentUniqueSeeds.erase(hash);
    }
    
    // Check if seed is in reads
    auto readIt = state.seedFreqInReads.find(hash);
    bool inReads = (readIt != state.seedFreqInReads.end());
    int64_t readCount = inReads ? readIt->second : 0;
    
    if (inReads) {
        // Seed is (or was) in both reads and genome
        
        // Jaccard numerator: if becoming absent, subtract read frequency
        if (wasInGenome && !nowInGenome) {
            genomeState.currentJaccardNumerator -= readCount;
            genomeState.currentPresenceIntersectionCount--;
        }
        
        // Jaccard denominator: no change for read-matching seeds
        // (already counted in read part of union)
        
        // Weighted Jaccard numerator: update min(read, genome)
        int64_t oldMin = std::min(readCount, oldGenomeCount);
        int64_t newMin = std::min(readCount, newGenomeCount);
        int64_t wJaccNumDelta = (newMin - oldMin);
        genomeState.currentWeightedJaccardNumerator += wJaccNumDelta;
        
        // Weighted Jaccard denominator: update sum of max(read, genome)
        int64_t oldMax = (oldGenomeCount == 0) ? readCount : std::max(readCount, oldGenomeCount);
        int64_t newMax = (newGenomeCount == 0) ? readCount : std::max(readCount, newGenomeCount);
        genomeState.currentWeightedJaccardDenominator += (newMax - oldMax);
        
        // Cosine: update dot product
        double oldDot = static_cast<double>(readCount * oldGenomeCount);
        double newDot = static_cast<double>(readCount * newGenomeCount);
        double cosNumDelta = (newDot - oldDot);
        genomeState.currentCosineNumerator += cosNumDelta;
        
        // Debug logging for node_15
        if (nodeId == "node_15") {
            logging::debug("      [DEL-READ] readCount={}, oldMin={}, newMin={}, wJaccDelta={} ({} -> {}), oldDot={}, newDot={}, cosDelta={} ({} -> {})",
                         readCount, oldMin, newMin, wJaccNumDelta,
                         oldWJaccNum, genomeState.currentWeightedJaccardNumerator,
                         oldDot, newDot, cosNumDelta,
                         oldCosNum, genomeState.currentCosineNumerator);
        }
        
        // Presence union: no change (seed still in reads even if removed from genome)
    } else {
        // Genome-only seed (not in reads)
        
        // Jaccard denominator: genome-only seeds contribute
        if (wasInGenome && !nowInGenome) {
            // Seed completely removed from genome: subtract from denominator
            genomeState.currentJaccardDenominator -= oldGenomeCount;
        } else {
            // Seed still in genome: update with frequency delta
            genomeState.currentJaccardDenominator += (newGenomeCount - oldGenomeCount);
        }
        
        // Weighted Jaccard denominator: genome-only seeds contribute their count
        if (newGenomeCount == 0) {
            genomeState.currentWeightedJaccardDenominator -= oldGenomeCount;
        } else {
            genomeState.currentWeightedJaccardDenominator += (newGenomeCount - oldGenomeCount);
        }
        
        // Presence union: if removed from genome, decrement union
        if (wasInGenome && !nowInGenome) {
            genomeState.currentPresenceUnionCount--;
        }
    }
    
    // Update genome magnitude for cosine: maintain squared magnitude to avoid rounding errors
    double oldMagContribution = static_cast<double>(oldGenomeCount * oldGenomeCount);
    double newMagContribution = static_cast<double>(newGenomeCount * newGenomeCount);
    double magnitudeSquaredDelta = newMagContribution - oldMagContribution;
    
    // Update squared magnitude directly (no sqrt/square roundtrip)
    genomeState.currentGenomeMagnitudeSquared += magnitudeSquaredDelta;
    genomeState.currentGenomeMagnitude = (genomeState.currentGenomeMagnitudeSquared > 0.0) ? 
        std::sqrt(genomeState.currentGenomeMagnitudeSquared) : 0.0;
}

// Toggle for verification mode: compute scores from scratch instead of using incremental updates
// Can be enabled via command line --verify-scores flag
void calculateMgsrStyleScores(const std::string& nodeId,
                             const placement::MgsrGenomeState& genomeState,
                             const placement::PlacementGlobalState& state,
                             placement::PlacementResult& result,
                             const placement::TraversalParams& params) {
    
    // Calculate final scores from current accumulated metrics (fully dynamic - no iteration)
    double jaccardScore = 0.0;
    double cosineScore = 0.0;
    double weightedJaccardScore = 0.0;
    double presenceJaccardScore = 0.0;
    
    // Verification path: recompute from scratch for validation
    const bool VERIFY_INCREMENTAL_SCORES = params.verify_scores;
    if (VERIFY_INCREMENTAL_SCORES) {
        // CRITICAL: Use getStringFromReference() to get the TRUE genome sequence from tree structure
        // NOT from the incrementally-maintained positionMap which may be corrupted!
        if (state.fullTree == nullptr) {
            logging::critical("VERIFICATION ERROR: fullTree is nullptr but verification is enabled!");
            std::exit(1);
        }
        
        // Get the true genome sequence from the tree
        std::string genomeSequence = state.fullTree->getStringFromReference(nodeId, false, true);
        logging::debug("VERIFICATION: Got genome sequence of length {} for node {}", 
                      genomeSequence.length(), nodeId);
        
        // Extract syncmer seeds from the TRUE genome sequence using rollingSyncmers
        // CRITICAL: Use returnAll=false to match how MGSR indexing extracts seeds
        auto syncmers = seeding::rollingSyncmers(genomeSequence, params.k, params.s, params.open, params.t, false);
        
        // Build complete seed frequency map from extracted seeds
        // With returnAll=false, all returned items are syncmers (isSyncmer is always true)
        absl::flat_hash_map<size_t, int64_t> verifyGenomeSeedCounts;
        for (const auto& [hash, isReverse, isSyncmer, endPos] : syncmers) {
            verifyGenomeSeedCounts[hash]++;
        }
        
        // Check for duplicate seeds (genomeCount > 1)
        int duplicateSeedCount = 0;
        int64_t totalDuplicateFrequency = 0;
        for (const auto& [hash, count] : verifyGenomeSeedCounts) {
            if (count > 1) {
                duplicateSeedCount++;
                totalDuplicateFrequency += count;
                if (duplicateSeedCount <= 5) {
                    auto readIt = state.seedFreqInReads.find(hash);
                    int64_t readCount = (readIt != state.seedFreqInReads.end()) ? readIt->second : 0;
                    logging::critical("    DUPLICATE SEED: hash={}, genomeCount={}, readCount={}", hash, count, readCount);
                }
            }
        }
        
        logging::debug("VERIFICATION: Extracted {} syncmers ({} unique) from true genome sequence", 
                      syncmers.size(), verifyGenomeSeedCounts.size());
        logging::critical("VERIFICATION: Found {} unique seeds with duplicates (total freq={})", 
                         duplicateSeedCount, totalDuplicateFrequency);
        
        // DEBUG: If this is node_2, compare verified seeds vs incremental seeds for the problem hashes
        if (nodeId == "node_2") {
            std::vector<size_t> problemHashes = {7376186388207194583ULL, 6864244690904188740ULL};
            for (size_t problemHash : problemHashes) {
                auto verifyIt = verifyGenomeSeedCounts.find(problemHash);
                auto incrIt = genomeState.currentSeedCounts.find(problemHash);
                int64_t verifyCount = (verifyIt != verifyGenomeSeedCounts.end()) ? verifyIt->second : 0;
                int64_t incrCount = (incrIt != genomeState.currentSeedCounts.end()) ? incrIt->second : 0;
                
                logging::critical("DEBUG PROBLEM HASH {}: verifyCount={}, incrCount={}", 
                                 problemHash, verifyCount, incrCount);
                
                // Find which positions have this hash in incremental state
                auto hashPosIt = genomeState.hashToPositionMap.find(problemHash);
                if (hashPosIt != genomeState.hashToPositionMap.end()) {
                    logging::critical("  Incremental: hash appears at {} position(s):", hashPosIt->second.size());
                    for (auto posMapIt : hashPosIt->second) {
                        logging::critical("    Position: {}, seedIndex: {}", posMapIt->first, posMapIt->second);
                    }
                }
            }
        }
        
        // Recompute all metrics from scratch
        int64_t verifyJaccardNum = 0;
        int64_t verifyJaccardDenom = 0;
        int64_t verifyWeightedJaccNum = 0;
        int64_t verifyWeightedJaccDenom = 0;
        double verifyCosineNum = 0.0;
        double verifyGenomeMagSq = 0.0;
        size_t verifyPresenceIntersection = 0;
        size_t verifyPresenceUnionGenomeOnly = 0;
        
        // Iterate through all read seeds
        for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
            auto genomeIt = verifyGenomeSeedCounts.find(seedHash);
            int64_t genomeCount = (genomeIt != verifyGenomeSeedCounts.end()) ? genomeIt->second : 0;
            
            if (genomeCount > 0) {
                // Seed in both read and genome
                verifyJaccardNum += readCount;
                verifyWeightedJaccNum += std::min(readCount, genomeCount);
                verifyWeightedJaccDenom += std::max(readCount, genomeCount);
                verifyCosineNum += static_cast<double>(readCount) * genomeCount;
                verifyPresenceIntersection++;
            } else {
                // Read-only seed
                verifyWeightedJaccDenom += readCount;
            }
        }
        
        // Add genome-only seeds to denominators
        for (const auto& [seedHash, genomeCount] : verifyGenomeSeedCounts) {
            verifyGenomeMagSq += static_cast<double>(genomeCount) * genomeCount;
            
            auto readIt = state.seedFreqInReads.find(seedHash);
            if (readIt == state.seedFreqInReads.end()) {
                // Genome-only seed
                verifyJaccardDenom += genomeCount;
                verifyWeightedJaccDenom += genomeCount;
                verifyPresenceUnionGenomeOnly++;
            }
        }
        
        verifyJaccardDenom += state.totalReadSeedFrequency;
        double verifyGenomeMag = std::sqrt(verifyGenomeMagSq);
        
        // Compute final scores
        double verifyJaccardScore = (verifyJaccardDenom > 0) ? 
            static_cast<double>(verifyJaccardNum) / verifyJaccardDenom : 0.0;
        double verifyWeightedJaccardScore = (verifyWeightedJaccDenom > 0) ? 
            static_cast<double>(verifyWeightedJaccNum) / verifyWeightedJaccDenom : 0.0;
        double verifyCosineScore = (state.readMagnitude > 0.0 && verifyGenomeMag > 0.0) ? 
            verifyCosineNum / (state.readMagnitude * verifyGenomeMag) : 0.0;
        double verifyPresenceScore = (state.seedFreqInReads.size() + verifyPresenceUnionGenomeOnly > 0) ?
            static_cast<double>(verifyPresenceIntersection) / (state.seedFreqInReads.size() + verifyPresenceUnionGenomeOnly) : 0.0;
        
        // Compare with incremental values
        int64_t incrJaccardDenom = state.totalReadSeedFrequency + genomeState.currentJaccardDenominator;
        double incrJaccardScore = (incrJaccardDenom > 0) ? 
            static_cast<double>(genomeState.currentJaccardNumerator) / incrJaccardDenom : 0.0;
        double incrWeightedJaccardScore = (genomeState.currentWeightedJaccardDenominator > 0) ? 
            static_cast<double>(genomeState.currentWeightedJaccardNumerator) / genomeState.currentWeightedJaccardDenominator : 0.0;
        double incrCosineScore = (state.readMagnitude > 0.0 && genomeState.currentGenomeMagnitude > 0.0) ? 
            genomeState.currentCosineNumerator / (state.readMagnitude * genomeState.currentGenomeMagnitude) : 0.0;
        
        // Check for discrepancies - BREAK ON FIRST MISMATCH FOR INVESTIGATION
        const double TOLERANCE = 1e-6;
        bool jacDiff = std::abs(verifyJaccardScore - incrJaccardScore) > TOLERANCE;
        bool wJacDiff = std::abs(verifyWeightedJaccardScore - incrWeightedJaccardScore) > TOLERANCE;
        bool cosDiff = std::abs(verifyCosineScore - incrCosineScore) > TOLERANCE;
        
        if (jacDiff || wJacDiff || cosDiff) {
            logging::critical("========================================");
            logging::critical("VERIFICATION MISMATCH DETECTED at node: {}", nodeId);
            logging::critical("========================================");
            
            if (jacDiff) {
                logging::critical("JACCARD MISMATCH:");
                logging::critical("  Verified:    {:.12f} (numerator={}, denominator={})", 
                                verifyJaccardScore, verifyJaccardNum, verifyJaccardDenom);
                logging::critical("  Incremental: {:.12f} (numerator={}, denominator={})",
                                incrJaccardScore, genomeState.currentJaccardNumerator, incrJaccardDenom);
                logging::critical("  Difference:  {:.12f}", std::abs(verifyJaccardScore - incrJaccardScore));
                logging::critical("  Numerator diff: {}", verifyJaccardNum - genomeState.currentJaccardNumerator);
                logging::critical("  Denominator diff: {}", verifyJaccardDenom - incrJaccardDenom);
            }
            
            if (wJacDiff) {
                logging::critical("WEIGHTED JACCARD MISMATCH:");
                logging::critical("  Verified:    {:.12f} (numerator={}, denominator={})",
                                verifyWeightedJaccardScore, verifyWeightedJaccNum, verifyWeightedJaccDenom);
                logging::critical("  Incremental: {:.12f} (numerator={}, denominator={})",
                                incrWeightedJaccardScore, genomeState.currentWeightedJaccardNumerator, 
                                genomeState.currentWeightedJaccardDenominator);
                logging::critical("  Difference:  {:.12f}", std::abs(verifyWeightedJaccardScore - incrWeightedJaccardScore));
                logging::critical("  Numerator diff: {}", verifyWeightedJaccNum - genomeState.currentWeightedJaccardNumerator);
                logging::critical("  Denominator diff: {}", verifyWeightedJaccDenom - genomeState.currentWeightedJaccardDenominator);
            }
            
            if (cosDiff) {
                logging::critical("COSINE MISMATCH:");
                logging::critical("  Verified:    {:.12f} (numerator={:.6f}, magnitude={:.6f})",
                                verifyCosineScore, verifyCosineNum, verifyGenomeMag);
                logging::critical("  Incremental: {:.12f} (numerator={:.6f}, magnitude={:.6f})",
                                incrCosineScore, genomeState.currentCosineNumerator, genomeState.currentGenomeMagnitude);
                logging::critical("  Difference:  {:.12f}", std::abs(verifyCosineScore - incrCosineScore));
                logging::critical("  Numerator diff: {:.6f}", verifyCosineNum - genomeState.currentCosineNumerator);
                logging::critical("  Magnitude diff: {:.6f}", verifyGenomeMag - genomeState.currentGenomeMagnitude);
            }
            
            // Print genome state details for investigation
            logging::critical("Current genome state:");
            logging::critical("  Position map size: {}", genomeState.positionMap.size());
            logging::critical("  Unique seeds: {}", genomeState.currentUniqueSeeds.size());
            logging::critical("  Seed counts map size: {}", genomeState.currentSeedCounts.size());
            
            // Print first few seeds in genome for debugging
            logging::critical("First 10 seeds in verified genome:");
            int seedCount = 0;
            for (const auto& [hash, count] : verifyGenomeSeedCounts) {
                if (seedCount++ >= 10) break;
                auto readIt = state.seedFreqInReads.find(hash);
                int64_t readCount = (readIt != state.seedFreqInReads.end()) ? readIt->second : 0;
                logging::critical("    hash={}, genomeCount={}, readCount={}", hash, count, readCount);
            }
            
            logging::critical("========================================");
            logging::critical("STOPPING at first mismatch for investigation");
            logging::critical("========================================");
            
            // Exit immediately to allow investigation
            std::exit(1);
        }
        
        // Use verified scores
        jaccardScore = verifyJaccardScore;
        weightedJaccardScore = verifyWeightedJaccardScore;
        cosineScore = verifyCosineScore;
        presenceJaccardScore = verifyPresenceScore;
    } else {
        // Normal path: use incremental scores (much faster)
        
        // Regular Jaccard: numerator and denominator both cached
        // Denominator = total_read_freq + genome_only_freq (cached incrementally)
        int64_t jaccardDenominator = state.totalReadSeedFrequency + genomeState.currentJaccardDenominator;
        if (jaccardDenominator > 0) {
            jaccardScore = static_cast<double>(genomeState.currentJaccardNumerator) / jaccardDenominator;
        }

        // Weighted Jaccard: numerator and denominator both cached
        int64_t weightedJaccardDenominator = genomeState.currentWeightedJaccardDenominator;
        if (weightedJaccardDenominator > 0) {
            weightedJaccardScore = static_cast<double>(genomeState.currentWeightedJaccardNumerator) / weightedJaccardDenominator;
        }
        
        // Debug logging for weighted Jaccard
        if (weightedJaccardScore > 1.0 || weightedJaccardScore < 0.0) {
            logging::warn("Node {}: INVALID weighted Jaccard = {:.6f} (numerator={}, denominator={})",
                         nodeId, weightedJaccardScore, genomeState.currentWeightedJaccardNumerator, 
                         weightedJaccardDenominator);
        }
        
        // Cosine similarity: numerator and magnitude both cached
        if (state.readMagnitude > 0.0 && genomeState.currentGenomeMagnitude > 0.0) {
            cosineScore = genomeState.currentCosineNumerator / (state.readMagnitude * genomeState.currentGenomeMagnitude);
        }
        
        // Presence/absence Jaccard: intersection and union both cached
        size_t intersectionCount = genomeState.currentPresenceIntersectionCount;
        size_t unionSize = state.seedFreqInReads.size() + genomeState.currentPresenceUnionCount;
        if (unionSize > 0) {
            presenceJaccardScore = static_cast<double>(intersectionCount) / unionSize;
        }
    }
    
    // Debug: Log first few nodes with detailed numeric diagnostics (AFTER computing scores)
    // Increase precision and add a sanity check for zero/NaN values so we can detect ordering/type bugs.
    static int debugCount = 0;
    if (debugCount < 10) {
        // Get current denominator for logging (whether from verify or incremental path)
        int64_t currentWeightedJaccDenom = VERIFY_INCREMENTAL_SCORES ? 
            0 : genomeState.currentWeightedJaccardDenominator; // For verify path, we already logged above
        
        // Log with higher precision
        logging::info("DEBUG Node {}: wJacc={:.9f}, cosine={:.9f}, jaccard={:.9f}{}",
                     nodeId,
                     weightedJaccardScore,
                     cosineScore,
                     jaccardScore,
                     VERIFY_INCREMENTAL_SCORES ? " [VERIFIED]" : "");

        debugCount++;
    }
    
    // Combined weighted score (scale = 1.0 means use weighted Jaccard as the combined metric)
    double combinedWeightedScore = weightedJaccardScore;
    
    // Get intersection count (same in both paths)
    size_t intersectionCount = genomeState.currentPresenceIntersectionCount;
    
    // Special tracking for true source node (read from environment variable)
    static std::string expectedNode;
    static bool checkedEnv = false;
    if (!checkedEnv) {
        const char* envNode = std::getenv("EXPECTED_NODE");
        if (envNode) {
            expectedNode = std::string(envNode);
            logging::info("Tracking expected node: {}", expectedNode);
        }
        checkedEnv = true;
    }
    
    if (!expectedNode.empty() && nodeId == expectedNode) {
        logging::info("*** FOUND EXPECTED SOURCE NODE {} ***", nodeId);
        logging::info("    Jaccard: {:.9f}", jaccardScore);
        logging::info("    Cosine: {:.9f}", cosineScore);
        logging::info("    Weighted Jaccard: {:.9f}", weightedJaccardScore);
        logging::info("    Raw matches: {}", genomeState.currentJaccardNumerator);
        logging::info("    Current best Jaccard: {:.9f} (node: {})", result.bestJaccardScore, result.bestJaccardNodeId);
        logging::info("    Current best Cosine: {:.9f} (node: {})", result.bestCosineScore, result.bestCosineNodeId);
        logging::info("    Current best Weighted: {:.9f} (node: {})", result.bestWeightedJaccardScore, result.bestWeightedJaccardNodeId);
    }
    
    // Update best scores in result
    result.updateJaccardScore(nodeId, jaccardScore);
    result.updateWeightedJaccardScore(nodeId, weightedJaccardScore);
    result.updateCosineScore(nodeId, cosineScore);
    result.updateJaccardPresenceScore(nodeId, presenceJaccardScore);
    result.updateRawSeedMatchScore(nodeId, genomeState.currentJaccardNumerator);  // Raw = sum of read frequencies for matched seeds
    result.updateHitsScore(nodeId, intersectionCount);  // Hits = count of unique seeds matched
    result.updateWeightedScore(nodeId, combinedWeightedScore, 1.0);  // Combined weighted score
    
    logging::debug("Node {}: Jaccard={:.6f}, WeightedJaccard={:.6f}, Cosine={:.6f}, PresenceJaccard={:.6f}, RawMatches={}, Hits={}", 
                  nodeId, jaccardScore, weightedJaccardScore, cosineScore, presenceJaccardScore, 
                  genomeState.currentJaccardNumerator, intersectionCount);
}

// Update metrics after backtracking  
void updateMetricsAfterBacktrack(placement::MgsrGenomeState& genomeState,
                                const placement::PlacementGlobalState& state) {
    // After backtracking, we may need to recalculate some metrics
    // For efficiency, we could track changes more precisely, but for correctness
    // a simple approach is sufficient since this happens after score calculation
}

namespace placement {

using ::panmanUtils::Node;
using ::panmanUtils::Tree;

// Forward declarations
class PlacementResult;


void dumpKmerDebugData(
    const std::vector<std::string>& readSequences,
    int k,
    const std::string& outputFilename);

// Add this new declaration
void dumpSyncmerDetails(const std::string& filename, 
                       const std::string& label,
                       const absl::flat_hash_map<size_t, int64_t>& seedMap,
                       const absl::flat_hash_map<size_t, std::string>& kmerMap);

// Function to load dictionary entries from index
void loadGlobalDictionary(
    Index::Reader& indexReader,
    PlacementGlobalState& state,
    int k_param, 
    int s_param  
);

// Shared progress state for UI updates
std::shared_ptr<PlacementProgressState> progress_state;

// Maximum number of nodes for which to dump recomputation range details
const size_t MAX_NODES_TO_DUMP = 5;

// Helper function to read a packed MGSR index file
std::unique_ptr<::capnp::MessageReader> loadMgsrIndexFromFile(const std::string& indexPath) {
    logging::debug("Loading MGSR index from file: {}", indexPath);
    
    // Normalize path to ensure consistent resolution
    boost::filesystem::path normalizedPath = boost::filesystem::absolute(indexPath);
    std::string resolvedPath = normalizedPath.string();
    logging::debug("DEBUG-PLACE: Using normalized path: {}", resolvedPath);
    
    // First check if the file exists
    if (!boost::filesystem::exists(resolvedPath)) {
        throw std::runtime_error("Index file not found: " + resolvedPath);
    }
    
    // Get the file size for validation
    uintmax_t fileSize = boost::filesystem::file_size(resolvedPath);
    if (fileSize == 0) {
        throw std::runtime_error("Index file is empty: " + resolvedPath);
    }
    logging::debug("Index file exists and has size {} bytes", fileSize);
    
    // Open the file
    int fd = open(resolvedPath.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("Failed to open index file: " + resolvedPath);
    }
    
    try {
        // Configure reader options
        ::capnp::ReaderOptions opts;
        opts.traversalLimitInWords = std::numeric_limits<uint64_t>::max();
        opts.nestingLimit = 1024;
        
        // Use PackedFdMessageReader directly
        auto reader = std::make_unique<::capnp::PackedFdMessageReader>(fd, opts);
        
        // Validate the MGSR index immediately to ensure it's usable
        auto mgsrIndex = reader->getRoot<MGSRIndex>();
        uint32_t k = mgsrIndex.getK();
        uint32_t s = mgsrIndex.getS();
        size_t seedCount = mgsrIndex.getSeedInfo().size();
        bool useRawSeeds = mgsrIndex.getUseRawSeeds();
        
        if (k <= 0 || seedCount <= 0) {
            close(fd); // Close the fd since we're not returning the reader
            throw std::runtime_error(
                "Invalid MGSR index data: k=" + std::to_string(k) + 
                ", seeds=" + std::to_string(seedCount) + 
                ". The index appears to be corrupted or incomplete.");
        }
        
        logging::debug("Successfully loaded MGSR index with k={}, s={}, seeds={}, useRawSeeds={}",
                    k, s, seedCount, useRawSeeds);
        
        // Return the MessageReader
        return reader;
    } catch (const ::kj::Exception& e) {
        // Clean up fd on exception
        close(fd);
        throw std::runtime_error("Cap'n Proto error: " + std::string(e.getDescription().cStr()));
    } catch (const std::exception& e) {
        // Clean up fd on exception
        close(fd);
        throw std::runtime_error("Error reading MGSR index: " + std::string(e.what()));
    }
}


// Helper functions for processNodeMutations

// Process seed deletion (quaternary value 1)
void processSeedDeletion(
    ::panmanUtils::Node* node,
    int64_t pos,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int k) {
    
    // Get identifier, using a fallback only if node is null
    std::string nodeId = node ? node->identifier : (result.bestWeightedNodeId.empty() ? "Unknown" : result.bestWeightedNodeId);
    
    // Check if this seed exists in our map
    if (node && result.nodeSeedMap.find(nodeId) != result.nodeSeedMap.end() && 
        result.nodeSeedMap[nodeId].find(pos) != result.nodeSeedMap[nodeId].end()) {
        
        // Get the seed from our map
        const seeding::seed_t& seed = result.nodeSeedMap[nodeId][pos];
        const size_t seedHash = seed.hash;
        
        // Remove it from our map
        result.nodeSeedMap[nodeId].erase(pos);
        
        // Update scores if this seed was in reads
        auto readIt = state.seedFreqInReads.find(seedHash);
        if (readIt != state.seedFreqInReads.end()) {
            const int64_t readCount = readIt->second;
            result.hitsInThisGenome -= readCount;
            result.currentJaccardNumerator -= readCount;
            
            // Handle weighted Jaccard numerator update 
            auto genomeSeedIt = result.currentAllGenomeSeedCounts.find(seedHash);
            int64_t oldGenomeCount = 0;
            int64_t newGenomeCount = 0;
            
            if (genomeSeedIt != result.currentAllGenomeSeedCounts.end()) {
                oldGenomeCount = genomeSeedIt->second;
                newGenomeCount = oldGenomeCount - 1;
                
                // Update weighted Jaccard numerator: subtract old contribution, add new contribution
                int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                int64_t newMinCount = (newGenomeCount > 0) ? std::min(readCount, newGenomeCount) : 0;
                result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                
                if (newGenomeCount <= 0) { 
                    result.currentAllGenomeSeedCounts.erase(genomeSeedIt);
                } else {
                    genomeSeedIt->second = newGenomeCount;
                }
            }
            
            // Update cosine similarity components (numerator: readCount * genomeCount)
            double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount); // Before deletion
            double newCosineTerm = static_cast<double>(readCount * newGenomeCount); // After deletion
            result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
            
            uniqueSeedHashes.erase(seedHash);
            
            logging::info("Deleted seed from scores: pos={}, node={}, hash={}, readCount={}", 
                        pos, nodeId, seedHash, readCount);
        }
    } else {
        logging::debug("No existing seed found at position {} to delete for node {}", pos, nodeId);
    }
}

// Process existing seed for modification (part of quaternary value 3)
void processExistingSeedRemoval(
    ::panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    int k) {
    
    if (!node) {
        logging::warn("Cannot process seed removal: node is null");
        return;
    }
    
    std::string nodeId = node->identifier;
    std::optional<seeding::seed_t> existingSeedOpt;
    
    // Check if this seed exists in our map
    if (result.nodeSeedMap.find(nodeId) != result.nodeSeedMap.end() && 
        result.nodeSeedMap[nodeId].find(pos) != result.nodeSeedMap[nodeId].end()) {
        // Found the seed in our map
        existingSeedOpt = result.nodeSeedMap[nodeId][pos];
        
        // Process the found seed
        if (existingSeedOpt) {
            const seeding::seed_t& existingSeed = existingSeedOpt.value();
            const size_t seedHash = existingSeed.hash;
            
            // Remove it from our map
            result.nodeSeedMap[nodeId].erase(pos);
            
            // Update scores if this seed was in reads
            auto readIt = state.seedFreqInReads.find(seedHash);
            if (readIt != state.seedFreqInReads.end()) {
                int64_t readCount = readIt->second;
                result.hitsInThisGenome -= readCount;
                result.currentJaccardNumerator -= readCount;
                
                // Handle weighted Jaccard numerator update
                auto genomeSeedIt = result.currentAllGenomeSeedCounts.find(seedHash);
                int64_t oldGenomeCount = 0;
                int64_t newGenomeCount = 0;
                
                if (genomeSeedIt != result.currentAllGenomeSeedCounts.end()) {
                    oldGenomeCount = genomeSeedIt->second;
                    newGenomeCount = oldGenomeCount - 1;
                    
                    // Update weighted Jaccard numerator: subtract old contribution, add new contribution
                    int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                    int64_t newMinCount = (newGenomeCount > 0) ? std::min(readCount, newGenomeCount) : 0;
                    result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                    
                    if (newGenomeCount <= 0) {
                        result.currentAllGenomeSeedCounts.erase(genomeSeedIt);
                    } else {
                        genomeSeedIt->second = newGenomeCount;
                    }
                }
                
                // Update cosine similarity components (numerator: readCount * genomeCount)
                double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount); // Before removal
                double newCosineTerm = static_cast<double>(readCount * newGenomeCount); // After removal
                result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
                
                uniqueSeedHashes.erase(seedHash);
                
                logging::info("Removed seed from scores: pos={}, node={}, hash={}, readCount={}", 
                            pos, nodeId, seedHash, readCount);
            }
        }
    } else {
        logging::debug("No existing seed found at position {} to remove", pos);
    }
}

// Add new seed and update scores
void addSeedAndUpdateScores(
    ::panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    const seeding::seed_t& newSeed) {
    

    if (!node) {
        logging::warn("Cannot add seed: node is null");
        return;
    }
    
    // Add hash to uniqueSeedHashes
    uniqueSeedHashes.insert(newSeed.hash);
    
    // Add seed to our direct map for this node
    if (result.nodeSeedMap.find(node->identifier) == result.nodeSeedMap.end()) {
        result.nodeSeedMap[node->identifier] = {};
    }
    result.nodeSeedMap[node->identifier][pos] = newSeed;
    
    // Update scores if seed exists in reads
    auto readIt = state.seedFreqInReads.find(newSeed.hash);
    bool inReads = readIt != state.seedFreqInReads.end();


    
    if (inReads) {
        int64_t readCount = readIt->second;
        result.hitsInThisGenome += readCount;
        result.currentJaccardNumerator += readCount;
        
        // Safely increment count in currentAllGenomeSeedCounts
        auto [it, inserted] = result.currentAllGenomeSeedCounts.try_emplace(newSeed.hash, 0);
        int64_t oldGenomeCount = it->second;  // Should be 0 for new seeds
        it->second++; // Increment the count
        int64_t newGenomeCount = it->second;  // Should be 1 for new seeds
        
        // Update weighted Jaccard numerator: min(readCount, genomeCount)
        int64_t oldMinCount = std::min(readCount, oldGenomeCount);  // Should be 0
        int64_t newMinCount = std::min(readCount, newGenomeCount);  // Should be min(readCount, 1)
        result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
        
        // Update cosine similarity components (numerator: readCount * genomeCount)
        double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount); // Before addition
        double newCosineTerm = static_cast<double>(readCount * newGenomeCount); // After addition
        result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
        
        // Log successful match (only if seed is in reads)
        logging::info("Added seed at pos {} for node {} with hash {}, read count: {}", 
                     pos, node->identifier, newSeed.hash, readCount);
    }
}

// Process new seed addition (quaternary value 2 or part of 3)
// --- HELPER: Move the getNodeIndex function up here ---

/**
 * @brief Helper function to find the node's index in the Cap'n Proto arrays
 * 
 * @param nodeId The node identifier
 * @param state The global placement state containing node path info
 * @param nodeIndex Output parameter to store the found index
 * @return True if node was found, false otherwise
 */
bool getNodeIndex(const std::string& nodeId, 
                  const PlacementGlobalState& state, 
                  uint64_t& nodeIndex) {
    // MGSR format doesn't use nodePathInfo lookup - DFS traversal handles indexing
    logging::debug("MGSR mode: getNodeIndex called for node '{}' but MGSR uses DFS traversal", nodeId);
    return false;
}


/**
 * @brief Update Jaccard similarity score for a node and track the best score
 * 
 * @param nodeId The node ID being evaluated
 * @param jaccardScore The Jaccard similarity score
 */
void PlacementResult::updateJaccardScore(const std::string& nodeId, double jaccardScore) {
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (jaccardScore > bestJaccardScore + TIED_THRESHOLD) {
        // Found a new best Jaccard score
        bestJaccardScore = jaccardScore;
        bestJaccardNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedJaccardNodeIds.clear();
        tiedJaccardNodeIds.push_back(nodeId);
        
        logging::debug("New best Jaccard node: {} with score {:.6f}", 
                     nodeId, jaccardScore);
                     
    } else if (std::abs(jaccardScore - bestJaccardScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedJaccardNodeIds.push_back(nodeId);
        logging::debug("Tied Jaccard node: {} with score {:.6f}", 
                     nodeId, jaccardScore);
    }
}

/**
 * Update the weighted Jaccard similarity score tracking
 * @param nodeId The node ID being evaluated  
 * @param weightedJaccardScore The weighted Jaccard similarity score
 */
void PlacementResult::updateWeightedJaccardScore(const std::string& nodeId, double weightedJaccardScore) {
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (weightedJaccardScore > bestWeightedJaccardScore + TIED_THRESHOLD) {
        // Found a new best weighted Jaccard score
        bestWeightedJaccardScore = weightedJaccardScore;
        bestWeightedJaccardNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedWeightedJaccardNodeIds.clear();
        tiedWeightedJaccardNodeIds.push_back(nodeId);
        
        logging::debug("New best weighted Jaccard node: {} with score {:.6f}", 
                     nodeId, weightedJaccardScore);
                     
    } else if (std::abs(weightedJaccardScore - bestWeightedJaccardScore) <= TIED_THRESHOLD) {
        // Add to tied nodes if not already present
        if (std::find(tiedWeightedJaccardNodeIds.begin(), tiedWeightedJaccardNodeIds.end(), nodeId) 
            == tiedWeightedJaccardNodeIds.end()) {
            tiedWeightedJaccardNodeIds.push_back(nodeId);
            logging::debug("Tied weighted Jaccard node: {} with score {:.6f}", 
                         nodeId, weightedJaccardScore);
        }
    }
}

/**
 * @brief Update Raw Seed Match score for a node and track the best score
 * 
 * @param nodeId The node ID being evaluated
 * @param score The raw seed match score (sum of read frequencies for matched seeds)
 */
void PlacementResult::updateRawSeedMatchScore(const std::string& nodeId, int64_t score) {
    // Check if score is better than current best
    if (score > bestRawSeedMatchScore) {
        // Found a new best raw seed match score
        bestRawSeedMatchScore = score;
        bestRawSeedMatchNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedRawSeedMatchNodeIds.clear();
        tiedRawSeedMatchNodeIds.push_back(nodeId);
        
        logging::debug("New best Raw Seed Match node: {} with score {}", 
                     nodeId, score);
                     
    } else if (score == bestRawSeedMatchScore) {
        // Add to tied nodes
        tiedRawSeedMatchNodeIds.push_back(nodeId);
        logging::debug("Tied Raw Seed Match node: {} with score {}", 
                     nodeId, score);
    }
}

/**
 * @brief Update hits score for a node and track the best score
 * 
 * @param nodeId The node ID being evaluated
 * @param hits The number of hits for this node
 */
void PlacementResult::updateHitsScore(const std::string& nodeId, int64_t hits) {
    // Check if hits is better than current best
    if (hits > maxHitsInAnyGenome) {
        // Found a new best hits score
        maxHitsInAnyGenome = hits;
        maxHitsNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedMaxHitsNodeIds.clear();
        tiedMaxHitsNodeIds.push_back(nodeId);
        
        logging::debug("New best hits node: {} with hits {}", 
                     nodeId, hits);
                     
    } else if (hits == maxHitsInAnyGenome) {
        // Add to tied nodes
        tiedMaxHitsNodeIds.push_back(nodeId);
        logging::debug("Tied hits node: {} with hits {}", 
                     nodeId, hits);
    }
}

/**
 * @brief Update Jaccard (Presence/Absence) score for a node and track the best score
 * 
 * @param nodeId The node ID being evaluated
 * @param score The Jaccard score based on presence/absence of seeds
 */
void PlacementResult::updateJaccardPresenceScore(const std::string& nodeId, double score) {
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (score > bestJaccardPresenceScore + TIED_THRESHOLD) {
        // Found a new best Jaccard (Presence/Absence) score
        bestJaccardPresenceScore = score;
        bestJaccardPresenceNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedJaccardPresenceNodeIds.clear();
        tiedJaccardPresenceNodeIds.push_back(nodeId);
        
        logging::debug("New best Jaccard (Presence) node: {} with score {:.6f}", 
                     nodeId, score);
                     
    } else if (std::abs(score - bestJaccardPresenceScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedJaccardPresenceNodeIds.push_back(nodeId);
        logging::debug("Tied Jaccard (Presence) node: {} with score {:.6f}", 
                     nodeId, score);
    }
}

/**
 * @brief Update cosine similarity score for a node and track the best score
 * 
 * @param nodeId The node ID being evaluated
 * @param cosineScore The cosine similarity score
 */
void PlacementResult::updateCosineScore(const std::string& nodeId, double cosineScore) {
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (cosineScore > bestCosineScore + TIED_THRESHOLD) {
        // Found a new best cosine score
        bestCosineScore = cosineScore;
        bestCosineNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedCosineNodeIds.clear();
        tiedCosineNodeIds.push_back(nodeId);
        
        logging::debug("New best cosine node: {} with score {:.6f}", 
                     nodeId, cosineScore);
                     
    } else if (std::abs(cosineScore - bestCosineScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedCosineNodeIds.push_back(nodeId);
        logging::debug("Tied cosine node: {} with score {:.6f}", 
                     nodeId, cosineScore);
    }
}

/**
 * @brief Update weighted score (combined Jaccard and cosine) for a node
 * 
 * @param nodeId The node ID being evaluated
 * @param weightedScore The weighted similarity score
 * @param scoreScale The weight for Jaccard in the combined score (1-scoreScale is used for cosine)
 */
void PlacementResult::updateWeightedScore(const std::string& nodeId, double weightedScore, double scoreScale) {
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (weightedScore > bestWeightedScore + TIED_THRESHOLD) {
        // Found a new best weighted score
        bestWeightedScore = weightedScore;
        bestWeightedNodeId = nodeId;
        
        // Reset tied nodes and add this one
        tiedWeightedNodeIds.clear();
        tiedWeightedNodeIds.push_back(nodeId);
        
        logging::debug("New best weighted node: {} with score {:.6f} (scale={:.2f})", 
                     nodeId, weightedScore, scoreScale);
                     
    } else if (std::abs(weightedScore - bestWeightedScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedWeightedNodeIds.push_back(nodeId);
        logging::debug("Tied weighted node: {} with score {:.6f}", 
                     nodeId, weightedScore);
    }
}



// LiteTree-specific scoring function (avoids StateManager dependency)
std::tuple<double, double, size_t> calculateAndUpdateScoresLite(
    const std::string& nodeId,
    PlacementGlobalState& state,
    PlacementResult& result,
    const TraversalParams& params,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    uint32_t nodeChangeIndex) {
    
    // Process current node's changes if available
    if (nodeChangeIndex < state.perNodeChanges.size()) {
        auto nodeChanges = state.perNodeChanges[nodeChangeIndex];
        
        // Process seed deltas (combined insertions and deletions)
        auto seedDeltas = nodeChanges.getSeedDeltas();
        for (auto delta : seedDeltas) {
            uint32_t seedIdx = delta.getSeedIndex();
            bool isDeleted = delta.getIsDeleted();
            
            if (seedIdx < state.seedInfo.size()) {
                auto seedInfo = state.seedInfo[seedIdx];
                size_t seedHash = seedInfo.getHash();
                
                if (isDeleted) {
                    // Process deletion
                    // Remove from unique seed hashes
                    uniqueSeedHashes.erase(seedHash);
                    
                    // Update ALL genome seed counts and magnitude
                    auto allGenomeIt = result.currentAllGenomeSeedCounts.find(seedHash);
                    if (allGenomeIt != result.currentAllGenomeSeedCounts.end()) {
                        int64_t oldGenomeCount = allGenomeIt->second;
                        int64_t newGenomeCount = oldGenomeCount - 1;
                        
                        // Update genome magnitude squared
                        result.currentGenomeMagnitudeSquared += (static_cast<double>(newGenomeCount * newGenomeCount) - 
                                                                  static_cast<double>(oldGenomeCount * oldGenomeCount));
                        
                        if (newGenomeCount <= 0) {
                            result.currentAllGenomeSeedCounts.erase(allGenomeIt);
                        } else {
                            allGenomeIt->second = newGenomeCount;
                        }
                        
                        // Update metrics if this seed exists in reads
                        auto readIt = state.seedFreqInReads.find(seedHash);
                        if (readIt != state.seedFreqInReads.end()) {
                            int64_t readCount = readIt->second;
                            
                            // Update weighted Jaccard numerator
                            int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                            int64_t newMinCount = std::min(readCount, newGenomeCount);
                            result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                            
                            // Update regular Jaccard numerator (if seed becomes absent, decrease by readCount)
                            if (oldGenomeCount > 0 && newGenomeCount == 0) {
                                result.currentJaccardNumerator -= readCount;
                            }
                            
                            // Update raw match score
                            result.hitsInThisGenome += (readCount * newGenomeCount) - (readCount * oldGenomeCount);
                            
                            // Update cosine similarity numerator
                            double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount);
                            double newCosineTerm = static_cast<double>(readCount * newGenomeCount);
                            result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
                        }
                    }
                } else {
                    // Process insertion/substitution
                    // Add to unique seed hashes
                    uniqueSeedHashes.insert(seedHash);
                    
                    // Update ALL genome seed counts and magnitude
                    auto [allGenomeIt, inserted] = result.currentAllGenomeSeedCounts.try_emplace(seedHash, 0);
                    int64_t oldGenomeCount = allGenomeIt->second;
                    int64_t newGenomeCount = oldGenomeCount + 1;
                    allGenomeIt->second = newGenomeCount;
                    
                    // Update genome magnitude squared
                    result.currentGenomeMagnitudeSquared += (static_cast<double>(newGenomeCount * newGenomeCount) - 
                                                              static_cast<double>(oldGenomeCount * oldGenomeCount));
                    
                    // Update metrics if this seed exists in reads
                    auto readIt = state.seedFreqInReads.find(seedHash);
                    if (readIt != state.seedFreqInReads.end()) {
                        int64_t readCount = readIt->second;
                        
                        // Update weighted Jaccard numerator
                        int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                        int64_t newMinCount = std::min(readCount, newGenomeCount);
                        result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                        
                        // Update regular Jaccard numerator (if seed becomes present, increase by readCount)
                        if (oldGenomeCount == 0 && newGenomeCount > 0) {
                            result.currentJaccardNumerator += readCount;
                        }
                        
                        // Update raw match score
                        result.hitsInThisGenome += (readCount * newGenomeCount) - (readCount * oldGenomeCount);
                        
                        // Update cosine similarity numerator
                        double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount);
                        double newCosineTerm = static_cast<double>(readCount * newGenomeCount);
                        result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
                    }
                }
            }
        }
    }
    
    // Calculate final scores from current metric components
    double jaccardScore = 0.0;
    double cosineScore = 0.0;
    size_t intersectionCount = 0;
    
    // Count intersection for presence/absence Jaccard
    for (const size_t seedHash : uniqueSeedHashes) {
        if (state.seedFreqInReads.find(seedHash) != state.seedFreqInReads.end()) {
            intersectionCount++;
        }
    }
    
    // Calculate regular Jaccard similarity (presence/absence with frequency weighting)
    // Denominator for regular Jaccard: sum of read frequencies for seeds that are present in genome
    int64_t regularJaccardDenominator = 0;
    
    // Add contributions from read seeds that are present in genome (intersection + read-only seeds)
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        auto genomeIt = result.currentAllGenomeSeedCounts.find(seedHash);
        if (genomeIt != result.currentAllGenomeSeedCounts.end() && genomeIt->second > 0) {
            // Seed is in intersection - already counted in numerator
            regularJaccardDenominator += readCount;
        } else {
            // Seed is read-only - add to denominator but not in numerator  
            regularJaccardDenominator += readCount;
        }
    }
    
    if (regularJaccardDenominator > 0) {
        jaccardScore = static_cast<double>(result.currentJaccardNumerator) / static_cast<double>(regularJaccardDenominator);
    }
    
    // Calculate weighted Jaccard similarity 
    // Denominator: sum of max(read_freq, genome_freq) for each seed in union
    int64_t weightedJaccardDenominator = 0;
    
    // Add contributions from all read seeds
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        auto genomeIt = result.currentAllGenomeSeedCounts.find(seedHash);
        int64_t genomeCount = (genomeIt != result.currentAllGenomeSeedCounts.end()) ? genomeIt->second : 0;
        weightedJaccardDenominator += std::max(readCount, genomeCount);
    }
    
    // Add contributions from genome seeds not in reads
    for (const auto& [seedHash, genomeCount] : result.currentAllGenomeSeedCounts) {
        if (state.seedFreqInReads.find(seedHash) == state.seedFreqInReads.end()) {
            weightedJaccardDenominator += genomeCount;
        }
    }
    
    double weightedJaccardScore = 0.0;
    if (weightedJaccardDenominator > 0) {
        weightedJaccardScore = static_cast<double>(result.currentWeightedJaccardNumerator) / static_cast<double>(weightedJaccardDenominator);
    }
    
    // Calculate cosine similarity using tracked genome magnitude
    double genomeMagnitude = std::sqrt(result.currentGenomeMagnitudeSquared);
    
    // Use precomputed read magnitude and tracked genome magnitude
    if (state.readMagnitude > 0.0 && genomeMagnitude > 0.0) {
        cosineScore = result.currentCosineNumerator / (state.readMagnitude * genomeMagnitude);
    }
    
    // Calculate presence/absence Jaccard
    double jaccardPresenceScore = 0.0;
    size_t totalReadSeeds = state.seedFreqInReads.size();
    size_t genomeUniqueSeeds = uniqueSeedHashes.size();
    size_t unionSize = totalReadSeeds + genomeUniqueSeeds - intersectionCount;
    
    if (unionSize > 0) {
        jaccardPresenceScore = static_cast<double>(intersectionCount) / static_cast<double>(unionSize);
    }
    
    // Update PlacementResult with computed scores
    result.updateJaccardScore(nodeId, jaccardScore);
    result.updateCosineScore(nodeId, cosineScore);
    result.updateRawSeedMatchScore(nodeId, result.hitsInThisGenome);
    result.updateJaccardPresenceScore(nodeId, jaccardPresenceScore);
    result.updateHitsScore(nodeId, result.hitsInThisGenome);
    result.updateWeightedScore(nodeId, weightedJaccardScore, 1.0);  // Use weighted Jaccard as weighted score
    
    logging::debug("Node {}: Jaccard={:.6f}, WeightedJaccard={:.6f}, Cosine={:.6f}, JaccardPresence={:.6f}, RawMatch={}, "
                  "intersectionCount={}", 
                  nodeId, jaccardScore, weightedJaccardScore, cosineScore, jaccardPresenceScore, 
                  result.hitsInThisGenome, intersectionCount);
    
    return {jaccardScore, cosineScore, intersectionCount};
}

void placeLiteHelper(panmapUtils::LiteNode* node, 
                    placement::MgsrGenomeState& genomeState,
                    uint32_t& dfsIndex,
                    const placement::PlacementGlobalState& state,
                    placement::PlacementResult& result,
                    const placement::TraversalParams& params) {

    if(dfsIndex % 1000 == 0) {
        std::cout << "\rplace dfsIndex: " << dfsIndex << std::flush;
    }
    
    if (!node) return;
    
    const std::string& nodeId = node->identifier;
    
    // TODO
    }
    
void placeLite(placement::PlacementResult &result, 
                       panmapUtils::LiteTree *liteTree,
                       ::MGSRIndex::Reader &mgsrIndex, 
                       const std::string &reads1,
                       const std::string &reads2,
                       std::vector<std::vector<seeding::seed_t>>& readSeeds,
                       std::vector<std::string>& readSequences,
                       std::vector<std::string>& readNames,
                       std::vector<std::string>& readQuals,
                       std::string &outputPath,
                       const std::string &indexPath,
                       const std::string &debug_node_id_param,
                       bool verify_scores,
                       panmanUtils::Tree *fullTree) {
    
    logging::info("Starting MGSR-style recursive placement");
    
    // Store verify_scores flag in a place accessible to calculateMgsrStyleScores
    if (verify_scores) {
        logging::info("VERIFICATION MODE: Will recompute all scores from scratch at each node");
        if (fullTree == nullptr) {
            logging::err("VERIFICATION MODE requires fullTree pointer but got nullptr!");
            throw std::runtime_error("Verification mode requires full Tree to be loaded");
        }
    }
    
    auto time_init_start = std::chrono::high_resolution_clock::now();
    
    // Initialize placement global state
    PlacementGlobalState state;
    state.seedInfo = mgsrIndex.getSeedInfo();
    state.perNodeChanges = mgsrIndex.getPerNodeChanges();
    state.kmerSize = mgsrIndex.getK();
    state.fullTree = fullTree;  // Store full tree pointer for verification mode
    
    // Set traversal parameters
    TraversalParams params;
    params.k = mgsrIndex.getK();
    params.s = mgsrIndex.getS();
    params.t = mgsrIndex.getT();
    params.open = mgsrIndex.getOpen();
    params.debug_node_id = debug_node_id_param;
    params.verify_scores = verify_scores; // Pass flag through params
    
    auto time_init_end = std::chrono::high_resolution_clock::now();
    auto duration_init = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_init_end - time_init_start);
    logging::info("[TIME] Placement initialization: {}ms", duration_init.count());
    
    // Process reads and extract seeds using MGSR method
    auto time_read_start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> allReadSequences;
    
    if (!reads1.empty()) {
        mgsr::extractReadSequences(reads1, reads2, allReadSequences);
    } else if (!readSequences.empty()) {
        allReadSequences = readSequences;
    }
    
    mgsr::ThreadsManager threadsManager(liteTree, 1);  
    threadsManager.initializeMGSRIndex(mgsrIndex);
    
    uint32_t l = mgsrIndex.getL();
    logging::info("Index parameters: k={}, s={}, t={}, l={}", 
                  params.k, params.s, params.t, l);
    
    // If l=0, use raw syncmers instead of k-minimizers (parallelized with deduplication)
    if (l == 0) {
        logging::info("l=0: Using raw syncmers instead of k-minimizers");
        
        // Step 1: Sort reads to group duplicates together
        std::vector<std::pair<std::string, size_t>> sortedReads;
        sortedReads.reserve(allReadSequences.size());
        for (size_t i = 0; i < allReadSequences.size(); i++) {
            sortedReads.emplace_back(allReadSequences[i], i);
        }
        tbb::parallel_sort(sortedReads.begin(), sortedReads.end(), 
            [](const auto& a, const auto& b) { return a.first < b.first; });
        
        // Step 2: Build dupReadsIndex - map unique sequences to their duplicate indices
        std::vector<std::pair<std::string, std::vector<size_t>>> dupReadsIndex;
        for (size_t i = 0; i < sortedReads.size(); ) {
            std::vector<size_t> duplicates;
            const std::string& uniqueSeq = sortedReads[i].first;
            while (i < sortedReads.size() && sortedReads[i].first == uniqueSeq) {
                duplicates.push_back(sortedReads[i].second);
                i++;
            }
            dupReadsIndex.emplace_back(uniqueSeq, std::move(duplicates));
        }
        
        logging::info("Deduplication: {} reads → {} unique sequences", 
                      allReadSequences.size(), dupReadsIndex.size());
        
        // Step 3: Extract syncmers from unique sequences in parallel
        size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
        std::vector<absl::flat_hash_map<size_t, int64_t>> threadLocalMaps(num_cpus);
        
        tbb::parallel_for(tbb::blocked_range<size_t>(0, dupReadsIndex.size()), 
            [&](const tbb::blocked_range<size_t>& range) {
                size_t threadId = tbb::this_task_arena::current_thread_index();
                auto& localMap = threadLocalMaps[threadId];
                
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const auto& [seq, duplicates] = dupReadsIndex[i];
                    int64_t multiplicity = duplicates.size();
                    
                    const auto& syncmers = seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false);
                    for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : syncmers) {
                        if (!isSyncmer) continue;
                        localMap[kmerHash] += multiplicity;
                    }
                }
            });
        
        // Step 4: Merge thread-local maps into global map
        for (const auto& localMap : threadLocalMaps) {
            for (const auto& [hash, count] : localMap) {
                state.seedFreqInReads[hash] += count;
            }
        }
        
        logging::info("Extracted {} unique syncmers from {} reads ({} unique patterns)", 
                      state.seedFreqInReads.size(), allReadSequences.size(), dupReadsIndex.size());
    } else {
        // Use MGSR's k-minimizer processing
        threadsManager.initializeQueryData(allReadSequences, false);
        
        // Debug: Check what we got from MGSR processing
        logging::info("MGSR processed {} unique read patterns from {} total reads", 
                      threadsManager.reads.size(), allReadSequences.size());
        
        size_t readsWithSeeds = 0;
        size_t totalKminmers = 0;
        for (const auto& read : threadsManager.reads) {
            if (!read.uniqueSeedmers.empty()) {
                readsWithSeeds++;
                for (const auto& [hash, positions] : read.uniqueSeedmers) {
                    totalKminmers += positions.size();
                }
            }
        }
        logging::info("Reads with k-minimizers: {}/{}, total k-minimizers: {}", 
                      readsWithSeeds, threadsManager.reads.size(), totalKminmers);
        
        for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
            const auto& read = threadsManager.reads[i];
            for (const auto& [kminmerHash, positions] : read.uniqueSeedmers) {
                // Count each unique k-minimizer occurrence
                state.seedFreqInReads[kminmerHash] += positions.size() * threadsManager.readSeedmersDuplicatesIndex[i].size();
            }
        }
        
        logging::info("Extracted {} unique k-minimizers from {} reads ({} unique read patterns)", 
                      state.seedFreqInReads.size(), allReadSequences.size(), threadsManager.reads.size());
    }
    
    auto time_read_end = std::chrono::high_resolution_clock::now();
    auto duration_read = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_read_end - time_read_start);
    logging::info("[TIME] Read processing & seed extraction: {}ms", duration_read.count());
    
    // Precompute read magnitude for cosine similarity and total read seed frequency
    auto time_precomp_start = std::chrono::high_resolution_clock::now();
    state.totalReadSeedFrequency = 0;
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        state.readMagnitude += static_cast<double>(readCount * readCount);
        state.totalReadSeedFrequency += readCount;
    }
    state.readMagnitude = std::sqrt(state.readMagnitude);
    logging::info("Precomputed read magnitude: {:.6f}, total read seed frequency: {}", 
                  state.readMagnitude, state.totalReadSeedFrequency);
    
    auto time_precomp_end = std::chrono::high_resolution_clock::now();
    auto duration_precomp = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_precomp_end - time_precomp_start);
    if (duration_precomp.count() > 0) {
        logging::info("[TIME] Magnitude precomputation: {}ms", duration_precomp.count());
    }
    
    // Set root pointer in state
    state.root = liteTree->root;
    
    // Initialize MGSR-style global genome state (single instance, modified in-place)
    // Start with empty state - root will be initialized by placeLiteHelper
    placement::MgsrGenomeState globalGenomeState;
    
    // Initialize denominators for empty genome state:
    // - Weighted Jaccard denominator = sum of all read frequencies (max when genome is empty)
    // - Jaccard denominator starts at 0 (no genome-only seeds yet)
    // - Presence union count tracks GENOME-ONLY seeds (starts at 0, read seeds added via formula)
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        globalGenomeState.currentWeightedJaccardDenominator += readCount;
    }
    globalGenomeState.currentPresenceUnionCount = 0;  // Genome-only seeds (none initially)
    globalGenomeState.currentJaccardDenominator = 0;  // No genome-only seeds initially
    
    logging::debug("Initialized: weighted Jaccard denom={}, presence union genome-only={}, Jaccard denom={}", 
                  globalGenomeState.currentWeightedJaccardDenominator,
                  globalGenomeState.currentPresenceUnionCount,
                  globalGenomeState.currentJaccardDenominator);
    
    uint32_t currentDfsIndex = 0;
    
    logging::info("Starting MGSR-style recursive placement traversal with {} tree nodes", 
                  state.perNodeChanges.size());
    
    // Start recursive traversal from root (dfsIndex=0)
    // placeLiteHelper will apply root's changes and calculate its scores
    auto time_traversal_start = std::chrono::high_resolution_clock::now();
    if (liteTree && liteTree->root) {
        placeLiteHelper(liteTree->root, globalGenomeState, currentDfsIndex, state, result, params);
    } else {
        logging::err("LiteTree or root is null, cannot perform placement");
        return;
    }
    
    auto time_traversal_end = std::chrono::high_resolution_clock::now();
    auto duration_traversal = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_traversal_end - time_traversal_start);
    logging::info("[TIME] Tree traversal & scoring: {}ms", duration_traversal.count());
    
    // Performance metrics
    result.totalReadsProcessed = allReadSequences.size();
    
    // Write placements.tsv file with all metrics
    // outputPath is actually the prefix directory, not a filename
    std::string placementsFilePath = outputPath;
    std::ofstream placementsFile(placementsFilePath);
    if (placementsFile.is_open()) {
        // Write header
        placementsFile << "metric\tscore\thits\tnodes\n";
        
        // Write raw matches (sum of read frequencies for matched seeds)
        placementsFile << "raw\t" << result.bestRawSeedMatchScore 
                      << "\t" << result.maxHitsInAnyGenome
                      << "\t" << result.bestRawSeedMatchNodeId << "\n";
        
        // Write jaccard (regular Jaccard with frequency weighting)
        placementsFile << "jaccard\t" << std::fixed << std::setprecision(6) << result.bestJaccardScore
                      << "\t\t" << result.bestJaccardNodeId << "\n";
        
        // Write cosine similarity
        placementsFile << "cosine\t" << std::fixed << std::setprecision(6) << result.bestCosineScore
                      << "\t\t" << result.bestCosineNodeId << "\n";
        
        // Write weighted Jaccard
        placementsFile << "weighted_jaccard\t" << std::fixed << std::setprecision(6) << result.bestWeightedJaccardScore
                      << "\t\t" << result.bestWeightedJaccardNodeId << "\n";
        
        placementsFile.close();
        logging::info("Wrote placement results to {}", placementsFilePath);
    } else {
        logging::err("Failed to open placements file: {}", placementsFilePath);
    }
    
    logging::info("MGSR-style placement completed. Best Jaccard score: {} (node: {})", 
                 result.bestJaccardScore, result.bestJaccardNodeId);
    logging::info("Best Cosine score: {} (node: {})", 
                 result.bestCosineScore, result.bestCosineNodeId);
    logging::info("Best Weighted score: {} (node: {})", 
                 result.bestWeightedScore, result.bestWeightedNodeId);
    logging::info("Best Raw matches: {} (node: {})", 
                 result.bestRawSeedMatchScore, result.bestRawSeedMatchNodeId);
}

} // namespace placement


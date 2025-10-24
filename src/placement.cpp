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
                               const placement::PlacementGlobalState& state);

void updateSimilarityMetricsDel(size_t hash,
                               placement::MgsrGenomeState& genomeState,
                               const placement::PlacementGlobalState& state);

void calculateMgsrStyleScores(const std::string& nodeId,
                             const placement::MgsrGenomeState& genomeState,
                             const placement::PlacementGlobalState& state,
                             placement::PlacementResult& result,
                             const std::unordered_set<size_t>& affectedHashes);

void updateMetricsAfterBacktrack(placement::MgsrGenomeState& genomeState,
                                const placement::PlacementGlobalState& state,
                                const std::unordered_set<size_t>& affectedHashes);

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
    // Process seed deletions
    auto deletions = nodeChanges.getSeedDeletions();
    for (uint32_t seedIdx : deletions) {
        if (seedIdx < state.seedInfo.size()) {
            auto seedReader = state.seedInfo[seedIdx];
            size_t seedHash = seedReader.getHash();
            auto it = genomeState.seedCounts.find(seedHash);
            if (it != genomeState.seedCounts.end()) {
                it->second--;
                if (it->second <= 0) {
                    genomeState.seedCounts.erase(it);
                    genomeState.uniqueSeeds.erase(seedHash);
                }
            }
        }
    }
    
    // Process seed additions (called seedInsubIndices in Cap'n Proto schema)
    auto additions = nodeChanges.getSeedInsubIndices();
    for (uint32_t seedIdx : additions) {
        if (seedIdx < state.seedInfo.size()) {
            auto seedReader = state.seedInfo[seedIdx];
            size_t seedHash = seedReader.getHash();
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
                         std::unordered_set<size_t>& affectedHashes) {
    
    auto seedReader = state.seedInfo[seedIndex];
    uint64_t pos = seedReader.getStartPos();
    size_t hash = seedReader.getHash();
    
    auto posIt = genomeState.positionMap.find(pos);
    if (posIt != genomeState.positionMap.end()) {
        // Position exists - substitution
        uint64_t oldSeedIndex = posIt->second;
        auto oldSeedReader = state.seedInfo[oldSeedIndex];
        size_t oldHash = oldSeedReader.getHash();
        
        backtrack.emplace_back(oldSeedIndex, 3);  // SUB
        affectedHashes.insert(oldHash);
        affectedHashes.insert(hash);
        
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
        affectedHashes.insert(hash);
        genomeState.hashToPositionMap[hash].push_back(posIt);
    }
    
    // Update similarity metrics
    updateSimilarityMetricsAdd(hash, genomeState, state);
}

// MGSR-style seed deletion (matches mgsr::mgsrPlacer::delSeedAtPosition)
void delSeedFromGenomeState(uint64_t pos,
                           placement::MgsrGenomeState& genomeState,
                           const placement::PlacementGlobalState& state,
                           std::vector<std::pair<uint64_t, uint8_t>>& backtrack,
                           std::unordered_set<size_t>& affectedHashes) {
    
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
    affectedHashes.insert(hash);
    
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
    updateSimilarityMetricsDel(hash, genomeState, state);
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
                               const placement::PlacementGlobalState& state) {    // Update seed counts
    auto [it, inserted] = genomeState.currentSeedCounts.try_emplace(hash, 0);
    int64_t oldGenomeCount = it->second;
    it->second++;
    int64_t newGenomeCount = it->second;
    
    genomeState.currentUniqueSeeds.insert(hash);
    
    // Update similarity metrics if seed is in reads
    auto readIt = state.seedFreqInReads.find(hash);
    if (readIt != state.seedFreqInReads.end()) {
        int64_t readCount = readIt->second;
        
        // Jaccard: if becoming present, add read frequency
        if (oldGenomeCount == 0 && newGenomeCount > 0) {
            genomeState.currentJaccardNumerator += readCount;
        }
        
        // Weighted Jaccard: update min(read, genome)
        int64_t oldMin = std::min(readCount, oldGenomeCount);
        int64_t newMin = std::min(readCount, newGenomeCount);
        genomeState.currentWeightedJaccardNumerator += (newMin - oldMin);
        
        // Cosine: update dot product
        double oldDot = static_cast<double>(readCount * oldGenomeCount);
        double newDot = static_cast<double>(readCount * newGenomeCount);
        genomeState.currentCosineNumerator += (newDot - oldDot);
    }
}

// Update similarity metrics when deleting a seed
void updateSimilarityMetricsDel(size_t hash,
                               placement::MgsrGenomeState& genomeState,
                               const placement::PlacementGlobalState& state) {
    
    auto countIt = genomeState.currentSeedCounts.find(hash);
    if (countIt == genomeState.currentSeedCounts.end()) return;
    
    int64_t oldGenomeCount = countIt->second;
    countIt->second--;
    int64_t newGenomeCount = countIt->second;
    
    if (newGenomeCount <= 0) {
        genomeState.currentSeedCounts.erase(countIt);
        genomeState.currentUniqueSeeds.erase(hash);
    }
    
    // Update similarity metrics if seed is in reads
    auto readIt = state.seedFreqInReads.find(hash);
    if (readIt != state.seedFreqInReads.end()) {
        int64_t readCount = readIt->second;
        
        // Jaccard: if becoming absent, subtract read frequency
        if (oldGenomeCount > 0 && newGenomeCount == 0) {
            genomeState.currentJaccardNumerator -= readCount;
        }
        
        // Weighted Jaccard: update min(read, genome)
        int64_t oldMin = std::min(readCount, oldGenomeCount);
        int64_t newMin = std::min(readCount, newGenomeCount);
        genomeState.currentWeightedJaccardNumerator += (newMin - oldMin);
        
        // Cosine: update dot product
        double oldDot = static_cast<double>(readCount * oldGenomeCount);
        double newDot = static_cast<double>(readCount * newGenomeCount);
        genomeState.currentCosineNumerator += (newDot - oldDot);
    }
}

// Calculate final similarity scores using MGSR-style approach
void calculateMgsrStyleScores(const std::string& nodeId,
                             const placement::MgsrGenomeState& genomeState,
                             const placement::PlacementGlobalState& state,
                             placement::PlacementResult& result,
                             const std::unordered_set<size_t>& affectedHashes) {
    
    // Calculate final scores from current accumulated metrics
    double jaccardScore = 0.0;
    double cosineScore = 0.0;
    double weightedJaccardScore = 0.0;
    double presenceJaccardScore = 0.0;
    
    // Regular Jaccard
    int64_t jaccardDenominator = 0;
    for (const auto& [hash, readCount] : state.seedFreqInReads) {
        jaccardDenominator += readCount;  // Union includes all read seeds
    }
    if (jaccardDenominator > 0) {
        jaccardScore = static_cast<double>(genomeState.currentJaccardNumerator) / jaccardDenominator;
    }
    
    // Weighted Jaccard  
    int64_t weightedJaccardDenominator = 0;
    for (const auto& [hash, readCount] : state.seedFreqInReads) {
        auto genomeIt = genomeState.currentSeedCounts.find(hash);
        int64_t genomeCount = (genomeIt != genomeState.currentSeedCounts.end()) ? genomeIt->second : 0;
        weightedJaccardDenominator += std::max(readCount, genomeCount);
    }
    // Add genome-only seeds
    for (const auto& [hash, genomeCount] : genomeState.currentSeedCounts) {
        if (state.seedFreqInReads.find(hash) == state.seedFreqInReads.end()) {
            weightedJaccardDenominator += genomeCount;
        }
    }
    if (weightedJaccardDenominator > 0) {
        weightedJaccardScore = static_cast<double>(genomeState.currentWeightedJaccardNumerator) / weightedJaccardDenominator;
    }
    
    // Cosine similarity
    double genomeMagnitude = 0.0;
    for (const auto& [hash, genomeCount] : genomeState.currentSeedCounts) {
        genomeMagnitude += static_cast<double>(genomeCount * genomeCount);
    }
    genomeMagnitude = std::sqrt(genomeMagnitude);
    
    if (state.readMagnitude > 0.0 && genomeMagnitude > 0.0) {
        cosineScore = genomeState.currentCosineNumerator / (state.readMagnitude * genomeMagnitude);
    }
    
    // Presence/absence Jaccard
    size_t intersectionCount = 0;
    for (const auto& [hash, _] : state.seedFreqInReads) {
        if (genomeState.currentUniqueSeeds.find(hash) != genomeState.currentUniqueSeeds.end()) {
            intersectionCount++;
        }
    }
    size_t unionSize = state.seedFreqInReads.size() + genomeState.currentUniqueSeeds.size() - intersectionCount;
    if (unionSize > 0) {
        presenceJaccardScore = static_cast<double>(intersectionCount) / unionSize;
    }
    
    // Update best scores in result
    result.updateJaccardScore(nodeId, jaccardScore);
    result.updateWeightedJaccardScore(nodeId, weightedJaccardScore);
    result.updateCosineScore(nodeId, cosineScore);
    result.updateJaccardPresenceScore(nodeId, presenceJaccardScore);
    result.updateRawSeedMatchScore(nodeId, genomeState.currentJaccardNumerator);
    result.updateHitsScore(nodeId, intersectionCount);
    
    logging::debug("Node {}: Jaccard={:.6f}, WeightedJaccard={:.6f}, Cosine={:.6f}, PresenceJaccard={:.6f}, intersect={}", 
                  nodeId, jaccardScore, weightedJaccardScore, cosineScore, presenceJaccardScore, intersectionCount);
}

// Update metrics after backtracking  
void updateMetricsAfterBacktrack(placement::MgsrGenomeState& genomeState,
                                const placement::PlacementGlobalState& state,
                                const std::unordered_set<size_t>& affectedHashes) {
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
            auto genomeSeedIt = result.currentGenomeSeedCounts.find(seedHash);
            int64_t oldGenomeCount = 0;
            int64_t newGenomeCount = 0;
            
            if (genomeSeedIt != result.currentGenomeSeedCounts.end()) {
                oldGenomeCount = genomeSeedIt->second;
                newGenomeCount = oldGenomeCount - 1;
                
                // Update weighted Jaccard numerator: subtract old contribution, add new contribution
                int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                int64_t newMinCount = (newGenomeCount > 0) ? std::min(readCount, newGenomeCount) : 0;
                result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                
                if (newGenomeCount <= 0) { 
                    result.currentGenomeSeedCounts.erase(genomeSeedIt);
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
                auto genomeSeedIt = result.currentGenomeSeedCounts.find(seedHash);
                int64_t oldGenomeCount = 0;
                int64_t newGenomeCount = 0;
                
                if (genomeSeedIt != result.currentGenomeSeedCounts.end()) {
                    oldGenomeCount = genomeSeedIt->second;
                    newGenomeCount = oldGenomeCount - 1;
                    
                    // Update weighted Jaccard numerator: subtract old contribution, add new contribution
                    int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                    int64_t newMinCount = (newGenomeCount > 0) ? std::min(readCount, newGenomeCount) : 0;
                    result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                    
                    if (newGenomeCount <= 0) {
                        result.currentGenomeSeedCounts.erase(genomeSeedIt);
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
        
        // Safely increment count in currentGenomeSeedCounts
        auto [it, inserted] = result.currentGenomeSeedCounts.try_emplace(newSeed.hash, 0);
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
                    
                    // Update metrics if this seed was in reads
                    auto readIt = state.seedFreqInReads.find(seedHash);
                    if (readIt != state.seedFreqInReads.end()) {
                        int64_t readCount = readIt->second;
                        
                        // Get current genome count before modification
                        auto genomeIt = result.currentGenomeSeedCounts.find(seedHash);
                        int64_t oldGenomeCount = (genomeIt != result.currentGenomeSeedCounts.end()) ? genomeIt->second : 0;
                        int64_t newGenomeCount = std::max(0L, oldGenomeCount - 1);
                        
                        // Update weighted Jaccard numerator: subtract old contribution, add new contribution
                        int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                        int64_t newMinCount = std::min(readCount, newGenomeCount);
                        result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                        
                        // Update regular Jaccard numerator (if seed becomes absent, decrease by readCount)
                        if (oldGenomeCount > 0 && newGenomeCount == 0) {
                            result.currentJaccardNumerator -= readCount;
                        }
                        
                        // Update raw match score: subtract old contribution, add new contribution
                        result.hitsInThisGenome += (readCount * newGenomeCount) - (readCount * oldGenomeCount);
                        
                        // Update cosine similarity components (numerator: readCount * genomeCount)
                        double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount);
                        double newCosineTerm = static_cast<double>(readCount * newGenomeCount);
                        result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
                        
                        // Update genome seed counts
                        if (genomeIt != result.currentGenomeSeedCounts.end()) {
                            if (--(genomeIt->second) <= 0) {
                                result.currentGenomeSeedCounts.erase(genomeIt);
                            }
                        }
                    }
                } else {
                    // Process insertion/substitution
                    // Add to unique seed hashes
                    uniqueSeedHashes.insert(seedHash);
                    
                    // Update metrics if this seed exists in reads
                    auto readIt = state.seedFreqInReads.find(seedHash);
                    if (readIt != state.seedFreqInReads.end()) {
                        int64_t readCount = readIt->second;
                        
                        // Get current genome count before modification
                        auto genomeIt = result.currentGenomeSeedCounts.find(seedHash);
                        int64_t oldGenomeCount = (genomeIt != result.currentGenomeSeedCounts.end()) ? genomeIt->second : 0;
                        int64_t newGenomeCount = oldGenomeCount + 1;
                        
                        // Update weighted Jaccard numerator: subtract old contribution, add new contribution
                        int64_t oldMinCount = std::min(readCount, oldGenomeCount);
                        int64_t newMinCount = std::min(readCount, newGenomeCount);
                        result.currentWeightedJaccardNumerator += (newMinCount - oldMinCount);
                        
                        // Update regular Jaccard numerator (if seed becomes present, increase by readCount)
                        if (oldGenomeCount == 0 && newGenomeCount > 0) {
                            result.currentJaccardNumerator += readCount;
                        }
                        
                        // Update raw match score: subtract old contribution, add new contribution
                        result.hitsInThisGenome += (readCount * newGenomeCount) - (readCount * oldGenomeCount);
                        
                        // Update cosine similarity components (numerator: readCount * genomeCount)
                        double oldCosineTerm = static_cast<double>(readCount * oldGenomeCount);
                        double newCosineTerm = static_cast<double>(readCount * newGenomeCount);
                        result.currentCosineNumerator += (newCosineTerm - oldCosineTerm);
                        
                        // Update genome seed counts
                        auto [it, inserted] = result.currentGenomeSeedCounts.try_emplace(seedHash, 0);
                        it->second++;
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
        auto genomeIt = result.currentGenomeSeedCounts.find(seedHash);
        if (genomeIt != result.currentGenomeSeedCounts.end() && genomeIt->second > 0) {
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
        auto genomeIt = result.currentGenomeSeedCounts.find(seedHash);
        int64_t genomeCount = (genomeIt != result.currentGenomeSeedCounts.end()) ? genomeIt->second : 0;
        weightedJaccardDenominator += std::max(readCount, genomeCount);
    }
    
    // Add contributions from genome seeds not in reads
    for (const auto& [seedHash, genomeCount] : result.currentGenomeSeedCounts) {
        if (state.seedFreqInReads.find(seedHash) == state.seedFreqInReads.end()) {
            weightedJaccardDenominator += genomeCount;
        }
    }
    
    double weightedJaccardScore = 0.0;
    if (weightedJaccardDenominator > 0) {
        weightedJaccardScore = static_cast<double>(result.currentWeightedJaccardNumerator) / static_cast<double>(weightedJaccardDenominator);
    }
    
    // Calculate cosine similarity
    // Compute current genome magnitude from genomeSeedCounts
    double genomeMagnitude = 0.0;
    for (const auto& [seedHash, genomeCount] : result.currentGenomeSeedCounts) {
        genomeMagnitude += static_cast<double>(genomeCount * genomeCount);
    }
    genomeMagnitude = std::sqrt(genomeMagnitude);
    
    // Use precomputed read magnitude
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

// Calculate MGSR-style scores and update result if better
void calculateMgsrStyleScores(const std::string& nodeId,
                             const placement::MgsrGenomeState& genomeState,
                             const placement::PlacementGlobalState& state,
                             placement::PlacementResult& result,
                             const std::unordered_set<size_t>& affectedHashes) {
    
    // Calculate scores (similar to original but using MgsrGenomeState)
    double jaccardScore = 0.0;
    double cosineScore = 0.0;
    size_t intersectionCount = 0;
    
    if (!state.seedFreqInReads.empty()) {
        // Calculate intersection
        for (const auto& [hash, genomeCount] : genomeState.currentSeedCounts) {
            auto readIt = state.seedFreqInReads.find(hash);
            if (readIt != state.seedFreqInReads.end()) {
                intersectionCount += std::min(genomeCount, readIt->second);
            }
        }
        
        // Calculate union for Jaccard
        std::unordered_set<size_t> unionHashes;
        for (const auto& [hash, count] : genomeState.currentSeedCounts) {
            unionHashes.insert(hash);
        }
        for (const auto& [hash, count] : state.seedFreqInReads) {
            unionHashes.insert(hash);
        }
        
        if (!unionHashes.empty()) {
            jaccardScore = static_cast<double>(intersectionCount) / unionHashes.size();
        }
        
        // Calculate cosine similarity
        double genomeMagnitude = 0.0;
        for (const auto& [hash, count] : genomeState.currentSeedCounts) {
            genomeMagnitude += count * count;
        }
        genomeMagnitude = std::sqrt(genomeMagnitude);
        
        if (genomeMagnitude > 0.0 && state.readMagnitude > 0.0) {
            double dotProduct = 0.0;
            for (const auto& [hash, genomeCount] : genomeState.currentSeedCounts) {
                auto readIt = state.seedFreqInReads.find(hash);
                if (readIt != state.seedFreqInReads.end()) {
                    dotProduct += genomeCount * readIt->second;
                }
            }
            cosineScore = dotProduct / (genomeMagnitude * state.readMagnitude);
        }
    }
    
    // Update best scores
    if (jaccardScore > result.bestJaccardScore) {
        result.bestJaccardScore = jaccardScore;
        result.bestJaccardNodeId = nodeId;
    }
    
    if (cosineScore > result.bestCosineScore) {
        result.bestCosineScore = cosineScore;
        result.bestCosineNodeId = nodeId;
    }
}

// MGSR-style recursive placement helper (mimics MGSR's placeReadsHelper)
void placeLiteHelper(panmapUtils::LiteNode* node, 
                    placement::MgsrGenomeState& genomeState,
                    uint32_t& dfsIndex,
                    const placement::PlacementGlobalState& state,
                    placement::PlacementResult& result,
                    const placement::TraversalParams& params) {
    
    if (!node) return;
    
    const std::string& nodeId = node->identifier;
    
    // Debug logging
    if (!params.debug_node_id.empty() && nodeId == params.debug_node_id) {
        logging::info("DEBUG: Processing debug node: {} at DFS index {}", nodeId, dfsIndex);
    }
    
    // *** FORWARD STEP: Apply mutations for this node ***
    std::vector<std::pair<uint64_t, uint8_t>> seedBacktracks;  // (seedIndex or position, changeType)
    std::unordered_set<size_t> affectedSeedHashes;
    
    if (dfsIndex < state.perNodeChanges.size()) {
        auto nodeChanges = state.perNodeChanges[dfsIndex];
        
        // Apply seed deltas
        auto seedDeltas = nodeChanges.getSeedDeltas();
        for (auto delta : seedDeltas) {
            if (delta.getIsDeleted()) {
                // Apply seed deletions
                uint32_t deletionPos = delta.getSeedIndex();
                delSeedFromGenomeState(deletionPos, genomeState, state, seedBacktracks, affectedSeedHashes);
            } else {
                // Apply seed additions
                uint32_t seedIdx = delta.getSeedIndex();
                if (seedIdx < state.seedInfo.size()) {
                    addSeedToGenomeState(seedIdx, genomeState, state, seedBacktracks, affectedSeedHashes);
                }
            }
        }
    }
    
    // *** CALCULATE SCORES for current node ***
    calculateMgsrStyleScores(nodeId, genomeState, state, result, affectedSeedHashes);
    
    // *** RECURSE TO CHILDREN ***
    for (panmapUtils::LiteNode* child : node->children) {
        dfsIndex++;
        placeLiteHelper(child, genomeState, dfsIndex, state, result, params);
    }
    
    // *** BACKWARD STEP: Backtrack mutations (restore parent state) ***
    for (auto it = seedBacktracks.rbegin(); it != seedBacktracks.rend(); ++it) {
        const auto& [seedIndexOrPos, changeType] = *it;
        
        if (changeType == 1) {  // Was addition, now remove
            auto seedReader = state.seedInfo[seedIndexOrPos];
            uint64_t pos = seedReader.getStartPos();
            delSeedFromGenomeStateBacktrack(pos, genomeState, state);
        } else if (changeType == 2) {  // Was deletion, now add back
            addSeedToGenomeStateBacktrack(seedIndexOrPos, genomeState, state);
        } else if (changeType == 3) {  // Was substitution, restore old seed
            // For substitution backtrack, seedIndexOrPos is the old seed index
            addSeedToGenomeStateBacktrack(seedIndexOrPos, genomeState, state);
        }
    }
}// MGSR-style placement with recursive DFS and backtracking
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
                       const std::string &debug_node_id_param) {
    
    logging::info("Starting MGSR-style recursive placement");
    
    auto time_init_start = std::chrono::high_resolution_clock::now();
    
    // Initialize placement global state
    PlacementGlobalState state;
    state.seedInfo = mgsrIndex.getSeedInfo();
    state.perNodeChanges = mgsrIndex.getPerNodeChanges();
    state.kmerSize = mgsrIndex.getK();
    
    // Set traversal parameters
    TraversalParams params;
    params.k = mgsrIndex.getK();
    params.s = mgsrIndex.getS();
    params.t = mgsrIndex.getT();
    params.open = mgsrIndex.getOpen();
    params.debug_node_id = debug_node_id_param;
    
    auto time_init_end = std::chrono::high_resolution_clock::now();
    auto duration_init = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_init_end - time_init_start);
    logging::info("[TIME] Placement initialization: {}ms", duration_init.count());
    
    // Process reads and extract seeds (same as before)
    auto time_read_start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> allReadSequences;
    
    if (!reads1.empty()) {
        mgsr::extractReadSequences(reads1, reads2, allReadSequences);
        
        absl::flat_hash_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
        std::vector<std::string> tempReadQuals;
        std::vector<std::string> tempReadNames;
        std::vector<std::vector<seeding::seed_t>> tempReadSeeds;
        std::vector<std::vector<std::string>> tempReadSeedSeqs;
        
        seeding::seedsFromFastq(params.k, params.s, params.t, params.open, 0, 
                      readSeedCounts, allReadSequences, tempReadQuals, 
                      tempReadNames, tempReadSeeds, tempReadSeedSeqs, reads1, reads2);
        
        for (const auto& [seedHash, counts] : readSeedCounts) {
            state.seedFreqInReads[seedHash] = counts.first + counts.second;
        }
        
        logging::info("Extracted {} unique seeds from reads", state.seedFreqInReads.size());
    } else if (!readSequences.empty()) {
        allReadSequences = readSequences;
        
        for (const std::string& seq : allReadSequences) {
            for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : seeding::rollingSyncmers(seq, params.k, params.s, params.open, params.t, false)) {
                if (!isSyncmer) continue;
                state.seedFreqInReads[kmerHash]++;
            }
        }
        
        logging::info("Extracted {} unique seeds from provided sequences", state.seedFreqInReads.size());
    }
    
    auto time_read_end = std::chrono::high_resolution_clock::now();
    auto duration_read = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_read_end - time_read_start);
    logging::info("[TIME] Read processing & seed extraction: {}ms", duration_read.count());
    
    // Precompute read magnitude for cosine similarity
    auto time_precomp_start = std::chrono::high_resolution_clock::now();
    for (const auto& [seedHash, readCount] : state.seedFreqInReads) {
        state.readMagnitude += static_cast<double>(readCount * readCount);
    }
    state.readMagnitude = std::sqrt(state.readMagnitude);
    logging::info("Precomputed read magnitude: {:.6f}", state.readMagnitude);
    
    auto time_precomp_end = std::chrono::high_resolution_clock::now();
    auto duration_precomp = std::chrono::duration_cast<std::chrono::milliseconds>(
        time_precomp_end - time_precomp_start);
    if (duration_precomp.count() > 0) {
        logging::info("[TIME] Magnitude precomputation: {}ms", duration_precomp.count());
    }
    
    // Set root pointer in state
    state.root = liteTree->root;
    
    // Initialize MGSR-style global genome state (single instance, modified in-place)
    placement::MgsrGenomeState globalGenomeState;
    uint32_t currentDfsIndex = 0;
    
    // Initialize root state if root has mutations (DFS index 0)
    if (liteTree && liteTree->root && state.perNodeChanges.size() > 0) {
        auto rootChanges = state.perNodeChanges[0];
        
        // Add root seeds (no backtracking needed for root)
        std::vector<std::pair<uint64_t, uint8_t>> dummy_backtrack;
        std::unordered_set<size_t> dummy_affected;
        
        auto seedDeltas = rootChanges.getSeedDeltas();
        for (auto delta : seedDeltas) {
            if (!delta.getIsDeleted()) {  // If not deleted, it's an addition
                uint32_t seedIdx = delta.getSeedIndex();
                if (seedIdx < state.seedInfo.size()) {
                    addSeedToGenomeState(seedIdx, globalGenomeState, state, dummy_backtrack, dummy_affected);
                }
            }
        }
        
        logging::info("Initialized root '{}' with {} seeds, {} unique hashes", 
                     liteTree->root->identifier,
                     globalGenomeState.positionMap.size(),
                     globalGenomeState.currentSeedCounts.size());
    }
    
    logging::info("Starting MGSR-style recursive placement traversal with {} tree nodes", 
                  state.perNodeChanges.size());
    
    // Start recursive traversal from root
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
    
    logging::info("MGSR-style placement completed. Best Jaccard score: {} (node: {})", 
                 result.bestJaccardScore, result.bestJaccardNodeId);
    logging::info("Best Cosine score: {} (node: {})", 
                 result.bestCosineScore, result.bestCosineNodeId);
    logging::info("Best Weighted score: {} (node: {})", 
                 result.bestWeightedScore, result.bestWeightedNodeId);
}

} // namespace placement


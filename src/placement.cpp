#include "panman.hpp"
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include "placement.hpp"
#include "caller_logging.hpp"
#include "index.capnp.h"
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

#include "progress_state.hpp"
#include "indexing.hpp"

using namespace coordinates;
using namespace logging;
using namespace state;

namespace placement {

using panmanUtils::Node;
using panmanUtils::Tree;

// Forward declarations
class PlacementResult;

// Shared progress state for UI updates
std::shared_ptr<PlacementProgressState> progress_state;

// Maximum number of nodes for which to dump recomputation range details
const size_t MAX_NODES_TO_DUMP = 5;

// Helper function for calculating cosine similarity delta
std::pair<double, double> getCosineDelta(
    bool isRemoval, 
    bool isAddition,
    size_t seedHash, 
    int64_t count,
    const std::unordered_map<size_t, size_t>& readSeedCounts,
    const std::unordered_map<size_t, int64_t>& genomeSeedCounts) {

  // Default values if seed not found in reads
  double numeratorDelta = 0.0;
  double sumOfSquaresDelta = 0.0;

  // Check if this seed exists in the read set
  const auto readIt = readSeedCounts.find(seedHash);
  if (readIt != readSeedCounts.end()) {
      const auto readCount = readIt->second;
      
      if (isRemoval) {
          // Seed is being removed - decrease similarity
          numeratorDelta = -static_cast<double>(readCount);
          sumOfSquaresDelta = -1.0;  // One seed removed from genome set
      } else if (isAddition) {
          // Seed is being added - increase similarity
          numeratorDelta = static_cast<double>(readCount);
          sumOfSquaresDelta = 1.0;   // One seed added to genome set
      }
  }

  return {numeratorDelta, sumOfSquaresDelta};
}

// Process mutations for a single node and update scores
void processNodeMutations(
    panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    const TraversalParams& params) {
    
    if (!node) return;
    
    // Get node's seed mutations from index using DFS index
    const int64_t dfsIndex = stateManager.getDfsIndex(node->identifier);
    if (dfsIndex < 0 || dfsIndex >= static_cast<int64_t>(state.perNodeSeedMutations.size())) {
        throw std::runtime_error("No mutations found for node: " + node->identifier);
    }
    
    // Get seed mutations for this node
    const auto seedMutation = state.perNodeSeedMutations[dfsIndex];
    
    // Access base positions and masks
    const auto basePositions = seedMutation.getBasePositions();
    const auto perPosMasks = seedMutation.getPerPosMasks();
    
    // Check for dictionary-based k-mers
    bool hasDictionaryKmers = seedMutation.hasKmerDictionaryIds() && \
                             seedMutation.hasKmerPositions() && \
                             !state.kmerDictionary.empty();
    
    // Get dictionary data if available
    auto kmerDictionaryIds = hasDictionaryKmers ? \
                           seedMutation.getKmerDictionaryIds() : \
                           ::capnp::List<uint32_t>::Reader();
    auto kmerPositions = hasDictionaryKmers ? \
                       seedMutation.getKmerPositions() : \
                       ::capnp::List<int64_t>::Reader();
    
    // Stats tracking
    size_t seedAdditions = 0;
    size_t seedDeletions = 0;
    
    // Group positions by block for better memory locality and cache efficiency
    std::unordered_map<int32_t, std::vector<int64_t>> positionsByBlock;
    std::vector<std::pair<int64_t, uint64_t>> positionsWithValues;
    
    // Constants for quaternary encoding
    constexpr uint8_t BITS_PER_VALUE = 2;
    constexpr uint8_t POSITIONS_PER_MASK = 32;
    constexpr uint8_t VALUE_MASK = 0x3;
    
    // First pass: decode quaternary masks to individual position-value pairs
    for (size_t i = 0; i < basePositions.size(); i++) {
        const int64_t basePos = basePositions[i];
        const uint64_t bitMask = perPosMasks[i];
        
        // Skip empty masks
        if (bitMask == 0) continue;
        
        // Decode up to 32 positions from each mask (2 bits per position)
        for (uint8_t offset = 0; offset < POSITIONS_PER_MASK; offset++) {
            // Extract quaternary value (2 bits)
            const uint8_t value = (bitMask >> (offset * BITS_PER_VALUE)) & VALUE_MASK;
            
            // Skip if no change (value 0 = unchanged seed status)
            if (value == 0) continue;
            
            // Calculate actual position
            const int64_t pos = basePos - offset;
            
            // Store position and value for processing
            positionsWithValues.emplace_back(pos, value);
            
            // Get the block this position belongs to
            auto blockCoord = stateManager.mapGlobalToBlockCoords(node->identifier, pos);
            if (blockCoord.blockId >= 0) {
                positionsByBlock[blockCoord.blockId].push_back(pos);
            }
        }
    }
    
    // Sort positions within each block for better sequential access
    for (auto& [blockId, positions] : positionsByBlock) {
        std::sort(positions.begin(), positions.end());
    }
    
    // Helper function to find position index in extracted sequence
    auto findPositionIndex = [](const std::vector<int64_t>& positions, int64_t targetPos) -> int {
        for (size_t i = 0; i < positions.size(); ++i) {
            if (positions[i] == targetPos) return static_cast<int>(i);
        }
        return -1;
    };
    
    // Helper function to create a new seed at a position
    auto createSeedAt = [&](int64_t pos, 
                           int posIndex, 
                           const std::string& sequence, 
                           const std::vector<int64_t>& positions, 
                           int32_t blockId) -> std::optional<size_t> {
        
        // First check if we have this position in the k-mer dictionary
        if (hasDictionaryKmers) {
            for (size_t i = 0; i < kmerPositions.size(); i++) {
                if (kmerPositions[i] == pos) {
                    // Found the k-mer in dictionary, get its ID
                    uint32_t kmerId = kmerDictionaryIds[i];
                    
                    // Look up k-mer sequence in dictionary
                    auto dictIt = state.kmerDictionary.find(kmerId);
                    if (dictIt != state.kmerDictionary.end()) {
                        // Use the k-mer from dictionary directly
                        const std::string& kmerStr = dictIt->second;
                        
                        // Hash the k-mer
                        auto hashResult = seeding::hashSeq(kmerStr, params.k);
                        
                        // Create and store the seed
                        seeding::seed_t newSeed{
                            .hash = hashResult.first,
                            .reversed = static_cast<bool>(hashResult.second),
                            .endPos = pos + params.k - 1 // Approximate end position
                        };
                        
                        stateManager.setSeedAtPosition(pos, newSeed);
                        stateManager.addSeedToBlock(blockId, pos);
                        seedAdditions++;
                        
                        return hashResult.first;
                    }
                }
            }
        }
        
        // Skip if position is invalid or k-mer won't fit
        if (posIndex < 0 || posIndex + params.k > sequence.length()) {
            return std::nullopt;
        }
        
        // Fall back to extracting k-mer from sequence
        std::string kmerStr = sequence.substr(posIndex, params.k);
        
        // Hash the k-mer
        auto hashResult = seeding::hashSeq(kmerStr, params.k);
        
        // Create and store the seed
        seeding::seed_t newSeed{
            .hash = hashResult.first,
            .reversed = static_cast<bool>(hashResult.second),
            .endPos = positions[posIndex + params.k - 1]
        };
        
        stateManager.setSeedAtPosition(pos, newSeed);
        stateManager.addSeedToBlock(blockId, pos);
        seedAdditions++;
        
        return hashResult.first;
    };
    
    // Process each block's mutations
    for (const auto& [blockId, positions] : positionsByBlock) {
        // Skip empty blocks
        if (positions.empty()) continue;
        
        // Get block range and validate
        const auto blockRange = stateManager.getBlockRange(blockId);
        if (blockRange.start < 0 || blockRange.end <= blockRange.start) continue;
        
        // Determine sequence extraction range with padding for k-mer context
        const int64_t minPos = positions.front();
        const int64_t maxPos = positions.back();
        const int64_t startPos = std::max(int64_t(0), minPos - params.k);
        const int64_t endPos = maxPos;
        const coordinates::CoordRange extractRange{startPos, endPos};
        
        try {
            // Extract sequence once for the entire block
            auto [blockSequence, blockPositions] = stateManager.extractSequence(
                node->identifier, extractRange, true);
            
            // Process each position-value pair in this block
            for (const auto& [pos, value] : positionsWithValues) {
                // Skip positions that don't belong to current block
                auto blockCoord = stateManager.mapGlobalToBlockCoords(node->identifier, pos);
                if (blockCoord.blockId != blockId) continue;
                
                // Find position in the extracted sequence
                const int posIndex = findPositionIndex(blockPositions, pos);
                
                // Handle based on quaternary value
                switch (value) {
                    case 1: {  // Deletion (on → off)
                        auto existingSeed = stateManager.getSeedAtPosition(pos);
                        if (!existingSeed) continue;  // Skip if no seed to delete
                        
                        const size_t seedHash = existingSeed->hash;
                        
                        // Update scores if seed exists in reads
                        auto readIt = state.seedFreqInReads.find(seedHash);
                        if (readIt != state.seedFreqInReads.end()) {
                            const int64_t readCount = readIt->second;
                            
                            // Update metrics: hits, Jaccard, genome seeds
                            result.hitsInThisGenome -= readCount;
                            result.currentJaccardNumerator -= readCount;
                            
                            auto& seedCount = result.currentGenomeSeedCounts[seedHash];
                            seedCount--;
                            if (seedCount <= 0) {
                                result.currentGenomeSeedCounts.erase(seedHash);
                            }
                            
                            // Update cosine similarity components
                            auto [numeratorDelta, denomDelta] = getCosineDelta(
                                true, false, seedHash, 0, 
                                state.seedFreqInReads, 
                                result.currentGenomeSeedCounts);
                            
                            result.currentCosineNumerator += numeratorDelta;
                            result.currentCosineDenominator += denomDelta;
                        }
                        
                        // Remove the seed
                        stateManager.clearSeedAtPosition(pos);
                        seedDeletions++;
                        break;
                    }
                    
                    case 2: {  // Addition (off → on)
                        // Create new seed and get its hash
                        auto seedHashOpt = createSeedAt(pos, posIndex, blockSequence, blockPositions, blockId);
                        if (!seedHashOpt) continue;
                        
                        const size_t seedHash = *seedHashOpt;
                        
                        // Update scores if seed exists in reads
                        auto readIt = state.seedFreqInReads.find(seedHash);
                        if (readIt != state.seedFreqInReads.end()) {
                            const int64_t readCount = readIt->second;
                            
                            // Update metrics: hits, Jaccard, genome seeds
                            result.hitsInThisGenome += readCount;
                            result.currentJaccardNumerator += readCount;
                            result.currentGenomeSeedCounts[seedHash]++;
                            
                            // Update cosine similarity components
                            auto [numeratorDelta, denomDelta] = getCosineDelta(
                                false, true, seedHash, 1, 
                                state.seedFreqInReads, 
                                result.currentGenomeSeedCounts);
                            
                            result.currentCosineNumerator += numeratorDelta;
                            result.currentCosineDenominator += denomDelta;
                        }
                        break;
                    }
                    
                    case 3: {  // Modification (on → on')
                        // Check for existing seed
                        auto existingSeed = stateManager.getSeedAtPosition(pos);
                        if (!existingSeed) {
                            // Fall back to treating as addition if seed missing
                            auto seedHashOpt = createSeedAt(pos, posIndex, blockSequence, blockPositions, blockId);
                            if (!seedHashOpt) continue;
                            
                            const size_t seedHash = *seedHashOpt;
                            
                            // Update scores for addition
                            auto readIt = state.seedFreqInReads.find(seedHash);
                            if (readIt != state.seedFreqInReads.end()) {
                                const int64_t readCount = readIt->second;
                                result.hitsInThisGenome += readCount;
                                result.currentJaccardNumerator += readCount;
                                result.currentGenomeSeedCounts[seedHash]++;
                                
                                // Update cosine similarity components
                                auto [numeratorDelta, denomDelta] = getCosineDelta(
                                    false, true, seedHash, 1, 
                                    state.seedFreqInReads, 
                                    result.currentGenomeSeedCounts);
                                
                                result.currentCosineNumerator += numeratorDelta;
                                result.currentCosineDenominator += denomDelta;
                            }
                            continue;
                        }
                        
                        // Get hash of existing seed for removal
                        const size_t oldSeedHash = existingSeed->hash;
                        
                        // Remove old seed from scores
                        auto oldReadIt = state.seedFreqInReads.find(oldSeedHash);
                        if (oldReadIt != state.seedFreqInReads.end()) {
                            const int64_t oldReadCount = oldReadIt->second;
                            
                            // Update hit count for removal
                            result.hitsInThisGenome -= oldReadCount;
                            result.currentJaccardNumerator -= oldReadCount;
                            
                            // Update genome seed counts
                            auto& seedCount = result.currentGenomeSeedCounts[oldSeedHash];
                            seedCount--;
                            if (seedCount <= 0) {
                                result.currentGenomeSeedCounts.erase(oldSeedHash);
                            }
                        }
                        
                        // Create new seed and get its hash
                        auto seedHashOpt = createSeedAt(pos, posIndex, blockSequence, blockPositions, blockId);
                        if (!seedHashOpt) continue;
                        
                        const size_t newSeedHash = *seedHashOpt;
                        
                        // Update scores for new seed
                        auto newReadIt = state.seedFreqInReads.find(newSeedHash);
                        if (newReadIt != state.seedFreqInReads.end()) {
                            const int64_t newReadCount = newReadIt->second;
                            
                            // Update hit count for addition
                            result.hitsInThisGenome += newReadCount;
                            result.currentJaccardNumerator += newReadCount;
                            
                            // Update genome seed counts
                            result.currentGenomeSeedCounts[newSeedHash]++;
                        }
                        break;
                    }
                    
                    default:
                        // Should never happen due to bitmasking to 2 bits
                        break;
                }
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Error extracting sequence for node " + 
                                  node->identifier + ": " + e.what());
        }
    }
    
    // Update Jaccard denominator after all changes
    result.currentJaccardDenominator = state.jaccardDenominator + 
        result.currentGenomeSeedCounts.size() - result.currentJaccardNumerator;
    
    // Update score tracking for the node
    result.updateHitsScore(node, result.hitsInThisGenome);
    
    // Calculate and update Jaccard score
    double jaccardScore = 0.0;
    if (result.currentJaccardDenominator > 0) {
        jaccardScore = static_cast<double>(result.currentJaccardNumerator) / 
                     result.currentJaccardDenominator;
        result.updateJaccardScore(node, jaccardScore);
    }
    
    // Calculate and update cosine similarity
    double cosineScore = 0.0;
    if (result.currentCosineDenominator > 0 && state.totalReadSeedCount > 0) {
        const double normGenome = std::sqrt(std::abs(result.currentCosineDenominator));
        const double normReads = std::sqrt(state.totalReadSeedCount);
        cosineScore = result.currentCosineNumerator / (normGenome * normReads);
        result.updateCosineScore(node, cosineScore);
    }
    
    // Calculate weighted score based on parameters
    const double weightedScore = params.scoreScale * jaccardScore + 
                               (1.0 - params.scoreScale) * cosineScore;
    result.updateWeightedScore(node, weightedScore, params.scoreScale);
}

// Main placement traversal function
void placementTraversal(
    state::StateManager& stateManager,
    PlacementResult& result,
    panmanUtils::Tree* T,
    PlacementGlobalState& state,
    const TraversalParams& params) {
    
    if (!T || !T->root) {
        throw std::invalid_argument("Invalid tree or root node for placement traversal");
    }
    
    // Setup progress tracking
    std::atomic<size_t> nodesProcessed{0};
    const size_t totalNodes = T->allNodes.size();
    const auto startTime = std::chrono::high_resolution_clock::now();
    
    logging::info("Pre-computing paths for all nodes...");
    // Pre-compute paths for all nodes using common function from state
    const auto nodePaths = state::computeNodePaths(T, T->root);
    
    logging::info("Grouping nodes by level...");
    // Group nodes by level using common function from state
    const auto nodesByLevel = state::groupNodesByLevel(T, T->root);
    
    // Track processed nodes and previous node for common ancestor optimization
    std::unordered_set<std::string> processedNodes;
    panmanUtils::Node* previousNode = nullptr;
    std::mutex previousNodeMutex; // Protect access to previousNode
    std::mutex processedNodesMutex; // Protect access to processedNodes
    
    // Track errors during parallel processing
    std::atomic<bool> encounteredErrors{false};
    std::vector<std::string> errorMessages;
    std::mutex errorMessagesMutex;
    
    logging::info("Starting placement traversal with {} nodes across {} levels", 
                 totalNodes, nodesByLevel.size());
    
    // Process each level in order (breadth-first)
    for (size_t level = 0; level < nodesByLevel.size(); ++level) {
        // Check if we should abort traversal due to errors
        if (encounteredErrors) {
            logging::warn("Aborting traversal due to errors in previous level");
            break;
        }
        
        const auto& currentLevelNodes = nodesByLevel[level];
        logging::info("Processing level {}: {} nodes", level, currentLevelNodes.size());
        
        // Clear caches periodically to manage memory usage
        if (level > 0 && level % 5 == 0) {
            stateManager.clearCaches();
        }
        
        // Define the node processing function to be run in parallel 
        auto processNode = [&](panmanUtils::Node* currentNode) {
            try {
                // Skip if already processed
                {
                    std::lock_guard<std::mutex> lock(processedNodesMutex);
                    if (processedNodes.contains(currentNode->identifier)) {
                        return;
                    }
                }
                
                // Get path from root to current node
                const auto pathIt = nodePaths.find(currentNode->identifier);
                if (pathIt == nodePaths.end() || pathIt->second.empty()) {
                    throw std::runtime_error("No path found for node " + currentNode->identifier);
                }
                const auto& pathFromRoot = pathIt->second;
                
                // Find common ancestor with previous node for optimization
                panmanUtils::Node* commonAncestor = nullptr;
                {
                    std::lock_guard<std::mutex> lock(previousNodeMutex);
                    if (previousNode) {
                        const auto prevPathIt = nodePaths.find(previousNode->identifier);
                        if (prevPathIt != nodePaths.end()) {
                            const auto& prevPath = prevPathIt->second;
                            const size_t minPathSize = std::min(pathFromRoot.size(), prevPath.size());
                            
                            // Find the point where paths diverge
                            for (size_t i = 0; i < minPathSize; i++) {
                                if (pathFromRoot[i] != prevPath[i]) break;
                                commonAncestor = pathFromRoot[i];
                            }
                        }
                    }
                }
                
                // Initialize state for this node
                stateManager.initializeNode(currentNode->identifier);
                
                // If no common ancestor, propagate state from parent
                if (currentNode != T->root && currentNode->parent && commonAncestor == nullptr) {
                    stateManager.propagateState(currentNode->parent->identifier, currentNode->identifier);
                }
                
                // Process seed mutations for this node
                processNodeMutations(currentNode, stateManager, state, result, params);
                
                // Mark as processed and update previous node
                {
                    std::lock_guard<std::mutex> lock(processedNodesMutex);
                    processedNodes.insert(currentNode->identifier);
                    previousNode = currentNode;
                    
                    // Update progress tracking
                    nodesProcessed.fetch_add(1, std::memory_order_relaxed);
                }
            }
            catch (const std::exception& e) {
                // Log error and continue with other nodes
                const std::string errorMsg = "Error processing node " + 
                    currentNode->identifier + ": " + e.what();
                logging::err(errorMsg);
                
                {
                    std::lock_guard<std::mutex> lock(errorMessagesMutex);
                    errorMessages.push_back(errorMsg);
                }
                
                // Set error flag for non-critical errors
                if (std::string(e.what()).find("CRITICAL") != std::string::npos) {
                    encounteredErrors.store(true, std::memory_order_relaxed);
                }
            }
        };
        
        // Process nodes at this level, choosing sequential or parallel based on count
        constexpr size_t PARALLEL_THRESHOLD = 4;
        if (currentLevelNodes.size() <= PARALLEL_THRESHOLD) {
            // Process sequentially for small levels
            for (auto* node : currentLevelNodes) {
                processNode(node);
                
                // Check for critical errors after each node in sequential mode
                if (encounteredErrors) break;
            }
        } else {
            // Use parallel processing for larger levels
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, currentLevelNodes.size()),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        // Skip if critical errors encountered
                        if (encounteredErrors) break;
                        
                        processNode(currentLevelNodes[i]);
                    }
                },
                tbb::auto_partitioner());
        }
        
        // If we encountered errors, log them all at once
        if (!errorMessages.empty()) {
            std::lock_guard<std::mutex> lock(errorMessagesMutex);
            logging::warn("Encountered {} errors during level {} processing", 
                        errorMessages.size(), level);
            for (const auto& msg : errorMessages) {
                logging::debug("  Error: {}", msg);
            }
            errorMessages.clear();
        }
    }
    
    // Set final result stats
    result.totalReadsProcessed = state.totalReadSeedCount;
    result.totalTimeSeconds = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::high_resolution_clock::now() - startTime).count();
    
    logging::info("Successfully processed {} of {} nodes in placement traversal", 
                 nodesProcessed.load(), totalNodes);
    
    // If we didn't process all nodes due to errors, log a warning
    if (nodesProcessed.load() < totalNodes) {
        logging::warn("Placement traversal did not complete successfully for all nodes");
    }
}

// Main placement function - simplified and consolidated
void place(
    PlacementResult& result,
    panmanUtils::Tree* T,
    ::Index::Reader& index,
    const std::string& reads1Path,
    const std::string& reads2Path,
    std::string& placementFileName) {
    
    // Initialize node tracking log
    logging::initNodeTracking();
    
    // Initialize progress state
    progress_state = std::make_shared<PlacementProgressState>();
    progress_state->startTime = std::chrono::high_resolution_clock::now();
    progress_state->running = true;
    
    // Set up parameters
    TraversalParams params;
    params.k = index.getK();
    params.s = index.getS();
    params.t = index.getT();
    params.open = index.getOpen();
    params.scoreScale = 0.5; // Default to equal weighting
    
    logging::info("Placement parameters: k={}, s={}, t={}, open={}, scoreScale={}",
                params.k, params.s, params.t, params.open ? "true" : "false", params.scoreScale);
    
    if (!T || !T->root) {
        throw std::runtime_error("Invalid tree or root node in place() function");
    }
    
    // Initialize state manager using the lighter, optimized version
    auto stateManager = indexing::initializeStateManagerLight(T, T->root, params.k, params.s);
    
    // Load k-mer dictionary from the index if available
    std::unordered_map<uint32_t, std::string> idToKmer;
    if (index.hasKmerDictionary()) {
        auto kmerDictionary = index.getKmerDictionary();
        for (uint32_t i = 0; i < kmerDictionary.size(); i++) {
            auto entry = kmerDictionary[i];
            std::string kmerSeq = entry.getSequence();
            idToKmer[i] = kmerSeq;
        }
        std::cout << "Loaded " << idToKmer.size() << " k-mers from dictionary" << std::endl;
    }
    
    // Initialize global state
    PlacementGlobalState state;
    state.perNodeSeedMutations = index.getPerNodeSeedMutations();
    state.perNodeGapMutations = index.getPerNodeGapMutations();
    state.nodePathInfo = index.getNodePathInfo();
    state.blockInfo = index.getBlockInfo();
    state.ancestorMatrix = index.getAncestorMatrix();
    
    // Add the k-mer dictionary to the state
    state.kmerDictionary = idToKmer;
    
    // Process reads to get seed counts
    if (reads1Path.empty() && reads2Path.empty()) {
        logging::warn("No read files provided, just loading the index. Exiting.");
        return;
    }
    
    logging::info("Processing reads from {}{}", 
                std::string(reads1Path),
                reads2Path.empty() ? "" : " and " + std::string(reads2Path));
    
    auto readStart = std::chrono::high_resolution_clock::now();
    
    // Process reads and extract seeds
    std::unordered_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::vector<std::vector<seeding::seed_t>> readSeeds;
    
    try {
        // Use the seeding function to extract seeds from FASTQ files
        seeding::seedsFromFastq(
            params.k, params.s, params.t, params.open, 0,  // 0 for l parameter
            readSeedCounts, readSequences, readQuals, readNames, readSeeds,
            std::string(reads1Path), std::string(reads2Path)
        );
        
        auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - readStart).count();
            
        logging::info("Processed {} reads in {}ms, found {} unique seeds", 
            readSequences.size(), readDuration, readSeedCounts.size());
        
        // Convert read seed counts to our required format
        state.seedFreqInReads.clear();
        state.totalReadSeedCount = 0;
        
        for (const auto& [seedHash, countPair] : readSeedCounts) {
            // Sum forward and reverse counts
            size_t totalCount = countPair.first + countPair.second;
            state.seedFreqInReads[seedHash] = totalCount;
            state.totalReadSeedCount += totalCount;
        }
        
        // Calculate Jaccard denominator
        state.jaccardDenominator = readSeedCounts.size();
        
        logging::info("Total seed occurrences in reads: {}", state.totalReadSeedCount);
    } 
    catch (const std::exception& e) {
        throw std::runtime_error("Error processing reads: " + std::string(e.what()));
    }
    
    if (state.totalReadSeedCount == 0) {
        logging::warn("No seeds found in reads. Check read format and k/s parameters.");
        return;
    }
    
    // Run the traversal
    auto traversalStart = std::chrono::high_resolution_clock::now();
    logging::info("Starting placement traversal");
    
    try {
        placementTraversal(*stateManager, result, T, state, params);
        
        auto traversalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - traversalStart).count();
            
        logging::info("Placement traversal completed in {}ms", traversalDuration);
        
        // Write results to file if a filename is provided
        if (!placementFileName.empty()) {
            std::ofstream outFile(placementFileName);
            if (outFile.is_open()) {
                outFile << "Placement Results:\n";
                
                outFile << "\nBest hit count: " << result.maxHitsInAnyGenome;
                if (result.maxHitsNode) {
                    outFile << " in node " << result.maxHitsNode->identifier << "\n";
                    if (result.tiedMaxHitsNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedMaxHitsNodes.size() << "): ";
                        for (auto* node : result.tiedMaxHitsNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\n";
                    }
                }
                
                outFile << "\nBest Jaccard similarity: " << result.bestJaccardScore;
                if (result.bestJaccardNode) {
                    outFile << " in node " << result.bestJaccardNode->identifier << "\n";
                    if (result.tiedJaccardNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedJaccardNodes.size() << "): ";
                        for (auto* node : result.tiedJaccardNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\n";
                    }
                }
                
                outFile << "\nBest Cosine similarity: " << result.bestCosineScore;
                if (result.bestCosineNode) {
                    outFile << " in node " << result.bestCosineNode->identifier << "\n";
                    if (result.tiedCosineNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedCosineNodes.size() << "): ";
                        for (auto* node : result.tiedCosineNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\n";
                    }
                }
                
                outFile << "\nBest Weighted score: " << result.bestWeightedScore;
                if (result.bestWeightedNode) {
                    outFile << " in node " << result.bestWeightedNode->identifier << "\n";
                    if (result.tiedWeightedNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedWeightedNodes.size() << "): ";
                        for (auto* node : result.tiedWeightedNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\n";
                    }
                }
                
                outFile << "\nPerformance metrics:";
                outFile << "\nTotal reads processed: " << result.totalReadsProcessed;
                outFile << "\nTotal time: " << result.totalTimeSeconds << " seconds\n";
                
                outFile.close();
            }
        }
    }
    catch (const std::exception& e) {
        logging::err("Error during placement: {}", e.what());
        throw;
    }
    
    if (progress_state) {
        progress_state->running = false;
    }
}

// Simplified placeBatch function
void placeBatch(
    panmanUtils::Tree* T, 
    ::Index::Reader& index,
    const std::string& batchFilePath,
    std::string prefixBase, 
    std::string refFileNameBase,
    std::string samFileNameBase, 
    std::string bamFileNameBase,
    std::string mpileupFileNameBase, 
    std::string vcfFileNameBase,
    std::string aligner, 
    const std::string& refNode,
    const bool& save_jaccard, 
    const bool& show_time,
    const float& score_proportion, 
    const int& max_tied_nodes) {
    
    logging::info("Starting batch placement with file: {}", batchFilePath);
    
    if (!boost::filesystem::exists(batchFilePath)) {
        throw std::runtime_error("Batch file not found: " + batchFilePath);
    }
    
    // Read batch file line by line using C-style I/O
    FILE* batchFilePtr = fopen(batchFilePath.c_str(), "r");
    if (!batchFilePtr) {
        throw std::runtime_error("Failed to open batch file: " + batchFilePath);
    }
    
    char lineBuffer[4096];
    while (fgets(lineBuffer, sizeof(lineBuffer), batchFilePtr)) {
        // Convert to std::string for easier handling
        std::string line(lineBuffer);
        
        // Remove trailing newline if present
        if (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
            line.pop_back();
        }
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Parse line - format is: sample_name,reads1[,reads2]
        std::vector<std::string> parts;
        boost::split(parts, line, boost::is_any_of(","));
        
        if (parts.size() < 2) {
            logging::warn("Invalid batch entry: {}", line);
            continue;
        }
        
        std::string sampleName = parts[0];
        std::string reads1Path = parts[1];
        std::string reads2Path = parts.size() > 2 ? parts[2] : "";
        
        logging::info("Processing sample: {}", sampleName);
        
        // Create result object for this sample
        PlacementResult result;
        
        // Build output file path
        std::string outputFile = std::string(prefixBase) + "_" + sampleName + ".placement";
        
        try {
            // Call the main place function
            place(result, T, index, reads1Path, reads2Path, outputFile);
            
            // Log successful placement
            logging::info("Completed placement for sample: {}", sampleName);
            logging::info("Best hit node: {}", 
                        result.bestWeightedNode ? result.bestWeightedNode->identifier : "None");
        }
        catch (const std::exception& e) {
            logging::err("Error placing sample {}: {}", sampleName, e.what());
        }
    }
    
    logging::info("Batch placement completed");
    
    // Close the batch file
    if (batchFilePtr) {
        fclose(batchFilePtr);
    }
}

// Update the hits score for a node
void PlacementResult::updateHitsScore(panmanUtils::Node* node, int64_t hits) {
    if (!node) return;
    
    // Update node-specific hits
    hitsInThisGenome = hits;
    
    // Update global max hits tracking
    if (hits > maxHitsInAnyGenome) {
        // New maximum - clear previous tied nodes and set new max
        maxHitsInAnyGenome = hits;
        maxHitsNode = node;
        tiedMaxHitsNodes.clear();
        tiedMaxHitsNodes.push_back(node);
    } else if (hits == maxHitsInAnyGenome) {
        // Tied for maximum - add to tied nodes list
        tiedMaxHitsNodes.push_back(node);
    }
}

// Update the Jaccard score for a node
void PlacementResult::updateJaccardScore(panmanUtils::Node* node, double score) {
    // Validate input
    if (!node) return;
    if (score < 0.0 || score > 1.0) {
        logging::warn("Invalid Jaccard score {} (must be 0-1), ignoring", score);
        return;
    }
    
    // Update global max Jaccard tracking
    if (score > bestJaccardScore) {
        // New maximum - clear previous tied nodes and set new max
        bestJaccardScore = score;
        bestJaccardNode = node;
        tiedJaccardNodes.clear();
        tiedJaccardNodes.push_back(node);
    } else if (std::abs(score - bestJaccardScore) < 1e-9) {
        // Tied for maximum (using small epsilon for floating point comparison)
        tiedJaccardNodes.push_back(node);
    }
}

// Update the cosine score for a node
void PlacementResult::updateCosineScore(panmanUtils::Node* node, double score) {
    // Validate input
    if (!node) return;
    if (score < 0.0 || score > 1.0) {
        logging::warn("Invalid cosine score {} (must be 0-1), ignoring", score);
        return;
    }
    
    // Update global max cosine tracking
    if (score > bestCosineScore) {
        // New maximum - clear previous tied nodes and set new max
        bestCosineScore = score;
        bestCosineNode = node;
        tiedCosineNodes.clear();
        tiedCosineNodes.push_back(node);
    } else if (std::abs(score - bestCosineScore) < 1e-9) {
        // Tied for maximum (using small epsilon for floating point comparison)
        tiedCosineNodes.push_back(node);
    }
}

// Update the weighted score for a node
void PlacementResult::updateWeightedScore(panmanUtils::Node* node, double score, double scale) {
    // Validate input
    if (!node) return;
    if (score < 0.0 || score > 1.0) {
        logging::warn("Invalid weighted score {} (must be 0-1), ignoring", score);
        return;
    }
    
    // Update global max weighted score tracking
    if (score > bestWeightedScore) {
        // New maximum - clear previous tied nodes and set new max
        bestWeightedScore = score;
        bestWeightedNode = node;
        tiedWeightedNodes.clear();
        tiedWeightedNodes.push_back(node);
    } else if (std::abs(score - bestWeightedScore) < 1e-9) {
        // Tied for maximum (using small epsilon for floating point comparison)
        tiedWeightedNodes.push_back(node);
    }
}

} // namespace placement


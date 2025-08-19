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

#include "progress_state.hpp"
#include "indexing.hpp"
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>

using namespace coordinates;
using namespace logging;
using namespace state;

namespace placement {

using panmanUtils::Node;
using panmanUtils::Tree;

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

// Helper function to process reads from a FASTQ file
void processReadsFromFastq(
    const std::string& fastqPath,
    int k,
    int s,
    PlacementGlobalState& state,
    std::vector<std::string>& readSequences);

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

// Helper function for calculating cosine similarity delta
std::pair<double, double> getCosineDelta(
    bool isRemoval, 
    bool isAddition,
    size_t seedHash, 
    const absl::flat_hash_map<size_t, int64_t>& readSeedCounts,
    const absl::flat_hash_map<size_t, int64_t>& genomeSeedCounts) {

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

// Helper functions for processNodeMutations

// Process seed deletion (quaternary value 1)
void processSeedDeletion(
    panmanUtils::Node* node,
    int64_t pos,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int k) {
    
    // Get identifier, using a fallback only if node is null
    std::string nodeId = node ? node->identifier : (result.bestWeightedNode ? result.bestWeightedNode->identifier : "Unknown");
    
    // Check if this seed exists in our map
    if (node && result.nodeSeedMap.find(nodeId) != result.nodeSeedMap.end() && 
        result.nodeSeedMap[nodeId].find(pos) != result.nodeSeedMap[nodeId].end()) {
        
        // Get the seed from our map
        const seed_t& seed = result.nodeSeedMap[nodeId][pos];
        const size_t seedHash = seed.hash;
        
        // Remove it from our map
        result.nodeSeedMap[nodeId].erase(pos);
        
        // Update scores if this seed was in reads
        auto readIt = state.seedFreqInReads.find(seedHash);
        if (readIt != state.seedFreqInReads.end()) {
            const int64_t readCount = readIt->second;
            result.hitsInThisGenome -= readCount;
            result.currentJaccardNumerator -= readCount;
            
            
            auto genomeSeedIt = result.currentGenomeSeedCounts.find(seedHash);
            if (genomeSeedIt != result.currentGenomeSeedCounts.end()) {
                if (--(genomeSeedIt->second) <= 0) { 
                    result.currentGenomeSeedCounts.erase(genomeSeedIt);
                }
            }
            
            auto [numeratorDelta, denomDelta] = getCosineDelta(
                true, false, seedHash, 
                state.seedFreqInReads, 
                result.currentGenomeSeedCounts);
            
            result.currentCosineNumerator += numeratorDelta;
            result.currentCosineDenominator += denomDelta;
            
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
    panmanUtils::Node* node,
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
    std::optional<seed_t> existingSeedOpt;
    
    // Check if this seed exists in our map
    if (result.nodeSeedMap.find(nodeId) != result.nodeSeedMap.end() && 
        result.nodeSeedMap[nodeId].find(pos) != result.nodeSeedMap[nodeId].end()) {
        // Found the seed in our map
        existingSeedOpt = result.nodeSeedMap[nodeId][pos];
        
        // Process the found seed
        if (existingSeedOpt) {
            const seed_t& existingSeed = existingSeedOpt.value();
            const size_t seedHash = existingSeed.hash;
            
            // Remove it from our map
            result.nodeSeedMap[nodeId].erase(pos);
            
            // Update scores if this seed was in reads
            auto readIt = state.seedFreqInReads.find(seedHash);
            if (readIt != state.seedFreqInReads.end()) {
                int64_t readCount = readIt->second;
                result.hitsInThisGenome -= readCount;
                result.currentJaccardNumerator -= readCount;
                
                auto genomeSeedIt = result.currentGenomeSeedCounts.find(seedHash);
                if (genomeSeedIt != result.currentGenomeSeedCounts.end()) {
                    if (--(genomeSeedIt->second) <= 0) {
                        result.currentGenomeSeedCounts.erase(genomeSeedIt);
                    }
                }
                
                auto [numeratorDelta, denomDelta] = getCosineDelta(
                    true, false, seedHash,
                    state.seedFreqInReads, 
                    result.currentGenomeSeedCounts);
                    
                result.currentCosineNumerator += numeratorDelta;
                result.currentCosineDenominator += denomDelta;
                
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
    panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    const seed_t& newSeed) {
    

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
        it->second++; // Increment the count
        
        // Update cosine similarity
        auto [numeratorDelta, denomDelta] = getCosineDelta(
            false, true, newSeed.hash,
            state.seedFreqInReads, 
            result.currentGenomeSeedCounts);
        
        result.currentCosineNumerator += numeratorDelta;
        result.currentCosineDenominator += denomDelta;
        
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
 * @param node The node being evaluated
 * @param jaccardScore The Jaccard similarity score
 */
void PlacementResult::updateJaccardScore(panmanUtils::Node* node, double jaccardScore) {
    if (!node) return;
    
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (jaccardScore > bestJaccardScore + TIED_THRESHOLD) {
        // Found a new best Jaccard score
        bestJaccardScore = jaccardScore;
        bestJaccardNode = node;
        
        // Reset tied nodes and add this one
        tiedJaccardNodes.clear();
        tiedJaccardNodes.push_back(node);
        
        logging::debug("New best Jaccard node: {} with score {:.6f}", 
                     node->identifier, jaccardScore);
                     
    } else if (std::abs(jaccardScore - bestJaccardScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedJaccardNodes.push_back(node);
        logging::debug("Tied Jaccard node: {} with score {:.6f}", 
                     node->identifier, jaccardScore);
    }
}

/**
 * @brief Update Raw Seed Match score for a node and track the best score
 * 
 * @param node The node being evaluated
 * @param score The raw seed match score (sum of read frequencies for matched seeds)
 */
void PlacementResult::updateRawSeedMatchScore(panmanUtils::Node* node, int64_t score) {
    if (!node) return;
    
    // Check if score is better than current best
    if (score > bestRawSeedMatchScore) {
        // Found a new best raw seed match score
        bestRawSeedMatchScore = score;
        bestRawSeedMatchNode = node;
        
        // Reset tied nodes and add this one
        tiedRawSeedMatchNodes.clear();
        tiedRawSeedMatchNodes.push_back(node);
        
        logging::debug("New best Raw Seed Match node: {} with score {}", 
                     node->identifier, score);
                     
    } else if (score == bestRawSeedMatchScore) {
        // Add to tied nodes
        tiedRawSeedMatchNodes.push_back(node);
        logging::debug("Tied Raw Seed Match node: {} with score {}", 
                     node->identifier, score);
    }
}

/**
 * @brief Update hits score for a node and track the best score
 * 
 * @param node The node being evaluated
 * @param hits The number of hits for this node
 */
void PlacementResult::updateHitsScore(panmanUtils::Node* node, int64_t hits) {
    if (!node) return;
    
    // Check if hits is better than current best
    if (hits > maxHitsInAnyGenome) {
        // Found a new best hits score
        maxHitsInAnyGenome = hits;
        maxHitsNode = node;
        
        // Reset tied nodes and add this one
        tiedMaxHitsNodes.clear();
        tiedMaxHitsNodes.push_back(node);
        
        logging::debug("New best hits node: {} with hits {}", 
                     node->identifier, hits);
                     
    } else if (hits == maxHitsInAnyGenome) {
        // Add to tied nodes
        tiedMaxHitsNodes.push_back(node);
        logging::debug("Tied hits node: {} with hits {}", 
                     node->identifier, hits);
    }
}

/**
 * @brief Update Jaccard (Presence/Absence) score for a node and track the best score
 * 
 * @param node The node being evaluated
 * @param score The Jaccard score based on presence/absence of seeds
 */
void PlacementResult::updateJaccardPresenceScore(panmanUtils::Node* node, double score) {
    if (!node) return;
    
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (score > bestJaccardPresenceScore + TIED_THRESHOLD) {
        // Found a new best Jaccard (Presence/Absence) score
        bestJaccardPresenceScore = score;
        bestJaccardPresenceNode = node;
        
        // Reset tied nodes and add this one
        tiedJaccardPresenceNodes.clear();
        tiedJaccardPresenceNodes.push_back(node);
        
        logging::debug("New best Jaccard (Presence) node: {} with score {:.6f}", 
                     node->identifier, score);
                     
    } else if (std::abs(score - bestJaccardPresenceScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedJaccardPresenceNodes.push_back(node);
        logging::debug("Tied Jaccard (Presence) node: {} with score {:.6f}", 
                     node->identifier, score);
    }
}

/**
 * @brief Update cosine similarity score for a node and track the best score
 * 
 * @param node The node being evaluated
 * @param cosineScore The cosine similarity score
 */
void PlacementResult::updateCosineScore(panmanUtils::Node* node, double cosineScore) {
    if (!node) return;
    
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (cosineScore > bestCosineScore + TIED_THRESHOLD) {
        // Found a new best cosine score
        bestCosineScore = cosineScore;
        bestCosineNode = node;
        
        // Reset tied nodes and add this one
        tiedCosineNodes.clear();
        tiedCosineNodes.push_back(node);
        
        logging::debug("New best cosine node: {} with score {:.6f}", 
                     node->identifier, cosineScore);
                     
    } else if (std::abs(cosineScore - bestCosineScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedCosineNodes.push_back(node);
        logging::debug("Tied cosine node: {} with score {:.6f}", 
                     node->identifier, cosineScore);
    }
}

/**
 * @brief Update weighted score (combined Jaccard and cosine) for a node
 * 
 * @param node The node being evaluated
 * @param weightedScore The weighted similarity score
 * @param scoreScale The weight for Jaccard in the combined score (1-scoreScale is used for cosine)
 */
void PlacementResult::updateWeightedScore(panmanUtils::Node* node, double weightedScore, double scoreScale) {
    if (!node) return;
    
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (weightedScore > bestWeightedScore + TIED_THRESHOLD) {
        // Found a new best weighted score
        bestWeightedScore = weightedScore;
        bestWeightedNode = node;
        
        // Reset tied nodes and add this one
        tiedWeightedNodes.clear();
        tiedWeightedNodes.push_back(node);
        
        logging::debug("New best weighted node: {} with score {:.6f} (scale={:.2f})", 
                     node->identifier, weightedScore, scoreScale);
                     
    } else if (std::abs(weightedScore - bestWeightedScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedWeightedNodes.push_back(node);
        logging::debug("Tied weighted node: {} with score {:.6f}", 
                     node->identifier, weightedScore);
    }
}


// Main placement function - simplified and consolidated
void place(
    PlacementResult& result,
    panmanUtils::Tree* T,
    ::MGSRIndex::Reader& mgsrIndex,
    const std::string& reads1Path,
    const std::string& reads2Path,
    std::vector<std::vector<seed_t>>& readSeeds,
    std::vector<std::string>& readSequences,
    std::vector<std::string>& readNames,
    std::vector<std::string>& readQuals,
    std::string& placementFileName,
    const std::string& indexPath,
    const std::string& debug_node_id_param) {  
    
    // FIXED: Store the message reader to manage its lifetime
    // This ensures the underlying data stays alive throughout the function
    std::shared_ptr<::capnp::MessageReader> messageHolder;
    
    // Only use the explicit indexPath if it\'s provided and not empty
    if (!indexPath.empty()) {
        logging::debug("Loading MGSR index directly from file: {}", indexPath);
        
        try {
            // Load and keep the message reader alive
            messageHolder = std::shared_ptr<::capnp::MessageReader>(loadMgsrIndexFromFile(indexPath).release());
            
            // Now safely get a reference from the message holder
            mgsrIndex = messageHolder->getRoot<MGSRIndex>();
            logging::debug("Successfully loaded index from: {}", indexPath);
        } catch (const std::exception& e) {
            logging::err("Failed to load index from file: {}", e.what());
            throw std::runtime_error("Failed to load index from file: " + std::string(e.what()));
        }
    }
    
    try {
        // Verify the MGSR index message has valid data
        uint32_t k = mgsrIndex.getK();
        uint32_t s = mgsrIndex.getS();
        size_t seedCount = mgsrIndex.getSeedInfo().size();
        bool useRawSeeds = mgsrIndex.getUseRawSeeds();
        
        if (k == 0 || s == 0 || seedCount == 0) {
            throw std::runtime_error(
                "Invalid MGSR index: k=" + std::to_string(k) + 
                ", s=" + std::to_string(s) + 
                ", seeds=" + std::to_string(seedCount) + 
                ". The index appears to be corrupted or incomplete.");
        }
        
        logging::debug("Verified MGSR index has valid root pointer with k={}, s={}, seeds={}, useRawSeeds={}", 
                     k, s, seedCount, useRawSeeds);
    } catch (const ::kj::Exception& e) {
        throw std::runtime_error("Cap\'n Proto error: Message did not contain a valid root pointer. The index file appears to be corrupted or not properly initialized.");
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error validating MGSR index: ") + e.what());
    }
    
    // Initialize node tracking log
    logging::initNodeTracking();
    
    // Initialize progress state
    progress_state = std::make_shared<PlacementProgressState>();
    progress_state->startTime = std::chrono::high_resolution_clock::now();
    progress_state->running = true;
    
    // Set up parameters from MGSR index
    TraversalParams params;
    params.k = mgsrIndex.getK();
    params.s = mgsrIndex.getS();
    params.t = mgsrIndex.getT();
    params.open = mgsrIndex.getL(); // MGSR uses L parameter instead of open
    params.scoreScale = 0.5; // Default to equal weighting
    params.debug_node_id = debug_node_id_param; // <-- SET THE PARAM IN STRUCT
    
    logging::debug("Placement parameters: k={}, s={}, t={}, open={}, scoreScale={}, debug_node_id='{}'",
                params.k, params.s, params.t, params.open ? "true" : "false", params.scoreScale, params.debug_node_id);
    
    if (!T || !T->root) {
        throw std::runtime_error("Invalid tree or root node in place() function");
    }
    
    // Initialize state manager using the lighter, optimized version
    logging::info("Initializing StateManager (light) with k={}, s={}", params.k, params.s);
    auto stateManager = indexing::initializeStateManagerLight(T, T->root, params.k, params.s);

    // Initialize simplified MGSR placement state (no legacy dictionary)
    PlacementGlobalState state;
    state.kmerSize = params.k;
    auto seedInfos = mgsrIndex.getSeedInfo();
    auto perNodeChanges = mgsrIndex.getPerNodeChanges();
    logging::info("Initializing MGSR placement: {} seeds; {} node change records", seedInfos.size(), perNodeChanges.size());
    
    // DEBUG: Detailed index information
    logging::info("=== MGSR INDEX DETAILS ===");
    logging::info("Index parameters: k={}, s={}, t={}, l={}, open={}, useRawSeeds={}", 
                  mgsrIndex.getK(), mgsrIndex.getS(), mgsrIndex.getT(), 
                  mgsrIndex.getL(), mgsrIndex.getOpen(), mgsrIndex.getUseRawSeeds());
    
    // Sample first few seeds from index
    logging::info("Sample seeds from index (first 10):");
    for (size_t i = 0; i < std::min(static_cast<size_t>(10), static_cast<size_t>(seedInfos.size())); i++) {
        auto seed = seedInfos[i];
        logging::info("  Seed[{}]: hash={}, startPos={}, endPos={}, isReverse={}", 
                      i, seed.getHash(), seed.getStartPos(), seed.getEndPos(), seed.getIsReverse());
    }
    
    // Sample first few node changes
    // Minimal StateManager usage retained for potential future extension
    stateManager->initializeSeedStorage();
    
    // Pre-initialize nodes in breadth-first order to ensure proper parent-child relationship
    if (T->root) {
        // First, ensure the root node is initialized
        stateManager->initializeNode(T->root->identifier);
        
        // Then initialize the rest level by level (breadth-first)
        std::vector<panmanUtils::Node*> currentLevel = {T->root};
        std::vector<panmanUtils::Node*> nextLevel;
        
        while (!currentLevel.empty()) {
            for (auto* node : currentLevel) {
                // Add all children to the next level
                for (auto* child : node->children) {
                    if (child) {
                        nextLevel.push_back(child);
                        
                        // Initialize this child node
                        try {
                            stateManager->initializeNode(child->identifier);
                        } catch (const std::exception& e) {
                            logging::debug("Initialization for node {} failed: {}", child->identifier, e.what());
                        }
                    }
                }
            }
            
            // Move to next level
            currentLevel = std::move(nextLevel);
            nextLevel.clear();
        }
    } else {
        throw std::runtime_error("Invalid tree or root node in place() function");
    }
    
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
    absl::flat_hash_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
    
    try {
        std::vector<std::vector<std::string>> readSeedSeqs;
        // Use the seeding function to extract seeds from FASTQ files
        seeding::seedsFromFastq(
            params.k, params.s, params.t, false, 0,  // 0 for l parameter
            readSeedCounts, readSequences, readQuals, readNames, readSeeds,
            readSeedSeqs,
            std::string(reads1Path), std::string(reads2Path)
        );
        
        auto readDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - readStart).count();
            
        logging::info("Processed {} reads in {}ms, found {} unique seed hashes", 
            readSequences.size(), readDuration, readSeedCounts.size());
        
        // Show top 10 most frequent seed hashes
        std::vector<std::pair<size_t, std::pair<size_t, size_t>>> seedFreqVec(readSeedCounts.begin(), readSeedCounts.end());
        std::sort(seedFreqVec.begin(), seedFreqVec.end(), 
                 [](const auto& a, const auto& b) { 
                     return (a.second.first + a.second.second) > (b.second.first + b.second.second); 
                 });
        
        logging::debug("Top 10 most frequent seed hashes in reads:");
        for (size_t i = 0; i < std::min(static_cast<size_t>(10), seedFreqVec.size()); i++) {
            auto hash = seedFreqVec[i].first;
            auto forwardCount = seedFreqVec[i].second.first;
            auto reverseCount = seedFreqVec[i].second.second;
            logging::debug("  Hash[{}]: forward={}, reverse={}, total={}", 
                          hash, forwardCount, reverseCount, forwardCount + reverseCount);
        }
        
        // Convert read seed counts to our required format
        state.seedFreqInReads.clear();
        state.totalReadSeedCount = 0;
        
        logging::debug("=== SEED FREQUENCY CONVERSION ===");
        for (const auto& [seedHash, countPair] : readSeedCounts) {
            // Sum forward and reverse counts
            size_t totalCount = countPair.first + countPair.second;
            state.seedFreqInReads[seedHash] = totalCount;
            state.totalReadSeedCount += totalCount;
            
            // Debug first few conversions
            if (state.seedFreqInReads.size() <= 5) {
                logging::debug("  Converting hash[{}]: forward={}, reverse={}, total={}", 
                              seedHash, countPair.first, countPair.second, totalCount);
            }
        }
        
        // Calculate Jaccard denominator
        state.jaccardDenominator = readSeedCounts.size();
        state.readUniqueSeedCount = readSeedCounts.size();
        
        logging::info("Seed processing complete: {} unique seeds, {} total occurrences", 
                      state.seedFreqInReads.size(), state.totalReadSeedCount);
        
        // Populate hashToKmer map from read seeds for debugging
        for (int i = 0; i < readSeeds.size(); i++) {
            for (int j = 0; j < readSeeds[i].size(); j++) {
                const auto& seed = readSeeds[i][j];
                const std::string& kmer_str = readSeedSeqs[i][j];
                if (state.hashToKmer.find(seed.hash) == state.hashToKmer.end()) {
                    state.hashToKmer[seed.hash] = kmer_str;
                }
            }
        }

    } 
    catch (const std::exception& e) {
        logging::err("Error processing reads: {}", e.what());
        throw std::runtime_error("Error processing reads: " + std::string(e.what()));
    }
    
    if (state.totalReadSeedCount == 0) {
        logging::warn("No seeds found in reads. Check read format and k/s parameters.");
        return;
    }
    
    // Run the traversal
    auto traversalStart = std::chrono::high_resolution_clock::now();
    logging::info("=== STARTING PLACEMENT TRAVERSAL ===");
    logging::info("Total read seeds: {}, Unique hashes: {}", state.totalReadSeedCount, state.seedFreqInReads.size());
    
    try {
        // MGSR DFS traversal applying per-node seed changes
        auto seedInfos = mgsrIndex.getSeedInfo();
        auto perNodeChanges = mgsrIndex.getPerNodeChanges();
        
        // Create position-to-seed-index mapping for deletion operations
        // In MGSR, seedDeletions contains positions, but we need seed indices for activeSeedIndices
        std::unordered_map<uint32_t, uint32_t> positionToSeedIndex;
        for (size_t i = 0; i < seedInfos.size(); i++) {
            uint32_t startPos = seedInfos[i].getStartPos();
            positionToSeedIndex[startPos] = i;
        }
        
        // Initialize reference seed state (MGSR starts with empty seed set)
        logging::info("Initialized empty reference state (MGSR starts with no seeds)");
        
        std::function<void(panmanUtils::Node*, uint64_t&)> dfs;
        dfs = [&](panmanUtils::Node* node, uint64_t& dfsIndex) {
            if (!node) return;
            uint64_t myIndex = dfsIndex;
            
            // Debug logging for first few nodes
            bool debugNode = myIndex < 5;
            
            // Apply node-specific seed changes - track operations for proper backtracking
            struct Operation {
                enum Type { ADD, REMOVE } type;
                uint32_t seedIdx;
                uint32_t originalCount;  // For REMOVE operations, store what the count was before removal
            };
            std::vector<Operation> operations;
            operations.reserve(24);
            
            if (myIndex < perNodeChanges.size()) {
                auto ch = perNodeChanges[myIndex];
                
                
                // Process seed deletions FIRST (seedDeletions contains positions, map to indices)
                for (auto pos : ch.getSeedDeletions()) { 
                    auto it = positionToSeedIndex.find(pos);
                    if (it != positionToSeedIndex.end()) {
                        uint32_t seedIdx = it->second;
                        auto activeIt = state.activeSeedIndices.find(seedIdx);
                        if (activeIt != state.activeSeedIndices.end()) { 
                            // Store the original count for backtracking
                            uint32_t originalCount = activeIt->second;
                            
                            if (activeIt->second == 1) {
                                state.activeSeedIndices.erase(activeIt); 
                            } else {
                                activeIt->second--; 
                            }
                            
                            // Record this removal operation
                            operations.push_back({Operation::REMOVE, seedIdx, originalCount});
                        }
                    }
                }
                
                // Process seed additions SECOND (seedInsubIndices contains seed indices)
                for (auto idx : ch.getSeedInsubIndices()) { 
                    // Check if this is a substitution (new seed at position where old seed exists)
                    uint32_t newSeedPos = seedInfos[idx].getStartPos();
                    auto posIt = positionToSeedIndex.find(newSeedPos);
                    
                    if (posIt != positionToSeedIndex.end()) {
                        uint32_t oldSeedIdx = posIt->second;
                        // Check if the old seed is currently active and is different from the new seed
                        auto activeIt = state.activeSeedIndices.find(oldSeedIdx);
                        if (activeIt != state.activeSeedIndices.end() && oldSeedIdx != idx) {
                            // This is a substitution - remove the old seed first
                            uint32_t originalCount = activeIt->second;
                                                        
                            if (activeIt->second == 1) {
                                state.activeSeedIndices.erase(activeIt);
                            } else {
                                activeIt->second--;
                            }
                            
                            // Record the substitution removal
                            operations.push_back({Operation::REMOVE, oldSeedIdx, originalCount});
                        }
                    }
                    
                    // Add the new seed - track original count for proper backtracking
                    uint32_t addOriginalCount = state.activeSeedIndices[idx]; // 0 if new entry
                    uint32_t newCount = ++state.activeSeedIndices[idx];
                    operations.push_back({Operation::ADD, idx, addOriginalCount});
                    
                }
            }
            
            // Calculate scores for this node
            int64_t hits = 0;
            std::set<size_t> matchingHashes; // Track unique hashes that match
            std::set<size_t> genomeHashes;   // Track all hashes in current genome state
            
            // First pass: collect all unique hashes in current genome state
            for (auto &kv : state.activeSeedIndices) {
                size_t h = seedInfos[kv.first].getHash();
                genomeHashes.insert(h);
            }
            
            // Second pass: calculate hits only once per unique hash
            for (size_t h : genomeHashes) {
                auto rit = state.seedFreqInReads.find(h);
                if (rit != state.seedFreqInReads.end()) { 
                    hits += rit->second; 
                    matchingHashes.insert(h);
                }
            }
            
            // // VERIFICATION: Calculate expected hits by direct comparison
            // int64_t expectedHits = 0;
            // std::set<size_t> expectedMatchingHashes;
            // for (const auto& [readHash, readCount] : state.seedFreqInReads) {
            //     if (genomeHashes.count(readHash)) {
            //         expectedHits += readCount;
            //         expectedMatchingHashes.insert(readHash);
            //     }
            // }
            
            // // Debug verification for selected nodes
            // logging::info("VERIFY Node {}: MGSR_hits={}, Expected_hits={}, Active_seeds={}, Matching_hashes={}/{}", 
            //                 node->identifier, hits, expectedHits, state.activeSeedIndices.size(), 
            //                 matchingHashes.size(), state.seedFreqInReads.size());
            // if (hits != expectedHits) {
            //     logging::warn("  MISMATCH: MGSR calculation differs from expected!");
            // }
        
            
            // Update result scores
            result.updateHitsScore(node, hits);
            result.updateRawSeedMatchScore(node, hits);
            
            // Calculate proper Jaccard presence score
            if (state.readUniqueSeedCount > 0) {
                // Jaccard = |intersection| / |union|
                // intersection = matchingHashes.size()
                // union = |read_hashes| + |genome_hashes| - |intersection|
                size_t intersection = matchingHashes.size();
                size_t union_size = state.readUniqueSeedCount + genomeHashes.size() - intersection;
                double jPresence = union_size > 0 ? (double)intersection / (double)union_size : 0.0;
                result.updateJaccardPresenceScore(node, jPresence);
                
                // Calculate standard Jaccard similarity (intersection over union)
                double jaccardSim = union_size > 0 ? (double)intersection / (double)union_size : 0.0;
                result.updateJaccardScore(node, jaccardSim);
                
                // Calculate cosine similarity: intersection / sqrt(|read| * |genome|)
                double cosineSim = 0.0;
                if (state.readUniqueSeedCount > 0 && genomeHashes.size() > 0) {
                    cosineSim = (double)intersection / sqrt((double)state.readUniqueSeedCount * (double)genomeHashes.size());
                }
                result.updateCosineScore(node, cosineSim);
                
                // Calculate weighted score (combine Jaccard and cosine with equal weighting)
                double weightedScore = 0.5 * jaccardSim + 0.5 * cosineSim;
                result.updateWeightedScore(node, weightedScore, 0.5);
            }
            
            // Recursively process children with correct DFS indexing
            for (auto *child : node->children) {
                dfsIndex++;
                dfs(child, dfsIndex);
            }
            
            for (auto it = operations.rbegin(); it != operations.rend(); ++it) {
                const Operation& op = *it;
                
                if (op.type == Operation::ADD) {
                    // Undo addition: restore the seed to its original count
                    if (op.originalCount == 0) {
                        // Seed didn't exist before, remove it completely
                        state.activeSeedIndices.erase(op.seedIdx);
                    } else {
                        // Seed existed before, restore original count
                        state.activeSeedIndices[op.seedIdx] = op.originalCount;
                    }
                } else { // Operation::REMOVE
                    // Undo removal: restore the seed with its original count
                    state.activeSeedIndices[op.seedIdx] = op.originalCount;
                }
            }
        };
        
        logging::info("Starting DFS traversal from root: {}", T->root->identifier);
        uint64_t startIndex = 0; dfs(T->root, startIndex);
        
        auto traversalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - traversalStart).count();
            
        logging::info("Placement traversal completed in {}ms", traversalDuration);
        
        // Log compact placement results
        logging::info("=== PLACEMENT RESULTS ===");
        
        // Raw hits score
        if (result.maxHitsNode) {
            logging::info("Raw Hits: {} (best: {})", 
                         result.maxHitsInAnyGenome,
                         result.maxHitsNode->identifier);
            if (result.tiedMaxHitsNodes.size() > 1) {
                std::string tiedNodes;
                for (size_t i = 0; i < result.tiedMaxHitsNodes.size(); i++) {
                    if (i > 0) tiedNodes += ", ";
                    tiedNodes += result.tiedMaxHitsNodes[i]->identifier;
                }
                logging::info("  Tied nodes ({}): {}", result.tiedMaxHitsNodes.size(), tiedNodes);
            }
        } else {
            logging::info("Raw Hits: {} (no node found)", result.maxHitsInAnyGenome);
        }
        
        // Jaccard presence score  
        if (result.bestJaccardPresenceNode) {
            logging::info("Jaccard Presence: {:.6f} (best: {})",
                         result.bestJaccardPresenceScore,
                         result.bestJaccardPresenceNode->identifier);
            if (result.tiedJaccardPresenceNodes.size() > 1) {
                std::string tiedNodes;
                for (size_t i = 0; i < result.tiedJaccardPresenceNodes.size(); i++) {
                    if (i > 0) tiedNodes += ", ";
                    tiedNodes += result.tiedJaccardPresenceNodes[i]->identifier;
                }
                logging::info("  Tied nodes ({}): {}", result.tiedJaccardPresenceNodes.size(), tiedNodes);
            }
        } else {
            logging::info("Jaccard Presence: {:.6f} (no node found)", result.bestJaccardPresenceScore);
        }
        
        // Jaccard similarity score
        if (result.bestJaccardNode) {
            logging::info("Jaccard Similarity: {:.6f} (best: {})",
                         result.bestJaccardScore,
                         result.bestJaccardNode->identifier);
            if (result.tiedJaccardNodes.size() > 1) {
                std::string tiedNodes;
                for (size_t i = 0; i < result.tiedJaccardNodes.size(); i++) {
                    if (i > 0) tiedNodes += ", ";
                    tiedNodes += result.tiedJaccardNodes[i]->identifier;
                }
                logging::info("  Tied nodes ({}): {}", result.tiedJaccardNodes.size(), tiedNodes);
            }
        } else {
            logging::info("Jaccard Similarity: {:.6f} (no node found)", result.bestJaccardScore);
        }
        
        // Cosine similarity score
        if (result.bestCosineNode) {
            logging::info("Cosine Similarity: {:.6f} (best: {})",
                         result.bestCosineScore,
                         result.bestCosineNode->identifier);
            if (result.tiedCosineNodes.size() > 1) {
                std::string tiedNodes;
                for (size_t i = 0; i < result.tiedCosineNodes.size(); i++) {
                    if (i > 0) tiedNodes += ", ";
                    tiedNodes += result.tiedCosineNodes[i]->identifier;
                }
                logging::info("  Tied nodes ({}): {}", result.tiedCosineNodes.size(), tiedNodes);
            }
        } else {
            logging::info("Cosine Similarity: {:.6f} (no node found)", result.bestCosineScore);
        }
        
        // Weighted score
        if (result.bestWeightedNode) {
            logging::info("Weighted Score: {:.6f} (best: {})",
                         result.bestWeightedScore,
                         result.bestWeightedNode->identifier);
            if (result.tiedWeightedNodes.size() > 1) {
                std::string tiedNodes;
                for (size_t i = 0; i < result.tiedWeightedNodes.size(); i++) {
                    if (i > 0) tiedNodes += ", ";
                    tiedNodes += result.tiedWeightedNodes[i]->identifier;
                }
                logging::info("  Tied nodes ({}): {}", result.tiedWeightedNodes.size(), tiedNodes);
            }
        } else {
            logging::info("Weighted Score: {:.6f} (no node found)", result.bestWeightedScore);
        }
        
        // Summary stats
        logging::info("Total unique read seeds: {}, Total seed occurrences: {}", 
                     state.readUniqueSeedCount,
                     state.totalReadSeedCount);
        
        // Write results to file if a filename is provided
        if (!placementFileName.empty()) {
            std::ofstream outFile(placementFileName);
            if (outFile.is_open()) {
                outFile << "Placement Results:\\n";
                
                outFile << "\\nBest hit count: " << result.maxHitsInAnyGenome;
                if (result.maxHitsNode) {
                    outFile << " in node " << result.maxHitsNode->identifier << "\\n";
                    if (result.tiedMaxHitsNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedMaxHitsNodes.size() << "): ";
                        for (auto* node : result.tiedMaxHitsNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\\n";
                    }
                }
                
                outFile << "\\nBest Jaccard similarity: " << result.bestJaccardScore;
                if (result.bestJaccardNode) {
                    outFile << " in node " << result.bestJaccardNode->identifier << "\\n";
                    if (result.tiedJaccardNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedJaccardNodes.size() << "): ";
                        for (auto* node : result.tiedJaccardNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\\n";
                    }
                }
                
                outFile << "\\nBest Cosine similarity: " << result.bestCosineScore;
                if (result.bestCosineNode) {
                    outFile << " in node " << result.bestCosineNode->identifier << "\\n";
                    if (result.tiedCosineNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedCosineNodes.size() << "): ";
                        for (auto* node : result.tiedCosineNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\\n";
                    }
                }
                
                outFile << "\\nBest Weighted score: " << result.bestWeightedScore;
                if (result.bestWeightedNode) {
                    outFile << " in node " << result.bestWeightedNode->identifier << "\\n";
                    if (result.tiedWeightedNodes.size() > 1) {
                        outFile << "Tied nodes (" << result.tiedWeightedNodes.size() << "): ";
                        for (auto* node : result.tiedWeightedNodes) {
                            outFile << node->identifier << " ";
                        }
                        outFile << "\\n";
                    }
                }
                
                outFile << "\\nPerformance metrics:";
                outFile << "\\nTotal reads processed: " << result.totalReadsProcessed;
                outFile << "\\nTotal time: " << result.totalTimeSeconds << " seconds\\n";
                
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

    // Log read seed hashes for debugging
    std::ofstream read_seeds_log_file("read_seeds_debug.log");
    if (read_seeds_log_file.is_open()) {
        read_seeds_log_file << "# Read Seed Hashes and Frequencies (Total Unique: " << state.seedFreqInReads.size() << ", Total Occurrences: " << state.totalReadSeedCount << ")\n";
        read_seeds_log_file << "# Hash\tFrequency\tKmer (if found in state.hashToKmer)\n";
        for (const auto& entry : state.seedFreqInReads) {
            std::string kmer_str = "<kmer_not_in_read_fallback_map>";
            if (state.hashToKmer.count(entry.first)) {
                kmer_str = state.hashToKmer.at(entry.first);
            }
            read_seeds_log_file << entry.first << "\t" << entry.second << "\t" << kmer_str << "\n"; 
        }
        read_seeds_log_file.close(); 
        logging::info("Dumped read seed frequencies to read_seeds_debug.log");
    } else {
        logging::warn("Could not open read_seeds_debug.log for writing.");
    }

    if (state.totalReadSeedCount == 0) {
        logging::warn("No seeds found in reads. Check read format and k/s parameters.");
        return;
    }
} // Closing brace for void place(...)

// Simplified placeBatch function
void placeBatch(
    panmanUtils::Tree* T, 
    ::MGSRIndex::Reader& mgsrIndex,
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
    const int& max_tied_nodes,
    const std::string& indexPath,
    const std::string& debug_node_id_param) {  // Add indexPath parameter
    
    logging::debug("Starting batch placement with file: {}", batchFilePath);
    
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
        std::string line_str(lineBuffer); // Renamed to avoid conflict with any other 'line' identifier
        
        // Remove trailing newline if present
        if (!line_str.empty() && (line_str.back() == '\n' || line_str.back() == '\r')) { // Corrected single quotes for char literals
            line_str.pop_back();
        }
        
        // Skip empty lines and comments
        if (line_str.empty() || line_str[0] == '#') continue; // Corrected single quotes for char literal
        
        // Parse line - format is: sample_name,reads1[,reads2]
        std::vector<std::string> parts;
        boost::split(parts, line_str, boost::is_any_of(",")); // Ensured comma is part of string literal
        
        if (parts.size() < 2) {
            logging::warn("Invalid batch entry: {}", line_str);
            continue;
        }
        
        std::string sampleName = parts[0];
        std::string reads1Path = parts[1];
        std::string reads2Path = parts.size() > 2 ? parts[2] : "";
        
        logging::debug("Processing sample: {}", sampleName);
        
        // Create result object for this sample
        PlacementResult result;
        
        // Build output file path
        std::string outputFile = std::string(prefixBase) + "_" + sampleName + ".placement";
        
        std::vector<std::vector<seed_t>> readSeeds;
        std::vector<std::string> readSequences;
        std::vector<std::string> readNames;
        std::vector<std::string> readQuals;

        try {
            // Call the main place function with correct parameters
            place(result, T, mgsrIndex, reads1Path, reads2Path, readSeeds, readSequences, readNames, readQuals, outputFile, indexPath, debug_node_id_param);
            
            // Log successful placement
            logging::debug("Completed placement for sample: {}", sampleName);
            logging::debug("Best hit node: {}", 
                        result.bestWeightedNode ? result.bestWeightedNode->identifier : "None");
        }
        catch (const std::exception& e) {
            logging::err("Error placing sample {}: {}", sampleName, e.what());
        }
    }
    
    logging::debug("Batch placement completed");
    
    // Close the batch file
    if (batchFilePtr) {
        fclose(batchFilePtr);
    }
}


} // namespace placement


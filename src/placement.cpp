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

// Function declarations
void dumpPlacementDebugData(
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    const std::vector<std::string>& nodesToDump,
    size_t maxNodes,
    const TraversalParams& params,
    const std::string& outputFilename,
    const std::vector<std::string>& readSequences);

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

// Helper function to read a packed index file
std::unique_ptr<::capnp::MessageReader> loadIndexFromFile(const std::string& indexPath) {
    logging::debug("Loading index from file: {}", indexPath);
    
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
        
        // Use PackedFdMessageReader directly instead of the wrapper
        auto reader = std::make_unique<::capnp::PackedFdMessageReader>(fd, opts);
        
        // Validate the index immediately to ensure it's usable
        auto indexRoot = reader->getRoot<Index>();
        uint32_t k = indexRoot.getK();
        uint32_t s = indexRoot.getS();
        size_t nodeCount = indexRoot.getPerNodeSeedMutations().size();
        
        if (k <= 0 || nodeCount <= 0) {
            close(fd); // Close the fd since we're not returning the reader
            throw std::runtime_error(
                "Invalid index data: k=" + std::to_string(k) + 
                ", nodes=" + std::to_string(nodeCount) + 
                ". The index appears to be corrupted or incomplete.");
        }
        
        logging::debug("Successfully loaded packed index with k={}, s={}, nodes={}",
                    k, s, nodeCount);
        
        // Return the MessageReader
        return reader;
    } catch (const ::kj::Exception& e) {
        // Clean up fd on exception
        close(fd);
        throw std::runtime_error("Cap'n Proto error: " + std::string(e.getDescription().cStr()));
    } catch (const std::exception& e) {
        // Clean up fd on exception
        close(fd);
        throw std::runtime_error("Error reading index: " + std::string(e.what()));
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

// Create seed from dictionary lookup
seeding::seed_t createSeedFromDictionary(
    const std::string& kmerStr,
    int64_t pos,
    int64_t endPos,
    int params_k,
    int params_s) {
    
    // Normalize k-mer to uppercase for consistent hashing
    std::string upperKmer = kmerStr;
    std::transform(upperKmer.begin(), upperKmer.end(), upperKmer.begin(),
                  [](unsigned char c){ return std::toupper(c); });
    
    auto syncmerResults = seeding::rollingSyncmers(upperKmer, params_k, params_s, false, 0, false);
    
    // Create default seed
    seeding::seed_t seed;
    
    // Only if a valid syncmer was found
    if (!syncmerResults.empty()) {
        auto [hash, isReversed, isSyncmer, startPos] = syncmerResults[0];
        
        // Only use the result if it's actually a syncmer
        if (isSyncmer) {
            // Fill in the seed fields
            seed.hash = hash;
            seed.reversed = isReversed;
            seed.endPos = endPos;
            
            // Return immediately when we have a valid seed
            return seed;
        }
    }
    
    // If we reach here, no valid syncmer was found, return default seed
    seed.hash = 0;
    seed.reversed = false;
    seed.endPos = 0;
    return seed;
}

// Add new seed and update scores
void addSeedAndUpdateScores(
    panmanUtils::Node* node,
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
void processNewSeedAddition(
    panmanUtils::Node* node,
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    PlacementResult& result,
    absl::flat_hash_set<size_t>& uniqueSeedHashes,
    int64_t pos,
    const TraversalParams& params,
    absl::flat_hash_map<int64_t, uint32_t>& positionToDictId,
    absl::flat_hash_map<int64_t, uint32_t>& positionToEndOffset,
    size_t& seedAdditions,
    size_t& dictionaryLookups) {
    
    bool dictLookupOK = false;
    bool kmerFoundInDict = false;
    std::string kmerStr = "";
    seeding::seed_t newSeed = {}; // Initialize

    // ---> TRACE: Log dictionary lookup attempt <---
    logging::verbose("TRACE_ADD: Attempting lookup for Node={}, Pos={}", node->identifier, pos);

    // CRITICAL FIX: Add detailed diagnostic logging for node_3689
    bool isTargetNode = (node->identifier == "node_3689");
    if (isTargetNode) {
        logging::info("NODE-3689: Dictionary lookup at position {}", pos);
    }

    // First try to use dictionary lookup if available
    if (positionToDictId.count(pos) > 0 && positionToEndOffset.count(pos) > 0) {
        uint32_t dictId = positionToDictId[pos];
        
        // ---> TRACE: Log dictionary ID found <---
        logging::verbose("TRACE_ADD_DICT: Node={}, Pos={}, DictID={}", node->identifier, pos, dictId);
        // ---> END TRACE <---
        
        // ---> DEBUG: Log found Dict ID <---
        logging::debug("DEBUG_ADD: Node={}, Pos={}, Found DictID={}, EndOffset={}", 
                     node->identifier, pos, dictId, positionToEndOffset[pos]);

        if (state.kmerDictionary.count(dictId) > 0) {
            kmerStr = state.kmerDictionary.at(dictId);
            kmerFoundInDict = true;
            int64_t endPos = pos + positionToEndOffset[pos];
            
            // ---> TRACE: Log k-mer found in dictionary <---
            logging::verbose("TRACE_ADD_DICT_KMER: Node={}, Pos={}, DictID={}, Kmer='{}'", node->identifier, pos, dictId, kmerStr);
            // ---> END TRACE <---

            // ---> DEBUG: Log Kmer found in dictionary <---
            logging::debug("DEBUG_ADD: Node={}, Pos={}, DictID={}, Found Kmer='{}'", 
                         node->identifier, pos, dictId, kmerStr);


            // Create seed using dictionary kmer - function handles normalization now
            newSeed = createSeedFromDictionary(kmerStr, pos, endPos, params.k, params.s);
            dictLookupOK = true; // Mark lookup as successful
            
            addSeedAndUpdateScores(node, stateManager, state, result, uniqueSeedHashes, pos, newSeed);
            seedAdditions++;
            dictionaryLookups++;
            
            // ---> TRACE: Log seed added (from dict) and read lookup status <---
            bool inReads = state.seedFreqInReads.count(newSeed.hash) > 0;
            logging::verbose("TRACE_ADD_DICT_DONE: Node={}, Pos={}, Hash={}, InReads={}", 
                           node->identifier, pos, newSeed.hash, inReads);
            // ---> END TRACE <---
            
        } else {
            // Invalid Dictionary ID - log error
            logging::warn("Invalid dictionary ID {} for position {} in node {}", 
                       dictId, pos, node->identifier);
        }
    } else {
        // Missing Position/Offset Info - log issue
        logging::warn("No dictionary position/offset info for pos {} in node {}", 
                   pos, node->identifier);
    }
    
    // ---> TRACE: Log final outcome details <---
    logging::verbose("TRACE_ADD_FINAL: Node={}, Pos={}, Hash={}, Kmer='{}', DictOK={}, KmerInDict={}", 
                 node->identifier, pos, (dictLookupOK ? newSeed.hash : 0), kmerStr, dictLookupOK, kmerFoundInDict);
    // ---> END TRACE <---

    // ---> DEBUG: Log final outcome (only reachable on success now) <---
    if (dictLookupOK) {
        logging::debug("DEBUG_ADD_SUCCESS: Node={}, Pos={}, Kmer='{}', NewHash={}", 
                     node->identifier, pos, kmerStr, newSeed.hash);
    } else {
        // Log failure without throwing
        logging::warn("Failed to add seed at position {} in node {}", pos, node->identifier);
    }
    // ---> END DEBUG <---
}

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
    // First check if nodePathInfo is empty - this is a critical error
    if (state.nodePathInfo.size() == 0) {
        logging::err("ERROR: nodePathInfo structure is empty! This means the index was built without node path information.");
        logging::err("Please rebuild the index with the -f/--reindex flag to fix this issue.");
        return false;
    }
    
    // Log the search attempt for debugging
    logging::debug("Searching for node '{}' in index with {} nodePathInfo entries", 
                 nodeId, state.nodePathInfo.size());
                 
    // If this is the first search, log some sample node IDs from the index
    static bool first_search = true;
    if (first_search && state.nodePathInfo.size() > 0) {
        first_search = false;
        // logging::info("Sample node IDs in index:");
        for (size_t i = 0; i < std::min<size_t>(5, state.nodePathInfo.size()); i++) {
            auto indexNodeId = state.nodePathInfo[i].getNodeId();
            std::string indexNodeIdStr(indexNodeId.begin(), indexNodeId.end());
            // logging::info("  [{}]: '{}'", i, indexNodeIdStr);
        }
    }
    
    // Search through nodePathInfo for this node with proper string conversion
    for (size_t i = 0; i < state.nodePathInfo.size(); i++) {
        auto indexNodeId = state.nodePathInfo[i].getNodeId();
        // Convert Text::Reader to std::string for proper comparison
        std::string indexNodeIdStr(indexNodeId.begin(), indexNodeId.end());
        
        if (indexNodeIdStr == nodeId) {
            logging::debug("Found node '{}' at index {}", nodeId, i);
            nodeIndex = i;
            return true;
        }
    }
    
    logging::debug("Node '{}' not found in index (checked {} entries), using default scoring", 
                 nodeId, state.nodePathInfo.size());
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
 * @brief Update Jaccard (Presence/Absence) score for a node and track the best score
 * 
 * @param node The node being evaluated
 * @param score The Jaccard score based on presence/absence of seeds
 */
void PlacementResult::updateJaccardPresenceScore(panmanUtils::Node* node, double score, int64_t rawCount) {
    if (!node) return;
    
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Check if score is better than current best
    if (score > bestJaccardPresenceScore + TIED_THRESHOLD) {
        // Found a new best Jaccard (Presence/Absence) score
        bestJaccardPresenceScore = score;
        bestJaccardPresenceCount = rawCount;
        bestJaccardPresenceNode = node;
        
        // Reset tied nodes and add this one
        tiedJaccardPresenceNodes.clear();
        tiedJaccardPresenceNodes.push_back(node);
        
        logging::debug("New best Jaccard (Presence) node: {} with score {:.6f} ({} raw matches)", 
                     node->identifier, score, rawCount);
                     
    } else if (std::abs(score - bestJaccardPresenceScore) <= TIED_THRESHOLD) {
        // Add to tied nodes
        tiedJaccardPresenceNodes.push_back(node);
        logging::debug("Tied Jaccard (Presence) node: {} with score {:.6f} ({} raw matches)", 
                     node->identifier, score, rawCount);
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

// After dumpPlacementDebugData - Add a new utility function for debugging k-mer hashing

/**
 * @brief Dumps detailed k-mer extraction information from reads to help debug hashing issues
 * 
 * @param readSequences Vector of read sequences
 * @param k k-mer size
 * @param outputFilename File to write the debug data to
 */
void dumpKmerDebugData(
    const std::vector<std::string>& readSequences,
    int k,
    const std::string& outputFilename) {
    
    std::ofstream outFile(outputFilename);
    if (!outFile.is_open()) {
        logging::err("Could not open k-mer debug file: {}", outputFilename);
        return;
    }
    
    logging::info("Dumping detailed k-mer debug data to {}", outputFilename);
    
    outFile << "# K-mer Debug Data (k=" << k << ")\n\n";
    
    // Process a few reads and k-mers for detailed inspection
    size_t readCount = std::min(static_cast<size_t>(10), readSequences.size());
    
    for (size_t readIdx = 0; readIdx < readCount; readIdx++) {
        const std::string& readSeq = readSequences[readIdx];
        
        outFile << "## Read " << readIdx << " (length: " << readSeq.length() << ")\n\n";
        
        if (readSeq.length() < static_cast<size_t>(k)) {
            outFile << "Read too short to extract k-mers\n\n";
            continue;
        }
        
        // Create all-uppercase version of the read for comparison
        std::string upperSeq = readSeq;
        std::transform(upperSeq.begin(), upperSeq.end(), upperSeq.begin(), 
                      [](unsigned char c){ return std::toupper(c); });
                      
        // Check if conversion was needed
        bool isAlreadyUpper = (upperSeq == readSeq);
        outFile << "Read is already all uppercase: " << (isAlreadyUpper ? "Yes" : "No") << "\n\n";
        
        outFile << "| Position | Original K-mer | OrigFwd Hash | OrigRC Hash | Original Canonical | Uppercase K-mer | UpperFwd Hash | UpperRC Hash | Upper Canonical | Hash Match? |\n";
        outFile << "| -------- | -------------- | ------------ | ----------- | ----------------- | --------------- | ------------- | ------------ | --------------- | ----------- |\n";
        
        // Process first 10 k-mers in this read
        size_t kmerCount = std::min(static_cast<size_t>(10), readSeq.length() - k + 1);
        
        
        outFile << "\n";
        
        // Add extra section to display character frequencies in the read
        outFile << "### Character Frequencies\n\n";
        std::unordered_map<char, size_t> charCounts;
        
        for (char c : readSeq) {
            charCounts[c]++;
        }
        
        outFile << "| Character | Count | Percentage |\n";
        outFile << "| --------- | ----- | ---------- |\n";
        
        for (char c = 'A'; c <= 'z'; c++) {
            if (charCounts.find(c) != charCounts.end()) {
                double percent = (static_cast<double>(charCounts[c]) / readSeq.length()) * 100.0;
                outFile << "| '" << c << "' | " << charCounts[c] << " | " 
                       << std::fixed << std::setprecision(2) << percent << "% |\n";
            }
        }
        
        // Check for other characters (non-alphabetic)
        for (const auto& [c, count] : charCounts) {
            if (!std::isalpha(c)) {
                double percent = (static_cast<double>(count) / readSeq.length()) * 100.0;
                outFile << "| '" << c << "' | " << count << " | " 
                       << std::fixed << std::setprecision(2) << percent << "% |\n";
            }
        }
        
        outFile << "\n";
    }
    
    // Add a section to compare original vs uppercase for each k-mer in the index dictionary
    outFile << "## Case Sensitivity Analysis for Index Dictionary\n\n";
    
    // We need to access the PlacementGlobalState to get the dictionary
    // This is a limitation of this debug function, but we'll note it in the output
    outFile << "Note: This section would normally display index dictionary entries, but requires access to PlacementGlobalState.\n";
    outFile << "Check the debug logs for detailed dictionary information.\n\n";
    
    outFile.close();
    logging::info("K-mer debug data written to {}", outputFilename);
}

void placementTraversal(
    state::StateManager& stateManager,
    PlacementResult& result,
    panmanUtils::Tree* T, 
    PlacementGlobalState& state,
    const TraversalParams& params) {
    
    if (!T || !T->root) {
        throw std::invalid_argument("Invalid tree or root node for placement traversal");
    }
    
    const double TIED_THRESHOLD = 0.0001; // Define threshold for considering scores tied
    
    // Use a comprehensive approach with block-aware traversal and full seed processing
    logging::info("Performing placement traversal using quaternary-encoded seed mutations");
    
    // Setup progress tracking
    std::atomic<size_t> nodesProcessed{0};
    const size_t totalNodes = T->allNodes.size();
    const auto startTime = std::chrono::high_resolution_clock::now();
    
    
    
    logging::debug("Grouping nodes by level...");
    // Group nodes by level using common function from state
    const auto nodesByLevel = state::groupNodesByLevel(T, T->root);
    
    // Create a mutex for updating the best results
    std::mutex resultMutex;
    
    // Create a task arena for parallel processing
    int numThreads = std::min(static_cast<int>(std::thread::hardware_concurrency()), 
                             params.t > 0 ? params.t : 4);
    tbb::task_arena arena(numThreads);
    
    logging::info("Processing tree with {} nodes using {} threads", totalNodes, numThreads);
    
    std::ofstream countsFile("node_seed_counts.tsv");
    countsFile << "node_id\tseed_count\tdeletions\tinsertions\tmodifications\ttotal_seed_count\tseed_matches\tweighted_seed_matches\n";

    // Process each level in breadth-first order
    for (const auto& level : nodesByLevel) {
        logging::debug("Processing level with {} nodes", level.size());
        
        arena.execute([&]() {
            // OPTIMIZED: Use better grain size calculation for TBB parallel_for
            // Calculate optimal grain size based on level size and thread count
            const size_t levelSize = level.size();
            const size_t optimalGrainSize = std::max(1UL, 
                std::min(levelSize / (numThreads * 4), 64UL)); // Prevent too small or too large grains
            
            // Process nodes at this level in parallel with optimized partitioning
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, levelSize, optimalGrainSize),
                [&](const tbb::blocked_range<size_t>& r) {
                    // Thread-local result for best scores (existing ones)
                    double localBestJaccardScore = 0.0; // This is for weighted Jaccard
                    double localBestCosineScore = 0.0;
                    double localBestWeightedScore = 0.0;
                    panmanUtils::Node* localBestJaccardNode = nullptr;
                    panmanUtils::Node* localBestCosineNode = nullptr;
                    panmanUtils::Node* localBestWeightedNode = nullptr;
                    std::vector<panmanUtils::Node*> localTiedJaccardNodes;
                    std::vector<panmanUtils::Node*> localTiedCosineNodes;
                    std::vector<panmanUtils::Node*> localTiedWeightedNodes;

                    // Thread-local result for new scores
                    int64_t localBestRawSeedMatchScore = 0;
                    panmanUtils::Node* localBestRawSeedMatchNode = nullptr;
                    std::vector<panmanUtils::Node*> localTiedRawSeedMatchNodes;

                    double localBestJaccardPresenceScore = 0.0;
                    int64_t localBestJaccardPresenceCount = 0;
                    panmanUtils::Node* localBestJaccardPresenceNode = nullptr;
                    std::vector<panmanUtils::Node*> localTiedJaccardPresenceNodes;

                    // Thread-local collectors to reduce mutex contention
                    std::vector<std::pair<std::string, std::unordered_map<int64_t, seeding::seed_t>>> localNodeSeedMaps;
                    std::vector<std::string> localTsvOutputLines;
                    
                    // Reserve space for better performance
                    localNodeSeedMaps.reserve(r.end() - r.begin());
                    localTsvOutputLines.reserve(r.end() - r.begin());
                    

                    for (size_t i = r.begin(); i < r.end(); i++) {
                        panmanUtils::Node* node = level[i];
                        if (!node) continue;
                        
                        try {
                            // ---> ADDED: Log processed node ID <---
                            static std::mutex processed_nodes_log_mutex;
                            std::lock_guard<std::mutex> pn_lock(processed_nodes_log_mutex);
                            std::ofstream pn_log_file("processed_nodes_in_traversal.log", std::ios_base::app);
                            if (pn_log_file.is_open()) {
                                pn_log_file << node->identifier << "\n";
                                pn_log_file.close();
                            }
                            // --- END LOG ---

                            PlacementNodeScore nodeScore;
                            nodeScore.hitsInThisGenome = 0;
                            nodeScore.currentJaccardNumerator = 0;
                            nodeScore.currentCosineNumerator = 0.0;
                            nodeScore.currentCosineDenominator = 0.0;
                            nodeScore.kmerSeedMap.clear();
                            nodeScore.currentGenomeSeedCounts.clear();
                            nodeScore.rawSeedMatchScore = 0; // Initialize new field
                            nodeScore.jaccardPresenceNumerator = 0; // Initialize new field
                            nodeScore.currentGenomeUniqueSeedHashes.clear(); // Initialize new field
                            
                            absl::flat_hash_set<size_t> uniqueSeedHashes; // This is for weighted Jaccard denominator, based on currentGenomeSeedCounts keys
                            uniqueSeedHashes.clear();


                            if (node->parent != nullptr) {
                                std::lock_guard<std::mutex> lock(resultMutex); // Protects result.nodeSeedMap
                                if (result.nodeSeedMap.count(node->parent->identifier)) {
                                    const auto& parent_map = result.nodeSeedMap.at(node->parent->identifier);
                                    nodeScore.kmerSeedMap.reserve(parent_map.size());
                                    
                                    // Inherit seeds from parent
                                    for (const auto& entry : parent_map) {
                                        nodeScore.kmerSeedMap.insert(entry);
                                    }
                                    
                                    // Initialize scores from inherited seeds
                                    for (const auto& pair : nodeScore.kmerSeedMap) {
                                        const seeding::seed_t& inherited_seed = pair.second;
                                        size_t inherited_hash = inherited_seed.hash;
                                        
                                        nodeScore.currentGenomeUniqueSeedHashes.insert(inherited_hash); // For presence/absence Jaccard
                                        nodeScore.currentGenomeSeedCounts[inherited_hash]++; // For weighted Jaccard and Cosine
                                        uniqueSeedHashes.insert(inherited_hash); // For weighted Jaccard denominator

                                        auto readIt = state.seedFreqInReads.find(inherited_hash);
                                        if (readIt != state.seedFreqInReads.end()) {
                                            int64_t readCount = readIt->second;
                                            nodeScore.hitsInThisGenome += readCount; // Used for weighted Jaccard numerator
                                            nodeScore.currentJaccardNumerator += readCount; // Used for weighted Jaccard numerator
                                            nodeScore.rawSeedMatchScore += readCount; // For raw seed match score
                                            nodeScore.jaccardPresenceNumerator++; // For presence/absence Jaccard numerator
                                            
                                            auto [num_delta, den_delta] = getCosineDelta(false, true, inherited_hash, state.seedFreqInReads, nodeScore.currentGenomeSeedCounts);
                                            nodeScore.currentCosineNumerator += num_delta;
                                            nodeScore.currentCosineDenominator += den_delta;
                                        }
                                    }
                                }
                            }
                            
                            uint64_t nodeIndex;
                            bool nodeFound = getNodeIndex(node->identifier, state, nodeIndex);
                            
                            if (!nodeFound) {
                                logging::warn("Node {} not found in index, skipping", node->identifier);
                                continue;
                            }

                            bool shouldLog = node->identifier == "node_434";
                            
                            auto seedMutation = state.perNodeSeedMutations[nodeIndex];
                            auto basePositions = seedMutation.getBasePositions();
                            auto perPosMasks = seedMutation.getPerPosMasks();
                            
                            auto kmerDictionaryIds = seedMutation.getKmerDictionaryIds();
                            auto kmerPositions = seedMutation.getKmerPositions();
                            auto kmerEndOffsets = seedMutation.getKmerEndOffsets();
                            
                            absl::flat_hash_map<int64_t, uint32_t> positionToDictId;
                            absl::flat_hash_map<int64_t, uint32_t> positionToEndOffset;
                            
                            // Ensure all three lists kmerDictionaryIds, kmerPositions, kmerEndOffsets are accessed safely
                            // Assuming kmerPositions is the reference for actual count of dictionary-linked seeds
                            for (size_t j = 0; j < kmerPositions.size(); j++) {
                                if (j < kmerDictionaryIds.size() && j < kmerEndOffsets.size()) {
                                    positionToDictId[kmerPositions[j]] = kmerDictionaryIds[j];
                                    positionToEndOffset[kmerPositions[j]] = kmerEndOffsets[j];
                                } else {
                                    // Log if there's a mismatch, but potentially continue if non-critical
                                    logging::warn("Node {}: Mismatch in kmer dictionary-related array lengths at index {}. Pos: {}. DictSize: {}, PosSize: {}, OffsetSize: {}", 
                                                  node->identifier, j, kmerPositions[j], kmerDictionaryIds.size(), kmerPositions.size(), kmerEndOffsets.size());
                                    // Decide if this is a fatal error or if we can skip this entry
                                }
                            }

                           

                            // Initialize tracking variables for seed operations
                            size_t deletionCount = 0;
                            size_t insertionCount = 0;
                            size_t modificationCount = 0;

                            if (basePositions.size() > 0 && perPosMasks.size() > 0 && basePositions.size() == perPosMasks.size()) {
                                for (size_t j = 0; j < basePositions.size(); j++) {
                                        int64_t basePos = basePositions[j];
                                    uint64_t mask = perPosMasks[j];
                                        
                                        for (uint8_t offset = 0; offset < 32; offset++) {
                                            uint8_t value = (mask >> (offset * 2)) & 0x3;
                                                int64_t pos = basePos - offset;
                                    
                                        
                                        if (value == 0) continue; 

                                        seeding::seed_t current_seed_at_pos{}; // Default initialize
                                        bool existed_before_mutation = nodeScore.kmerSeedMap.count(pos);
                                        if (existed_before_mutation) {
                                            current_seed_at_pos = nodeScore.kmerSeedMap.at(pos);
                                        }

                                        if (value == 1) { // Delete
                                            deletionCount++;
                                            if (existed_before_mutation) {
                                                size_t deleted_hash = current_seed_at_pos.hash;
                                                
                                                
                                                nodeScore.kmerSeedMap.erase(pos);

                                                if (state.seedFreqInReads.count(deleted_hash)) {
                                                    int64_t delReadCount = state.seedFreqInReads.at(deleted_hash);
                                                    nodeScore.hitsInThisGenome -= delReadCount;
                                                    nodeScore.currentJaccardNumerator -= delReadCount;
                                                    nodeScore.rawSeedMatchScore -= delReadCount; // Update raw seed match score
                                                    nodeScore.jaccardPresenceNumerator--;    // Update presence/absence Jaccard numerator

                                                    auto [num_delta, den_delta] = getCosineDelta(true, false, deleted_hash, state.seedFreqInReads, nodeScore.currentGenomeSeedCounts);
                                                    nodeScore.currentCosineNumerator += num_delta;
                                                    nodeScore.currentCosineDenominator += den_delta;
                                                }
                                                
                                                // Safely decrement and erase from currentGenomeSeedCounts and uniqueSeedHashes
                                                auto it_genome_counts = nodeScore.currentGenomeSeedCounts.find(deleted_hash);
                                                if (it_genome_counts != nodeScore.currentGenomeSeedCounts.end()) {
                                                    it_genome_counts->second--;
                                                    if (it_genome_counts->second == 0) {
                                                        nodeScore.currentGenomeSeedCounts.erase(it_genome_counts);
                                                        uniqueSeedHashes.erase(deleted_hash); // For weighted Jaccard denominator
                                                        nodeScore.currentGenomeUniqueSeedHashes.erase(deleted_hash); // For presence/absence Jaccard
                                                    }
                                                }
                                            }
                                        } else { // Add (value == 2) or Modify (value == 3)
                                            // std::cerr << " ANd now this." << std::endl;
                                            if (!positionToDictId.count(pos)) {
                                                logging::warn("Node {}: Cannot find dictionary ID for position {} to add/modify seed.", node->identifier, pos);
                                                continue;
                                            }
                                            uint32_t dictId = positionToDictId.at(pos);
                                            if (!state.kmerDictionary.count(dictId)){
                                                logging::warn("Node {}: Dictionary ID {} (for pos {}) not found in global kmerDictionary.", node->identifier, dictId, pos);
                                                continue;
                                            }

                                            const std::string& kmerStr = state.kmerDictionary.at(dictId);
                                            int64_t end_offset_val = positionToEndOffset.count(pos) ? positionToEndOffset.at(pos) : params.k -1;
                                            seeding::seed_t new_seed = createSeedFromDictionary(kmerStr, pos, pos + end_offset_val, params.k, params.s);
                                            size_t new_hash = new_seed.hash;

                                            // Track operation type
                                            if (value == 2) {
                                                insertionCount++;
                                            } else if (value == 3) {
                                                modificationCount++;
                                            }

                                            if (value == 3 && existed_before_mutation) { 
                                                size_t old_hash = current_seed_at_pos.hash;
                                                if (old_hash != new_hash) { 
                                                    if (state.seedFreqInReads.count(old_hash)) {
                                                        int64_t oldReadCount = state.seedFreqInReads.at(old_hash);
                                                        nodeScore.hitsInThisGenome -= oldReadCount;
                                                        nodeScore.currentJaccardNumerator -= oldReadCount;
                                                        nodeScore.rawSeedMatchScore -= oldReadCount; // Update raw seed match score
                                                        nodeScore.jaccardPresenceNumerator--;    // Update presence/absence Jaccard numerator

                                                        auto [num_delta, den_delta] = getCosineDelta(true, false, old_hash, state.seedFreqInReads, nodeScore.currentGenomeSeedCounts);
                                                        nodeScore.currentCosineNumerator += num_delta;
                                                        nodeScore.currentCosineDenominator += den_delta;
                                                    }
                                                    auto it_genome_counts_old = nodeScore.currentGenomeSeedCounts.find(old_hash);
                                                    if (it_genome_counts_old != nodeScore.currentGenomeSeedCounts.end()) {
                                                        it_genome_counts_old->second--;
                                                        if (it_genome_counts_old->second == 0) {
                                                            nodeScore.currentGenomeSeedCounts.erase(it_genome_counts_old);
                                                            uniqueSeedHashes.erase(old_hash); // For weighted Jaccard denominator
                                                            nodeScore.currentGenomeUniqueSeedHashes.erase(old_hash); // For presence/absence Jaccard
                                                        }
                                                    }
                                                }
                                            }

                                            nodeScore.kmerSeedMap[pos] = new_seed;
                                            bool new_hash_was_absent_in_genome_counts = !nodeScore.currentGenomeSeedCounts.count(new_hash) || nodeScore.currentGenomeSeedCounts.at(new_hash) == 0;
                                            if (new_hash_was_absent_in_genome_counts) {
                                               uniqueSeedHashes.insert(new_hash); // For weighted Jaccard denominator
                                               nodeScore.currentGenomeUniqueSeedHashes.insert(new_hash); // For presence/absence Jaccard
                                            }
                                            nodeScore.currentGenomeSeedCounts[new_hash]++;
                                            
                                            if (state.seedFreqInReads.count(new_hash)) {
                                                int64_t addReadCount = state.seedFreqInReads.at(new_hash);
                                                bool new_hash_was_absent_in_read_matches = !nodeScore.currentGenomeUniqueSeedHashes.count(new_hash); // Check before insert if this is the first match for this hash

                                                // Only add to scores if the hash truly changed or if it's a new addition (value==2)
                                                // Or if it was modified from a different hash
                                                if (value == 2 || (value == 3 && (!existed_before_mutation || current_seed_at_pos.hash != new_hash))) {
                                                   nodeScore.hitsInThisGenome += addReadCount;
                                                   nodeScore.currentJaccardNumerator += addReadCount;
                                                   nodeScore.rawSeedMatchScore += addReadCount; // Update raw seed match score
                                                   // Update presence/absence Jaccard numerator only if it's a truly new hash match for reads
                                                   if (new_hash_was_absent_in_genome_counts || (value ==3 && current_seed_at_pos.hash != new_hash) || value == 2) {
                                                       if(state.seedFreqInReads.count(new_hash)) nodeScore.jaccardPresenceNumerator++; 
                                                   }

                                                   auto [num_delta, den_delta] = getCosineDelta(false, true, new_hash, state.seedFreqInReads, nodeScore.currentGenomeSeedCounts);
                                                   nodeScore.currentCosineNumerator += num_delta;
                                                   nodeScore.currentCosineDenominator += den_delta;
                                                } else if (value == 3 && existed_before_mutation && current_seed_at_pos.hash == new_hash) {
                                                   // Hash is the same, scores related to read counts should not change unless genome count was 0 before.
                                                   // The currentCosineDenominator might change if this is the first instance of this seed hash in genome.
                                                   // This is handled by getCosineDelta when genome count goes from 0 to 1 for this hash.
                                                   if (new_hash_was_absent_in_genome_counts) { // if it was 0, and now 1
                                                        auto [num_delta, den_delta] = getCosineDelta(false, true, new_hash, state.seedFreqInReads, nodeScore.currentGenomeSeedCounts);
                                                        nodeScore.currentCosineNumerator += num_delta; 
                                                        nodeScore.currentCosineDenominator += den_delta; 
                                                        // If this is the first time this specific hash is seen for this node (e.g. count went from 0 to 1)
                                                        // and it's in reads, increment jaccardPresenceNumerator
                                                        if(state.seedFreqInReads.count(new_hash)) nodeScore.jaccardPresenceNumerator++;
                                                   }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            
                            // Calculate Weighted Jaccard Score
                            double weightedJaccardScore = 0.0;
                            // Denominator for weighted Jaccard: (total unique seeds in reads) + (total unique seeds in current genome) - (weighted intersection)
                            // state.jaccardDenominator is total unique seeds in reads (state.seedFreqInReads.size())
                            // uniqueSeedHashes.size() is total unique seeds in current genome (keys of nodeScore.currentGenomeSeedCounts)
                            // nodeScore.currentJaccardNumerator is the sum of read frequencies of matched seeds (weighted intersection)
                            double weightedJaccardDenominator = state.jaccardDenominator + uniqueSeedHashes.size() - nodeScore.currentJaccardNumerator;
                            if (weightedJaccardDenominator > 0) {
                                weightedJaccardScore = static_cast<double>(nodeScore.currentJaccardNumerator) / weightedJaccardDenominator;
                            }

                            // Calculate Jaccard (Presence/Absence) Score
                            double jaccardPresenceScore = 0.0;
                            double jaccardPresenceDenominator_val = static_cast<double>(state.readUniqueSeedCount) + nodeScore.currentGenomeUniqueSeedHashes.size() - nodeScore.jaccardPresenceNumerator;
                            if (jaccardPresenceDenominator_val > 0) {
                                jaccardPresenceScore = static_cast<double>(nodeScore.jaccardPresenceNumerator) / jaccardPresenceDenominator_val;
                            }
                            
                            // Calculate cosine similarity
                            double cosineScore = 0.0;
                                        if (nodeScore.currentCosineDenominator > 0 && state.totalReadSeedCount > 0) {
                                            const double normGenome = std::sqrt(std::abs(nodeScore.currentCosineDenominator));
                        const double normReads = std::sqrt(state.totalReadSeedCount);
                                            cosineScore = nodeScore.currentCosineNumerator / (normGenome * normReads);
                                        }
                                        
                                        // Calculate weighted score
                                        double weightedScore = params.scoreScale * weightedJaccardScore + 
                                                              (1.0 - params.scoreScale) * cosineScore;
                                        
                                        // Check if this is one of our special debug nodes
                                        bool is_special_node = (node->identifier == "KU950624.1" || node->identifier == "DJ068250.1");
                                            
                                        // Only do debug logging for special nodes
                                        if (is_special_node) {
                                            std::stringstream debug_ss;
                                            debug_ss << "DEBUG_PLACEMENT_NODE: " << node->identifier << "\n";
                                            debug_ss << "  RawSeedMatchScore: " << nodeScore.rawSeedMatchScore << "\n";
                                            debug_ss << "  WeightedJaccard (adapted): " << weightedJaccardScore << " (num: " << nodeScore.currentJaccardNumerator << ", den_reads_unique: " << state.jaccardDenominator << ", den_genome_unique_hashes: " << uniqueSeedHashes.size() << ")\n";
                                            debug_ss << "  JaccardPresence: " << jaccardPresenceScore << " (num: " << nodeScore.jaccardPresenceNumerator << ", den_reads_unique: " << state.readUniqueSeedCount << ", den_genome_unique: " << nodeScore.currentGenomeUniqueSeedHashes.size() << ")\n";
                                            debug_ss << "  CosineScore: " << cosineScore << " (num: " << nodeScore.currentCosineNumerator << ", den_genome_sqrt: " << std::sqrt(std::abs(nodeScore.currentCosineDenominator)) << ", den_reads_sqrt: " << std::sqrt(state.totalReadSeedCount) << ")\n";
                                            debug_ss << "  OverallWeightedScore: " << weightedScore << "\n";
                                            debug_ss << "  Total seeds in this node's kmerSeedMap (after muts): " << nodeScore.kmerSeedMap.size() << "\n";
                                            debug_ss << "  Contributing seeds to RawSeedMatchScore (hash: kmer_from_dict (read_freq)) :\n";
                                            
                                            int seeds_logged_count = 0;
                                            std::vector<std::pair<size_t, int64_t>> contributing_seeds; // hash, read_freq

                                            for (const auto& entry : nodeScore.kmerSeedMap) { // entry is pair<int64_t_pos, seeding::seed_t>
                                                const seeding::seed_t& seed_on_node = entry.second;
                                                auto read_freq_it = state.seedFreqInReads.find(seed_on_node.hash);
                                                if (read_freq_it != state.seedFreqInReads.end() && read_freq_it->second > 0) {
                                                    contributing_seeds.push_back({seed_on_node.hash, read_freq_it->second});
                                                }
                                            }
                                            // Sort by read frequency descending to see most impactful
                                            std::sort(contributing_seeds.rbegin(), contributing_seeds.rend(), [](const auto&a, const auto&b){
                                                return a.second < b.second; // Sorts descending by frequency
                                            });

                                            for(const auto& p : contributing_seeds){
                                                std::string kmer_str = "<kmer_not_found_in_dict>";
                                                // Try to find kmer string for this hash via any kmerPosition that might have this hash
                                                // This is indirect; ideally, we'd have hash->dictId or kmerSeedMap would store dictId
                                                // For now, we iterate kmerPositions to find one that yields this hash.
                                                bool kmer_found_for_hash = false;
                                                // positionToDictId is local to this node's processing in the loop, derived from its seedMutation entry
                                                for(const auto& kmer_pos_original_entry : positionToDictId){ // positionToDictId is absl::flat_hash_map<int64_t, uint32_t>
                                                    int64_t kmer_pos_original = kmer_pos_original_entry.first;
                                                    uint32_t dict_id_original = kmer_pos_original_entry.second;

                                                    if(kmer_pos_original < 0) continue; 
                                                    auto kmer_map_it = state.kmerDictionary.find(dict_id_original);
                                                    if(kmer_map_it != state.kmerDictionary.end()){
                                                        auto end_offset_it = positionToEndOffset.find(kmer_pos_original);
                                                        if(end_offset_it != positionToEndOffset.end()){
                                                            seeding::seed_t temp_seed = createSeedFromDictionary(kmer_map_it->second, kmer_pos_original, kmer_pos_original + end_offset_it->second, params.k, params.s);
                                                            if(temp_seed.hash == p.first){
                                                                kmer_str = kmer_map_it->second;
                                                                kmer_found_for_hash = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                if(!kmer_found_for_hash && state.hashToKmer.count(p.first)){ // Fallback to read kmer if populated
                                                     kmer_str = state.hashToKmer.at(p.first) + " (from read_hashToKmer)";
                                                }

                                                debug_ss << "    - " << p.first << ": " << kmer_str << " (read_freq: " << p.second << ")\n";
                                                seeds_logged_count++;
                                            }

                                            if (seeds_logged_count == 0) {
                                                debug_ss << "    - (No seeds on this node matched seeds found in reads and populated in kmerSeedMap)\n";
                                            }
                                            
                                            // --- > NEW: Log all seeds in final kmerSeedMap (up to 50) < ---
                                            debug_ss << "  All seeds in this node's final kmerSeedMap (pos: hash kmer_str read_found? read_freq) - max 50 shown:\n";
                                            int all_seeds_log_count = 0;
                                            // Sort kmerSeedMap by position for consistent logging order
                                            std::vector<std::pair<int64_t, seeding::seed_t>> sorted_kmer_map(nodeScore.kmerSeedMap.begin(), nodeScore.kmerSeedMap.end());
                                            std::sort(sorted_kmer_map.begin(), sorted_kmer_map.end(), [](const auto&a, const auto&b){
                                                return a.first < b.first;
                                            });

                                            for (const auto& map_entry : sorted_kmer_map) {
                                                if (all_seeds_log_count >= 50) break;
                                                int64_t pos = map_entry.first;
                                                const seeding::seed_t& seed = map_entry.second;
                                                std::string kmer_s = "<kmer_lookup_failed>";
                                                if (state.hashToKmer.count(seed.hash)) {
                                                    kmer_s = state.hashToKmer.at(seed.hash);
                                                }
                                                bool found_in_reads = state.seedFreqInReads.count(seed.hash);
                                                int64_t read_freq = found_in_reads ? state.seedFreqInReads.at(seed.hash) : 0;
                                                debug_ss << "    - Pos: " << pos << ", Hash: " << seed.hash 
                                                         << ", Kmer: " << kmer_s 
                                                         << ", InReads: " << (found_in_reads ? "Yes" : "No") 
                                                         << ", ReadFreq: " << read_freq << "\n";
                                                all_seeds_log_count++;
                                            }
                                            if (all_seeds_log_count == 0 && !nodeScore.kmerSeedMap.empty()) {
                                                debug_ss << "    - (kmerSeedMap is not empty, but failed to log entries - check limits or map content)\n";
                                            } else if (nodeScore.kmerSeedMap.empty()) {
                                                debug_ss << "    - (kmerSeedMap is empty for this node)\n";
                                            }
                                            // --- > END NEW LOG SECTION < ---
                                            
                                            // Write to the debug file
                                            static std::mutex placement_debug_log_mutex;
                                            std::lock_guard<std::mutex> lock(placement_debug_log_mutex);
                                            
                                            std::ofstream debug_file("placement_debug_scores.log", std::ios_base::app);
                                            if (debug_file.is_open()) {
                                                debug_file << debug_ss.str() << std::endl; // Add an extra newline for readability between entries
                                                debug_file.close();
                                            } else {
                                                // Fallback or error logging if file can't be opened
                                                std::cerr << "ERROR: Could not open placement_debug_scores.log for writing node " << node->identifier << std::endl;
                                            }
                                            debug_ss << "\n==== DETAILED DEBUG FOR SPECIAL NODE: " << node->identifier << " ====\n";
                                            
                                            // 1. Print ALL seeds in the node (not just contributing ones)
                                            debug_ss << "\nALL SEEDS IN NODE (sorted by read frequency):\n";
                                            
                                            // Collect all seeds first for sorting
                                            std::vector<std::tuple<int64_t, size_t, std::string, int64_t, bool>> all_seeds;
                                            
                                            for (const auto& entry : nodeScore.kmerSeedMap) {
                                                int64_t pos = entry.first;
                                                const seeding::seed_t& seed = entry.second;
                                                size_t hash = seed.hash;
                                                std::string kmer_str = "<unknown>";
                                                
                                                // Try to get kmer string from hashToKmer map
                                                if (state.hashToKmer.count(hash)) {
                                                    kmer_str = state.hashToKmer.at(hash);
                                                }
                                                
                                                // Get read frequency if it exists
                                                int64_t read_count = 0;
                                                bool in_reads = state.seedFreqInReads.count(hash) > 0;
                                                if (in_reads) {
                                                    read_count = state.seedFreqInReads.at(hash);
                                                }
                                                
                                                all_seeds.emplace_back(pos, hash, kmer_str, read_count, in_reads);
                                            }
                                            
                                            // Sort seeds by read frequency (descending), then by position
                                            std::sort(all_seeds.begin(), all_seeds.end(),
                                                [](const auto& a, const auto& b) {
                                                    // Primary sort by read frequency (descending)
                                                    if (std::get<3>(a) != std::get<3>(b)) {
                                                        return std::get<3>(a) > std::get<3>(b);
                                                    }
                                                    // Secondary sort by presence in reads
                                                    if (std::get<4>(a) != std::get<4>(b)) {
                                                        return std::get<4>(a); // true comes before false
                                                    }
                                                    // Tertiary sort by position
                                                    return std::get<0>(a) < std::get<0>(b);
                                                });
                                            
                                            // Print all seeds in sorted order
                                            debug_ss << "  Total seeds: " << all_seeds.size() << "\n\n";
                                            
                                            for (const auto& [pos, hash, kmer_str, read_count, in_reads] : all_seeds) {
                                                debug_ss << "  Pos: " << pos << ", Hash: " << hash 
                                                         << ", Kmer: " << kmer_str 
                                                         << ", ReadCount: " << read_count 
                                                         << ", InReads: " << (in_reads ? "YES" : "no") << "\n";
                                            }
                                            
                                            // 2. Print all seeds that matched with reads and their details
                                            debug_ss << "\nREAD-MATCHING SEEDS (sorted by frequency):\n";
                                            std::vector<std::tuple<size_t, std::string, int64_t, int64_t>> matching_seeds;
                                            
                                            for (const auto& entry : nodeScore.kmerSeedMap) {
                                                int64_t pos = entry.first;
                                                const seeding::seed_t& seed = entry.second;
                                                size_t hash = seed.hash;
                                                
                                                if (state.seedFreqInReads.count(hash) > 0) {
                                                    int64_t read_count = state.seedFreqInReads.at(hash);
                                                    std::string kmer_str = "<unknown>";
                                                    if (state.hashToKmer.count(hash)) {
                                                        kmer_str = state.hashToKmer.at(hash);
                                                    }
                                                    matching_seeds.emplace_back(hash, kmer_str, read_count, pos);
                                                }
                                            }
                                            
                                            // Sort by read frequency (descending)
                                            std::sort(matching_seeds.begin(), matching_seeds.end(), 
                                                [](const auto& a, const auto& b) {
                                                    return std::get<2>(a) > std::get<2>(b);
                                                });
                                            
                                            for (const auto& [hash, kmer, freq, pos] : matching_seeds) {
                                                debug_ss << "  Hash: " << hash 
                                                         << ", Kmer: " << kmer 
                                                         << ", ReadFreq: " << freq 
                                                         << ", Position: " << pos << "\n";
                                            }
                                            
                                            // 3. Add information about seed mutations specific to this node
                                            debug_ss << "\nSEED MUTATIONS FOR THIS NODE:\n";
                                            uint64_t nodeIndex;
                                            if (getNodeIndex(node->identifier, state, nodeIndex)) {
                                                auto seedMutation = state.perNodeSeedMutations[nodeIndex];
                                                auto basePositions = seedMutation.getBasePositions();
                                                auto perPosMasks = seedMutation.getPerPosMasks();
                                                
                                                debug_ss << "  Total mutation entries: " << basePositions.size() << "\n";
                                                
                                                // Print detailed mutation info (print all entries)
                                                size_t mut_count = 0;
                                                
                                                // First collect all mutations for sorting
                                                std::vector<std::tuple<int64_t, std::string, std::string, int64_t>> all_mutations;
                                                
                                                for (size_t j = 0; j < basePositions.size(); j++) {
                                                    int64_t basePos = basePositions[j];
                                                    uint64_t mask = perPosMasks[j];
                                                    
                                                    for (uint8_t offset = 0; offset < 32; offset++) {
                                                        uint8_t value = (mask >> (offset * 2)) & 0x3;
                                                        int64_t pos = basePos - offset;
                                                        
                                                        if (value == 0) continue; // Skip no-op
                                                        
                                                        mut_count++;
                                                        std::string op_type = 
                                                            (value == 1) ? "DELETE" : 
                                                            (value == 2) ? "ADD" : "MODIFY";
                                                        
                                                        // Try to find kmer info for this position
                                                        std::string kmer_info = "<unknown>";
                                                        if (positionToDictId.count(pos) && 
                                                            state.kmerDictionary.count(positionToDictId[pos])) {
                                                            kmer_info = state.kmerDictionary.at(positionToDictId[pos]);
                                                        }
                                                        
                                                        // Find read frequency for this kmer
                                                        int64_t read_freq = 0;
                                                        if (!kmer_info.empty() && kmer_info != "<unknown>") {
                                                            // Compute hash for this kmer
                                                            auto syncmers = seeding::rollingSyncmers(kmer_info, params.k, params.s, false, 0, false);
                                                            if (!syncmers.empty()) {
                                                                size_t hash = std::get<0>(syncmers[0]);
                                                                if (state.seedFreqInReads.count(hash) > 0) {
                                                                    read_freq = state.seedFreqInReads.at(hash);
                                                                }
                                                            }
                                                        }
                                                        
                                                        all_mutations.emplace_back(pos, op_type, kmer_info, read_freq);
                                                    }
                                                }
                                                
                                                // Sort mutations by read frequency (if kmer matched)
                                                std::sort(all_mutations.begin(), all_mutations.end(),
                                                    [](const auto& a, const auto& b) {
                                                        // Primary sort by read frequency (descending)
                                                        if (std::get<3>(a) != std::get<3>(b)) {
                                                            return std::get<3>(a) > std::get<3>(b);
                                                        }
                                                        // Secondary sort by position
                                                        return std::get<0>(a) < std::get<0>(b);
                                                    });
                                                
                                                // Print all mutations
                                                for (const auto& [pos, op_type, kmer_info, read_freq] : all_mutations) {
                                                    debug_ss << "  " << op_type << " at pos " << pos 
                                                             << ", Kmer: " << kmer_info;
                                                    
                                                    if (read_freq > 0) {
                                                        debug_ss << ", ReadFreq: " << read_freq << " (MATCHES READS)";
                                                    }
                                                    
                                                    debug_ss << "\n";
                                                }
                                                
                                                if (mut_count == 0) {
                                                    debug_ss << "  No active mutations found\n";
                                                }
                                            } else {
                                                debug_ss << "  Node not found in index by getNodeIndex\n";
                                            }
                                            
                                            debug_ss << "==== END DETAILED DEBUG ====\n";
                                        }
                                        
                                        // END DEBUG LOGGING
                                        
                                        // Update local best scores
                                        // const double TIED_THRESHOLD = 0.0001; // Moved to function scope
                                        
                                        // Update best Jaccard (this is the weighted Jaccard)
                                        if (weightedJaccardScore > localBestJaccardScore + TIED_THRESHOLD) {
                                            localBestJaccardScore = weightedJaccardScore;
                                            localBestJaccardNode = node;
                                            localTiedJaccardNodes.clear();
                                            localTiedJaccardNodes.push_back(node);
                                        } else if (std::abs(weightedJaccardScore - localBestJaccardScore) <= TIED_THRESHOLD) {
                                            localTiedJaccardNodes.push_back(node);
                                        }

                                        // Update best Raw Seed Match Score (Local)
                                        if (nodeScore.rawSeedMatchScore > localBestRawSeedMatchScore) {
                                            localBestRawSeedMatchScore = nodeScore.rawSeedMatchScore;
                                            localBestRawSeedMatchNode = node;
                                            localTiedRawSeedMatchNodes.clear();
                                            localTiedRawSeedMatchNodes.push_back(node);
                                        } else if (nodeScore.rawSeedMatchScore == localBestRawSeedMatchScore) {
                                            localTiedRawSeedMatchNodes.push_back(node);
                                        }

                                        // Update best Jaccard (Presence/Absence) (Local)
                                         if (jaccardPresenceScore > localBestJaccardPresenceScore + TIED_THRESHOLD) {
                                            localBestJaccardPresenceScore = jaccardPresenceScore;
                                            localBestJaccardPresenceCount = nodeScore.jaccardPresenceNumerator;
                                            localBestJaccardPresenceNode = node;
                                            localTiedJaccardPresenceNodes.clear();
                                            localTiedJaccardPresenceNodes.push_back(node);
                                        } else if (std::abs(jaccardPresenceScore - localBestJaccardPresenceScore) <= TIED_THRESHOLD) {
                                            localTiedJaccardPresenceNodes.push_back(node);
                                        }
                                        
                                        // Update best Cosine
                                        if (cosineScore > localBestCosineScore + TIED_THRESHOLD) {
                                            localBestCosineScore = cosineScore;
                                            localBestCosineNode = node;
                                            localTiedCosineNodes.clear();
                                            localTiedCosineNodes.push_back(node);
                                        } else if (std::abs(cosineScore - localBestCosineScore) <= TIED_THRESHOLD) {
                                            localTiedCosineNodes.push_back(node);
                                        }
                                        
                                        // Update best Weighted
                                        if (weightedScore > localBestWeightedScore + TIED_THRESHOLD) {
                                            localBestWeightedScore = weightedScore;
                                            localBestWeightedNode = node;
                                            localTiedWeightedNodes.clear();
                                            localTiedWeightedNodes.push_back(node);
                                        } else if (std::abs(weightedScore - localBestWeightedScore) <= TIED_THRESHOLD) {
                                            localTiedWeightedNodes.push_back(node);
                                        }
                                        
                                        // Collect node seed map for batch processing (thread-local)
                                        std::unordered_map<int64_t, seeding::seed_t> temp_map;
                                        temp_map.reserve(nodeScore.kmerSeedMap.size());
                                        for (const auto& entry : nodeScore.kmerSeedMap) {
                                            temp_map.insert(entry);
                                        }
                                        localNodeSeedMaps.emplace_back(node->identifier, std::move(temp_map));
                                        
                                        // Prepare TSV output for batch writing (thread-local)
                                        // Calculate additional metrics
                                        int64_t totalSeedCount = nodeScore.kmerSeedMap.size(); // Use actual seed map size
                                        
                                        int64_t seedMatches = nodeScore.jaccardPresenceNumerator;
                                        int64_t weightedSeedMatches = nodeScore.rawSeedMatchScore;
                                        
                                        // Collect TSV line for batch writing
                                        std::ostringstream tsvLine;
                                        tsvLine << node->identifier << "\t" 
                                               << nodeScore.kmerSeedMap.size() << "\t"
                                               << deletionCount << "\t"
                                               << insertionCount << "\t"
                                               << modificationCount << "\t"
                                               << totalSeedCount << "\t"
                                               << seedMatches << "\t"
                                               << weightedSeedMatches << "\n";
                                        localTsvOutputLines.push_back(tsvLine.str());
                                        
                                        nodesProcessed++;
                                        
                        } catch (const std::exception& e) {
                            logging::err("Error processing node {}: {}", node->identifier, e.what());
                        }


                        // std::cerr << node->identifier << "\t" << localBestWeightedScore << "\t" << localBestJaccardScore << "\t" << localBestRawSeedMatchScore << "\t" << localBestCosineScore << "\t" << localBestJaccardPresenceScore << std::endl;
                    }
                    
                    // Batch process thread-local collections to reduce mutex contention
                    {
                        std::lock_guard<std::mutex> lock(resultMutex);
                        
                        // Batch insert node seed maps
                        for (const auto& [nodeId, seedMap] : localNodeSeedMaps) {
                            result.nodeSeedMap[nodeId] = seedMap;
                        }
                        
                        // Batch write TSV output lines
                        for (const std::string& line : localTsvOutputLines) {
                            countsFile << line;
                        }
                    }
                    
                    // Merge results with global best
                    {
                        std::lock_guard<std::mutex> lock(resultMutex);
                        
                        // Update best Jaccard (Weighted)
                        if (localBestJaccardScore > result.bestJaccardScore + TIED_THRESHOLD) { 
                            result.bestJaccardScore = localBestJaccardScore;
                            result.bestJaccardNode = localBestJaccardNode;
                            result.tiedJaccardNodes = localTiedJaccardNodes;
                        } else if (std::abs(localBestJaccardScore - result.bestJaccardScore) <= TIED_THRESHOLD) { 
                            result.tiedJaccardNodes.insert(
                                result.tiedJaccardNodes.end(),
                                localTiedJaccardNodes.begin(),
                                localTiedJaccardNodes.end());
                        }

                        // Update Raw Seed Match Score (Global)
                        if (localBestRawSeedMatchNode) {
                            // The main best node from this thread for this metric
                            result.updateRawSeedMatchScore(localBestRawSeedMatchNode, localBestRawSeedMatchScore);
                            // Process other nodes from this thread that tied for this thread's local best raw seed score
                            for (panmanUtils::Node* tied_node : localTiedRawSeedMatchNodes) {
                                // updateRawSeedMatchScore will correctly add to global ties if the score matches the global best
                                // No need to check if tied_node == localBestRawSeedMatchNode, update function handles it.
                                result.updateRawSeedMatchScore(tied_node, localBestRawSeedMatchScore);
                            }
                        }

                        // Update Jaccard (Presence/Absence) Score (Global)
                        if (localBestJaccardPresenceNode) {
                            result.updateJaccardPresenceScore(localBestJaccardPresenceNode, localBestJaccardPresenceScore, localBestJaccardPresenceCount);
                            for (panmanUtils::Node* tied_node : localTiedJaccardPresenceNodes) {
                                result.updateJaccardPresenceScore(tied_node, localBestJaccardPresenceScore, localBestJaccardPresenceCount);
                            }
                        }
                        
                        // Update best Cosine
                        if (localBestCosineScore > result.bestCosineScore + TIED_THRESHOLD) {
                            result.bestCosineScore = localBestCosineScore;
                            result.bestCosineNode = localBestCosineNode;
                            result.tiedCosineNodes = localTiedCosineNodes;
                        } else if (std::abs(localBestCosineScore - result.bestCosineScore) <= TIED_THRESHOLD) {
                            result.tiedCosineNodes.insert(
                                result.tiedCosineNodes.end(),
                                localTiedCosineNodes.begin(),
                                localTiedCosineNodes.end());
                        }
                        
                        // Update best Weighted
                        if (localBestWeightedScore > result.bestWeightedScore + TIED_THRESHOLD) {
                            result.bestWeightedScore = localBestWeightedScore;
                            result.bestWeightedNode = localBestWeightedNode;
                            result.tiedWeightedNodes = localTiedWeightedNodes;
                        } else if (std::abs(localBestWeightedScore - result.bestWeightedScore) <= TIED_THRESHOLD) {
                            result.tiedWeightedNodes.insert(
                                result.tiedWeightedNodes.end(),
                                localTiedWeightedNodes.begin(),
                                localTiedWeightedNodes.end());
                        }
                    }

                }
            );
        });
    }
    
    // Set performance metrics
    auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::high_resolution_clock::now() - startTime).count();
    result.totalTimeSeconds = totalTime;
    
    // Write simple placement result files
    try {
        // Write raw seed match results
        std::ofstream rawFile("placement.raw.txt");
        if (rawFile.is_open()) {
            if (result.bestRawSeedMatchNode) {
                rawFile << result.bestRawSeedMatchNode->identifier << "\t" << result.bestRawSeedMatchScore << std::endl;
            } else {
                rawFile << "None\t0" << std::endl;
            }
            rawFile.close();
        }
        
        // Write weighted Jaccard results  
        std::ofstream weightedFile("placement.weighted.txt");
        if (weightedFile.is_open()) {
            if (result.bestJaccardNode) {
                weightedFile << result.bestJaccardNode->identifier << "\t" << result.bestJaccardScore << std::endl;
            } else {
                weightedFile << "None\t0.0" << std::endl;
            }
            weightedFile.close();
        }
        
        // Write cosine similarity results
        std::ofstream cosineFile("placement.cosine.txt");
        if (cosineFile.is_open()) {
            if (result.bestCosineNode) {
                cosineFile << result.bestCosineNode->identifier << "\t" << result.bestCosineScore << std::endl;
            } else {
                cosineFile << "None\t0.0" << std::endl;
            }
            cosineFile.close();
        }
        
        logging::debug("Written placement result files: placement.raw.txt, placement.weighted.txt, placement.cosine.txt");
    } catch (const std::exception& e) {
        logging::warn("Failed to write placement result files: {}", e.what());
    }
    
    // Log final results
    /* logging::info("Placement completed in {}s, processed {} nodes", 
                  static_cast<long long>(totalTime), 
                  static_cast<long long>(nodesProcessed.load())); */
    // logging::info("Best node by weighted score: {}", 
    //             result.bestWeightedNode ? result.bestWeightedNode->identifier : "None");
} // End of placementTraversal

// Add implementation of dumpSyncmerDetails function here
void dumpSyncmerDetails(const std::string& filename, 
                      const std::string& label,
                      const absl::flat_hash_map<size_t, int64_t>& seedMap,
                      const absl::flat_hash_map<size_t, std::string>& kmerMap) {
    std::ofstream outFile(filename, std::ios_base::app); // Append mode
    if (!outFile.is_open()) {
        logging::warn("Failed to open syncmer debug file: {}", filename);
        return;
    }
    
    outFile << "\n=== " << label << " Syncmers ===\n";
    outFile << "Hash\tSequence\tCount\n";
    
    for (const auto& [hash, count] : seedMap) {
        std::string kmerSeq = "<unknown>";
        if (kmerMap.count(hash) > 0) {
            kmerSeq = kmerMap.at(hash);
        }
        
        outFile << hash << "\t" << kmerSeq << "\t" << count << "\n";
    }
    
    outFile.close();
    logging::info("Wrote {} syncmers to {}", seedMap.size(), filename);
}

// Dump placement debug data to file
void dumpPlacementDebugData(
    state::StateManager& stateManager,
    PlacementGlobalState& state,
    const std::vector<std::string>& nodesToDump,
    size_t maxNodes,
    const TraversalParams& params,
    const std::string& outputFilename,
    const std::vector<std::string>& readSequences) {
    
    logging::info("Dumping placement debug data to: {}", outputFilename);
    
    // Write debug info to file
    std::ofstream outFile(outputFilename);
    if (!outFile.is_open()) {
        logging::warn("Failed to open debug output file: {}", outputFilename);
        return;
    }
    
    // Write basic placement parameters
    outFile << "# Placement Debug Information\n\n";
    outFile << "Parameters:\n";
    outFile << "  k: " << params.k << "\n";
    outFile << "  s: " << params.s << "\n";
    outFile << "  scoreScale: " << params.scoreScale << "\n";
    
    // Write read information
    outFile << "\nRead Information:\n";
    outFile << "  Total reads: " << readSequences.size() << "\n";
    outFile << "  Unique read seeds: " << state.seedFreqInReads.size() << "\n";
    outFile << "  Total read seed count: " << state.totalReadSeedCount << "\n";
    outFile << "  Jaccard denominator: " << state.jaccardDenominator << "\n";
    
    // Sample read k-mers if available
    if (!state.seedFreqInReads.empty()) {
        outFile << "\nSample Read Seeds:\n";
        size_t count = 0;
        for (const auto& [hash, freq] : state.seedFreqInReads) {
            if (count++ >= 10) break;
            std::string kmer = "<unknown>";
            if (state.hashToKmer.count(hash) > 0) {
                kmer = state.hashToKmer.at(hash);
            }
            outFile << "  Hash: " << hash << ", Frequency: " << freq << ", K-mer: " << kmer << "\n";
        }
    }
    
    // Write placement index information
    outFile << "\nPlacement Index Information:\n";
    outFile << "  Nodes in index: " << state.nodePathInfo.size() << "\n";
    outFile << "  Dictionary entries: " << state.kmerDictionary.size() << "\n";
    
    // Sample some node information from the index
    if (state.nodePathInfo.size() > 0) {
        outFile << "\nSample Nodes from Index:\n";
        const size_t nodesToShow = (maxNodes < state.nodePathInfo.size()) ? maxNodes : state.nodePathInfo.size();
        for (size_t i = 0; i < nodesToShow; i++) {
            auto nodeInfo = state.nodePathInfo[i];
            std::string nodeId(nodeInfo.getNodeId().begin(), nodeInfo.getNodeId().end());
            std::string parentId = "<none>";
            if (nodeInfo.hasParentId()) {
                parentId = std::string(nodeInfo.getParentId().begin(), nodeInfo.getParentId().end());
            }
            
            outFile << "  Node: " << nodeId << "\n";
            outFile << "    Level: " << nodeInfo.getLevel() << "\n";
            outFile << "    Parent: " << parentId << "\n";
            outFile << "    Active Blocks: " << nodeInfo.getActiveBlocks().size() << "\n";
        }
    }
    
    // Add specific nodes to dump if requested
    if (!nodesToDump.empty()) {
        outFile << "\nRequested Nodes to Dump:\n";
        for (const auto& nodeId : nodesToDump) {
            outFile << "  Node: " << nodeId << "\n";
            
            // Try to find this node in the index
            bool foundInIndex = false;
            for (size_t i = 0; i < state.nodePathInfo.size(); i++) {
                auto indexNodeId = state.nodePathInfo[i].getNodeId();
                std::string indexNodeIdStr(indexNodeId.begin(), indexNodeId.end());
                if (indexNodeIdStr == nodeId) {
                    foundInIndex = true;
                    auto nodeInfo = state.nodePathInfo[i];
                    outFile << "    Found in index at position " << i << "\n";
                    outFile << "    Level: " << nodeInfo.getLevel() << "\n";
                    outFile << "    Active Blocks: " << nodeInfo.getActiveBlocks().size() << "\n";
                    break;
                }
            }
            
            if (!foundInIndex) {
                outFile << "    Not found in index\n";
            }
        }
    }
    
    outFile.close();
    logging::info("Debug data dump complete");
}

/**
 * @brief Load seed information from index into the StateManager for placement
 * 
 * This function extracts seed info from the index but doesn\'t activate blocks.
 * It\'s a lightweight approach focused only on seed hash collection for comparison.
 * 
 * @param stateManager StateManager to populate with seed data
 * @param index The index to read seed data from
 */
void loadSeedsFromIndex(state::StateManager& stateManager, const ::Index::Reader& index) {
    auto perNodeSeedMutations = index.getPerNodeSeedMutations();
    logging::info("Loading seeds from index for {} nodes...", perNodeSeedMutations.size());
    
    // FIXED: Defensive check for empty seed mutations
    if (perNodeSeedMutations.size() == 0) {
        logging::warn("Index contains no seed mutations - this will likely cause issues with placement");
        return;
    }
    
    // Get block information from the index
    auto blockInfoList = index.getBlockInfo();
    logging::info("Populating block ranges in StateManager from index with {} blocks defined in BlockInfo.", blockInfoList.size());

    // Iterate through the blocks defined in the index and set their ranges in StateManager
    for (auto block : blockInfoList) {
        int32_t blockId = block.getBlockId();
        // Directly get rangeStart and rangeEnd. They default to 0 if not set in the index.
        uint64_t startCoord = block.getRangeStart();
        uint64_t endCoord = block.getRangeEnd();

        if (startCoord == 0 && endCoord == 0 && blockId != 0) { // blockId 0 might legitimately be [0,0] if empty
             logging::debug("Block {} from index has range [0, 0]. This might be an uninitialized or genuinely empty block.", blockId);
        }
        
        coordinates::CoordRange range{ (int64_t)startCoord, (int64_t)endCoord };
        stateManager.setBlockRange(blockId, range);
        logging::debug("Set block {} range from index: [{}, {}]", blockId, range.start, range.end);
    }
    
    // The old loop that tried to initialize block ranges has been replaced by the loop above.
    
    // First check if there are any seed mutations in the index
    size_t totalMutationEntries = 0;
    size_t nonEmptyMutationNodes = 0;
    
    logging::debug("DEBUG_SEEDS: Starting detailed seed diagnostics...");
    
    // Check various fields in the index
    try {
        logging::debug("DEBUG_SEEDS: Index k={}, s={}, nodes in mutation list={}", 
                      index.getK(), index.getS(), perNodeSeedMutations.size());
        
        // Check if the k-mer dictionary exists
        if (index.hasKmerDictionary()) {
            auto kmerDict = index.getKmerDictionary();
            logging::debug("DEBUG_SEEDS: K-mer dictionary has {} entries", kmerDict.size());
            
            // Log a few entries if they exist
            if (kmerDict.size() > 0) {
                for (size_t i = 0; i < (kmerDict.size() < 5 ? kmerDict.size() : 5); i++) {
                    // Use to_string to avoid formatter issues
                    std::string seqStr = std::string(kmerDict[i].getSequence());
                    logging::debug("DEBUG_SEEDS: Dictionary entry {}: seq={}, canonical={}",
                                 i, seqStr, kmerDict[i].getCanonicalForm());
                }
            }
        } else {
            logging::warn("DEBUG_SEEDS: K-mer dictionary is not present in the index!");
        }
        
        // Check for node path info
        if (index.hasNodePathInfo()) {
            auto nodePathInfo = index.getNodePathInfo();
            logging::debug("DEBUG_SEEDS: NodePathInfo has {} entries", nodePathInfo.size());
            
            // Log a few entries if they exist
            if (nodePathInfo.size() > 0) {
                for (size_t i = 0; i < (nodePathInfo.size() < 5 ? nodePathInfo.size() : 5); i++) {
                    std::string nodeIdStr = std::string(nodePathInfo[i].getNodeId());
                    logging::debug("DEBUG_SEEDS: NodePathInfo entry {}: nodeId={}, level={}, activeBlocks={}",
                                 i, nodeIdStr, nodePathInfo[i].getLevel(), 
                                 nodePathInfo[i].getActiveBlocks().size());
                }
            }
        } else {
            logging::warn("DEBUG_SEEDS: NodePathInfo is not present in the index!");
        }
    } catch (const std::exception& e) {
        logging::err("DEBUG_SEEDS: Error inspecting index: {}", e.what());
    }
    
    // Check each node for seed mutations
    for (size_t dfsIndex = 0; dfsIndex < perNodeSeedMutations.size(); dfsIndex++) {
        auto mutations = perNodeSeedMutations[dfsIndex];
        auto basePositions = mutations.getBasePositions();
        auto perPosMasks = mutations.getPerPosMasks();
        auto kmerPositions = mutations.getKmerPositions();
        
        totalMutationEntries += basePositions.size();
        
        if (basePositions.size() > 0 || kmerPositions.size() > 0) {
            nonEmptyMutationNodes++;
            
            // Log first few for diagnostics
            if (nonEmptyMutationNodes <= 5) {
                logging::info("DEBUG_SEEDS: Found node at DFS index {} with {} base positions and {} kmer positions",
                             dfsIndex, basePositions.size(), kmerPositions.size());
                
                // Log detailed content for the first few entries
                if (basePositions.size() > 0) {
                    for (size_t i = 0; i < (basePositions.size() < 3 ? basePositions.size() : 3); i++) {
                        logging::debug("DEBUG_SEEDS: Base position {}: pos={}, mask={}",
                                     i, basePositions[i], 
                                     i < perPosMasks.size() ? perPosMasks[i] : 0);
                    }
                }
                
                if (kmerPositions.size() > 0) {
                    for (size_t i = 0; i < (kmerPositions.size() < 3 ? kmerPositions.size() : 3); i++) {
                        logging::debug("DEBUG_SEEDS: KmerPosition {}: pos={}",
                                     i, kmerPositions[i]);
                    }
                }
                
                // Find node ID for this DFS index
                try {
                    if (dfsIndex < stateManager.nodeIdsByDfsIndex.size()) {
                        std::string nodeId = stateManager.nodeIdsByDfsIndex[dfsIndex];
                        logging::debug("DEBUG_SEEDS: DFS index {} corresponds to node ID: {}", 
                                     dfsIndex, nodeId);
                        
                        // Check active blocks for this node
                        try {
                            const auto& activeBlocks = stateManager.getActiveBlocks(nodeId);
                            logging::debug("DEBUG_SEEDS: Node {} has {} active blocks", 
                                         nodeId, activeBlocks.size());
                        } catch (const std::exception& e) {
                            logging::warn("DEBUG_SEEDS: Cannot get active blocks for node {}: {}", 
                                        nodeId, e.what());
                        }
                    } else {
                        logging::warn("DEBUG_SEEDS: DFS index {} is out of bounds for nodeIdsByDfsIndex (size {})", 
                                    dfsIndex, stateManager.nodeIdsByDfsIndex.size());
                    }
                } catch (const std::exception& e) {
                    logging::warn("DEBUG_SEEDS: Error mapping DFS index {} to node ID: {}", 
                                dfsIndex, e.what());
                }
            }
        }
    }
    
    logging::info("MUTATION_SUMMARY: Index contains {} total mutation entries across {} nodes with mutations", 
                 totalMutationEntries, nonEmptyMutationNodes);
    
    if (nonEmptyMutationNodes == 0) {
        logging::warn("NO MUTATIONS FOUND IN INDEX! This explains why no seeds are being loaded.");
        
        // Check if we have a root node or if DFS indices are initialized
        logging::debug("DEBUG_SEEDS: StateManager has {} nodes in nodeIdsByDfsIndex", 
                     stateManager.nodeIdsByDfsIndex.size());
        
        // Check the first few nodes in the index by DFS
        for (size_t i = 0; i < (stateManager.nodeIdsByDfsIndex.size() < 5 ? stateManager.nodeIdsByDfsIndex.size() : 5); i++) {
            std::string nodeId = stateManager.nodeIdsByDfsIndex[i];
            logging::debug("DEBUG_SEEDS: Node at DFS index {}: {}", i, nodeId);
            
            try {
                // Check if there are active blocks
                const auto& activeBlocks = stateManager.getActiveBlocks(nodeId);
                logging::debug("DEBUG_SEEDS: Node {} has {} active blocks", nodeId, activeBlocks.size());
                
                // Check if it has any mutations by direct check
                try {
                    auto& nodeState = stateManager.getNodeState(nodeId);
                    size_t blockMutations = 0;
                    size_t nucMutations = 0;
                    
                    // Log nucleotide changes if any
                    if (!nodeState.nucleotideChanges.empty()) {
                        nucMutations = nodeState.nucleotideChanges.size();
                        logging::debug("DEBUG_SEEDS: Node {} has {} nucleotide changes in its state", 
                                     nodeId, nucMutations);
                    }
                    
                    // Log if there are recomputation ranges
                    if (!nodeState.recompRanges.empty()) {
                        logging::debug("DEBUG_SEEDS: Node {} has {} recomputation ranges", 
                                     nodeId, nodeState.recompRanges.size());
                    }
                    
                    // Check seed changes directly
                    if (!nodeState.seedChangeBasePositions.empty() && !nodeState.seedChangeBitMasks.empty()) {
                        logging::debug("DEBUG_SEEDS: Node {} has {} seed change base positions and {} bit masks", 
                                     nodeId, nodeState.seedChangeBasePositions.size(), 
                                     nodeState.seedChangeBitMasks.size());
                    } else {
                        logging::debug("DEBUG_SEEDS: Node {} has NO seed changes in its state", nodeId);
                    }
                    
                } catch (const std::exception& e) {
                    logging::warn("DEBUG_SEEDS: Error checking node state for {}: {}", nodeId, e.what());
                }
            } catch (const std::exception& e) {
                logging::warn("DEBUG_SEEDS: Error getting active blocks for node {}: {}", nodeId, e.what());
            }
        }
        
        return;
    }
    
    // We don\'t need to actually load or store the seeds here
    // They will be extracted directly during placement in collectExistingSeeds
    logging::info("Seed info verified in index - will be extracted during placement");
    
    // Show some detailed statistics about seeds in the index
    logging::debug("DEBUG_SEEDS: Index contains {} total base position entries across {} nodes with mutations",
                 totalMutationEntries, nonEmptyMutationNodes);
}

// Main placement function - simplified and consolidated
void place(
    PlacementResult& result,
    panmanUtils::Tree* T,
    ::Index::Reader& index,
    const std::string& reads1Path,
    const std::string& reads2Path,
    std::vector<std::vector<seeding::seed_t>>& readSeeds,
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
        logging::debug("Loading index directly from file: {}", indexPath);
        
        try {
            // Load and keep the message reader alive
            messageHolder = std::shared_ptr<::capnp::MessageReader>(loadIndexFromFile(indexPath).release());
            
            // Now safely get a reference from the message holder
            index = messageHolder->getRoot<Index>();
            logging::debug("Successfully loaded index from: {}", indexPath);
        } catch (const std::exception& e) {
            logging::err("Failed to load index from file: {}", e.what());
            throw std::runtime_error("Failed to load index from file: " + std::string(e.what()));
        }
    }
    
    // CRITICAL FIX: Validate the index message has a root pointer
    try {
        // Verify the index message has valid data
        uint32_t k = index.getK();
        uint32_t s = index.getS();
        size_t nodeCount = index.getPerNodeSeedMutations().size();
        
        if (k == 0 || s == 0 || nodeCount == 0) {
            throw std::runtime_error(
                "Invalid index: k=" + std::to_string(k) + 
                ", s=" + std::to_string(s) + 
                ", nodes=" + std::to_string(nodeCount) + 
                ". The index appears to be corrupted or incomplete.");
        }
        
        logging::debug("Verified index has valid root pointer with k={}, s={}, nodes={}", 
                     k, s, nodeCount);
    } catch (const ::kj::Exception& e) {
        throw std::runtime_error("Cap\'n Proto error: Message did not contain a valid root pointer. The index file appears to be corrupted or not properly initialized.");
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error validating index: ") + e.what());
    }
    
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
    params.debug_node_id = debug_node_id_param; // <-- SET THE PARAM IN STRUCT
    
    logging::debug("Placement parameters: k={}, s={}, t={}, open={}, scoreScale={}, debug_node_id='{}'",
                params.k, params.s, params.t, params.open ? "true" : "false", params.scoreScale, params.debug_node_id);
    
    if (!T || !T->root) {
        throw std::runtime_error("Invalid tree or root node in place() function");
    }
    
    // Initialize state manager using the lighter, optimized version
    // logging::info("Initializing StateManager (light) with k={}, s={}", params.k, params.s);
    auto stateManager = indexing::initializeStateManagerLight(T, T->root, params.k, params.s);

    // Initialize global placement state
    PlacementGlobalState state;
    state.perNodeSeedMutations = index.getPerNodeSeedMutations();
    state.perNodeGapMutations = index.getPerNodeGapMutations();
    state.nodePathInfo = index.getNodePathInfo();
    state.blockInfo = index.getBlockInfo();
    state.ancestorMatrix = index.getAncestorMatrix();

    // ADDED: Explicitly initialize hierarchical seed storage
    logging::info("PLACEMENT: Explicitly initializing hierarchical seed storage");
    stateManager->initializeSeedStorage();

    // ADDED: Load seeds from index into the StateManager (now simplified)
    loadSeedsFromIndex(*stateManager, index);

    // Initialize k-mer dictionary
    if (index.hasKmerDictionary()) {
        auto kmerDictionary = index.getKmerDictionary();
        for (uint32_t i = 0; i < kmerDictionary.size(); i++) {
            auto entry = kmerDictionary[i];
            std::string kmerSeq = entry.getSequence();
            state.kmerDictionary[i] = kmerSeq;
        }
        logging::info("Loaded {} k-mers from dictionary", kmerDictionary.size());
        logging::info("Loaded {} k-mers from dictionary into state.kmerDictionary", state.kmerDictionary.size());
    }
    // --- > NEW: Populate state.hashToKmer from state.kmerDictionary < ---
    state.hashToKmer.clear();
    // std::cerr << "This is the kmer table built from the genome index: " << std::endl;
    if (!state.kmerDictionary.empty() && params.k > 0 && params.s > 0 && params.k > params.s) {
        state.hashToKmer.reserve(state.kmerDictionary.size());
        logging::info("Populating state.hashToKmer from state.kmerDictionary ({} entries) using k={}, s={}...", 
                      state.kmerDictionary.size(), params.k, params.s);
        int processed_for_hash_map = 0;
        for (const auto& dict_entry : state.kmerDictionary) { // dict_entry is pair<uint32_t, std::string>
            const std::string& kmerSeqStr = dict_entry.second;
            if (!kmerSeqStr.empty()) {
                std::string upperKmerSeqStr = kmerSeqStr;
                // Normalize k-mer to uppercase for consistent hashing
                std::transform(upperKmerSeqStr.begin(), upperKmerSeqStr.end(), upperKmerSeqStr.begin(),
                              [](unsigned char c){ return std::toupper(c); });
                
                auto syncmers = seeding::rollingSyncmers(upperKmerSeqStr, params.k, params.s, false, 0, false);
                if (!syncmers.empty()) {
                    size_t hash = std::get<0>(syncmers[0]);
                    bool is_syncmer = std::get<2>(syncmers[0]);
                    if (is_syncmer) {
                        state.hashToKmer[hash] = kmerSeqStr; // Store original k-mer from dict
                        // std::cerr << "Hash: " << hash << ", Kmer: " << kmerSeqStr << std::endl;
                    }
                }
            }
            processed_for_hash_map++;
            if (processed_for_hash_map % 200000 == 0 || processed_for_hash_map == state.kmerDictionary.size()) {
                 logging::debug("Processed {}/{} entries for state.hashToKmer.", processed_for_hash_map, state.kmerDictionary.size());
            }
        }
        logging::info("Populated state.hashToKmer with {} unique hash-to-kmer mappings.", state.hashToKmer.size());
    } else {
        if (state.kmerDictionary.empty()) {
            logging::warn("state.kmerDictionary is empty, cannot populate state.hashToKmer.");
        } else {
            logging::warn("Invalid k/s params (k={}, s={}), cannot populate state.hashToKmer.", params.k, params.s);
        }
    }
    // --- > END POPULATE state.hashToKmer < ---
    
    // Pre-initialize nodes in breadth-first order to ensure proper parent-child relationship
    logging::debug("Pre-initializing nodes in breadth-first order");

    if (T->root) {
        // First, ensure the root node is initialized
        stateManager->initializeNode(T->root->identifier);
        logging::debug("Initialized root node: {}", T->root->identifier);
        
        // Then initialize the rest level by level (breadth-first)
        // This ensures parent nodes are always initialized before children
        std::vector<panmanUtils::Node*> currentLevel = {T->root};
        std::vector<panmanUtils::Node*> nextLevel;
        
        while (!currentLevel.empty()) {
            for (auto* node : currentLevel) {
                // Add all children to the next level
                for (auto* child : node->children) {
                    if (child) {
                        nextLevel.push_back(child);
                        
                        // Initialize this child node
                        // Parent is guaranteed to be initialized already
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
    
    // ---> DEBUG: Verify Index Load <---
    logging::debug("DEBUG: Loaded index has {} kmer dictionary entries.", state.kmerDictionary.size());
    logging::debug("DEBUG: Loaded index has seed mutations for {} nodes.", state.perNodeSeedMutations.size());
    // ---> END DEBUG <---

    // Process reads to get seed counts
    if (reads1Path.empty() && reads2Path.empty()) {
        logging::warn("No read files provided, just loading the index. Exiting.");
        return;
    }
    
    logging::debug("Processing reads from {}{}", 
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
            
        logging::debug("Processed {} reads in {}ms, found {} unique seeds", 
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
        state.readUniqueSeedCount = readSeedCounts.size(); // Initialize new field for presence/absence Jaccard
        
        // std::cerr << "These are the seeds from the reads: " << std::endl;
        // // print all read seeds, hashes and their kmers
        // for (int i = 0; i < readSeeds.size(); i++) {
        //     for (int j = 0; j < readSeeds[i].size(); j++) {
        //         std::string kmer_str = readSeedSeqs[i][j];
        //         std::cerr << "Seed Hash: " << readSeeds[i][j].hash << ", Kmer: " << kmer_str << std::endl;
        //     }
        // }

        logging::debug("Total seed occurrences in reads: {}", state.totalReadSeedCount);
        
        // ---> DEBUG: Verify Read Seeds <---
        logging::debug("DEBUG: Populated state.seedFreqInReads with {} unique seed hashes.", state.seedFreqInReads.size());
        if (!state.seedFreqInReads.empty()) {
            int print_count = 0;
            logging::debug("DEBUG: First few read seed hashes and counts:");
            for (const auto& [hash, count] : state.seedFreqInReads) {
                logging::debug("  Hash: {}, Count: {}", hash, count);
                if (++print_count >= 5) break;
            }
        }
        // ---> END DEBUG <---

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
    logging::debug("Starting placement traversal");
    
    try {
        placementTraversal(*stateManager, result, T, state, params);
        
        auto traversalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - traversalStart).count();
            
        logging::debug("Placement traversal completed in {}ms", traversalDuration);
        
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

    // --- > DUMP READ SEED HASHES < ---  // THIS IS THE CORRECT BLOCK TO KEEP
    std::ofstream read_seeds_log_file_correct("read_seeds_debug.log"); // Variable is read_seeds_log_file_correct
    if (read_seeds_log_file_correct.is_open()) { // CORRECTED
        std::ofstream read_seeds_log_file_correct("read_seeds_debug.log"); // Renamed variable to avoid conflict if any old one remains
        if (read_seeds_log_file_correct.is_open()) {
            read_seeds_log_file_correct << "# Read Seed Hashes and Frequencies (Total Unique: " << state.seedFreqInReads.size() << ", Total Occurrences: " << state.totalReadSeedCount << ")\n";
            read_seeds_log_file_correct << "# Hash\tFrequency\tKmer (if found in state.hashToKmer)\n";
            for (const auto& entry : state.seedFreqInReads) { // entry is std::pair<const size_t, int64_t>
                std::string kmer_str = "<kmer_not_in_read_fallback_map>";
                if (state.hashToKmer.count(entry.first)) {
                    kmer_str = state.hashToKmer.at(entry.first);
                }
                read_seeds_log_file_correct << entry.first << "\t" << entry.second << "\t" << kmer_str << "\n"; 
            }
            read_seeds_log_file_correct.close(); 
            // logging::info("Dumped read seed frequencies to read_seeds_debug.log");
        } else {
            logging::warn("Could not open read_seeds_debug.log for writing.");
        }
        // --- > END DUMP < ---

        if (state.totalReadSeedCount == 0) {
            logging::warn("No seeds found in reads. Check read format and k/s parameters.");
            return;
        }
    }
} // Closing brace for void place(...)

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
    const int& max_tied_nodes,
    const std::string& indexPath,
    const std::string& debug_node_id_param) {  
    
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
        std::string line_str(lineBuffer);
        
        // Remove trailing newline if present
        if (!line_str.empty() && (line_str.back() == '\n' || line_str.back() == '\r')) {
            line_str.pop_back();
        }
        
        // Skip empty lines and comments
        if (line_str.empty() || line_str[0] == '#') continue;
        
        // Parse line - format is: sample_name,reads1[,reads2]
        std::vector<std::string> parts;
        boost::split(parts, line_str, boost::is_any_of(","));
        
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
        
        std::vector<std::vector<seeding::seed_t>> readSeeds;
        std::vector<std::string> readSequences;
        std::vector<std::string> readNames;
        std::vector<std::string> readQuals;

        try {
            // Call the main place function with correct parameters
            place(result, T, index, reads1Path, reads2Path, readSeeds, readSequences, readNames, readQuals, outputFile, indexPath, debug_node_id_param);
            
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

/**
 * @brief Dumps a summary of the placement results to a Markdown file.
 *
 * @param result The PlacementResult object containing the results.
 * @param outputFilename The path to the output Markdown file.
 */
void dumpPlacementSummary(const placement::PlacementResult& result, const std::string& outputFilename) {
    std::ofstream outFile(outputFilename);
    if (!outFile.is_open()) {
        logging::warn("Could not open placement summary file for writing: {}", outputFilename);
        return;
    }

    logging::debug("Dumping placement summary to: {}", outputFilename);

    outFile << "# Panmap Placement Summary\\n\\n";
    outFile << "| Metric                  | Value                 | Best Node(s)                     |\\n";
    outFile << "| :---------------------- | :-------------------- | :------------------------------- |\\n";

    // Helper lambda to format node lists
    auto formatNodes = [](const panmanUtils::Node* bestNode, const std::vector<panmanUtils::Node*>& tiedNodes) -> std::string {
        if (!bestNode) {
            return "N/A";
        }
        std::string nodeStr = bestNode->identifier;
        if (tiedNodes.size() > 1) {
            nodeStr += " (Tied: ";
            int count = 0;
            for (const auto* node : tiedNodes) {
                if (node != bestNode) { // Avoid repeating the best node
                    if (count > 0) nodeStr += ", ";
                    nodeStr += node->identifier;
                    if (++count >= 3) { // Limit tied nodes shown
                       nodeStr += ", ...";
                    break;
                }
            }
            }
            nodeStr += ")";
        }
        return nodeStr;
    };

    // --- Placement Scores ---
    outFile << "| Raw Seed Matches        | " << result.bestRawSeedMatchScore << " | "
            << formatNodes(result.bestRawSeedMatchNode, result.tiedRawSeedMatchNodes) << " |\\n";
    outFile << "| Best Jaccard Score      | " << std::fixed << std::setprecision(5) << result.bestJaccardScore << " | " 
            << formatNodes(result.bestJaccardNode, result.tiedJaccardNodes) << " |\\n";
    outFile << "| Jaccard (Presence)      | " << std::fixed << std::setprecision(5) << result.bestJaccardPresenceScore << " | "
            << formatNodes(result.bestJaccardPresenceNode, result.tiedJaccardPresenceNodes) << " |\\n";
    outFile << "| Best Cosine Score       | " << std::fixed << std::setprecision(5) << result.bestCosineScore << " | " 
            << formatNodes(result.bestCosineNode, result.tiedCosineNodes) << " |\\n";
    outFile << "| Best Weighted Score     | " << std::fixed << std::setprecision(5) << result.bestWeightedScore << " | " 
            << formatNodes(result.bestWeightedNode, result.tiedWeightedNodes) << " |\\n";

    // --- Performance Metrics ---
    outFile << "\\n## Performance Metrics\\n\\n";
    outFile << "| Metric                  | Value                 |\\n";
    outFile << "| :---------------------- | :-------------------- |\\n";
    outFile << "| Total Reads Processed   | " << result.totalReadsProcessed << " |\\n";
    outFile << "| Total Placement Time    | " << result.totalTimeSeconds << " seconds |\\n";
    
    outFile.close();
    logging::debug("Placement summary written successfully.");
}

} // namespace placement


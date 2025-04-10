#include "state.hpp"
#include "caller_logging.hpp"
#include "coordinates.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "panman.hpp"
#include "seeding.hpp"
#include "tbb/enumerable_thread_specific.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iterator>
#include <memory>
#include <mutex> // For std::mutex and std::lock_guard
#include <optional>
#include <ostream>
#include <queue>
#include <shared_mutex> // For std::shared_mutex
#include <stdexcept>
#include <string>
#include <string_view>
#include <tbb/tbb.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <sstream>
#include <iomanip>

namespace state {

/*
 * Character Data Access Pattern:
 * 
 * Hierarchical delta-based approach to storing character data:
 * 
 * 1. Initial sequences are stored in the global blockSequences map
 * 2. Each node only stores mutations relative to its parent
 * 3. Ancestor relationship is maintained via parentId reference
 * 4. When accessing characters:
 *    a. First check current node's characterData for mutations
 *    b. If not found, recursively check each ancestor node
 *        - When inversion status differs, transform coordinates before lookup
 *        - Position mirroring: pos → (length-1-pos) 
 *        - Gap position handling: reverse traversal order
 *    c. If not found in any ancestor, get from global blockSequences
 *    d. Apply appropriate inversion transformations (position and nucleotide)
 * 
 * This approach minimizes memory usage by:
 * - Storing only reference sequences once
 * - Avoiding duplication of mutations across the tree
 * - Each node only stores its specific differences
 */

using ::coordinates::CoordRange;
using ::coordinates::GapUpdate;
using ::gap_map::GapRange;


StateManager::StateManager() : numCoords(0), numBlocks(0), kmerSize(0) {
  // We'll update the block coordinate mappings when block ranges are populated
}

StateManager::StateManager(size_t numCoordinates)
    : numCoords(numCoordinates), numBlocks(0), kmerSize(0) {
  positionSeeds.resize(numCoordinates);
  
  // Initialize thread-local caches
  // Pre-allocate a reasonable number of entries for each cache
  for (auto& cache : gapListLengthCaches) {
    // Cache is initialized by default constructor, just ensure it's ready
    bool dummy = false;
    cache.get(0, 0, dummy); // Touch the cache to ensure it's initialized
  }
  
  // We'll update the block coordinate mappings when block ranges are populated
}

void StateManager::clearCaches() {
  sequenceCache.clear();
  // Removed gapSkipCache reference

  logging::debug("Cleared StateManager caches");
}

// Initialize DFS indices for traversal
void StateManager::initializeDfsIndices(panmanUtils::Tree *tree) {
  if (!tree || !tree->root) {
    logging::err("Invalid tree in initializeDfsIndices");
    return;
  }

  logging::info("Starting DFS index initialization with {} nodes in tree", tree->allNodes.size());
  logging::info("Root node: {}", tree->root->identifier);

  
  int64_t dfsIndex = 0;
  
  // First assign index to root node
  nodeDfsIndices[tree->root->identifier] = dfsIndex++;
  logging::info("Assigned DFS index 0 to root node {}", tree->root->identifier);
  
  // Use DFS traversal for the rest
  initializeDfsIndices(tree, tree->root, dfsIndex);
  
  // Check for nodes without indices and assign them
  logging::info("Initial DFS traversal assigned {} indices", dfsIndex);
  
  int missingCount = 0;
  for (const auto& [nodeId, node] : tree->allNodes) {
    if (nodeDfsIndices.find(nodeId) == nodeDfsIndices.end()) {
      nodeDfsIndices[nodeId] = dfsIndex++;
      missingCount++;
      
      if (missingCount <= 5) {
        logging::info("Assigned fallback DFS index {} to node {}", dfsIndex-1, nodeId);
      }
    }
  }
  
  if (missingCount > 0) {
    logging::info("Assigned fallback DFS indices to {} nodes not covered by standard traversal", missingCount);
  }
  
  logging::info("DFS indexing complete, assigned {} indices total", dfsIndex);
}

// Recursive helper for DFS index initialization
void StateManager::initializeDfsIndices(panmanUtils::Tree *tree,
                                        panmanUtils::Node *node,
                                        int64_t &dfsIndex) {
  if (!node) {
    logging::warn("Null node encountered in DFS traversal");
    return;
  }

  // Root node already assigned in the calling function
  if (node != tree->root) {
    // Assign current DFS index to this node
    nodeDfsIndices[node->identifier] = dfsIndex++;
    
    // Log first few assignments
    if (dfsIndex <= 5) {
      logging::info("Assigned DFS index {} to node {}", dfsIndex-1, node->identifier);
    }
  }

  // Recursively process children
  for (auto *child : node->children) {
    if (!child) {
      logging::warn("Null child encountered in node {}", node->identifier);
      continue;
    }
    initializeDfsIndices(tree, child, dfsIndex);
  }
}

// Get the DFS index for a node
int64_t StateManager::getDfsIndex(const std::string &nodeId) const {
  auto it = nodeDfsIndices.find(nodeId);
  if (it != nodeDfsIndices.end()) {
    return it->second;
  }
  return -1; // Invalid index
}

// Set the k-mer size for recomputation range calculations
void StateManager::setKmerSize(int16_t k) {
  if (k > 0) {
    kmerSize = k;
  }
}

// Get the current k-mer size
int16_t StateManager::getKmerSize() const { return kmerSize; }

// Get the total number of blocks
size_t StateManager::getNumBlocks() const { return numBlocks; }

// Set the total number of blocks
void StateManager::setNumBlocks(size_t blockCount) { 
  numBlocks = blockCount; 
  
  // Reserve space in block-related maps
  blockToSeeds.reserve(blockCount);
  
  logging::debug("Set number of blocks to {}", blockCount);
}

// Initialize a node's state
void StateManager::initializeNode(const std::string &nodeId) {
  // First check if node is already initialized (without locking)
  if (nodeStates.find(nodeId) != nodeStates.end()) {
    return;
  }
  
  // Track nodes being initialized to avoid recursion loops
  static thread_local std::unordered_set<std::string> nodesInProgress;
  if (nodesInProgress.find(nodeId) != nodesInProgress.end()) {
    logging::debug("Avoiding recursion: node {} is already being initialized", nodeId);
    return;
  }
  
  // Mark this node as being initialized
  nodesInProgress.insert(nodeId);
  
  // Use RAII guard to ensure we remove the node from inProgress on exit
  struct InitGuard {
    std::string id;
    InitGuard(const std::string& nodeId) : id(nodeId) {}
    ~InitGuard() { nodesInProgress.erase(id); }
  } initGuard(nodeId);

  // Check again with lock
  {
    std::lock_guard<std::mutex> lock(nodeStatesMutex);
    if (nodeStates.find(nodeId) != nodeStates.end()) {
      logging::debug("Node {} already initialized (checked with lock)", nodeId);
      return;
    }
    
    // Create a new node state
    nodeStates.emplace(nodeId, NodeState());
    logging::debug("Created new node state for node {}", nodeId);
  }
  
  // Ensure we have a valid root gap map (lazy initialization)
  if (!rootGapMap) {
    std::lock_guard<std::mutex> lock(nodeStatesMutex);
    if (!rootGapMap) {
      logging::info("Creating root gap map during initialization of node {}", nodeId);
      rootGapMap = std::make_shared<HierarchicalGapMap>();
    }
  }

  // Find parent node in hierarchy
  std::string parentId;
  {
    auto hierarchyIt = nodeHierarchy.find(nodeId);
    if (hierarchyIt != nodeHierarchy.end()) {
      parentId = hierarchyIt->second.parentId;
    }
  }

  // If this is not the root, propagate from parent
  if (!parentId.empty()) {
    logging::debug("Node {} has parent {}, initializing parent first", nodeId, parentId);
    
    // Make sure parent is initialized FIRST (without holding our lock)
    initializeNode(parentId);
    
    // Now propagate state from parent to child
    try {
      propagateState(parentId, nodeId);
    } catch (const std::exception& e) {
      logging::err("Error propagating from {} to {}: {}", parentId, nodeId, e.what());
    }
  } else {
    // Root node - set up root gap map
    std::lock_guard<std::mutex> lock(nodeStatesMutex);
    auto& nodeState = nodeStates[nodeId];
    
    if (!rootGapMap) {
      rootGapMap = std::make_shared<HierarchicalGapMap>();
    }
    
    nodeState.gapMap = rootGapMap;
    logging::debug("Assigned root gap map to root node {}", nodeId);
  }

  // Final verification and setup (with lock)
  {
    std::lock_guard<std::mutex> lock(nodeStatesMutex);
    
    // Ensure node has a valid gap map
    if (nodeStates[nodeId].gapMap == nullptr) {
      if (rootGapMap) {
        nodeStates[nodeId].gapMap = std::make_shared<HierarchicalGapMap>(rootGapMap);
      } else {
        nodeStates[nodeId].gapMap = std::make_shared<HierarchicalGapMap>();
      }
    }
    
    // Reset cache for this node
    resetNodeCache(nodeId);
  }
}

// Track which block is used by which node
void StateManager::trackBlockUsage(int32_t blockId, const std::string &nodeId) {
  if (blockId >= 0 && blockId < static_cast<int32_t>(numBlocks)) {
    // Ensure the block-to-seeds map has an entry for this block
    blockToSeeds[blockId]; // Create entry if doesn't exist
    std::cout << "Block " << blockId << " used by node " << nodeId << std::endl;
  }
}

// Get block status
bool StateManager::isBlockOn(std::string_view nodeId, int32_t blockId) const {
  // Find the node state directly with string_view
  std::lock_guard<std::mutex> lock(nodeStatesMutex);
  auto it = nodeStates.find(std::string(nodeId));
  if (it == nodeStates.end()) {
    return false;
  }
  return it->second.isBlockOn(blockId);
}

bool StateManager::isBlockInverted(std::string_view nodeId,
                                   int32_t blockId) const {
  if (auto stateIt = nodeStates.find(std::string(nodeId));
      stateIt != nodeStates.end()) {
    return stateIt->second.isBlockInverted(blockId);
  }
  return false;
}

// Set block status
void StateManager::setBlockOn(std::string_view nodeId, int32_t blockId, bool on) {
  std::lock_guard<std::mutex> lock(nodeStatesMutex);
  
  // Find the node state directly
  auto it = nodeStates.find(std::string(nodeId));
  if (it == nodeStates.end()) {
    return;
  }
  
  // Get reference to node state to modify it
  auto& nodeState = it->second;
  
  // Use NodeState's built-in method to change block status
  nodeState.setBlockOn(blockId, on);
  
  // Update tracking metadata
  if (on) {
    // Block was turned on, update node cache
    nodeBlockRangeCache.erase(std::string(nodeId));
  }
}

void StateManager::setBlockForward(std::string_view nodeId, int32_t blockId,
                                   bool forward) {
  auto &nodeState = nodeStates[std::string(nodeId)];
  nodeState.setBlockForward(blockId, forward);

  // Reset cache since block orientation changed
  resetNodeCache(nodeId);
}

void StateManager::setBlockInverted(std::string_view nodeId, int32_t blockId,
                                    bool inverted) {
  bool forward = !inverted;
  setBlockForward(nodeId, blockId, forward);
}

// Get character at position
char StateManager::getCharAtPosition(std::string_view nodeId, int32_t blockId,
                                     int32_t nucPos, int32_t gapPos) const {
    
    // Validate node exists
    auto stateIt = nodeStates.find(std::string(nodeId));
    if (stateIt == nodeStates.end()) {
        std::stringstream ss;
        ss << "Node not found: " << nodeId;
        logging::err("{}", ss.str());
        throw std::runtime_error(ss.str());
    }
    
    // Check if block is active for this node - return gap character if not
    bool isBlockActive = isBlockOn(nodeId, blockId);
    if (!isBlockActive) {
        return '-'; // Return gap character for inactive blocks - matches PanMAN
    }
    
    const auto& nodeState = stateIt->second;
    
    // Create position key for character lookup
    PositionKey key{blockId, nucPos, gapPos};
    
    // STEP 1: Check node's direct character data (most specific)
    auto charIt = nodeState.characterData.find(key);
    if (charIt != nodeState.characterData.end()) {
        return charIt->second;
    }
    
    // STEP 2: Check node's mutation index to find which node has this mutation
    auto mutIt = nodeState.mutationIndex.find(key);
    if (mutIt != nodeState.mutationIndex.end() && mutIt->second != std::string(nodeId)) {
        std::string sourceNodeId = mutIt->second;
        
        try {
            // Verify the source node exists
            auto sourceNodeIt = nodeStates.find(sourceNodeId);
            if (sourceNodeIt != nodeStates.end()) {
                // Check if block is active in the source node
                if (!isBlockOn(sourceNodeId, blockId)) {
                    return '-';
                }
                
                // Get character directly from source node's character data
                auto sourceCharIt = sourceNodeIt->second.characterData.find(key);
                if (sourceCharIt != sourceNodeIt->second.characterData.end()) {
                    char result = sourceCharIt->second;
                    
                    return result;
                }
            }
        } catch (const std::exception& e) {
            // Continue to global lookup on error
        }
    }
    
    // STEP 3: Check global mutation index as a fallback
    auto globalMutIt = mutationIndex.find(key);
    if (globalMutIt != mutationIndex.end() && globalMutIt->second != std::string(nodeId)) {
        std::string sourceNodeId = globalMutIt->second;
        
        try {
            // Verify the source node exists
            auto sourceNodeIt = nodeStates.find(sourceNodeId);
            if (sourceNodeIt != nodeStates.end()) {
                // Check if block is active in the source node
                if (!isBlockOn(sourceNodeId, blockId)) {
                    return '-';
                }
                
                // Get character directly from source node's character data
                auto sourceCharIt = sourceNodeIt->second.characterData.find(key);
                if (sourceCharIt != sourceNodeIt->second.characterData.end()) {
                    char result = sourceCharIt->second;
                    
                    return result;
                }
            }
        } catch (const std::exception& e) {
            // Continue to reference sequence on error
        }
    }
    
    // STEP 4: Fall back to reference sequence for this block
    auto seqIt = blockSequences.find(blockId);
    if (seqIt != blockSequences.end()) {
        const auto& sequence = seqIt->second;
        
        // Try using our pre-computed mapping first
        int32_t stringIndex = -1;
        
        // Get gap list length for this position
        size_t gapListLength = 0;
        try {
            gapListLength = getGapListLength(blockId, nucPos);
        } catch (const std::exception& e) {
          // Error logged inside getGapListLength if needed
        }
        
        // Check if this should be a valid lookup
        if (gapPos >= 0 && gapPos >= static_cast<int32_t>(gapListLength)) {
            // Gap position is out of bounds for the defined gap list
        }

        // Check if we have a direct mapping for this position
        auto indexIt = blockCoordToStringIndex.find(key);
        if (indexIt != blockCoordToStringIndex.end()) {
            // Use pre-computed mapping - this is in FORWARD ordering
            stringIndex = indexIt->second;
            
            // Get block inversion status for this node
            bool isInverted = isBlockInverted(nodeId, blockId);
            
            // For inverted blocks, convert the forward string index to inverted index
            if (isInverted) {
                int32_t seqLength = static_cast<int32_t>(sequence.length());
                if (stringIndex >= 0 && stringIndex < seqLength) {
                    // Invert the index: (length - 1) - forwardIndex
                    stringIndex = (seqLength - 1) - stringIndex;
                }
            }
        } else {
            logging::err("No mapping found for position: {}:{}:{}", blockId, nucPos, gapPos);
            throw std::runtime_error("No mapping found for position: " + std::to_string(blockId) + ":" + 
                                    std::to_string(nucPos) + ":" + std::to_string(gapPos));
        }
        
        // Verify the calculated index is valid
        if (stringIndex >= 0 && stringIndex < static_cast<int32_t>(sequence.length())) {
            char result = sequence[stringIndex];
            
            // Check if we need to complement the nucleotide for inverted blocks
            bool isInverted = isBlockInverted(nodeId, blockId);
            if (isInverted && isNonGapChar(result)) {
                // Complement the nucleotide (A<->T, C<->G)
                switch (result) {
                    case 'A': result = 'T'; break;
                    case 'T': result = 'A'; break;
                    case 'C': result = 'G'; break;
                    case 'G': result = 'C'; break;
                    case 'a': result = 't'; break;
                    case 't': result = 'a'; break;
                    case 'c': result = 'g'; break;
                    case 'g': result = 'c'; break;
                    // Other characters remain unchanged
                    default: break;
                }
            }

            return result;
        } else {
            logging::err("Invalid string index {} for sequence of length {}", stringIndex, sequence.length());
            throw std::runtime_error("Invalid string index: " + std::to_string(stringIndex) + 
                                  " for sequence of length " + std::to_string(sequence.length()));
        }
    }
    logging::err("Block {} sequence not found for node {}", blockId, nodeId);
    throw std::runtime_error("Block " + std::to_string(blockId) + " sequence not found for node " + std::string(nodeId));
}

// Set character at (blockId, nucPos, gapPos) in nodeId
bool StateManager::setCharAtPosition(std::string_view nodeId, int32_t blockId,
                                     int32_t nucPos, int32_t gapPos, char c) {
  try {
    auto &nodeState = getNodeState(std::string(nodeId));

   

    // Store the character
    PositionKey key{blockId, nucPos, gapPos};
    nodeState.characterData[key] = c;

    // Update the mutation index for this position
    updateMutationIndex(key, std::string(nodeId));


    return true;
  } catch (const std::exception& e) {
    logging::warn("Error in setCharAtPosition({}:{}:{}:{}): {}", 
                std::string(nodeId), blockId, nucPos, gapPos, e.what());
    return false;
  }
}

// Get range in scalar global coordinates for a block
CoordRange StateManager::getBlockRange(int32_t blockId) const {
  auto blockRangeIt = blockRanges.find(blockId);
  if (blockRangeIt != blockRanges.end()) {
    return blockRangeIt->second;
  }
  throw std::runtime_error("Block range not found for block ID: " + std::to_string(blockId));
}

// Set range in scalar global coordinates for a block
void StateManager::setBlockRange(int32_t blockId, const CoordRange &range) {
  // Store the block range
  blockRanges[blockId] = range;
  
  // Update global position cache
  if (blockId >= 0 && blockId < static_cast<int32_t>(globalPosCache.blockStartOffsets.size())) {
    globalPosCache.blockStartOffsets[blockId] = range.start;
  }

  // Log the block range setting
  logging::debug("Set block {} range to [{}, {}), size: {}", blockId,
                 range.start, range.end, range.end - range.start);
}

// Calculate recomputation range based on mutation type with efficient initial 4*k guess
std::optional<CoordRange>
StateManager::calculateRecompRange(std::string_view nodeId, int32_t blockId,
                                   int32_t pos, int32_t len, bool isBlockMutation,
                                   bool isBlockDeactivation) {
    if (kmerSize <= 0) kmerSize = getKmerSize();
    
    bool blockActive = isBlockOn(nodeId, blockId);
    
    auto blockRangeIt = blockRanges.find(blockId);
    if (blockRangeIt == blockRanges.end()) {
      throw std::runtime_error("Block " + std::to_string(blockId) + " not found in blockRanges");
    }
    
    auto blockRange = blockRangeIt->second;
    bool isInverted = blockActive && isBlockInverted(nodeId, blockId);
    
    CoordRange result = coordinates::RangeOperations::calculateRecompRange(
        blockRange, pos, len, kmerSize, isBlockMutation, isBlockDeactivation,
        numCoords, isInverted);
    
    // Constrain to block boundaries for initial range
    result.start = std::max(result.start, blockRange.start);
    result.end = std::min(result.end, blockRange.end);
    
    if (isBlockMutation || isBlockDeactivation) {
      auto sortedBlocks = getSortedActiveBlocks(nodeId);
      
      int currentPos = -1;
      for (size_t i = 0; i < sortedBlocks.size(); ++i) {
        if (sortedBlocks[i].first == blockId) {
          currentPos = static_cast<int>(i);
          break;
        }
      }
      
      if (currentPos >= 0 || isBlockDeactivation) {
        if (isBlockDeactivation) {
          result.start = std::max(result.start, blockRange.start);
          result.end = std::min(result.end, blockRange.end);
        }
      }
    }
    
    
    return result;
}

bool StateManager::applyBlockMutation(std::string_view nodeId, int32_t blockId,
                                      bool isInsertion, bool isInversion) {
  std::string strNodeId(nodeId); // Convert string_view to string for map keys

  // Get node state - ensure it exists or is initialized
  NodeState& nodeState = getNodeState(strNodeId);

  // Determine if this is a deactivation (turning block off) or activation/inversion change
  bool wasOn = isBlockOn(strNodeId, blockId);
  bool wasInverted = wasOn && isBlockInverted(strNodeId, blockId);
  bool isDeactivation = !isInsertion && wasOn;
  bool inversionChanged = false; // Initialize inversionChanged

  // Determine activation change
  bool activationChanged = isInsertion && !wasOn; // Turning block ON

  // Find block range - throw if not found
  auto blockRangeIt = blockRanges.find(blockId);
  if (blockRangeIt == blockRanges.end()) {
    throw std::runtime_error("Block " + std::to_string(blockId) +
                           " range not found for block mutation in node " + strNodeId);
  }
  auto blockRange = blockRangeIt->second;

  // Apply the block mutation state changes
  if (isInsertion) {
    nodeState.setBlockOn(blockId, true);
    nodeState.setBlockForward(blockId, !isInversion); // setBlockForward takes true for forward
    // Check if inversion status actually changed upon activation
    inversionChanged = (wasInverted != isInversion);
  } else {
    // Handle actual deletion (not inversion toggle)
    if (!isInversion) {
        nodeState.setBlockOn(blockId, false);
        nodeState.setBlockForward(blockId, true); // Reset strand on deletion
    } else {
        // Handle inversion toggle (mutInfo=0, inversion=true)
        bool newInvertedState = !wasInverted;
        nodeState.setBlockForward(blockId, !newInvertedState);
        inversionChanged = true; // Inversion always changes if toggled
    }
  }

  // --- Recomputation Range Logic ---
  auto &recompRanges = nodeState.recompRanges;
  
  // Get k-mer size for proper expansion
  int k = getKmerSize();
  if (k <= 0) k = 32; // Default k-mer size if not set
  
  // Helper to add a range to the node's recomp ranges
  auto addRange = [&recompRanges](const coordinates::CoordRange& range) {
    if (range.start < range.end) {
      recompRanges.push_back(range);
    }
  };
  
  // CASE 1: Block activation or inversion change
  if (activationChanged || inversionChanged) {
    // For off->on or inversion change: include the entire block
    CoordRange range = {blockRange.start, blockRange.end};
    std::cout << "NODE_MUT_DETAIL: BlockMut " << blockId << " adding FullBlockRecompRange=[" << range.start << ", " << range.end << ")" << std::endl;
    addRange(range);
  }
  // CASE 2: Block deactivation
  else if (isDeactivation) {
    // For on->off: add k positions at start and end of block
    
    // Start boundary range (first k positions)
       CoordRange startBoundary = {blockRange.start, std::min(blockRange.start + k, blockRange.end)};
       std::cout << "NODE_MUT_DETAIL: BlockMut " << blockId << " adding StartBoundaryRecompRange=[" << startBoundary.start << ", " << startBoundary.end << ")" << std::endl;
       addRange(startBoundary);
    
    // End boundary range (last k positions)
    CoordRange endBoundary = {std::max(blockRange.end - k, blockRange.start), blockRange.end};
       std::cout << "NODE_MUT_DETAIL: BlockMut " << blockId << " adding EndBoundaryRecompRange=[" << endBoundary.start << ", " << endBoundary.end << ")" << std::endl;
       addRange(endBoundary);
 
    // Also need to update ranges for adjacent blocks (but don't cross block boundaries)
       int32_t prevBlockId = -1;
       int32_t nextBlockId = -1;

    // Find closest active preceding block
       for (int32_t checkId = blockId - 1; checkId >= 0; --checkId) {
           if (isBlockOn(strNodeId, checkId)) {
               prevBlockId = checkId;
               break;
           }
       }

    // Find closest active succeeding block
       for (int32_t checkId = blockId + 1; checkId < static_cast<int32_t>(getNumBlocks()); ++checkId) {
          if (isBlockOn(strNodeId, checkId)) {
              nextBlockId = checkId;
              break;
          }
       }

       // Add boundary range for the end of the previous active block
       if (prevBlockId != -1) {
           try {
               auto prevBlockRange = getBlockRange(prevBlockId);
        bool isPrevInverted = isBlockInverted(strNodeId, prevBlockId);
        
        // For previous block, add the end boundary (or start if inverted)
        CoordRange prevEndBoundary;
        if (isPrevInverted) {
          // For inverted blocks, add k positions at the start
          prevEndBoundary = {prevBlockRange.start, std::min(prevBlockRange.start + k, prevBlockRange.end)};
        } else {
          // For normal blocks, add k positions at the end
          prevEndBoundary = {std::max(prevBlockRange.end - k, prevBlockRange.start), prevBlockRange.end};
        }
        
        std::cout << "NODE_MUT_DETAIL: BlockMut " << blockId << " adding PrevAdjacentRecompRange=[" 
                 << prevEndBoundary.start << ", " << prevEndBoundary.end << ") (from Block " << prevBlockId << ")" << std::endl;
               addRange(prevEndBoundary);
           } catch (const std::exception& e) {
        std::cout << "NODE_MUT_DETAIL: WARN - Could not get range for prev block " << prevBlockId << ": " << e.what() << std::endl;
           }
       }

       // Add boundary range for the start of the next active block
       if (nextBlockId != -1) {
            try {
               auto nextBlockRange = getBlockRange(nextBlockId);
        bool isNextInverted = isBlockInverted(strNodeId, nextBlockId);
        
        // For next block, add the start boundary (or end if inverted)
        CoordRange nextStartBoundary;
        if (isNextInverted) {
          // For inverted blocks, add k positions at the end
          nextStartBoundary = {std::max(nextBlockRange.end - k, nextBlockRange.start), nextBlockRange.end};
        } else {
          // For normal blocks, add k positions at the start
          nextStartBoundary = {nextBlockRange.start, std::min(nextBlockRange.start + k, nextBlockRange.end)};
        }
        
        std::cout << "NODE_MUT_DETAIL: BlockMut " << blockId << " adding NextAdjacentRecompRange=[" 
                 << nextStartBoundary.start << ", " << nextStartBoundary.end << ") (from Block " << nextBlockId << ")" << std::endl;
               addRange(nextStartBoundary);
           } catch (const std::exception& e) {
               std::cout << "NODE_MUT_DETAIL: WARN - Could not get range for next block " << nextBlockId << ": " << e.what() << std::endl;
           }
       }
  }

  // Reset node-specific caches because block state changed
  resetNodeCache(nodeId);

  return true;
}

// Helper method to propagate state
void StateManager::propagateState(const std::string &fromNode, const std::string &toNode) {
  if (nodeStates.find(fromNode) == nodeStates.end()) {
    logging::err("Parent node {} not found for propagation", fromNode);
    throw std::runtime_error("Parent node " + fromNode + " not found for propagation");
  }
  
  // Get existing child state or create a new one
  auto childStateIt = nodeStates.find(toNode);
  if (childStateIt == nodeStates.end()) {
    // Create new child state
    nodeStates.emplace(toNode, NodeState());
  }
  
  // Get parent and child node states
  auto &parentState = getNodeState(fromNode);
  auto &childState = getNodeState(toNode);

  // Set parent relationship for hierarchical lookups
  childState.parentId = fromNode;

  // Track problematic positions for debugging
  std::vector<int64_t> problematicPositions;
  std::unordered_map<int32_t, bool> blockPositionChecked;
  
  // Log initial state before propagation
  logging::info("BLOCK_PROPAGATION: Parent '{}' to child '{}' starting, parent has {} active blocks",
             fromNode, toNode, parentState.activeBlocks.size());
  
  // Ensure child node inherits parent's block state
  for (int32_t blockId : parentState.activeBlocks) {
    // Check if parent has this block
    if (parentState.isBlockOn(blockId)) {
      
      // Propagate block activation status
      childState.setBlockOn(blockId, true);
      
      // Also propagate block orientation (forward/inverted)
      bool isInverted = parentState.isBlockInverted(blockId);
      childState.setBlockForward(blockId, !isInverted);
    }
  }

  // Set up hierarchical gap map - ALWAYS ensure a gap map is set
  if (parentState.gapMap) {
    // Best case - inherit from parent's gap map
    childState.gapMap = std::make_shared<HierarchicalGapMap>(parentState.gapMap);
  } else {
    // Fallback to root gap map if parent has none
    // Happens when, e.g. node_1, the first node, has no mutations => node_2 will inherit the root gap map
    childState.gapMap = std::make_shared<HierarchicalGapMap>(rootGapMap);
  }
  
  // FIXED: Copy the mutation index entries from parent to child
  // These tell the child which ancestor nodes have mutations at each position
  for (const auto& [posKey, ancestorId] : parentState.mutationIndex) {
    // Simply copy the entry - this tells the child which ancestor has the mutation
    childState.mutationIndex[posKey] = ancestorId;
  }
  
  // Reset cache for child node
  resetNodeCache(toNode);
}

void StateManager::applyNucleotideMutation(const std::string& nodeId, NodeState& nodeState,
                                         int32_t blockId, int32_t nucPos, int32_t gapPos, char value) {
  // Add validation to check if the position is valid for the block
  auto blockSeqIt = blockSequences.find(blockId);
  if (blockSeqIt == blockSequences.end()) {
    throw std::runtime_error("Block " + std::to_string(blockId) + " not found for mutation in node " + nodeId);
  }
  
  // Check if position exceeds block sequence length
  if (nucPos >= static_cast<int32_t>(blockSeqIt->second.length())) {
    throw std::runtime_error("Mutation position " + std::to_string(nucPos) + 
                           " exceeds block " + std::to_string(blockId) + 
                           " length (" + std::to_string(blockSeqIt->second.length()) + ")");
  }

  // Get gap list length for this position
  size_t gapListLength = 0;
  try {
    gapListLength = getGapListLength(blockId, nucPos);
  } catch (const std::exception& e) {
    std::cout << ", ERROR: Failed to get gap list length: " << e.what();
  }
  
  PositionKey key{blockId, nucPos, gapPos};
  nodeState.characterData[key] = value;
  
  updateMutationIndex(key, nodeId);
  
}

// Helper to set gap list length
// This sets the number of positions in the gap list for a specific nucleotide
// Note: This is NOT related to gap characters ('-') in the sequence tracked by
// gap maps. The gap list contains optional nucleotides that:
// - In normal blocks: come BEFORE the main nucleotide in traversal order
// - In inverted blocks: come AFTER the main nucleotide in traversal order (and
// in reverse sequence)
void StateManager::setGapListLength(int32_t blockId,
                                    int32_t nucPos, size_t length) {
  // Create key for lookup
  GapListKey key = GapListKey::create(blockId, nucPos);

  // Set in global gap list lengths (legacy)
  globalGapListLengths[key] = length;

  // Set in new optimized array if initialized
  if (blockId >= 0 && static_cast<size_t>(blockId) < gapListLengthArray.size()) {
    auto& blockArray = gapListLengthArray[blockId];
    if (nucPos >= 0 && static_cast<size_t>(nucPos) < blockArray.size()) {
      blockArray[nucPos] = length;
      
      // Also update thread-local cache if entry exists
      bool found = false;
      gapListLengthCaches.local().get(blockId, nucPos, found);
      if (found) {
        gapListLengthCaches.local().put(blockId, nucPos, length);
      }
    }
  }
}

// Get gap list length for a position - optimized version
size_t StateManager::getGapListLength(int32_t blockId, int32_t nucPos) const {
  // Check thread-local cache first (fastest path)
  bool found = false;
  size_t cachedLength = gapListLengthCaches.local().get(blockId, nucPos, found);
  if (found) {
    return cachedLength;
  }
  
  // Next check optimized array (second fastest)
  if (blockId >= 0 && static_cast<size_t>(blockId) < gapListLengthArray.size()) {
    const auto& blockArray = gapListLengthArray[blockId];
    if (nucPos >= 0 && static_cast<size_t>(nucPos) < blockArray.size()) {
      size_t length = blockArray[nucPos];
      
      // Update cache for future lookups
      gapListLengthCaches.local().put(blockId, nucPos, length);
      
      return length;
    }
  }
  
  // Fallback to legacy hash table lookup (slowest path)
  GapListKey key = GapListKey::create(blockId, nucPos);
  auto it = globalGapListLengths.find(key);
  if (it != globalGapListLengths.end()) {
    size_t length = it->second;
    
    // Update cache for future lookups
    gapListLengthCaches.local().put(blockId, nucPos, length);
    
    return length;
  }
  
  // Default to 0 length if no gap list exists
  return 0;
}

// Initialize the optimized gap list length array
void StateManager::initializeGapListLengthArray(size_t numBlocks, size_t maxNucPos) {
  // Initialize with appropriate size
  gapListLengthArray.clear();
  gapListLengthArray.resize(numBlocks);
  
  // Find the maximum nucPos for each block directly from gap lists
  std::vector<size_t> maxNucPosPerBlock(numBlocks, 0);
  
  // First pass - determine the maximum nucPos for each block from gap lists
  for (const auto& [key, length] : globalGapListLengths) {
    if (key.blockId >= 0 && static_cast<size_t>(key.blockId) < numBlocks) {
      maxNucPosPerBlock[key.blockId] = std::max(maxNucPosPerBlock[key.blockId], 
                                              static_cast<size_t>(key.nucPos + 1));
    }
  }
  
  // Calculate the overall max nucPos for fallback purposes
  size_t overallMaxNucPos = 0;
  for (size_t blockId = 0; blockId < numBlocks; ++blockId) {
    overallMaxNucPos = std::max(overallMaxNucPos, maxNucPosPerBlock[blockId]);
  }
  
  // Use the passed maxNucPos value if it's valid and larger than what we calculated
  if (maxNucPos > overallMaxNucPos) {
    overallMaxNucPos = maxNucPos;
  }
  
  // If we still don't have a reasonable nucPos value, use a default
  if (overallMaxNucPos == 0) {
    overallMaxNucPos = 10000; // Conservative default
  }
  
  // Allocate and initialize arrays
  for (size_t blockId = 0; blockId < numBlocks; ++blockId) {
    size_t blockMaxPos = maxNucPosPerBlock[blockId];
    if (blockMaxPos == 0) {
      // Use the overall max as a reasonable fallback if this block has no gap lists
      blockMaxPos = overallMaxNucPos;
    }
    
    // Ensure minimum reasonable size
    blockMaxPos = std::max(blockMaxPos, static_cast<size_t>(1000));
    
    gapListLengthArray[blockId].resize(blockMaxPos, 0);
  }
  
  // Populate from existing hash table
  for (const auto& [key, length] : globalGapListLengths) {
    if (key.blockId >= 0 && static_cast<size_t>(key.blockId) < numBlocks) {
      auto& blockArray = gapListLengthArray[key.blockId];
      if (key.nucPos >= 0 && static_cast<size_t>(key.nucPos) < blockArray.size()) {
        blockArray[key.nucPos] = length;
      }
    }
  }
  
  logging::info("Initialized gap list length array with {} blocks and max nucPos {}", 
                numBlocks, overallMaxNucPos);
}

// Helper method to get a sorted list of active blocks
std::vector<std::pair<int32_t, coordinates::CoordRange>>
StateManager::getSortedActiveBlocks(std::string_view nodeId) const {

  std::vector<std::pair<int32_t, coordinates::CoordRange>> sortedBlocks;
  const auto &state = getNodeState(std::string(nodeId));

  for (int32_t blockId : state.activeBlocks) {
    auto rangeIt = blockRanges.find(blockId);
    if (rangeIt != blockRanges.end()) {
      sortedBlocks.emplace_back(blockId, rangeIt->second);
    }
  }

  // Sort blocks by start position
  std::sort(sortedBlocks.begin(), sortedBlocks.end(),
            [](const auto &a, const auto &b) {
              return a.second.start < b.second.start;
            });

  // We can't cache the result in a const method
  // nodeBlockRangeCache[std::string(nodeId)] = sortedBlocks;

  return sortedBlocks;
}

// Helper to get node gap map
std::shared_ptr<HierarchicalGapMap>
StateManager::getNodeGapMap(std::string_view nodeId) const {
  auto it = nodeStates.find(std::string(nodeId));
  if (it != nodeStates.end() && it->second.gapMap) {
    return it->second.gapMap;
  }

  logging::warn("No gap map found for node {}, returning root gap map", nodeId);
  
  // Ensure rootGapMap is valid before returning
  if (!rootGapMap) {
    logging::err("Root gap map is null when trying to get gap map for node {}", nodeId);
    // This is a const method, so we can't modify rootGapMap directly
    // Return a new temporary gap map as a last resort
    return std::make_shared<HierarchicalGapMap>();
  }
  
  // Fallback to root gap map if node not found
  return rootGapMap;
}

// Get active blocks for a node
const std::unordered_set<int32_t> &
StateManager::getActiveBlocks(std::string_view nodeId) const {
  static const std::unordered_set<int32_t> emptySet;
  
  auto it = nodeStates.find(std::string(nodeId));
  if (it != nodeStates.end()) {
    return it->second.activeBlocks;
  }
  
  // Return empty set if node not found
  return emptySet;
}


// Helper method to merge a new range with existing ranges
void StateManager::mergeRangeWithExisting(std::vector<CoordRange> &ranges,
                                          const CoordRange &newRange) {
  // Use the consolidated RangeOperations function
  coordinates::RangeOperations::mergeRangeWithExisting(ranges, newRange);
}


void StateManager::backtrackNode(const std::string &nodeId) {
  // Remove the node state to free memory
  nodeStates.erase(nodeId);


  // Clear sequence cache for this node
  auto it = sequenceCache.begin();
  while (it != sequenceCache.end()) {
    // Extract node ID from cache key (format is "nodeId:start:end:skipGaps")
    size_t firstColon = it->first.find(':');
    if (firstColon != std::string::npos) {
      std::string cacheNodeId = it->first.substr(0, firstColon);
      if (cacheNodeId == nodeId) {
        it = sequenceCache.erase(it);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }

  logging::debug("Backtracked node {}", nodeId);
}

// Get node state reference - use structured bindings for clarity
NodeState &StateManager::getNodeState(const std::string &nodeId) {
  if (auto it = nodeStates.find(nodeId); it != nodeStates.end()) {
    return it->second;
  }
  // Call initializeNode if the node doesn't exist
  initializeNode(nodeId);
  
  auto newIt = nodeStates.find(nodeId);
  if (newIt == nodeStates.end()) {
    logging::err("Failed to initialize node state for {}", nodeId);
    throw std::runtime_error("Failed to initialize node state for " + nodeId);
  }
  return newIt->second;
}

// Method to get const reference to a node's state
const NodeState &StateManager::getNodeState(const std::string &nodeId) const {
  auto it = nodeStates.find(nodeId);
  if (it == nodeStates.end()) {
    std::stringstream ss;
    ss << "Node state not found for node: " << nodeId
       << ". Total states available: " << nodeStates.size()
       << ". This error indicates that the node was not properly initialized.";
    logging::err("{}", ss.str());
    throw std::runtime_error(ss.str());
  }
  return it->second;
}

// Seed management methods (from TraversalGlobalState)

// Optimized seed management using hierarchical storage
std::optional<seeding::seed_t> StateManager::getSeedAtPosition(int64_t pos) const {
  if (pos < 0) {
    return std::nullopt;
  }

  // First try to get from hierarchical store for the current node
  // This would require context of which node we're working with

  // Fallback to global position seeds
  if (pos < static_cast<int64_t>(positionSeeds.size())) {
    return positionSeeds[pos];
  }

  return std::nullopt;
}

void StateManager::setSeedAtPosition(int64_t pos, const seeding::seed_t &seed) {
  if (pos < 0) {
    return;
  }

  // Resize if needed
  if (pos >= static_cast<int64_t>(positionSeeds.size())) {
    positionSeeds.resize(pos + 1);
  }

  // Update global position seeds
  positionSeeds[pos].emplace(seed);
}

// Efficiently process seeds in a range to avoid redundant sequence extraction
void StateManager::processSeedsInRange(
    std::string_view nodeId, const CoordRange &range, int k, int s,
    const std::function<void(int64_t startPos, int64_t endPos, std::optional<seeding::seed_t> &seed, std::string_view kmerView)> &processFn) {

  // Get active blocks for this node for block-aware processing
  std::vector<std::pair<int32_t, CoordRange>> relevantBlocks;
  auto sortedBlocks = getSortedActiveBlocks(nodeId);
  
  // Find blocks that overlap with our target range
  for (const auto& [blockId, blockRange] : sortedBlocks) {
    // Check if this block overlaps with our range
    if (blockRange.end <= range.start || blockRange.start >= range.end) {
      continue; // No overlap
    }
    
    // Calculate overlap between block and target range
    CoordRange overlapRange = {
      std::max(blockRange.start, range.start),
      std::min(blockRange.end, range.end)
    };
    
    relevantBlocks.emplace_back(blockId, overlapRange);
  }
  
  // If no relevant blocks, nothing to process
  if (relevantBlocks.empty()) {
    return;
  }
  
  logging::debug("Processing seeds in range [{}, {}) - found {} overlapping blocks", 
               range.start, range.end, relevantBlocks.size());
  
  // Process each block separately for better cache locality
  for (const auto& [blockId, blockRange] : relevantBlocks) {
    // Check cache first before extracting sequence
    std::string cacheKey = std::string(nodeId) + ":block:" + std::to_string(blockId);
    
    // Extract sequence for this block range
    std::pair<std::string, std::vector<int64_t>> result;
    try {
      // Extract with skipGaps=true for k-mer processing
      // This correctly uses node-specific gap map handling through extractSequence
      logging::info("EXTRACT_SEQ_CALL (processSeedsInRange): Block Processing - Node {} BlockId {} Range [{},{})", 
                   nodeId, blockId, blockRange.start, blockRange.end);
      // Log caller before extractSequence
      LOG_CALLER("extractSequence", nodeId, blockRange.start, blockRange.end);
      result = extractSequence(nodeId, blockRange, true);
    } catch (const std::exception& e) {
      logging::err("Error extracting sequence for block {} range [{}, {}): {}", 
                  blockId, blockRange.start, blockRange.end, e.what());
      throw std::runtime_error(e.what());
    }
    
    auto& sequence = result.first;
    auto& positions = result.second;
    
    // Ensure we have enough sequence for at least one k-mer
    if (sequence.length() < k) {
      continue;
    }
    
    // Track block-local start index to avoid rescanning from beginning each time
    size_t startIndex = 0;
    
    // Find first position inside our target range
    while (startIndex < positions.size() && positions[startIndex] < range.start) {
      startIndex++;
    }
    
    // Process each k-mer within this block (sequence is already gapless)
    for (size_t i = startIndex; i + k <= sequence.length(); ++i) {
      // Get global coordinates for start and end of this k-mer
      int64_t startPos = positions[i];
      
      // Make sure we're still in our target range
      if (startPos >= range.end) {
        break;
      }
      
      // Ensure we have enough positions for the entire k-mer
      if (i + k - 1 >= positions.size()) {
        break;
      }
      
      int64_t endPos = positions[i + k - 1];
      
      // Check if this position is already processed
      auto existingSeed = getSeedAtPosition(startPos);
      
      // Create view of the k-mer for efficient processing
      std::string_view kmerView(sequence.c_str() + i, k);
      
      // Call the processing function
      processFn(startPos, endPos, existingSeed, kmerView);
    }
  }
}

// Helper to check if a position is in a gap
bool StateManager::isGapPosition(std::shared_ptr<HierarchicalGapMap> nodeGapMap, int64_t pos) const {
  // Validate the provided position
  if (pos < 0) {
    logging::warn("Negative position {} provided to isGapPosition", pos);
    return false;
  }
  
  // Null check on gap map
  if (!nodeGapMap) {
    return false;
  }
  
  // Directly check if position is in a gap using the provided gap map
  try {
    return nodeGapMap->isGap(pos);
  } catch (const std::exception& e) {
    // Enhance error message with coordinate system context
    logging::debug("Failed to check gap status for position {}: {}", pos, e.what());
    return false;
  }
}

// Centralized gap map update function that handles all types of updates
void StateManager::consolidatedGapUpdate(const std::string& nodeId, 
                                      const std::vector<coordinates::GapUpdate>& updates) {
  if (updates.empty()) {
    return;
  }
  
  // Use lock for thread safety
  std::lock_guard<std::mutex> lock(nodeStatesMutex);
  
  auto nodeStateIt = nodeStates.find(nodeId);
  if (nodeStateIt == nodeStates.end()) {
    logging::err("Node {} not found for gap update", nodeId);
    return;
  }
  
  auto& nodeState = nodeStateIt->second;
  if (!nodeState.gapMap) {
    logging::err("Node {} has no gap map for update", nodeId);
    return;
  }
  
  // Convert to gap_map::GapUpdate format directly
  std::vector<gap_map::GapUpdate> gapMapUpdates;
  gapMapUpdates.reserve(updates.size());
  
  for (const auto& update : updates) {
    // Create gap_map::GapUpdate (isRemoval, {start, end})
    gap_map::GapUpdate gapMapUpdate(!update.isGapAddition, 
                                  {update.pos, update.pos + update.length - 1});
    gapMapUpdates.push_back(gapMapUpdate);
  }
  
  // Apply updates to node's gap map
  size_t updateCount = nodeState.gapMap->applyUpdates(gapMapUpdates);
  logging::debug("Applied {} gap updates to node {}", updateCount, nodeId);
  
  // Get all direct children (no need for recursive updates since we use hierarchical gap maps)
  std::vector<std::string> directChildren;
  for (const auto& [childId, childState] : nodeStates) {
    if (childState.parentId == nodeId) {
      directChildren.push_back(childId);
    }
  }
  
  // Reset cache for the current node
  resetNodeCache(nodeId);
  
  // Update only direct children's gap maps - don't re-apply updates recursively
  for (const auto& childId : directChildren) {
    auto childStateIt = nodeStates.find(childId);
    if (childStateIt != nodeStates.end() && childStateIt->second.gapMap) {
      // Update the child's gap map to inherit from parent's changed map
      childStateIt->second.gapMap = std::make_shared<HierarchicalGapMap>(nodeState.gapMap);
      
      // Reset the child's cache
      resetNodeCache(childId);
    }
  }
}

// Initialize the state manager with a tree
void StateManager::initialize(panmanUtils::Tree *tree) {
  if (!tree || !tree->root) {
    logging::err("Invalid tree in initialize");
    return;
  }
  
  logging::info("Initializing state manager with tree");
  
  // Set up the node hierarchy first
  initializeNodeHierarchy(tree, tree->root);
  
  // Initialize DFS indices for traversal
  initializeDfsIndices(tree);
  
  // Initialize mutation index from existing character data
  initializeMutationIndex();
  
  // Initialize optimized gap list length array 
  // We don't need to pre-calculate maxNucPos - the function will determine it
  initializeGapListLengthArray(numBlocks, 0);
  
  // Ensure all coordinate mappings are properly flushed
  flushCoordinateMappings();
  
  // Initialize active block boundaries for all nodes
  for (const auto& [nodeId, _] : tree->allNodes) {
    initializeBlockBoundaries(nodeId);
  }
  
  // Initialize global position cache for fast lookups 
  if (!globalPosCacheInitialized && numCoords > 0) {
    globalPosCache.posToBlockId.resize(numCoords, -1);
    globalPosCache.blockStartOffsets.resize(numBlocks, -1);
    globalPosCache.maxPosition = numCoords;
    
    // Fill in the cache with block IDs
    for (const auto& [blockId, range] : blockRanges) {
      globalPosCache.blockStartOffsets[blockId] = range.start;
      for (int64_t pos = range.start; pos < range.end; pos++) {
        if (pos >= 0 && pos < numCoords) {
          globalPosCache.posToBlockId[pos] = blockId;
        }
      }
    }
    
    globalPosCacheInitialized = true;
    logging::info("Initialized global position cache with {} coordinates", numCoords);
  }
  
  logging::info("StateManager initialization complete");
}

// Helper to initialize block boundaries
void StateManager::initializeBlockBoundaries(const std::string& nodeId) {
  auto& nodeState = getNodeState(nodeId);
  
  // Clear current boundaries
  nodeState.firstActiveBlockPos.reset();
  nodeState.lastActiveBlockPos.reset();
  
  // Calculate and set new boundaries if there are active blocks
  if (!nodeState.activeBlocks.empty()) {
    for (int32_t blockId : nodeState.activeBlocks) {
      try {
        auto blockRange = getBlockRange(blockId);
        
        // Update first active block position
        if (!nodeState.firstActiveBlockPos.has_value() || 
            blockRange.start < nodeState.firstActiveBlockPos.value()) {
          nodeState.firstActiveBlockPos = blockRange.start;
        }
        
        // Update last active block position
        if (!nodeState.lastActiveBlockPos.has_value() || 
            blockRange.end > nodeState.lastActiveBlockPos.value()) {
          nodeState.lastActiveBlockPos = blockRange.end;
        }
      } catch (const std::exception& e) {
        // Log but continue with other blocks
        logging::warn("Could not get range for block {} in node {}: {}", 
                    blockId, nodeId, e.what());
      }
    }
    
    logging::debug("Initialized block boundaries for node {}: first={}, last={}", 
                 nodeId, nodeState.firstActiveBlockPos.value_or(-1), 
                 nodeState.lastActiveBlockPos.value_or(-1));
  } else {
    logging::debug("No active blocks for node {}, no boundaries set", nodeId);
  }
}

// Helper to check if a character is a non-gap character
bool StateManager::isNonGapChar(char c) const {
    return c != '-' && c != 'x';
}

/**
 * @brief Set the sequence for a specific block
 * @param blockId The ID of the block
 * @param sequence The sequence to set for the block
 */
void StateManager::setBlockSequence(int32_t blockId, const std::string& sequence) {
    blockSequences[blockId] = sequence;
}

// Method to initialize all coordinate mappings for a block
void StateManager::initializeBlockMappings(int32_t blockId) {
  static std::mutex initMutex;
  static std::unordered_set<int32_t> initializedBlocks;
  
  // Check if this block has already been initialized
  {
    std::lock_guard<std::mutex> lock(initMutex);
    if (initializedBlocks.find(blockId) != initializedBlocks.end()) {
      // Already initialized
      return;
    }
  }
  
  // Check if the block sequence exists
  auto seqIt = blockSequences.find(blockId);
  if (seqIt == blockSequences.end()) {
    logging::warn("Cannot initialize mappings for block {}: sequence not found", blockId);
    return;
  }
  
  // Get block range
  auto rangeIt = blockRanges.find(blockId);
  if (rangeIt == blockRanges.end()) {
    logging::warn("Cannot initialize mappings for block {}: range not found", blockId);
    return;
  }
  
  const auto& sequence = seqIt->second;
  const auto& range = rangeIt->second;
  
  logging::debug("Initializing coordinate mappings for block {} with sequence length {}", 
                blockId, sequence.length());
  
  // Calculate all mappings in one pass for efficiency
  std::vector<std::tuple<int32_t, int32_t, int64_t>> mappings;
  mappings.reserve(sequence.length() * 500); // Estimate including some gaps
  
  int64_t pos = range.start;
  int32_t strIndex = 0;
  
  // Process each nucleotide position
  for (int32_t nucPos = 0; nucPos < static_cast<int32_t>(sequence.length()); nucPos++) {
    // First register gap positions for this nucPos
    size_t gapListLen = getGapListLength(blockId, nucPos);
    
    // Gap positions come before main position in the 3D coordinate system
    for (int32_t gapPos = static_cast<int32_t>(gapListLen) - 1; gapPos >= 0; gapPos--) {
      // Skip registration for positions that would be negative
      if (pos + gapPos >= 0) {
        mappings.emplace_back(nucPos, gapPos, pos + gapPos);
      }
    }
    
    // Register main nucleotide position last (follows proper coordinate order)
    mappings.emplace_back(nucPos, -1, pos + gapListLen);
    
    // Move to next position group
    pos += gapListLen + 1;
  }
  
  // Now batch register all the mappings at once
  for (const auto& [nPos, gPos, globalPos] : mappings) {
    registerCoordinateMapping(blockId, nPos, gPos, globalPos);
  }
  
  // Mark as initialized
  {
    std::lock_guard<std::mutex> lock(initMutex);
    initializedBlocks.insert(blockId);
  }
  
  logging::debug("Initialized coordinate mappings for block {}", blockId);
  
  // Ensure all changes are flushed
  flushCoordinateMappings();
}

// Method to map 3D coordinates to global position
int64_t StateManager::mapToGlobalCoordinate(int32_t blockId, int32_t nucPos, int32_t gapPos, bool inverted) const {
  // Create position key for precise lookup
  PositionKey key = PositionKey::create(blockId, nucPos, gapPos);
  
  // Direct hash lookup - fastest path
  auto it = blockCoordToGlobalPos.find(key);
  if (it != blockCoordToGlobalPos.end()) {
    int64_t globalPos = it->second;
    
    // If the caller indicates this is an inverted coordinate, transform it
    if (inverted) {
      auto rangeIt = blockRanges.find(blockId);
      if (rangeIt != blockRanges.end()) {
        const auto& range = rangeIt->second;
        int64_t blockStart = range.start;
        int64_t blockEnd = range.end;
        
        // Mirror the position around the block's center
        // For a block with range [100, 200), pos 120 becomes 180
        globalPos = blockStart + (blockEnd - blockStart) - (globalPos - blockStart) - 1;
      }
    }
    
    return globalPos;
  }
  
  // If no direct mapping found, compute it using block range and gap list length
  auto rangeIt = blockRanges.find(blockId);
  if (rangeIt != blockRanges.end()) {
    const auto& range = rangeIt->second;
    
    // Get gap list length
    size_t gapListLength = getGapListLength(blockId, nucPos);
    
    // Compute global position in coordinate space
    int64_t globalPos = range.start + nucPos + gapListLength;
    
    // Adjust for specific gap position if provided
    if (gapPos >= 0 && gapPos < static_cast<int32_t>(gapListLength)) {
      globalPos = range.start + nucPos + gapPos;
    }
    
    // If inverted, mirror around the block's center
    if (inverted) {
      int64_t blockStart = range.start;
      int64_t blockEnd = range.end;
      globalPos = blockStart + (blockEnd - blockStart) - (globalPos - blockStart) - 1;
    }
    
    return globalPos;
  }
  
  // Fallback for error case
  return -1;
}

// Register mapping from block coordinates to global position
void StateManager::registerCoordinateMapping(int32_t blockId, int32_t nucPos, int32_t gapPos, int64_t globalPos) {
  // Create position key
  PositionKey key = PositionKey::create(blockId, nucPos, gapPos);
  
  // Store mapping to global position
  blockCoordToGlobalPos[key] = globalPos;
  
  // Also compute and store string index if needed
  if (nucPos >= 0) {
    int32_t stringIndex = nucPos;
    blockCoordToStringIndex[key] = stringIndex;
  }
}

// Flush coordinate mappings to ensure all are properly stored
void StateManager::flushCoordinateMappings() {
  // No specific flush operation needed in current implementation
  // This method exists for compatibility with code that calls it
  logging::debug("Flushed coordinate mappings");
}

// Reset cache for a specific node
void StateManager::resetNodeCache(std::basic_string_view<char, std::char_traits<char>> nodeId) {
  // Clear cache entries related to this node
  nodeBlockRangeCache.erase(std::string(nodeId));
  
  // Clear sequence cache entries for this node
  auto it = sequenceCache.begin();
  while (it != sequenceCache.end()) {
    // Extract node ID from cache key (format is "nodeId:start:end:skipGaps")
    size_t firstColon = it->first.find(':');
    if (firstColon != std::string::npos) {
      std::string cacheNodeId = it->first.substr(0, firstColon);
      if (cacheNodeId == std::string(nodeId)) {
        it = sequenceCache.erase(it);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }
  
  logging::debug("Reset cache for node {}", nodeId);
}

// Update the mutation index for a position
void StateManager::updateMutationIndex(const PositionKey& key, const std::string& nodeId) {
  // Add/update entry in global mutation index
  mutationIndex[key] = nodeId;
  
  // Also update node-specific mutation index if needed
  auto nodeIt = nodeStates.find(nodeId);
  if (nodeIt != nodeStates.end()) {
    nodeIt->second.mutationIndex[key] = nodeId;
  }
}

// Set the s-mer size for distance calculations
void StateManager::setSmerSize(int s) {
  if (s > 0) {
    smerSize = s;
  }
}

// Get recomputation ranges for a node
std::vector<coordinates::CoordRange> StateManager::getRecompRanges(std::basic_string_view<char, std::char_traits<char>> nodeId) const {
  auto it = nodeStates.find(std::string(nodeId));
  if (it != nodeStates.end()) {
    return it->second.recompRanges;
  }
  
  // Return empty vector if node not found
  return {};
}

// Expand recomputation ranges to ensure complete k-mers
std::vector<coordinates::CoordRange> StateManager::expandRecompRanges(
    std::basic_string_view<char, std::char_traits<char>> nodeId,
    const std::vector<coordinates::CoordRange>& ranges, 
    int k) {
  
  if (ranges.empty()) {
    return {};
  }
  
  std::vector<coordinates::CoordRange> expandedRanges;
  expandedRanges.reserve(ranges.size());
  
  // Expand each range to ensure it contains complete k-mers
  for (const auto& range : ranges) {
    // Expand range by (k-1) in each direction to ensure complete k-mers
    coordinates::CoordRange expandedRange{
      std::max<int64_t>(0, range.start - (k - 1)),
      std::min<int64_t>(numCoords, range.end + (k - 1))
    };
    
    // Add the expanded range
    expandedRanges.push_back(expandedRange);
  }
  
  // Now merge overlapping ranges
  std::vector<coordinates::CoordRange> mergedRanges;
  
  // Sort expanded ranges by start position
  std::sort(expandedRanges.begin(), expandedRanges.end(),
            [](const coordinates::CoordRange& a, const coordinates::CoordRange& b) {
              return a.start < b.start;
            });
  
  // Merge overlapping ranges
  for (const auto& range : expandedRanges) {
    if (mergedRanges.empty() || mergedRanges.back().end < range.start) {
      // No overlap with previous range, add new range
      mergedRanges.push_back(range);
    } else {
      // Overlap with previous range, extend it
      mergedRanges.back().end = std::max(mergedRanges.back().end, range.end);
    }
  }
  
  return mergedRanges;
}

// Add a seed to a block for tracking
void StateManager::addSeedToBlock(int32_t blockId, int64_t seedPos) {
  if (blockId >= 0 && blockId < static_cast<int32_t>(blockToSeeds.size())) {
    blockToSeeds[blockId].insert(seedPos);
  }
}

// Clear a seed at a specific position
void StateManager::clearSeedAtPosition(int64_t pos) {
  if (pos >= 0 && pos < static_cast<int64_t>(positionSeeds.size())) {
    positionSeeds[pos].reset();
  }
}

// Map global position to block coordinates
coordinates::BlockCoordinate StateManager::mapGlobalToBlockCoords(
    std::basic_string_view<char, std::char_traits<char>> nodeId,
    int64_t globalPos) const {
  
  // Use cached data for fast lookup if available
  if (globalPosCacheInitialized && 
      globalPos >= 0 && 
      globalPos < static_cast<int64_t>(globalPosCache.posToBlockId.size())) {
    
    int32_t blockId = globalPosCache.posToBlockId[globalPos];
    if (blockId >= 0) {
      // Compute offset within block
      int64_t blockStart = globalPosCache.blockStartOffsets[blockId];
      int64_t offset = globalPos - blockStart;
      
      return coordinates::BlockCoordinate{
        blockId,             // blockId
        static_cast<int32_t>(offset), // blockPos (relative position)
        globalPos            // globalPos
      };
    }
  }
  
  // Fallback method: search through all blocks
  const auto& sortedBlocks = getSortedActiveBlocks(nodeId);
  
  for (const auto& [blockId, range] : sortedBlocks) {
    if (globalPos >= range.start && globalPos < range.end) {
      // Found the block containing this position
      int64_t offset = globalPos - range.start;
      
      return coordinates::BlockCoordinate{
        blockId,             // blockId
        static_cast<int32_t>(offset), // blockPos (relative position)
        globalPos            // globalPos
      };
    }
  }
  
  // Position not found in any block
  return coordinates::BlockCoordinate{
    -1,   // Invalid blockId
    -1,   // Invalid blockPos
    globalPos  // Original global position
  };
}

// Extract a sequence from a range of positions
std::pair<std::string, std::vector<int64_t>> StateManager::extractSequence(
    std::basic_string_view<char, std::char_traits<char>> nodeId, 
    const coordinates::CoordRange& range, 
    bool skipGaps) {
  
  // Check cache first
  std::string cacheKey = std::string(nodeId) + ":" + 
                        std::to_string(range.start) + ":" + 
                        std::to_string(range.end) + ":" + 
                        (skipGaps ? "1" : "0");
  
  auto cacheIt = sequenceCache.find(cacheKey);
  if (cacheIt != sequenceCache.end()) {
    return cacheIt->second;
  }
  
  // Get gap map for this node
  auto nodeGapMap = getNodeGapMap(nodeId);
  
  // Prepare result
  std::string sequence;
  std::vector<int64_t> positions;
  
  // Reserve reasonable capacity
  int64_t rangeSize = range.end - range.start;
  sequence.reserve(rangeSize);
  positions.reserve(rangeSize);
  
  // Get sorted active blocks for this node
  const auto& sortedBlocks = getSortedActiveBlocks(nodeId);
  
  // Extract sequence for each overlapping block
  for (const auto& [blockId, blockRange] : sortedBlocks) {
    // Check if this block overlaps with our target range
    if (blockRange.end <= range.start || blockRange.start >= range.end) {
      continue; // No overlap
    }
    
    // Calculate overlap
    int64_t overlapStart = std::max(blockRange.start, range.start);
    int64_t overlapEnd = std::min(blockRange.end, range.end);
    
    // Check if block is active and inverted
    bool isActive = isBlockOn(nodeId, blockId);
    if (!isActive) {
      continue; // Skip inactive blocks
    }
    
    bool isInverted = isBlockInverted(nodeId, blockId);
    
    // Get block sequence
    auto seqIt = blockSequences.find(blockId);
    if (seqIt == blockSequences.end()) {
      continue; // Skip blocks with no sequence
    }
    
    const auto& blockSeq = seqIt->second;
    
    // Extract nucleotides from overlap range
    for (int64_t pos = overlapStart; pos < overlapEnd; pos++) {
      // Skip gaps if requested
      if (skipGaps && isGapPosition(nodeGapMap, pos)) {
        continue;
      }
      
      // Find which block and position this belongs to
      try {
        // Map to block coordinates
        auto blockCoord = mapGlobalToBlockCoords(nodeId, pos);
        
        if (blockCoord.blockId == blockId) {
          // Get character at this position
          try {
            int32_t nucPos = blockCoord.blockPos;
            // For the root sequence or inactive blocks, access the block sequence directly
            if (nucPos >= 0 && nucPos < static_cast<int32_t>(blockSeq.length())) {
              char c = blockSeq[nucPos];
              
              // Complement nucleotide if block is inverted
              if (isInverted && isNonGapChar(c)) {
                switch (c) {
                  case 'A': c = 'T'; break;
                  case 'T': c = 'A'; break;
                  case 'C': c = 'G'; break;
                  case 'G': c = 'C'; break;
                  case 'a': c = 't'; break;
                  case 't': c = 'a'; break;
                  case 'c': c = 'g'; break;
                  case 'g': c = 'c'; break;
                  default: break;
                }
              }
              
              sequence.push_back(c);
              positions.push_back(pos);
            }
          } catch (const std::exception& e) {
            logging::debug("Error getting character at pos {} in block {}: {}", 
                          pos, blockId, e.what());
          }
        }
      } catch (const std::exception& e) {
        logging::debug("Error mapping position {}: {}", pos, e.what());
      }
    }
  }
  
  // Cache the result
  auto result = std::make_pair(sequence, positions);
  sequenceCache[cacheKey] = result;
  
  return result;
}

} // namespace state


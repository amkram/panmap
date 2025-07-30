#include "state.hpp"
#include "coordinates.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "panman.hpp" // Assuming panmanUtils types are here
#include "seeding.hpp"
#include "seq_utils.hpp"
#include <algorithm>
#include <atomic> // Added
#include <cassert>
#include <cstddef>
#include <cstdint> // Added
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <exception> // Added
#include <functional> // Added
#include <iterator>
#include <memory> // Added
#include <mutex> // For std::mutex and std::lock_guard
#include <optional> // Added
#include <shared_mutex> // For std::shared_mutex
#include <sstream> // Added
#include <stdexcept> // Added
#include <string> // Added explicitly

// OPTIMIZATION: Branch prediction hints for performance
#ifdef __GNUC__
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define LIKELY(x)   (x)
#define UNLIKELY(x) (x)
#endif
#include <string_view>
#include <tbb/tbb.h>
#include <tuple> // Added
#include <unordered_map> // Added for consistency
#include <unordered_set> // Added explicitly
#include <utility> // Added
#include <vector> // Added explicitly
#include "absl/container/flat_hash_set.h"
#include <fstream>

namespace state {
    // Define the global counters that are causing linker errors
    std::atomic<size_t> totalKmerAttempts{0};
    std::atomic<size_t> successfulKmerExtractions{0};
    std::atomic<size_t> kmerFailedDueToLength{0};
    std::atomic<size_t> kmerFailedDueToOnlyGaps{0};
}

namespace state {

// Forward declaration for the diagnoseBlockStatus function
void diagnoseBlockStatus(const NodeState& nodeState, const std::string& nodeId, int32_t blockId);

// Define the thread-local cache storage
thread_local StateManager::NodeStateCache StateManager::tlsNodeCache;

// Define the thread-local cache storage
// Use the CacheType alias defined in the header to match the declaration
thread_local HierarchicalStore<char>::CacheType HierarchicalStore<char>::tlsCache;

// Implementation of HierarchicalStore<char>::get method
std::optional<char> HierarchicalStore<char>::get(const PositionKey& key) const {
  // First, check local store
  {
    std::lock_guard<std::mutex> lock(storeMutex);
    auto it = localValues.find(key);
    if (it != localValues.end()) {
      localCacheHits.fetch_add(1, std::memory_order_relaxed);
      return it->second;
    }
  }
  
  // If not found locally, check parent
  auto parentPtr = parent.lock();
  if (parentPtr) {
    parentTraversals.fetch_add(1, std::memory_order_relaxed);
    return parentPtr->get(key);
  }
  
  // Not found anywhere
  return std::nullopt;
}

/*
 * Character Data Access Pattern:
 * 
 * Hierarchical delta-based approach to storing character data:
 * 
 * 1. Initial sequences are stored in the global blockSequences map
 * 2. Each node only stores mutations relative to its parent
 * 3. Ancestor relationship is maintained via parentId reference
 * 4. When accessing characters:
 *    a. First check current node's characterStore (hierarchical) for mutations
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

// Gap run tracking is now part of NodeState struct in the header file

// Global variables for initialization tracking
std::unordered_set<std::string> nodesInProgress;
// std::recursive_mutex initializationMutex; // REMOVE THIS LINE

// Error reporting mutex
std::mutex outputMutex;
// Delete this line: mutable std::shared_mutex nodeMutex; // This comment is wrong, nodeMutex IS used
std::mutex logMutex;

std::unordered_map<std::string, int64_t> nodeDfsIndices;
std::vector<std::string> nodeIdsByDfsIndex; // Added to map index back to ID

StateManager::StateManager() : numCoords(0), numBlocks(0), kmerSize(0) {
  // No cache counters initialization needed; per-node caches are removed
}

StateManager::StateManager(size_t numCoordinates)
    : numCoords(numCoordinates), numBlocks(0), kmerSize(0) {
  // Initialize position seeds array
  positionSeeds.resize(numCoordinates);

  // Early initialization of empty gap list array
  // Will be properly sized when setNumBlocks is called
  gapListLengthArray.clear();
}

void StateManager::clearCaches() {
    // Node-specific caches have been removed, so this is now just a placeholder
    // that logs a debug message for clarity
    std::cerr << "[DEBUG] clearCaches called - node-specific caches removed" << std::endl;
}

// Initialize DFS indices for traversal
void StateManager::initializeDfsIndices(panmanUtils::Tree *tree) {
  if (!tree || !tree->root) {
    std::cerr << "[ERROR] Invalid tree in initializeDfsIndices" << std::endl;
    return;
  }

  std::cerr << "[DEBUG] Starting DFS index initialization with " << tree->allNodes.size() << " nodes in tree" << std::endl;
  std::cerr << "[DEBUG] Root node: " << tree->root->identifier << std::endl;
  
  // Resize the reverse mapping vector
  nodeIdsByDfsIndex.clear();
  nodeIdsByDfsIndex.resize(tree->allNodes.size()); // Exact size
  nodeDfsIndices.clear(); // Clear forward map too
  
  int64_t dfsIndex = 0;
  
  // First assign index to root node
  nodeDfsIndices[tree->root->identifier] = dfsIndex;
  nodeIdsByDfsIndex[dfsIndex] = tree->root->identifier; // Store reverse mapping
  dfsIndex++;
  std::cerr << "[DEBUG] Assigned DFS index 0 to root node " << tree->root->identifier << std::endl;
  
  // Use DFS traversal for the rest
  initializeDfsIndices(tree, tree->root, dfsIndex); // Pass vector by ref if needed by recursive call
  
  // Check for nodes without indices and assign them
  std::cerr << "[INFO] Initial DFS traversal assigned " << dfsIndex << " indices" << std::endl;
  
  int missingCount = 0;
  for (const auto& [nodeId, node] : tree->allNodes) {
    if (nodeDfsIndices.find(nodeId) == nodeDfsIndices.end()) {
      // Ensure vector is large enough if fallback assignment happens
      if (static_cast<size_t>(dfsIndex) >= nodeIdsByDfsIndex.size()) {
           nodeIdsByDfsIndex.resize(dfsIndex + 1);
      }
      nodeDfsIndices[nodeId] = dfsIndex;
      nodeIdsByDfsIndex[dfsIndex] = nodeId; // Store reverse mapping
      dfsIndex++;
      
      if (missingCount <= 5) {
        std::cerr << "[INFO] Assigned fallback DFS index " << (dfsIndex-1) << " to node " << nodeId << std::endl;
      }
    }
  }
  
  if (missingCount > 0) {
    std::cerr << "[INFO] Assigned fallback DFS indices to " << missingCount << " nodes not covered by standard traversal" << std::endl;
  }
  
  // Final resize to exact count if fallbacks occurred
  nodeIdsByDfsIndex.resize(dfsIndex);
  
  std::cerr << "[INFO] DFS indexing complete, assigned " << dfsIndex << " indices total" << std::endl;
}

// Recursive helper for DFS index initialization
void StateManager::initializeDfsIndices(panmanUtils::Tree *tree,
                                        panmanUtils::Node *node,
                                        int64_t &dfsIndex) {
  if (!node) {
    std::cerr << "[WARN] Null node encountered in DFS traversal" << std::endl;
    return;
  }

  // Root node already assigned in the calling function
  if (node != tree->root) {
    // Assign current DFS index to this node
    nodeDfsIndices[node->identifier] = dfsIndex;
    // Ensure vector is large enough
    if (static_cast<size_t>(dfsIndex) >= nodeIdsByDfsIndex.size()) {
         nodeIdsByDfsIndex.resize(dfsIndex + 1);
    }
    nodeIdsByDfsIndex[dfsIndex] = node->identifier; // Store reverse mapping
    dfsIndex++;
    
    // Log first few assignments
    if (dfsIndex <= 5) {
      std::cerr << "[INFO] Assigned DFS index " << (dfsIndex-1) << " to node " << node->identifier << std::endl;
    }
  }

  // Recursively process children
  for (auto *child : node->children) {
    if (!child) {
      std::cerr << "[WARN] Null child encountered in node " << node->identifier << std::endl;
      continue;
    }
    initializeDfsIndices(tree, child, dfsIndex); // Pass vector by ref if needed
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
  
  // Initialize gap list array with correct number of blocks
  // This happens early, before any setGapListLength calls
  if (gapListLengthArray.size() != blockCount) {
    gapListLengthArray.resize(blockCount);
    std::cerr << "[DEBUG] Initialized gap list array with " << blockCount << " blocks" << std::endl;
  }
  
  std::cerr << "[DEBUG] Set number of blocks to " << blockCount << std::endl;
}

// Initialize a node's state
void StateManager::initializeNode(const std::string &nodeId) {
  // Acquire the lock at the beginning of the function
  std::unique_lock<std::shared_mutex> lock(nodeMutex);

  // Get the static verbosity control variable
  static bool& verbose_logging = getInitLoggingRef();

  // Check if node is already initialized
  if (nodeStates.find(nodeId) != nodeStates.end()) {
    if (verbose_logging) {
      std::cerr << "[DEBUG] Node " << nodeId << " already initialized" << std::endl;
    }
    return;
  }

  // Create a new node state
  nodeStates.emplace(nodeId, NodeState());
  if (verbose_logging) {
    std::cerr << "[DEBUG] Created new node state for node " << nodeId << std::endl;
  }

  // Lazy initialization of rootGapMap
  if (!rootGapMap) {
    rootGapMap = std::make_shared<HierarchicalGapMap>();
  }

  // Lazy initialization of rootCharacterStore
  static std::shared_ptr<HierarchicalCharacterStore> rootCharacterStore;
  if (!rootCharacterStore) {
    rootCharacterStore = std::make_shared<HierarchicalCharacterStore>();
  }
  
  // Get or create root seed store for initialization
  static std::shared_ptr<HierarchicalSeedStore> rootSeedStore;
  if (!rootSeedStore) {
    rootSeedStore = std::make_shared<HierarchicalSeedStore>();
    if (verbose_logging) {
      std::cerr << "[DEBUG] Created root seed store for initialization" << std::endl;
    }
  }

  // Find parent node in hierarchy
  std::string parentId;
  auto hierarchyIt = nodeHierarchy.find(nodeId);
  if (hierarchyIt != nodeHierarchy.end()) {
    parentId = hierarchyIt->second.parentId;
    
    // Set up hierarchical structure inheritance from nodeHierarchy directly
    auto& nodeState = nodeStates[nodeId];
    nodeState.parentId = parentId;
  }
  
  // Get reference to the node state
  auto& nodeState = nodeStates[nodeId];
    
  // Initialize gap map and character store from nodeHierarchy if available
  if (hierarchyIt != nodeHierarchy.end()) {
    if (hierarchyIt->second.gapMap) {
      nodeState.gapMap = std::make_shared<HierarchicalGapMap>(hierarchyIt->second.gapMap);
      if (verbose_logging) {
        std::cerr << "[DEBUG] Node " << nodeId << " inherited gap map from hierarchy" << std::endl;
      }
    }
    
    if (hierarchyIt->second.characterStore) {
      nodeState.characterStore = std::make_shared<HierarchicalCharacterStore>(hierarchyIt->second.characterStore);
      if (verbose_logging) {
        std::cerr << "[DEBUG] Node " << nodeId << " inherited character store from hierarchy" << std::endl;
      }
    }
  }
   

  // If this is not the root and we have a parent,
  // propagate state from parent to this node
  if (!parentId.empty()) {
    // With topological initialization, the parent should already be initialized
    // No need to recursively initialize it
    
    // Check if parent exists in nodeStates - it should if we're using topological order
    if (nodeStates.find(parentId) == nodeStates.end()) {
      if (verbose_logging) {
        std::cerr << "[WARN] Parent " << parentId << " for node " << nodeId << " not found in nodeStates - initialization order issue" << std::endl;
      }
      
      // No need to initialize the parent - it should have been initialized before
      // in the topological order. If not, that's an error in our sorting or the tree structure.
    } else {
      // OPTIMIZATION: Store parent reference for later bulk materialization
      // Don't materialize here to avoid quadratic behavior during initialization
      if (verbose_logging) {
        std::cerr << "[DEBUG] Node " << nodeId << " parent reference set to " << parentId 
                  << " (materialization deferred)" << std::endl;
      }
    }
  } else {
    // Root node setup
    auto& nodeState = nodeStates[nodeId];
    
    // Ensure gapMap and characterStore are assigned
    if (!nodeState.gapMap) nodeState.gapMap = rootGapMap; 
    if (!nodeState.characterStore) nodeState.characterStore = rootCharacterStore;
    
    // OPTIMIZATION: Root node starts with empty materialized block state
    nodeState.materializeBlockState(nullptr);
    if (verbose_logging) {
      std::cerr << "[DEBUG] Root node " << nodeId << " initialized with empty materialized block state" << std::endl;
    }
  }

  // Final verification and setup
  auto& finalNodeState = nodeStates[nodeId];
  
  // Ensure node has a valid gap map
  if (finalNodeState.gapMap == nullptr) {
    std::cerr << "[WARN] Node " << nodeId << " gapMap still null after init, assigning root." << std::endl;
    finalNodeState.gapMap = rootGapMap; // Assign root as fallback
  }
  if (finalNodeState.characterStore == nullptr) {
    std::cerr << "[WARN] Node " << nodeId << " characterStore still null after init, assigning root." << std::endl;
    finalNodeState.characterStore = rootCharacterStore; // Assign root as fallback
  }
  
  // IMPORTANT FIX: Clear any inherited recomputation ranges
  // This ensures ranges used for recomputation are always valid for this node's coordinate space
  finalNodeState.recompRanges.clear();
  if (verbose_logging) {
    std::cerr << "[DEBUG] Cleared inherited recomputation ranges for node " << nodeId << std::endl;
  }
  
  // Materialize character data for performance optimization
  if (!parentId.empty()) {
    try {
      auto parentStateIt = nodeStates.find(parentId);
      if (parentStateIt != nodeStates.end()) {
        // Get parent's materialized characters and local character store
        const NodeState* parentState = &parentStateIt->second;
        const HierarchicalCharacterStore* characterStore = finalNodeState.characterStore.get();
        
        // Note: Materialization is deferred until explicitly called via materializeNodeState()
        // This ensures proper parent-child inheritance order
        if (verbose_logging) {
          std::cerr << "[DEBUG] Initialization complete for node " << nodeId << " (materialization deferred)" << std::endl;
        }
      }
    } catch (const std::exception& e) {
      logging::warn("Failed to materialize character data for node {}: {}", nodeId, e.what());
    }
  }
  
  // Reset cache for this node
  resetNodeCache(nodeId);
  
  if (verbose_logging) {
    std::cerr << "[DEBUG] Node " << nodeId << " initialized with " << finalNodeState.activeBlocks.size() << " active blocks" << std::endl;
  }
}

// OPTIMIZATION: Bulk materialize all node block states for maximum performance
void StateManager::materializeAllNodeBlockStates() {
  std::unique_lock<std::shared_mutex> lock(nodeMutex);
  
  // Build topological order for efficient materialization
  std::vector<std::string> topoOrder;
  std::unordered_set<std::string> visited;
  std::unordered_set<std::string> processing;
  
  // Helper function for DFS topological sort
  std::function<void(const std::string&)> dfsVisit = [&](const std::string& nodeId) {
    if (processing.count(nodeId)) {
      std::cerr << "[WARN] Cycle detected in node hierarchy at " << nodeId << std::endl;
      return;
    }
    if (visited.count(nodeId)) return;
    
    processing.insert(nodeId);
    
    auto it = nodeStates.find(nodeId);
    if (it != nodeStates.end() && !it->second.parentId.empty()) {
      dfsVisit(it->second.parentId);
    }
    
    processing.erase(nodeId);
    visited.insert(nodeId);
    topoOrder.push_back(nodeId);
  };
  
  // Visit all nodes to build topological order
  for (const auto& [nodeId, nodeState] : nodeStates) {
    if (!visited.count(nodeId)) {
      dfsVisit(nodeId);
    }
  }
  
  // Materialize in topological order (parents before children)
  size_t materializedCount = 0;
  for (const std::string& nodeId : topoOrder) {
    auto it = nodeStates.find(nodeId);
    if (it == nodeStates.end()) continue;
    
    auto& nodeState = it->second;
    if (!nodeState.blockStateMaterialized) {
      // Find parent state
      const NodeState* parentState = nullptr;
      if (!nodeState.parentId.empty()) {
        auto parentIt = nodeStates.find(nodeState.parentId);
        if (parentIt != nodeStates.end()) {
          parentState = &parentIt->second;
        }
      }
      
      // Materialize this node's block state
      nodeState.materializeBlockState(parentState);
      materializedCount++;
    }
    
    // OPTIMIZATION: Cache parent state pointer for ultra-fast access
    nodeState.cacheParentState(nodeStates);
  }
  
  std::cerr << "[INFO] Bulk materialized block states for " << materializedCount << " nodes" << std::endl;
}

// Track which block is used by which node
void StateManager::trackBlockUsage(int32_t blockId, const std::string &nodeId) {
  if (blockId >= 0 && blockId < static_cast<int32_t>(numBlocks)) {
    blockToSeeds[blockId]; // Create entry if doesn't exist
    std::cerr << "[DEBUG] Block " << blockId << " used by node " << nodeId << std::endl; // Replaced trace with debug
  }
}

// Get block status - inheritance-based approach
bool StateManager::isBlockOn(std::string_view nodeId, int32_t blockId) const {
    std::shared_lock<std::shared_mutex> lock(nodeMutex);
    return isBlockOn_unsafe(nodeId, blockId);
}

bool StateManager::isBlockOn_unsafe(std::string_view nodeId, int32_t blockId) const {
    // CRITICAL OPTIMIZATION: Avoid string construction in hot path
    // Use string_view directly with a specialized cache lookup
    const NodeState* nodeState = tlsNodeCache.getCachedStateByView(nodeId, nodeStates);
    
    if (LIKELY(nodeState)) {
        // OPTIMIZATION: Fast path for already materialized state (most common case)
        if (LIKELY(nodeState->blockStateMaterialized)) {
            return nodeState->isBlockActive(blockId);
        }
        
        // OPTIMIZATION: Check local mutations first - use likely/unlikely hints
        if (UNLIKELY(nodeState->localBlockActivations.contains(blockId))) {
            return true;
        }
        if (UNLIKELY(nodeState->localBlockDeactivations.contains(blockId))) {
            return false;
        }
        
        // OPTIMIZATION: Cache parent node pointer to avoid repeated string conversions
        const NodeState* currentState = nodeState;
        while (UNLIKELY(!currentState->parentId.empty())) {
            // Use cached parent pointer if available
            if (currentState->cachedParentState) {
                currentState = currentState->cachedParentState;
            } else {
                // Fallback to string lookup (should be rare after warmup)
                currentState = tlsNodeCache.getCachedStateByView(currentState->parentId, nodeStates);
                if (!currentState) break;
            }
            
            // If parent has materialized state, use it
            if (LIKELY(currentState->blockStateMaterialized)) {
                return currentState->isBlockActive(blockId);
            }
            
            // Check parent's local mutations
            if (UNLIKELY(currentState->localBlockActivations.contains(blockId))) {
                return true;
            }
            if (UNLIKELY(currentState->localBlockDeactivations.contains(blockId))) {
                return false;
            }
        }
        
        // Fallback: use legacy activeBlocks for compatibility
        return nodeState->activeBlocks.contains(blockId);
    }
    return false; // Node not found
}

bool StateManager::isBlockInverted(std::string_view nodeId, int32_t blockId) const {
    std::shared_lock<std::shared_mutex> lock(nodeMutex);
    return isBlockInverted_unsafe(nodeId, blockId);
}

bool StateManager::isBlockInverted_unsafe(std::string_view nodeId, int32_t blockId) const {
    auto it = nodeStates.find(std::string(nodeId));
    if (it != nodeStates.end()) {
        // Check for explicit inversion setting at this node level only
        auto orientIt = it->second.blockOrientation.find(blockId);
        if (orientIt != it->second.blockOrientation.end() && orientIt->second.has_value()) {
            return !orientIt->second.value(); // false in map means inverted
        }
    }
    return false; // Default to forward orientation
}

// Unsafe version of getActiveBlockRanges - caller must hold the lock
std::vector<std::pair<int32_t, coordinates::CoordRange>>
StateManager::getActiveBlockRanges(std::string_view nodeId) const {
  std::vector<std::pair<int32_t, coordinates::CoordRange>> activeBlocks;
  
  // NOTE: This is the unsafe version - caller must hold the lock
  
  try {
    // SIMPLIFIED: Direct lookup without cache overhead - cache was 37.74% bottleneck
    std::string nodeIdStr(nodeId); // Convert once, reuse
    auto it = nodeStates.find(nodeIdStr);
    if (it == nodeStates.end()) {
      return activeBlocks; // Empty result for non-existent node
    }
    const NodeState* nodeState = &it->second;
    
    // SIMPLIFIED: Always use materialized state path if available
    if (nodeState->blockStateMaterialized) {
      const size_t totalBlocks = numBlocks;
      
      // OPTIMIZATION: Use vectorized batch block status checking
      // Process blocks in larger chunks for better cache locality
      const int32_t batchSize = 128; // Increased from 64 for better vectorization
      activeBlocks.reserve(nodeState->materializedActiveBlocks.size());
      
      // Pre-allocate batch status vector to avoid repeated allocations
      // NOTE: Using std::vector<char> instead of std::vector<bool> because std::vector<bool>
      // is specialized and doesn't provide a usable .data() method
      thread_local std::vector<char> batchStatus;
      batchStatus.clear();
      batchStatus.reserve(batchSize);
      
      for (int32_t startBlock = 0; startBlock < static_cast<int32_t>(totalBlocks); startBlock += batchSize) {
        int32_t remainingBlocks = std::min(batchSize, static_cast<int32_t>(totalBlocks) - startBlock);
        
        // Get batch status for this range - use a temporary bool array
        batchStatus.resize(remainingBlocks);
        std::unique_ptr<bool[]> tempBoolArray(new bool[remainingBlocks]);
        nodeState->getBlockStatusBatchFast(startBlock, remainingBlocks, tempBoolArray.get());
        
        // Copy results to our char vector for processing
        for (int32_t i = 0; i < remainingBlocks; ++i) {
          batchStatus[i] = tempBoolArray[i] ? 1 : 0;
        }
        
        // Process results and add active blocks - unroll inner loop for better performance
        for (int32_t i = 0; i < remainingBlocks; ++i) {
          if (LIKELY(!batchStatus[i])) continue;
          
          int32_t blockId = startBlock + i;
          auto rangeIt = blockRanges.find(blockId);
          if (LIKELY(rangeIt != blockRanges.end())) {
            activeBlocks.emplace_back(blockId, rangeIt->second);
          }
        }
      }
    } else {
      // OPTIMIZATION: Fallback path - use more efficient collection
      const size_t totalBlocks = numBlocks;
      
      // Pre-allocate based on typical active block ratio
      activeBlocks.reserve(totalBlocks / 4); // Assume ~25% blocks are active
      
      // Use batch checking even for inheritance model
      for (int32_t blockId = 0; blockId < static_cast<int32_t>(totalBlocks); ++blockId) {
        if (isBlockOn_unsafe(nodeId, blockId)) {
          auto rangeIt = blockRanges.find(blockId);
          if (LIKELY(rangeIt != blockRanges.end())) {
            activeBlocks.emplace_back(blockId, rangeIt->second);
          }
        }
      }
    }
    
    // Sort by global position for consistent ordering
    std::sort(activeBlocks.begin(), activeBlocks.end(), 
              [](const auto& a, const auto& b) {
                return a.second.start < b.second.start;
              });
              
  } catch (const std::exception& e) {
    logging::err("Exception in getActiveBlockRanges for node {}: {}", std::string(nodeId), e.what());
  }
  
  return activeBlocks;
}
// Set block status using inheritance-based approach
void StateManager::setBlockOn(std::string_view nodeId, int32_t blockId, bool on) {
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // Use unique lock for modifications
  
  // Find the node state directly
  std::string strNodeId(nodeId); // Convert once
  auto it = nodeStates.find(strNodeId);
  if (it == nodeStates.end()) {
    return;
  }
  
  // Get reference to node state to modify it
  auto& nodeState = it->second;
  
  // NEW: Set local block mutation instead of maintaining complete activeBlocks
  nodeState.setLocalBlockState(blockId, on);
  
  // Re-materialize this node's block state
  auto parentIt = nodeStates.find(nodeState.parentId);
  const NodeState* parentState = (parentIt != nodeStates.end()) ? &parentIt->second : nullptr;
  nodeState.materializeBlockState(parentState);
  
  // Invalidate materialized state for all children (they need to re-inherit)
  for (const auto& [childId, hierarchy] : nodeHierarchy) {
    if (hierarchy.parentId == strNodeId) {
      auto childStateIt = nodeStates.find(childId);
      if (childStateIt != nodeStates.end()) {
        childStateIt->second.blockStateMaterialized = false;
        // Children will re-materialize lazily when accessed
      }
    }
  }
  
  // Legacy: Also update old activeBlocks for backward compatibility during transition
  if (on) {
    nodeState.activeBlocks.insert(blockId);
  } else {
    nodeState.activeBlocks.erase(blockId);
  }
  
  // Reset node cache since block status changed
  resetNodeCache(nodeId);
}

void StateManager::setBlockForward(std::string_view nodeId, int32_t blockId,
                                   bool forward) {
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // ADD THIS: Use unique lock for modifications
  std::string strNodeId(nodeId); // Convert once
  
  // Ensure node exists before modifying it
  auto nodeIt = nodeStates.find(strNodeId);
  if (nodeIt == nodeStates.end()) {
    lock.unlock(); // Release lock temporarily
    getNodeState(strNodeId); // This will initialize the node properly
    lock.lock(); // Re-acquire lock
    nodeIt = nodeStates.find(strNodeId);
    if (nodeIt == nodeStates.end()) {
      std::cerr << "ERROR: Failed to initialize node " << strNodeId << std::endl;
      return;
    }
  }
  
  auto &nodeState = nodeIt->second;
  nodeState.setBlockForward(blockId, forward);

  // Reset cache since block orientation changed
  resetNodeCache(strNodeId);
}

void StateManager::setBlockInverted(std::string_view nodeId, int32_t blockId,
                                    bool inverted) {
  bool forward = !inverted;
  setBlockForward(nodeId, blockId, forward); // Calls the one above, conversion happens there.
}

// Get character at position using hierarchical lookup with node's mutation index
char StateManager::getCharAtPosition(std::string_view nodeId_sv, int32_t blockId,
                                     int32_t nucPos, int32_t gapPos) const {
  // OPTIMIZATION: Avoid string construction in hot path
  try {
    // Create position key once
    state::PositionKey posKey{blockId, nucPos, gapPos};
    
    // STEP 1: Check if block is OFF for this node → return '-'
    bool isTargetNodeBlockActive = isBlockOn_unsafe(nodeId_sv, blockId);
    if (!isTargetNodeBlockActive) {
      return '-';
    }
    
    // STEP 2: Check current node's character data (accumulated mutations)
    const NodeState* nodeState = tlsNodeCache.getCachedStateByView(nodeId_sv, nodeStates);
    if (nodeState != nullptr && nodeState->characterStore) {
      // Try materialized state first for O(1) lookup
      std::optional<char> materializedChar = nodeState->getMaterializedCharacter(posKey);
      if (materializedChar.has_value()) {
        char canonicalChar = materializedChar.value();
        
        // Apply inversion if needed before returning
        bool isBlockInverted = isBlockInverted_unsafe(nodeId_sv, blockId);
        if (isBlockInverted) {
          return seq_utils::getComplementCharacter(canonicalChar);
        } else {
          return canonicalChar;
        }
      }
      
      // Check node's local store (no hierarchical lookup)
      std::optional<char> localCharOpt = nodeState->characterStore->getLocal(posKey);
      if (localCharOpt.has_value()) {
        char canonicalChar = localCharOpt.value();
        
        // Apply inversion if needed before returning
        bool isBlockInverted = isBlockInverted_unsafe(nodeId_sv, blockId);
        if (isBlockInverted) {
          return seq_utils::getComplementCharacter(canonicalChar);
        } else {
          return canonicalChar;
        }
      }
    }
    
    // STEP 3: Reference sequence fallback (set when block first turns ON)
    auto blockSeqIt = blockSequences.find(blockId);
    if (blockSeqIt != blockSequences.end()) {
      state::PositionKey structuralKey_for_fallback = state::PositionKey::create(blockId, nucPos, gapPos);
      int64_t flatIndex = -1;
      auto blockMapIt = blockRootCharFlatIndices.find(blockId);
      if (blockMapIt != blockRootCharFlatIndices.end()) {
        auto& blockSpecificMap = blockMapIt->second; 
        auto flatIndexIt = blockSpecificMap.find(structuralKey_for_fallback);
        if (flatIndexIt != blockSpecificMap.end()) {
          flatIndex = flatIndexIt->second;
        }
      }

      if (flatIndex >= 0 && flatIndex < static_cast<int64_t>(blockSeqIt->second.size())) {
        char canonicalChar = blockSeqIt->second[flatIndex];
        
        // Apply inversion if needed before returning
        bool isBlockInverted = isBlockInverted_unsafe(nodeId_sv, blockId);
        if (isBlockInverted) {
          return seq_utils::getComplementCharacter(canonicalChar);
        } else {
          return canonicalChar;
        }
      }
    }
    
    // If no reference found, return gap
    return '-';
    
  } catch (const std::exception &e) {
    logging::warn("Exception in getCharAtPosition for node {}, block {}, position ({}, {}): {}", 
                std::string(nodeId_sv), blockId, nucPos, gapPos, e.what());
    return '-';
  }
}

// PERFORMANCE OPTIMIZATION: Bulk character extraction to avoid individual StateManager calls
std::vector<char> StateManager::getCharactersAtBlock(std::string_view nodeId, int32_t blockId,
                                                    int32_t startNucPos, int32_t startGapPos, int32_t length) const {
  std::vector<char> result;
  result.reserve(length);
  
  if (length <= 0) {
    return result;
  }
  
  std::string nodeId_str(nodeId);
  
  try {
    // Cache common lookups for efficiency
    const NodeState* nodeState = nullptr;
    {
      auto it = nodeStates.find(nodeId_str);
      if (it != nodeStates.end()) {
        nodeState = &it->second;
      }
    }
    
    // Cache block sequence lookup
    const std::string* blockSequence = nullptr;
    auto blockSeqIt = blockSequences.find(blockId);
    if (blockSeqIt != blockSequences.end()) {
      blockSequence = &blockSeqIt->second;
    }
    
    // Cache block orientation
    bool isBlockInverted = isBlockInverted_unsafe(nodeId, blockId);
    bool isBlockActive = isBlockOn_unsafe(nodeId, blockId);
    
    if (!isBlockActive) {
      // If block is not active, all positions return gap character
      result.assign(length, '-');
      return result;
    }
    
    // Extract characters in sequence with minimal overhead
    for (int32_t i = 0; i < length; ++i) {
      int32_t currentNucPos = startNucPos;
      int32_t currentGapPos = startGapPos;
      
      if (startGapPos == -1) {
        currentNucPos = startNucPos + i;
        currentGapPos = -1;
      } else {
        currentNucPos = startNucPos;
        currentGapPos = startGapPos + i;
      }
      
      char canonicalChar = '?';
      bool foundInStore = false;
      
      // Quick check in materialized state first
      if (nodeState && nodeState->characterStore) {
        state::PositionKey posKey{blockId, currentNucPos, currentGapPos};
        std::optional<char> materializedChar = nodeState->getMaterializedCharacter(posKey);
        if (materializedChar.has_value()) {
          canonicalChar = materializedChar.value();
          foundInStore = true;
        } else {
          // Check local character store
          std::optional<char> localCharOpt = nodeState->characterStore->getLocal(posKey);
          if (localCharOpt.has_value()) {
            canonicalChar = localCharOpt.value();
            foundInStore = true;
          }
        }
      }
      
      // Fallback to block sequences if not found in stores
      if (!foundInStore && blockSequence) {
        state::PositionKey structuralKey = state::PositionKey::create(blockId, currentNucPos, currentGapPos);
        auto blockMapIt = blockRootCharFlatIndices.find(blockId);
        if (blockMapIt != blockRootCharFlatIndices.end()) {
          auto flatIndexIt = blockMapIt->second.find(structuralKey);
          if (flatIndexIt != blockMapIt->second.end()) {
            int64_t flatIndex = flatIndexIt->second;
            if (flatIndex >= 0 && flatIndex < static_cast<int64_t>(blockSequence->size())) {
              canonicalChar = (*blockSequence)[flatIndex];
              foundInStore = true;
            }
          }
        }
      }
      
      if (!foundInStore) {
        canonicalChar = '-';
      }
      
      // Apply orientation if needed
      char finalChar = canonicalChar;
      if (isBlockInverted && canonicalChar != '-' && canonicalChar != 'x') {
        switch (canonicalChar) {
          case 'A': finalChar = 'T'; break; case 'T': finalChar = 'A'; break;
          case 'C': finalChar = 'G'; break; case 'G': finalChar = 'C'; break;
          case 'a': finalChar = 't'; break; case 't': finalChar = 'a'; break;
          case 'c': finalChar = 'g'; break; case 'g': finalChar = 'c'; break;
          default: finalChar = canonicalChar; break;
        }
      }
      
      result.push_back(finalChar);
    }
  } catch (const std::exception& e) {
    logging::warn("Exception in getCharactersAtBlock for node {}, block {}, position ({}, {}), length {}: {}", 
                  std::string(nodeId), blockId, startNucPos, startGapPos, length, e.what());
    // Return partial result up to the point of failure
  }
  
  return result;
}

// Set character at (blockId, nucPos, gapPos) in nodeId
bool StateManager::setCharAtPosition(std::string_view nodeId, int32_t blockId,
                                     int32_t nucPos, int32_t gapPos, char c) {
  try {
    // OPTIMIZATION: Use string_view as much as possible, convert only when necessary
    const NodeState* cachedNodeState = tlsNodeCache.getCachedStateByView(nodeId, nodeStates);
    if (!cachedNodeState) {
      // Need to create node state, convert string_view to string only now
      std::string nodeIdStr(nodeId);
      auto &nodeState = getNodeState(nodeIdStr);
      return setCharAtPositionImpl(nodeIdStr, nodeState, blockId, nucPos, gapPos, c);
    }
    
    // We have cached state, but need mutable reference for modification
    std::string nodeIdStr(nodeId);
    auto &nodeState = getNodeState(nodeIdStr);
    return setCharAtPositionImpl(nodeIdStr, nodeState, blockId, nucPos, gapPos, c);
    
  } catch (const std::exception& e) {
    logging::warn("Error in setCharAtPosition({}:{}:{}:{}): {}", 
                std::string(nodeId), blockId, nucPos, gapPos, e.what());
    return false;
  }
}

// Helper method to reduce code duplication
bool StateManager::setCharAtPositionImpl(const std::string& nodeIdStr, NodeState& nodeState,
                                        int32_t blockId, int32_t nucPos, int32_t gapPos, char c) {
    bool is_debug_node = (nodeIdStr == "node_2");
    
    // Create position key
    state::PositionKey posKey{blockId, nucPos, gapPos};
    // Create combined key for caching
    NodePositionKey cacheKey{nodeIdStr, posKey};
    
    // Make sure characterStore exists
    if (!nodeState.characterStore) {
      // Critical issue - characterStore doesn't exist for this node!
      logging::warn("Creating missing characterStore for node {}", nodeIdStr);
      
      // Look up parent's characterStore to inherit from
      std::string parentId = nodeState.parentId;
      std::shared_ptr<HierarchicalCharacterStore> parentCharStore = nullptr;
      
      if (!parentId.empty()) {
        try {
          const auto& parentState = getNodeState(parentId);
          if (parentState.characterStore) {
            parentCharStore = parentState.characterStore;
          }
        } catch (const std::exception& e) {
          logging::err("Error accessing parent state during character store creation: {}", e.what());
        }
      }
      
      // Create a new characterStore, inheriting from parent if possible
      if (parentCharStore != nullptr) {
        nodeState.characterStore = std::make_shared<HierarchicalCharacterStore>(parentCharStore);
        logging::warn("Created characterStore for node {} with parent inheritance", nodeIdStr);
        
        // Note: Materialization is deferred until explicitly called via materializeNodeState()
        // This ensures proper parent-child inheritance order
      } else {
        nodeState.characterStore = std::make_shared<HierarchicalCharacterStore>();
        logging::warn("Created characterStore for node {} without parent", nodeIdStr);
      }
    }
    
    // Store character in the node's local store and directly in nodeHierarchy
    bool storeSuccess = false;
    if (nodeState.characterStore) {
      // Store in node's local store with position key
      storeSuccess = nodeState.characterStore->set(posKey, c);
      
      // Also store in nodeHierarchy to ensure proper lookups
      auto hierIt = nodeHierarchy.find(nodeIdStr);
      if (hierIt != nodeHierarchy.end()) {
        if (!hierIt->second.characterStore) {
          hierIt->second.characterStore = std::make_shared<HierarchicalCharacterStore>();
          logging::debug("Created missing hierarchy store for node {}", nodeIdStr);
        }
        hierIt->second.characterStore->set(posKey, c);
      }
      
      if (is_debug_node) {
        // Verify we can retrieve the character directly from the local store
        auto localValue = nodeState.characterStore->getLocal(posKey);
        if (localValue) {
          if (*localValue == c) {
            logging::debug("Verified local character storage for node_2 at {}:{}:{}: '{}'", 
                         blockId, nucPos, gapPos, c);
          } else {
            logging::warn("Local character mismatch after storage in node_2! Stored '{}', got '{}'", 
                         c, *localValue);
          }
        } else {
          logging::warn("Failed to verify local character storage in node_2 at {}:{}:{}", 
                       blockId, nucPos, gapPos);
        }
        
        // Also verify the hierarchical lookup works correctly
        auto hierarchicalValue = nodeState.characterStore->get(posKey);
        if (!hierarchicalValue || *hierarchicalValue != c) {
          logging::warn("Hierarchical lookup failed in node_2 at {}:{}:{}. Expected '{}', got '{}'",
                       blockId, nucPos, gapPos, c, 
                       hierarchicalValue ? *hierarchicalValue : '?');
        }
        
        // Check hierarchy store if available
        if (hierIt != nodeHierarchy.end() && hierIt->second.characterStore) {
          auto hierValue = hierIt->second.characterStore->get(posKey);
          if (!hierValue || *hierValue != c) {
            logging::warn("Hierarchy store lookup failed for node_2 at {}:{}:{}. Expected '{}', got '{}'",
                        blockId, nucPos, gapPos, c,
                        hierValue ? *hierValue : '?');
          }
        }
      }
    } else {
      logging::err("Failed to store character in characterStore for node {} - it's still NULL!", nodeIdStr);
      return false;
    }
    
    if (!storeSuccess && is_debug_node) {
      logging::warn("Failed to store character '{}' at position {}:{}:{} in node_2", 
                  c, blockId, nucPos, gapPos);
    }
    
    return true;
}

// Get range in scalar global coordinates for a block
CoordRange StateManager::getBlockRange(int32_t blockId) const {
  auto blockRangeIt = blockRanges.find(blockId);
  if (blockRangeIt != blockRanges.end()) {
    return blockRangeIt->second;
  }
  
  // Add detailed error logging to debug why we can't find block ranges
  std::stringstream ss;
  ss << "Block range not found for block ID: " << blockId << ". ";
  ss << "Total blocks in map: " << blockRanges.size() << ". ";
  
  // Log some of the existing block IDs for context
  ss << "Available blocks: ";
  int i = 0;
  for (const auto& [id, range] : blockRanges) {
    ss << id << " ";
    if (++i >= 5) { // Limit logging to first 5 blocks
      ss << "...";
      break;
    }
  }
  
  // Add initialization information to help diagnose the issue
  ss << "\nIMPORTANT: Block ranges are populated during initialization. ";
  ss << "If this error occurs, the most likely cause is that block ranges have not been fully ";
  ss << "initialized yet, or that the initialization sequence is incorrect. ";
  ss << "Make sure that setBlockRange has been called for all blocks before attempting to access them.";
  
  // Log the error with details
  logging::err("{}", ss.str());
  
  throw std::runtime_error(ss.str());
}

// Set range in scalar global coordinates for a block
void StateManager::setBlockRange(int32_t blockId, const CoordRange &range) {
  // Static counter to limit logging globally
  static std::atomic<int> set_range_logs_count = 0;
  const int MAX_SET_RANGE_LOGS = 5;

  // Store the block range
  blockRanges[blockId] = range;
  
  // Update global position cache
  if (blockId >= 0 && blockId < static_cast<int32_t>(globalPosCache.blockStartOffsets.size())) {
    globalPosCache.blockStartOffsets[blockId] = range.start;
  }

  // Log the block range setting, but limit frequency
  bool log_this_set = (set_range_logs_count.load() < MAX_SET_RANGE_LOGS);
  if (log_this_set) {
      logging::debug("Set block {} range to [{}, {}), size: {}", blockId,
                     range.start, range.end, range.end - range.start);
      set_range_logs_count++;
  }
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
    
    // For nucleotide mutations, expand upstream by k non-gap positions
    if (!isBlockMutation && kmerSize > 0) {
      try {
        std::string strNodeId(nodeId);
        
        // Get node state for fast gap checking
        const NodeState* nodeStatePtr = nullptr;
        {
          auto it = nodeStates.find(strNodeId);
          if (it != nodeStates.end()) {
            nodeStatePtr = &it->second;
          }
        }
        
        if (nodeStatePtr) {
          // Calculate global position of the mutation
          int64_t globalMutationPos = blockRange.start + pos;
          
          // Find k non-gap positions upstream from the mutation
          int64_t targetStart = globalMutationPos - kmerSize;
          int64_t currentPos = globalMutationPos - 1;
          int positionsNeeded = kmerSize;
          int iterations = 0;
          
          // Work backwards from the mutation position to account for gaps
          while (positionsNeeded > 0 && currentPos >= blockRange.start && iterations < kmerSize * 100) {
            iterations++;
            
            // Check if current position is a gap using fast gap run checking
            bool isGap = false;
            try {
              isGap = nodeStatePtr->isGapRun(blockId, currentPos - blockRange.start);
            } catch (const std::exception& e) {
              isGap = false; // If we can't check for gaps, assume it's not a gap
            }
            
            if (!isGap) {
              // This position counts toward our k positions
              positionsNeeded--;
              targetStart = currentPos;
            }
            // If it's a gap, we need to go further back but don't count it
            
            currentPos--;
          }
          
          // Ensure we don't go below the block start
          targetStart = std::max(targetStart, blockRange.start);
          
          // Expand the result range upstream
          result.start = std::min(result.start, targetStart);
        }
      } catch (const std::exception& e) {
        // Fallback to naive expansion if gap-aware expansion fails
        logging::debug("Gap-aware expansion failed for nucleotide mutation, using naive expansion: {}", e.what());
      }
    }
    
    if (isBlockMutation) {
      auto activeBlocks = getActiveBlockRanges(nodeId);
      
      int currentPos = -1;
      for (size_t i = 0; i < activeBlocks.size(); ++i) {
        if (activeBlocks[i].first == blockId) {
          currentPos = static_cast<int>(i);
          break;
        }
      }
      
      // For block inversions, we need to recompute the entire block
      // since k-mers can be affected throughout the block
      bool isInversion = isBlockMutation && !isBlockDeactivation && blockActive;
      if (isInversion) {
        // Use the full block range for inversions
        result.start = blockRange.start;
        result.end = blockRange.end;
        
        // Log the expanded range for inversions
        logging::debug("Expanded recomp range for block inversion: [{}, {})",
                     result.start, result.end);
      }
      else if (currentPos >= 0 || isBlockDeactivation) {
        if (isBlockDeactivation) {
          result.start = std::max(result.start, blockRange.start);
          result.end = std::min(result.end, blockRange.end);
        }
      }
    }
    
    // Properly wrap the result in std::optional before returning
    return {result}; // Use brace initialization
}

bool StateManager::applyBlockMutation(std::string_view nodeId, int32_t blockId,
                                      bool isInsertion, bool isInversion_flag) {
  std::string strNodeId(nodeId);
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // Acquire lock
  
  // CRITICAL FIX: Verify that blockRanges has entries before proceeding
  if (blockRanges.empty()) {
    logging::err("CRITICAL ERROR: blockRanges map is empty during applyBlockMutation for node {} block {}!", 
                strNodeId, blockId);
    return false;
  }

  // Targetted debugging for node_1004 and block 703
  bool is_target_debug_node = (strNodeId == "node_1004"); // MODIFIED
  bool is_target_debug_block = (blockId == 703); // MODIFIED
  bool enable_specific_debug = is_target_debug_node && is_target_debug_block;

  if (enable_specific_debug) {
    logging::debug("APPLY_BM_DEBUG node_1004 (block {}): ENTER. isInsertion={}, isInversion_flag={}", blockId, isInsertion, isInversion_flag);
  }

  auto nodeStateIt = nodeStates.find(strNodeId);
  if (nodeStateIt == nodeStates.end()) {
    logging::err("APPLY_BM_DEBUG: Node {} not found in nodeStates.", strNodeId);
    // lock is released by RAII
    return false; // Or throw, depending on desired error handling
  }
  NodeState& nodeState = nodeStateIt->second;

  bool wasOn = isBlockOn_unsafe(strNodeId, blockId); // Use _unsafe
  bool wasInverted = wasOn && isBlockInverted_unsafe(strNodeId, blockId); // Use _unsafe

  if (enable_specific_debug) {
    logging::debug("APPLY_BM_DEBUG node_1004 (block {}): wasOn={}, wasInverted={}", blockId, wasOn, wasInverted);
  }

  // Determine the nature of the operation based on input flags and current state
  bool intendsDeactivation = !isInsertion && !isInversion_flag;
  bool intendsActivationOrOrientationChangeViaInsert = isInsertion;
  bool intendsOrientationToggle = !isInsertion && isInversion_flag;

  bool effectiveActivation = intendsActivationOrOrientationChangeViaInsert && !wasOn;
  bool effectiveDeactivation = intendsDeactivation && wasOn;
  bool effectiveOrientationChange = 
      (intendsActivationOrOrientationChangeViaInsert && wasOn && wasInverted != isInversion_flag) || 
      (intendsOrientationToggle && wasOn);

  if (enable_specific_debug) {
    logging::debug("APPLY_BM_DEBUG node_1004 (block {}): intendsDeactivation={}, effectiveDeactivation={}", blockId, intendsDeactivation, effectiveDeactivation);
  }
  
  auto blockRangeIt = blockRanges.find(blockId);
  if (blockRangeIt == blockRanges.end()) {
    logging::err("Block {} range not found for block mutation in node {}", blockId, strNodeId);
    // lock is released by RAII
    return false; // Or throw
  }
  auto blockRange = blockRangeIt->second;
  int16_t k = getKmerSize(); 

  // Determine and add recomputation ranges to nodeState.recompRanges
  if (effectiveActivation || effectiveOrientationChange) {
    // Block OFF to ON: recompute whole new block range (no expansion needed here)
    // The upstream k non-gap logic will be handled separately below
    if (blockRange.start < blockRange.end) {
      mergeRangeWithExisting(nodeState.recompRanges, blockRange);
    }
  } else if (effectiveDeactivation) {
    // Block ON to OFF: recompute k non-gap range at end of previous ON block
    // This is entirely handled by the upstream recomputation logic below
    // No ranges from the deleted block itself
  }
  

  // Add upstream recomputation range for block insertions and deletions
  // This implements the critical logic: k non-gap characters upstream in previous ON block
  if (effectiveActivation || effectiveDeactivation || effectiveOrientationChange) {
    // DEADLOCK FIX: Cannot call getActiveBlockRanges here because it will try to acquire the lock again
    // Instead, manually collect active blocks since we already hold the lock
    
    // CRITICAL FIX: We need the PARENT state to determine upstream blocks, not the current node state
    // The current node state may have already been modified by previous mutations
    std::vector<std::pair<int32_t, coordinates::CoordRange>> activeBlocks;
    
    // Get parent's active blocks (before mutations) for upstream recomputation logic
    std::string parentId = nodeState.parentId;
    if (!parentId.empty()) {
      // Use parent's block state to determine which blocks are upstream
      for (int32_t testBlockId = 0; testBlockId < static_cast<int32_t>(numBlocks); ++testBlockId) {
        if (isBlockOn_unsafe(parentId, testBlockId)) {
          auto rangeIt = blockRanges.find(testBlockId);
          if (rangeIt != blockRanges.end()) {
            activeBlocks.emplace_back(testBlockId, rangeIt->second);
          }
        }
      }
    }
    
    // Sort by global position for consistent ordering
    std::sort(activeBlocks.begin(), activeBlocks.end(), 
              [](const auto& a, const auto& b) {
                return a.second.start < b.second.start;
              });
    
    // Find the previous active block (upstream from current block)
    int32_t previousActiveBlockId = -1;
    coordinates::CoordRange previousActiveRange;
    
    for (const auto& [activeBlockId, activeRange] : activeBlocks) {
      if (activeRange.end <= blockRange.start) {
        // This block is upstream from our target block
        if (previousActiveBlockId == -1 || activeRange.end > previousActiveRange.end) {
          // This is either the first upstream block we found, or it's closer to our target
          previousActiveBlockId = activeBlockId;
          previousActiveRange = activeRange;
        }
      }
    }
    
    // Add upstream recomputation range if we found a previous active block
    if (previousActiveBlockId != -1) {
      // Lock-safe gap-aware expansion: find the last k non-gap characters of the previous ON block
      // This is a simplified version of expandRecompRanges that doesn't acquire additional locks
      
      // Get the gap map for this node (this doesn't acquire locks)
      const NodeState* nodeStatePtr = nullptr;
      {
        auto it = nodeStates.find(strNodeId);
        if (it != nodeStates.end()) {
          nodeStatePtr = &it->second;
        }
      }
      
      if (nodeStatePtr) {
        // Simplified gap-aware expansion: expand from the end of the block backwards by k positions,
        // but account for gaps by expanding further when gaps are encountered
        int64_t targetStart = previousActiveRange.end - k;
        int64_t currentPos = previousActiveRange.end - 1;
        int positionsNeeded = k;
        int iterations = 0;
        
        // Work backwards from the end of the block to account for gaps
        while (positionsNeeded > 0 && currentPos >= previousActiveRange.start && iterations < k * 100) {
          iterations++;
          
          // Check if current position is a gap using fast gap run checking
          bool isGap = false;
          try {
            isGap = nodeStatePtr->isGapRun(previousActiveBlockId, currentPos - previousActiveRange.start);
          } catch (const std::exception& e) {
            // If we can't check for gaps, assume it's not a gap
            isGap = false;
          }
          
          if (!isGap) {
            // This position counts toward our k positions
            positionsNeeded--;
            targetStart = currentPos;
          }
          // If it's a gap, we need to go further back but don't count it
          
          currentPos--;
        }
        
        // Ensure we don't go below the block start
        targetStart = std::max(targetStart, previousActiveRange.start);
        
        // Create the upstream recomputation range
        coordinates::CoordRange upstreamRange;
        upstreamRange.start = targetStart;
        upstreamRange.end = previousActiveRange.end;
        
        if (upstreamRange.start < upstreamRange.end) {
          mergeRangeWithExisting(nodeState.recompRanges, upstreamRange);
          
          if (enable_specific_debug) {
            logging::debug("APPLY_BM_DEBUG: Added upstream recomputation range [{}, {}) from previous block {} for block {} mutation (gap-aware expansion)",
                           upstreamRange.start, upstreamRange.end, previousActiveBlockId, blockId);
          }
        }
      } else {
        // Fallback to naive expansion when nodeState is not available
        logging::warn("APPLY_BM_DEBUG: Could not get node state for node {}: Using naive k-position expansion.", strNodeId);
        coordinates::CoordRange naiveUpstreamRange;
        naiveUpstreamRange.start = std::max(previousActiveRange.end - k, previousActiveRange.start);
        naiveUpstreamRange.end = previousActiveRange.end;
        if (naiveUpstreamRange.start < naiveUpstreamRange.end) {
          mergeRangeWithExisting(nodeState.recompRanges, naiveUpstreamRange);
        }
      }
    } else if (enable_specific_debug) {
      logging::debug("APPLY_BM_DEBUG: No upstream active block found for block {} mutation", blockId);
    }
  }

  // Apply the actual block mutation state changes - using inheritance approach
  if (isInsertion) { 
    // Always set an explicit state for this block to ensure proper inheritance
    nodeState.setExplicitBlockState(blockId, true);
    nodeState.setBlockForward(blockId, !isInversion_flag); 
    if (enable_specific_debug) {
        logging::debug("APPLY_BM_DEBUG node_1004 (block {}): Called setExplicitBlockState(true), setBlockForward({})", blockId, !isInversion_flag);
    }
  } else { 
    if (intendsDeactivation) { 
      // For block deletions, always set an explicit state, even if the block was already off
      if (enable_specific_debug) {
        std::cerr << "[DEBUG] APPLY_BM_DEBUG node_1004 (block " << blockId << "): Intends deactivation. Clearing seeds and calling setExplicitBlockState(false)." << std::endl;
      }
      auto it = blockToSeeds.find(blockId);
      if (it != blockToSeeds.end()) {
        it->second.clear();
      }
      
      // CRITICAL: Always set explicit OFF state to ensure proper inheritance
      nodeState.setExplicitBlockState(blockId, false);
      
      if (enable_specific_debug) {
        bool checkOn = isBlockOn_unsafe(strNodeId, blockId); // Re-check with _unsafe
        std::cerr << "[DEBUG] APPLY_BM_DEBUG node_1004 (block " << blockId << "): After setExplicitBlockState(false), isBlockOn_unsafe is " << (checkOn ? "true" : "false") << std::endl;
      }
    } else if (intendsOrientationToggle) { 
      if (wasOn) { 
        bool newInvertedState = !wasInverted;
        // Set explicit block state to true to ensure we're not relying on inheritance for this block
        nodeState.setExplicitBlockState(blockId, true);
        nodeState.setBlockForward(blockId, !newInvertedState);
        if (enable_specific_debug) {
            logging::debug("APPLY_BM_DEBUG node_1004 (block {}): Toggled orientation. New inverted state: {}", blockId, newInvertedState);
        }
      }
    }
  }

  // Update gap runs when block state changes for fast gap detection
  if (effectiveActivation) {
    // Block turned ON: need to update gap runs for this block
    // For now, clear any existing gap runs for this block and they will be rebuilt on demand
    nodeState.gapRuns.erase(blockId);
  } else if (effectiveDeactivation) {
    // Block turned OFF: clear gap runs for this block since it's no longer active
    nodeState.gapRuns.erase(blockId);
  }

  resetNodeCache(nodeId); // Stays as is, uses string_view
  // lock is released by RAII
  return true;
}

// Helper method to initialize hierarchical connections without explicit propagation
void StateManager::propagateState(const std::string &fromNode, const std::string &toNode) {
  std::unique_lock<std::shared_mutex> lock(nodeMutex); 
  static bool& verbose_logging = getPropLoggingRef(); // Ensure this is accessible or pass as param
  
  
  // Check if parent node exists
  if (nodeStates.find(fromNode) == nodeStates.end()) {
    logging::err("Parent node {} not found for child {}", fromNode, toNode);
    throw std::runtime_error("Parent node " + fromNode + " not found");
  }
  
  // Create child node if it doesn't exist
  auto childStateIt = nodeStates.find(toNode);
  if (childStateIt == nodeStates.end()) {
    nodeStates.emplace(toNode, NodeState());
    childStateIt = nodeStates.find(toNode); // Re-find after emplace
    
  }
  
  auto &parentState = nodeStates.at(fromNode); 
  auto &childState = childStateIt->second;   

  // Set the parent-child relationship for inheritance
  childState.parentId = fromNode;
  
  // OPTIMIZATION: Materialize block state from parent inheritance
  childState.materializeBlockState(&parentState);
  
  // Legacy: Copy parent's active blocks to child (inheritance) for backward compatibility
  childState.activeBlocks = parentState.activeBlocks;
  
  // Inherit parent's total seed count
  childState.initializeSeedCountFromParent(parentState.getTotalSeedCount());
  
  // Set up hierarchical gap map
  if (parentState.gapMap) {
    childState.gapMap = std::make_shared<HierarchicalGapMap>(parentState.gapMap);
  } else {
    childState.gapMap = std::make_shared<HierarchicalGapMap>(rootGapMap);
  }
  
  // Set up hierarchical character store  
  if (parentState.characterStore) {
    childState.characterStore = std::make_shared<HierarchicalCharacterStore>(parentState.characterStore);
  } else {
    static std::shared_ptr<HierarchicalCharacterStore> localRootCharacterStore;
    if (!localRootCharacterStore) { localRootCharacterStore = std::make_shared<HierarchicalCharacterStore>(); }
    childState.characterStore = std::make_shared<HierarchicalCharacterStore>(localRootCharacterStore);
  }
  
  // Inherit parent's gap runs for fast gap detection
  childState.gapRuns = parentState.gapRuns;
  
  // Set up hierarchical seed store
  auto hierarchyIt_child = nodeHierarchy.find(toNode);
  if (hierarchyIt_child != nodeHierarchy.end()) {
    auto parentHierIt = nodeHierarchy.find(fromNode);
    if (parentHierIt != nodeHierarchy.end() && parentHierIt->second.seedStore) {
      hierarchyIt_child->second.seedStore = std::make_shared<HierarchicalSeedStore>(parentHierIt->second.seedStore);
    } else {
      static std::shared_ptr<HierarchicalSeedStore> localRootSeedStore;
      if (!localRootSeedStore) { localRootSeedStore = std::make_shared<HierarchicalSeedStore>(); }
      hierarchyIt_child->second.seedStore = std::make_shared<HierarchicalSeedStore>(localRootSeedStore);
    }
  }
  
  // Materialize character and seed data for fast indexing access
  try {
    const NodeState* parentStatePtr = &parentState;
    const HierarchicalCharacterStore* childCharStore = childState.characterStore.get();
    
    // Use empty root characters as fallback
    static const absl::flat_hash_map<PositionKey, char, PositionKey::Hash> emptyRootChars;
    
    // Note: Materialization is deferred until explicitly called via materializeNodeState()
    // This ensures proper parent-child inheritance order and avoids premature materialization
  } catch (const std::exception& e) {
    logging::warn("Failed to materialize character data during propagation from {} to {}: {}", 
                  fromNode, toNode, e.what());
  }
  
  resetNodeCache(toNode);
  
}

void StateManager::applyNucleotideMutation(const std::string& nodeId, NodeState& nodeState,
                                         int32_t blockId, int32_t nucPos, int32_t gapPos, char value) {

  

  
  // Store the original value directly. Complementation will be handled by getCharAtPosition.
  char storeValue = value; 

  // Create position key for lookup
  PositionKey posKey{blockId, nucPos, gapPos};
  // Create combined key for caching
  NodePositionKey cacheKey{nodeId, posKey};
  
  // Make sure the character store exists
  if (!nodeState.characterStore) {
    // Critical issue - characterStore doesn't exist for this node!
    logging::warn("Creating missing characterStore during applyNucleotideMutation for node {}", nodeId);
    
    
    // Look up parent's characterStore to inherit from
    std::string parentId = nodeState.parentId;
    std::shared_ptr<HierarchicalCharacterStore> parentCharStore = nullptr;
    
    if (!parentId.empty()) {
      try {
        const auto& parentState = getNodeState(parentId);
        if (parentState.characterStore) {
          parentCharStore = parentState.characterStore;
          
        }
      } catch (const std::exception& e) {
        logging::err("Error accessing parent state during character store creation: {}", e.what());
        
      }
    }
    
    // Create a new characterStore, inheriting from parent if possible
    if (parentCharStore != nullptr) {
      nodeState.characterStore = std::make_shared<HierarchicalCharacterStore>(parentCharStore);
      logging::warn("Created characterStore for node {} with parent inheritance", nodeId);
      
     
    } else {
      nodeState.characterStore = std::make_shared<HierarchicalCharacterStore>();
      logging::warn("Created characterStore for node {} without parent", nodeId);
      
    }
  }
  
  // CRITICAL FIX: Directly store mutations in the local character store and hierarchy store
  // This ensures they're found when later looking them up
  bool directStoreSuccess = false;
  bool hierStoreSuccess = false;
  
  // First store in the node's local store
  if (nodeState.characterStore) {
    directStoreSuccess = nodeState.characterStore->set(posKey, storeValue);
    
    // Update materialized state for performance
    nodeState.setMaterializedCharacter(posKey, storeValue);
  }
  
  // Also store in hierarchy store to ensure lookups work properly
  auto hierIt = nodeHierarchy.find(nodeId);
  if (hierIt != nodeHierarchy.end()) {
    if (!hierIt->second.characterStore) {
      hierIt->second.characterStore = std::make_shared<HierarchicalCharacterStore>();
      
    }
    
    hierStoreSuccess = hierIt->second.characterStore->set(posKey, storeValue);
    
  }
  
  // Then also use setCharAtPosition to ensure consistent behavior
  bool setResult = setCharAtPosition(nodeId, blockId, nucPos, gapPos, storeValue);
  
  if (!setResult) {
    logging::err("Failed to set character at position for node {} at block {}:{}:{}", 
                nodeId, blockId, nucPos, gapPos);
                
   
  }
}

// Helper to set gap list length
// This sets the number of positions in the gap list for a specific nucleotide
// Note: This is NOT related to gap characters ('-') in the sequence tracked by
// gap maps. The gap list contains optional nucleotides that:
// - In normal blocks: come BEFORE the main nucleotide in traversal order
// - In inverted blocks: come AFTER the main nucleotide in traversal order (and
// in reverse sequence)
//
// IMPORTANT: Gap list lengths are a fixed structural property of the tree defined
// by the reference sequence. They are the same for all nodes in the tree and
// should only be set during initialization, never modified per-node.
void StateManager::setGapListLength(int32_t blockId, int32_t nucPos, size_t length) {
  // Static counters to limit logging globally
  static std::atomic<int> resize_logs_count = 0;
  static std::atomic<int> set_logs_count = 0;
  const int MAX_LOGS_PER_TYPE = 5;

  // Check for valid block ID
  if (blockId < 0) {
    logging::warn("Invalid block ID {} in setGapListLength", blockId);
    return;
  }
  
  // Ensure the outer array has enough space
  if (static_cast<size_t>(blockId) >= gapListLengthArray.size()) {
    // Log this resize only once, as it's less frequent and indicates overall sizing
    static std::once_flag outer_resize_flag;
    std::call_once(outer_resize_flag, [&](){
        logging::debug("Resizing outer gapListLengthArray once from {} to at least {}", 
                      gapListLengthArray.size(), blockId + 1);
    });
    gapListLengthArray.resize(blockId + 1);
  }
  
  // Only resize the inner array if necessary - this is key for memory efficiency
  if (static_cast<size_t>(nucPos) >= gapListLengthArray[blockId].size()) {
    // If we need to store a non-zero value, resize just this block's array
    if (length > 0) {
      // Limit logging for inner array resize
      bool log_this_resize = (resize_logs_count.load() < MAX_LOGS_PER_TYPE);
      if (log_this_resize) {
          logging::debug("Resizing inner array for block {} from {} to {}", 
                        blockId, gapListLengthArray[blockId].size(), nucPos + 1);
          resize_logs_count++;
      }
      gapListLengthArray[blockId].resize(nucPos + 1, 0);
    } else {
      // If length is 0, we don't need to resize since the default is 0
      return;
    }
  }
  
  // Store the value
  gapListLengthArray[blockId][nucPos] = length;
  
  // Log significant changes to help with debugging, but limit frequency
  if (length > 0) {
    bool log_this_set = (set_logs_count.load() < MAX_LOGS_PER_TYPE);
    if (log_this_set) {
        logging::debug("Set gap list length for block {} pos {} to {}", 
                      blockId, nucPos, length);
        set_logs_count++;
    }
  }
}

// Get gap list length for a position
size_t StateManager::getGapListLength(int32_t blockId, int32_t nucPos) const {
  // Check for valid block ID
  if (blockId < 0 || static_cast<size_t>(blockId) >= gapListLengthArray.size()) {
    return 0; // No gap list for invalid block
  }
  
  // Efficient early return if the nucPos is outside the allocated range
  if (static_cast<size_t>(nucPos) >= gapListLengthArray[blockId].size()) {
    return 0; // Position outside the allocated range
  }
  
  // Return the stored value (which will be 0 if not set)
  return gapListLengthArray[blockId][nucPos];
}

// Check if a position is tracked in the gap list for a specific block
// Note: This does NOT check if the character at this position is a gap ("-"),
// but rather if the position is tracked in the block's gap list data structure
bool StateManager::isGapAt(int32_t blockId, int64_t blockPos) const {
  if (blockId < 0) {
    logging::err("Negative block ID {} in isGapAt", blockId);
    return false;
  }
  
  // Look up the gap list for this block
  auto it = blockIdToGapList.find(blockId);
  if (it == blockIdToGapList.end()) {
    return false; // No gap list for this block
  }
  
  // Check if the position exists in the gap list
  auto& gapList = it->second;
  return std::binary_search(gapList.begin(), gapList.end(), blockPos);
}

void StateManager::addGap(int32_t blockId, int64_t blockPos) {
  if (blockId < 0) {
    logging::err("Negative block ID {} in addGap", blockId);
    return;
  }
  
  // Get or create the gap list for this block
  auto& gapList = blockIdToGapList[blockId];
  
  // Check if the gap is already registered
  auto it = std::lower_bound(gapList.begin(), gapList.end(), blockPos);
  if (it != gapList.end() && *it == blockPos) {
    logging::debug("Gap already exists at block {} position {}", blockId, blockPos);
    return;
  }
  
  // Insert the gap position (maintaining sorted order)
  gapList.insert(it, blockPos);
  
  // Update the gap list length - use 3-parameter version
  setGapListLength(blockId, 0, gapList.size());
  
  logging::debug("Added gap at block {} position {}, total gaps: {}", 
               blockId, blockPos, gapList.size());
}

void StateManager::removeGap(int32_t blockId, int64_t blockPos) {
  if (blockId < 0) {
    logging::err("Negative block ID {} in removeGap", blockId);
    return;
  }
  
  // Look up the gap list for this block
  auto it = blockIdToGapList.find(blockId);
  if (it == blockIdToGapList.end()) {
    logging::debug("No gaps registered for block {} in removeGap", blockId);
    return;
  }
  
  // Find and remove the gap position
  auto& gapList = it->second;
  auto pos = std::lower_bound(gapList.begin(), gapList.end(), blockPos);
  if (pos != gapList.end() && *pos == blockPos) {
    gapList.erase(pos);
    
    // Update the gap list length - use 3-parameter version
    setGapListLength(blockId, 0, gapList.size());
    
    logging::debug("Removed gap at block {} position {}, total gaps: {}", 
                  blockId, blockPos, gapList.size());
  } else {
    logging::debug("Gap at block {} position {} not found in removeGap", blockId, blockPos);
  }
}

// Initialize the optimized gap list length array - now just resize inner arrays
void StateManager::initializeGapListLengthArray(size_t numBlocks, size_t providedMaxNucPos) {
  // Validate inputs
  if (numBlocks == 0) {
    logging::err("Cannot initialize gap list length array with 0 blocks");
    return;
  }
  
  // Log initial array state
  logging::debug("GAPLIST: Initializing gap list array with {} blocks, max nucPos hint {}", 
              numBlocks, providedMaxNucPos);
  logging::debug("GAPLIST: Current array state - outer size: {}", gapListLengthArray.size());
  
  if (!gapListLengthArray.empty()) {
    for (size_t i = 0; i < std::min(size_t(5), gapListLengthArray.size()); i++) {
      logging::debug("GAPLIST: Block {} inner array size: {}", i, gapListLengthArray[i].size());
    }
  }
  
  // Ensure outer array has the right size (should already be done in setNumBlocks)
  if (gapListLengthArray.size() != numBlocks) {
    logging::debug("GAPLIST: Resizing outer array from {} to {} blocks", 
                gapListLengthArray.size(), numBlocks);
    gapListLengthArray.resize(numBlocks);
  }
  
  // Instead of resizing all inner arrays upfront, we'll keep track of which ones 
  // actually have non-zero values so far
  size_t nonEmptyInnerArrays = 0;
  size_t totalNonZeroEntries = 0;
  
  // Check current state instead of blindly resizing everything
  for (size_t blockId = 0; blockId < numBlocks; ++blockId) {
    bool hasNonZeroValues = false;
    for (size_t j = 0; j < gapListLengthArray[blockId].size(); ++j) {
      if (gapListLengthArray[blockId][j] > 0) {
        hasNonZeroValues = true;
        totalNonZeroEntries++;
      }
    }
    
    if (hasNonZeroValues) {
      nonEmptyInnerArrays++;
    }
  }
  
  logging::debug("GAPLIST: Found {} blocks with non-zero gap lists containing {} total entries", 
              nonEmptyInnerArrays, totalNonZeroEntries);
  
  // No need to resize anything here - we'll do it on-demand in setGapListLength
  // This significantly reduces memory usage by not allocating large arrays for blocks
  // that don't need them
  
  // Final state logging
  logging::debug("GAPLIST: Completed gap list array initialization with lazy allocation strategy");
  
  if (!gapListLengthArray.empty()) {
    for (size_t i = 0; i < std::min(size_t(5), gapListLengthArray.size()); i++) {
      // Log some of the gap list values to check for non-zero entries
      std::stringstream ss;
      ss << "GAPLIST: Block " << i << " sample values: ";
      int value_count = 0; // Counter for logged values
      for (size_t j = 0; j < gapListLengthArray[i].size(); j++) {
        if (gapListLengthArray[i][j] > 0) {
          ss << "[" << j << "]=" << gapListLengthArray[i][j] << " ";
          if (++value_count >= 5) { // Limit logging to first 5 non-zero values
              ss << "...";
              break;
          }
        }
      }
      logging::debug("{}", ss.str());
    }
  }
  
  logging::debug("Successfully initialized gap list outer array with {} blocks", numBlocks);
}

// Internal helper method to get active block ranges efficiently - RENAMED
std::vector<std::pair<int32_t, coordinates::CoordRange>>
StateManager::getActiveBlockRangesImpl(std::string_view nodeId) const {
  // Convert string_view to string once and reuse
  std::string nodeIdStr(nodeId);
  std::vector<std::pair<int32_t, coordinates::CoordRange>> activeBlocks;
  
  // Create a shared lock for thread safety
  std::shared_lock<std::shared_mutex> lock(nodeMutex);
  
  try {
    // Collect all active blocks using inheritance model
    std::unordered_set<int32_t> inheritedActiveBlocks;
    
    // Get total available blocks
    const size_t totalBlocks = numBlocks;
    
    // Check each possible block ID using inheritance-based isBlockOn_unsafe
    for (int32_t blockId = 0; blockId < static_cast<int32_t>(totalBlocks); ++blockId) {
      if (isBlockOn_unsafe(nodeId, blockId)) {
        inheritedActiveBlocks.insert(blockId);
      }
    }
    
    // Log active blocks (limit logging to reduce spam)
    static thread_local int callCount = 0;
    bool enableDebugLogging = (callCount++ < 10) || inheritedActiveBlocks.empty();
    
    if (enableDebugLogging) {
      logging::debug("Node {} has {} active blocks (with inheritance)", 
                    nodeIdStr, inheritedActiveBlocks.size());
      
      // Log some active blocks for debugging
      if (!inheritedActiveBlocks.empty()) {
        std::stringstream ss;
        ss << "Active blocks for node " << nodeIdStr << ": ";
        int count = 0;
        for (int32_t blockId : inheritedActiveBlocks) {
          ss << blockId << " ";
          if (++count >= 5) { // Limit logging to first 5 blocks
            ss << "...";
            break;
          }
        }
        logging::debug("{}", ss.str());
      }
      
      // Count total available block ranges for debugging
      logging::debug("Total blocks in map: {}", blockRanges.size());
    }
    
    // Gather active blocks with their ranges
    for (int32_t blockId : inheritedActiveBlocks) {
      auto rangeIt = blockRanges.find(blockId);
      if (rangeIt != blockRanges.end()) {
        activeBlocks.emplace_back(blockId, rangeIt->second);
        
        if (enableDebugLogging && activeBlocks.size() <= 5) {
          logging::debug("Block {} has range [{}, {})", blockId, 
                       rangeIt->second.start, rangeIt->second.end);
        }
      } else {
        logging::warn("Block range not found for block ID: {}. Total blocks in map: {}.", 
                     blockId, blockRanges.size());
      }
    }
    
    // CRITICAL: Sort active blocks by their coordinate position for correct sequence order
    // This ensures that blocks appear in proper genomic order, not arbitrary ID order
    std::sort(activeBlocks.begin(), activeBlocks.end(),
              [](const auto& a, const auto& b) {
                return a.second.start < b.second.start;
              });
    
    if (enableDebugLogging && activeBlocks.size() > 0) {
      logging::debug("Sorted {} active blocks by position for node {}", 
                   activeBlocks.size(), nodeIdStr);
    }
  } catch (const std::exception &e) {
    logging::err("Error getting active block ranges for node {}: {}", nodeIdStr, e.what());
    throw;
  }

  return activeBlocks;
}

// Public method to get active block ranges (non-const version)
std::vector<std::pair<int32_t, coordinates::CoordRange>>
StateManager::getActiveBlockRanges(std::string_view nodeId) {
  // Now calls the renamed public implementation method
  // Need to cast away constness for the non-const call signature, 
  // but the underlying implementation is const-safe.
  return const_cast<StateManager*>(this)->getActiveBlockRangesImpl(nodeId); 
}

// Helper to get node gap map
std::shared_ptr<HierarchicalGapMap>
StateManager::getNodeGapMap(std::string_view nodeId) const {
  // Convert once and reuse
  std::string nodeIdStr(nodeId);
  
  auto it = nodeStates.find(nodeIdStr);
  if (it != nodeStates.end() && it->second.gapMap) {
    return it->second.gapMap;
  }

  // Throw instead of warning and falling back to root gap map
  throw std::runtime_error("No gap map found for node " + nodeIdStr);
}

// Get active blocks for a node
const absl::flat_hash_set<int32_t> &
StateManager::getActiveBlocks(std::string_view nodeId) const {
  // Lock node states for reading
  std::shared_lock<std::shared_mutex> lock(nodeMutex);

  std::string nodeIdStr(nodeId);
  auto it = nodeStates.find(nodeIdStr);
  if (it == nodeStates.end()) {
    // This should not happen if initializeNode was called correctly
    // Log detailed error and return a static empty set
    logging::err("Node state not found for {} in getActiveBlocks. Returning empty set.", nodeIdStr);
    static const absl::flat_hash_set<int32_t> empty_set;
    return empty_set; // Return reference to a static empty set
  }
  return it->second.activeBlocks;
}


// Helper method to merge a new range with existing ranges
void StateManager::mergeRangeWithExisting(std::vector<CoordRange> &ranges,
                                          const CoordRange &newRange) {
  // Use the consolidated RangeOperations function
  coordinates::RangeOperations::mergeRangeWithExisting(ranges, newRange);
}


void StateManager::backtrackNode(const std::string &nodeId) {
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // ADD THIS: Use unique lock for modifications
  // Remove the node state to free memory
  nodeStates.erase(nodeId);

  // No sequence cache to clear anymore - removed

  logging::debug("Backtracked node {}", nodeId);
}

// Get node state reference - use structured bindings for clarity
NodeState &StateManager::getNodeState(const std::string &nodeId) {
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // ADD THIS: Use unique lock for potential modifications
  if (auto it = nodeStates.find(nodeId); it != nodeStates.end()) {
    return it->second;
  }
  // Call initializeNode if the node doesn't exist
  // IMPORTANT: initializeNode must handle its own locking internally, 
  // AND the sequential pre-initialization pass must prevent concurrent calls here.
  // --> Release lock before calling initializeNode <--
  lock.unlock();
  initializeNode(nodeId);
  lock.lock(); // --> Re-acquire lock <--
  
  auto newIt = nodeStates.find(nodeId);
  if (newIt == nodeStates.end()) {
    logging::err("Failed to initialize node state for {}", nodeId);
    throw std::runtime_error("Failed to initialize node state for " + nodeId);
  }
  return newIt->second;
}

// Method to get const reference to a node's state
const NodeState &StateManager::getNodeState(const std::string &nodeId) const {
  std::shared_lock<std::shared_mutex> lock(nodeMutex); // ADD THIS: Use shared lock for read access
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

// Optimized seed management using materialized state only
std::optional<seeding::seed_t> StateManager::getSeedAtPosition(std::string_view nodeId, int64_t pos) const {
    if (pos < 0) {
        logging::warn("Invalid position {} in getSeedAtPosition for node {}", pos, nodeId);
        return std::nullopt;
    }
    
    std::shared_lock<std::shared_mutex> readLock(nodeMutex);
    std::string strNodeId(nodeId);
    
    // Only check the node's materialized seeds (no hierarchical lookups)
    auto nodeStateIt = nodeStates.find(strNodeId);
    if (nodeStateIt == nodeStates.end()) {
        throw std::runtime_error("Node state not found for " + strNodeId);
    }
    
    // Debug: Track inheritance for specific nodes
    static thread_local int debug_call_count = 0;
    // All inherited seeds should already be materialized before mutations are applied
    auto result = nodeStateIt->second.getMaterializedSeed(pos);
    return result; // Return the result directly (might be nullopt if no seed at this position)
}

// Efficiently process seeds in a range to avoid redundant sequence extraction
void StateManager::processSeedsInRange(
    std::string_view nodeId, const CoordRange &range, int k, int s,
    const std::function<void(int64_t startPos, int64_t endPos, std::optional<seeding::seed_t> &seed, std::string_view kmerView)> &processFn) {

  // Extract the entire range at once instead of character by character
  auto [sequence, positions, gaps, endPositions] = extractSequence(nodeId, range, true);
  
  if (sequence.length() < k || positions.size() < k) {
    return; // Not enough characters for a k-mer
  }
  
  // Create a string_view to avoid copying
  std::string_view seqView(sequence);
  
  // Process each potential k-mer position in a single pass
  for (size_t i = 0; i + k <= sequence.length() && i + k <= positions.size(); i++) {
      int64_t startPos = positions[i];
      int64_t endPos = positions[i + k - 1];
      
    // Skip if positions aren't contiguous (might contain gaps we skipped)
    if (endPos - startPos + 1 != k) continue;
    
    // Get or create seed reference - create a temporary that lives through the function call
    std::optional<seeding::seed_t> seedTemp;
    std::optional<seeding::seed_t>* seedPtr;
    
    if (startPos < static_cast<int64_t>(positionSeeds.size())) {
      seedPtr = &positionSeeds[startPos];
    } else {
      seedTemp = std::optional<seeding::seed_t>();
      seedPtr = &seedTemp;
    }
    
    // Process with the callback using a substring view (no copy)
    std::string_view kmerView = seqView.substr(i, k);
    processFn(startPos, endPos, *seedPtr, kmerView);
  }
}

// Helper to check if a position is in a gap
bool StateManager::isGapPosition(std::shared_ptr<HierarchicalGapMap> nodeGapMap, int64_t pos) const {
  // Validate the provided position
  if (pos < 0) {
    throw std::invalid_argument("Negative position " + std::to_string(pos) + " provided to isGapPosition");
  }
  
  // Null check on gap map 
  if (nodeGapMap == nullptr) { // Check against nullptr
    throw std::invalid_argument("Null gap map provided to isGapPosition");
  }
  
  // Directly check if position is in a gap using the provided gap map
  try {
    return nodeGapMap->isGap(pos); // Assumes HierarchicalGapMap has isGap member
  } catch (const std::exception& e) {
    // Rethrow with enhanced error context
    throw std::runtime_error("Failed to check gap status for position " + std::to_string(pos) + 
                            ": " + e.what());
  }
}

// Centralized gap map update function that handles all types of updates
void StateManager::batchUpdateGapMap(const std::string& nodeId, 
                                      const std::vector<coordinates::GapUpdate>& updates) {
  if (updates.empty()) {
    return;
  }
  
  // Use lock for thread safety
  // std::lock_guard<std::mutex> lock(nodeStatesMutex); // REMOVE THIS LINE
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // ADD THIS: Use unique lock for modifications

    auto nodeStateIt = nodeStates.find(nodeId);
  if (nodeStateIt == nodeStates.end()) {
    throw std::runtime_error("Node not found during consolidatedGapUpdate: " + nodeId);
  }
  
  auto& nodeState = nodeStateIt->second;
  if (!nodeState.gapMap) {
    throw std::runtime_error("Node has no gap map during consolidatedGapUpdate: " + nodeId);
  }
  
  // Convert to gap_map::GapUpdate format directly
  std::vector<gap_map::GapUpdate> gapMapUpdates;
  gapMapUpdates.reserve(updates.size());
  
  for (const auto& update : updates) {
    // Create gap_map::GapUpdate (isRemoval, {start, end})
    bool isRemoval = !update.isGapAddition; // REVERTED: Original logic was correct
    int64_t start = update.pos;
    int64_t end = update.pos + update.length - 1;
    gap_map::GapUpdate gapMapUpdate(isRemoval, {start, end});
    gapMapUpdates.push_back(gapMapUpdate);
  }
  
  // Apply updates to node's hierarchical gap map
  size_t updateCount = nodeState.gapMap->applyUpdates(gapMapUpdates);
  logging::debug("Applied {} hierarchical gap updates to node {}", updateCount, nodeId);
  
  // Get all direct children (no need for recursive updates since we use hierarchical gap maps)
  std::vector<std::string> directChildren;
  for (const auto& [childId, childState] : nodeStates) {
    if (childState.parentId == nodeId) {
      directChildren.push_back(childId);
    }
  }
  
  // Reset cache for the current node AND ITS CHILDREN (since gap map changes affect them)
  resetNodeCache(nodeId); // Set resetChildCaches to false - optimization
}

// Initialize the state manager with a tree
void StateManager::initialize(panmanUtils::Tree *tree, size_t maxNucPosHint) {
  if (!tree || !tree->root) {
    logging::err("Invalid tree in initialize");
    return;
  }
  
  // CRITICAL FIX: Log the current state of blockRanges before initialization
  logging::info("Starting StateManager initialization with {} existing block ranges", blockRanges.size());
  
  logging::debug("Initializing state manager with tree");
  
  // Set up the node hierarchy first
  initializeNodeHierarchy(tree, tree->root);
  
  // Initialize DFS indices for traversal
  initializeDfsIndices(tree);
  
  
  // Initialize seed storage hierarchies
  logging::info("Initializing hierarchical seed storage");
  initializeSeedStorage();
  
  // Verify seed stores were properly initialized
  size_t nodesWithSeedStores = 0;
  size_t nodesWithoutSeedStores = 0;
  std::vector<std::string> firstFewNodesWithoutStores;
  
  for (const auto& [nodeId, hierarchy] : nodeHierarchy) {
    if (hierarchy.seedStore) {
      nodesWithSeedStores++;
    } else {
      nodesWithoutSeedStores++;
      if (firstFewNodesWithoutStores.size() < 5) {
        firstFewNodesWithoutStores.push_back(nodeId);
      }
    }
  }
  
  // Log verification results
  if (nodesWithoutSeedStores > 0) {
    logging::err("SEED_VERIFY: Found {} nodes without seed stores after initialization!", 
                nodesWithoutSeedStores);
    for (const auto& nodeId : firstFewNodesWithoutStores) {
      logging::err("SEED_VERIFY: Node without seed store: {}", nodeId);
    }
  } else {
    logging::info("SEED_VERIFY: All {} nodes have seed stores properly initialized", 
                 nodesWithSeedStores);
  }
  
  // Ensure gap list inner arrays are properly sized
  logging::debug("Ensuring gap list arrays are properly sized with hint {}", maxNucPosHint);
  initializeGapListLengthArray(numBlocks, maxNucPosHint);
  
  // After initial sizing, optimize memory usage
  optimizeGapListMemory();
  
  // Get topologically sorted nodes (parents before children)
  std::vector<std::string> sortedNodes = buildTopologicalOrder(tree);
  
  // Initialize nodes in topological order to avoid deadlocks
  int nodesInitialized = 0;
  logging::info("Initializing {} nodes in topological order", sortedNodes.size());
  
  for (const std::string& nodeId : sortedNodes) {
    try {
      // Initialize node - this will not use recursion now
      initializeNode(nodeId);
      nodesInitialized++;
      
      // Log progress periodically
      if (nodesInitialized % 1000 == 0 || nodesInitialized == 1) {
        logging::info("Initialized {}/{} nodes", nodesInitialized, sortedNodes.size());
      }
    } catch (const std::exception& e) {
      logging::err("Failed to initialize node {}: {}", nodeId, e.what());
    }
  }
  
  // globalPosCache will be initialized by a separate call to initializeGlobalCache()
  // after blockRanges and blockCoordToGlobalPos are fully populated by the caller (indexing::initializeStateManager)
  
  logging::info("StateManager structural initialization complete - initialized {} nodes. Global position cache to be built separately.", nodesInitialized);
}

void StateManager::initializeGlobalCache() {
  std::unique_lock<std::shared_mutex> lock(nodeMutex); // Ensure thread safety for cache init
  if (!globalPosCacheInitialized.load(std::memory_order_acquire) && numCoords > 0 && !blockRanges.empty()) {
    logging::debug("Initializing global position cache with {} coordinates and {} blocks", 
                 numCoords, blockRanges.size());

    globalPosCache.posToCoords.resize(numCoords, std::nullopt);
    if (numBlocks > 0) { // Only resize if numBlocks is positive
        globalPosCache.blockStartOffsets.resize(numBlocks, -1);
    } else {
        // If numBlocks is 0, ensure blockStartOffsets is empty or handle appropriately
        globalPosCache.blockStartOffsets.clear();
        logging::debug("Number of blocks is 0, blockStartOffsets will be empty.");
    }
    globalPosCache.maxPosition = numCoords;

    size_t blocksFilled = 0;
    size_t positionsMapped = 0;

    std::unordered_map<int64_t, std::tuple<int32_t, int32_t, int32_t>> tempGlobalToLocalMap;
    tempGlobalToLocalMap.reserve(blockCoordToGlobalPos.size()); // Pre-allocate
    for (const auto& [key, globalPos] : blockCoordToGlobalPos) {
        if (globalPos >= 0 && globalPos < static_cast<int64_t>(numCoords)) { // Ensure globalPos is within bounds
            tempGlobalToLocalMap[globalPos] = std::make_tuple(key.blockId, key.nucPos, key.gapPos);
        }
    }
    
    logging::debug("Collected {} valid initial mappings from blockCoordToGlobalPos for global cache", tempGlobalToLocalMap.size());

    for (int64_t pos = 0; pos < static_cast<int64_t>(numCoords); ++pos) {
        auto it = tempGlobalToLocalMap.find(pos);
        if (it != tempGlobalToLocalMap.end()) {
            globalPosCache.posToCoords[pos] = it->second;
            positionsMapped++;
        }
    }

    for (const auto& [blockId, range] : blockRanges) {
        if (blockId >= 0 && static_cast<size_t>(blockId) < globalPosCache.blockStartOffsets.size()) {
            globalPosCache.blockStartOffsets[blockId] = range.start;
            blocksFilled++;
        } else {
            logging::warn("Block ID {} out of range for blockStartOffsets (size {}) during cache init.", blockId, globalPosCache.blockStartOffsets.size());
        }
    }

    // Ensure all cache data is written before setting the initialized flag
    std::atomic_thread_fence(std::memory_order_release);
    globalPosCacheInitialized.store(true, std::memory_order_release);
    logging::debug("Global position cache initialized: {} block starts, {} positions mapped to coords. Cache range: [0, {})", 
                 blocksFilled, positionsMapped, globalPosCache.maxPosition);
  } else {
    if (globalPosCacheInitialized.load(std::memory_order_acquire)) {
        logging::debug("Global position cache already initialized.");
    } else if (numCoords == 0) {
        logging::warn("Global position cache not initialized: numCoords is 0.");
    } else if (blockRanges.empty()) {
        logging::warn("Global position cache not initialized: blockRanges is empty.");
    }
  }
}

bool StateManager::isGlobalCacheInitialized() const {
  return globalPosCacheInitialized.load(std::memory_order_acquire);
}

// Helper to initialize block boundaries
void StateManager::initializeBlockBoundaries(const std::string& nodeId) {
  auto& nodeState = getNodeState(nodeId);
  
  // Log the number of active blocks at the start
  logging::debug("Initializing block boundaries for node {}, starting with {} active blocks", 
               nodeId, nodeState.activeBlocks.size());
  
  // Clear current boundaries
  nodeState.firstActiveBlockPos.reset();
  nodeState.lastActiveBlockPos.reset();
  
  // Calculate and set new boundaries if there are active blocks
  if (!nodeState.activeBlocks.empty()) {
    // Log first few block IDs for debugging
    std::stringstream ss;
    ss << "Active blocks for node " << nodeId << ": ";
    int count = 0;
    for (int32_t blockId : nodeState.activeBlocks) {
      ss << blockId << " ";
      if (++count >= 5) { // Limit logging to first 5 blocks
        ss << "...";
        break;
      }
    }
    logging::debug("{}", ss.str());

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
        
        // Log successful range retrieval
        logging::debug("Retrieved range [{}, {}) for block {} in node {}", 
                     blockRange.start, blockRange.end, blockId, nodeId);
      } catch (const std::exception& e) {
        // Log but continue with other blocks
        logging::warn("Could not get range for block {} in node {}: {}", 
                    blockId, nodeId, e.what());
      }
    }
    
    // Log the boundaries
    if (nodeState.firstActiveBlockPos.has_value() && nodeState.lastActiveBlockPos.has_value()) {
      logging::debug("Node {} boundaries: first active block position: {}, last active block position: {}", 
                   nodeId, nodeState.firstActiveBlockPos.value(), nodeState.lastActiveBlockPos.value());
    } else {
      logging::warn("Failed to determine boundaries for node {}", nodeId);
    }
  } else {
    logging::debug("Node {} has no active blocks, skipping boundary initialization", nodeId);
  }
}

// Helper to check if a character is a non-gap character
bool StateManager::isNonGapChar(char c) const {
    // Simple counter to track character distribution
    static std::atomic<size_t> totalCalls{0};
    static std::atomic<size_t> dashCount{0};
    static std::atomic<size_t> nCount{0};
    static std::atomic<size_t> otherCount{0};
    
    // Only count/log at intervals to reduce overhead
    size_t count = totalCalls.fetch_add(1, std::memory_order_relaxed);
    bool isReportInterval = (count > 0 && count % 10000000 == 0);
    
    // Check if the character is a gap or not - simplified to ONLY check for '-'
    bool isGap = (c == '-' || c == 'x');
    
    // Update stats
    if (isGap) {
        dashCount.fetch_add(1, std::memory_order_relaxed);
    } else if (c == 'N' || c == 'n') {
        nCount.fetch_add(1, std::memory_order_relaxed);
    } else {
        otherCount.fetch_add(1, std::memory_order_relaxed);
    }
    
    // Log stats at reporting intervals
    if (isReportInterval) {
        size_t total = totalCalls.load(std::memory_order_relaxed);
        size_t dashes = dashCount.load(std::memory_order_relaxed);
        size_t ns = nCount.load(std::memory_order_relaxed);
        size_t others = otherCount.load(std::memory_order_relaxed);
        
        double dashPct = 100.0 * dashes / total;
        double nPct = 100.0 * ns / total;
        double otherPct = 100.0 * others / total;
        
        fprintf(stderr, "\nCHARACTER DISTRIBUTION after %zu checks: '-': %.1f%%, 'N/n': %.1f%%, other: %.1f%%\n",
               total, dashPct, nPct, otherPct);
    }
    
    return !isGap;
}

/**
 * @brief Set the sequence for a specific block
 * @param blockId The ID of the block
 * @param sequence The sequence to set for the block
 */
void StateManager::setBlockSequence(int32_t blockId, const std::string& sequence) {
    // Log this operation for debugging
    static thread_local int seq_set_log_count = 0;
    if (seq_set_log_count < 5) {
        logging::debug("Setting block {} sequence with length {}", blockId, sequence.length());
        if (!sequence.empty()) {
            logging::debug("Sample: '{}'", sequence.substr(0, std::min(size_t(10), sequence.length())));
        }
        seq_set_log_count++;
    }
    
    blockSequences[blockId] = sequence;
}

// Method to get access to the block sequences (root/reference character data)
const absl::flat_hash_map<int32_t, std::string>& StateManager::getBlockSequences() const {
    return blockSequences;
}

// Method to initialize all coordinate mappings for a block
void StateManager::initializeBlockMappings(int32_t blockId) {
  static std::mutex initMutex;
  static std::unordered_set<int32_t> initializedBlocks;
  static std::atomic<int> blocks_logged_count = 0; // Atomic counter for logging limit
  const int MAX_BLOCKS_TO_LOG = 5;
  
  // Fast check without lock first
  if (initializedBlocks.count(blockId) > 0) {
    return; // Already initialized
  }
  
  // Acquire lock for initialization check
  std::lock_guard<std::mutex> lock(initMutex);
  
  // Check again with lock held
  if (initializedBlocks.find(blockId) != initializedBlocks.end()) {
    return; // Another thread initialized it while we were waiting
  }
  
  // --- REDUNDANT MAPPING CALCULATION REMOVED ---
  // The correct mappings are registered directly in StateManager::initialize
  // This function now only serves to mark the block as initialized.
  
  // Get block sequence and range (still needed for logging if kept)
  auto seqIt = blockSequences.find(blockId);
  auto rangeIt = blockRanges.find(blockId);
  
  if (seqIt == blockSequences.end() || rangeIt == blockRanges.end()) {
    // Log warning but don't return early, still mark as initialized if possible?
    // Or maybe return, as something is wrong if sequence/range are missing here.
    logging::warn("State inconsistency: sequence or range not found for block {} during final initialization marking.", blockId);
    // Let's still mark it initialized to prevent repeated attempts, but log clearly.
    // return; // Decide if returning here is better. Let's not for now.
  }
  
  // Limit logging (Optional: Can remove this logging block entirely)
  bool should_log = (blocks_logged_count.load() < MAX_BLOCKS_TO_LOG);
  if (should_log) {
      logging::debug("Marking block {} mappings as initialized (mappings generated during StateManager::initialize).", blockId);
  }
  // --- END REDUNDANT MAPPING CALCULATION ---
  
  // Mark as initialized
  initializedBlocks.insert(blockId);
  
  // Only log completion and increment counter if we logged the start
  if (should_log) {
      // logging::debug("Marked block {} mappings as initialized.", blockId); // Log adjusted
      blocks_logged_count++; // Increment only after successfully logging start+end for a block
  }
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

// Reset cache for a specific node
void StateManager::resetNodeCache(std::string_view nodeId) {
    // Reset thread-local character cache when node state changes significantly
    HierarchicalStore<char>::clearThreadLocalCache();

    // Invalidate memoized results for this node - REMOVED references to old cache members
    std::string nodeIdStr(nodeId);
    logging::debug("Reset cache for node {}", nodeId);
}

// Set the s-mer size for distance calculations
void StateManager::setSmerSize(int s) {
  if (s > 0) {
    smerSize = s;
  }
}
int StateManager::getSmerSize() const {
  return smerSize;
}

// Get recomputation ranges for a node
const std::vector<coordinates::CoordRange>& StateManager::getRecompRanges(std::string_view nodeId) const {
  auto it = nodeStates.find(std::string(nodeId));
  if (it != nodeStates.end()) {
    // logging::info("STATE_DEBUG: getRecompRanges found {} ranges for node {}", it->second.recompRanges.size(), nodeId); // Modified to direct call
    logging::debug("STATE_DEBUG: getRecompRanges found {} ranges for node {}", it->second.recompRanges.size(), nodeId);
    return it->second.recompRanges;
  }
  
  // Return empty vector if node not found
  static const std::vector<coordinates::CoordRange> emptyRanges;
  // logging::info("STATE_DEBUG: getRecompRanges found 0 ranges (node {} not found or no ranges)", nodeId); // Modified to direct call
  logging::debug("STATE_DEBUG: getRecompRanges found 0 ranges (node {} not found or no ranges)", nodeId);
  return emptyRanges;
}

// Helper function to find the next active block's starting/ending position,
// skipping potential inactive blocks between them.
std::optional<std::pair<int64_t, int32_t>> findAdjacentActiveBlockPosition(
    const std::vector<std::pair<int32_t, coordinates::CoordRange>>& activeBlocks, // Node's currently active blocks, sorted by start coord
    int32_t sourceBlockId, // The block we are trying to move from (might be inactive for the current node)
    const coordinates::CoordRange& sourceBlockGlobalRange, // Global range of the sourceBlockId
    bool searchForward) {

    if (activeBlocks.empty()) {
        logging::debug("findAdjacentActiveBlockPosition: No active blocks for node to jump to from block {}.", sourceBlockId); // Corrected: Added sourceBlockId for {}
        return std::nullopt; // No active blocks to jump to
    }

    if (searchForward) {
        // We are at the end of sourceBlockGlobalRange (e.g., sourceBlockGlobalRange.end - 1).
        // Find the first active block that starts at or after sourceBlockGlobalRange.end.
        int64_t bestNextBlockStart = -1;
        int32_t bestNextBlockId = -1;

        for (const auto& activeBlockPair : activeBlocks) {
            if (activeBlockPair.second.start >= sourceBlockGlobalRange.end) {
                // This active block starts after or at the end of the source block.
                if (bestNextBlockId == -1 || activeBlockPair.second.start < bestNextBlockStart) {
                    bestNextBlockStart = activeBlockPair.second.start;
                    bestNextBlockId = activeBlockPair.first;
                }
            }
        }
        if (bestNextBlockId != -1) {
            logging::debug("findAdjacentActiveBlockPosition: From block {} (ends {}), jumping FORWARD to active block {} (starts {}).", 
                         sourceBlockId, sourceBlockGlobalRange.end, bestNextBlockId, bestNextBlockStart);
            return std::make_pair(bestNextBlockStart, bestNextBlockId);
        }
    } else { // Search backward
        // We are at the start of sourceBlockGlobalRange.
        // Find the first active block that ends at or before sourceBlockGlobalRange.start.
        // Iterate in reverse to find the closest one that ends before current starts.
        int64_t bestPrevBlockEnd = -1; 
        int32_t bestPrevBlockId = -1;

        for (auto it = activeBlocks.rbegin(); it != activeBlocks.rend(); ++it) {
            const auto& activeBlockPair = *it;
            if (activeBlockPair.second.end <= sourceBlockGlobalRange.start) {
                // This active block ends before or at the start of the source block.
                if (bestPrevBlockId == -1 || activeBlockPair.second.end > bestPrevBlockEnd) { 
                    bestPrevBlockEnd = activeBlockPair.second.end;
                    bestPrevBlockId = activeBlockPair.first;
                }
            }
        }
        if (bestPrevBlockId != -1) {
            logging::debug("findAdjacentActiveBlockPosition: From block {} (starts {}), jumping BACKWARD to active block {} (ends {}). Returning {}.", // Corrected: Added one more {} for return value 
                         sourceBlockId, sourceBlockGlobalRange.start, bestPrevBlockId, bestPrevBlockEnd, bestPrevBlockEnd -1);
            return std::make_pair(bestPrevBlockEnd - 1, bestPrevBlockId); 
        }
    }
    logging::debug("findAdjacentActiveBlockPosition: No adjacent active block found for block {} in {} direction.", sourceBlockId, (searchForward ? "FORWARD" : "BACKWARD")); // Corrected: Cast bool to string for logging
    return std::nullopt;
}

// Expand recomputation ranges to ensure complete k-mers
std::vector<coordinates::CoordRange> StateManager::expandRecompRanges(
    std::string_view nodeId_sv, // Corrected to nodeId_sv as in actual code
    const std::vector<coordinates::CoordRange>& ranges, 
    int k_param) { // k_param is kmerSize
  
  std::string nodeId(nodeId_sv); // Convert once

  if (ranges.empty() || k_param <= 0) { // Added check for k_param <= 0
    return {};
  }
  
  // Early check: ensure global cache is initialized before proceeding
  if (!globalPosCacheInitialized.load(std::memory_order_acquire)) {
    logging::err("expandRecompRanges: Global position cache not initialized. Cannot expand ranges for node '{}'. Returning original ranges.", nodeId);
    return ranges;
  }
  
  std::vector<coordinates::CoordRange> expandedRanges;
  expandedRanges.reserve(ranges.size());
  
  std::vector<std::pair<int32_t, coordinates::CoordRange>> activeBlocksForNode;
  std::shared_ptr<HierarchicalGapMap> nodeGapMap;
  try {
      activeBlocksForNode = getActiveBlockRanges(nodeId); // State AFTER mutation application by caller
      nodeGapMap = getNodeGapMap(nodeId); 
  } catch (const std::exception& e) {
      logging::err("expandRecompRanges: Failed to get active blocks or gap map for node '{}': {}. Returning original ranges.", nodeId, e.what());
      return ranges; 
  }

  // Define overall bounds for expansion for this node
  int64_t minNodeExpansionBound = 0;
  int64_t maxNodeExpansionBound = static_cast<int64_t>(numCoords); // Use full coordinate range

  if (!activeBlocksForNode.empty()) {
      minNodeExpansionBound = activeBlocksForNode.front().second.start;
      maxNodeExpansionBound = activeBlocksForNode.back().second.end;
      logging::debug("expandRecompRanges: Node '{}' active region for expansion bounds: [{}, {})", nodeId, minNodeExpansionBound, maxNodeExpansionBound);
  } else if (!this->blockRanges.empty()) { // Node has no active blocks, but global blocks exist
      // If a node has no active blocks, expansion should still be tried if ranges are provided,
      // typically for a block that was just deactivated. Use pangenome limits.
      logging::warn("expandRecompRanges: Node '{}' has no active blocks. Expansion bounds set to pangenome range [{}, {}).", nodeId, minNodeExpansionBound, maxNodeExpansionBound);
  } else {
       logging::warn("expandRecompRanges: Node '{}' has no active blocks and no global blocks defined. Cannot expand. Returning original ranges.", nodeId);
    return ranges; 
  }

  auto ensureExactlyKNonGapChars = 
      [&](int64_t boundaryPos, bool isForward, 
          std::shared_ptr<HierarchicalGapMap> gapMap) -> int64_t { 
    int nonGapCount = 0;
    int64_t currentPos = boundaryPos;
    int64_t lastValidPos = boundaryPos; 
    int iterations = 0;
    int traversalDir = isForward ? 1 : -1;

    int32_t cachedIterationBlockId = -1;
    coordinates::CoordRange cachedIterationBlockRange;
    bool cachedIterationBlockInverted = false;

    while (nonGapCount < k_param) {
        iterations++;
        if (iterations > k_param * 2000 && k_param > 0) { // Increased safety break slightly
            logging::warn("expandRecompRanges: Exceeded max iterations ({}) for node '{}' expanding from boundary {}. Accumulated {} non-gap chars.", iterations, nodeId, boundaryPos, nonGapCount);
             break; 
        }

        // Overall boundary checks for the node/pangenome
        if (isForward && currentPos >= maxNodeExpansionBound) {
            logging::debug("expandRecompRanges: currentPos {} hit maxNodeExpansionBound {}. Stopping.", currentPos, maxNodeExpansionBound);
            break;
        }
        if (!isForward && currentPos < minNodeExpansionBound) {
            logging::debug("expandRecompRanges: currentPos {} hit minNodeExpansionBound {}. Stopping.", currentPos, minNodeExpansionBound);
            break;
        }
        
        if (gapMap) {
            try {
                if (gapMap->isGap(currentPos)) {
                    int64_t nextNonGapPos = gapMap->skipGap(currentPos, isForward, maxNodeExpansionBound); // Use maxNodeExpansionBound as limit for skip
                    if (nextNonGapPos != currentPos && 
                        ((isForward && nextNonGapPos > currentPos) || (!isForward && nextNonGapPos < currentPos))) {
                        if (isForward && nextNonGapPos >= maxNodeExpansionBound) { currentPos = maxNodeExpansionBound; break; }
                        if (!isForward && nextNonGapPos < minNodeExpansionBound) { currentPos = minNodeExpansionBound -1 ; break; }
                    currentPos = nextNonGapPos; 
                    continue; 
                } 
            }
        } catch (const std::exception& e) {
                logging::warn("expandRecompRanges: Error checking/skipping gap at pos {} for node '{}': {}. Proceeding to char check.", currentPos, nodeId, e.what());
            }
        }

        int32_t currentSegmentBlockId = -1;
        coordinates::CoordRange currentSegmentBlockRange;
        bool currentSegmentBlockInverted = false;

        if (cachedIterationBlockId != -1 && currentPos >= cachedIterationBlockRange.start && currentPos < cachedIterationBlockRange.end) {
            currentSegmentBlockId = cachedIterationBlockId;
            currentSegmentBlockRange = cachedIterationBlockRange;
            currentSegmentBlockInverted = cachedIterationBlockInverted;
        } else {
            cachedIterationBlockId = -1; 
            for (const auto& [id, range_iter] : this->blockRanges) { // Use this->blockRanges
                if (currentPos >= range_iter.start && currentPos < range_iter.end) {
                    currentSegmentBlockId = id;
                    currentSegmentBlockRange = range_iter;
                    try { currentSegmentBlockInverted = isBlockInverted(nodeId, id); } 
                    catch (const std::exception& e) { 
                        logging::warn("expandRecompRanges: Error getting inversion status for block {} node '{}': {}. Assuming not inverted.", id, nodeId, e.what());
                        currentSegmentBlockInverted = false; 
                    }
                    cachedIterationBlockId = id; 
                    cachedIterationBlockRange = range_iter; 
                    cachedIterationBlockInverted = currentSegmentBlockInverted;
                    break;
                }
            }
        }
        
        if (currentSegmentBlockId == -1) {
            logging::debug("expandRecompRanges: Pos {} for node '{}' is outside any known global block definition. Stopping expansion {} from {}.", currentPos, nodeId, (isForward ? "forward" : "backward"), boundaryPos);
            break; 
        }

        char c = '-'; 
            try {
                // Bounds check before calling fastMapGlobalToLocal
                if (currentPos < 0 || currentPos >= static_cast<int64_t>(numCoords)) {
                    logging::warn("expandRecompRanges: Position {} outside valid coordinate range [0, {}) for node '{}'. Using fallback '-'.", currentPos, numCoords, nodeId);
                    c = '-';
                } else {
                    auto coord3D = fastMapGlobalToLocal(currentPos);
                    if (coord3D) {
                        int32_t mappedBlockId, structuralNucPos, structuralGapPos;
                        std::tie(mappedBlockId, structuralNucPos, structuralGapPos) = *coord3D;
                        if (mappedBlockId == currentSegmentBlockId) {
                            c = getCharAtPosition(nodeId, mappedBlockId, structuralNucPos, structuralGapPos);
                        } else {
                            logging::warn("expandRecompRanges: Mismatch! currentPos {} in currentSegmentBlockId {} but fastMapGlobalToLocal returned mappedBlockId {}. Using fallback '-'.", currentPos, currentSegmentBlockId, mappedBlockId);
                        }
                    } else {
                        logging::warn("expandRecompRanges: Failed to map globalPos {} to local for node '{}'. Using fallback '-'.", currentPos, nodeId);
                    }
                }
            } catch (const std::exception& e) { 
                 logging::warn("expandRecompRanges: Exception getting char at pos {} (node '{}'): {}. Using fallback '-'.", currentPos, nodeId, e.what());
            }

        if (c != '-' && c != 'x') { 
            nonGapCount++;
            lastValidPos = currentPos;
        }

        if (nonGapCount >= k_param) break;

        bool effectivelyForwardDirection = (isForward && !currentSegmentBlockInverted) || (!isForward && currentSegmentBlockInverted);
        bool atCurrentBlockBoundary;
        if (effectivelyForwardDirection) {
            atCurrentBlockBoundary = (currentPos == currentSegmentBlockRange.end - 1);
        } else {
            atCurrentBlockBoundary = (currentPos == currentSegmentBlockRange.start);
        }

        if (atCurrentBlockBoundary) { 
            // Pass currentSegmentBlockId AND currentSegmentBlockRange to the modified findAdjacentActiveBlockPosition
            auto nextBlockJumpInfo = findAdjacentActiveBlockPosition(activeBlocksForNode, currentSegmentBlockId, currentSegmentBlockRange, effectivelyForwardDirection);
            if (nextBlockJumpInfo) {
                currentPos = nextBlockJumpInfo->first; 
                cachedIterationBlockId = -1; // Force re-evaluation of block for currentPos
            } else {
                logging::debug("expandRecompRanges: Hit boundary of block {} and no adjacent active block to jump to for node '{}'. Stopping.", currentSegmentBlockId, nodeId);
                break; 
            }
        } else {
            currentPos += traversalDir; 
        }
    } // End while

    if (isForward) {
        return (nonGapCount >= k_param || boundaryPos == lastValidPos) ? lastValidPos + 1 : boundaryPos; 
    } else { 
        return (nonGapCount >= k_param || boundaryPos == lastValidPos) ? lastValidPos : boundaryPos;
    }
  }; 
  
  auto expandSingleRange = 
      [&](const coordinates::CoordRange& range) -> coordinates::CoordRange {
    // (Lambda definition body remains the same, calls ensureExactlyKNonGapChars)
    int64_t clampedStart = std::clamp(range.start, minNodeExpansionBound, maxNodeExpansionBound);
    int64_t clampedEnd = std::clamp(range.end, minNodeExpansionBound, maxNodeExpansionBound);
    
    // Only reject ranges that are truly invalid (negative size after clamping)
    // Empty ranges (start == end) are valid and should be processed
    if (clampedEnd < clampedStart) {
         logging::warn("expandRecompRanges: Initial range [{}, {}) invalid or outside active bounds [{}, {}) node '{}'. Skipping.",
                   range.start, range.end, minNodeExpansionBound, maxNodeExpansionBound, nodeId);
         return range; 
    }
    int64_t expandedStart = ensureExactlyKNonGapChars(clampedStart, false, nodeGapMap); 
    int64_t expandedEnd = clampedEnd;
    coordinates::CoordRange expandedRange{expandedStart, expandedEnd};
    if (expandedRange.end < expandedRange.start) {
         logging::err("expandRecompRanges: Invalid expanded range [{}, {}) node '{}'. Original=[{}, {}), Clamped=[{}, {}).", expandedRange.start, expandedRange.end, nodeId, range.start, range.end, clampedStart, clampedEnd);
         return {clampedStart, clampedEnd};
    }
    return expandedRange;
  }; // --- End lambda expandSingleRange ---
  
  // Process each range using expandSingleRange
  for (const auto& range : ranges) {
    try {
      coordinates::CoordRange expanded = expandSingleRange(range);
      if (expanded.start < expanded.end) {
         expandedRanges.push_back(expanded);
      } else {
         logging::warn("expandRecompRanges: Skipping invalid expanded range [{}, {}) for original [{}, {}) node '{}'", expanded.start, expanded.end, range.start, range.end, nodeId);
      }
    } catch (const std::exception& e) {
      logging::err("expandRecompRanges: Exception expanding range [{}, {}) node '{}': {}. Keeping original.", range.start, range.end, nodeId, e.what());
       if (range.start < range.end) { 
          expandedRanges.push_back(range);
       }
    }
  }
  
  // Merge overlapping ranges
  size_t preMergeCount = expandedRanges.size();
  if (expandedRanges.size() > 1) {
    std::sort(expandedRanges.begin(), expandedRanges.end(), 
             [](const auto& a, const auto& b) { return a.start < b.start; });
    std::vector<coordinates::CoordRange> mergedRanges;
    if (!expandedRanges.empty()) { 
      mergedRanges.push_back(expandedRanges[0]);
      for (size_t i = 1; i < expandedRanges.size(); ++i) {
        auto& currentRange = expandedRanges[i];
        auto& lastRange = mergedRanges.back();
        if (currentRange.start <= lastRange.end) {
          lastRange.end = std::max(lastRange.end, currentRange.end);
        } else {
          mergedRanges.push_back(currentRange);
        }
      }
      expandedRanges = std::move(mergedRanges);
    }
  }
  
  logging::info("expandRecompRanges: Node '{}' Original ranges: {}, Expanded: {}, Merged: {}",
              nodeId, ranges.size(), preMergeCount, expandedRanges.size());
  
  return expandedRanges;
} 

// Add a seed to a block for tracking
void StateManager::addSeedToBlock(int32_t blockId, int64_t seedPos) {
  if (blockId >= 0 && blockId < static_cast<int32_t>(numBlocks)) {
    blockToSeeds[blockId].insert(seedPos);
  } else {
    logging::warn("addSeedToBlock: Invalid blockId {} (must be in range [0, {})), seed position {}", 
                 blockId, numBlocks, seedPos);
  }
}

// Remove a seed from a block's tracking
void StateManager::removeSeedFromBlock(int32_t blockId, int64_t seedPos) {
  if (blockId >= 0 && blockId < static_cast<int32_t>(numBlocks)) {
    auto it = blockToSeeds.find(blockId);
    if (it != blockToSeeds.end()) {
      it->second.erase(seedPos);
    }
  }
}

// Map global position to block coordinates
coordinates::BlockCoordinate StateManager::mapGlobalToBlockCoords(
    std::basic_string_view<char, std::char_traits<char>> nodeId,
    int64_t globalPos) const {
  
  // Validate input position
  if (globalPos < 0 || globalPos >= static_cast<int64_t>(numCoords)) {
    return {-1, -1, -1}; // Invalid position
  }
  
  // First, determine which block contains this position
  int32_t blockId = -1;
  
  // Use global position cache for fast lookup if initialized
  if (globalPosCacheInitialized.load(std::memory_order_acquire) && globalPos < globalPosCache.maxPosition) {
    blockId = globalPosCache.getBlockId(globalPos); // NEW - Use helper method
  } else {
    // Fallback to linear search through block ranges
    for (const auto& [id, range] : blockRanges) {
      if (globalPos >= range.start && globalPos < range.end) {
        blockId = id;
        break;
      }
    }
  }
  
  if (blockId < 0) {
    // Position not in any block
    return {-1, -1, -1};
  }
  
  // Check if block is active for this node
  if (!isBlockOn(nodeId, blockId)) {
    // Block exists but is not active for this node
    return {blockId, -1, -1}; // Return blockId but invalid nuc/gap pos
  }
  
  // Get if block is inverted
  bool inverted = isBlockInverted(nodeId, blockId);
  
  // Get block range
  const auto& range = blockRanges.at(blockId);
  
  // Calculate offset within block
  int64_t blockOffset = globalPos - range.start;
  
  // Map to 3D coordinates using fastMapGlobalToLocal
  auto coordsOpt = fastMapGlobalToLocal(globalPos);
  if (!coordsOpt) {
    return {blockId, -1, -1}; // Mapping failed
  }
  
  auto [mappedBlockId, nucPos, gapPos] = *coordsOpt;
  
  // Double-check block ID matches what we expect
  if (mappedBlockId != blockId) {
    // This is a serious inconsistency that shouldn't happen
    logging::warn("Block ID mismatch during coordinate mapping: expected {}, got {}",
                 blockId, mappedBlockId);
    return {blockId, -1, -1};
  }
  
  // Return the final coordinates, handling inversion if needed
  if (inverted) {
    // For inverted blocks, invert the nucPos within the block
    // The exact inversion logic depends on block-specific details
    // (implemented within the fastMapGlobalToLocal)
    return {blockId, nucPos, gapPos};
  } else {
    return {blockId, nucPos, gapPos};
  }
}

// Extract a sequence from a range of positions
std::tuple<std::string, std::vector<int64_t>, std::vector<bool>, std::vector<int64_t>> StateManager::extractSequence(
    std::string_view nodeId_sv, 
    const coordinates::CoordRange& range, 
    bool skipGaps) const {
    
    std::string nodeId(nodeId_sv);
    
    // Add specific debug logging for problem blocks
    bool debug_problem_blocks = (nodeId == "node_2");
    std::ofstream debug_extract;
    
    
    CharacterBuffer& charBuffer = charBuffers.local(); 
    charBuffer.clear();
    
    std::vector<std::pair<int32_t, coordinates::CoordRange>> active_blocks_in_node;
    std::shared_ptr<HierarchicalGapMap> nodeGapMap;
    const NodeState* nodeState = nullptr; // Declare in broader scope for gap run checking

    try {
        std::shared_lock<std::shared_mutex> lock(nodeMutex); // Lock for reading nodeStates 
        auto nodeStateIt = nodeStates.find(nodeId);
        if (nodeStateIt == nodeStates.end()) {
            logging::err("extractSequence: Node state not found for '{}'. Returning empty.", nodeId);
            if (debug_problem_blocks && debug_extract.is_open()) {
                debug_extract << "  -> ERROR: Node state not found, returning empty" << std::endl;
                debug_extract.close();
            }
            return {"", {}, {}, {}};
        }
        
        // Get node state pointer for fast gap run checking
        nodeState = &nodeStateIt->second;
        
        // nodeStatePtr = &nodeStateIt->second; // Assign if needed later, but gapMap is key here
        if (!nodeStateIt->second.gapMap) { // Check directly from iterator
             logging::err("extractSequence: GapMap is null for node '{}'. Returning empty.", nodeId);
             if (debug_problem_blocks && debug_extract.is_open()) {
                debug_extract << "  -> ERROR: GapMap is null, returning empty" << std::endl;
                debug_extract.close();
             }
             return {"", {}, {}, {}};
        }
        nodeGapMap = nodeStateIt->second.gapMap;
        // Active blocks can be fetched outside the lock if getActiveBlockRanges is thread-safe or re-locks
        // However, to ensure consistency of state used (gapMap + activeBlocks from same logical point for nodeId)
        // it might be better to fetch active_blocks_in_node here too if it doesn't cause deadlocks.
        // For now, assuming getActiveBlockRanges handles its own locking or is safe after this read lock.
        lock.unlock();
        
        active_blocks_in_node = getActiveBlockRanges(nodeId); 
    } catch (const std::exception& e) {
        logging::err("extractSequence: Failed to get active blocks or gap map for node '{}': {}. Returning empty.", nodeId, e.what());
        if (debug_problem_blocks && debug_extract.is_open()) {
            debug_extract << "  -> EXCEPTION getting active blocks/gap map: " << e.what() << std::endl;
            debug_extract.close();
        }
        return {"", {}, {}, {}};
    }
    
    if (active_blocks_in_node.empty()) {
      logging::debug("extractSequence: No active blocks found for node '{}' to process for range [{}, {})", nodeId, range.start, range.end);
      return {"", {}, {}, {}};
    }
    
    int totalCharsProcessed_debug = 0;
    int nonGapCharsFound_debug = 0;

    for (const auto& [blockId, block_global_range] : active_blocks_in_node) {
        coordinates::CoordRange processing_range_for_this_block{
            std::max(block_global_range.start, range.start),
            std::min(block_global_range.end, range.end)
        };

        if (processing_range_for_this_block.start >= processing_range_for_this_block.end) {
            continue;
        }
        
        // Special debug for problem blocks 617 and 941
        bool is_problem_block = (debug_problem_blocks && (blockId == 617 || blockId == 941));
        

        bool blockIsInverted = isBlockInverted(nodeId, blockId); // Uses isBlockInverted_unsafe via public wrapper
        

        int32_t max_nuc_pos_for_block = -1; 
        auto flat_indices_map_it = blockRootCharFlatIndices.find(blockId);
        if (flat_indices_map_it != blockRootCharFlatIndices.end()){
            for(const auto& [pos_key, flat_idx] : flat_indices_map_it->second){
                max_nuc_pos_for_block = std::max(max_nuc_pos_for_block, pos_key.nucPos);
            }
        }
        
        

        std::vector<char> temp_chars_for_this_block;
        std::vector<int64_t> temp_positions_for_this_block;
        std::vector<bool> temp_gaps_for_this_block;

        // OPTIMIZATION: Use bulk character extraction for the entire block range
        int64_t blockProcessingStart = processing_range_for_this_block.start;
        int64_t blockProcessingEnd = processing_range_for_this_block.end;
        
        std::vector<char> blockBulkChars;
        std::vector<int64_t> blockBulkPositions;
        std::vector<bool> blockBulkGaps;
        
        try {
          // Extract all characters for this block range at once
          for (int64_t globalPos = blockProcessingStart; globalPos < blockProcessingEnd; ++globalPos) {
            auto coordsOpt = fastMapGlobalToLocal(globalPos);
            if (!coordsOpt) continue;
            
            auto [blockIdMapped, nucPos, gapPos] = *coordsOpt;
            if (blockIdMapped != blockId) continue; // Only process characters from current block
            
            char char_to_add;
            if (nodeState && nodeState->isGapRun(blockId, nucPos)) {
              char_to_add = '-';
            } else {
              char_to_add = getCharAtPosition(nodeId_sv, blockId, nucPos, gapPos);
            }
            
            bool isGap = (char_to_add == '-' || char_to_add == 'x');
            
            if (!(skipGaps && isGap)) {
              blockBulkChars.push_back(char_to_add);
              blockBulkPositions.push_back(globalPos);
              blockBulkGaps.push_back(isGap);
              if (isNonGapChar(char_to_add)) nonGapCharsFound_debug++;
            }
            totalCharsProcessed_debug++;
          }
          
          // Apply block inversion if needed
          if (blockIsInverted) {
            if (is_problem_block && debug_extract.is_open()) {
              debug_extract << "     Reversing " << blockBulkChars.size() << " characters due to block inversion" << std::endl;
            }
            std::reverse(blockBulkChars.begin(), blockBulkChars.end());
            std::reverse(blockBulkPositions.begin(), blockBulkPositions.end());
            std::reverse(blockBulkGaps.begin(), blockBulkGaps.end());
          }
          
          // Add to result buffers
          charBuffer.buffer.insert(charBuffer.buffer.end(), blockBulkChars.begin(), blockBulkChars.end());
          charBuffer.positions.insert(charBuffer.positions.end(), blockBulkPositions.begin(), blockBulkPositions.end());
          charBuffer.gaps.insert(charBuffer.gaps.end(), blockBulkGaps.begin(), blockBulkGaps.end());
          
          if (is_problem_block && debug_extract.is_open()) {
            debug_extract << "     Final chars for block " << blockId << ": '" 
                         << std::string(blockBulkChars.begin(), blockBulkChars.end()) << "'" << std::endl;
            debug_extract << "     Length: " << blockBulkChars.size() << " characters" << std::endl;
          }
          
        } catch (const std::exception& e) {
        } catch (const std::exception& e) {
          // Fallback to original individual character extraction if bulk fails
          logging::warn("Bulk extraction failed for block {}, falling back to individual extraction: {}", blockId, e.what());
          
          std::vector<char> temp_chars_for_this_block;
          std::vector<int64_t> temp_positions_for_this_block;
          std::vector<bool> temp_gaps_for_this_block;

          for (int32_t nuc_p = 0; nuc_p <= max_nuc_pos_for_block; ++nuc_p) {
              size_t gap_list_len = getGapListLength(blockId, nuc_p);
              
              for (size_t gap_idx = 0; gap_idx < gap_list_len; ++gap_idx) {
                  totalCharsProcessed_debug++;
                  int64_t global_pos_gap_char = fastMapLocalToGlobal(blockId, nuc_p, static_cast<int32_t>(gap_idx));
                  if (global_pos_gap_char >= processing_range_for_this_block.start && global_pos_gap_char < processing_range_for_this_block.end) {
                      char char_to_add;
                      if (nodeState && nodeState->isGapRun(blockId, nuc_p)) {
                          char_to_add = '-';
                      } else {
                          char_to_add = getCharAtPosition(nodeId_sv, blockId, nuc_p, static_cast<int32_t>(gap_idx));
                      }
                      
                      if (!(skipGaps && (char_to_add == '-' || char_to_add == 'x'))) {
                          temp_chars_for_this_block.push_back(char_to_add);
                          temp_positions_for_this_block.push_back(global_pos_gap_char);
                          temp_gaps_for_this_block.push_back(char_to_add == '-' || char_to_add == 'x');
                          if (isNonGapChar(char_to_add)) nonGapCharsFound_debug++;
                      }
                  }
              }

              totalCharsProcessed_debug++;
              int64_t global_pos_main_nuc = fastMapLocalToGlobal(blockId, nuc_p, -1);
              if (global_pos_main_nuc >= processing_range_for_this_block.start && global_pos_main_nuc < processing_range_for_this_block.end) {
                  char char_to_add;
                  if (nodeState && nodeState->isGapRun(blockId, nuc_p)) {
                      char_to_add = '-';
                  } else {
                      char_to_add = getCharAtPosition(nodeId_sv, blockId, nuc_p, -1);
                  }
                  
                  if (!(skipGaps && (char_to_add == '-' || char_to_add == 'x'))) {
                      temp_chars_for_this_block.push_back(char_to_add);
                      temp_positions_for_this_block.push_back(global_pos_main_nuc);
                      temp_gaps_for_this_block.push_back(char_to_add == '-' || char_to_add == 'x');
                      if (isNonGapChar(char_to_add)) nonGapCharsFound_debug++;
                  }
              }
          }

          if (blockIsInverted) {
              std::reverse(temp_chars_for_this_block.begin(), temp_chars_for_this_block.end());
              std::reverse(temp_positions_for_this_block.begin(), temp_positions_for_this_block.end());
              std::reverse(temp_gaps_for_this_block.begin(), temp_gaps_for_this_block.end());
          }
          
          charBuffer.buffer.insert(charBuffer.buffer.end(), temp_chars_for_this_block.begin(), temp_chars_for_this_block.end());
          charBuffer.positions.insert(charBuffer.positions.end(), temp_positions_for_this_block.begin(), temp_positions_for_this_block.end());
          charBuffer.gaps.insert(charBuffer.gaps.end(), temp_gaps_for_this_block.begin(), temp_gaps_for_this_block.end());
        }
    }
    
    std::string result_str(charBuffer.buffer.begin(), charBuffer.buffer.end());
    std::vector<int64_t> result_pos = charBuffer.positions;
    std::vector<bool> result_gaps = charBuffer.gaps;
    
    // Calculate end positions for each potential k-mer starting position
    std::vector<int64_t> result_end_positions;
    result_end_positions.reserve(result_pos.size());
    
    // Get k value from StateManager
    int k = getKmerSize();
    
    for (size_t i = 0; i < result_pos.size(); i++) {
        int64_t endPos = -1; // Default end position if no k-mer found
        
        // Only calculate end position for non-gap starting positions
        if (i < result_gaps.size() && !result_gaps[i]) {
            size_t nonGapCount = 0;
            size_t currentIdx = i;
            
            // Count k non-gap bases starting from position i
            while (currentIdx < result_gaps.size() && nonGapCount < k) {
                if (!result_gaps[currentIdx]) {
                    nonGapCount++;
                }
                // If we've found k non-gap bases, the end position is the position after this character
                if (nonGapCount == k && currentIdx < result_pos.size()) {
                    endPos = result_pos[currentIdx] + 1; // Position after the k-th character
                    break;
                }
                currentIdx++;
            }
        }
        
        result_end_positions.push_back(endPos);
    }
    
    if (result_str.empty() && range.start < range.end && !active_blocks_in_node.empty()) {
      // Rate limit this warning to avoid spam
      static std::atomic<int> extract_warning_count{0};
      if (extract_warning_count.fetch_add(1) < 10) {
        logging::warn("extractSequence: Failed to extract any characters for node '{}' from range [{}:{}) using structural iteration, though active blocks overlap.",
                    nodeId, range.start, range.end);
      }
      if (debug_problem_blocks && debug_extract.is_open()) {
        debug_extract << "  -> WARNING: Failed to extract any characters despite active blocks" << std::endl;
      }
    }
    
    if (debug_problem_blocks && debug_extract.is_open()) {
        debug_extract << "  -> Final result: \"" << result_str << "\"" << std::endl;
        debug_extract << "  -> Length: " << result_str.size() << " characters" << std::endl;
        debug_extract.close();
    }
    
    logging::debug("extractSequence (structural): Extracted {} chars ({} non-gap) for node '{}' from range [{}:{}), total processed from blocks: {}",
                 result_str.size(), nonGapCharsFound_debug, nodeId, range.start, range.end, totalCharsProcessed_debug);
    
    return {result_str, result_pos, result_gaps, result_end_positions};
}

// Helper method to get block ID from a global position
int32_t StateManager::getBlockIdFromPosition(int64_t pos) const {
  // Validate position
  if (pos < 0 || pos >= static_cast<int64_t>(numCoords)) {
    throw std::out_of_range("Position " + std::to_string(pos) + 
                           " is out of range (must be between 0 and " + 
                           std::to_string(numCoords - 1) + ")");
  }
  
  // Use the global position cache for efficient lookup if initialized
  if (globalPosCacheInitialized.load(std::memory_order_acquire) && 
      pos < static_cast<int64_t>(globalPosCache.maxPosition)) { // Check against maxPosition
    // int32_t blockId = globalPosCache.posToBlockId[pos]; // OLD - Incorrect access
    int32_t blockId = globalPosCache.getBlockId(pos); // NEW - Use helper method
    if (blockId >= 0) {
      return blockId;
    }
    // Fall through to linear search if cache doesn't have a mapping
  } else if (!globalPosCacheInitialized.load(std::memory_order_acquire)) {
    throw std::runtime_error("getBlockIdFromPosition called before globalPosCache was initialized");
  }
  
  // Fallback to linear search through block ranges
  for (const auto& [blockId, range] : blockRanges) {
    if (pos >= range.start && pos < range.end) {
      return blockId;
    }
  }
  
  // If we get here, the position isn't in any known block
  throw std::runtime_error("Position " + std::to_string(pos) + 
                          " is not contained in any known block");
}

// Initialize node hierarchy tree structure
void StateManager::initializeNodeHierarchy(panmanUtils::Tree *tree, panmanUtils::Node *rootNode) {
  if (!tree || !rootNode) {
    logging::err("Invalid tree or root node in initializeNodeHierarchy");
    return;
  }

  logging::info("Initializing node hierarchy for {} nodes", tree->allNodes.size());
  
  // CRITICAL FIX: Before clearing the hierarchy, preserve the block ranges if they exist
  absl::flat_hash_map<int32_t, coordinates::CoordRange> preservedBlockRanges;
  if (!blockRanges.empty()) {
    logging::info("Preserving {} existing block ranges before reinitializing hierarchy", blockRanges.size());
    // Copy entries individually
    for (const auto& [blockId, range] : blockRanges) {
      preservedBlockRanges[blockId] = range;
    }
  }
  
  // Clear existing hierarchy
  nodeHierarchy.clear();
  
  // Initialize the root node first
  std::string rootId = rootNode->identifier;
  nodeHierarchy[rootId] = NodeHierarchy{};
  nodeHierarchy[rootId].parentId = ""; // Root has no parent
  
  // Process all nodes in the tree
  std::function<void(panmanUtils::Node*, const std::string&)> processNode;
  processNode = [&](panmanUtils::Node* node, const std::string& parentId) {
    if (!node) return;
    
    std::string nodeId = node->identifier;
    
    // Create or get the node hierarchy entry
    auto& hierarchy = nodeHierarchy[nodeId];
    
    // Set parent relationship
    hierarchy.parentId = parentId;
    
    // Record parent-child relationship for both parent and child
    if (!parentId.empty()) {
      nodeHierarchy[parentId].childrenIds.push_back(nodeId);
    }
    
    // Recursively process all children
    for (auto* child : node->children) {
      if (child) {
        processNode(child, nodeId);
      }
    }
  };
  
  // Start processing from the root node
  processNode(rootNode, "");
  
  // CRITICAL FIX: Restore preserved block ranges after hierarchy initialization
  if (!preservedBlockRanges.empty()) {
    logging::info("Restoring {} preserved block ranges after hierarchy initialization", preservedBlockRanges.size());
    
    // Copy entries individually
    for (const auto& [blockId, range] : preservedBlockRanges) {
      blockRanges[blockId] = range;
    }
    
    // Extra verification for critical blocks 617 and 941
    auto it617 = blockRanges.find(617);
    auto it941 = blockRanges.find(941);
    if (it617 != blockRanges.end()) {
      logging::info("Verified block 617 range preserved: [{}, {})", it617->second.start, it617->second.end);
    }
    if (it941 != blockRanges.end()) {
      logging::info("Verified block 941 range preserved: [{}, {})", it941->second.start, it941->second.end);
    }
  }
  
  // Ensure rootGapMap is initialized
  if (!rootGapMap) {
    rootGapMap = std::make_shared<HierarchicalGapMap>();
    logging::info("Created root gap map");
  }
  
  // Create root character store if needed
  static std::shared_ptr<HierarchicalCharacterStore> rootCharacterStore;
  if (!rootCharacterStore) {
    rootCharacterStore = std::make_shared<HierarchicalCharacterStore>();
    logging::info("Created root character store");
  }
  
  // Create root seed store if needed (parallel to character store)
  static std::shared_ptr<HierarchicalSeedStore> rootSeedStore;
  if (!rootSeedStore) {
    rootSeedStore = std::make_shared<HierarchicalSeedStore>();
    logging::info("Created root seed store");
  }
  
  // Set up hierarchical structures for all nodes
  std::vector<std::string> processedNodes;
  processedNodes.reserve(nodeHierarchy.size());
  
  // First pass - set up the root node
  nodeHierarchy[rootId].gapMap = rootGapMap;
  nodeHierarchy[rootId].characterStore = rootCharacterStore;
  nodeHierarchy[rootId].seedStore = rootSeedStore; // Initialize root seed store
  processedNodes.push_back(rootId);
  
  // Second pass - process all other nodes in topological order
  auto sortedNodes = buildTopologicalOrder(tree);
  
  // Skip the root node as we already processed it
  for (const auto& nodeId : sortedNodes) {
    if (nodeId == rootId) continue;
    
    auto& hierarchy = nodeHierarchy[nodeId];
    auto& parentId = hierarchy.parentId;
    
    // Get parent's hierarchical structures
    if (!parentId.empty() && nodeHierarchy.find(parentId) != nodeHierarchy.end()) {
      auto& parentHierarchy = nodeHierarchy[parentId];
      
      // Set up gap map inheritance
      if (parentHierarchy.gapMap) {
        hierarchy.gapMap = std::make_shared<HierarchicalGapMap>(parentHierarchy.gapMap);
      } else {
        hierarchy.gapMap = std::make_shared<HierarchicalGapMap>(rootGapMap);
        logging::warn("Parent {} has no gap map during hierarchy setup for {}", parentId, nodeId);
      }
      
      // Set up character store inheritance
      if (parentHierarchy.characterStore) {
        hierarchy.characterStore = std::make_shared<HierarchicalCharacterStore>(parentHierarchy.characterStore);
      } else {
        hierarchy.characterStore = std::make_shared<HierarchicalCharacterStore>(rootCharacterStore);
        logging::warn("Parent {} has no character store during hierarchy setup for {}", parentId, nodeId);
      }
      
      // Set up seed store inheritance (parallel to character store)
      if (parentHierarchy.seedStore) {
        hierarchy.seedStore = std::make_shared<HierarchicalSeedStore>(parentHierarchy.seedStore);
      } else {
        hierarchy.seedStore = std::make_shared<HierarchicalSeedStore>(rootSeedStore);
        logging::warn("Parent {} has no seed store during hierarchy setup for {}", parentId, nodeId);
      }
    } else {
      // No parent found, use root structures
      hierarchy.gapMap = std::make_shared<HierarchicalGapMap>(rootGapMap);
      hierarchy.characterStore = std::make_shared<HierarchicalCharacterStore>(rootCharacterStore);
      hierarchy.seedStore = std::make_shared<HierarchicalSeedStore>(rootSeedStore);
      logging::warn("Node {} has no valid parent in hierarchy, using root structures", nodeId);
    }
    
    // Add to processed nodes
    processedNodes.push_back(nodeId);
  }
  
  // Verify that all nodes in the hierarchy have proper structures
  size_t missingGapMaps = 0;
  size_t missingCharacterStores = 0;
  size_t missingSeedStores = 0;
  
  for (const auto& [nodeId, hierarchy] : nodeHierarchy) {
    if (!hierarchy.gapMap) {
      missingGapMaps++;
      logging::warn("Node {} missing gap map after hierarchy initialization", nodeId);
    }
    
    if (!hierarchy.characterStore) {
      missingCharacterStores++;
      logging::warn("Node {} missing character store after hierarchy initialization", nodeId);
    }
    
    if (!hierarchy.seedStore) {
      missingSeedStores++;
      logging::warn("Node {} missing seed store after hierarchy initialization", nodeId);
    }
  }
  
  if (missingGapMaps > 0 || missingCharacterStores > 0 || missingSeedStores > 0) {
    logging::err("Hierarchy initialization incomplete: {} nodes missing gap maps, {} missing character stores, {} missing seed stores",
               missingGapMaps, missingCharacterStores, missingSeedStores);
  }
  
  // Log summary
  size_t nodeCount = nodeHierarchy.size();
  logging::info("Initialized hierarchy for {} nodes, {} processed successfully", 
               nodeCount, processedNodes.size());
}

// Set logging verbosity level for state operations
void StateManager::setVerboseLogging(bool verbose) {
  // Access the static variables through static methods to ensure
  // they're properly initialized and shared across the code
  static bool& init_logging = getInitLoggingRef();
  static bool& prop_logging = getPropLoggingRef();
  
  init_logging = verbose;
  prop_logging = verbose;
  
  logging::info("State manager verbose logging set to {}", verbose);
}

// Helper methods to access the static variables safely
bool& StateManager::getInitLoggingRef() {
  static bool verbose_init_logging = false;
  return verbose_init_logging;
}

bool& StateManager::getPropLoggingRef() {
  static bool verbose_prop_logging = false;
  return verbose_prop_logging;
}

// Fast lookup from global position to (blockId, nucPos, gapPos)
std::optional<std::tuple<int32_t, int32_t, int32_t>>
StateManager::fastMapGlobalToLocal(int64_t globalPos) const {
  // Check if the cache is initialized and the position is valid
  if (!globalPosCacheInitialized.load(std::memory_order_acquire) || !globalPosCache.isValid(globalPos)) {
    // Keep minimal logging for critical issues like uninitialized cache or out-of-bounds access
    if (!globalPosCacheInitialized.load(std::memory_order_acquire)) {
        logging::warn("fastMapGlobalToLocal: Cache not initialized for coordinate mapping (pos: {}, numCoords: {}, blockRanges.size: {})", 
                      globalPos, numCoords, blockRanges.size());
    } else { // Must be invalid position
        logging::warn("fastMapGlobalToLocal: Position {} outside valid cache range [0, {}) (cache maxPos: {})", 
                      globalPos, globalPosCache.maxPosition, globalPosCache.maxPosition);
    }
    return std::nullopt;
  }

  // Direct lookup using the precomputed cache
  const auto& cachedCoords = globalPosCache.posToCoords[globalPos];

  // Log if the mapping is missing in the cache (e.g., position between blocks)
  // Limit this logging to reduce spam
  if (!cachedCoords.has_value()) {
    static thread_local int map_miss_log_count = 0;
    const int MAX_MAP_MISS_LOGS = 5;
    if (map_miss_log_count < MAX_MAP_MISS_LOGS) {
        logging::debug("fastMapGlobalToLocal: Position {} not mapped in cache (likely between blocks).", globalPos);
        map_miss_log_count++;
    }
  }
  
  // Return the cached value (which might be std::nullopt)
  return cachedCoords; 
}

// Fast lookup from (blockId, nucPos, gapPos) to global position
int64_t StateManager::fastMapLocalToGlobal(int32_t blockId, int32_t nucPos, int32_t gapPos) const {
    PositionKey key = PositionKey::create(blockId, nucPos, gapPos);
    auto it = blockCoordToGlobalPos.find(key);
    if (it != blockCoordToGlobalPos.end()) {
        return it->second;
    }
    return -1; // Return -1 to indicate the mapping was not found.
}

// Extract a k-mer starting from a specific position (scanning until finding k non-gap chars)
std::pair<std::string, std::vector<int64_t>> StateManager::extractKmer(
    std::string_view nodeId, int64_t pos, int k, bool reverse) const {
  
  // Declare tracking variables for diagnostics with correct global scope
  // These are in global scope (no namespace), not in state::
  extern ::std::atomic<size_t> totalKmerAttempts;
  extern ::std::atomic<size_t> successfulKmerExtractions;
  extern ::std::atomic<size_t> kmerFailedDueToLength;
  extern ::std::atomic<size_t> kmerFailedDueToOnlyGaps;
  
  // Increment attempt counter if available
  try {
    totalKmerAttempts.fetch_add(1, std::memory_order_relaxed);
  } catch (...) {
    // Ignore if counter isn't available
  }
  
  // Main kmer extraction implementation
  std::string rawKmer; // Store raw sequence with gaps
  std::string ungappedKmer; // Store sequence with gaps removed
  std::vector<int64_t> positions; // non-gapped positions of each k-mer character
  int nonGapCount = 0;

  rawKmer.reserve(k * 2); // Reserve more to account for gaps
  ungappedKmer.reserve(k); // Reserve exact k size for ungapped result
  positions.reserve(k * 2);

  // OPTIMIZATION: Use bulk character extraction to avoid individual getCharAtPosition calls
  // Estimate characters needed (accounting for gaps)
  int64_t estimatedLength = k * 2; // Heuristic: assume up to 50% gaps
  int64_t startPos = reverse ? std::max(pos - estimatedLength + 1, static_cast<int64_t>(0)) : pos;
  int64_t endPos = reverse ? pos + 1 : std::min(pos + estimatedLength, static_cast<int64_t>(numCoords));
  
  // Adjust bounds if necessary
  if (startPos < 0) startPos = 0;
  if (endPos > static_cast<int64_t>(numCoords)) endPos = numCoords;
  
  // Extract bulk characters using coordinate mapping
  std::vector<char> bulkChars;
  std::vector<int64_t> bulkPositions;
  
  try {
    std::shared_ptr<HierarchicalGapMap> nodeGapMap;
    try {
        nodeGapMap = getNodeGapMap(nodeId); 
    } catch (const std::exception& e) {
        logging::warn("StateManager::extractKmer: Could not get gap map for node {}: {}", nodeId, e.what());
    }
    
    // Extract characters for the range
    for (int64_t currentPos = startPos; currentPos < endPos; ++currentPos) {
      auto coordsOpt = fastMapGlobalToLocal(currentPos);
      if (!coordsOpt) continue;
      
      auto [blockId, nucPos, gapPos] = *coordsOpt;
      char c = getCharAtPosition(nodeId, blockId, nucPos, gapPos);
      
      bulkChars.push_back(c);
      bulkPositions.push_back(currentPos);
    }
  } catch (const std::exception& e) {
    // Fallback to individual character extraction if bulk fails
    int64_t currentPos = pos;
    int direction = reverse ? -1 : 1;
    
    while (nonGapCount < k) { 
      if (currentPos < 0 || currentPos >= static_cast<int64_t>(numCoords)) {
        try {
          kmerFailedDueToLength.fetch_add(1, std::memory_order_relaxed);
        } catch (...) {}
        return {"", {}};
      }

      auto coordsOpt = fastMapGlobalToLocal(currentPos);
      if (!coordsOpt) {
        throw std::runtime_error("Failed to map position " + std::to_string(currentPos) + 
                             " to block coordinates for node " + std::string(nodeId));
      }

      auto [blockId, nucPos, gapPos] = *coordsOpt;
      char c = getCharAtPosition(nodeId, blockId, nucPos, gapPos);

      if (reverse) {
        rawKmer.insert(0, 1, c);
        if (c != '-' && c != 'x') {
          ungappedKmer.insert(0, 1, c);
          positions.insert(positions.begin(), currentPos);
          nonGapCount++;
        }
      } else {
        rawKmer.push_back(c);
        if (c != '-' && c != 'x') {
          ungappedKmer.push_back(c);
          positions.push_back(currentPos);
          nonGapCount++;
        }
      }
      currentPos += direction;
    }
    return {ungappedKmer, positions};
  }
  
  // Process bulk characters to extract k-mer
  if (reverse) {
    // Process from end to beginning for reverse mode
    for (int i = bulkChars.size() - 1; i >= 0 && nonGapCount < k; --i) {
      char c = bulkChars[i];
      rawKmer.insert(0, 1, c);
      if (c != '-' && c != 'x') {
        ungappedKmer.insert(0, 1, c);
        positions.insert(positions.begin(), bulkPositions[i]);
        nonGapCount++;
      }
    }
  } else {
    // Process from beginning for forward mode
    for (size_t i = 0; i < bulkChars.size() && nonGapCount < k; ++i) {
      char c = bulkChars[i];
      rawKmer.push_back(c);
      if (c != '-' && c != 'x') {
        ungappedKmer.push_back(c);
        positions.push_back(bulkPositions[i]);
        nonGapCount++;
      }
    }
  }

  if (nonGapCount < k) {
    // Failed to get a full k-mer
    try {
      if (rawKmer.length() > 0 && nonGapCount == 0) {
        kmerFailedDueToOnlyGaps.fetch_add(1, std::memory_order_relaxed);
      } else {
        kmerFailedDueToLength.fetch_add(1, std::memory_order_relaxed);
      }
    } catch (...) {}
    
    return {"", {}};
  }

  // Count successful extraction
  try {
    successfulKmerExtractions.fetch_add(1, std::memory_order_relaxed);
  } catch (...) {}

  // Return the ungapped k-mer but preserve the original positions vector
  // which allows correct tracking of gapped end positions
  return {ungappedKmer, positions};
}

void StateManager::initializeBlockRangeMappings() {
  std::unique_lock<std::shared_mutex> lock(nodeMutex);
  if (blockRangeMappingsInitialized) return;
  
  blockRangeMappings.clear();
  blockRangeMappings.reserve(blockRanges.size());
  
  // Create mappings directly (no need to sort as blocks are already sorted by start position)
  for (const auto& [blockId, range] : blockRanges) {
    blockRangeMappings[blockId] = range;
  }
  
  blockRangeMappingsInitialized = true;
  std::cout << "Initialized " << blockRangeMappings.size() << " block range mappings" << std::endl;
}

std::string StateManager::extractSequence(const std::string& nodeId, int64_t start, 
                                 int64_t end, bool skipGaps, bool needLock) const {
    
  // Get a read lock if needed
  std::shared_lock<std::shared_mutex> lock(nodeMutex, std::defer_lock);
  if (needLock) {
    lock.lock();
  }

  // Validate the range
  if (start < 0 || start >= end) {
    throw std::invalid_argument("Invalid range [" + std::to_string(start) + 
                               ", " + std::to_string(end) + ") for node " + nodeId);
  }

  int64_t extractionLength = end - start;
  std::string result;
  result.reserve(extractionLength); // initial reservation, might need more for gaps
  
  try {
    // Log the extraction request with minimal info
    logging::debug("Extracting: node={} range=[{},{})", nodeId, start, end);
    
    // Get the node state to check which blocks are active
    const auto& nodeState = getNodeState(nodeId);
    
    // Get ALL block ranges that overlap with our extraction range, regardless of active status
    std::vector<std::pair<int32_t, coordinates::CoordRange>> overlappingBlocks;
    
    for (const auto& [blockId, blockRange] : blockRanges) {
      // Check if this block overlaps with the requested range
      if (blockRange.end <= start || blockRange.start >= end) {
        continue; // No overlap
      }
      
      // Calculate overlap between block and target range
      coordinates::CoordRange overlapRange{
        std::max(blockRange.start, start),
        std::min(blockRange.end, end)
      };
      
      overlappingBlocks.emplace_back(blockId, overlapRange);
    }
    
    if (overlappingBlocks.empty()) {
      throw std::runtime_error("No blocks found for extraction range [" + 
                              std::to_string(start) + ", " + std::to_string(end) + ")");
    }
    
    // Sort blocks by start position
    std::sort(overlappingBlocks.begin(), overlappingBlocks.end(),
           [](const auto& a, const auto& b) {
             return a.second.start < b.second.start;
           });
           
    // Verify blocks cover the entire extraction range without gaps
    int64_t covered = start;
    for (const auto& [blockId, blockRange] : overlappingBlocks) {
      if (blockRange.start > covered) {
        logging::err("Block coverage gap: pos={} to block_start={}", covered, blockRange.start);
        throw std::runtime_error("Gap in block coverage for extraction range");
      }
      covered = std::max(covered, blockRange.end);
    }
    
    if (covered < end) {
      logging::err("Incomplete coverage: blocks_end={}, extraction_end={}", covered, end);
      throw std::runtime_error("Blocks don't cover entire extraction range");
    }

    // Count statistics to help diagnose character issues
    int gapCount = 0;
    int charCount = 0;
    int inactiveBlockGaps = 0;
    
    // Process each position in the range
    for (int64_t pos = start; pos < end; ++pos) {
      try {
        // Map the global position to local coordinates
        auto maybeCoords = fastMapGlobalToLocal(pos);
        
        if (!maybeCoords) {
          logging::err("Position {} mapping failed", pos);
          throw std::runtime_error("Position mapping failed");
        }
        
        auto [blockId, nucPos, gapPos] = *maybeCoords;
        
        // Check if this block is active for the node
        bool blockActive = isBlockOn_unsafe(nodeId, blockId);
        
        // If block is inactive, add a gap character and continue
        if (!blockActive && !skipGaps) {
          result.push_back('-');
          gapCount++;
          inactiveBlockGaps++;
          continue;
        } else if (!blockActive && skipGaps) { // If block is inactive and we skip gaps, do nothing for this pos
            gapCount++; // Still count it as a conceptual gap encountered
            continue;
        }
        
        // Get the character at this position (for active blocks)
        char c = getCharAtPosition(nodeId, blockId, nucPos, gapPos);

        if (skipGaps && (c == '-' || c == 'x')) {
          gapCount++;
          continue;
        } else if (c == '-' || c == 'x') { // Not skipping gaps, but it is a gap char
            gapCount++;
        }

        result.push_back(c);
        if (c != '-' && c != 'x') { // Count non-gap characters added to result
            charCount++;
        }
        
      } catch (const std::exception& e_inner) { // Changed variable name
        logging::err("Process error at pos={}: {}", pos, e_inner.what()); // Used e_inner.what()
        throw;
      }
    }
    
    // Log summary statistics to help diagnose gap issues
    double gapPercentage = (extractionLength > 0) ? (100.0 * gapCount / extractionLength) : 0.0;
    logging::info("Extraction complete: len={}, gaps={} ({:.1f}%), chars={}, inactive_block_gaps={}",
                 extractionLength, gapCount, gapPercentage, charCount, inactiveBlockGaps);
    
  } catch (const std::exception& e_outer) { // Changed variable name
    logging::err("Extraction failed: {}", e_outer.what()); // Used e_outer.what()
    throw;
  }
  
  return result;
}

// Optimize memory usage by shrinking gap list arrays that are mostly zeros
void StateManager::optimizeGapListMemory() {
  size_t shrunkArrays = 0;
  size_t bytesFreed = 0;
  
  for (size_t blockId = 0; blockId < gapListLengthArray.size(); ++blockId) {
    auto& innerArray = gapListLengthArray[blockId];
    
    // Skip small arrays - not worth optimizing
    if (innerArray.size() < 100) {
            continue;
        }
        
    // Find the last non-zero element position
    int64_t lastNonZeroPos = -1;
    for (int64_t i = static_cast<int64_t>(innerArray.size()) - 1; i >= 0; --i) {
      if (innerArray[i] > 0) {
        lastNonZeroPos = i;
        break;
      }
    }
    
    // If we found no non-zero elements, we can shrink to minimum size
    if (lastNonZeroPos == -1) {
      size_t oldSize = innerArray.size();
      size_t newSize = 1; // Minimum size to avoid complete deallocation
      bytesFreed += (oldSize - newSize) * sizeof(size_t);
      
      // Create a new vector with the smaller size and swap
      std::vector<size_t> newVector(newSize, 0);
      innerArray.swap(newVector);
      
      shrunkArrays++;
            continue;
        }
        
    // If the last non-zero is less than 50% of the capacity, shrink the array
    if (static_cast<size_t>(lastNonZeroPos + 1) < innerArray.size() / 2) {
      size_t oldSize = innerArray.size();
      size_t newSize = static_cast<size_t>(lastNonZeroPos + 1 + 16); // Keep a small buffer
      bytesFreed += (oldSize - newSize) * sizeof(size_t);
      
      // Resize to just beyond the last non-zero element
      innerArray.resize(newSize);
      innerArray.shrink_to_fit();
      
      shrunkArrays++;
    }
  }
  
  if (shrunkArrays > 0) {
    // Convert to MB for human-readable output
    double mbFreed = static_cast<double>(bytesFreed) / (1024.0 * 1024.0);
    logging::info("Optimized gap list memory usage: shrunk {} arrays, freed approximately {:.2f} MB", 
                 shrunkArrays, mbFreed);
  }
}

// Method to flush any pending coordinate mappings
void StateManager::flushCoordinateMappings() {
  // TODO: Implement if necessary - currently, mappings are registered directly.
  // This function might be a remnant of a previous batching approach.
  // For now, it does nothing to satisfy the linker.
  logging::debug("StateManager::flushCoordinateMappings() called - currently no-op.");
}


// Helper to build a topologically sorted list of nodes (parents before children)
std::vector<std::string> StateManager::buildTopologicalOrder(panmanUtils::Tree *tree) {
  std::vector<std::string> sortedNodes;
  
  if (!tree || !tree->root) {
    logging::err("Invalid tree in buildTopologicalOrder");
    return sortedNodes;
  }
  
  // Map to track visited status (0=unvisited, 1=visiting, 2=visited)
  std::unordered_map<std::string, int> visitStatus;
  
  // Pre-populate the map with all nodes to handle disconnected nodes
  for (const auto& [nodeId, _] : tree->allNodes) {
    visitStatus[nodeId] = 0; // Mark all as unvisited initially
  }
  
  // Stack for non-recursive DFS (faster and safer than recursion)
  std::vector<std::pair<std::string, bool>> stack; // (nodeId, processed)
  stack.emplace_back(tree->root->identifier, false);
  
  // Process nodes in DFS order
  while (!stack.empty()) {
    auto [nodeId, processed] = stack.back();
    stack.pop_back();
    
    if (processed) {
      // Node is fully processed, add to sorted list
      sortedNodes.push_back(nodeId);
                continue;
            }
            
    // Mark as being processed
    visitStatus[nodeId] = 1;
    
    // Add node again but marked as processed (will be added to result later)
    stack.emplace_back(nodeId, true);
    
    // Find children through the hierarchy map
    auto hierarchyIt = nodeHierarchy.find(nodeId);
    if (hierarchyIt != nodeHierarchy.end()) {
      const auto& children = hierarchyIt->second.childrenIds;
      
      // Process children in reverse order (to maintain original order in final result)
      for (auto it = children.rbegin(); it != children.rend(); ++it) {
        const auto& childId = *it;
        
        // Skip if already processed or being processed
        if (visitStatus[childId] == 2) {
          continue; // Already fully processed
        }
        
        if (visitStatus[childId] == 1) {
          // Being processed - indicates a cycle which shouldn't happen in a tree
          logging::warn("Cycle detected in node hierarchy during topological sort: {} -> {}", nodeId, childId);
          continue; // Skip this child to avoid infinite loop
        }
        
        // Add child to stack
        stack.emplace_back(childId, false);
      }
    }
  }
  
  // Reverse the result to get parents before children (topological order)
  std::reverse(sortedNodes.begin(), sortedNodes.end());
  
  logging::info("Built topological order for {} nodes", sortedNodes.size());
  return sortedNodes;
}

// Get all seeds in a specific block
const absl::flat_hash_set<int64_t>& StateManager::getBlockSeeds(int32_t blockId) const {
  static const absl::flat_hash_set<int64_t> emptySet;
  auto it = blockToSeeds.find(blockId);
  if (it != blockToSeeds.end()) {
    return it->second;
  }
  return emptySet;
}

// ... existing code ...
void diagnoseBlockStatus(const NodeState& nodeState, const std::string& nodeId, int32_t blockId) {
    logging::debug("Block status diagnosis for block {} in node {}:", blockId, nodeId);
    
    // Check local state
    bool locallyActive = nodeState.activeBlocks.contains(blockId);
    logging::debug("  - Locally active: {}", locallyActive ? "YES" : "NO");
    
    // Check orientation if available
        auto orientIt = nodeState.blockOrientation.find(blockId);
        if (orientIt != nodeState.blockOrientation.end()) {
            if (orientIt->second.has_value()) {
            logging::debug("  - Local orientation: {}", orientIt->second.value() ? "FORWARD" : "INVERTED");
            } else {
            logging::debug("  - Local orientation: UNSPECIFIED (NULL)");
            }
        } else {
        logging::debug("  - Local orientation: NOT SET");
    }
    
    // Log parent ID
    if (!nodeState.parentId.empty()) {
        logging::debug("  - Parent node: {}", nodeState.parentId);
        } else {
        logging::debug("  - No parent (root node)");
    }
}

// Restore the seed storage initialization method with enhanced diagnostics
void StateManager::initializeSeedStorage() {
    std::unique_lock<std::shared_mutex> lock(nodeMutex);
    
    // Log some diagnostics about the node hierarchy
    logging::info("SEED_INIT: Initializing seed storage for {} nodes", nodeHierarchy.size());
    
    // Track successful initializations
    size_t createdRootStores = 0;
    size_t createdChildStores = 0;
    std::vector<std::string> nodesWithoutStores;
    
    // Initialize hierarchical seed storage for all nodes
    for (auto& [nodeId, hierarchy] : nodeHierarchy) {
        // Create a new seed store if it doesn't exist
        if (!hierarchy.seedStore) {
            // If node has a parent, use parent's seed store as base
            if (!hierarchy.parentId.empty()) {
                auto parentIt = nodeHierarchy.find(hierarchy.parentId);
                if (parentIt != nodeHierarchy.end() && parentIt->second.seedStore) {
                    hierarchy.seedStore = std::make_shared<HierarchicalSeedStore>(parentIt->second.seedStore);
                    logging::info("SEED_INIT: Created child seed store for node {} with parent {}", 
                                 nodeId, hierarchy.parentId);
                    createdChildStores++;
    } else {
                    // Parent exists but doesn't have a seed store yet
                    // Create a root-level store for this node instead
                    hierarchy.seedStore = std::make_shared<HierarchicalSeedStore>();
                    logging::warn("SEED_INIT: Parent {} of node {} has no seed store, created root store instead",
                                hierarchy.parentId, nodeId);
                    createdRootStores++;
                }
            } else {
                // For root nodes (or if parent ID is empty), create a root-level store
                hierarchy.seedStore = std::make_shared<HierarchicalSeedStore>();
                logging::info("SEED_INIT: Created root seed store for node {} (no parent)", nodeId);
                createdRootStores++;
            }
        } else {
            logging::debug("SEED_INIT: Seed store already exists for node {}", nodeId);
        }
        
        // Verify the seed store was created successfully
        if (!hierarchy.seedStore) {
            nodesWithoutStores.push_back(nodeId);
        }
    }
    
    // Print summary
    logging::info("SEED_INIT: Successfully initialized seed stores: {} root stores, {} child stores",
                createdRootStores, createdChildStores);
    
    // Log any failures
    if (!nodesWithoutStores.empty()) {
        logging::err("SEED_INIT: Failed to create seed stores for {} nodes:", nodesWithoutStores.size());
        for (size_t i = 0; i < std::min(nodesWithoutStores.size(), size_t(5)); i++) {
            logging::err("  - {}", nodesWithoutStores[i]);
        }
        if (nodesWithoutStores.size() > 5) {
            logging::err("  - ... and {} more", nodesWithoutStores.size() - 5);
        }
    }
}

// Add node-specific setSeedAtPosition method
void StateManager::setSeedAtPosition(std::string_view nodeId_sv, int64_t pos, const seeding::seed_t& seed) {
    if (pos < 0) {
        logging::warn("Invalid position {} in setSeedAtPosition for node {}", pos, nodeId_sv);
        return;
    }
    
    std::string strNodeId(nodeId_sv);
    
    // THREAD SAFETY FIX: Use atomic node initialization with proper synchronization
    // This ensures only one thread initializes a node and prevents race conditions
    
    // Use a static mutex specifically for node initialization to avoid deadlocks
    static std::mutex initMutex;
    
    // Check if node exists (with proper locking)
    {
        std::shared_lock<std::shared_mutex> readLock(nodeMutex);
        auto hierIt = nodeHierarchy.find(strNodeId);
        if (hierIt != nodeHierarchy.end()) {
            // Node exists - fast path, set seed directly
            if (hierIt->second.seedStore) {
                hierIt->second.seedStore->set(pos, seed);
                logging::debug("Set seed at pos={} in existing node {}", pos, strNodeId);
                return;
            }
        }
    }
    
    // Node doesn't exist or needs initialization - use init mutex to serialize this
    {
        std::lock_guard<std::mutex> initLock(initMutex);
        
        // Double-check pattern: verify node still needs initialization
        {
            std::shared_lock<std::shared_mutex> readLock(nodeMutex);
            auto hierIt = nodeHierarchy.find(strNodeId);
            if (hierIt != nodeHierarchy.end() && hierIt->second.seedStore) {
                // Another thread initialized it - set seed and return
                hierIt->second.seedStore->set(pos, seed);
                logging::debug("Set seed at pos={} in node {} (initialized by another thread)", pos, strNodeId);
                return;
            }
        }
        
        // Safe to initialize node
        try {
            initializeNode(strNodeId);
            logging::debug("Initialized node {} during setSeedAtPosition", strNodeId);
        } catch (const std::exception& e) {
            logging::warn("Failed to initialize node {} during setSeedAtPosition: {}", strNodeId, e.what());
            return;
        }
    }
    
    // Now set the seed with proper write lock
    std::shared_lock<std::shared_mutex> readLock(nodeMutex);
    
    // Get the node's hierarchical seed store
    auto hierIt = nodeHierarchy.find(strNodeId);
    if (hierIt == nodeHierarchy.end()) {
        logging::warn("Node {} not found in hierarchy after initialization", strNodeId);
        return;
    }
    
    // Ensure seed store exists - this should be safe now that initialization is atomic
    if (!hierIt->second.seedStore) {
        // Need to upgrade to write lock to create seed store
        readLock.unlock();
        std::unique_lock<std::shared_mutex> writeLock(nodeMutex);
        
        // Double-check pattern again
        hierIt = nodeHierarchy.find(strNodeId);
        if (hierIt == nodeHierarchy.end()) {
            logging::warn("Node {} disappeared during seed store creation", strNodeId);
            return;
        }
        
        if (!hierIt->second.seedStore) {
            // Create new seed store, inheriting from parent if possible
            std::string parentId = hierIt->second.parentId;
            if (!parentId.empty()) {
                auto parentIt = nodeHierarchy.find(parentId);
                if (parentIt != nodeHierarchy.end() && parentIt->second.seedStore) {
                    hierIt->second.seedStore = std::make_shared<HierarchicalSeedStore>(parentIt->second.seedStore);
                    logging::debug("Created seed store for node {} using parent {}", strNodeId, parentId);
                } else {
                    hierIt->second.seedStore = std::make_shared<HierarchicalSeedStore>();
                    logging::debug("Created root-level seed store for node {}", strNodeId);
                }
            } else {
                hierIt->second.seedStore = std::make_shared<HierarchicalSeedStore>();
                logging::debug("Created root-level seed store for node {}", strNodeId);
            }
        }
        
        // Downgrade to read lock for setting the seed
        writeLock.unlock();
        readLock.lock();
        hierIt = nodeHierarchy.find(strNodeId);
    }
    
    // Now set the seed in the node's store
    if (hierIt != nodeHierarchy.end() && hierIt->second.seedStore) {
        hierIt->second.seedStore->set(pos, seed);
        
        // Also update materialized state if it exists
        auto nodeStateIt = nodeStates.find(strNodeId);
        if (nodeStateIt != nodeStates.end()) {
            nodeStateIt->second.setMaterializedSeed(pos, seed);
        }
        
        logging::debug("Set seed at pos={} in node {}", pos, strNodeId);
    } else {
        logging::warn("Failed to set seed at pos={} in node {} - no seed store", pos, strNodeId);
    }
}

// Add node-specific clearSeedAtPosition method
void StateManager::clearSeedAtPosition(std::string_view nodeId, int64_t pos) {
    if (pos < 0) {
        logging::warn("Invalid position {} in clearSeedAtPosition for node {}", pos, nodeId);
        return;
    }
    
    std::unique_lock<std::shared_mutex> writeLock(nodeMutex);
    std::string strNodeId(nodeId);
    
    // Get the node's hierarchical seed store
    auto hierIt = nodeHierarchy.find(strNodeId);
    if (hierIt != nodeHierarchy.end() && hierIt->second.seedStore) {
        // Use the remove method of the hierarchical seed store
        hierIt->second.seedStore->remove(pos);
        
        // Also update materialized state if it exists
        auto nodeStateIt = nodeStates.find(strNodeId);
        if (nodeStateIt != nodeStates.end()) {
            // Remove lock to prevent deadlock - this is just for cleanup
            nodeStateIt->second.materializedSeeds.erase(pos);
        }
        
        logging::debug("SEED_CLEAR: Cleared seed at pos={} in node {}", pos, strNodeId);
        return;
    } else {
        if (hierIt == nodeHierarchy.end()) {
            logging::warn("SEED_CLEAR: Node {} not found in hierarchy when clearing seed", strNodeId);
        } else {
            logging::warn("SEED_CLEAR: Node {} has no seed store", strNodeId);
        }
    }
    
    logging::warn("Failed to clear seed at position {} for node {}: no valid seed store", pos, strNodeId);
}

// Materialize node state for efficient lookups
void StateManager::materializeNodeState(const std::string& nodeId) {
    std::unique_lock<std::shared_mutex> lock(nodeMutex);
    
    auto nodeStateIt = nodeStates.find(nodeId);
    if (nodeStateIt == nodeStates.end()) {
        logging::warn("materializeNodeState: Node {} not found in nodeStates", nodeId);
        return;
    }
    
    auto& nodeState = nodeStateIt->second;
    
    // Check if already materialized
    if (nodeState.materializedStateComputed) {
        return;
    }
    
    // Get parent state to copy from (parent should already be materialized in topological order)
    const NodeState* parentState = nullptr;
    
    // Debug: Check if node exists in hierarchy
    auto hierIt = nodeHierarchy.find(nodeId);
    if (hierIt != nodeHierarchy.end()) {
        if (hierIt->second.parentId != nodeState.parentId) {
            // FIX: Copy the correct parentId from hierarchy to nodeState
            nodeState.parentId = hierIt->second.parentId;
        }
    }
    
    if (!nodeState.parentId.empty()) {
        auto parentIt = nodeStates.find(nodeState.parentId);
        if (parentIt != nodeStates.end()) {
            parentState = &parentIt->second;
            
            // Parent should already be materialized due to topological initialization order
            if (!parentState->materializedStateComputed) {
                // Recursively materialize parent first
                lock.unlock();
                materializeNodeState(nodeState.parentId);
                lock.lock();
                // Re-get the parent state after potential changes
                auto parentIt2 = nodeStates.find(nodeState.parentId);
                if (parentIt2 != nodeStates.end()) {
                    parentState = &parentIt2->second;
                }
            }
        }
    }
    
    // Materialize character state
    if (nodeState.characterStore) {
        static const absl::flat_hash_map<PositionKey, char, PositionKey::Hash> emptyRootChars;
        nodeState.materializeCharacterState(parentState, nodeState.characterStore.get(), emptyRootChars);
    }
    
    // Materialize seed state
    if (hierIt != nodeHierarchy.end() && hierIt->second.seedStore) {
        nodeState.materializeSeedState(parentState, hierIt->second.seedStore.get());
    } else {
        // Even if node has no local seedStore, it should inherit parent's materialized seeds
        nodeState.materializeSeedState(parentState, nullptr);
    }
}

// Consolidated gap update method for efficient batch operations
void StateManager::consolidatedGapUpdate(const std::string& nodeId, 
                                        const std::vector<coordinates::GapUpdate>& updates) {
  if (updates.empty()) {
    logging::debug("consolidatedGapUpdate: No updates to apply for node {}", nodeId);
    return;
  }
  
  std::unique_lock<std::shared_mutex> lock(nodeMutex);
  
  // Get the node's gap map
  auto nodeStateIt = nodeStates.find(nodeId);
  if (nodeStateIt == nodeStates.end()) {
    logging::err("consolidatedGapUpdate: Node state not found for '{}'", nodeId);
    return;
  }
  
  auto& nodeState = nodeStateIt->second;
  if (!nodeState.gapMap) {
    logging::err("consolidatedGapUpdate: Gap map is null for node '{}'", nodeId);
    return;
  }
  
  logging::debug("consolidatedGapUpdate: Applying {} gap updates for node {}", updates.size(), nodeId);
  
  // Convert coordinates::GapUpdate to gap_map::GapUpdate and apply
  for (const auto& update : updates) {
    try {
      // Convert coordinates::GapUpdate to gap_map::GapUpdate
      // coordinates::GapUpdate has: pos, length, isGapAddition
      // gap_map::GapUpdate is: std::pair<bool, std::pair<int64_t, int64_t>>
      bool isRemoval = !update.isGapAddition;
      int64_t start = update.pos;
      int64_t end = update.pos + update.length - 1;
      gap_map::GapUpdate gapMapUpdate(isRemoval, {start, end});
      
      nodeState.gapMap->applyUpdate(gapMapUpdate);
    } catch (const std::exception& e) {
      logging::err("consolidatedGapUpdate: Failed to apply update at position {} for node {}: {}", 
                  update.pos, nodeId, e.what());
    }
  }
  
  // Invalidate any cached data that depends on gap structure
  resetNodeCache(nodeId);
  
  logging::debug("consolidatedGapUpdate: Successfully applied {} updates for node {}", updates.size(), nodeId);
}

// String overload for getRecompRanges
std::vector<coordinates::CoordRange> StateManager::getRecompRanges(const std::string& nodeId) const {
  return getRecompRanges(std::string_view(nodeId));
}

// String overload for mapGlobalToBlockCoords
coordinates::BlockCoordinate StateManager::mapGlobalToBlockCoords(const std::string& nodeId, int64_t globalPos) const {
  return mapGlobalToBlockCoords(std::string_view(nodeId), globalPos);
}

void StateManager::updateMaterializedSeedsAfterRecomputation(const std::string& nodeId, const std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>& detailedSeedChanges) {
    std::unique_lock<std::shared_mutex> lock(nodeMutex);
    
    auto nodeStateIt = nodeStates.find(nodeId);
    if (nodeStateIt == nodeStates.end()) {
        logging::warn("updateMaterializedSeedsAfterRecomputation: Node {} not found", nodeId);
        return;
    }
    
    auto& nodeState = nodeStateIt->second;
    
    // Apply the detailed seed changes computed during recomputation
    nodeState.updateMaterializedSeedsAfterRecomputation(detailedSeedChanges);
}

} // namespace state

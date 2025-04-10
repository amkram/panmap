#pragma once

#include "coordinates.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "panman.hpp"
#include "seeding.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include <string>
#include <string_view>
#include <tbb/enumerable_thread_specific.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <cctype>
#include <boost/pool/object_pool.hpp>

/**
 * @file state.hpp
 * @brief State management system for PanMAN traversal and decoding
 * 
 * This file defines the state management system that tracks the state of each node
 * in the phylogenetic tree, manages mutations, and provides access to
 * sequence data.
 * 
 * The character data system uses a hierarchical delta-based approach where:
 * 
 * 1. Each node only stores its specific mutations relative to its parent
 * 2. Parent-child relationships are tracked with NodeState.parentId
 * 3. When accessing a character:
 *    - First check the node itself
 *    - Then check ancestors recursively
 *    - Finally fall back to reference sequences
 */

namespace state {

// Forward declarations
class StateManager;
struct NodeState;
struct SequenceWithGaps;
struct SharedBlockData;
struct NodeDelta;
struct Mutation;
class LazySequenceView;
struct GapListKey;

// Define position_key_t 
using position_key_t = std::tuple<int, int, int>;

/**
 * @struct SequenceWithGaps
 * @brief Represents a sequence with gap information
 *
 * This structure stores a sequence and its associated gap map for efficient
 * sequence manipulation that accounts for insertions and deletions.
 */
struct SequenceWithGaps {
  /** @brief The actual sequence data */
  std::string sequence;
  
  /** @brief Gap map associated with this sequence */
  std::shared_ptr<gap_map::GapMap> gapMap;
  
  /** @brief Default constructor */
  SequenceWithGaps() : gapMap(std::make_shared<gap_map::GapMap>()) {}
  
  /** @brief Constructor with sequence */
  SequenceWithGaps(const std::string& seq) 
    : sequence(seq), gapMap(std::make_shared<gap_map::GapMap>()) {}
  
  /** @brief Constructor with sequence and gap map */
  SequenceWithGaps(const std::string& seq, std::shared_ptr<gap_map::GapMap> gaps)
    : sequence(seq), gapMap(gaps) {}
    
  /**
   * @brief Remove a gap at the specified position
   * 
   * @param position Position where the gap starts
   * @param length Length of the gap to remove
   * @return true if successful, false otherwise
   */
  bool removeGap(int64_t position, int64_t length) {
    if (position < 0 || length <= 0) {
      return false;
    }

    if (!gapMap) {
      gapMap = std::make_shared<gap_map::GapMap>();
    }
    
    // Create and apply a gap removal update
    gap_map::GapUpdate update(true, std::make_pair(position, length));
    gap_map::applyGapUpdate(*gapMap, update);
    return true;
  }

  /**
   * @brief Insert a gap at the specified position
   * 
   * @param position Position where the gap should be inserted
   * @param length Length of the gap to insert
   * @return true if successful, false otherwise
   */
  bool insertGap(int64_t position, int64_t length) {
    if (position < 0 || length <= 0) {
      return false;
    }
    
    if (!gapMap) {
      gapMap = std::make_shared<gap_map::GapMap>();
    }
    
    // Create and apply a gap insertion update
    gap_map::GapUpdate update(false, std::make_pair(position, length));
    gap_map::applyGapUpdate(*gapMap, update);
    return true;
  }
};


// Removed PackedPositionKey forward declaration

/**
 * @struct PositionKey
 * @brief Represents a 3D coordinate in the block-nucleotide-gap coordinate system
 * - blockId: Identifies the block
 * - nucPos: "Main nucleotide" position within the nucleotide sequence of the block
 * - gapPos: Gap position relative to nucPos (-1 for main nucleotide, 0+ for gap positions)
 * -- Values > 0 index into the optional gap list for the main nucleotide
 */
struct PositionKey {
  int32_t blockId;
  int32_t nucPos;
  int32_t gapPos;

  bool operator==(const PositionKey &other) const {
    return blockId == other.blockId && nucPos == other.nucPos &&
           gapPos == other.gapPos;
  }

  // Helper for creating a position key
  static PositionKey create(int32_t blockId, int32_t nucPos, int32_t gapPos) {
    return {blockId, nucPos, gapPos};
  }
  
  // Removed toPacked method

  /**
   * @struct Hash
   * @brief Hash function for PositionKey to enable use in unordered containers
   */
  struct Hash {
    std::size_t operator()(const PositionKey &k) const {
      // Use standard C++ hashing approach
      std::size_t h1 = std::hash<int32_t>{}(k.blockId);
      std::size_t h2 = std::hash<int32_t>{}(k.nucPos);
      std::size_t h3 = std::hash<int32_t>{}(k.gapPos);
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };
};


/**
 * @struct GapListKey
 * @brief Simplified 2D coordinate for gap list position lookup
 * Maps directly to the blockId and nucPos components of the 3D coordinate system
 */
struct GapListKey {
  int32_t blockId;
  int32_t nucPos;

  // Factory method to create a GapListKey
  static GapListKey create(int32_t blockId, int32_t nucPos) {
    GapListKey key;
    key.blockId = blockId;
    key.nucPos = nucPos;
    return key;
  }
  
  // Equality comparison for hash table usage
  bool operator==(const GapListKey &other) const {
    return blockId == other.blockId && nucPos == other.nucPos;
  }

  // Efficient hash implementation
  struct Hash {
    std::size_t operator()(const GapListKey &k) const {
      // Standard hash combination technique using XOR and bit shift
      return std::hash<int32_t>()(k.blockId) ^ 
             (std::hash<int32_t>()(k.nucPos) << 1);
    }
  };
};

/**
 * @struct CharacterBuffer
 * @brief Thread-local buffer for raw sequence extraction
 * 
 * A buffer for character data and corresponding position information
 * during sequence extraction. It helps avoid frequent memory allocations
 * during high-performance sequence operations.
 */
struct CharacterBuffer {
  std::vector<char> buffer;
  
  /** @brief Corresponding positions for each character */
  std::vector<int64_t> positions;

  CharacterBuffer() {
    // Preallocate with reasonable size
    buffer.reserve(4096);
    positions.reserve(4096);
  }

  void clear() {
    buffer.clear();
    positions.clear();
  }
};

/**
 * @struct UpdateBatch
 * @brief Simplified structure for batched updates
 * 
 * This lightweight structure replaces the more complex UpdateTransaction class,
 * providing the essential functionality without unnecessary abstraction.
 * 
 * @tparam UpdateType The type of update to be batched (e.g., GapUpdate)
 */
template <typename UpdateType>
struct UpdateBatch {
  // Node ID this batch applies to
  std::string nodeId;
  
  // Collection of updates
  std::vector<UpdateType> updates;
  
  // Simple constructor
  UpdateBatch(const std::string &id) : nodeId(id) {}
  
  // Core functionality
  void addUpdate(const UpdateType &update) { updates.push_back(update); }
  void addUpdates(const std::vector<UpdateType> &newUpdates) {
    updates.insert(updates.end(), newUpdates.begin(), newUpdates.end());
  }

  // Basic properties
  size_t size() const { return updates.size(); }
  bool empty() const { return updates.empty(); }
  const std::string &getNodeId() const { return nodeId; }
  const std::vector<UpdateType> &getUpdates() const { return updates; }
};

/**
 * @class HierarchicalGapMap
 * @brief Hierarchical data structure for efficient gap management
 * 
 * A hierarchical gap-run map. Each node has its own gap map
 * that inherits from its parent, enabling local modifications while
 * reusing most of the structure.
 */
class HierarchicalGapMap {
private:
  /** @brief Type alias for shared gap map pointer */
  using GapMapPtr = std::shared_ptr<gap_map::GapMap>;

  /** @brief Parent gap map (weak_ptr to avoid circular references) */
  std::weak_ptr<HierarchicalGapMap> parent;

  /** @brief This node's local gap map representation */
  GapMapPtr localMap;

  /** @brief Flag indicating if this map has been modified */
  bool modified = false;

  /** @brief Backtracking information for this level only */
  std::vector<gap_map::GapUpdate> localBacktrack;

  /** @brief Mutex for thread-safe access to this node's gap map */
  mutable std::mutex mapMutex;
  
  /** @brief Cache for materialized map to avoid recomputation */
  mutable GapMapPtr materializedMapCache;
  mutable bool materializedMapCacheValid = false;

  /**
   * @brief Materialize a complete gap map by merging with ancestors
   * 
   * This method creates a complete map by applying the hierarchical override pattern:
   * 1. Start with the local map
   * 2. Recursively incorporate all parent maps
   * 3. Apply a strict precedence rule: local definitions override parent definitions
   * 
   * This is a deterministic process equivalent to linear inheritance with overriding.
   * 
   * @return Shared pointer to a materialized gap map
   */
  GapMapPtr materializeMap() const {
    std::lock_guard<std::mutex> lock(mapMutex);
    
    // Check if we have a valid cached materialized map
    if (materializedMapCacheValid && materializedMapCache) {
      return materializedMapCache;
    }
    
    // Start with our local map
    auto result = std::make_shared<gap_map::GapMap>(*localMap);

    // If we have a parent, merge with their map
    if (auto parentPtr = parent.lock()) {
      // Get materialized parent map
      auto parentMap = parentPtr->materializeMap();

      // Merge parent entries that don't conflict with our local entries
      for (const auto &entry : *parentMap) {
        if (result->find(entry.first) == result->end()) {
          (*result)[entry.first] = entry.second;
        }
      }
    } else if (!parent.expired()) {
      // Parent pointer exists but cannot be locked (being deleted)
      logging::debug("Parent gap map is being deleted during materialization");
    }
    
    // Cache the result
    materializedMapCache = result;
    materializedMapCacheValid = true;

    return result;
  }

public:
  // Constructor for root map
  HierarchicalGapMap() : localMap(std::make_shared<gap_map::GapMap>()) {}

  // Constructor for child maps
  HierarchicalGapMap(std::shared_ptr<HierarchicalGapMap> parentMap)
      : parent(parentMap), localMap(std::make_shared<gap_map::GapMap>()) {
    // if (!parentMap) {
    //   logging::warn("HierarchicalGapMap created with null parent");
    // }
  }

  /**
   * @brief Check if a position contains a gap
   * 
   * Checks if the specified position falls within any gap range.
   * First checks the local map, then falls back to the parent map if needed.
   * 
   * @param pos Position to check
   * @return true if the position is a gap, false otherwise
   */
  bool isGap(int64_t pos) const {
    if (pos < 0) {
      logging::warn("Negative position {} passed to isGap", pos);
      return false;
    }
    
    // First check our local map
    {
      std::lock_guard<std::mutex> lock(mapMutex);

      auto it = localMap->upper_bound(pos);
      if (it != localMap->begin()) {
        it = std::prev(it);
        int64_t gapStart = it->first;
        int64_t gapLength = it->second;
        
        // Validate gap length
        if (gapLength <= 0) {
          logging::warn("Invalid gap length {} at position {}", gapLength, gapStart);
          return false;
        }
        
        int64_t gapEnd = gapStart + gapLength - 1;

        if (pos >= gapStart && pos <= gapEnd) {
          return true;
        }
      }
    }

    // If not found in local map, check parent map
    if (auto parentPtr = parent.lock()) {
      return parentPtr->isGap(pos);
    }

    return false;
  }

  /**
   * @brief Apply a single gap update
   * 
   * Applies a single gap update to the local map and stores
   * backtracking information.
   * 
   * @param update Gap update to apply
   * @return true if update was successfully applied, false otherwise
   */
  bool applyUpdate(const gap_map::GapUpdate &update) {
    // Validate update
    if (update.second.first < 0) { // position is second.first
      logging::warn("Invalid gap update position: {}", update.second.first);
      return false;
    }
    
    if (update.second.second <= 0 && update.first) { // length is second.second, isGapAddition is first
      logging::warn("Invalid gap update length: {}", update.second.second);
      return false;
    }
    
    std::lock_guard<std::mutex> lock(mapMutex);

    // Mark as modified
    modified = true;
    materializedMapCacheValid = false;

    // Store backtracking info
    localBacktrack.push_back(update);

    // Apply to local map
    gap_map::applyGapUpdate(*localMap, update);
    return true;
  }

  /**
   * @brief Apply multiple gap updates
   * 
   * Applies a batch of gap updates to the local map and stores
   * backtracking information for each one.
   * 
   * @param updates Vector of gap updates to apply
   * @return Number of successfully applied updates
   */
  size_t applyUpdates(const std::vector<gap_map::GapUpdate> &updates) {
    if (updates.empty())
      return 0;

    size_t successCount = 0;
    std::vector<gap_map::GapUpdate> validUpdates;
    
    // Validate updates first
    for (const auto &update : updates) {
      if (update.second.first < 0) { // position is second.first
        logging::warn("Skipping invalid gap update position: {}", update.second.first);
        continue;
      }
      
      if (update.second.second <= 0 && update.first) { // length is second.second, isGapAddition is first
        logging::warn("Skipping invalid gap update length: {}", update.second.second);
        continue;
      }
      
      validUpdates.push_back(update);
      successCount++;
    }
    
    if (validUpdates.empty()) {
      return 0;
    }

    std::lock_guard<std::mutex> lock(mapMutex);

    // Mark as modified
    modified = true;
    materializedMapCacheValid = false;

    // Store backtracking info for valid updates
    localBacktrack.insert(localBacktrack.end(), validUpdates.begin(), validUpdates.end());

    // Apply valid updates to local map
    for (const auto &update : validUpdates) {
      gap_map::applyGapUpdate(*localMap, update);
    }
    
    return successCount;
  }

  /**
   * @brief Find next non-gap position
   * 
   * Skips any gaps and returns the next non-gap position in the
   * specified direction.
   * 
   * @param pos Starting position
   * @param forward Direction (true for forward, false for backward)
   * @param totalCoords Total number of coordinates in the sequence
   * @return Next non-gap position
   */
  int64_t skipGap(int64_t pos, bool forward, int64_t totalCoords) const {
    // Validate inputs
    if (totalCoords <= 0) {
      logging::warn("Invalid totalCoords: {}", totalCoords);
      return 0;
    }
    
    // Clamp the starting position to valid range
    pos = std::max(int64_t(0), std::min(pos, totalCoords - 1));
    
    // Get a materialized map once
    auto map = materializeMap();

    if (forward) {
      // Forward search: find next non-gap position
      while (pos < totalCoords) {
      auto it = map->upper_bound(pos);
        bool inGap = false;
        
      if (it != map->begin()) {
          auto prev = std::prev(it);
          int64_t gapStart = prev->first;
          int64_t gapLength = prev->second;
          
          if (gapLength <= 0) {
            logging::warn("Invalid gap length {} at position {}", gapLength, gapStart);
            break;
          }
          
        int64_t gapEnd = gapStart + gapLength - 1;

          // If we're in a gap, skip to the end of the gap
        if (pos >= gapStart && pos <= gapEnd) {
          pos = gapEnd + 1;
            inGap = true;
          }
        }
        
        // If we're not in a gap, we've found our position
        if (!inGap) break;
      }
      
      // Ensure we don't exceed bounds
      return std::min(pos, totalCoords - 1);
    } else {
      // Backward search: find previous non-gap position
      while (pos > 0) {
      auto it = map->upper_bound(pos);
        bool inGap = false;
        
      if (it != map->begin()) {
          auto prev = std::prev(it);
          int64_t gapStart = prev->first;
          int64_t gapLength = prev->second;
          
          if (gapLength <= 0) {
            logging::warn("Invalid gap length {} at position {}", gapLength, gapStart);
            break;
          }
          
        int64_t gapEnd = gapStart + gapLength - 1;

          // If we're in a gap, skip to before the gap
        if (pos >= gapStart && pos <= gapEnd) {
          pos = gapStart - 1;
            inGap = true;
          }
        }
        
        // If we're not in a gap, we've found our position
        if (!inGap) break;
    }

      // Ensure we don't go below zero
      return std::max(int64_t(0), pos);
    }
  }

  // Get a complete copy of the gap map for this node
  gap_map::GapMap getFullMap() const { return *materializeMap(); }

  // Get the local changes only
  gap_map::GapMap getLocalMap() const {
    std::lock_guard<std::mutex> lock(mapMutex);
    return *localMap;
  }

  // Get backtracking information for this level only
  std::vector<gap_map::GapUpdate> getLocalBacktrack() const {
    std::lock_guard<std::mutex> lock(mapMutex);
    return localBacktrack;
  }

  // Check if this map has local modifications
  bool isModified() const {
    std::lock_guard<std::mutex> lock(mapMutex);
    return modified;
  }

  // Clear local modifications
  void clearModifications() {
    std::lock_guard<std::mutex> lock(mapMutex);
    localMap->clear();
    localBacktrack.clear();
    modified = false;
    materializedMapCacheValid = false;
  }
};

// Add this new template after HierarchicalGapMap
/**
 * @class HierarchicalStore
 * @brief Generic hierarchical storage for any type of data
 * 
 * Provides efficient inheritance of values from parent to child with
 * local overrides at each level.
 * 
 * @tparam T The type of values stored
 */
template <typename T> class HierarchicalStore {
private:
  /** @brief Parent store (weak_ptr to avoid circular references) */
  std::weak_ptr<HierarchicalStore<T>> parent;
  
  /** @brief Local values that override parent values */
  std::unordered_map<int64_t, T> localValues;
  
  /** @brief Mutex for thread-safe access */
  mutable std::mutex storeMutex;
  
  /** @brief Flag indicating if this store has been modified */
  bool modified = false;

public:
  /**
   * @brief Constructor for root store
   */
  HierarchicalStore() {}

  /**
   * @brief Constructor for child stores
   * @param parentStore Parent store to inherit from
   */
  HierarchicalStore(std::shared_ptr<HierarchicalStore<T>> parentStore)
      : parent(parentStore) {
    if (!parentStore) {
      logging::warn("HierarchicalStore created with null parent");
    }
  }

  /**
   * @brief Get value, checking parent if not in local store
   * 
   * @param key The key to look up
   * @return The value if found, std::nullopt otherwise
   */
  std::optional<T> get(int64_t key) const {
    if (key < 0) {
      logging::warn("Negative key {} passed to HierarchicalStore::get", key);
      return std::nullopt;
    }
    
    // First check local values
    {
      std::lock_guard<std::mutex> lock(storeMutex);
      auto it = localValues.find(key);
      if (it != localValues.end()) {
        return it->second;
      }
    }

    // Then check parent
    if (auto parentPtr = parent.lock()) {
      return parentPtr->get(key);
    } else if (!parent.expired()) {
      logging::debug("Parent store is being deleted during get operation");
    }

    return std::nullopt;
  }

  /**
   * @brief Set value only in local store
   * 
   * @param key The key to set
   * @param value The value to store
   * @return true if successful, false otherwise
   */
  bool set(int64_t key, const T &value) {
    if (key < 0) {
      logging::warn("Negative key {} passed to HierarchicalStore::set", key);
      return false;
    }
    
    std::lock_guard<std::mutex> lock(storeMutex);
    localValues[key] = value;
    modified = true;
    return true;
  }

  /**
   * @brief Check if this store has local modifications
   * @return true if local modifications exist, false otherwise
   */
  bool isModified() const {
    std::lock_guard<std::mutex> lock(storeMutex);
    return modified;
  }

  /**
   * @brief Get all local values
   * @return Map of all locally stored key-value pairs
   */
  std::unordered_map<int64_t, T> getLocalValues() const {
    std::lock_guard<std::mutex> lock(storeMutex);
    return localValues;
  }

  /**
   * @brief Apply multiple updates efficiently
   * 
   * @param updates Vector of key-value pairs to update
   * @return Number of successfully applied updates
   */
  size_t batchUpdate(const std::vector<std::pair<int64_t, T>> &updates) {
    if (updates.empty())
      return 0;

    std::lock_guard<std::mutex> lock(storeMutex);
    
    size_t successCount = 0;
    for (const auto &[key, value] : updates) {
      if (key < 0) {
        logging::warn("Skipping negative key {} in batchUpdate", key);
        continue;
      }
      
      localValues[key] = value;
      successCount++;
    }
    
    if (successCount > 0) {
    modified = true;
  }

    return successCount;
  }

  /**
   * @brief Clear local modifications, reverting to parent state
   */
  void clearModifications() {
    std::lock_guard<std::mutex> lock(storeMutex);
    localValues.clear();
    modified = false;
  }
  
  /**
   * @brief Check if a key exists in this store or any parent
   * 
   * @param key The key to check for
   * @return true if the key exists, false otherwise
   */
  bool contains(int64_t key) const {
    if (key < 0) {
      logging::warn("Negative key {} passed to HierarchicalStore::contains", key);
      return false;
    }
    
    // Check local values first
    {
      std::lock_guard<std::mutex> lock(storeMutex);
      if (localValues.find(key) != localValues.end()) {
        return true;
      }
    }
    
    // Then check parent
    if (auto parentPtr = parent.lock()) {
      return parentPtr->contains(key);
    }
    
    return false;
  }
  
  /**
   * @brief Get the number of local entries
   * 
   * @return Count of local entries
   */
  size_t localSize() const {
    std::lock_guard<std::mutex> lock(storeMutex);
    return localValues.size();
  }
};

/**
 * @struct SharedBlockData
 * @brief Shared data for a block across multiple nodes
 *
 * This structure contains shared information for a specific block
 * that can be accessed by multiple nodes in the tree.
 */
struct SharedBlockData {
  /** @brief Reference sequence for this block */
  std::string refSequence;
  
  /** @brief Hierarchical gap map for this block */
  std::shared_ptr<HierarchicalGapMap> gapMap;
  
  /** @brief Global starting position of this block */
  int64_t globalStartPos;
  
  /** @brief Length of this block */
  int64_t length;
  
  /** @brief Whether this block is active */
  bool isActive = true;
  
  /** @brief Whether this block is inverted */
  bool isInverted = false;
  
  /** @brief Mutex for thread-safe access */
  mutable std::mutex mutex;
  
  /** @brief Default constructor */
  SharedBlockData() : gapMap(std::make_shared<HierarchicalGapMap>()),
                      globalStartPos(0), length(0) {}
  
  /** @brief Constructor with reference sequence */
  SharedBlockData(const std::string& seq, int64_t startPos = 0) 
    : refSequence(seq), 
      gapMap(std::make_shared<HierarchicalGapMap>()),
      globalStartPos(startPos),
      length(seq.length()) {}
};

// Define the NodeState structure
struct NodeState {
  // Default constructor required for std::unordered_map::operator[]
  NodeState() = default;

  // Parent reference for hierarchical delta lookup
  std::string parentId;

  std::unordered_set<int32_t> activeBlocks;
  std::unordered_map<int32_t, bool>
      blockStrands; // true = forward, false = inverted
  std::vector<std::tuple<int32_t, int32_t, int32_t, char, char>>
      nucleotideChanges;
  std::vector<coordinates::CoordRange> recompRanges;
  
  // Track boundaries of active blocks for efficient range expansion
  std::optional<int64_t> firstActiveBlockPos;  // First position of any active block
  std::optional<int64_t> lastActiveBlockPos;   // Last position of any active block
  
  // Recorded mutations for debugging
  std::vector<std::pair<int64_t, uint8_t>> recordedMutations;

  // Use efficient key structs instead of string keys
  std::unordered_map<state::PositionKey, char, state::PositionKey::Hash> characterData;
  std::unordered_map<int32_t, bool> blockOnStatus;

  // Fields for gap map updates
  std::vector<coordinates::GapUpdate> gapMapUpdates;

  // Mutation index for direct position lookup
  // Maps position keys to the node ID that last mutated this position
  std::unordered_map<PositionKey, std::string, PositionKey::Hash> mutationIndex;
  
  // Node-specific hierarchical gap map
  std::shared_ptr<HierarchicalGapMap> gapMap;

  // Integrated sequence and gap representation for better cache locality
  std::unordered_map<int32_t, state::SequenceWithGaps> blockSequences;
  mutable std::shared_mutex
      blockSequencesMutex; // For thread-safe access with concurrent reads

  // Compressed gap runs for this node
  std::vector<gap_map::CompressedGapRun> compressedGapRuns;

  // Delta-based storage
  std::shared_ptr<state::NodeDelta> delta; // Delta from parent to this node

  // Shared immutable block data
  std::unordered_map<int32_t, std::shared_ptr<SharedBlockData>> sharedBlockData;

  // Cache for extracted sequences to avoid redundant computation
  mutable std::unordered_map<int64_t, std::string> extractedKmers;
  mutable std::mutex kmerCacheMutex;

  // Record the last 100 position lookups for cache locality
  static constexpr size_t RECENT_LOOKUPS_SIZE = 100;
  static thread_local std::vector<PositionKey> recentLookups;
  
  // Quaternary-encoded seed changes (for efficient storage and transmission)
  std::vector<int64_t> seedChangeBasePositions;  // Base positions for seed changes
  std::vector<uint64_t> seedChangeBitMasks;      // Quaternary-encoded masks
  
  /**
   * @brief Add seed changes in quaternary-encoded format
   * 
   * @param basePositions Base positions for seed changes
   * @param bitMasks Quaternary-encoded masks
   */
  void addSeedChanges(const std::vector<int64_t>& basePositions, 
                     const std::vector<uint64_t>& bitMasks) {
    // Store the base positions and bit masks                 
    seedChangeBasePositions = basePositions;
    seedChangeBitMasks = bitMasks;
  }
  
  /**
   * @brief Get seed changes in decoded form
   * 
   * @return Vector of (position, wasOn, isOn) tuples representing seed changes
   */
  std::vector<std::tuple<int64_t, bool, bool>> getSeedChanges() const {
    // This function depends on an external decoder function, add a forward declaration
    // or direct implementation if needed
    if (seedChangeBasePositions.empty() || seedChangeBitMasks.empty()) {
      return {};
    }
    
    std::vector<std::tuple<int64_t, bool, bool>> decodedChanges;
    decodedChanges.reserve(seedChangeBasePositions.size() * 4); // Rough estimate
    
    for (size_t i = 0; i < seedChangeBasePositions.size(); ++i) {
      int64_t basePos = seedChangeBasePositions[i];
      uint64_t mask = seedChangeBitMasks[i];
      
      // Decode up to 32 positions from each mask
      for (int k = 0; k < 32; ++k) {
        // Extract the quaternary value for this position
        uint8_t value = (mask >> (k * 2)) & 0x3;
        
        // Skip if no change (should be rare in practice)
        if (value == 0) continue;
        
        // Calculate position
        int64_t pos = basePos - k;
        
        // Decode was-on and is-on states from quaternary value
        bool wasOn = (value == 1 || value == 3);  // Deleted or Modified
        bool isOn = (value == 2 || value == 3);   // Added or Modified
        
        // Add decoded change
        decodedChanges.emplace_back(pos, wasOn, isOn);
      }
    }
    
    // Sort by position
    std::sort(decodedChanges.begin(), decodedChanges.end(),
              [](const auto& a, const auto& b) {
                return std::get<0>(a) < std::get<0>(b);
              });
              
    return decodedChanges;
  }
  
  // Helper methods for block status
  bool isBlockOn(int32_t blockId) const {
    auto it = blockOnStatus.find(blockId);
    return it != blockOnStatus.end() && it->second;
  }

  bool isBlockInverted(int32_t blockId) const {
    auto it = blockStrands.find(blockId);
    return it == blockStrands.end() || !it->second;
  }

  void setBlockOn(int32_t blockId, bool on) {
    if (on) {
      activeBlocks.insert(blockId);
    } else {
      activeBlocks.erase(blockId);
    }
    blockOnStatus[blockId] = on;
  }

  void setBlockForward(int32_t blockId, bool forward) {
    blockStrands[blockId] = forward;
  }

  // Convert updates to compressed gap runs for storage efficiency
  void compressGapRuns() {
    compressedGapRuns.clear();

    // Convert from gapMapUpdates to compressed format
    for (const auto &update : gapMapUpdates) {
      if (update.isGapAddition) {
        compressedGapRuns.push_back(
            gap_map::CompressedGapRun::create(update.pos, update.length));
      }
    }
  }

  // Add a gap update to the node's update list
  void addGapUpdate(const coordinates::GapUpdate &update) {
    gapMapUpdates.push_back(update);
  }

  // Apply pending gap updates directly to hierarchical gap map
  void applyGapUpdates() {
    if (!gapMap || gapMapUpdates.empty())
      return;

    std::vector<gap_map::GapUpdate> updates;
    updates.reserve(gapMapUpdates.size());

    for (const auto &update : gapMapUpdates) {
      // For gap_map, true means deletion, false means addition
      bool isRemoval = !update.isGapAddition;
      int64_t start = update.pos;
      int64_t end = update.pos + update.length - 1;

      gap_map::GapRange range = {start, end};
      updates.emplace_back(isRemoval, range);

      // Also update the integrated SequenceWithGaps structures
      // Find which block this gap belongs to
      for (auto &[blockId, seqWithGaps] : blockSequences) {
        // Get the block range
        // This is just a placeholder - in the real implementation,
        // you would need to access the actual block range information
        auto blockRange = sharedBlockData.find(blockId);
        if (blockRange != sharedBlockData.end()) {
          int64_t blockStart = blockRange->second->globalStartPos;
          int64_t blockEnd = blockStart + blockRange->second->length - 1;

          // Check if this gap is within this block
          if (start >= blockStart && start <= blockEnd) {
            // Convert to block-local position
            int64_t localStart = start - blockStart;
            int64_t length = update.length;

            // Update the SequenceWithGaps
            if (isRemoval) {
              seqWithGaps.removeGap(localStart, length);
            } else {
              seqWithGaps.insertGap(localStart, length);
            }
          }
        }
      }
    }

    // Apply all updates at once
    gapMap->applyUpdates(updates);

    // Clear the updates now that they're applied
    gapMapUpdates.clear();
  }

  // Custom move constructor
  NodeState(NodeState &&other) noexcept
      : activeBlocks(std::move(other.activeBlocks)),
        blockStrands(std::move(other.blockStrands)),
        nucleotideChanges(std::move(other.nucleotideChanges)),
        recompRanges(std::move(other.recompRanges)),
        characterData(std::move(other.characterData)),
        blockOnStatus(std::move(other.blockOnStatus)),
        gapMapUpdates(std::move(other.gapMapUpdates)),
        gapMap(std::move(other.gapMap)),
        blockSequences(std::move(other.blockSequences)),
        compressedGapRuns(std::move(other.compressedGapRuns)),
        delta(std::move(other.delta)),
        sharedBlockData(std::move(other.sharedBlockData)),
        extractedKmers(std::move(other.extractedKmers)),
        seedChangeBasePositions(std::move(other.seedChangeBasePositions)),
        seedChangeBitMasks(std::move(other.seedChangeBitMasks)) {
    // Do not move blockSequencesMutex and kmerCacheMutex
  }

  // Custom move assignment operator
  NodeState &operator=(NodeState &&other) noexcept {
    if (this != &other) {
      activeBlocks = std::move(other.activeBlocks);
      blockStrands = std::move(other.blockStrands);
      nucleotideChanges = std::move(other.nucleotideChanges);
      recompRanges = std::move(other.recompRanges);
      characterData = std::move(other.characterData);
      blockOnStatus = std::move(other.blockOnStatus);
      gapMapUpdates = std::move(other.gapMapUpdates);
      gapMap = std::move(other.gapMap);
      blockSequences = std::move(other.blockSequences);
      compressedGapRuns = std::move(other.compressedGapRuns);
      delta = std::move(other.delta);
      sharedBlockData = std::move(other.sharedBlockData);
      extractedKmers = std::move(other.extractedKmers);
      seedChangeBasePositions = std::move(other.seedChangeBasePositions);
      seedChangeBitMasks = std::move(other.seedChangeBitMasks);
      // Do not move blockSequencesMutex and kmerCacheMutex
    }
    return *this;
  }

  NodeState(const NodeState &) = delete;
  NodeState &operator=(const NodeState &) = delete;

  // Create a clone of this NodeState with fresh mutexes
  NodeState clone() const {
    NodeState newState;

    // Copy all the data (but not the mutexes)
    newState.activeBlocks = activeBlocks;
    newState.blockStrands = blockStrands;
    newState.nucleotideChanges = nucleotideChanges;
    newState.recompRanges = recompRanges;
    newState.characterData = characterData;
    newState.blockOnStatus = blockOnStatus;
    newState.gapMapUpdates = gapMapUpdates;
    newState.gapMap = gapMap;
    newState.blockSequences = blockSequences;
    newState.compressedGapRuns = compressedGapRuns;
    newState.delta = delta;
    newState.sharedBlockData = sharedBlockData;
    newState.extractedKmers = extractedKmers;
    newState.seedChangeBasePositions = seedChangeBasePositions;
    newState.seedChangeBitMasks = seedChangeBitMasks;

    return newState;
  }
};

class StateManager {
private:
  // Node states for each processed node
  std::unordered_map<std::string, NodeState> nodeStates;
  
  // Mutex for thread-safe access to nodeStates
  mutable std::mutex nodeStatesMutex;
  
  // Mutex for thread-safe access to mutation index
  mutable std::mutex mutationIndexMutex;

  // Add mutationIndex here to match where it was originally declared
  // Mutation index for direct position lookup
  // Maps position keys to the node ID that last mutated this position
  std::unordered_map<PositionKey, std::string, PositionKey::Hash> mutationIndex;

  // Streamlined hierarchical node structure for efficient traversal
  struct NodeHierarchy {
    // Core hierarchy relationships
    std::string parentId;
    std::vector<std::string> childrenIds;
    
    // Essential data stores for the 3D coordinate system
    std::shared_ptr<HierarchicalGapMap> gapMap;         // Critical for gap position handling
    std::shared_ptr<HierarchicalStore<char>> characterStore; // For nucleotide data
    
    // Seed storage for optimized lookups
    std::shared_ptr<HierarchicalStore<seeding::seed_t>> seedStore;
    int64_t dfsIndex = -1; // DFS index for ancestor/descendant checks
    int64_t subtreeSize = 0; // Size of the subtree rooted at this node
  };

  // Hierarchy management
  std::unordered_map<std::string, NodeHierarchy> nodeHierarchy;

  // Reference block sequences
  std::unordered_map<int32_t, std::string> blockSequences;

  // Block ranges
  std::unordered_map<int32_t, coordinates::CoordRange> blockRanges;
  
  // Replace hash table with a 2D array for faster lookups
  std::unordered_map<GapListKey, size_t, GapListKey::Hash> globalGapListLengths;
  
  // New efficient data structure for gap list lengths
  std::vector<std::vector<size_t>> gapListLengthArray;
  
  // Cache for frequently accessed gap list lengths
  struct GapListLengthCache {
    static constexpr const size_t CACHE_SIZE = 1024 * 5;
    struct CacheEntry {
      int32_t blockId;
      int32_t nucPos;
      size_t length;
      uint64_t timestamp;
    };
    
    CacheEntry entries[CACHE_SIZE];
    uint64_t currentTimestamp = 0;
    
    GapListLengthCache() {
      for (size_t i = 0; i < CACHE_SIZE; i++) {
        entries[i].blockId = -1; // Invalid blockId
      }
    }
    
    void put(int32_t blockId, int32_t nucPos, size_t length) {
      size_t index = (static_cast<size_t>(blockId) * 31 + static_cast<size_t>(nucPos)) % CACHE_SIZE;
      entries[index].blockId = blockId;
      entries[index].nucPos = nucPos;
      entries[index].length = length;
      entries[index].timestamp = currentTimestamp++;
    }
    
    size_t get(int32_t blockId, int32_t nucPos, bool& found) {
      size_t index = (static_cast<size_t>(blockId) * 31 + static_cast<size_t>(nucPos)) % CACHE_SIZE;
      if (entries[index].blockId == blockId && entries[index].nucPos == nucPos) {
        entries[index].timestamp = currentTimestamp++; // Update timestamp on access
        found = true;
        return entries[index].length;
      }
      found = false;
      return 0;
    }
  };
  
  // Thread-local cache for gap list lengths
  mutable tbb::enumerable_thread_specific<GapListLengthCache> gapListLengthCaches;
  
  // Global gap map
  gap_map::GapMap globalGapMap;
  std::vector<gap_map::GapUpdate> gapMapBacktrack;
  mutable std::shared_mutex globalGapMapMutex;

  // Thread-local gap maps for efficient batch processing
  tbb::enumerable_thread_specific<gap_map::GapMap> threadLocalGapMaps;

  // Cache for sorted active blocks
  std::unordered_map<std::string,
                     std::vector<std::pair<int32_t, coordinates::CoordRange>>>
      nodeBlockRangeCache;

  // Number of blocks and coordinates
  size_t numBlocks = 0;
  size_t numCoords = 0;

  // Flag indicating if cache is initialized
  bool globalPosCacheInitialized = false;

  // Global position cache for fast mapping
  struct GlobalPosCache {
    std::vector<int32_t> posToBlockId;      // Maps global position to block ID
    std::vector<int64_t> blockStartOffsets; // Starting offset of each block
    int64_t maxPosition = 0;                // Maximum global position

    void clear() {
      posToBlockId.clear();
      blockStartOffsets.clear();
      maxPosition = 0;
    }

    bool isValid(int64_t pos) const { return pos >= 0 && pos < maxPosition; }

    int32_t getBlockId(int64_t pos) const {
      if (!isValid(pos) || static_cast<size_t>(pos) >= posToBlockId.size()) {
        return -1;
      }
      return posToBlockId[pos];
    }
  } globalPosCache;

  // Node DFS indices for ancestral nodes discovery
  std::unordered_map<std::string, int64_t> nodeDfsIndices;

  // Syncmer k-mer size
  int16_t kmerSize = 0;

  // Syncmer s-mer size
  int smerSize = 0; // Default s-mer size

  // Sequence cache for expensive operations
  static constexpr size_t MAX_SEQUENCE_CACHE_SIZE =
      100000; // Maximum number of cached sequences
  mutable std::shared_mutex
      sequenceCacheMutex; // Mutex for thread-safe cache access

  // Seed at each position (if any)
  std::vector<std::optional<seeding::seed_t>> positionSeeds;

  // Mapping from blocks to seed positions
  std::unordered_map<int32_t, std::unordered_set<int64_t>> blockToSeeds;

  // Cache for extracted sequences
  std::unordered_map<std::string, std::pair<std::string, std::vector<int64_t>>>
      sequenceCache;

  // Thread-local character buffers
  tbb::enumerable_thread_specific<CharacterBuffer> charBuffers;

  // 7. Shared Immutable Block Data
  std::unordered_map<int32_t, std::shared_ptr<SharedBlockData>>
      sharedBlockDataMap;

  // Root gap map (shared by all nodes)
  std::shared_ptr<HierarchicalGapMap> rootGapMap;

  // Node relationship tracking for gap map inheritance
  std::unordered_map<std::string, std::string> nodeParents;

  // Transaction management
  std::unordered_map<std::string,
                     std::shared_ptr<UpdateBatch<coordinates::GapUpdate>>> activeTransactions;
  std::mutex transactionMutex;

  // We use the hierarchical gap map for efficient sequence traversal

  // Gap statistics for adaptive optimization
  gap_map::GapStatistics gapStats;
  
  // Coordinate mapping between 3D coordinates and scalar global positions
  std::unordered_map<PositionKey, int64_t, PositionKey::Hash> blockCoordToGlobalPos;
  // Maps blockId, nucPos, gapPos to string index in sequences
  std::unordered_map<PositionKey, int32_t, PositionKey::Hash> blockCoordToStringIndex;
  
  // Method to update block coordinate mappings using block ranges as the source of truth
  void updateBlockCoordMapping();

  // Helper method to get a sorted list of active blocks
  std::vector<std::pair<int32_t, coordinates::CoordRange>>
  getSortedActiveBlocks(std::string_view nodeId) const;

  // Helper method to merge overlapping ranges
  std::vector<coordinates::CoordRange>
  mergeRanges(std::vector<coordinates::CoordRange> &ranges);


  // Helper for DFS index computation
  void computeDfsIndices(const std::string &nodeId, int64_t &index);

  // Helper to check if one node is a descendant of another
  bool isDescendant(const std::string &ancestorId,
                    const std::string &nodeId) const;

  // Private methods
  // Helper method to log a single position's mapping (used by logGapListsForBlock)
  void logPositionMapping(const std::string& nodeId, int32_t blockId, 
                         int32_t nucPos, const std::string& sequence) const;
                         
  // Helper method to find and skip over gap runs efficiently
  int64_t skipGapRun(std::string_view nodeId, int64_t pos, bool forward) const;
  
  // Helper to diagnose block boundary issues for recomputation ranges
  void diagnoseBoundaryBlocks(std::string_view nodeId, int32_t blockId, 
                             const coordinates::CoordRange& range, 
                             const std::string& context);

  // Helper to initialize block boundaries
  void initializeBlockBoundaries(const std::string& nodeId);

  // Helper function to update block gap mappings after gap changes
  void updateBlockGapMappings(const std::string &nodeId);

  // Efficiently process seeds in a range to avoid redundant sequence extraction
  void processSeedsInRange(
      std::string_view nodeId, const coordinates::CoordRange &range, int k, int s,
      const std::function<void(int64_t startPos, int64_t endPos, std::optional<seeding::seed_t> &seed, std::string_view kmerView)>
          &processFn);

public:
  // Constructor
  StateManager();
  StateManager(size_t numCoordinates);
  
  // Initialize the state manager with a tree
  void initialize(panmanUtils::Tree *tree);

  // Initialize node hierarchy tree structure
  void initializeNodeHierarchy(panmanUtils::Tree *tree,
                               panmanUtils::Node *rootNode);
  

  // Initialize node using stored DFS indices
  void initializeNode(const std::string &nodeId);

  // Initialize sequence data
  void initializeSequenceData(panmanUtils::Tree *tree,
                              panmanUtils::Node *rootNode);

  // Get character at position
  char getCharAtPosition(std::string_view nodeId, int32_t blockId,
                        int32_t nucPos, int32_t gapPos) const;

  // Get active blocks for a node
  const std::unordered_set<int32_t> &
  getActiveBlocks(std::string_view nodeId) const;

  // Set character at position, returns true if successful
  bool setCharAtPosition(std::string_view nodeId, int32_t blockId,
                         int32_t nucPos, int32_t gapPos, char c);

  // Methods to check block status
  bool isBlockOn(std::string_view nodeId, int32_t blockId) const;
  bool isBlockInverted(std::string_view nodeId, int32_t blockId) const;
  void setBlockOn(std::string_view nodeId, int32_t blockId, bool on);
  void setBlockForward(std::string_view nodeId, int32_t blockId, bool forward);
  void setBlockInverted(std::string_view nodeId, int32_t blockId,
                        bool inverted);


  // Access to node state
  NodeState &getNodeState(const std::string &nodeId);
  const NodeState &getNodeState(const std::string &nodeId) const;

  // Get recomputation ranges for a node
  const std::vector<coordinates::CoordRange>&
  getRecompRanges(std::string_view nodeId) const;

  // Methods to handle block mutations
  bool applyBlockMutation(std::string_view nodeId, int32_t blockId,
                          bool isInsertion, bool isInversion);
  
  // Overloaded version that accepts NodeState directly for better performance
  void applyNucleotideMutation(const std::string& nodeId, NodeState& nodeState,
                               int32_t blockId, int32_t nucPos, int32_t gapPos, char newChar);
  
  // Process multi-position nucleotide mutations
  void processMultiPositionMutation(std::string_view nodeId, int32_t blockId,
                                  int32_t nucPos, int32_t gapPos,
                                  int mutationType, int length,
                                  uint32_t packedNucs);
  
  // Check if ranges overlap
  bool hasOverlap(const coordinates::CoordRange &range1,
                  const coordinates::CoordRange &range2) const;

  // Get block range
  coordinates::CoordRange getBlockRange(int32_t blockId) const;

  // Set block range
  void setBlockRange(int32_t blockId, const coordinates::CoordRange &range);

  // Get gap list length for a position - faster implementation
  size_t getGapListLength(int32_t blockId, int32_t nucPos) const;

  // Extract sequence from a node
  std::pair<std::string, std::vector<int64_t>>
  extractSequence(std::string_view nodeId, const coordinates::CoordRange &range,
                  bool skipGaps = true);

  // Helper method to find common ancestor
  std::string findLowestCommonAncestor(const std::string &nodeA,
                                       const std::string &nodeB);

  // Helper method to find path to node
  std::vector<std::string> findPathToNode(const std::string &fromNode,
                                          const std::string &toNode);

  // Helper method to propagate state
  void propagateState(const std::string &fromNode, const std::string &toNode);

  // Fast lookup from global position to (blockId, nucPos, gapPos)
  // Note: This method does NOT handle inversions - caller must adjust nucPos if
  // needed This is by design to optimize for the fixed global coordinate system
  std::optional<std::tuple<int32_t, int32_t, int32_t>>
  fastMapGlobalToLocal(int64_t globalPos) const;

  // Fast lookup from (blockId, nucPos, gapPos) to global position
  int64_t fastMapLocalToGlobal(int32_t blockId, int32_t nucPos,
                               int32_t gapPos) const;

  // Initialize DFS indices for traversal
  void initializeDfsIndices(panmanUtils::Tree *tree);

  // Recursive helper for DFS index initialization
  void initializeDfsIndices(panmanUtils::Tree *tree, panmanUtils::Node *node,
                            int64_t &dfsIndex);

  // Get the DFS index for a node
  int64_t getDfsIndex(const std::string &nodeId) const;

  void setKmerSize(int16_t k);

  void setSmerSize(int s);

  int16_t getKmerSize() const;

  int getSmerSize() const;

  // Add a recomputation range for a node
  void addRecompRange(const std::string &nodeId,
                      const coordinates::CoordRange &range);

  // Expand recomputation ranges to include k before and after
  void expandRecompRanges(const std::string &nodeId);

  // Backtrack node state to ancestor nodes
  void backtrackNode(const std::string &nodeId);

  // Get number of blocks
  size_t getNumBlocks() const;

  // Set the number of blocks
  void setNumBlocks(size_t blockCount);

  // Seed management methods
  std::optional<seeding::seed_t> getSeedAtPosition(int64_t pos) const;
  void setSeedAtPosition(int64_t pos, const seeding::seed_t &seed);
  void clearSeedAtPosition(int64_t pos);
  void addSeedToBlock(int32_t blockId, int64_t pos);
  void removeSeedFromBlock(int32_t blockId, int64_t pos);
  const std::unordered_set<int64_t> &getBlockSeeds(int32_t blockId) const;
  const std::vector<std::optional<seeding::seed_t>> &getAllSeeds() const;

  // Updated gap map methods
  void updateGapMap(std::string_view nodeId, int64_t pos, int64_t length,
                   bool isRemoval);
  void updateGapMapFromNodeUpdates(const std::string &nodeId);
  void updateGapMapFromInversion(const coordinates::CoordRange &blockRange,
                                 bool inverted);
  bool isGapPosition(std::shared_ptr<HierarchicalGapMap> nodeGapMap, int64_t pos) const;

  // Transaction-based gap map operations
  std::shared_ptr<UpdateBatch<coordinates::GapUpdate>>
  beginGapTransaction(const std::string &nodeId);
  bool commitGapTransaction(
      std::shared_ptr<UpdateBatch<coordinates::GapUpdate>> transaction);
  bool rollbackGapTransaction(
      std::shared_ptr<UpdateBatch<coordinates::GapUpdate>> transaction);

  // Helper method to find non-gap positions efficiently
  int64_t findNonGapPosition(std::string_view nodeId, int64_t startPos,
                             int count, bool forward) const;

  // Debug info
  void printGapMapInfo(std::string_view nodeId) const;
  void printGapMapInfo() const;
  void printGapCacheStats() const;
  void printNodeMutations(const std::string &nodeId) const;

  // 1. Lazy Sequence View - returns a view instead of materializing the full
  // sequence
  LazySequenceView getLazySequenceView(std::string_view nodeId,
                                       const coordinates::CoordRange &range,
                                       bool skipGaps = true) const;

  // 4. Region-Based Mutation Application
  void applyMutationsForRange(const std::string &nodeId,
                              const coordinates::CoordRange &range,
                              const std::vector<Mutation> &mutations);

  // 5. On-Demand Block Activation
  bool ensureBlockLoaded(const std::string &nodeId, int32_t blockId);

  // 6. Streaming K-mer Generation
  void streamKmers(
      const std::string &nodeId, const coordinates::CoordRange &range, int k,
      std::function<void(size_t hash, bool isReverse, int64_t pos)> processor);

  // 7. Range Batching
  std::vector<coordinates::CoordRange>
  optimizeBatches(const std::vector<coordinates::CoordRange> &ranges,
                  size_t optimalBatchSize = 4096);

  // 7. Shared Immutable Block Data
  std::shared_ptr<SharedBlockData> getSharedBlockData(int32_t blockId);
  void createSharedBlockData(int32_t blockId, const std::string &nodeId);

  // 8. Delta-Based Node Representation
  void propagateWithDelta(const std::string &fromNode,
                          const std::string &toNode);

  void clearAllSeeds(); // New method for placement

  // Optimized mapping from global position to block coordinates
  // This method handles block inversions and returns the inverted flag
  // Use this instead of fastMapGlobalToLocal when you need complete coordinate
  // mapping
  coordinates::BlockCoordinate
  mapGlobalToBlockCoords(std::string_view nodeId, int64_t globalPos) const;

  // Optimized method to find next/previous position accounting for block
  // boundaries
  int64_t movePosition(std::string_view nodeId, int64_t globalPos,
                       bool forward);

  // Expand ranges to include k non-gap characters in each direction
  std::vector<coordinates::CoordRange>
  expandRecompRanges(std::string_view nodeId,
                     const std::vector<coordinates::CoordRange> &ranges, int k);

  // Reset node cache
  void resetNodeCache(std::string_view nodeId);

  // Clear caches
  void clearCaches();

  // Helper method to clear caches for a specific node
  void clearCacheForNode(const std::string& nodeId);

  // Helper method to track block usage
  void trackBlockUsage(int32_t blockId, const std::string &nodeId);

  // Helper method to merge ranges
  std::vector<coordinates::CoordRange>
  mergeRanges(const std::vector<coordinates::CoordRange> &ranges) const;

  // Helper for recomp range calculation
  std::optional<coordinates::CoordRange>
  calculateRecompRange(std::string_view nodeId, int32_t blockId, int32_t pos,
                     int32_t len, bool isBlockMutation = false,
                     bool isBlockDeactivation = false);

  // Helper to merge ranges
  void mergeRangeWithExisting(std::vector<coordinates::CoordRange> &ranges,
                              const coordinates::CoordRange &newRange);

  // Helper to check if a character is a non-gap character
  bool isNonGapChar(char c) const;
  
  // Helper to consolidate gap map updates across nodes
  void consolidatedGapUpdate(const std::string& nodeId,
                            const std::vector<coordinates::GapUpdate>& updates);

  // Helper to initialize gap lists from the PanMAN tree structure
  void initializeGapLists(panmanUtils::Tree* tree);

  // Helper to set gap list length - update both legacy hash table and new array
  void setGapListLength(int32_t blockId, int32_t nucPos, size_t length);

  // Helper for batch processing gap updates
  void batchProcessGapUpdates(const std::string &nodeId);

  // Batch gap map operations
  void batchUpdateGapMap(const std::vector<coordinates::GapUpdate>& updates,
                        const std::string& nodeId);

  // Initialize the mutation index from existing character data
  void initializeMutationIndex();

  // Helper to get node gap map
  std::shared_ptr<HierarchicalGapMap>
  getNodeGapMap(std::string_view nodeId) const;

  // Update mutation index directly for a specific position and node
  void updateMutationIndex(const PositionKey& key, const std::string& nodeId);

  // Set sequence for a specific block
  void setBlockSequence(int32_t blockId, const std::string &sequence);
  
  // Access block sequences for debugging
  const std::unordered_map<int32_t, std::string>& getBlockSequences() const {
    return blockSequences;
  }

  // Method to map 3D coordinates to global position
  int64_t mapToGlobalCoordinate(int32_t blockId, int32_t nucPos, int32_t gapPos, bool inverted = false) const;

  // Method to explicitly register a mapping from block coordinates to global position
  void registerCoordinateMapping(int32_t blockId, int32_t nucPos, int32_t gapPos, int64_t globalPos);

  // Method to initialize all coordinate mappings for a block
  void initializeBlockMappings(int32_t blockId);

  // Method to flush any pending coordinate mappings
  void flushCoordinateMappings();

  // Debug functions for coordinate mapping
  void dumpBlockCoordinateMappings(const std::string& nodeId) const;
  void logGapListsForBlock(const std::string& nodeId, int32_t blockId) const;

  // Print comprehensive mapping information for all coordinates
  void dumpAllCoordinateMappings() const;
  
  // Extract and dump actual sequence data and coordinates for critical blocks
  void dumpExtractedSequences() const;

  // Initialize gap list length array
  void initializeGapListLengthArray(size_t numBlocks, size_t maxNucPos);
};

// 1. Lazy Sequence View - avoids materializing entire sequences
class LazySequenceView {
private:
  std::string_view nodeId;
  coordinates::CoordRange range;
  const StateManager &stateManager;
  bool skipGaps;

  // Cached position mapping for faster access
  mutable std::vector<std::tuple<int32_t, int32_t, int32_t, bool>>
      positionMapping;
  mutable bool positionsMapped = false;

public:
  LazySequenceView(std::string_view nodeId,
                   const coordinates::CoordRange &range,
                   const StateManager &stateManager, bool skipGaps = true);

  // Size is unknown until traversed
  size_t size() const;

  // Get character at index (relative to view start)
  char operator[](size_t index) const;

  // Get global position for index
  int64_t globalPosition(size_t index) const;

  // Get all position mappings - builds cache
  const std::vector<std::tuple<int32_t, int32_t, int32_t, bool>> &
  getPositionMappings() const;

  // Create a materialized string from this view
  std::string toString() const;

  // Get a pair of sequence string and positions (compatible with current API)
  std::pair<std::string, std::vector<int64_t>> toSequenceAndPositions() const;

  // Iterator support
  class iterator {
  private:
    const LazySequenceView *view;
    size_t index;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = char;
    using difference_type = std::ptrdiff_t;
    using pointer = const char *;
    using reference = const char &;

    iterator(const LazySequenceView *view, size_t index)
        : view(view), index(index) {}

    char operator*() const { return (*view)[index]; }
    iterator &operator++() {
      ++index;
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++index;
      return tmp;
    }
    bool operator==(const iterator &other) const {
      return index == other.index && view == other.view;
    }
    bool operator!=(const iterator &other) const { return !(*this == other); }
  };

  iterator begin() const { return iterator(this, 0); }
  iterator end() const { return iterator(this, size()); }
};

// 2. Circular Buffer for streaming k-mer generation
template <typename T> class CircularBuffer {
private:
  std::vector<T> buffer;
  size_t capacity;
  size_t head = 0;
  size_t count = 0;

public:
  CircularBuffer(size_t capacity) : capacity(capacity), buffer(capacity) {}

  void push(const T &item) {
    size_t insertPos = (head + count) % capacity;
    if (count == capacity) {
      head = (head + 1) % capacity;
    } else {
      count++;
    }
    buffer[insertPos] = item;
  }

  T &operator[](size_t index) { return buffer[(head + index) % capacity]; }

  const T &operator[](size_t index) const {
    return buffer[(head + index) % capacity];
  }

  bool isFull() const { return count == capacity; }

  bool isValid() const {
    // For genomic data, check if any position is a gap
    for (size_t i = 0; i < count; ++i) {
      if (buffer[(head + i) % capacity] == '-') {
        return false;
      }
    }
    return isFull();
  }

  void clear() {
    head = 0;
    count = 0;
  }

  size_t size() const { return count; }
};

// 3. Delta-based Node Representation
struct NodeDelta {
  // Parent node ID
  std::string parentNodeId;

  // Block mutations (blockId, isOn, isInverted)
  std::vector<std::tuple<int32_t, bool, bool>> blockMutations;

  // Character mutations
  std::unordered_map<PositionKey, char, PositionKey::Hash> characterMutations;

  // Gap list mutations (blockId, nucPos, length)
  std::vector<std::tuple<int32_t, int32_t, size_t>> gapListMutations;

  // Recomputation ranges
  std::vector<coordinates::CoordRange> recompRanges;

  // Methods to apply delta to a state
  void applyToNodeState(NodeState &state,
                        const StateManager &stateManager) const;

  // Create delta from difference between two node states
  static NodeDelta createDelta(const std::string &parentNodeId,
                               const NodeState &parentState,
                               const NodeState &childState);
};

// Mutation type for region-based mutation application
struct Mutation {
  enum class Type { Block, Nucleotide };

  Type type;
  int32_t blockId;
  int32_t nucPos;
  int32_t gapPos;
  char newChar;     // For nucleotide mutations
  bool isInsertion; // For block mutations
  bool isInversion; // For block mutations

  // Constructors for different mutation types
  static Mutation createBlockMutation(int32_t blockId, bool isInsertion,
                                      bool isInversion) {
    Mutation m;
    m.type = Type::Block;
    m.blockId = blockId;
    m.isInsertion = isInsertion;
    m.isInversion = isInversion;
    return m;
  }

  static Mutation createNucleotideMutation(int32_t blockId, int32_t nucPos,
                                           int32_t gapPos, char newChar) {
    Mutation m;
    m.type = Type::Nucleotide;
    m.blockId = blockId;
    m.nucPos = nucPos;
    m.gapPos = gapPos;
    m.newChar = newChar;
    return m;
  }
};

// Helper function declarations for tree processing utilities
template <typename T, typename ETS>
std::vector<T> mergeThreadLocalVectors(ETS& threadLocalVectors);

/**
 * Compute paths from root to each node in the tree
 * 
 * @param tree Pointer to the tree structure
 * @param rootNode Pointer to the root node
 * @return A map of node identifiers to their paths from root
 */
std::unordered_map<std::string, std::vector<panmanUtils::Node*>> 
computeNodePaths(panmanUtils::Tree* tree, panmanUtils::Node* rootNode);

/**
 * Group nodes by their level in the tree
 * 
 * @param tree Pointer to the tree structure
 * @param rootNode Pointer to the root node
 * @return A vector of node vectors, where each inner vector contains nodes at the same level
 */
std::vector<std::vector<panmanUtils::Node*>> 
groupNodesByLevel(panmanUtils::Tree* tree, panmanUtils::Node* rootNode);

} // namespace state
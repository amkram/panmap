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
#include <atomic>
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <boost/functional/hash.hpp>

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

// Define position_key_t 
using position_key_t = std::tuple<int, int, int>;

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
   * @brief Optimized hash function for PositionKey to enable use in unordered containers
   */
  struct Hash {
    std::size_t operator()(const PositionKey &k) const {
      // Optimized hash function using bit manipulation for better performance
      // This is faster than boost::hash_combine and should have good distribution
      std::size_t h = static_cast<std::size_t>(k.blockId);
      h = (h << 16) ^ static_cast<std::size_t>(k.nucPos);
      h = (h << 8) ^ static_cast<std::size_t>(k.gapPos);
      
      // Apply additional mixing for better distribution
      h ^= h >> 16;
      h *= 0x85ebca6b;
      h ^= h >> 13;
      h *= 0xc2b2ae35;
      h ^= h >> 16;
      
      return h;
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
  
  /** @brief Indicates if each position contains a gap character */
  std::vector<bool> gaps;
  
  /** @brief Pre-calculated end positions for k-mer seeds starting at each position */
  std::vector<int64_t> endPositions;

  CharacterBuffer() {
    // Preallocate with reasonable size
    buffer.reserve(4096);
    positions.reserve(4096);
    gaps.reserve(4096);
    endPositions.reserve(4096);
  }

  void clear() {
    buffer.clear();
    positions.clear();
    gaps.clear();
    endPositions.clear();
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

  /** @brief Shared mutex for thread-safe access to this node's gap map */
  mutable std::shared_mutex mapMutex;
  
  /** @brief Cache for materialized map to avoid recomputation */
  mutable GapMapPtr materializedMapCache;
  mutable bool materializedMapCacheValid = false;

  /** @brief Cache statistics for performance analysis */
  mutable std::atomic<uint64_t> materializedCacheHits{0};
  mutable std::atomic<uint64_t> materializedCacheMisses{0};

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
    // First, try to get cached result if available
    {
      std::shared_lock<std::shared_mutex> lock(mapMutex);
      if (materializedMapCacheValid && materializedMapCache) {
        return materializedMapCache;
      }
    }
    
    // Start with our local map (make copy outside lock to reduce contention)
    auto result = std::make_shared<gap_map::GapMap>();
    GapMapPtr parentMap = nullptr;
    
    // Get parent map WITHOUT holding our own lock to prevent deadlock
    if (auto parentPtr = parent.lock()) {
      try {
        // Recursive call with no locks held - safe for hierarchy traversal
        parentMap = parentPtr->materializeMap();
      } catch (const std::exception& e) {
        logging::warn("Failed to materialize parent gap map: {}", e.what());
        // Continue with local map only
      }
    }
    
    // Now acquire our lock once to do all the work atomically
    std::unique_lock<std::shared_mutex> lock(mapMutex);
    
    // Double-check cache after acquiring lock (another thread might have updated it)
    if (materializedMapCacheValid && materializedMapCache) {
      return materializedMapCache;
    }
    
    // Copy our local map
    *result = *localMap;

    // Merge parent entries that don't conflict with our local entries
    if (parentMap) {
      for (const auto &entry : *parentMap) {
        if (result->find(entry.first) == result->end()) {
          (*result)[entry.first] = entry.second;
        }
      }
    }
    
    // Cache the result atomically
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
   * @brief Check if a position contains a gap (optimized cache-first version)
   * 
   * Checks if the specified position falls within any gap range.
   * Optimized to use materialized cache when available for maximum performance.
   * 
   * @param pos Position to check
   * @return true if the position is a gap, false otherwise
   */
  bool isGap(int64_t pos) const {
    if (pos < 0) {
      return false;
    }
    
    // Try to use materialized cache first for optimal performance
    {
      std::shared_lock<std::shared_mutex> lock(mapMutex);
      if (materializedMapCacheValid && materializedMapCache) {
        materializedCacheHits.fetch_add(1, std::memory_order_relaxed);
        auto it = materializedMapCache->upper_bound(pos);
        if (it != materializedMapCache->begin()) {
          it = std::prev(it);
          int64_t gapStart = it->first;
          int64_t gapLength = it->second;
          
          if (gapLength > 0) {
            int64_t gapEnd = gapStart + gapLength - 1;
            if (pos >= gapStart && pos <= gapEnd) {
              return true;
            }
          }
        }
        return false;
      } else {
        materializedCacheMisses.fetch_add(1, std::memory_order_relaxed);
      }
    }
    
    // Fallback: use hierarchy traversal if cache not available
    const HierarchicalGapMap* current = this;
    
    while (current) {
      // Check this level's local map with shared lock for reads
      {
        std::shared_lock<std::shared_mutex> lock(current->mapMutex);
        
        auto it = current->localMap->upper_bound(pos);
        if (it != current->localMap->begin()) {
          it = std::prev(it);
          int64_t gapStart = it->first;
          int64_t gapLength = it->second;
          
          if (gapLength > 0) {
            int64_t gapEnd = gapStart + gapLength - 1;
            if (pos >= gapStart && pos <= gapEnd) {
              return true;
            }
          }
        }
      }
      
      // Move to parent - only one shared_ptr lock per level
      auto parentPtr = current->parent.lock();
      current = parentPtr.get();
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
    
    std::unique_lock<std::shared_mutex> lock(mapMutex);

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

    std::unique_lock<std::shared_mutex> lock(mapMutex);

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
    std::shared_lock<std::shared_mutex> lock(mapMutex);
    return *localMap;
  }

  // Get backtracking information for this level only
  std::vector<gap_map::GapUpdate> getLocalBacktrack() const {
    std::shared_lock<std::shared_mutex> lock(mapMutex);
    return localBacktrack;
  }

  // Check if this map has local modifications
  bool isModified() const {
    std::shared_lock<std::shared_mutex> lock(mapMutex);
    return modified;
  }

  // Clear local modifications
  void clearModifications() {
    std::unique_lock<std::shared_mutex> lock(mapMutex);
    localMap->clear();
    localBacktrack.clear();
    modified = false;
    materializedMapCacheValid = false;
  }

  /**
   * @brief Eagerly materialize the cache for performance optimization
   * 
   * This method forces materialization of the cache, which can improve
   * performance for frequently accessed gap maps by avoiding repeated
   * hierarchy traversals.
   */
  void ensureMaterializedCache() const {
    std::shared_lock<std::shared_mutex> lock(mapMutex);
    if (!materializedMapCacheValid || !materializedMapCache) {
      lock.unlock(); // Release shared lock before acquiring unique lock
      materializeMap(); // This will populate the cache
    }
  }

private:
  /** @brief Thread-local cache for gap position queries */
  struct GapQueryCache {
    int64_t lastPos = -1;
    bool lastResult = false;
    const HierarchicalGapMap* lastMap = nullptr;
    
    void clear() {
      lastPos = -1;
      lastResult = false;
      lastMap = nullptr;
    }
  };
  static thread_local GapQueryCache tlsGapCache;

  /**
   * @struct CacheStatistics
   * @brief Statistics for cache performance monitoring
   */
  struct CacheStatistics {
    std::atomic<uint64_t> totalQueries{0};
    std::atomic<uint64_t> materializedCacheHits{0};
    std::atomic<uint64_t> materializedCacheMisses{0};
    std::atomic<uint64_t> hierarchyTraversals{0};
    std::atomic<uint64_t> totalHierarchyLevels{0};
    std::atomic<uint64_t> threadLocalCacheHits{0};
    std::atomic<uint64_t> threadLocalCacheMisses{0};
    
    // Get cache hit rate percentage
    double getMaterializedHitRate() const {
      uint64_t total = materializedCacheHits + materializedCacheMisses;
      return total > 0 ? (100.0 * materializedCacheHits) / total : 0.0;
    }
    
    double getThreadLocalHitRate() const {
      uint64_t total = threadLocalCacheHits + threadLocalCacheMisses;
      return total > 0 ? (100.0 * threadLocalCacheHits) / total : 0.0;
    }
    
    double getAverageHierarchyDepth() const {
      return hierarchyTraversals > 0 ? (double)totalHierarchyLevels / hierarchyTraversals : 0.0;
    }
    
    void reset() {
      totalQueries = 0;
      materializedCacheHits = 0;
      materializedCacheMisses = 0;
      hierarchyTraversals = 0;
      totalHierarchyLevels = 0;
      threadLocalCacheHits = 0;
      threadLocalCacheMisses = 0;
    }
  };

  /** @brief Global cache statistics (static for all instances) */
  static CacheStatistics cacheStats;

  /** @brief Enable/disable cache statistics (can be toggled for performance) */
  static std::atomic<bool> enableCacheStats;
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
  absl::flat_hash_map<int64_t, T> localValues;
  
  /** @brief Mutex for thread-safe access */
  mutable std::mutex storeMutex;
  
  /** @brief Flag indicating if this store has been modified */
  bool modified = false;

  /** @brief Cache statistics for performance analysis */
  mutable std::atomic<uint64_t> localCacheHits{0};
  mutable std::atomic<uint64_t> parentTraversals{0};

public:
  /**
   * @brief Constructor for root store
   */
  HierarchicalStore() {
  }

  /**
   * @brief Constructor for child stores
   * @param parentStore Parent store to inherit from
   */
  HierarchicalStore(std::shared_ptr<HierarchicalStore<T>> parentStore)
      : parent(parentStore) {
  }

  /**
   * @brief Get a value from the store or its ancestors (optimized iterative version)
   * @param key The key to look up
   * @return The value if found, or nullopt if not found
   */
  std::optional<T> get(int64_t key) const {
    // Start with current node
    const HierarchicalStore<T>* current = this;
    bool firstLevel = true;
    
    while (current) {
      // Check local values at current level with lock scoped to this level only
      {
        std::lock_guard<std::mutex> lock(current->storeMutex);
        
        // Look for key in local values
        auto it = current->localValues.find(key);
        if (it != current->localValues.end()) {
          if (firstLevel) {
            localCacheHits.fetch_add(1, std::memory_order_relaxed);
          } else {
            parentTraversals.fetch_add(1, std::memory_order_relaxed);
          }
          return it->second;
        }
        
        // If not found locally, get parent pointer before releasing lock
        if (auto parentShared = current->parent.lock()) {
          // Move up the hierarchy
          current = parentShared.get();
          firstLevel = false;
          continue;
        }
      }
      
      // No parent or couldn't lock parent, end traversal
      break;
    }
    
    // Not found anywhere
    return std::nullopt;
  }

  /**
   * @brief Get a value from the local store only.
   * @param key The key to look up.
   * @return The value if found, or std::nullopt if not found in the local store.
   */
  std::optional<char> getLocal(const PositionKey& key) const {
    std::lock_guard<std::mutex> lock(storeMutex);
    auto it = localValues.find(key);
    if (it != localValues.end()) {
      return it->second;
    }
    return std::nullopt;
  }

  /**
   * @brief Get a value from the local store only.
   * @param key The key to look up.
   * @return The value if found, or std::nullopt if not found in the local store.
   */
  std::optional<T> getLocal(int64_t key) const {
    std::lock_guard<std::mutex> lock(storeMutex);
    auto it = localValues.find(key);
    if (it != localValues.end()) {
      return it->second;
    }
    return std::nullopt;
  }
  
  /**
   * @brief Set a value in the store
   * @param key The key to set
   * @param value The value to set
   */
  void set(int64_t key, const T& value) {
    std::lock_guard<std::mutex> lock(storeMutex);
    localValues[key] = value;
    modified = true;
  }

  /**
   * @brief Remove a value from the store
   * @param key The key to remove
   */
  void remove(int64_t key) {
    std::lock_guard<std::mutex> lock(storeMutex);
    
    // Add an entry but with a nullopt value
    localValues.erase(key);
    modified = true;
  }

  /**
   * @brief Check if the store has been modified
   * @return True if the store has been modified
   */
  bool isModified() const {
    std::lock_guard<std::mutex> lock(storeMutex);
    return modified;
  }

  /**
   * @brief Get all local values
   * @return A copy of the local values map
   */
  absl::flat_hash_map<int64_t, T> getLocalValues() const {
    std::lock_guard<std::mutex> lock(storeMutex);
    return localValues;
  }
};

// Define HierarchicalCharacterStore as an alias for HierarchicalStore<char>
// This will fix type compatibility issues between the two types
using HierarchicalCharacterStore = HierarchicalStore<char>;

// Define type alias for hierarchical seed store
using HierarchicalSeedStore = HierarchicalStore<seeding::seed_t>;

/**
 * @struct NodePositionKey
 * @brief Combines Node ID and PositionKey for caching
 */
struct NodePositionKey {
    std::string nodeId;
    PositionKey posKey;

    bool operator==(const NodePositionKey& other) const {
        return nodeId == other.nodeId && posKey == other.posKey;
    }

    struct Hash {
        std::size_t operator()(const NodePositionKey& k) const {
            std::size_t h1 = std::hash<std::string>{}(k.nodeId);
            std::size_t h2 = PositionKey::Hash{}(k.posKey);
            // Combine hashes - simple XOR for now, consider better mixing if collisions occur
            return h1 ^ (h2 << 1);
        }
    };
};

/**
 * @class HierarchicalCharacterStore
 * @brief Specialized version of HierarchicalStore for characters with PositionKey
 * 
 * This specialized version allows using PositionKey instead of int64_t as keys
 */
template<>
class HierarchicalStore<char> {
private:
  /** @brief Parent store (weak_ptr to avoid circular references) */
  std::weak_ptr<HierarchicalStore<char>> parent;
  
  /** @brief Local values that override parent values - USE absl::flat_hash_map */
  absl::flat_hash_map<PositionKey, char, PositionKey::Hash> localValues;
  
  /** @brief Mutex for thread-safe access */
  mutable std::mutex storeMutex;
  
  /** @brief Flag indicating if this store has been modified */
  bool modified = false;

  /** @brief Cache statistics for performance analysis */
  mutable std::atomic<uint64_t> localCacheHits{0};
  mutable std::atomic<uint64_t> parentTraversals{0};
  mutable std::atomic<uint64_t> tlsCacheHits{0};
  mutable std::atomic<uint64_t> tlsCacheMisses{0};

  // Define the cache type for this specialization
  using CacheType = absl::flat_hash_map<NodePositionKey, char, NodePositionKey::Hash>;
  // Declare the thread-local cache (definition in .cpp file)
  static thread_local CacheType tlsCache;

public:
  /**
   * @brief Constructor for root store
   */
  HierarchicalStore() {}

  /**
   * @brief Constructor for child stores
   * @param parentStore Parent store to inherit from
   */
  HierarchicalStore(std::shared_ptr<HierarchicalStore<char>> parentStore)
      : parent(parentStore) {
    // Remove this warning - root node intentionally has null parent
    // if (!parentStore) {
    //   logging::warn("HierarchicalStore created with null parent");
    // }
  }

  /**
   * @brief Get value, checking parent if not in local store
   * 
   * @param key The key to look up
   * @return The value if found, std::nullopt otherwise
   */
  std::optional<char> get(const PositionKey& key) const {
    // Optimized iterative traversal to minimize shared_ptr overhead
    const HierarchicalStore<char>* current = this;
    bool firstLevel = true;
    
    while (current) {
        // Check this level with minimal lock time
        {
            std::lock_guard<std::mutex> lock(current->storeMutex);
            auto it = current->localValues.find(key);
            if (it != current->localValues.end()) {
                if (firstLevel) {
                    localCacheHits.fetch_add(1, std::memory_order_relaxed);
                } else {
                    parentTraversals.fetch_add(1, std::memory_order_relaxed);
                }
                return it->second;
            }
        }
        
        // Move to parent - only one shared_ptr lock per level
        auto parentPtr = current->parent.lock();
        current = parentPtr.get();
        firstLevel = false;
    }
    
    return std::nullopt;
  }

  /**
   * @brief Get a value from the local store only.
   * @param key The key to look up.
   * @return The value if found, or std::nullopt if not found in the local store.
   */
  std::optional<char> getLocal(const PositionKey& key) const {
    std::lock_guard<std::mutex> lock(storeMutex);
    auto it = localValues.find(key);
    if (it != localValues.end()) {
      return it->second;
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
  bool set(const PositionKey& key, const char& value) {
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
  absl::flat_hash_map<PositionKey, char, PositionKey::Hash> getLocalValues() const {
    std::lock_guard<std::mutex> lock(storeMutex);
    return localValues;
  }

  /**
   * @brief Clears the thread-local cache for the current thread.
   */
  static void clearThreadLocalCache() {
    tlsCache.clear();
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

  absl::flat_hash_set<int32_t> activeBlocks;
  // std::unordered_map<int32_t, bool> blockStrands; // OLD: true = forward, false = inverted
  // Use optional<bool> to represent orientation: nullopt=OFF, true=ON_FORWARD, false=ON_INVERTED
  absl::flat_hash_map<int32_t, std::optional<bool>> blockOrientation;

  std::vector<std::tuple<int32_t, int32_t, int32_t, char, char>>
      nucleotideChanges;
  std::vector<coordinates::CoordRange> recompRanges;
  
  // Track boundaries of active blocks for efficient range expansion
  std::optional<int64_t> firstActiveBlockPos;  // First position of any active block
  std::optional<int64_t> lastActiveBlockPos;   // Last position of any active block
  
  // Recorded mutations for debugging
  std::vector<std::pair<int64_t, uint8_t>> recordedMutations;

  // Fields for gap map updates
  std::vector<coordinates::GapUpdate> gapMapUpdates;

  // Node-specific hierarchical gap map
  std::shared_ptr<HierarchicalGapMap> gapMap;

  // Node-specific hierarchical character store
  std::shared_ptr<HierarchicalCharacterStore> characterStore;

  // Compressed gap runs for this node
  std::vector<gap_map::CompressedGapRun> compressedGapRuns;

  // Quaternary-encoded seed changes (for efficient storage and transmission)
  std::vector<int64_t> seedChangeBasePositions;  // Base positions for seed changes
  std::vector<uint64_t> seedChangeBitMasks;      // Quaternary-encoded masks
  
  // Efficient seed count tracking to avoid recomputation
  int64_t totalSeedCount = 0;        // Total seeds at this node
  int64_t seedDeletions = 0;         // Seeds deleted from parent
  int64_t seedInsertions = 0;        // Seeds added at this node
  int64_t seedModifications = 0;     // Seeds modified from parent
  
  // NEW: Track inheritance vs local changes
  int64_t inheritedSeedCount = 0;    // Seeds inherited from parent
  int64_t localSeedChanges = 0;      // Net change from local mutations (+/-)
  
  // Helper methods for seed tracking
  void initializeSeedCountFromParent(int64_t parentCount) {
    inheritedSeedCount = parentCount;
    localSeedChanges = 0;
    totalSeedCount = inheritedSeedCount + localSeedChanges;
  }
  
  void updateSeedCounts(int64_t cleared, int64_t added, int64_t changed) {
    seedDeletions += cleared;
    seedInsertions += added;
    seedModifications += changed;
    localSeedChanges += (added - cleared);
    totalSeedCount = inheritedSeedCount + localSeedChanges;
  }
  
  int64_t getTotalSeedCount() const {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    if (materializedStateComputed && !materializedSeeds.empty()) {
      return materializedSeeds.size();
    }
    return inheritedSeedCount + localSeedChanges;
  }
  
  // Materialized inheritance storage for performance optimization
  // These store the complete state (parent + local changes) for O(1) lookups
  absl::flat_hash_map<PositionKey, char, PositionKey::Hash> materializedCharacters;
  absl::flat_hash_map<int64_t, seeding::seed_t> materializedSeeds;
  
  // Flag to indicate if materialized state has been computed
  bool materializedStateComputed = false;
  
  // Mutex for thread-safe access to materialized state
  mutable std::mutex materializedStateMutex;
  
  /**
   * @brief Add seed changes in quaternary-encoded format
   * 
   * @param basePositions Base positions for seed changes
   * @param bitMasks Quaternary-encoded masks
   */
  void addSeedChanges(const std::vector<int64_t>& basePositions, 
                     const std::vector<uint64_t>& bitMasks) {
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
  

  
  /**
   * @brief Materialize complete character state by inheriting from parent + applying local changes
   * This method should be called during node initialization for optimal performance.
   * 
   * @param parentState Parent node state to inherit from (nullptr for root)
   * @param characterStore Local character store with node-specific mutations
   * @param rootCharacters Fallback to root/reference sequences for unmutated positions
   */
  void materializeCharacterState(const NodeState* parentState, 
                                 const HierarchicalCharacterStore* characterStore,
                                 const absl::flat_hash_map<PositionKey, char, PositionKey::Hash>& rootCharacters) {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    
    // Clear any existing materialized state
    materializedCharacters.clear();
    
    // Step 1: Inherit from parent state
    if (parentState && parentState->materializedStateComputed) {
      std::lock_guard<std::mutex> parentLock(parentState->materializedStateMutex);
      materializedCharacters = parentState->materializedCharacters;
    }
    
    // Step 2: Apply local mutations
    if (characterStore) {
      auto localValues = characterStore->getLocalValues();
      for (const auto& [posKey, character] : localValues) {
        materializedCharacters[posKey] = character;
      }
    }
    
    materializedStateComputed = true;
  }
  
  /**
   * @brief Materialize complete seed state by inheriting from parent + applying local changes
   * 
   * @param parentState Parent node state to inherit from (nullptr for root)
   * @param seedStore Local seed store with node-specific mutations
   */
  void materializeSeedState(const NodeState* parentState,
                           const HierarchicalSeedStore* seedStore) {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    
    materializedSeeds.clear();
    
    // Step 1: Inherit from parent seed state
    if (parentState && parentState->materializedStateComputed) {
      std::lock_guard<std::mutex> parentLock(parentState->materializedStateMutex);
      materializedSeeds = parentState->materializedSeeds;
    }
    
    // Step 2: Apply local seed mutations
    if (seedStore) {
      auto localValues = seedStore->getLocalValues();
      for (const auto& [position, seed] : localValues) {
        materializedSeeds[position] = seed;
      }
    }
    
    materializedStateComputed = true;
  }
  
  /**
   * @brief Fast O(1) lookup for characters using materialized state
   * 
   * @param posKey Position key to look up
   * @return Character if found, nullopt otherwise
   */
  std::optional<char> getMaterializedCharacter(const PositionKey& posKey) const {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    
    if (!materializedStateComputed) {
      return std::nullopt;  // Materialized state not ready
    }
    
    auto it = materializedCharacters.find(posKey);
    if (it != materializedCharacters.end()) {
      return it->second;
    }
    
    return std::nullopt;
  }
  
  /**
   * @brief Fast O(1) lookup for seeds using materialized state
   * 
   * @param position Position to look up
   * @return Seed if found, nullopt otherwise
   */
  std::optional<seeding::seed_t> getMaterializedSeed(int64_t position) const {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    
    if (!materializedStateComputed) {
      return std::nullopt;  // Materialized state not ready
    }
    
    auto it = materializedSeeds.find(position);
    if (it != materializedSeeds.end()) {
      return it->second;
    }
    
    return std::nullopt;
  }
  
  /**
   * @brief Add a character mutation to materialized state
   * This should be called when mutations are applied to maintain consistency
   * 
   * @param posKey Position key
   * @param character Character value
   */
  void setMaterializedCharacter(const PositionKey& posKey, char character) {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    materializedCharacters[posKey] = character;
  }
  
  /**
   * @brief Add a seed mutation to materialized state
   * 
   * @param position Position
   * @param seed Seed value
   */
  void setMaterializedSeed(int64_t position, const seeding::seed_t& seed) {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    materializedSeeds[position] = seed;
  }
  
  /**
   * @brief Optimized bulk materialization for processing multiple child nodes
   * Allows for shared parent state to be reused efficiently across siblings
   * 
   * @param parentState Parent state to inherit from
   * @param characterStore Local character mutations
   * @param seedStore Local seed mutations  
   * @param rootCharacters Root/reference fallback
   * @param enableMoveOptimization If true, attempts to use move semantics where possible
   */
  void materializeBulkState(const NodeState* parentState,
                           const HierarchicalCharacterStore* characterStore,
                           const HierarchicalSeedStore* seedStore,
                           const absl::flat_hash_map<PositionKey, char, PositionKey::Hash>& rootCharacters,
                           bool enableMoveOptimization = true) {
    std::lock_guard<std::mutex> lock(materializedStateMutex);
    
    // Clear existing state
    materializedCharacters.clear();
    materializedSeeds.clear();
    
    // Materialize characters
    if (parentState && parentState->materializedStateComputed) {
      std::lock_guard<std::mutex> parentLock(parentState->materializedStateMutex);
      
      // Calculate sizes for optimal allocation
      size_t charParentSize = parentState->materializedCharacters.size();
      size_t charLocalSize = characterStore ? characterStore->getLocalValues().size() : 0;
      materializedCharacters.reserve(charParentSize + charLocalSize);
      
      size_t seedParentSize = parentState->materializedSeeds.size();
      size_t seedLocalSize = seedStore ? seedStore->getLocalValues().size() : 0;
      materializedSeeds.reserve(seedParentSize + seedLocalSize);
      
      // Bulk copy parent state
      materializedCharacters = parentState->materializedCharacters;
      materializedSeeds = parentState->materializedSeeds;
    }
    
    // Apply local character mutations
    if (characterStore) {
      auto localValues = characterStore->getLocalValues();
      for (const auto& [posKey, character] : localValues) {
        materializedCharacters.insert_or_assign(posKey, character);
      }
    }
    
    // Apply local seed mutations
    if (seedStore) {
      auto localValues = seedStore->getLocalValues();
      for (const auto& [position, seed] : localValues) {
        materializedSeeds.insert_or_assign(position, seed);
      }
    }
    
    materializedStateComputed = true;
  }
  
  // Helper method for LOCAL block status only (doesn't implement inheritance)
  // The StateManager::isBlockOn_unsafe method handles true inheritance by checking parent nodes
  bool isBlockOn(int32_t blockId) const {
    // First check if we have an explicit setting for this block
    if (hasExplicitBlockState(blockId)) {
      return isBlockExplicitlyOn(blockId);
    }
    
    // Fallback: When no explicit state is set, rely ONLY on activeBlocks (without inheritance)
    // This should generally not be used for inheritance-sensitive code
    return activeBlocks.contains(blockId);
  }
  
  // Method to "capture" an inherited state by setting it explicitly in this node
  void captureInheritedBlockState(int32_t blockId, bool isOn) {
    // Set the explicit state to match what was inherited using activeBlocks
    if (isOn) {
      activeBlocks.insert(blockId);
    } else {
      activeBlocks.erase(blockId);
    }
  }

  // Helper method for LOCAL block inversion status (doesn't implement inheritance)
  // The StateManager::isBlockInverted_unsafe method handles true inheritance
  bool isBlockInverted(int32_t blockId) const {
    // Keep using blockOrientation map for inversion status, but only if block is active locally
    if (!isBlockOn(blockId)) {
      return false; // Not meaningful if block isn't active
    }
    
    auto it = blockOrientation.find(blockId);
    // Block is INVERTED if it exists and the optional value is present and false
    return it != blockOrientation.end() && it->second.has_value() && !it->second.value();
  }

  void setBlockOn(int32_t blockId, bool on) {
    if (on) {
        activeBlocks.insert(blockId); 
        blockOrientation.try_emplace(blockId, true); // Default to ON-FORWARD if not already specifically set
    } else {
        activeBlocks.erase(blockId);            
        blockOrientation[blockId] = std::nullopt; // Explicitly mark as OFF at this node level
    }
  }

  void setBlockForward(int32_t blockId, bool forward) {
    // Only set orientation if the block is effectively ON at this node level.
    // An explicit entry in blockOrientation (even if nullopt for OFF) means this node has a say.
    // If not in blockOrientation, but in activeBlocks, it's ON (e.g. inherited/root).
    auto it = blockOrientation.find(blockId);
    if (it != blockOrientation.end()) {
        if (it->second.has_value()) { // It's explicitly ON (true or false value)
            it->second = forward;
        }
        // If it was std::nullopt (explicitly OFF), do nothing - orientation is irrelevant.
    } else if (activeBlocks.count(blockId)) {
        // It's in activeBlocks (so considered ON), but no specific orientation override yet.
        blockOrientation[blockId] = forward;
    }
    // If not in blockOrientation and not in activeBlocks, it's OFF from this node's perspective.
    // Do not activate it here; applyBlockMutation should use setBlockOn first.
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

  // Custom move constructor
  NodeState(NodeState &&other) noexcept
      : activeBlocks(std::move(other.activeBlocks)),
        blockOrientation(std::move(other.blockOrientation)),
        nucleotideChanges(std::move(other.nucleotideChanges)),
        recompRanges(std::move(other.recompRanges)),
        gapMapUpdates(std::move(other.gapMapUpdates)),
        gapMap(std::move(other.gapMap)),
        characterStore(std::move(other.characterStore)),
        compressedGapRuns(std::move(other.compressedGapRuns)),
        seedChangeBasePositions(std::move(other.seedChangeBasePositions)),
        seedChangeBitMasks(std::move(other.seedChangeBitMasks)) {
    // No mutexes to worry about anymore
  }

  // Custom move assignment operator
  NodeState &operator=(NodeState &&other) noexcept {
    if (this != &other) {
      activeBlocks = std::move(other.activeBlocks);
      blockOrientation = std::move(other.blockOrientation);
      nucleotideChanges = std::move(other.nucleotideChanges);
      recompRanges = std::move(other.recompRanges);
      gapMapUpdates = std::move(other.gapMapUpdates);
      gapMap = std::move(other.gapMap);
      characterStore = std::move(other.characterStore);
      compressedGapRuns = std::move(other.compressedGapRuns);
      seedChangeBasePositions = std::move(other.seedChangeBasePositions);
      seedChangeBitMasks = std::move(other.seedChangeBitMasks);
      // No mutexes to worry about anymore
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
    newState.blockOrientation = blockOrientation;
    newState.nucleotideChanges = nucleotideChanges;
    newState.recompRanges = recompRanges;
    newState.gapMapUpdates = gapMapUpdates;
    newState.gapMap = gapMap;
    newState.characterStore = characterStore;
    newState.compressedGapRuns = compressedGapRuns;
    newState.seedChangeBasePositions = seedChangeBasePositions;
    newState.seedChangeBitMasks = seedChangeBitMasks;

    return newState;
  }

  // Check if block has an explicit setting at this node
  bool hasExplicitBlockState(int32_t blockId) const {
    // Use activeBlocks as the source of truth for block state
    return activeBlocks.contains(blockId);
  }
  
  // Check if block is explicitly ON at this node
  bool isBlockExplicitlyOn(int32_t blockId) const {
    return activeBlocks.contains(blockId);
  }
  
  // Set explicit block state for inheritance-based approach
  void setExplicitBlockState(int32_t blockId, bool on) {
    // Use activeBlocks as the primary block state tracker
    if (on) {
      activeBlocks.insert(blockId);
    } else {
      activeBlocks.erase(blockId);
    }
    
    // Important: If turning off, clear any orientation data since it's irrelevant
    if (!on) {
      auto it = blockOrientation.find(blockId);
      if (it != blockOrientation.end()) {
        it->second = std::nullopt; // Mark as explicitly OFF
      }
    }
  }

  // Recomputation ranges
  std::optional<coordinates::CoordRange> calculateRecompRange(
      std::string_view nodeId, int32_t blockId, int32_t pos, int32_t len,
      bool isBlockMutation, bool isBlockDeactivation);
};

class StateManager {
private:
  // Basic State
  size_t numCoords;  // Total number of coordinates in the system
  size_t numBlocks;  // Total number of blocks
  // Thread-safe access
  mutable std::mutex nodeStatesMutex;

  // Shared mutex for thread safety in node state access
  mutable std::shared_mutex nodeMutex;
  
  // k-mer size for recomputation range calculations
  int16_t kmerSize;
  int16_t smerSize;  // s-mer size for seeding

  // Mapping from node ID to node state
  mutable absl::flat_hash_map<std::string, NodeState> nodeStates;

  // Streamlined hierarchical node structure for efficient traversal
  struct NodeHierarchy {
    // Core hierarchy relationships
    std::string parentId;
    std::vector<std::string> childrenIds;
    
    // Essential data stores for the 3D coordinate system
    std::shared_ptr<HierarchicalGapMap> gapMap;         // Critical for gap position handling
    std::shared_ptr<HierarchicalCharacterStore> characterStore; // For nucleotide data
    
    // Seed storage for optimized lookups
    std::shared_ptr<HierarchicalSeedStore> seedStore;
    int64_t dfsIndex = -1; // DFS index for ancestor/descendant checks
    int64_t subtreeSize = 0; // Size of the subtree rooted at this node
  };

  // Hierarchy management
  absl::flat_hash_map<std::string, NodeHierarchy> nodeHierarchy;

  // Reference block sequences
  absl::flat_hash_map<int32_t, std::string> blockSequences;

  // Block ranges and sequences are global properties (not node-specific)
  // All nodes reference the same blocks, they just have different active sets
  absl::flat_hash_map<int32_t, coordinates::CoordRange> blockRanges;
  
  // 2D array for fast gap list length lookup
  std::vector<std::vector<size_t>> gapListLengthArray;
  
  // Mapping from block ID to sorted list of gap positions
  absl::flat_hash_map<int32_t, std::vector<int64_t>> blockIdToGapList;
  
  // Global gap map
  gap_map::GapMap globalGapMap;
  std::vector<gap_map::GapUpdate> gapMapBacktrack;
  mutable std::shared_mutex globalGapMapMutex;

  // Flag indicating if cache is initialized
  std::atomic<bool> globalPosCacheInitialized{false};

  // Global position cache for fast mapping
  struct GlobalPosCache {
    // std::vector<int32_t> posToBlockId;      // Maps global position to block ID - REPLACED
    std::vector<int64_t> blockStartOffsets; // Starting offset of each block
    int64_t maxPosition = 0;                // Maximum global position

    // NEW: Direct mapping from global position to 3D coordinates
    std::vector<std::optional<std::tuple<int32_t, int32_t, int32_t>>> posToCoords; 

    void clear() {
      // posToBlockId.clear(); // REPLACED
      posToCoords.clear();
      blockStartOffsets.clear();
      maxPosition = 0;
    }

    bool isValid(int64_t pos) const { return pos >= 0 && pos < maxPosition; }

    // Replaced getBlockId with getCoords
    std::optional<std::tuple<int32_t, int32_t, int32_t>> getCoords(int64_t pos) const {
      if (!isValid(pos) || static_cast<size_t>(pos) >= posToCoords.size()) {
        return std::nullopt;
      }
      return posToCoords[pos];
    }

    // Optional: Keep a way to get block ID if needed elsewhere, derived from coords
    int32_t getBlockId(int64_t pos) const {
      auto coordsOpt = getCoords(pos);
      if (coordsOpt.has_value()) {
          return std::get<0>(coordsOpt.value());
      }
      return -1; 
    }

  } globalPosCache;

  // NOTE: Seeds are now managed through hierarchical seed storage per node
  // No global positionSeeds array needed
  // Legacy positionSeeds for backward compatibility (should be removed eventually)
  std::vector<std::optional<seeding::seed_t>> positionSeeds;

  // Mapping from blocks to seed positions
  absl::flat_hash_map<int32_t, absl::flat_hash_set<int64_t>> blockToSeeds;

  // Thread-local character buffers
  mutable tbb::enumerable_thread_specific<CharacterBuffer> charBuffers;

  // 7. Shared Immutable Block Data
  absl::flat_hash_map<int32_t, std::shared_ptr<SharedBlockData>>
      sharedBlockDataMap;

  // Root gap map (shared by all nodes)
  std::shared_ptr<HierarchicalGapMap> rootGapMap;

  // Transaction management
  absl::flat_hash_map<std::string,
                     std::shared_ptr<UpdateBatch<coordinates::GapUpdate>>> activeTransactions;
  std::mutex transactionMutex;

  // Gap statistics for adaptive optimization
  gap_map::GapStatistics gapStats;
  
  // Cached root node ID for faster lookups
  mutable std::string cachedRootNodeId;
  
  // Coordinate mapping between 3D coordinates and scalar global positions
  absl::flat_hash_map<PositionKey, int64_t, PositionKey::Hash> blockCoordToGlobalPos;
  // Maps blockId, nucPos, gapPos to string index in sequences
  absl::flat_hash_map<PositionKey, int32_t, PositionKey::Hash> blockCoordToStringIndex;
  
  // Block range mappings
  absl::flat_hash_map<int32_t, coordinates::CoordRange> blockRangeMappings;
  bool blockRangeMappingsInitialized = false;
  
  // Method to update block coordinate mappings using block ranges as the source of truth
  void updateBlockCoordMapping();

  // Helper method to merge overlapping ranges
  std::vector<coordinates::CoordRange>
  mergeRanges(std::vector<coordinates::CoordRange> &ranges);

  // Helper for DFS index computation
  void computeDfsIndices(const std::string &nodeId, int64_t &index);

  // Helper to check if one node is a descendant of another
  bool isDescendant(const std::string &ancestorId,
                    const std::string &nodeId) const;
  
  // Helper to get root/reference character for fallback lookup
  std::optional<char> getRootCharacter(const PositionKey& posKey) const;

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

  // Internal helper methods
  std::vector<std::pair<int32_t, coordinates::CoordRange>> getActiveBlockRangesImpl(std::string_view nodeId) const;
  std::shared_ptr<HierarchicalGapMap> getNodeGapMap(std::string_view nodeId) const;
  
  // Logging reference functions
  static bool& getInitLoggingRef();
  static bool& getPropLoggingRef();

public:

  // Constructors
  StateManager();
  StateManager(size_t numCoordinates);

  // Initialization and Setup
  void initialize(panmanUtils::Tree* T, size_t numBlocks);
  void initializeNodeHierarchy(panmanUtils::Tree* T, panmanUtils::Node* rootNode);

  // Sequence and K-mer Extraction
  std::pair<std::string, std::vector<int64_t>> extractKmer(std::string_view nodeId, int64_t pos, int k, bool reverse = false) const;
  std::tuple<std::string, std::vector<int64_t>, std::vector<bool>, std::vector<int64_t>> extractSequence(std::string_view nodeId, const coordinates::CoordRange& range, bool includeGaps = false) const;
  
  // Additional extractSequence overloads for backward compatibility
  std::string extractSequence(const std::string& nodeId, int64_t start, int64_t end, bool includeGaps, bool unused) const;
  std::string extractSequence(const std::string& nodeId, int64_t start, int64_t end, bool includeGaps);
  std::string extractSequence(const std::string& nodeId, int64_t start, int64_t end);

  // Block and Node Information
  const absl::flat_hash_set<int32_t>& getActiveBlocks(std::string_view nodeId) const;
  coordinates::CoordRange getBlockRange(int32_t blockId) const;


  // Public method to get active block ranges (non-const version)
  std::vector<std::pair<int32_t, coordinates::CoordRange>>
  getActiveBlockRanges(std::string_view nodeId);

  // Public method to get active block ranges (const version)
  std::vector<std::pair<int32_t, coordinates::CoordRange>>
  getActiveBlockRanges(std::string_view nodeId) const;


  void setBlockOn(std::string_view nodeId, int32_t blockId, bool on);
  void setBlockForward(std::string_view nodeId, int32_t blockId, bool forward);
  void setBlockInverted(std::string_view nodeId, int32_t blockId, bool inverted);
  bool isBlockOn(std::string_view nodeId, int32_t blockId) const;
  bool isBlockInverted(std::string_view nodeId, int32_t blockId) const;
  bool applyBlockMutation(std::string_view nodeId, int32_t blockId, bool isInsertion, bool isInversion_flag);
  void initializeNode(const std::string &nodeId);
  
  // Materialize node state for optimal performance
  void materializeNodeState(const std::string& nodeId);

  int16_t getKmerSize() const;
  void setKmerSize(int16_t k);

  // Coordinate mapping
  std::optional<std::tuple<int32_t, int32_t, int32_t>>
  fastMapGlobalToLocal(int64_t globalPos) const;

  // Fast lookup from (blockId, nucPos, gapPos) to global position
  int64_t fastMapLocalToGlobal(int32_t blockId, int32_t nucPos, int32_t gapPos) const;

  // Data structures for k-mer sequences and offsets
  std::vector<std::string> nodeIdsByDfsIndex; // Maps DFS indices to node IDs
  absl::flat_hash_map<std::string, absl::flat_hash_map<int64_t, std::string>> nodeKmerSequences; // Maps node IDs to position-indexed k-mer sequences
  absl::flat_hash_map<std::string, absl::flat_hash_map<int64_t, uint32_t>> nodeKmerEndOffsets; // Maps node IDs to position-indexed k-mer end offsets
  
  // Block root character flat indices for fast lookup
  absl::flat_hash_map<int32_t, absl::flat_hash_map<PositionKey, int64_t, PositionKey::Hash>> blockRootCharFlatIndices;
  
  // Character and Seed Access
  char getCharAtPosition(std::string_view nodeId, int32_t blockId, int32_t nucPos, int32_t gapPos) const;
  bool setCharAtPosition(std::string_view nodeId, int32_t blockId, int32_t nucPos, int32_t gapPos, char newChar);
  std::optional<seeding::seed_t> getSeedAtPosition(std::string_view nodeId, int64_t pos) const;
  void setSeedAtPosition(std::string_view nodeId_sv, int64_t pos, const seeding::seed_t& seed);
  void clearSeedAtPosition(std::string_view nodeId, int64_t pos);
  
  // Block State and Mutations
  bool isBlockOn_unsafe(std::string_view nodeId, int32_t blockId) const;
  bool isBlockInverted_unsafe(std::string_view nodeId, int32_t blockId) const;
  bool isGapAt(int32_t blockId, int64_t blockPos) const;
  bool isGapPosition(std::shared_ptr<HierarchicalGapMap> nodeGapMap, int64_t pos) const;
  bool isNonGapChar(char c) const;
  
  // Coordinate Mapping
  int64_t getDfsIndex(const std::string &nodeId) const;
  int64_t mapToGlobalCoordinate(int32_t blockId, int32_t nucPos, int32_t gapPos, bool inverted) const;
  int32_t getBlockIdFromPosition(int64_t pos) const;
  coordinates::BlockCoordinate mapGlobalToBlockCoords(int64_t globalPos) const;
  coordinates::BlockCoordinate mapGlobalToBlockCoords(const std::string& nodeId, int64_t globalPos) const;
  coordinates::BlockCoordinate mapGlobalToBlockCoords(std::string_view nodeId, int64_t globalPos) const;
  
  // Gap Management
  void addGap(int32_t blockId, int64_t blockPos);
  void removeGap(int32_t blockId, int64_t blockPos);
  void batchUpdateGapMap(const std::string& nodeId, const std::vector<coordinates::GapUpdate>& updates);
  void consolidatedGapUpdate(const std::string& nodeId, const std::vector<coordinates::GapUpdate>& updates);
  size_t getGapListLength(int32_t blockId, int32_t nucPos) const;
  
  // Seed Management
  void addSeedToBlock(int32_t blockId, int64_t seedPos);
  void removeSeedFromBlock(int32_t blockId, int64_t seedPos);
  const absl::flat_hash_set<int64_t>& getBlockSeeds(int32_t blockId) const;
  int getSmerSize() const;
  void setSmerSize(int s);
  
  // Node and Block Management
  void resetNodeCache(std::string_view nodeId);
  void propagateState(const std::string &fromNode, const std::string &toNode);
  void applyNucleotideMutation(const std::string& nodeId, NodeState& nodeState, int32_t blockId, int32_t nucPos, int32_t gapPos, char value);
  void trackBlockUsage(int32_t blockId, const std::string &nodeId);
  void backtrackNode(const std::string &nodeId);
  void setBlockSequence(int32_t blockId, const std::string& sequence);
  const absl::flat_hash_map<int32_t, std::string>& getBlockSequences() const;
  
  // Initialization and Setup
  void initializeDfsIndices(panmanUtils::Tree *tree);
  void initializeDfsIndices(panmanUtils::Tree *tree, panmanUtils::Node *node, int64_t &dfsIndex);
  void initializeGapListLengthArray(size_t numBlocks, size_t providedMaxNucPos);
  void initializeBlockBoundaries(const std::string& nodeId);
  void initializeBlockMappings(int32_t blockId);
  void initializeBlockRangeMappings();
  void initializeGlobalCache();
  bool isGlobalCacheInitialized() const;
  void initializeSeedStorage();
  
  // Utility Functions
  void flushCoordinateMappings();
  void registerCoordinateMapping(int32_t blockId, int32_t nucPos, int32_t gapPos, int64_t globalPos);
  void setVerboseLogging(bool verbose);
  void optimizeGapListMemory();
  std::vector<std::string> buildTopologicalOrder(panmanUtils::Tree *tree);
  std::vector<coordinates::CoordRange> expandRecompRanges(std::string_view nodeId, const std::vector<coordinates::CoordRange>& ranges, int k_param);
  std::vector<coordinates::CoordRange> getRecompRanges(const std::string& nodeId) const;
  const std::vector<coordinates::CoordRange>& getRecompRanges(std::string_view nodeId) const;
  void mergeRangeWithExisting(std::vector<coordinates::CoordRange> &ranges, const coordinates::CoordRange &newRange);
  void processSeedsInRange(const std::string& nodeId, const coordinates::CoordRange& range, bool includeGaps);
  void processSeedsInRange(std::string_view nodeId, const coordinates::CoordRange& range, int k, int s,
                          const std::function<void(int64_t, int64_t, std::optional<seeding::seed_t>&, std::string_view)>& processFn);

  // Node state access
  NodeState &getNodeState(const std::string &nodeId);
  const NodeState &getNodeState(const std::string &nodeId) const;

  // Utility
  size_t getNumCoords() const { return numCoords; }
  size_t getNumBlocks() const;
  void setNumBlocks(size_t blockCount);
  void setBlockRange(int32_t blockId, const coordinates::CoordRange &range);
  void setGapListLength(int32_t blockId, int32_t nucPos, size_t length);
  void clearCaches();

  // Recomputation range calculation
  std::optional<coordinates::CoordRange> calculateRecompRange(std::string_view nodeId, int32_t blockId,
                                                            int32_t pos, int32_t len, bool isBlockMutation,
                                                            bool isBlockDeactivation);
};

// Namespace-level utility functions (implemented in state_helpers.cpp)
std::vector<std::vector<panmanUtils::Node*>> groupNodesByLevel(panmanUtils::Tree* tree, panmanUtils::Node* rootNode);
std::unordered_map<std::string, std::vector<panmanUtils::Node*>> computeNodePaths(panmanUtils::Tree* tree, panmanUtils::Node* rootNode);

} // namespace state
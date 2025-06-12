#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

#include "panman.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <new>
#include <queue>
#include <stack>
#include <stdlib.h>
#include <unordered_set>
#include <utility>
#include <vector>

/**
 * @brief Performance utilities for threading and memory management
 */

namespace threading {

struct GroupInfo {
  std::shared_ptr<panmanUtils::Node> startNode; // Where this group starts
  std::shared_ptr<panmanUtils::Node> stopNode;  // Where this group ends
  std::unordered_set<std::shared_ptr<panmanUtils::Node>>
      groupNodes; // Hash set for O(1) membership testing
};

/**
 * @brief Computes a DFS traversal order of nodes starting from root
 * @param root The root node to start traversal from
 * @return Vector of nodes in DFS order
 */
inline std::vector<std::shared_ptr<panmanUtils::Node>>
computeDFSOrder(std::shared_ptr<panmanUtils::Node> root) {
  std::vector<std::shared_ptr<panmanUtils::Node>> dfsOrder;
  if (!root)
    return dfsOrder;

  std::stack<std::shared_ptr<panmanUtils::Node>> stack;
  stack.push(root);

  while (!stack.empty()) {
    auto current = stack.top();
    stack.pop();
    dfsOrder.push_back(current);

    // Push children in reverse order to process leftmost child first
    for (auto it = current->children.rbegin(); it != current->children.rend();
         ++it) {
      // Create shared_ptr from raw pointer to maintain ownership semantics
      stack.push(std::make_shared<panmanUtils::Node>(**it));
    }
  }

  return dfsOrder;
}

/**
 * @brief Splits a DFS ordered list of nodes into balanced groups for parallel
 * processing
 * @param dfsOrder Vector of nodes in DFS order
 * @param numGroups Number of groups to split into
 * @return Vector of GroupInfo containing group boundaries and membership
 */
inline std::vector<GroupInfo> computeBalancedGroups(
    const std::vector<std::shared_ptr<panmanUtils::Node>> &dfsOrder,
    size_t numGroups) {
  std::vector<GroupInfo> groups;
  if (dfsOrder.empty() || numGroups == 0)
    return groups;

  // Ensure we don't create more groups than nodes
  numGroups = std::min(numGroups, dfsOrder.size());

  // Calculate approximate group size
  size_t nodesPerGroup = dfsOrder.size() / numGroups;
  size_t remainder = dfsOrder.size() % numGroups;

  groups.reserve(numGroups);

  size_t startIdx = 0;
  for (size_t g = 0; g < numGroups; g++) {
    GroupInfo group;

    // Calculate end index for this group, distributing the remainder
    size_t groupSize = nodesPerGroup + (g < remainder ? 1 : 0);
    size_t endIdx = std::min(startIdx + groupSize, dfsOrder.size());

    // Set start and stop nodes
    group.startNode = dfsOrder[startIdx];
    group.stopNode = dfsOrder[endIdx - 1];

    // Add nodes to the group
    for (size_t i = startIdx; i < endIdx; i++) {
      group.groupNodes.insert(dfsOrder[i]);
    }

    groups.push_back(std::move(group));
    startIdx = endIdx;
  }

  return groups;
}

} // namespace threading

namespace memory {

/**
 * @brief Aligned memory allocation for SIMD operations
 */
template <typename T>
inline T *alignedAlloc(size_t size, size_t alignment = 64) {
  void *ptr = nullptr;
  if (posix_memalign(&ptr, alignment, size * sizeof(T)) != 0) {
    throw std::bad_alloc();
  }
  return static_cast<T *>(ptr);
}

/**
 * @brief Cache-aligned optimized k-mer buffer for SIMD processing
 */
class KmerBufferPool {
private:
  // Hardware constants
  static constexpr size_t CACHE_LINE_SIZE = 64;
  static constexpr size_t MAX_KMER_LENGTH = 1024;

  // Thread-local storage for bulk processing
  static thread_local char *nucleotide_buffer;
  static thread_local size_t *hash_buffer;
  static thread_local bool *orientation_buffer;
  static thread_local int64_t *position_buffer;
  static thread_local size_t buffer_size;
  
  // Initialization protection
  static thread_local bool is_initialized;
  static inline std::mutex initialization_mutex;

public:
  /**
   * @brief Initialize buffer pool with appropriate size
   * @param k k-mer length
   * @param batch_size Number of k-mers to process in a batch
   */
  static void initialize(size_t k, size_t batch_size = 1024) {
    // Protect against race conditions during initialization
    if (is_initialized) return;
    
    std::lock_guard<std::mutex> lock(initialization_mutex);
    if (is_initialized) return; // Double-check after acquiring lock
    
    // Align buffers to cache lines for optimal SIMD
    size_t aligned_k =
        ((k + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;

    // Free any existing buffers
    cleanup();

    // Allocate aligned memory for all buffers
    nucleotide_buffer = alignedAlloc<char>(batch_size * aligned_k);
    hash_buffer = alignedAlloc<size_t>(batch_size);
    orientation_buffer = alignedAlloc<bool>(batch_size);
    position_buffer = alignedAlloc<int64_t>(batch_size);
    buffer_size = batch_size;

    // Prefetch buffers into L1 cache
    for (size_t i = 0; i < batch_size; i += CACHE_LINE_SIZE / sizeof(size_t)) {
      __builtin_prefetch(&hash_buffer[i], 1, 3);
      __builtin_prefetch(&position_buffer[i], 1, 3);
    }
    
    is_initialized = true;
  }

  /**
   * @brief Get buffer pointers for batch processing
   */
  static void getBatchBuffers(char **nucBuffer, size_t **hashBuffer,
                              bool **orientBuffer, int64_t **posBuffer,
                              size_t k) {

    if (!nucleotide_buffer) {
      initialize(k);
    }

    *nucBuffer = nucleotide_buffer;
    *hashBuffer = hash_buffer;
    *orientBuffer = orientation_buffer;
    *posBuffer = position_buffer;
  }

  /**
   * @brief Clean up buffers
   */
  static void cleanup() {
    if (nucleotide_buffer) {
      free(nucleotide_buffer);
      nucleotide_buffer = nullptr;
    }
    if (hash_buffer) {
      free(hash_buffer);
      hash_buffer = nullptr;
    }
    if (orientation_buffer) {
      free(orientation_buffer);
      orientation_buffer = nullptr;
    }
    if (position_buffer) {
      free(position_buffer);
      position_buffer = nullptr;
    }
    is_initialized = false;
  }

  /**
   * @brief Get the current buffer size
   */
  static size_t getBufferSize() { return buffer_size; }
};

/**
 * @brief Reusable buffer for sequence data
 */
class Buffer {
public:
  std::vector<char> seqBuffer;
  std::vector<int> gapsBuffer;
  std::vector<int> coordsBuffer;
  std::vector<int> deadBlocksBuffer;

  /**
   * @brief Resizes all internal buffers
   * @param size The size to resize to
   */
  inline void resize(size_t size) {
    seqBuffer.resize(size);
    gapsBuffer.resize(size);
    coordsBuffer.resize(size);
    deadBlocksBuffer.resize(std::min(size, size_t(1024)));
  }
};

/**
 * @brief Thread-safe pool of reusable memory buffers
 */
class BufferPool {
public:
  static constexpr size_t MAX_BUFFER_SIZE = 1024 * 1024; // 1MB
  static constexpr size_t MAX_POOL_SIZE = 32;

  /**
   * @brief Acquires a buffer of requested size from the pool
   * @param requested_size Desired buffer size
   * @return Pointer to acquired buffer
   */
  inline Buffer *acquire(size_t requested_size) {
    std::lock_guard<std::mutex> lock(mutex_);

    if (!available_buffers_.empty()) {
      Buffer *buffer = available_buffers_.front();
      available_buffers_.pop();
      buffer->resize(requested_size);
      return buffer;
    }

    Buffer *new_buffer = new Buffer();
    new_buffer->resize(requested_size);
    return new_buffer;
  }

  /**
   * @brief Returns a buffer to the pool
   * @param buffer Pointer to buffer to release
   */
  inline void release(Buffer *buffer) {
    if (!buffer)
      return;

    std::lock_guard<std::mutex> lock(mutex_);
    if (available_buffers_.size() < MAX_POOL_SIZE) {
      available_buffers_.push(buffer);
    } else {
      delete buffer;
    }
  }

  inline ~BufferPool() {
    while (!available_buffers_.empty()) {
      Buffer *buffer = available_buffers_.front();
      available_buffers_.pop();
      delete buffer;
    }
  }

private:
  std::queue<Buffer *> available_buffers_;
  std::mutex mutex_;
};

} // namespace memory

#endif // PERFORMANCE_HPP
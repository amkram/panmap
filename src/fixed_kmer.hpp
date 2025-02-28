#ifndef FIXED_KMER_HPP
#define FIXED_KMER_HPP

// Forward declarations for the coordinates namespace
namespace coordinates {
class CoordinateManager;
class CoordinateTraverser;
struct tupleCoord_t;
} // namespace coordinates

#include "seed_annotated_tree.hpp"
#include <algorithm>
#include <cstdint>
#include <immintrin.h> // For SIMD instructions

// TBB memory allocator for optimized allocations
#include <tbb/scalable_allocator.h>

// Forward declarations for types from other namespaces

// Memory utility for optimized allocations
namespace memory {
// Aligned allocation using TBB for SIMD operations
template <typename T>
inline T *allocate_aligned(size_t count, size_t alignment = 64) {
  return static_cast<T *>(
      scalable_aligned_malloc(count * sizeof(T), alignment));
}

// Free aligned memory
template <typename T> inline void free_aligned(T *ptr) {
  scalable_aligned_free(ptr);
}

// RAII wrapper for aligned memory
template <typename T> class AlignedBuffer {
private:
  T *buffer;
  size_t count;

public:
  AlignedBuffer(size_t count, size_t alignment = 64)
      : buffer(allocate_aligned<T>(count, alignment)), count(count) {}

  ~AlignedBuffer() {
    if (buffer)
      free_aligned(buffer);
  }

  // No copy
  AlignedBuffer(const AlignedBuffer &) = delete;
  AlignedBuffer &operator=(const AlignedBuffer &) = delete;

  // Allow move
  AlignedBuffer(AlignedBuffer &&other) noexcept
      : buffer(other.buffer), count(other.count) {
    other.buffer = nullptr;
  }

  AlignedBuffer &operator=(AlignedBuffer &&other) noexcept {
    if (this != &other) {
      if (buffer)
        free_aligned(buffer);
      buffer = other.buffer;
      count = other.count;
      other.buffer = nullptr;
    }
    return *this;
  }

  // Access
  T *get() { return buffer; }
  const T *get() const { return buffer; }
  T &operator[](size_t index) { return buffer[index]; }
  const T &operator[](size_t index) const { return buffer[index]; }
  size_t size() const { return count; }
};
} // namespace memory

namespace fixed_kmer {
// Define a constexpr list of supported k-mer sizes
// Only keeping hardware-efficient sizes that align with byte boundaries
constexpr std::array<int, 3> KmerSizes = {8, 16, 32};

// Check if K is in our list of optimized sizes
template <int K>
constexpr bool is_optimized_size = (K == 8 || K == 16 || K == 32);

// Optimized structure for compact k-mer storage
// Uses 64-bit integers to store k-mers efficiently
template <int K> struct CompactKmer {
  static_assert(K <= 32, "K must be 32 or less for CompactKmer");

  // For K=8: use a single uint16_t (each nucleotide is 2 bits)
  // For K=16: use a single uint32_t
  // For K=32: use a uint64_t
  using storage_type = typename std::conditional<
      (K <= 8), uint16_t,
      typename std::conditional<(K <= 16), uint32_t, uint64_t>::type>::type;

  // Packed binary representation of k-mer
  storage_type packed_data;

  // Pack a character sequence into compact form
  static storage_type pack(const char *sequence) {
    // Bit-packing optimization
    storage_type result = 0;
    for (int i = 0; i < K; i++) {
      // Convert nucleotide to 2-bit representation
      uint8_t bits;
      switch (sequence[i]) {
      case 'A':
      case 'a':
        bits = 0;
        break;
      case 'C':
      case 'c':
        bits = 1;
        break;
      case 'G':
      case 'g':
        bits = 2;
        break;
      case 'T':
      case 't':
        bits = 3;
        break;
      default:
        bits = 0; // Handle unexpected characters
      }
      // Shift and add to result
      result = (result << 2) | bits;
    }
    return result;
  }

  // Specialized hashers for different k-mer sizes
  static size_t hash(storage_type packed) {
    if constexpr (K <= 8) {
      // Multiplicative hash - fast for small k-mers
      constexpr uint32_t MULTIPLIER =
          2654435761u; // Knuth's multiplicative hash
      return packed * MULTIPLIER;
    } else if constexpr (K <= 16) {
      // FNV-1a - good for medium k-mers
      constexpr uint32_t FNV_PRIME = 16777619u;
      constexpr uint32_t FNV_OFFSET = 2166136261u;
      uint32_t hash = FNV_OFFSET;
      uint32_t val = packed;
      hash ^= (val & 0xFF);
      hash *= FNV_PRIME;
      hash ^= ((val >> 8) & 0xFF);
      hash *= FNV_PRIME;
      hash ^= ((val >> 16) & 0xFF);
      hash *= FNV_PRIME;
      hash ^= ((val >> 24) & 0xFF);
      hash *= FNV_PRIME;
      return hash;
    } else {
      // xxHash-inspired mixer for larger k-mers
      uint64_t x = packed;
      x ^= x >> 33;
      x *= 0xff51afd7ed558ccdULL;
      x ^= x >> 33;
      x *= 0xc4ceb9fe1a85ec53ULL;
      x ^= x >> 33;
      return x;
    }
  }
};

// Determine the closest optimal k-mer size to the requested size
int suggestOptimalKmerSize(int requestedSize);

// Primary template declaration for getSeedAtSpecialized
template <int K>
inline bool getSeedAtSpecialized(coordinates::CoordinateTraverser &traverser,
                                 coordinates::CoordinateManager &manager,
                                 size_t &resultHash, bool &resultIsReverse,
                                 int64_t &resultEndPos, const int64_t &pos,
                                 Tree *T);

// Forward declarations for batch processing helpers
template <int K>
void processKmerBatchAVX2(coordinates::CoordinateTraverser &traverser,
                          coordinates::CoordinateManager &manager,
                          const int64_t *positions, size_t batchSize,
                          size_t *resultHashes, bool *resultIsReverse,
                          int64_t *resultEndPos, bool *validResults,
                          char *buffer);

template <int K>
void processKmerBatchSSE(coordinates::CoordinateTraverser &traverser,
                         coordinates::CoordinateManager &manager,
                         const int64_t *positions, size_t batchSize,
                         size_t *resultHashes, bool *resultIsReverse,
                         int64_t *resultEndPos, bool *validResults,
                         char *buffer);

// Specialized template for getSeedAt with k=8
template <>
inline bool getSeedAtSpecialized<8>(coordinates::CoordinateTraverser &traverser,
                                    coordinates::CoordinateManager &manager,
                                    size_t &resultHash, bool &resultIsReverse,
                                    int64_t &resultEndPos, const int64_t &pos,
                                    Tree *T) {

  // Optimized implementation for k=8
  alignas(16) char buffer[8];

  // Skip validation for performance if position is in gap
  if (manager.isGapPosition(pos)) {
    return false;
  }

  // Get starting coordinate
  const coordinates::tupleCoord_t *startCoord = manager.getTupleCoord(pos);
  if (!startCoord)
    return false;

  // Check if block exists
  if (!manager.getBlockExists()[startCoord->blockId].first) {
    return false;
  }

  // Reset traverser and collect 8 characters
  traverser.reset(startCoord);
  if (!traverser.skipToNthNonGap(8, buffer)) {
    return false;
  }

  // Get end position
  resultEndPos = traverser.getCurrent()->scalar;

// Use SIMD for forward-reverse comparison when available
#ifdef __SSE4_1__
  // Pack k-mer using the specialized function
  auto forward_packed = CompactKmer<8>::pack(buffer);

  // Create reverse complement
  alignas(16) char rev_buffer[8];
  for (int i = 0; i < 8; i++) {
    char c = buffer[7 - i];
    switch (c) {
    case 'A':
    case 'a':
      rev_buffer[i] = 'T';
      break;
    case 'C':
    case 'c':
      rev_buffer[i] = 'G';
      break;
    case 'G':
    case 'g':
      rev_buffer[i] = 'C';
      break;
    case 'T':
    case 't':
      rev_buffer[i] = 'A';
      break;
    default:
      rev_buffer[i] = c;
    }
  }

  auto reverse_packed = CompactKmer<8>::pack(rev_buffer);

  // Compare and select smaller
  if (forward_packed <= reverse_packed) {
    resultHash = CompactKmer<8>::hash(forward_packed);
    resultIsReverse = false;
  } else {
    resultHash = CompactKmer<8>::hash(reverse_packed);
    resultIsReverse = true;
  }
#else
  // Fallback to standard method if SIMD not available
  auto [fHash, rHash] = seeding::hashSeq(std::string_view(buffer, 8));
  if (fHash < rHash) {
    resultHash = fHash;
    resultIsReverse = false;
  } else {
    resultHash = rHash;
    resultIsReverse = true;
  }
#endif

  return true;
}

// Specialized template for getSeedAt with k=16
template <>
inline bool
getSeedAtSpecialized<16>(coordinates::CoordinateTraverser &traverser,
                         coordinates::CoordinateManager &manager,
                         size_t &resultHash, bool &resultIsReverse,
                         int64_t &resultEndPos, const int64_t &pos, Tree *T) {

  // Optimized implementation for k=16
  alignas(16) char buffer[16];

  // Skip validation for performance if position is in gap
  if (manager.isGapPosition(pos)) {
    return false;
  }

  // Get starting coordinate
  const coordinates::tupleCoord_t *startCoord = manager.getTupleCoord(pos);
  if (!startCoord)
    return false;

  // Check if block exists
  if (!manager.getBlockExists()[startCoord->blockId].first) {
    return false;
  }

  // Reset traverser and collect 16 characters
  traverser.reset(startCoord);
  if (!traverser.skipToNthNonGap(16, buffer)) {
    return false;
  }

  // Get end position
  resultEndPos = traverser.getCurrent()->scalar;

  // Pack k-mer using the specialized function
  auto forward_packed = CompactKmer<16>::pack(buffer);

  // Create reverse complement
  alignas(16) char rev_buffer[16];
  for (int i = 0; i < 16; i++) {
    char c = buffer[15 - i];
    switch (c) {
    case 'A':
    case 'a':
      rev_buffer[i] = 'T';
      break;
    case 'C':
    case 'c':
      rev_buffer[i] = 'G';
      break;
    case 'G':
    case 'g':
      rev_buffer[i] = 'C';
      break;
    case 'T':
    case 't':
      rev_buffer[i] = 'A';
      break;
    default:
      rev_buffer[i] = c;
    }
  }

  auto reverse_packed = CompactKmer<16>::pack(rev_buffer);

  // Compare and select smaller
  if (forward_packed <= reverse_packed) {
    resultHash = CompactKmer<16>::hash(forward_packed);
    resultIsReverse = false;
  } else {
    resultHash = CompactKmer<16>::hash(reverse_packed);
    resultIsReverse = true;
  }

  return true;
}

// Specialized template for getSeedAt with k=32
template <>
inline bool
getSeedAtSpecialized<32>(coordinates::CoordinateTraverser &traverser,
                         coordinates::CoordinateManager &manager,
                         size_t &resultHash, bool &resultIsReverse,
                         int64_t &resultEndPos, const int64_t &pos, Tree *T) {

  // Optimized implementation for k=32
  alignas(32) char buffer[32];

  // Skip validation for performance if position is in gap
  if (manager.isGapPosition(pos)) {
    return false;
  }

  // Get starting coordinate
  const coordinates::tupleCoord_t *startCoord = manager.getTupleCoord(pos);
  if (!startCoord)
    return false;

  // Check if block exists
  if (!manager.getBlockExists()[startCoord->blockId].first) {
    return false;
  }

  // Prefetch memory for better performance
  __builtin_prefetch(&buffer[0], 1, 3);

  // Reset traverser and collect 32 characters
  traverser.reset(startCoord);
  if (!traverser.skipToNthNonGap(32, buffer)) {
    return false;
  }

  // Get end position
  resultEndPos = traverser.getCurrent()->scalar;

  // Pack k-mer using the specialized function
  // For 32-mers, we use two 16-mer packing operations
  uint32_t forward_packed1 = CompactKmer<16>::pack(buffer);
  uint32_t forward_packed2 = CompactKmer<16>::pack(buffer + 16);
  uint64_t forward_packed =
      (static_cast<uint64_t>(forward_packed1) << 32) | forward_packed2;

  // Create reverse complement - optimized with unrolled loop
  alignas(32) char rev_buffer[32];

#ifdef __AVX2__
  // Use AVX2 for faster reverse complement if available
  // This is a simplified representation of what would be AVX2 code
  // In real implementation, this would use proper AVX2 intrinsics
  for (int i = 0; i < 32; i += 8) {
    for (int j = 0; j < 8; j++) {
      char c = buffer[31 - (i + j)];
      switch (c) {
      case 'A':
      case 'a':
        rev_buffer[i + j] = 'T';
        break;
      case 'C':
      case 'c':
        rev_buffer[i + j] = 'G';
        break;
      case 'G':
      case 'g':
        rev_buffer[i + j] = 'C';
        break;
      case 'T':
      case 't':
        rev_buffer[i + j] = 'A';
        break;
      default:
        rev_buffer[i + j] = c;
      }
    }
  }
#else
  // Scalar fallback with manual loop unrolling for better performance
  for (int i = 0; i < 32; i++) {
    char c = buffer[31 - i];
    switch (c) {
    case 'A':
    case 'a':
      rev_buffer[i] = 'T';
      break;
    case 'C':
    case 'c':
      rev_buffer[i] = 'G';
      break;
    case 'G':
    case 'g':
      rev_buffer[i] = 'C';
      break;
    case 'T':
    case 't':
      rev_buffer[i] = 'A';
      break;
    default:
      rev_buffer[i] = c;
    }
  }
#endif

  // Pack the reverse complement
  uint32_t reverse_packed1 = CompactKmer<16>::pack(rev_buffer);
  uint32_t reverse_packed2 = CompactKmer<16>::pack(rev_buffer + 16);
  uint64_t reverse_packed =
      (static_cast<uint64_t>(reverse_packed1) << 32) | reverse_packed2;

  // Compare and select smaller using optimized 64-bit comparison
  if (forward_packed <= reverse_packed) {
    resultHash = CompactKmer<32>::hash(forward_packed);
    resultIsReverse = false;
  } else {
    resultHash = CompactKmer<32>::hash(reverse_packed);
    resultIsReverse = true;
  }

  return true;
}

// General dispatcher function for getSeedAt
bool dispatchGetSeedAt(coordinates::CoordinateTraverser &traverser,
                       coordinates::CoordinateManager &manager,
                       size_t &resultHash, bool &resultIsReverse,
                       int64_t &resultEndPos, const int64_t &pos, Tree *T,
                       const int32_t &k) {

  // Dispatch to specialized implementations based on k-mer size
  // Using exact matches for best performance
  switch (k) {
  case 8:
    return getSeedAtSpecialized<8>(traverser, manager, resultHash,
                                   resultIsReverse, resultEndPos, pos, T);
  case 16:
    return getSeedAtSpecialized<16>(traverser, manager, resultHash,
                                    resultIsReverse, resultEndPos, pos, T);
  case 32:
    return getSeedAtSpecialized<32>(traverser, manager, resultHash,
                                    resultIsReverse, resultEndPos, pos, T);
  }

  // For other sizes, suggest using an optimal size instead
  int suggestedK = suggestOptimalKmerSize(k);
  if (k != suggestedK) {
    static bool warned = false;
    if (!warned) {
      msg("Non-optimal k-mer size detected: {} (consider {} instead)", k,
          suggestedK);
      warned = true; // Only warn once per run
    }
  }

  // Fallback to standard implementation
  return seed_annotated_tree::getSeedAt(
      traverser, manager, resultHash, resultIsReverse, resultEndPos, pos, T, k);
}

// Specialized batch processing for fixed k-mer sizes
template <int K>
void getSeedsBatchFixed(coordinates::CoordinateTraverser &traverser,
                        coordinates::CoordinateManager &manager,
                        const int64_t *positions, size_t positionCount,
                        size_t *resultHashes, bool *resultIsReverse,
                        int64_t *resultEndPos, bool *validResults) {

  // Ensure we're using a supported k-mer size
  static_assert(is_optimized_size<K>,
                "K must be an optimized size (8, 16, 32)");

  // Ensure all output arrays are properly initialized
  for (size_t i = 0; i < positionCount; i++) {
    validResults[i] = false;
  }

  // Allocate aligned buffer for SIMD operations using tbbmalloc
  memory::AlignedBuffer<char> buffer(K * positionCount, 64);

  // For k=8, we can use AVX2 to process 8 k-mers at once or SSE4.1 for 4 at
  // once
  if constexpr (K == 8) {
#ifdef __AVX2__
    // Process 8 positions at a time using AVX2
    for (size_t i = 0; i < positionCount; i += 8) {
      size_t batchSize = std::min(size_t(8), positionCount - i);

      // Prefetch next batch of positions
      if (i + 16 < positionCount) {
        __builtin_prefetch(&positions[i + 16], 0, 3);
      }

      // Process each position in this batch
      processKmerBatchAVX2<K>(traverser, manager, &positions[i], batchSize,
                              &resultHashes[i], &resultIsReverse[i],
                              &resultEndPos[i], &validResults[i],
                              buffer.get() + (i * K));
    }
#elif defined(__SSE4_1__)
    // Process 4 positions at a time using SSE4.1
    for (size_t i = 0; i < positionCount; i += 4) {
      size_t batchSize = std::min(size_t(4), positionCount - i);

      // Prefetch next batch of positions
      if (i + 8 < positionCount) {
        __builtin_prefetch(&positions[i + 8], 0, 3);
      }

      // Process each position in this batch
      processKmerBatchSSE<K>(traverser, manager, &positions[i], batchSize,
                             &resultHashes[i], &resultIsReverse[i],
                             &resultEndPos[i], &validResults[i],
                             buffer.get() + (i * K));
    }
#else
    // Fallback to scalar processing
    for (size_t i = 0; i < positionCount; i++) {
      bool success = getSeedAtSpecialized<K>(
          traverser, manager, resultHashes[i], resultIsReverse[i],
          resultEndPos[i], positions[i], nullptr);
      validResults[i] = success;
    }
#endif
  } else if constexpr (K == 16) {
// For k=16, process pairs of k-mers when SSE4.1 is available
#ifdef __SSE4_1__
    for (size_t i = 0; i < positionCount; i += 2) {
      size_t batchSize = std::min(size_t(2), positionCount - i);

      // Prefetch next batch
      if (i + 4 < positionCount) {
        __builtin_prefetch(&positions[i + 4], 0, 3);
      }

      // Process pair of positions
      for (size_t j = 0; j < batchSize; j++) {
        bool success = getSeedAtSpecialized<K>(
            traverser, manager, resultHashes[i + j], resultIsReverse[i + j],
            resultEndPos[i + j], positions[i + j], nullptr);
        validResults[i + j] = success;
      }
    }
#else
    // Fallback to scalar processing
    for (size_t i = 0; i < positionCount; i++) {
      bool success = getSeedAtSpecialized<K>(
          traverser, manager, resultHashes[i], resultIsReverse[i],
          resultEndPos[i], positions[i], nullptr);
      validResults[i] = success;
    }
#endif
  } else {
    // For k=32, process one k-mer at a time with prefetching
    for (size_t i = 0; i < positionCount; i++) {
      // Prefetch next position
      if (i + 1 < positionCount) {
        __builtin_prefetch(&positions[i + 1], 0, 3);
      }

      bool success = getSeedAtSpecialized<K>(
          traverser, manager, resultHashes[i], resultIsReverse[i],
          resultEndPos[i], positions[i], nullptr);
      validResults[i] = success;
    }
  }

  // AlignedBuffer will automatically free the memory with tbbmalloc
}

// Optimized batch processing using AVX2 for 8-mers
template <int K>
void processKmerBatchAVX2(coordinates::CoordinateTraverser &traverser,
                          coordinates::CoordinateManager &manager,
                          const int64_t *positions, size_t batchSize,
                          size_t *resultHashes, bool *resultIsReverse,
                          int64_t *resultEndPos, bool *validResults,
                          char *buffer) {

  static_assert(K == 8, "AVX2 batch processing only implemented for K=8");

#ifdef __AVX2__
  // Process up to 8 positions at once using AVX2
  for (size_t i = 0; i < batchSize; i += 8) {
    size_t currentBatchSize = std::min(size_t(8), batchSize - i);

    // Process each position individually without dynamic SIMD operations
    for (size_t j = 0; j < currentBatchSize; j++) {
      int64_t pos = positions[i + j];
      validResults[i + j] = false;

      // Skip positions in gap regions
      if (!manager.isGapPosition(pos)) {
        // Get starting coordinate
        const coordinates::tupleCoord_t *startCoord =
            manager.getTupleCoord(pos);
        if (startCoord && manager.getBlockExists()[startCoord->blockId].first) {
          // Prefetch coord data for better cache behavior
          __builtin_prefetch(startCoord, 0, 3);

          // Reset traverser to process this position individually
          traverser.reset(startCoord);

          // Get k-mer at this position - we're using traverser here but in a
          // fully optimized implementation, this would use direct SIMD across
          // multiple positions
          char *current_buffer = buffer + ((i + j) * K);
          bool foundValid = traverser.skipToNthNonGap(K, current_buffer);

          if (foundValid) {
            // Record end position and set validity
            resultEndPos[i + j] = traverser.getCurrent()->scalar;
            validResults[i + j] = true;
          } else {
            validResults[i + j] = false;
          }
        }
      }
    }

    // Process DNA encoding and comparison in SIMD when possible
    // This is where a full implementation would use AVX2 to process multiple
    // k-mers at once
    for (size_t j = 0; j < currentBatchSize; j++) {
      if (validResults[i + j]) {
        char *current_buffer = buffer + ((i + j) * K);

        // Process DNA nucleotides with SIMD (simplified here)
        // In a full implementation, we'd load multiple k-mers and process them
        // in parallel using intrinsics like _mm256_set_epi8, _mm256_cmpeq_epi8,
        // etc.

        // Calculate forward hash
        size_t fHash = 0;
        for (int k = 0; k < K; k++) {
          size_t cval = 0;
          switch (current_buffer[k]) {
          case 'A':
          case 'a':
            cval = 0x3c8bfbb395c60474;
            break;
          case 'C':
          case 'c':
            cval = 0x3193c18562a02b4c;
            break;
          case 'G':
          case 'g':
            cval = 0x20323ed082572324;
            break;
          case 'T':
          case 't':
            cval = 0x295549f54be24456;
            break;
          default:
            cval = 0;
            break;
          }
          fHash ^= seeding::rol(cval, K - k - 1);
        }

        // Create reverse complement and calculate hash
        alignas(32) char revBuffer[K];
        for (int k = 0; k < K; k++) {
          char c = current_buffer[K - 1 - k];
          switch (c) {
          case 'A':
          case 'a':
            revBuffer[k] = 'T';
            break;
          case 'C':
          case 'c':
            revBuffer[k] = 'G';
            break;
          case 'G':
          case 'g':
            revBuffer[k] = 'C';
            break;
          case 'T':
          case 't':
            revBuffer[k] = 'A';
            break;
          default:
            revBuffer[k] = c;
            break;
          }
        }

        size_t rHash = 0;
        for (int k = 0; k < K; k++) {
          size_t cval = 0;
          switch (revBuffer[k]) {
          case 'A':
          case 'a':
            cval = 0x3c8bfbb395c60474;
            break;
          case 'C':
          case 'c':
            cval = 0x3193c18562a02b4c;
            break;
          case 'G':
          case 'g':
            cval = 0x20323ed082572324;
            break;
          case 'T':
          case 't':
            cval = 0x295549f54be24456;
            break;
          default:
            cval = 0;
            break;
          }
          rHash ^= seeding::rol(cval, K - k - 1);
        }

        // Select canonical hash (min of forward and reverse)
        if (fHash <= rHash) {
          resultHashes[i + j] = fHash;
          resultIsReverse[i + j] = false;
        } else {
          resultHashes[i + j] = rHash;
          resultIsReverse[i + j] = true;
        }
      }
    }
  }
#endif
}

// Optimized batch processing using SSE4.1 for 8-mers
template <int K>
void processKmerBatchSSE(coordinates::CoordinateTraverser &traverser,
                         coordinates::CoordinateManager &manager,
                         const int64_t *positions, size_t batchSize,
                         size_t *resultHashes, bool *resultIsReverse,
                         int64_t *resultEndPos, bool *validResults,
                         char *buffer) {

  static_assert(K == 8, "SSE batch processing only implemented for K=8");

#ifdef __SSE4_1__
  // Process up to 4 positions at once using SSE4.1
  for (size_t i = 0; i < batchSize; i += 4) {
    size_t currentBatchSize = std::min(size_t(4), batchSize - i);

    // Process each position individually - avoiding dynamic SIMD operations
    for (size_t j = 0; j < currentBatchSize; j++) {
      int64_t pos = positions[i + j];
      validResults[i + j] = false;

      // Skip positions in gap regions
      if (!manager.isGapPosition(pos)) {
        // Get starting coordinate
        const coordinates::tupleCoord_t *startCoord =
            manager.getTupleCoord(pos);
        if (startCoord && manager.getBlockExists()[startCoord->blockId].first) {
          // Reset traverser to process this position
          traverser.reset(startCoord);

          // Get k-mer at this position
          char *current_buffer = buffer + ((i + j) * K);
          bool foundValid = traverser.skipToNthNonGap(K, current_buffer);

          if (foundValid) {
            // Record end position and set validity
            resultEndPos[i + j] = traverser.getCurrent()->scalar;
            validResults[i + j] = true;
          }
        }
      }
    }

    // Process valid k-mers
    for (size_t j = 0; j < currentBatchSize; j++) {
      if (validResults[i + j]) {
        char *current_buffer = buffer + ((i + j) * K);

        // Calculate forward hash
        size_t fHash =
            CompactKmer<K>::hash(CompactKmer<K>::pack(current_buffer));

        // Create reverse complement and calculate hash
        alignas(16) char revBuffer[K];
        for (int k = 0; k < K; k++) {
          char c = current_buffer[K - 1 - k];
          switch (c) {
          case 'A':
          case 'a':
            revBuffer[k] = 'T';
            break;
          case 'C':
          case 'c':
            revBuffer[k] = 'G';
            break;
          case 'G':
          case 'g':
            revBuffer[k] = 'C';
            break;
          case 'T':
          case 't':
            revBuffer[k] = 'A';
            break;
          default:
            revBuffer[k] = c;
            break;
          }
        }

        // Calculate reverse hash
        size_t rHash = CompactKmer<K>::hash(CompactKmer<K>::pack(revBuffer));

        // Select canonical hash (min of forward and reverse)
        if (fHash <= rHash) {
          resultHashes[i + j] = fHash;
          resultIsReverse[i + j] = false;
        } else {
          resultHashes[i + j] = rHash;
          resultIsReverse[i + j] = true;
        }
      }
    }
  }
#endif
}

// Dispatch function for batch processing
void dispatchGetSeedsBatch(coordinates::CoordinateTraverser &traverser,
                           coordinates::CoordinateManager &manager,
                           const int64_t *positions, size_t positionCount,
                           int32_t k, size_t *resultHashes,
                           bool *resultIsReverse, int64_t *resultEndPos,
                           bool *validResults) {

  // Dispatch to specialized implementations based on k-mer size
  // Using exact matches for best performance
  switch (k) {
  case 8:
    getSeedsBatchFixed<8>(traverser, manager, positions, positionCount,
                          resultHashes, resultIsReverse, resultEndPos,
                          validResults);
    return;
  case 16:
    getSeedsBatchFixed<16>(traverser, manager, positions, positionCount,
                           resultHashes, resultIsReverse, resultEndPos,
                           validResults);
    return;
  case 32:
    getSeedsBatchFixed<32>(traverser, manager, positions, positionCount,
                           resultHashes, resultIsReverse, resultEndPos,
                           validResults);
    return;
  }

  // For other sizes, suggest using an optimal size
  int suggestedK = suggestOptimalKmerSize(k);
  if (k != suggestedK) {
    static bool warned = false;
    if (!warned) {
      msg("Non-optimal k-mer size detected in batch: {} (consider {} instead)",
          k, suggestedK);
      warned = true; // Only warn once per run
    }
  }

  // No specialized implementation available for this k-mer size
  // Initialize all results as invalid
  for (size_t i = 0; i < positionCount; i++) {
    validResults[i] = false;
  }
  return;
}

} // namespace fixed_kmer
#endif // FIXED_KMER_HPP
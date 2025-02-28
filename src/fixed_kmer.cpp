#include "fixed_kmer.hpp"
#include "coordinates.hpp"
#include "seed_annotated_tree.hpp"

// This file contains explicit instantiations of template functions
// from fixed_kmer.hpp

namespace fixed_kmer {
// Implementation of suggestOptimalKmerSize
int suggestOptimalKmerSize(int requestedSize) {
  // Early return if the requested size is already optimal
  if (requestedSize <= 0)
    return 8; // Default to smallest optimal size

  if (requestedSize == 8 || requestedSize == 16 || requestedSize == 32) {
    return requestedSize;
  }

  // Find the closest optimal size
  if (requestedSize < 12)
    return 8;
  if (requestedSize < 24)
    return 16;
  return 32; // Default to largest size for anything >= 24
}

// Explicit instantiation of each template for the specific sizes we need
template void getSeedsBatchFixed<8>(coordinates::CoordinateTraverser &,
                                    coordinates::CoordinateManager &,
                                    const int64_t *, size_t, size_t *, bool *,
                                    int64_t *, bool *);

template void getSeedsBatchFixed<16>(coordinates::CoordinateTraverser &,
                                     coordinates::CoordinateManager &,
                                     const int64_t *, size_t, size_t *, bool *,
                                     int64_t *, bool *);

template void getSeedsBatchFixed<32>(coordinates::CoordinateTraverser &,
                                     coordinates::CoordinateManager &,
                                     const int64_t *, size_t, size_t *, bool *,
                                     int64_t *, bool *);

// Implementation of processKmerBatchAVX2 for K=8
template <>
void processKmerBatchAVX2<8>(coordinates::CoordinateTraverser &traverser,
                             coordinates::CoordinateManager &manager,
                             const int64_t *positions, size_t batchSize,
                             size_t *resultHashes, bool *resultIsReverse,
                             int64_t *resultEndPos, bool *validResults,
                             char *buffer) {

#ifdef __AVX2__
  // Process up to 8 positions at once using AVX2
  for (size_t i = 0; i < batchSize; i += 8) {
    size_t currentBatchSize = std::min(size_t(8), batchSize - i);

    // Process each position individually - avoid SIMD mask operations with
    // variable indices
    for (size_t j = 0; j < currentBatchSize; j++) {
      int64_t pos = positions[i + j];
      validResults[i + j] = false;

      // Skip validation for performance if position is in gap
      if (!manager.isGapPosition(pos)) {
        // Get starting coordinate
        const coordinates::tupleCoord_t *startCoord =
            manager.getTupleCoord(pos);
        if (startCoord && manager.getBlockExists()[startCoord->blockId].first) {
          // Reset traverser to process this position
          traverser.reset(startCoord);

          // Get k-mer at this position
          char *current_buffer = buffer + ((i + j) * 8);
          bool foundValid = traverser.skipToNthNonGap(8, current_buffer);

          if (foundValid) {
            // Record end position and set validity
            resultEndPos[i + j] = traverser.getCurrent()->scalar;
            validResults[i + j] = true;

            // Pack k-mer using the specialized function
            auto forward_packed = CompactKmer<8>::pack(current_buffer);

            // Calculate forward hash
            size_t fHash = CompactKmer<8>::hash(forward_packed);

            // Create reverse complement
            char rev_buffer[8];
            for (int k = 0; k < 8; k++) {
              char c = current_buffer[7 - k];
              switch (c) {
              case 'A':
              case 'a':
                rev_buffer[k] = 'T';
                break;
              case 'C':
              case 'c':
                rev_buffer[k] = 'G';
                break;
              case 'G':
              case 'g':
                rev_buffer[k] = 'C';
                break;
              case 'T':
              case 't':
                rev_buffer[k] = 'A';
                break;
              default:
                rev_buffer[k] = c;
                break;
              }
            }

            // Pack reverse complement
            auto reverse_packed = CompactKmer<8>::pack(rev_buffer);

            // Calculate reverse hash
            size_t rHash = CompactKmer<8>::hash(reverse_packed);

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
    }
  }
#else
  // Fallback implementation when AVX2 is not available
  for (size_t i = 0; i < batchSize; i++) {
    validResults[i] = false;
    int64_t pos = positions[i];

    // Skip validation for performance if position is in gap
    if (!manager.isGapPosition(pos)) {
      // Get starting coordinate
      const coordinates::tupleCoord_t *startCoord = manager.getTupleCoord(pos);
      if (startCoord && manager.getBlockExists()[startCoord->blockId].first) {
        // Reset traverser to process this position
        traverser.reset(startCoord);

        // Get k-mer at this position
        char *current_buffer = buffer + (i * 8);
        bool foundValid = traverser.skipToNthNonGap(8, current_buffer);

        if (foundValid) {
          // Record end position and set validity
          resultEndPos[i] = traverser.getCurrent()->scalar;
          validResults[i] = true;

          // Simple hash calculation without SIMD
          size_t fHash = 0;
          size_t rHash = 0;

          // Calculate forward hash
          auto forward_packed = CompactKmer<8>::pack(current_buffer);
          fHash = CompactKmer<8>::hash(forward_packed);

          // Create reverse complement
          char rev_buffer[8];
          for (int k = 0; k < 8; k++) {
            char c = current_buffer[7 - k];
            switch (c) {
            case 'A':
            case 'a':
              rev_buffer[k] = 'T';
              break;
            case 'C':
            case 'c':
              rev_buffer[k] = 'G';
              break;
            case 'G':
            case 'g':
              rev_buffer[k] = 'C';
              break;
            case 'T':
            case 't':
              rev_buffer[k] = 'A';
              break;
            default:
              rev_buffer[k] = c;
              break;
            }
          }

          // Calculate reverse hash
          auto reverse_packed = CompactKmer<8>::pack(rev_buffer);
          rHash = CompactKmer<8>::hash(reverse_packed);

          // Select canonical hash
          if (fHash <= rHash) {
            resultHashes[i] = fHash;
            resultIsReverse[i] = false;
          } else {
            resultHashes[i] = rHash;
            resultIsReverse[i] = true;
          }
        }
      }
    }
  }
#endif
}

// Implementation of processKmerBatchSSE for K=8
template <>
void processKmerBatchSSE<8>(coordinates::CoordinateTraverser &traverser,
                            coordinates::CoordinateManager &manager,
                            const int64_t *positions, size_t batchSize,
                            size_t *resultHashes, bool *resultIsReverse,
                            int64_t *resultEndPos, bool *validResults,
                            char *buffer) {

#ifdef __SSE4_1__
  // Process up to 4 positions at once using SSE4.1
  for (size_t i = 0; i < batchSize; i += 4) {
    size_t currentBatchSize = std::min(size_t(4), batchSize - i);

    // Process each position individually
    for (size_t j = 0; j < currentBatchSize; j++) {
      int64_t pos = positions[i + j];
      validResults[i + j] = false;

      // Skip validation for performance if position is in gap
      if (!manager.isGapPosition(pos)) {
        // Get starting coordinate
        const coordinates::tupleCoord_t *startCoord =
            manager.getTupleCoord(pos);
        if (startCoord && manager.getBlockExists()[startCoord->blockId].first) {
          // Reset traverser to process this position
          traverser.reset(startCoord);

          // Get k-mer at this position
          char *current_buffer = buffer + ((i + j) * 8);
          bool foundValid = traverser.skipToNthNonGap(8, current_buffer);

          if (foundValid) {
            // Record end position and set validity
            resultEndPos[i + j] = traverser.getCurrent()->scalar;
            validResults[i + j] = true;

            // Pack k-mer using the specialized function
            auto forward_packed = CompactKmer<8>::pack(current_buffer);

            // Calculate forward hash
            size_t fHash = CompactKmer<8>::hash(forward_packed);

            // Create reverse complement
            char rev_buffer[8];
            for (int k = 0; k < 8; k++) {
              char c = current_buffer[7 - k];
              switch (c) {
              case 'A':
              case 'a':
                rev_buffer[k] = 'T';
                break;
              case 'C':
              case 'c':
                rev_buffer[k] = 'G';
                break;
              case 'G':
              case 'g':
                rev_buffer[k] = 'C';
                break;
              case 'T':
              case 't':
                rev_buffer[k] = 'A';
                break;
              default:
                rev_buffer[k] = c;
                break;
              }
            }

            // Pack reverse complement
            auto reverse_packed = CompactKmer<8>::pack(rev_buffer);

            // Calculate reverse hash
            size_t rHash = CompactKmer<8>::hash(reverse_packed);

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
    }
  }
#else
  // Fallback implementation when SSE4.1 is not available
  for (size_t i = 0; i < batchSize; i++) {
    validResults[i] = false;
    int64_t pos = positions[i];

    // Skip validation for performance if position is in gap
    if (!manager.isGapPosition(pos)) {
      // Get starting coordinate
      const coordinates::tupleCoord_t *startCoord = manager.getTupleCoord(pos);
      if (startCoord && manager.getBlockExists()[startCoord->blockId].first) {
        // Reset traverser to process this position
        traverser.reset(startCoord);

        // Get k-mer at this position
        char *current_buffer = buffer + (i * 8);
        bool foundValid = traverser.skipToNthNonGap(8, current_buffer);

        if (foundValid) {
          // Record end position and set validity
          resultEndPos[i] = traverser.getCurrent()->scalar;
          validResults[i] = true;

          // Simple hash calculation without SIMD
          size_t fHash = 0;
          size_t rHash = 0;

          // Calculate forward hash
          auto forward_packed = CompactKmer<8>::pack(current_buffer);
          fHash = CompactKmer<8>::hash(forward_packed);

          // Create reverse complement
          char rev_buffer[8];
          for (int k = 0; k < 8; k++) {
            char c = current_buffer[7 - k];
            switch (c) {
            case 'A':
            case 'a':
              rev_buffer[k] = 'T';
              break;
            case 'C':
            case 'c':
              rev_buffer[k] = 'G';
              break;
            case 'G':
            case 'g':
              rev_buffer[k] = 'C';
              break;
            case 'T':
            case 't':
              rev_buffer[k] = 'A';
              break;
            default:
              rev_buffer[k] = c;
              break;
            }
          }

          // Calculate reverse hash
          auto reverse_packed = CompactKmer<8>::pack(rev_buffer);
          rHash = CompactKmer<8>::hash(reverse_packed);

          // Select canonical hash
          if (fHash <= rHash) {
            resultHashes[i] = fHash;
            resultIsReverse[i] = false;
          } else {
            resultHashes[i] = rHash;
            resultIsReverse[i] = true;
          }
        }
      }
    }
  }
#endif
}
} // namespace fixed_kmer
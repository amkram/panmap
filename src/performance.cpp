#include "performance.hpp"
#include <cmath> // For std::abs

// Include the template definitions fully
// The error occurs because the compiler can't see the complete template
// definitions when trying to instantiate them

namespace memory {
// Define thread_local variables for KmerBufferPool
thread_local char *KmerBufferPool::nucleotide_buffer = nullptr;
thread_local size_t *KmerBufferPool::hash_buffer = nullptr;
thread_local bool *KmerBufferPool::orientation_buffer = nullptr;
thread_local int64_t *KmerBufferPool::position_buffer = nullptr;
thread_local size_t KmerBufferPool::buffer_size = 0;
} // namespace memory

// Add implementation of suggestOptimalKmerSize that was in fixed_kmer.cpp
namespace fixed_kmer {
    // This function suggests an optimal k-mer size for hardware efficiency
    // Previously in fixed_kmer.cpp, now moved here to resolve linking issues
    int32_t suggestOptimalKmerSize(int32_t requestedK) {
        // Common optimal k-mer sizes for hardware efficiency
        // Prefer sizes that are multiples of 8, 16, or 32 for better memory alignment
        constexpr int32_t optimalSizes[] = {8, 16, 24, 32, 48, 64};
        
        // Find the closest optimal size to the requested k
        int32_t bestK = optimalSizes[0];
        int32_t minDiff = std::abs(requestedK - bestK);
        
        for (auto size : optimalSizes) {
            int32_t diff = std::abs(requestedK - size);
            if (diff < minDiff) {
                minDiff = diff;
                bestK = size;
            }
            // If we're already past the requested k by a significant margin, stop
            if (size > requestedK && diff > minDiff) {
                break;
            }
        }
        
        return bestK;
    }
} // namespace fixed_kmer

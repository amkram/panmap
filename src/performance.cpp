#include "performance.hpp"
#include "fixed_kmer.hpp"
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

// Implementation of non-template functions in the fixed_kmer namespace
namespace fixed_kmer {
// This file no longer needs the explicit template instantiations
// since they're now fully implemented in fixed_kmer.cpp
}

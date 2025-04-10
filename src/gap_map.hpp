#pragma once

#include <boost/icl/interval_set.hpp>
#include <cstddef>
#include <utility>
#include <vector>
#include <map>
#include <algorithm>
#include <immintrin.h> // For SIMD instructions
#include <cstdint>

// Include full definitions from coordinates.hpp instead of forward declarations
#include "coordinates.hpp"

namespace gap_map {

// Define the gap map type
using GapMap = std::map<int64_t, int64_t>;
using GapRange = std::pair<int64_t, int64_t>;
using GapUpdate = std::pair<bool, GapRange>;

// Move namespace boost::icl usage outside of any namespace
using namespace boost::icl;

// Compressed Run-Length Encoding for Gap Runs
struct CompressedGapRun {
    uint32_t position : 30;  // Position (30 bits allows up to 1B positions)
    uint32_t type : 2;       // Type of run (0=short, 1=medium, 2=long, 3=very long)
    
    union {
        uint8_t shortLength;     // For runs <= 255
        uint16_t mediumLength;   // For runs <= 65535
        uint32_t longLength;     // For runs > 65535
    };
    
    // Helper methods for efficient size determination
    static CompressedGapRun create(int64_t start, int64_t length) {
        CompressedGapRun run;
        run.position = start;
        
        if (length <= 255) {
            run.type = 0;
            run.shortLength = length;
        } else if (length <= 65535) {
            run.type = 1;
            run.mediumLength = length;
        } else {
            run.type = 2;
            run.longLength = length;
        }
        
        return run;
    }
    
    int64_t getLength() const {
        switch (type) {
            case 0: return shortLength;
            case 1: return mediumLength;
            case 2: return longLength;
            default: return 0;
        }
    }
    
    int64_t getEndPosition() const {
        return position + getLength() - 1;
    }
    
    bool contains(int64_t pos) const {
        return pos >= position && pos <= getEndPosition();
    }
};

// Statistics for gap maps to enable adaptive optimization
struct GapStatistics {
    size_t totalGapRuns = 0;
    size_t totalGapPositions = 0;
    size_t maxGapRunLength = 0;
    double avgGapRunLength = 0.0;
    size_t smallGapRuns = 0;  // < 10 positions
    size_t mediumGapRuns = 0; // 10-100 positions
    size_t largeGapRuns = 0;  // > 100 positions
    
    void update(const GapMap& gapMap) {
        totalGapRuns = gapMap.size();
        totalGapPositions = 0;
        maxGapRunLength = 0;
        smallGapRuns = 0;
        mediumGapRuns = 0;
        largeGapRuns = 0;
        
        for (const auto& [pos, length] : gapMap) {
            totalGapPositions += length;
            maxGapRunLength = std::max(maxGapRunLength, static_cast<size_t>(length));
            
            if (length < 10) smallGapRuns++;
            else if (length < 100) mediumGapRuns++;
            else largeGapRuns++;
        }
        
        avgGapRunLength = totalGapRuns > 0 ? 
            static_cast<double>(totalGapPositions) / totalGapRuns : 0.0;
    }
    
    bool shouldUseCompression() const {
        // Use compression if we have many gap runs or large gaps
        return totalGapRuns > 1000 || maxGapRunLength > 10000;
    }
    
    bool shouldUseArray() const {
        // Use array-based representation if we have mostly small gaps
        return smallGapRuns > 0.8 * totalGapRuns;
    }
};

// SIMD-accelerated gap detection for large blocks of sequence data
inline bool anyGaps_SIMD(const char* data, unsigned long size) {
    if (size < 16) {
        // Fallback for small blocks
        return std::any_of(data, data + size, [](char c) { return c == '-'; });
    }
    
    // Process 16 characters at once with SSE instructions
    __m128i gap_char = _mm_set1_epi8('-');
    for (size_t i = 0; i <= size - 16; i += 16) {
        __m128i block = _mm_loadu_si128(reinterpret_cast<const __m128i*>(data + i));
        __m128i cmp = _mm_cmpeq_epi8(block, gap_char);
        if (_mm_movemask_epi8(cmp) != 0) {
            return true;
        }
    }
    
    // Check any remaining characters
    for (size_t i = (size / 16) * 16; i < size; i++) {
        if (data[i] == '-') return true;
    }
    
    return false;
}

// Count gaps using SIMD for improved performance
inline bool countGaps_SIMD(const char* data, unsigned long size) {
    if (size < 16) {
        // Fallback for small blocks
        return std::count(data, data + size, '-');
    }
    
    size_t gapCount = 0;
    __m128i gap_char = _mm_set1_epi8('-');
    
    for (size_t i = 0; i <= size - 16; i += 16) {
        __m128i block = _mm_loadu_si128(reinterpret_cast<const __m128i*>(data + i));
        __m128i cmp = _mm_cmpeq_epi8(block, gap_char);
        unsigned mask = _mm_movemask_epi8(cmp);
        gapCount += __builtin_popcount(mask);
    }
    
    // Count any remaining characters
    for (size_t i = (size / 16) * 16; i < size; i++) {
        if (data[i] == '-') gapCount++;
    }
    
    return gapCount > 0;
}

// Find contiguous gap runs with SIMD for improved performance
// Returns vector of (start, length) pairs
inline bool findGapRuns_SIMD(const char* data, unsigned long size) {
    std::vector<std::pair<size_t, size_t>> gapRuns;
    
    if (size < 16) {
        // Manual processing for small sequences
        size_t i = 0;
        while (i < size) {
            // Find start of a gap run
            while (i < size && data[i] != '-') i++;
            if (i >= size) break;
            
            size_t start = i;
            // Find end of the gap run
            while (i < size && data[i] == '-') i++;
            
            gapRuns.emplace_back(start, i - start);
        }
        return !gapRuns.empty();
    }
    
    // SIMD-accelerated gap run detection
    __m128i gap_char = _mm_set1_epi8('-');
    size_t runStart = SIZE_MAX;  // Invalid position to indicate no run in progress
    
    for (size_t i = 0; i <= size - 16; i += 16) {
        __m128i block = _mm_loadu_si128(reinterpret_cast<const __m128i*>(data + i));
        __m128i cmp = _mm_cmpeq_epi8(block, gap_char);
        unsigned mask = _mm_movemask_epi8(cmp);
        
        // Process each bit in the mask (1 = gap, 0 = non-gap)
        for (unsigned j = 0; j < 16; j++) {
            bool isGap = (mask & (1 << j)) != 0;
            
            if (isGap && runStart == SIZE_MAX) {
                // Start of a new gap run
                runStart = i + j;
            }
            else if (!isGap && runStart != SIZE_MAX) {
                // End of a gap run
                gapRuns.emplace_back(runStart, i + j - runStart);
                runStart = SIZE_MAX;
            }
        }
    }
    
    // Process any remaining characters
    for (size_t i = (size / 16) * 16; i < size; i++) {
        bool isGap = data[i] == '-';
        
        if (isGap && runStart == SIZE_MAX) {
            // Start of a new gap run
            runStart = i;
        }
        else if (!isGap && runStart != SIZE_MAX) {
            // End of a gap run
            gapRuns.emplace_back(runStart, i - runStart);
            runStart = SIZE_MAX;
        }
    }
    
    // Handle case where a gap run extends to the end of the sequence
    if (runStart != SIZE_MAX) {
        gapRuns.emplace_back(runStart, size - runStart);
    }
    
    return !gapRuns.empty();
}

/**
 * @brief Convert a gap map to a nucleotide run set
 * @param gapMap Gap map to be converted
 * @param blockRanges Block range information
 * @return Nucleotide run set
 */
interval_set<int64_t>
gapMapToNucRunSet(const GapMap &gapMap,
                  const std::vector<std::pair<int64_t, int64_t>> &blockRanges);

/**
 * @brief Update a gap map with a single step
 * @param gapMap Gap map to be updated
 * @param update Update to be applied (deletion flag and range)
 * @param backtrack Vector to store backtrack information
 * @param gapMapUpdates Vector to store all updates
 * @param recordGapMapUpdates Whether to record updates
 * @param totalCoordinateCount Total number of coordinates for bounds checking
 */
void updateGapMapStep(GapMap &gapMap, const GapUpdate &update,
                      std::vector<GapUpdate> &backtrack,
                      std::vector<GapUpdate> &gapMapUpdates,
                      bool recordGapMapUpdates = true,
                      int64_t totalCoordinateCount = -1);

/**
 * @brief Update a gap map with multiple updates
 * @param gapMap Gap map to be updated
 * @param updates Vector of updates to be applied
 * @param backtrack Vector to store backtrack information
 * @param gapMapUpdates Vector to store all updates
 * @param totalCoordinateCount Total number of coordinates for bounds checking
 */
void updateGapMap(GapMap &gapMap, const std::vector<GapUpdate> &updates,
                  std::vector<GapUpdate> &backtrack,
                  std::vector<GapUpdate> &gapMapUpdates,
                  int64_t totalCoordinateCount = -1);

/**
 * @brief Simplified gap map update for single update
 * @param gapMap Gap map to be updated
 * @param update Single update to apply
 * @param backtrack Vector to store backtrack information
 * @param totalCoordinateCount Total number of coordinates for bounds checking
 */
void updateGapMap(GapMap &gapMap, const GapUpdate &update,
                  std::vector<GapUpdate> &backtrack,
                  int64_t totalCoordinateCount = -1);

/**
 * @brief Validate and fix issues in a gap map
 * @param gapMap Gap map to validate and fix
 * @param totalCoordinateCount Total number of coordinates for bounds checking
 * @param fixErrors Whether to fix errors or just report them
 * @return True if the gap map is valid (or was fixed), false otherwise
 */
bool validateAndFixGapMap(GapMap &gapMap, int64_t totalCoordinateCount,
                         bool fixErrors = true);

/**
 * @brief Invert ranges in a gap map
 * @param nucRanges Ranges to be inverted
 * @param invertRange Range within which to perform inversion
 * @return Vector of inverted ranges
 */
std::vector<GapRange> invertRanges(const std::vector<GapRange> &nucRanges,
                                   const GapRange &invertRange);

/**
 * @brief Invert gap map for a range
 * @param gapMap Gap map to be inverted
 * @param invertRange Range within which to perform inversion
 * @param backtrack Vector to store backtrack information
 */
void invertGapMap(GapMap &gapMap, const GapRange &invertRange,
                  std::vector<GapUpdate> &backtrack);

/**
 * @brief Invert gap map for a range and record updates
 * @param gapMap Gap map to be inverted
 * @param invertRange Range within which to perform inversion
 * @param backtrack Vector to store backtrack information
 * @param gapMapUpdates Vector to store all updates made
 */
void invertGapMap(GapMap &gapMap, const GapRange &invertRange,
                  std::vector<GapUpdate> &backtrack,
                  std::vector<GapUpdate> &gapMapUpdates);

/**
 * @brief Create coordinate indices for a gap map
 * @param degapCoordIndex Map for degapped coordinate lookup
 * @param regapCoordIndex Map for regapped coordinate lookup
 * @param gapMap Source gap map
 * @param blockRanges Block range information
 */
void makeCoordIndex(GapMap &degapCoordIndex, GapMap &regapCoordIndex,
                    const GapMap &gapMap,
                    const std::vector<GapRange> &blockRanges);

/**
 * @brief Detect and remove suspiciously large gaps from a gap map
 * @param gapMap Gap map to check and fix
 * @param totalSequenceSize Total sequence size for percentage calculation
 * @param threshold Percentage threshold (default 0.9 or 90% of sequence)
 * @return The number of suspicious gaps removed
 */
int removeSuspiciousGaps(GapMap &gapMap, int64_t totalSequenceSize, double threshold = 0.9);

/**
 * @brief Prevent overlaps in a gap map by merging adjacent/overlapping gaps
 * @param gapMap Gap map to check and fix
 * @return True if changes were made, false if no overlaps found
 */
bool preventOverlaps(GapMap &gapMap);

/**
 * @brief Batch process gap updates for better performance
 * @param gapMap Gap map to be updated
 * @param updates Vector of updates to be processed
 * @param backtrack Vector to store backtrack information
 * @param totalCoordinateCount Total number of coordinates for bounds checking
 * @return Number of actual updates applied after merging
 */
size_t batchProcessGapUpdates(GapMap &gapMap, 
                            const std::vector<GapUpdate> &updates,
                            std::vector<GapUpdate> &backtrack,
                            int64_t totalCoordinateCount = -1);

/**
 * @brief Calculate gap statistics for a gap map
 * @param gapMap Gap map to analyze
 * @return GapStatistics object with detailed statistics
 */
GapStatistics calculateGapStatistics(const GapMap &gapMap);

// Add a gap to the map
void addGap(GapMap& gapMap, int64_t pos, int64_t length);

// Remove a gap from the map
void removeGap(GapMap& gapMap, int64_t pos, int64_t length);

// Update gap map with a single update
void updateGapMap(GapMap& gapMap, const GapUpdate& update, std::vector<GapUpdate>& backtrack, int64_t totalCoords);

// Invert a section of the gap map
void invertGapMap(GapMap& gapMap, const GapRange& range, std::vector<GapUpdate>& backtrack);

// Prevent overlapping gaps by merging
bool preventOverlaps(GapMap& gapMap);

// Apply a single gap update to a map
inline void applyGapUpdate(GapMap& gapMap, const GapUpdate& update) {
    const auto& [isRemoval, range] = update;
    const auto& [start, end] = range;
    int64_t length = end - start + 1;
    
    if (isRemoval) {
        // Remove gap
        removeGap(gapMap, start, length);
    } else {
        // Add gap
        addGap(gapMap, start, length);
    }
}

} // namespace gap_map
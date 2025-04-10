#pragma once

#include "logging.hpp" // Add logging header
#include <cctype>
#include <cstdint>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "panman.hpp"


// Forward declare gap_map types to break circular dependency
namespace gap_map {
    using GapRange = std::pair<int64_t, int64_t>;
    using GapUpdate = std::pair<bool, GapRange>;
}

// Define position_key_t outside of any namespace
using position_key_t = std::tuple<int, int, int>;

// Define the hash function specialization for position_key_t
namespace std {
template <> struct hash<position_key_t> {
  inline size_t operator()(const position_key_t &key) const {
    // Combine the hash of the three integers
    auto hash1 = std::hash<int>{}(std::get<0>(key));
    auto hash2 = std::hash<int>{}(std::get<1>(key));
    auto hash3 = std::hash<int>{}(std::get<2>(key));

    // Use the same hash function as in panman's tuple_hash
    return hash1 ^ hash2 ^ hash3;
  }
};
} // namespace std


namespace panmanUtils {

typedef std::vector<
    std::pair<std::vector<std::pair<char, std::vector<char>>>,
              std::vector<std::vector<std::pair<char, std::vector<char>>>>>>
    sequence_t;

    typedef std::vector<std::pair<char, std::vector<char>>> block_t;
    typedef std::vector<std::pair<bool, std::vector<bool>>> blockExists_t;
    typedef std::vector<std::pair<bool, std::vector<bool>>> blockStrand_t;

    // Forward declarations for PanmanUtils types
    struct Block;
    struct BlockGapList;
    struct GapList;
} // namespace panmanUtils

namespace coordinates {

// Define the GapUpdate structure in the coordinates namespace before it's used
struct GapUpdate {
  int64_t pos;        // Position of the gap
  int64_t length;     // Length of the gap run
  bool isGapAddition; // True if this adds a gap, false if it removes a gap

  /**
   * @brief Construct a new Gap Update object
   *
   * @param p Position of the gap
   * @param l Length of the gap run
   * @param add True if this adds a gap, false if it removes a gap
   */
  inline GapUpdate(int64_t p, int64_t l, bool add)
      : pos(p), length(l), isGapAddition(add) {}
};

// Forward declarations
class CoordinateManager;

/**
 * @brief Helper method for validating scalars with consistent behavior
 * 
 * @param scalar Scalar coordinate to validate
 * @param total_size Total size of coordinate space
 * @param throwOnError Whether to throw an exception (true) or return false (false)
 * @return bool True if scalar is valid, false otherwise (if not throwing)
 */
inline bool validateScalar(int64_t scalar, int64_t total_size, bool throwOnError = false) {
    if (scalar < 0 || scalar >= total_size) {
        if (throwOnError) {
            std::string error_msg = "Invalid scalar: ";
            error_msg += std::to_string(scalar);
            error_msg += " (total_size: ";
            error_msg += std::to_string(total_size);
            error_msg += ")";
            throw std::runtime_error(error_msg);
        }
        return false;
    }
    return true;
}

/**
 * @brief Represents a position within a block in the genome
 * 
 * BlockPosition is a 3D coordinate system used to uniquely identify
 * locations within the genomic data, handling both gaps and nucleotides.
 */
struct BlockPosition {
  int32_t blockId; // Block identifier
  int32_t nucPos;  // Nucleotide position within the block
  int32_t gapPos;  // Gap position (-1 if not a gap)

  // Default constructor creates an invalid position
  inline BlockPosition() : blockId(-1), nucPos(-1), gapPos(-1) {}

  // Constructor with all parameters
  inline BlockPosition(int32_t b, int32_t n, int32_t g)
        : blockId(b), nucPos(n), gapPos(g) {}
    
    // Check if position is valid
  inline bool isValid() const {
        return blockId >= 0 && nucPos >= 0 && (gapPos >= 0 || gapPos == -1);
    }
    
  // Equality comparison
  inline bool operator==(const BlockPosition &other) const {
    return blockId == other.blockId && nucPos == other.nucPos &&
           gapPos == other.gapPos;
  }

  // Inequality comparison
  inline bool operator!=(const BlockPosition &other) const {
        return !(*this == other);
    }

  // String representation for debugging
  inline std::string toString() const {
    std::stringstream ss;
    ss << "BlockPosition{blockId=" << blockId << ", nucPos=" << nucPos
       << ", gapPos=" << gapPos << "}";
    return ss.str();
    }
};

struct CoordRange {
  int64_t start;
  int64_t end;
  
  bool operator==(const CoordRange& other) const {
    return start == other.start && end == other.end;
  }
  
  size_t size() const {
    return static_cast<size_t>(std::max<int64_t>(0, end - start));
  }
  
  bool contains(int64_t pos) const {
    return pos >= start && pos < end;
  }
  
  bool overlaps(const CoordRange& other) const {
    return start < other.end && end > other.start;
  }
  
  bool adjacent(const CoordRange& other) const {
    return end == other.start || start == other.end;
  }
};


extern bool debug;

// These functions are part of the coordinates namespace but outside the class
bool canMergeBlocks(int blockId1, int blockId2,
                  const panmanUtils::blockExists_t &blockExists);

std::pair<int, int>
findAdjacentLiveBlocks(int blockId,
                      const panmanUtils::blockExists_t &blockExists);

bool isValidRange(const std::pair<int64_t, int64_t> &range,
                const std::vector<std::pair<int64_t, int64_t>> &blockRanges);

std::vector<GapUpdate>
convertFromGapMapUpdates(const std::vector<gap_map::GapUpdate> &updates);

std::string scalarToString(int64_t scalar, const CoordinateManager &manager);

/**
 * @namespace RangeOperations
 * @brief Utility functions for range manipulation
 *
 * This namespace contains functions for working with CoordRange objects,
 * including merging, expansion, and verification.
 */
namespace RangeOperations {
  /**
   * @brief Merge overlapping or adjacent ranges
   * 
   * Sorts ranges by start position and merges any that overlap or are adjacent.
   * 
   * @param ranges Vector of ranges to merge
   * @return New vector containing merged ranges
   */
  inline std::vector<CoordRange> mergeRanges(std::vector<CoordRange> ranges) {
    if (ranges.empty()) {
      return ranges;
    }

    // Sort ranges by start position
    std::sort(ranges.begin(), ranges.end(),
              [](const auto &a, const auto &b) { return a.start < b.start; });

    // Merge overlapping or adjacent ranges
    std::vector<CoordRange> mergedRanges;
    mergedRanges.reserve(ranges.size());

    // Start with the first range
    mergedRanges.push_back(ranges.front());

    // Use range-based for loop with auto references
    for (size_t i = 1; i < ranges.size(); ++i) {
      auto &last = mergedRanges.back();
      const auto &current = ranges[i];

      // Check if current range overlaps or is adjacent to the last merged range
      if (current.start <= last.end + 1) {
        // Extend the last range if current range extends beyond it
        last.end = std::max(last.end, current.end);
      } else {
        // Otherwise, add as a new range
        mergedRanges.push_back(current);
      }
    }

    return mergedRanges;
  }
  
  /**
   * @brief Merge a new range with existing ranges
   * 
   * Adds a new range to a vector of ranges and merges if needed.
   * 
   * @param ranges Vector of existing ranges to merge with
   * @param newRange New range to add
   */
  inline void mergeRangeWithExisting(std::vector<CoordRange> &ranges,
                                const CoordRange &newRange) {
    // Validate the new range - ensure start <= end
    if (newRange.start > newRange.end) {
      return; // Invalid range
    }

    // If this is the first range, just add it
    if (ranges.empty()) {
      ranges.push_back(newRange);
      return;
    }

    // Add the new range and use the mergeRanges utility to efficiently merge
    // overlapping ranges
    ranges.push_back(newRange);
    ranges = mergeRanges(std::move(ranges));
  }
  
  /**
   * @brief Calculate recomputation range based on mutation parameters
   * 
   * Determines the range of coordinates that needs to be recomputed
   * after a mutation occurs.
   * 
   * @param blockRange Range of the affected block
   * @param pos Position within the block
   * @param len Length of the mutation
   * @param kmerSize Size of k-mers for ensuring coverage
   * @param isBlockMutation Whether this is a block-level mutation
   * @param isBlockDeactivation Whether this is a block deactivation
   * @param maxCoord Maximum coordinate value allowed
   * @param isInverted Whether the block is inverted
   * @return Recomputation range
   */
  inline CoordRange calculateRecompRange(
      const CoordRange& blockRange,
      int32_t pos,
      int32_t len,
      size_t kmerSize,
      bool isBlockMutation,
      bool isBlockDeactivation,
      int64_t maxCoord,
      bool isInverted) {
      
    CoordRange result;
    
    // Different handling based on mutation type
    if (isBlockMutation) {
      if (isBlockDeactivation) {
        // On->Off block mutation (deactivation)
        // Only include k positions at block boundaries
        // Start boundary
        result.start = blockRange.start;
        result.end = std::min(blockRange.start + static_cast<int64_t>(kmerSize), blockRange.end);
        
        // This function returns one range, so we'll modify the caller to handle the second range
      } else {
        // Off->On or inversion change
        // Include the entire block
        result.start = blockRange.start;
        result.end = blockRange.end;
      }
    } else {
      // Nucleotide mutation
      // Add k positions on either side, but don't cross block boundaries
      int64_t mutationPos;
      if (isInverted) {
        // For inverted blocks, the position is mirrored
        int32_t blockLength = blockRange.end - blockRange.start;
        mutationPos = blockRange.start + (blockLength - 1 - pos);
      } else {
        // For normal blocks, simply add the position to the block start
        mutationPos = blockRange.start + pos;
      }
      
      // Add k positions in each direction, constrained by block boundaries
      result.start = std::max(blockRange.start, mutationPos - static_cast<int64_t>(kmerSize));
      result.end = std::min(blockRange.end, mutationPos + len + static_cast<int64_t>(kmerSize));
    }
    
    // Ensure result is within valid bounds
    result.start = std::max<int64_t>(0, result.start);
    result.end = std::min<int64_t>(maxCoord, result.end);
    
    return result;
  }
  
  /**
   * @brief Expand ranges to ensure they contain enough non-gap characters
   * 
   * Verifies that ranges have at least k non-gap characters at boundaries,
   * and expands them if needed.
   * 
   * @param ranges Ranges to expand
   * @param k K-mer size requiring non-gap characters
   * @param findNonGapFn Function that finds positions with non-gap characters
   * @param maxCoord Maximum coordinate value allowed
   * @return Expanded ranges that have been merged if they overlap
   */
  template<typename FindNonGapFn>
  inline std::vector<CoordRange> expandRanges(
      const std::vector<CoordRange>& ranges,
      int k,
      FindNonGapFn findNonGapFn,
      int64_t maxCoord) {
      
    // Early return for empty ranges or k=0
    if (ranges.empty() || k <= 0) {
      return ranges;
    }
    
    // First merge the input ranges to avoid redundant processing
    auto mergedRanges = mergeRanges(std::vector<CoordRange>(ranges));
    
    // Create result vector
    std::vector<CoordRange> verifiedRanges;
    verifiedRanges.reserve(mergedRanges.size());
    
    // Process each merged range
    for (const auto& range : mergedRanges) {
      // Our initial guess should already be generous
      // But we need to verify there are at least k non-gap characters at each boundary
      CoordRange verifiedRange = range;
      
      // Expand start if needed
      int64_t expandedStart = findNonGapFn(range.start, k, false);
      verifiedRange.start = std::max<int64_t>(0, expandedStart);
      
      // Expand end if needed
      int64_t expandedEnd = findNonGapFn(range.end, k, true);
      verifiedRange.end = std::min<int64_t>(maxCoord - 1, expandedEnd);
      
      verifiedRanges.push_back(verifiedRange);
    }
    
    // Final merge in case expansions caused overlaps
    return mergeRanges(verifiedRanges);
  }
  
  /**
   * @brief Convert from GapRange (inclusive end) to CoordRange (exclusive end)
   * 
   * @param gapRanges Vector of GapRange objects with inclusive end
   * @return Vector of CoordRange objects with exclusive end
   */
  template<typename GapRange>
  inline std::vector<CoordRange> gapRangesToCoordRanges(
      const std::vector<GapRange>& gapRanges) {
    std::vector<CoordRange> coordRanges;
    coordRanges.reserve(gapRanges.size());
    
    for (const auto& range : gapRanges) {
      // Extract start and end (inclusive) from the range pair
      auto start = range.first;
      auto endInclusive = (range.second > 0) ? range.second : (start + range.second - 1);
      
      // Convert to CoordRange (exclusive end)
      coordRanges.emplace_back(start, endInclusive + 1);
    }
    
    return coordRanges;
  }
  
  /**
   * @brief Convert from CoordRange (exclusive end) to GapRange (inclusive end)
   * 
   * @param coordRanges Vector of CoordRange objects with exclusive end
   * @return Vector of pairs representing GapRange objects with inclusive end
   */
  template<typename GapRangePair>
  inline std::vector<GapRangePair> coordRangesToGapRanges(
      const std::vector<CoordRange>& coordRanges) {
    std::vector<GapRangePair> gapRanges;
    gapRanges.reserve(coordRanges.size());
    
    for (const auto& range : coordRanges) {
      // Convert from exclusive end to inclusive end
      auto start = range.start;
      auto length = range.end - range.start;
      
      // Create gap range pair (start, length)
      gapRanges.emplace_back(start, length);
    }
    
    return gapRanges;
  }
}

/**
 * @struct BlockCoordinate
 * @brief Represents a 3D coordinate in block space
 * 
 */
struct BlockCoordinate {
  int32_t blockId;  // ID of the block
  int32_t nucPos;   // Position within the block's nucleotide sequence
  int32_t gapPos;   // Gap position (-1 for main nucleotide, >= 0 for gap list positions)
  
  // Default constructor
  BlockCoordinate() : blockId(-1), nucPos(-1), gapPos(-1) {}
  
  // Full constructor
  BlockCoordinate(int32_t bid, int32_t np, int32_t gp) 
    : blockId(bid), nucPos(np), gapPos(gp) {}
    
  // Conversion from tuple
  BlockCoordinate(const std::tuple<int32_t, int32_t, int32_t>& t)
    : blockId(std::get<0>(t)), nucPos(std::get<1>(t)), 
      gapPos(std::get<2>(t)) {}
  
  // Conversion to tuple (for backward compatibility)
  operator std::tuple<int32_t, int32_t, int32_t>() const {
    return std::make_tuple(blockId, nucPos, gapPos);
  }
  
  // Check if this is a valid coordinate
  bool isValid() const {
    return blockId >= 0 && nucPos >= 0;
  }
  
  // Equality operator
  bool operator==(const BlockCoordinate& other) const {
    return blockId == other.blockId && 
           nucPos == other.nucPos && 
           gapPos == other.gapPos;
  }
};

} // namespace coordinates

// Add hash specialization in std namespace
namespace std {
    template<>
    struct hash<coordinates::BlockPosition> {
        std::size_t operator()(const coordinates::BlockPosition& pos) const {
            // Combine the hash of the three integers
            std::size_t h1 = std::hash<int32_t>{}(pos.blockId);
            std::size_t h2 = std::hash<int32_t>{}(pos.nucPos);
            std::size_t h3 = std::hash<int32_t>{}(pos.gapPos);
            
            // Simple combining function - can be improved if needed
            return ((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1);
        }
    };
} // namespace std
#ifndef COORDINATES_HPP
#define COORDINATES_HPP

#include "logging.hpp"
#include "seq_utils.hpp"
#include <chrono>
#include <map>
#include <memory>
#include <panmanUtils.hpp>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// Move kj-related includes and declarations outside the coordinates namespace
#include <capnp/common.h>
#include <kj/array.h>
#include <kj/memory.h>
#include <kj/string.h>

#include <tbb/scalable_allocator.h>

// Forward declare fixed_kmer functions needed by coordinates
namespace fixed_kmer {
int suggestOptimalKmerSize(int requestedSize);
}

namespace coordinates {

using namespace logging;

// Core types
typedef std::vector<
    std::pair<std::vector<std::pair<char, std::vector<char>>>,
              std::vector<std::vector<std::pair<char, std::vector<char>>>>>>
    sequence_t;

typedef std::vector<std::pair<bool, std::vector<bool>>> blockExists_t;
typedef std::vector<std::pair<bool, std::vector<bool>>> blockStrand_t;
typedef std::vector<std::tuple<int32_t, int32_t, bool, bool, bool, bool>>
    blockMutationInfo_t;

typedef std::vector<std::pair<
    std::vector<std::pair<int64_t, std::vector<int64_t>>>,
    std::vector<std::vector<std::pair<int64_t, std::vector<int64_t>>>>>>
    globalCoords_t;

// Gap-related types
using GapRange = std::pair<int64_t, int64_t>;
using GapUpdate = std::pair<bool, GapRange>;
using GapMap = std::map<int64_t, int64_t>;

// Forward declarations
class CoordinateManager;
class CoordinateTraverser;
struct tupleCoord_t;
struct CoordRange;

// First define tupleCoord_t
struct tupleCoord_t {
  int blockId;
  int nucPos;
  int nucGapPos;
  int64_t scalar;
  tupleCoord_t *next;
  tupleCoord_t *prev;

  constexpr tupleCoord_t() noexcept
      : blockId(-1), nucPos(-1), nucGapPos(-1), scalar(-1), next(nullptr),
        prev(nullptr) {}
  constexpr tupleCoord_t(int b, int n, int g, int64_t s = -1) noexcept
      : blockId(b), nucPos(n), nucGapPos(g), scalar(s), next(nullptr),
        prev(nullptr) {}

  // Core validation
  bool isValidBasic() const noexcept {
    // Basic range checks
    if (blockId < 0)
      return false;
    if (nucPos < 0)
      return false;
    if (nucGapPos < -1)
      return false; // -1 is valid for non-gap positions
    if (scalar < 0)
      return false;

    return true;
  }

  // Core comparison operators
  bool operator==(const tupleCoord_t &other) const noexcept {
    return scalar == other.scalar;
  }

  bool operator!=(const tupleCoord_t &other) const noexcept {
    return !(*this == other);
  }

  bool operator<(const tupleCoord_t &other) const noexcept {
    return scalar < other.scalar;
  }

  bool operator>(const tupleCoord_t &other) const noexcept {
    return other < *this;
  }

  bool operator<=(const tupleCoord_t &other) const noexcept {
    return !(other < *this);
  }

  bool operator>=(const tupleCoord_t &other) const noexcept {
    return !(*this < other);
  }

  // Simple pointer setup without scalar propagation
  void setNext(tupleCoord_t *n) noexcept {
    next = n;
    if (n) {
      n->prev = this;
    }
  }

  void setPrev(tupleCoord_t *p) noexcept {
    prev = p;
    if (p) {
      p->next = this;
    }
  }

  // Safe pointer access
  tupleCoord_t *getNext() const noexcept { return next; }

  tupleCoord_t *getPrev() const noexcept { return prev; }
};

// Define CoordinateManager before CoordinateTraverser
class CoordinateManager {
private:
  std::vector<tupleCoord_t> coords;
  std::vector<size_t> block_sizes;
  std::vector<std::pair<int64_t, int64_t>> block_ranges;
  std::vector<int64_t> scalarCoordToBlockId;
  std::map<int64_t, int64_t> gap_map;
  std::vector<int> BlockSizes;
  int64_t max_block_id = 0;
  int64_t max_nuc_pos = 0;
  int64_t max_gap_pos = 0;
  int64_t total_size = 0;
  int64_t max_scalar = 0;
  bool isCircular = false;  // Track if sequence is circular
  bool initialized = false; // Track if coordinates are initialized

  sequence_t *sequence;       // Non-owned pointer
  blockExists_t *blockExists; // Non-owned pointer
  blockStrand_t *blockStrand; // Non-owned pointer

  // Function declarations
  bool validateCoordinateState(const tupleCoord_t &coord,
                               const sequence_t &sequence) const noexcept;
  void setupCoordinates();

public:
  // Constructor now takes pointers to existing data
  CoordinateManager(sequence_t *seq, blockExists_t *exists,
                    blockStrand_t *strands)
      : sequence(seq), blockExists(exists), blockStrand(strands) {
    if (!seq) {
      throw std::runtime_error(
          "Null sequence pointer in CoordinateManager constructor");
    }

    // Initialize block existence and strand vectors
    blockExists->resize(seq->size(), {true, {}});
    blockStrand->resize(seq->size(), {true, {}});

    // Calculate maximum positions and total size in one pass
    total_size = 0;
    max_nuc_pos = 0;
    max_gap_pos = 0;
    max_block_id = seq->size() - 1;

    block_sizes.resize(seq->size(), 0);
    block_ranges.resize(seq->size(), {-1, -1});

    // msg("Scanning blocks for size calculation...");
    for (size_t blockId = 0; blockId < seq->size(); blockId++) {
      if (!seq->at(blockId).first.empty()) {
        size_t block_total = 0;
        for (const auto &nuc : seq->at(blockId).first) {
          block_total += nuc.second.size() + 1; // Count gaps + nucleotide
          block_sizes[blockId]++;
          if (!nuc.second.empty()) {
            max_gap_pos = std::max(max_gap_pos,
                                   static_cast<int64_t>(nuc.second.size() - 1));
          }
        }
        total_size += block_total;
        max_nuc_pos =
            std::max(max_nuc_pos,
                     static_cast<int64_t>(seq->at(blockId).first.size() - 1));
      }
    }

    // msg("Block scan complete:");
    // msg("  Total coordinates: {}", total_size);
    // msg("  Max block ID: {}", max_block_id);
    // msg("  Max nucleotide position: {}", max_nuc_pos);
    // msg("  Max gap position: {}", max_gap_pos);

    if (total_size <= 0) {
      throw std::runtime_error(
          "Invalid total_size calculated in CoordinateManager constructor");
    }

    // Initialize scalarCoordToBlockId with proper size
    scalarCoordToBlockId.resize(total_size);

    // Reserve space for coordinates to prevent reallocation
    coords.reserve(total_size);

    // Initialize coordinate vector and gap map
    gap_map[0] = total_size - 1;

    // Set up all coordinates
    setupCoordinates();

    // Mark as initialized
    initialized = true;
  }

  // Add a method to check if coordinates are initialized
  bool isInitialized() const { return initialized; }

  const int getNumBlocks() const { return sequence->size(); }
  const int getNumCoords() const {
    if (coords.empty()) {
      throw std::runtime_error("Coordinate vector is empty in getNumCoords()");
    }
    return coords.size();
  }

  const int64_t getScalarCoord(const tupleCoord_t &coord) const {
    return coord.scalar;
  }
  const int64_t getScalarCoord(const tupleCoord_t *coord) const {
    return coord->scalar;
  }

  const int64_t getBlockIdOfScalarCoord(int64_t scalar) const {
    return scalarCoordToBlockId[scalar];
  }
  const tupleCoord_t *getFirstCoordinateInBlock(int32_t blockId) const {
    if (!blockStrand->at(blockId).first) {
      // Inverted block - first coordinate is rightmost
      return &coords[block_ranges[blockId].second];
    }
    // Forward block - first coordinate is leftmost
    return &coords[block_ranges[blockId].first];
  }
  const tupleCoord_t *getLastCoordinateInBlock(int32_t blockId) const {
    if (!blockStrand->at(blockId).first) {
      // Inverted block - last coordinate is leftmost
      return &coords[block_ranges[blockId].first];
    }
    // Forward block - last coordinate is rightmost
    return &coords[block_ranges[blockId].second];
  }

  // Fast method to get the first coordinate in an "on" block
  const tupleCoord_t *getFirstCoordinateInOnBlock() const {
    const auto &blockExistVec = *blockExists;

    // Start from the lowest block ID and find the first "on" block
    for (int blockId = 0; blockId < block_ranges.size(); ++blockId) {
      if (blockExistVec[blockId].first) {
        // Found an "on" block, return its first coordinate
        return getFirstCoordinateInBlock(blockId);
      }
    }
    return nullptr; // No "on" blocks found
  }

  // Fast method to get the last coordinate in an "on" block
  // This avoids traversing the entire sequence
  const tupleCoord_t *getLastCoordinateInOnBlock() const {
    const auto &blockExistVec = *blockExists;

    // Start from the highest block ID and find the first "on" block
    for (int blockId = block_ranges.size() - 1; blockId >= 0; --blockId) {
      if (blockExistVec[blockId].first) {
        // Found an "on" block, return its last coordinate
        return getLastCoordinateInBlock(blockId);
      }
    }
    return nullptr; // No "on" blocks found
  }

  const tupleCoord_t *getLeftmostCoordinate() const { return &coords.front(); }

  const tupleCoord_t *getRightmostCoordinate() const {
    if (coords.size() < 2) {
      return nullptr;
    }
    return &coords[coords.size() - 1]; // skip sentinel
  }

  // Thread-safe block range updates
  void updateBlockRange(int blockId, int64_t start, int64_t end) {

    if (blockId < 0 || blockId >= block_ranges.size()) {
      err("Invalid block ID {} for range update (max: {})", blockId,
          block_ranges.size() - 1);
      return;
    }

    if (start < 0 || end < start || end >= coords.size()) {
      err("Invalid range [{},{}] for block {} (max scalar: {})", start, end,
          blockId, coords.size() - 1);
      return;
    }

    block_ranges[blockId] = {start, end};
  }

  // Thread-safe block range access with validation
  std::pair<int64_t, int64_t> getBlockRange(size_t blockId) const {
    if (blockId >= block_ranges.size())
      return {0, 0};
    return block_ranges[blockId];
  }

  const blockExists_t &getBlockExists() const { return *blockExists; }
  blockExists_t &getBlockExists() { return *blockExists; }

  const blockStrand_t &getBlockStrand() const { return *blockStrand; }
  blockStrand_t &getBlockStrand() { return *blockStrand; }

  const sequence_t &getSequence() const { return *sequence; }
  sequence_t &getSequence() { return *sequence; }

  void updateGapMapStep(
      const std::pair<bool, std::pair<int64_t, int64_t>> &update,
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &backtrack,
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapMapUpdates,
      bool recordGapMapUpdates = true) {
    if (update.first == true && update.second.first == update.second.second)
      return;

    bool toGap = update.first;
    int64_t start = update.second.first;
    int64_t end = update.second.second;

    // Safety check for invalid gap range
    if (start > end) {

      warn("Invalid gap range: start ({}) > end ({})", start, end);
      return;
    }

    auto rightIt = gap_map.upper_bound(start);
    auto leftIt =
        (rightIt == gap_map.begin()) ? gap_map.end() : std::prev(rightIt);

    bool rightItExists = rightIt != gap_map.end();
    bool leftItExists = leftIt != gap_map.end();

    if (toGap) {
      // add gap range
      if (gap_map.empty()) {
        gap_map[start] = end;
        backtrack.emplace_back(true, std::make_pair(start, end));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(start, end));
        }
        return;
      }

      decltype(rightIt) curIt;

      // curIt starts outside of any range
      if (!leftItExists || (!rightItExists && start > leftIt->second) ||
          (leftItExists && start > leftIt->second && rightItExists &&
           start < rightIt->first)) {
        if (leftItExists && start == leftIt->second + 1) {
          // 1 base after left range and merge with left
          curIt = leftIt;
          backtrack.emplace_back(false,
                                 std::make_pair(curIt->first, curIt->second));
          curIt->second = end;
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(
                false, std::make_pair(curIt->first, curIt->second));
          }
        } else {
          // insert new range
          auto tmpIt = gap_map.emplace(start, end);
          curIt = tmpIt.first;
          backtrack.emplace_back(true,
                                 std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(
                false, std::make_pair(curIt->first, curIt->second));
          }
        }
      } else {
        curIt = leftIt;
        if (end <= curIt->second) {
          return;
        }
        backtrack.emplace_back(false,
                               std::make_pair(curIt->first, curIt->second));
        curIt->second = end;
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(
              false, std::make_pair(curIt->first, curIt->second));
        }
      }

      auto nextIt = std::next(curIt);
      while (true) {
        if (nextIt == gap_map.end()) {
          break;
        }

        if (nextIt->second <= curIt->second) {
          auto tmpIt = nextIt;
          nextIt = std::next(nextIt);
          backtrack.emplace_back(false,
                                 std::make_pair(tmpIt->first, tmpIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(
                true, std::make_pair(tmpIt->first, tmpIt->second));
          }
          gap_map.erase(tmpIt);
        } else if (nextIt->first <= end + 1) {
          backtrack.emplace_back(false,
                                 std::make_pair(curIt->first, curIt->second));
          curIt->second = nextIt->second;
          backtrack.emplace_back(false,
                                 std::make_pair(nextIt->first, nextIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(
                false, std::make_pair(curIt->first, curIt->second));
            gapMapUpdates.emplace_back(
                true, std::make_pair(nextIt->first, nextIt->second));
          }
          gap_map.erase(nextIt);
          break;
        } else {
          break;
        }
      }
    } else {
      // remove gap range
      if (gap_map.empty() || (!leftItExists && end < rightIt->first) ||
          (!rightItExists && start > leftIt->second)) {
        return;
      }

      decltype(rightIt) curIt;
      decltype(rightIt) nextIt;
      if (!leftItExists || (leftItExists && start > leftIt->second &&
                            rightItExists && start < rightIt->first)) {
        // curIt starts outside of any range
        curIt = rightIt;

        if (end < curIt->first) {
          return;
        }

        // ends within the curIt range
        if (end <= curIt->second) {
          if (end == curIt->second) {
            backtrack.emplace_back(false,
                                   std::make_pair(curIt->first, curIt->second));
            if (recordGapMapUpdates) {
              gapMapUpdates.emplace_back(
                  true, std::make_pair(curIt->first, curIt->second));
            }
            gap_map.erase(curIt);
          } else {
            gap_map[end + 1] = curIt->second;
            backtrack.emplace_back(true,
                                   std::make_pair(end + 1, curIt->second));
            backtrack.emplace_back(false,
                                   std::make_pair(curIt->first, curIt->second));
            if (recordGapMapUpdates) {
              gapMapUpdates.emplace_back(
                  false, std::make_pair(end + 1, curIt->second));
              gapMapUpdates.emplace_back(
                  true, std::make_pair(curIt->first, curIt->second));
            }
            gap_map.erase(curIt);
          }
          return;
        } else {
          // curIt starts inside of a range
          curIt = leftIt;

          if (end <= curIt->second) {
            // contained in the curIt range
            if (start == curIt->first && end == curIt->second) {
              backtrack.emplace_back(
                  false, std::make_pair(curIt->first, curIt->second));
              if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(
                    true, std::make_pair(curIt->first, curIt->second));
              }
              gap_map.erase(curIt);
            } else if (start == curIt->first) {
              gap_map[end + 1] = curIt->second;
              backtrack.emplace_back(true,
                                     std::make_pair(end + 1, curIt->second));
              backtrack.emplace_back(
                  false, std::make_pair(curIt->first, curIt->second));
              if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(
                    false, std::make_pair(end + 1, curIt->second));
                gapMapUpdates.emplace_back(
                    true, std::make_pair(curIt->first, curIt->second));
              }
              gap_map.erase(curIt);
            } else if (end == curIt->second) {
              backtrack.emplace_back(
                  false, std::make_pair(curIt->first, curIt->second));
              curIt->second = start - 1;
              if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(
                    false, std::make_pair(curIt->first, curIt->second));
              }
            } else {
              gap_map[end + 1] = curIt->second;
              backtrack.emplace_back(true,
                                     std::make_pair(end + 1, curIt->second));
              backtrack.emplace_back(
                  false, std::make_pair(curIt->first, curIt->second));
              if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(
                    false, std::make_pair(end + 1, curIt->second));
                gapMapUpdates.emplace_back(
                    false, std::make_pair(curIt->first, start - 1));
              }
              curIt->second = start - 1;
            }
            return;
          } else {
            if (start == curIt->first) {
              nextIt = std::next(curIt);
              backtrack.emplace_back(
                  false, std::make_pair(curIt->first, curIt->second));
              if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(
                    true, std::make_pair(curIt->first, curIt->second));
              }
              gap_map.erase(curIt);
            } else {
              backtrack.emplace_back(
                  false, std::make_pair(curIt->first, curIt->second));
              curIt->second = start - 1;
              if (recordGapMapUpdates) {
                gapMapUpdates.emplace_back(
                    false, std::make_pair(curIt->first, curIt->second));
              }
              nextIt = std::next(curIt);
            }
          }
        }

        while (true) {
          if (nextIt == gap_map.end()) {
            break;
          }

          if (nextIt->first > end) {
            break;
          } else if (nextIt->second <= end) {
            auto tmpIt = nextIt;
            nextIt = std::next(nextIt);
            backtrack.emplace_back(false,
                                   std::make_pair(tmpIt->first, tmpIt->second));
            if (recordGapMapUpdates) {
              gapMapUpdates.emplace_back(
                  true, std::make_pair(tmpIt->first, tmpIt->second));
            }
            gap_map.erase(tmpIt);
          } else {
            gap_map[end + 1] = nextIt->second;
            backtrack.emplace_back(true,
                                   std::make_pair(end + 1, nextIt->second));
            backtrack.emplace_back(
                false, std::make_pair(nextIt->first, nextIt->second));
            if (recordGapMapUpdates) {
              gapMapUpdates.emplace_back(
                  false, std::make_pair(end + 1, nextIt->second));
              gapMapUpdates.emplace_back(
                  true, std::make_pair(nextIt->first, nextIt->second));
            }
            gap_map.erase(nextIt);
            break;
          }
        }
      }
    }
  }

  void updateGapMap(
      const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &updates,
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &backtrack,
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>
          &gapMapUpdates) {
    auto &gapMap = gap_map;
    for (const auto &update : updates) {
      updateGapMapStep(update, backtrack, gapMapUpdates);
    }
  }

  const std::map<int64_t, int64_t> &getGapMap() const { return gap_map; }

  void setGapMap(int64_t pos, int64_t val) { gap_map[pos] = val; }

  void removeFromGapMap(int64_t pos) { gap_map.erase(pos); }

  std::vector<std::pair<int64_t, int64_t>>
  invertRanges(const std::vector<std::pair<int64_t, int64_t>> &nucRanges,
               const std::pair<int64_t, int64_t> &invertRange) {
    std::vector<std::pair<int64_t, int64_t>> invertedRanges;

    auto [start, end] = invertRange;

    for (auto it = nucRanges.rbegin(); it != nucRanges.rend(); ++it) {
      const auto &[curStart, curEnd] = *it;
      invertedRanges.emplace_back(start + end - curEnd, start + end - curStart);
    }

    return invertedRanges;
  }

  void invertGapMap(
      const std::pair<int64_t, int64_t> &invertRange,
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &backtrack,
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>
          &gapMapUpdates) {
    auto &gapMap = gap_map;
    const auto &[start, end] = invertRange;

    // Validate range before proceeding
    if (start < 0 || end < start || end >= coords.size()) {
      warn("Invalid inversion range [{}:{}]. Skipping gap map inversion.",
           start, end);
      return;
    }

    // First, collect all existing gap runs within this range
    std::vector<std::pair<int64_t, int64_t>> existingGaps;
    for (auto it = gap_map.lower_bound(start);
         it != gap_map.end() && it->first <= end; ++it) {
      // Don't allow gaps to extend beyond the inversion range
      int64_t gapEnd = std::min(it->first + it->second - 1, end);
      if (gapEnd >= it->first) { // Sanity check
        existingGaps.push_back({it->first, gapEnd});
        // Record backtracking info for existing gap
        backtrack.push_back({false, {it->first, it->second}});
        // Mark this gap for removal
        gapMapUpdates.push_back({true, {it->first, it->second}});
      }
    }

    // Remove all existing gaps in this range
    for (const auto &[pos, _] : existingGaps) {
      gap_map.erase(pos);
    }

    // Now create inverted gaps
    for (const auto &[gapStart, gapEnd] : existingGaps) {
      // Calculate inverted positions
      int64_t invertedEnd = start + end - gapStart;
      int64_t invertedStart = start + end - gapEnd;

      // Ensure invertedStart <= invertedEnd
      if (invertedStart > invertedEnd) {
        std::swap(invertedStart, invertedEnd);
      }

      // Add new inverted gap
      int64_t length = invertedEnd - invertedStart + 1;
      gap_map[invertedStart] = length;

      // Record this new gap for backtracking
      backtrack.push_back({true, {invertedStart, 0}});
      gapMapUpdates.push_back({false, {invertedStart, length}});
    }

    // Log summary
    msg("Inverted {} gap runs in range [{}:{}]", existingGaps.size(), start,
        end);
  }

  void removeGapMapEntry(int64_t pos) { gap_map.erase(pos); }

  bool validateCoordinate(const tupleCoord_t &coord) const noexcept {

    // Essential range checks
    if (coord.blockId < 0 || coord.blockId >= sequence->size())
      return false;

    const auto &block = sequence->at(coord.blockId);
    if (coord.nucPos < 0 || coord.nucPos >= block.first.size())
      return false;

    // Gap position validation
    if (coord.nucGapPos >= 0) {
      if (coord.nucGapPos >= block.first[coord.nucPos].second.size())
        return false;
    }

    return true;
  }

  size_t size() const { return coords.size(); }

  // Update block pointers when existence or strand changes
  bool reinitializeBlockCoordinates(int32_t blockId) {
    // msg("Reinitializing block coordinates for blockId: {}", blockId);
    // Validate blockId
    if (blockId < 0 || blockId >= getNumBlocks()) {
      err("Invalid blockId: {} (valid range: 0-{})", blockId,
          getNumBlocks() - 1);
      return false;
    }

    // Check if block is inverted
    bool isInverted = !blockStrand->at(blockId).first;
    // msg("Block {} is{} inverted", blockId, isInverted ? "" : " not");

    // Get the range of this block
    int64_t start_scalar = block_ranges[blockId].first;
    int64_t end_scalar = block_ranges[blockId].second;

    if (start_scalar < 0 || end_scalar < 0 || start_scalar > end_scalar) {
      err("Invalid block range for blockId {}: [{}, {}]", blockId, start_scalar,
          end_scalar);
      return false;
    }

    bool isLiveBlock = blockExists->at(blockId).first;

    // First, reset all pointers in this block
    for (int64_t i = start_scalar; i <= end_scalar; i++) {
      coords[i].next = nullptr;
      coords[i].prev = nullptr;
    }

    if (isLiveBlock) {
      if (isInverted) {
        // For inverted blocks, connect in REVERSE scalar order (traversal goes
        // from high to low scalar) The actual traversal will be: end_scalar ->
        // end_scalar-1 -> ... -> start_scalar
        for (int64_t i = end_scalar; i > start_scalar; --i) {
          tupleCoord_t *current = &coords[i];
          tupleCoord_t *next = &coords[i - 1];
          current->next = next;
          next->prev = current;
        }

        // The first coordinate in traversal order for inverted blocks is at
        // end_scalar
        tupleCoord_t *first_in_traversal = &coords[end_scalar];

        // The last coordinate in traversal order for inverted blocks is at
        // start_scalar
        tupleCoord_t *last_in_traversal = &coords[start_scalar];

        // Find previous live block
        int prevLiveBlockId = -1;
        for (int i = blockId - 1; i >= 0; i--) {
          if (blockExists->at(i).first) {
            prevLiveBlockId = i;
            break;
          }
        }

        // Connect to previous live block if it exists, otherwise explicitly set
        // prev to nullptr
        if (prevLiveBlockId >= 0) {
          int64_t prev_block_end = block_ranges[prevLiveBlockId].second;
          tupleCoord_t *prev_block_last = &coords[prev_block_end];

          first_in_traversal->prev = prev_block_last;
          prev_block_last->next = first_in_traversal;
        } else {
          first_in_traversal->prev = nullptr;
        }

        // Find next live block
        int nextLiveBlockId = -1;
        for (int i = blockId + 1; i < getNumBlocks(); i++) {
          if (blockExists->at(i).first) {
            nextLiveBlockId = i;
            break;
          }
        }

        // Connect to next live block if it exists, otherwise explicitly set
        // next to nullptr
        if (nextLiveBlockId >= 0) {
          int64_t next_block_start = block_ranges[nextLiveBlockId].first;
          tupleCoord_t *next_block_first = &coords[next_block_start];

          last_in_traversal->next = next_block_first;
          next_block_first->prev = last_in_traversal;
        } else {
          last_in_traversal->next = nullptr;
        }
      } else {
        // For forward blocks, connect in FORWARD scalar order (traversal
        // matches scalar order) The actual traversal will be: start_scalar ->
        // start_scalar+1 -> ... -> end_scalar
        for (int64_t i = start_scalar; i < end_scalar; ++i) {
          tupleCoord_t *current = &coords[i];
          tupleCoord_t *next = &coords[i + 1];
          current->next = next;
          next->prev = current;
        }

        // The first coordinate in traversal order for forward blocks is at
        // start_scalar
        tupleCoord_t *first_in_traversal = &coords[start_scalar];

        // The last coordinate in traversal order for forward blocks is at
        // end_scalar
        tupleCoord_t *last_in_traversal = &coords[end_scalar];

        // Find previous live block
        int prevLiveBlockId = -1;
        for (int i = blockId - 1; i >= 0; i--) {
          if (blockExists->at(i).first) {
            prevLiveBlockId = i;
            break;
          }
        }

        // Connect to previous live block if it exists, otherwise explicitly set
        // prev to nullptr
        if (prevLiveBlockId >= 0) {
          int64_t prev_block_end = block_ranges[prevLiveBlockId].second;
          tupleCoord_t *prev_block_last = &coords[prev_block_end];

          first_in_traversal->prev = prev_block_last;
          prev_block_last->next = first_in_traversal;
        } else {
          first_in_traversal->prev = nullptr;
        }

        // Find next live block
        int nextLiveBlockId = -1;
        for (int i = blockId + 1; i < getNumBlocks(); i++) {
          if (blockExists->at(i).first) {
            nextLiveBlockId = i;
            break;
          }
        }

        // Connect to next live block if it exists, otherwise explicitly set
        // next to nullptr
        if (nextLiveBlockId >= 0) {
          int64_t next_block_start = block_ranges[nextLiveBlockId].first;
          tupleCoord_t *next_block_first = &coords[next_block_start];

          last_in_traversal->next = next_block_first;
          next_block_first->prev = last_in_traversal;
        } else {
          last_in_traversal->next = nullptr;
        }
      }
    }
    // We don't do anything for dead blocks - they're simply skipped in the
    // traversal

    // For dead blocks, connect the previous live block to the next live block
    if (!isLiveBlock) {
      // Find previous live block
      int prevLiveBlockId = -1;
      for (int i = blockId - 1; i >= 0; i--) {
        if (blockExists->at(i).first) {
          prevLiveBlockId = i;
          break;
        }
      }

      // Find next live block
      int nextLiveBlockId = -1;
      for (int i = blockId + 1; i < getNumBlocks(); i++) {
        if (blockExists->at(i).first) {
          nextLiveBlockId = i;
          break;
        }
      }

      // Connect previous live block to next live block if both exist
      if (prevLiveBlockId >= 0 && nextLiveBlockId >= 0) {
        int64_t prev_block_end = block_ranges[prevLiveBlockId].second;
        int64_t next_block_start = block_ranges[nextLiveBlockId].first;

        tupleCoord_t *prev_block_last = &coords[prev_block_end];
        tupleCoord_t *next_block_first = &coords[next_block_start];

        prev_block_last->next = next_block_first;
        next_block_first->prev = prev_block_last;

        // Add gap information for the dead block
        // Calculate the size of the dead block
        int64_t dead_block_size = end_scalar - start_scalar + 1;

        // Create a gap entry in the gap map to represent the dead block
        // The gap starts at the end of the previous live block and extends to
        // the start of the next live block
        std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> backtrack;
        std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapMapUpdates;

        // Add a gap from prev_block_end to next_block_start
        // The first parameter is true to indicate adding a gap (not removing)
        // The second parameter is the range of the gap
        // Safety check to ensure valid gap range
        if (prev_block_end < next_block_start - 1) {
          updateGapMapStep({true, {prev_block_end + 1, next_block_start - 1}},
                           backtrack, gapMapUpdates);
        } else {
          // Log warning if gap range is invalid
          warn("Invalid gap range for dead block {}: [{}, {}]", blockId,
               prev_block_end + 1, next_block_start - 1);
        }
      }
    }

    // After reinitialization, let's verify the block's traversal
    if (isLiveBlock) {
      const tupleCoord_t *firstCoord = getFirstCoordinateInBlock(blockId);
      if (firstCoord) {
        std::string traversalSeq;
        const tupleCoord_t *curr = firstCoord;
        int maxSteps = 50; // Limit to avoid potential infinite loops
        int step = 0;

        while (curr && curr->blockId == blockId && step < maxSteps) {
          char c;
          if (curr->nucGapPos == -1) {
            c = sequence->at(curr->blockId).first[curr->nucPos].first;
          } else {
            c = sequence->at(curr->blockId)
                    .first[curr->nucPos]
                    .second[curr->nucGapPos];
          }

          // Apply complementation for inverted blocks
          if (isInverted && c != '-' && c != 'x') {
            c = seq_utils::getComplementCharacter(c);
          }

          traversalSeq += c;
          curr = curr->next;
          step++;
        }

        // msg("Block {} traversal after reinitialization: '{}'", blockId,
        // traversalSeq);
      }
    }

    return true;
  }

  // Update block strand state
  void updateBlockStrand(int32_t blockId, bool isForward) {
    if (blockId < 0 || blockId >= blockStrand->size()) {
      err("Invalid block ID for strand update: {}", blockId);
      return;
    }
    blockStrand->at(blockId).first = isForward;
  }

  // Update block exists state
  void updateBlockExists(int32_t blockId, bool exists) {
    if (blockId < 0 || blockId >= blockExists->size()) {
      err("Invalid block ID for exists update: {}", blockId);
      return;
    }
    blockExists->at(blockId).first = exists;
  }

  void updateSequence(int32_t blockId, int32_t nucPos, int32_t nucGapPos,
                      const char &nuc) {
    if (blockId < 0 || blockId >= sequence->size()) {
      err("Invalid block ID for sequence update: {}", blockId);
      return;
    }
    if (nucGapPos == -1) {
      sequence->at(blockId).first[nucPos].first = nuc;
    } else {
      sequence->at(blockId).first[nucPos].second[nucGapPos] = nuc;
    }
  }

  // Compare coordinates respecting traversal order (inversion-aware)
  bool compareCoordsByTraversalOrder(const tupleCoord_t *a,
                                     const tupleCoord_t *b) const {
    // Handle null pointers
    if (!a || !b)
      return false;

    // If coordinates are in different blocks
    if (a->blockId != b->blockId) {
      return a->blockId < b->blockId;
    }

    // For coordinates in the same block
    if (!blockStrand->at(a->blockId).first) {
      // Block is inverted, reverse comparison
      return a->scalar > b->scalar;
    }

    // Normal forward comparison
    return a->scalar < b->scalar;
  }

  // Fast vector lookup of tupleCoord by scalar
  const tupleCoord_t *getTupleCoord(int64_t scalar) const {
    if (scalar < 0 || scalar >= coords.size())
      return nullptr;
    return &coords[scalar];
  }

  // Find a coordinate based on block ID, nucleotide position, and gap position
  // using std::lower_bound
  const tupleCoord_t *findTupleCoord(const tupleCoord_t &coord) const {
    // Validate input coordinate

    // Get block range
    auto range = getBlockRange(coord.blockId);
    if (range.first < 0 || range.second < range.first ||
        range.second >= coords.size()) {
      err("findTupleCoord: invalid block range [{},{}] for block {}",
          range.first, range.second, coord.blockId);
      return nullptr;
    }

    // Define comparator for coordinate ordering
    auto comp = [](const tupleCoord_t &a, const tupleCoord_t &b) {
      if (a.blockId != b.blockId)
        return a.blockId < b.blockId;
      if (a.nucPos != b.nucPos)
        return a.nucPos < b.nucPos;
      if (a.nucGapPos != b.nucGapPos) {
        if (a.nucGapPos == -1)
          return false;
        if (b.nucGapPos == -1)
          return true;

        return a.nucGapPos < b.nucGapPos;
      }
      return false;
    };

    // Use lower_bound to find first element not less than coord
    auto it = std::lower_bound(coords.begin() + range.first,
                               coords.begin() + range.second + 1, coord, comp);

    return &(*it);
  }

  // Gap map operations - enhanced versions

  // Get the next non-gap position after the provided scalar coordinate
  const tupleCoord_t *getNextNonGapPosition(int64_t scalar) const {
    if (scalar < 0 || scalar >= coords.size()) {
      return nullptr;
    }

    // Check if this position is in a gap run
    auto it = gap_map.find(scalar);
    if (it != gap_map.end()) {
      // Found a gap run starting at this position
      // Skip ahead by the gap length
      int64_t next_pos = scalar + it->second;

      // Validate the next position is actually after the current one
      if (next_pos <= scalar) {
        err("getNextNonGapPosition: invalid gap map entry at position {} - gap "
            "length {} is not positive",
            scalar, it->second);
        return nullptr;
      }

      if (next_pos < coords.size()) {
        return &coords[next_pos];
      } else {
        return nullptr; // Beyond sequence end
      }
    }

    // If not the start of a gap run, check if we're in the middle of one
    auto upper_it = gap_map.upper_bound(scalar);
    if (upper_it != gap_map.begin()) {
      auto lower_it = std::prev(upper_it);
      if (scalar >= lower_it->first &&
          scalar < lower_it->first + lower_it->second) {
        // We're in the middle of a gap run
        // Skip to the end of this run
        int64_t next_pos = lower_it->first + lower_it->second;

        // Validate the next position is actually after the current one
        if (next_pos <= scalar) {
          err("getNextNonGapPosition: invalid gap map entry at position {} - "
              "calculated next position {} is not after current position",
              scalar, next_pos);
          return nullptr;
        }

        if (next_pos < coords.size()) {
          return &coords[next_pos];
        } else {
          return nullptr; // Beyond sequence end
        }
      }
    }

    // If we're not in any gap run, return the next position
    if (scalar + 1 < coords.size()) {
      return &coords[scalar + 1];
    }

    return nullptr; // No next position
  }

  // Get the previous non-gap position before the provided scalar coordinate
  const tupleCoord_t *getPreviousNonGapPosition(int64_t scalar) const {
    if (scalar <= 0 || scalar >= coords.size()) {
      return nullptr;
    }

    // Start by checking the previous position
    int64_t prev_pos = scalar - 1;

    // Check if that position is in a gap run
    if (isGapPosition(prev_pos)) {
      // We're in a gap run, need to find its beginning
      auto upper_it = gap_map.upper_bound(prev_pos);

      // If we find a relevant gap run entry, we'll use it to calculate
      if (upper_it != gap_map.begin()) {
        auto lower_it = std::prev(upper_it);

        // Check if we're in this gap run
        if (prev_pos >= lower_it->first &&
            prev_pos < lower_it->first + lower_it->second) {
          // We're in this gap run, jump to position before it starts
          prev_pos = lower_it->first - 1;

          // If that's still valid
          if (prev_pos >= 0) {
            return &coords[prev_pos];
          } else {
            return nullptr; // No previous position
          }
        }
      }

      // If we didn't find the gap run entry (or it wasn't relevant),
      // we'll search backward position by position, but this should be rare
      while (prev_pos >= 0 && isGapPosition(prev_pos)) {
        prev_pos--;
      }

      if (prev_pos >= 0) {
        return &coords[prev_pos];
      } else {
        return nullptr; // No previous position
      }
    } else {
      // Previous position isn't a gap, return it directly
      return &coords[prev_pos];
    }
  }

  // Check if a position is within a gap
  bool isGapPosition(int64_t scalar) const {
    if (scalar < 0 || scalar >= coords.size()) {
      return false;
    }

    // Check if it's the start of a gap run
    if (gap_map.find(scalar) != gap_map.end()) {
      return true;
    }

    // Check if it's in the middle of a gap run
    auto it = gap_map.upper_bound(scalar);
    if (it != gap_map.begin()) {
      auto prev_it = std::prev(it);
      if (scalar >= prev_it->first &&
          scalar < prev_it->first + prev_it->second) {
        return true;
      }
    }

    return false;
  }

  // Get the gap run length at a position
  int64_t getGapRunLength(int64_t scalar) const {
    if (scalar < 0 || scalar >= coords.size()) {
      return 0;
    }

    // Check if it's the start of a gap run
    auto it = gap_map.find(scalar);
    if (it != gap_map.end()) {
      return it->second;
    }

    // Check if it's in the middle of a gap run
    auto upper_it = gap_map.upper_bound(scalar);
    if (upper_it != gap_map.begin()) {
      auto lower_it = std::prev(upper_it);
      if (scalar >= lower_it->first &&
          scalar < lower_it->first + lower_it->second) {
        // We're in the middle of a gap run
        return lower_it->first + lower_it->second - scalar;
      }
    }

    return 0; // Not in a gap run
  }

  // Optimize gap map by merging adjacent runs and fixing issues
  void optimizeGapMap() {
    if (gap_map.empty()) {
      return;
    }

    std::map<int64_t, int64_t> optimized_gap_map;
    std::vector<std::pair<int64_t, int64_t>> to_process;

    // First pass: collect valid entries and fix out-of-bounds issues
    for (const auto &[pos, length] : gap_map) {
      // Skip negative positions
      if (pos < 0) {
        err("optimizeGapMap: Skipping negative position entry at {}", pos);
        continue;
      }

      // Skip non-positive lengths
      if (length <= 0) {
        err("optimizeGapMap: Skipping non-positive length entry at position {}",
            length);
        continue;
      }

      // Validate end position doesn't exceed sequence length
      int64_t end_pos = pos + length - 1; // Inclusive end position
      if (end_pos >= coords.size()) {
        err("optimizeGapMap: Truncating out-of-bounds entry at position {} "
            "from length {} to {}",
            pos, length, coords.size() - pos);
        to_process.push_back({pos, coords.size() - pos});
      } else {
        to_process.push_back({pos, length});
      }
    }

    // Log the number of collected entries
    msg("optimizeGapMap: Processing {} gap runs", to_process.size());

    // Sort entries by position for easier processing
    std::sort(to_process.begin(), to_process.end());

    // Merge overlapping or adjacent runs
    if (!to_process.empty()) {
      int64_t current_start = to_process[0].first;
      int64_t current_end = to_process[0].first + to_process[0].second -
                            1; // Convert to inclusive end

      for (size_t i = 1; i < to_process.size(); ++i) {
        int64_t next_start = to_process[i].first;
        int64_t next_end = to_process[i].first + to_process[i].second -
                           1; // Convert to inclusive end

        // Check if runs overlap or are adjacent
        if (next_start <= current_end + 1) {
          // Found overlapping or adjacent run - extend current run if needed
          current_end = std::max(current_end, next_end);
          msg("optimizeGapMap: Merging overlapping/adjacent gap runs at {} and "
              "{}",
              current_start, next_start);
        } else {
          // No overlap - add current run to optimized map and start a new one
          optimized_gap_map[current_start] =
              current_end - current_start + 1; // Convert back to length
          current_start = next_start;
          current_end = next_end;
        }
      }

      // Add the last run
      optimized_gap_map[current_start] =
          current_end - current_start + 1; // Convert back to length
    }

    // Final validation - check for any remaining overlaps
    bool has_overlaps = false;
    for (auto it = optimized_gap_map.begin(); it != optimized_gap_map.end();
         ++it) {
      auto next_it = std::next(it);
      if (next_it != optimized_gap_map.end()) {
        int64_t current_end = it->first + it->second - 1;
        if (current_end >= next_it->first) {
          err("optimizeGapMap: CRITICAL - Still have overlap between entries "
              "at {} (length {}) and {} (length {})",
              it->first, it->second, next_it->first, next_it->second);
          has_overlaps = true;
        }
      }
    }

    if (has_overlaps) {
      err("optimizeGapMap: Failed to eliminate all overlaps. Gap map may still "
          "be corrupted.");
    }

    // Replace the original gap map with the optimized one
    msg("Gap map optimization: reduced from {} to {} entries", gap_map.size(),
        optimized_gap_map.size());
    gap_map = std::move(optimized_gap_map);
  }

  // Get gap count between two positions
  int64_t getGapCount(int64_t start_scalar, int64_t end_scalar) const {
    if (start_scalar >= end_scalar || start_scalar < 0 ||
        end_scalar >= coords.size()) {
      return 0;
    }

    int64_t total_gaps = 0;

    // Find all gap runs that intersect with our range
    auto it = gap_map.lower_bound(start_scalar);

    // If we're not at the first gap run, check if we're inside a previous one
    if (it != gap_map.begin()) {
      auto prev_it = std::prev(it);
      if (start_scalar >= prev_it->first &&
          start_scalar < prev_it->first + prev_it->second) {
        // We start inside a gap run
        int64_t run_end = prev_it->first + prev_it->second;
        total_gaps += std::min(run_end, end_scalar) - start_scalar;

        // If the run covers our entire range, we're done
        if (run_end >= end_scalar) {
          return total_gaps;
        }
      }
    }

    // Process gap runs that start within our range
    while (it != gap_map.end() && it->first < end_scalar) {
      int64_t run_start = std::max(it->first, start_scalar);
      int64_t run_end = std::min(it->first + it->second, end_scalar);

      if (run_end > run_start) {
        total_gaps += run_end - run_start;
      }

      ++it;
    }

    return total_gaps;
  }

  // Restore gap map from backtracking information
  void restoreGapMapFromBacktracks(
      const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>
          &backtrackInfo) {
    msg("Restoring gap map with {} backtrack entries", backtrackInfo.size());

    // First pass: Identify and validate all entries that will be part of the
    // gap map after restoration
    std::map<int64_t, int64_t>
        targetGapMap; // What the gap map should look like after restoration
    std::vector<int64_t>
        positionsToRemove; // Positions that should be removed from gap map

    // Start with current gap map
    targetGapMap = gap_map;

    // Apply backtracking in reverse order
    for (auto it = backtrackInfo.rbegin(); it != backtrackInfo.rend(); ++it) {
      const auto &[is_new_entry, range] = *it;
      const auto &[pos, value] = range;

      // Basic validation of the backtrack entry
      if (pos < 0) {
        warn("Invalid position {} in backtrack entry", pos);
        continue;
      }

      if (is_new_entry) {
        // This was a new entry added during the mutation, so remove it
        positionsToRemove.push_back(pos);
        targetGapMap.erase(pos);
      } else {
        // This was an existing entry that was modified, restore its original
        // value
        if (value <= 0) {
          warn("Invalid value {} for position {} in backtrack entry", value,
               pos);
          positionsToRemove.push_back(pos);
          targetGapMap.erase(pos);
        } else if (pos + value > coords.size()) {
          warn(
              "Entry at position {} with length {} exceeds coordinate count {}",
              pos, value, coords.size());
          // Keep entry but trim to valid range
          targetGapMap[pos] = coords.size() - pos;
        } else {
          // Valid entry, restore it
          targetGapMap[pos] = value;
        }
      }
    }

    // Second pass: Check for overlapping entries and resolve conflicts
    bool has_overlaps = false;
    for (auto it = targetGapMap.begin(); it != targetGapMap.end(); ++it) {
      auto next_it = std::next(it);
      if (next_it != targetGapMap.end()) {
        int64_t curr_end = it->first + it->second - 1;
        if (curr_end >= next_it->first) {
          has_overlaps = true;
          warn("Detected overlap between gap entries at positions {} and {}",
               it->first, next_it->first);
          break;
        }
      }
    }

    // If overlaps are detected or there are out-of-bounds entries, run
    // optimization
    if (has_overlaps) {
      warn("Overlapping entries detected in target gap map, optimizing");
      // Create temp map with proper boundaries
      std::map<int64_t, int64_t> optimizedMap;
      std::vector<std::pair<int64_t, int64_t>> entries;

      // Extract valid entries
      for (const auto &[pos, len] : targetGapMap) {
        if (pos >= 0 && len > 0 && pos + len <= coords.size()) {
          entries.push_back({pos, len});
        }
      }

      // Sort by position
      std::sort(entries.begin(), entries.end());

      // Merge overlapping or adjacent entries
      if (!entries.empty()) {
        int64_t curr_start = entries[0].first;
        int64_t curr_end = curr_start + entries[0].second - 1;

        for (size_t i = 1; i < entries.size(); ++i) {
          int64_t next_start = entries[i].first;
          int64_t next_end = next_start + entries[i].second - 1;

          if (next_start <= curr_end + 1) {
            // Merge with current entry
            curr_end = std::max(curr_end, next_end);
          } else {
            // Add current entry and start a new one
            optimizedMap[curr_start] = curr_end - curr_start + 1;
            curr_start = next_start;
            curr_end = next_end;
          }
        }

        // Add the last entry
        optimizedMap[curr_start] = curr_end - curr_start + 1;
      }

      // Use optimized map
      targetGapMap = optimizedMap;
    }

    // Apply changes to actual gap map
    // First, remove entries that should be removed
    for (int64_t pos : positionsToRemove) {
      gap_map.erase(pos);
    }

    // Then update all entries in the target map
    for (const auto &[pos, len] : targetGapMap) {
      // Skip invalid entries as a final safety check
      if (pos < 0 || len <= 0 || pos + len > coords.size()) {
        warn("Skipping invalid entry: pos={}, len={}", pos, len);
        continue;
      }

      // Update or add entry
      gap_map[pos] = len;
    }

    // Final validation
    for (const auto &[pos, len] : gap_map) {
      if (pos < 0 || len <= 0 || pos + len > coords.size()) {
        warn("Invalid entry after restoration: pos={}, len={}", pos, len);
      }
    }

    msg("Gap map restoration complete, final size: {}", gap_map.size());
  }

  // Method to update or remove a gap map entry
  void updateGapMapEntry(int64_t pos, int64_t val, bool erase) {
    if (erase) {
      gap_map.erase(pos);
    } else {
      gap_map[pos] = val;
    }
  }

  // Method to clear the entire gap map
  void clearGapMap() { gap_map.clear(); }

  // Validate the integrity of the gap map
  bool validateGapMap(bool fix_issues = false) const {
    bool isValid = true;
    std::vector<std::pair<int64_t, int64_t>> overlappingEntries;
    std::vector<std::pair<int64_t, int64_t>> outOfBoundsEntries;
    std::vector<std::pair<int64_t, int64_t>> negativeEntries;
    std::vector<std::pair<int64_t, int64_t>> invertedRangeEntries;
    std::vector<std::pair<int64_t, int64_t>> suspiciouslyLargeEntries;

    // Calculate average and maximum block sizes for reference
    int64_t totalBlockSize = 0;
    int64_t maxBlockSize = 0;
    int blockCount = 0;

    for (size_t i = 0; i < block_ranges.size(); i++) {
      if (i < blockExists->size() && blockExists->at(i).first) {
        // This is a live block
        int64_t size = block_ranges[i].second - block_ranges[i].first + 1;
        totalBlockSize += size;
        maxBlockSize = std::max(maxBlockSize, size);
        blockCount++;
      }
    }

    // Calculate average block size and reasonable maximum gap size
    int64_t avgBlockSize = blockCount > 0 ? totalBlockSize / blockCount : 1000;
    int64_t reasonableMaxGapSize =
        std::max(maxBlockSize * 3, avgBlockSize * 10);

    // Check each gap map entry
    for (const auto &[pos, length] : gap_map) {
      // Validate position is in range
      if (pos < 0 || pos >= coords.size()) {
        err("validateGapMap: gap map entry at position {} is out of range "
            "[0,{}]",
            pos, coords.size() - 1);
        outOfBoundsEntries.push_back({pos, length});
        isValid = false;
        continue;
      }

      // Validate length is positive
      if (length <= 0) {
        err("validateGapMap: gap map entry at position {} has invalid length "
            "{}",
            pos, length);
        negativeEntries.push_back({pos, length});
        isValid = false;
        continue;
      }

      // Validate end position is in range
      int64_t end_pos = pos + length - 1; // Inclusive end position
      if (end_pos >= coords.size()) {
        err("validateGapMap: gap map entry at position {} with length {} "
            "extends beyond coordinate count {}",
            pos, length, coords.size());
        outOfBoundsEntries.push_back({pos, length});
        isValid = false;
        continue;
      }

      // Validate there are no overlapping gap runs
      auto next_it = gap_map.upper_bound(pos);
      if (next_it != gap_map.end() && next_it->first <= end_pos) {
        err("validateGapMap: gap map entry at position {} with length {} "
            "overlaps with entry at position {}",
            pos, length, next_it->first);
        overlappingEntries.push_back({pos, length});
        isValid = false;
      }
    }

    // If issues were found, print detailed diagnostics
    if (!isValid) {
      std::cout << "Gap map validation failed with "
                << overlappingEntries.size() << " overlapping entries, "
                << outOfBoundsEntries.size() << " out-of-bounds entries, "
                << negativeEntries.size() << " negative length entries, and "
                << invertedRangeEntries.size() << " inverted range entries."
                << std::endl;

      if (!overlappingEntries.empty()) {
        err("Overlapping entries:");
        int count = 0;
        for (const auto &[pos, length] : overlappingEntries) {
          if (count++ >= 5) {
            err("  ... and {} more", overlappingEntries.size() - 5);
            break;
          }
          err("  Position: {}, Length: {}, End: {}", pos, length,
              pos + length - 1);

          // Get block information for this position
          int64_t blockId = pos < scalarCoordToBlockId.size()
                                ? scalarCoordToBlockId[pos]
                                : -1;
          if (blockId >= 0 && blockId < blockExists->size()) {
            bool blockIsLive = blockExists->at(blockId).first;
            err("    In block {}, which is {}", blockId,
                blockIsLive ? "LIVE" : "DEAD");
          }
        }
      }

      // If fix_issues is true, try to optimize the gap map
      if (fix_issues) {
        const_cast<CoordinateManager *>(this)->optimizeGapMap();
      }
    }

    return isValid;
  }

  ~CoordinateManager() {}
};

class CoordinateTraverser {
private:
  const tupleCoord_t *current;
  const tupleCoord_t *end;
  CoordinateManager *coordManager;
  bool done;

public:
  CoordinateTraverser(const tupleCoord_t *start, CoordinateManager *manager);
  CoordinateTraverser();
  CoordinateTraverser(const CoordinateTraverser &other) = default;
  CoordinateTraverser &operator=(const CoordinateTraverser &other) = default;

  bool validateState() const noexcept;
  bool isDone() const noexcept { return done; }
  void next() noexcept;
  void prev() noexcept;
  void skipGaps() noexcept;
  void skipGapsBackward() noexcept;
  bool skipToNthNonGap(int n, char *buffer) noexcept;
  void nextSkipGaps() noexcept;
  void reset(const tupleCoord_t *start) noexcept;
  void reset(const tupleCoord_t *start, const tupleCoord_t *end) noexcept;
  bool isInverted() const noexcept;
  char getChar(bool complement = true) const noexcept;
  void moveToNextBlock() noexcept;
  void moveToPreviousBlock() noexcept;
  bool isBlockInverted(int blockId) const noexcept;
  bool isValidNucleotide(char c) const noexcept;

  const tupleCoord_t *getCurrent() const noexcept { return current; }
  const tupleCoord_t *getEnd() const noexcept { return end; }
  CoordinateManager &getCoordManager() const noexcept { return *coordManager; }

  // Thread-safe aligned memory allocation helper (for SIMD operations)
  template <typename T>
  T *allocateAligned(size_t count, size_t alignment = 64) const noexcept {
    return static_cast<T *>(
        scalable_aligned_malloc(count * sizeof(T), alignment));
  }

  // Thread-safe aligned memory deallocation
  template <typename T> void freeAligned(T *ptr) const noexcept {
    scalable_aligned_free(ptr);
  }
};

// Now define CoordRange after all dependencies are defined
struct CoordRange {
  int64_t start_scalar;
  int64_t stop_scalar;

  CoordRange() : start_scalar(-1), stop_scalar(-1) {}
  CoordRange(const tupleCoord_t *s, const tupleCoord_t *e)
      : start_scalar(s ? s->scalar : -1), stop_scalar(e ? e->scalar : -1) {
    if (s) {
      debug("CoordRange: Creating range with start blockId={}, nucPos={}, "
            "nucGapPos={}, scalar={}",
            s->blockId, s->nucPos, s->nucGapPos, s->scalar);
    }
    if (e) {
      debug("CoordRange: Creating range with end blockId={}, nucPos={}, "
            "nucGapPos={}, scalar={}",
            e->blockId, e->nucPos, e->nucGapPos, e->scalar);
    }
    normalize();
  }

  // Ensure start_scalar <= stop_scalar
  void normalize() {
    if (start_scalar >= 0 && stop_scalar >= 0 && start_scalar > stop_scalar) {

      // Normal inverted range case - swap values
      if (start_scalar > stop_scalar) {
        warn("Normalizing inverted range: start ({}) > stop ({})", start_scalar,
             stop_scalar);
        std::swap(start_scalar, stop_scalar);
        debug("AFTER Normalization: start={}, stop={}", start_scalar,
              stop_scalar);
      }
    }
  }

  bool isValid() const {
    return start_scalar >= 0 && stop_scalar >= 0 && start_scalar <= stop_scalar;
  }

  const tupleCoord_t *getStart(const CoordinateTraverser &traverser) const {
    if (!isValid())
      return nullptr;
    return traverser.getCoordManager().getTupleCoord(start_scalar);
  }

  const tupleCoord_t *getStop(const CoordinateTraverser &traverser) const {
    if (!isValid())
      return nullptr;
    return traverser.getCoordManager().getTupleCoord(stop_scalar);
  }

  void setStart(const int64_t startScalar) {
    start_scalar = startScalar;
    normalize();
  }

  void setStop(const int64_t stopScalar) {
    stop_scalar = stopScalar;
    normalize();
  }

  // Get length of the range (number of positions)
  int64_t length() const {
    if (!isValid())
      return 0;
    return stop_scalar - start_scalar + 1;
  }

  CoordRange merge(const CoordRange &other) const {
    // Create a new range that encompasses both this range and the other range
    if (!isValid())
      return other;
    if (!other.isValid())
      return *this;

    // Take the minimum of the start positions and maximum of the end positions
    CoordRange result;
    result.start_scalar = std::min(start_scalar, other.start_scalar);
    result.stop_scalar = std::max(stop_scalar, other.stop_scalar);
    result.normalize();

    return result;
  }

  // Note: This operator< should ONLY be used for storage/sorting when traversal
  // order doesn't matter For traversal-aware comparisons, use
  // compareCoordsByTraversalOrder through the CoordinateManager
  bool operator<(const CoordRange &other) const noexcept {
    if (!isValid() || !other.isValid())
      return false;
    // Compare by scalar since that's fixed and unique
    return start_scalar < other.start_scalar;
  }
};

// Go upstream (against traversal direction) until neededNongap nucleotides are
// seen
const tupleCoord_t *expandUpstream(CoordinateTraverser &traverser,
                                   const tupleCoord_t *coord, int neededNongap,
                                   blockExists_t &blockExists,
                                   blockStrand_t &blockStrand,
                                   const tupleCoord_t *stop_coord);

// Go downstream (following traversal direction) until neededNongap nucleotides
// are seen
const tupleCoord_t *expandDownstream(CoordinateTraverser &traverser,
                                     const tupleCoord_t *coord,
                                     int neededNongap,
                                     blockExists_t &blockExists,
                                     blockStrand_t &blockStrand);

void expandRanges(std::vector<CoordRange> &ranges,
                  CoordinateTraverser &traverser, int k);

void mergeRanges(std::vector<CoordRange> &ranges,
                 CoordinateTraverser &traverser);

} // namespace coordinates

#endif // COORDINATES_HPP

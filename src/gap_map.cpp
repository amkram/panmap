#include "gap_map.hpp"
#include "coordinates.hpp"
#include "logging.hpp"
#include <algorithm>
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_set.hpp>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <utility>
#include <vector>

using namespace boost::icl;

namespace gap_map {

void updateGapMapStep(GapMap &gapMap, const GapUpdate &update,
                      std::vector<GapUpdate> &backtrack,
                      std::vector<GapUpdate> &gapMapUpdates,
                      bool recordGapMapUpdates, int64_t totalCoordinateCount) {

  const auto &[del, range] = update;
  const auto &[start, end] = range;

  // Validate the update has a valid range
  if (start < 0 || end < start) {
    std::cerr << "Warning: Invalid gap update range: [" << start << ", " << end
              << "]" << std::endl;
    return;
  }

  // Always validate coordinates, totalCoordinateCount == 0 means no validation
  if (totalCoordinateCount > 0) {
    if (start < 0 || end < 0 || start >= totalCoordinateCount ||
        end >= totalCoordinateCount) {
      std::cerr << "Warning: Gap update (" << start << "," << end
                << ") out of bounds (0," << totalCoordinateCount - 1 << ")"
                << std::endl;
      // Skip or truncate as appropriate
      if (start < 0 || start >= totalCoordinateCount)
        return;

      int64_t validEnd = std::min(end, totalCoordinateCount - 1);
      if (validEnd < start)
        return;

      // Continue with adjusted range
      GapUpdate adjustedUpdate{del, {start, validEnd}};
      return updateGapMapStep(gapMap, adjustedUpdate, backtrack, gapMapUpdates,
                              recordGapMapUpdates);
    }
  }

  if (del) {
    // Record backtrack info before erasing
    auto it = gapMap.find(start);
    if (it != gapMap.end()) {
      backtrack.emplace_back(false, std::make_pair(start, it->second));
      if (recordGapMapUpdates) {
        gapMapUpdates.emplace_back(true, std::make_pair(start, it->second));
      }
      gapMap.erase(it);
    }
  } else {
    // Calculate length - to be used for storage
    int64_t newLength = end - start + 1;

    // Skip invalid entries
    if (newLength <= 0) {
      std::cerr << "Warning: Skipping gap update with non-positive length: "
                << newLength << " at position " << start << std::endl;
      return;
    }

    // SIMPLIFIED APPROACH: Each entry simply indicates the number of positions
    // to the next non-gap, with no need for overlapping entry detection or
    // merging

    // Check if there's already an entry at this position
    auto existingEntry = gapMap.find(start);
    if (existingEntry != gapMap.end()) {
      // Record the existing entry for backtracking
      backtrack.emplace_back(false,
                             std::make_pair(start, existingEntry->second));
      if (recordGapMapUpdates) {
        gapMapUpdates.emplace_back(
            true, std::make_pair(start, existingEntry->second));
      }
    } else {
      // For new entries, record a "null" entry for backtracking
      backtrack.emplace_back(true, std::make_pair(start, 0));
    }

    // Set the new length directly - no merging, no complex processing
    gapMap[start] = newLength;

    // Record the update
    if (recordGapMapUpdates) {
      gapMapUpdates.emplace_back(false, std::make_pair(start, newLength));
    }
  }
}

void updateGapMap(GapMap &gapMap, const std::vector<GapUpdate> &updates,
                  std::vector<GapUpdate> &backtrack,
                  std::vector<GapUpdate> &gapMapUpdates,
                  int64_t totalCoordinateCount) {
  // Optimize by reserving space for backtrack vector
  backtrack.reserve(backtrack.size() + updates.size() * 2);

  if (!updates.empty()) {
    logging::msg("Processing {} gap updates", updates.size());

    // Sort updates by position first to make processing more efficient
    std::vector<GapUpdate> sortedUpdates = updates;
    logging::msg("Sorting gap updates by position");
    std::sort(sortedUpdates.begin(), sortedUpdates.end(),
              [](const GapUpdate &a, const GapUpdate &b) {
                return a.second.first < b.second.first;
              });

    // Process updates in order of position
    int processedCount = 0;
    int largeGapCount = 0;
    for (const auto &update : sortedUpdates) {
      // Fast path for large gaps - skip expensive duplicate checking
      const auto &[del, range] = update;
      const auto &[start, length] = range;

      // Log progress every 1000 updates
      if (++processedCount % 1000 == 0) {
        logging::msg("Processed {}/{} gap updates", processedCount,
                     sortedUpdates.size());
      }

      if (!del && length > 10000) {
        // For very large gaps, just store directly without extensive checking
        // Each entry indicates positions to next non-gap position
        largeGapCount++;
        if (largeGapCount % 10 == 0) {
          logging::msg(
              "Processed {} large gaps (>10000), current: start={}, length={}",
              largeGapCount, start, length);
        }

        gapMap[start] = length;

        // Record for backtracking
        backtrack.emplace_back(true, std::make_pair(start, 0));

        // Record the update if needed
        gapMapUpdates.emplace_back(false, std::make_pair(start, length));
      } else {
        // Normal path for smaller gaps
        updateGapMapStep(gapMap, update, backtrack, gapMapUpdates, true,
                         totalCoordinateCount);
      }
    }

    logging::msg("Finished processing all {} gap updates ({} large gaps)",
                 sortedUpdates.size(), largeGapCount);
  }
}

void updateGapMap(GapMap &gapMap, const GapUpdate &update,
                  std::vector<GapUpdate> &backtrack,
                  int64_t totalCoordinateCount) {

  std::vector<GapUpdate> gapMapUpdates;
  updateGapMapStep(gapMap, update, backtrack, gapMapUpdates, false,
                   totalCoordinateCount);
}

bool validateAndFixGapMap(GapMap &gapMap, int64_t totalCoordinateCount,
                          bool fixErrors) {
  bool isValid = true;
  std::vector<int64_t> keysToRemove;
  std::map<int64_t, int64_t> keysToUpdate;
  int invalidCount = 0;
  int overlapCount = 0;
  int outOfBoundsCount = 0;
  int suspiciousCount = 0;

  // Check for negative or inverted ranges
  for (const auto &[start, length] : gapMap) {
    // Calculate actual end position
    int64_t end = start + length - 1;

    // Check negative coordinates
    if (start < 0 || length <= 0) {
      std::cerr << "Error: Invalid gap coordinates/length at position " << start
                << " (length: " << length << ")" << std::endl;
      isValid = false;
      invalidCount++;
      if (fixErrors)
        keysToRemove.push_back(start);
    }

    // Check bounds
    if (totalCoordinateCount > 0 && end >= totalCoordinateCount) {
      std::cerr << "Error: Gap extends beyond coordinates at position " << start
                << " with length " << length << ", end=" << end
                << " (max coordinate: " << totalCoordinateCount - 1 << ")"
                << std::endl;
      isValid = false;
      outOfBoundsCount++;
      if (fixErrors) {
        if (start < totalCoordinateCount) {
          // Adjust length to remain within bounds
          keysToUpdate[start] = totalCoordinateCount - start;
        } else {
          keysToRemove.push_back(start);
        }
      }
    }

    // Check for suspiciously large gaps
    if (totalCoordinateCount > 0 && length > 0.9 * totalCoordinateCount) {
      std::cerr << "Warning: Suspiciously large gap detected at position "
                << start << " with length " << length << " ("
                << (length * 100.0 / totalCoordinateCount)
                << "% of total sequence)" << std::endl;
      isValid = false;
      suspiciousCount++;
      if (fixErrors) {
        keysToRemove.push_back(start);
      }
    }
  }

  // Check for overlaps - in gap map, values are lengths, not end positions
  auto it = gapMap.begin();
  while (it != gapMap.end() && std::next(it) != gapMap.end()) {
    auto nextIt = std::next(it);
    int64_t current_end =
        it->first + it->second - 1; // Calculate real end position

    if (current_end >= nextIt->first) {
      std::cerr << "Error: Overlapping gap ranges at positions " << it->first
                << " (end: " << current_end << ") and " << nextIt->first
                << " (overlap: " << (current_end - nextIt->first + 1)
                << " bases)" << std::endl;
      isValid = false;
      overlapCount++;
      if (fixErrors) {
        // Calculate merged length considering the true end positions
        int64_t next_end = nextIt->first + nextIt->second - 1;
        int64_t merged_end = std::max(current_end, next_end);
        int64_t merged_length = merged_end - it->first + 1;

        // Merge the entries
        keysToUpdate[it->first] = merged_length;
        keysToRemove.push_back(nextIt->first);
      }
    }
    ++it;
  }

  // Provide a summary of all issues
  if (!isValid) {
    std::cerr << "Gap map validation summary: " << std::endl;
    std::cerr << "  - " << invalidCount
              << " gaps with invalid coordinates/lengths" << std::endl;
    std::cerr << "  - " << outOfBoundsCount
              << " gaps extending beyond coordinate bounds" << std::endl;
    std::cerr << "  - " << overlapCount << " overlapping gap pairs"
              << std::endl;
    std::cerr << "  - " << suspiciousCount
              << " suspiciously large gaps (>90% of sequence)" << std::endl;

    if (fixErrors) {
      std::cerr << "Applying fixes: " << std::endl;
      std::cerr << "  - Removing " << keysToRemove.size() << " invalid gaps"
                << std::endl;
      std::cerr << "  - Adjusting " << keysToUpdate.size() << " gaps"
                << std::endl;
    }
  }

  // Fix errors if requested
  if (fixErrors) {
    // Remove invalid entries
    for (int64_t key : keysToRemove) {
      gapMap.erase(key);
    }

    // Update truncated entries
    for (const auto &[key, newLength] : keysToUpdate) {
      gapMap[key] = newLength;
    }

    // Verify that all issues are fixed
    if (!keysToRemove.empty() || !keysToUpdate.empty()) {
      bool stillHasIssues = false;

      // Check for overlaps again after fixes
      it = gapMap.begin();
      while (it != gapMap.end() && std::next(it) != gapMap.end()) {
        auto nextIt = std::next(it);
        int64_t current_end = it->first + it->second - 1;

        if (current_end >= nextIt->first) {
          std::cerr << "Error: Overlapping gap ranges still exist after fixes!"
                    << std::endl;
          stillHasIssues = true;
          break;
        }
        ++it;
      }

      if (!stillHasIssues) {
        std::cerr << "All gap map issues successfully fixed." << std::endl;
        return true;
      }
    }
  }

  return isValid;
}

interval_set<int64_t>
gapMapToNucRunSet(const GapMap &gapMap,
                  const std::vector<std::pair<int64_t, int64_t>> &blockRanges) {

  interval_set<int64_t> curNucRunSet;
  int64_t start = -1;
  for (const auto &[curStart, length] : gapMap) {
    // Calculate actual end position
    int64_t curEnd = curStart + length - 1;

    if (start == -1) {
      if (curStart != 0) {
        curNucRunSet.add(interval<int64_t>::closed(0, curStart - 1));
      }
    } else {
      curNucRunSet.add(interval<int64_t>::closed(start, curStart - 1));
    }
    start = curEnd + 1;
  }
  if (start <= blockRanges.back().second) {
    curNucRunSet.add(
        interval<int64_t>::closed(start, blockRanges.back().second));
  }
  return curNucRunSet;
}

std::vector<GapRange> invertRanges(const std::vector<GapRange> &nucRanges,
                                   const GapRange &invertRange) {
  std::vector<GapRange> result;
  result.reserve(nucRanges.size());

  const auto &[blockStart, blockEnd] = invertRange;

  // Quick validation of invert range
  if (blockStart < 0 || blockEnd < 0 || blockStart > blockEnd) {
    return result; // Return empty for invalid invert range
  }

  for (const auto &range : nucRanges) {
    const auto &[rangeStart, rangeEnd] = range;

    // Skip invalid ranges
    if (rangeStart < 0 || rangeEnd < 0 || rangeStart > rangeEnd ||
        rangeStart < blockStart || rangeEnd > blockEnd) {
      continue;
    }

    // Create the inverted range
    int64_t invertedStart = blockEnd - (rangeEnd - blockStart);
    int64_t invertedEnd = blockEnd - (rangeStart - blockStart);

    // Skip invalid results (should never happen with valid inputs)
    if (invertedStart < 0 || invertedEnd < 0 || invertedStart > invertedEnd) {
      continue;
    }

    result.emplace_back(invertedStart, invertedEnd);
  }

  std::sort(result.begin(), result.end());
  return result;
}

void invertGapMap(GapMap &gapMap, const GapRange &invertRange,
                  std::vector<GapUpdate> &backtrack,
                  std::vector<GapUpdate> &gapMapUpdates) {

  const auto &[blockStart, blockEnd] = invertRange;
  bool recordGapMapUpdates = !gapMapUpdates.empty();

  // Validate the inversion range - return early if invalid
  if (blockStart < 0 || blockEnd < blockStart) {
    std::cerr << "Warning: Invalid inversion range [" << blockStart << ":"
              << blockEnd << "]. Skipping gap map inversion." << std::endl;
    return;
  }

  // 1. Backup all gaps in affected range for restoration
  std::map<int64_t, int64_t> originalGaps;
  for (auto it = gapMap.lower_bound(blockStart);
       it != gapMap.end() && it->first <= blockEnd; ++it) {
    originalGaps[it->first] = it->second;

    // Record backtracking info
    backtrack.push_back({false, {it->first, it->second}});
    if (recordGapMapUpdates) {
      gapMapUpdates.push_back({true, {it->first, it->second}});
    }
  }

  // 2. Handle gaps that extend from before into our range
  if (blockStart > 0) {
    auto prevIt = gapMap.upper_bound(blockStart);
    if (prevIt != gapMap.begin()) {
      prevIt = std::prev(prevIt);
      int64_t prevEnd = prevIt->first + prevIt->second - 1;

      if (prevEnd >= blockStart) {
        // Record original for backtracking
        backtrack.push_back({false, {prevIt->first, prevIt->second}});
        if (recordGapMapUpdates) {
          gapMapUpdates.push_back({true, {prevIt->first, prevIt->second}});
        }

        // Truncate the portion that overlaps with our range
        int64_t newLength = blockStart - prevIt->first;
        if (newLength > 0) {
          gapMap[prevIt->first] = newLength;
          if (recordGapMapUpdates) {
            gapMapUpdates.push_back({false, {prevIt->first, newLength}});
          }
        } else {
          gapMap.erase(prevIt->first);
        }
      }
    }
  }

  // 3. Remove all gaps in the range - we'll recompute them
  for (auto it = gapMap.lower_bound(blockStart);
       it != gapMap.end() && it->first <= blockEnd;) {
    it = gapMap.erase(it);
  }

  // 4. Extract gap positions from original gaps
  std::vector<GapRange> nucRanges;
  for (const auto &[gapStart, gapLength] : originalGaps) {
    int64_t gapEnd = gapStart + gapLength - 1;
    if (gapStart >= blockStart && gapEnd <= blockEnd) {
      nucRanges.emplace_back(gapStart, gapEnd);
    }
  }

  // 5. Invert the ranges using our invertRanges function
  auto invertedRanges = invertRanges(nucRanges, invertRange);

  // 6. Convert to CoordRange, merge, and convert back using RangeOperations
  auto coordRanges =
      coordinates::RangeOperations::gapRangesToCoordRanges(invertedRanges);
  auto mergedCoordRanges =
      coordinates::RangeOperations::mergeRanges(coordRanges);
  auto mergedRanges =
      coordinates::RangeOperations::coordRangesToGapRanges<GapRange>(
          mergedCoordRanges);

  // 7. Add the merged inverted ranges back to the gap map
  for (const auto &range : mergedRanges) {
    int64_t length = range.second - range.first + 1;
    if (length > 0) {
      gapMap[range.first] = length;

      // Record backtracking info
      backtrack.push_back({true, {range.first, 0}});
      if (recordGapMapUpdates) {
        gapMapUpdates.push_back({false, {range.first, length}});
      }
    }
  }

  // 8. Apply preventOverlaps to fix any overlaps
  preventOverlaps(gapMap);

  // 9. Check for and remove suspiciously large gaps
  int64_t totalSize = blockEnd - blockStart + 1;
  for (auto it = gapMap.lower_bound(blockStart);
       it != gapMap.end() && it->first <= blockEnd;) {
    double coverage = static_cast<double>(it->second) / totalSize;
    if (coverage > 0.9) {
      // Record removal for backtracking
      backtrack.push_back({true, {it->first, 0}});
      if (recordGapMapUpdates) {
        gapMapUpdates.push_back({true, {it->first, it->second}});
      }

      // Remove the suspicious gap
      it = gapMap.erase(it);
    } else {
      ++it;
    }
  }

  // 10. Ensure backtracking info is complete
  std::set<int64_t> backtrackedPositions;
  for (const auto &update : backtrack) {
    backtrackedPositions.insert(update.second.first);
  }

  // Check that we have backtracking entries for all current gaps in the range
  for (auto it = gapMap.lower_bound(blockStart);
       it != gapMap.end() && it->first <= blockEnd; ++it) {
    if (backtrackedPositions.find(it->first) == backtrackedPositions.end()) {
      backtrack.push_back({true, {it->first, 0}});
      if (recordGapMapUpdates) {
        gapMapUpdates.push_back({false, {it->first, it->second}});
      }
    }
  }
}

void invertGapMap(GapMap &gapMap, const GapRange &invertRange,
                  std::vector<GapUpdate> &backtrack) {

  std::vector<GapUpdate> gapMapUpdates;
  invertGapMap(gapMap, invertRange, backtrack, gapMapUpdates);
}

void makeCoordIndex(GapMap &degapCoordIndex, GapMap &regapCoordIndex,
                    const GapMap &gapMap,
                    const std::vector<GapRange> &blockRanges) {

  degapCoordIndex.clear();
  regapCoordIndex.clear();

  int64_t degapCtr = 0;
  int64_t lastEnd = -1;

  for (const auto &[start, length] : gapMap) {
    // Calculate actual end position
    int64_t end = start + length - 1;

    if (lastEnd != -1) {
      for (int64_t i = lastEnd + 1; i < start; ++i) {
        degapCoordIndex[i] = degapCtr;
        regapCoordIndex[degapCtr] = i;
        ++degapCtr;
      }
    }
    lastEnd = end;
  }

  // Handle any remaining coordinates after the last gap
  if (!blockRanges.empty() && lastEnd != -1) {
    const auto &[_, blockEnd] = blockRanges.back();
    for (int64_t i = lastEnd + 1; i <= blockEnd; ++i) {
      degapCoordIndex[i] = degapCtr;
      regapCoordIndex[degapCtr] = i;
      ++degapCtr;
    }
  }
}

/**
 * Implementation of suspicious gap removal function
 */
int removeSuspiciousGaps(GapMap &gapMap, int64_t totalSequenceSize,
                         double threshold) {
  // Validate inputs
  if (totalSequenceSize <= 0) {
    std::cerr << "Error: Invalid sequence size for suspicious gap detection: "
              << totalSequenceSize << std::endl;
    return 0;
  }

  std::vector<int64_t> keysToRemove;

  // First pass: identify suspicious gaps
  for (const auto &[start, length] : gapMap) {
    // Check if this gap covers a suspiciously large portion of the sequence
    double coverage = static_cast<double>(length) / totalSequenceSize;

    if (coverage > threshold) {
      keysToRemove.push_back(start);
      std::cerr << "Warning: Found suspicious gap at position " << start
                << " with length " << length << " (" << (coverage * 100.0)
                << "% of sequence), will be removed." << std::endl;
    }
  }

  // Second pass: remove identified gaps
  for (const auto &key : keysToRemove) {
    gapMap.erase(key);
  }

  // Return number of removed gaps
  if (!keysToRemove.empty()) {
    std::cerr << "Removed " << keysToRemove.size()
              << " suspicious gaps that covered more than "
              << (threshold * 100.0) << "% of the sequence." << std::endl;
  }

  return keysToRemove.size();
}

/**
 * Helper function to merge adjacent or overlapping range-based entries
 *
 * @param entries Vector of pairs with (start, length) for each range
 * @return Vector of merged entries with no overlaps
 */
std::vector<std::pair<int64_t, int64_t>>
mergeAdjacentRanges(const std::vector<std::pair<int64_t, int64_t>> &entries) {

  if (entries.empty()) {
    return {};
  }

  // Sort by start position
  std::vector<std::pair<int64_t, int64_t>> sorted = entries;
  std::sort(sorted.begin(), sorted.end());

  std::vector<std::pair<int64_t, int64_t>> result;
  result.reserve(sorted.size());

  // Add first entry as current
  auto current = sorted.front();

  // Process and merge subsequent entries
  for (size_t i = 1; i < sorted.size(); ++i) {
    int64_t curEnd = current.first + current.second - 1;
    int64_t nextStart = sorted[i].first;

    // If current overlaps or is adjacent to next, merge them
    if (curEnd + 1 >= nextStart) {
      int64_t nextEnd = sorted[i].first + sorted[i].second - 1;
      current.second = std::max(curEnd, nextEnd) - current.first + 1;
    } else {
      // No overlap, add current to result and move to next
      result.push_back(current);
      current = sorted[i];
    }
  }

  // Add final merged entry
  result.push_back(current);

  return result;
}

/**
 * Implementation of overlap prevention
 */
bool preventOverlaps(GapMap &gapMap) {
  bool changesWereMade = false;
  bool debugMode = false; // Set to true to enable verbose logging

  // First pass: remove invalid entries (negative positions or non-positive
  // lengths)
  for (auto it = gapMap.begin(); it != gapMap.end();) {
    if (it->first < 0 || it->second <= 0) {
      if (debugMode) {
        std::cerr << "Warning: Removing invalid gap entry at position "
                  << it->first << " with length " << it->second << std::endl;
      }
      it = gapMap.erase(it);
      changesWereMade = true;
    } else {
      ++it;
    }
  }

  // Check for suspiciously large gaps (covering almost the entire sequence)
  // Get the largest position + length to estimate sequence size
  int64_t maxPosition = 0;
  for (const auto &[pos, len] : gapMap) {
    maxPosition = std::max(maxPosition, pos + len);
  }

  if (maxPosition > 0) {
    for (auto it = gapMap.begin(); it != gapMap.end();) {
      // If a gap starts at or near the beginning and covers more than 90% of
      // the sequence
      if (it->first <= 10 && it->second > maxPosition * 0.9) {
        if (debugMode) {
          std::cerr
              << "Warning: Removing suspicious full-sequence gap from position "
              << it->first << " with length " << it->second << " (covers "
              << (100.0 * it->second / maxPosition) << "% of the sequence)"
              << std::endl;
        }
        it = gapMap.erase(it);
        changesWereMade = true;
      } else {
        ++it;
      }
    }
  }

  // Second pass: merge overlapping or adjacent gaps using the shared helper
  std::vector<std::pair<int64_t, int64_t>> gapEntries;
  gapEntries.reserve(gapMap.size());

  // Convert map to vector of entries
  for (const auto &[pos, len] : gapMap) {
    gapEntries.emplace_back(pos, len);
  }

  // Skip if no entries to merge
  if (!gapEntries.empty()) {
    // Merge overlapping/adjacent entries
    auto mergedEntries = mergeAdjacentRanges(gapEntries);

    // If the size changed, we need to update the map
    if (mergedEntries.size() != gapEntries.size()) {
      changesWereMade = true;

      // Clear and rebuild the map with merged entries
      gapMap.clear();
      for (const auto &[pos, len] : mergedEntries) {
        gapMap[pos] = len;
      }
    }
  }

  // Final verification - ensure the gap map is strictly ordered with no
  // overlaps
  bool is_valid = true;
  int64_t prev_end =
      -2; // Start at -2 so adjacent gaps at position 0 can be detected

  for (const auto &[pos, len] : gapMap) {
    int64_t end = pos + len - 1;

    // Check that this gap doesn't overlap with previous gap
    if (pos <= prev_end + 1) {
      if (debugMode) {
        std::cerr << "Warning: Gap map still has overlaps after prevention. "
                  << "Previous gap ended at " << prev_end
                  << ", but new gap starts at " << pos << std::endl;
      }
      is_valid = false;
    }

    // Update previous end for next iteration
    prev_end = end;
  }

  if (!is_valid && debugMode) {
    std::cerr << "Warning: Unable to fully resolve all overlaps. The gap map "
                 "may still be inconsistent."
              << std::endl;
  }

  return changesWereMade;
}

size_t batchProcessGapUpdates(GapMap &gapMap,
                              const std::vector<GapUpdate> &updates,
                              std::vector<GapUpdate> &backtrack,
                              int64_t totalCoordinateCount) {
  if (updates.empty()) {
    return 0;
  }

  logging::debug("Batch processing {} gap updates", updates.size());

  // Separate updates by type
  std::vector<std::pair<int64_t, int64_t>> additions, removals;
  additions.reserve(updates.size());
  removals.reserve(updates.size());

  for (const auto &update : updates) {
    if (update.first) { // Deletion
      removals.push_back(update.second);
    } else { // Addition
      additions.push_back(update.second);
    }
  }

  // Use the shared helper to merge ranges of same type
  auto mergedAdditions = mergeAdjacentRanges(additions);
  auto mergedRemovals = mergeAdjacentRanges(removals);

  // Create merged updates list with removals first
  std::vector<GapUpdate> mergedUpdates;
  mergedUpdates.reserve(mergedAdditions.size() + mergedRemovals.size());

  // Add removals first (they should be processed first)
  for (const auto &range : mergedRemovals) {
    mergedUpdates.emplace_back(true, range);
  }

  // Then add additions
  for (const auto &range : mergedAdditions) {
    mergedUpdates.emplace_back(false, range);
  }

  // Log if we achieved reduction
  if (mergedUpdates.size() < updates.size()) {
    logging::debug("Merged into {} updates ({}% reduction)",
                   mergedUpdates.size(),
                   100.0 * (1.0 - static_cast<double>(mergedUpdates.size()) /
                                      updates.size()));
  }

  // Apply the merged updates
  std::vector<GapUpdate>
      dummyUpdates; // Empty vector for updates we don't need to track
  for (const auto &update : mergedUpdates) {
    updateGapMapStep(gapMap, update, backtrack, dummyUpdates, false,
                     totalCoordinateCount);
  }

  return mergedUpdates.size();
}

GapStatistics calculateGapStatistics(const GapMap &gapMap) {
  GapStatistics stats;
  stats.update(gapMap);

  // Log statistics
  logging::info("Gap Statistics: {} runs, {} positions, max length: {}, avg "
                "length: {:.2f}",
                stats.totalGapRuns, stats.totalGapPositions,
                stats.maxGapRunLength, stats.avgGapRunLength);
  logging::debug(
      "Gap run sizes: small: {} ({}%), medium: {} ({}%), large: {} ({}%)",
      stats.smallGapRuns,
      stats.totalGapRuns > 0 ? 100.0 * stats.smallGapRuns / stats.totalGapRuns
                             : 0.0,
      stats.mediumGapRuns,
      stats.totalGapRuns > 0 ? 100.0 * stats.mediumGapRuns / stats.totalGapRuns
                             : 0.0,
      stats.largeGapRuns,
      stats.totalGapRuns > 0 ? 100.0 * stats.largeGapRuns / stats.totalGapRuns
                             : 0.0);

  return stats;
}

// Implementation of the missing functions that are causing linker errors
void addGap(GapMap &gapMap, int64_t pos, int64_t length) {
  if (length <= 0)
    return;

  // Find the insertion point
  auto it = gapMap.upper_bound(pos);

  // Check if we can merge with previous gap
  if (it != gapMap.begin()) {
    auto prev = std::prev(it);
    int64_t prevStart = prev->first;
    int64_t prevLength = prev->second;
    int64_t prevEnd = prevStart + prevLength - 1;

    // If previous gap overlaps or is adjacent to new gap, merge
    if (prevEnd >= pos - 1) {
      // Extend previous gap
      prev->second = std::max(prevEnd, pos + length - 1) - prevStart + 1;

      // Check if we need to merge with subsequent gaps too
      auto next = it;
      while (next != gapMap.end()) {
        int64_t nextStart = next->first;
        int64_t nextLength = next->second;
        int64_t nextEnd = nextStart + nextLength - 1;

        // If this gap overlaps with the extended previous gap, merge
        if (nextStart <= prevStart + prev->second) {
          // Update length of previous gap to include this one
          prev->second =
              std::max(prevStart + prev->second - 1, nextEnd) - prevStart + 1;
          next = gapMap.erase(next);
        } else {
          break;
        }
      }

      return;
    }
  }

  // Check if we can merge with next gap
  if (it != gapMap.end()) {
    int64_t nextStart = it->first;
    int64_t nextEnd = nextStart + it->second - 1;

    // If new gap overlaps or is adjacent to next gap, merge
    if (pos + length >= nextStart - 1) {
      // Create merged gap
      int64_t mergedStart = pos;
      int64_t mergedEnd = std::max(pos + length - 1, nextEnd);
      int64_t mergedLength = mergedEnd - mergedStart + 1;

      // Replace next gap with merged one - instead of modifying first, erase
      // and insert
      int64_t oldLength = it->second;
      gapMap.erase(it);
      it = gapMap.emplace(mergedStart, mergedLength).first;

      // Check if we need to merge with subsequent gaps
      auto next = std::next(it);
      while (next != gapMap.end()) {
        int64_t nextStart = next->first;
        int64_t nextLength = next->second;
        int64_t nextEnd = nextStart + nextLength - 1;

        // If this gap overlaps with the merged gap, merge
        if (nextStart <= mergedStart + it->second) {
          // Update merged gap to include this one
          int64_t newEnd = std::max(mergedStart + it->second - 1, nextEnd);
          int64_t newLength = newEnd - mergedStart + 1;

          // Update the merged gap - erase and insert again
          gapMap.erase(it);
          it = gapMap.emplace(mergedStart, newLength).first;

          next = gapMap.erase(next);
        } else {
          break;
        }
      }

      return;
    }
  }

  // If we can't merge, insert new gap
  gapMap.emplace(pos, length);
}

void removeGap(GapMap &gapMap, int64_t pos, int64_t length) {
  if (length <= 0)
    return;

  int64_t end = pos + length - 1;

  // Find all gaps that overlap with the removal range
  auto it = gapMap.lower_bound(pos);
  if (it != gapMap.begin()) {
    // Check if previous gap overlaps with removal range
    auto prev = std::prev(it);
    int64_t prevStart = prev->first;
    int64_t prevLength = prev->second;
    int64_t prevEnd = prevStart + prevLength - 1;

    if (prevEnd >= pos) {
      // Previous gap overlaps with removal range

      // Case 1: Removal completely within gap
      if (pos > prevStart && end < prevEnd) {
        // Split the gap
        int64_t firstLength = pos - prevStart;
        int64_t secondStart = end + 1;
        int64_t secondLength = prevEnd - end;

        // Update previous gap
        prev->second = firstLength;

        // Insert new gap for second part
        gapMap.emplace(secondStart, secondLength);
        return;
      }

      // Case 2: Removal at start of gap
      if (pos <= prevStart && end < prevEnd) {
        // Shorten gap at start
        int64_t newStart = end + 1;
        int64_t newLength = prevEnd - end;

        // Replace old gap
        gapMap.erase(prev);
        gapMap.emplace(newStart, newLength);
        it = gapMap.lower_bound(newStart);
      }

      // Case 3: Removal at end of gap
      else if (pos > prevStart && end >= prevEnd) {
        // Shorten gap at end
        prev->second = pos - prevStart;
      }

      // Case 4: Removal completely covers gap
      else if (pos <= prevStart && end >= prevEnd) {
        // Remove the gap entirely
        it = gapMap.erase(prev);
      }
    }
  }

  // Process any remaining gaps that overlap with the removal range
  while (it != gapMap.end() && it->first <= end) {
    int64_t gapStart = it->first;
    int64_t gapLength = it->second;
    int64_t gapEnd = gapStart + gapLength - 1;

    // Case 1: Removal completely within gap
    if (pos > gapStart && end < gapEnd) {
      // Split the gap
      int64_t firstLength = pos - gapStart;
      int64_t secondStart = end + 1;
      int64_t secondLength = gapEnd - end;

      // Update current gap
      it->second = firstLength;

      // Insert new gap for second part
      gapMap.emplace(secondStart, secondLength);
      break;
    }

    // Case 2: Removal at start of gap
    if (pos <= gapStart && end < gapEnd) {
      // Shorten gap at start
      int64_t newStart = end + 1;
      int64_t newLength = gapEnd - end;

      // Replace old gap
      auto next = std::next(it);
      gapMap.erase(it);
      gapMap.emplace(newStart, newLength);
      it = next;
      continue;
    }

    // Case 3: Removal at end of gap
    if (pos > gapStart && end >= gapEnd) {
      // Shorten gap at end
      it->second = pos - gapStart;
      it++;
      continue;
    }

    // Case 4: Removal completely covers gap
    if (pos <= gapStart && end >= gapEnd) {
      // Remove the gap entirely
      it = gapMap.erase(it);
      continue;
    }

    it++;
  }
}

} // namespace gap_map

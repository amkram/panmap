#include "gap_map.hpp"
#include "timing.hpp"
#include <algorithm>
#include <boost/icl/interval_set.hpp>
#include <iostream>

using namespace boost::icl;

namespace gap_map {

void updateGapMapStep(GapMap &gapMap, const GapUpdate &update,
                      std::vector<GapUpdate> &backwtrack,
                      std::vector<GapUpdate> &gapMapUpdates,
                      bool recordGapMapUpdates) {

  TIME_FUNCTION;
  const auto &[del, range] = update;
  const auto &[start, end] = range;

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
    // Record backtrack info for any existing mapping
    auto it = gapMap.find(start);
    if (it != gapMap.end()) {
      backtrack.emplace_back(false, std::make_pair(start, it->second));
    } else {
      backtrack.emplace_back(true, std::make_pair(start, end));
    }

    if (recordGapMapUpdates) {
      gapMapUpdates.emplace_back(false, std::make_pair(start, end));
    }
    gapMap[start] = end;
  }
}

void updateGapMap(GapMap &gapMap, const std::vector<GapUpdate> &updates,
                  std::vector<GapUpdate> &backtrack,
                  std::vector<GapUpdate> &gapMapUpdates) {

  TIME_FUNCTION;
  for (const auto &update : updates) {
    updateGapMapStep(gapMap, update, backtrack, gapMapUpdates, true);
  }
}

void updateGapMap(GapMap &gapMap, const GapUpdate &update,
                  std::vector<GapUpdate> &backtrack) {

  std::vector<GapUpdate> gapMapUpdates;
  updateGapMapStep(gapMap, update, backtrack, gapMapUpdates, false);
}

interval_set<int64_t>
gapMapToNucRunSet(const GapMap &gapMap,
                  const std::vector<std::pair<int64_t, int64_t>> &blockRanges) {

  interval_set<int64_t> curNucRunSet;
  int64_t start = -1;
  for (const auto &[curStart, curEnd] : gapMap) {
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

  TIME_FUNCTION;
  std::vector<GapRange> result;
  result.reserve(nucRanges.size());

  const auto &[blockStart, blockEnd] = invertRange;

  // Validate the invert range itself
  if (blockStart < 0 || blockEnd < 0 || blockStart > blockEnd) {
    std::cerr << "Warning: Invalid invert range: [" << blockStart << ", "
              << blockEnd << "]" << std::endl;
    return result; // Return empty result
  }

  for (const auto &range : nucRanges) {
    // Validate input range
    const auto &[rangeStart, rangeEnd] = range;
    if (rangeStart < 0 || rangeEnd < 0 || rangeStart > rangeEnd) {
      std::cerr << "Warning: Skipping invalid range: [" << rangeStart << ", "
                << rangeEnd << "]" << std::endl;
      continue;
    }

    // Ensure the range is within the block range
    if (rangeStart < blockStart || rangeEnd > blockEnd) {
      std::cerr << "Warning: Range [" << rangeStart << ", " << rangeEnd
                << "] outside invert range [" << blockStart << ", " << blockEnd
                << "]" << std::endl;
      continue;
    }

    // Create the inverted range
    int64_t invertedStart = blockEnd - (rangeEnd - blockStart);
    int64_t invertedEnd = blockEnd - (rangeStart - blockStart);

    // Final validation
    if (invertedStart < 0 || invertedEnd < 0 || invertedStart > invertedEnd) {
      std::cerr << "Warning: Invalid inverted range computed: ["
                << invertedStart << ", " << invertedEnd << "]" << std::endl;
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

  TIME_FUNCTION;
  const auto &[blockStart, blockEnd] = invertRange;

  std::vector<GapRange> nucRanges;
  auto it = gapMap.lower_bound(blockStart);

  // Collect ranges that need to be inverted
  while (it != gapMap.end() && it->first <= blockEnd) {
    if (it->second <= blockEnd) {
      nucRanges.emplace_back(it->first, it->second);
    }
    ++it;
  }

  // Remove old ranges
  for (const auto &range : nucRanges) {
    updateGapMapStep(gapMap, {true, range}, backtrack, gapMapUpdates, true);
  }

  // Add inverted ranges
  auto invertedRanges = invertRanges(nucRanges, invertRange);
  for (const auto &range : invertedRanges) {
    updateGapMapStep(gapMap, {false, range}, backtrack, gapMapUpdates, true);
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

  TIME_FUNCTION;
  degapCoordIndex.clear();
  regapCoordIndex.clear();

  int64_t degapCtr = 0;
  int64_t lastEnd = -1;

  for (const auto &[start, end] : gapMap) {
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

} // namespace gap_map

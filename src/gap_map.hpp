#ifndef GAP_MAP_HPP
#define GAP_MAP_HPP

#include "coordinates.hpp"
#include <boost/icl/interval_set.hpp>
#include <vector>

namespace gap_map {

using namespace boost::icl;
using coordinates::GapMap;
using coordinates::GapRange;
using coordinates::GapUpdate;

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
 */
void updateGapMapStep(GapMap &gapMap, const GapUpdate &update,
                      std::vector<GapUpdate> &backtrack,
                      std::vector<GapUpdate> &gapMapUpdates,
                      bool recordGapMapUpdates = true);

/**
 * @brief Update a gap map with multiple updates
 * @param gapMap Gap map to be updated
 * @param updates Vector of updates to be applied
 * @param backtrack Vector to store backtrack information
 * @param gapMapUpdates Vector to store all updates
 */
void updateGapMap(GapMap &gapMap, const std::vector<GapUpdate> &updates,
                  std::vector<GapUpdate> &backtrack,
                  std::vector<GapUpdate> &gapMapUpdates);

/**
 * @brief Simplified gap map update for single update
 * @param gapMap Gap map to be updated
 * @param update Single update to apply
 * @param backtrack Vector to store backtrack information
 */
void updateGapMap(GapMap &gapMap, const GapUpdate &update,
                  std::vector<GapUpdate> &backtrack);

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
 * @brief Create coordinate indices for a gap map
 * @param degapCoordIndex Map for degapped coordinate lookup
 * @param regapCoordIndex Map for regapped coordinate lookup
 * @param gapMap Source gap map
 * @param blockRanges Block range information
 */
void makeCoordIndex(GapMap &degapCoordIndex, GapMap &regapCoordIndex,
                    const GapMap &gapMap,
                    const std::vector<GapRange> &blockRanges);

} // namespace gap_map

#endif // GAP_MAP_HPP
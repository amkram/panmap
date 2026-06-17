#pragma once

/**
 * @file traversal.hpp
 * @brief Reconstruct a node's genome seed set by accumulating index deltas along the
 *        root..node path. Absorbs the copy-pasted delta-traversal + nodeChangeOffsets
 *        snippets. Single concern: produce the seed multiset; metric checks are done
 *        elsewhere against the live PlacementGlobalState.
 */

#include "metrics_oracle.hpp"  // indexUtils::SeedCountMap
#include "panmap_utils.hpp"
#include "test_index.hpp"

#include <vector>

namespace ts {

// Apply seed-change deltas (childCount is the absolute count; <=0 means deleted)
// for each node along the path, returning the final hash -> genome count map.
indexUtils::SeedCountMap reconstructGenomeSeeds(const std::vector<panmapUtils::LiteNode*>& path, const IndexData& idx);

}  // namespace ts

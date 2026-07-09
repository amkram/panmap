#pragma once

/**
 * @file tree_helpers.hpp
 * @brief LiteTree lookup helpers + the shared RSV panman fixture.
 */

#include "panman.hpp"
#include "panmanUtils.hpp"
#include "panmap_utils.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace ts {

// Find a node's DFS index by identifier; -1 if not present.
int32_t findNodeIndex(const panmapUtils::LiteTree& tree, const std::string& id);

// Path from root down to the target node (inclusive), in root..target order.
std::vector<panmapUtils::LiteNode*> pathToRoot(const panmapUtils::LiteTree& tree, int32_t targetDfsIndex);

// Loads src/test/data/rsv_4K.panman once; exposes the tree and on-demand genome
// reconstruction so tests don't need committed per-node FASTA files.
struct RSVPanmanFixture {
    RSVPanmanFixture();

    panmanUtils::Tree* tree() const { return T; }

    std::string genomeOf(const std::string& nodeId) const;
    std::vector<std::string> nodeIds() const;

    std::unique_ptr<panmanUtils::TreeGroup> TG;
    panmanUtils::Tree* T = nullptr;
};

}  // namespace ts

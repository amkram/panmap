#include "traversal.hpp"

namespace ts {

indexUtils::SeedCountMap reconstructGenomeSeeds(const std::vector<panmapUtils::LiteNode*>& path, const IndexData& idx) {
    indexUtils::SeedCountMap seeds;
    for (const auto* node : path) {
        const uint32_t nodeIdx = node->nodeIndex;
        const uint64_t start = idx.nodeChangeOffset(nodeIdx);
        const uint64_t end = idx.nodeChangeOffset(nodeIdx + 1);
        for (uint64_t j = start; j < end; ++j) {
            const uint64_t hash = idx.hashAt(j);
            const int64_t childCount = idx.childCountAt(j);
            if (childCount > 0) {
                seeds[hash] = childCount;
            } else {
                seeds.erase(hash);
            }
        }
    }
    return seeds;
}

}  // namespace ts

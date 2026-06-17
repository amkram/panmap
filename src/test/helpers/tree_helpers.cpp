#include "tree_helpers.hpp"

#include "paths.hpp"

#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <algorithm>
#include <fstream>
#include <functional>
#include <stdexcept>

namespace ts {

int32_t findNodeIndex(const panmapUtils::LiteTree& tree, const std::string& id) {
    for (size_t i = 0; i < tree.dfsIndexToNode.size(); ++i) {
        if (tree.dfsIndexToNode[i] && tree.dfsIndexToNode[i]->identifier == id) {
            return static_cast<int32_t>(i);
        }
    }
    return -1;
}

std::vector<panmapUtils::LiteNode*> pathToRoot(const panmapUtils::LiteTree& tree, int32_t targetDfsIndex) {
    std::vector<panmapUtils::LiteNode*> path;
    if (targetDfsIndex < 0 || static_cast<size_t>(targetDfsIndex) >= tree.dfsIndexToNode.size()) {
        return path;
    }
    for (panmapUtils::LiteNode* cur = tree.dfsIndexToNode[targetDfsIndex]; cur != nullptr; cur = cur->parent) {
        path.push_back(cur);
    }
    std::reverse(path.begin(), path.end());
    return path;
}

RSVPanmanFixture::RSVPanmanFixture() {
    const std::string panmanPath = dataPath("rsv_4K.panman");
    std::ifstream inputFile(panmanPath, std::ios::binary);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Failed to open test panman: " + panmanPath);
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inBuffer;
    inBuffer.push(boost::iostreams::lzma_decompressor());
    inBuffer.push(inputFile);

    std::istream inputStream(&inBuffer);
    TG = std::make_unique<panmanUtils::TreeGroup>(inputStream);
    if (TG->trees.empty()) {
        throw std::runtime_error("Test panman has no trees: " + panmanPath);
    }
    T = &TG->trees[0];
}

std::string RSVPanmanFixture::genomeOf(const std::string& nodeId) const {
    return T->getStringFromReference(nodeId, false);
}

std::vector<std::string> RSVPanmanFixture::nodeIds() const {
    std::vector<std::string> ids;
    std::function<void(panmanUtils::Node*)> collect = [&](panmanUtils::Node* node) {
        if (!node) return;
        ids.push_back(node->identifier);
        for (auto* child : node->children) collect(child);
    };
    collect(T->root);
    return ids;
}

}  // namespace ts

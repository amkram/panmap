
#include "panmap_utils.hpp"
#include <queue>

namespace panmapUtils {

uint32_t LiteTree::getBlockStartScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].first;
}

uint32_t LiteTree::getBlockEndScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].second;
}

void LiteTree::initialize(::LiteTree::Reader liteTreeReader) {
  // initialize blockScalarRanges
  auto blockScalarRangesReader = liteTreeReader.getBlockRanges();
  blockScalarRanges.resize(blockScalarRangesReader.size());
  for (size_t i = 0; i < blockScalarRangesReader.size(); i++) {
    blockScalarRanges[i] = {blockScalarRangesReader[i].getRangeBeg(), blockScalarRangesReader[i].getRangeEnd()};
  }

  // initialize allLiteNodes
  auto liteNodesReader = liteTreeReader.getLiteNodes();
  size_t numNodes = liteNodesReader.size();
  
  // Pre-allocate dfsIndexToNode vector
  dfsIndexToNode.resize(numNodes, nullptr);
  
  for (size_t i = 0; i < numNodes; i++) {
    const auto liteNodeReader = liteNodesReader[i];
    const auto& nodeIdentifier = liteNodeReader.getId();
    const auto parentIndex = liteNodeReader.getParentIndex();
    nodeToDfsIndex.emplace(nodeIdentifier, i);
    auto [it, inserted] = allLiteNodes.emplace(nodeIdentifier, new LiteNode(nodeIdentifier, nullptr, {}));
    
    // Store pointer in dfsIndexToNode for index-based access
    it->second->nodeIndex = i;
    dfsIndexToNode[i] = it->second;
    
    if (i == 0) continue;
    const auto parentNodeReader = liteNodesReader[parentIndex];
    const auto& parentNodeId = parentNodeReader.getId();
    it->second->parent = allLiteNodes[parentNodeId];
    allLiteNodes[parentNodeId]->children.push_back(it->second);
  }

  root = allLiteNodes[liteNodesReader[0].getId()];
}

//end of namespace panmapUtils
}
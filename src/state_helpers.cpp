#include "state.hpp"
#include <queue>
#include <unordered_set>
#include <functional>

namespace state {

// Implementation of helper functions for tree processing utilities
template <typename T, typename ETS>
std::vector<T> mergeThreadLocalVectors(ETS& threadLocalVectors) {
  std::vector<T> mergedResult;
  
  // Calculate total size needed
  size_t totalSize = 0;
  for (const auto& localVec : threadLocalVectors) {
    totalSize += localVec.size();
  }
  
  // Pre-allocate result vector
  mergedResult.reserve(totalSize);
  
  // Merge all thread-local vectors
  for (const auto& localVec : threadLocalVectors) {
    mergedResult.insert(mergedResult.end(), localVec.begin(), localVec.end());
  }
  
  return mergedResult;
}

// Explicitly instantiate the template for the types used in the project
template std::vector<std::tuple<long, unsigned long, bool, long, long, bool, long>> 
mergeThreadLocalVectors<std::tuple<long, unsigned long, bool, long, long, bool, long>, 
                       tbb::enumerable_thread_specific<std::vector<std::tuple<long, unsigned long, bool, long, long, bool, long>>>>(
    tbb::enumerable_thread_specific<std::vector<std::tuple<long, unsigned long, bool, long, long, bool, long>>>& threadLocalVectors);

/**
 * Compute paths from root to each node in the tree
 * 
 * @param tree Pointer to the tree structure
 * @param rootNode Pointer to the root node
 * @return A map of node identifiers to their paths from root
 */
std::unordered_map<std::string, std::vector<panmanUtils::Node*>> 
computeNodePaths(panmanUtils::Tree* tree, panmanUtils::Node* rootNode) {
  std::unordered_map<std::string, std::vector<panmanUtils::Node*>> nodePaths;
  
  if (!tree || !rootNode) {
    return nodePaths;
  }

  // Start with root node's path (just itself)
  nodePaths[rootNode->identifier] = {rootNode};
  
  // Function to recursively compute paths
  std::function<void(panmanUtils::Node*, const std::vector<panmanUtils::Node*>&)> 
  computePaths = [&](panmanUtils::Node* node, const std::vector<panmanUtils::Node*>& parentPath) {
    if (!node) return;
    
    // Current node's path is parent's path plus itself
    std::vector<panmanUtils::Node*> currentPath = parentPath;
    currentPath.push_back(node);
    
    // Store the path
    nodePaths[node->identifier] = currentPath;
    
    // Compute paths for all children
    for (auto* child : node->children) {
      computePaths(child, currentPath);
    }
  };
  
  // Compute paths for all root's children
  for (auto* child : rootNode->children) {
    std::vector<panmanUtils::Node*> rootPath = {rootNode};
    computePaths(child, rootPath);
  }
  
  return nodePaths;
}

/**
 * Group nodes by their level in the tree
 * 
 * @param tree Pointer to the tree structure
 * @param rootNode Pointer to the root node
 * @return A vector of node vectors, where each inner vector contains nodes at the same level
 */
std::vector<std::vector<panmanUtils::Node*>> 
groupNodesByLevel(panmanUtils::Tree* tree, panmanUtils::Node* rootNode) {
  std::vector<std::vector<panmanUtils::Node*>> nodesByLevel;
  
  if (!tree || !rootNode) {
    return nodesByLevel;
  }
  
  // BFS queue for level traversal
  std::queue<std::pair<panmanUtils::Node*, int>> nodeQueue;
  std::unordered_set<std::string> visitedNodes;
  
  // Start with root at level 0
  nodeQueue.push({rootNode, 0});
  visitedNodes.insert(rootNode->identifier);
  
  while (!nodeQueue.empty()) {
    auto [node, level] = nodeQueue.front();
    nodeQueue.pop();
    
    // Ensure we have enough levels
    if (nodesByLevel.size() <= static_cast<size_t>(level)) {
      nodesByLevel.resize(level + 1);
    }
    
    // Add node to its level
    nodesByLevel[level].push_back(node);
    
    // Queue all unvisited children at next level
    for (auto* child : node->children) {
      if (child && visitedNodes.find(child->identifier) == visitedNodes.end()) {
        nodeQueue.push({child, level + 1});
        visitedNodes.insert(child->identifier);
      }
    }
  }
  
  return nodesByLevel;
}

} // namespace state 
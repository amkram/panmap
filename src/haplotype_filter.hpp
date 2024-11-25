#ifndef __HAPLOTYPE_FILTER_HPP
#define __HAPLOTYPE_FILTER_HPP

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <tbb/concurrent_vector.h>
#include <boost/icl/interval_map.hpp>
#include <boost/icl/split_interval_map.hpp>
#include "index.capnp.h"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Dense>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/combinable.h>
#include <tbb/task_group.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>


using namespace panmanUtils;

namespace haplotype_filter {

void noFilter(
  std::vector<std::string>& nodes, Eigen::MatrixXd& probs, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
  const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::vector<bool>& lowScoreReads, const size_t& numLowScoreReads,
  const std::string& excludeNode, const std::vector<bool>& excludeReads
) {
  std::cerr << "No preEM filtering" << std::endl;

  size_t numExcludedReads = std::count(excludeReads.begin(), excludeReads.end(), true);

  std::cerr << "Excluding " << numExcludedReads << " reads with high number of duplicates" << std::endl;
  
  if (!excludeNode.empty()) {
    probs.resize(allScores.begin()->second.size() - numExcludedReads, allScores.size() - leastRecentIdenticalAncestors.size() - 1);
  } else {
    probs.resize(allScores.begin()->second.size() - numExcludedReads, allScores.size() - leastRecentIdenticalAncestors.size());
  }
  size_t colIndex = 0;


  for (const auto& node : allScores) {
    if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
    if (!excludeNode.empty() && node.first == excludeNode) continue;
    const auto& curNodeScores = node.second;
    size_t rowIndex = 0;
    for (size_t i = 0; i < curNodeScores.size(); ++i) {
      if (excludeReads[i]) continue;
      probs(rowIndex, colIndex) = curNodeScores[i].second;
      ++rowIndex;
    }
    nodes.push_back(node.first);
    ++colIndex;
  }
  std::cerr << "Finished noFilter: " << nodes.size() << " nodes" << std::endl;
}

std::unordered_set<std::string> get_nth_order_neighbors(Tree *T, const std::string& node, int n_order) {
  std::unordered_set<std::string> result;

  if (T->allNodes.find(node) == T->allNodes.end()) {
    return result;
  }

  std::queue<std::pair<std::string, int>> bfs_queue;
  std::unordered_set<std::string> visited;

  bfs_queue.push({node, 0});
  visited.insert(node);

  while (!bfs_queue.empty()) {
    auto [current_node, level] = bfs_queue.front();
    bfs_queue.pop();

    // If we've reached the desired level, add the node to the result
    if (level <= n_order) {
      result.insert(current_node);

      if (level < n_order) {
        // Add children to the queue
        for (Node* child : T->allNodes[current_node]->children) {
          if (visited.find(child->identifier) == visited.end()) {
            bfs_queue.push({child->identifier, level + 1});
            visited.insert(child->identifier);
          }
        }
      }
      // Add parent to the queue (if exists)
      if (T->allNodes[current_node]->parent && visited.find(T->allNodes[current_node]->parent->identifier) == visited.end()) {
        bfs_queue.push({T->allNodes[current_node]->parent->identifier, level + 1});
        visited.insert(T->allNodes[current_node]->parent->identifier);
      }
    }
  }

  return result;
  
}

void filter_method_1(
  std::vector<std::string>& nodes, Eigen::MatrixXd& probs, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
  const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
  const std::vector<bool>& lowScoreReads, const size_t& numLowScoreReads, const std::string& excludeNode, const std::vector<bool>& excludeReads, Tree *T, int n_order
) {
  std::cerr << "Filter method 1: filter out haplotypes that do not have a unique best read score" << std::endl;

  std::unordered_set<std::string> filteredNodes;
  size_t numExcludedReads = std::count(excludeReads.begin(), excludeReads.end(), true);

  std::cerr << "Excluding " << numExcludedReads << " reads with high number of duplicates" << std::endl;

  std::vector<std::vector<int32_t>> scoresMatrix(allScores.begin()->second.size() - numExcludedReads);
  for (auto& vec : scoresMatrix) vec.resize(allScores.size() - leastRecentIdenticalAncestors.size());

  std::vector<std::string> allNodes(allScores.size() - leastRecentIdenticalAncestors.size());
  size_t colIndex = 0;
  for (const auto& node : allScores) {
    if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
    allNodes[colIndex] = node.first;
    const auto& curNodeScores = node.second;
    size_t rowIndex = 0;
    for (size_t j = 0; j < curNodeScores.size(); ++j) {
      if (excludeReads[j]) continue;
      scoresMatrix[rowIndex][colIndex] = curNodeScores[j].first;
      ++rowIndex;
    }
    ++colIndex;
  }

  for (size_t i = 0; i < scoresMatrix.size(); ++i) {
    const auto& curReadScores = scoresMatrix[i];
    int32_t highestScore = 0;
    std::string highestScoringNode = "";
    for (size_t j = 0; j < curReadScores.size(); ++j) {
      if (curReadScores[j] == 0) continue;
      if (curReadScores[j] > highestScore) {
        highestScore = curReadScores[j];
        highestScoringNode = allNodes[j];
      } else if (curReadScores[j] == highestScore) {
        highestScoringNode = "";
      }
    }
    if (highestScoringNode != "") filteredNodes.insert(highestScoringNode);
  }

  std::unordered_set<std::string> nodes_seen;
  for (const auto& filteredNode : filteredNodes) {
    std::unordered_set<std::string> identicalSet;
    if (identicalSets.find(filteredNode) != identicalSets.end()) {
      identicalSet = identicalSets.at(filteredNode);
    } else {
      identicalSet.insert(filteredNode);
    }
    for (const auto& node : identicalSet) {
      std::unordered_set<std::string> nth_order_neighbors = get_nth_order_neighbors(T, node, n_order);
      for (const auto& neighbor : nth_order_neighbors) {
        std::string neighbor_leastRecentIdenticalAncestor = neighbor;
        if (leastRecentIdenticalAncestors.find(neighbor) != leastRecentIdenticalAncestors.end()) {
          neighbor_leastRecentIdenticalAncestor = leastRecentIdenticalAncestors.at(neighbor);
        }
        if (nodes_seen.find(neighbor_leastRecentIdenticalAncestor) == nodes_seen.end() &&
            leastRecentIdenticalAncestors.find(neighbor_leastRecentIdenticalAncestor) == leastRecentIdenticalAncestors.end()) {
          nodes.push_back(neighbor_leastRecentIdenticalAncestor);
          nodes_seen.insert(neighbor_leastRecentIdenticalAncestor);
        }
      }
    }
  } 

  probs.resize(scoresMatrix.size(), nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) {
    const auto& curNodeScores = allScores.find(nodes[i])->second;
    for (size_t j = 0; j < curNodeScores.size(); ++j) {
      probs(j, i) = curNodeScores[j].second;
    }
  }

  std::cerr << "Finished filter method 1: " << nodes.size() << " nodes" << std::endl;
}

void filter_method_2(
  std::vector<std::string>& nodes, Eigen::MatrixXd& probs, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
  const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::vector<bool>& lowScoreReads, const size_t& numLowScoreReads,
  const std::string& excludeNode, const std::vector<bool>& excludeReads, size_t minHighestScoreCount
) {
  std::cerr << "Filter method 2: filter out haplotypes that do not have a best read score" << std::endl;

  std::unordered_map<std::string, size_t> highestScoreCounts;
  size_t numExcludedReads = std::count(excludeReads.begin(), excludeReads.end(), true);

  std::cerr << "Excluding " << numExcludedReads << " reads with high number of duplicates" << std::endl;

  std::vector<std::vector<int32_t>> scoresMatrix(allScores.begin()->second.size() - numExcludedReads);
  for (auto& vec : scoresMatrix) vec.resize(allScores.size() - leastRecentIdenticalAncestors.size());

  std::vector<std::string> allNodes(allScores.size() - leastRecentIdenticalAncestors.size());
  size_t colIndex = 0;
  for (const auto& node : allScores) {
    if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
    allNodes[colIndex] = node.first;
    const auto& curNodeScores = node.second;
    size_t rowIndex = 0;
    for (size_t j = 0; j < curNodeScores.size(); ++j) {
      if (excludeReads[j]) continue;
      scoresMatrix[rowIndex][colIndex] = curNodeScores[j].first;
      ++rowIndex;
    }
    ++colIndex;
  }

  
  for (size_t i = 0; i < scoresMatrix.size(); ++i) {
    const auto& curReadScores = scoresMatrix[i];
    int32_t highestScore = 0;
    std::vector<std::string> highestScoringNodes;
    for (size_t j = 0; j < curReadScores.size(); ++j) {
      if (curReadScores[j] == 0) continue;
      if (curReadScores[j] > highestScore) {
        highestScore = curReadScores[j];
        highestScoringNodes.clear();
        highestScoringNodes.push_back(allNodes[j]);
      } else if (curReadScores[j] == highestScore) {
        highestScoringNodes.push_back(allNodes[j]);
        }
      }
      for (const auto& node : highestScoringNodes) {
        ++highestScoreCounts[node];
    }
  }

  for (const auto& node : highestScoreCounts) {
    if (node.second >= minHighestScoreCount) {
      nodes.push_back(node.first);
    }
  }

  probs.resize(scoresMatrix.size(), nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) {
    const auto& curNodeScores = allScores.find(nodes[i])->second;
    for (size_t j = 0; j < curNodeScores.size(); ++j) {
      probs(j, i) = curNodeScores[j].second;
    }
  }

  std::cerr << "Finished filter method 2: " << nodes.size() << " nodes" << std::endl;
}


void filter_method_3(
  std::vector<std::string>& nodes, Eigen::MatrixXd& probs, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
  const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::vector<bool>& lowScoreReads, const size_t& numLowScoreReads,
  const std::string& excludeNode, const std::vector<bool>& excludeReads, Tree *T, int n_order
) {
  std::cerr << "Filter method 3: select nodes that are local optima and their n-th order neighbors" << std::endl;

  std::unordered_set<std::string> filteredNodes;
  size_t numExcludedReads = std::count(excludeReads.begin(), excludeReads.end(), true);

  std::cerr << "Excluding " << numExcludedReads << " reads with high number of duplicates" << std::endl;

  std::vector<std::vector<int32_t>> scoresMatrix(allScores.begin()->second.size() - numExcludedReads);
  for (auto& vec : scoresMatrix) vec.resize(allScores.size());

  std::vector<std::string> allNodes(allScores.size());
  size_t colIndex = 0;
  for (const auto& node : allScores) {
    allNodes[colIndex] = node.first;
    const auto& curNodeScores = node.second;
    size_t rowIndex = 0;
    for (size_t j = 0; j < curNodeScores.size(); ++j) {
      if (excludeReads[j]) continue;
      scoresMatrix[rowIndex][colIndex] = curNodeScores[j].first;
      ++rowIndex;
    }
    ++colIndex;
  }

  std::unordered_map<std::string, size_t> highestScoreReadCounts;
  for (size_t i = 0; i < scoresMatrix.size(); ++i) {
    const auto& curReadScores = scoresMatrix[i];
    int32_t highestScore = 0;
    std::vector<size_t> highestScoringNodeIndices;
    for (size_t j = 0; j < curReadScores.size(); ++j) {
      if (curReadScores[j] == 0) continue;
      if (curReadScores[j] > highestScore) {
        highestScoringNodeIndices.clear();
        highestScoringNodeIndices.push_back(j);
        highestScore = curReadScores[j];
      } else if (curReadScores[j] == highestScore) {
        highestScoringNodeIndices.push_back(j);
      }
    }
    for (const auto& index : highestScoringNodeIndices) {
      ++highestScoreReadCounts[allNodes[index]];
    }
  }

  std::vector<std::string> localOptima;
  for (const auto& node : allNodes) {

  }

  probs.resize(scoresMatrix.size(), nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) {
    const auto& curNodeScores = allScores.find(nodes[i])->second;
    for (size_t j = 0; j < curNodeScores.size(); ++j) {
      probs(j, i) = curNodeScores[j].second;
    }
  }

  std::cerr << "Finished filter method 3: " << nodes.size() << " nodes" << std::endl;
}

}



#endif
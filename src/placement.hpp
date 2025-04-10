#pragma once

#include "capnp/list.h"
#include "index.capnp.h"
#include "panman.hpp"
#include "progress_state.hpp"
#include "state.hpp"
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tbb/concurrent_vector.h>
#include <unordered_map>
#include <utility>
#include <vector>

// Forward declaration for indexing 
namespace indexing {
void processNodesByLevel(
    const std::vector<panmanUtils::Node*>& nodes,
    std::function<void(panmanUtils::Node*)> processFunction);
}

namespace placement {

// Forward declarations
struct PlacementResult;

// External declarations
extern std::shared_ptr<PlacementProgressState> progress_state;

// Parameters for traversal - consolidated
struct TraversalParams {
  int k = 32;              // k-mer size
  int s = 16;              // syncmer parameter s
  int t = 0;               // t-syncmer parameter
  bool open = true;        // Whether to use open syncmers
  double scoreScale = 1.0; // Scaling factor for scores
};

// Track global state during placement
struct PlacementGlobalState {
  // Reference to seed mutations from the index
  ::capnp::List<::SeedMutations>::Reader perNodeSeedMutations;

  // Reference to gap mutations from the index
  ::capnp::List<::GapMutations>::Reader perNodeGapMutations;

  // Reference to node path information from the index
  ::capnp::List<::NodePathInfo>::Reader nodePathInfo;

  // Reference to block information from the index
  ::capnp::List<::BlockInfo>::Reader blockInfo;

  // Reference to ancestor matrix from the index
  ::capnp::List<::capnp::List<bool>>::Reader ancestorMatrix;

  // K-mer dictionary mapping from IDs to sequences
  std::unordered_map<uint32_t, std::string> kmerDictionary;

  // Seed frequency in reads (hash -> count)
  std::unordered_map<size_t, int64_t> seedFreqInReads;

  // Total number of seed occurrences in reads
  size_t totalReadSeedCount = 0;

  // Jaccard denominator (total unique seeds in reads)
  size_t jaccardDenominator = 0;

  // Flag to stop the traversal
  std::atomic<bool> stopTraversal{false};

  // Helper methods for node relationships
  bool isAncestor(const std::string &potentialAncestorId, const std::string &nodeId) const;
  uint32_t getNodeLevel(const std::string &nodeId) const;
  std::vector<int32_t> getNodeActiveBlocks(const std::string &nodeId) const;
};

// Consolidated score tracking and results
struct PlacementResult {
  // Hits-based results
  int64_t maxHitsInAnyGenome = 0;
  panmanUtils::Node *maxHitsNode = nullptr;
  std::vector<panmanUtils::Node *> tiedMaxHitsNodes;

  // Jaccard-based results
  double bestJaccardScore = 0.0;
  panmanUtils::Node *bestJaccardNode = nullptr;
  std::vector<panmanUtils::Node *> tiedJaccardNodes;

  // Cosine-based results
  double bestCosineScore = 0.0;
  panmanUtils::Node *bestCosineNode = nullptr;
  std::vector<panmanUtils::Node *> tiedCosineNodes;

  // Weighted result
  double bestWeightedScore = 0.0;
  panmanUtils::Node *bestWeightedNode = nullptr;
  std::vector<panmanUtils::Node *> tiedWeightedNodes;

  // Current tracking values
  int64_t hitsInThisGenome = 0;
  int64_t currentJaccardNumerator = 0;
  int64_t currentJaccardDenominator = 0;
  double currentCosineNumerator = 0.0;
  double currentCosineDenominator = 0.0;
  std::unordered_map<size_t, int64_t> currentGenomeSeedCounts;

  // Performance metrics
  int64_t totalReadsProcessed = 0;
  double totalTimeSeconds = 0.0;
  
  // Helper methods
  void updateHitsScore(panmanUtils::Node* node, int64_t hits);
  void updateJaccardScore(panmanUtils::Node* node, double score);
  void updateCosineScore(panmanUtils::Node* node, double score);
  void updateWeightedScore(panmanUtils::Node* node, double score, double scale);
};

// Core functions for placement
void processNodeMutations(panmanUtils::Node *node,
                          state::StateManager &stateManager,
                          PlacementGlobalState &state,
                          PlacementResult &result,
                          const TraversalParams &params);

void placementTraversal(state::StateManager &stateManager,
                        PlacementResult &result,
                        panmanUtils::Tree *T, 
                        PlacementGlobalState &state,
                        const TraversalParams &params);

void place(PlacementResult &result, panmanUtils::Tree *T,
           ::Index::Reader &index, const std::string &reads1Path,
           const std::string &reads2Path,
           std::string &placementFileName);

void placeBatch(panmanUtils::Tree *T, ::Index::Reader &index,
                const std::string &batchFilePath,
                std::string prefixBase, std::string refFileNameBase,
                std::string samFileNameBase, std::string bamFileNameBase,
                std::string mpileupFileNameBase, std::string vcfFileNameBase,
                std::string aligner, const std::string &refNode,
                const bool &save_jaccard, const bool &show_time,
                const float &score_proportion, const int &max_tied_nodes);

// Helper function
std::pair<double, double> getCosineDelta(
    bool isRemoval, bool isAddition, size_t seedHash, int64_t count,
    const std::unordered_map<size_t, size_t>& readSeedCounts,
    const std::unordered_map<size_t, int64_t>& genomeSeedCounts);

} // namespace placement
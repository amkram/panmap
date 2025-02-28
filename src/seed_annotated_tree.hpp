#ifndef __TREE_HPP
#define __TREE_HPP

#include "coordinates.hpp"
#include "gap_map.hpp"
#include "index.capnp.h"
#include "panmanUtils.hpp"
#include "performance.hpp"
#include "seeding.hpp"
#include "timing.hpp"
#include <atomic>
#include <capnp/common.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <capnp/serialize.h>
#include <iostream>
#include <mutex>
#include <tuple>
#include <unordered_map>
#include <vector>

double time_stamp();

using namespace seeding;

// Forward declare types from coordinates namespace
namespace coordinates {
struct tupleCoord_t;
struct CoordRange;
using GapRange = std::pair<int64_t, int64_t>;
using GapUpdate = std::pair<bool, GapRange>;
using GapMap = std::map<int64_t, int64_t>;
using sequence_t = std::vector<
    std::pair<std::vector<std::pair<char, std::vector<char>>>,
              std::vector<std::vector<std::pair<char, std::vector<char>>>>>>;
} // namespace coordinates

// Now use the coordinate types
using coordinates::CoordRange;
using coordinates::GapMap;
using coordinates::GapRange;
using coordinates::GapUpdate;
using coordinates::sequence_t;
using coordinates::tupleCoord_t;

// Parameters for indexing and placement
struct PanmapParams {
  int32_t k; // k-mer length
  int32_t s; // syncmer length
  int32_t t; // threshold
  bool open; // open/closed syncmers
  int32_t l; // minimum length
};

inline auto seed_cmp = [](const std::pair<int32_t, std::string> &a,
                          const std::pair<int32_t, std::string> &b) {
  if (a.second != b.second) {
    return a.second < b.second;
  }
  return a.first < b.first;
};

typedef std::unordered_map<
    std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>>
    seedmerIndex_t;

using namespace panmanUtils;

typedef std::vector<std::tuple<int32_t, int32_t, bool, bool, bool, bool>>
    blockMutationInfo_t;
typedef std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, char, char>>
    mutationInfo_t;

namespace seed_annotated_tree {

// Helper method that uses a hybrid approach - scalar logic within blocks,
// traversal order at boundaries
int getValidNucleotidesEfficiently(coordinates::CoordinateManager &manager,
                                   int64_t startPos, char *buffer, int maxChars,
                                   int64_t &resultEndPos,
                                   const blockExists_t &blockExists,
                                   const blockStrand_t &blockStrand,
                                   const sequence_t &sequence);

inline void
fillDfsIndexes(panmanUtils::Tree *T, panmanUtils::Node *node, int64_t &dfsIndex,
               std::unordered_map<std::string, int64_t> &dfsIndexes) {
  TIME_FUNCTION;
  dfsIndexes[node->identifier] = dfsIndex;
  dfsIndex++;
  // enforce an ordering of children (alphabetical)
  std::set<Node *, std::function<bool(const Node *, const Node *)>>
      orderedChildren(node->children.begin(), node->children.end(),
                      [](const Node *a, const Node *b) {
                        return a->identifier < b->identifier;
                      });
  for (panmanUtils::Node *child : orderedChildren) {
    fillDfsIndexes(T, child, dfsIndex, dfsIndexes);
  }
}

// Struct to hold seed information
struct SeedInfo {
  size_t hash;
  int64_t endPos;
  bool isReverse;
};

// Type for seed change: {pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse,
// newIsReverse, oldEndPos, newEndPos}
using SeedChange =
    std::tuple<int64_t, bool, bool, std::optional<size_t>,
               std::optional<size_t>, std::optional<bool>, std::optional<bool>,
               std::optional<int64_t>, std::optional<int64_t>>;

struct SeedChangeComparator { // for set of seed changes, ordered by scalar
                              // position
  bool operator()(const SeedChange &a, const SeedChange &b) const {
    return std::get<0>(a) < std::get<0>(b);
  }
};

// Common state tracking for both indexing and placement
struct CommonTraversalState {
  std::vector<std::optional<seeding::onSeedsHash>> onSeedsHash;
  std::vector<std::unordered_set<int64_t>> BlocksToSeeds;
  std::unordered_map<std::string, int64_t> dfsIndexes;
  std::unordered_set<int> inverseBlockIds;

  CommonTraversalState(panmanUtils::Tree *T, size_t numBlocks,
                       size_t numCoords) {
    initialize(T, numBlocks, numCoords);
  }

  void initialize(panmanUtils::Tree *T, size_t numBlocks, size_t numCoords) {
    onSeedsHash.clear();
    onSeedsHash.resize(numCoords);
    BlocksToSeeds.clear();
    BlocksToSeeds.resize(numBlocks);
    dfsIndexes.clear();
    inverseBlockIds.clear();
    int64_t dfsIndex = 0;
    fillDfsIndexes(T, T->root, dfsIndex, dfsIndexes);
  }
};

// Common local state used during node processing
struct CommonNodeState {
  // Mutation tracking
  coordinates::blockExists_t oldBlockExists;
  coordinates::blockStrand_t oldBlockStrand;

  // Pre-allocated mutation buffers
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t nucleotideMutationInfo;
  std::vector<coordinates::CoordRange> recompRanges;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapMapUpdates;
  std::vector<std::pair<bool, int>> inverseBlockIdsBacktrack;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>
      gapRunBlocksBacktracks;

  CommonNodeState(size_t numBlocks, size_t numCoords) {
    // Initialize mutation tracking vectors
    oldBlockExists.resize(numBlocks, {false, {}});
    oldBlockStrand.resize(numBlocks, {true, {}});

    // Reserve space for mutation buffers
    blockMutationInfo.reserve(numBlocks);
    nucleotideMutationInfo.reserve(numCoords);
    recompRanges.reserve(numBlocks);
    gapRunUpdates.reserve(numBlocks);
    gapRunBacktracks.reserve(numBlocks);
    gapMapUpdates.reserve(numBlocks);
    inverseBlockIdsBacktrack.reserve(numBlocks);
    gapRunBlocksBacktracks.reserve(numBlocks);
  }
};

// Indexing-specific state
struct TraversalGlobalState : public CommonTraversalState {
  using CommonTraversalState::CommonTraversalState;
};

struct NodeLocalState : public CommonNodeState {
  NodeLocalState(size_t numBlocks, size_t numCoords)
      : CommonNodeState(numBlocks, numCoords) {}
};

// Placement-specific state
struct PlacementGlobalState : public CommonTraversalState {
  // Additional placement-specific fields
  int64_t maxHitsInAnyGenome = 0;
  double bestJaccardScore = 0.0;
  double bestCosineScore = 0.0;

  PlacementGlobalState(panmanUtils::Tree *T, size_t numBlocks, size_t numCoords)
      : CommonTraversalState(T, numBlocks, numCoords) {}
};

struct PlacementNodeState : public CommonNodeState {
  // Additional placement-specific fields
  int64_t hitsInThisGenome = 0;
  double cosineNumerator = 0.0;
  double cosineSumOfSquares = 0.0;

  PlacementNodeState(size_t numBlocks, size_t numCoords)
      : CommonNodeState(numBlocks, numCoords) {}
};

void applyMutations(blockMutationInfo_t &blockMutationInfo,
                    std::vector<coordinates::CoordRange> &recompRanges,
                    mutationInfo_t &mutationInfo,
                    std::vector<gap_map::GapUpdate> &gapRunUpdates,
                    std::vector<gap_map::GapUpdate> &gapRunBacktracks,
                    std::vector<gap_map::GapUpdate> &gapMapUpdates,
                    panmanUtils::Tree *T, panmanUtils::Node *node,
                    coordinates::CoordinateTraverser &traverser,
                    bool isPlacement, std::unordered_set<int> &inverseBlockIds,
                    std::vector<std::pair<bool, int>> &inverseBlockIdsBacktrack,
                    int k);

struct mutationMatrices {
  // Store mutation matrices
  std::vector<std::vector<double>> submat; // 4 x 4 substitution rate matrix
  std::unordered_map<int64_t, double>
      insmat; // 1 x N insertion rate by length matrix
  std::unordered_map<int64_t, double>
      delmat; // 1 x N deletion rate by length matrix

  bool filled = false;

  double maxInsLogProb = 40;
  double maxDelLogProb = 40;

  mutationMatrices() {
    // initialize mutationMatrices object and intialize the correct size for
    // substitution amtrix
    // total_submuts.resize(4);
    submat.resize(4);
    for (size_t i = 0; i < 4; ++i) {
      submat[i].resize(4);
    }
  }
};

/* Interface */

void setupIndexing(sequence_t &sequence, blockExists_t &blockExists,
                   blockStrand_t &blockStrand, const Tree *T);

void setupPlacement(
    std::vector<std::optional<seeding::onSeedsHash>> &onSeedsHash,
    sequence_t &sequence, blockExists_t &blockExists,
    blockStrand_t &blockStrand, const Tree *T);

bool getSeedAt(coordinates::CoordinateTraverser &traverser,
               coordinates::CoordinateManager &manager, size_t &resultHash,
               bool &resultIsReverse, int64_t &resultEndPos, const int64_t &pos,
               Tree *T, const int32_t &k);

/**
 * @brief Process multiple seed positions simultaneously to improve throughput
 * @param traverser Coordinate traverser
 * @param manager Coordinate manager
 * @param positions Array of positions to process
 * @param positionCount Number of positions in the array
 * @param k K-mer length
 * @param resultHashes Output hashes
 * @param resultIsReverse Output orientation flags
 * @param resultEndPos Output end positions
 * @param validResults Output validity flags
 */
void getSeedsBatch(coordinates::CoordinateTraverser &traverser,
                   coordinates::CoordinateManager &manager,
                   const int64_t *positions, // Input positions
                   size_t positionCount,     // Number of positions
                   int32_t k,                // k-mer length
                   size_t *resultHashes,     // Output hashes
                   bool *resultIsReverse,    // Output orientation flags
                   int64_t *resultEndPos,    // Output end positions
                   bool *validResults        // Output validity flags
);

// Fill mutation matrices from tree or file
std::pair<size_t, size_t> getMaskCoorsForMutmat(const std::string &s1,
                                                const std::string &s2,
                                                size_t window,
                                                double threshold);
// void fillMutationMatricesFromTree(mutationMatrices &mutMat, Tree *T,
//                                   size_t window, double threshold);
void fillMutationMatricesFromTree_test(mutationMatrices &mutMat, Tree *T,
                                       const std::string &path);
void fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf);

// Build mutation matrices by traversing through all parent-child pairs
void writeMutationMatrices(const mutationMatrices &mutMat,
                           std::ofstream &mmfout);

void getNucleotideSequenceFromBlockCoordinates(
    std::string &seq, std::vector<int64_t> &coords, std::vector<int64_t> &gaps,
    std::vector<int32_t> &deadBlocks, coordinates::CoordinateManager &manager,
    int64_t start_scalar, int64_t stop_scalar, const Tree *T, const Node *node);

std::string getStringAtNode(Node *node, Tree *T, bool aligned);

void undoMutations(Tree *T, const Node *node,
                   const blockMutationInfo_t &blockMutationInfo,
                   const mutationInfo_t &mutationInfo,
                   coordinates::CoordinateTraverser &traverser,
                   std::vector<gap_map::GapUpdate> &gapRunBacktracks);

// Common seed processing functions
inline void processSeedChange(
    CommonTraversalState &state, int64_t pos, bool oldVal, bool newVal,
    std::optional<size_t> oldSeed, std::optional<size_t> newSeed,
    std::optional<bool> oldIsReverse, std::optional<bool> newIsReverse,
    std::optional<int64_t> oldEndPos, std::optional<int64_t> newEndPos,
    coordinates::CoordinateTraverser &traverser) {

  if (oldVal && newVal && oldSeed != newSeed) {
    // Update existing seed
    state.onSeedsHash[pos] = {newSeed.value(), newEndPos.value(),
                              newIsReverse.value()};
  } else if (oldVal && !newVal) {
    // Remove seed
    state.onSeedsHash[pos].reset();
    int blockId = traverser.getCoordManager().getBlockIdOfScalarCoord(pos);
    state.BlocksToSeeds[blockId].erase(pos);
  } else if (!oldVal && newVal) {
    // Add new seed
    state.onSeedsHash[pos] = {newSeed.value(), newEndPos.value(),
                              newIsReverse.value()};
    int blockId = traverser.getCoordManager().getBlockIdOfScalarCoord(pos);
    state.BlocksToSeeds[blockId].insert(pos);
  }
}

inline void undoSeedChange(
    CommonTraversalState &state, int64_t pos, bool oldVal, bool newVal,
    std::optional<size_t> oldSeed, std::optional<size_t> newSeed,
    std::optional<bool> oldIsReverse, std::optional<bool> newIsReverse,
    std::optional<int64_t> oldEndPos, std::optional<int64_t> newEndPos,
    coordinates::CoordinateTraverser &traverser) {

  if (oldVal && newVal) {
    // Restore previous seed state
    state.onSeedsHash[pos] = {oldSeed.value(), oldEndPos.value(),
                              oldIsReverse.value()};
  } else if (oldVal && !newVal) {
    // Restore deleted seed
    state.onSeedsHash[pos] = {oldSeed.value(), oldEndPos.value(),
                              oldIsReverse.value()};
    int blockId = traverser.getCoordManager().getBlockIdOfScalarCoord(pos);
    state.BlocksToSeeds[blockId].insert(pos);
  } else if (!oldVal && newVal) {
    // Remove added seed
    state.onSeedsHash[pos].reset();
    int blockId = traverser.getCoordManager().getBlockIdOfScalarCoord(pos);
    state.BlocksToSeeds[blockId].erase(pos);
  }
}
} // namespace seed_annotated_tree

using namespace seed_annotated_tree;
/**
 * @brief Template specializations for fixed k-mer sizes
 * These optimized implementations handle specific k-mer sizes more efficiently
 */
namespace fixed_kmer {
// Specialized function for small k-mers (k <= 16)
template <int K>
inline bool getSeedAtFixed(coordinates::CoordinateTraverser &traverser,
                           coordinates::CoordinateManager &manager,
                           size_t &resultHash, bool &resultIsReverse,
                           int64_t &resultEndPos, const int64_t &pos, Tree *T);

// Specialized batch processing for small k-mers
template <int K>
inline void getSeedsBatchFixed(coordinates::CoordinateTraverser &traverser,
                               coordinates::CoordinateManager &manager,
                               const int64_t *positions, size_t positionCount,
                               size_t *resultHashes, bool *resultIsReverse,
                               int64_t *resultEndPos, bool *validResults);

// Helper for selecting the correct specialized implementation based on k-mer
// size
inline bool dispatchGetSeedAt(coordinates::CoordinateTraverser &traverser,
                              coordinates::CoordinateManager &manager,
                              size_t &resultHash, bool &resultIsReverse,
                              int64_t &resultEndPos, const int64_t &pos,
                              Tree *T, const int32_t &k);

// Helper for selecting the correct specialized batch implementation
inline void dispatchGetSeedsBatch(coordinates::CoordinateTraverser &traverser,
                                  coordinates::CoordinateManager &manager,
                                  const int64_t *positions,
                                  size_t positionCount, int32_t k,
                                  size_t *resultHashes, bool *resultIsReverse,
                                  int64_t *resultEndPos, bool *validResults);

} // namespace fixed_kmer
#endif

#ifndef __TREE_HPP
#define __TREE_HPP

#include "coordinates.hpp"  // Include this file for full definition of CoordinateTraverser
#include "gap_map.hpp"
#include "index.capnp.h"
#include "seeding.hpp"
#include <capnp/common.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <capnp/serialize.h>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <unordered_set>

using namespace seeding;
using namespace panmanUtils;

// Forward declare types from coordinates namespace
namespace coordinates {
class CoordinateManager;
class CoordinateTraverser;
struct tupleCoord_t;
struct CoordRange;
using GapRange = std::pair<int64_t, int64_t>;
using GapUpdate = std::pair<bool, GapRange>;
using GapMap = std::map<int64_t, int64_t>;
using sequence_t = std::vector<
    std::pair<std::vector<std::pair<char, std::vector<char>>>,
              std::vector<std::vector<std::pair<char, std::vector<char>>>>>>;
using blockExists_t = std::vector<std::pair<bool, std::vector<bool>>>;
using blockStrand_t = std::vector<std::pair<bool, std::vector<bool>>>;
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


typedef std::vector<std::tuple<int32_t, int32_t, bool, bool, bool, bool>>
    blockMutationInfo_t;
typedef std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, char, char>>
    mutationInfo_t;

namespace seed_annotated_tree {

// Global debug flag to control verbose output
extern bool debug;

// Helper method that uses a hybrid approach - scalar logic within blocks,
// traversal order at boundaries
int getValidNucleotidesEfficiently(coordinates::CoordinateManager &manager,
                                   int64_t startPos, char *buffer, int maxChars,
                                   int64_t &resultEndPos,
                                   const coordinates::blockExists_t &blockExists,
                                   const coordinates::blockStrand_t &blockStrand,
                                   const sequence_t &sequence);

inline void
fillDfsIndexes(panmanUtils::Tree *T, panmanUtils::Node *node, int64_t &dfsIndex,
               std::unordered_map<std::string, int64_t> &dfsIndexes) {

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

class TraversalNodeState {
public:
  blockExists_t oldBlockExists;
  blockStrand_t oldBlockStrand;
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t nucleotideMutationInfo;
  std::vector<coordinates::GapUpdate> gapRunUpdates;
  std::vector<coordinates::GapUpdate> gapRunBacktracks;
  std::vector<coordinates::GapUpdate> gapMapUpdates;
  std::vector<coordinates::GapUpdate> gapRunBlocksBacktracks;
  std::vector<coordinates::CoordRange> recompRanges;
  std::vector<std::pair<bool, int>> inverseBlockIdsBacktrack;

  TraversalNodeState(int num_blocks, int64_t num_scalars) {
    blockMutationInfo.clear();
    nucleotideMutationInfo.clear();
    gapRunUpdates.clear();
    gapRunBacktracks.clear();
    gapMapUpdates.clear();
    gapRunBlocksBacktracks.clear();
    recompRanges.clear();
    inverseBlockIdsBacktrack.clear();
  }

  TraversalNodeState() {}
};

// For indexing - same as the base for now
using IndexingNodeState = TraversalNodeState;

// NodeLocalState should be compatible with or derive from TraversalNodeState
using NodeLocalState = TraversalNodeState;

// For placement - extends the base with seed change tracking
class PlacementNodeState : public TraversalNodeState {
public:
  // Add seed changes tracking for placement
  std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>> seedChanges;

  // Flag to track if node has been modified
  bool isDirty = false;
  
  // Container for tie-breaking information
  std::unique_ptr<std::unordered_map<std::string, double>> tieBreakingInfo;

  PlacementNodeState(int num_blocks, int64_t num_scalars) 
    : TraversalNodeState(num_blocks, num_scalars) {
    seedChanges.clear();
    tieBreakingInfo = std::make_unique<std::unordered_map<std::string, double>>();
  }

  PlacementNodeState() : TraversalNodeState() {}

  // Copy constructor
  PlacementNodeState(const PlacementNodeState &other)
      : TraversalNodeState(other) {
    isDirty = other.isDirty;
    if (other.tieBreakingInfo) {
      tieBreakingInfo = std::make_unique<std::unordered_map<std::string, double>>(*other.tieBreakingInfo);
    }
  }

  // Move constructor
  PlacementNodeState(PlacementNodeState &&other) noexcept
      : TraversalNodeState(std::move(other)) {
    isDirty = other.isDirty;
    tieBreakingInfo = std::move(other.tieBreakingInfo);
  }

  // Virtual destructor
  virtual ~PlacementNodeState() {}
};

// Indexing-specific state
struct TraversalGlobalState : public CommonTraversalState {
  using CommonTraversalState::CommonTraversalState;
};

class PlacementGlobalState : public TraversalGlobalState {
public:
  // Collections for seed tracking during placement
  std::unordered_map<size_t, std::pair<size_t, size_t>> perNodeReadSeedCounts;
  std::unordered_map<size_t, int64_t> perNodeGenomeSeedCounts;
  size_t totalReadSeedCount = 0;
  
  // Access to index data
  ::capnp::List<GapMutations>::Reader perNodeGapMutations;
  ::capnp::List<SeedMutations>::Reader perNodeSeedMutations;
  
  // Inverse block IDs for handling inversions
  std::unordered_set<int32_t> inverseBlockIds;
  
  // K-mer size for seed extraction
  int32_t kmerSize = 31; // Default value

  PlacementGlobalState(panmanUtils::Tree *tree, int32_t numBlocks,
                      int64_t numCoords)
      : TraversalGlobalState(tree, numBlocks, numCoords) {
    // Initialize placement-specific data structures
    perNodeReadSeedCounts.clear();
    perNodeGenomeSeedCounts.clear();
    totalReadSeedCount = 0;
    inverseBlockIds.clear();
  }

  // Method to set index data references
  void setIndexData(const ::capnp::List<GapMutations>::Reader &gapMutations,
                   const ::capnp::List<SeedMutations>::Reader &seedMutations) {
    perNodeGapMutations = gapMutations;
    perNodeSeedMutations = seedMutations;
  }
};

void applyMutations(blockMutationInfo_t &blockMutationInfo,
                    std::vector<coordinates::CoordRange> &recompRanges,
                    mutationInfo_t &mutationInfo,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapMapUpdates,
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
                   std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks);

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

// Add declaration for buildMutationMatricesHelper_test
void buildMutationMatricesHelper_test(
    mutationMatrices &mutMat, panmanUtils::Tree *T, panmanUtils::Node *node,
    std::map<int64_t, int64_t> &gapMap,
    coordinates::CoordinateTraverser &traverser,
    std::vector<int64_t> &scalarCoordToBlockId,
    std::vector<std::unordered_set<int>> &BlocksToSeeds,
    std::vector<int> &BlockSizes,
    const std::vector<std::pair<int64_t, int64_t>> &blockRanges,
    std::vector<int64_t> &parentBaseCounts,
    std::vector<int64_t> &totalBaseCounts,
    std::vector<std::vector<int64_t>> &subCount,
    std::unordered_map<int64_t, int64_t> &insCount,
    std::unordered_map<int64_t, int64_t> &delCount);

// Declare shared functions for mutation application and backtracking
void applyTreeNodeMutations(
  TraversalNodeState &nodeState,
  coordinates::CoordinateTraverser &traverser,
  panmanUtils::Tree *T, 
  panmanUtils::Node *node,
  bool isPlacement,
  std::unordered_set<int32_t> &inverseBlockIds,
  int k);

void backtrackNodeState(
  panmanUtils::Tree *T, 
  panmanUtils::Node *node, 
  TraversalNodeState &nodeState,
  coordinates::CoordinateTraverser &traverser,
  bool undoSeedChanges = false,
  std::unordered_map<size_t, int64_t> *currentGenomeSeedCounts = nullptr,
  std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>> *seedChanges = nullptr);

void processGapMutationsForPlacement(
  const ::capnp::List<GapMutations>::Reader &gapMutationsList,
  PlacementNodeState &nodeState,
  coordinates::CoordinateManager &manager);

// Helper for processing seed mutations in placement
void processSeedMutationsForPlacement(
  PlacementGlobalState &state,
  PlacementNodeState &nodeState,
  coordinates::CoordinateTraverser &traverser,
  const ::capnp::List<SeedMutations>::Reader &seedIndex,
  int dfsIndex,
  std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
  std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
  int64_t &hitsInThisGenome);

// Process a single gap mutation record
void processGapMutationForPlacement(
    const GapMutations::Reader &gapMutation,
    PlacementNodeState &nodeState,
    coordinates::CoordinateManager &coordManager);

// Process a single seed mutation record
void processSeedMutationForPlacement(
    PlacementGlobalState &state,
    PlacementNodeState &nodeState,
    coordinates::CoordinateTraverser &coordTraverser,
    const SeedMutations::Reader &seedMutation,
    int64_t dfsIndex,
    std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> &seedToNodeMap,
    std::unordered_map<uint64_t, int64_t> &seedToDfsMap,
    int64_t &seedCount);

// Shared utility functions for gap map operations
void initializeGapMapFromSequence(
  coordinates::CoordinateTraverser &traverser,
  coordinates::CoordinateManager &manager);
  
bool isSuspiciousGap(
  int64_t position, 
  int64_t length, 
  int64_t sequenceSize);
  
bool validateAndFixGapMap(coordinates::CoordinateManager &manager, 
                         const std::string &context);

// Process runs of blocks with the same state (on/off/inverted) efficiently
void processBlockRun(
    std::vector<CoordRange> &recompRanges,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
    std::vector<std::pair<int64_t, int64_t>> &gapRunBacktracks,
    std::vector<std::pair<size_t, std::pair<bool, int64_t>>> &gapMapUpdates,
    coordinates::CoordinateTraverser &traverser,
    int32_t startBlockId,
    int32_t endBlockId,
    bool isPlacement);

const int DEFAULT_SHORT_RECOMP_SIZE = 21;
const int DEFAULT_GAP_RECOMP_SIZE = 5;
const int MIN_BLOCK_RUN_SIZE = 2; // Minimum number of blocks to consider as a "run" for optimization

// Add declaration for validateGapMapAfterMutations
bool validateGapMapAfterMutations(
  coordinates::CoordinateTraverser &traverser,
  panmanUtils::Node *node);

} // namespace seed_annotated_tree
#endif

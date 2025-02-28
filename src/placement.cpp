#include "placement.hpp"
#include "alignment.hpp"
#include "fixed_kmer.hpp"
#include "genotyping.hpp"
#include "indexing.hpp"
#include "logging.hpp"
#include "performance.hpp"
#include "seed_annotated_tree.hpp"
#include "timing.hpp"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_group.h>

// Add FTXUI includes for dynamic progress display
#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>

namespace fs = boost::filesystem;

using namespace seed_annotated_tree;
using namespace coordinates; // For coordinate types and traversal
using namespace logging;

namespace placement {

// Structure to hold the state of the UI
struct PlacementProgressState {
  std::atomic<bool> running{true};
  std::mutex mtx;
  std::string currentNodeId = "none";
  std::string maxNodeId = "none";
  int64_t maxHits = 0;
  int64_t nodesVisited = 0;
  std::vector<std::string> tiedNodes;
};

// Shared state for tracking progress
std::shared_ptr<PlacementProgressState> progress_state;

// Function to render the placement progress UI
ftxui::Component CreatePlacementProgressUI() {
  auto renderer = ftxui::Renderer([&] {
    std::lock_guard<std::mutex> lock(progress_state->mtx);

    // Create elements for tied nodes if any
    ftxui::Elements tied_elements;
    if (!progress_state->tiedNodes.empty()) {
      for (const auto &node : progress_state->tiedNodes) {
        tied_elements.push_back(ftxui::text(" • " + node));
      }
    } else {
      tied_elements.push_back(ftxui::text(" • none") | ftxui::dim);
    }

    return ftxui::vbox({
               ftxui::hbox(
                   {ftxui::text("Panmap Placement Progress") | ftxui::bold |
                        ftxui::color(ftxui::Color::BlueLight),
                    ftxui::filler(),
                    ftxui::text("Nodes visited: ") | ftxui::dim,
                    ftxui::text(std::to_string(progress_state->nodesVisited)) |
                        ftxui::bold}),
               ftxui::separator(),
               ftxui::hbox({ftxui::text("Current node: "),
                            ftxui::text(progress_state->currentNodeId) |
                                ftxui::color(ftxui::Color::Yellow)}),
               ftxui::hbox({ftxui::text("Best node:    "),
                            ftxui::text(progress_state->maxNodeId) |
                                ftxui::bold |
                                ftxui::color(ftxui::Color::GreenLight)}),
               ftxui::hbox(
                   {ftxui::text("Max hits:     "),
                    ftxui::text(std::to_string(progress_state->maxHits)) |
                        ftxui::bold}),
               ftxui::separator(),
               ftxui::text("Tied nodes:") | ftxui::bold,
               ftxui::vbox(tied_elements),
           }) |
           ftxui::border | ftxui::color(ftxui::Color::White);
  });

  return renderer;
}

std::pair<double, double> getCosineDelta(
    bool oldVal, bool newVal, size_t oldSeedVal, size_t newSeedVal,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts) {
  TIME_FUNCTION;
  auto getFreqs = [&](size_t seedVal) {
    double freqInReads = 0.0;
    auto readIt = readSeedCounts.find(seedVal);
    if (readIt != readSeedCounts.end()) {
      freqInReads = std::log2(1.0 + readIt->second.first +
                              readIt->second.second); // Set default frequency
    }

    double freqInCurrGenome = 0.0;
    auto currGenomeIt = currentGenomeSeedCounts.find(seedVal);
    if (currGenomeIt != currentGenomeSeedCounts.end()) {
      freqInCurrGenome =
          std::log2(1.0 + currGenomeIt->second); // Initialize with actual count
    }

    return std::make_pair(freqInReads, freqInCurrGenome);
  };

  double numeratorDelta = 0.0;
  double sumOfSquaresDelta = 0.0;

  if (!oldVal && newVal) { // Adding a new seed
    auto [readWeight, genomeWeight] = getFreqs(newSeedVal);
    numeratorDelta = readWeight * genomeWeight;
    sumOfSquaresDelta = std::pow(genomeWeight, 2);
  } else if (oldVal && !newVal) { // Removing a seed
    auto [readWeight, genomeWeight] = getFreqs(oldSeedVal);
    numeratorDelta = -readWeight * genomeWeight;
    sumOfSquaresDelta = -std::pow(genomeWeight, 2);
  } else if (oldVal && newVal) { // Changing a seed
    auto [oldReadWeight, oldGenomeWeight] = getFreqs(oldSeedVal);
    auto [newReadWeight, newGenomeWeight] = getFreqs(newSeedVal);

    numeratorDelta =
        newReadWeight * newGenomeWeight - oldReadWeight * oldGenomeWeight;
    sumOfSquaresDelta =
        std::pow(newGenomeWeight, 2) - std::pow(oldGenomeWeight, 2);
  }

  return {numeratorDelta, sumOfSquaresDelta};
}

// Specialized template version that's instantiated for specific k-mer sizes
template <int K>
void placementTraversalSpecialized(
    CoordinateTraverser &traverser, int64_t hitsInThisGenome,
    int64_t &maxHitsInAnyGenome, Node *&maxHitsNode,
    int64_t jaccardNumeratorThisGenome, int64_t jaccardDenominatorThisGenome,
    double &bestJaccardScore, Node *&bestJaccardNode,
    double cosineNumeratorThisGenome, double cosineSumOfSquaresThisGenome,
    double &bestCosineScore, Node *&bestCosineNode, Node *parent, Node *current,
    PlacementGlobalState &state, PlacementResult &result,
    ::capnp::List<SeedMutations>::Reader &seedIndex,
    ::capnp::List<GapMutations>::Reader &gapIndex, int seedK, int seedS,
    int seedT, bool open, int seedL, Tree *T,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
    const size_t &totalReadSeedCount, const size_t &numGenomes,
    blockExists_t &oldBlockExists, coordinates::blockStrand_t &oldBlockStrand,
    const std::string &true_node_id, const std::string &species,
    const int &read_count, const int &mutation_count, bool firstCall,
    bool &stopReached, Node *startNode, Node *stopNode,
    const std::unordered_set<Node *> &groupNodes, int64_t &searchCount) {

  TIME_FUNCTION;

  // Cache manager reference at start
  auto &manager = traverser.getCoordManager();
  auto &sequence = manager.getSequence();
  auto &blockExists = manager.getBlockExists();
  auto &blockStrand = manager.getBlockStrand();

  // Check if current node is in the group
  if (groupNodes.find(current) == groupNodes.end()) {
    msg("[DFS] Node {} not in group, stopping traversal branch",
        current->identifier);
    stopReached = true;
    return;
  }

  int64_t hitsBefore = hitsInThisGenome;
  int64_t dfsIndex = state.dfsIndexes[current->identifier];

  // Initialize block ranges if needed
  size_t ranges_initialized = 0;
  for (size_t blockId = 0; blockId < sequence.size(); blockId++) {
    auto range = manager.getBlockRange(blockId);
    if (range.first == 0 && range.second == 0 && blockExists[blockId].first) {
      ranges_initialized++;

      // Scan coordinates to find actual range
      int64_t start = -1, end = -1;
      for (int64_t i = 0; i < manager.size(); i++) {
        const auto *coord = manager.getTupleCoord(i);
        if (coord && coord->blockId == blockId) {
          if (start == -1)
            start = i;
          end = i;
        }
      }
      if (start != -1 && end != -1) {
        manager.updateBlockRange(blockId, start, end);
      }
    }
  }

  PlacementNodeState nodeState(traverser.getCoordManager().getNumBlocks(),
                               traverser.getCoordManager().getNumCoords());

  // Initialize mutation state
  nodeState.oldBlockExists = manager.getBlockExists();
  nodeState.oldBlockStrand = manager.getBlockStrand();
  nodeState.gapRunBacktracks.clear();

  // Apply mutations
  seed_annotated_tree::applyMutations(
      nodeState.blockMutationInfo, nodeState.recompRanges,
      nodeState.nucleotideMutationInfo, nodeState.gapRunUpdates,
      nodeState.gapRunBacktracks, nodeState.gapMapUpdates, T, current,
      traverser, true, state.inverseBlockIds,
      nodeState.inverseBlockIdsBacktrack,
      K); // Note: K is now known at compile time, not seedK

  // Update gap map with current mutations
  manager.updateGapMap(nodeState.gapRunUpdates, nodeState.gapRunBacktracks,
                       nodeState.gapMapUpdates);

  // Process gap mutations
  auto gapMutationsList = gapIndex[dfsIndex].getDeltas();
  for (const auto &gapMutation : gapMutationsList) {
    int32_t pos = gapMutation.getPos();
    auto maybeValue = gapMutation.getMaybeValue();
    if (maybeValue.isValue()) {
      auto it = manager.getGapMap().find(pos);
      if (it != manager.getGapMap().end()) {
        nodeState.gapRunBacktracks.emplace_back(
            false, std::make_pair(pos, it->second));
      } else {
        nodeState.gapRunBacktracks.emplace_back(
            true, std::make_pair(pos, maybeValue.getValue()));
      }
      manager.setGapMap(pos, maybeValue.getValue());
    } else {
      auto it = manager.getGapMap().find(pos);
      nodeState.gapRunBacktracks.emplace_back(false,
                                              std::make_pair(pos, it->second));
      manager.removeFromGapMap(pos);
    }
  }

  // <seedPos, seedHash, isReverse, endPos, seedCount>
  std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>> seedChanges;
  seedChanges.reserve(100); // Reasonable estimate of changes per node

  // Process seed mutations directly instead of storing in localSeedChanges
  auto currBasePositions = seedIndex[dfsIndex].getBasePositions();
  auto currPerPosMasks = seedIndex[dfsIndex].getPerPosMasks();

  // Using fixed-size buffer for better performance
  alignas(64) char seedBuffer[K];

  // Create batch processing variables for K-mer positions
  constexpr size_t BATCH_SIZE = 32;
  alignas(64) int64_t batchPositions[BATCH_SIZE];
  alignas(64) size_t batchHashes[BATCH_SIZE];
  alignas(64) bool batchIsReverse[BATCH_SIZE];
  alignas(64) int64_t batchEndPos[BATCH_SIZE];
  alignas(64) bool batchValid[BATCH_SIZE];

  for (int i = 0; i < currBasePositions.size(); ++i) {
    int64_t pos = currBasePositions[i];
    uint64_t tritMask = currPerPosMasks[i];

    // For each batch of positions to process
    size_t batchIdx = 0;

    for (int j = 0; j < 32; ++j) {
      uint8_t ternaryNumber = (tritMask >> (j * 2)) & 0x3;
      int64_t seedPos = pos - j;

      if (ternaryNumber == 1) { // on -> off
        if (!state.onSeedsHash[seedPos].has_value()) {
          continue;
        }

        // Get current values before changing
        auto [oldSeed, oldEndPos, oldIsReverse] =
            state.onSeedsHash[seedPos].value();

        // Track seed count changes directly
        size_t oldSeedVal = oldSeed;
        size_t oldSeedCount = currentGenomeSeedCounts[oldSeedVal];

        seedChanges.emplace_back(seedPos, oldSeedVal, oldIsReverse, oldEndPos,
                                 oldSeedCount);

        if (currentGenomeSeedCounts.find(oldSeedVal) ==
            currentGenomeSeedCounts.end()) {
          currentGenomeSeedCounts[oldSeedVal] = 0;
        }
        if (currentGenomeSeedCounts[oldSeedVal] > 0) {
          currentGenomeSeedCounts[oldSeedVal]--;
        }

        // Update metrics directly
        auto readIt = readSeedCounts.find(oldSeedVal);
        if (readIt != readSeedCounts.end() && oldSeedCount > 0) {
          // Decrement hits since we're removing a seed
          hitsInThisGenome--;
        }

        // Apply state change directly
        state.onSeedsHash[seedPos].reset();
        int blockId = manager.getBlockIdOfScalarCoord(seedPos);
        state.BlocksToSeeds[blockId].erase(seedPos);

      } else if (ternaryNumber == 2) { // off -> on or on -> on (change)
        // Add this position to the batch
        batchPositions[batchIdx++] = seedPos;

        // Process the batch if it's full or we're at the end
        if (batchIdx == BATCH_SIZE || j == 31) {
          // Process the batch of positions
          fixed_kmer::getSeedsBatchFixed<K>(
              traverser, manager, batchPositions, batchIdx, batchHashes,
              batchIsReverse, batchEndPos, batchValid);

          // Process the results
          for (size_t idx = 0; idx < batchIdx; idx++) {
            if (!batchValid[idx])
              continue;

            int64_t seedPos = batchPositions[idx];
            size_t newSeedHash = batchHashes[idx];
            bool newIsReverse = batchIsReverse[idx];
            int64_t newEndPos = batchEndPos[idx];

            if (!state.onSeedsHash[seedPos].has_value()) {
              seedChanges.emplace_back(seedPos, newSeedHash, newIsReverse,
                                       newEndPos, 0);

              // off -> on case
              // Track seed count changes directly
              size_t newSeedVal = newSeedHash;
              if (currentGenomeSeedCounts.find(newSeedVal) ==
                  currentGenomeSeedCounts.end()) {
                currentGenomeSeedCounts[newSeedVal] = 0;
              }
              currentGenomeSeedCounts[newSeedVal]++;

              // Update metrics directly
              size_t newSeedCount = currentGenomeSeedCounts[newSeedVal];
              auto readIt = readSeedCounts.find(newSeedVal);
              if (readIt != readSeedCounts.end() && newSeedCount > 0) {
                // Increment hits since we're adding a seed
                hitsInThisGenome++;
              }

              // Apply state change directly
              state.onSeedsHash[seedPos] = {newSeedHash, newEndPos,
                                            newIsReverse};
              int blockId =
                  traverser.getCoordManager().getBlockIdOfScalarCoord(seedPos);
              state.BlocksToSeeds[blockId].insert(seedPos);
            } else {
              // on -> on (change) case
              auto [oldSeed, oldEndPos, oldIsReverse] =
                  state.onSeedsHash[seedPos].value();

              // Skip if no actual change
              if (oldSeed == newSeedHash && oldIsReverse == newIsReverse &&
                  oldEndPos == newEndPos) {
                continue;
              }

              // Track seed count changes directly
              size_t oldSeedVal = oldSeed;
              size_t newSeedVal = newSeedHash;
              size_t oldSeedCount = currentGenomeSeedCounts[oldSeedVal];

              seedChanges.emplace_back(seedPos, oldSeedVal, oldIsReverse,
                                       oldEndPos, oldSeedCount);

              if (currentGenomeSeedCounts.find(oldSeedVal) ==
                  currentGenomeSeedCounts.end()) {
                currentGenomeSeedCounts[oldSeedVal] = 0;
              }
              if (currentGenomeSeedCounts[oldSeedVal] > 0) {
                currentGenomeSeedCounts[oldSeedVal]--;
              }

              if (currentGenomeSeedCounts.find(newSeedVal) ==
                  currentGenomeSeedCounts.end()) {
                currentGenomeSeedCounts[newSeedVal] = 0;
              }
              currentGenomeSeedCounts[newSeedVal]++;

              // Update metrics directly - for replacement, remove old and add
              // new
              auto readItOld = readSeedCounts.find(oldSeedVal);
              if (readItOld != readSeedCounts.end() && oldSeedCount == 0) {
                // Decrement hits for the old seed value
                hitsInThisGenome--;
              }

              auto readItNew = readSeedCounts.find(newSeedVal);
              if (readItNew != readSeedCounts.end() &&
                  currentGenomeSeedCounts[newSeedVal] == 1) {
                // Increment hits for the new seed value
                hitsInThisGenome++;
              }

              // Apply state change directly
              state.onSeedsHash[seedPos] = {newSeedHash, newEndPos,
                                            newIsReverse};
            }
          }

          // Reset batch index
          batchIdx = 0;
        }
      }
    }
  }

  // Check if this node has a better score
  if (hitsInThisGenome > maxHitsInAnyGenome) {
    maxHitsInAnyGenome = hitsInThisGenome;
    maxHitsNode = current;

    // Update progress state with new max node
    if (progress_state) {
      std::lock_guard<std::mutex> lock(progress_state->mtx);
      progress_state->maxNodeId = current->identifier;
      progress_state->maxHits = maxHitsInAnyGenome;
      progress_state->tiedNodes.clear();
      progress_state->tiedNodes.push_back(current->identifier);
    }
  }
  // Track nodes that have the same score (tied)
  else if (hitsInThisGenome == maxHitsInAnyGenome && maxHitsInAnyGenome > 0) {
    // Update progress state with tied node
    if (progress_state) {
      std::lock_guard<std::mutex> lock(progress_state->mtx);
      // Check if this node is already in the tied nodes list
      if (std::find(progress_state->tiedNodes.begin(),
                    progress_state->tiedNodes.end(),
                    current->identifier) == progress_state->tiedNodes.end()) {
        progress_state->tiedNodes.push_back(current->identifier);
      }
    }
  }

  nodeState.recompRanges.clear();
  nodeState.gapRunUpdates.clear();
  nodeState.gapMapUpdates.clear();

  // Recursively process children
  for (Node *child : current->children) {
    if (stopReached) {
      err("BREAKING PLACEMENT TRAVERSAL");
      break;
    }
    if (groupNodes.find(child) == groupNodes.end()) {
      err("BREAKING PLACEMENT TRAVERSAL");
      stopReached = true;
      break;
    }

    // Process child
    placementTraversalSpecialized<K>(
        traverser, hitsInThisGenome, maxHitsInAnyGenome, maxHitsNode,
        jaccardNumeratorThisGenome, jaccardDenominatorThisGenome,
        bestJaccardScore, bestJaccardNode, cosineNumeratorThisGenome,
        cosineSumOfSquaresThisGenome, bestCosineScore, bestCosineNode, current,
        child, state, result, seedIndex, gapIndex,
        K, // Use compile-time K instead of seedK
        seedS, seedT, open, seedL, T, readSeedCounts, currentGenomeSeedCounts,
        totalReadSeedCount, numGenomes, oldBlockExists, oldBlockStrand,
        true_node_id, species, read_count, mutation_count, false, stopReached,
        startNode, stopNode, groupNodes, searchCount);
  }

  backtrackNode(state, traverser, T, current, nodeState,
                currentGenomeSeedCounts, seedChanges,
                nodeState.gapRunBacktracks);

  // Update progress state with current node info
  if (progress_state) {
    std::lock_guard<std::mutex> lock(progress_state->mtx);
    progress_state->currentNodeId = current->identifier;
    progress_state->nodesVisited++;
  }
}

// Modify the main placementTraversal function to dispatch to specialized
// versions
void placementTraversal(
    CoordinateTraverser &traverser, int64_t hitsInThisGenome,
    int64_t &maxHitsInAnyGenome, Node *&maxHitsNode,
    int64_t jaccardNumeratorThisGenome, int64_t jaccardDenominatorThisGenome,
    double &bestJaccardScore, Node *&bestJaccardNode,
    double cosineNumeratorThisGenome, double cosineSumOfSquaresThisGenome,
    double &bestCosineScore, Node *&bestCosineNode, Node *parent, Node *current,
    PlacementGlobalState &state, PlacementResult &result,
    ::capnp::List<SeedMutations>::Reader &seedIndex,
    ::capnp::List<GapMutations>::Reader &gapIndex, int seedK, int seedS,
    int seedT, bool open, int seedL, Tree *T,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
    const size_t &totalReadSeedCount, const size_t &numGenomes,
    blockExists_t &oldBlockExists, coordinates::blockStrand_t &oldBlockStrand,
    const std::string &true_node_id, const std::string &species,
    const int &read_count, const int &mutation_count, bool firstCall,
    bool &stopReached, Node *startNode, Node *stopNode,
    const std::unordered_set<Node *> &groupNodes, int64_t &searchCount) {

  // Dispatch to specialized template implementations based on k-mer size
  switch (seedK) {
  case 8:
    return placementTraversalSpecialized<8>(
        traverser, hitsInThisGenome, maxHitsInAnyGenome, maxHitsNode,
        jaccardNumeratorThisGenome, jaccardDenominatorThisGenome,
        bestJaccardScore, bestJaccardNode, cosineNumeratorThisGenome,
        cosineSumOfSquaresThisGenome, bestCosineScore, bestCosineNode, parent,
        current, state, result, seedIndex, gapIndex, seedK, seedS, seedT, open,
        seedL, T, readSeedCounts, currentGenomeSeedCounts, totalReadSeedCount,
        numGenomes, oldBlockExists, oldBlockStrand, true_node_id, species,
        read_count, mutation_count, firstCall, stopReached, startNode, stopNode,
        groupNodes, searchCount);
  case 16:
    return placementTraversalSpecialized<16>(
        traverser, hitsInThisGenome, maxHitsInAnyGenome, maxHitsNode,
        jaccardNumeratorThisGenome, jaccardDenominatorThisGenome,
        bestJaccardScore, bestJaccardNode, cosineNumeratorThisGenome,
        cosineSumOfSquaresThisGenome, bestCosineScore, bestCosineNode, parent,
        current, state, result, seedIndex, gapIndex, seedK, seedS, seedT, open,
        seedL, T, readSeedCounts, currentGenomeSeedCounts, totalReadSeedCount,
        numGenomes, oldBlockExists, oldBlockStrand, true_node_id, species,
        read_count, mutation_count, firstCall, stopReached, startNode, stopNode,
        groupNodes, searchCount);
  case 32:
    return placementTraversalSpecialized<32>(
        traverser, hitsInThisGenome, maxHitsInAnyGenome, maxHitsNode,
        jaccardNumeratorThisGenome, jaccardDenominatorThisGenome,
        bestJaccardScore, bestJaccardNode, cosineNumeratorThisGenome,
        cosineSumOfSquaresThisGenome, bestCosineScore, bestCosineNode, parent,
        current, state, result, seedIndex, gapIndex, seedK, seedS, seedT, open,
        seedL, T, readSeedCounts, currentGenomeSeedCounts, totalReadSeedCount,
        numGenomes, oldBlockExists, oldBlockStrand, true_node_id, species,
        read_count, mutation_count, firstCall, stopReached, startNode, stopNode,
        groupNodes, searchCount);
  }
  // If no optimized implementation is available, log a message and suggest an
  // optimal size
  if (firstCall) {
    int32_t suggestedK = fixed_kmer::suggestOptimalKmerSize(seedK);
    msg("Using default implementation for k={} (not optimized). Consider using "
        "k={} instead.",
        seedK, suggestedK);
  }

  // The original implementation follows
  TIME_FUNCTION;

  // Cache manager reference at start
  auto &manager = traverser.getCoordManager();
  auto &sequence = manager.getSequence();
  auto &blockExists = manager.getBlockExists();
  auto &blockStrand = manager.getBlockStrand();

  // Check if current node is in the group
  if (groupNodes.find(current) == groupNodes.end()) {
    msg("[DFS] Node {} not in group, stopping traversal branch",
        current->identifier);
    stopReached = true;
    return;
  }

  int64_t hitsBefore = hitsInThisGenome;
  int64_t dfsIndex = state.dfsIndexes[current->identifier];

  // Initialize block ranges if needed
  size_t ranges_initialized = 0;
  for (size_t blockId = 0; blockId < sequence.size(); blockId++) {
    auto range = manager.getBlockRange(blockId);
    if (range.first == 0 && range.second == 0 && blockExists[blockId].first) {
      ranges_initialized++;

      // Scan coordinates to find actual range
      int64_t start = -1, end = -1;
      for (int64_t i = 0; i < manager.size(); i++) {
        const auto *coord = manager.getTupleCoord(i);
        if (coord && coord->blockId == blockId) {
          if (start == -1)
            start = i;
          end = i;
        }
      }
      if (start != -1 && end != -1) {
        manager.updateBlockRange(blockId, start, end);
      }
    }
  }

  PlacementNodeState nodeState(traverser.getCoordManager().getNumBlocks(),
                               traverser.getCoordManager().getNumCoords());

  // Initialize mutation state
  nodeState.oldBlockExists = manager.getBlockExists();
  nodeState.oldBlockStrand = manager.getBlockStrand();
  nodeState.gapRunBacktracks.clear();

  // Apply mutations
  seed_annotated_tree::applyMutations(
      nodeState.blockMutationInfo, nodeState.recompRanges,
      nodeState.nucleotideMutationInfo, nodeState.gapRunUpdates,
      nodeState.gapRunBacktracks, nodeState.gapMapUpdates, T, current,
      traverser, true, state.inverseBlockIds,
      nodeState.inverseBlockIdsBacktrack, seedK);

  // Update gap map with current mutations
  manager.updateGapMap(nodeState.gapRunUpdates, nodeState.gapRunBacktracks,
                       nodeState.gapMapUpdates);

  // Process gap mutations
  auto gapMutationsList = gapIndex[dfsIndex].getDeltas();
  for (const auto &gapMutation : gapMutationsList) {
    int32_t pos = gapMutation.getPos();
    auto maybeValue = gapMutation.getMaybeValue();
    if (maybeValue.isValue()) {
      auto it = manager.getGapMap().find(pos);
      if (it != manager.getGapMap().end()) {
        nodeState.gapRunBacktracks.emplace_back(
            false, std::make_pair(pos, it->second));
      } else {
        nodeState.gapRunBacktracks.emplace_back(
            true, std::make_pair(pos, maybeValue.getValue()));
      }
      manager.setGapMap(pos, maybeValue.getValue());
    } else {
      auto it = manager.getGapMap().find(pos);
      nodeState.gapRunBacktracks.emplace_back(false,
                                              std::make_pair(pos, it->second));
      manager.removeFromGapMap(pos);
    }
  }

  // <seedPos, seedHash, isReverse, endPos, seedCount>
  std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>> seedChanges;
  seedChanges.reserve(100); // Reasonable estimate of changes per node

  // Process seed mutations directly instead of storing in localSeedChanges
  auto currBasePositions = seedIndex[dfsIndex].getBasePositions();
  auto currPerPosMasks = seedIndex[dfsIndex].getPerPosMasks();

  // Using fixed-size buffer for better performance
  alignas(64) char seedBuffer[seedK];

  for (int i = 0; i < currBasePositions.size(); ++i) {
    int64_t pos = currBasePositions[i];
    uint64_t tritMask = currPerPosMasks[i];

    for (int j = 0; j < 32; ++j) {
      uint8_t ternaryNumber = (tritMask >> (j * 2)) & 0x3;
      int64_t seedPos = pos - j;

      if (ternaryNumber == 1) { // on -> off
        if (!state.onSeedsHash[seedPos].has_value()) {
          continue;
        }

        // Get current values before changing
        auto [oldSeed, oldEndPos, oldIsReverse] =
            state.onSeedsHash[seedPos].value();

        // Track seed count changes directly
        size_t oldSeedVal = oldSeed;
        size_t oldSeedCount = currentGenomeSeedCounts[oldSeedVal];

        seedChanges.emplace_back(seedPos, oldSeedVal, oldIsReverse, oldEndPos,
                                 oldSeedCount);

        if (currentGenomeSeedCounts.find(oldSeedVal) ==
            currentGenomeSeedCounts.end()) {
          currentGenomeSeedCounts[oldSeedVal] = 0;
        }
        if (currentGenomeSeedCounts[oldSeedVal] > 0) {
          currentGenomeSeedCounts[oldSeedVal]--;
        }

        // Update metrics directly
        auto readIt = readSeedCounts.find(oldSeedVal);
        if (readIt != readSeedCounts.end() && oldSeedCount > 0) {
          // Decrement hits since we're removing a seed
          hitsInThisGenome--;
        }

        // Apply state change directly
        state.onSeedsHash[seedPos].reset();
        int blockId = manager.getBlockIdOfScalarCoord(seedPos);
        state.BlocksToSeeds[blockId].erase(seedPos);

      } else if (ternaryNumber == 2) { // off -> on or on -> on (change)
        size_t newSeedHash;            // set by getSeedAt
        int64_t newEndPos;             // set by getSeedAt
        bool newIsReverse;             // set by getSeedAt

        if (!seed_annotated_tree::getSeedAt(traverser, manager, newSeedHash,
                                            newIsReverse, newEndPos, seedPos, T,
                                            seedK)) {
          continue;
        }

        if (!state.onSeedsHash[seedPos].has_value()) {
          seedChanges.emplace_back(seedPos, newSeedHash, newIsReverse,
                                   newEndPos, 0);

          // off -> on case
          // Track seed count changes directly
          size_t newSeedVal = newSeedHash;
          if (currentGenomeSeedCounts.find(newSeedVal) ==
              currentGenomeSeedCounts.end()) {
            currentGenomeSeedCounts[newSeedVal] = 0;
          }
          currentGenomeSeedCounts[newSeedVal]++;

          // Update metrics directly
          size_t newSeedCount = currentGenomeSeedCounts[newSeedVal];
          auto readIt = readSeedCounts.find(newSeedVal);
          if (readIt != readSeedCounts.end() && newSeedCount > 0) {
            // Increment hits since we're adding a seed
            hitsInThisGenome++;
          }

          // Apply state change directly
          state.onSeedsHash[seedPos] = {newSeedHash, newEndPos, newIsReverse};
          int blockId =
              traverser.getCoordManager().getBlockIdOfScalarCoord(seedPos);
          state.BlocksToSeeds[blockId].insert(seedPos);
        } else {
          // on -> on (change) case
          auto [oldSeed, oldEndPos, oldIsReverse] =
              state.onSeedsHash[seedPos].value();

          // Skip if no actual change
          if (oldSeed == newSeedHash && oldIsReverse == newIsReverse &&
              oldEndPos == newEndPos) {
            continue;
          }

          // Track seed count changes directly
          size_t oldSeedVal = oldSeed;
          size_t newSeedVal = newSeedHash;
          size_t oldSeedCount = currentGenomeSeedCounts[oldSeedVal];

          seedChanges.emplace_back(seedPos, oldSeedVal, oldIsReverse, oldEndPos,
                                   oldSeedCount);

          if (currentGenomeSeedCounts.find(oldSeedVal) ==
              currentGenomeSeedCounts.end()) {
            currentGenomeSeedCounts[oldSeedVal] = 0;
          }
          if (currentGenomeSeedCounts[oldSeedVal] > 0) {
            currentGenomeSeedCounts[oldSeedVal]--;
          }

          if (currentGenomeSeedCounts.find(newSeedVal) ==
              currentGenomeSeedCounts.end()) {
            currentGenomeSeedCounts[newSeedVal] = 0;
          }
          currentGenomeSeedCounts[newSeedVal]++;

          // Update metrics directly - for replacement, remove old and add new
          auto readItOld = readSeedCounts.find(oldSeedVal);
          if (readItOld != readSeedCounts.end() && oldSeedCount == 0) {
            // Decrement hits for the old seed value
            hitsInThisGenome--;
          }

          auto readItNew = readSeedCounts.find(newSeedVal);
          if (readItNew != readSeedCounts.end() &&
              currentGenomeSeedCounts[newSeedVal] == 1) {
            // Increment hits for the new seed value
            hitsInThisGenome++;
          }

          // Apply state change directly
          state.onSeedsHash[seedPos] = {newSeedHash, newEndPos, newIsReverse};
        }
      }
    }
  }

  // Check if this node has a better score
  if (hitsInThisGenome > maxHitsInAnyGenome) {
    maxHitsInAnyGenome = hitsInThisGenome;
    maxHitsNode = current;

    // Update progress state with new max node
    if (progress_state) {
      std::lock_guard<std::mutex> lock(progress_state->mtx);
      progress_state->maxNodeId = current->identifier;
      progress_state->maxHits = maxHitsInAnyGenome;
      progress_state->tiedNodes.clear();
      progress_state->tiedNodes.push_back(current->identifier);
    }
  }
  // Track nodes that have the same score (tied)
  else if (hitsInThisGenome == maxHitsInAnyGenome && maxHitsInAnyGenome > 0) {
    // Update progress state with tied node
    if (progress_state) {
      std::lock_guard<std::mutex> lock(progress_state->mtx);
      // Check if this node is already in the tied nodes list
      if (std::find(progress_state->tiedNodes.begin(),
                    progress_state->tiedNodes.end(),
                    current->identifier) == progress_state->tiedNodes.end()) {
        progress_state->tiedNodes.push_back(current->identifier);
      }
    }
  }

  nodeState.recompRanges.clear();
  nodeState.gapRunUpdates.clear();
  nodeState.gapMapUpdates.clear();

  // Recursively process children
  for (Node *child : current->children) {
    if (stopReached) {
      err("BREAKING PLACEMENT TRAVERSAL");
      break;
    }
    if (groupNodes.find(child) == groupNodes.end()) {
      err("BREAKING PLACEMENT TRAVERSAL");
      stopReached = true;
      break;
    }

    // Process child
    placementTraversal(
        traverser, hitsInThisGenome, maxHitsInAnyGenome, maxHitsNode,
        jaccardNumeratorThisGenome, jaccardDenominatorThisGenome,
        bestJaccardScore, bestJaccardNode, cosineNumeratorThisGenome,
        cosineSumOfSquaresThisGenome, bestCosineScore, bestCosineNode, current,
        child, state, result, seedIndex, gapIndex, seedK, seedS, seedT, open,
        seedL, T, readSeedCounts, currentGenomeSeedCounts, totalReadSeedCount,
        numGenomes, oldBlockExists, oldBlockStrand, true_node_id, species,
        read_count, mutation_count, false, stopReached, startNode, stopNode,
        groupNodes, searchCount);
  }

  backtrackNode(state, traverser, T, current, nodeState,
                currentGenomeSeedCounts, seedChanges,
                nodeState.gapRunBacktracks);

  // Update progress state with current node info
  if (progress_state) {
    std::lock_guard<std::mutex> lock(progress_state->mtx);
    progress_state->currentNodeId = current->identifier;
    progress_state->nodesVisited++;
  }
}

void backtrackNode(
    PlacementGlobalState &state, CoordinateTraverser &traverser, Tree *T,
    Node *current, PlacementNodeState &nodeState,
    std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
    std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>>
        &seedChanges,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>
        &gapRunBacktracks) {

  // Undo mutations
  undoMutations(T, current, nodeState.blockMutationInfo,
                nodeState.nucleotideMutationInfo, traverser, gapRunBacktracks);

  // Now restore the seed state for all changed positions
  for (const auto &[pos, oldSeedVal, oldIsReverse, oldEndPos, oldSeedCount] :
       seedChanges) {
    state.onSeedsHash[pos] = {oldSeedVal, oldEndPos, oldIsReverse};
    currentGenomeSeedCounts[oldSeedVal] = oldSeedCount;
  }
}

void place(int64_t &maxHitsInAnyGenome, Node *&maxHitsNode,
           double &bestJaccardScore, Node *&bestJaccardNode,
           double &bestCosineScore, Node *&bestCosineNode,
           PlacementResult &result, Tree *T, Index::Reader &index,
           const std::string &reads1Path, const std::string &reads2Path,
           seed_annotated_tree::mutationMatrices &mutMat, std::string prefix,
           std::string refFileName, std::string samFileName,
           std::string bamFileName, std::string mpileupFileName,
           std::string vcfFileName, std::string aligner,
           const std::string &refNode, const bool &save_jaccard,
           const bool &show_time, const float &score_proportion,
           const int &max_tied_nodes, const std::string &true_node_id,
           const std::string &species, const int &read_count,
           const int &mutation_count, const std::string &placementFileName) {
  auto start = std::chrono::high_resolution_clock::now();

  // Replace Timer constructors with static method calls
  timing::Timer::start("read");
  timing::Timer::start("seedIndex");
  timing::Timer::start("traverse");
  timing::Timer::start("write");

  // Initialize progress state
  progress_state = std::make_shared<PlacementProgressState>();

  // Skip UI initialization - use console logs instead
  bool use_ui = false;

  /* Sets up structures used in DFS */
  int32_t k = 0;
  int32_t s = 0;
  int32_t t = 0;
  bool open = false;
  int32_t l = 0;

  coordinates::blockExists_t rootBlockExists;
  coordinates::blockStrand_t rootBlockStrand;

  // Use scoped smart pointers that will auto-cleanup
  {
    auto sequence = std::make_unique<sequence_t>();
    auto blockExists = std::make_unique<blockExists_t>();
    auto blockStrand = std::make_unique<blockStrand_t>();

    // The CoordinateManager uses but doesn't own these resources
    auto coordManager = std::make_unique<CoordinateManager>(
        sequence.get(), blockExists.get(), blockStrand.get());

    // Initialize mutable tree data
    std::vector<std::optional<seeding::onSeedsHash>> onSeedsHash;

    // Initialize data structures
    seed_annotated_tree::setupPlacement(onSeedsHash, *sequence, *blockExists,
                                        *blockStrand, T);

    // Initialize traverser
    coordinates::CoordinateTraverser traverser(
        coordManager->getLeftmostCoordinate(), coordManager.get());

    msg("Starting placement.");

    {
      TIME_BLOCK("setup");
      // PlacementObjects prePlacementObjects(T);

      // rootBlockExists = prePlacementObjects.data.blockExists;
      // rootBlockStrand = prePlacementObjects.data.blockStrand;

      // err("Reading index parameters:");
      k = index.getK();
      // err("  k: {}", k);
      s = index.getS();
      // err("  s: {}", s);
      t = index.getT();
      // err("  t: {}", t);
      open = index.getOpen();
      // err("  open: {}", open);
      l = index.getL();
      // err("  l: {}", l);

      // Check if k is an optimal size, and suggest a better one
      int32_t suggestedK = fixed_kmer::suggestOptimalKmerSize(k);
      if (k != suggestedK) {
        msg("Notice: k={} is not optimized for hardware efficiency", k);
        msg("For better performance, using k={} instead", suggestedK);

        // Automatically use optimized k-mer size
        k = suggestedK;
        msg("Using optimized k-mer size: k={}", k);
      } else {
        msg("Using optimized k-mer size: k={}", k);
      }

      if (k <= 0 || s <= 0) {
        throw std::runtime_error(
            fmt::format("Invalid seed parameters: k={}, s={}, t={}, open={}, "
                        "l={}. These must be positive values.",
                        k, s, t, open, l));
      }
    }

    /* Process read seeds */
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::vector<std::vector<seed>> readSeeds;
    std::unordered_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
    size_t totalReadSeedCount = 0;
    size_t numGenomes = T->allNodes.size();

    {
      TIME_BLOCK("seedsFromFastq");

      seedsFromFastq(k, s, t, open, l, readSeedCounts, readSequences, readQuals,
                     readNames, readSeeds, reads1Path, reads2Path);

      msg("Found {} seeds in {} reads", readSeedCounts.size(),
          readSequences.size());
    }

    int numThreads = 1;

    bool pairedEndReads = reads2Path.size();

    // Step 1: Compute the DFS order
    std::vector<Node *> dfsOrder;
    {
      TIME_BLOCK("computeDFSOrder");
      dfsOrder = threading::computeDFSOrder(T->root);
    }

    // Step 2: Split into groups
    std::vector<threading::GroupInfo> groups;
    {
      TIME_BLOCK("splitDFSIntoGroups");
      groups = threading::splitDFSIntoGroups(dfsOrder, numThreads);
    }

    // err("Split nodes into {} groups", groups.size());

    // Initialize shared state
    tbb::concurrent_vector<ScoreGlobal> scoreGlobalsJaccard;
    tbb::concurrent_vector<ScoreGlobal> scoreGlobalsCosine;
    ::capnp::List<GapMutations>::Reader perNodeGapMutations_Reader =
        index.getPerNodeGapMutations();
    ::capnp::List<SeedMutations>::Reader perNodeSeedMutations_Reader =
        index.getPerNodeSeedMutations();

    seed_annotated_tree::PlacementGlobalState state(
        T, traverser.getCoordManager().getNumBlocks(),
        traverser.getCoordManager().getNumCoords());

    // Process each group
    for (const auto &group : groups) {

      // err("Processing group starting at node {}",
      // group.startNode->identifier); Create new objects for this group

      // err("Created group placement objects");

      rootBlockExists = traverser.getCoordManager().getBlockExists();
      rootBlockStrand = traverser.getCoordManager().getBlockStrand();

      // Initialize group state
      int64_t hitsInThisGenome = 0;
      int64_t jaccardNumeratorThisGenome = 0;
      double cosineNumeratorThisGenome = 0;
      double cosineSumOfSquaresThisGenome = 0;

      std::unordered_map<size_t, int64_t> currentGenomeSeedCounts;
      bool stopReached = false;
      int64_t searchCount = 0;

      int64_t jaccardDenominatorThisGenome = [&]() {
        int64_t sum = 0;
        for (const auto &[seed, counts] : readSeedCounts) {
          sum += counts.first + counts.second;
        }
        return sum;
      }();

      // Initialize node state and buffers

      // err("Starting placement DFS for group");
      placementTraversal(
          traverser, hitsInThisGenome, maxHitsInAnyGenome, maxHitsNode,
          jaccardNumeratorThisGenome, jaccardDenominatorThisGenome,
          bestJaccardScore, bestJaccardNode, cosineNumeratorThisGenome,
          cosineSumOfSquaresThisGenome, bestCosineScore, bestCosineNode,
          group.startNode, group.startNode, state, result,
          perNodeSeedMutations_Reader, perNodeGapMutations_Reader, k, s, t,
          open, l, T, readSeedCounts, currentGenomeSeedCounts,
          totalReadSeedCount, numGenomes, rootBlockExists, rootBlockStrand,
          true_node_id, species, read_count, mutation_count, true, stopReached,
          group.startNode, group.stopNode, group.groupNodes, searchCount);

      // err("Completed processing group starting at node {}",
      // group.startNode->identifier);
    }

    // LinkedNode* curr = maxHitsNode;
    // LinkedNode* bestNode = nullptr;
    std::cout << "[Placement result] Most hits: " << maxHitsNode->identifier
              << " (" << maxHitsInAnyGenome << ")" << std::endl;
    double bestNodeMapQ = 0;
    std::string bestNodeSequence;
    std::unordered_map<size_t,
                       std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
        bestSeedToRefPositions;
    auto start_time = std::chrono::high_resolution_clock::now();

    // while (curr != nullptr) {
    //     // std::cout << "[mapq] processing node: " << curr->node->identifier
    //     << std::endl;
    //     // std::cout << "[candidate] " << curr->node->identifier <<
    //     std::endl; struct Chain {
    //         int32_t score;
    //         int32_t query_start;  // min y coordinate
    //         int32_t query_end;    // max y coordinate + span
    //         std::vector<size_t> anchor_ids;
    //         bool is_primary;
    //         std::vector<Chain *> secondary_chains;
    //     };

    //     std::vector<Chain> node_chains;
    //     const float chn_pen_gap = 0.01f;  // Gap penalty coefficient
    //     const float chn_pen_skip = 0.01f;
    //     const int max_dist_x = 5000;      // Max reference distance
    //     const int max_dist_y = 5000;      // Max query distance
    //     const int bw = 500;               // Bandwidth

    //     // Get node sequence and create degap mapping
    //     std::string currNodeSequence = "";
    //     std::string gappedSeq =
    //     T->getStringFromReference(curr->node->identifier, true, true);
    //     std::vector<int32_t> degap;
    //     for (int32_t i = 0; i < gappedSeq.size(); i++) {
    //         char &c = gappedSeq[i];
    //         degap.push_back(currNodeSequence.size());
    //         if (c != '-') {
    //             currNodeSequence += c;
    //         }
    //     }

    //     // Map seeds to reference positions
    //     std::unordered_map<size_t, std::pair<std::vector<uint32_t>,
    //     std::vector<uint32_t>>> currSeedToRefPositions; for(int i = 0; i <
    //     curr->seeds.size(); i++) {
    //         if(curr->seeds[i].has_value()) {
    //             size_t seed = curr->seeds[i].value().hash;
    //             bool reversed = curr->seeds[i].value().isReverse;
    //             int pos = degap[i];

    //             if (currSeedToRefPositions.find(seed) ==
    //             currSeedToRefPositions.end()) {
    //                 std::vector<uint32_t> a;
    //                 std::vector<uint32_t> b;
    //                 currSeedToRefPositions[seed] = std::make_pair(a,b);
    //             }

    //             if(reversed) {
    //                 currSeedToRefPositions[seed].second.push_back(pos);
    //             } else {
    //                 currSeedToRefPositions[seed].first.push_back(pos);
    //             }
    //         }
    //     }

    //     // std::cout << "[getAnchorsmapq] got " <<
    //     currSeedToRefPositions.size() << " seeds for node: " <<
    //     curr->node->identifier << std::endl;

    //     // Get anchors from seeds
    //     std::vector<std::tuple<int64_t, int32_t, int>> anchors;
    //     alignment::getAnchors(const_cast<std::vector<std::tuple<int64_t,
    //     int32_t, int>> &>(anchors),
    //     const_cast<std::vector<std::vector<seeding::seed>> &>(readSeeds),
    //     const_cast<std::vector<std::string> &>(readSequences),
    //     const_cast<std::unordered_map<size_t,
    //     std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>
    //     &>(currSeedToRefPositions), k);
    //     // std::cout << "[mapq] got " << anchors.size() << " anchors for
    //     node: " << curr->node->identifier << std::endl;

    //     if (anchors.empty()) {
    //         // std::cout << "[DEBUG] No anchors found for node: " <<
    //         curr->node->identifier << std::endl; continue;
    //     }

    //     // Sort anchors by reference position
    //     std::sort(anchors.begin(), anchors.end(), [](const auto& a, const
    //     auto& b) {
    //         return std::get<0>(a) < std::get<0>(b);
    //     });

    //     // Initialize DP arrays with negative values to ensure proper chain
    //     formation std::vector<int32_t> f(anchors.size(),
    //     std::numeric_limits<int32_t>::min());  // Chaining scores
    //     std::vector<int64_t> p(anchors.size(), -1);  // Backtrack pointers
    //     std::vector<int32_t> v(anchors.size(),
    //     std::numeric_limits<int32_t>::min());  // Peak scores

    //     // Fill DP arrays using minimap2's algorithm
    //     for (size_t i = 0; i < anchors.size(); i++) {
    //         const auto& [xi, yi, wi] = anchors[i];
    //         f[i] = wi;  // Initialize with current anchor span
    //         v[i] = wi;  // Initialize peak score with current span

    //         // std::cout << "[DEBUG] Processing anchor " << i << " - pos: ("
    //         << xi << "," << yi << "), span: " << wi << std::endl;

    //         // Find best predecessor
    //         int64_t max_j = -1;
    //         int64_t st = std::max(0L, static_cast<int64_t>(i) - 50);  //
    //         Limit search to last 50 anchors like minimap2

    //         for (int64_t j = i-1; j >= st; j--) {
    //             const auto& [xj, yj, wj] = anchors[j];

    //             // Compute distances like minimap2
    //             int32_t dq = yi - yj;  // Query distance
    //             int32_t dr = xi - xj;  // Reference distance

    //             // Skip if distances too large
    //             if (dq <= 0 || dq > max_dist_x) continue;
    //             if (dr <= 0 || dr > max_dist_y) continue;

    //             // Compute diagonal difference and gap
    //             int32_t dd = dr > dq ? dr - dq : dq - dr;
    //             if (dd > bw) continue;  // Skip if outside bandwidth

    //             // Compute minimum gap and score
    //             int32_t dg = dr < dq? dr : dq;  // Minimum gap
    //             int32_t sc = std::min(wi, dg);  // Score using current
    //             anchor's span

    //             // Apply gap cost when diagonal difference exists or gap
    //             larger than span if (dd || dg > wi) {  // Compare with
    //             current span
    //                 float lin_pen = chn_pen_gap * dd + chn_pen_skip * dg;
    //                 float log_pen = dd >= 1? std::log2(dd + 1) : 0.0f;
    //                 sc -= static_cast<int32_t>(lin_pen + 0.5f * log_pen);
    //             }

    //             // Add previous chain score
    //             sc += f[j];

    //             if (sc > f[i]) {
    //                 f[i] = sc;
    //                 p[i] = j;
    //                 // std::cout << "[DEBUG] Updated chain - anchor: " << i
    //                 << ", prev: " << j << ", score: " << sc << std::endl;
    //             }
    //         }

    //         // Update peak score - ensure we're not using uninitialized
    //         values if (p[i] >= 0) {
    //             v[i] = std::max(f[i], v[p[i]]);
    //         } else {
    //             v[i] = f[i];
    //         }

    //         // std::cout << "[DEBUG] Final scores for anchor " << i << " - f:
    //         " << f[i] << ", v: " << v[i] << std::endl;
    //     }

    //     // After computing f[] and p[], collect chains by backtracking
    //     std::vector<bool> used(anchors.size(), false);
    //     int total_chains_attempted = 0;
    //     int chains_filtered_by_size = 0;
    //     int chains_filtered_by_bounds = 0;
    //     int chains_filtered_by_coords = 0;

    //     for(size_t i = 0; i < anchors.size(); i++) {
    //         if(used[i] || f[i] < 0) {
    //             // if (f[i] < 0) std::cout << "[DEBUG] Skipping anchor " << i
    //             << " due to negative score: " << f[i] << std::endl; continue;
    //         }
    //         total_chains_attempted++;

    //         Chain chain;
    //         chain.score = f[i];
    //         chain.query_start = std::numeric_limits<int32_t>::max();
    //         chain.query_end = std::numeric_limits<int32_t>::min();

    //         // Backtrack to collect anchors
    //         size_t curr_anchor = i;
    //         int32_t steps = 0;
    //         const int32_t max_steps = anchors.size();
    //         // std::cout << "[DEBUG] Before while loop - curr_anchor: " <<
    //         curr_anchor
    //         //           << ", used[curr_anchor]: " << (used[curr_anchor] ?
    //         "true" : "false")
    //         //           << ", steps: " << steps
    //         //           << ", max_steps: " << max_steps << std::endl;

    //         while(curr_anchor != (size_t)-1 && steps++ < max_steps) {
    //             // std::cout << "[DEBUG] Inside while loop - iteration: " <<
    //             steps
    //             //           << ", curr_anchor: " << curr_anchor
    //             //           << ", p[curr_anchor]: " << p[curr_anchor] <<
    //             std::endl;

    //             const auto& [x, y, w] = anchors[curr_anchor];
    //             // std::cout << "[DEBUG] Anchor details - x: " << x << ", y:
    //             " << y << ", w: " << w << std::endl;

    //             if (y > std::numeric_limits<int32_t>::max() - w) {
    //                 // std::cout << "[WARN:OVERFLOW] y: " << y << " w: " << w
    //                 << " for node: " << curr->node->identifier << std::endl;
    //                 continue;
    //             }
    //         chain.anchor_ids.push_back(curr_anchor);

    //         // Track query coordinates
    //         chain.query_end = std::max(chain.query_end, y + w);
    //         chain.query_start = std::min(chain.query_start, y);
    //         int64_t next_anchor = p[curr_anchor];
    //         if (next_anchor == -1) {
    //           break;
    //         }
    //         curr_anchor = static_cast<size_t>(next_anchor);
    //         // std::cout << "[DEBUG] Updated curr_anchor to: " << curr_anchor
    //         << std::endl;
    //     }

    //     // std::cout << "[DEBUG] After while loop - final curr_anchor: " <<
    //     curr_anchor
    //     //           << ", final steps: " << steps << std::endl;

    //     // std::cout << "[DEBUG] Chain details - size: " <<
    //     chain.anchor_ids.size()
    //     //           << ", start: " << chain.query_start
    //     //           << ", end: " << chain.query_end
    //     //           << ", score: " << chain.score << std::endl;

    //     if(chain.anchor_ids.size() < 3) {
    //         chains_filtered_by_size++;
    //         continue;
    //     }

    //     if(chain.query_start > chain.query_end) {
    //         chains_filtered_by_bounds++;
    //         continue;
    //     }

    //     if(chain.query_start == std::numeric_limits<int32_t>::max() ||
    //        chain.query_end == std::numeric_limits<int32_t>::min()) {
    //         chains_filtered_by_coords++;
    //         continue;
    //     }
    //     if (chain.query_start <= chain.query_end &&
    //         chain.query_start != std::numeric_limits<int32_t>::max() &&
    //         chain.query_end != std::numeric_limits<int32_t>::min()) {
    //           // std::cout << "[mapq] adding chain with " <<
    //           chain.anchor_ids.size() << " anchors for node: " <<
    //           curr->node->identifier << std::endl; for (size_t anchor_id :
    //           chain.anchor_ids) {
    //             used[anchor_id] = true;
    //           }
    //           node_chains.push_back(chain);
    //     }
    // }

    //     // std::cout << "[DEBUG] Chain filtering stats for node " <<
    //     curr->node->identifier << ":" << std::endl
    //     //           << "  Total chains attempted: " <<
    //     total_chains_attempted << std::endl
    //     //           << "  Filtered by size (<3): " <<
    //     chains_filtered_by_size << std::endl
    //     //           << "  Filtered by bounds (start>end): " <<
    //     chains_filtered_by_bounds << std::endl
    //     //           << "  Filtered by invalid coords: " <<
    //     chains_filtered_by_coords << std::endl
    //     //           << "  Chains remaining: " << node_chains.size() <<
    //     std::endl;

    //     // Sort all chains by score
    //     std::sort(node_chains.begin(), node_chains.end(),
    //         [](const Chain& a, const Chain& b) { return a.score > b.score;
    //         });

    //     // Find primary chains
    //     std::vector<Chain*> primary_chains;
    //     for(auto& chain : node_chains) {
    //         bool is_secondary = false;
    //         for(auto* pc : primary_chains) {
    //             // Check overlap on query coordinates
    //             int32_t overlap_start = std::max(chain.query_start,
    //             pc->query_start); int32_t overlap_end =
    //             std::min(chain.query_end, pc->query_end); if(overlap_end >
    //             overlap_start) {
    //                 int32_t overlap_len = overlap_end - overlap_start;
    //                 int32_t chain_len = chain.query_end - chain.query_start;
    //                 int32_t pc_len = pc->query_end - pc->query_start;
    //                 float overlap_frac = static_cast<float>(overlap_len) /
    //                 std::min(chain_len, pc_len);

    //                 // std::cout << "[DEBUG] Chain overlaps with primary
    //                 chain - "
    //                 //           << "overlap: " << overlap_len
    //                 //           << "/" << std::min(chain_len, pc_len)
    //                 //           << " (" << (overlap_frac * 100) << "%)"
    //                 //           << ", primary span: [" << pc->query_start <<
    //                 "," << pc->query_end << "]"
    //                 //           << ", primary score: " << pc->score <<
    //                 std::endl;

    //                 if(overlap_frac >= 0.5f && chain.score < pc->score) {
    //                     is_secondary = true;
    //                     // std::cout << "[DEBUG] Marked as secondary due to "
    //                     //           << (overlap_frac * 100) << "% overlap
    //                     with higher scoring primary chain" << std::endl;
    //                     Chain *secondary = new Chain(chain);
    //                     pc->secondary_chains.push_back(secondary);
    //                     break;
    //                 }
    //         }
    //         if(!is_secondary) {
    //             chain.is_primary = true;
    //             Chain *primary = new Chain(chain);
    //             primary_chains.push_back(primary);
    //             // std::cout << "[DEBUG] Added to primary chains" <<
    //             std::endl;
    //         }
    //     }
    //     // std::cout << "[mapq] got " << primary_chains.size() << " primary
    //     chains for node: " << curr->node->identifier << std::endl;
    //     // Process all primary chains and compute mapQ scores
    //     if(!primary_chains.empty()) {
    //         double totalMapQ = 0;

    //         double numPrimary = static_cast<double>(primary_chains.size());
    //         double numTotal = static_cast<double>(node_chains.size());

    //         for(const auto* primary_chain : primary_chains) {
    //             // std::cout << "[mapq] processing primary chain: " <<
    //             primary_chain->score << " for node: " <<
    //             curr->node->identifier << std::endl;
    //             // Find best secondary chain that overlaps with this primary
    //             chain double secondaryScore = 0.0; for(const auto&
    //             secondary_chain : primary_chain->secondary_chains) {
    //               if (secondary_chain->score > secondaryScore) {
    //                 secondaryScore = secondary_chain->score;
    //               }
    //             }
    //             double m = primary_chain->anchor_ids.size();
    //             double f1 = primary_chain->score;
    //             double f2 = secondaryScore;

    //             // Compute mapQ for this primary chain
    //             double mapQ = (1.0 - f2/f1) * log(f1 + m + numPrimary +
    //             numTotal);
    //             // std::cout << "m=" << m << " f1=" << f1 << " f2=" << f2 <<
    //             " mapQ=" << mapQ << " numPrimary=" << numPrimary << "
    //             numTotal=" << numTotal << std::endl; totalMapQ += mapQ;
    //         }
    //         // std::cout << "totalMapQ=" << totalMapQ << std::endl;
    //         // std::cout << "totalMapQ(withlogfactor)=" << totalMapQ *
    //         log2(numPrimary) << std::endl; for (auto* chain : primary_chains)
    //         {
    //           for (auto* secondary : chain->secondary_chains) {
    //             delete secondary;
    //           }
    //           chain->secondary_chains.clear();
    //           delete chain;
    //         }
    //         primary_chains.clear();
    //         // std::cout << "thisNodeMapQ: " << totalMapQ << " vs.
    //         bestNodeMapQ: " << bestNodeMapQ << std::endl; if (totalMapQ >
    //         bestNodeMapQ) {
    //             // std::cout << "[DEBUG] New best MAPQ: " << totalMapQ << "
    //             for node: " << (curr->node != nullptr ?
    //             curr->node->identifier : "nullptr") << std::endl;
    //             bestNodeMapQ = totalMapQ;
    //             bestNode = curr;
    //             bestNodeSequence = currNodeSequence;
    //             bestSeedToRefPositions = currSeedToRefPositions;
    //         }
    //     }

    //     curr = curr->next;
    // }

    // std::cout << "[Placement result] After chain scoring: " <<
    // bestNode->node->identifier << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time);
    // std::cout << "[TIME] candidate processing took: " << duration.count() <<
    // "ms" << std::endl;

    // std::cout << "[DEBUG] Best MAPQ: " << bestNodeMapQ << " for node: " <<
    // bestNode->node->identifier << std::endl;
    std::string bestMatchSequence;
    if (!refNode.empty()) {
      bestMatchSequence = T->getStringFromReference(refNode, false, true);
    } else {
      bestMatchSequence = bestNodeSequence;
    }

    // Write placement information to file if a filename was provided
    if (!placementFileName.empty()) {
      try {
        std::ofstream placementFile(placementFileName);
        if (placementFile.is_open()) {
          placementFile << "PlacementNode\tHits\tTiedNodes\n";

          // Write the primary placement node
          if (maxHitsNode) {
            placementFile << maxHitsNode->identifier << "\t"
                          << maxHitsInAnyGenome << "\t";

            // Count and write any tied nodes
            std::vector<std::string> tiedNodes;

            // Write tied nodes as a comma-separated list
            if (!tiedNodes.empty()) {
              for (size_t i = 0; i < tiedNodes.size(); ++i) {
                placementFile << tiedNodes[i];
                if (i < tiedNodes.size() - 1) {
                  placementFile << ",";
                }
              }
            } else {
              placementFile << "none";
            }
            placementFile << "\n";
          } else {
            placementFile << "No placement found\t0\tnone\n";
          }

          placementFile.close();
          msg("Placement results written to {}", placementFileName);
        } else {
          err("Failed to open placement file for writing: {}",
              placementFileName);
        }
      } catch (const std::exception &e) {
        err("Error writing placement file: {}", e.what());
      }
    }

    if (show_time) {
      timing::Timer::report();
    }

    // After the traversal is complete, gather tied nodes
    std::vector<std::string> tiedNodes;
    if (maxHitsNode != nullptr) {
      // Find nodes with the same max score
      for (auto &nodePair : result.placementScoresJaccard) {
        if (nodePair.second == maxHitsInAnyGenome) {
          tiedNodes.push_back(nodePair.first->identifier);
        }
      }
    }

    // Update UI with final tied nodes
    if (progress_state) {
      std::lock_guard<std::mutex> lock(progress_state->mtx);
      progress_state->tiedNodes = tiedNodes;
    }

    // Clean up the UI thread
    if (use_ui && progress_state) {
      progress_state->running = false;
    }

    // Final message with placement results
    if (maxHitsNode) {
      std::vector<std::string> tiedNodeIds;
      for (const auto &nodeId : progress_state->tiedNodes) {
        tiedNodeIds.push_back(nodeId);
      }

      std::string tiedNodesStr =
          tiedNodeIds.size() > 1 ? fmt::format(" (tied with {} other nodes)",
                                               tiedNodeIds.size() - 1)
                                 : "";

      msg("Final placement: {} with {} hits{}", maxHitsNode->identifier,
          maxHitsInAnyGenome, tiedNodesStr);

      if (!tiedNodeIds.empty() && tiedNodeIds.size() > 1) {
        msg("Tied nodes: {}", fmt::join(tiedNodeIds, ", "));
      }
    } else {
      msg("No placement found");
    }

    auto placement_duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
    msg("Placement completed in {}ms", placement_duration.count());
  }

  // After the placement is complete, reset progress state
  progress_state.reset();
}

void placeBatch(Tree *T, Index::Reader &index, const std::string &batchFilePath,
                seed_annotated_tree::mutationMatrices &mutMat,
                std::string prefixBase, std::string refFileNameBase,
                std::string samFileNameBase, std::string bamFileNameBase,
                std::string mpileupFileNameBase, std::string vcfFileNameBase,
                std::string aligner, const std::string &refNode,
                const bool &save_jaccard, const bool &show_time,
                const float &score_proportion, const int &max_tied_nodes) {

  // Count total samples in batch file first
  std::ifstream countFile(batchFilePath);
  size_t totalSamples = 0;
  std::string line;
  while (std::getline(countFile, line)) {
    if (!line.empty())
      totalSamples++;
  }
  countFile.close();

  msg("Starting batch processing of {} samples", totalSamples);

  // Read k-mer size from index
  int32_t k = index.getK();

  // Check if k is an optimal size, and suggest a better one
  int32_t suggestedK = fixed_kmer::suggestOptimalKmerSize(k);
  if (k != suggestedK) {
    msg("Notice: k={} is not optimized for hardware efficiency", k);
    msg("For better performance, using k={} instead", suggestedK);

    // Automatically use optimized k-mer size
    k = suggestedK;
    msg("Using optimized k-mer size: k={}", k);
  } else {
    msg("Using optimized k-mer size: k={}", k);
  }

  // Create a master progress state for the batch
  auto batch_progress = std::make_shared<PlacementProgressState>();
  batch_progress->running = true;

  // Process each sample
  std::ifstream batchFile(batchFilePath);
  size_t sampleCount = 0;

  while (std::getline(batchFile, line)) {
    if (line.empty())
      continue;

    sampleCount++;
    std::stringstream ss(line);
    std::string sampleName, reads1, reads2;

    // Parse the line - format is either "sample_name\tread1_path" or
    // "sample_name\tread1_path\tread2_path"
    ss >> sampleName >> reads1;
    ss >> reads2; // This will be empty if there is no second read file

    std::string progressMsg = fmt::format(
        "Processing sample {}/{}: {}", sampleCount, totalSamples, sampleName);
    std::string separator(progressMsg.length(), '=');

    msg("\n{}\n{}\n{}", separator, progressMsg, separator);

    try {
      // Setup output file paths for this sample
      std::string prefix = prefixBase;
      if (!prefix.empty() && prefix.back() != '.')
        prefix += ".";
      prefix += sampleName;

      std::string refFileName =
          refFileNameBase.empty() ? "" : prefix + ".reference.fa";
      std::string samFileName = samFileNameBase.empty() ? "" : prefix + ".sam";
      std::string bamFileName = bamFileNameBase.empty() ? "" : prefix + ".bam";
      std::string mpileupFileName =
          mpileupFileNameBase.empty() ? "" : prefix + ".mpileup";
      std::string vcfFileName = vcfFileNameBase.empty() ? "" : prefix + ".vcf";
      std::string placementFileName = prefix + ".placement.tsv";

      // Initialize for this sample
      int64_t maxHitsInAnyGenome = 0;
      double bestJaccardScore = 0.0;
      double bestCosineScore = 0.0;
      Node *maxHitsNode = nullptr;
      Node *bestJaccardNode = nullptr;
      Node *bestCosineNode = nullptr;
      PlacementResult result;

      // Update batch progress state with sample info
      if (batch_progress) {
        std::lock_guard<std::mutex> lock(batch_progress->mtx);
        batch_progress->currentNodeId = fmt::format(
            "Sample {}/{}: {}", sampleCount, totalSamples, sampleName);
      }

      // Process this sample
      place(maxHitsInAnyGenome, maxHitsNode, bestJaccardScore, bestJaccardNode,
            bestCosineScore, bestCosineNode, result, T, index, reads1, reads2,
            mutMat, prefix, refFileName, samFileName, bamFileName,
            mpileupFileName, vcfFileName, aligner, refNode, save_jaccard,
            show_time, score_proportion, max_tied_nodes, "", "", 0, 0,
            placementFileName);

      // Collect and report tied nodes
      std::vector<std::string> tiedNodes;
      if (maxHitsNode != nullptr) {
        for (auto &nodePair : result.placementScoresJaccard) {
          if (nodePair.second == maxHitsInAnyGenome) {
            tiedNodes.push_back(nodePair.first->identifier);
          }
        }

        std::string tiedNodesStr =
            tiedNodes.size() > 1 ? fmt::format(" (tied with {} other nodes)",
                                               tiedNodes.size() - 1)
                                 : "";

        msg("Sample {} placed at {} with {} hits{}", sampleName,
            maxHitsNode->identifier, maxHitsInAnyGenome, tiedNodesStr);

        if (!tiedNodes.empty() && tiedNodes.size() > 1) {
          msg("Found {} nodes tied for best placement: {}", tiedNodes.size(),
              fmt::join(tiedNodes, ", "));
        }
      } else {
        msg("Sample {} placement failed or produced no results", sampleName);
      }

      // Clean up
      delete maxHitsNode;
      delete bestJaccardNode;
      delete bestCosineNode;

    } catch (const std::exception &e) {
      err("ERROR placing sample {}: {}", sampleName, e.what());
    }
  }

  // Clean up batch progress
  if (batch_progress) {
    batch_progress->running = false;
  }

  msg("Batch processing completed: {} samples processed", sampleCount);
}

} // namespace placement

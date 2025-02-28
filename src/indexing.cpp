#include "indexing.hpp"
#include "coordinates.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "seed_annotated_tree.hpp"
#include <ranges>

using namespace coordinates; // For coordinate types and traversal
using namespace logging;     // For logging utilities

namespace indexing {

bool debug = false;

void indexingTraversal(
    seed_annotated_tree::TraversalGlobalState &state,
    ::capnp::List<SeedMutations>::Builder &perNodeSeedMutations,
    ::capnp::List<GapMutations>::Builder &perNodeGapMutations,
    const PanmapParams &params, panmanUtils::Tree *T, panmanUtils::Node *node,
    coordinates::CoordinateTraverser &traverser) {

  // Track timing for performance monitoring
  auto start_time = std::chrono::high_resolution_clock::now();

  // Cache manager reference at start
  auto &manager = traverser.getCoordManager();
  auto &sequence = manager.getSequence();
  auto &blockExists = manager.getBlockExists();
  auto &blockStrand = manager.getBlockStrand();

  // Create local state for this node
  seed_annotated_tree::NodeLocalState nodeState(manager.getNumBlocks(),
                                                manager.getNumCoords());
  nodeState.oldBlockExists = blockExists;
  nodeState.oldBlockStrand = blockStrand;

  // Get DFS index for current node
  int64_t dfsIndex = state.dfsIndexes[node->identifier];
  msg("Processing node {} of {} (Node ID: {})", dfsIndex, T->allNodes.size(),
      node->identifier);

  try {
    // ===== STEP 1: Apply all block mutations =====
    msg("Applying block mutations for node {}", node->identifier);

    // Apply mutations - this function handles both block and nucleotide
    // mutations
    seed_annotated_tree::applyMutations(
        nodeState.blockMutationInfo, // Records block mutations for backtracking
        nodeState.recompRanges,      // Collects recomputation ranges
        nodeState.nucleotideMutationInfo, // Records nucleotide mutations for
                                          // backtracking
        nodeState.gapRunUpdates,          // Gap run updates to be applied
        nodeState.gapRunBacktracks,       // For backtracking gap changes
        nodeState.gapMapUpdates,          // Tracks gap map updates
        T, node, traverser,
        false,                              // Not placement mode
        state.inverseBlockIds,              // Track inverted blocks
        nodeState.inverseBlockIdsBacktrack, // For backtracking inversions
        params.k); // K-mer length needed for range computation

    // After block mutations, update gap map
    msg("Updating gap map after mutations for node {}", node->identifier);
    manager.updateGapMap(nodeState.gapRunUpdates, nodeState.gapRunBacktracks,
                         nodeState.gapMapUpdates);

    // ===== STEP 2: Process block mutations to update the gap map =====
    // Handle special cases for block mutations (on->off, off->on, inversions)
    std::vector<int64_t> invertedBlocks;
    for (int i = 0; i < blockExists.size(); i++) {
      const bool &oldExists = nodeState.oldBlockExists[i].first;
      const bool &newExists = blockExists[i].first;

      // Track inverted blocks for gap map updates
      if (newExists && !blockStrand[i].first) {
        invertedBlocks.push_back(i);
      }

      // Block turned off - convert to gaps
      if (oldExists && !newExists) {
        const auto &[start, end] = manager.getBlockRange(i);

        // Validate range before updating
        if (start < 0 || end < start || end >= manager.size()) {
          warn("Invalid block range for blockId {} when turning off: [{}:{}]. "
               "Skipping.",
               i, start, end);
          continue;
        }

        manager.updateGapMapStep({true, {start, end}},
                                 nodeState.gapRunBacktracks,
                                 nodeState.gapMapUpdates);
      }
      // Block turned on - convert from gaps
      else if (!oldExists && newExists) {
        // Get block's coordinate range
        const tupleCoord_t *coord = manager.getFirstCoordinateInBlock(i);
        const tupleCoord_t *end = manager.getLastCoordinateInBlock(i);

        if (!coord || !end) {
          warn("Could not get valid coordinates for blockId {} when turning "
               "on. Skipping.",
               i);
          continue;
        }

        // Invert coordinates if block is inverted
        if (!blockStrand[i].first) {
          std::swap(coord, end);
        }

        // Find non-gap positions in this block
        std::pair<int64_t, int64_t> curNucRange = {-1, -1};
        std::vector<std::pair<int64_t, int64_t>> nucRanges;

        // Use block traversal to find nucleotide ranges
        traverser.reset(coord, end);

        // Loop safely
        int safety_counter = 0;
        int max_iterations = 1000000;

        while (!traverser.isDone() && safety_counter < max_iterations) {
          coord = traverser.getCurrent();
          if (!coord) {
            warn("Null coordinate encountered during block {} traversal. "
                 "Breaking.",
                 i);
            break;
          }

          char c = traverser.getChar();
          int64_t scalar = coord->scalar;

          // Track non-gap positions
          if (c != '-' && c != 'x') {
            if (curNucRange.first != -1 && curNucRange.second + 1 == scalar) {
              ++curNucRange.second;
            } else {
              if (curNucRange.first != -1) {
                nucRanges.emplace_back(curNucRange);
              }
              curNucRange = {scalar, scalar};
            }
          }

          traverser.next();
          safety_counter++;

          if (safety_counter >= max_iterations) {
            warn("Reached maximum iterations ({}) traversing blockId {}. "
                 "Breaking.",
                 max_iterations, i);
            break;
          }
        }

        // Add last range if exists
        if (curNucRange.first != -1) {
          nucRanges.emplace_back(curNucRange);
        }

        // Apply discovered ranges to gap map
        if (blockStrand[i].first) {
          // Forward block - use ranges directly
          for (const auto &range : nucRanges) {
            // Validate range
            if (range.first < 0 || range.second < range.first ||
                range.second >= manager.size()) {
              warn("Invalid nucleotide range [{}:{}] for block {}. Skipping.",
                   range.first, range.second, i);
              continue;
            }
            // Remove gaps at nucleotide positions
            manager.updateGapMapStep({false, range}, nodeState.gapRunBacktracks,
                                     nodeState.gapMapUpdates);
          }
        } else {
          // Inverted block - invert ranges first
          std::vector<std::pair<int64_t, int64_t>> invertedRanges =
              manager.invertRanges(nucRanges, manager.getBlockRange(i));

          for (const auto &range : invertedRanges) {
            // Validate range
            if (range.first < 0 || range.second < range.first ||
                range.second >= manager.size()) {
              warn("Invalid inverted range [{}:{}] for block {}. Skipping.",
                   range.first, range.second, i);
              continue;
            }
            // Remove gaps at nucleotide positions
            manager.updateGapMapStep({false, range}, nodeState.gapRunBacktracks,
                                     nodeState.gapMapUpdates);
          }
        }
      }
    }

    // Handle inverted blocks specially
    msg("Processing {} inverted blocks for node {}", invertedBlocks.size(),
        node->identifier);
    for (const auto &blockId : invertedBlocks) {
      auto blockRange = manager.getBlockRange(blockId);

      // Validate range before inversion
      if (blockRange.first < 0 || blockRange.second < blockRange.first ||
          blockRange.second >= manager.size()) {
        warn("Invalid block range for inversion [{}:{}]. Skipping.",
             blockRange.first, blockRange.second);
        continue;
      }

      // Invert gap map for this block range
      manager.invertGapMap(blockRange, nodeState.gapRunBlocksBacktracks,
                           nodeState.gapMapUpdates);
    }

    // Process block backtracking info
    for (auto it = nodeState.gapRunBlocksBacktracks.rbegin();
         it != nodeState.gapRunBlocksBacktracks.rend(); ++it) {
      const auto &[del, range] = *it;
      if (del) {
        manager.removeFromGapMap(range.first);
      } else {
        manager.setGapMap(range.first, range.second);
      }
    }

    // ===== STEP 3: Gather ranges that need recomputation after mutations =====
    msg("Processing recomputation ranges for node {}", node->identifier);

    // Merge ranges before expanding to reduce redundancy
    coordinates::mergeRanges(nodeState.recompRanges, traverser);

    // Expand ranges to include k-mer context
    msg("Expanding ranges to include {} flanking nucleotides", params.k);
    coordinates::expandRanges(nodeState.recompRanges, traverser, params.k);

    // Merge again after expanding
    coordinates::mergeRanges(nodeState.recompRanges, traverser);
    msg("After merging and expanding: {} ranges to process",
        nodeState.recompRanges.size());

    // ===== STEP 4: Process recomp ranges to get seed changes =====
    // First determine global boundary for ranges
    const coordinates::tupleCoord_t *globalEnd =
        manager.getLastCoordinateInOnBlock();
    if (!globalEnd) {
      warn("Could not find last coordinate in any ON block. Using rightmost "
           "coordinate instead.");
      globalEnd = manager.getRightmostCoordinate();
    }

    const int32_t lastOnBlock =
        globalEnd ? globalEnd->blockId : sequence.size();
    const tupleCoord_t *lastOnCoord =
        globalEnd ? globalEnd : manager.getRightmostCoordinate();

    if (!lastOnCoord) {
      throw std::runtime_error(
          "Could not determine last coordinate for processing");
    }

    // Process seed changes for each recomputation range
    msg("Processing seed changes for node {}", node->identifier);

    // Set for collecting seed changes
    std::set<seed_annotated_tree::SeedChange,
             seed_annotated_tree::SeedChangeComparator>
        seedChanges;

    // Process ranges to collect seed changes
    processSeedChanges(state, nodeState, traverser, params, T, node,
                       lastOnBlock, lastOnCoord, seedChanges);

    // ===== STEP 5: Apply seed changes =====
    // Note: processSeedChanges function already applies the changes

    // ===== STEP 6: Store seed changes in index =====
    msg("Storing {} seed changes for node {} in index", seedChanges.size(),
        node->identifier);
    storeSeedChanges(perNodeSeedMutations, perNodeGapMutations, dfsIndex,
                     nodeState, seedChanges);

    // ===== STEP 7: Recursive step - process children =====
    msg("Processing {} children of node {}", node->children.size(),
        node->identifier);
    for (panmanUtils::Node *child : node->children) {
      try {
        // Process child nodes with current state
        indexingTraversal(state, perNodeSeedMutations, perNodeGapMutations,
                          params, T, child, traverser);
      } catch (const std::exception &e) {
        warn("Error processing child node {}: {}. Continuing with next child.",
             child->identifier, e.what());
      }
    }

    // ===== STEP 8: Undo seed changes after all children processed =====
    msg("Undoing seed changes for node {}", node->identifier);
    undoSeedChanges(state, nodeState, traverser, seedChanges);

    // ===== STEP 9: Undo all mutations =====
    msg("Undoing mutations for node {}", node->identifier);

    // CRITICAL FIX: Properly reinitialize all block coordinates BEFORE undoing
    // mutations This is essential to ensure coordinate system consistency
    msg("Reinitializing block coordinates prior to undo");

    // First validate gap map and attempt to fix issues if found
    bool gap_map_valid = manager.validateGapMap(true);
    if (!gap_map_valid) {
      warn("Gap map had issues before undoing mutations - attempted to fix "
           "automatically");

      // Check if optimization succeeded
      if (!manager.validateGapMap(false)) {
        err("Gap map still has issues after optimization. Proceeding with "
            "caution.");
      }
    }

    // First reset all block coordinates
    for (int32_t blockId = 0; blockId < blockExists.size(); blockId++) {
      // Note: We reinitialize ALL blocks, not just mutated ones, to ensure
      // proper consistency throughout the coordinate system
      try {
        manager.reinitializeBlockCoordinates(blockId);
      } catch (const std::exception &e) {
        warn("Exception reinitializing block {}: {}. Continuing with next "
             "block.",
             blockId, e.what());
      }
    }

    // Now undo all mutations
    try {
      seed_annotated_tree::undoMutations(T, node, nodeState.blockMutationInfo,
                                         nodeState.nucleotideMutationInfo,
                                         traverser, nodeState.gapRunBacktracks);

      // Validate gap map again after undoing mutations
      if (!manager.validateGapMap(true)) {
        warn("Gap map had issues after undoing mutations - attempted to fix "
             "automatically");
      }
    } catch (const std::exception &e) {
      err("Error during undoMutations: {}. Will attempt recovery.", e.what());

      // Recovery: Try to optimize gap map and ensure blocks are properly
      // initialized
      try {
        manager.optimizeGapMap();

        // Reinitialize again after optimization
        for (int32_t blockId = 0; blockId < blockExists.size(); blockId++) {
          if (blockExists[blockId]
                  .first) { // Only reinitialize live blocks during recovery
            manager.reinitializeBlockCoordinates(blockId);
          }
        }
      } catch (...) {
        err("Critical error during recovery. Coordinate system may be "
            "inconsistent.");
      }
    }

    // Calculate and report time taken
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time)
                        .count();

    msg("Node {} processing completed in {}ms", node->identifier, duration);

  } catch (const std::exception &e) {
    err("Critical error processing node {}: {}", node->identifier, e.what());
    msg("Attempting recovery to continue with siblings...");

    // Try to restore state as much as possible before returning
    try {
      // First reset seed changes
      undoSeedChanges(state, nodeState, traverser,
                      std::set<seed_annotated_tree::SeedChange,
                               seed_annotated_tree::SeedChangeComparator>());

      // Then undo mutations
      seed_annotated_tree::undoMutations(T, node, nodeState.blockMutationInfo,
                                         nodeState.nucleotideMutationInfo,
                                         traverser, nodeState.gapRunBacktracks);
    } catch (...) {
      err("Error during recovery after exception. State may be inconsistent.");
    }

    // Rethrow if this is a critical issue that should stop processing
    // throw; // Comment this out if you want to try to continue despite
    // failures
  }
}

void index(panmanUtils::Tree *T, Index::Builder &index, int k, int s) {
  // msg("Starting indexing");

  if (!T) {
    throw std::runtime_error("Null tree pointer provided to indexing");
  }

  // Create owned data structures
  auto sequence = std::make_unique<sequence_t>();
  auto blockExists = std::make_unique<blockExists_t>();
  auto blockStrand = std::make_unique<blockStrand_t>();

  // msg("Setting up sequence data");
  seed_annotated_tree::setupIndexing(*sequence, *blockExists, *blockStrand, T);

  // Create coordinate manager
  auto coordManager = std::make_unique<CoordinateManager>(
      sequence.get(), blockExists.get(), blockStrand.get());

  // msg("Creating traverser");
  coordinates::CoordinateTraverser traverser(
      coordManager->getLeftmostCoordinate(), coordManager.get());

  if (!traverser.validateState()) {
    throw std::runtime_error("Invalid traverser state after initialization");
  }

  // Initialize the gap map by scanning the entire sequence
  // This creates the initial state of gap runs
  msg("Initializing gap map");
  auto start_time = std::chrono::high_resolution_clock::now();

  // Reset traverser to the beginning
  traverser.reset(coordManager->getLeftmostCoordinate());

  int64_t current_gap_start = -1;
  int64_t current_run_length = 0;
  std::map<int64_t, int64_t> gap_runs;

  // Scan through the sequence to identify all gap runs
  while (!traverser.isDone()) {
    const tupleCoord_t *coord = traverser.getCurrent();
    int64_t scalar = coord->scalar;
    bool is_gap = (traverser.getChar() == '-');

    if (is_gap) {
      // Found a gap position
      if (current_gap_start == -1) {
        // Start a new gap run
        current_gap_start = scalar;
        current_run_length = 1;
      } else {
        // Continue the current gap run
        current_run_length++;
      }
    } else {
      // Found a non-gap position
      if (current_gap_start != -1) {
        // End the current gap run
        gap_runs[current_gap_start] = current_run_length;
        current_gap_start = -1;
        current_run_length = 0;
      }
    }

    traverser.next();
  }

  // Don't forget to add the last gap run if it extends to the end
  if (current_gap_start != -1) {
    gap_runs[current_gap_start] = current_run_length;
  }

  // Update the coordinate manager's gap map with our findings
  auto &gap_map = coordManager->getGapMap();
  // Don't try to clear the const map directly
  // gap_map.clear(); // Clear any existing data

  // Use our new method to clear the gap map
  coordManager->clearGapMap();

  // Add our new gap runs
  for (const auto &[start, length] : gap_runs) {
    coordManager->setGapMap(start, length);
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                      end_time - start_time)
                      .count();
  msg("Gap map initialized with {} gap runs in {}ms", gap_runs.size(),
      duration);

  // Reset traverser to the beginning
  traverser.reset(coordManager->getLeftmostCoordinate());

  // msg("Initializing indexing structures");
  seed_annotated_tree::TraversalGlobalState state(
      T, coordManager->getNumBlocks(), coordManager->getNumCoords());

  // Initialize builders
  auto perNodeSeedMutations_Builder =
      index.initPerNodeSeedMutations(T->allNodes.size());
  auto perNodeGapMutations_Builder =
      index.initPerNodeGapMutations(T->allNodes.size());

  // Set up parameters
  PanmapParams params;
  params.k = k;
  params.s = s;
  params.open = false;
  params.t = 1;
  params.l = 1;

  // Store parameters in index
  index.setK(params.k);
  index.setS(params.s);
  index.setT(params.t);
  index.setOpen(params.open);
  index.setL(params.l);

  msg("Starting indexing");

  auto traversal_start = std::chrono::high_resolution_clock::now();

  // Perform indexing traversal starting from root
  indexingTraversal(state, perNodeSeedMutations_Builder,
                    perNodeGapMutations_Builder, params, T, T->root, traverser);

  auto traversal_end = std::chrono::high_resolution_clock::now();
  auto indexing_duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(traversal_end -
                                                            traversal_start)
          .count();

  msg("Indexing completed in {}ms", indexing_duration);
}

void processSeedChanges(
    seed_annotated_tree::TraversalGlobalState &state,
    seed_annotated_tree::NodeLocalState &nodeState,
    coordinates::CoordinateTraverser &traverser, const PanmapParams &params,
    panmanUtils::Tree *T, panmanUtils::Node *node, int32_t lastOnBlock,
    const tupleCoord_t *lastOnCoord,
    std::set<seed_annotated_tree::SeedChange,
             seed_annotated_tree::SeedChangeComparator> &seedChanges) {

  // Clear any existing seed changes
  seedChanges.clear();

  // Pre-allocate a vector for batch collection of seed changes before inserting
  // into set This reduces individual insertions and tree rebalancing operations
  std::vector<seed_annotated_tree::SeedChange> batchChanges;

  // Estimate the typical number of changes to pre-allocate memory
  // This avoids frequent reallocations
  constexpr size_t TYPICAL_CHANGES_PER_RANGE = 100;
  batchChanges.reserve(nodeState.recompRanges.size() *
                       TYPICAL_CHANGES_PER_RANGE);

  auto &manager = traverser.getCoordManager();
  auto &gap_map = manager.getGapMap();
  auto &blockExists = manager.getBlockExists();
  auto &blockStrand = manager.getBlockStrand();

  // EARLY DETECTION: Log processing info
  std::cout << "Processing " << nodeState.recompRanges.size()
            << " ranges for node " << node->identifier << std::endl;

  // Process each recomputation range
  for (auto &range : nodeState.recompRanges) {
    std::cout << "Processing range: " << range.start_scalar << " -> "
              << range.stop_scalar << std::endl;
    bool atGlobalEnd = false;

    // EARLY DETECTION: Validate range before proceeding
    if (range.start_scalar < 0 || range.stop_scalar < 0 ||
        range.start_scalar >= manager.size() ||
        range.stop_scalar >= manager.size()) {
      std::cout << "ERROR: Invalid scalar range [" << range.start_scalar << ", "
                << range.stop_scalar << "] exceeds valid bounds [0, "
                << manager.size() - 1 << "]. Skipping." << std::endl;
      continue;
    }

    if (range.getStop(traverser)->blockId >= lastOnBlock) {
      std::cout << " => global end." << std::endl;
      atGlobalEnd = true;
      range.setStop(lastOnCoord->scalar);
    }

    // Ensure the range is valid (start <= stop)
    if (range.start_scalar > range.stop_scalar) {
      std::cout
          << "Warning: Invalid range (start > stop), swapping start and stop."
          << std::endl;
      std::swap(range.start_scalar, range.stop_scalar);
    }

    // Additional validation to ensure minimal valid range
    if (range.start_scalar == range.stop_scalar) {
      std::cout << "Warning: Empty range detected, skipping." << std::endl;
      continue;
    }

    // EARLY DETECTION: Validate start and stop coordinates
    const tupleCoord_t *start_coord = manager.getTupleCoord(range.start_scalar);
    const tupleCoord_t *stop_coord = manager.getTupleCoord(range.stop_scalar);

    if (!start_coord || !stop_coord) {
      std::cout << "ERROR: Failed to get valid tuple coordinates for range ["
                << range.start_scalar << ", " << range.stop_scalar
                << "]. Skipping." << std::endl;
      continue;
    }

    if (!manager.validateCoordinate(*start_coord) ||
        !manager.validateCoordinate(*stop_coord)) {
      std::cout << "ERROR: Range has invalid coordinates. Start: "
                << start_coord->blockId << "," << start_coord->nucPos << ","
                << start_coord->nucGapPos << " Stop: " << stop_coord->blockId
                << "," << stop_coord->nucPos << "," << stop_coord->nucGapPos
                << ". Skipping." << std::endl;
      continue;
    }

    // Get mutated sequence to recalculate seeds
    std::string seq;
    std::vector<int64_t> coords;
    std::vector<int64_t> gaps;
    std::vector<int32_t> deadBlocks;
    std::cout << "Calling getNucleotideSequenceFromBlockCoordinates"
              << std::endl;

    // Additional safety - log range values for debugging
    std::cout << "Range values pre-call: start_scalar=" << range.start_scalar
              << ", stop_scalar=" << range.stop_scalar << std::endl;

    seed_annotated_tree::getNucleotideSequenceFromBlockCoordinates(
        seq, coords, gaps, deadBlocks, manager, range.start_scalar,
        range.stop_scalar, T, node);

    // EARLY DETECTION: Verify sequence extraction results
    if (seq.empty() && !gaps.empty()) {
      std::cout << "Warning: Extracted sequence is empty but found "
                << gaps.size() << " gaps." << std::endl;
    }

    if (seq.empty() && coords.empty() && gaps.empty()) {
      std::cout << "Warning: No content extracted from range. Skipping."
                << std::endl;
      continue;
    }

    // Enhanced gap tracking: Track gap runs in this range
    std::map<int64_t, int64_t> local_gap_runs;
    int64_t current_gap_start = -1;
    int64_t current_run_length = 0;

    // Determine traversal direction for gaps
    bool isForwardTraversal = true;
    if (!coords.empty() && coords.size() > 1) {
      isForwardTraversal = coords[1] > coords[0];
    }

    // Process and track gaps
    if (!gaps.empty() && !coords.empty()) {
      // Sort coords and gaps for easier sequential processing
      std::sort(coords.begin(), coords.end());
      std::sort(gaps.begin(), gaps.end());

      // Find gap runs by examining all positions
      std::vector<int64_t> all_positions;
      all_positions.reserve(coords.size() + gaps.size());
      all_positions.insert(all_positions.end(), coords.begin(), coords.end());
      all_positions.insert(all_positions.end(), gaps.begin(), gaps.end());
      std::sort(all_positions.begin(), all_positions.end());

      // Identify gap runs
      for (size_t i = 0; i < all_positions.size(); i++) {
        int64_t scalar = all_positions[i];
        bool is_gap = std::binary_search(gaps.begin(), gaps.end(), scalar);

        if (is_gap) {
          // Found a gap position
          if (current_gap_start == -1) {
            // Start a new gap run
            current_gap_start = scalar;
            current_run_length = 1;
          } else {
            // Continue the current gap run
            current_run_length++;
          }
        } else {
          // Found a non-gap position
          if (current_gap_start != -1) {
            // End the current gap run
            local_gap_runs[current_gap_start] = current_run_length;
            current_gap_start = -1;
            current_run_length = 0;
          }
        }
      }

      // Don't forget to add the last gap run if it extends to the end
      if (current_gap_start != -1) {
        local_gap_runs[current_gap_start] = current_run_length;
      }

      // Update the gap map with our local findings
      for (const auto &[gap_start, gap_length] : local_gap_runs) {
        // EARLY DETECTION: Validate gap position and length
        if (gap_start < 0 || gap_start >= manager.size()) {
          std::cout << "ERROR: Gap start position " << gap_start
                    << " is out of valid range [0, " << manager.size() - 1
                    << "]. Skipping this gap run." << std::endl;
          continue;
        }

        if (gap_length <= 0) {
          std::cout << "ERROR: Invalid gap length " << gap_length
                    << " at position " << gap_start << ". Skipping."
                    << std::endl;
          continue;
        }

        // Create gap map update
        auto existing_it = gap_map.find(gap_start);
        if (existing_it != gap_map.end()) {
          // Record existing value for backtracking
          nodeState.gapRunBacktracks.emplace_back(
              false, std::make_pair(gap_start, existing_it->second));

          // Update if different
          if (existing_it->second != gap_length) {
            // Update gap map and record the change
            manager.setGapMap(gap_start, gap_length);
            nodeState.gapMapUpdates.emplace_back(
                false, std::make_pair(gap_start, gap_length));
          }
        } else {
          // Record that this is a new entry
          nodeState.gapRunBacktracks.emplace_back(true,
                                                  std::make_pair(gap_start, 0));

          // Update gap map and record the change
          manager.setGapMap(gap_start, gap_length);
          nodeState.gapMapUpdates.emplace_back(
              false, std::make_pair(gap_start, gap_length));
        }
      }

      // Look for gap runs in the existing map that should be removed
      // (because they're now nucleotides or part of modified gap runs)
      std::vector<int64_t> to_remove;
      for (auto it = gap_map.lower_bound(range.start_scalar);
           it != gap_map.end() && it->first <= range.stop_scalar; ++it) {
        if (local_gap_runs.find(it->first) == local_gap_runs.end()) {
          to_remove.push_back(it->first);
          // Record this removal for backtracking
          nodeState.gapRunBacktracks.emplace_back(
              false, std::make_pair(it->first, it->second));
          // Record the gap map update
          nodeState.gapMapUpdates.emplace_back(
              true, std::make_pair(it->first, it->second));
        }
      }

      // Remove obsolete gap runs
      for (auto pos : to_remove) {
        manager.removeFromGapMap(pos);
      }
    }

    // Remove seeds in dead blocks - THIS IS THE BOTTLENECK SECTION
    // Optimize by batching the creation of seed changes and only processing
    // positions with active seeds
    for (int deadBlock : deadBlocks) {
      // EARLY DETECTION: Validate block ID
      if (deadBlock < 0 || deadBlock >= state.BlocksToSeeds.size()) {
        std::cout << "ERROR: Invalid dead block ID " << deadBlock
                  << " (max: " << state.BlocksToSeeds.size() - 1
                  << "). Skipping." << std::endl;
        continue;
      }

      // Get all positions affected by this dead block in one pass
      const auto &blockPositions = state.BlocksToSeeds[deadBlock];

      // Create a filtered list of only positions with active seeds
      // This avoids checking positions that don't have seeds later
      std::vector<int> activePositions;
      activePositions.reserve(blockPositions.size() / 2); // Rough estimate

      for (auto &pos : blockPositions) {
        // EARLY DETECTION: Validate position
        if (pos < 0 || pos >= state.onSeedsHash.size()) {
          std::cout << "ERROR: Invalid position " << pos
                    << " for seed in block " << deadBlock << ". Skipping."
                    << std::endl;
          continue;
        }

        if (state.onSeedsHash[pos].has_value()) {
          activePositions.push_back(pos);
        }
      }

      // Skip this block if no active seeds
      if (activePositions.empty())
        continue;

      // Process only positions with active seeds
      for (int pos : activePositions) {
        auto [oldSeed, oldEndPos, oldIsReverse] =
            state.onSeedsHash[pos].value();

        // Create seed change object directly in the vector
        batchChanges.emplace_back(pos,
                                  true,  // old seed on
                                  false, // new seed off
                                  oldSeed, std::nullopt, oldIsReverse,
                                  std::nullopt, oldEndPos, std::nullopt);
      }
    }

    // Remove seeds that now start as gaps
    for (int gapPos : gaps) {
      // EARLY DETECTION: Validate position
      if (gapPos < 0 || gapPos >= state.onSeedsHash.size()) {
        std::cout << "ERROR: Invalid gap position " << gapPos
                  << " (max: " << state.onSeedsHash.size() - 1 << "). Skipping."
                  << std::endl;
        continue;
      }

      if (state.onSeedsHash[gapPos].has_value()) {
        auto [oldSeed, oldEndPos, oldIsReverse] =
            state.onSeedsHash[gapPos].value();
        batchChanges.emplace_back(gapPos,
                                  true,  // old seed on
                                  false, // new seed off
                                  oldSeed, std::nullopt, oldIsReverse,
                                  std::nullopt, oldEndPos, std::nullopt);
      }
    }

    // Process sequence and build new seeds if long enough
    if (seq.size() >= params.k) {
      std::cout << "Processing kmers: seq.size()=" << seq.size()
                << ", coords.size()=" << coords.size() << std::endl;
      try {
        std::vector<std::tuple<size_t, bool, bool, int64_t>> kmers =
            seeding::rollingSyncmers(seq, params.k, params.s, params.open,
                                     params.t, true);
        std::cout << "Generated " << kmers.size() << " kmers" << std::endl;

        // EARLY DETECTION: Check for reasonable kmer count
        if (kmers.empty() && seq.size() > params.k) {
          std::cout << "Warning: No kmers generated from sequence of length "
                    << seq.size() << " with k=" << params.k << std::endl;
        }

        if (kmers.size() > seq.size()) {
          std::cout << "Warning: More kmers (" << kmers.size()
                    << ") than sequence positions (" << seq.size()
                    << "). This may indicate an issue." << std::endl;
        }

        for (int64_t i = 0; i < kmers.size(); ++i) {
          const auto &[hash, isReverse, isSeed, startPos] = kmers[i];

          // Extensive bounds checking
          if (startPos < 0 || startPos >= seq.size()) {
            std::cout << "Warning: kmer startPos " << startPos
                      << " out of bounds for seq size " << seq.size()
                      << std::endl;
            continue;
          }

          if (startPos >= coords.size()) {
            std::cout << "Warning: kmer startPos " << startPos
                      << " out of bounds for coords size " << coords.size()
                      << std::endl;
            continue;
          }

          // Check if coords[startPos] is a valid index for onSeedsHash
          if (coords[startPos] < 0 ||
              coords[startPos] >= state.onSeedsHash.size()) {
            std::cout << "Warning: coords[startPos] = " << coords[startPos]
                      << " out of bounds for onSeedsHash size "
                      << state.onSeedsHash.size() << std::endl;
            continue;
          }

          // Check endpoint coordinates too
          if (i + params.k - 1 >= coords.size()) {
            std::cout << "Warning: end coordinate index " << (i + params.k - 1)
                      << " out of bounds for coords size " << coords.size()
                      << std::endl;
            continue;
          }

          // Now safely proceed with the operation
          bool inMap = state.onSeedsHash[coords[startPos]].has_value();

          if (!inMap && isSeed) {
            // Add new seed - with bounds check
            if (i + params.k - 1 < coords.size()) {
              batchChanges.emplace_back(coords[startPos],
                                        false, // old seed off
                                        true,  // new seed on
                                        std::nullopt, hash, std::nullopt,
                                        isReverse, std::nullopt,
                                        coords[i + params.k - 1]);
            } else {
              std::cout << "Warning: Skipping new seed at " << startPos
                        << " due to coords bounds" << std::endl;
            }
          } else if (inMap && !isSeed) {
            // Remove existing seed
            auto [oldSeed, oldEndPos, oldIsReverse] =
                state.onSeedsHash[coords[startPos]].value();
            batchChanges.emplace_back(coords[startPos],
                                      true,  // old seed on
                                      false, // new seed off
                                      oldSeed, std::nullopt, oldIsReverse,
                                      std::nullopt, oldEndPos, std::nullopt);
          } else if (inMap && isSeed) {
            // Replace seed if changed - with bounds check
            if (i + params.k - 1 < coords.size()) {
              auto [oldSeed, oldEndPos, oldIsReverse] =
                  state.onSeedsHash[coords[startPos]].value();
              if (hash != oldSeed || isReverse != oldIsReverse ||
                  coords[i + params.k - 1] != oldEndPos) {
                batchChanges.emplace_back(coords[startPos],
                                          true, // old seed on
                                          true, // new seed on
                                          oldSeed, hash, oldIsReverse,
                                          isReverse, oldEndPos,
                                          coords[i + params.k - 1]);
              }
            } else {
              std::cout << "Warning: Skipping seed update at " << startPos
                        << " due to coords bounds" << std::endl;
            }
          }
        }
      } catch (const std::exception &e) {
        std::cout << "Exception in kmer processing: " << e.what() << std::endl;
        // Continue processing other ranges even if this one fails
      }
    }
  }

  // Now we've built all the changes in a vector, sort them first to match the
  // set's ordering
  std::sort(batchChanges.begin(), batchChanges.end(),
            seed_annotated_tree::SeedChangeComparator());

  // EARLY DETECTION: Log total number of changes before adding to set
  std::cout << "Total seed changes to process: " << batchChanges.size()
            << std::endl;

  // Insert all changes into the set at once
  seedChanges.insert(batchChanges.begin(), batchChanges.end());

  // Process all seed changes (sorted in set)
  for (const auto &change : seedChanges) {
    const auto &[pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse,
                 newIsReverse, oldEndPos, newEndPos] = change;
    processSeedChange(state, pos, oldVal, newVal, oldSeed, newSeed,
                      oldIsReverse, newIsReverse, oldEndPos, newEndPos,
                      traverser);
  }
}

void undoSeedChanges(
    seed_annotated_tree::TraversalGlobalState &state,
    seed_annotated_tree::NodeLocalState &nodeState,
    coordinates::CoordinateTraverser &traverser,
    const std::set<seed_annotated_tree::SeedChange,
                   seed_annotated_tree::SeedChangeComparator> &seedChanges) {

  // Undo seed changes in reverse order using common function
  for (auto it = seedChanges.rbegin(); it != seedChanges.rend(); ++it) {
    const auto &[pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse,
                 newIsReverse, oldEndPos, newEndPos] = *it;
    seed_annotated_tree::undoSeedChange(state, pos, oldVal, newVal, oldSeed,
                                        newSeed, oldIsReverse, newIsReverse,
                                        oldEndPos, newEndPos, traverser);
  }
}

void storeSeedChanges(
    ::capnp::List<SeedMutations>::Builder &perNodeSeedMutations,
    ::capnp::List<GapMutations>::Builder &perNodeGapMutations, int64_t dfsIndex,
    const seed_annotated_tree::NodeLocalState &nodeState,
    const std::set<seed_annotated_tree::SeedChange,
                   seed_annotated_tree::SeedChangeComparator> &seedChanges) {

  // Process seed changes into ternary encoding
  std::vector<int64_t> basePositions;
  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> masks_all;

  // Work with iterators since we can't index directly into a set
  auto it = seedChanges.rbegin();

  while (it != seedChanges.rend()) {
    // Get current position as our anchor
    const auto &[curPos, curOldVal, curNewVal, curOldSeed, curNewSeed,
                 curOldIsReverse, curIsReverse, curOldEndPos, curEndPos] = *it;
    basePositions.push_back(curPos);
    std::vector<std::pair<uint64_t, uint64_t>> masks;

    // Process all changes near this position (within 32 positions)
    auto maskIt = it;

    while (maskIt != seedChanges.rend()) {
      const auto &[maskPos, maskOldVal, maskNewVal, maskOldSeed, maskNewSeed,
                   maskOldIsReverse, maskIsReverse, maskOldEndPos, maskEndPos] =
          *maskIt;

      // Break if we're too far from current position
      if (curPos - maskPos >= 32) {
        break;
      }

      // Determine ternary code for this change
      int8_t ternaryNumber;
      if (maskOldVal && maskNewVal)
        ternaryNumber = 2; // changed/inserted
      else if (!maskOldVal && maskNewVal)
        ternaryNumber = 2; // changed/inserted
      else if (maskOldVal && !maskNewVal)
        ternaryNumber = 1; // deleted
      else
        ternaryNumber = 0; // same

      // Add to masks with position offset
      masks.emplace_back(std::make_pair(ternaryNumber, curPos - maskPos));

      // Move to next change
      ++maskIt;
    }

    // Save all masks for this position
    masks_all.push_back(masks);

    // Move main iterator to where maskIt left off
    it = maskIt;
  }

  // Store seed changes in index
  auto basePositionsBuilder =
      perNodeSeedMutations[dfsIndex].initBasePositions(basePositions.size());
  auto perPosMasksBuilder =
      perNodeSeedMutations[dfsIndex].initPerPosMasks(masks_all.size());

  for (int i = 0; i < masks_all.size(); i++) {
    const auto &masks = masks_all[i];
    uint64_t tritMask = 0;
    for (const auto &[ternaryNumber, offset] : masks) {
      // Since we already validated the offset when building the masks,
      // this is just an additional safeguard
      if (offset < 32) {
        tritMask |= (ternaryNumber & 0x3) << (offset * 2);
      }
    }
    perPosMasksBuilder.set(masks_all.size() - i - 1, tritMask);
    basePositionsBuilder.set(masks_all.size() - i - 1, basePositions[i]);
  }

  // Store gap mutations
  auto nodeGapBuilder = perNodeGapMutations[dfsIndex];
  auto gapMutationsBuilder =
      nodeGapBuilder.initDeltas(nodeState.gapMapUpdates.size());

  for (int i = 0; i < nodeState.gapMapUpdates.size(); ++i) {
    const auto &[erase, range] = nodeState.gapMapUpdates[i];
    if (erase) {
      gapMutationsBuilder[i].setPos(range.first);
      gapMutationsBuilder[i].initMaybeValue().setNone();
    } else {
      gapMutationsBuilder[i].setPos(range.first);
      gapMutationsBuilder[i].initMaybeValue().setValue(range.second);
    }
  }
}

} // namespace indexing
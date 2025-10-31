#include "seed_annotated_tree.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "seq_utils.hpp"
#include "coordinates.hpp"
#include "visualization.hpp"
#include <chrono>
#include <spdlog/fmt/bundled/format.h>

namespace seed_annotated_tree {

using namespace logging; // Add this to use logging functions directly
using namespace coordinates;

// Define the global debug variable that's declared as extern in the header
bool debug = false;

// Use the typedefs from coordinates namespace
using blockExists_t = blockExists_t;
using blockStrand_t = blockStrand_t;
using sequence_t = sequence_t;
using tupleCoord_t = tupleCoord_t;

// Time tracking functions
static std::chrono::time_point<std::chrono::high_resolution_clock> global_timer 
    = std::chrono::high_resolution_clock::now();

// Change from bool debug to a function
bool isDebugEnabled() {
    return debug;  // Use the namespace global variable
}

void setDebug(bool enabled) {
    debug = enabled;
}

// Helper function for debug messages
template <typename... Args>
void debug_msg(const std::string& fmt, Args&&... args) {
    if (isDebugEnabled()) {
        logging::debug(fmt, std::forward<Args>(args)...);
    }
}

// Helper method that uses a hybrid approach - scalar logic within blocks,
// traversal order at boundaries
int getValidNucleotidesEfficiently(
    CoordinateManager &manager, int64_t startPos, char *buffer,
    int maxChars, int64_t &resultEndPos, const blockExists_t &blockExists,
    const blockStrand_t &blockStrand, const sequence_t &sequence) {
  // Get starting coordinate
  const tupleCoord_t *startCoord = manager.getTupleCoord(startPos);
  if (!startCoord)
    return 0;

  // Validate block existence
  if (!blockExists[startCoord->blockId].first) {
    return 0;
  }

  // Use gap map to check if starting in a gap
  if (manager.isGapPosition(startPos)) {
    // Check if we're in an inverted block to determine traversal direction
    bool isBlockInverted = !blockStrand[startCoord->blockId].first;

    // Starting in a gap - use appropriate function based on block orientation
    const tupleCoord_t *nextNonGap =
        isBlockInverted ? manager.getPreviousNonGapPosition(startPos)
                        : // Inverted blocks move backward in scalar space
            manager.getNextNonGapPosition(
                startPos); // Normal blocks move forward in scalar space

    if (!nextNonGap)
      return 0; // No more non-gap positions
    startCoord = nextNonGap;
    startPos = nextNonGap->scalar;
  }

  int count = 0;
  int64_t maxScalar = manager.size() - 1;

  // Start with current coordinate
  const tupleCoord_t *coord = startCoord;

  // Main loop to collect nucleotides
  while (coord && count < maxChars && coord->scalar <= maxScalar) {
    // Cache current block info for multiple uses
    int32_t currentBlockId = coord->blockId;
    bool isBlockInverted = !blockStrand[currentBlockId].first;

    // Skip if block doesn't exist
    if (!blockExists[currentBlockId].first) {
      // Find next block according to traversal order
      int32_t nextBlockId =
          isBlockInverted ? currentBlockId - 1 : currentBlockId + 1;
      if (nextBlockId < 0 || static_cast<size_t>(nextBlockId) >= blockExists.size())
        break;

      // Move to correct position in next block
      const tupleCoord_t *nextBlockCoord =
          isBlockInverted ? manager.getLastCoordinateInBlock(nextBlockId)
                          : manager.getFirstCoordinateInBlock(nextBlockId);

      if (!nextBlockCoord)
        break;
      coord = nextBlockCoord;
      continue;
    }

    // Get the range of this block for fast within-block processing
    std::pair<int64_t, int64_t> blockRange =
        manager.getBlockRange(currentBlockId);

    // Fast path: process all positions within this block using scalar logic
    // This avoids pointer chasing when we're inside a block
    int64_t currentPos = coord->scalar;
    int64_t blockEnd = isBlockInverted ? blockRange.first : blockRange.second;

    // Determine direction for scalar increments based on block orientation
    int scalarStep = isBlockInverted ? -1 : 1;

    // Process all positions in this block with scalar logic
    while (count < maxChars && ((isBlockInverted && currentPos >= blockEnd) ||
                                (!isBlockInverted && currentPos <= blockEnd))) {

      // Check if current position is in a gap using the gap map
      if (manager.isGapPosition(currentPos)) {
        // Skip this gap run efficiently based on block orientation
        if (isBlockInverted) {
          // In an inverted block, use getPreviousNonGapPosition for consistent
          // behavior
          const tupleCoord_t *prevNonGap =
              manager.getPreviousNonGapPosition(currentPos);
          if (prevNonGap && prevNonGap->scalar >= blockEnd) {
            currentPos = prevNonGap->scalar;
            continue;
          } else {
            // Either no previous non-gap or it's outside our block
            break;
          }
        } else {
          // In a forward block, use getNextNonGapPosition for consistent
          // behavior
          const tupleCoord_t *nextNonGap =
              manager.getNextNonGapPosition(currentPos);
          if (nextNonGap && nextNonGap->scalar <= blockEnd) {
            currentPos = nextNonGap->scalar;
            continue;
          } else {
            // Either no next non-gap or it's outside our block
            break;
          }
        }
      }

      // Get coordinate at current position
      const tupleCoord_t *currentCoord = manager.getTupleCoord(currentPos);
      if (!currentCoord)
        break;

      // Sanity check - make sure we're still in the same block
      if (currentCoord->blockId != currentBlockId)
        break;

      // Get character at this position
      char c;
      if (currentCoord->nucGapPos == -1) {
        c = sequence[currentBlockId].first[currentCoord->nucPos].first;
      } else {
        c = sequence[currentBlockId]
                .first[currentCoord->nucPos]
                .second[currentCoord->nucGapPos];
      }

      // Apply complementation for inverted blocks
      if (isBlockInverted) {
        c = seq_utils::getComplementCharacter(c);
      }

      // Only collect valid nucleotides
      if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' ||
          c == 'c' || c == 'g' || c == 't') {
        buffer[count++] = c;
        resultEndPos = currentPos;
      }

      // Move to next position within block using scalar logic
      currentPos += scalarStep;
    }

    // After processing the block, we need traversal logic to move to the next
    // block Find the next/prev block based on traversal order
    int32_t nextBlockId =
        isBlockInverted ? currentBlockId - 1 : currentBlockId + 1;
    if (nextBlockId < 0 || static_cast<size_t>(nextBlockId) >= blockExists.size())
      break;

    // Skip non-existent blocks
    while (nextBlockId >= 0 && static_cast<size_t>(nextBlockId) < blockExists.size() &&
           !blockExists[nextBlockId].first) {
      nextBlockId = isBlockInverted ? nextBlockId - 1 : nextBlockId + 1;
      if (nextBlockId < 0 || static_cast<size_t>(nextBlockId) >= blockExists.size()) {
        nextBlockId = -1;
        break;
      }
    }

    if (nextBlockId < 0 || static_cast<size_t>(nextBlockId) >= blockExists.size())
      break;

    // Move to the appropriate position in the next block
    bool isNextBlockInverted = !blockStrand[nextBlockId].first;
    const tupleCoord_t *nextBlockCoord =
        isNextBlockInverted ? manager.getLastCoordinateInBlock(nextBlockId)
                            : manager.getFirstCoordinateInBlock(nextBlockId);

    if (!nextBlockCoord)
      break;
    coord = nextBlockCoord;
  }

  return count;
}

bool getSeedAt(coordinates::CoordinateTraverser &traverser,
                                    coordinates::CoordinateManager &manager,
                                    size_t &resultHash, bool &resultIsReverse,
                                    int64_t &resultEndPos, const int64_t &pos,
                                    panmanUtils::Tree *T, const int32_t &k) {
  
  // Use const references to avoid copies and clarify non-ownership
  const auto &sequence = manager.getSequence();
  const auto &blockExists = manager.getBlockExists();
  const auto &blockStrand = manager.getBlockStrand();
  const auto &gap_map = manager.getGapMap();

  // expect more non-gaps that gaps -- should also try reverse
  if (__builtin_expect(manager.isGapPosition(pos), 0)) {
    return false;
  }

  // Preallocate a fixed buffer for performance
  alignas(64) static thread_local char
      buffer[256]; // Align to cache line and make thread local

  // Reset traverser to start at our coordinate
  const tupleCoord_t *startCoord = manager.getTupleCoord(pos);
  if (__builtin_expect(!startCoord, 0))
    throw std::runtime_error("Invalid coordinate");

  // expect more blocks to exist than not
  if (__builtin_expect(!blockExists[startCoord->blockId].first, 0)) {
    return false;
  }

  // Check if we're starting in an inverted block
  bool startInInvertedBlock = !blockStrand[startCoord->blockId].first;

  // If we're in an inverted block, ensure we have enough room to collect k
  // bases This is important because inverted blocks read in reverse scalar
  // order
  if (startInInvertedBlock) {
    // In inverted blocks, we need to check if we have k positions available in
    // the traversal direction (which is decreasing scalar for inverted blocks)
    int availablePositions = 0;
    const tupleCoord_t *checkCoord = startCoord;

    // Count available non-gap positions in traversal direction
    while (checkCoord && availablePositions < k) {
      if (manager.isGapPosition(checkCoord->scalar)) {
        checkCoord = manager.getPreviousNonGapPosition(checkCoord->scalar);
      } else {
        availablePositions++;
        // Get previous position respecting block boundaries and inversions
        if (checkCoord->prev) {
          checkCoord = checkCoord->prev;
        } else {
          break;
        }
      }
    }

    if (availablePositions < k) {
      debug_msg("Not enough positions ({}) available in inverted block from scalar "
            "{} to extract k={} bases",
            availablePositions, startCoord->scalar, k);
      return false;
    }
  }

  traverser.reset(startCoord);

  // Get k non-gap characters and store in buffer directly in one pass
  // This now correctly handles both forward and inverted blocks
  bool foundEnough = traverser.skipToNthNonGap(k, buffer);

  // Check if we found enough characters
  if (__builtin_expect(!foundEnough, 0)) {
    return false;
  }

  // Get the end position
  resultEndPos = traverser.getCurrent()->scalar;

  // Calculate hash without string construction (performance critical)
  try {
    // Avoid creating std::string - work directly with the buffer
    auto [fHash, rHash] = seeding::hashSeq(std::string_view(buffer, k));
    if (fHash < rHash) {
      resultHash = fHash;
      resultIsReverse = false;
    } else {
      resultHash = rHash;
      resultIsReverse = true;
    }
    return true;
  } catch (...) {
    return false;
  }
}

// coords (blockId, nucPosition, nucGapPosition)
// Sequence includes both boundary coordinates
// Sequence, scalarCoords of sequence, scalarCoords of Gaps, dead blocks
void getNucleotideSequenceFromBlockCoordinates(
    std::string &seq, std::vector<int64_t> &coords, std::vector<int64_t> &gaps,
    std::vector<int32_t> &deadBlocks, coordinates::CoordinateManager &manager,
    int64_t start_scalar, int64_t stop_scalar, const Tree *T,
    const Node *node) {

  // Clear output vectors
  seq.clear();
  coords.clear();
  gaps.clear();
  deadBlocks.clear();

  // Validate input scalars
  if (start_scalar < 0 || stop_scalar < 0) {
    std::cout << "Error: Invalid negative scalar values: start_scalar="
              << start_scalar << ", stop_scalar=" << stop_scalar << std::endl;
    return;
  }

  // Normalize scalar order if needed - this is for handling inverted blocks
  // where traversal goes backward
  if (start_scalar > stop_scalar) {
    warn("Normalizing inverted range: start ({}) > stop ({})", 
         start_scalar, stop_scalar);

    // Visualize block structure and gap map before normalization
    std::cout << "Before normalization:" << std::endl;
    std::string rangeInfo = "Range_" + std::to_string(start_scalar) + "_" + std::to_string(stop_scalar);
    visualization::printNodeVisualization(
        rangeInfo, 
        manager.getBlockExists(), 
        manager.getBlockStrand(),
        manager.getGapMap(),
        manager.getNumCoords());

    // Explicit swap for clarity - this ensures we always process from lower to
    // higher scalar This is needed because we're iterating over the coords
    // array in scalar order
    int64_t temp = start_scalar;
    start_scalar = stop_scalar;
    stop_scalar = temp;
    
    // Visualize after normalization
    std::cout << "After normalization:" << std::endl;
    rangeInfo = "Range_" + std::to_string(start_scalar) + "_" + std::to_string(stop_scalar);
    std::cout << "Normalized range: " << rangeInfo << std::endl;
  }

  


  // Get references to block existence and strand information
  const blockExists_t &blockExists = manager.getBlockExists();
  const blockStrand_t &blockStrand = manager.getBlockStrand();

  // Get start and stop coordinates
  const tupleCoord_t *start = manager.getTupleCoord(start_scalar);
  const tupleCoord_t *stop = manager.getTupleCoord(stop_scalar);

  if (!start || !stop) {
    std::cout << "Error: Invalid start or stop coordinate" << std::endl;
    return;
  }

  // Validate start and stop coordinates
  if (!manager.validateCoordinate(*start) ||
      !manager.validateCoordinate(*stop)) {
    std::cout << "Error: Start or stop coordinate validation failed"
              << std::endl;
    return;
  }

  // Create traverser
  coordinates::CoordinateTraverser traverser(start, &manager);

  // Set end point for traverser
  traverser.reset(start, stop);

  // Track seen dead blocks to avoid duplicates
  std::unordered_set<int32_t> seenDeadBlocks;

  // Calculate traversal length for buffer allocation
  int64_t traversal_length = stop_scalar - start_scalar + 1;

  // Handle inverted range
  if (traversal_length <= 0) {
    // std::cout << "Warning: Inverted or zero-length range detected"
    //           << " (start_scalar=" << start_scalar
    //           << ", stop_scalar=" << stop_scalar << ")" << std::endl;
    traversal_length = 1; // Set minimum valid length
  }

  // Estimate a reasonable initial capacity (assuming about 20% of positions are
  // gaps)
  seq.reserve(traversal_length);
  coords.reserve(traversal_length);
  gaps.reserve(traversal_length / 5);
  deadBlocks.reserve(10); 

  // CRITICAL FIX: Create a new visit counter map for each function call
  // This ensures coordinates visited in previous calls don't trigger false
  // cycles
  std::unordered_map<int64_t, int> visitCount;
  int maxVisits = 3; // Maximum allowed visits to same coordinate

  // Also track positions seen to detect direct back-and-forth cycling
  const tupleCoord_t *prevCoord = nullptr;

  while (!traverser.isDone()) {

    const tupleCoord_t *currCoord = traverser.getCurrent();

    // std::cout << "currCoord = " << currCoord->blockId << "," <<
    // currCoord->nucPos << "," << currCoord->nucGapPos << " = " <<
    // currCoord->scalar << std::endl;
    if (!currCoord) {
      std::cout << "Error: Null coordinate encountered during traversal"
                << std::endl;
      break;
    }

    // Enhanced cycle detection - using the non-static visit counter
    visitCount[currCoord->scalar]++;
    if (visitCount[currCoord->scalar] > maxVisits) {
      std::cout << "Error: Detected cycle in coordinate traversal at position "
                << currCoord->blockId << "," << currCoord->nucPos << ","
                << currCoord->nucGapPos << " (scalar=" << currCoord->scalar
                << "). Breaking out of loop." << std::endl;
      break;
    }

    // Additional check for direct back-and-forth cycles (ping-pong)
    if (prevCoord && prevCoord->next == currCoord &&
        currCoord->next == prevCoord) {
      std::cout
          << "Error: Detected direct back-and-forth cycle between positions "
          << prevCoord->scalar << " and " << currCoord->scalar
          << ". Breaking out." << std::endl;
      break;
    }
    prevCoord = currCoord;

    int64_t scalar = currCoord->scalar;
    if (scalar > manager.getNumCoords()) {
      std::cout << "Error: Scalar value " << scalar
                << " exceeds coordinate count " << manager.getNumCoords()
                << std::endl;
      return;
    }

    // Get the character at the current position
    char c = traverser.getChar();

    // Check for end marker 'x'
    if (c == 'x') {
      // std::cout << "End of block reached" << std::endl;
      traverser.next();
      continue;
    }

    // Determine if we're in an inverted block
    bool isInvertedBlock = currCoord->blockId < blockStrand.size() &&
                           !blockStrand[currCoord->blockId].first;

    // Check for valid nucleotide characters
    if (c != '-' && c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N' &&
        c != 'a' && c != 'c' && c != 'g' && c != 't' && c != 'n') {
      // std::cout << "Warning: Invalid character '" << c
      //           << "' encountered at position " << currCoord->blockId << ","
      //           << currCoord->nucPos << "," << currCoord->nucGapPos
      //           << ". Skipping." << std::endl;
      traverser.next();
      continue;
    }

    // Check for consecutive gap positions with same blockId and nucPos
    if (currCoord->nucGapPos > 0) {
      // We're in a gap position within a nucleotide
      int currentBlockId = currCoord->blockId;
      int currentNucPos = currCoord->nucPos;

      // Check if the next position continues the pattern
      const tupleCoord_t *nextCoord = currCoord->next;
      if (nextCoord && nextCoord->blockId == currentBlockId &&
          nextCoord->nucPos == currentNucPos &&
          nextCoord->nucGapPos == currCoord->nucGapPos + 1) {

        // We're in a consecutive gap run - find the end of it
        const tupleCoord_t *endOfRun = nextCoord;
        int skipCount = 1;

        while (endOfRun->next && endOfRun->next->blockId == currentBlockId &&
               endOfRun->next->nucPos == currentNucPos &&
               endOfRun->next->nucGapPos == endOfRun->nucGapPos + 1) {
          endOfRun = endOfRun->next;
          skipCount++;
        }

        // If we found a significant run, skip to the end
        if (skipCount > 10) {
          // std::cout << "Skipping " << skipCount
          //           << " consecutive gap positions at " << currentBlockId << ","
          //           << currentNucPos << std::endl;

          // Add the current gap position
          gaps.push_back(scalar);

          // For handling inverted blocks correctly
          int64_t endScalar = endOfRun->scalar;

          // Add the end position of the run to gaps
          gaps.push_back(endScalar);

          // Log the gap range - useful for debugging
          debug_msg(
              "Recording gap run from scalar={} to scalar={} in block {} ({})",
              scalar, endScalar, currentBlockId,
              isInvertedBlock ? "inverted" : "forward");


          // Skip to the position after the run
          if (endOfRun->next && endOfRun->next->scalar <= stop_scalar) {
            debug_msg("Skipping to position after gap run: blockId={}, nucPos={}, "
                  "nucGapPos={}, scalar={}",
                  endOfRun->next->blockId, endOfRun->next->nucPos,
                  endOfRun->next->nucGapPos, endOfRun->next->scalar);
            traverser.reset(endOfRun->next, stop);
            continue;
          } else {
            // We've reached the end of the sequence or the stop position
            debug_msg("Reached end of sequence or stop position after gap run");
            break;
          }
        }
      }
    }

    // std::cout <<"Curr: " << currCoord->blockId << "," << currCoord->nucPos <<
    // "," << currCoord->nucGapPos << " = " << currCoord->scalar << " -> " <<
    // traverser.getChar() << std::endl;

    // Check if block is alive
    if (blockExists[currCoord->blockId].first) {
      // Live block - process character
      char c = traverser.getChar();

      // Track gap positions
      if (c == '-') {
        gaps.push_back(scalar);

        debug_msg("Found gap at position: blockId={}, nucPos={}, nucGapPos={}, "
              "scalar={}",
              currCoord->blockId, currCoord->nucPos, currCoord->nucGapPos,
              scalar);

        // Check if this is the start of a gap run
        int64_t gap_run_length = manager.getGapRunLength(scalar);
        if (gap_run_length > 1) {
          debug_msg("Detected gap run of length {} at scalar={}", gap_run_length,
                scalar);

          // We're at a gap run - more efficient handling

          // Add all positions in this gap run to the gaps vector
          // But only up to the stop_scalar boundary
          // Limit the number of individual positions we add
          int64_t max_positions_to_add =
              std::min(gap_run_length - 1, static_cast<int64_t>(100));
          for (int64_t i = 1;
               i <= max_positions_to_add && scalar + i <= stop_scalar; ++i) {
            gaps.push_back(scalar + i);
          }

          // Log details about the gap run and our next move
          debug_msg("Gap run: adding {} gap positions and attempting to skip to "
                "position after run",
                max_positions_to_add);

          // Skip directly to the end of this gap run, accounting for whether
          // we're in an inverted block
          bool isInvertedBlock = currCoord->blockId < blockStrand.size() &&
                                 !blockStrand[currCoord->blockId].first;

          const tupleCoord_t *next_non_gap =
              isInvertedBlock ? manager.getPreviousNonGapPosition(scalar)
                              :                          // For inverted blocks
                  manager.getNextNonGapPosition(scalar); // For normal blocks

          if (next_non_gap) {
            debug_msg("Next non-gap position found: blockId={}, nucPos={}, "
                  "nucGapPos={}, scalar={}, in {} block",
                  next_non_gap->blockId, next_non_gap->nucPos,
                  next_non_gap->nucGapPos, next_non_gap->scalar,
                  isInvertedBlock ? "inverted" : "forward");

            if ((isInvertedBlock && next_non_gap->scalar >= start_scalar) ||
                (!isInvertedBlock && next_non_gap->scalar <= stop_scalar)) {
              debug_msg("Skipping to next non-gap position within boundary");
              traverser.reset(next_non_gap, stop);
              continue; // Skip to next iteration with the new position
            } else {
              debug_msg("Next non-gap position is outside boundary, ending "
                    "traversal");
              break;
            }
          } else {
            debug_msg("No next non-gap position found, ending traversal");
            // Either we're at the end of the sequence or past stop_scalar
            // In either case, we're done
            break;
          }
        }
      } else {
        // Only add non-gap characters to the sequence
        // Check for valid nucleotide characters
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't' && c != 'n') {
          std::cout << "Warning: Skipping invalid character '" << c
                    << "' at position " << currCoord->blockId << ","
                    << currCoord->nucPos << "," << currCoord->nucGapPos
                    << std::endl;
        } else {
          seq.push_back(c);
          coords.push_back(scalar);
        }
      }
    } else {
      // Dead block - don't add to sequence but track it
      // Add to deadBlocks if we haven't seen this block already
      if (seenDeadBlocks.find(currCoord->blockId) == seenDeadBlocks.end()) {
        deadBlocks.push_back(currCoord->blockId);
        seenDeadBlocks.insert(currCoord->blockId);
      }

      // Skip to the next block efficiently
      traverser.moveToNextBlock();
      continue;
    }

    // Move to next coordinate, using the skipGaps optimization
    if (traverser.getChar() == '-') {
      // Properly account for inversion when skipping gaps
      traverser
          .skipGaps(); // skipGaps() now internally handles inversions correctly
    } else {
      // Move to next position in traversal order
      traverser.next(); // next() follows the correct pointer which respects the
                        // block orientation
    }
  }

  // If we extracted nothing, log a message
  if (seq.empty() && !gaps.empty()) {
    // std::cout << "Warning: Extracted sequence is empty but found "
    //           << gaps.size() << " gaps" << std::endl;
  }

  // Report statistics
  // std::cout << "Sequence extraction complete: " << seq.size()
  //           << " nucleotides, " << gaps.size() << " gaps, " << deadBlocks.size()
  //           << " dead blocks" << std::endl;

  // Debug visualization of gap map
  visualization::printNodeVisualization(
      "getNucleotideSequenceFromBlockCoordinates",
      manager.getBlockExists(),
      manager.getBlockStrand(),
      manager.getGapMap(),
      manager.size());
}

void setupSequence(const Tree *T, sequence_t &sequence,
                   blockExists_t &blockExists, blockStrand_t &blockStrand) {
  // Validate tree input
  if (!T) {
    throw std::runtime_error("Null tree pointer provided to setupSequence");
  }

  // Get references to key tree structures with validation
  const BlockGapList &blockGaps = T->blockGaps;
  const std::vector<GapList> &gaps = T->gaps;
  const std::vector<Block> &blocks = T->blocks;

  // Validate key components
  if (blocks.empty()) {
    throw std::runtime_error("Tree has no blocks");
  }

  // First find maxBlockId
  int32_t maxBlockId = 0;
  for (size_t i = 0; i < blocks.size(); i++) {
    int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
    maxBlockId = std::max(maxBlockId, primaryBlockId);
  }

  // Resize all structures once to final size
  size_t finalSize = maxBlockId + 1;
  msg("Resizing sequence structures to {} blocks", finalSize);
  
  try {
    sequence.resize(finalSize);
    blockExists.resize(finalSize, {false, {}});
    blockStrand.resize(finalSize, {true, {}});
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        "Failed to resize sequence structures to {} blocks: {}", finalSize, e.what()));
  }

  // Assigning block gaps
  int validBlockGaps = 0;
  int skippedBlockGaps = 0;
  
  for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
    int32_t pos = blockGaps.blockPosition[i];
    if (pos >= finalSize) {
      err("Block gap position {} exceeds max block id {}", pos, maxBlockId);
      skippedBlockGaps++;
      continue;
    }
    try {
      sequence[pos].second.resize(blockGaps.blockGapLength[i]);
      blockExists[pos].second.resize(blockGaps.blockGapLength[i], false);
      blockStrand[pos].second.resize(blockGaps.blockGapLength[i], true);
      validBlockGaps++;
    } catch (const std::exception &e) {
      err("Error resizing for block gap at position {}: {}", pos, e.what());
      skippedBlockGaps++;
    }
  }
  
  if (validBlockGaps > 0) {
    msg("Processed {} block gaps ({} skipped)", validBlockGaps, skippedBlockGaps);
  }

  // Process blocks
  int validBlocks = 0;
  int skippedBlocks = 0;
  int emptyBlocks = 0;
  
  for (size_t i = 0; i < blocks.size(); i++) {
    int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
    int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

    if (primaryBlockId >= finalSize) {
      err("Primary block ID {} exceeds max block id {}", primaryBlockId,
          maxBlockId);
      skippedBlocks++;
      continue;
    }

    bool hasContent = false;
    for (size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
      bool endFlag = false;
      for (size_t k = 0; k < 8; k++) {
        const int nucCode =
            (((blocks[i].consensusSeq[j]) >> (4 * (7 - k))) & 15);

        if (nucCode == 0) {
          endFlag = true;
          break;
        }
        const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

        try {
          if (secondaryBlockId != -1) {
            // Validate secondary block ID is within the resized structure
            if (secondaryBlockId < 0 || secondaryBlockId >= (int32_t)sequence[primaryBlockId].second.size()) {
              err("Invalid secondary block ID {} for primary block {}", 
                  secondaryBlockId, primaryBlockId);
              break;
            }
            sequence[primaryBlockId].second[secondaryBlockId].push_back(
                {nucleotide, {}});
          } else {
            sequence[primaryBlockId].first.push_back({nucleotide, {}});
          }
          hasContent = true;
        } catch (const std::exception &e) {
          err("Error processing nucleotide for block {},{}: {}", 
              primaryBlockId, secondaryBlockId, e.what());
        }
      }
      if (endFlag) {
        break;
      }
    }

    try {
      // Add the end marker
      sequence[primaryBlockId].first.push_back({'x', {}});
      
      // Update tracking
      if (hasContent) {
        validBlocks++;
      } else {
        emptyBlocks++;
      }
    } catch (const std::exception &e) {
      err("Error adding end marker to block {}: {}", primaryBlockId, e.what());
      skippedBlocks++;
    }
  }
  
  msg("Processed {} blocks ({} with content, {} empty, {} skipped)", 
      validBlocks + emptyBlocks, validBlocks, emptyBlocks, skippedBlocks);

  // Assigning nucleotide gaps
  int validGaps = 0;
  int skippedGaps = 0;
  
  for (size_t i = 0; i < gaps.size(); i++) {
    int32_t primaryBId = (gaps[i].primaryBlockId);
    int32_t secondaryBId = (gaps[i].secondaryBlockId);

    if (primaryBId >= finalSize) {
      err("Gap primary block ID {} exceeds max block id {}", primaryBId,
          maxBlockId);
      skippedGaps++;
      continue;
    }

    for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
      int len = gaps[i].nucGapLength[j];
      int pos = gaps[i].nucPosition[j];

      try {
        if (secondaryBId != -1) {
          // Validate secondary block ID and position are within bounds
          if (secondaryBId < 0 || secondaryBId >= (int32_t)sequence[primaryBId].second.size()) {
            err("Invalid secondary block ID {} for gap in primary block {}", 
                secondaryBId, primaryBId);
            skippedGaps++;
            continue;
          }
          
          if (pos < 0 || pos >= (int32_t)sequence[primaryBId].second[secondaryBId].size()) {
            err("Invalid position {} for gap in block {},{}", 
                pos, primaryBId, secondaryBId);
            skippedGaps++;
            continue;
          }
          
          // Resize the gap vector if needed
          if (sequence[primaryBId].second[secondaryBId][pos].second.size() < (size_t)len) {
            sequence[primaryBId].second[secondaryBId][pos].second.resize(len, 'N');
          }
          
          validGaps++;
        } else {
          // Handle primary block gaps
          if (pos < 0 || pos >= (int32_t)sequence[primaryBId].first.size()) {
            err("Invalid position {} for gap in primary block {}", pos, primaryBId);
            skippedGaps++;
            continue;
          }
          
          // Resize the gap vector if needed
          if (sequence[primaryBId].first[pos].second.size() < (size_t)len) {
            sequence[primaryBId].first[pos].second.resize(len, 'N');
          }
          
          validGaps++;
        }
      } catch (const std::exception &e) {
        err("Error processing gap at position {},{},{}: {}", 
            primaryBId, secondaryBId, pos, e.what());
        skippedGaps++;
      }
    }
  }
  
  if (validGaps > 0) {
    msg("Processed {} nucleotide gaps ({} skipped)", validGaps, skippedGaps);
  }
}

static size_t getStart(const std::string &s1, const std::string &s2,
                      size_t window, double threshold) {
  if (s1.empty()) {
    return 0;
  }

  size_t numAlign = 0;
  size_t numMatch = 0;
  size_t beg = 0;
  size_t idx = 0;
  std::queue<size_t> begs;

  while (idx < s1.size()) {
    if (beg == 0) {
      beg = idx;
    } else {
      begs.push(idx);
    }
    if (s1[idx] == s2[idx]) {
      ++numMatch;
    }
    ++numAlign;
    ++idx;

    if (numAlign == window) {
      double pcid = static_cast<double>(numMatch) / static_cast<double>(window);
      if (pcid >= threshold && s1[beg] == s2[beg]) {
        return beg;
      }

      if (s1[beg] == s2[beg]) {
        --numMatch;
      }
      --numAlign;
      beg = begs.front();
      begs.pop();
    }
  }

  return s1.size() - 1;
}

static size_t getEnd(const std::string &s1, const std::string &s2,
                     size_t window, double threshold) {

  if (s1.empty()) {
    return 0;
  }

  size_t numAlign = 0;
  size_t numMatch = 0;
  size_t end = s1.size();
  size_t idx = s1.size() - 1;
  std::queue<size_t> ends;

  while (true) {
    if (s1[idx] == '-' && s2[idx] == '-') {
      if (idx == 0) {
        break;
      }
      --idx;
      continue;
    }
    if (end == s1.size()) {
      end = idx;
    } else {
      ends.push(idx);
    }
    if (s1[idx] == s2[idx]) {
      ++numMatch;
    }
    ++numAlign;

    if (idx == 0) {
      break;
    }
    --idx;

    if (numAlign == window) {
      double pcid = static_cast<double>(numMatch) / static_cast<double>(window);
      if (pcid >= threshold && s1[end] == s2[end]) {
        return end;
      }

      if (s1[end] == s2[end]) {
        --numMatch;
      }
      --numAlign;
      end = ends.front();
      ends.pop();
    }
  }

  return 0;
}

std::pair<size_t, size_t>
getMaskCoorsForMutmat(const std::string &s1,
                      const std::string &s2, size_t window,
                      double threshold) {

  assert(s1.size() == s2.size());
  if (window == 0 || threshold == 0.0) {
    return std::make_pair<size_t, size_t>(0, s1.size() - 1);
  }
  return std::make_pair<size_t, size_t>(getStart(s1, s2, window, threshold),
                                        getEnd(s1, s2, window, threshold));
}

void writeMutationMatrices(const mutationMatrices &mutMat,
                                                std::ofstream &mmfout) {

  for (const std::vector<double> &row : mutMat.submat) {
    for (const double &prob : row) {
      mmfout << prob << " ";
    }
    mmfout << "\n";
  }

  for (const auto &[size, count] : mutMat.insmat) {
    mmfout << size << ":" << count << " ";
  }
  mmfout << "\n";

  for (const auto &[size, count] : mutMat.delmat) {
    mmfout << size << ":" << count << " ";
  }
  mmfout << "\n";
}

void applyMutations(
    blockMutationInfo_t &blockMutationInfo,
    std::vector<coordinates::CoordRange> &recompRanges,
    mutationInfo_t &mutationInfo,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapMapUpdates, 
    panmanUtils::Tree *T,
    panmanUtils::Node *node, coordinates::CoordinateTraverser &traverser,
    bool isPlacement, std::unordered_set<int> &inverseBlockIds,
    std::vector<std::pair<bool, int>> &inverseBlockIdsBacktrack, int k) {
  auto &manager = traverser.getCoordManager();
  auto &blockExists = manager.getBlockExists();
  auto &blockStrand = manager.getBlockStrand();
  auto &sequence = manager.getSequence();

  auto &blocks = T->blocks;
  auto &gaps = T->gaps;
  auto &blockGaps = T->blockGaps;
  auto &sequenceInverted = T->sequenceInverted;
  auto &circularSequences = T->circularSequences;

  // Helper function to normalize gap ranges
  auto normalizeGapRange = [&](int64_t &start, int64_t &end) {
    if (start > end) {
      warn("Normalizing invalid gap range: original start ({}) > end ({})",
           start, end);
      std::swap(start, end);
    }
  };

  // Helper function to validate mutation coordinates before access
  auto isValidMutation = [&](int32_t blockId, int32_t secondaryBlockId,
                             int32_t nucPos, int32_t nucGapPos) -> bool {
    // Check primary block ID
    if (blockId < 0 || blockId >= sequence.size()) {
      warn("Invalid blockId {} when applying mutation (max: {})", blockId,
           sequence.size() - 1);
      return false;
    }

    if (secondaryBlockId != -1) {
      // Check secondary block ID
      if (secondaryBlockId < 0 ||
          secondaryBlockId >= sequence[blockId].second.size()) {
        warn(
            "Invalid secondaryBlockId {} for blockId {} when applying mutation",
            secondaryBlockId, blockId);
        return false;
      }

      // Check nucleotide position for secondary block
      if (nucPos < 0 ||
          nucPos >= sequence[blockId].second[secondaryBlockId].size()) {
        warn("Invalid nucPos {} for blockId {}, secondaryBlockId {} when "
             "applying mutation",
             nucPos, blockId, secondaryBlockId);
        return false;
      }

      // Check gap position for secondary block
      if (nucGapPos >= 0 && nucGapPos >= sequence[blockId]
                                             .second[secondaryBlockId][nucPos]
                                             .second.size()) {
        warn("Invalid nucGapPos {} for blockId {}, secondaryBlockId {}, nucPos "
             "{} when applying mutation",
             nucGapPos, blockId, secondaryBlockId, nucPos);
        return false;
      }
    } else {
      // Check nucleotide position for primary block
      if (nucPos < 0 || nucPos >= sequence[blockId].first.size()) {
        warn("Invalid nucPos {} for blockId {} when applying mutation", nucPos,
             blockId);
        return false;
      }

      // Check gap position for primary block
      if (nucGapPos >= 0 &&
          nucGapPos >= sequence[blockId].first[nucPos].second.size()) {
        warn("Invalid nucGapPos {} for blockId {}, nucPos {} when applying "
             "mutation",
             nucGapPos, blockId, nucPos);
        return false;
      }
    }

    return true;
  };

  // Block Mutations
  std::unordered_map<int32_t, bool> blockMutated;

  // First collect block mutation info
  for (auto mutation : node->blockMutation) {
    int32_t primaryBlockId = mutation.primaryBlockId;
    int32_t secondaryBlockId = mutation.secondaryBlockId;
    bool type = mutation.blockMutInfo;
    bool inversion = mutation.inversion;

    // Validate primary block ID
    if (primaryBlockId < 0 || primaryBlockId >= sequence.size()) {
      warn("Skipping block mutation with invalid primaryBlockId {}",
           primaryBlockId);
      continue;
    }

    // Validate secondary block ID if relevant
    if (secondaryBlockId != -1) {
      if (secondaryBlockId < 0 ||
          secondaryBlockId >= sequence[primaryBlockId].second.size()) {
        warn("Skipping block mutation with invalid secondaryBlockId {} for "
             "blockId {}",
             secondaryBlockId, primaryBlockId);
        continue;
      }
    }

    bool startingStrand = blockStrand[primaryBlockId].first;
    bool oldMut;
    bool oldStrand;
    // Process block mutation and update state
    if (type == 1) {
      // insertion
      if (secondaryBlockId != -1) {
        oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
        oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
        blockExists[primaryBlockId].second[secondaryBlockId] = true;
        blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
      } else {
        oldStrand = blockStrand[primaryBlockId].first;
        oldMut = blockExists[primaryBlockId].first;
        blockExists[primaryBlockId].first = true;
        blockStrand[primaryBlockId].first = !inversion;
      }
      blockMutationInfo.emplace_back(
          std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId,
                          oldMut, oldStrand, true, !inversion));
    } else {
      if (inversion) {
        // This means that this is not a deletion, but instead an inversion
        if (secondaryBlockId != -1) {
          oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
          oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
          blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
        } else {
          oldStrand = blockStrand[primaryBlockId].first;
          oldMut = blockExists[primaryBlockId].first;
          blockStrand[primaryBlockId].first = !oldStrand;
        }
        if (oldMut != true) {
          warn("There was a problem in PanMAT generation. Please Report.");
        }
        blockMutationInfo.emplace_back(
            std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId,
                            oldMut, oldStrand, oldMut, !oldStrand));
      } else {
        // Actually a deletion
        if (secondaryBlockId != -1) {
          oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
          oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
          blockExists[primaryBlockId].second[secondaryBlockId] = false;
          blockStrand[primaryBlockId].second[secondaryBlockId] = true;
        } else {
          oldStrand = blockStrand[primaryBlockId].first;
          oldMut = blockExists[primaryBlockId].first;
          blockExists[primaryBlockId].first = false;
          blockStrand[primaryBlockId].first = true;
        }
        blockMutationInfo.emplace_back(
            std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId,
                            oldMut, oldStrand, false, true));
      }
    }

    // Track strand changes
    if (startingStrand != blockStrand[primaryBlockId].first) {
      if (blockStrand[primaryBlockId].first) {
        inverseBlockIds.erase(primaryBlockId);
        inverseBlockIdsBacktrack.emplace_back(
            std::make_pair(false, primaryBlockId));
      } else {
        inverseBlockIds.insert(primaryBlockId);
        inverseBlockIdsBacktrack.emplace_back(
            std::make_pair(true, primaryBlockId));
      }

      // If this is an inversion (rather than simply turning on a block), invert
      // the gap map
      if (oldMut) { // Only invert if the block was already active
        // Get the range of the block for gap map inversion
        auto range = manager.getBlockRange(primaryBlockId);
        if (range.first >= 0 && range.second >= range.first) {
          debug_msg("Inverting gap map for block {} with range [{},{}] due to "
                "strand change",
                primaryBlockId, range.first, range.second);

          // Create backtrack vectors for the inversion operation
          std::vector<gap_map::GapUpdate> inversionBacktracks;
          std::vector<gap_map::GapUpdate> inversionUpdates;

          // Invert the gap map for this block
          manager.invertGapMap(range, inversionBacktracks, inversionUpdates);
          
          // Validate the gap map after inversion
          auto &gap_map = const_cast<coordinates::GapMap&>(manager.getGapMap());
          if (!visualization::validateGapMap(gap_map, manager.size())) {
            msg("Gap map validation failed after inverting range [{}:{}]",
                     range.first, range.second);
            
            // Apply gap map fixes if needed
            if (gap_map::preventOverlaps(gap_map)) {
              msg("Fixed overlapping gaps after inversion in applyMutations");
            }
          }

          // Add backtracking entries to main lists
          gapRunBacktracks.insert(gapRunBacktracks.end(),
                                  inversionBacktracks.begin(),
                                  inversionBacktracks.end());
          gapMapUpdates.insert(gapMapUpdates.end(), inversionUpdates.begin(),
                               inversionUpdates.end());
        }
      }
    }

    // For indexing, track affected block ranges
    if (!isPlacement) {
      blockMutated[primaryBlockId] = true; // Mark block as needing recomp range
    }
  }

  // Reinitialize block coordinates after all block mutations
  for (const auto &key : blockMutated) {
    manager.reinitializeBlockCoordinates(key.first);
  }

  if (!isPlacement) {
    // Process block mutations for indexing
    for (const auto &key : blockMutated) {
      try {
        const tupleCoord_t *start =
            manager.getFirstCoordinateInBlock(key.first);
        const tupleCoord_t *stop = manager.getLastCoordinateInBlock(key.first);

        if (start && stop) {
          recompRanges.emplace_back(start, stop);
        }
      } catch (const std::exception &e) {
        err("Failed to create block range: {}", e.what());
      }
    }
  }

  // Process nucleotide mutations
  std::vector<int64_t> mutationLengths;
  // Store minimal info needed for recomp range creation: blockId, start pos,
  // gap pos, length
  struct RecompInfo {
    int32_t blockId;
    int32_t startPos;
    int32_t gapPos;
    int32_t length;
    bool isSNP;
    int32_t nucPos;     // Added missing field
    int32_t nucGapPos;  // Added missing field
    int32_t len;        // Added missing field (alias for length)
  };
  std::vector<RecompInfo> recompMutations;

  // First collect all nucleotide mutation info
  for (int64_t i = 0; i < node->nucMutation.size(); i++) {
    try {
      int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
      int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;
      int32_t nucPosition = node->nucMutation[i].nucPosition;
      int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
      uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
      int len = (type < 3) ? (node->nucMutation[i].mutInfo >> 4) : 0;

      // Store minimal mutation info for later recomp range creation
      if (!isPlacement && blockExists[primaryBlockId].first &&
          !blockMutated[primaryBlockId]) {
        recompMutations.push_back({
            primaryBlockId, nucPosition, nucGapPosition, len,
            type >= 3 // isSNP
        });
      }

      // Apply the mutation - with bounds checking
      if (type < 3) {
        if (type == panmanUtils::NucMutationType::NS) {
          // Substitution
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition, nucGapPosition + j)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .second[secondaryBlockId][nucPosition]
                                  .second[nucGapPosition + j];
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition,
                    nucGapPosition + j, oldVal, newVal));
              }
            } else {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition + j, nucGapPosition)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .second[secondaryBlockId][nucPosition + j]
                                  .first;
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, newVal));
              }
            }
          } else {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition, nucGapPosition + j)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .first[nucPosition]
                                  .second[nucGapPosition + j];
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition,
                    nucGapPosition + j, oldVal, newVal));
              }
            } else {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition + j, nucGapPosition)) {
                  continue;
                }

                char oldVal =
                    sequence[primaryBlockId].first[nucPosition + j].first;
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId].first[nucPosition + j].first = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, newVal));
              }
            }
          }
        } else if (type == panmanUtils::NucMutationType::NI) {
          // Insertion
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition, nucGapPosition + j)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .second[secondaryBlockId][nucPosition]
                                  .second[nucGapPosition + j];
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition,
                    nucGapPosition + j, oldVal, newVal));
              }
            } else {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition + j, nucGapPosition)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .second[secondaryBlockId][nucPosition + j]
                                  .first;
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, newVal));
              }
            }
          } else {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition, nucGapPosition + j)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .first[nucPosition]
                                  .second[nucGapPosition + j];
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, newVal));
              }
            } else {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition + j, nucGapPosition)) {
                  continue;
                }

                char oldVal =
                    sequence[primaryBlockId].first[nucPosition + j].first;
                char newVal = panmanUtils::getNucleotideFromCode(
                    ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId].first[nucPosition + j].first = newVal;
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, newVal));
              }
            }
          }
        } else if (type == panmanUtils::NucMutationType::ND) {
          // Deletion
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition, nucGapPosition + j)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .second[secondaryBlockId][nucPosition]
                                  .second[nucGapPosition + j];
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = '-';
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition,
                    nucGapPosition + j, oldVal, '-'));
              }
            } else {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition + j, nucGapPosition)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .second[secondaryBlockId][nucPosition + j]
                                  .first;
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = '-';
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, '-'));
              }
            }
          } else {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition, nucGapPosition + j)) {
                  continue;
                }

                char oldVal = sequence[primaryBlockId]
                                  .first[nucPosition]
                                  .second[nucGapPosition + j];
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = '-';
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition,
                    nucGapPosition + j, oldVal, '-'));
              }
            } else {
              for (int j = 0; j < len; j++) {
                if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                     nucPosition + j, nucGapPosition)) {
                  continue;
                }

                char oldVal =
                    sequence[primaryBlockId].first[nucPosition + j].first;
                sequence[primaryBlockId].first[nucPosition + j].first = '-';
                mutationInfo.emplace_back(std::make_tuple(
                    primaryBlockId, secondaryBlockId, nucPosition + j,
                    nucGapPosition, oldVal, '-'));
              }
            }
          }
        }
      } else {
        len = 0;
        if (type == panmanUtils::NucMutationType::NSNPS) {
          // SNP Substitution
          char newVal = panmanUtils::getNucleotideFromCode(
              ((node->nucMutation[i].nucs) >> 20) & 0xF);
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .second[secondaryBlockId][nucPosition]
                                .second[nucGapPosition];
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            } else {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .second[secondaryBlockId][nucPosition]
                                .first;
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            }
          } else {
            if (nucGapPosition != -1) {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .first[nucPosition]
                                .second[nucGapPosition];
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            } else {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId].first[nucPosition].first;
              sequence[primaryBlockId].first[nucPosition].first = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            }
          }
        } else if (type == panmanUtils::NucMutationType::NSNPI) {
          // SNP Insertion
          char newVal = panmanUtils::getNucleotideFromCode(
              ((node->nucMutation[i].nucs) >> 20) & 0xF);
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .second[secondaryBlockId][nucPosition]
                                .second[nucGapPosition];
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            } else {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .second[secondaryBlockId][nucPosition]
                                .first;
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            }
          } else {
            if (nucGapPosition != -1) {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .first[nucPosition]
                                .second[nucGapPosition];
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            } else {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId].first[nucPosition].first;
              sequence[primaryBlockId].first[nucPosition].first = newVal;
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, newVal));
            }
          }
        } else if (type == panmanUtils::NucMutationType::NSNPD) {
          // SNP Deletion
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .second[secondaryBlockId][nucPosition]
                                .second[nucGapPosition];
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = '-';
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, '-'));
            } else {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .second[secondaryBlockId][nucPosition]
                                .first;
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = '-';
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, '-'));
            }
          } else {
            if (nucGapPosition != -1) {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId]
                                .first[nucPosition]
                                .second[nucGapPosition];
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = '-';
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, '-'));
            } else {
              if (!isValidMutation(primaryBlockId, secondaryBlockId,
                                   nucPosition, nucGapPosition)) {
                continue;
              }

              char oldVal = sequence[primaryBlockId].first[nucPosition].first;
              sequence[primaryBlockId].first[nucPosition].first = '-';
              mutationInfo.emplace_back(
                  std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition,
                                  nucGapPosition, oldVal, '-'));
            }
          }
        }
      }
    } catch (const std::exception &e) {
      err("Exception in applyMutations for nucleotide mutation at index {}: {}",
          i, e.what());
      continue; // Skip this mutation and continue with the rest
    }
  }

  // Create recomp ranges after all mutations are applied
  if (!isPlacement && !recompMutations.empty()) {
    std::sort(recompMutations.begin(), recompMutations.end(), 
      [](const RecompInfo& a, const RecompInfo& b) {
        if (a.blockId != b.blockId) return a.blockId < b.blockId;
        if (a.startPos != b.startPos) return a.startPos < b.startPos;
        // Special handling for gapPos: -1 should come after positive values
        if (a.gapPos != b.gapPos) {
          if (a.gapPos == -1) return false; // a.gapPos = -1 comes after b.gapPos
          if (b.gapPos == -1) return true;  // b.gapPos = -1 comes after a.gapPos
          return a.gapPos < b.gapPos;       // Normal comparison for positive values
        }
        return a.nucPos < b.nucPos;
      });
    for (size_t i = 0; i < recompMutations.size(); i++) {
      auto &mut = recompMutations[i];
      int32_t blockId = mut.blockId;
      int32_t nucPos = mut.nucPos;
      int32_t nucGapPos = mut.nucGapPos;

      auto range = manager.getBlockRange(blockId);
      int64_t start_scalar = range.first;
      int64_t end_scalar = range.second;

      if (start_scalar < 0 || end_scalar < 0) {
        warn("Invalid block range for blockId {}: [{}, {}]", blockId,
             start_scalar, end_scalar);
        continue;
      }

      // Create range for this mutation
      // For SNPs, we only need a small range, otherwise use k + 1 positions
      int len = mut.len;
      int recompSize = mut.isSNP ? DEFAULT_SHORT_RECOMP_SIZE : (k + 1);

      int64_t start_pos = -1;
      const auto *firstCoord = manager.getFirstCoordinateInBlock(blockId);
      if (firstCoord) {
        start_pos = firstCoord->scalar;
      }

      const tupleCoord_t *startCoord = nullptr;
      const tupleCoord_t *endCoord = nullptr;

      // Get the coordinates for the recomputation range
      traverser.reset(manager.getLeftmostCoordinate());
      int64_t counter = 0;
      int64_t remaining = start_scalar + nucPos;
      while (!traverser.isDone() && counter < remaining) {
        traverser.next();
        counter++;
      }
      if (!traverser.isDone()) {
        startCoord = traverser.getCurrent();
        counter = 0;
        int maxRemaining = std::min<int64_t>(recompSize, end_scalar - remaining);
        while (!traverser.isDone() && counter < maxRemaining) {
          traverser.next();
          counter++;
        }
        endCoord = traverser.getCurrent();
      }

      if (startCoord && endCoord) {
        recompRanges.push_back(CoordRange(startCoord, endCoord));
      } else {
        warn("Failed to create coordinates for recomp range");
      }
    }
  }

  // Detect and handle continuous runs of off blocks
  // This optimization is particularly useful for large genomes with many blocks
  std::vector<std::pair<int32_t, int32_t>> offBlockRuns;
  
  // Find runs of consecutive off blocks
  int32_t runStart = -1;
  for (int32_t blockId = 0; blockId < blockExists.size(); ++blockId) {
    bool isOff = !blockExists[blockId].first;
    
    if (isOff) {
      // Start a new run or continue the current one
      if (runStart == -1) {
        runStart = blockId;
      }
    } else {
      // End the current run if there was one
      if (runStart != -1) {
        int32_t runEnd = blockId - 1;
        if (runEnd - runStart >= MIN_BLOCK_RUN_SIZE) {
          // Only store runs of at least MIN_BLOCK_RUN_SIZE blocks
          offBlockRuns.emplace_back(runStart, runEnd);
        }
        runStart = -1;
      }
    }
  }
  
  // Don't forget the last run if it extends to the end
  if (runStart != -1) {
    int32_t runEnd = blockExists.size() - 1;
    if (runEnd - runStart >= MIN_BLOCK_RUN_SIZE) {
      offBlockRuns.emplace_back(runStart, runEnd);
    }
  }
  
  // Process each run of off blocks
  for (const auto &[startBlockId, endBlockId] : offBlockRuns) {
    // Create temporary vectors of the correct type
    std::vector<std::pair<int64_t, int64_t>> tempGapRunBacktracks;
    std::vector<std::pair<size_t, std::pair<bool, int64_t>>> tempGapMapUpdates;
    
    processBlockRun(
      recompRanges,
      gapRunUpdates, 
      tempGapRunBacktracks,  // Use the temporary vector of the correct type
      tempGapMapUpdates,     // Use the temporary vector of the correct type
      traverser, 
      startBlockId, 
      endBlockId, 
      isPlacement);
    
    if (isDebugEnabled()) {
      msg("Processed run of {} off blocks from {} to {}", 
          endBlockId - startBlockId + 1, startBlockId, endBlockId);
    }
  }
  
  // Now also detect runs of blocks with the same orientation
  std::vector<std::pair<int32_t, int32_t>> orientationRuns;
  
  // Find runs of consecutive blocks with the same orientation (forward or inverted)
  runStart = -1;
  bool runOrientation = true; // true = forward, false = inverted
  
  for (int32_t blockId = 0; blockId < blockExists.size(); ++blockId) {
    // Skip off blocks
    if (!blockExists[blockId].first) {
      // End any current run
      if (runStart != -1) {
        int32_t runEnd = blockId - 1;
        if (runEnd - runStart >= MIN_BLOCK_RUN_SIZE) {
          orientationRuns.emplace_back(runStart, runEnd);
        }
        runStart = -1;
      }
      continue;
    }
    
    bool isForward = blockStrand[blockId].first;
    
    if (runStart == -1) {
      // Start a new run
      runStart = blockId;
      runOrientation = isForward;
    } else if (isForward != runOrientation) {
      // End the current run due to orientation change
      int32_t runEnd = blockId - 1;
      if (runEnd - runStart >= MIN_BLOCK_RUN_SIZE) {
        orientationRuns.emplace_back(runStart, runEnd);
      }
      // Start a new run
      runStart = blockId;
      runOrientation = isForward;
    }
  }
  
  // Don't forget the last run
  if (runStart != -1) {
    int32_t runEnd = blockExists.size() - 1;
    if (runEnd - runStart >= MIN_BLOCK_RUN_SIZE) {
      orientationRuns.emplace_back(runStart, runEnd);
    }
  }
  
  // Process each run of blocks with the same orientation
  for (const auto &[startBlockId, endBlockId] : orientationRuns) {
    bool isInverted = !blockStrand[startBlockId].first;
    if (isDebugEnabled()) {
      msg("Processing run of {} {} blocks from {} to {}", 
          endBlockId - startBlockId + 1, 
          isInverted ? "inverted" : "forward", 
          startBlockId, endBlockId);
    }
    
    // Create temporary vectors of the correct type
    std::vector<std::pair<int64_t, int64_t>> tempGapRunBacktracks;
    std::vector<std::pair<size_t, std::pair<bool, int64_t>>> tempGapMapUpdates;
    
    processBlockRun(
      recompRanges,
      gapRunUpdates, 
      tempGapRunBacktracks,  // Use the temporary vector of the correct type
      tempGapMapUpdates,     // Use the temporary vector of the correct type
      traverser, 
      startBlockId, 
      endBlockId, 
      isPlacement);
  }

  if (!isPlacement) {
    // Process block mutations for indexing
    for (const auto &key : blockMutated) {
      try {
        const tupleCoord_t *start =
            manager.getFirstCoordinateInBlock(key.first);
        const tupleCoord_t *stop = manager.getLastCoordinateInBlock(key.first);

        if (start && stop) {
          recompRanges.emplace_back(start, stop);
        }
      } catch (const std::exception &e) {
        err("Failed to create block range: {}", e.what());
      }
    }
  }

  // Update gap map with current mutations
  traverser.getCoordManager().updateGapMap(
    gapRunUpdates, 
    gapRunBacktracks, 
    gapMapUpdates);
}

void undoMutations(
    Tree *T, const Node *node, const blockMutationInfo_t &blockMutationInfo,
    const mutationInfo_t &mutationInfo,
    coordinates::CoordinateTraverser &traverser,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks) {
  auto &sequence = traverser.getCoordManager().getSequence();
  auto &blockExists = traverser.getCoordManager().getBlockExists();
  auto &blockStrand = traverser.getCoordManager().getBlockStrand();

  // Counters for tracking validity issues
  int invalid_blockId = 0;
  int invalid_secondaryBlockId = 0;
  int invalid_nucPos = 0;
  int invalid_nucGapPos = 0;

  // For debugging
  // msg("Undoing {} block mutations and {} nucleotide mutations for node {}",
  //     blockMutationInfo.size(), mutationInfo.size(), node->identifier);


  // STEP 2: Restore the gap map to its previous state using backtracking
  // information
  // msg("Restoring gap map state from {} backtrack entries",
  //     gapRunBacktracks.size());

  // Enhanced restoration - apply backtracking entries more safely
  try {
    // Backup current gap map in case restoration fails
    auto originalGapMap = traverser.getCoordManager().getGapMap();

    // Try to restore the gap map
    traverser.getCoordManager().restoreGapMapFromBacktracks(gapRunBacktracks);

    // Verify the restored gap map is valid
    auto &gap_map = traverser.getCoordManager().getGapMap();

    // Validate and fix the gap map if needed after restoration
    bool is_valid = visualization::validateGapMap(gap_map, traverser.getCoordManager().size());
    if (!is_valid) {
      msg("Gap map validation failed after restoration - attempting to fix overlaps");
      if (gap_map::preventOverlaps(const_cast<coordinates::GapMap&>(gap_map))) {
        msg("Fixed overlapping gaps in gap map after restoration by merging");
      }
    }
  } catch (const std::exception &e) {
    err("Exception during gap map restoration: {}. Continuing with current "
        "state.",
        e.what());
  }

  // STEP 3: Take a snapshot of the current state for validation
  // This is useful for detecting state changes during the undo process
  auto sequenceSnapshot = sequence;
  auto blockExistsSnapshot = blockExists;
  auto blockStrandSnapshot = blockStrand;

  // STEP 4: First pass - Identify blocks with strand changes
  std::vector<int32_t> strandChangedBlocks;
  for (auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend();
       it++) {
    auto mutation = *it;
    try {
      int32_t blockId = std::get<0>(mutation);
      int32_t secondaryBlockId = std::get<1>(mutation);
      bool oldStrand = std::get<3>(mutation);

      // If primary block's strand is changing
      if (secondaryBlockId == -1 && blockStrand[blockId].first != oldStrand) {
        // Record for special handling
        strandChangedBlocks.push_back(blockId);
      }
      // If secondary block's strand is changing
      else if (secondaryBlockId != -1 && blockId < blockStrand.size() &&
               secondaryBlockId < blockStrand[blockId].second.size() &&
               blockStrand[blockId].second[secondaryBlockId] != oldStrand) {
        // Secondary blocks with strand changes need special handling too
        msg("Secondary block strand change detected: blockId={}, "
            "secondaryBlockId={}",
            blockId, secondaryBlockId);
        // Currently we don't handle secondary blocks specially - this would
        // need implementation
      }
    } catch (const std::exception &e) {
      warn("Exception identifying strand changes: {}. Skipping.", e.what());
    }
  }

  if (!strandChangedBlocks.empty()) {
    // msg("Detected {} blocks with strand changes that need special handling",
    //     strandChangedBlocks.size());
  }

  // STEP 5: Apply the block mutations in reverse order
  // msg("Undoing block mutations");
  for (auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend();
       it++) {
    auto mutation = *it;
    try {
      int32_t blockId = std::get<0>(mutation);
      int32_t secondaryBlockId = std::get<1>(mutation);
      bool oldExists = std::get<2>(mutation);
      bool oldStrand = std::get<3>(mutation);

      // Validate blockId
      if (blockId < 0 || blockId >= blockExists.size()) {
        warn("Invalid blockId {} during undo (max: {}). Skipping.", blockId,
             blockExists.size() - 1);
        invalid_blockId++;
        continue;
      }

      // Process secondary block mutations if needed
      if (secondaryBlockId != -1) {
        // Validate secondaryBlockId
        if (secondaryBlockId < 0 ||
            secondaryBlockId >= blockExists[blockId].second.size()) {
          warn("Invalid secondaryBlockId {} for blockId {} during undo. "
               "Skipping.",
               secondaryBlockId, blockId);
          invalid_secondaryBlockId++;
          continue;
        }

        // Restore old block existence state
        blockExists[blockId].second[secondaryBlockId] = oldExists;

        // Restore old strand state if different
        if (blockStrand[blockId].second[secondaryBlockId] != oldStrand) {
          blockStrand[blockId].second[secondaryBlockId] = oldStrand;
        }
      } else {
        // Restore old block existence state
        blockExists[blockId].first = oldExists;

        // Restore old strand state if different
        if (blockStrand[blockId].first != oldStrand) {
          blockStrand[blockId].first = oldStrand;
        }
      }

      // Crucial step: After each block mutation undo, reinitialize the block
      // coordinates This ensures correct traversal for subsequent operations
      try {
        traverser.getCoordManager().reinitializeBlockCoordinates(blockId);
        debug_msg("Successfully reinitialized block {}", blockId);
      } catch (const std::exception &e) {
        warn("Error reinitializing block {}: {}", blockId, e.what());
      }

    } catch (const std::exception &e) {
      warn("Exception in undoMutations (block mutation): {}. Skipping.",
           e.what());
      continue; // Skip this mutation and continue processing
    }
  }

  // STEP 6: Special handling for blocks that had strand changes
  // For each block that changed strand, re-invert its gap map
  std::vector<gap_map::GapUpdate> strandChangeBacktracks;
  std::vector<gap_map::GapUpdate> strandChangeUpdates;

  for (int32_t blockId : strandChangedBlocks) {
    // Get block range
    auto range = traverser.getCoordManager().getBlockRange(blockId);
    if (range.first < 0 || range.second < range.first) {
      warn("Invalid block range [{}, {}] for block {} with strand change. "
           "Skipping.",
           range.first, range.second, blockId);
      continue;
    }

    msg("Re-inverting gap map for block {} with range [{}, {}] after strand "
        "change",
        blockId, range.first, range.second);

    // Re-invert gap runs in this range
    traverser.getCoordManager().invertGapMap(range, strandChangeBacktracks,
                                             strandChangeUpdates);
  }

  if (!strandChangeBacktracks.empty()) {
    // msg("Created {} backtrack entries while re-inverting gaps for strand "
    //     "changes",
    //     strandChangeBacktracks.size());
  }

  // STEP 7: Undo nucleotide mutations in reverse order from application
  // msg("Undoing nucleotide mutations");
  for (auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++) {
    auto mutation = *it;
    try {
      // Extract mutation components
      int32_t blockId = std::get<0>(mutation);
      int32_t secondaryBlockId = std::get<1>(mutation);
      int32_t nucPos = std::get<2>(mutation);
      int32_t nucGapPos = std::get<3>(mutation);
      char oldValue = std::get<4>(mutation);
      char newValue = std::get<5>(mutation);

      // Skip if invalid blockId
      if (blockId < 0 || blockId >= sequence.size()) {
        invalid_blockId++;
        if (invalid_blockId <= 5) { // Limit excessive logging
          // warn("Invalid blockId {} when undoing mutation (max: {})", blockId,
          //     sequence.size() - 1);
        }
        continue;
      }

      // Skip if block doesn't exist
      if (!blockExists[blockId].first) {
        debug_msg("Skipping nucleotide mutation for dead block {}", blockId);
        continue;
      }

      // Separate handling for secondary blocks vs primary blocks
      if (secondaryBlockId != -1) {
        // Check if secondaryBlockId is valid
        if (secondaryBlockId < 0 ||
            secondaryBlockId >= sequence[blockId].second.size()) {
          invalid_secondaryBlockId++;
          if (invalid_secondaryBlockId <= 5) {
            // warn("Invalid secondaryBlockId {} for blockId {}", secondaryBlockId,
            //     blockId);
          }
          continue;
        }

        // Handle gap positions vs non-gap positions
        if (nucGapPos != -1) {
          // Validate nucPos
          if (nucPos < 0 ||
              nucPos >= sequence[blockId].second[secondaryBlockId].size()) {
            invalid_nucPos++;
            if (invalid_nucPos <= 5) {
              // warn("Invalid nucPos {} for blockId {}, secondaryBlockId {}",
              //     nucPos, blockId, secondaryBlockId);
            }
            continue;
          }

          // Validate nucGapPos
          if (nucGapPos < 0 ||
              nucGapPos >= sequence[blockId]
                               .second[secondaryBlockId][nucPos]
                               .second.size()) {
            invalid_nucGapPos++;
            if (invalid_nucGapPos <= 5) {
              // warn("Invalid nucGapPos {} for blockId {}, secondaryBlockId {}, "
              //     "nucPos {}",
              //     nucGapPos, blockId, secondaryBlockId, nucPos);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue = sequence[blockId]
                                  .second[secondaryBlockId][nucPos]
                                  .second[nucGapPos];
          if (currentValue != newValue && newValue != '\0') {
            // warn("Value mismatch during undo at [{},{},{},{}]: expected '{}', "
            //      "found '{}'",
            //      blockId, secondaryBlockId, nucPos, nucGapPos, newValue,
            //      currentValue);
          }
          sequence[blockId].second[secondaryBlockId][nucPos].second[nucGapPos] =
              oldValue;

        } else {
          // Validate nucPos
          if (nucPos < 0 ||
              nucPos >= sequence[blockId].second[secondaryBlockId].size()) {
            invalid_nucPos++;
            if (invalid_nucPos <= 5) {
              // warn("Invalid nucPos {} for blockId {}, secondaryBlockId {}",
              //     nucPos, blockId, secondaryBlockId);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue =
              sequence[blockId].second[secondaryBlockId][nucPos].first;
          if (currentValue != newValue && newValue != '\0') {
            // warn("Value mismatch during undo at [{},{},{}]: expected '{}', "
            //      "found '{}'",
            //      blockId, secondaryBlockId, nucPos, newValue, currentValue);
          }
          sequence[blockId].second[secondaryBlockId][nucPos].first = oldValue;
        }
      } else {
        // Handle primary block mutations
        if (nucGapPos != -1) {
          // Validate nucPos
          if (nucPos < 0 || nucPos >= sequence[blockId].first.size()) {
            invalid_nucPos++;
            if (invalid_nucPos <= 5) {
              // warn("Invalid nucPos {} for blockId {}", nucPos, blockId);
            }
            continue;
          }

          // Validate nucGapPos
          if (nucGapPos < 0 ||
              nucGapPos >= sequence[blockId].first[nucPos].second.size()) {
            invalid_nucGapPos++;
            if (invalid_nucGapPos <= 5) {
              // warn("Invalid nucGapPos {} for blockId {}, nucPos {}",
              //     nucGapPos, blockId, nucPos);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue = sequence[blockId].first[nucPos].second[nucGapPos];
          if (currentValue != newValue && newValue != '\0') {
            // warn("Value mismatch during undo at [{},{},{}]: expected '{}', "
            //      "found '{}'",
            //      blockId, nucPos, nucGapPos, newValue, currentValue);
          }
          sequence[blockId].first[nucPos].second[nucGapPos] = oldValue;

        } else {
          // Validate nucPos
          if (nucPos < 0 || nucPos >= sequence[blockId].first.size()) {
            invalid_nucPos++;
            if (invalid_nucPos <= 5) {
              // warn("Invalid nucPos {} for blockId {}", nucPos, blockId);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue = sequence[blockId].first[nucPos].first;
          if (currentValue != newValue && newValue != '\0') {
            // warn("Value mismatch during undo at [{},{}]: expected '{}', found "
            //      "'{}'",
            //      blockId, nucPos, newValue, currentValue);
          }
          sequence[blockId].first[nucPos].first = oldValue;
        }
      }
    } catch (const std::exception &e) {
      warn("Exception in undoMutations: {}. Continuing with next mutation.",
           e.what());
      continue;
    }
  }

  // Log summary of skipped invalid positions
  if (invalid_blockId > 0 || invalid_secondaryBlockId > 0 ||
      invalid_nucPos > 0 || invalid_nucGapPos > 0) {
    // warn("Skipped invalid positions during undo - "
    //      "blockId: {}, secondaryBlockId: {}, nucPos: {}, nucGapPos: {}",
    //      invalid_blockId, invalid_secondaryBlockId, invalid_nucPos,
    //      invalid_nucGapPos);
  } else {
    // msg("All nucleotide mutations successfully undone");
  }


  // STEP 9: One final pass to ensure block coordinates are properly linked
  try {
    // msg("Performing final coordinate system cleanup");
    bool prev_linked = false;
    int32_t prev_live_blockId = -1;

    // Link all blocks in sequence
    for (int32_t blockId = 0; blockId < blockExists.size(); blockId++) {
      if (blockExists[blockId].first) {
        // This is a live block
        traverser.getCoordManager().reinitializeBlockCoordinates(blockId);

        // Connect this block to previous live block if it exists
        if (prev_linked && prev_live_blockId >= 0) {
          const tupleCoord_t *prev_last =
              traverser.getCoordManager().getLastCoordinateInBlock(
                  prev_live_blockId);
          const tupleCoord_t *curr_first =
              traverser.getCoordManager().getFirstCoordinateInBlock(blockId);

          if (prev_last && curr_first) {
            // Link previous block's last coordinate to this block's first
            // We cast away const because we need to update the pointers
            tupleCoord_t *mutable_prev_last =
                const_cast<tupleCoord_t *>(prev_last);
            tupleCoord_t *mutable_curr_first =
                const_cast<tupleCoord_t *>(curr_first);

            mutable_prev_last->next = mutable_curr_first;
            mutable_curr_first->prev = mutable_prev_last;
          }
        }

        prev_linked = true;
        prev_live_blockId = blockId;
      }
    }
  } catch (const std::exception &e) {
    // warn("Exception during final coordinate cleanup: {}", e.what());
  }

  // msg("Mutation undo complete for node {}", node->identifier);
}


void fillMutationMatricesFromTree_test(
    mutationMatrices &mutMat, Tree *T, const std::string &path) {
  return;
}

void fillMutationMatricesFromFile(mutationMatrices &mutMat,
                                                       std::ifstream &inf) {

  std::string line;
  int idx = 0;
  while (getline(inf, line)) {
    std::vector<double> probs;
    std::vector<std::string> fields;
    stringSplit(line, ' ', fields);
    for (const auto &f : fields) {
      probs.push_back(std::stod(f));
    }
    if (probs.size() == 0) {
      break;
    }

    if (idx < 4) {
      if (probs.size() != 4) {
        throw std::invalid_argument(
            "Received invalid mutation matrix (.mm) file");
      }

      std::vector<double> probs;
      for (const auto &f : fields) {
        probs.push_back(std::stod(f));
      }
      mutMat.submat[idx] = std::move(probs);
    } else if (idx == 4) {
      if (probs.size() < 1) {
        throw std::invalid_argument(
            "Received invalid mutation matrix (.mm) file");
      }

      for (const auto &f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.insmat[size] = prob;
      }
    } else if (idx == 5) {
      if (probs.size() < 1) {
        throw std::invalid_argument(
            "Received invalid mutation matrix (.mm) file");
      }

      for (const auto &f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.delmat[size] = prob;
      }
    }
    idx++;
  }

  if (idx != 6) {
    throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
  }
  mutMat.filled = true;
}

void buildMutationMatricesHelper_test(
    seed_annotated_tree::mutationMatrices &mutMat, panmanUtils::Tree *T, panmanUtils::Node *node,
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
    std::unordered_map<int64_t, int64_t> &delCount) {

  return;
}

void setupIndexing(sequence_t &sequence, blockExists_t &blockExists,
                   blockStrand_t &blockStrand, const Tree *T) {
  // Call the common setupSequence function
  setupSequence(T, sequence, blockExists, blockStrand);
}

void setupPlacement(
    std::vector<std::optional<seeding::onSeedsHash>> &onSeedsHash,
    sequence_t &sequence, blockExists_t &blockExists,
    blockStrand_t &blockStrand, const Tree *T) {
  // Call the common setupSequence function
  setupSequence(T, sequence, blockExists, blockStrand);
  
  // Initialize onSeedsHash vector to the appropriate size
  // This will be populated later during the placement process
  if (T && !sequence.empty()) {
    int64_t totalCoords = 0;
    
    try {
      // Correctly count all coordinates including gap positions
      // Ignore secondary block IDs since they're deprecated
      for (size_t blockId = 0; blockId < sequence.size(); blockId++) {
        // Count primary block coordinates
        for (size_t nucPos = 0; nucPos < sequence[blockId].first.size(); nucPos++) {
          // Count the base coordinate (nucGapPos = -1)
          totalCoords++;
          
          // Count all gap positions for this nucleotide
          totalCoords += sequence[blockId].first[nucPos].second.size();
        }
      }
      
      msg("Initializing onSeedsHash with {} total coordinates", totalCoords);
      
      
      // Resize with error handling
      try {
        onSeedsHash.clear();
        onSeedsHash.resize(totalCoords);
        msg("Successfully resized onSeedsHash to {} elements", totalCoords);
      } catch (const std::bad_alloc& e) {
        err("Memory allocation failed when resizing onSeedsHash to {} elements: {}", totalCoords, e.what());
        throw std::runtime_error(fmt::format("Failed to allocate memory for {} coordinates", totalCoords));
      } catch (const std::exception& e) {
        err("Exception when resizing onSeedsHash: {}", e.what());
        throw;
      }
    } catch (const std::exception& e) {
      err("Exception in setupPlacement while counting coordinates: {}", e.what());
      throw;
    }
  } else {
    err("Warning: Empty sequence or null tree in setupPlacement");
    onSeedsHash.clear();
  }
}

// Implementation of shared functions for mutation application and backtracking
void applyTreeNodeMutations(
  seed_annotated_tree::TraversalNodeState &nodeState,
  coordinates::CoordinateTraverser &traverser,
  panmanUtils::Tree *T, 
  panmanUtils::Node *node,
  bool isPlacement,
  std::unordered_set<int32_t> &inverseBlockIds,
  int k) {
  
  // Initialize mutation state
  nodeState.oldBlockExists = traverser.getCoordManager().getBlockExists();
  nodeState.oldBlockStrand = traverser.getCoordManager().getBlockStrand();
  nodeState.gapRunBacktracks.clear();

  // Apply mutations - core logic shared between indexing and placement
  applyMutations(
    nodeState.blockMutationInfo, 
    nodeState.recompRanges,
    nodeState.nucleotideMutationInfo, 
    nodeState.gapRunUpdates,
    nodeState.gapRunBacktracks, 
    nodeState.gapMapUpdates, 
    T, node, traverser, 
    isPlacement, 
    inverseBlockIds,
    nodeState.inverseBlockIdsBacktrack, 
    k);

  // Update gap map with current mutations
  traverser.getCoordManager().updateGapMap(
    nodeState.gapRunUpdates, 
    nodeState.gapRunBacktracks,
    nodeState.gapMapUpdates);
}

bool validateGapMapAfterMutations(
  coordinates::CoordinateTraverser &traverser,
  panmanUtils::Node *node) {
  
  auto &manager = traverser.getCoordManager();
  // Use the global debug variable instead of defining a local one
  // bool debug = true; // Could make this configurable
  
  bool gap_map_valid = manager.validateGapMap(true);
  if (!gap_map_valid) {
    msg("Gap map needed optimization after updates for node {}", node->identifier);
  }
  
  // Perform comprehensive validation after significant modifications
  bool state_valid = manager.validateFullState();
  if (!state_valid) {
    msg("Warning: Full state validation failed after updating gap map for node {}", 
        node->identifier);
    
    // This would be a serious error, but we'll try to continue
    msg("Attempting to proceed despite validation failure");
  }
  
  return gap_map_valid && state_valid;
}

void backtrackNodeState(
  panmanUtils::Tree *T, 
  panmanUtils::Node *node, 
  seed_annotated_tree::TraversalNodeState &nodeState,
  coordinates::CoordinateTraverser &traverser,
  bool undoSeedChanges,
  std::unordered_map<size_t, int64_t> *currentGenomeSeedCounts,
  std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t>> *seedChanges) {
  
  // Undo mutations using the stored mutation information
  undoMutations(
    T, node, 
    nodeState.blockMutationInfo,
    nodeState.nucleotideMutationInfo, 
    traverser, 
    nodeState.gapRunBacktracks);
  
  // Handle seed changes if applicable (placement-specific)
  if (undoSeedChanges && currentGenomeSeedCounts && seedChanges) {
    for (const auto &seedChange : *seedChanges) {
      int64_t pos = std::get<0>(seedChange);
      size_t oldSeedVal = std::get<1>(seedChange);
      bool oldIsReverse = std::get<2>(seedChange);
      int64_t oldEndPos = std::get<3>(seedChange);
      int64_t oldSeedCount = std::get<4>(seedChange);
      
      // Restore the old seed count
      (*currentGenomeSeedCounts)[oldSeedVal] = oldSeedCount;
    }
  }
}

void processGapMutationsForPlacement(
  const ::capnp::List<GapMutations>::Reader &gapMutationsList,
  PlacementNodeState &nodeState,
  coordinates::CoordinateManager &manager) {
    
  for (const auto &gapMutation : gapMutationsList) {
    // Process each delta in the list
    for (const auto &delta : gapMutation.getDeltas()) {
      // Access the position and value
      int32_t pos = delta.getPos();
      auto maybeValue = delta.getMaybeValue();
      
      if (maybeValue.isValue()) {
        // This is a set or update operation
        auto it = manager.getGapMap().find(pos);
        if (it != manager.getGapMap().end()) {
          // Found existing value at this position - store for backtracking
          nodeState.gapRunBacktracks.emplace_back(
              std::make_pair(true, std::make_pair(pos, it->second)));
        } else {
          // No existing value - store deletion for backtracking
          nodeState.gapRunBacktracks.emplace_back(
              std::make_pair(false, std::make_pair(pos, 0)));
        }
        
        // Apply the update
        int64_t newValue = maybeValue.getValue();
        nodeState.gapRunUpdates.emplace_back(pos, newValue, true);
      } else {
        // This is a deletion operation
        auto it = manager.getGapMap().find(pos);
        if (it != manager.getGapMap().end()) {
          // Found existing value - store for backtracking
          nodeState.gapRunBacktracks.emplace_back(
              std::make_pair(true, std::make_pair(pos, it->second)));
          
          // Apply the deletion
          nodeState.gapRunUpdates.emplace_back(
              std::make_pair(false, std::make_pair(pos, 0)));
        }
      }
    }
  }
}

// Implementation of helper for processing seed mutations in placement
void processSeedMutationsForPlacement(
  PlacementGlobalState &state,
  PlacementNodeState &nodeState,
  coordinates::CoordinateTraverser &traverser,
  const ::capnp::List<SeedMutations>::Reader &seedIndex,
  int dfsIndex,
  std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
  std::unordered_map<size_t, int64_t> &currentGenomeSeedCounts,
  int64_t &hitsInThisGenome) {
  
  auto &manager = traverser.getCoordManager();
  auto currBasePositions = seedIndex[dfsIndex].getBasePositions();
  auto currPerPosMasks = seedIndex[dfsIndex].getPerPosMasks();
  
  // Process the seed mutations
  for (int i = 0; i < currBasePositions.size(); ++i) {
    int64_t pos = currBasePositions[i];
    uint64_t tritMask = currPerPosMasks[i];

    for (int j = 0; j < 32; ++j) {
      uint8_t ternaryNumber = (tritMask >> (j * 2)) & 0x3;
      
      // Fix: Use kmerSize from state instead of state.k
      int k = state.kmerSize; // Use the appropriate field name for k-mer size
      
      if (ternaryNumber == 1) { // on -> off
        if (!state.onSeedsHash[pos].has_value()) {
          continue;
        }

        // Get current values before changing
        auto [oldSeed, oldEndPos, oldIsReverse] =
            state.onSeedsHash[pos].value();

        // Track seed count changes for backtracking
        int64_t oldSeedCount = 0;
        auto it = currentGenomeSeedCounts.find(oldSeed);
        if (it != currentGenomeSeedCounts.end()) {
            oldSeedCount = it->second;
            it->second--;
            
            // If seed count is 0, this seed is no longer present
            if (it->second == 0) {
                // Decrement the counter if the seed exists in the reads
                auto readIt = readSeedCounts.find(oldSeed);
                if (readIt != readSeedCounts.end()) {
                    hitsInThisGenome--;
                }
                
                // Remove seeds with zero count to save memory
                currentGenomeSeedCounts.erase(oldSeed);
            }
        }

        // Record change for backtracking
        nodeState.seedChanges.emplace_back(pos, oldSeed, oldIsReverse, 
                                          oldEndPos, oldSeedCount);

        // Turn off seed
        state.onSeedsHash[pos].reset();
        int blockId = manager.getBlockIdOfScalarCoord(pos);
        state.BlocksToSeeds[blockId].erase(pos);
      }
      else if (ternaryNumber == 2) { // off -> on or update
        bool wasOn = state.onSeedsHash[pos].has_value();
        size_t oldSeed = 0;
        int64_t oldEndPos = 0;
        bool oldIsReverse = false;
        int64_t oldSeedCount = 0;
        
        // Get current values if seed exists
        if (wasOn) {
            auto [seed, endPos, isReverse] = state.onSeedsHash[pos].value();
            oldSeed = seed;
            oldEndPos = endPos;
            oldIsReverse = isReverse;
            
            // Track current count
            auto it = currentGenomeSeedCounts.find(oldSeed);
            if (it != currentGenomeSeedCounts.end()) {
                oldSeedCount = it->second;
                it->second--;
                
                // If seed count is 0, this seed is no longer present
                if (it->second == 0) {
                    // Decrement the counter if the seed exists in the read
                    auto readIt = readSeedCounts.find(oldSeed);
                    if (readIt != readSeedCounts.end()) {
                        hitsInThisGenome--;
                    }
                    
                    // Remove seeds with zero count
                    currentGenomeSeedCounts.erase(oldSeed);
                }
            }
        }

        // Determine new seed at this position
        size_t resultHash = 0;
        int64_t resultEndPos = 0;
        bool resultIsReverse = false;
        
        if (getSeedAt(traverser, manager, resultHash, resultIsReverse, 
                     resultEndPos, pos, nullptr, k)) {
            
            // Record change for backtracking
            nodeState.seedChanges.emplace_back(pos, oldSeed, oldIsReverse, 
                                              oldEndPos, oldSeedCount);
            
            // Update seed state
            state.onSeedsHash[pos] = {resultHash, resultEndPos, resultIsReverse};
            int blockId = manager.getBlockIdOfScalarCoord(pos);
            state.BlocksToSeeds[blockId].insert(pos);
            
            // Update seed count for this genome
            currentGenomeSeedCounts[resultHash]++;
            
            // If this is the first occurrence of this seed in the genome
            // and it exists in the read, increment the hit counter
            if (currentGenomeSeedCounts[resultHash] == 1) {
                auto readIt = readSeedCounts.find(resultHash);
                if (readIt != readSeedCounts.end()) {
                    hitsInThisGenome++;
                }
            }
        }
        else if (wasOn) {
            // Failed to get a new seed, but need to remove the old one
            state.onSeedsHash[pos].reset();
            int blockId = manager.getBlockIdOfScalarCoord(pos);
            state.BlocksToSeeds[blockId].erase(pos);
        }
      }
    }
  }
}

// Shared utility function for initializing gap map from sequence traversal
void initializeGapMapFromSequence(
    coordinates::CoordinateTraverser &traverser,
    coordinates::CoordinateManager &manager) {
    
  // Reset traverser to the beginning
  traverser.reset(manager.getLeftmostCoordinate());

  int64_t current_gap_start = -1;
  int64_t current_run_length = 0;
  std::map<int64_t, int64_t> gap_runs;

  // Scan through the sequence to identify all gap runs
  while (!traverser.isDone()) {
    const tupleCoord_t *coord = traverser.getCurrent();
    int64_t scalar = coord->scalar;
    bool is_gap = (traverser.getChar() == '-');
    
    // Check if we're in an inverted block - important for consistent gap handling
    int32_t blockId = coord->blockId;
    bool isInverted = false;
    
    if (blockId >= 0 && blockId < manager.getBlockStrand().size()) {
      isInverted = !manager.getBlockStrand()[blockId].first;
    }

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
        // Always store gap runs with start < end regardless of block orientation
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

  // Filter out suspicious gaps and add valid ones to the map
  int gap_runs_filtered = 0;
  
  // Clear the coordinate manager's gap map
  manager.clearGapMap();

  // Add the gap runs to the map
  for (const auto &[start, length] : gap_runs) {
    // Skip any suspicious full-sequence gaps
    if (isSuspiciousGap(start, length, manager.size())) {
      warn("Skipping suspicious full-sequence gap: position={}, length={}", start, length);
      gap_runs_filtered++;
      continue;
    }
    
    // Validate the gap boundaries to ensure they're within sequence bounds
    if (start < 0 || start >= manager.size() || length <= 0) {
      warn("Skipping invalid gap: position={}, length={}", start, length);
      gap_runs_filtered++;
      continue;
    }
    
    // Handle case where the gap would extend beyond sequence bounds
    if (start + length > manager.size()) {
      int64_t adjusted_length = manager.size() - start;
      warn("Adjusting gap length from {} to {} at position {}", length, adjusted_length, start);
      
      if (adjusted_length <= 0) {
        gap_runs_filtered++;
        continue;
      }
      
      manager.setGapMap(start, adjusted_length);
    } else {
      manager.setGapMap(start, length);
    }
  }
  
  if (gap_runs_filtered > 0) {
    msg("Filtered out {} suspicious or invalid gap runs", gap_runs_filtered);
  }
  
  msg("Gap map initialized with {} gap runs", manager.getGapMap().size());
  
  // Ensure no overlaps in the gap map after initialization
  auto &gap_map = const_cast<coordinates::GapMap&>(manager.getGapMap());
  if (gap_map::preventOverlaps(gap_map)) {
    msg("Fixed overlapping gaps during gap map initialization");
  }
}

// Check if a gap is suspiciously large (covering most of the sequence)
bool isSuspiciousGap(int64_t position, int64_t length, int64_t sequenceSize) {
  // A gap is suspicious if it starts at the beginning and covers >90% of the sequence
  return (position == 0 && length > sequenceSize * 0.9);
}

// Helper function to log gap map validation errors with consistent format
void logGapMapValidationError(const std::string &context, const std::string &errorMessage) {
  if (context.empty()) {
    err("Gap map validation error: {}", errorMessage);
  } else {
    err("Gap map validation error during {}: {}", context, errorMessage);
  }
}

bool validateAndFixGapMap(coordinates::CoordinateManager &manager,
                                              const std::string &context) {
  // Get gap map reference
  auto &gap_map = const_cast<coordinates::GapMap&>(manager.getGapMap());
  int64_t totalCoordinateCount = manager.size();
  
  if (gap_map.empty()) {
    // Empty map is valid, nothing to fix
    return true;
  }
  
  bool valid = true;
  bool fixesApplied = false;
  
  // Check for basic validity conditions
  for (const auto &[pos, length] : gap_map) {
    if (pos < 0) {
      logGapMapValidationError(context, fmt::format("Negative position: {}", pos));
      valid = false;
    }
    
    if (length <= 0) {
      logGapMapValidationError(context, fmt::format("Non-positive length: {} at position {}", length, pos));
      valid = false;
    }
    
    if (pos + length > totalCoordinateCount) {
      logGapMapValidationError(context, 
          fmt::format("Gap at position {} with length {} extends beyond coordinate count {}", 
                     pos, length, totalCoordinateCount));
      valid = false;
    }
  }
  
  // Check for overlaps
  auto prevIt = gap_map.begin();
  if (prevIt != gap_map.end()) {
    auto it = std::next(prevIt);
    while (it != gap_map.end()) {
      int64_t prevEnd = prevIt->first + prevIt->second - 1;
      
      if (prevEnd >= it->first) {
        logGapMapValidationError(context, 
            fmt::format("Overlapping gaps: [{},{}] and [{},{}]", 
                       prevIt->first, prevEnd, it->first, it->first + it->second - 1));
        valid = false;
      }
      
      prevIt = it;
      ++it;
    }
  }
  
  // If invalid, fix the gap map
  if (!valid) {
    // First remove any suspicious gaps
    auto suspiciousGapsCount = gap_map::removeSuspiciousGaps(gap_map, totalCoordinateCount);
    if (suspiciousGapsCount > 0) {
      debug_msg("Removed {} suspicious gaps", suspiciousGapsCount);
      fixesApplied = true;
    }
    
    // Fix overlaps by merging gaps
    if (gap_map::preventOverlaps(gap_map)) {
      debug_msg("Fixed overlapping gaps");
      fixesApplied = true;
    }
    
    // Remove any remaining invalid entries
    std::vector<int64_t> invalidPositions;
    for (const auto &[pos, length] : gap_map) {
      if (pos < 0 || length <= 0 || pos + length > totalCoordinateCount) {
        invalidPositions.push_back(pos);
      }
    }
    
    for (int64_t pos : invalidPositions) {
      gap_map.erase(pos);
      fixesApplied = true;
    }
    
    if (!invalidPositions.empty()) {
      debug_msg("Removed {} invalid gap entries", invalidPositions.size());
    }
    
    // Perform a final validation check
    valid = true;
    for (const auto &[pos, length] : gap_map) {
      if (pos < 0 || length <= 0 || pos + length > totalCoordinateCount) {
        valid = false;
        break;
      }
    }
    
    if (!valid) {
      throw std::runtime_error(fmt::format("Failed to fix gap map during {}", context));
    }
    
    if (fixesApplied) {
      debug_msg("Gap map successfully fixed during {}", context);
    }
  }
  
  return valid;
}

// Process runs of blocks with the same state (on/off/inverted) efficiently
void processBlockRun(
    std::vector<coordinates::CoordRange> &recompRanges,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
    std::vector<std::pair<int64_t, int64_t>> &gapRunBacktracks,
    std::vector<std::pair<size_t, std::pair<bool, int64_t>>> &gapMapUpdates,
    coordinates::CoordinateTraverser &traverser,
    int32_t startBlockId, int32_t endBlockId, 
    bool isPlacement) {
    
  auto &manager = traverser.getCoordManager();
  
  // Validate block range
  if (startBlockId < 0 || endBlockId < startBlockId || 
      endBlockId >= manager.getNumBlocks()) {
    err("Invalid block run range: {}-{}", startBlockId, endBlockId);
    return;
  }
  
  // See if all blocks in the run have the same on/off state
  auto &blockExists = manager.getBlockExists();
  bool allBlocksOff = true;
  bool allBlocksOn = true;
  
  for (int32_t blockId = startBlockId; blockId <= endBlockId; ++blockId) {
    bool isOn = blockExists[blockId].first;
    if (isOn) {
      allBlocksOff = false;
    } else {
      allBlocksOn = false;
    }
  }
  
  // Check if all blocks have the same orientation
  auto &blockStrand = manager.getBlockStrand();
  bool allBlocksForward = true;
  bool allBlocksInverted = true;
  
  for (int32_t blockId = startBlockId; blockId <= endBlockId; ++blockId) {
    if (!blockExists[blockId].first) continue; // Skip off blocks
    
    bool isForward = blockStrand[blockId].first;
    if (isForward) {
      allBlocksInverted = false;
    } else {
      allBlocksForward = false;
    }
  }
  
  // Get the scalar range for the entire run
  auto [startScalar, endScalar] = manager.getBlockRunScalarRange(startBlockId, endBlockId);
  
  if (startScalar == -1 || endScalar == -1) {
    err("Invalid scalar range for block run {}-{}", startBlockId, endBlockId);
    return;
  }
  
  msg("Processing block run from {}-{}, scalars {}-{}, all off: {}, all on: {}, all forward: {}, all inverted: {}",
      startBlockId, endBlockId, startScalar, endScalar, 
      allBlocksOff, allBlocksOn, allBlocksForward, allBlocksInverted);
  
  // Handle different cases based on block states
  if (allBlocksOff) {
    // For off blocks, we can handle the entire run as a single gap region
    if (!isPlacement) {
      // During indexing, mark the entire region as a gap
      
      // Add the range to updates without scanning each position
      gapRunUpdates.emplace_back(startScalar, endScalar - startScalar + 1, false);
      msg("Added entire off-block run as gap: start={}, length={}", 
          startScalar, endScalar - startScalar + 1);
    }
  } else if (allBlocksOn) {
    // For on blocks with same orientation, we can do optimized traversal
    if (allBlocksForward || allBlocksInverted) {
      // Use standard traversal instead of the template-based traverseOrientationRun
      // to avoid lambda template issues
      traverser.reset(manager.getTupleCoord(startScalar));
      
      // Create buffers to store extracted data
      std::vector<char> sequenceBuffer;
      std::vector<int64_t> scalarBuffer;
      sequenceBuffer.reserve(endScalar - startScalar + 1);
      scalarBuffer.reserve(endScalar - startScalar + 1);
      
      // Simple manual traversal loop
      while (!traverser.isDone() && traverser.getCurrent()->scalar <= endScalar) {
        const coordinates::tupleCoord_t* coord = traverser.getCurrent();
        char c = traverser.getChar(true);  // Get character with proper complementation
        
        // Process the character
        if (coord && traverser.isValidNucleotide(c)) {
          if (isPlacement) {
            // For placement, track extracted characters
            debug_msg("Placement: extracted seed character '{}' at scalar {}", c, coord->scalar);
          } else {
            // For indexing, collect the characters
            sequenceBuffer.push_back(c);
            scalarBuffer.push_back(coord->scalar);
          }
        }
        
        // Move to next position, skipping gaps efficiently
        if (traverser.getChar() == '-') {
            traverser.skipGaps();
        } else {
            traverser.next();
        }
      }
      
      // Process extracted data for indexing
      if (!isPlacement && !sequenceBuffer.empty()) {
        debug_msg("Indexing: extracted {} characters from block run {}-{}", 
              sequenceBuffer.size(), startBlockId, endBlockId);
              
        // Add to recomp ranges if we extracted usable sequence
        if (!sequenceBuffer.empty()) {
          recompRanges.emplace_back(
            manager.getTupleCoord(scalarBuffer.front()),
            manager.getTupleCoord(scalarBuffer.back())
          );
        }
      }
      
      msg("Used standard traversal for blocks {}-{}", 
          startBlockId, endBlockId);
    } else {
      // Mixed orientation blocks - use standard traversal
      msg("Mixed orientation blocks {}-{}, using standard traversal", 
          startBlockId, endBlockId);
    }
  }
}

// Enhanced version of getNucleotideSequenceFromBlockCoordinates that uses orientation run optimization
void getNucleotideSequenceFromBlockCoordinatesOptimized(
    std::string &sequence,
    std::vector<int64_t> &blockIds,
    std::vector<int64_t> &blockIdPos,
    std::vector<int> &runLengths,
    coordinates::CoordinateManager &manager,
    int64_t start_scalar, int64_t stop_scalar, const Tree *T,
    const Node *node) {
  // Function implementation should be here
}

// Process a single gap mutation record
void processGapMutationForPlacement(
    const GapMutations::Reader &gapMutation,
    PlacementNodeState &nodeState,
    coordinates::CoordinateManager &manager) {
   
  // Process the individual gap mutations from this specific record
  auto deltas = gapMutation.getDeltas();
  for (auto delta : deltas) {
    int32_t pos = delta.getPos();
    auto maybeValue = delta.getMaybeValue();
    
    if (maybeValue.isValue()) {
      // Add gap
      int64_t length = maybeValue.getValue();
      
      // Skip adding any suspicious gaps
      if (isSuspiciousGap(pos, length, manager.size())) {
        logging::warn("Skipping suspicious gap: position={}, length={}", pos, length);
        continue;
      }
      
      manager.setGapMap(pos, length);
    } else {
      // Remove gap
      manager.removeFromGapMap(pos);
    }
  }
}

// Process a single seed mutation record
void processSeedMutationForPlacement(
    PlacementGlobalState &state,
    PlacementNodeState &nodeState,
    coordinates::CoordinateTraverser &traverser,
    const SeedMutations::Reader &seedMutation,
    int64_t dfsIndex,
    std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> &readSeedCounts,
    std::unordered_map<uint64_t, int64_t> &genomeSeedCounts,
    int64_t &hitsInThisGenome) {
   
  // Process the seed mutations directly
  auto basePositions = seedMutation.getBasePositions();
  auto perPosMasks = seedMutation.getPerPosMasks();
   
  // Ensure we have both arrays and they're the same size
  if (basePositions.size() > 0 && basePositions.size() == perPosMasks.size()) {
    for (size_t i = 0; i < basePositions.size(); i++) {
      int64_t basePos = basePositions[i];
      uint64_t mask = perPosMasks[i];
      
      // Process each seed position with its mask
      for (int j = 0; j < 32; j++) {
        // Extract the ternary value for this position
        uint8_t ternaryCode = (mask >> (j * 2)) & 0x3;
        
        if (ternaryCode == 0) {
          // No change
          continue;
        }
        
        int64_t pos = basePos - j;
        
        // Skip positions outside the valid range
        if (pos < 0 || pos >= traverser.getCoordManager().size()) {
          continue;
        }
        
        if (ternaryCode == 1) {
          // Seed deleted
          if (state.onSeedsHash[pos].has_value()) {
            auto [oldSeed, oldEndPos, oldIsReverse] = state.onSeedsHash[pos].value();
            // Update hit counts by removing this seed
            auto found = genomeSeedCounts.find(oldSeed);
            if (found != genomeSeedCounts.end()) {
              found->second--;
              hitsInThisGenome--;
            }
            
            // Remove the seed from hash
            state.onSeedsHash[pos] = std::nullopt;
          }
        } 
        else if (ternaryCode == 2) {
          // Seed added or changed - handled by the existing counts
          // already processed and stored in perNodeReadSeedCounts/perNodeGenomeSeedCounts
        }
      }
    }
  }
   
  // Update hit counts based on processed seeds
  for (const auto& [seed, count] : genomeSeedCounts) {
    auto readCountIter = readSeedCounts.find(seed);
    if (readCountIter != readSeedCounts.end() && count > 0) {
      hitsInThisGenome++;
    }
  }
}
}
#include "seed_annotated_tree.hpp"
#include "gap_map.hpp"
#include "logging.hpp"
#include "performance.hpp"
#include "seeding.hpp"
#include "seq_utils.hpp"
#include <boost/icl/interval_set.hpp>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>

#include "fixed_kmer.hpp"

using namespace boost::icl;
using namespace panman;
using namespace seed_annotated_tree;
using namespace coordinates;
using namespace logging;

std::chrono::time_point<std::chrono::high_resolution_clock> global_timer =
    std::chrono::high_resolution_clock::now();

double time_stamp() {
  std::chrono::time_point<std::chrono::high_resolution_clock> newtime =
      std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> duration = newtime - global_timer;
  global_timer = newtime;
  return duration.count();
}

// Helper method that uses a hybrid approach - scalar logic within blocks,
// traversal order at boundaries
int seed_annotated_tree::getValidNucleotidesEfficiently(
    coordinates::CoordinateManager &manager, int64_t startPos, char *buffer,
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
    // Starting in a gap - use nextNonGapPosition to efficiently skip
    const tupleCoord_t *nextNonGap = manager.getNextNonGapPosition(startPos);
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
      if (nextBlockId < 0 || nextBlockId >= blockExists.size())
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
        // Skip this gap run efficiently
        int64_t gapLength = manager.getGapRunLength(currentPos);
        if (gapLength > 0) {
          // Skip the entire gap run, respecting block direction
          if (isBlockInverted) {
            // In an inverted block, we're going backwards (decreasing scalar)
            // Since gap runs are defined by start position and length in
            // forward direction, we need to find where this gap run starts when
            // going backwards
            auto upper_it = manager.getGapMap().upper_bound(currentPos);
            if (upper_it != manager.getGapMap().begin()) {
              auto run_start_it = std::prev(upper_it);
              if (currentPos >= run_start_it->first &&
                  currentPos < run_start_it->first + run_start_it->second) {
                // We're in the middle of a gap run
                // Skip to the start of this gap run (in backwards direction)
                currentPos = run_start_it->first - 1;
                continue;
              }
            }
            // If we couldn't find a gap run, just skip one position
            currentPos -= 1;
          } else {
            // In a forward block, just move forward by the gap length
            currentPos += gapLength;
          }
          continue;
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
        c = panmanUtils::getComplementCharacter(c);
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
    if (nextBlockId < 0 || nextBlockId >= blockExists.size())
      break;

    // Skip non-existent blocks
    while (nextBlockId >= 0 && nextBlockId < blockExists.size() &&
           !blockExists[nextBlockId].first) {
      nextBlockId = isBlockInverted ? nextBlockId - 1 : nextBlockId + 1;
      if (nextBlockId < 0 || nextBlockId >= blockExists.size()) {
        nextBlockId = -1;
        break;
      }
    }

    if (nextBlockId < 0 || nextBlockId >= blockExists.size())
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

bool seed_annotated_tree::getSeedAt(coordinates::CoordinateTraverser &traverser,
                                    CoordinateManager &manager,
                                    size_t &resultHash, bool &resultIsReverse,
                                    int64_t &resultEndPos, const int64_t &pos,
                                    Tree *T, const int32_t &k) {
  // When optimization is enabled, try to dispatch to specialized
  // implementations
  return fixed_kmer::dispatchGetSeedAt(
      traverser, manager, resultHash, resultIsReverse, resultEndPos, pos, T, k);
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

  traverser.reset(startCoord);

  // Get k non-gap characters and store in buffer directly in one pass
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
void seed_annotated_tree::getNucleotideSequenceFromBlockCoordinates(
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

  // Normalize scalar order if needed
  if (start_scalar > stop_scalar) {
    std::cout << "Warning: Normalizing inverted range: start_scalar="
              << start_scalar << ", stop_scalar=" << stop_scalar << std::endl;
    std::swap(start_scalar, stop_scalar);
  }

  // Check for extremely large values that might indicate overflow
  if (start_scalar > 1000000000 || stop_scalar > 1000000000) {
    std::cout << "Error: Suspiciously large scalar values: start_scalar="
              << start_scalar << ", stop_scalar=" << stop_scalar << std::endl;
    return;
  }

  // Validate gap map integrity
  if (!manager.validateGapMap()) {
    std::cout
        << "Error: Gap map validation failed. Fixing gap map before proceeding."
        << std::endl;
    // Attempt to fix the gap map by optimizing it
    manager.optimizeGapMap();

    // Check again after optimization
    if (!manager.validateGapMap()) {
      std::cout << "Error: Gap map still invalid after optimization. Traversal "
                   "may be unreliable."
                << std::endl;
    }
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
    std::cout << "Warning: Inverted or zero-length range detected"
              << " (start_scalar=" << start_scalar
              << ", stop_scalar=" << stop_scalar << ")" << std::endl;
    traversal_length = 1; // Set minimum valid length
  }

  // Estimate a reasonable initial capacity (assuming about 20% of positions are
  // gaps)
  seq.reserve(traversal_length);
  coords.reserve(traversal_length);
  gaps.reserve(traversal_length / 5);
  deadBlocks.reserve(10); // Usually there aren't many dead blocks

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
      std::cout << "End of block reached" << std::endl;
      traverser.next();
      continue;
    }

    // Determine if we're in an inverted block
    bool isInvertedBlock = currCoord->blockId < blockStrand.size() &&
                           !blockStrand[currCoord->blockId].first;

    // Check for valid nucleotide characters
    if (c != '-' && c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N' &&
        c != 'a' && c != 'c' && c != 'g' && c != 't' && c != 'n') {
      std::cout << "Warning: Invalid character '" << c
                << "' encountered at position " << currCoord->blockId << ","
                << currCoord->nucPos << "," << currCoord->nucGapPos
                << ". Skipping." << std::endl;
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
          std::cout << "Skipping " << skipCount
                    << " consecutive gap positions at " << currentBlockId << ","
                    << currentNucPos << std::endl;

          // Add the current gap position
          gaps.push_back(scalar);

          // For handling inverted blocks correctly
          int64_t endScalar = endOfRun->scalar;

          // Add the end position of the run to gaps
          gaps.push_back(endScalar);

          // Log the gap range - useful for debugging
          debug(
              "Recording gap run from scalar={} to scalar={} in block {} ({})",
              scalar, endScalar, currentBlockId,
              isInvertedBlock ? "inverted" : "forward");

          // Safety check for extremely large gap runs
          if (skipCount > 100000) {
            std::cout
                << "Warning: Extremely large gap run detected (" << skipCount
                << " positions). This may indicate a problem with the data."
                << std::endl;
          }

          // Skip to the position after the run
          if (endOfRun->next && endOfRun->next->scalar <= stop_scalar) {
            debug("Skipping to position after gap run: blockId={}, nucPos={}, "
                  "nucGapPos={}, scalar={}",
                  endOfRun->next->blockId, endOfRun->next->nucPos,
                  endOfRun->next->nucGapPos, endOfRun->next->scalar);
            traverser.reset(endOfRun->next, stop);
            continue;
          } else {
            // We've reached the end of the sequence or the stop position
            debug("Reached end of sequence or stop position after gap run");
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

        debug("Found gap at position: blockId={}, nucPos={}, nucGapPos={}, "
              "scalar={}",
              currCoord->blockId, currCoord->nucPos, currCoord->nucGapPos,
              scalar);

        // Check if this is the start of a gap run
        int64_t gap_run_length = manager.getGapRunLength(scalar);
        if (gap_run_length > 1) {
          debug("Detected gap run of length {} at scalar={}", gap_run_length,
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
          debug("Gap run: adding {} gap positions and attempting to skip to "
                "position after run",
                max_positions_to_add);

          // Skip directly to the end of this gap run
          const tupleCoord_t *next_non_gap =
              manager.getNextNonGapPosition(scalar);
          if (next_non_gap) {
            debug("Next non-gap position found: blockId={}, nucPos={}, "
                  "nucGapPos={}, scalar={}",
                  next_non_gap->blockId, next_non_gap->nucPos,
                  next_non_gap->nucGapPos, next_non_gap->scalar);

            if (next_non_gap->scalar <= stop_scalar) {
              debug("Skipping to next non-gap position within stop boundary");
              traverser.reset(next_non_gap, stop);
              continue; // Skip to next iteration with the new position
            } else {
              debug("Next non-gap position is beyond stop_scalar ({}), ending "
                    "traversal",
                    stop_scalar);
              break;
            }
          } else {
            debug("No next non-gap position found, ending traversal");
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
      traverser.skipGaps();
    } else {
      traverser.next();
    }
  }

  // If we extracted nothing, log a message
  if (seq.empty() && !gaps.empty()) {
    std::cout << "Warning: Extracted sequence is empty but found "
              << gaps.size() << " gaps" << std::endl;
  }

  // Report statistics
  std::cout << "Sequence extraction complete: " << seq.size()
            << " nucleotides, " << gaps.size() << " gaps, " << deadBlocks.size()
            << " dead blocks" << std::endl;
}

void setupSequence(const Tree *T, sequence_t &sequence,
                   blockExists_t &blockExists, blockStrand_t &blockStrand) {
  const BlockGapList &blockGaps = T->blockGaps;
  const std::vector<GapList> &gaps = T->gaps;
  const std::vector<Block> &blocks = T->blocks;

  // First find maxBlockId
  int32_t maxBlockId = 0;
  for (size_t i = 0; i < blocks.size(); i++) {
    int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
    maxBlockId = std::max(maxBlockId, primaryBlockId);
  }

  // Resize all structures once to final size
  size_t finalSize = maxBlockId + 1;
  // msg("Resizing sequence structures to {} blocks", finalSize);
  sequence.resize(finalSize);
  blockExists.resize(finalSize, {false, {}});
  blockStrand.resize(finalSize, {true, {}});

  // Assigning block gaps
  for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
    int32_t pos = blockGaps.blockPosition[i];
    if (pos >= finalSize) {
      err("Block gap position {} exceeds max block id {}", pos, maxBlockId);
      continue;
    }
    sequence[pos].second.resize(blockGaps.blockGapLength[i]);
    blockExists[pos].second.resize(blockGaps.blockGapLength[i], false);
    blockStrand[pos].second.resize(blockGaps.blockGapLength[i], true);
  }

  // Process blocks
  for (size_t i = 0; i < blocks.size(); i++) {
    int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
    int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

    if (primaryBlockId >= finalSize) {
      err("Primary block ID {} exceeds max block id {}", primaryBlockId,
          maxBlockId);
      continue;
    }

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

        if (secondaryBlockId != -1) {
          sequence[primaryBlockId].second[secondaryBlockId].push_back(
              {nucleotide, {}});
        } else {
          sequence[primaryBlockId].first.push_back({nucleotide, {}});
        }
      }
      if (endFlag) {
        break;
      }
    }

    sequence[primaryBlockId].first.push_back({'x', {}});
  }

  // Assigning nucleotide gaps
  for (size_t i = 0; i < gaps.size(); i++) {
    int32_t primaryBId = (gaps[i].primaryBlockId);
    int32_t secondaryBId = (gaps[i].secondaryBlockId);

    if (primaryBId >= finalSize) {
      err("Gap primary block ID {} exceeds max block id {}", primaryBId,
          maxBlockId);
      continue;
    }

    for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
      int len = gaps[i].nucGapLength[j];
      int pos = gaps[i].nucPosition[j];

      if (secondaryBId != -1) {
        sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
      } else {
        sequence[primaryBId].first[pos].second.resize(len, '-');
      }
    }
  }
}

void seed_annotated_tree::setupPlacement(
    std::vector<std::optional<seeding::onSeedsHash>> &onSeedsHash,
    sequence_t &sequence, blockExists_t &blockExists,
    blockStrand_t &blockStrand, const Tree *T) {

  setupSequence(T, sequence, blockExists, blockStrand);

  onSeedsHash.resize(sequence.size());
}

void seed_annotated_tree::setupIndexing(sequence_t &sequence,
                                        blockExists_t &blockExists,
                                        blockStrand_t &blockStrand,
                                        const Tree *T) {

  setupSequence(T, sequence, blockExists, blockStrand);
}

static int getIndexFromNucleotide(char nuc) {

  switch (nuc) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  case '*':
    return 4;
  default:
    return 5;
  }
  return 5;
}
static size_t getBeg(const std::string &s1, const std::string &s2,
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
    if (s1[idx] == '-' && s2[idx] == '-') {
      ++idx;
      continue;
    }
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
seed_annotated_tree::getMaskCoorsForMutmat(const std::string &s1,
                                           const std::string &s2, size_t window,
                                           double threshold) {

  assert(s1.size() == s2.size());
  if (window == 0 || threshold == 0.0) {
    return std::make_pair<size_t, size_t>(0, s1.size() - 1);
  }
  return std::make_pair<size_t, size_t>(getBeg(s1, s2, window, threshold),
                                        getEnd(s1, s2, window, threshold));
}

void seed_annotated_tree::writeMutationMatrices(const mutationMatrices &mutMat,
                                                std::ofstream &mmfout) {
  TIME_FUNCTION;
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

void seed_annotated_tree::applyMutations(
    blockMutationInfo_t &blockMutationInfo,
    std::vector<coordinates::CoordRange> &recompRanges,
    mutationInfo_t &mutationInfo,
    std::vector<gap_map::GapUpdate> &gapRunUpdates,
    std::vector<gap_map::GapUpdate> &gapRunBacktracks,
    std::vector<gap_map::GapUpdate> &gapMapUpdates, panmanUtils::Tree *T,
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

    // Process block mutation and update state
    if (type == 1) {
      // insertion
      bool oldStrand;
      bool oldMut;
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
      bool oldMut;
      bool oldStrand;
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
  if (!isPlacement) {
    // Then handle nucleotide mutations
    for (const auto &mut : recompMutations) {
      if (blockMutated[mut.blockId]) {
        continue; // already created recomp range for this block
      }

      try {
        // Create start coordinate
        tupleCoord_t startCoord{mut.blockId, mut.startPos, mut.gapPos};
        const tupleCoord_t *start = manager.findTupleCoord(startCoord);

        traverser.reset(start);
        traverser.prev();
        start = traverser.getCurrent();

        if (!start)
          continue;

        // Create stop coordinate
        tupleCoord_t stopCoord;

        if (mut.gapPos != -1) {
          stopCoord.blockId = mut.blockId;
          stopCoord.nucPos = mut.startPos;
          stopCoord.nucGapPos = std::min(
              mut.gapPos + mut.length + k,
              (int)sequence[mut.blockId].first[mut.startPos].second.size() - 1);
        } else {
          stopCoord.blockId = mut.blockId;
          stopCoord.nucPos =
              std::min(mut.startPos + mut.length + k,
                       (int)sequence[mut.blockId].first.size() - 1);
          stopCoord.nucGapPos = -1;
        }

        // vector lookup, fast
        const tupleCoord_t *stop = manager.findTupleCoord(stopCoord);

        traverser.reset(stop);
        traverser.next();
        stop = traverser.getCurrent();

        if (stop) {
          recompRanges.emplace_back(start, stop);
        } else {
          err("Bad stop lookup {},{},{} for block {}", stopCoord.blockId,
              stopCoord.nucPos, stopCoord.nucGapPos, mut.blockId);
          throw std::runtime_error("Bad stop");
        }
      } catch (const std::exception &e) {
        err("Failed to create mutation range: {}", e.what());
      }
    }
  }

  // Process gap updates and normalize ranges before updating gap map
  for (auto &update : gapRunUpdates) {
    auto &[erase, range] = update;
    auto &[pos, length] = range;

    // Skip invalid updates
    if (pos < 0 || pos >= manager.size()) {
      warn("Skipping gap update with invalid position {} (max: {})", pos,
           manager.size() - 1);
      continue;
    }

    // For non-erase operations, ensure length is positive
    if (!erase && length <= 0) {
      warn("Skipping gap update with invalid length {} at position {}", length,
           pos);
      continue;
    }
  }

  // Update gap map with current mutations
  manager.updateGapMap(gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  // Process gap mutations for proper backtracking
  for (auto &mutation : gapRunUpdates) {
    const auto &[erase, range] = mutation;
    const auto &[pos, length] = range;

    if (erase) {
      // This is removing a gap run
      auto it = manager.getGapMap().find(pos);
      if (it != manager.getGapMap().end()) {
        // Record this for backtracking - we're modifying an existing entry
        gapRunBacktracks.emplace_back(false, std::make_pair(pos, it->second));
        // Record the update for later use
        gapMapUpdates.emplace_back(true, std::make_pair(pos, it->second));
      }
    } else {
      // This is adding or modifying a gap run
      auto it = manager.getGapMap().find(pos);
      if (it != manager.getGapMap().end()) {
        // Record this for backtracking - we're modifying an existing entry
        gapRunBacktracks.emplace_back(false, std::make_pair(pos, it->second));
      } else {
        // Record this for backtracking - we're adding a new entry
        gapRunBacktracks.emplace_back(true, std::make_pair(pos, 0));
      }
      // Record the update for later use
      gapMapUpdates.emplace_back(false, std::make_pair(pos, length));
    }
  }
}

void seed_annotated_tree::undoMutations(
    Tree *T, const Node *node, const blockMutationInfo_t &blockMutationInfo,
    const mutationInfo_t &mutationInfo,
    coordinates::CoordinateTraverser &traverser,
    std::vector<gap_map::GapUpdate> &gapRunBacktracks) {
  auto &sequence = traverser.getCoordManager().getSequence();
  auto &blockExists = traverser.getCoordManager().getBlockExists();
  auto &blockStrand = traverser.getCoordManager().getBlockStrand();

  // Counters for tracking validity issues
  int invalid_blockId = 0;
  int invalid_secondaryBlockId = 0;
  int invalid_nucPos = 0;
  int invalid_nucGapPos = 0;

  // For debugging
  msg("Undoing {} block mutations and {} nucleotide mutations for node {}",
      blockMutationInfo.size(), mutationInfo.size(), node->identifier);

  // STEP 1: Verify gap map integrity before restoring
  if (!traverser.getCoordManager().validateGapMap()) {
    warn("Gap map validation failed before undo. Optimizing gap map.");
    traverser.getCoordManager().optimizeGapMap();
  }

  // STEP 2: Restore the gap map to its previous state using backtracking
  // information
  msg("Restoring gap map state from {} backtrack entries",
      gapRunBacktracks.size());

  // Enhanced restoration - apply backtracking entries more safely
  try {
    // Backup current gap map in case restoration fails
    auto originalGapMap = traverser.getCoordManager().getGapMap();

    // Try to restore the gap map
    traverser.getCoordManager().restoreGapMapFromBacktracks(gapRunBacktracks);

    // Verify the restored gap map is valid
    if (!traverser.getCoordManager().validateGapMap(false)) {
      warn("Gap map validation failed after restoring previous state. Will "
           "attempt repair.");
      traverser.getCoordManager().optimizeGapMap();
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
    msg("Detected {} blocks with strand changes that need special handling",
        strandChangedBlocks.size());
  }

  // STEP 5: Apply the block mutations in reverse order
  msg("Undoing block mutations");
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
        debug("Successfully reinitialized block {}", blockId);
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
    msg("Created {} backtrack entries while re-inverting gaps for strand "
        "changes",
        strandChangeBacktracks.size());
  }

  // STEP 7: Undo nucleotide mutations in reverse order from application
  msg("Undoing nucleotide mutations");
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
          warn("Invalid blockId {} when undoing mutation (max: {})", blockId,
               sequence.size() - 1);
        }
        continue;
      }

      // Skip if block doesn't exist
      if (!blockExists[blockId].first) {
        debug("Skipping nucleotide mutation for dead block {}", blockId);
        continue;
      }

      // Separate handling for secondary blocks vs primary blocks
      if (secondaryBlockId != -1) {
        // Check if secondaryBlockId is valid
        if (secondaryBlockId < 0 ||
            secondaryBlockId >= sequence[blockId].second.size()) {
          invalid_secondaryBlockId++;
          if (invalid_secondaryBlockId <= 5) {
            warn("Invalid secondaryBlockId {} for blockId {}", secondaryBlockId,
                 blockId);
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
              warn("Invalid nucPos {} for blockId {}, secondaryBlockId {}",
                   nucPos, blockId, secondaryBlockId);
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
              warn("Invalid nucGapPos {} for blockId {}, secondaryBlockId {}, "
                   "nucPos {}",
                   nucGapPos, blockId, secondaryBlockId, nucPos);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue = sequence[blockId]
                                  .second[secondaryBlockId][nucPos]
                                  .second[nucGapPos];
          if (currentValue != newValue && newValue != '\0') {
            warn("Value mismatch during undo at [{},{},{},{}]: expected '{}', "
                 "found '{}'",
                 blockId, secondaryBlockId, nucPos, nucGapPos, newValue,
                 currentValue);
          }
          sequence[blockId].second[secondaryBlockId][nucPos].second[nucGapPos] =
              oldValue;

        } else {
          // Validate nucPos
          if (nucPos < 0 ||
              nucPos >= sequence[blockId].second[secondaryBlockId].size()) {
            invalid_nucPos++;
            if (invalid_nucPos <= 5) {
              warn("Invalid nucPos {} for blockId {}, secondaryBlockId {}",
                   nucPos, blockId, secondaryBlockId);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue =
              sequence[blockId].second[secondaryBlockId][nucPos].first;
          if (currentValue != newValue && newValue != '\0') {
            warn("Value mismatch during undo at [{},{},{}]: expected '{}', "
                 "found '{}'",
                 blockId, secondaryBlockId, nucPos, newValue, currentValue);
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
              warn("Invalid nucPos {} for blockId {}", nucPos, blockId);
            }
            continue;
          }

          // Validate nucGapPos
          if (nucGapPos < 0 ||
              nucGapPos >= sequence[blockId].first[nucPos].second.size()) {
            invalid_nucGapPos++;
            if (invalid_nucGapPos <= 5) {
              warn("Invalid nucGapPos {} for blockId {}, nucPos {}", nucGapPos,
                   blockId, nucPos);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue = sequence[blockId].first[nucPos].second[nucGapPos];
          if (currentValue != newValue && newValue != '\0') {
            warn("Value mismatch during undo at [{},{},{}]: expected '{}', "
                 "found '{}'",
                 blockId, nucPos, nucGapPos, newValue, currentValue);
          }
          sequence[blockId].first[nucPos].second[nucGapPos] = oldValue;

        } else {
          // Validate nucPos
          if (nucPos < 0 || nucPos >= sequence[blockId].first.size()) {
            invalid_nucPos++;
            if (invalid_nucPos <= 5) {
              warn("Invalid nucPos {} for blockId {}", nucPos, blockId);
            }
            continue;
          }

          // Safe access - restore original value
          char currentValue = sequence[blockId].first[nucPos].first;
          if (currentValue != newValue && newValue != '\0') {
            warn("Value mismatch during undo at [{},{}]: expected '{}', found "
                 "'{}'",
                 blockId, nucPos, newValue, currentValue);
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
    warn("Skipped invalid positions during undo - "
         "blockId: {}, secondaryBlockId: {}, nucPos: {}, nucGapPos: {}",
         invalid_blockId, invalid_secondaryBlockId, invalid_nucPos,
         invalid_nucGapPos);
  } else {
    msg("All nucleotide mutations successfully undone");
  }

  // STEP 8: Final verification step and gap map sanity check
  if (!traverser.getCoordManager().validateGapMap(true)) {
    warn("Gap map has issues after undoing all mutations - attempting to "
         "optimize");
    traverser.getCoordManager().optimizeGapMap();
  }

  // STEP 9: One final pass to ensure block coordinates are properly linked
  try {
    msg("Performing final coordinate system cleanup");
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
    warn("Exception during final coordinate cleanup: {}", e.what());
  }

  msg("Mutation undo complete for node {}", node->identifier);
}

void buildMutationMatricesHelper_test(
    mutationMatrices &mutMat, Tree *T, Node *node,
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

void seed_annotated_tree::fillMutationMatricesFromTree_test(
    mutationMatrices &mutMat, Tree *T, const std::string &path) {
  return;
}

void seed_annotated_tree::fillMutationMatricesFromFile(mutationMatrices &mutMat,
                                                       std::ifstream &inf) {
  TIME_FUNCTION;
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
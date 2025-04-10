#include "coordinates.hpp"
#include "logging.hpp"
#include "seed_annotated_tree.hpp"
#include <algorithm>
#include <sstream>
#include <stdexcept>

using namespace logging;

namespace coordinates {

// Define the debug variable
bool debug = false;

bool CoordinateManager::validateCoordinateState(
    const tupleCoord_t &coord, const sequence_t &sequence) const noexcept {
  if (!coord.isValidBasic())
    return false;
  if (coord.blockId > max_block_id)
    return false;
  if (coord.nucPos > max_nuc_pos)
    return false;
  if (coord.nucGapPos > max_gap_pos)
    return false;
  return true;
}

// Add comprehensive validation method
bool CoordinateManager::validateFullState() const {
  // Validate all coordinates
  for (size_t i = 0; i < coords.size(); i++) {
    if (!validateCoordinate(coords[i])) {
      err("validateFullState: Coordinate {} failed validation", i);
      return false;
    }
  }
  
  // Validate gap map consistency with actual sequence
  for (const auto& [start, length] : gap_map) {
    // Check that the gap run is within bounds
    int64_t end = start + length - 1;
    if (start < 0 || end >= size()) {
      err("validateFullState: Gap run [{}, {}] (length {}) exceeds coordinates range [0, {}]", 
          start, end, length, size() - 1);
      return false;
    }
    
    // Sample some positions to verify they are gaps (checking all would be too slow)
    // Check first, last, and some middle positions
    std::vector<int64_t> positions_to_check = {start, end};
    if (length > 2) {
      positions_to_check.push_back(start + length/2); // middle
    }
    
    for (int64_t pos : positions_to_check) {
      const tupleCoord_t* coord = getTupleCoord(pos);
      if (!coord) {
        err("validateFullState: No coordinate at scalar position {}", pos);
        return false;
      }
      
      // Get character at this position
      char c;
      if (coord->nucGapPos == -1) {
        c = sequence->at(coord->blockId).first[coord->nucPos].first;
      } else {
        c = sequence->at(coord->blockId).first[coord->nucPos].second[coord->nucGapPos];
      }
      
      if (c != '-') {
        err("validateFullState: Position {} in gap run [{},{}] is not a gap ('{}')", 
            pos, start, end, c);
        return false;
      }
    }
  }
  
  // Validate block ranges
  for (size_t blockId = 0; blockId < block_ranges.size(); blockId++) {
    const auto& [start, end] = block_ranges[blockId];
    if (start < 0 || end < start || end >= coords.size()) {
      err("validateFullState: Invalid block range for block {}: [{}, {}]", 
          blockId, start, end);
      return false;
    }
  }
  
  // Validate traversal pointers
  for (size_t i = 0; i < coords.size(); i++) {
    const tupleCoord_t& coord = coords[i];
    
    // Check next pointer
    if (coord.next) {
      if (coord.next->prev != &coord) {
        err("validateFullState: Inconsistent next/prev pointers at position {}", i);
        return false;
      }
      
      // Check if next pointer points to a valid coordinate
      if (!validateCoordinate(*coord.next)) {
        err("validateFullState: Next pointer at position {} points to invalid coordinate", i);
        return false;
      }
    }
    
    // Check prev pointer
    if (coord.prev) {
      if (coord.prev->next != &coord) {
        err("validateFullState: Inconsistent prev/next pointers at position {}", i);
        return false;
      }
      
      // Check if prev pointer points to a valid coordinate
      if (!validateCoordinate(*coord.prev)) {
        err("validateFullState: Prev pointer at position {} points to invalid coordinate", i);
        return false;
      }
    }
  }
  
  // All checks passed
  return true;
}

void CoordinateManager::setupCoordinates() {
  // msg("Starting coordinate setup");
  int64_t ctr = 0;

  // For gap map initialization
  int64_t current_gap_start = -1;
  int64_t current_gap_length = 0;

  // First pass: Initialize all coordinates with scalars
  for (size_t blockId = 0; blockId < sequence->size(); blockId++) {
    const auto &block = sequence->at(blockId);
    if (block.first.empty()) {
      // msg("Block {} is empty, skipping", blockId);
      continue;
    }

    // Initialize coordinates for this block
    int64_t blockStartScalar = ctr;
    // msg("Initializing block {}: start_scalar={}, positions={}",
    //     blockId, blockStartScalar, block.first.size());

    // Process each position in order
    for (size_t nucPos = 0; nucPos < block.first.size(); nucPos++) {
      const auto &nucEntry = block.first[nucPos];
      const auto &gapVec = nucEntry.second;

      // First add gap positions if they exist
      for (size_t gapPos = 0; gapPos < gapVec.size(); gapPos++) {
        if (ctr >= total_size) {
          throw std::runtime_error(
              fmt::format("Coordinate counter exceeded total_size at block {} "
                          "pos {} gap {} (ctr={}, total={})",
                          blockId, nucPos, gapPos, ctr, total_size));
        }

        tupleCoord_t coord(blockId, nucPos, gapPos, ctr);
        if (!validateCoordinateState(coord, *sequence)) {
          throw std::runtime_error(
              fmt::format("Invalid gap coordinate state at block {} pos {} gap "
                          "{} scalar {}",
                          blockId, nucPos, gapPos, ctr));
        }

        // Track gap runs for gap map
        if (gapVec[gapPos] == '-') {
          if (current_gap_start == -1) {
            current_gap_start = ctr;
            current_gap_length = 1;
          } else {
            current_gap_length++;
          }
        } else if (current_gap_start != -1) {
          // End of a gap run
          gap_map[current_gap_start] = current_gap_length;
          current_gap_start = -1;
          current_gap_length = 0;
        }

        coords.push_back(coord);
        scalarCoordToBlockId[ctr] = blockId;
        ctr++;
      }

      // Then add nucleotide position
      if (ctr >= total_size) {
        throw std::runtime_error(
            fmt::format("Coordinate counter exceeded total_size at block {} "
                        "pos {} (ctr={}, total={})",
                        blockId, nucPos, ctr, total_size));
      }

      tupleCoord_t coord(blockId, nucPos, -1, ctr);
      if (!validateCoordinateState(coord, *sequence)) {
        throw std::runtime_error(fmt::format(
            "Invalid nucleotide coordinate state at block {} pos {} scalar {}",
            blockId, nucPos, ctr));
      }

      // Check if this nucleotide position is a gap
      if (nucEntry.first == '-') {
        if (current_gap_start == -1) {
          current_gap_start = ctr;
          current_gap_length = 1;
        } else {
          current_gap_length++;
        }
      } else if (current_gap_start != -1) {
        // End of a gap run
        gap_map[current_gap_start] = current_gap_length;
        current_gap_start = -1;
        current_gap_length = 0;
      }

      coords.push_back(coord);
      scalarCoordToBlockId[ctr] = blockId;
      ctr++;
    }

    // Update block range
    if (blockStartScalar >= 0 && ctr > blockStartScalar) {
      block_ranges[blockId] = {blockStartScalar, ctr - 1};
      // msg("Block {} range initialized: [{},{}], positions={},
      // total_coords={}",
      //     blockId, blockStartScalar, ctr-1, block.first.size(), ctr -
      //     blockStartScalar);
    } else {
      err("Invalid block range computed: start={}, end={}", blockStartScalar,
          ctr - 1);
    }
  }

  // If we ended with an active gap run, add it to the map
  if (current_gap_start != -1) {
    gap_map[current_gap_start] = current_gap_length;
  }

  max_scalar = ctr - 1;

  // msg("Initial coordinate setup complete: {} coordinates created",
  // coords.size());

  // Initialize pointers between coordinates
  for (size_t i = 0; i < coords.size() - 1; i++) {
    coords[i].setNext(&coords[i + 1]);
    coords[i + 1].setPrev(&coords[i]);
  }

  // Validate all coordinates and ranges
  for (size_t i = 0; i < coords.size(); i++) {
    if (!coords[i].isValidBasic()) {
      throw std::runtime_error(
          fmt::format("Invalid coordinate at position {}: block={}, pos={}, "
                      "gap={}, scalar={}",
                      i, coords[i].blockId, coords[i].nucPos,
                      coords[i].nucGapPos, coords[i].scalar));
    }

    if (coords[i].blockId >= sequence->size()) {
      throw std::runtime_error(
          fmt::format("Block ID {} exceeds sequence size {} at position {}",
                      coords[i].blockId, sequence->size(), i));
    }

    const auto &block = sequence->at(coords[i].blockId);
    if (coords[i].nucPos >= block.first.size()) {
      throw std::runtime_error(fmt::format(
          "Nucleotide position {} exceeds block size {} at position {}",
          coords[i].nucPos, block.first.size(), i));
    }

    if (coords[i].nucGapPos >= 0 &&
        coords[i].nucGapPos >= block.first[coords[i].nucPos].second.size()) {
      throw std::runtime_error(fmt::format(
          "Gap position {} exceeds gap vector size {} at position {}",
          coords[i].nucGapPos, block.first[coords[i].nucPos].second.size(), i));
    }
  }

  for (size_t blockId = 0; blockId < sequence->size(); blockId++) {
    reinitializeBlockCoordinates(blockId);
  }

}

CoordinateTraverser::CoordinateTraverser(const tupleCoord_t *start,
                                         CoordinateManager *manager)
    : current(start), end(nullptr), coordManager(manager), done(false) {
  if (!start || !manager) {
    err("CoordinateTraverser: null start coordinate or manager");
    done = true;
    return;
  }

  // Ensure coordinates are initialized
  if (!manager->isInitialized()) {
    err("CoordinateTraverser: coordinate manager not initialized");
    done = true;
    return;
  }

  // Validate start coordinate
  if (!start->isValidBasic()) {
    err("CoordinateTraverser: invalid start coordinate: block={}, pos={},{}, "
        "scalar={}",
        start->blockId, start->nucPos, start->nucGapPos, start->scalar);
    done = true;
    return;
  }

  // Validate block ID is in range
  if (start->blockId < 0 || start->blockId >= manager->getNumBlocks()) {
    err("CoordinateTraverser: start coordinate block ID {} out of range [0,{}]",
        start->blockId, manager->getNumBlocks() - 1);
    done = true;
    return;
  }

  // Validate scalar is in range
  if (start->scalar < 0 || start->scalar >= manager->size()) {
    err("CoordinateTraverser: start coordinate scalar {} out of range [0,{}]",
        start->scalar, manager->size() - 1);
    done = true;
    return;
  }

  // Validate coordinate is in the manager's vector
  const tupleCoord_t *stored = manager->getTupleCoord(start->scalar);
  if (!stored || stored != start) {
    err("CoordinateTraverser: start coordinate not found in manager's storage");
    done = true;
    return;
  }

  if (!manager->validateCoordinate(*start)) {
    err("CoordinateTraverser: coordinate validation failed");
    done = true;
    return;
  }
}

CoordinateTraverser::CoordinateTraverser()
    : current(nullptr), end(nullptr), coordManager(nullptr), done(true) {}

// Go upstream (against traversal direction) until neededNongap nucleotides are
// seen
const tupleCoord_t *expandUpstream(CoordinateTraverser &traverser,
                                   const tupleCoord_t *coord, int neededNongap,
                                   blockExists_t &blockExists,
                                   blockStrand_t &blockStrand,
                                   const tupleCoord_t *stop_coord) {

  auto &manager = traverser.getCoordManager();
  if (!coord)
    return nullptr;

  const tupleCoord_t *result = coord;
  int count = 0;

  // Reset traverser
  if (stop_coord) {
    traverser.reset(coord, stop_coord);
  } else {
    traverser.reset(coord);
  }

  // Safety check for boundary conditions
  if (coord->scalar <= 0 || coord->scalar >= manager.size()) {
    warn("expandUpstream: coord scalar {} is out of bounds [0,{}]", 
        coord->scalar, manager.size() - 1);
    return result;
  }

  // Check if we're in an inverted block
  bool isInverted = traverser.isInverted();
  debug_msg("expandUpstream: Starting at scalar {} in {} block", 
       coord->scalar, isInverted ? "inverted" : "normal");

  // Skip any gaps at the current position using the gap map
  if (manager.isGapPosition(coord->scalar)) {
    debug_msg("expandUpstream: Starting position is a gap, skipping gaps");
    // Use the appropriate gap-skipping method based on block orientation
    traverser.skipGapsBackward(); // This now handles inversions internally
    if (traverser.isDone()) {
      return result;
    }
    result = traverser.getCurrent();
  }

  // Safety counter to prevent infinite loops
  int safety_counter = 0;
  int max_iterations = 10000;

  while (count < neededNongap && !traverser.isDone() && safety_counter < max_iterations) {
    const tupleCoord_t *pos = traverser.getCurrent();
    if (!pos) {
      err("expandUpstream: got null position from traverser");
      break;
    }

    safety_counter++;
    if (safety_counter >= max_iterations) {
      err("expandUpstream: reached maximum iterations ({})", max_iterations);
      break;
    }

    // Validate block ID before access
    if (pos->blockId < 0 || pos->blockId >= blockExists.size()) {
      err("expandUpstream: invalid block ID {} (max: {})", pos->blockId,
          blockExists.size() - 1);
      break;
    }

    // Check if current position is out of bounds
    if (pos->scalar < 0 || pos->scalar >= manager.size()) {
      err("expandUpstream: scalar {} is out of bounds [0,{}]", 
          pos->scalar, manager.size() - 1);
      break;
    }

    // Now safe to check block existence
    if (!blockExists[pos->blockId].first) {
      debug_msg("expandUpstream: Block {} is off, moving to previous block", pos->blockId);
      // Block is off, move to previous block
      if (pos->blockId == 0) {
        return result;
      }
      traverser.moveToPreviousBlock();
      continue;
    }

    // Stop if we've reached stop_coord
    if (stop_coord && pos->scalar == stop_coord->scalar) {
      debug_msg("expandUpstream: Reached stop coordinate at scalar {}", pos->scalar);
      break;
    }

    // Count non-gap positions
    char c = traverser.getChar();
    if (c == 'x') {
      debug_msg("expandUpstream: Found end marker 'x' at scalar {}", pos->scalar);
      return result;
    }
    if (c != '-') {
      count++;
      result = pos;
      debug_msg("expandUpstream: Found non-gap character '{}' at scalar {}, count={}/{}", 
            c, pos->scalar, count, neededNongap);
    } else {
      debug_msg("expandUpstream: Found gap at scalar {}, using skipGapsBackward", pos->scalar);
      // Use efficient gap skipping that properly handles inversions
      traverser.skipGapsBackward(); 
      continue;
    }

    traverser.prev();
  }

  debug_msg("expandUpstream: Returning result at scalar {}, found {} non-gap positions",
        result->scalar, count);
  return result;
}

// Go downstream (following traversal direction) until neededNongap nucleotides
// are seen
const tupleCoord_t *expandDownstream(CoordinateTraverser &traverser,
                                     const tupleCoord_t *coord,
                                     int neededNongap,
                                     blockExists_t &blockExists,
                                     blockStrand_t &blockStrand) {

  if (!coord) {
    err("expandDownstream: null coordinate");
    return nullptr;
  }

  auto &manager = traverser.getCoordManager();
  
  // Safety check for boundary conditions
  if (coord->scalar < 0 || coord->scalar >= manager.size()) {
    warn("expandDownstream: coord scalar {} is out of bounds [0,{}]", 
        coord->scalar, manager.size() - 1);
    return coord;
  }

  int count = 0;
  const tupleCoord_t *result = coord;

  traverser.reset(coord);

  // Check if we're in an inverted block
  bool isInverted = traverser.isInverted();
  debug_msg("expandDownstream: Starting at scalar {} in {} block", 
       coord->scalar, isInverted ? "inverted" : "normal");

  // Skip any gaps at the current position using the gap map
  if (manager.isGapPosition(coord->scalar)) {
    debug_msg("expandDownstream: Starting position is a gap, skipping gaps");
    traverser.skipGaps(); // This now handles inversions internally
    if (traverser.isDone()) {
      return result;
    }
    result = traverser.getCurrent();
  }

  // Safety counter to prevent infinite loops
  int safety_counter = 0;
  int max_iterations = 10000;

  while (count < neededNongap && !traverser.isDone() && safety_counter < max_iterations) {
    const tupleCoord_t *pos = traverser.getCurrent();
    if (!pos) {
      err("expandDownstream: null position from traverser");
      break;
    }
    
    safety_counter++;
    if (safety_counter >= max_iterations) {
      err("expandDownstream: reached maximum iterations ({})", max_iterations);
      break;
    }

    // Check if current position is out of bounds
    if (pos->scalar < 0 || pos->scalar >= manager.size()) {
      err("expandDownstream: scalar {} is out of bounds [0,{}]", 
          pos->scalar, manager.size() - 1);
      break;
    }

    // Validate block ID
    if (pos->blockId < 0 || pos->blockId >= blockExists.size()) {
      err("expandDownstream: invalid block ID {} (max: {})", pos->blockId,
          blockExists.size() - 1);
      break;
    }

    // Handle dead blocks
    if (!blockExists[pos->blockId].first) {
      debug_msg("expandDownstream: Block {} is off, moving to next block", pos->blockId);
      if (pos->blockId >= traverser.getCoordManager().getNumBlocks() - 1) {
        return result;
      }
      traverser.moveToNextBlock();
      continue;
    }

    // Count non-gap positions
    char c = traverser.getChar();
    if (c == 'x') {
      debug_msg("expandDownstream: Found end marker 'x' at scalar {}", pos->scalar);
      return result;
    }
    if (c != '-') {
      count++;
      result = pos;
      debug_msg("expandDownstream: Found non-gap character '{}' at scalar {}, count={}/{}", 
            c, pos->scalar, count, neededNongap);
    } else {
      debug_msg("expandDownstream: Found gap at scalar {}, using skipGaps", pos->scalar);
      // Use the skip gaps optimization which now properly handles inversions
      traverser.skipGaps();
      continue;
    }

    traverser.next();
  }

  debug_msg("expandDownstream: Returning result at scalar {}, found {} non-gap positions",
        result->scalar, count);
  return result;
}

void expandRanges(std::vector<CoordRange> &ranges,
                  CoordinateTraverser &traverser, int k) {
  // msg("Expanding {} ranges with k={}", ranges.size(), k);

  // Remove invalid ranges
  ranges.erase(
      std::remove_if(ranges.begin(), ranges.end(),
                     [&](const CoordRange &range) {
                       if (!range.isValid()) {
                         err("Removing invalid range before expansion");
                         return true;
                       }
                       return false;
                     }),
      ranges.end());

  std::vector<CoordRange> validRanges;
  validRanges.reserve(ranges.size());

  for (const auto &range : ranges) {
    const tupleCoord_t *start = range.getStart(traverser);
    const tupleCoord_t *stop = range.getStop(traverser);

    if (!start || !stop) {
      err("Invalid range pointers - skipping range");
      continue;
    }

    const tupleCoord_t *newStart = expandUpstream(
        traverser, start, k, traverser.getCoordManager().getBlockExists(),
        traverser.getCoordManager().getBlockStrand(), stop);

    if (!newStart) {
      err("Failed to expand range upstream - skipping range");
      continue;
    }

    const tupleCoord_t *newStop = expandDownstream(
        traverser, stop, k, traverser.getCoordManager().getBlockExists(),
        traverser.getCoordManager().getBlockStrand());

    if (!newStop) {
      err("Failed to expand range downstream - skipping range");
      continue;
    }

    validRanges.emplace_back(newStart, newStop);
    // msg("Successfully expanded range: block {} [{},{}] -> block {} [{},{}]",
    //     start->blockId, start->scalar, stop->scalar,
    //     newStart->blockId, newStart->scalar, newStop->scalar);
  }

  ranges = std::move(validRanges);
  // msg("Expansion complete - {} valid ranges remain", ranges.size());
}

void mergeRanges(std::vector<CoordRange> &ranges,
                 CoordinateTraverser &traverser) {
  if (ranges.empty())
    return;

  // msg("Merging {} ranges", ranges.size());

  // Sort ranges by start scalar
  std::sort(ranges.begin(), ranges.end());

  // Index for where to place the next merged range
  size_t mergedIdx = 0;

  // Iterate through remaining ranges, merging as we go
  for (size_t i = 1; i < ranges.size(); ++i) {
    auto &current = ranges[mergedIdx];
    const auto &next = ranges[i];

    if (!current.isValid() || !next.isValid()) {
      err("Invalid range encountered during merge at index {}", i);
      continue;
    }

    // msg("Comparing ranges - current: [{},{}], next: [{},{}]",
    //     current.start_scalar, current.stop_scalar,
    //     next.start_scalar, next.stop_scalar);

    // Check if ranges overlap or are adjacent
    if (current.stop_scalar >= next.start_scalar - 1) {
      // Extend current range if needed
      if (current.stop_scalar < next.stop_scalar) {
        // msg("Extending current range with next range");
        current.stop_scalar = next.stop_scalar;
      }
    } else {
      // Ranges don't overlap, start a new merged range
      mergedIdx++;
      if (mergedIdx < i) {
        // msg("Starting new range at index {}", mergedIdx);
        ranges[mergedIdx] = next;
      }
    }
  }

  // Resize the vector to remove any unused ranges
  size_t oldSize = ranges.size();
  ranges.resize(mergedIdx + 1);
  // msg("Merged {} ranges into {} ranges", oldSize, ranges.size());
}

bool canMergeBlocks(int blockId1, int blockId2,
                    const blockExists_t &blockExists) {
  if (blockId1 == blockId2)
    return true;

  // Ensure blocks are in order
  if (blockId1 > blockId2) {
    std::swap(blockId1, blockId2);
  }

  // Find next existing block after blockId1
  int nextLiveBlock = blockId1;
  while (nextLiveBlock < blockExists.size()) {
    nextLiveBlock++;
    if (nextLiveBlock >= blockExists.size())
      break;
    if (blockExists[nextLiveBlock].first) {
      // Found next live block
      if (nextLiveBlock == blockId2)
        return true; // Direct connection
      if (nextLiveBlock > blockId2)
        return false; // Overshot - there's a live block in between
      break;          // Found a different live block - can't merge
    }
  }
  return false;
}

// Helper function to find next/prev live blocks
std::pair<int, int> findAdjacentLiveBlocks(int blockId,
                                           const blockExists_t &blockExists) {
  int prevBlock = blockId - 1;
  while (prevBlock >= 0 && !blockExists[prevBlock].first) {
    prevBlock--;
  }

  int nextBlock = blockId + 1;
  while (nextBlock < blockExists.size() && !blockExists[nextBlock].first) {
    nextBlock++;
  }

  return {prevBlock, nextBlock};
}

CoordRange createBlockRange(CoordinateManager *manager, int32_t blockId,
                            const blockStrand_t &blockStrand) {
  const tupleCoord_t *start = manager->getFirstCoordinateInBlock(blockId);
  const tupleCoord_t *stop = manager->getLastCoordinateInBlock(blockId);

  if (!start || !stop) {
    throw std::runtime_error("Failed to create block range: null coordinates");
  }

  if (!blockStrand[blockId].first) {
    std::swap(start, stop);
  }

  return CoordRange{start, stop};
}

bool isValidRange(const std::pair<int64_t, int64_t> &range,
                  const std::vector<std::pair<int64_t, int64_t>> &blockRanges) {

  
  if (range.first > range.second)
    return false;

  for (const auto &blockRange : blockRanges) {
    if (range.first >= blockRange.first && range.second <= blockRange.second) {
      return true;
    }
  }

  return false;
}

std::string coordinateToString(const tupleCoord_t &coord) {
  std::stringstream ss;
  ss << coord.blockId << ":" << coord.nucPos;
  if (coord.nucGapPos >= 0) {
    ss << ":" << coord.nucGapPos;
  }
  return ss.str();
}

tupleCoord_t stringToCoordinate(const std::string &str) {
  std::stringstream ss(str);
  std::string token;
  std::vector<int> parts;

  while (std::getline(ss, token, ':')) {
    parts.push_back(std::stoi(token));
  }

  if (parts.size() < 2 || parts.size() > 3) {
    throw std::runtime_error("Invalid coordinate string format");
  }

  tupleCoord_t coord;
  coord.blockId = parts[0];
  coord.nucPos = parts[1];
  coord.nucGapPos = parts.size() == 3 ? parts[2] : -1;

  return coord;
}

bool CoordinateTraverser::isInverted() const noexcept {
  return current->blockId >= 0 &&
         current->blockId < coordManager->getBlockStrand().size() &&
         !coordManager->getBlockStrand().at(current->blockId).first;
}

char CoordinateTraverser::getChar(bool complement) const noexcept {
  if (done || !current)
    return 'x'; // done

  // Check if block exists - if not, return gap character
  if (!coordManager->getBlockExists().at(current->blockId).first) {
    return '-';
  }

  // Get the block strand orientation
  bool isForwardStrand = !isBlockInverted(current->blockId);

  // Get character from sequence
  char c;
  if (current->nucGapPos == -1) {
    c = coordManager->getSequence()
            .at(current->blockId)
            .first[current->nucPos]
            .first;
  } else {
    c = coordManager->getSequence()
            .at(current->blockId)
            .first[current->nucPos]
            .second[current->nucGapPos];
  }

  // Early return for gap or end marker
  if (c == 'x' || c == '-') {
    return c;
  }

  // Apply complementation for inverted blocks
  if (complement && !isForwardStrand) {
    return seq_utils::getComplementCharacter(c);
  }

  return c;
}

void CoordinateTraverser::moveToNextBlock() noexcept {
  if (done)
    return;

  int currentBlock = current->blockId;
  while (!done && current->blockId == currentBlock) {
    next();
  }
}

void CoordinateTraverser::moveToPreviousBlock() noexcept {
  if (done)
    return;

  int currentBlock = current->blockId;
  while (!done && current->blockId == currentBlock) {
    prev();
  }
}

bool CoordinateTraverser::isBlockInverted(int blockId) const noexcept {
  return blockId >= 0 && blockId < coordManager->getBlockStrand().size() &&
         !coordManager->getBlockStrand().at(blockId).first;
}

bool CoordinateTraverser::isValidNucleotide(char c) const noexcept {
  c = std::toupper(c);
  return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

bool CoordinateTraverser::validateState() const noexcept {
  if (!coordManager)
    return false;
  if (!current)
    return false;

  // Validate current pointer
  if (!current->isValidBasic()) {
    err("Traverser validation failed: invalid current coordinate");
    return false;
  }

  // Validate scalar is in range
  if (current->scalar < 0 || current->scalar >= coordManager->size()) {
    err("Traverser validation failed: current coordinate scalar {} out of "
        "range [0,{}]",
        current->scalar, coordManager->size() - 1);
    return false;
  }

  // Validate current coordinate is in the manager's vector
  const tupleCoord_t *stored = coordManager->getTupleCoord(current->scalar);
  if (!stored || stored != current) {
    err("Traverser validation failed: current coordinate not found in "
        "manager's storage");
    return false;
  }

  return true;
}

void CoordinateTraverser::skipGaps() noexcept {
  if (done || !current)
    return;

  // If the current position is not a gap, nothing to do
  char c = getChar();
  if (c != '-')
    return;

  // Store the current scalar to detect if we're not making progress
  int64_t starting_scalar = current->scalar;
  int starting_block_id = current->blockId;

  // Check if we're in an inverted block
  bool in_inverted_block = isInverted();
  
  debug_msg("skipGaps: Starting at scalar {} in {} block", 
       current->scalar, in_inverted_block ? "inverted" : "normal");

  // Use the gap map to find the next non-gap position
  const tupleCoord_t *next_non_gap = in_inverted_block ?
      coordManager->getPreviousNonGapPosition(current->scalar) :
      coordManager->getNextNonGapPosition(current->scalar);

  // If we found a valid next position, jump to it
  if (next_non_gap) {
    // Safety check - make sure the position is valid
    if (next_non_gap->scalar < 0 || next_non_gap->scalar >= coordManager->size()) {
      err("skipGaps: next non-gap position {} is out of bounds [0,{}]",
          next_non_gap->scalar, coordManager->size() - 1);
      done = true;
      return;
    }
    
    // Safety check - make sure we're actually moving in the right direction
    if ((in_inverted_block && next_non_gap->scalar >= starting_scalar) ||
        (!in_inverted_block && next_non_gap->scalar <= starting_scalar)) {
      err("skipGaps: detected invalid gap map entry - next position ({}) is "
          "not in correct direction from current position ({})",
          next_non_gap->scalar, starting_scalar);
      done = true;
      return;
    }

    // Check if we're crossing a block boundary
    if (next_non_gap->blockId != starting_block_id) {
      debug_msg("skipGaps: Crossing block boundary from {} to {}", 
           starting_block_id, next_non_gap->blockId);
            
      // Check if the target block exists
      if (!coordManager->getBlockExists().at(next_non_gap->blockId).first) {
        // Need to skip further to the next live block
        debug_msg("skipGaps: Target block {} is dead, finding next live block");
        
        // Find the next live block in traversal order
        int nextLiveBlockId = -1;
        if (in_inverted_block) {
          // When inverted, we need to go to lower block IDs
          for (int i = next_non_gap->blockId; i >= 0; --i) {
            if (coordManager->getBlockExists().at(i).first) {
              nextLiveBlockId = i;
              break;
            }
          }
        } else {
          // When normal, we go to higher block IDs
          for (int i = next_non_gap->blockId; i < coordManager->getNumBlocks(); ++i) {
            if (coordManager->getBlockExists().at(i).first) {
              nextLiveBlockId = i;
              break;
            }
          }
        }
        
        // If we found a live block, get its first/last coordinate based on its orientation
        if (nextLiveBlockId >= 0) {
          bool nextBlockInverted = !coordManager->getBlockStrand().at(nextLiveBlockId).first;
          
          // Check for traversal direction consistency when changing blocks
          bool directionChange = in_inverted_block != nextBlockInverted;
          if (directionChange) {
            debug_msg("skipGaps: Direction change detected when moving from block {} ({}) to {} ({})",
                  starting_block_id, in_inverted_block ? "inverted" : "normal",
                  nextLiveBlockId, nextBlockInverted ? "inverted" : "normal");
          }
          
          next_non_gap = nextBlockInverted ? 
              coordManager->getLastCoordinateInBlock(nextLiveBlockId) :
              coordManager->getFirstCoordinateInBlock(nextLiveBlockId);
              
          if (!next_non_gap) {
            err("skipGaps: Failed to get valid coordinate for next live block {}", nextLiveBlockId);
            done = true;
            return;
          }
        } else {
          // No more live blocks in traversal direction
          debug_msg("skipGaps: No more live blocks in traversal direction, marking as done");
          done = true;
          return;
        }
      }
    }

    current = next_non_gap;
    debug_msg("skipGaps: Jumped to scalar {} in block {}", 
         current->scalar, current->blockId);

    // Check if we've reached the end coordinate
    if (end && current == end) {
      done = true;
    }
  } else {
    // If there's no valid next position, we're done
    debug_msg("skipGaps: No next non-gap position found, traversal complete");
    done = true;
  }
}

void CoordinateTraverser::skipGapsBackward() noexcept {
  if (done || !current)
    return;

  // If the current position is not a gap, nothing to do
  char c = getChar();
  if (c != '-')
    return;

  // Store the current information to detect issues
  int64_t starting_scalar = current->scalar;
  int starting_block_id = current->blockId;
  
  // Check if we're in an inverted block
  bool in_inverted_block = isInverted();
  
  debug_msg("skipGapsBackward: Starting at scalar {} in {} block", 
       current->scalar, in_inverted_block ? "inverted" : "normal");
  
  // Use the gap map to find the previous non-gap position - direction is reversed in inverted blocks
  const tupleCoord_t *prev_non_gap = in_inverted_block ?
      coordManager->getNextNonGapPosition(current->scalar) :
      coordManager->getPreviousNonGapPosition(current->scalar);

  // If we found a valid previous position, jump to it
  if (prev_non_gap) {
    // Safety check - make sure the position is valid
    if (prev_non_gap->scalar < 0 || prev_non_gap->scalar >= coordManager->size()) {
      err("skipGapsBackward: previous non-gap position {} is out of bounds [0,{}]",
          prev_non_gap->scalar, coordManager->size() - 1);
      done = true;
      return;
    }
    
    // Safety check - make sure we're moving in the expected direction
    if ((in_inverted_block && prev_non_gap->scalar <= starting_scalar) ||
        (!in_inverted_block && prev_non_gap->scalar >= starting_scalar)) {
      err("skipGapsBackward: detected invalid movement - prev position ({}) is "
          "not in correct direction from current position ({})",
          prev_non_gap->scalar, starting_scalar);
      done = true;
      return;
    }

    // Check if we're crossing a block boundary
    if (prev_non_gap->blockId != starting_block_id) {
      debug_msg("skipGapsBackward: Crossing block boundary from {} to {}", 
           starting_block_id, prev_non_gap->blockId);
            
      // Check if the target block exists
      if (!coordManager->getBlockExists().at(prev_non_gap->blockId).first) {
        // Need to skip further to the previous live block
        debug_msg("skipGapsBackward: Target block {} is dead, finding previous live block");
        
        // Find the previous live block in traversal order
        int prevLiveBlockId = -1;
        if (in_inverted_block) {
          // When inverted, going backward means higher block IDs
          for (int i = prev_non_gap->blockId; i < coordManager->getNumBlocks(); ++i) {
            if (coordManager->getBlockExists().at(i).first) {
              prevLiveBlockId = i;
              break;
            }
          }
        } else {
          // When normal, going backward means lower block IDs
          for (int i = prev_non_gap->blockId; i >= 0; --i) {
            if (coordManager->getBlockExists().at(i).first) {
              prevLiveBlockId = i;
              break;
            }
          }
        }
        
        // If we found a live block, get its first/last coordinate based on its orientation
        if (prevLiveBlockId >= 0) {
          bool prevBlockInverted = !coordManager->getBlockStrand().at(prevLiveBlockId).first;
          
          // Check for traversal direction consistency when changing blocks
          bool directionChange = in_inverted_block != prevBlockInverted;
          if (directionChange) {
            debug_msg("skipGapsBackward: Direction change detected when moving from block {} ({}) to {} ({})",
                  starting_block_id, in_inverted_block ? "inverted" : "normal",
                  prevLiveBlockId, prevBlockInverted ? "inverted" : "normal");
          }
          
          prev_non_gap = prevBlockInverted ? 
              coordManager->getFirstCoordinateInBlock(prevLiveBlockId) :
              coordManager->getLastCoordinateInBlock(prevLiveBlockId);
              
          if (!prev_non_gap) {
            err("skipGapsBackward: Failed to get valid coordinate for previous live block {}", prevLiveBlockId);
            done = true;
            return;
          }
        } else {
          // No more live blocks in reverse traversal direction
          debug_msg("skipGapsBackward: No more live blocks in reverse direction, marking as done");
          done = true;
          return;
        }
      }
    }

    current = prev_non_gap;
    debug_msg("skipGapsBackward: Jumped to scalar {} in block {}", 
         current->scalar, current->blockId);

    // Check if we've reached the end coordinate
    if (end && current == end) {
      done = true;
    }
  } else {
    // If there's no valid previous position, we're done
    debug_msg("skipGapsBackward: No previous non-gap position found, traversal complete");
    done = true;
  }
}

bool CoordinateTraverser::skipToNthNonGap(int n, char *buffer) noexcept {
  if (done || !current || n <= 0)
    return false;

  int found_count = 0;

  // Main loop to collect valid nucleotides
  while (found_count < n && !done) {
    // Get current character
    char c = getChar();

    // If it's a valid nucleotide (not a gap), add to buffer
    if (c != '-' && c != 'x') {
      if (buffer) { // Only store if buffer is provided
        buffer[found_count] = c;
      }
      found_count++;

      // If we've found enough, we're done
      if (found_count >= n)
        break;
    }

    // Move to next position efficiently by skipping gaps
    if (c == '-') {
      skipGaps(); // Use gap map to efficiently skip runs of gaps - already handles inversions now
    } else {
      next(); // Move to next position in traversal order (handles inversions via pointer structure)
    }

    // Check if we're done
    if (done && found_count < n)
      return false;
  }
  return true;
}

void CoordinateTraverser::next() noexcept {
  if (done)
    return;

  const tupleCoord_t *next = current->next;
  if (!next) {
    debug_msg("next: No next coordinate, traversal complete");
    done = true;
    return;
  }

  // Check if we're crossing a block boundary
  if (next->blockId != current->blockId) {
    debug_msg("next: Crossing block boundary from {} to {}", current->blockId, next->blockId);

    // Check if the next block is off - if so, we can skip it
    if (!coordManager->getBlockExists()[next->blockId].first) {
      debug_msg("next: Detected off block {}, will attempt to skip the run");
      
      // Save current position
      const tupleCoord_t *savedCurrent = current;
      
      // Try to skip the run of off blocks
      skipOffBlockRun();
      
      // If skipOffBlockRun set done=true or didn't change current, we're done
      if (done || current == savedCurrent) {
        return;
      }
      
      // If we get here, we've successfully skipped to the next live block
      debug_msg("next: Skipped off block run to block {}", current->blockId);
      return;
    } else {
      // Check if we're about to enter a run of blocks with the same orientation
      bool nextBlockInverted = !coordManager->getBlockStrand().at(next->blockId).first;
      
      // Detect potential run of blocks with same orientation
      auto [runStart, runEnd] = coordManager->findBlockRunWithSameOrientation(next->blockId);
      
      if (runStart != -1 && runEnd > runStart && runEnd - runStart > 1) {
        // We have at least 3 consecutive blocks with same orientation
        debug_msg("next: Detected a run of {} blocks with same orientation ({})",
             runEnd - runStart + 1, nextBlockInverted ? "inverted" : "forward");
             
        // Continue with normal traversal - to optimize this path, the caller can use traverseOrientationRun
      }
    }
  }

  // Validate next coordinate is in the manager's vector
  const tupleCoord_t *stored = coordManager->getTupleCoord(next->scalar);
  if (!stored || stored != next) {
    err("Traverser next: next coordinate not found in manager's storage");
    done = true;
    return;
  }

  // Check if we've reached the end coordinate
  if (end && next == end) {
    current = next;
    done = true;
    return;
  }
  
  current = next;
  debug_msg("next: Moved to scalar {} in block {}", current->scalar, current->blockId);
}

// New method to skip over a run of consecutive off blocks
void CoordinateTraverser::skipOffBlockRun() noexcept {
  if (done || !current)
    return;
    
  // Must be at the boundary of an off block
  int32_t currentBlockId = current->blockId;
  int32_t nextBlockId = -1;
  
  if (current->next) {
    nextBlockId = current->next->blockId;
  } else {
    // No next position, can't skip
    return;
  }
  
  // If next block is not off, nothing to do
  if (coordManager->getBlockExists().at(nextBlockId).first) {
    return;
  }
  
  // Find the run of off blocks
  auto [runStart, runEnd] = coordManager->findOffBlockRun(nextBlockId);
  
  if (runStart == -1 || runEnd == -1) {
    // No valid run found
    return;
  }
  
  // Find the next live block after the run
  int32_t nextLiveBlockId = coordManager->findNextLiveBlockInDirection(runEnd, true);
  
  if (nextLiveBlockId == -1) {
    // No more live blocks, we're done
    done = true;
    return;
  }
  
  // Get the first/last coordinate of the next live block based on its orientation
  bool nextBlockInverted = !coordManager->getBlockStrand().at(nextLiveBlockId).first;
  const tupleCoord_t *next = nextBlockInverted ? 
      coordManager->getLastCoordinateInBlock(nextLiveBlockId) :
      coordManager->getFirstCoordinateInBlock(nextLiveBlockId);
      
  if (!next) {
    err("skipOffBlockRun: Failed to get coordinate for next live block {}", nextLiveBlockId);
    done = true;
    return;
  }
  
  // Jump to the new position
  current = next;
  
  // Check if we've reached the end
  if (end && current == end) {
    done = true;
  }
}

// Define the templated implementation for traverseOrientationRun
template<typename SeedCallback>
void CoordinateTraverser::traverseOrientationRunImpl(SeedCallback callback, bool extractSeeds) noexcept {
  if (!current || done) {
    warn("traverseOrientationRun: Traverser already at end or invalid state");
    return;
  }

  // Save current position to restore later if needed
  const tupleCoord_t* startCoord = current;
  
  // Determine if this is an inverted run
  bool isInverted = !coordManager->getBlockStrand().at(current->blockId).first;
  
  // Get the end of the run (where orientation changes or block ends)
  auto blockRunRange = coordManager->findBlockRunWithSameOrientation(current->blockId);
  int32_t startBlockId = blockRunRange.first;
  int32_t endBlockId = blockRunRange.second;
  
  // Get scalar range for this run
  auto scalarRange = coordManager->getBlockRunScalarRange(startBlockId, endBlockId);
  int64_t startScalar = scalarRange.first;
  int64_t endScalar = scalarRange.second;
  
  if (startScalar < 0 || endScalar < 0 || startScalar > endScalar) {
    warn("traverseOrientationRun: Invalid scalar range ({}, {})", startScalar, endScalar);
    return;
  }
  
  // Implementation depends on what we're extracting
  if (extractSeeds) {
    // Common traversal logic for both inverted and forward runs
    const tupleCoord_t *endCoord = coordManager->getTupleCoord(endScalar);
    debug_msg("traverseOrientationRun: Traversing {} run from scalar {} to {}", 
          isInverted ? "inverted" : "forward", current->scalar, endScalar);
    
    // Buffer for collecting seeds if needed (for batch processing)
    constexpr size_t MAX_BUFFER_SIZE = 4096;
    char nucleotideBuffer[MAX_BUFFER_SIZE] = {0};
    size_t bufferPos = 0;
    std::vector<const tupleCoord_t*> coordBuffer;
    coordBuffer.reserve(MAX_BUFFER_SIZE);
    
    // Determines if we're in batch mode or callback-per-nucleotide mode
    bool batchMode = false;
    
    // Use a single pass traversal for efficiency
    while (current && !done) {
      // Check if we've reached the end of the run
      if (current->scalar > endScalar) {
        debug_msg("traverseOrientationRun: Reached end of run at scalar {}", current->scalar);
        break;
      }
      
      // Extract the nucleotide at this position
      char c = getChar(true);  // Get character with proper complementation
      
      // For seed extraction, we need to check if this is a valid nucleotide
      if (isValidNucleotide(c)) {
        if (batchMode) {
          // Store in buffer for batch processing
          if (bufferPos < MAX_BUFFER_SIZE - 1) {
            nucleotideBuffer[bufferPos] = c;
            coordBuffer.push_back(current);
            bufferPos++;
          }
        } else {
          // Call the callback function with the character and coordinate
          debug_msg("traverseOrientationRun: Extracting seed at position {} ({})", current->scalar, c);
          // Always call with consistent parameters - use a direct match to the callback signature
          callback(c, current, 1);
        }
      }
      
      // Move to the next position
      next();
      if (done) break;
    }
    
    // If in batch mode, process the buffered data
    if (batchMode && bufferPos > 0) {
      // Only call batch-style callback if there's actual data and the callback supports it
      // This requires checking the callback function signature for compatibility
      try {
        // Try to call the callback with the buffer data
        // This is a simplified approach - we'll need to handle this better later
        if (bufferPos > 0) {
          for (size_t i = 0; i < bufferPos; i++) {
            callback(nucleotideBuffer[i], coordBuffer[i], 1);
          }
        }
      } catch (...) {
        // If calling the callback fails, we'll just log an error and continue
        warn("traverseOrientationRun: Failed to call batch callback");
      }
    }
  } else {
    // Just traversal without extraction - we can optimize by jumping to the end
    // of the run directly if no processing is needed for each position
    debug_msg("traverseOrientationRun: Fast traversal to end of run (scalar {})", endScalar);
    
    const tupleCoord_t *targetPos = coordManager->getTupleCoord(endScalar);
    if (targetPos) {
      current = targetPos;
      if (end && current == end) {
        done = true;
      }
    } else {
      // If we can't find the target position, traverse normally
      while (current && !done && current->scalar <= endScalar) {
        next();
      }
    }
  }
}

// Explicit instantiations for common callback types
// This allows the compiler to generate the template code for these specific types
// making them available to other compilation units without including the implementation details
template void CoordinateTraverser::traverseOrientationRunImpl<std::function<void(char, const tupleCoord_t*, size_t)>>(
    std::function<void(char, const tupleCoord_t*, size_t)> callback, bool extractSeeds) noexcept;

// Add wrapper function template for traverseOrientationRun to use the implementation
template<typename SeedCallback>
void CoordinateTraverser::traverseOrientationRun(SeedCallback callback, bool extractSeeds) noexcept {
  traverseOrientationRunImpl(callback, extractSeeds);
}

// Explicit instantiations for the wrapper function
template void CoordinateTraverser::traverseOrientationRun<std::function<void(char, const tupleCoord_t*, size_t)>>(
    std::function<void(char, const tupleCoord_t*, size_t)> callback, bool extractSeeds) noexcept;

// No-op version for backward compatibility
void CoordinateTraverser::traverseOrientationRun(bool extractSeeds) noexcept {
  auto noOpCallback = [](char c, const tupleCoord_t* coord, size_t) {};
  traverseOrientationRunImpl(noOpCallback, extractSeeds);
}

void CoordinateTraverser::nextSkipGaps() noexcept {
  if (done)
    return;

  // First move to the next position
  next();

  // Then skip any gaps
  if (!done) {
    char c = getChar();
    if (c == '-') {
      skipGaps();
    }
  }
}

void CoordinateTraverser::prev() noexcept {
  if (done)
    return;

  const tupleCoord_t *prev = current->prev;
  if (!prev) {
    debug_msg("prev: No previous coordinate, traversal complete");
    done = true;
    return;
  }

  // Detect circular reference
  if (prev == current || prev->prev == current) {
    err("Traverser prev: detected circular reference in coordinate chain at "
        "scalar {}",
        current->scalar);
    done = true;
    return;
  }

  // Validate prev pointer
  if (!prev->isValidBasic()) {
    err("Traverser prev: invalid prev coordinate");
    done = true;
    return;
  }

  // Validate scalar is in range
  if (prev->scalar < 0 || prev->scalar >= coordManager->size()) {
    err("Traverser prev: prev coordinate scalar {} out of range [0,{}]",
        prev->scalar, coordManager->size() - 1);
    done = true;
    return;
  }

  // Check if we're crossing a block boundary
  if (prev->blockId != current->blockId) {
    debug_msg("prev: Crossing block boundary from {} to {}", 
         current->blockId, prev->blockId);
         
    // Check if the previous block is alive
    if (!coordManager->getBlockExists().at(prev->blockId).first) {
      debug_msg("prev: Block {} is dead, looking for previous live block");
      
      // Find previous live block
      int prevLiveBlockId = -1;
      for (int i = prev->blockId; i >= 0; i--) {
        if (coordManager->getBlockExists().at(i).first) {
          prevLiveBlockId = i;
          break;
        }
      }
      
      if (prevLiveBlockId >= 0) {
        // Get the last coordinate of the previous live block based on its orientation
        bool prevBlockInverted = !coordManager->getBlockStrand().at(prevLiveBlockId).first;
        prev = prevBlockInverted ? 
            coordManager->getFirstCoordinateInBlock(prevLiveBlockId) : 
            coordManager->getLastCoordinateInBlock(prevLiveBlockId);
            
        if (!prev) {
          err("prev: Failed to get coordinate for previous live block {}", prevLiveBlockId);
          done = true;
          return;
        }
      } else {
        // No more live blocks behind
        debug_msg("prev: No more live blocks behind, traversal complete");
        done = true; 
        return;
      }
    }
  }

  // Validate prev coordinate is in the manager's vector
  const tupleCoord_t *stored = coordManager->getTupleCoord(prev->scalar);
  if (!stored || stored != prev) {
    err("Traverser prev: prev coordinate not found in manager's storage");
    done = true;
    return;
  }

  if (prev == end) {
    current = prev;
    done = true;
    return;
  }
  
  current = prev;
  debug_msg("prev: Moved to scalar {} in block {}", current->scalar, current->blockId);
}

void CoordinateTraverser::reset(const tupleCoord_t *start) noexcept {
  current = start;
  end = nullptr;
  done = (current == nullptr);
}

void CoordinateTraverser::reset(const tupleCoord_t *start, const tupleCoord_t *end) noexcept {
  current = start;
  this->end = end;
  done = (current == nullptr);
}

} // namespace coordinates

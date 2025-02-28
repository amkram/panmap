#include "coordinates.hpp"
#include "logging.hpp"
#include "seed_annotated_tree.hpp"
#include "timing.hpp"
#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace logging;

namespace coordinates {

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

  // msg("Coordinate validation complete");

  // Optimize the gap map by merging adjacent runs
  optimizeGapMap();
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

  TIME_FUNCTION;
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

  // Skip any gaps at the current position using the gap map
  if (manager.isGapPosition(coord->scalar)) {
    traverser.skipGapsBackward();
    if (traverser.isDone()) {
      return result;
    }
    result = traverser.getCurrent();
  }

  while (count < neededNongap && !traverser.isDone()) {
    const tupleCoord_t *pos = traverser.getCurrent();
    if (!pos) {
      err("expandUpstream: got null position from traverser");
      break;
    }

    // Validate block ID before access
    if (pos->blockId < 0 || pos->blockId >= blockExists.size()) {
      err("expandUpstream: invalid block ID {} (max: {})", pos->blockId,
          blockExists.size() - 1);
      break;
    }

    // Now safe to check block existence
    if (!blockExists[pos->blockId].first) {
      // Block is off, move to previous block
      if (pos->blockId == 0) {
        return result;
      }
      traverser.moveToPreviousBlock();
      continue;
    }

    // Stop if we've reached stop_coord
    if (stop_coord && pos == stop_coord) {
      break;
    }

    // Count non-gap positions
    char c = traverser.getChar();
    if (c == 'x') {
      return result;
    }
    if (c != '-') {
      count++;
      result = pos;
    } else {
      // Now we can use skipGapsBackward for more efficient traversal
      traverser.skipGapsBackward();
      continue;
    }

    traverser.prev();
  }

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
    throw std::runtime_error("expandDownstream: null coordinate");
  }

  auto &manager = traverser.getCoordManager();

  int count = 0;
  const tupleCoord_t *result = coord;

  traverser.reset(coord);

  // Skip any gaps at the current position using the gap map
  if (manager.isGapPosition(coord->scalar)) {
    traverser.skipGaps();
    if (traverser.isDone()) {
      return result;
    }
    result = traverser.getCurrent();
  }

  while (count < neededNongap && !traverser.isDone()) {
    const tupleCoord_t *pos = traverser.getCurrent();

    // Handle dead blocks
    if (!blockExists[pos->blockId].first) {
      if (pos->blockId >= traverser.getCoordManager().getNumBlocks() - 1) {
        return result;
      }
      traverser.moveToNextBlock();
      continue;
    }

    // Count non-gap positions
    char c = traverser.getChar();
    if (c == 'x') {
      return result;
    }
    if (c != '-') {
      count++;
      result = pos;
    } else {
      // Use the skip gaps optimization when going forward
      traverser.skipGaps();
      continue;
    }

    traverser.next();
  }

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

  TIME_FUNCTION;
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

  // Use the gap map to find the next non-gap position
  const tupleCoord_t *next_non_gap =
      coordManager->getNextNonGapPosition(current->scalar);

  // If we found a valid next position, jump to it
  if (next_non_gap) {
    // Safety check - make sure we're actually moving forward
    if (next_non_gap->scalar <= starting_scalar) {
      err("skipGaps: detected invalid gap map entry - next position ({}) is "
          "not after current position ({})",
          next_non_gap->scalar, starting_scalar);
      done = true;
      return;
    }

    current = next_non_gap;

    // Check if we've reached the end coordinate
    if (end && current == end) {
      done = true;
    }
  } else {
    // If there's no valid next position, we're done
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

  // Use the gap map to find the previous non-gap position
  const tupleCoord_t *prev_non_gap =
      coordManager->getPreviousNonGapPosition(current->scalar);

  // If we found a valid previous position, jump to it
  if (prev_non_gap) {
    current = prev_non_gap;

    // Check if we've reached the end coordinate (somewhat unlikely in backward
    // direction)
    if (end && current == end) {
      done = true;
    }
  } else {
    // If there's no valid previous position, we're done
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
      skipGaps(); // Use gap map to efficiently skip runs of gaps
    } else {
      next(); // Move to next position
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
    done = true;
    return;
  }

  // Detect circular reference - if next points back to current or itself
  if (next == current || next->next == current) {
    err("Traverser next: detected circular reference in coordinate chain at "
        "scalar {}",
        current->scalar);
    done = true;
    return;
  }

  // Validate next pointer
  if (!next->isValidBasic()) {
    err("Traverser next: invalid next coordinate");
    done = true;
    return;
  }

  // Validate scalar is in range
  if (next->scalar < 0 || next->scalar >= coordManager->size()) {
    err("Traverser next: next coordinate scalar {} out of range [0,{}]",
        next->scalar, coordManager->size() - 1);
    done = true;
    return;
  }

  // Validate next coordinate is in the manager's vector
  const tupleCoord_t *stored = coordManager->getTupleCoord(next->scalar);
  if (!stored || stored != next) {
    err("Traverser next: next coordinate not found in manager's storage");
    done = true;
    return;
  }

  if (next == end) {
    done = true;
    return;
  }
  current = next;
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

  // Validate prev coordinate is in the manager's vector
  const tupleCoord_t *stored = coordManager->getTupleCoord(prev->scalar);
  if (!stored || stored != prev) {
    err("Traverser prev: prev coordinate not found in manager's storage");
    done = true;
    return;
  }

  if (prev == end) {
    done = true;
    return;
  }
  current = prev;
}

void CoordinateTraverser::reset(const tupleCoord_t *start,
                                const tupleCoord_t *end) noexcept {
  if (!start) {
    current = nullptr;
    this->end = nullptr;
    done = true;
    return;
  }

  current = start;
  this->end = end;
  done = false;
}

void CoordinateTraverser::reset(const tupleCoord_t *start) noexcept {
  reset(start, nullptr);
}

} // namespace coordinates

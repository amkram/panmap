#pragma once

#include "mgsr.hpp"
#include "panmap_utils.hpp"

static void compareBruteForce(
  panmanUtils::Tree *T,
  panmanUtils::Node *node,
  const panmapUtils::BlockSequences& blockSequences,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& degapCoordIndex,
  const std::map<uint64_t, uint64_t>& regapCoordIndex
) {
  bool printCorrect = false;
  std::cout << "Checking " << node->identifier << " states with brute force..." << std::endl;

  // check sequence object
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequenceBruteForce;
  std::vector<bool> blockExistsBruteForce;
  std::vector<bool> blockStrandBruteForce;
  std::unordered_map<int, int> blockLengthsBruteForce;
  panmapUtils::getSequenceFromReference(T, sequenceBruteForce, blockExistsBruteForce, blockStrandBruteForce, blockLengthsBruteForce, node->identifier);


  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequenceDynamic = blockSequences.sequence;
  const std::vector<bool>& blockExistsDynamic = blockSequences.blockExists;
  const std::vector<bool>& blockStrandDynamic = blockSequences.blockStrand;
  if (sequenceDynamic.size() != sequenceBruteForce.size()) {
    std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic.size() << " != brute force " << sequenceBruteForce.size() << std::endl;
    std::exit(1);
  } else {
    if (printCorrect) std::cout << "\tIdentical sequence size... passed" << std::endl;
  }

  uint64_t localScalarCoordBruteForce = 0;
  uint64_t globalScalarCoord = 0;
  for (int blockId = 0; blockId < blockSequences.numBlocks(); blockId++) {
    if (!blockExistsDynamic[blockId]) {
      globalScalarCoord += blockLengthsBruteForce[blockId];
      continue;
    }

    if (globalScalarCoord != globalCoords.getBlockStartScalar(blockId)) {
      std::cerr << "Global scalar coord block " << blockId << " start mismatch: dynamic " << globalScalarCoord << " != brute force " << globalCoords.getBlockStartScalar(blockId) << std::endl;
      std::exit(1);
    } else {
      if (printCorrect) std::cout << "\tIdentical block " << blockId << " start scalar coord... passed" << std::endl;
    }

    if (blockSequences.blockExists[blockId] != blockExistsBruteForce[blockId]) {
      std::cerr << "Block " << blockId << " exists state mismatch: dynamic " << blockSequences.blockExists[blockId] << " != brute force " << blockExistsBruteForce[blockId] << std::endl;
      std::exit(1);
    } else {
      if (printCorrect) std::cout << "\tIdentical block " << blockId << " exists state... passed" << std::endl;
    }

    if (blockSequences.blockStrand[blockId] != blockStrandBruteForce[blockId]) {
      std::cerr << "Block " << blockId << " strand state mismatch: dynamic " << blockSequences.blockStrand[blockId] << " != brute force " << blockStrandBruteForce[blockId] << std::endl;
      std::exit(1);
    } else {
      if (printCorrect) std::cout << "\tIdentical block " << blockId << " strand state... passed" << std::endl;
    }

    for (int i = 0; i < sequenceDynamic[blockId].size(); i++) {      
      if (sequenceDynamic[blockId][i].second.size() != sequenceBruteForce[blockId][i].second.size()) {
        std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic[blockId][i].second.size() << " != brute force " << sequenceBruteForce[blockId][i].second.size() << std::endl;
        std::exit(1);
      } else {
        if (printCorrect) std::cout << "\tIdentical gap nuc size at (" << blockId << ", " << i << ")... passed" << std::endl;
      }

      for (int j = 0; j < sequenceDynamic[blockId][i].second.size(); j++) {
        if (sequenceDynamic[blockId][i].second[j] != sequenceBruteForce[blockId][i].second[j]) {
          std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", " << j << "): dynamic " << sequenceDynamic[blockId][i].second[j] << " != brute force " << sequenceBruteForce[blockId][i].second[j] << std::endl;
          std::exit(1);
        } else {
          if (printCorrect) std::cout << "\tIdentical gap nuc at (" << blockId << ", " << i << ", " << j << ")... passed" << std::endl;
        }

        if (sequenceDynamic[blockId][i].second[j] != '-') {
          if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
            std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") | " << globalScalarCoord << ": dynamic " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != brute force " << localScalarCoordBruteForce << std::endl;
            std::exit(1);
          } else {
            if (printCorrect) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed" << std::endl;
          }
          ++localScalarCoordBruteForce;
        }
        ++globalScalarCoord;
      }

      if (sequenceDynamic[blockId][i].first != sequenceBruteForce[blockId][i].first) {
        std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", -1): dynamic " << sequenceDynamic[blockId][i].first << " != brute force " << sequenceBruteForce[blockId][i].first << std::endl;
        std::exit(1);
      } else {
        if (printCorrect) std::cout << "\tIdentical main nuc at (" << blockId << ", " << i << ", -1)... passed" << std::endl;
      }

      if (sequenceDynamic[blockId][i].first != 'x') {
        if (sequenceDynamic[blockId][i].first != '-') {
          if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
            std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) | " << globalScalarCoord << ": dynamic " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != brute force " << localScalarCoordBruteForce << std::endl;
            std::exit(1);
          } else {
            if (printCorrect) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed" << std::endl;
          }
          ++localScalarCoordBruteForce;
        }
        ++globalScalarCoord;
      }
    }
    if (globalScalarCoord - 1 != globalCoords.getBlockEndScalar(blockId)) {
      std::cerr << "Global scalar coord block " << blockId << " end mismatch: dynamic " << globalScalarCoord - 1 << " != brute force " << globalCoords.getBlockEndScalar(blockId) << std::endl;
      std::exit(1);
    } else {
      if (printCorrect) std::cout << "\tIdentical block " << blockId << " end scalar coord... passed" << std::endl;
    }
  } 

  std::cout << "         " << node->identifier << " states passed brute force check" << std::endl;
  
}

static void applyMutations (
  panmanUtils::Node *node,
  panmapUtils::BlockSequences &blockSequences,
  std::unordered_set<uint64_t>& invertedBlocks,
  panmapUtils::GlobalCoords& globalCoords,
  std::vector<std::tuple<uint32_t, bool, bool, bool, bool>>& blockMutationRecord,
  std::vector<std::tuple<panmapUtils::Coordinate, char, char>>& nucMutationRecord,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapRunUpdates,
  std::vector<std::pair<uint64_t, bool>>& invertedBlocksBacktracks,
  const std::vector<bool>& oldBlockExists,
  const std::vector<bool>& oldBlockStrand
) {
  std::vector<bool>& blockExists = blockSequences.blockExists;
  std::vector<bool>& blockStrand = blockSequences.blockStrand;
  std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;

  // process block mutations
  for (const auto& blockMutation : node->blockMutation) {
    int32_t blockId = blockMutation.primaryBlockId;
    bool isInsertion = blockMutation.blockMutInfo;
    bool isInversion = blockMutation.inversion;
    bool oldExists = blockExists[blockId];
    bool oldStrand = blockStrand[blockId];
    if (isInsertion) {
      blockExists[blockId] = true;
      blockStrand[blockId] = !isInversion;
      if (isInversion) {
        std::cout << "\t\t" << node->identifier << " inserted inversion block " << blockId << std::endl;
      }
      if (!blockStrand[blockId]) {
        invertedBlocks.insert(blockId);
        invertedBlocksBacktracks.emplace_back(blockId, true);
      }
    } else if (isInversion) {
      blockStrand[blockId] = !blockStrand[blockId];
      std::cout << "\t\t" << node->identifier << " inverted block " << blockId << " from " << oldStrand << " to " << blockStrand[blockId] << std::endl;
      if (!blockStrand[blockId]) {
        invertedBlocks.insert(blockId);
        invertedBlocksBacktracks.emplace_back(blockId, true);
      } else {
        invertedBlocks.erase(blockId);
        invertedBlocksBacktracks.emplace_back(blockId, false);
      }
    } else {
      blockExists[blockId] = false;
      blockStrand[blockId] = true;
      if (!oldStrand) {
        invertedBlocks.erase(blockId);
        invertedBlocksBacktracks.emplace_back(blockId, false);
      }
    }
    blockMutationRecord.emplace_back(std::make_tuple(blockId, oldExists, oldStrand, blockExists[blockId], blockStrand[blockId]));
  }

  // process nuc mutations
  for (const auto& nucMutation : node->nucMutation) {
    int length = nucMutation.mutInfo >> 4;
    for (int i = 0; i < length; i++) {
      panmapUtils::Coordinate pos = panmapUtils::Coordinate(nucMutation, i);
      char oldNuc = blockSequences.getSequenceBase(pos);
      int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
      char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
      blockSequences.setSequenceBase(pos, newNuc);
      nucMutationRecord.emplace_back(std::make_tuple(pos, oldNuc, newNuc));

      if (oldBlockExists[pos.primaryBlockId] && blockExists[pos.primaryBlockId]) {
        int64_t scalarCoord = globalCoords.getScalarFromCoord(pos);
        if (oldNuc != '-' && newNuc == '-') {
          // nuc to gap
          if (!gapRunUpdates.empty() && gapRunUpdates.back().first == true && gapRunUpdates.back().second.second + 1 == scalarCoord) {
            ++(gapRunUpdates.back().second.second);
          }
          else {
            gapRunUpdates.emplace_back(true, std::make_pair(scalarCoord, scalarCoord)); 
          }
        } else if (oldNuc == '-' && newNuc != '-') {
          // gap to nuc
          if (!gapRunUpdates.empty() && gapRunUpdates.back().first == false && gapRunUpdates.back().second.second + 1 == scalarCoord) {
            ++(gapRunUpdates.back().second.second);
          } else {
            gapRunUpdates.emplace_back(false, std::make_pair(scalarCoord, scalarCoord));
          }
        }
      }
    }
  }
}

uint64_t mgsr::degapGlobal(const uint64_t& globalCoord, const std::map<uint64_t, uint64_t>& degapCoordsIndex) {
  auto coordIt = degapCoordsIndex.upper_bound(globalCoord);
  if (coordIt == degapCoordsIndex.begin()) {
      return 0;
  }
  return globalCoord - std::prev(coordIt)->second;
}

uint64_t mgsr::regapGlobal(const uint64_t& localCoord, const std::map<uint64_t, uint64_t>& regapCoordsIndex) {
  auto coordIt = regapCoordsIndex.upper_bound(localCoord);
  if (coordIt == regapCoordsIndex.begin()) {
      return 0;
  }
  return localCoord + std::prev(coordIt)->second;
}

void mgsr::updateGapMapStep(
  std::map<uint64_t, uint64_t>& gapMap,
  const std::pair<bool, std::pair<uint64_t, uint64_t>>& update,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates,
  bool recordGapMapUpdates
) {
  bool toGap = update.first;
  uint64_t start = update.second.first;
  uint64_t end = update.second.second;

  auto rightIt = gapMap.upper_bound(start);
  auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

  bool rightItExists = rightIt != gapMap.end();
  bool leftItExists = leftIt != gapMap.end();

  if (toGap) {
    // add gap range
    if (gapMap.empty()) {
      gapMap[start] = end;
      backtrack.emplace_back(true, std::make_pair(start, end));
      if (recordGapMapUpdates) {
        gapMapUpdates.emplace_back(false, std::make_pair(start, end));
      }
      return;
    }
    
    decltype(rightIt) curIt;

    // curIt starts outside of any range
    if (!leftItExists || (!rightItExists && start > leftIt->second) || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
      if (leftItExists && start == leftIt->second + 1) {
        // 1 base after left range and merge with left
        curIt = leftIt;
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        curIt->second = end;
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        }
      } else {
        // insert new range
        auto tmpIt = gapMap.emplace(start, end);
        curIt = tmpIt.first;
        backtrack.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        }
      }
    } else {
      curIt = leftIt;
      if (end <= curIt->second) {
        return;
      }
      backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
      curIt->second = end;
      if (recordGapMapUpdates) {
        gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
      }
    }

    auto nextIt = std::next(curIt);
    while (true) {
      if (nextIt == gapMap.end()) {
        break;
      }

      if (nextIt->second <= curIt->second) {
        auto tmpIt = nextIt;
        nextIt = std::next(nextIt);
        backtrack.emplace_back(false, std::make_pair(tmpIt->first, tmpIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(true, std::make_pair(tmpIt->first, tmpIt->second));
        }
        gapMap.erase(tmpIt);
      } else if (nextIt->first <= end + 1) {
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        curIt->second = nextIt->second;
        backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          gapMapUpdates.emplace_back(true, std::make_pair(nextIt->first, nextIt->second));
        }
        gapMap.erase(nextIt);
        break;
      } else {
        break;
      }
    }
  } else {
    // remove gap range
    if (gapMap.empty() || (!leftItExists && end < rightIt->first) || (!rightItExists && start > leftIt->second)) {
      return;
    }

    decltype(rightIt) curIt;
    decltype(rightIt) nextIt;
    if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
      // curIt starts outside of any range
      curIt = rightIt;

      if (end < curIt->first) {
        return;
      }

      // ends within the curIt range
      if (end <= curIt->second) {
        if (end == curIt->second) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else {
          gapMap[end+1] = curIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        }
        return;
      } else {
        nextIt = std::next(curIt);
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        }
        gapMap.erase(curIt);
      }
      
    } else {
      // curIt starts inside of a range
      curIt = leftIt;
      
      if (end <= curIt->second) {
        // contained in the curIt range
        if (start == curIt->first && end == curIt->second) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else if (start == curIt->first) {
          gapMap[end + 1] = curIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else if (end == curIt->second) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          curIt->second = start - 1;
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          }
        } else {
          gapMap[end + 1] = curIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
            gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, start-1));
          }
          curIt->second = start - 1;
        }
        return;
      } else {
        if (start == curIt->first) {
          nextIt = std::next(curIt);
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          curIt->second = start - 1;
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          }
          nextIt = std::next(curIt);
        }
      }
    }

    
    while (true) {
      if (nextIt == gapMap.end()) {
        break;
      }

      if (nextIt->first > end) {
        break;
      } else if (nextIt->second <= end) {
        auto tmpIt = nextIt;
        nextIt = std::next(nextIt);
        backtrack.emplace_back(false, std::make_pair(tmpIt->first, tmpIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(true, std::make_pair(tmpIt->first, tmpIt->second));
        }
        gapMap.erase(tmpIt);
      } else {
        gapMap[end + 1] = nextIt->second;
        backtrack.emplace_back(true, std::make_pair(end+1, nextIt->second));
        backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(end+1, nextIt->second));
          gapMapUpdates.emplace_back(true, std::make_pair(nextIt->first, nextIt->second));
        }
        gapMap.erase(nextIt);
        break;
      }
    }
  }
}

void mgsr::updateGapMap(
  std::map<uint64_t, uint64_t>& gapMap,
  const std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& updates,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates
) {
  for (const auto& update : updates) {
    updateGapMapStep(gapMap, update, backtrack, gapMapUpdates, true);
  }
}

void mgsr::invertGapMap(
  std::map<uint64_t, uint64_t>& gapMap,
  const std::pair<uint64_t, uint64_t>& invertRange,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates
) {
  const auto& [start, end] = invertRange;

  auto rightIt = gapMap.upper_bound(start);
  auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

  bool rightItExists = rightIt != gapMap.end();
  bool leftItExists = leftIt != gapMap.end();

  // completely inside or outside a gap range -> do nothing
  if (
    gapMap.empty() || // empty gap map
    (!leftItExists && end < rightIt->first) || // completely left of first gap range
    (!rightItExists && start > leftIt->second) // completely right of last gap range
    // (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first) || // completely between two gap ranges
    // (leftItExists && start >= leftIt->first && end <= leftIt->second) // completely inside a gap range
  ) {
    return;
  }
  
  //                    gaps            beg      end
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> blockRuns;
  if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
    // start outside of a range
    auto curIt = rightIt;
    blockRuns.emplace_back(false, std::make_pair(start, curIt->first - 1));
    if (end <= curIt->second) {
      blockRuns.emplace_back(true, std::make_pair(blockRuns.back().second.second + 1, end));
    } else {
      blockRuns.emplace_back(true, std::make_pair(blockRuns.back().second.second + 1, curIt->second));
      curIt = std::next(curIt);
      while (true) {
        if (curIt == gapMap.end()) {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
          break;
        } else if (end < curIt->first) {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
          break;
        } else if (end > curIt->second) {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
          blockRuns.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          curIt = std::next(curIt);
        } else {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
          blockRuns.emplace_back(true, std::make_pair(curIt->first, end));
          break;
        }
      }
    }
  } else {
    // start inside of a range
    auto curIt = leftIt;
    blockRuns.emplace_back(true, std::make_pair(start, curIt->second));
    curIt = std::next(curIt);
    while (true) {
      if (curIt == gapMap.end()) {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
        break;
      } else if (end < curIt->first) {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
        break;
      } else if (end > curIt->second) {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
        blockRuns.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        curIt = std::next(curIt);
      } else {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
        blockRuns.emplace_back(true, std::make_pair(curIt->first, end));
        break;
      }
    }
  }

  int64_t curBeg = blockRuns.front().second.first;
  for (auto it = blockRuns.rbegin(); it != blockRuns.rend(); ++it) {
    int64_t curEnd = curBeg + (it->second.second - it->second.first);
    updateGapMapStep(gapMap, {it->first, {curBeg, curEnd}}, backtrack, gapMapUpdates, false);
    curBeg = curEnd + 1;
  }
}

void mgsr::makeCoordIndex(
  std::map<uint64_t, uint64_t>& degapCoordIndex,
  std::map<uint64_t, uint64_t>& regapCoordIndex,
  const std::map<uint64_t, uint64_t>& gapMap,
  const panmapUtils::GlobalCoords& globalCoords
) {
  uint64_t totalGapSize = 0;
  if (gapMap.empty() || gapMap.begin()->first > 0) {
    degapCoordIndex[0] == totalGapSize;
    regapCoordIndex[0] == totalGapSize;
  }
  for (auto &gap : gapMap) {
    uint64_t gapStart = gap.first;
    uint64_t gapEnd = gap.second;
    uint64_t gapSize = gapEnd - gapStart + 1;
    if (gapEnd == (uint64_t)globalCoords.lastScalarCoord) break;
    totalGapSize += gapSize;
    degapCoordIndex[gapEnd+1] = totalGapSize;
    regapCoordIndex[gapEnd+1-totalGapSize] = totalGapSize;
  }
}

void mgsr::mgsrIndexBuilder::buildIndexHelper(
  panmanUtils::Node *node,
  panmapUtils::BlockSequences &blockSequences,
  std::vector<bool> &blockExistsDelayed,
  std::vector<bool> &blockStrandDelayed,
  panmapUtils::GlobalCoords &globalCoords,
  std::map<uint64_t, uint64_t> &gapMap,
  std::unordered_set<uint64_t> &invertedBlocks
) {
  // record old and new block and nuc states for backtracking and calculating gaps
  std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> blockMutationRecord;
  std::vector<std::tuple<panmapUtils::Coordinate, char, char>> nucMutationRecord;

  // for building gap map
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapUpdates;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBlocksBacktracks;
  std::vector<std::pair<uint64_t, bool>> invertedBlocksBacktracks;

  applyMutations(node, blockSequences, invertedBlocks, globalCoords, blockMutationRecord, nucMutationRecord, gapRunUpdates, invertedBlocksBacktracks, blockExistsDelayed, blockStrandDelayed);
  
  updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      // on to off -> block range to all gaps
      uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
      uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
      updateGapMapStep(gapMap, {true, {beg, end}}, gapRunBacktracks, gapMapUpdates, true);
    } else if (!oldExists && newExists) {
      // off to on -> recompute across entire block
      panmapUtils::Coordinate coord = globalCoords.getBlockStartCoord(blockId);
      panmapUtils::Coordinate end = globalCoords.getBlockEndCoord(blockId);
      std::pair<int64_t, int64_t> curNucRange = {-1, -1};
      std::vector<std::pair<int64_t, int64_t>> nucRanges;
      while (true) {
        char nuc = blockSequences.getSequenceBase(coord);
        nuc = nuc == 'x' ? '-' : nuc;
        int64_t scalar = globalCoords.getScalarFromCoord(coord);
        if (nuc != '-') {
          if (curNucRange.first != -1 && curNucRange.second + 1 == scalar) {
            ++curNucRange.second;
          } else {
            if (curNucRange.first != -1) {
              if (node->identifier == "MG027860.1" && blockId == 1459) {
                std::cout << "At node " << node->identifier << " block " << blockId << " updating gap map from " << curNucRange.first << " to " << curNucRange.second << std::endl;
              }
              updateGapMapStep(gapMap, {false, {(uint64_t)curNucRange.first, (uint64_t)curNucRange.second}}, gapRunBacktracks, gapMapUpdates, true);
            }
            curNucRange = {scalar, scalar};
          }
        }

        if (coord == end) break;
        coord = globalCoords.stepRight(coord);
      }
      if (curNucRange.first != -1) {
        if (node->identifier == "MG027860.1" && blockId == 1459) {
          std::cout << "At node " << node->identifier << " block " << blockId << " updating gap map from " << curNucRange.first << " to " << curNucRange.second << std::endl;
        }
        updateGapMapStep(gapMap, {false, {(uint64_t)curNucRange.first, (uint64_t)curNucRange.second}}, gapRunBacktracks, gapMapUpdates, true);
      }
    }
  }

  for (const auto& blockId : invertedBlocks) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
    if (node->identifier == "MG027860.1" && blockId == 1459) {
      std::cout << "At node " << node->identifier << " block " << blockId << " inverting gap map from " << beg << " to " << end << std::endl;
    }
    invertGapMap(gapMap, {beg, end}, gapRunBlocksBacktracks, gapMapUpdates);
  }

  std::map<uint64_t, uint64_t> degapCoordIndex;
  std::map<uint64_t, uint64_t> regapCoordIndex;
  makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, globalCoords);

  for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  // update delayed block states
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExistsDelayed[blockId] = newExists;
    blockStrandDelayed[blockId] = newStrand;
  }

  // compare with brute force for debugging
  compareBruteForce(T, node, blockSequences, globalCoords, degapCoordIndex, regapCoordIndex);

  // update delayed inverted blocks
  for (panmanUtils::Node *child : node->children) {
    buildIndexHelper(child, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks);
  }

  // backtrack
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockSequences.blockExists[blockId] = oldExists;
    blockSequences.blockStrand[blockId] = oldStrand;
    blockExistsDelayed[blockId] = oldExists;
    blockStrandDelayed[blockId] = oldStrand;
  }

  for (const auto& [coord, oldNuc, newNuc] : nucMutationRecord) {
    blockSequences.setSequenceBase(coord, oldNuc);
  }

  for (const auto& [del, blockId] : invertedBlocksBacktracks) {
    if (del) {
      invertedBlocks.erase(blockId);
    } else {
      invertedBlocks.insert(blockId);
    }
  }

  for (auto it = gapRunBacktracks.rbegin(); it != gapRunBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }
}

void mgsr::mgsrIndexBuilder::buildIndex() {
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences.sequence);

  std::vector<bool> blockExistsDelayed = blockSequences.blockExists;
  std::vector<bool> blockStrandDelayed = blockSequences.blockStrand;

  std::map<uint64_t, uint64_t> gapMap{{0, globalCoords.lastScalarCoord}};
  std::unordered_set<uint64_t> invertedBlocks;

  buildIndexHelper(T->root, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks);
}


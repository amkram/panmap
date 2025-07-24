#pragma once

#include "mgsr.hpp"
#include "panmap_utils.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"

static void compareBruteForce(
  panmanUtils::Tree *T,
  panmanUtils::Node *node,
  const panmapUtils::BlockSequences& blockSequences,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  const std::map<uint64_t, uint64_t>& degapCoordIndex,
  const std::map<uint64_t, uint64_t>& regapCoordIndex,
  const std::vector<std::optional<seeding::rsyncmer_t>>& refOnSyncmers,
  const std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers,
  int k, int s, int t, int l, bool open
) {
  bool printCorrect = false;
  std::string nodeToDebug = "KJ627733.1";
  std::cout << "Checking " << node->identifier << " states with brute force..." << std::endl;

  // check sequence object
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequenceBruteForce;
  std::vector<bool> blockExistsBruteForce;
  std::vector<bool> blockStrandBruteForce;
  std::unordered_map<int, int> blockLengthsBruteForce;
  panmapUtils::getSequenceFromReference(T, sequenceBruteForce, blockExistsBruteForce, blockStrandBruteForce, blockLengthsBruteForce, node->identifier);

  std::string gappedSequenceBruteForce = panmapUtils::getStringFromSequence(sequenceBruteForce, blockLengthsBruteForce, blockExistsBruteForce, blockStrandBruteForce, true);
  std::vector<std::pair<uint64_t, uint64_t>> gapMapBruteForce;
  // checking gap map with brute force
  for (uint64_t i = 0; i < gappedSequenceBruteForce.size(); i++) {
    char nuc = gappedSequenceBruteForce[i];
    if (nuc == '-') {
      if (!gapMapBruteForce.empty() && gapMapBruteForce.back().second + 1 == i) {
        ++gapMapBruteForce.back().second;
      } else {
        gapMapBruteForce.emplace_back(i, i);
      }
    }
  }
  
  if (gapMapBruteForce.size() != gapMap.size()) {
    std::cout << "Gap map size mismatch: dynamic " << gapMap.size() << " != brute force " << gapMapBruteForce.size() << std::endl;
    std::cout << "Dynamic gap map:";
    for (const auto& [a, b] : gapMap) {
      std::cout << "(" << a << "," << b << ") ";
    }
    std::cout << std::endl;
    std::cout << "Brute force gap map:";
    for (const auto& [a, b] : gapMapBruteForce) {
      std::cout << "(" << a << "," << b << ") ";
    }
    std::cout << std::endl;
    std::exit(1);
  }

  size_t gapMapIndex = 0;
  for (const auto& [a, b] : gapMap) {
    if (a == gapMapBruteForce[gapMapIndex].first && b == gapMapBruteForce[gapMapIndex].second) {
      ++gapMapIndex;
    } else {
      std::cout << "Gap map mismatch at gap map index " << gapMapIndex << " dynamic " << a << " " << b << " != brute force " << gapMapBruteForce[gapMapIndex].first << " " << gapMapBruteForce[gapMapIndex].second << std::endl;
      std::cout << "Dynamic gap map:";
      for (const auto& [a, b] : gapMap) {
        std::cout << "(" << a << "," << b << ") ";
      }
      std::cout << std::endl;
      std::cout << "Brute force gap map:";
      for (const auto& [a, b] : gapMapBruteForce) {
        std::cout << "(" << a << "," << b << ") ";
      }
      std::cout << std::endl;
      std::exit(1);
    }
  }


  // check block sequence objects and coordinates
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequenceDynamic = blockSequences.sequence;
  const std::vector<bool>& blockExistsDynamic = blockSequences.blockExists;
  const std::vector<bool>& blockStrandDynamic = blockSequences.blockStrand;
  if (sequenceDynamic.size() != sequenceBruteForce.size()) {
    std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic.size() << " != brute force " << sequenceBruteForce.size() << std::endl;
    std::exit(1);
  } else {
    if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical sequence size... passed: " << sequenceDynamic.size() << " == " << sequenceBruteForce.size() << std::endl;
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
      if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " start scalar coord... passed: " << globalScalarCoord << " == " << globalCoords.getBlockStartScalar(blockId) << std::endl;
    }

    if (blockSequences.blockExists[blockId] != blockExistsBruteForce[blockId]) {
      std::cerr << "Block " << blockId << " exists state mismatch: dynamic " << blockSequences.blockExists[blockId] << " != brute force " << blockExistsBruteForce[blockId] << std::endl;
      std::exit(1);
    } else {
      if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " exists state... passed: " << blockSequences.blockExists[blockId] << " == " << blockExistsBruteForce[blockId] << std::endl;
    }

    if (blockSequences.blockStrand[blockId] != blockStrandBruteForce[blockId]) {
      std::cerr << "Block " << blockId << " strand state mismatch: dynamic " << blockSequences.blockStrand[blockId] << " != brute force " << blockStrandBruteForce[blockId] << std::endl;
      std::exit(1);
    } else {
      if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " strand state... passed: " << blockSequences.blockStrand[blockId] << " == " << blockStrandBruteForce[blockId] << std::endl;
    }

    if (blockStrandDynamic[blockId]) {
      for (int i = 0; i < sequenceDynamic[blockId].size(); i++) {      
        if (sequenceDynamic[blockId][i].second.size() != sequenceBruteForce[blockId][i].second.size()) {
          std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic[blockId][i].second.size() << " != brute force " << sequenceBruteForce[blockId][i].second.size() << std::endl;
          std::exit(1);
        } else {
          if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc size at (" << blockId << ", " << i << ")... passed: " << sequenceDynamic[blockId][i].second.size() << " == " << sequenceBruteForce[blockId][i].second.size() << std::endl;
        }
  
        for (int j = 0; j < sequenceDynamic[blockId][i].second.size(); j++) {
          if (sequenceDynamic[blockId][i].second[j] != sequenceBruteForce[blockId][i].second[j]) {
            std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", " << j << "): dynamic " << sequenceDynamic[blockId][i].second[j] << " != brute force " << sequenceBruteForce[blockId][i].second[j] << std::endl;
            std::exit(1);
          } else {
            if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc at (" << blockId << ", " << i << ", " << j << ")... passed: " << sequenceDynamic[blockId][i].second[j] << " == " << sequenceBruteForce[blockId][i].second[j] << std::endl;
          }
  
          if (sequenceDynamic[blockId][i].second[j] != '-') {
            if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") and global scalar coord  " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].second[j] << " ?= " << sequenceBruteForce[blockId][i].second[j] << std::endl;
              std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
            }
            ++localScalarCoordBruteForce;
          }
          ++globalScalarCoord;
        }
  
        if (sequenceDynamic[blockId][i].first != sequenceBruteForce[blockId][i].first) {
          std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", -1): dynamic " << sequenceDynamic[blockId][i].first << " != brute force " << sequenceBruteForce[blockId][i].first << std::endl;
          std::exit(1);
        } else {
          if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical main nuc at (" << blockId << ", " << i << ", -1)... passed: " << sequenceDynamic[blockId][i].first << " == " << sequenceBruteForce[blockId][i].first << std::endl;
        }
  
        if (sequenceDynamic[blockId][i].first != 'x') {
          if (sequenceDynamic[blockId][i].first != '-') {
            if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) and global scalar coord " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].first << " ?= " << sequenceBruteForce[blockId][i].first << std::endl;
              std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
            }
            ++localScalarCoordBruteForce;
          }
          ++globalScalarCoord;
        }
      }
    } else {
      for (int i = sequenceDynamic[blockId].size() - 1; i >= 0; i--) {
        if (sequenceDynamic[blockId][i].first != sequenceBruteForce[blockId][i].first) {
          std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", -1): dynamic " << sequenceDynamic[blockId][i].first << " != brute force " << sequenceBruteForce[blockId][i].first << std::endl;
          std::exit(1);
        } else {
          if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical main nuc at (" << blockId << ", " << i << ", -1)... passed: " << sequenceDynamic[blockId][i].first << " == " << sequenceBruteForce[blockId][i].first << std::endl;
        }
  
        if (sequenceDynamic[blockId][i].first != 'x') {
          if (sequenceDynamic[blockId][i].first != '-') {
            if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) and global scalar coord " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].first << " ?= " << sequenceBruteForce[blockId][i].first << std::endl;
              std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
            }
            ++localScalarCoordBruteForce;
          }
          ++globalScalarCoord;
        }

        if (sequenceDynamic[blockId][i].second.size() != sequenceBruteForce[blockId][i].second.size()) {
          std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic[blockId][i].second.size() << " != brute force " << sequenceBruteForce[blockId][i].second.size() << std::endl;
          std::exit(1);
        } else {
          if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc size at (" << blockId << ", " << i << ")... passed: " << sequenceDynamic[blockId][i].second.size() << " == " << sequenceBruteForce[blockId][i].second.size() << std::endl;
        }
  
        for (int j = sequenceDynamic[blockId][i].second.size() - 1; j >= 0; j--) {
          if (sequenceDynamic[blockId][i].second[j] != sequenceBruteForce[blockId][i].second[j]) {
            std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", " << j << "): dynamic " << sequenceDynamic[blockId][i].second[j] << " != brute force " << sequenceBruteForce[blockId][i].second[j] << std::endl;
            std::exit(1);
          } else {
            if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc at (" << blockId << ", " << i << ", " << j << ")... passed: " << sequenceDynamic[blockId][i].second[j] << " == " << sequenceBruteForce[blockId][i].second[j] << std::endl;
          }
  
          if (sequenceDynamic[blockId][i].second[j] != '-') {
            if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") and global scalar coord " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].second[j] << " ?= " << sequenceBruteForce[blockId][i].second[j] << std::endl;
              std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
            }
            ++localScalarCoordBruteForce;
          }
          ++globalScalarCoord;
        }

      }
    }

    if (globalScalarCoord - 1 != globalCoords.getBlockEndScalar(blockId)) {
      std::cerr << "Global scalar coord block " << blockId << " end mismatch: dynamic " << globalScalarCoord - 1 << " != brute force " << globalCoords.getBlockEndScalar(blockId) << std::endl;
      std::exit(1);
    } else {
      if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " end scalar coord... passed: " << globalScalarCoord - 1 << " == " << globalCoords.getBlockEndScalar(blockId) << std::endl;
    }
  }

  // check syncmers
  std::string ungappedSequence = panmapUtils::getStringFromSequence(sequenceBruteForce, blockLengthsBruteForce, blockExistsBruteForce, blockStrandBruteForce, false);
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersBruteForce = seeding::rollingSyncmers(ungappedSequence, k, s, open, t, false);
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersDynamic;
  for (size_t i = 0; i < refOnSyncmers.size(); i++) {
    if (refOnSyncmers[i].has_value()) {
      const auto& [hash, endPos, isReverse] = refOnSyncmers[i].value();
      syncmersDynamic.emplace_back(std::make_tuple(hash, isReverse, true, mgsr::degapGlobal(i, degapCoordIndex)));
    }
  }



  // check all syncmers
  if (syncmersDynamic.size() != syncmersBruteForce.size()) {
    std::cout << "Syncmer count mismatch: dynamic " << syncmersDynamic.size() << " != brute force " << syncmersBruteForce.size() << std::endl;
    std::cout << "Dynamic syncmers: ";
    for (const auto& syncmer : syncmersDynamic) {
      std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << mgsr::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
    }
    std::cout << std::endl;
    std::cout << "Brute force syncmers: ";
    for (const auto& syncmer : syncmersBruteForce) {
      std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << mgsr::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
    }
    std::cout << std::endl;
    std::exit(1);
  } else if (printCorrect && node->identifier == nodeToDebug) {
    std::cout << "Identical syncmer count... passed: " << syncmersDynamic.size() << " == " << syncmersBruteForce.size() << std::endl;
  }

  for (size_t i = 0; i < syncmersDynamic.size(); i++) {
    const auto& [hash, isReverse, isSeed, startPos] = syncmersDynamic[i];
    const auto& [hashBruteForce, isReverseBruteForce, isSeedBruteForce, startPosBruteForce] = syncmersBruteForce[i];
    if (hash != hashBruteForce || isReverse != isReverseBruteForce || startPos != startPosBruteForce) {
      std::cout << "Syncmer mismatch at " << i << "th syncmer: dynamic (" << hash << ", " << startPos << ", " << isReverse << ") != brute force (" << hashBruteForce << ", " << startPosBruteForce << ", " << isReverseBruteForce << ")" << std::endl;
      std::exit(1);
    }
  }


  std::cout << "         " << node->identifier << " states passed brute force check" << std::endl;
  
}

static void applyMutations (
  panmanUtils::Node *node,
  panmapUtils::BlockSequences &blockSequences,
  std::unordered_set<uint64_t>& invertedBlocks,
  panmapUtils::GlobalCoords& globalCoords,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
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
      if (!blockStrand[blockId]) {
        invertedBlocks.insert(blockId);
        invertedBlocksBacktracks.emplace_back(blockId, true);
      }
    } else if (isInversion) {
      blockStrand[blockId] = !blockStrand[blockId];
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


    if (blockStrand[blockId]) {
      // forward strand
      localMutationRanges.emplace_back(std::make_pair(globalCoords.getBlockStartCoord(blockId), globalCoords.getBlockEndCoord(blockId)));
    } else {
      // reversed strand
      localMutationRanges.emplace_back(std::make_pair(globalCoords.getBlockEndCoord(blockId), globalCoords.getBlockStartCoord(blockId)));
    }
  }

  // process nuc mutations
  for (const auto& nucMutation : node->nucMutation) {
    int length = nucMutation.mutInfo >> 4;
    int blockId;
    for (int i = 0; i < length; i++) {
      panmapUtils::Coordinate pos = panmapUtils::Coordinate(nucMutation, i);
      blockId = pos.primaryBlockId;
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
    if (blockExists[blockId] && oldBlockExists[blockId] && blockStrand[blockId] == oldBlockStrand[blockId]) {
      if (blockStrand[blockId]) {
        localMutationRanges.emplace_back(std::make_pair(panmapUtils::Coordinate(nucMutation, 0), panmapUtils::Coordinate(nucMutation, length - 1)));
      } else {
        localMutationRanges.emplace_back(std::make_pair(panmapUtils::Coordinate(nucMutation, length - 1), panmapUtils::Coordinate(nucMutation, 0)));
      }
    }
  }

  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      // on to off -> block range to all gaps
      uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
      uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
      gapRunUpdates.emplace_back(true, std::make_pair(beg, end));
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
              gapRunUpdates.emplace_back(false, std::make_pair((uint64_t)curNucRange.first, (uint64_t)curNucRange.second));
            }
            curNucRange = {scalar, scalar};
          }
        }

        if (coord == end) break;
        coord = globalCoords.stepRightCoordinate(coord);
      }
      if (curNucRange.first != -1) {
        gapRunUpdates.emplace_back(false, std::make_pair((uint64_t)curNucRange.first, (uint64_t)curNucRange.second));
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
  panmanUtils::Node *node,
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

    if (end < rightIt->first) {
      // completely completely between two ranges... all nucs are on.. do nothing
      return;
    }

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

std::vector<panmapUtils::NewSyncmerRange> mgsr::mgsrIndexBuilder::computeNewSyncmerRanges(
  panmanUtils::Node* node,
  const panmapUtils::BlockSequences& blockSequences,
  const panmapUtils::GlobalCoords& globalCoords,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
  std::vector<std::tuple<uint64_t, uint64_t, bool>>& blockOnSyncmersBacktracks
) {
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges;
  if (localMutationRanges.empty()) {
    return newSyncmerRanges;
  }

  const std::vector<bool>& blockExists = blockSequences.blockExists;
  const std::vector<bool>& blockStrand = blockSequences.blockStrand;
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;

  std::sort(localMutationRanges.begin(), localMutationRanges.end(), [&globalCoords, &blockStrand](const auto& a, const auto& b) {
    return globalCoords.getScalarFromCoord(a.first, blockStrand[a.first.primaryBlockId]) < globalCoords.getScalarFromCoord(b.first, blockStrand[b.first.primaryBlockId]);
  });


  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> mergedLocalMutationRanges{localMutationRanges.front()};
  for (size_t i = 1; i < localMutationRanges.size(); ++i) {
    const auto& [curBeg, curEnd] = mergedLocalMutationRanges.back();
    const auto& [nextBeg, nextEnd] = localMutationRanges[i];
    
    // check if the current range and the next range are adjacent on their global scalar coordinates
    if (globalCoords.getScalarFromCoord(curEnd, blockStrand[curBeg.primaryBlockId]) + 1 >= globalCoords.getScalarFromCoord(nextBeg, blockStrand[nextBeg.primaryBlockId])) {
      if (globalCoords.getScalarFromCoord(nextEnd, blockStrand[nextBeg.primaryBlockId]) > globalCoords.getScalarFromCoord(curEnd, blockStrand[curBeg.primaryBlockId])) {
        mergedLocalMutationRanges.back().second = nextEnd;
      }
    } else {
      mergedLocalMutationRanges.emplace_back(nextBeg, nextEnd);
    }
  }

  int k = indexBuilder.getK();
  size_t localMutationRangeIndex = 0;
  int offsetsToDelete = -1;
  int endOffset = -1;
  panmapUtils::NewSyncmerRange curSyncmerRange;
  while (localMutationRangeIndex < mergedLocalMutationRanges.size()) {
    auto [curBegCoord, curEndCoord] = mergedLocalMutationRanges[localMutationRangeIndex];
    auto syncmerRangeBegCoord = curBegCoord;
    auto syncmerRangeEndCoord = curEndCoord;
    // expand to the left... if reach newSyncmerRanges.back(), merge
    bool reachedEnd = false;
    uint32_t offset = 0;
    while (offset < k - 1) {
      if (globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]) == 0) {
        break;
      }
      curBegCoord = globalCoords.stepBackwardScalar(curBegCoord, blockStrand);
      if (!blockExists[curBegCoord.primaryBlockId]) {
        curBegCoord = globalCoords.getBlockStartCoord(curBegCoord.primaryBlockId);
        continue;
      }
      if (!newSyncmerRanges.empty() 
          && globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]) <= globalCoords.getScalarFromCoord(curSyncmerRange.endCoord, blockStrand[curSyncmerRange.endCoord.primaryBlockId])
      ) {
        // reached current newSyncmerRange... merge
        curBegCoord = curSyncmerRange.begCoord;
        syncmerRangeBegCoord = curBegCoord;
        newSyncmerRanges.pop_back();
        break;
      }

      if (blockExists[curBegCoord.primaryBlockId] && blockSequences.getSequenceBase(curBegCoord) != '-') {
        offset++;
        syncmerRangeBegCoord = curBegCoord;
      }
    }

    // expand to the right... if reach mergedLocalMutationRanges[localMutationRangeIndex + 1], merge
    offset = 0;
    while (offset < k - 1) {
      if (globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]) == globalCoords.lastScalarCoord) {
        reachedEnd = true;
        break;
      }
      curEndCoord = globalCoords.stepForwardScalar(curEndCoord, blockStrand);
      if (!blockExists[curEndCoord.primaryBlockId]) {
        curEndCoord = globalCoords.getBlockEndCoord(curEndCoord.primaryBlockId);
        continue;
      }
      if (localMutationRangeIndex != mergedLocalMutationRanges.size() - 1
          && globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]) >= globalCoords.getScalarFromCoord(mergedLocalMutationRanges[localMutationRangeIndex + 1].first, blockStrand[mergedLocalMutationRanges[localMutationRangeIndex + 1].first.primaryBlockId])
      ) {
        // reached next mutation range... merge
        curEndCoord = mergedLocalMutationRanges[localMutationRangeIndex + 1].second;
        syncmerRangeEndCoord = curEndCoord;
        localMutationRangeIndex++;
        offset = 0;
        continue;
      }

      if (blockExists[curEndCoord.primaryBlockId] && blockSequences.getSequenceBase(curEndCoord) != '-') {
        offset++;
        syncmerRangeEndCoord = curEndCoord;
      }
    }

    newSyncmerRanges.emplace_back(syncmerRangeBegCoord, syncmerRangeEndCoord, "", std::vector<uint64_t>(), std::vector<uint64_t>());
    if (reachedEnd) {
      offsetsToDelete = k - offset - 1;
      endOffset = offset;
      break;
    }
    localMutationRangeIndex++;
  }

  for (size_t i = 0; i < newSyncmerRanges.size(); i++) {
    panmapUtils::NewSyncmerRange& syncmerRange = newSyncmerRanges[i];
    panmapUtils::Coordinate curCoord = syncmerRange.begCoord;
    panmapUtils::Coordinate curEndCoord = syncmerRange.endCoord;
    
    std::string& localRangeSeq = syncmerRange.localRangeSeq;
    std::vector<uint64_t>& localRangeCoordToGlobalScalarCoords = syncmerRange.localRangeCoordToGlobalScalarCoords;
    std::vector<uint64_t>& seedsToDelete = syncmerRange.seedsToDelete;
    std::vector<uint64_t>& localRangeCoordToBlockId = syncmerRange.localRangeCoordToBlockId;
    localRangeSeq = "";
    std::vector<uint64_t>().swap(localRangeCoordToGlobalScalarCoords);
    std::vector<uint64_t>().swap(seedsToDelete);
    std::vector<uint64_t>().swap(localRangeCoordToBlockId);
    while (true) {
      if (!blockExists[curCoord.primaryBlockId]) {
        if (curCoord.primaryBlockId == curEndCoord.primaryBlockId) {
          break;
        }
        curCoord = globalCoords.stepForwardScalar(globalCoords.getBlockEndCoord(curCoord.primaryBlockId), blockStrand);
        continue;
      }
      
      auto curScalarCoord = globalCoords.getScalarFromCoord(curCoord, blockStrand[curCoord.primaryBlockId]);
      char curNuc = blockSequences.getSequenceBase(curCoord);
      if (curNuc != '-') {
        if (!blockStrand[curCoord.primaryBlockId]) curNuc = panmanUtils::getComplementCharacter(curNuc);
        localRangeSeq += curNuc;
        localRangeCoordToGlobalScalarCoords.push_back(curScalarCoord);
        localRangeCoordToBlockId.push_back(curCoord.primaryBlockId);
      } else if (refOnSyncmers[curScalarCoord].has_value()) {
        blockOnSyncmers[curCoord.primaryBlockId].erase(curScalarCoord);
        if (blockOnSyncmers[curCoord.primaryBlockId].empty()) blockOnSyncmers.erase(curCoord.primaryBlockId);
        blockOnSyncmersBacktracks.emplace_back(curCoord.primaryBlockId, curScalarCoord, false);
        seedsToDelete.push_back(curScalarCoord);
      }
      if (curCoord == curEndCoord) break;
      curCoord = globalCoords.stepForwardScalar(curCoord, blockStrand);
    }

    if (i == newSyncmerRanges.size() - 1 && offsetsToDelete != -1) {
      for (size_t j = localRangeSeq.size() - endOffset - offsetsToDelete; j < localRangeSeq.size(); j++) {
        if (refOnSyncmers[localRangeCoordToGlobalScalarCoords[j]].has_value()) {
          auto curScalarCoord = localRangeCoordToGlobalScalarCoords[j];
          auto curBlockId = localRangeCoordToBlockId[j];
          blockOnSyncmers[curBlockId].erase(curScalarCoord);
          if (blockOnSyncmers[curBlockId].empty()) blockOnSyncmers.erase(curBlockId);
          blockOnSyncmersBacktracks.emplace_back(curBlockId, curScalarCoord, false);
          seedsToDelete.push_back(localRangeCoordToGlobalScalarCoords[j]);
        }
      }
    }
  }
  return newSyncmerRanges;
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
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBlockInversionBacktracks;
  std::vector<std::pair<uint64_t, bool>> invertedBlocksBacktracks;

  // for computing new syncmers
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> localMutationRanges;
  std::vector<std::tuple<uint64_t, bool, seeding::rsyncmer_t>> refOnSyncmersBacktracks;
  std::vector<std::tuple<uint64_t, uint64_t, bool>> blockOnSyncmersBacktracks;


  applyMutations(node, blockSequences, invertedBlocks, globalCoords, localMutationRanges, blockMutationRecord, nucMutationRecord, gapRunUpdates, invertedBlocksBacktracks, blockExistsDelayed, blockStrandDelayed);
  
  updateGapMap(node, gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapRunUpdates); // gapRunUpdates is no longer needed... clear memory

  std::vector<uint64_t> invertedBlocksVec(invertedBlocks.begin(), invertedBlocks.end());
  std::sort(invertedBlocksVec.begin(), invertedBlocksVec.end());
  for (const auto& blockId : invertedBlocksVec) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);

    invertGapMap(gapMap, {beg, end}, gapRunBlockInversionBacktracks, gapMapUpdates);
  }

  std::map<uint64_t, uint64_t> degapCoordIndex;
  std::map<uint64_t, uint64_t> regapCoordIndex;
  makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, globalCoords);


  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = computeNewSyncmerRanges(node, blockSequences, globalCoords, localMutationRanges, blockOnSyncmersBacktracks);

  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    if (localRangeSeq.size() >= indexBuilder.getK()) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(localRangeSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = localRangeCoordToGlobalScalarCoords[startPos];
        auto endPosGlobal = localRangeCoordToGlobalScalarCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = localRangeCoordToBlockId[startPos];
        bool wasSeed = refOnSyncmers[startPosGlobal].has_value();
        if (!wasSeed && isSeed) {
          refOnSyncmersBacktracks.emplace_back(startPosGlobal, true, seeding::rsyncmer_t());
          blockOnSyncmersBacktracks.emplace_back(curBlockId, startPosGlobal, true);
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
          blockOnSyncmers[curBlockId].insert(startPosGlobal);
        } else if (wasSeed && !isSeed) {
          refOnSyncmersBacktracks.emplace_back(startPosGlobal, false, refOnSyncmers[startPosGlobal].value());
          blockOnSyncmersBacktracks.emplace_back(curBlockId, startPosGlobal, false);
          refOnSyncmers[startPosGlobal] = std::nullopt;
          blockOnSyncmers[curBlockId].erase(startPosGlobal);
          if (blockOnSyncmers[localRangeCoordToBlockId[startPos]].empty()) {
            blockOnSyncmers.erase(localRangeCoordToBlockId[startPos]);
          }
        } else if (wasSeed && isSeed) {
          refOnSyncmersBacktracks.emplace_back(startPosGlobal, false, refOnSyncmers[startPosGlobal].value());
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }

    for (uint64_t pos : seedsToDelete) {
      if (!refOnSyncmers[pos].has_value()) {
        std::cerr << "Error: refOnSyncmers[" << pos << "] is null" << std::endl;
        std::exit(1);
      }
      refOnSyncmersBacktracks.emplace_back(pos, false, refOnSyncmers[pos].value());
      refOnSyncmers[pos] = std::nullopt;
      // blockOnSyncmers and blockOnSyncmersBacktracks are updated in computeNewSyncmerRanges()
    }
  }
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      if (blockOnSyncmers.find(blockId) != blockOnSyncmers.end()) {
        for (uint64_t pos : blockOnSyncmers[blockId]) {
          refOnSyncmersBacktracks.emplace_back(pos, false, refOnSyncmers[pos].value());
          blockOnSyncmersBacktracks.emplace_back(blockId, pos, false);
          refOnSyncmers[pos] = std::nullopt;
        }
        blockOnSyncmers.erase(blockId);
      }
    }
  }


  // compare with brute force for debugging
  if (node->children.size() == 0) {
    compareBruteForce(T, node, blockSequences, globalCoords, gapMap, degapCoordIndex, regapCoordIndex, refOnSyncmers, blockOnSyncmers, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getT(), indexBuilder.getL(), indexBuilder.getOpen());
  }


  for (auto it = gapRunBlockInversionBacktracks.rbegin(); it != gapRunBlockInversionBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapRunBlockInversionBacktracks); // gapRunBlockInversionBacktracks is no longer needed... clear memory

  // update delayed block states
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExistsDelayed[blockId] = newExists;
    blockStrandDelayed[blockId] = newStrand;
  }

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

  for (const auto& [blockId, del] : invertedBlocksBacktracks) {
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

  for (const auto& [pos, del, rsyncmer] : refOnSyncmersBacktracks) {
    if (del) {
      refOnSyncmers[pos] = std::nullopt;
    } else {
      refOnSyncmers[pos] = rsyncmer;
    }
  }
  
  for (const auto& [blockId, pos, del] : blockOnSyncmersBacktracks) {
    if (del) {
      blockOnSyncmers[blockId].erase(pos);
      if (blockOnSyncmers[blockId].empty()) {
        blockOnSyncmers.erase(blockId);
      }
    } else {
      blockOnSyncmers[blockId].insert(pos);
    }
  }
}

void mgsr::mgsrIndexBuilder::buildIndex() {
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences.sequence);
  refOnSyncmers.resize(globalCoords.lastScalarCoord + 1);

  std::vector<bool> blockExistsDelayed = blockSequences.blockExists;
  std::vector<bool> blockStrandDelayed = blockSequences.blockStrand;

  std::map<uint64_t, uint64_t> gapMap{{0, globalCoords.lastScalarCoord}};
  std::unordered_set<uint64_t> invertedBlocks;

  buildIndexHelper(T->root, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks);
  std::cout << "Finished building index!" << std::endl;
}


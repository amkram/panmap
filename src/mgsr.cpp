#pragma once

#include "mgsr.hpp"
#include "panmap_utils.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include <fcntl.h>
#include <unistd.h>

static void compareBruteForceBuild(
  panmanUtils::Tree *T,
  panmanUtils::Node *node,
  const panmapUtils::BlockSequences& blockSequences,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  const std::map<uint64_t, uint64_t>& degapCoordIndex,
  const std::map<uint64_t, uint64_t>& regapCoordIndex,
  const std::vector<std::optional<seeding::rsyncmer_t>>& refOnSyncmers,
  const std::set<uint64_t>& refOnSyncmersMap,
  const std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers,
  const std::vector<std::optional<uint64_t>>& refOnKminmers,
  const std::vector<seeding::uniqueKminmer_t>& uniqueKminmers,
  const std::unordered_map<seeding::uniqueKminmer_t, uint64_t>& kminmerToUniqueIndex,
  int k, int s, int t, int l, bool open
) {
  bool printCorrect = false;
  std::string nodeToDebug = "KJ627733.1";
  std::cout << "Checking " << node->identifier << " states with brute force... " << std::flush;

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

  std::cout << "sequence and coordinate objects passed... " << std::flush;

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
  
  auto curSyncmerOnMapIt = refOnSyncmersMap.begin();
  for (size_t i = 0; i < syncmersDynamic.size(); i++) {
    const auto& [hash, isReverse, isSeed, startPos] = syncmersDynamic[i];
    const auto& [hashBruteForce, isReverseBruteForce, isSeedBruteForce, startPosBruteForce] = syncmersBruteForce[i];
    if (hash != hashBruteForce || isReverse != isReverseBruteForce || startPos != startPosBruteForce) {
      std::cout << "Syncmer mismatch at " << i << "th syncmer: dynamic (" << hash << ", " << startPos << ", " << isReverse << ") != brute force (" << hashBruteForce << ", " << startPosBruteForce << ", " << isReverseBruteForce << ")" << std::endl;
      std::exit(1);
    }
    if (mgsr::regapGlobal(startPos, regapCoordIndex) != *curSyncmerOnMapIt) {
      std::cout << "Syncmer on map mismatch at " << i << "th syncmer: dynamic/bruteforce " << startPos << " != map " << *curSyncmerOnMapIt << std::endl;
      std::exit(1);
    }
    ++curSyncmerOnMapIt;
  }
  if (curSyncmerOnMapIt != refOnSyncmersMap.end()) {
    std::cout << "SyncmerOnMap has more elements than syncmers: " << refOnSyncmersMap.size() << " != " << syncmersDynamic.size() << std::endl;
    std::exit(1);
  }
  std::cout << "syncmers passed... " << std::flush;

  // check k-min-mers
  std::vector<std::tuple<size_t, size_t, size_t, bool>> kminmersBruteForce;
  if (syncmersBruteForce.size() >= l) {
    size_t forwardRolledHash = 0;
    size_t reverseRolledHash = 0;
  
    // first kminmer
    for (size_t i = 0; i < l; ++i) {
      forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ std::get<0>(syncmersBruteForce[i]);
      reverseRolledHash = seeding::rol(reverseRolledHash, k) ^ std::get<0>(syncmersBruteForce[l-i-1]);
    }
  
    if (forwardRolledHash != reverseRolledHash) {
      size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
      kminmersBruteForce.emplace_back(minHash, std::get<3>(syncmersBruteForce[0]), std::get<3>(syncmersBruteForce[l-1])+k-1, reverseRolledHash < forwardRolledHash);
    }
  
    // rest of kminmers
    for (uint64_t i = 1; i < syncmersBruteForce.size()-l+1; ++i) {
      if (!std::get<2>(syncmersBruteForce[i-1]) || !std::get<2>(syncmersBruteForce[i+l-1])) {
        std::cout << "invalid syncmer" << std::endl;
        exit(0);
      }
      const size_t& prevSyncmerHash = std::get<0>(syncmersBruteForce[i-1]);
      const size_t& nextSyncmerHash = std::get<0>(syncmersBruteForce[i+l-1]);
      forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ seeding::rol(prevSyncmerHash, k * l) ^ nextSyncmerHash;
      reverseRolledHash = seeding::ror(reverseRolledHash, k) ^ seeding::ror(prevSyncmerHash, k)     ^ seeding::rol(nextSyncmerHash, k * (l-1));
  
      if (forwardRolledHash != reverseRolledHash) {
        size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
        kminmersBruteForce.emplace_back(minHash, std::get<3>(syncmersBruteForce[i]), std::get<3>(syncmersBruteForce[i+l-1])+k-1, reverseRolledHash < forwardRolledHash);
      }
    }
  }


  std::vector<std::tuple<size_t, size_t, size_t, bool>> kminmersDynamic;
  for (size_t i = 0; i < refOnKminmers.size(); i++) {
    if (refOnKminmers[i].has_value()) {
      const auto& [startPos, endPos, hash, isReverse] = uniqueKminmers[refOnKminmers[i].value()];
      kminmersDynamic.emplace_back(std::make_tuple(hash, mgsr::degapGlobal(startPos, degapCoordIndex), mgsr::degapGlobal(endPos, degapCoordIndex), isReverse));
    }
  }

  if (kminmersDynamic.size() != kminmersBruteForce.size()) {
    std::cout << "K-min-mer count mismatch: dynamic " << kminmersDynamic.size() << " != brute force " << kminmersBruteForce.size() << std::endl;
    std::cout << "Dynamic k-min-mers: ";
    for (const auto& kminmer : kminmersDynamic) {
      std::cout << "(" << std::get<0>(kminmer) << ", " << std::get<1>(kminmer) << ", " << mgsr::regapGlobal(std::get<1>(kminmer), regapCoordIndex) << ", " << std::get<2>(kminmer) << ", " << mgsr::regapGlobal(std::get<2>(kminmer), regapCoordIndex) << ", " << std::get<3>(kminmer) << ") ";
    }
    std::cout << std::endl;
    std::cout << "Brute force k-min-mers: ";
    for (const auto& kminmer : kminmersBruteForce) {
      std::cout << "(" << std::get<0>(kminmer) << ", " << std::get<1>(kminmer) << ", " << mgsr::regapGlobal(std::get<1>(kminmer), regapCoordIndex) << ", " << std::get<2>(kminmer) << ", " << mgsr::regapGlobal(std::get<2>(kminmer), regapCoordIndex) << ", " << std::get<3>(kminmer) << ") ";
    }
    std::cout << std::endl;
    std::exit(1);
  } else if (printCorrect && node->identifier == nodeToDebug) {
    std::cout << "Identical k-min-mer count... passed: " << kminmersDynamic.size() << " == " << kminmersBruteForce.size() << std::endl;
  }

  for (size_t i = 0; i < kminmersDynamic.size(); i++) {
    const auto& [hash, startPos, endPos, isReverse] = kminmersDynamic[i];
    const auto& [hashBruteForce, startPosBruteForce, endPosBruteForce, isReverseBruteForce] = kminmersBruteForce[i];
    if (hash != hashBruteForce || startPos != startPosBruteForce || endPos != endPosBruteForce || isReverse != isReverseBruteForce) {
      std::cout << "K-min-mer mismatch at " << i << "th k-min-mer: dynamic (" << hash << ", " << startPos << ", " << endPos << ", " << isReverse << ") != brute force (" << hashBruteForce << ", " << startPosBruteForce << ", " << endPosBruteForce << ", " << isReverseBruteForce << ")" << std::endl;
      std::exit(1);
    }
  }
  std::cout << "k-min-mers passed... " << std::flush;

  std::cout << "         " << node->identifier << " states passed brute force check" << std::endl;
  
}

static void compareBruteForcePlace(
  panmanUtils::Tree *T,
  panmanUtils::Node *node,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  const std::map<uint64_t, uint64_t>& degapCoordIndex,
  const std::map<uint64_t, uint64_t>& regapCoordIndex,
  const std::map<uint64_t, uint64_t>& positionMap,
  const std::unordered_map<size_t, std::set<std::map<uint64_t, uint64_t>::iterator, mgsr::IteratorComparator>>& hashToPositionMap,
  const ::capnp::List<SeedInfo>::Reader& seedInfos,
  int k, int s, int t, int l, bool open
) {
  std::cout << "Checking " << node->identifier << " states with brute force at placement... " << std::flush;
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
  
  std::string ungappedSequence = panmapUtils::getStringFromSequence(sequenceBruteForce, blockLengthsBruteForce, blockExistsBruteForce, blockStrandBruteForce, false);
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersBruteForce = seeding::rollingSyncmers(ungappedSequence, k, s, open, t, false);
  // check k-min-mers
  std::vector<std::tuple<size_t, size_t, size_t, bool>> kminmersBruteForce;
  if (syncmersBruteForce.size() >= l) {
    size_t forwardRolledHash = 0;
    size_t reverseRolledHash = 0;
  
    // first kminmer
    for (size_t i = 0; i < l; ++i) {
      forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ std::get<0>(syncmersBruteForce[i]);
      reverseRolledHash = seeding::rol(reverseRolledHash, k) ^ std::get<0>(syncmersBruteForce[l-i-1]);
    }
  
    if (forwardRolledHash != reverseRolledHash) {
      size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
      kminmersBruteForce.emplace_back(minHash, std::get<3>(syncmersBruteForce[0]), std::get<3>(syncmersBruteForce[l-1])+k-1, reverseRolledHash < forwardRolledHash);
    }
  
    // rest of kminmers
    for (uint64_t i = 1; i < syncmersBruteForce.size()-l+1; ++i) {
      if (!std::get<2>(syncmersBruteForce[i-1]) || !std::get<2>(syncmersBruteForce[i+l-1])) {
        std::cout << "invalid syncmer" << std::endl;
        exit(0);
      }
      const size_t& prevSyncmerHash = std::get<0>(syncmersBruteForce[i-1]);
      const size_t& nextSyncmerHash = std::get<0>(syncmersBruteForce[i+l-1]);
      forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ seeding::rol(prevSyncmerHash, k * l) ^ nextSyncmerHash;
      reverseRolledHash = seeding::ror(reverseRolledHash, k) ^ seeding::ror(prevSyncmerHash, k)     ^ seeding::rol(nextSyncmerHash, k * (l-1));
  
      if (forwardRolledHash != reverseRolledHash) {
        size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
        kminmersBruteForce.emplace_back(minHash, std::get<3>(syncmersBruteForce[i]), std::get<3>(syncmersBruteForce[i+l-1])+k-1, reverseRolledHash < forwardRolledHash);
      }
    }
  }

  std::vector<std::tuple<size_t, size_t, size_t, bool>> kminmersBruteDynamic;
  for (const auto& [pos, infoIndex] : positionMap) {
    kminmersBruteDynamic.emplace_back(
      seedInfos[infoIndex].getHash(),
      mgsr::degapGlobal(seedInfos[infoIndex].getStartPos(), degapCoordIndex),
      mgsr::degapGlobal(seedInfos[infoIndex].getEndPos(), degapCoordIndex),
      seedInfos[infoIndex].getIsReverse()
    );
  }

  if (kminmersBruteDynamic.size() != kminmersBruteForce.size()) {
    std::cout << "K-min-mer count mismatch: dynamic " << kminmersBruteDynamic.size() << " != brute force " << kminmersBruteForce.size() << std::endl;
    std::exit(1);
  }

  for (size_t i = 0; i < kminmersBruteDynamic.size(); i++) {
    const auto& [hash, startPos, endPos, isReverse] = kminmersBruteDynamic[i];
    const auto& [hashBruteForce, startPosBruteForce, endPosBruteForce, isReverseBruteForce] = kminmersBruteForce[i];
    if (hash != hashBruteForce || startPos != startPosBruteForce || endPos != endPosBruteForce || isReverse != isReverseBruteForce) {
      std::cout << "K-min-mer mismatch at " << i << "th k-min-mer: dynamic (" << hash << ", " << startPos << ", " << endPos << ", " << isReverse << ") != brute force (" << hashBruteForce << ", " << startPosBruteForce << ", " << endPosBruteForce << ", " << isReverseBruteForce << ")" << std::endl;
      std::exit(1);
    }
  }

  std::cout << "passed" << std::endl;
}

static void applyMutations (
  panmanUtils::Node *node,
  size_t dfsIndex,
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
    int lastOffset = -1;
    for (int i = 0; i < length; i++) {
      panmapUtils::Coordinate pos = panmapUtils::Coordinate(nucMutation, i);
      if ((pos.nucPosition == sequence[pos.primaryBlockId].size() - 1 && pos.nucGapPosition == -1) ||
          (pos.nucPosition >= sequence[pos.primaryBlockId].size())) {
        continue;
      }
      lastOffset = i;
      blockId = pos.primaryBlockId;
      char oldNuc = blockSequences.getSequenceBase(pos);
      int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
      char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);

      if (oldNuc == newNuc) continue;
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
    if (lastOffset != -1 &&blockExists[blockId] && oldBlockExists[blockId] && blockStrand[blockId] == oldBlockStrand[blockId]) {
      if (blockStrand[blockId]) {
        localMutationRanges.emplace_back(std::make_pair(panmapUtils::Coordinate(nucMutation, 0), panmapUtils::Coordinate(nucMutation, lastOffset)));
      } else {
        localMutationRanges.emplace_back(std::make_pair(panmapUtils::Coordinate(nucMutation, lastOffset), panmapUtils::Coordinate(nucMutation, 0)));
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
  size_t dfsIndex,
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

void mgsr::revertGapMapInversions(
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBlocksBacktracks,
  std::map<uint64_t, uint64_t>& gapMap
) {
  for (auto it = gapMapBlocksBacktracks.rbegin(); it != gapMapBlocksBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }
}

void mgsr::makeCoordIndex(
  std::map<uint64_t, uint64_t>& degapCoordIndex,
  std::map<uint64_t, uint64_t>& regapCoordIndex,
  const std::map<uint64_t, uint64_t>& gapMap,
  uint64_t lastScalarCoord
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
    if (gapEnd == lastScalarCoord) break;
    totalGapSize += gapSize;
    degapCoordIndex[gapEnd+1] = totalGapSize;
    regapCoordIndex[gapEnd+1-totalGapSize] = totalGapSize;
  }
}

std::vector<panmapUtils::NewSyncmerRange> mgsr::mgsrIndexBuilder::computeNewSyncmerRanges(
  panmanUtils::Node* node,
  size_t dfsIndex,
  const panmapUtils::BlockSequences& blockSequences,
  const panmapUtils::GlobalCoords& globalCoords,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersChangeRecord
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
        blockOnSyncmersChangeRecord.emplace_back(curCoord.primaryBlockId, curScalarCoord, panmapUtils::seedChangeType::DEL);
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
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, curScalarCoord, panmapUtils::seedChangeType::DEL);
          seedsToDelete.push_back(localRangeCoordToGlobalScalarCoords[j]);
        }
      }
    }
  }
  return newSyncmerRanges;
}

std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> mgsr::mgsrIndexBuilder::computeNewKminmerRanges(
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord
) {
  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges;

  if (refOnSyncmersMap.size() < indexBuilder.getL()) {
    // no k-min-mers to compute, erase all k-min-mers
    return newKminmerRanges;
  }

  
  std::sort(refOnSyncmersChangeRecord.begin(), refOnSyncmersChangeRecord.end(), [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });


  int64_t syncmerChangeIndex = 0;
  while (syncmerChangeIndex < refOnSyncmersChangeRecord.size()) {
    const auto& [syncmerPos, changeType, rsyncmer] = refOnSyncmersChangeRecord[syncmerChangeIndex];
    std::set<uint64_t>::iterator curBegIt, curEndIt;
    if (changeType == panmapUtils::seedChangeType::DEL) {
      if (syncmerPos < *refOnSyncmersMap.begin()) {
        syncmerChangeIndex++;
        continue;
      } else {
        curEndIt = refOnSyncmersMap.upper_bound(syncmerPos);
        curBegIt = std::prev(curEndIt);
      }
    } else {
      auto it = refOnSyncmersMap.lower_bound(syncmerPos);
      curBegIt = syncmerPos == *refOnSyncmersMap.begin() ? refOnSyncmersMap.begin() : std::prev(it);
      curEndIt = syncmerPos == *refOnSyncmersMap.rbegin() ? refOnSyncmersMap.end() : std::next(it);
    }
    
    
    int offset = 1;
    // expand to the left
    while (offset < indexBuilder.getL() - 1 && curBegIt != refOnSyncmersMap.begin()) {
      if (!newKminmerRanges.empty() && *curBegIt <= *(newKminmerRanges.back().second)) {
        curBegIt = newKminmerRanges.back().first;
        newKminmerRanges.pop_back();
        break;
      }
      --curBegIt;
      offset++;
    }

    offset = 1;
    // expand to the right
    while (offset < indexBuilder.getL() - 1 && curEndIt != refOnSyncmersMap.end()) {
      if (syncmerChangeIndex != refOnSyncmersChangeRecord.size() - 1) {
        const auto& [nextSyncmerPos, nextChangeType, nextRsyncmer] = refOnSyncmersChangeRecord[syncmerChangeIndex + 1];
        if (nextChangeType != panmapUtils::seedChangeType::DEL) {
          if (*curEndIt >= nextSyncmerPos) {
            auto nextIt = refOnSyncmersMap.lower_bound(nextSyncmerPos);
            curEndIt = nextSyncmerPos == *refOnSyncmersMap.rbegin() ? refOnSyncmersMap.end() : std::next(nextIt);
            syncmerChangeIndex++;
            offset = 1;
            continue;
          }
        } else {
          if (nextSyncmerPos < *refOnSyncmersMap.rbegin()) {
            if (*curEndIt >= nextSyncmerPos) {
              curEndIt = std::next(curEndIt);
              syncmerChangeIndex++;
              offset = 1;
              continue;
            }
          } else {
            curEndIt = refOnSyncmersMap.end();
            break;
          }
        }
      }
      ++curEndIt;
      offset++;
    }

    newKminmerRanges.emplace_back(curBegIt, curEndIt);
    if (curEndIt == refOnSyncmersMap.end()) {
      break;
    }
    syncmerChangeIndex++;
  }
  return newKminmerRanges;
}

void mgsr::mgsrIndexBuilder::buildIndexHelper(
  panmanUtils::Node *node,
  panmapUtils::BlockSequences &blockSequences,
  std::vector<bool> &blockExistsDelayed,
  std::vector<bool> &blockStrandDelayed,
  panmapUtils::GlobalCoords &globalCoords,
  std::map<uint64_t, uint64_t> &gapMap,
  std::unordered_set<uint64_t> &invertedBlocks,
  uint64_t &dfsIndex
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
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>> refOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>> blockOnSyncmersChangeRecord;

  // for computing new k-min-mers
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, uint64_t>> refOnKminmersChangeRecord;

  applyMutations(node, dfsIndex, blockSequences, invertedBlocks, globalCoords, localMutationRanges, blockMutationRecord, nucMutationRecord, gapRunUpdates, invertedBlocksBacktracks, blockExistsDelayed, blockStrandDelayed);

  updateGapMap(node, dfsIndex, gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapRunUpdates); // gapRunUpdates is no longer needed... clear memory

  std::vector<uint64_t> invertedBlocksVec(invertedBlocks.begin(), invertedBlocks.end());
  std::sort(invertedBlocksVec.begin(), invertedBlocksVec.end());
  for (const auto& blockId : invertedBlocksVec) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);

    invertGapMap(gapMap, {beg, end}, gapRunBlockInversionBacktracks, gapMapUpdates);
  }

  // // not really needed for building the index... But do need to keep it for debugging and comparing with brute force
  // std::map<uint64_t, uint64_t> degapCoordIndex;
  // std::map<uint64_t, uint64_t> regapCoordIndex;
  // makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, (uint64_t)globalCoords.lastScalarCoord);



  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = computeNewSyncmerRanges(node, dfsIndex, blockSequences, globalCoords, localMutationRanges, blockOnSyncmersChangeRecord);

  // processing syncmers
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    if (localRangeSeq.size() >= indexBuilder.getK()) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(localRangeSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = localRangeCoordToGlobalScalarCoords[startPos];
        auto endPosGlobal = localRangeCoordToGlobalScalarCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = localRangeCoordToBlockId[startPos];
        bool wasSeed = refOnSyncmers[startPosGlobal].has_value();
        if (!wasSeed && isSeed) {
          auto it = refOnSyncmersMap.insert(startPosGlobal).first;
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::ADD, seeding::rsyncmer_t());
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::ADD);
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
          blockOnSyncmers[curBlockId].insert(startPosGlobal);
        } else if (wasSeed && !isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::DEL, refOnSyncmers[startPosGlobal].value());
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::DEL);
          refOnSyncmers[startPosGlobal] = std::nullopt;
          refOnSyncmersMap.erase(startPosGlobal);
          blockOnSyncmers[curBlockId].erase(startPosGlobal);
          if (blockOnSyncmers[localRangeCoordToBlockId[startPos]].empty()) {
            blockOnSyncmers.erase(localRangeCoordToBlockId[startPos]);
          }
        } else if (wasSeed && isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::SUB, refOnSyncmers[startPosGlobal].value());
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }

    for (uint64_t pos : seedsToDelete) {
      if (!refOnSyncmers[pos].has_value()) {
        std::cerr << "Error: refOnSyncmers[" << pos << "] is null" << std::endl;
        std::exit(1);
      }
      refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers[pos].value());
      refOnSyncmers[pos] = std::nullopt;
      refOnSyncmersMap.erase(pos);
      // blockOnSyncmers and blockOnSyncmersChangeRecord are updated in computeNewSyncmerRanges()
    }
  }

  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      if (blockOnSyncmers.find(blockId) != blockOnSyncmers.end()) {
        for (uint64_t pos : blockOnSyncmers[blockId]) {
          refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers[pos].value());
          blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);

          refOnSyncmers[pos] = std::nullopt;
          refOnSyncmersMap.erase(pos);
        }
        blockOnSyncmers.erase(blockId);
      }
    }
  }

  // processing k-min-mers
  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges = computeNewKminmerRanges(refOnSyncmersChangeRecord);

  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto beg = *newKminmerRanges[i].first;
    auto end = *newKminmerRanges[i].second;
    if (newKminmerRanges[i].second != refOnSyncmersMap.end() && beg > end) {
      std::cerr << "Error: beg (" << beg << ") > end (" << end << ") in node " << node->identifier << std::endl;
      std::exit(1);
    }
    if (i != 0 && *newKminmerRanges[i-1].second > beg) {
      std::cerr << "Error: newKminmerRanges[" << i-1 << "].second (" << *newKminmerRanges[i-1].second << ") > newKminmerRanges[" << i << "].first (" << *newKminmerRanges[i].first << ") in node " << node->identifier << std::endl;
      std::exit(1);
    }
  }

  auto k = indexBuilder.getK();
  auto l = indexBuilder.getL();
  std::vector<uint64_t> deletedSeedIndices;
  std::vector<uint64_t> addedSeedIndices;
  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto [curIt, endIt] = newKminmerRanges[i];
    auto indexingIt = curIt;
    size_t forwardHash = 0;
    size_t reverseHash = 0;
    bool shortRange = false;

    std::vector<size_t> startingSyncmerHashes;
    for (size_t j = 0; j < l; j++) {
      if (curIt == refOnSyncmersMap.end()) {
        break;
      } else if (endIt != refOnSyncmersMap.end() && *curIt == *endIt && j != l - 1) {
        break;
      }
      startingSyncmerHashes.push_back(refOnSyncmers[*curIt].value().hash);
      if (j != l - 1) ++curIt;
    }
    if (startingSyncmerHashes.size() < l) {
      continue;
    }
    for (size_t j = 0; j < l; j++) {
      forwardHash = seeding::rol(forwardHash, k) ^ startingSyncmerHashes[j];
      reverseHash = seeding::rol(reverseHash, k) ^ startingSyncmerHashes[l - j - 1];
    }
    if (forwardHash != reverseHash) {
      seeding::uniqueKminmer_t uniqueKminmer{*indexingIt, refOnSyncmers[*curIt].value().endPos, std::min(forwardHash, reverseHash), reverseHash < forwardHash};

      if (refOnKminmers[*indexingIt].has_value()) {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, refOnKminmers[*indexingIt].value());
      } else {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
      }

      auto uniqueKminmerIndexIt = kminmerToUniqueIndex.find(uniqueKminmer);
      if (uniqueKminmerIndexIt == kminmerToUniqueIndex.end()) {
        uniqueKminmers.emplace_back(uniqueKminmer);
        kminmerToUniqueIndex[uniqueKminmer] = uniqueKminmers.size() - 1;
        addedSeedIndices.push_back(uniqueKminmers.size() - 1);
        refOnKminmers[*indexingIt] = uniqueKminmers.size() - 1;
      } else {
        addedSeedIndices.push_back(uniqueKminmerIndexIt->second);
        refOnKminmers[*indexingIt] = uniqueKminmerIndexIt->second;
      }
      
    } else {
      if (refOnKminmers[*indexingIt].has_value()) {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, refOnKminmers[*indexingIt].value());
        refOnKminmers[*indexingIt] = std::nullopt;
        deletedSeedIndices.push_back(*indexingIt);
      }
    }

    if (curIt == endIt) continue;
    
    while (curIt != endIt) {
      ++curIt;
      if (curIt == refOnSyncmersMap.end()) break;
      forwardHash = seeding::rol(forwardHash, k) ^ seeding::rol(refOnSyncmers[*indexingIt].value().hash, k * l) ^ refOnSyncmers[*curIt].value().hash;
      reverseHash = seeding::ror(reverseHash, k) ^ seeding::ror(refOnSyncmers[*indexingIt].value().hash, k)     ^ seeding::rol(refOnSyncmers[*curIt].value().hash, k * (l-1));
      ++indexingIt;

      if (forwardHash != reverseHash) {
        seeding::uniqueKminmer_t uniqueKminmer{*indexingIt, refOnSyncmers[*curIt].value().endPos, std::min(forwardHash, reverseHash), reverseHash < forwardHash};
        if (refOnKminmers[*indexingIt].has_value()) {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, refOnKminmers[*indexingIt].value());
        } else {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
        }

        auto uniqueKminmerIndexIt = kminmerToUniqueIndex.find(uniqueKminmer);
        if (uniqueKminmerIndexIt == kminmerToUniqueIndex.end()) {
          uniqueKminmers.emplace_back(uniqueKminmer);
          kminmerToUniqueIndex[uniqueKminmer] = uniqueKminmers.size() - 1;
          addedSeedIndices.push_back(uniqueKminmers.size() - 1);
          refOnKminmers[*indexingIt] = uniqueKminmers.size() - 1;
        } else {
          addedSeedIndices.push_back(uniqueKminmerIndexIt->second);
          refOnKminmers[*indexingIt] = uniqueKminmerIndexIt->second;
        }
      } else {
        if (refOnKminmers[*indexingIt].has_value()) {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, refOnKminmers[*indexingIt].value());
          refOnKminmers[*indexingIt] = std::nullopt;
          deletedSeedIndices.push_back(*indexingIt);
        }
      }
    }
  }
  if (!newKminmerRanges.empty() && newKminmerRanges.back().second == refOnSyncmersMap.end()) {
    auto delIt = newKminmerRanges.back().second;
    for (size_t j = 0; j < l - 1; j++) {
      --delIt;
      if (refOnKminmers[*delIt].has_value()) {
        refOnKminmersChangeRecord.emplace_back(*delIt, panmapUtils::seedChangeType::DEL, refOnKminmers[*delIt].value());
        refOnKminmers[*delIt] = std::nullopt;
        deletedSeedIndices.push_back(*delIt);
      }
      if (delIt == refOnSyncmersMap.begin()) break;
    }
  }
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL && refOnKminmers[syncmerPos].has_value()) {
      refOnKminmersChangeRecord.emplace_back(syncmerPos, panmapUtils::seedChangeType::DEL, refOnKminmers[syncmerPos].value());
      refOnKminmers[syncmerPos] = std::nullopt;
      deletedSeedIndices.push_back(syncmerPos);
    }
  }


  
  //  Adding node changes to index
  NodeChanges::Builder curNodeChanges = perNodeChanges[dfsIndex];
  curNodeChanges.setNodeIndex(dfsIndex);

  // adding inserted/substituted seeds to index
  capnp::List<uint32_t>::Builder seedInsubBuilder = curNodeChanges.initSeedInsubIndices(addedSeedIndices.size());
  for (size_t i = 0; i < addedSeedIndices.size(); i++) {
    seedInsubBuilder.set(i, addedSeedIndices[i]);
  }

  // adding deleted seeds to index
  capnp::List<uint32_t>::Builder seedDeletionBuilder = curNodeChanges.initSeedDeletions(deletedSeedIndices.size());
  for (size_t i = 0; i < deletedSeedIndices.size(); i++) {
    seedDeletionBuilder.set(i, deletedSeedIndices[i]);
  }

  // adding coord deltas to index
  capnp::List<CoordDelta>::Builder coordDeltaBuilder = curNodeChanges.initCoordDeltas(gapMapUpdates.size());
  for (size_t i = 0; i < gapMapUpdates.size(); i++) {
    const auto& [del, range] = gapMapUpdates[i];
    coordDeltaBuilder[i].setPos(range.first);
    if (del) {
      coordDeltaBuilder[i].getEndPos().setNone();
    } else {
      coordDeltaBuilder[i].getEndPos().setValue(range.second);
    }
  }

  // adding inverted blocks to index
  capnp::List<uint32_t>::Builder invertedBlocksBuilder = curNodeChanges.initInvertedBlocks(invertedBlocksVec.size());
  for (size_t i = 0; i < invertedBlocksVec.size(); i++) {
    invertedBlocksBuilder.set(i, invertedBlocksVec[i]);
  }
  
  // compare with brute force for debugging
  // compareBruteForceBuild(T, node, blockSequences, globalCoords, gapMap, degapCoordIndex, regapCoordIndex, refOnSyncmers, refOnSyncmersMap, blockOnSyncmers, refOnKminmers, uniqueKminmers, kminmerToUniqueIndex, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getT(), indexBuilder.getL(), indexBuilder.getOpen());


  revertGapMapInversions(gapRunBlockInversionBacktracks, gapMap);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapRunBlockInversionBacktracks); // gapRunBlockInversionBacktracks is no longer needed... clear memory

  // update delayed block states
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExistsDelayed[blockId] = newExists;
    blockStrandDelayed[blockId] = newStrand;
  }

  std::cout << "\rdfsIndex: " << dfsIndex << std::flush;
  for (panmanUtils::Node *child : node->children) {
    dfsIndex++;
    buildIndexHelper(child, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks, dfsIndex);
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

  for (const auto& [pos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      // Was added... need to delete
      refOnSyncmers[pos] = std::nullopt;
      refOnSyncmersMap.erase(pos);
    } else {
      // was deleted or replaced... need to restore
      refOnSyncmers[pos] = rsyncmer;
      refOnSyncmersMap.insert(pos);
    }
  }
  
  for (const auto& [blockId, pos, changeType] : blockOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      // was added... need to delete
      blockOnSyncmers[blockId].erase(pos);
      if (blockOnSyncmers[blockId].empty()) {
        blockOnSyncmers.erase(blockId);
      }
    } else {
      blockOnSyncmers[blockId].insert(pos);
    }
  }

  for (const auto& [pos, changeType, kminmerIndex] : refOnKminmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      refOnKminmers[pos] = std::nullopt;
    } else {
      refOnKminmers[pos] = kminmerIndex;
    }
  }
}

void mgsr::mgsrIndexBuilder::buildIndex() {
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences.sequence);
  refOnSyncmers.resize(globalCoords.lastScalarCoord + 1);
  refOnKminmers.resize(globalCoords.lastScalarCoord + 1);

  std::vector<bool> blockExistsDelayed = blockSequences.blockExists;
  std::vector<bool> blockStrandDelayed = blockSequences.blockStrand;

  std::map<uint64_t, uint64_t> gapMap{{0, globalCoords.lastScalarCoord}};
  std::unordered_set<uint64_t> invertedBlocks;

  uint64_t dfsIndex = 0;
  buildIndexHelper(T->root, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks, dfsIndex);
  
  // Finally add unique k-min-mers to index
  capnp::List<SeedInfo>::Builder seedInfoBuilder = indexBuilder.initSeedInfo(uniqueKminmers.size());
  for (size_t i = 0; i < uniqueKminmers.size(); i++) {
    seedInfoBuilder[i].setHash(uniqueKminmers[i].hash);
    seedInfoBuilder[i].setStartPos(uniqueKminmers[i].startPos);
    seedInfoBuilder[i].setEndPos(uniqueKminmers[i].endPos);
    seedInfoBuilder[i].setIsReverse(uniqueKminmers[i].isReverse);
  }

  std::cout << "Finished building index!" << std::endl;
}

void mgsr::mgsrIndexBuilder::writeIndex(const std::string& path) {
  int fd = open(path.c_str(), O_RDWR | O_CREAT, 0644);
  if (fd == -1) {
    std::cerr << "Error: failed to open file " << path << std::endl;
    std::exit(1);
  }
  std::cout << "Writing index to " << path << std::endl;
  capnp::writePackedMessageToFd(fd, outMessage);
  close(fd);
  std::cout << "Index written to " << path << std::endl;
}

int mgsr::mgsrIndexReader::open_file(const std::string& path) {
  int fd = open(path.c_str(), O_RDONLY);
  if (fd == -1) {
    std::cerr << "Error: failed to open file " << path << std::endl;
    std::exit(1);
  }
  return fd;
}

void mgsr::mgsrPlacer::addSeedAtPosition(uint64_t uniqueKminmerIndex, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks) {
  uint32_t pos = indexReader.seedInfos[uniqueKminmerIndex].getStartPos();
  auto posMapIt = positionMap.find(pos);
  if (posMapIt != positionMap.end()) {
    seedBacktracks.emplace_back(posMapIt->second, panmapUtils::seedChangeType::SUB);
    size_t oldHash = indexReader.seedInfos[posMapIt->second].getHash();
    size_t newHash = indexReader.seedInfos[uniqueKminmerIndex].getHash();
    auto oldHashToPositionIt = hashToPositionMap.find(oldHash);
    oldHashToPositionIt->second.erase(posMapIt);
    if (oldHashToPositionIt->second.empty()) hashToPositionMap.erase(oldHashToPositionIt);
    posMapIt->second = uniqueKminmerIndex;
    hashToPositionMap[newHash].insert(posMapIt);
  } else {
    posMapIt = positionMap.emplace(pos, uniqueKminmerIndex).first;
    size_t hash = indexReader.seedInfos[uniqueKminmerIndex].getHash();
    hashToPositionMap[hash].insert(posMapIt);
    seedBacktracks.emplace_back(uniqueKminmerIndex, panmapUtils::seedChangeType::ADD);
  }
}

void mgsr::mgsrPlacer::addSeedAtPosition(uint64_t uniqueKminmerIndex) {
  uint32_t pos = indexReader.seedInfos[uniqueKminmerIndex].getStartPos();
  auto posMapIt = positionMap.find(pos);
  if (posMapIt != positionMap.end()) {
    size_t oldHash = indexReader.seedInfos[posMapIt->second].getHash();
    size_t newHash = indexReader.seedInfos[uniqueKminmerIndex].getHash();
    auto oldHashToPositionIt = hashToPositionMap.find(oldHash);
    oldHashToPositionIt->second.erase(posMapIt);
    if (oldHashToPositionIt->second.empty()) hashToPositionMap.erase(oldHashToPositionIt);
    posMapIt->second = uniqueKminmerIndex;
    hashToPositionMap[newHash].insert(posMapIt);
  } else {
    posMapIt = positionMap.emplace(pos, uniqueKminmerIndex).first;
    size_t hash = indexReader.seedInfos[uniqueKminmerIndex].getHash();
    hashToPositionMap[hash].insert(posMapIt);
  }
}


void mgsr::mgsrPlacer::delSeedAtPosition(uint64_t pos, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks) {
  auto posMapIt = positionMap.find(pos);
  seedBacktracks.emplace_back(posMapIt->second, panmapUtils::seedChangeType::DEL);
  size_t hash = indexReader.seedInfos[posMapIt->second].getHash();
  auto hashToPositionIt = hashToPositionMap.find(hash);
  hashToPositionIt->second.erase(posMapIt);
  if (hashToPositionIt->second.empty()) hashToPositionMap.erase(hashToPositionIt);
  positionMap.erase(posMapIt);
}

void mgsr::mgsrPlacer::delSeedAtPosition(uint64_t pos) {
  auto posMapIt = positionMap.find(pos);
  size_t hash = indexReader.seedInfos[posMapIt->second].getHash();
  auto hashToPositionIt = hashToPositionMap.find(hash);
  hashToPositionIt->second.erase(posMapIt);
  if (hashToPositionIt->second.empty()) hashToPositionMap.erase(hashToPositionIt);
  positionMap.erase(posMapIt);
}

void mgsr::mgsrPlacer::updateSeeds(uint64_t currentDfsIndex, std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>>& seedBacktracks) {
  for (uint32_t seedInsubSeedIndex : indexReader.perNodeChanges[currentDfsIndex].getSeedInsubIndices()) {
    addSeedAtPosition(seedInsubSeedIndex, seedBacktracks);
  }

  for (uint32_t deletedPos : indexReader.perNodeChanges[currentDfsIndex].getSeedDeletions()) {
    delSeedAtPosition(deletedPos, seedBacktracks);
  }
}

void mgsr::mgsrPlacer::updateGapMap(uint64_t currentDfsIndex, const panmapUtils::GlobalCoords& globalCoords, std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBacktracks, std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBlocksBacktracks) {
  for (const auto& coordDelta : indexReader.perNodeChanges[currentDfsIndex].getCoordDeltas()) {
    const uint32_t& pos = coordDelta.getPos();
    const auto& endPos = coordDelta.getEndPos();
    switch (endPos.which()) {
      case CoordDelta::EndPos::VALUE: {
        auto gapMapIt = gapMap.find(pos);
        if (gapMapIt != gapMap.end()) {
          gapMapBacktracks.emplace_back(false, std::make_pair(pos, gapMapIt->second));
          gapMapIt->second = endPos.getValue();
        } else {
          gapMapBacktracks.emplace_back(true, std::make_pair(pos, endPos.getValue()));
          gapMap[pos] = endPos.getValue();
        }
        break;
      }
      case CoordDelta::EndPos::NONE: {
        gapMapBacktracks.emplace_back(false, std::make_pair(pos, gapMap[pos]));
        gapMap.erase(pos);
        break;
      }
      default:
        std::cerr << "Error: unknown endPos type" << std::endl;
        std::exit(1);
    }
  }

  for (const auto& invertedBlock : indexReader.perNodeChanges[currentDfsIndex].getInvertedBlocks()) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(invertedBlock);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(invertedBlock);
    std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> tmpGapMapUpdates;
    invertGapMap(gapMap, {beg, end}, gapMapBlocksBacktracks, tmpGapMapUpdates);
  }
}


void mgsr::mgsrPlacer::placeReadsHelper(panmanUtils::Node* node, uint64_t& dfsIndex, const panmapUtils::GlobalCoords& globalCoords) {
  std::vector<std::pair<uint64_t, panmapUtils::seedChangeType>> seedBacktracks;
  updateSeeds(dfsIndex, seedBacktracks);

  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapBacktracks;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapBlocksBacktracks;
  updateGapMap(dfsIndex, globalCoords, gapMapBacktracks, gapMapBlocksBacktracks);

  std::map<uint64_t, uint64_t> degapCoordIndex;
  std::map<uint64_t, uint64_t> regapCoordIndex;
  makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, (uint64_t)globalCoords.lastScalarCoord);

  // compareBruteForcePlace(T, node, globalCoords, gapMap, degapCoordIndex, regapCoordIndex, positionMap, hashToPositionMap, indexReader.seedInfos, indexReader.indexReader.getK(), indexReader.indexReader.getS(), indexReader.indexReader.getT(), indexReader.indexReader.getL(), indexReader.indexReader.getOpen());

  revertGapMapInversions(gapMapBlocksBacktracks, gapMap);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapMapBlocksBacktracks); // gapMapBlocksBacktracks is no longer needed... clear memory
  
  std::cout << "\rdfsIndex: " << dfsIndex << std::flush;
  for (panmanUtils::Node *child : node->children) {
    dfsIndex++;
    placeReadsHelper(child, dfsIndex, globalCoords);
  }

  // backtrack
  for (const auto& [uniqueKminmerIndex, changeType] : seedBacktracks) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      delSeedAtPosition(indexReader.seedInfos[uniqueKminmerIndex].getStartPos());
    } else {
      addSeedAtPosition(uniqueKminmerIndex);
    }
  }

  for (auto it = gapMapBacktracks.rbegin(); it != gapMapBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }
}

void mgsr::mgsrPlacer::placeReads() {
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences.sequence);

  gapMap.clear();
  gapMap.insert(std::make_pair(0, globalCoords.lastScalarCoord));

  uint64_t dfsIndex = 0;
  placeReadsHelper(T->root, dfsIndex, globalCoords);
}
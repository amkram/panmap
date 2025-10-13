#pragma once

#include "mgsr.hpp"
#include "panmap_utils.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_for.h>
#include <tbb/global_control.h>
#include <random>
#include <numeric>
#include <algorithm>
#include <stack>
#include <bitset>

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
  std::vector<char> blockExistsBruteForce;
  std::vector<char> blockStrandBruteForce;
  std::unordered_map<int, int> blockLengthsBruteForce;
  panmapUtils::getSequenceFromReference(T, sequenceBruteForce, blockExistsBruteForce, blockStrandBruteForce, blockLengthsBruteForce, node->identifier);

  std::string gappedSequenceBruteForce = panmapUtils::getStringFromSequence(sequenceBruteForce, blockLengthsBruteForce, blockExistsBruteForce, blockStrandBruteForce, true);
  std::vector<std::pair<uint64_t, uint64_t>> gapMapBruteForce;
  // checking gap map with brute force
  uint64_t curScalar = 0;
  while (curScalar < gappedSequenceBruteForce.size()){
    const auto& curCoord = globalCoords.getCoordFromScalar(curScalar);
    if (curCoord == globalCoords.blockEdgeCoords[curCoord.primaryBlockId].start && !blockSequences.blockExists[curCoord.primaryBlockId]) {
      const auto curBlockLength = blockLengthsBruteForce[curCoord.primaryBlockId];
      if (!gapMapBruteForce.empty() && gapMapBruteForce.back().second + 1 == curScalar) {
        gapMapBruteForce.back().second += curBlockLength;
      } else {
        gapMapBruteForce.emplace_back(curScalar, curScalar + curBlockLength - 1);
      }
      curScalar += curBlockLength;
      continue;
    }
    char nuc = gappedSequenceBruteForce[curScalar];
    if (nuc == '-') {
      if (!gapMapBruteForce.empty() && gapMapBruteForce.back().second + 1 == curScalar) {
        ++gapMapBruteForce.back().second;
      } else {
        gapMapBruteForce.emplace_back(curScalar, curScalar);
      }
    }
    curScalar++;
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

  std::cout << "GapMap passed... " << std::flush;


  // // check block sequence objects and coordinates
  // const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequenceDynamic = blockSequences.sequence;
  // const std::vector<char>& blockExistsDynamic = blockSequences.blockExists;
  // const std::vector<char>& blockStrandDynamic = blockSequences.blockStrand;
  // if (sequenceDynamic.size() != sequenceBruteForce.size()) {
  //   std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic.size() << " != brute force " << sequenceBruteForce.size() << std::endl;
  //   std::exit(1);
  // } else {
  //   if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical sequence size... passed: " << sequenceDynamic.size() << " == " << sequenceBruteForce.size() << std::endl;
  // }

  // uint64_t localScalarCoordBruteForce = 0;
  // uint64_t globalScalarCoord = 0;
  // for (int blockId = 0; blockId < blockSequences.numBlocks(); blockId++) {
  //   if (!blockExistsDynamic[blockId]) {
  //     globalScalarCoord += blockLengthsBruteForce[blockId];
  //     continue;
  //   }

  //   if (globalScalarCoord != globalCoords.getBlockStartScalar(blockId)) {
  //     std::cerr << "Global scalar coord block " << blockId << " start mismatch: dynamic " << globalScalarCoord << " != brute force " << globalCoords.getBlockStartScalar(blockId) << std::endl;
  //     std::exit(1);
  //   } else {
  //     if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " start scalar coord... passed: " << globalScalarCoord << " == " << globalCoords.getBlockStartScalar(blockId) << std::endl;
  //   }

  //   if (blockSequences.blockExists[blockId] != blockExistsBruteForce[blockId]) {
  //     std::cerr << "Block " << blockId << " exists state mismatch: dynamic " << blockSequences.blockExists[blockId] << " != brute force " << blockExistsBruteForce[blockId] << std::endl;
  //     std::exit(1);
  //   } else {
  //     if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " exists state... passed: " << blockSequences.blockExists[blockId] << " == " << blockExistsBruteForce[blockId] << std::endl;
  //   }

  //   if (blockSequences.blockStrand[blockId] != blockStrandBruteForce[blockId]) {
  //     std::cerr << "Block " << blockId << " strand state mismatch: dynamic " << blockSequences.blockStrand[blockId] << " != brute force " << blockStrandBruteForce[blockId] << std::endl;
  //     std::exit(1);
  //   } else {
  //     if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " strand state... passed: " << blockSequences.blockStrand[blockId] << " == " << blockStrandBruteForce[blockId] << std::endl;
  //   }

  //   if (blockStrandDynamic[blockId]) {
  //     for (int i = 0; i < sequenceDynamic[blockId].size(); i++) {      
  //       if (sequenceDynamic[blockId][i].second.size() != sequenceBruteForce[blockId][i].second.size()) {
  //         std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic[blockId][i].second.size() << " != brute force " << sequenceBruteForce[blockId][i].second.size() << std::endl;
  //         std::exit(1);
  //       } else {
  //         if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc size at (" << blockId << ", " << i << ")... passed: " << sequenceDynamic[blockId][i].second.size() << " == " << sequenceBruteForce[blockId][i].second.size() << std::endl;
  //       }
  
  //       for (int j = 0; j < sequenceDynamic[blockId][i].second.size(); j++) {
  //         if (sequenceDynamic[blockId][i].second[j] != sequenceBruteForce[blockId][i].second[j]) {
  //           std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", " << j << "): dynamic " << sequenceDynamic[blockId][i].second[j] << " != brute force " << sequenceBruteForce[blockId][i].second[j] << std::endl;
  //           std::exit(1);
  //         } else {
  //           if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc at (" << blockId << ", " << i << ", " << j << ")... passed: " << sequenceDynamic[blockId][i].second[j] << " == " << sequenceBruteForce[blockId][i].second[j] << std::endl;
  //         }
  
  //         if (sequenceDynamic[blockId][i].second[j] != '-') {
  //           if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
  //             std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") and global scalar coord  " << globalScalarCoord << std::endl;
  //             std::cerr << "Nuc: " << sequenceDynamic[blockId][i].second[j] << " ?= " << sequenceBruteForce[blockId][i].second[j] << std::endl;
  //             std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
  //             std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
  //             std::exit(1);
  //           } else {
  //             if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
  //           }
  //           ++localScalarCoordBruteForce;
  //         }
  //         ++globalScalarCoord;
  //       }
  
  //       if (sequenceDynamic[blockId][i].first != sequenceBruteForce[blockId][i].first) {
  //         std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", -1): dynamic " << sequenceDynamic[blockId][i].first << " != brute force " << sequenceBruteForce[blockId][i].first << std::endl;
  //         std::exit(1);
  //       } else {
  //         if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical main nuc at (" << blockId << ", " << i << ", -1)... passed: " << sequenceDynamic[blockId][i].first << " == " << sequenceBruteForce[blockId][i].first << std::endl;
  //       }
  
  //       if (sequenceDynamic[blockId][i].first != 'x') {
  //         if (sequenceDynamic[blockId][i].first != '-') {
  //           if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
  //             std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) and global scalar coord " << globalScalarCoord << std::endl;
  //             std::cerr << "Nuc: " << sequenceDynamic[blockId][i].first << " ?= " << sequenceBruteForce[blockId][i].first << std::endl;
  //             std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
  //             std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
  //             std::exit(1);
  //           } else {
  //             if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
  //           }
  //           ++localScalarCoordBruteForce;
  //         }
  //         ++globalScalarCoord;
  //       }
  //     }
  //   } else {
  //     for (int i = sequenceDynamic[blockId].size() - 1; i >= 0; i--) {
  //       if (sequenceDynamic[blockId][i].first != sequenceBruteForce[blockId][i].first) {
  //         std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", -1): dynamic " << sequenceDynamic[blockId][i].first << " != brute force " << sequenceBruteForce[blockId][i].first << std::endl;
  //         std::exit(1);
  //       } else {
  //         if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical main nuc at (" << blockId << ", " << i << ", -1)... passed: " << sequenceDynamic[blockId][i].first << " == " << sequenceBruteForce[blockId][i].first << std::endl;
  //       }
  
  //       if (sequenceDynamic[blockId][i].first != 'x') {
  //         if (sequenceDynamic[blockId][i].first != '-') {
  //           if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
  //             std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) and global scalar coord " << globalScalarCoord << std::endl;
  //             std::cerr << "Nuc: " << sequenceDynamic[blockId][i].first << " ?= " << sequenceBruteForce[blockId][i].first << std::endl;
  //             std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
  //             std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
  //             std::exit(1);
  //           } else {
  //             if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
  //           }
  //           ++localScalarCoordBruteForce;
  //         }
  //         ++globalScalarCoord;
  //       }

  //       if (sequenceDynamic[blockId][i].second.size() != sequenceBruteForce[blockId][i].second.size()) {
  //         std::cerr << "Sequence size mismatch: dynamic " << sequenceDynamic[blockId][i].second.size() << " != brute force " << sequenceBruteForce[blockId][i].second.size() << std::endl;
  //         std::exit(1);
  //       } else {
  //         if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc size at (" << blockId << ", " << i << ")... passed: " << sequenceDynamic[blockId][i].second.size() << " == " << sequenceBruteForce[blockId][i].second.size() << std::endl;
  //       }
  
  //       for (int j = sequenceDynamic[blockId][i].second.size() - 1; j >= 0; j--) {
  //         if (sequenceDynamic[blockId][i].second[j] != sequenceBruteForce[blockId][i].second[j]) {
  //           std::cerr << "Nuc mismatch at coord (" << blockId << ", " << i << ", " << j << "): dynamic " << sequenceDynamic[blockId][i].second[j] << " != brute force " << sequenceBruteForce[blockId][i].second[j] << std::endl;
  //           std::exit(1);
  //         } else {
  //           if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical gap nuc at (" << blockId << ", " << i << ", " << j << ")... passed: " << sequenceDynamic[blockId][i].second[j] << " == " << sequenceBruteForce[blockId][i].second[j] << std::endl;
  //         }
  
  //         if (sequenceDynamic[blockId][i].second[j] != '-') {
  //           if (mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
  //             std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") and global scalar coord " << globalScalarCoord << std::endl;
  //             std::cerr << "Nuc: " << sequenceDynamic[blockId][i].second[j] << " ?= " << sequenceBruteForce[blockId][i].second[j] << std::endl;
  //             std::cerr << "Degapped scalar coord: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
  //             std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
  //             std::exit(1);
  //           } else {
  //             if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed: " << mgsr::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
  //           }
  //           ++localScalarCoordBruteForce;
  //         }
  //         ++globalScalarCoord;
  //       }

  //     }
  //   }

  //   if (globalScalarCoord - 1 != globalCoords.getBlockEndScalar(blockId)) {
  //     std::cerr << "Global scalar coord block " << blockId << " end mismatch: dynamic " << globalScalarCoord - 1 << " != brute force " << globalCoords.getBlockEndScalar(blockId) << std::endl;
  //     std::exit(1);
  //   } else {
  //     if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical block " << blockId << " end scalar coord... passed: " << globalScalarCoord - 1 << " == " << globalCoords.getBlockEndScalar(blockId) << std::endl;
  //   }
  // }

  // std::cout << "sequence and coordinate objects passed... " << std::flush;

  // // check syncmers
  // std::string ungappedSequence = panmapUtils::getStringFromSequence(sequenceBruteForce, blockLengthsBruteForce, blockExistsBruteForce, blockStrandBruteForce, false);
  // std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersBruteForce = seeding::rollingSyncmers(ungappedSequence, k, s, open, t, false);
  // std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersDynamic;
  // for (size_t i = 0; i < refOnSyncmers.size(); i++) {
  //   if (refOnSyncmers[i].has_value()) {
  //     const auto& [hash, endPos, isReverse] = refOnSyncmers[i].value();
  //     syncmersDynamic.emplace_back(std::make_tuple(hash, isReverse, true, mgsr::degapGlobal(i, degapCoordIndex)));
  //   }
  // }



  // // check all syncmers
  // if (syncmersDynamic.size() != syncmersBruteForce.size()) {
  //   std::cout << "Syncmer count mismatch: dynamic " << syncmersDynamic.size() << " != brute force " << syncmersBruteForce.size() << std::endl;
  //   std::cout << "Dynamic syncmers: ";
  //   for (const auto& syncmer : syncmersDynamic) {
  //     std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << mgsr::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << "Brute force syncmers: ";
  //   for (const auto& syncmer : syncmersBruteForce) {
  //     std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << mgsr::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
  //   }
  //   std::cout << std::endl;
  //   std::exit(1);
  // } else if (printCorrect && node->identifier == nodeToDebug) {
  //   std::cout << "Identical syncmer count... passed: " << syncmersDynamic.size() << " == " << syncmersBruteForce.size() << std::endl;
  // }
  
  // auto curSyncmerOnMapIt = refOnSyncmersMap.begin();
  // for (size_t i = 0; i < syncmersDynamic.size(); i++) {
  //   const auto& [hash, isReverse, isSeed, startPos] = syncmersDynamic[i];
  //   const auto& [hashBruteForce, isReverseBruteForce, isSeedBruteForce, startPosBruteForce] = syncmersBruteForce[i];
  //   if (hash != hashBruteForce || isReverse != isReverseBruteForce || startPos != startPosBruteForce) {
  //     std::cout << "Syncmer mismatch at " << i << "th syncmer: dynamic (" << hash << ", " << startPos << ", " << isReverse << ") != brute force (" << hashBruteForce << ", " << startPosBruteForce << ", " << isReverseBruteForce << ")" << std::endl;
  //     std::cout << "Dynamic syncmers: ";
  //     for (const auto& syncmer : syncmersDynamic) {
  //       std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << mgsr::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
  //     }
  //     std::cout << std::endl;
  //     std::cout << "Brute force syncmers: ";
  //     for (const auto& syncmer : syncmersBruteForce) {
  //       std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << mgsr::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
  //     }
  //     std::cout << std::endl;
  //     std::exit(1);
  //   }
  //   if (mgsr::regapGlobal(startPos, regapCoordIndex) != *curSyncmerOnMapIt) {
  //     std::cout << "Syncmer on map mismatch at " << i << "th syncmer: dynamic/bruteforce " << startPos << " != map " << *curSyncmerOnMapIt << std::endl;
  //     std::exit(1);
  //   }
  //   ++curSyncmerOnMapIt;
  // }
  // if (curSyncmerOnMapIt != refOnSyncmersMap.end()) {
  //   std::cout << "SyncmerOnMap has more elements than syncmers: " << refOnSyncmersMap.size() << " != " << syncmersDynamic.size() << std::endl;
  //   std::exit(1);
  // }
  // std::cout << "syncmers passed... " << std::flush;

  // // check k-min-mers
  // std::vector<std::tuple<size_t, size_t, size_t, bool>> kminmersBruteForce;
  // if (syncmersBruteForce.size() >= l) {
  //   size_t forwardRolledHash = 0;
  //   size_t reverseRolledHash = 0;
  
  //   // first kminmer
  //   for (size_t i = 0; i < l; ++i) {
  //     forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ std::get<0>(syncmersBruteForce[i]);
  //     reverseRolledHash = seeding::rol(reverseRolledHash, k) ^ std::get<0>(syncmersBruteForce[l-i-1]);
  //   }
  
  //   if (forwardRolledHash != reverseRolledHash) {
  //     size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
  //     kminmersBruteForce.emplace_back(minHash, std::get<3>(syncmersBruteForce[0]), std::get<3>(syncmersBruteForce[l-1])+k-1, reverseRolledHash < forwardRolledHash);
  //   }
  
  //   // rest of kminmers
  //   for (uint64_t i = 1; i < syncmersBruteForce.size()-l+1; ++i) {
  //     if (!std::get<2>(syncmersBruteForce[i-1]) || !std::get<2>(syncmersBruteForce[i+l-1])) {
  //       std::cout << "invalid syncmer" << std::endl;
  //       exit(0);
  //     }
  //     const size_t& prevSyncmerHash = std::get<0>(syncmersBruteForce[i-1]);
  //     const size_t& nextSyncmerHash = std::get<0>(syncmersBruteForce[i+l-1]);
  //     forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ seeding::rol(prevSyncmerHash, k * l) ^ nextSyncmerHash;
  //     reverseRolledHash = seeding::ror(reverseRolledHash, k) ^ seeding::ror(prevSyncmerHash, k)     ^ seeding::rol(nextSyncmerHash, k * (l-1));
  
  //     if (forwardRolledHash != reverseRolledHash) {
  //       size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
  //       kminmersBruteForce.emplace_back(minHash, std::get<3>(syncmersBruteForce[i]), std::get<3>(syncmersBruteForce[i+l-1])+k-1, reverseRolledHash < forwardRolledHash);
  //     }
  //   }
  // }


  // std::vector<std::tuple<size_t, size_t, size_t, bool>> kminmersDynamic;
  // for (size_t i = 0; i < refOnKminmers.size(); i++) {
  //   if (refOnKminmers[i].has_value()) {
  //     const auto& [startPos, endPos, hash, isReverse] = uniqueKminmers[refOnKminmers[i].value()];
  //     kminmersDynamic.emplace_back(std::make_tuple(hash, mgsr::degapGlobal(startPos, degapCoordIndex), mgsr::degapGlobal(endPos, degapCoordIndex), isReverse));
  //   }
  // }

  // if (kminmersDynamic.size() != kminmersBruteForce.size()) {
  //   std::cout << "K-min-mer count mismatch: dynamic " << kminmersDynamic.size() << " != brute force " << kminmersBruteForce.size() << std::endl;
  //   std::cout << "Dynamic k-min-mers: ";
  //   for (const auto& kminmer : kminmersDynamic) {
  //     std::cout << "(" << std::get<0>(kminmer) << ", " << std::get<1>(kminmer) << ", " << mgsr::regapGlobal(std::get<1>(kminmer), regapCoordIndex) << ", " << std::get<2>(kminmer) << ", " << mgsr::regapGlobal(std::get<2>(kminmer), regapCoordIndex) << ", " << std::get<3>(kminmer) << ") ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << "Brute force k-min-mers: ";
  //   for (const auto& kminmer : kminmersBruteForce) {
  //     std::cout << "(" << std::get<0>(kminmer) << ", " << std::get<1>(kminmer) << ", " << mgsr::regapGlobal(std::get<1>(kminmer), regapCoordIndex) << ", " << std::get<2>(kminmer) << ", " << mgsr::regapGlobal(std::get<2>(kminmer), regapCoordIndex) << ", " << std::get<3>(kminmer) << ") ";
  //   }
  //   std::cout << std::endl;
  //   std::exit(1);
  // } else if (printCorrect && node->identifier == nodeToDebug) {
  //   std::cout << "Identical k-min-mer count... passed: " << kminmersDynamic.size() << " == " << kminmersBruteForce.size() << std::endl;
  // }

  // for (size_t i = 0; i < kminmersDynamic.size(); i++) {
  //   const auto& [hash, startPos, endPos, isReverse] = kminmersDynamic[i];
  //   const auto& [hashBruteForce, startPosBruteForce, endPosBruteForce, isReverseBruteForce] = kminmersBruteForce[i];
  //   if (hash != hashBruteForce || startPos != startPosBruteForce || endPos != endPosBruteForce || isReverse != isReverseBruteForce) {
  //     std::cout << "K-min-mer mismatch at " << i << "th k-min-mer: dynamic (" << hash << ", " << startPos << ", " << endPos << ", " << isReverse << ") != brute force (" << hashBruteForce << ", " << startPosBruteForce << ", " << endPosBruteForce << ", " << isReverseBruteForce << ")" << std::endl;
  //     std::exit(1);
  //   }
  // }
  // std::cout << "k-min-mers passed... " << std::flush;

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
  const std::unordered_map<size_t, std::vector<std::map<uint64_t, uint64_t>::iterator>>& hashToPositionMap,
  const std::vector<seeding::uniqueKminmer_t>& seedInfos,
  int k, int s, int t, int l, bool open
) {
  std::cout << "Checking " << node->identifier << " states with brute force at placement... " << std::flush;
  // check sequence object
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequenceBruteForce;
  std::vector<char> blockExistsBruteForce;
  std::vector<char> blockStrandBruteForce;
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
      seedInfos[infoIndex].hash,
      mgsr::degapGlobal(seedInfos[infoIndex].startPos, degapCoordIndex),
      mgsr::degapGlobal(seedInfos[infoIndex].endPos, degapCoordIndex),
      seedInfos[infoIndex].isReverse
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

static void printTopN(const Eigen::VectorXd& myVector, size_t n) {
  std::vector<size_t> indices(myVector.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&myVector](size_t i, size_t j) {
    return myVector(i) > myVector(j);
  });
  for (size_t i = 0; i < std::min(indices.size(), n); ++i) {
    std::cout << myVector(indices[i]) << " ";
  }
  std::cout << std::endl;
}


void mgsr::MgsrLiteTree::cleanup() {
  for (auto& pair : allLiteNodes) {
    delete pair.second;
  }
  allLiteNodes.clear();
  blockScalarRanges.clear();
  root = nullptr;
}

uint32_t mgsr::MgsrLiteTree::getBlockStartScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].first;
}

uint32_t mgsr::MgsrLiteTree::getBlockEndScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].second;
}

void mgsr::MgsrLiteTree::initialize(MGSRIndex::Reader indexReader, size_t numThreads, bool lowMemory) {
  std::cerr << "Starting to initialize MgsrLiteTree from index..." << std::endl;
  this->numThreads = numThreads;
  this->lowMemory = lowMemory;

  k = indexReader.getK();
  s = indexReader.getS();
  t = indexReader.getT();
  l = indexReader.getL();
  openSyncmer = indexReader.getOpen();

  std::vector<std::vector<std::pair<uint32_t, bool>>> seedDeltas;
  std::vector<std::vector<std::tuple<uint32_t, uint32_t, bool>>> gapRunDeltas;
  std::vector<std::vector<uint32_t>> invertedBlocks;
  capnp::List<SeedInfo>::Reader seedInfosReader = indexReader.getSeedInfo();
  capnp::List<NodeChanges>::Reader perNodeChangesReader = indexReader.getPerNodeChanges();
  seedInfos.resize(seedInfosReader.size());

  size_t seedInfos_chunkSize = (seedInfosReader.size() + numThreads - 1) / numThreads;
  std::cerr << "Start to read in seedInfos" << std::endl;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, seedInfosReader.size(), seedInfos_chunkSize),
    [&](const tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i < range.end(); i++) {
        const auto& seedReader = seedInfosReader[i];
        auto& seed = seedInfos[i];
        seed.hash = seedReader.getHash();
        seed.startPos = seedReader.getStartPos();
        seed.endPos = seedReader.getEndPos();
        seed.isReverse = seedReader.getIsReverse();
      }
    }, tbb::simple_partitioner());
  std::cerr << "Successfully read in seedInfos, size: " << seedInfos.size() << std::endl;

  seedDeltas.resize(perNodeChangesReader.size());
  gapRunDeltas.resize(perNodeChangesReader.size());
  invertedBlocks.resize(perNodeChangesReader.size());
  std::cerr << "Starting to read in perNodeChanges" << std::endl;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, perNodeChangesReader.size()),
    [&](const tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i < range.end(); i++) {
        const auto& currentPerNodeChangeReader = perNodeChangesReader[i];
        auto& currentSeedDeltas = seedDeltas[i];
        auto& currentGapRunDeltas = gapRunDeltas[i];
        auto& currentInvertedBlocks = invertedBlocks[i];
        const auto& currentSeedDeltasReader = currentPerNodeChangeReader.getSeedDeltas();
        const auto& currentGapRunDeltasReader = currentPerNodeChangeReader.getGapRunDeltas();
        const auto& currentInvertedBlocksReader = currentPerNodeChangeReader.getInvertedBlocks();
        
        currentSeedDeltas.resize(currentSeedDeltasReader.size());
        currentGapRunDeltas.reserve(currentGapRunDeltasReader.size());
        currentInvertedBlocks.resize(currentInvertedBlocksReader.size());
        
        for (size_t j = 0; j < currentSeedDeltasReader.size(); j++) {
          currentSeedDeltas[j].first = currentSeedDeltasReader[j].getSeedIndex();
          currentSeedDeltas[j].second = currentSeedDeltasReader[j].getIsDeleted();
        }
        
        for (size_t j = 0; j < currentGapRunDeltasReader.size(); j++) {
          const auto& currentGapRunDeltaReader = currentGapRunDeltasReader[j];
          currentGapRunDeltas.push_back({currentGapRunDeltaReader.getStartPos(), 
                                          currentGapRunDeltaReader.getEndPos(), 
                                          currentGapRunDeltaReader.getToGap()});
        }
        
        for (size_t j = 0; j < currentInvertedBlocksReader.size(); j++) {
          currentInvertedBlocks[j] = currentInvertedBlocksReader[j];
        }
      }
    });
  std::cerr << "Successfully read in perNodeChanges, size: " << perNodeChangesReader.size() << std::endl;

  auto liteTreeReader = indexReader.getLiteTree();
  // initialize blockScalarRanges
  auto blockScalarRangesReader = liteTreeReader.getBlockRanges();
  blockScalarRanges.resize(blockScalarRangesReader.size());
  for (size_t i = 0; i < blockScalarRangesReader.size(); i++) {
    blockScalarRanges[i] = {blockScalarRangesReader[i].getRangeBeg(), blockScalarRangesReader[i].getRangeEnd()};
  }

  // initialize allLiteNodes
  auto liteNodesReader = liteTreeReader.getLiteNodes();
  for (size_t i = 0; i < liteNodesReader.size(); i++) {
    const auto liteNodeReader = liteNodesReader[i];
    const auto& nodeIdentifier = liteNodeReader.getId();
    const auto parentIndex = liteNodeReader.getParentIndex();
    auto [it, inserted] = allLiteNodes.emplace(nodeIdentifier, new MgsrLiteNode(nodeIdentifier, nullptr, {}, i));
    it->second->seedDeltas = std::move(seedDeltas[i]);
    it->second->gapRunDeltas = std::move(gapRunDeltas[i]);
    it->second->invertedBlocks = std::move(invertedBlocks[i]);
    if (lowMemory) {
      it->second->readScoreDeltasLowMemory.resize(numThreads);
    } else {
      it->second->readScoreDeltas.resize(numThreads);
    }
    if (i == 0) continue;
    const auto parentNodeReader = liteNodesReader[parentIndex];
    const auto& parentNodeId = parentNodeReader.getId();
    it->second->parent = allLiteNodes[parentNodeId];
    allLiteNodes[parentNodeId]->children.push_back(it->second);
  }

  root = allLiteNodes[liteNodesReader[0].getId()];
}

void mgsr::MgsrLiteTree::detachNode(mgsr::MgsrLiteNode* node) {
  if (node == root) return;
  auto& parent = node->parent;
  auto& children = node->children;
  parent->children.erase(std::remove(parent->children.begin(), parent->children.end(), node), parent->children.end());
  for (auto child : children) {
    child->parent = parent;
    parent->children.push_back(child);
  }

  detachedNodes.insert(node);
}

std::pair<std::unordered_map<mgsr::MgsrLiteNode*, mgsr::MgsrLiteNode*>, std::unordered_map<mgsr::MgsrLiteNode*, int>> mgsr::MgsrLiteTree::findClosestTargets(
  const mgsr::MgsrLiteTree& tree,
  const std::vector<MgsrLiteNode*>& targets
) {
  std::unordered_map<mgsr::MgsrLiteNode*, mgsr::MgsrLiteNode*> closestTarget;
  std::unordered_map<mgsr::MgsrLiteNode*, int> distance;
  
  std::queue<mgsr::MgsrLiteNode*> q;
  
  for (mgsr::MgsrLiteNode* target : targets) {
    closestTarget[target] = target;
    distance[target] = 0;
    q.push(target);
  }
  
  while (!q.empty()) {
    mgsr::MgsrLiteNode* u = q.front();
    q.pop();
    
    if (u->parent && distance.find(u->parent) == distance.end()) {
      int edgeWeight = u->seedDeltas.size();
      distance[u->parent] = distance[u] + edgeWeight;
      closestTarget[u->parent] = closestTarget[u];
      q.push(u->parent);
    }
    
    for (mgsr::MgsrLiteNode* child : u->children) {
      if (distance.find(child) == distance.end()) {
        int edgeWeight = child->seedDeltas.size();
        distance[child] = distance[u] + edgeWeight;
        closestTarget[child] = closestTarget[u];
        q.push(child);
      }
    }
  }
  
  return {closestTarget, distance};
}



int64_t mgsr::mgsrPlacer::getReadBruteForceScore(
  size_t readIndex, absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>& hashCoordInfoCacheTable
) {
  auto readCopy = reads[readIndex];
  initializeReadMinichains(readCopy);
  int32_t pseudoScore = getReadPseudoScore(readCopy);
  return pseudoScore;
}

static size_t getBruteForceKminmerOverlapCount(
  const std::unordered_map<size_t, std::vector<std::pair<uint32_t, std::vector<uint64_t>>>>& seedmerToReads,
  const std::unordered_map<size_t, std::vector<std::map<uint64_t, uint64_t>::iterator>>& hashToPositionMap
) {
  size_t binaryOverlapKminmerCount = 0;
  for (const auto& [hash, positionMap] : hashToPositionMap) {
    if (seedmerToReads.find(hash) != seedmerToReads.end()) {
     ++binaryOverlapKminmerCount;
    }
  }
  return binaryOverlapKminmerCount;
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
  std::vector<uint32_t>& potentialSyncmerDeletions,
  const std::vector<char>& oldBlockExists,
  const std::vector<char>& oldBlockStrand
) {
  std::vector<char>& blockExists = blockSequences.blockExists;
  std::vector<char>& blockStrand = blockSequences.blockStrand;
  std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;
  
  // process block mutations
  for (const auto& blockMutation : node->blockMutation) {
    const int32_t blockId = blockMutation.primaryBlockId;
    const bool isInsertion = blockMutation.blockMutInfo;
    const bool isInversion = blockMutation.inversion;
    const bool oldExists = blockExists[blockId];
    const bool oldStrand = blockStrand[blockId];

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
    blockMutationRecord.emplace_back(blockId, oldExists, oldStrand, blockExists[blockId], blockStrand[blockId]);


    const auto& curBlockEdgeCoords = globalCoords.blockEdgeCoords[blockId];
    if (blockStrand[blockId]) {
      // forward strand
      localMutationRanges.emplace_back(curBlockEdgeCoords.start, curBlockEdgeCoords.end);
    } else {
      // reversed strand
      localMutationRanges.emplace_back(curBlockEdgeCoords.end, curBlockEdgeCoords.start);
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
      const char oldNuc = blockSequences.getSequenceBase(pos);
      const int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
      const char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);

      if (oldNuc == newNuc) continue;
      blockSequences.setSequenceBase(pos, newNuc);
      nucMutationRecord.emplace_back(pos, oldNuc, newNuc);

      if (oldBlockExists[pos.primaryBlockId] && blockExists[pos.primaryBlockId]) {
        const int64_t scalarCoord = globalCoords.getScalarFromCoord(pos);
        if (newNuc == '-') {
          // nuc to gap
          if (!gapRunUpdates.empty() && gapRunUpdates.back().first == true && gapRunUpdates.back().second.second + 1 == scalarCoord) {
            ++(gapRunUpdates.back().second.second);
          }
          else {
            gapRunUpdates.emplace_back(true, std::make_pair(scalarCoord, scalarCoord)); 
          }
          if (blockExists[blockId] && oldBlockExists[blockId] && blockStrand[blockId] == oldBlockStrand[blockId]) {
            potentialSyncmerDeletions.push_back(globalCoords.getScalarFromCoord(pos, blockStrand[pos.primaryBlockId]));
          }
        } else if (oldNuc == '-') {
          // gap to nuc
          if (!gapRunUpdates.empty() && gapRunUpdates.back().first == false && gapRunUpdates.back().second.second + 1 == scalarCoord) {
            ++(gapRunUpdates.back().second.second);
          } else {
            gapRunUpdates.emplace_back(false, std::make_pair(scalarCoord, scalarCoord));
          }
        }
      }
    }
    if (lastOffset != -1 && blockExists[blockId] && oldBlockExists[blockId] && blockStrand[blockId] == oldBlockStrand[blockId]) {
      if (blockStrand[blockId]) {
        localMutationRanges.emplace_back(panmapUtils::Coordinate(nucMutation, 0), panmapUtils::Coordinate(nucMutation, lastOffset));
      } else {
        localMutationRanges.emplace_back(panmapUtils::Coordinate(nucMutation, lastOffset), panmapUtils::Coordinate(nucMutation, 0));
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
      panmapUtils::Coordinate coord = globalCoords.blockEdgeCoords[blockId].start;
      panmapUtils::Coordinate end = globalCoords.blockEdgeCoords[blockId].end;
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
        globalCoords.stepRightCoordinate(coord);
      }
      if (curNucRange.first != -1) {
        gapRunUpdates.emplace_back(false, std::make_pair((uint64_t)curNucRange.first, (uint64_t)curNucRange.second));
      }
    }
  }
}

template <typename T>
static void inline perfect_shuffle(std::vector<T>& v) {
  int n = v.size();

  std::vector<T> canvas(n);

  for (int i = 0; i < n / 2; i++) {
    canvas[i*2] = v[i];
    canvas[i*2+1] = v[i + n/2];
  }

  v = std::move(canvas);
}

void mgsr::extractReadSequences(const std::string& readPath1, const std::string& readPath2, std::vector<std::string>& readSequences) {
  FILE *fp;
  kseq_t *seq;
  fp = fopen(readPath1.c_str(), "r");
  if(!fp){
    std::cerr << "Error: File " << readPath1 << " not found" << std::endl;
    exit(0);
  }
  seq = kseq_init(fileno(fp));
  int line;
  while ((line = kseq_read(seq)) >= 0) {
    readSequences.push_back(seq->seq.s);
  }
  if (readPath2.size() > 0) {
    fp = fopen(readPath2.c_str(), "r");
    if(!fp){
      std::cerr << "Error: File " << readPath2 << " not found" << std::endl;
      exit(0);
    }
    seq = kseq_init(fileno(fp));

    line = 0;
    int forwardReads = readSequences.size();
    while ((line = kseq_read(seq)) >= 0) {
      readSequences.push_back(seq->seq.s);
    }

    if (readSequences.size() != forwardReads*2){
      std::cerr << "Error: File " << readPath2 << " does not contain the same number of reads as " << readPath1 << std::endl;
      exit(0);
    }
    
    //Shuffle reads together, so that pairs are next to eatch other
    perfect_shuffle(readSequences);
  }
}

mgsr::RDGNode* mgsr::ReadDebruijnGraph::makeNewNode(size_t hash) {
  auto newNode = std::make_unique<RDGNode>(hash);
  auto rawPtr = newNode.get();
  hashToNode[hash] = std::move(newNode);
  return rawPtr;
}

std::pair<mgsr::RDGNode*, bool> mgsr::ReadDebruijnGraph::tryMakeNewNode(size_t hash) {
  auto hashToNodeIt = hashToNode.find(hash);
  if (hashToNodeIt == hashToNode.end()) {
    return {makeNewNode(hash), true};
  }
  return {hashToNodeIt->second.get(), false};
}

void mgsr::ReadDebruijnGraph::linkNodes(RDGNode* node1, RDGNode* node2) {
  if (node1 == nullptr || node2 == nullptr) {
    return;
  }
  node1->neighbors.insert(node2);
  node2->neighbors.insert(node1);
}

void mgsr::ReadDebruijnGraph::identifyConnectedComponents() {
  std::unordered_set<RDGNode*> visited;
  for (const auto& [_, nodeUniquePtr] : hashToNode) {
    auto node = nodeUniquePtr.get();
    if (visited.find(node) != visited.end()) continue;
    std::vector<RDGNode*> currentComponent;
    searchComponentDFS(node, visited, currentComponent);
    std::sort(currentComponent.begin(), currentComponent.end(), [](const RDGNode* a, const RDGNode* b) {
      return a->hash < b->hash;
    });
    connectedComponents.push_back(std::move(currentComponent));
  }
}


void mgsr::ReadDebruijnGraph::searchComponentDFS(RDGNode* node, std::unordered_set<RDGNode*>& visited, std::vector<RDGNode*>& currentComponent) {
  std::stack<RDGNode*> stack;
  stack.push(node);
  while (!stack.empty()) {
    RDGNode* currentNode = stack.top();
    stack.pop();
    if (visited.find(currentNode) != visited.end()) continue;
    visited.insert(currentNode);
    currentComponent.push_back(currentNode);
    for (const auto& neighbor : currentNode->neighbors) {
      if (visited.find(neighbor) == visited.end()) {
        stack.push(neighbor);
      }
    }
  }
}

void mgsr::ReadDebruijnGraph::sortReads(
  std::vector<mgsr::Read>& reads,
  std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex
) {
  size_t curReadIndex = 0;
  std::vector<mgsr::Read> sortedReads(reads.size());
  std::vector<std::vector<size_t>> sortedReadSeedmersDuplicatesIndex(readSeedmersDuplicatesIndex.size());
  for (const auto& component : connectedComponents) {
    sortReadsHelper(component[0], curReadIndex, reads, sortedReads, readSeedmersDuplicatesIndex, sortedReadSeedmersDuplicatesIndex);
  }
  reads.swap(sortedReads);
  readSeedmersDuplicatesIndex.swap(sortedReadSeedmersDuplicatesIndex);
}

void mgsr::ReadDebruijnGraph::sortReadsHelper(
  RDGNode* node, size_t& curReadIndex,
  std::vector<mgsr::Read>& reads, std::vector<mgsr::Read>& sortedReads,
  std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex, std::vector<std::vector<size_t>>& sortedReadSeedmersDuplicatesIndex
) {
  std::stack<RDGNode*> stack;
  std::unordered_set<RDGNode*> visited;
  stack.push(node);
  while (!stack.empty()) {
    RDGNode* currentNode = stack.top();
    stack.pop();
    if (visited.find(currentNode) != visited.end()) continue;
    visited.insert(currentNode);
    for (const auto& readIndex : currentNode->readIndicesMid) {
      sortedReads[curReadIndex] = std::move(reads[readIndex]);
      sortedReadSeedmersDuplicatesIndex[curReadIndex] = std::move(readSeedmersDuplicatesIndex[readIndex]);
      ++curReadIndex;
    }
    for (const auto& neighbor : currentNode->neighbors) {
      if (visited.find(neighbor) == visited.end()) {
        stack.push(neighbor);
      }
    }
  }

  // std::queue<RDGNode*> queue;
  // std::unordered_set<RDGNode*> visited;
  // queue.push(node);
  // visited.insert(node);
  // while (!queue.empty()) {
  //   RDGNode* currentNode = queue.front();
  //   queue.pop();
  //   for (const auto& readIndex : currentNode->readIndicesMid) {
  //     sortedReads[curReadIndex] = std::move(reads[readIndex]);
  //     sortedReadSeedmersDuplicatesIndex[curReadIndex] = std::move(readSeedmersDuplicatesIndex[readIndex]);
  //     ++curReadIndex;
  //   }
  //   for (const auto& neighbor : currentNode->neighbors) {
  //     if (visited.find(neighbor) == visited.end()) {
  //       visited.insert(neighbor);
  //       queue.push(neighbor);
  //     }
  //   }
  // }
}

void mgsr::ReadDebruijnGraph::exportToGFA(const std::string& filename) {
    std::ofstream gfaFile(filename);
  if (!gfaFile.is_open()) {
    throw std::runtime_error("Cannot open GFA file for writing: " + filename);
  }
  
  // Header line
  gfaFile << "H\tVN:Z:1.0\n";
  
  // Export segments (nodes)
  for (const auto& [hash, node] : hashToNode) {
    gfaFile << "S\t" << hash << "\t" << node.get()->hash << "\n";
  }
  
  // Export links (edges)
  std::unordered_set<std::pair<size_t, size_t>, mgsr::PairHash> exportedEdges;
  
  for (const auto& [hash, node] : hashToNode) {
    for (const auto& neighbor : node->neighbors) {
      size_t neighborHash = neighbor->hash;
      
      // Avoid duplicate edges (undirected graph)
      std::pair<size_t, size_t> edge = {std::min(hash, neighborHash), std::max(hash, neighborHash)};
      if (exportedEdges.find(edge) != exportedEdges.end()) continue;
      exportedEdges.insert(edge);
      
      // Determine orientation and overlap
      int overlap = 1;
      gfaFile << "L\t" << hash << "\t+\t" << neighborHash << "\t+\t" << overlap << "M\n";
    }
  }
  
  // Export paths for each connected component
  for (size_t i = 0; i < connectedComponents.size(); ++i) {
    if (connectedComponents[i].empty()) continue;
    
    gfaFile << "P\tcomponent_" << i << "\t";
    for (size_t j = 0; j < connectedComponents[i].size(); ++j) {
      if (j > 0) gfaFile << ",";
      gfaFile << connectedComponents[i][j]->hash << "+";
    }
    gfaFile << "\t*\n";
  }
  
  gfaFile.close();
}

void mgsr::ReadDebruijnGraph::buildGraph(std::vector<Read>& reads) {
  for (size_t i = 0; i < reads.size(); ++i) {
    RDGNode* prevNode = nullptr;
    for (size_t j = 0; j < reads[i].seedmersList.size(); ++j) {
      const auto& seedmer = reads[i].seedmersList[j];
      auto [curNode, isNew] = tryMakeNewNode(seedmer.hash);

      linkNodes(prevNode, curNode);

      if (j == reads[i].seedmersList.size() / 2 + 1) {
        curNode->readIndicesMid.push_back(i);
      }
      curNode->readIndicesCovered.push_back(i);

      prevNode = curNode;
    } 
  }

  // find connected components
  identifyConnectedComponents();

  // sort components by size (largest first)
  std::sort(connectedComponents.begin(), connectedComponents.end(),
    [](const std::vector<RDGNode*>& a, const std::vector<RDGNode*>& b) {
      return a.size() > b.size();
  });
  
}

void mgsr::ThreadsManager::initializeMGSRIndex(MGSRIndex::Reader indexReader) {
  k = indexReader.getK();
  s = indexReader.getS();
  t = indexReader.getT();
  l = indexReader.getL();
  openSyncmer = indexReader.getOpen();
}

void mgsr::mgsrPlacer::initializeMGSRIndex(MGSRIndex::Reader indexReader) {
  k = indexReader.getK();
  s = indexReader.getS();
  t = indexReader.getT();
  l = indexReader.getL();
  openSyncmer = indexReader.getOpen();
}

void mgsr::ThreadsManager::initializeQueryData(std::span<const std::string> readSequences, bool fast_mode) {
  std::unordered_map<std::string_view, std::vector<size_t>> seqToIndex;
  for (size_t i = 0; i < readSequences.size(); ++i) {
    seqToIndex[readSequences[i]].push_back(i);
  }
  std::vector<std::pair<std::string_view, std::vector<size_t>>> seqToIndexVec(seqToIndex.size());
  size_t seqToIndexVecIndex = 0;
  for (auto& [seq, index] : seqToIndex) {
    seqToIndexVec[seqToIndexVecIndex].first = seq;
    seqToIndexVec[seqToIndexVecIndex].second = std::move(index);
    ++seqToIndexVecIndex;
  }

  // seedmers for each unique read sequence
  size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
  std::vector<mgsr::Read> uniqueReadSeedmers(seqToIndexVec.size());
  tbb::parallel_for(tbb::blocked_range<size_t>(0, seqToIndexVec.size(), seqToIndexVec.size() / num_cpus), [&](const tbb::blocked_range<size_t>& range){
    for (size_t i = range.begin(); i < range.end(); ++i) {
      const auto& seq = seqToIndexVec[i].first;
      const auto& syncmers = seeding::rollingSyncmers(seq, k, s, openSyncmer, t, false);
      mgsr::Read& curRead = uniqueReadSeedmers[i];
      if (syncmers.size() < l) continue;

      size_t forwardRolledHash = 0;
      size_t reverseRolledHash = 0;
      // first kminmer
      for (size_t i = 0; i < l; ++i) {
        forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ std::get<0>(syncmers[i]);
        reverseRolledHash = seeding::rol(reverseRolledHash, k) ^ std::get<0>(syncmers[l-i-1]);
      }

      uint32_t iorder = 0;
      if (forwardRolledHash != reverseRolledHash) {
        size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
        curRead.uniqueSeedmers.emplace(minHash, std::vector<uint32_t>{iorder});
        curRead.seedmersList.emplace_back(mgsr::readSeedmer{
          minHash, std::get<3>(syncmers[0]), std::get<3>(syncmers[l-1])+k-1, reverseRolledHash < forwardRolledHash, iorder});
        ++iorder;
      }

      // rest of kminmer
      for (uint64_t i = 1; i < syncmers.size()-l+1; ++i) {
        if (!std::get<2>(syncmers[i-1]) || !std::get<2>(syncmers[i+l-1])) {
          std::cout << "invalid syncmer" << std::endl;
          exit(0);
        }
        const size_t& prevSyncmerHash = std::get<0>(syncmers[i-1]);
        const size_t& nextSyncmerHash = std::get<0>(syncmers[i+l-1]);
        forwardRolledHash = seeding::rol(forwardRolledHash, k) ^ seeding::rol(prevSyncmerHash, k * l) ^ nextSyncmerHash;
        reverseRolledHash = seeding::ror(reverseRolledHash, k) ^ seeding::ror(prevSyncmerHash, k)     ^ seeding::rol(nextSyncmerHash, k * (l-1));

        if (forwardRolledHash != reverseRolledHash) {
          size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
          auto uniqueSeedmersIt = curRead.uniqueSeedmers.find(minHash);
          if (uniqueSeedmersIt == curRead.uniqueSeedmers.end()) {
            curRead.uniqueSeedmers.emplace(minHash, std::vector<uint32_t>{iorder});
            curRead.seedmersList.emplace_back(mgsr::readSeedmer{
              minHash, std::get<3>(syncmers[i]), std::get<3>(syncmers[i+l-1])+k-1, reverseRolledHash < forwardRolledHash, iorder});
            ++iorder;
          } else {
            uniqueSeedmersIt->second.push_back(iorder);
            curRead.seedmersList.emplace_back(mgsr::readSeedmer{
              uniqueSeedmersIt->first, std::get<3>(syncmers[i]), std::get<3>(syncmers[i+l-1])+k-1, reverseRolledHash < forwardRolledHash, iorder});
            ++iorder;
          }
        }
      }
    }
  });
  
  std::vector<size_t> sortedUniqueReadSeedmersIndices(uniqueReadSeedmers.size());
  for (size_t i = 0; i < uniqueReadSeedmers.size(); ++i) sortedUniqueReadSeedmersIndices[i] = i;
  tbb::parallel_sort(sortedUniqueReadSeedmersIndices.begin(), sortedUniqueReadSeedmersIndices.end(), [&uniqueReadSeedmers, fast_mode](size_t i1, size_t i2) {
    const auto& lhs = uniqueReadSeedmers[i1].seedmersList;
    const auto& rhs = uniqueReadSeedmers[i2].seedmersList;
    
    size_t minSize = std::min(lhs.size(), rhs.size());

    for (size_t i = 0; i < minSize; ++i) {
      if (lhs[i].hash != rhs[i].hash) {
        return lhs[i].hash < rhs[i].hash;
      }
    }

    if (lhs.size() != rhs.size()) {
      return lhs.size() < rhs.size();
    }
    
    return std::lexicographical_compare(
      lhs.begin(), lhs.end(),
      rhs.begin(), rhs.end(),
      [fast_mode](const mgsr::readSeedmer& a, const mgsr::readSeedmer& b) {
        if (!fast_mode) {
          if (a.begPos != b.begPos) return a.begPos < b.begPos;
          if (a.endPos != b.endPos) return a.endPos < b.endPos;
        }
        if (a.rev != b.rev) return a.rev < b.rev;
        return a.iorder < b.iorder;
      }
    );
  });

  reads.emplace_back(std::move(uniqueReadSeedmers[sortedUniqueReadSeedmersIndices[0]]));
  readSeedmersDuplicatesIndex.emplace_back(std::vector<size_t>());
  for (const auto& seqSortedIndex : seqToIndexVec[sortedUniqueReadSeedmersIndices[0]].second) {
    readSeedmersDuplicatesIndex.back().push_back(seqSortedIndex);
  }

  for (size_t i = 1; i < sortedUniqueReadSeedmersIndices.size(); ++i) {
    const auto& currSeedmers = uniqueReadSeedmers[sortedUniqueReadSeedmersIndices[i]];
    
    bool isDuplicate = false;
    
    if (currSeedmers.seedmersList.size() == reads.back().seedmersList.size()) {
      const auto& curr = currSeedmers.seedmersList;
      const auto& prev = reads.back().seedmersList;
      
      if (fast_mode) {
        isDuplicate = std::equal(curr.begin(), curr.end(), prev.begin(),
          [](const mgsr::readSeedmer& a, const mgsr::readSeedmer& b) {
            return a.hash == b.hash && a.rev == b.rev && a.iorder == b.iorder;
          });
      } else {
        isDuplicate = true;
        
        for (size_t j = 0; j < curr.size() && isDuplicate; ++j) {
          if (curr[j].hash != prev[j].hash || 
              curr[j].rev != prev[j].rev || 
              curr[j].iorder != prev[j].iorder) {
            isDuplicate = false;
            break;
          }
          
          size_t currLength = curr[j].endPos - curr[j].begPos;
          size_t prevLength = prev[j].endPos - prev[j].begPos;
          if (currLength != prevLength) {
            isDuplicate = false;
            break;
          }
          
          if (j > 0) {
            size_t currGap = curr[j].begPos - curr[j-1].endPos;
            size_t prevGap = prev[j].begPos - prev[j-1].endPos;
            if (currGap != prevGap) {
              isDuplicate = false;
              break;
            }
          }
        }
      }
    }
    
    if (!isDuplicate) {
      reads.emplace_back(std::move(uniqueReadSeedmers[sortedUniqueReadSeedmersIndices[i]]));
      readSeedmersDuplicatesIndex.emplace_back(std::vector<size_t>());
    }
    
    for (const auto& seqSortedIndex : seqToIndexVec[sortedUniqueReadSeedmersIndices[i]].second) {
      readSeedmersDuplicatesIndex.back().push_back(seqSortedIndex);
    }
  }

  if (lowMemory) {
    ReadDebruijnGraph readDebruijnGraph;
    readDebruijnGraph.buildGraph(reads);
    readDebruijnGraph.sortReads(reads, readSeedmersDuplicatesIndex);
  } else {
    std::random_device rd;
    std::mt19937 g(rd());
    
    // Create indices for shuffling
    std::vector<size_t> indices(reads.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    // Shuffle the indices
    std::shuffle(indices.begin(), indices.end(), g);
    
    // Apply the shuffle to both vectors
    std::vector<mgsr::Read> shuffledReads;
    std::vector<std::vector<size_t>> shuffledReadSeedmersDuplicatesIndex;
    shuffledReads.reserve(reads.size());
    shuffledReadSeedmersDuplicatesIndex.reserve(readSeedmersDuplicatesIndex.size());
    
    for (size_t idx : indices) {
      shuffledReads.emplace_back(std::move(reads[idx]));
      shuffledReadSeedmersDuplicatesIndex.emplace_back(std::move(readSeedmersDuplicatesIndex[idx]));
    }
    
    // Replace original vectors with shuffled ones
    reads = std::move(shuffledReads);
    readSeedmersDuplicatesIndex = std::move(shuffledReadSeedmersDuplicatesIndex);
  }





  numPassedReads = reads.size();
  if (skipSingleton) {
    absl::flat_hash_map<uint64_t, uint64_t> kminmerCounts;
    for (uint32_t i = 0; i < reads.size(); ++i) {
      reads[i].seedmerStates.resize(reads[i].seedmersList.size(), mgsr::SeedmerState{false, false, false});
      for (const auto& seedmer : reads[i].uniqueSeedmers) {
        kminmerCounts[seedmer.first] += seedmer.second.size();
      }
    }

    for (uint32_t i = 0; i < reads.size(); ++i) {
      for (const auto& seedmer : reads[i].seedmersList) {
        if (kminmerCounts.at(seedmer.hash) == 1) {
          reads[i].readType = mgsr::ReadType::CONTAINS_SINGLETON;
          ++numSingletonReads;
          --numPassedReads;
          break;
        } 
      }
      if (reads[i].readType != mgsr::ReadType::CONTAINS_SINGLETON) {
        size_t estimatedUpdates = std::numeric_limits<size_t>::max();
        for (const auto& seedmer : reads[i].seedmersList) {
          allSeedmerHashesSet.insert(seedmer.hash);
        }
      }
    }
  } else {
    for (uint32_t i = 0; i < reads.size(); ++i) {
      reads[i].seedmerStates.resize(reads[i].seedmersList.size(), mgsr::SeedmerState{false, false, false});
      size_t estimatedUpdates = std::numeric_limits<size_t>::max();
      for (const auto& seedmer : reads[i].seedmersList) {
        allSeedmerHashesSet.insert(seedmer.hash);
      }
    }
  }

  const size_t chunkSize = (reads.size() + numThreads - 1) / numThreads;
  for (size_t i = 0; i < numThreads; ++i) {
    size_t start = i * chunkSize;
    size_t end = (i == numThreads - 1) ? reads.size() : (i + 1) * chunkSize;
    if (start < reads.size()) {
      threadRanges[i].first = start;
      threadRanges[i].second = end;
    }
  }

  std::cerr << "Collapsed " << readSequences.size() << " raw reads to " << reads.size() << " sketched kminmer sets" << std::endl;
  if (skipSingleton) {
    std::cerr << "SkipSingleton turned on: " << numSingletonReads << " reads with singletons will be skipped during placement and EM... "
              << "Total reads to process: " << numPassedReads << std::endl;
  }
  for (size_t i = 0; i < numThreads; ++i) {
    std::cerr << "  Thread " << i << " will process " << threadRanges[i].second - threadRanges[i].first << " reads: " << threadRanges[i].first << " -> " << threadRanges[i].second << std::endl;
  }

}



void mgsr::mgsrPlacer::initializeQueryData(std::span<mgsr::Read> reads, bool fast_mode) {
  this->reads = reads;
  for (uint32_t i = 0; i < reads.size(); ++i) {
    for (const auto& seedmer : reads[i].uniqueSeedmers) {
      for (const auto& seedmerIndex : seedmer.second) {
        seedmerToReads[seedmer.first].emplace_back(i, seedmerIndex);
      }
    }
  }

  // initialize score index structures
  size_t numReads = reads.size();
  size_t numNodes = liteTree->allLiteNodes.size();
  // maxScoreNodeIndex.resize(numReads);
  if (fast_mode) {
    kminmerMatches.resize(numReads, std::make_pair(0, 0));
    perNodeKminmerMatchesDeltasIndex.resize(numNodes);
  } else {
    readScores.resize(numReads, 0);
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
  uint64_t start,
  uint64_t end,
  bool toGap,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates,
  bool recordGapMapUpdates
) {
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
    updateGapMapStep(gapMap, update.second.first, update.second.second, update.first, backtrack, gapMapUpdates, true);
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
    updateGapMapStep(gapMap, curBeg, curEnd, it->first, backtrack, gapMapUpdates, false);
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

std::vector<panmapUtils::NewSyncmerRange> mgsr::mgsrIndexBuilder::computeNewSyncmerRangesWalk(
  panmanUtils::Node* node,
  size_t dfsIndex,
  const panmapUtils::BlockSequences& blockSequences,
  const std::vector<char>& blockExistsDelayed,
  const std::vector<char>& blockStrandDelayed,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersChangeRecord
) {
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges;
  if (localMutationRanges.empty()) {
    return newSyncmerRanges;
  }

  const std::vector<char>& blockExists = blockSequences.blockExists;
  const std::vector<char>& blockStrand = blockSequences.blockStrand;
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;

  std::sort(localMutationRanges.begin(), localMutationRanges.end(), [&globalCoords, &blockStrand](const auto& a, const auto& b) {
    return globalCoords.getScalarFromCoord(a.first, blockStrand[a.first.primaryBlockId]) < globalCoords.getScalarFromCoord(b.first, blockStrand[b.first.primaryBlockId]);
  });

  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> mergedLocalMutationRanges{localMutationRanges.front()};
  for (size_t i = 1; i < localMutationRanges.size(); ++i) {
    const auto& [curBeg, curEnd] = mergedLocalMutationRanges.back();
    const auto& [nextBeg, nextEnd] = localMutationRanges[i];
    
    // check if the current range and the next range are adjacent on their global scalar coordinates
    if (globalCoords.getScalarFromCoord(curEnd, blockStrand[curEnd.primaryBlockId]) + 1 >= globalCoords.getScalarFromCoord(nextBeg, blockStrand[nextBeg.primaryBlockId])) {
      if (globalCoords.getScalarFromCoord(nextEnd, blockStrand[nextEnd.primaryBlockId]) > globalCoords.getScalarFromCoord(curEnd, blockStrand[curEnd.primaryBlockId])) {
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
    auto curBegScalarTest = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
    auto curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
    auto leftGapMapIt = gapMap.lower_bound(globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]));
    auto rightGapMapIt = gapMap.lower_bound(globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]));

    // expand to the left... if reach newSyncmerRanges.back(), merge
    bool reachedEnd = false;
    uint32_t offset = 0;
    leftGapMapIt = leftGapMapIt == gapMap.begin() ? gapMap.begin() : std::prev(leftGapMapIt);
    while (offset < k - 1) {
      auto curBegScalar = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
      if (curBegScalar == 0) {
        break;
      }
      if (leftGapMapIt == gapMap.begin()) {
        if (curBegScalar - 1 > leftGapMapIt->second) {
          globalCoords.stepBackwardScalar(curBegCoord, blockStrand);
          --curBegScalarTest;
        } else if (curBegScalar >= leftGapMapIt->first) {
          if (leftGapMapIt->first == 0) {
            break;
          } else {
            curBegScalar = leftGapMapIt->first - 1;
            curBegScalarTest = leftGapMapIt->first - 1;
            curBegCoord = globalCoords.getCoordFromScalar(curBegScalar);
            if (!blockStrand[curBegCoord.primaryBlockId]) {
              curBegCoord = globalCoords.getCoordFromScalar(curBegScalar, false);
            }
          }
        }
      } else if (curBegScalar - 1 > leftGapMapIt->second) {
        globalCoords.stepBackwardScalar(curBegCoord, blockStrand);
        --curBegScalarTest;
      } else {
        curBegScalar = leftGapMapIt->first - 1;
        curBegScalarTest = leftGapMapIt->first - 1;
        curBegCoord = globalCoords.getCoordFromScalar(curBegScalar);
        if (!blockStrand[curBegCoord.primaryBlockId]) {
          curBegCoord = globalCoords.getCoordFromScalar(curBegScalar, false);
        }
        leftGapMapIt = leftGapMapIt == gapMap.begin() ? gapMap.begin() : std::prev(leftGapMapIt);
      }
      

      if (!newSyncmerRanges.empty() 
          && globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]) <= globalCoords.getScalarFromCoord(curSyncmerRange.endCoord, blockStrand[curSyncmerRange.endCoord.primaryBlockId])
      ) {
        // reached current newSyncmerRange... merge
        curBegCoord = curSyncmerRange.begCoord;
        curBegScalarTest = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
        syncmerRangeBegCoord = curBegCoord;
        newSyncmerRanges.pop_back();
        break;
      }

      if (!blockExists[curBegCoord.primaryBlockId]) {
        curBegCoord = globalCoords.blockEdgeCoords[curBegCoord.primaryBlockId].start;
        curBegScalarTest = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
        continue;
      }

      if (blockExists[curBegCoord.primaryBlockId] && blockSequences.getSequenceBase(curBegCoord) != '-') {
        offset++;
        syncmerRangeBegCoord = curBegCoord;
      }
    }


    // expand to the right... if reach mergedLocalMutationRanges[localMutationRangeIndex + 1], merge
    offset = 0;
    if (rightGapMapIt == gapMap.begin()) {
      rightGapMapIt = gapMap.begin();
    } else if (rightGapMapIt == gapMap.end()) {
      rightGapMapIt = std::prev(rightGapMapIt);
    } else if (rightGapMapIt->first != curEndScalar && curEndScalar <= std::prev(rightGapMapIt)->second) {
      rightGapMapIt = std::prev(rightGapMapIt);
    }

    while (offset < k - 1) {
      if (curEndScalar == globalCoords.lastScalarCoord) {
        reachedEnd = true;
        break;
      }

      if (rightGapMapIt == gapMap.end()) {
        auto lastGapMapIt = std::prev(gapMap.end());
        if (curEndScalar <= lastGapMapIt->second && lastGapMapIt->second != globalCoords.lastScalarCoord) {
          curEndScalar = lastGapMapIt->second + 1;
          curEndCoord = globalCoords.getCoordFromScalar(curEndScalar);
          if (!blockStrand[curEndCoord.primaryBlockId]) {
            curEndCoord = globalCoords.getCoordFromScalar(curEndScalar, false);
          }
        }
      } else if ((curEndScalar >= rightGapMapIt->first && curEndScalar <= rightGapMapIt->second) || curEndScalar + 1 >= rightGapMapIt->first) {
        if (rightGapMapIt->second == globalCoords.lastScalarCoord) {
          if (localMutationRangeIndex == mergedLocalMutationRanges.size() - 1) {
            reachedEnd = true;
            break;
          } else {
            curEndCoord = mergedLocalMutationRanges[localMutationRangeIndex + 1].second;
            curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
            syncmerRangeEndCoord = curEndCoord;
            localMutationRangeIndex++;
            offset = 0;
            continue;
          }
        }
        curEndScalar = rightGapMapIt->second + 1;
        curEndCoord = globalCoords.getCoordFromScalar(curEndScalar);
        if (!blockStrand[curEndCoord.primaryBlockId]) {
          curEndCoord = globalCoords.getCoordFromScalar(curEndScalar, false);
        }
        rightGapMapIt = std::next(rightGapMapIt);
      } else {
        globalCoords.stepForwardScalar(curEndCoord, blockStrand);
        ++curEndScalar;
      }
      
      if (!blockExists[curEndCoord.primaryBlockId]) {
        curEndCoord = globalCoords.blockEdgeCoords[curEndCoord.primaryBlockId].end;
        curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
        continue;
      }
      if (localMutationRangeIndex != mergedLocalMutationRanges.size() - 1
          && globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]) >= globalCoords.getScalarFromCoord(mergedLocalMutationRanges[localMutationRangeIndex + 1].first, blockStrand[mergedLocalMutationRanges[localMutationRangeIndex + 1].first.primaryBlockId])
      ) {
        // reached next mutation range... merge
        curEndCoord = mergedLocalMutationRanges[localMutationRangeIndex + 1].second;
        curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
        rightGapMapIt = gapMap.lower_bound(curEndScalar);
        if (rightGapMapIt == gapMap.begin()) {
          rightGapMapIt = gapMap.begin();
        } else if (rightGapMapIt == gapMap.end()) {
          rightGapMapIt = std::prev(rightGapMapIt);
        } else if (rightGapMapIt->first != curEndScalar && curEndScalar <= std::prev(rightGapMapIt)->second) {
          rightGapMapIt = std::prev(rightGapMapIt);
        }
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
    curSyncmerRange = newSyncmerRanges.back();
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
        curCoord = globalCoords.stepForwardScalar(globalCoords.blockEdgeCoords[curCoord.primaryBlockId].end, blockStrand);
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
      globalCoords.stepForwardScalar(curCoord, blockStrand);
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

std::vector<panmapUtils::NewSyncmerRange> mgsr::mgsrIndexBuilder::computeNewSyncmerRangesJump(
  panmanUtils::Node* node,
  size_t dfsIndex,
  const panmapUtils::BlockSequences& blockSequences,
  const std::vector<char>& blockExistsDelayed,
  const std::vector<char>& blockStrandDelayed,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersChangeRecord
) {
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges;
  if (localMutationRanges.empty()) {
    return newSyncmerRanges;
  }

  const std::vector<char>& blockExists = blockSequences.blockExists;
  const std::vector<char>& blockStrand = blockSequences.blockStrand;
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;

  std::sort(localMutationRanges.begin(), localMutationRanges.end(), [&globalCoords, &blockStrand](const auto& a, const auto& b) {
    return globalCoords.getScalarFromCoord(a.first, blockStrand[a.first.primaryBlockId]) < globalCoords.getScalarFromCoord(b.first, blockStrand[b.first.primaryBlockId]);
  });

  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> mergedLocalMutationRanges{localMutationRanges.front()};
  for (size_t i = 1; i < localMutationRanges.size(); ++i) {
    const auto& [curBeg, curEnd] = mergedLocalMutationRanges.back();
    const auto& [nextBeg, nextEnd] = localMutationRanges[i];

    
    // check if the current range and the next range are adjacent on their global scalar coordinates
    if (globalCoords.getScalarFromCoord(curEnd, blockStrand[curEnd.primaryBlockId]) + 1 >= globalCoords.getScalarFromCoord(nextBeg, blockStrand[nextBeg.primaryBlockId])) {
      if (globalCoords.getScalarFromCoord(nextEnd, blockStrand[nextEnd.primaryBlockId]) > globalCoords.getScalarFromCoord(curEnd, blockStrand[curEnd.primaryBlockId])) {
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
    auto curBegScalarTest = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
    auto curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
    auto leftGapMapIt = gapMap.lower_bound(globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]));
    auto rightGapMapIt = gapMap.lower_bound(globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]));

    // expand to the left... if reach newSyncmerRanges.back(), merge
    bool reachedEnd = false;
    uint32_t offset = 0;
    leftGapMapIt = leftGapMapIt == gapMap.begin() ? gapMap.begin() : std::prev(leftGapMapIt);
    while (offset < k - 1) {
      auto curBegScalar = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
      if (curBegScalar == 0) {
        break;
      }
      if (leftGapMapIt == gapMap.begin()) {
        if (curBegScalar - 1 > leftGapMapIt->second) {
          globalCoords.stepBackwardScalar(curBegCoord, blockStrand);
          --curBegScalarTest;
        } else if (curBegScalar >= leftGapMapIt->first) {
          if (leftGapMapIt->first == 0) {
            break;
          } else {
            curBegScalar = leftGapMapIt->first - 1;
            curBegScalarTest = leftGapMapIt->first - 1;
            curBegCoord = globalCoords.getCoordFromScalar(curBegScalar);
            if (!blockStrand[curBegCoord.primaryBlockId]) {
              curBegCoord = globalCoords.getCoordFromScalar(curBegScalar, false);
            }
          }
        }
      } else if (curBegScalar - 1 > leftGapMapIt->second) {
        globalCoords.stepBackwardScalar(curBegCoord, blockStrand);
        --curBegScalarTest;
      } else {
        curBegScalar = leftGapMapIt->first - 1;
        curBegScalarTest = leftGapMapIt->first - 1;
        curBegCoord = globalCoords.getCoordFromScalar(curBegScalar);
        if (!blockStrand[curBegCoord.primaryBlockId]) {
          curBegCoord = globalCoords.getCoordFromScalar(curBegScalar, false);
        }
        leftGapMapIt = leftGapMapIt == gapMap.begin() ? gapMap.begin() : std::prev(leftGapMapIt);
      }


      if (!newSyncmerRanges.empty() 
          && globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]) <= globalCoords.getScalarFromCoord(curSyncmerRange.endCoord, blockStrand[curSyncmerRange.endCoord.primaryBlockId])
      ) {
        // reached current newSyncmerRange... merge
        curBegCoord = curSyncmerRange.begCoord;
        curBegScalarTest = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
        syncmerRangeBegCoord = curBegCoord;
        newSyncmerRanges.pop_back();
        break;
      }

      if (!blockExists[curBegCoord.primaryBlockId]) {
        curBegCoord = globalCoords.blockEdgeCoords[curBegCoord.primaryBlockId].start;
        curBegScalarTest = globalCoords.getScalarFromCoord(curBegCoord, blockStrand[curBegCoord.primaryBlockId]);
        continue;
      }

      if (blockExists[curBegCoord.primaryBlockId] && blockSequences.getSequenceBase(curBegCoord) != '-') {
        offset++;
        syncmerRangeBegCoord = curBegCoord;
      }
    }


    // expand to the right... if reach mergedLocalMutationRanges[localMutationRangeIndex + 1], merge
    offset = 0;
    if (rightGapMapIt == gapMap.begin()) {
      rightGapMapIt = gapMap.begin();
    } else if (rightGapMapIt == gapMap.end()) {
      rightGapMapIt = std::prev(rightGapMapIt);
    } else if (rightGapMapIt->first != curEndScalar && curEndScalar <= std::prev(rightGapMapIt)->second) {
      rightGapMapIt = std::prev(rightGapMapIt);
    }
    while (offset < k - 1) {
      if (curEndScalar == globalCoords.lastScalarCoord) {
        reachedEnd = true;
        break;
      }

      if (rightGapMapIt == gapMap.end()) {
        auto lastGapMapIt = std::prev(gapMap.end());
        if (curEndScalar <= lastGapMapIt->second) {
          if (lastGapMapIt->second != globalCoords.lastScalarCoord) {
            curEndScalar = lastGapMapIt->second + 1;
            curEndCoord = globalCoords.getCoordFromScalar(curEndScalar);
            if (!blockStrand[curEndCoord.primaryBlockId]) {
              curEndCoord = globalCoords.getCoordFromScalar(curEndScalar, false);
            }
          }
        } else {
          globalCoords.stepForwardScalar(curEndCoord, blockStrand);
          ++curEndScalar;
        }
      } else if (curEndScalar <= rightGapMapIt->second && (curEndScalar >= rightGapMapIt->first || curEndScalar + 1 >= rightGapMapIt->first)) {
        if (rightGapMapIt->second == globalCoords.lastScalarCoord) {
          if (localMutationRangeIndex == mergedLocalMutationRanges.size() - 1) {
            reachedEnd = true;
            break;
          } else {
            curEndCoord = mergedLocalMutationRanges[localMutationRangeIndex + 1].second;
            curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
            syncmerRangeEndCoord = curEndCoord;
            localMutationRangeIndex++;
            offset = 0;
            continue;
          }
        }
        curEndScalar = rightGapMapIt->second + 1;
        curEndCoord = globalCoords.getCoordFromScalar(curEndScalar);
        if (!blockStrand[curEndCoord.primaryBlockId]) {
          curEndCoord = globalCoords.getCoordFromScalar(curEndScalar, false);
        }
        rightGapMapIt = std::next(rightGapMapIt);
      } else {
        globalCoords.stepForwardScalar(curEndCoord, blockStrand);
        ++curEndScalar;
      }

      if (!blockExists[curEndCoord.primaryBlockId]) {
        curEndCoord = globalCoords.blockEdgeCoords[curEndCoord.primaryBlockId].end;
        curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
        continue;
      }
      if (localMutationRangeIndex != mergedLocalMutationRanges.size() - 1
          && globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]) >= globalCoords.getScalarFromCoord(mergedLocalMutationRanges[localMutationRangeIndex + 1].first, blockStrand[mergedLocalMutationRanges[localMutationRangeIndex + 1].first.primaryBlockId])
      ) {
        // reached next mutation range... merge
        curEndCoord = mergedLocalMutationRanges[localMutationRangeIndex + 1].second;
        curEndScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
        rightGapMapIt = gapMap.lower_bound(curEndScalar);
        if (rightGapMapIt == gapMap.begin()) {
          rightGapMapIt = gapMap.begin();
        } else if (rightGapMapIt == gapMap.end()) {
          rightGapMapIt = std::prev(rightGapMapIt);
        } else if (rightGapMapIt->first != curEndScalar && curEndScalar <= std::prev(rightGapMapIt)->second) {
          rightGapMapIt = std::prev(rightGapMapIt);
        }
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
    curSyncmerRange = newSyncmerRanges.back();
    if (reachedEnd) {
      offsetsToDelete = k - offset - 1;
      endOffset = offset;
      break;
    }
    localMutationRangeIndex++;
  }

  
  const auto lastScalarCoord = globalCoords.lastScalarCoord;
  for (size_t i = 0; i < newSyncmerRanges.size(); i++) {
    panmapUtils::NewSyncmerRange& syncmerRange = newSyncmerRanges[i];
    panmapUtils::Coordinate curBegCoord = syncmerRange.begCoord;
    panmapUtils::Coordinate curCoord = curBegCoord;
    auto curCoordScalar = globalCoords.getScalarFromCoord(curCoord, blockStrand[curCoord.primaryBlockId]);
    const panmapUtils::Coordinate curEndCoord = syncmerRange.endCoord;
    const auto curEndCoordScalar = globalCoords.getScalarFromCoord(curEndCoord, blockStrand[curEndCoord.primaryBlockId]);
    auto curCoordGapMapIt = gapMap.lower_bound(curCoordScalar);

    // if startChar == '-', it means the start position is in the first gap run group and if the first gap run extends to the end of the genome, then we can skip this syncmer range
    const char startChar = blockSequences.getSequenceBase(curCoord);
    if (startChar == '-') {
      const auto curBlockId = curCoord.primaryBlockId;
      if (gapMap.begin()->second == lastScalarCoord) {
        continue;
      } else if (blockStrandDelayed[curBlockId] == blockStrand[curBlockId] || !blockExistsDelayed[curBlockId] || !blockExists[curBlockId]) {
        curCoord = globalCoords.getCoordFromScalar(gapMap.begin()->second + 1);
        if (!blockStrand[curCoord.primaryBlockId]) {
          curCoord = globalCoords.getCoordFromScalar(gapMap.begin()->second + 1, false);
        }
        curCoordScalar = globalCoords.getScalarFromCoord(curCoord, blockStrand[curCoord.primaryBlockId]);
        curBegCoord = curCoord;
      }
    }
  

    std::string& localRangeSeq = syncmerRange.localRangeSeq;
    std::vector<uint64_t>& localRangeCoordToGlobalScalarCoords = syncmerRange.localRangeCoordToGlobalScalarCoords;
    std::vector<uint64_t>& seedsToDelete = syncmerRange.seedsToDelete;
    std::vector<uint64_t>& localRangeCoordToBlockId = syncmerRange.localRangeCoordToBlockId;
    localRangeSeq = "";
    std::vector<uint64_t>().swap(localRangeCoordToGlobalScalarCoords);
    std::vector<uint64_t>().swap(seedsToDelete);
    std::vector<uint64_t>().swap(localRangeCoordToBlockId);
    bool recomputeBlock = false;
    bool recomputeInProgress = false;
    while (true) {
      const auto curBlockId = curCoord.primaryBlockId;
      if (curCoordGapMapIt == gapMap.end()) {
        // do nothing... Current and all subsequent positions are non-gap
      } else if (recomputeBlock || (blockStrandDelayed[curBlockId] != blockStrand[curBlockId] && blockExistsDelayed[curBlockId] && blockExists[curBlockId])) {
        // need to recompute this whole block
        if (!recomputeInProgress) {
          if (curCoord != curBegCoord) {
            if (blockStrand[curBlockId]) {
              curCoord = globalCoords.blockEdgeCoords[curBlockId].start;
            } else {
              curCoord = globalCoords.blockEdgeCoords[curBlockId].end;
            }
            curCoordScalar = globalCoords.getScalarFromCoord(curCoord, blockStrand[curCoord.primaryBlockId]);
          }
          recomputeInProgress = true;
        }
      }
      
      if (!blockExists[curCoord.primaryBlockId]) {
        if (curCoord.primaryBlockId == curEndCoord.primaryBlockId) {
          break;
        }
        curCoord = globalCoords.stepForwardScalar(globalCoords.blockEdgeCoords[curCoord.primaryBlockId].end, blockStrand);
        curCoordScalar = globalCoords.getScalarFromCoord(curCoord, blockStrand[curCoord.primaryBlockId]);
        continue;
      }

      char curNuc = blockSequences.getSequenceBase(curCoord);

      if (curNuc != '-') {
        if (!blockStrand[curCoord.primaryBlockId]) curNuc = panmanUtils::getComplementCharacter(curNuc);
        localRangeSeq += curNuc;
        localRangeCoordToGlobalScalarCoords.push_back(curCoordScalar);
        localRangeCoordToBlockId.push_back(curCoord.primaryBlockId);
      } else if (refOnSyncmers[curCoordScalar].has_value()) {
        blockOnSyncmers[curCoord.primaryBlockId].erase(curCoordScalar);
        if (blockOnSyncmers[curCoord.primaryBlockId].empty()) blockOnSyncmers.erase(curCoord.primaryBlockId);
        blockOnSyncmersChangeRecord.emplace_back(curCoord.primaryBlockId, curCoordScalar, panmapUtils::seedChangeType::DEL);
        seedsToDelete.push_back(curCoordScalar);
      }
      if (curCoordScalar == curEndCoordScalar)  {
        break;
      }
      if (recomputeInProgress) {
        globalCoords.stepForwardScalar(curCoord, blockStrand);
        ++curCoordScalar;
        if (curBlockId != curCoord.primaryBlockId) {
          recomputeBlock = false;
          recomputeInProgress = false;
          if (curCoordGapMapIt != gapMap.end()) {
            while (curCoordGapMapIt != gapMap.end() && !((std::prev(curCoordGapMapIt)->second < curCoordScalar) && (curCoordScalar < curCoordGapMapIt->first))) {
              ++curCoordGapMapIt;
            }
          }
        }
      } else if (curCoordGapMapIt != gapMap.end() && curCoordScalar == curCoordGapMapIt->first - 1) {
        // step over gap run 
        if (curCoordGapMapIt->second == lastScalarCoord) {
          break;
        } else {
          curCoord = globalCoords.getCoordFromScalar(curCoordGapMapIt->second + 1);
          if (!blockStrand[curCoord.primaryBlockId]) {
            curCoord = globalCoords.getCoordFromScalar(curCoordGapMapIt->second + 1, false);
          }
          curCoordScalar = globalCoords.getScalarFromCoord(curCoord, blockStrand[curCoord.primaryBlockId]);
          if (recomputeBlock && curBlockId != curBegCoord.primaryBlockId) {
            recomputeBlock = false;
            recomputeInProgress = false;
          }
          ++curCoordGapMapIt;
        }
      } else {
        globalCoords.stepForwardScalar(curCoord, blockStrand);
        ++curCoordScalar;
        if (recomputeBlock && curBlockId != curBegCoord.primaryBlockId) {
          recomputeBlock = false;
          recomputeInProgress = false;
        }
      }
    }

    if (i == newSyncmerRanges.size() - 1 && offsetsToDelete != -1) {
      for (size_t j = localRangeSeq.size() - endOffset - offsetsToDelete; j < localRangeSeq.size(); j++) {
        if (refOnSyncmers[localRangeCoordToGlobalScalarCoords[j]].has_value()) {
          auto curGlobalScalarCoord = localRangeCoordToGlobalScalarCoords[j];
          auto curBlockId = localRangeCoordToBlockId[j];
          blockOnSyncmers[curBlockId].erase(curGlobalScalarCoord);
          if (blockOnSyncmers[curBlockId].empty()) blockOnSyncmers.erase(curBlockId);
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, curGlobalScalarCoord, panmapUtils::seedChangeType::DEL);
          seedsToDelete.push_back(localRangeCoordToGlobalScalarCoords[j]);
        }
      }
    }
  }

  return newSyncmerRanges;
}

std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> mgsr::mgsrIndexBuilder::computeNewKminmerRanges(
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
  const uint64_t dfsIndex
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
          }
          //  else {
          //   curEndIt = refOnSyncmersMap.end();
          //   break;
          // }
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
  std::vector<char> &blockExistsDelayed,
  std::vector<char> &blockStrandDelayed,
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

  // Nuc deletions on block without block mutation
  std::vector<uint32_t> potentialSyncmerDeletions;

  blockMutationRecord.reserve(node->blockMutation.size());
  nucMutationRecord.reserve(node->nucMutation.size() * 6);
  potentialSyncmerDeletions.reserve(node->nucMutation.size() * 6);
  localMutationRanges.reserve(node->blockMutation.size() + node->nucMutation.size() * 6);
  gapRunUpdates.reserve(node->nucMutation.size() * 6 + node->blockMutation.size() * 10);
  applyMutations(node, dfsIndex, blockSequences, invertedBlocks, globalCoords, localMutationRanges, blockMutationRecord, nucMutationRecord, gapRunUpdates, invertedBlocksBacktracks, potentialSyncmerDeletions, blockExistsDelayed, blockStrandDelayed);
  blockMutationRecord.shrink_to_fit();
  nucMutationRecord.shrink_to_fit();
  potentialSyncmerDeletions.shrink_to_fit();
  localMutationRanges.shrink_to_fit();
  gapRunUpdates.shrink_to_fit();

  std::sort(gapRunUpdates.begin(), gapRunUpdates.end(), [&](const auto& a, const auto& b) { return a.second.first < b.second.first; });
  updateGapMap(node, dfsIndex, gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

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
  // if (dfsIndex >= 0) {
  //   makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, (uint64_t)globalCoords.lastScalarCoord);
  // }

  bool useJump = true;
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges;
  if (useJump) {
    newSyncmerRanges = computeNewSyncmerRangesJump(node, dfsIndex, blockSequences,  blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord);
  } else {
    newSyncmerRanges = computeNewSyncmerRangesWalk(node, dfsIndex, blockSequences,  blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord);
  }

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

  if (useJump) {
    for (uint32_t pos : potentialSyncmerDeletions) {
      if (refOnSyncmers[pos].has_value()) {
        refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers[pos].value());
        refOnSyncmers[pos] = std::nullopt;
        const auto blockId = globalCoords.getBlockIdFromScalar(pos);
        blockOnSyncmers[blockId].erase(pos);
        if (blockOnSyncmers[blockId].empty()) blockOnSyncmers.erase(blockId);
        blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);
        refOnSyncmersMap.erase(pos);
      }
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
  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges = computeNewKminmerRanges(refOnSyncmersChangeRecord, dfsIndex);

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
  std::vector<std::pair<uint64_t, uint64_t>> substitutedSeedIndices;
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
    auto& curRefOnKminmer = refOnKminmers[*indexingIt];
    if (forwardHash != reverseHash) {
      seeding::uniqueKminmer_t uniqueKminmer{*indexingIt, refOnSyncmers[*curIt].value().endPos, std::min(forwardHash, reverseHash), reverseHash < forwardHash};
      bool substitution = curRefOnKminmer.has_value();
      if (substitution) {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, curRefOnKminmer.value());
      } else {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
      }

      auto uniqueKminmerIndexIt = kminmerToUniqueIndex.find(uniqueKminmer);
      if (uniqueKminmerIndexIt == kminmerToUniqueIndex.end()) {
        uniqueKminmers.emplace_back(uniqueKminmer);
        kminmerToUniqueIndex[uniqueKminmer] = uniqueKminmers.size() - 1;
        if (substitution) {
          substitutedSeedIndices.emplace_back(curRefOnKminmer.value(), uniqueKminmers.size() - 1);
        } else {
          addedSeedIndices.push_back(uniqueKminmers.size() - 1);
        }
        curRefOnKminmer = uniqueKminmers.size() - 1;
      } else {
        if (substitution) {
          substitutedSeedIndices.emplace_back(curRefOnKminmer.value(), uniqueKminmerIndexIt->second);
        } else {
          addedSeedIndices.push_back(uniqueKminmerIndexIt->second);
        }
        curRefOnKminmer = uniqueKminmerIndexIt->second;
      }
      
    } else {
      if (curRefOnKminmer.has_value()) {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, curRefOnKminmer.value());
        deletedSeedIndices.push_back(curRefOnKminmer.value());
        curRefOnKminmer = std::nullopt;
      }
    }

    if (curIt == endIt) continue;
    
    while (curIt != endIt) {
      ++curIt;
      if (curIt == refOnSyncmersMap.end()) break;
      forwardHash = seeding::rol(forwardHash, k) ^ seeding::rol(refOnSyncmers[*indexingIt].value().hash, k * l) ^ refOnSyncmers[*curIt].value().hash;
      reverseHash = seeding::ror(reverseHash, k) ^ seeding::ror(refOnSyncmers[*indexingIt].value().hash, k)     ^ seeding::rol(refOnSyncmers[*curIt].value().hash, k * (l-1));
      ++indexingIt;

      auto& curRefOnKminmer = refOnKminmers[*indexingIt];
      if (forwardHash != reverseHash) {
        seeding::uniqueKminmer_t uniqueKminmer{*indexingIt, refOnSyncmers[*curIt].value().endPos, std::min(forwardHash, reverseHash), reverseHash < forwardHash};
        bool substitution = curRefOnKminmer.has_value();
        if (substitution) {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, curRefOnKminmer.value());
        } else {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
        }

        auto uniqueKminmerIndexIt = kminmerToUniqueIndex.find(uniqueKminmer);
        if (uniqueKminmerIndexIt == kminmerToUniqueIndex.end()) {
          uniqueKminmers.emplace_back(uniqueKminmer);
          kminmerToUniqueIndex[uniqueKminmer] = uniqueKminmers.size() - 1;
          if (substitution) {
            substitutedSeedIndices.emplace_back(curRefOnKminmer.value(), uniqueKminmers.size() - 1);
          } else {
            addedSeedIndices.push_back(uniqueKminmers.size() - 1);
          }
          curRefOnKminmer = uniqueKminmers.size() - 1;
        } else {
          if (substitution) {
            substitutedSeedIndices.emplace_back(curRefOnKminmer.value(), uniqueKminmerIndexIt->second);
          } else {
            addedSeedIndices.push_back(uniqueKminmerIndexIt->second);
          }
          curRefOnKminmer = uniqueKminmerIndexIt->second;
        }
      } else {
        if (curRefOnKminmer.has_value()) {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, curRefOnKminmer.value());
          deletedSeedIndices.push_back(curRefOnKminmer.value());
          curRefOnKminmer = std::nullopt;
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
        deletedSeedIndices.push_back(refOnKminmers[*delIt].value());
        refOnKminmers[*delIt] = std::nullopt;
      }
      if (delIt == refOnSyncmersMap.begin()) break;
    }
  }
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL && refOnKminmers[syncmerPos].has_value()) {
      refOnKminmersChangeRecord.emplace_back(syncmerPos, panmapUtils::seedChangeType::DEL, refOnKminmers[syncmerPos].value());
      deletedSeedIndices.push_back(refOnKminmers[syncmerPos].value());
      refOnKminmers[syncmerPos] = std::nullopt;
    }
  }

  std::sort(deletedSeedIndices.begin(), deletedSeedIndices.end(), [&](const auto& a, const auto& b) { return uniqueKminmers[a].startPos < uniqueKminmers[b].startPos; });
 
  //  Adding node changes to index
  NodeChanges::Builder curNodeChanges = perNodeChanges[dfsIndex];
  curNodeChanges.setNodeIndex(dfsIndex);

  // adding inserted/substituted seeds to index
  capnp::List<SeedDelta>::Builder seedDeltasBuilder = curNodeChanges.initSeedDeltas(addedSeedIndices.size() + deletedSeedIndices.size() + substitutedSeedIndices.size() * 2);
  size_t deltaSeedIndicesIndex = 0, deletedIdx = 0, addedIdx = 0, substitutedIdx = 0;
  
  while (deletedIdx < deletedSeedIndices.size() && addedIdx < addedSeedIndices.size() && substitutedIdx < substitutedSeedIndices.size()) {
    auto deletedStartPos = uniqueKminmers[deletedSeedIndices[deletedIdx]].startPos;
    auto addedStartPos = uniqueKminmers[addedSeedIndices[addedIdx]].startPos;
    auto substitutedStartPos = uniqueKminmers[substitutedSeedIndices[substitutedIdx].first].startPos;
    if (deletedStartPos <= addedStartPos && deletedStartPos <= substitutedStartPos) {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(deletedSeedIndices[deletedIdx]);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
      deltaSeedIndicesIndex++;
      deletedIdx++;
    } else if (addedStartPos <= substitutedStartPos) {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(addedSeedIndices[addedIdx]);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
      deltaSeedIndicesIndex++;
      addedIdx++;
    } else {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].first);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
      deltaSeedIndicesIndex++;
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].second);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
      deltaSeedIndicesIndex++;
      substitutedIdx++;
    }
  }
  
  while (deletedIdx < deletedSeedIndices.size() && addedIdx < addedSeedIndices.size()) {
    auto deletedStartPos = uniqueKminmers[deletedSeedIndices[deletedIdx]].startPos;
    auto addedStartPos = uniqueKminmers[addedSeedIndices[addedIdx]].startPos;
    if (deletedStartPos <= addedStartPos) {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(deletedSeedIndices[deletedIdx]);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
      deltaSeedIndicesIndex++;
      deletedIdx++;
    } else {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(addedSeedIndices[addedIdx]);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
      deltaSeedIndicesIndex++;
      addedIdx++;
    }
  }
  while (deletedIdx < deletedSeedIndices.size() && substitutedIdx < substitutedSeedIndices.size()) {
    auto deletedStartPos = uniqueKminmers[deletedSeedIndices[deletedIdx]].startPos;
    auto substitutedStartPos = uniqueKminmers[substitutedSeedIndices[substitutedIdx].first].startPos;
    if (deletedStartPos <= substitutedStartPos) {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(deletedSeedIndices[deletedIdx]);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
      deltaSeedIndicesIndex++;
      deletedIdx++;
    } else {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].first);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
      deltaSeedIndicesIndex++;
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].second);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
      deltaSeedIndicesIndex++;
      substitutedIdx++;
    }
  }
  while (addedIdx < addedSeedIndices.size() && substitutedIdx < substitutedSeedIndices.size()) {
    auto addedStartPos = uniqueKminmers[addedSeedIndices[addedIdx]].startPos;
    auto substitutedStartPos = uniqueKminmers[substitutedSeedIndices[substitutedIdx].first].startPos;
    if (addedStartPos <= substitutedStartPos) {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(addedSeedIndices[addedIdx]);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
      deltaSeedIndicesIndex++;
      addedIdx++;
    } else {
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].first);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
      deltaSeedIndicesIndex++;
      seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].second);
      seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
      deltaSeedIndicesIndex++;
      substitutedIdx++;
    }
  }
  
  while (deletedIdx < deletedSeedIndices.size()) {
    seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(deletedSeedIndices[deletedIdx]);
    seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
    deltaSeedIndicesIndex++;
    deletedIdx++;
  }
  while (addedIdx < addedSeedIndices.size()) {
    seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(addedSeedIndices[addedIdx]);
    seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
    deltaSeedIndicesIndex++;
    addedIdx++;
  }
  while (substitutedIdx < substitutedSeedIndices.size()) {
    seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].first);
    seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(true);
    deltaSeedIndicesIndex++;
    seedDeltasBuilder[deltaSeedIndicesIndex].setSeedIndex(substitutedSeedIndices[substitutedIdx].second);
    seedDeltasBuilder[deltaSeedIndicesIndex].setIsDeleted(false);
    deltaSeedIndicesIndex++;
    substitutedIdx++;
  }
  

  // adding coord deltas to index
  capnp::List<GapRunDelta>::Builder gapRunDeltaBuilder = curNodeChanges.initGapRunDeltas(gapRunUpdates.size());
  for (size_t i = 0; i < gapRunUpdates.size(); i++) {
    const auto& [toGap, range] = gapRunUpdates[i];
    gapRunDeltaBuilder[i].setStartPos(range.first);
    gapRunDeltaBuilder[i].setEndPos(range.second);
    gapRunDeltaBuilder[i].setToGap(toGap);
  }
  

  // adding inverted blocks to index
  capnp::List<uint32_t>::Builder invertedBlocksBuilder = curNodeChanges.initInvertedBlocks(invertedBlocksVec.size());
  for (size_t i = 0; i < invertedBlocksVec.size(); i++) {
    invertedBlocksBuilder.set(i, invertedBlocksVec[i]);
  }


  // // compare with brute force for debugging
  // if (dfsIndex >= 0) {

  //   compareBruteForceBuild(T, node, blockSequences, globalCoords, gapMap, degapCoordIndex, regapCoordIndex, refOnSyncmers, refOnSyncmersMap, blockOnSyncmers, refOnKminmers, uniqueKminmers, kminmerToUniqueIndex, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getT(), indexBuilder.getL(), indexBuilder.getOpen());
  // }



  revertGapMapInversions(gapRunBlockInversionBacktracks, gapMap);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapRunBlockInversionBacktracks); // gapRunBlockInversionBacktracks is no longer needed... clear memory

  // update delayed block states
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExistsDelayed[blockId] = newExists;
    blockStrandDelayed[blockId] = newStrand;
  }

  nodeToDfsIndex[node->identifier] = dfsIndex;
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
  panmapUtils::GlobalCoords globalCoords(blockSequences);
  refOnSyncmers.resize(globalCoords.lastScalarCoord + 1);
  refOnKminmers.resize(globalCoords.lastScalarCoord + 1);

  std::vector<char> blockExistsDelayed = blockSequences.blockExists;
  std::vector<char> blockStrandDelayed = blockSequences.blockStrand;

  std::map<uint64_t, uint64_t> gapMap{{0, globalCoords.lastScalarCoord}};
  std::unordered_set<uint64_t> invertedBlocks;

  // add lite tree to index
  LiteTree::Builder liteTreeBuilder = indexBuilder.initLiteTree();
  uint64_t dfsIndex = 0;
  buildIndexHelper(T->root, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks, dfsIndex);


  
  // Add unique k-min-mers to index
  capnp::List<SeedInfo>::Builder seedInfoBuilder = indexBuilder.initSeedInfo(uniqueKminmers.size());
  for (size_t i = 0; i < uniqueKminmers.size(); i++) {
    seedInfoBuilder[i].setHash(uniqueKminmers[i].hash);
    seedInfoBuilder[i].setStartPos(uniqueKminmers[i].startPos);
    seedInfoBuilder[i].setEndPos(uniqueKminmers[i].endPos);
    seedInfoBuilder[i].setIsReverse(uniqueKminmers[i].isReverse);
  }

  // Add block infos to index
  capnp::List<BlockRange>::Builder blockRangesBuilder = liteTreeBuilder.initBlockRanges(globalCoords.globalCoords.size());
  for (size_t i = 0; i < globalCoords.globalCoords.size(); i++) {
    blockRangesBuilder[i].setRangeBeg(globalCoords.getBlockStartScalar(i));
    blockRangesBuilder[i].setRangeEnd(globalCoords.getBlockEndScalar(i));
  }

  // Add node to index
  capnp::List<LiteNode>::Builder nodesBuilder = liteTreeBuilder.initLiteNodes(T->allNodes.size());
  for (const auto& [nodeId, node] : T->allNodes) {
    auto dfsIndex = nodeToDfsIndex[nodeId];
    nodesBuilder[dfsIndex].setId(nodeId);
    if (node->parent == nullptr) {
      nodesBuilder[dfsIndex].setParentIndex(0);
    } else {
      nodesBuilder[dfsIndex].setParentIndex(nodeToDfsIndex[node->parent->identifier]);
    }
  }




  std::cout << "Finished building index!" << std::endl;
}

void mgsr::mgsrIndexBuilder::writeIndex(const std::string& path) {
  int fd = ::open(path.c_str(), O_RDWR | O_CREAT, 0644);
  if (fd == -1) {
    std::cerr << "Error: failed to open file " << path << std::endl;
    std::exit(1);
  }
  std::cout << "Writing index to " << path << std::endl;
  capnp::writePackedMessageToFd(fd, outMessage);
  close(fd);
  std::cout << "Index written to " << path << std::endl;
}

int mgsr::open_file(const std::string& path) {
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd == -1) {
    std::cerr << "Error: failed to open file " << path << std::endl;
    std::exit(1);
  }
  return fd;
}

void mgsr::mgsrPlacer::addSeedAtPosition(uint64_t newKminmerIndex, std::unordered_set<uint64_t>& affectedSeedmers) {
  const auto& newKminmer = liteTree->seedInfos[newKminmerIndex];
  const uint32_t pos = newKminmer.startPos;
  auto posMapIt = positionMap.emplace(pos, newKminmerIndex).first;
  size_t hash = newKminmer.hash;
  affectedSeedmers.insert(hash);
  hashToPositionMap[hash].push_back(posMapIt);
}

void mgsr::mgsrPlacer::addSeedAtPosition(uint64_t newKminmerIndex) {
  const auto& newKminmer = liteTree->seedInfos[newKminmerIndex];
  const uint32_t pos = newKminmer.startPos;
  auto posMapIt = positionMap.emplace(pos, newKminmerIndex).first;
  size_t hash = newKminmer.hash;
  auto& positions = hashToPositionMap[hash];
  positions.push_back(posMapIt);
  delayedRefSeedmerStatus[hash] = positions.size() == 1 ? mgsr::RefSeedmerExistStatus::EXIST_UNIQUE : mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE;
}

void mgsr::mgsrPlacer::subSeedAtPosition(uint64_t newKminmerIndex, std::unordered_set<uint64_t>& affectedSeedmers) {
  const auto& newKminmer = liteTree->seedInfos[newKminmerIndex];
  const uint32_t pos = newKminmer.startPos;
  auto posMapIt = positionMap.find(pos);
  uint64_t& oldKminmerIndex = posMapIt->second;
  size_t oldHash = liteTree->seedInfos[oldKminmerIndex].hash;
  size_t newHash = newKminmer.hash;
  affectedSeedmers.insert(oldHash);
  affectedSeedmers.insert(newHash);
  auto oldHashToPositionIt = hashToPositionMap.find(oldHash);
  if (oldHashToPositionIt->second.size() == 1) {
    hashToPositionMap.erase(oldHashToPositionIt);
  } else {
    auto eraseIt = std::find(oldHashToPositionIt->second.begin(), oldHashToPositionIt->second.end(), posMapIt);
    oldHashToPositionIt->second.erase(eraseIt);
    if (oldHashToPositionIt->second.empty()) hashToPositionMap.erase(oldHashToPositionIt);
  }
  oldKminmerIndex = newKminmerIndex;
  hashToPositionMap[newHash].push_back(posMapIt);
}

void mgsr::mgsrPlacer::subSeedAtPosition(uint64_t newKminmerIndex) {
  const auto& newKminmer = liteTree->seedInfos[newKminmerIndex];
  const uint32_t pos = newKminmer.startPos;
  auto posMapIt = positionMap.find(pos);
  uint64_t& oldKminmerIndex = posMapIt->second;
  size_t oldHash = liteTree->seedInfos[oldKminmerIndex].hash;
  size_t newHash = newKminmer.hash;
  auto oldHashToPositionIt = hashToPositionMap.find(oldHash);
  if (oldHashToPositionIt->second.size() == 1) {
    hashToPositionMap.erase(oldHashToPositionIt);
    delayedRefSeedmerStatus.erase(oldHash);
  } else {
    auto eraseIt = std::find(oldHashToPositionIt->second.begin(), oldHashToPositionIt->second.end(), posMapIt);
    oldHashToPositionIt->second.erase(eraseIt);
    delayedRefSeedmerStatus[oldHash] = oldHashToPositionIt->second.size() == 1 ? mgsr::RefSeedmerExistStatus::EXIST_UNIQUE : mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE;
  }
  oldKminmerIndex = newKminmerIndex;
  auto& positions = hashToPositionMap[newHash];
  positions.push_back(posMapIt);
  delayedRefSeedmerStatus[newHash] = positions.size() == 1 ? mgsr::RefSeedmerExistStatus::EXIST_UNIQUE : mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE;
}

void mgsr::mgsrPlacer::delSeedAtPosition(uint64_t kminmerIndex, std::unordered_set<uint64_t>& affectedSeedmers) {
  auto posMapIt = positionMap.find(liteTree->seedInfos[kminmerIndex].startPos);
  const uint64_t oldKminmerIndex = posMapIt->second;
  size_t hash = liteTree->seedInfos[oldKminmerIndex].hash;
  affectedSeedmers.insert(hash);
  auto hashToPositionIt = hashToPositionMap.find(hash);
  if (hashToPositionIt->second.size() == 1) {
    hashToPositionMap.erase(hashToPositionIt);
  } else {
    auto eraseIt = std::find(hashToPositionIt->second.begin(), hashToPositionIt->second.end(), posMapIt);
    hashToPositionIt->second.erase(eraseIt);
    if (hashToPositionIt->second.empty()) hashToPositionMap.erase(hashToPositionIt);
  }
  positionMap.erase(posMapIt);
}

void mgsr::mgsrPlacer::delSeedAtPosition(uint64_t kminmerIndex) {
  auto posMapIt = positionMap.find(liteTree->seedInfos[kminmerIndex].startPos);
  const uint64_t oldKminmerIndex = posMapIt->second;
  size_t hash = liteTree->seedInfos[oldKminmerIndex].hash;
  auto hashToPositionIt = hashToPositionMap.find(hash);
  if (hashToPositionIt->second.size() == 1) {
    hashToPositionMap.erase(hashToPositionIt);
    delayedRefSeedmerStatus.erase(hash);
  } else {
    auto eraseIt = std::find(hashToPositionIt->second.begin(), hashToPositionIt->second.end(), posMapIt);
    hashToPositionIt->second.erase(eraseIt);
    delayedRefSeedmerStatus[hash] = hashToPositionIt->second.size() == 1 ? mgsr::RefSeedmerExistStatus::EXIST_UNIQUE : mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE;
  }
  positionMap.erase(posMapIt);
}

mgsr::RefSeedmerExistStatus mgsr::mgsrPlacer::getCurrentRefSeedmerExistStatus(uint64_t hash) {
  auto hashToPositionIt = hashToPositionMap.find(hash);
  if (hashToPositionIt == hashToPositionMap.end()) {
    return mgsr::RefSeedmerExistStatus::NOT_EXIST;
  } else if (hashToPositionIt->second.size() == 1) {
    return mgsr::RefSeedmerExistStatus::EXIST_UNIQUE;
  } else if (hashToPositionIt->second.size() > 1) {
    return mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE;
  }

  std::cerr << "Error: getRefSeedmerExistStatus: hash " << hash << " exists in " << hashToPositionIt->second.size() << " positions" << std::endl;
  std::exit(1);
  return mgsr::RefSeedmerExistStatus::NOT_EXIST;
}

mgsr::RefSeedmerExistStatus mgsr::mgsrPlacer::getDelayedRefSeedmerExistStatus(uint64_t hash) {
  auto hashToPositionIt = delayedRefSeedmerStatus.find(hash);
  if (hashToPositionIt == delayedRefSeedmerStatus.end()) {
    return mgsr::RefSeedmerExistStatus::NOT_EXIST;
  } else {
    return hashToPositionIt->second;
  }
  return mgsr::RefSeedmerExistStatus::NOT_EXIST;
}

void mgsr::mgsrPlacer::updateSeeds(
  MgsrLiteNode* node,
  std::unordered_set<uint64_t>& affectedSeedmers
) {
  const auto& curSeedDeltas = node->seedDeltas;
  const auto& seedInfos = liteTree->seedInfos;
  for (size_t i = 0; i < curSeedDeltas.size(); i++) {
    const auto [seedIndex, toDelete] = curSeedDeltas[i];
    if (i != curSeedDeltas.size() - 1 && seedInfos[seedIndex].startPos == seedInfos[curSeedDeltas[i + 1].first].startPos) {
      subSeedAtPosition(curSeedDeltas[i + 1].first, affectedSeedmers);
      i++;
    } else {
      if (toDelete) {
        delSeedAtPosition(seedIndex, affectedSeedmers);
      } else {
        addSeedAtPosition(seedIndex, affectedSeedmers);
      }
    }
  }
}

void mgsr::mgsrPlacer::backtrackSeeds(
  MgsrLiteNode* node,
  uint64_t nodeDfsIndex
) {
  const auto& curSeedDeltas = node->seedDeltas;
  const auto& seedInfos = liteTree->seedInfos;
  for (size_t i = 0; i < curSeedDeltas.size(); i++) {
    const auto [seedIndex, toDelete] = curSeedDeltas[i];
    if (i != curSeedDeltas.size() - 1 && seedInfos[seedIndex].startPos == seedInfos[curSeedDeltas[i + 1].first].startPos) {
      subSeedAtPosition(seedIndex);
      i++;
    } else {
      if (toDelete) {
        addSeedAtPosition(seedIndex);
      } else {
        delSeedAtPosition(seedIndex);
      }
    }
  }
}

void mgsr::mgsrPlacer::updateGapMap(
  MgsrLiteNode* node,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBacktracks,
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapBlocksBacktracks
) {
  const auto& curGapRunDeltas = node->gapRunDeltas;
  const auto& invertedBlocks = node->invertedBlocks;
  gapMapBacktracks.reserve(curGapRunDeltas.size() * 2);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> dummyGapMapUpdates;
  for (const auto [startPos, endPos, toGap] : curGapRunDeltas) {
    updateGapMapStep(gapMap, startPos, endPos, toGap, gapMapBacktracks, dummyGapMapUpdates, false);
  }
  gapMapBacktracks.shrink_to_fit();

  gapMapBlocksBacktracks.reserve(invertedBlocks.size() * 128);
  for (const auto& invertedBlock : invertedBlocks) {
    uint64_t beg = (uint64_t)liteTree->getBlockStartScalar(invertedBlock);
    uint64_t end = (uint64_t)liteTree->getBlockEndScalar(invertedBlock);
    invertGapMap(gapMap, {beg, end}, gapMapBlocksBacktracks, dummyGapMapUpdates);
  }
  gapMapBlocksBacktracks.shrink_to_fit();
}

void mgsr::squareEM::updateProps1() {    
  denoms.noalias() = probs * props0;
  inverseDenoms.noalias() = denoms.cwiseInverse();
  for (size_t i = 0; i < numNodes; ++i) {
    props1(i) = (readDuplicates.array() * (probs.col(i).array() * props0[i] * inverseDenoms.array())).sum() * invTotalWeight;
  }
}

void mgsr::squareEM::updateProps2() {
  denoms.noalias() = probs * props1;
  inverseDenoms.noalias() = denoms.cwiseInverse();
  for (size_t i = 0; i < numNodes; ++i) {
    props2(i) = (readDuplicates.array() * (probs.col(i).array() * props1[i] * inverseDenoms.array())).sum() * invTotalWeight;
  }
}

void mgsr::squareEM::updateProps(const Eigen::VectorXd& original, Eigen::VectorXd& result) {
  if (result.size() != numNodes) {
    result.resize(numNodes);
  }
  denoms.noalias() = probs * original;
  inverseDenoms.noalias() = denoms.cwiseInverse();
  for (size_t i = 0; i < numNodes; ++i) {
    result(i) = (readDuplicates.array() * (probs.col(i).array() * original[i] * inverseDenoms.array())).sum() * invTotalWeight;
  }
}

void mgsr::squareEM::normalizeProps(Eigen::VectorXd& props) {
  for (size_t i = 0; i < props.size(); ++i) {
    if (props(i) <= 0) {
      // std::cout << "prop[" << i << "]: " << props(i) << std::endl;
      props(i) = 1e-12;
    }
  }
  double sum = props.sum();
  props /= sum;
}

double mgsr::squareEM::getExp(const Eigen::VectorXd& props) {
  double llh = (readDuplicates.array() * (probs * props).array().log()).sum();
  return llh;
}

void mgsr::squareEM::resetIteration() {
  curIteration = 0;
}

void mgsr::squareEM::runSquareEM(uint64_t maximumIterations) {
  resetIteration();
  invTotalWeight = 1.0 / readDuplicates.array().sum();
  while (curIteration < maximumIterations) {
    props0.swap(props);
    updateProps1();
    normalizeProps(props1);
    updateProps2();
    normalizeProps(props2);

    r.noalias() = props1 - props0;
    v.noalias() = (props2 - props1) - r;

    double alpha = -r.norm() / v.norm();

    propsSq.noalias() = props0 - 2.0 * alpha * r + alpha * alpha * v;
    normalizeProps(propsSq);

    double llh2 = getExp(props2);
    double llhSq = getExp(propsSq);

    double difference = 0;
    if (llhSq > llh2 - eta) {
      props.swap(propsSq);
      difference = llhSq - llh;
      llh = llhSq;
    } else {
      props.swap(props2);
      difference = llh2 - llh;
      llh = llh2;
    }

    // std::cout << "\rEM iteration " << curIteration << " difference: " << std::fixed << std::setprecision(8) << difference << std::flush;
    // if (abs(difference) < 0.00001 || curIteration >= maximumIterations) {
    //   break;
    // }
    double max_change = (props - props0).array().abs().maxCoeff();
    std::cout << "\rEM iteration " << curIteration << " max_change: " << std::fixed << std::setprecision(8) << max_change << std::flush;

    if (max_change < maxChangeThreshold || curIteration >= maximumIterations) {
      break;
    }

    ++curIteration;
  }
}


bool mgsr::squareEM::removeLowPropNodes() {
  std::vector<size_t> passedIndices;
  passedIndices.reserve(props.size());

  for (size_t i = 0; i < props.size(); ++i) {
    if (props(i) >= propThresholdToRemove) {
      passedIndices.push_back(i);
    }
  }

  if (passedIndices.size() == props.size()) {
    return false;
  }

  std::vector<std::string> passedNodes;
  Eigen::MatrixXd passedProbs(probs.rows(), passedIndices.size());
  passedNodes.reserve(passedIndices.size());
  for (size_t i = 0; i < passedIndices.size(); ++i) {
    passedNodes.push_back(nodes[passedIndices[i]]);
    passedProbs.col(i) = probs.col(passedIndices[i]);
  }

  nodes.swap(passedNodes);
  probs.swap(passedProbs);
  numNodes = nodes.size();

  const size_t propsChunkSize = (numNodes + numThreads - 1) / numThreads;
  threadsRangeByProps.resize(numThreads);
  for (size_t i = 0; i < numThreads; ++i) {
    size_t start = i * propsChunkSize;
    size_t end = (i == numThreads - 1) ? numNodes : (i + 1) * propsChunkSize;
    if (start < numNodes) {
      threadsRangeByProps[i].first = start;
      threadsRangeByProps[i].second = end;
    }
  }

  props = Eigen::VectorXd::Constant(numNodes, 1.0 / static_cast<double>(numNodes));
  props0 = Eigen::VectorXd::Zero(numNodes);
  props1 = Eigen::VectorXd::Zero(numNodes);
  props2 = Eigen::VectorXd::Zero(numNodes);
  propsSq = Eigen::VectorXd::Zero(numNodes);
  r = Eigen::VectorXd::Zero(numNodes);
  v = Eigen::VectorXd::Zero(numNodes);

  return true;
}

void mgsr::mgsrPlacer::setReadScore(size_t readIndex, const int32_t score) {
  readScores[readIndex] = score;
  if (score > reads[readIndex].maxScore) {
    reads[readIndex].maxScore = score;
  }
}

void mgsr::mgsrPlacer::updateRefSeedmerStatus (
  size_t hash, mgsr::RefSeedmerChangeType& seedmerChangeType,
  mgsr::RefSeedmerExistStatus refSeedmerOldStatus,
  mgsr::RefSeedmerExistStatus refSeedmerNewStatus
) {
  switch (refSeedmerOldStatus) {
    case mgsr::RefSeedmerExistStatus::NOT_EXIST: {
      switch (refSeedmerNewStatus) {
        case mgsr::RefSeedmerExistStatus::NOT_EXIST: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::NOT_EXIST_TO_NOT_EXIST;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
        case mgsr::RefSeedmerExistStatus::EXIST_UNIQUE: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::NOT_EXIST_TO_EXIST_UNIQUE;
          if (allSeedmerHashesSet->find(hash) != allSeedmerHashesSet->end()) ++binaryOverlapKminmerCount;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
        case mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::NOT_EXIST_TO_EXIST_DUPLICATE;
          if (allSeedmerHashesSet->find(hash) != allSeedmerHashesSet->end()) ++binaryOverlapKminmerCount;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
      }
      break;
    }
    case mgsr::RefSeedmerExistStatus::EXIST_UNIQUE: {
      switch (refSeedmerNewStatus) {
        case mgsr::RefSeedmerExistStatus::NOT_EXIST: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_NOT_EXIST;
          if (allSeedmerHashesSet->find(hash) != allSeedmerHashesSet->end()) --binaryOverlapKminmerCount;
          delayedRefSeedmerStatus.erase(hash);
          break;
        }
        case mgsr::RefSeedmerExistStatus::EXIST_UNIQUE: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_EXIST_UNIQUE;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
        case mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_EXIST_DUPLICATE;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
      }
      break;
    }
    case mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE: {
      switch (refSeedmerNewStatus) {
        case mgsr::RefSeedmerExistStatus::NOT_EXIST: {  
          seedmerChangeType = mgsr::RefSeedmerChangeType::EXIST_DUPLICATE_TO_NOT_EXIST;
          if (allSeedmerHashesSet->find(hash) != allSeedmerHashesSet->end()) --binaryOverlapKminmerCount;
          delayedRefSeedmerStatus.erase(hash);
          break;
        }
        case mgsr::RefSeedmerExistStatus::EXIST_UNIQUE: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::EXIST_DUPLICATE_TO_EXIST_UNIQUE;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
        case mgsr::RefSeedmerExistStatus::EXIST_DUPLICATE: {
          seedmerChangeType = mgsr::RefSeedmerChangeType::EXIST_DUPLICATE_TO_EXIST_DUPLICATE;
          delayedRefSeedmerStatus[hash] = refSeedmerNewStatus;
          break;
        }
      }
      break;
    }
  }
}

void mgsr::mgsrPlacer::updateSeedmerChangesTypeFlag(mgsr::RefSeedmerChangeType seedmerChangeType, std::pair<bool, bool>& flags) {
  static const bool allUniqueToNonUnique[] = {
      false, // EXIST_UNIQUE_TO_EXIST_UNIQUE
      true,  // EXIST_UNIQUE_TO_EXIST_DUPLICATE 
      true,  // EXIST_UNIQUE_TO_NOT_EXIST
      false, // EXIST_DUPLICATE_TO_EXIST_UNIQUE
      false, // EXIST_DUPLICATE_TO_EXIST_DUPLICATE
      false, // EXIST_DUPLICATE_TO_NOT_EXIST
      false, // NOT_EXIST_TO_EXIST_UNIQUE
      false, // NOT_EXIST_TO_EXIST_DUPLICATE
      false  // NOT_EXIST_TO_NOT_EXIST
  };
  
  static const bool allNonUniqueToUnique[] = {
      false, // EXIST_UNIQUE_TO_EXIST_UNIQUE
      false, // EXIST_UNIQUE_TO_EXIST_DUPLICATE
      false, // EXIST_UNIQUE_TO_NOT_EXIST
      true,  // EXIST_DUPLICATE_TO_EXIST_UNIQUE
      false, // EXIST_DUPLICATE_TO_EXIST_DUPLICATE
      false, // EXIST_DUPLICATE_TO_NOT_EXIST
      true,  // NOT_EXIST_TO_EXIST_UNIQUE
      false, // NOT_EXIST_TO_EXIST_DUPLICATE
      false  // NOT_EXIST_TO_NOT_EXIST
  };

  // Use enum value as index (ensure your enum values are sequential)
  int index = static_cast<int>(seedmerChangeType);
  flags.first &= allUniqueToNonUnique[index];
  flags.second &= allNonUniqueToUnique[index];
}

inline uint64_t mgsr::mgsrPlacer::decodeBegFromMinichain(uint64_t minichain) {
  return (minichain >> 1) & 0x7FFFFFFF;
}

inline uint64_t mgsr::mgsrPlacer::decodeEndFromMinichain(uint64_t minichain) {
  return (minichain >> 32) & 0x7FFFFFFF;
}

inline bool mgsr::mgsrPlacer::decodeRevFromMinichain(uint64_t minichain) {
  return minichain & 1;
}

uint64_t mgsr::mgsrPlacer::extendMinichain(
  std::map<uint64_t, uint64_t>::const_iterator refPositionIt, const mgsr::Read& curRead, 
  uint64_t& curEnd, bool rev, uint64_t qidx, uint64_t c
) {
  const auto& seedInfos = liteTree->seedInfos;
  const auto& curSeedmerList = curRead.seedmersList;
  uint64_t currentQidx = qidx;
  uint64_t chainLength = c;
  auto currentRefPositionIt = refPositionIt;
  
  while (currentQidx < curSeedmerList.size() - 1) {
    // Get next seedmer
    const auto [nqhash, nqbeg, nqend, nqrev, nqidx] = curSeedmerList[currentQidx + 1];
    
    // Check if next hash exists and is unique
    auto nextHashToPositionIt = hashToPositionMap.find(nqhash);
    if (nextHashToPositionIt == hashToPositionMap.end() || 
        nextHashToPositionIt->second.size() != 1) {
        break;
    }
    
    auto nextRefPositionIt = *(nextHashToPositionIt->second.begin());
    bool nrev = nqrev != seedInfos[nextRefPositionIt->second].isReverse;
    
    // Check if orientation matches
    if (rev != nrev) {
        break;
    }
    
    // Check if positions are adjacent in the correct direction
    bool isAdjacent = false;
    if (rev && nextRefPositionIt == std::prev(currentRefPositionIt)) {
        isAdjacent = true;
    } else if (!rev && nextRefPositionIt == std::next(currentRefPositionIt)) {
        isAdjacent = true;
    }
    
    if (!isAdjacent) {
        break;
    }
    
    // Extend the chain
    curEnd = nqidx;
    currentQidx++;
    chainLength++;
    currentRefPositionIt = nextRefPositionIt;
  }
  
  return chainLength;
}


void mgsr::mgsrPlacer::initializeReadMinichains(size_t readIndex) {
  mgsr::Read& curRead = reads[readIndex];
  initializeReadMinichains(curRead);
}

void mgsr::mgsrPlacer::initializeReadMinichains(mgsr::Read& curRead) {
  auto& curMinichains = curRead.minichains;
  auto& curSeedmerList = curRead.seedmersList;
  curMinichains.clear();
  uint64_t i = 0;
  while (i < curSeedmerList.size()) {
    const auto& [hash, qbeg, qend, qrev, qidx] = curSeedmerList[i];
    uint64_t c = 1;
    auto hashToPositionIt = hashToPositionMap.find(hash);
    if (hashToPositionIt != hashToPositionMap.end()) {
      if (hashToPositionIt->second.size() < 2) {
        if (hashToPositionIt->second.size() != 1) {
          std::cerr << "Error: hash " << hash << " has multiple positions" << std::endl;
          exit(1);
        }
        auto curRefPositionIt = *(hashToPositionIt->second.begin());
        uint64_t curEnd = i;
        bool rev = qrev != liteTree->seedInfos[curRefPositionIt->second].isReverse;
        c = extendMinichain(curRefPositionIt, curRead, curEnd, rev, qidx, c);
        curMinichains.push_back({i, curEnd, rev});
      } else {
        // duplicates.insert(i);
      }
    }
    i += c; 
  }
  ++readMinichainsInitialized;
}

void mgsr::mgsrPlacer::extendChainRemoval(uint64_t& c, uint64_t& curEnd, const std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos, uint32_t lastSeedmerIndex) {
  while (curEnd != lastSeedmerIndex && c < affectedSeedmerInfos.size()) {
    uint64_t nextIndexOnSeedmerList = curEnd + 1;
    const auto& affectedSeedmerInfo = affectedSeedmerInfos[c];
    uint32_t affectedSeedmerIndex = affectedSeedmerInfo.index;
    
    if (nextIndexOnSeedmerList != affectedSeedmerIndex) {
        break;
    }
    
    mgsr::RefSeedmerChangeType seedmerChangeType = affectedSeedmerInfo.refSeedmerChangeType;
    
    if (seedmerChangeType == mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_EXIST_DUPLICATE ||
        seedmerChangeType == mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_NOT_EXIST
    ) {
      ++curEnd;
      ++c;
      // Continue the loop instead of recursive call
    } else {
      break;
    }
  }
}

void mgsr::mgsrPlacer::extendChainAddition(
  uint64_t& c, uint64_t& curEnd, const std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos,
  bool chainRev, std::map<uint64_t, uint64_t>::const_iterator refPositionIt, uint64_t readIndex
) {
  auto& curRead = reads[readIndex];
  const auto& curSeedmerList = curRead.seedmersList;
  
  while (curEnd != curSeedmerList.size() - 1 && c < affectedSeedmerInfos.size()) {
    const auto& affectedSeedmerInfo = affectedSeedmerInfos[c];
    uint32_t affectedSeedmerIndex = affectedSeedmerInfo.index;
    mgsr::RefSeedmerChangeType SeedmerChangeType = affectedSeedmerInfo.refSeedmerChangeType;
    bool refRev = affectedSeedmerInfo.rev;
    
    if (curEnd + 1 != affectedSeedmerIndex) {
      break;
    }
    
    if (SeedmerChangeType == mgsr::RefSeedmerChangeType::NOT_EXIST_TO_EXIST_UNIQUE ||
      SeedmerChangeType == mgsr::RefSeedmerChangeType::EXIST_DUPLICATE_TO_EXIST_UNIQUE
    ) {
      bool nextRev = refRev != curSeedmerList[affectedSeedmerIndex].rev;
      if (nextRev != chainRev) {
        break;
      }
      
      auto curRefPositionIt = hashToPositionMap.find(curSeedmerList[affectedSeedmerIndex].hash)->second.front();
      auto nextRefPositionIt = chainRev ? std::prev(refPositionIt) : std::next(refPositionIt);
      
      if (curRefPositionIt == nextRefPositionIt) {
        ++c;
        ++curEnd;
        refPositionIt = nextRefPositionIt; // Update refPositionIt for next iteration
        // Continue the loop instead of recursive call
      } else {
        break;
      }
    } else {
      break;
    }
  }
}

void mgsr::mgsrPlacer::extendChainUpdate(
  uint64_t& c, uint64_t& curEnd, const std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos,
  bool chainRev, std::map<uint64_t, uint64_t>::const_iterator refPositionIt, uint64_t readIndex
) {
  auto& curRead = reads[readIndex];
  const auto& curSeedmerList = curRead.seedmersList;
  
  while (curEnd != curSeedmerList.size() - 1 && c < affectedSeedmerInfos.size()) {
    const auto& affectedSeedmerInfo = affectedSeedmerInfos[c];
    uint32_t affectedSeedmerIndex = affectedSeedmerInfo.index;
    mgsr::RefSeedmerChangeType SeedmerChangeType = affectedSeedmerInfo.refSeedmerChangeType;
    bool refRev = affectedSeedmerInfo.rev;
    
    if (curEnd + 1 != affectedSeedmerIndex) {
      break;
    }
    
    if (SeedmerChangeType == mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_EXIST_UNIQUE) {
      bool nextRev = refRev != curSeedmerList[affectedSeedmerIndex].rev;
      if (nextRev != chainRev) {
        break;
      }
      
      auto curRefPositionIt = hashToPositionMap.find(curSeedmerList[affectedSeedmerIndex].hash)->second.front();
      auto nextRefPositionIt = chainRev ? std::prev(refPositionIt) : std::next(refPositionIt);
      
      if (curRefPositionIt == nextRefPositionIt) {
        ++c;
        ++curEnd;
        refPositionIt = nextRefPositionIt; // Update refPositionIt for next iteration
        // Continue the loop instead of recursive call
      } else {
        break;
      }
    } else {
      break;
    }
  }
}


bool mgsr::mgsrPlacer::colinearAdjacent(std::map<uint64_t, uint64_t>::const_iterator fromIterator, std::map<uint64_t, uint64_t>::const_iterator toIterator, bool fromRev, bool toRev) {
  if (fromRev != toRev) return false;
  auto nextIterator = fromRev ? std::prev(fromIterator) : std::next(fromIterator);
  return nextIterator == toIterator;
}

bool mgsr::mgsrPlacer::colinearAdjacent(std::map<uint64_t, uint64_t>::const_iterator fromIterator, std::map<uint64_t, uint64_t>::const_iterator toIterator, bool fromRev) {
  auto nextIterator = fromRev ? std::prev(fromIterator) : std::next(fromIterator);
  return nextIterator == toIterator;
}

void mgsr::mgsrPlacer::addToMinichains(const std::vector<readSeedmer>& curSeedmerList, std::vector<mgsr::Minichain>& curMinichains, mgsr::Minichain minichain) {
  // add and/or merge minichains
  bool addRangeRev = minichain.rev;
  uint64_t addRangeBeg = minichain.begIndex;
  uint64_t addRangeEnd = minichain.endIndex;
  if (curMinichains.empty()) {
    curMinichains.push_back(minichain);
  } else if (curMinichains.size() == 1) {
    ++readMinichainsAddedToSingleton;
    auto& originalMinichain = curMinichains[0];
    uint64_t originalMinichainBeg = originalMinichain.begIndex;
    uint64_t originalMinichainEnd = originalMinichain.endIndex;
    bool originalMinichainRev = originalMinichain.rev;
    
    if (addRangeEnd == originalMinichainBeg - 1 && originalMinichainBeg != 0) {
      // attempting to merge to the right
      if (addRangeRev != originalMinichainRev) {
        // different direction -> add without merging
        curMinichains.insert(curMinichains.begin(), minichain);
      } else {
        // same direction
        auto originalRangeBegKminmerPositionIt = hashToPositionMap.find(curSeedmerList[originalMinichainBeg].hash)->second.front();
        auto newRangeEndKminmerPositionIt = hashToPositionMap.find(curSeedmerList[addRangeEnd].hash)->second.front();
        if (colinearAdjacent(newRangeEndKminmerPositionIt, originalRangeBegKminmerPositionIt, addRangeRev)) {
          originalMinichain.begIndex = addRangeBeg;
        } else {
          curMinichains.insert(curMinichains.begin(), minichain);
        }
      }  
    } else if (addRangeBeg == originalMinichainEnd + 1) {
      // attempting to merge to the left
      if (addRangeRev != originalMinichainRev) {
        // different direction -> add without merging
        curMinichains.push_back(minichain);
      } else {
        // same direction
        auto originalRangeEndKminmerPositionIt = hashToPositionMap.find(curSeedmerList[originalMinichainEnd].hash)->second.front();
        auto newRangeBegKminmerPositionIt = hashToPositionMap.find(curSeedmerList[addRangeBeg].hash)->second.front();
        if (colinearAdjacent(originalRangeEndKminmerPositionIt, newRangeBegKminmerPositionIt, originalMinichainRev)) {
          originalMinichain.endIndex = addRangeEnd;
        } else {
          curMinichains.push_back(minichain);
        }
      }  
    } else {
      // add without merging
      if (addRangeEnd < originalMinichainBeg) {
        curMinichains.insert(curMinichains.begin(), minichain);
      } else {
        curMinichains.push_back(minichain);
      }
    }
  } else {
    ++readMinichainsAddedToMultiple;
    bool leftMinichainExists = false;
    bool rightMinichainExists = false;
    auto rightMinichainIt = std::upper_bound(curMinichains.begin(), curMinichains.end(), minichain, mgsr::compareMinichainByBeg);
    auto leftMinichainIt = curMinichains.end();
    if (rightMinichainIt != curMinichains.end()) {
      rightMinichainExists = true;
    }
    if (rightMinichainIt != curMinichains.begin()) {
      leftMinichainIt = std::prev(rightMinichainIt);
      leftMinichainExists = true;
    }

    bool mergeLeft = false;
    bool mergeRight = false;
    uint64_t leftMinichainBeg;
    uint64_t leftMinichainEnd;
    bool leftMinichainRev;
    uint64_t rightMinichainBeg;
    uint64_t rightMinichainEnd;
    bool rightMinichainRev;

    if (leftMinichainExists) {
      mgsr::Minichain& leftMinichain = *leftMinichainIt;
      leftMinichainBeg = leftMinichain.begIndex;
      leftMinichainEnd = leftMinichain.endIndex;
      leftMinichainRev = leftMinichain.rev;
      if (addRangeRev == leftMinichainRev && leftMinichainEnd + 1 == addRangeBeg) {
        auto newRangeBegKminmerPositionIt = hashToPositionMap.find(curSeedmerList[addRangeBeg].hash)->second.front();
        auto leftRangeEndKminmerPositionIt = hashToPositionMap.find(curSeedmerList[leftMinichainEnd].hash)->second.front();
        if (colinearAdjacent(leftRangeEndKminmerPositionIt, newRangeBegKminmerPositionIt, leftMinichainRev)) {
          mergeLeft = true;
        }
      }
    }  

    if (rightMinichainExists) {
      mgsr::Minichain& rightMinichain = *rightMinichainIt;
      rightMinichainBeg = rightMinichain.begIndex;
      rightMinichainEnd = rightMinichain.endIndex;
      rightMinichainRev = rightMinichain.rev;
      if (addRangeRev == rightMinichainRev && addRangeEnd + 1 == rightMinichainBeg) {
        auto newRangeEndKminmerPositionIt = hashToPositionMap.find(curSeedmerList[addRangeEnd].hash)->second.front();
        auto rightRangeBegKminmerPositionIt = hashToPositionMap.find(curSeedmerList[rightMinichainBeg].hash)->second.front();
        if (colinearAdjacent(newRangeEndKminmerPositionIt, rightRangeBegKminmerPositionIt, addRangeRev)) {
          mergeRight = true;
        }
      }
    }

    if (mergeLeft && mergeRight) {
      leftMinichainIt->endIndex = rightMinichainEnd;
      curMinichains.erase(rightMinichainIt);
    } else if (mergeLeft) {
      leftMinichainIt->endIndex = addRangeEnd;
    } else if (mergeRight) {
      rightMinichainIt->begIndex = addRangeBeg;
    } else {
      if (!leftMinichainExists) {
        curMinichains.insert(curMinichains.begin(), minichain);
      } else if (!rightMinichainExists) {
        curMinichains.push_back(minichain);
      } else {
        curMinichains.insert(leftMinichainIt+1, minichain);
      }
    }
  }
  return;
}


void mgsr::mgsrPlacer::removeFromMinichains(std::vector<mgsr::Minichain>& curMinichains, mgsr::Minichain minichain) {
  // remove from existing minichains
  uint64_t removeRangeBeg = minichain.begIndex;
  uint64_t removeRangeEnd = minichain.endIndex;
  if (curMinichains.size() == 1) {
    ++readMinichainsRemovedInplace;
    mgsr::Minichain& originalMinichain = curMinichains[0];
    uint64_t originalMinichainBeg = originalMinichain.begIndex;
    uint64_t originalMinichainEnd = originalMinichain.endIndex;
    bool originalMinichainRev = originalMinichain.rev;
    if (originalMinichainBeg == removeRangeBeg) {
      if (originalMinichainEnd == removeRangeEnd) {
        // remove the minichain
        curMinichains.clear();
      } else {
        // edit the beg
        originalMinichain.begIndex = removeRangeEnd + 1;
      }
    } else if (originalMinichainEnd == removeRangeEnd) {
      if (originalMinichainBeg == removeRangeBeg) {
        // remove the minichain
        curMinichains.clear();
      } else {
        // edit the end
        originalMinichain.endIndex = removeRangeBeg - 1;
      }
    } else {
      // need to split a minichain
      originalMinichain.endIndex = removeRangeBeg - 1;
      curMinichains.emplace_back(removeRangeEnd + 1, originalMinichainEnd, originalMinichainRev);
    }
  } else {
    ++readMinichainsRemovedFromMultiple;
    auto currentOriginalMinichainIt = std::upper_bound(curMinichains.begin(), curMinichains.end(), minichain, mgsr::compareMinichainByBeg);
    --currentOriginalMinichainIt;
    auto& currentOriginalMinichain = *currentOriginalMinichainIt;
    uint64_t currentOriginalMinichainBeg = currentOriginalMinichain.begIndex;
    uint64_t currentOriginalMinichainEnd = currentOriginalMinichain.endIndex;
    bool currentOriginalMinichainRev = currentOriginalMinichain.rev;
    
    if (removeRangeEnd > currentOriginalMinichainEnd) {
      uint32_t numToErase = 0;
      auto it = currentOriginalMinichainIt;
      if (currentOriginalMinichainBeg == removeRangeBeg) {
        // remove the minichain
        ++numToErase;
        ++currentOriginalMinichainIt;
      } else {
        // edit the end then step right
        currentOriginalMinichainIt->endIndex = removeRangeBeg - 1;
        ++currentOriginalMinichainIt;
        ++it;
      }

      currentOriginalMinichainEnd = currentOriginalMinichainIt->endIndex;
      while (currentOriginalMinichainIt != curMinichains.end() && currentOriginalMinichainEnd <= removeRangeEnd) {
        ++numToErase;
        ++currentOriginalMinichainIt;
        currentOriginalMinichainEnd = currentOriginalMinichainIt->endIndex;
      }

      currentOriginalMinichainBeg = currentOriginalMinichainIt->begIndex;
      if (currentOriginalMinichainIt != curMinichains.end() && currentOriginalMinichainBeg <= removeRangeEnd) {
        currentOriginalMinichainIt->begIndex = removeRangeEnd + 1;
      }

      for (size_t i = 0; i < numToErase; ++i) {
        auto curit = it+i;
        uint64_t curitMinichainBeg = curit->begIndex;
        uint64_t curitMinichainEnd = curit->endIndex;
        bool curitMinichainRev = curit->rev;
      }
      curMinichains.erase(it, it+numToErase);
    } else {
      // same as if there is one minichain but instead of clear, use erase
      if (currentOriginalMinichainBeg == removeRangeBeg) {
        if (currentOriginalMinichainEnd == removeRangeEnd) {
          // erase the minichain
          curMinichains.erase(currentOriginalMinichainIt);
        } else {
          // edit the beg
          currentOriginalMinichainIt->begIndex = removeRangeEnd + 1;
        }
      } else if (currentOriginalMinichainEnd  == removeRangeEnd) {
        if (currentOriginalMinichainBeg == removeRangeBeg) {
          // erase the minichain
          curMinichains.erase(currentOriginalMinichainIt);
        } else {
          // edit the end
          currentOriginalMinichainIt->endIndex = removeRangeBeg - 1;
        }
      } else {
        currentOriginalMinichainIt->endIndex = removeRangeBeg - 1;
        curMinichains.insert(currentOriginalMinichainIt+1, {removeRangeEnd + 1, currentOriginalMinichainEnd, currentOriginalMinichainRev});
      }
    }
  }
}

void mgsr::mgsrPlacer::updateMinichainsMixed(size_t readIndex, const std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos) {
  mgsr::Read& curRead = reads[readIndex];
  const auto& curSeedmerList = curRead.seedmersList;
  auto& curMinichains = curRead.minichains;

  if (affectedSeedmerInfos.size() > minichainsToUpdate.size()) minichainsToUpdate.resize(affectedSeedmerInfos.size());
  if (affectedSeedmerInfos.size() > minichainsToRemove.size()) minichainsToRemove.resize(affectedSeedmerInfos.size());
  if (affectedSeedmerInfos.size() > minichainsToAdd.size())    minichainsToAdd.resize(affectedSeedmerInfos.size());

  uint32_t i = 0;
  uint32_t curMinichainsToUpdateIndex = 0;
  uint32_t curMinichainsToRemoveIndex = 0;
  uint32_t curMinichainsToAddIndex = 0;

  while (i < affectedSeedmerInfos.size()) {
    auto [affectedSeedmerIndex, SeedmerChangeType, refRev] = affectedSeedmerInfos[i];
    uint64_t c = i + 1;
    uint64_t curEnd = affectedSeedmerIndex;

    switch (SeedmerChangeType) {
      case mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_EXIST_DUPLICATE:
      case mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_NOT_EXIST:
        // match to no match -> remove from minichains
        {
          extendChainRemoval(c, curEnd, affectedSeedmerInfos, curSeedmerList.size() - 1);
          mgsr::Minichain minichain = {affectedSeedmerIndex, curEnd, false};
          auto& curMinichainsToRemove = minichainsToRemove[curMinichainsToRemoveIndex];
          curMinichainsToRemove.first = minichain;
          ++curMinichainsToRemoveIndex;
          i += (curEnd - affectedSeedmerIndex + 1);
          break;
        }
      case mgsr::RefSeedmerChangeType::EXIST_DUPLICATE_TO_EXIST_UNIQUE:
      case mgsr::RefSeedmerChangeType::NOT_EXIST_TO_EXIST_UNIQUE:
        // no match to match -> create minichains
        {
          bool rev = refRev != curSeedmerList[affectedSeedmerIndex].rev;
          auto positionItFromCurrentHash = hashToPositionMap.at(curSeedmerList[affectedSeedmerIndex].hash).front();
          extendChainAddition(c, curEnd, affectedSeedmerInfos, rev, positionItFromCurrentHash, readIndex);
          // encode minichain_t
          mgsr::Minichain minichain = {affectedSeedmerIndex, curEnd, rev};
          auto& curMinichainToAdd = minichainsToAdd[curMinichainsToAddIndex];
          curMinichainToAdd.first = minichain;
          ++curMinichainsToAddIndex;
          i += (curEnd - affectedSeedmerIndex + 1);
          break;
        }
      case mgsr::RefSeedmerChangeType::EXIST_UNIQUE_TO_EXIST_UNIQUE:
        // match to match -> update minichains
        {
          bool rev = refRev != curSeedmerList[affectedSeedmerIndex].rev;
          auto positionItFromCurrentHash = hashToPositionMap.at(curSeedmerList[affectedSeedmerIndex].hash).front();
          extendChainUpdate(c, curEnd, affectedSeedmerInfos, rev, positionItFromCurrentHash, readIndex);
          // encode minichain_t
          mgsr::Minichain minichain = {affectedSeedmerIndex, curEnd, rev};
          auto& curMinichainToUpdate = minichainsToUpdate[curMinichainsToUpdateIndex];
          curMinichainToUpdate.first = minichain;
          ++curMinichainsToUpdateIndex;
          i += (curEnd - affectedSeedmerIndex + 1);
          break;
        }
        break;
      default:
        // no match to no match -> do nothing
        ++i;
        break;
    }
  }

  // remove minichains first for unique to non-unique changes
  for (size_t i = 0; i < curMinichainsToRemoveIndex; ++i) {
    const auto& [minichain, _] = minichainsToRemove[i];
    removeFromMinichains(curMinichains, minichain);
  }

  // remove minichains then add minichains for unique to unique changes
  for (size_t i = 0; i < curMinichainsToUpdateIndex; ++i) {
    const auto& [minichain, _] = minichainsToUpdate[i];
    removeFromMinichains(curMinichains, minichain);
  }

  for (size_t i = 0; i < curMinichainsToUpdateIndex; ++i) {
    const auto& [minichain, _] = minichainsToUpdate[i];
    addToMinichains(curSeedmerList, curMinichains, minichain);
  }

  // finally add minichains for non-unique to unique changes
  for (size_t i = 0; i < curMinichainsToAddIndex; ++i) {
    const auto& [minichain, _] = minichainsToAdd[i];
    addToMinichains(curSeedmerList, curMinichains, minichain);
  }
}

void mgsr::mgsrPlacer::updateMinichains(size_t readIndex, const std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos, bool allUniqueToNonUnique, bool allNonUniqueToUnique) {
  mgsr::Read& curRead = reads[readIndex];
  const auto& curSeedmerList = curRead.seedmersList;
  auto& curMinichains = curRead.minichains;
  // decode affectedSeedmerIndexCode:
  // 0th bit is 1 -> reversed
  // 1st to 8th bit is the seedmerChangeType
  // 9th to 40th bit is the affectedSeedmerIndex
  if (affectedSeedmerInfos.size() > minichainsToUpdate.size()) {
    minichainsToUpdate.resize(affectedSeedmerInfos.size());
  }
  uint32_t i = 0;
  uint32_t curMinichainsToUpdateIndex = 0;
  while (i < affectedSeedmerInfos.size()) {
    const auto& affectedSeedmerInfo = affectedSeedmerInfos[i];
    auto [affectedSeedmerIndex, SeedmerChangeType, refRev] = affectedSeedmerInfo;

    uint64_t c = i + 1;
    uint64_t curEnd = affectedSeedmerIndex;
    
    if (allUniqueToNonUnique) {
      // match to no match -> remove from minichains
      /*chaining debug*/if (false) std::cout << "\tRead is allUniqueToNonUnique... removing minichains" << std::endl;
      extendChainRemoval(c, curEnd, affectedSeedmerInfos, curSeedmerList.size() - 1);
      // encode minichain_t
      mgsr::Minichain minichain = {affectedSeedmerIndex, curEnd, false};
      /*chaining debug*/if (false) std::cout << "\t\tRemoving minichain: " << affectedSeedmerIndex << " " << curEnd << " " << 0ULL << "... encoded to " << minichain.begIndex << " " << minichain.endIndex << " " << minichain.rev << std::endl;
      auto& curMinichainToUpdate = minichainsToUpdate[curMinichainsToUpdateIndex];
      curMinichainToUpdate.first = minichain;
      curMinichainToUpdate.second = false;
      ++curMinichainsToUpdateIndex;
      ++readMinichainsRemoved;
    } else if (allNonUniqueToUnique) {
      // no match to match -> create minichains
      /*chaining debug*/if (false) std::cout << "\tRead is allNonUniqueToUnique... adding minichains" << std::endl;
      bool rev = refRev != curSeedmerList[affectedSeedmerIndex].rev;
      auto positionItFromCurrentHash = hashToPositionMap.at(curSeedmerList[affectedSeedmerIndex].hash).front();
      extendChainAddition(c, curEnd, affectedSeedmerInfos, rev, positionItFromCurrentHash, readIndex);
      // encode minichain_t
      mgsr::Minichain minichain = {affectedSeedmerIndex, curEnd, rev};
      /*chaining debug*/if (false) std::cout << "\t\tAdding minichain: " << affectedSeedmerIndex << " " << curEnd << " " << rev << "... encoded to " << minichain.begIndex << " " << minichain.endIndex << " " << minichain.rev << std::endl;
      auto& curMinichainToUpdate = minichainsToUpdate[curMinichainsToUpdateIndex];
      curMinichainToUpdate.first = minichain;
      curMinichainToUpdate.second = true;
      ++curMinichainsToUpdateIndex;
      ++readMinichainsAdded;
    } else {
      std::cerr << "Error: allUniqueToNonUnique and allNonUniqueToUnique are both false" << std::endl;
      exit(1);
    }
    i += (curEnd - affectedSeedmerIndex + 1);
  }
  /*chaining debug*/if (false) std::cout << "\tGenerated " << curMinichainsToUpdateIndex << " minichains to " << (allUniqueToNonUnique ? "remove" : "add") << std::endl;

  if (allUniqueToNonUnique) {
    for (uint32_t i = 0; i < curMinichainsToUpdateIndex; ++i) {
      const auto& [minichain, toAdd] = minichainsToUpdate[i];
      /*chaining debug*/if (false) std::cout << "\t\tRemoving minichain: " << minichain.begIndex << " " << minichain.endIndex << " " << minichain.rev << std::endl;
      removeFromMinichains(curMinichains, minichain);
      /*chaining debug*/if (false) {
        std::cout << "\t\tNow minichains are: ";
        for (const auto& minichain : curMinichains) {
          std::cout << "(" << minichain.begIndex << ", " << minichain.endIndex << ", " << minichain.rev << ") ";
        }
        std::cout << std::endl;
      }
    }
  } else if (allNonUniqueToUnique) {
    bool wasEmpty = curMinichains.empty();
    for (uint32_t i = 0; i < curMinichainsToUpdateIndex; ++i) {
      const auto& [minichain, toAdd] = minichainsToUpdate[i];
      /*chaining debug*/if (false) std::cout << "\t\tAdding minichain: " << minichain.begIndex << " " << minichain.endIndex << " " << minichain.rev << std::endl;
      if (wasEmpty) {
        curMinichains.push_back(minichain);
        ++readMinichainsAddedToEmpty;
      } else {
        addToMinichains(curSeedmerList, curMinichains, minichain);
      }
      /*chaining debug*/if (false) {
        std::cout << "\t\tNow minichains are: ";
        for (const auto& minichain : curMinichains) {
          std::cout << "(" << minichain.begIndex << ", " << minichain.endIndex << ", " << minichain.rev << ") ";
        }
        std::cout << std::endl;
      }
    }
  }
}

void mgsr::mgsrPlacer::fillReadToAffectedSeedmerIndex(
  absl::flat_hash_map<uint32_t, std::pair<std::vector<mgsr::affectedSeedmerInfo>, std::pair<bool, bool>>>& readToAffectedSeedmerIndex,
  const std::unordered_set<uint64_t>& affectedSeedmers
) {
  absl::flat_hash_map<uint32_t, uint32_t> readAffectedSeedmerCount;
  for (const auto& hash : affectedSeedmers) {
    auto affectedSeedmerToReads = seedmerToReads.find(hash);
    if (affectedSeedmerToReads == seedmerToReads.end()) continue;
    for (const auto& [readIndex, _] : affectedSeedmerToReads->second) {
      ++readAffectedSeedmerCount[readIndex];
    }
  }
  readToAffectedSeedmerIndex.reserve(readAffectedSeedmerCount.size());
  for (const auto& hash : affectedSeedmers) {
    mgsr::RefSeedmerChangeType seedmerChangeType;
    mgsr::RefSeedmerExistStatus refSeedmerOldStatus = getDelayedRefSeedmerExistStatus(hash);
    mgsr::RefSeedmerExistStatus refSeedmerNewStatus = getCurrentRefSeedmerExistStatus(hash);
    updateRefSeedmerStatus(hash, seedmerChangeType, refSeedmerOldStatus, refSeedmerNewStatus);

    bool refRev = false;
    if (refSeedmerNewStatus == mgsr::RefSeedmerExistStatus::EXIST_UNIQUE) {
      refRev = liteTree->seedInfos[hashToPositionMap[hash].front()->second].isReverse;
    }

    auto affectedseedmerToReads = seedmerToReads.find(hash);
    if (affectedseedmerToReads == seedmerToReads.end()) continue;
    for (const auto& [readIndex, affectedSeedmerIndex] : affectedseedmerToReads->second) {
      if (reads[readIndex].readType != mgsr::ReadType::PASS) continue;
      auto [affectedSeedmerIndexVectorIt, inserted] = readToAffectedSeedmerIndex.try_emplace(readIndex, std::vector<mgsr::affectedSeedmerInfo>{}, std::pair<bool, bool>{});
      auto& affectedSeedmerIndexVector = affectedSeedmerIndexVectorIt->second;
      if (affectedSeedmerIndexVector.first.empty()) {
        affectedSeedmerIndexVector.second.first = true;
        affectedSeedmerIndexVector.second.second = true;
        affectedSeedmerIndexVector.first.reserve(readAffectedSeedmerCount[readIndex]);
      }
      affectedSeedmerIndexVector.first.emplace_back(affectedSeedmerIndex, seedmerChangeType, refRev);
      updateSeedmerChangesTypeFlag(seedmerChangeType, affectedSeedmerIndexVector.second);
    }
  }
}





uint64_t mgsr::mgsrPlacer::getRefSeedmerBegFromHash(const size_t hash) const {
  auto hashToPositionIt = hashToPositionMap.find(hash);
  return hashToPositionIt->second.front()->first;
}


uint64_t mgsr::mgsrPlacer::getRefSeedmerEndFromHash(const size_t hash) const {
  auto hashToPositionIt = hashToPositionMap.find(hash);
  auto positionIt = hashToPositionIt->second.front();
  return liteTree->seedInfos[positionIt->second].endPos;
}

inline int32_t absDifference(const uint32_t a, const uint32_t b) {
  return a > b ? a - b : b - a;
}

int32_t mgsr::mgsrPlacer::getLocalGap(const uint32_t a, const uint32_t b) const {
  if (a == b) return 0;
  uint32_t left = a < b ? a : b;
  uint32_t right = a > b ? a : b;
  auto curIt = gapMap.upper_bound(left);

  if (curIt == gapMap.end() || curIt->first > right) {
    return right - left;
  }

  uint32_t gap = curIt->first - left - 1;
  auto prevIt = curIt;
  ++curIt;
  while (true) {
    if (curIt == gapMap.end()) {
      gap += right - prevIt->second;
      break;
    } else if (right > curIt->first) {
      gap += curIt->first - prevIt->second - 1;
    } else if (right < curIt->first) {
      gap += right - prevIt->second;
      break;
    } else if (right == curIt->first) {
      std::cerr << "Error: right == curIt->first" << std::endl;
      exit(1);
    }
    prevIt = curIt;
    ++curIt;
  }
  return gap;
}

bool mgsr::mgsrPlacer::isColinearFromMinichains(
  mgsr::Read& curRead, const mgsr::Minichain& minichain1, const mgsr::Minichain& minichain2
) {
  const bool minichainRev = minichain1.rev;
  const auto& seedmersList = curRead.seedmersList;
  const auto& beg1seedmer = seedmersList[minichain1.begIndex];
  const auto& end1seedmer = seedmersList[minichain1.endIndex];
  const auto& beg2seedmer = seedmersList[minichain2.begIndex];
  const auto& end2seedmer = seedmersList[minichain2.endIndex];

  const auto qbeg2 = beg2seedmer.begPos; // qbeg2: minichain2 -> query -> startKminmer
  const auto qend1 = end1seedmer.endPos; // qend1: minichain1 -> query -> endKminmer

  if (!minichainRev) {
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg1CoordInfoCacheIt = hashCoordInfoCacheTable.find(beg1seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg2CoordInfoCacheIt = hashCoordInfoCacheTable.find(beg2seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rend1CoordInfoCacheIt = hashCoordInfoCacheTable.find(end1seedmer.hash);
   
    // rbeg1: minichain1 -> ref -> startKminmer
    auto& rglobalbeg1 = rbeg1CoordInfoCacheIt->second.rGlobalBeg;
    if (rbeg1CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg1CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg1 = getRefSeedmerBegFromHash(beg1seedmer.hash);
    }
  
    // rbeg2: minichain2 -> ref -> startKminmer
    auto& rglobalbeg2 = rbeg2CoordInfoCacheIt->second.rGlobalBeg;
    if (rbeg2CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg2CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg2 = getRefSeedmerBegFromHash(beg2seedmer.hash);
    }
    // rend1: minichain1 -> ref -> endKminmer
    auto& rglobalend1 = rend1CoordInfoCacheIt->second.rGlobalEnd;
    if (rend1CoordInfoCacheIt->second.endDfsIndex != curDfsIndex) {
      rend1CoordInfoCacheIt->second.endDfsIndex = curDfsIndex;
      rglobalend1 = getRefSeedmerEndFromHash(end1seedmer.hash);
    }

    int32_t qgap = absDifference(qbeg2, qend1);
    int32_t rgap = getLocalGap(rglobalbeg2, rglobalend1);
    if (rglobalbeg1 < rglobalbeg2 && absDifference(qgap, rgap) < maximumGap) return true;
  } else {
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg1CoordInfoCacheIt = hashCoordInfoCacheTable.find(end1seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg2CoordInfoCacheIt = hashCoordInfoCacheTable.find(end2seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rend2CoordInfoCacheIt = hashCoordInfoCacheTable.find(beg2seedmer.hash);
   
    auto& rglobalbeg1 = rbeg1CoordInfoCacheIt->second.rGlobalBeg;
    if (rbeg1CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg1CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg1 = getRefSeedmerBegFromHash(end1seedmer.hash);
    }

    auto& rglobalbeg2 = rbeg2CoordInfoCacheIt->second.rGlobalBeg;
    if (rbeg2CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg2CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg2 = getRefSeedmerBegFromHash(end2seedmer.hash);
    }

    auto& rglobalend2 = rend2CoordInfoCacheIt->second.rGlobalEnd;
    if (rend2CoordInfoCacheIt->second.endDfsIndex != curDfsIndex) {
      rend2CoordInfoCacheIt->second.endDfsIndex = curDfsIndex;
      rglobalend2 = getRefSeedmerEndFromHash(beg2seedmer.hash);
    }

    int32_t qgap = absDifference(qbeg2, qend1);
    int32_t rgap = getLocalGap(rglobalbeg1, rglobalend2);
    if (rglobalbeg2 < rglobalbeg1 && absDifference(qgap, rgap) < maximumGap) return true;
  }
  return false;
}

bool mgsr::mgsrPlacer::isColinearFromMinichains(
  mgsr::Read& curRead, const mgsr::Minichain& minichain1, const mgsr::Minichain& minichain2,
  const std::map<uint64_t, uint64_t>& degapCoordIndex,
  const std::map<uint64_t, uint64_t>& regapCoordIndex
) {
  const bool minichainRev = minichain1.rev;
  const auto& seedmersList = curRead.seedmersList;
  const auto& beg1seedmer = seedmersList[minichain1.begIndex];
  const auto& end1seedmer = seedmersList[minichain1.endIndex];
  const auto& beg2seedmer = seedmersList[minichain2.begIndex];
  const auto& end2seedmer = seedmersList[minichain2.endIndex];

  const auto qbeg2 = beg2seedmer.begPos; // qbeg2: minichain2 -> query -> startKminmer
  const auto qend1 = end1seedmer.endPos; // qend1: minichain1 -> query -> endKminmer

  if (!minichainRev) {
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg1CoordInfoCacheIt = hashCoordInfoCacheTable.find(beg1seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg2CoordInfoCacheIt = hashCoordInfoCacheTable.find(beg2seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rend1CoordInfoCacheIt = hashCoordInfoCacheTable.find(end1seedmer.hash);
   
    // rbeg1: minichain1 -> ref -> startKminmer
    auto& rglobalbeg1 = rbeg1CoordInfoCacheIt->second.rGlobalBeg;
    auto& rbeg1       = rbeg1CoordInfoCacheIt->second.rLocalBeg;
    if (rbeg1CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg1CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg1 = getRefSeedmerBegFromHash(beg1seedmer.hash);
      rbeg1 = mgsr::degapGlobal(rglobalbeg1, degapCoordIndex);
    }
  
    // rbeg2: minichain2 -> ref -> startKminmer
    auto& rglobalbeg2 = rbeg2CoordInfoCacheIt->second.rGlobalBeg;
    auto& rbeg2       = rbeg2CoordInfoCacheIt->second.rLocalBeg;
    if (rbeg2CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg2CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg2 = getRefSeedmerBegFromHash(beg2seedmer.hash);
      rbeg2 = mgsr::degapGlobal(rglobalbeg2, degapCoordIndex);
    }
    // rend1: minichain1 -> ref -> endKminmer
    auto& rglobalend1 = rend1CoordInfoCacheIt->second.rGlobalEnd;
    auto& rend1       = rend1CoordInfoCacheIt->second.rLocalEnd;
    if (rend1CoordInfoCacheIt->second.endDfsIndex != curDfsIndex) {
      rend1CoordInfoCacheIt->second.endDfsIndex = curDfsIndex;
      rglobalend1 = getRefSeedmerEndFromHash(end1seedmer.hash);
      rend1 = mgsr::degapGlobal(rglobalend1, degapCoordIndex);
    }

    int32_t qgap = absDifference(qbeg2, qend1);
    int32_t rgap = absDifference(rbeg2, rend1);
    if (rbeg1 < rbeg2 && absDifference(qgap, rgap) < maximumGap) return true;
  } else {
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg1CoordInfoCacheIt = hashCoordInfoCacheTable.find(end1seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rbeg2CoordInfoCacheIt = hashCoordInfoCacheTable.find(end2seedmer.hash);
    absl::flat_hash_map<size_t, mgsr::hashCoordInfoCache>::iterator rend2CoordInfoCacheIt = hashCoordInfoCacheTable.find(beg2seedmer.hash);
   
    auto& rglobalbeg1 = rbeg1CoordInfoCacheIt->second.rGlobalBeg;
    auto& rbeg1       = rbeg1CoordInfoCacheIt->second.rLocalBeg;
    if (rbeg1CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg1CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg1 = getRefSeedmerBegFromHash(end1seedmer.hash);
      rbeg1 = mgsr::degapGlobal(rglobalbeg1, degapCoordIndex);
    }

    auto& rglobalbeg2 = rbeg2CoordInfoCacheIt->second.rGlobalBeg;
    auto& rbeg2       = rbeg2CoordInfoCacheIt->second.rLocalBeg;
    if (rbeg2CoordInfoCacheIt->second.begDfsIndex != curDfsIndex) {
      rbeg2CoordInfoCacheIt->second.begDfsIndex = curDfsIndex;
      rglobalbeg2 = getRefSeedmerBegFromHash(end2seedmer.hash);
      rbeg2 = mgsr::degapGlobal(rglobalbeg2, degapCoordIndex);
    }

    auto& rglobalend2 = rend2CoordInfoCacheIt->second.rGlobalEnd;
    auto& rend2       = rend2CoordInfoCacheIt->second.rLocalEnd;
    if (rend2CoordInfoCacheIt->second.endDfsIndex != curDfsIndex) {
      rend2CoordInfoCacheIt->second.endDfsIndex = curDfsIndex;
      rglobalend2 = getRefSeedmerEndFromHash(beg2seedmer.hash);
      rend2 = mgsr::degapGlobal(rglobalend2, degapCoordIndex);
    }


    int32_t qgap = absDifference(qbeg2, qend1);
    int32_t rgap = absDifference(rbeg1, rend2);
    if (rbeg2 < rbeg1 && absDifference(qgap, rgap) < maximumGap) return true;
  }
  return false;
}

int32_t mgsr::mgsrPlacer::getReadPseudoScore(mgsr::Read& curRead) {
  const auto& minichains = curRead.minichains;

  int32_t pseudoChainScore = 0;
  if (minichains.empty()) {
    return 0;
  } else if (minichains.size() == 1) {
    return  minichains[0].getLength();
  } else {
    // find longest minichain
    uint64_t longestMinichainLength = 0;
    int32_t longestMinichainIndex = -1;
    for (int32_t i = 0; i < minichains.size(); ++i) {
      const auto currentMinichainLength = minichains[i].getLength();
      if (currentMinichainLength > longestMinichainLength) {
        longestMinichainIndex = i;
        longestMinichainLength = currentMinichainLength;
      }
    }

    const auto& longestMinichain = minichains[longestMinichainIndex];
    bool longestMinichainIsReversed = longestMinichain.rev;

    for (int32_t i = 0; i < minichains.size(); ++i) {
      const mgsr::Minichain& currentMinichain = minichains[i];
      bool currentMinichainIsReversed = currentMinichain.rev;
      if (i == longestMinichainIndex) {
        pseudoChainScore += longestMinichain.getLength();
        continue;
      }

      if (currentMinichainIsReversed != longestMinichainIsReversed) {
        continue;
      }

      if (longestMinichainIndex < i) {
        if (isColinearFromMinichains(curRead, longestMinichain, currentMinichain)) {
          pseudoChainScore += currentMinichain.getLength();
        }
      } else if (longestMinichainIndex > i) {
        if (isColinearFromMinichains(curRead, currentMinichain, longestMinichain)) {
          pseudoChainScore += currentMinichain.getLength();
        }
      }
    }
  }
  return pseudoChainScore;
}


std::vector<uint32_t> mgsr::mgsrPlacer::getScoresAtNode(const std::string& nodeId) const {
  MgsrLiteNode* currentNode = liteTree->allLiteNodes[nodeId];

  // Get the path from root the current node
  std::vector<MgsrLiteNode*> nodePath;
  while (currentNode->parent != nullptr) {
    nodePath.push_back(currentNode);
    currentNode = currentNode->parent;
  }
  nodePath.push_back(currentNode);
  std::reverse(nodePath.begin(), nodePath.end());

  std::vector<uint32_t> curNodeScores(readScores.size(), 0);
  for (const auto& node : nodePath) {
    for (const auto& scoreDelta : node->readScoreDeltas[threadId]) {
      curNodeScores[scoreDelta.readIndex] = scoreDelta.scoreDelta;
    }
  }
  return curNodeScores;
}

void mgsr::mgsrPlacer::getScoresAtNode(const std::string& nodeId, std::vector<uint32_t>& curNodeScores) const {
  MgsrLiteNode* currentNode = liteTree->allLiteNodes[nodeId];

  // Get the path from root the current node
  std::vector<MgsrLiteNode*> nodePath;
  while (currentNode->parent != nullptr) {
    nodePath.push_back(currentNode);
    currentNode = currentNode->parent;
  }
  nodePath.push_back(currentNode);
  std::reverse(nodePath.begin(), nodePath.end());

  for (const auto& node : nodePath) {
    for (const auto& scoreDelta : node->readScoreDeltas[threadId]) {
      curNodeScores[scoreDelta.readIndex] = scoreDelta.scoreDelta;
    }
  }

}

bool mgsr::mgsrPlacer::identicalReadScores(const std::string& node1, const std::string& node2, bool fast_mode) const {
  MgsrLiteNode* currentNode1 = liteTree->allLiteNodes[node1];
  MgsrLiteNode* currentNode2 = liteTree->allLiteNodes[node2];
  std::vector<MgsrLiteNode*> nodePath1;
  std::vector<MgsrLiteNode*> nodePath2;

  while (currentNode1->parent != nullptr) {
    nodePath1.push_back(currentNode1);
    currentNode1 = currentNode1->parent;
  }
  nodePath1.push_back(currentNode1);
  std::reverse(nodePath1.begin(), nodePath1.end());

  while (currentNode2->parent != nullptr) {
    nodePath2.push_back(currentNode2);
    currentNode2 = currentNode2->parent;
  }
  nodePath2.push_back(currentNode2);
  std::reverse(nodePath2.begin(), nodePath2.end());

  size_t lcaIndex = 0;
  while (lcaIndex < nodePath1.size() && lcaIndex < nodePath2.size() && nodePath1[lcaIndex] == nodePath2[lcaIndex]) {
    ++lcaIndex;
  }
  lcaIndex -= 1;

  std::vector<uint32_t> lcaScores = getScoresAtNode(nodePath1[lcaIndex]->identifier);

  std::unordered_map<size_t, std::pair<uint32_t, uint32_t>> changedReads;
  for (size_t i = lcaIndex + 1; i < nodePath1.size(); ++i) {
    MgsrLiteNode* currNode = nodePath1[i];
    for (const auto& scoreDelta : currNode->readScoreDeltas[threadId]) {
      changedReads[scoreDelta.readIndex] = {scoreDelta.scoreDelta, lcaScores[scoreDelta.readIndex]};
    }
  }

  for (size_t i = lcaIndex + 1; i < nodePath2.size(); ++i) {
    MgsrLiteNode* currNode = nodePath2[i];
    for (const auto& scoreDelta : currNode->readScoreDeltas[threadId]) {
      if (changedReads.find(scoreDelta.readIndex) != changedReads.end()) {
        changedReads[scoreDelta.readIndex].second = scoreDelta.scoreDelta;
      } else {
        changedReads[scoreDelta.readIndex] = {lcaScores[scoreDelta.readIndex], scoreDelta.scoreDelta};
      }
    }
  }

  for (const auto& [readIndex, scoreDeltaInfo] : changedReads) {
    auto [score1, score2] = scoreDeltaInfo;
    if (score1 != score2) {
      return false;
    }
  }
  return true;
}

void mgsr::mgsrPlacer::mergeIdenticalNodes(const std::unordered_set<std::string>& identicalGroup) {
  std::unordered_set<std::string> seenNodes;
  std::unordered_set<std::string> unseenNodes = identicalGroup;
  for (const std::string& currNode : identicalGroup) {
    if (seenNodes.find(currNode) != seenNodes.end()) continue;
    seenNodes.insert(currNode);
    unseenNodes.erase(currNode);
    std::unordered_set<std::string> identicals;

    for (const std::string& potentialIdentical : unseenNodes) {
      if (identicalReadScores(currNode, potentialIdentical)) {
        identicals.insert(potentialIdentical);
        identicalGroups[currNode].push_back(potentialIdentical);
        identicalNodeToGroup[potentialIdentical] = currNode;
        if (identicalGroups.find(potentialIdentical) != identicalGroups.end()) {
          for (const auto& otherIdentical : identicalGroups[potentialIdentical]) {
            identicalGroups[currNode].push_back(otherIdentical);
            identicalNodeToGroup[otherIdentical] = currNode;
          }
          identicalGroups.erase(potentialIdentical);
        }
      }
    }

    for (const auto& identical : identicals) {
      seenNodes.insert(identical);
      unseenNodes.erase(identical);
    }
  }
}


static void insertionSortAffectedSeedmerInfos(std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos) {
  auto first = affectedSeedmerInfos.begin();
  auto last = affectedSeedmerInfos.end();
  for (auto it = first + 1; it != last; ++it) {
    auto key = std::move(*it);
    auto hole = it;
    
    while (hole != first && (hole - 1)->index > key.index) {
      *hole = std::move(*(hole - 1));
      --hole;
    }
    *hole = std::move(key);
  }
}

static void sortAffectedSeedmerInfos(std::vector<mgsr::affectedSeedmerInfo>& affectedSeedmerInfos) {
  if (affectedSeedmerInfos.size() <= 20) {
    insertionSortAffectedSeedmerInfos(affectedSeedmerInfos);
  } else {
    std::sort(affectedSeedmerInfos.begin(), affectedSeedmerInfos.end(), [](const auto& a, const auto& b) {
      return a.index < b.index;
    });
  }
}

void mgsr::mgsrPlacer::setProgressTracker(ProgressTracker* tracker, size_t tid) {
  progressTracker = tracker;
  threadId = tid;
}

void mgsr::mgsrPlacer::computeOverlapCoefficientsHelper(
  MgsrLiteNode* node,
  const absl::flat_hash_set<size_t>& allSeedmerHashesSet,
  std::vector<std::pair<std::string, double>>& overlapCoefficients,
  std::vector<std::optional<size_t>>& kminmerOnRef,
  std::unordered_map<size_t, size_t>& kminmerOnRefCount
) {
  if (curDfsIndex % 1000 == 0) {
    std::cout << "\r" << curDfsIndex << " / " << liteTree->allLiteNodes.size() << std::flush;
  }

  // **** Update seeds ****
  size_t binaryOverlapKminmerCountBacktrack = binaryOverlapKminmerCount;

  const auto& seedInfos = liteTree->seedInfos;
  const auto& curSeedDeltas = node->seedDeltas;
  for (size_t i = 0; i < curSeedDeltas.size(); i++) {
    const auto [seedIndex, toDelete] = curSeedDeltas[i];
    if (i != curSeedDeltas.size() - 1 && seedInfos[seedIndex].startPos == seedInfos[curSeedDeltas[i + 1].first].startPos) {
      // substitution
      const auto& newKminmer = seedInfos[curSeedDeltas[i + 1].first];
      const auto startPos = newKminmer.startPos;
      auto& kminmerOnRefState = kminmerOnRef[startPos];
      auto oldKminmer = kminmerOnRefState.value();
      auto oldKminmerCountIt = kminmerOnRefCount.find(oldKminmer);
      oldKminmerCountIt->second -= 1;
      if (oldKminmerCountIt->second == 0) {
        kminmerOnRefCount.erase(oldKminmerCountIt);
        if (allSeedmerHashesSet.find(oldKminmer) != allSeedmerHashesSet.end()) --binaryOverlapKminmerCount;
      }
      auto [insertedKminmerCountIt, newKminmerCountInserted] = kminmerOnRefCount.emplace(newKminmer.hash, 1);
      if (!newKminmerCountInserted) {
        insertedKminmerCountIt->second++;
      } else {
        if (allSeedmerHashesSet.find(newKminmer.hash) != allSeedmerHashesSet.end()) ++binaryOverlapKminmerCount;
      }
      kminmerOnRefState = newKminmer.hash;
      i++;
    } else {
      if (toDelete) {
        // deletion
        uint32_t deletePos = seedInfos[seedIndex].startPos;
        auto oldKminmer = kminmerOnRef[deletePos].value();
        auto oldKminmerCountIt = kminmerOnRefCount.find(oldKminmer);
        oldKminmerCountIt->second -= 1;
        if (oldKminmerCountIt->second == 0) {
          kminmerOnRefCount.erase(oldKminmerCountIt);
          if (allSeedmerHashesSet.find(oldKminmer) != allSeedmerHashesSet.end()) --binaryOverlapKminmerCount;
        }
        kminmerOnRef[deletePos] = std::nullopt;
      } else {
        // insertion
        const auto& newKminmer = seedInfos[seedIndex];
        const auto startPos = newKminmer.startPos;
        auto& kminmerOnRefState = kminmerOnRef[startPos];
        auto [insertedKminmerCountIt, newKminmerCountInserted] = kminmerOnRefCount.emplace(newKminmer.hash, 1);
        if (!newKminmerCountInserted) {
          insertedKminmerCountIt->second++;
        } else {
          if (allSeedmerHashesSet.find(newKminmer.hash) != allSeedmerHashesSet.end()) ++binaryOverlapKminmerCount;
        }
        kminmerOnRefState = newKminmer.hash;
      }
    }
  }

  overlapCoefficients[curDfsIndex].first = node->identifier;
  overlapCoefficients[curDfsIndex].second = static_cast<double>(binaryOverlapKminmerCount) / static_cast<double>(kminmerOnRefCount.size());

  for (MgsrLiteNode *child : node->children) {
    ++curDfsIndex;
    computeOverlapCoefficientsHelper(child, allSeedmerHashesSet, overlapCoefficients, kminmerOnRef, kminmerOnRefCount);
  }

  // Backtrack seeds
  for (size_t i = 0; i < curSeedDeltas.size(); i++) {
    const auto [seedIndex, toDelete] = curSeedDeltas[i];
    if (i != curSeedDeltas.size() - 1 && seedInfos[seedIndex].startPos == seedInfos[curSeedDeltas[i + 1].first].startPos) {
      // substitution
      const auto& newKminmer = seedInfos[seedIndex];
      const auto startPos = newKminmer.startPos;
      auto& kminmerOnRefState = kminmerOnRef[startPos];
      auto oldKminmer = kminmerOnRefState.value();
      auto oldKminmerCountIt = kminmerOnRefCount.find(oldKminmer);
      oldKminmerCountIt->second -= 1;
      if (oldKminmerCountIt->second == 0) {
        kminmerOnRefCount.erase(oldKminmerCountIt);
        if (allSeedmerHashesSet.find(oldKminmer) != allSeedmerHashesSet.end()) --binaryOverlapKminmerCount;
      }
      auto [insertedKminmerCountIt, newKminmerCountInserted] = kminmerOnRefCount.emplace(newKminmer.hash, 1);
      if (!newKminmerCountInserted) {
        insertedKminmerCountIt->second++;
      } else {
        if (allSeedmerHashesSet.find(newKminmer.hash) != allSeedmerHashesSet.end()) ++binaryOverlapKminmerCount;
      }
      kminmerOnRefState = newKminmer.hash;
      i++;
    } else {
      if (toDelete) {
        // insert deleted kminmer
        const auto& newKminmer = seedInfos[seedIndex];
        const auto startPos = newKminmer.startPos;
        auto& kminmerOnRefState = kminmerOnRef[startPos];
        auto [insertedKminmerCountIt, newKminmerCountInserted] = kminmerOnRefCount.emplace(newKminmer.hash, 1);
        if (!newKminmerCountInserted) {
          insertedKminmerCountIt->second++;
        } else {
          if (allSeedmerHashesSet.find(newKminmer.hash) != allSeedmerHashesSet.end()) ++binaryOverlapKminmerCount;
        }
        kminmerOnRefState = newKminmer.hash;

      } else {
        // delete inserted kminmer
        uint32_t deletePos = seedInfos[seedIndex].startPos;
        auto oldKminmer = kminmerOnRef[deletePos].value();
        auto oldKminmerCountIt = kminmerOnRefCount.find(oldKminmer);
        oldKminmerCountIt->second -= 1;
        if (oldKminmerCountIt->second == 0) {
          kminmerOnRefCount.erase(oldKminmerCountIt);
          if (allSeedmerHashesSet.find(oldKminmer) != allSeedmerHashesSet.end()) --binaryOverlapKminmerCount;
        }
        kminmerOnRef[deletePos] = std::nullopt;
      }
    }
  }

  binaryOverlapKminmerCount = binaryOverlapKminmerCountBacktrack;
}

std::vector<std::pair<std::string, double>> mgsr::mgsrPlacer::computeOverlapCoefficients(const absl::flat_hash_set<size_t>& allSeedmerHashesSet) {
  curDfsIndex = 0;
  size_t maxIndex = 0;
  const auto& seedInfos = liteTree->seedInfos;
  for (size_t i = 0; i < seedInfos.size(); ++i) {
    if (seedInfos[i].startPos > maxIndex) {
      maxIndex = seedInfos[i].startPos;
    }
  }
  std::vector<std::optional<size_t>> kminmerOnRef(maxIndex + 1, std::nullopt);
  std::unordered_map<size_t, size_t> kminmerOnRefCount;
  std::vector<std::pair<std::string, double>> overlapCoefficients(liteTree->allLiteNodes.size(), std::make_pair("", 0.0));
  computeOverlapCoefficientsHelper(liteTree->root, allSeedmerHashesSet, overlapCoefficients, kminmerOnRef, kminmerOnRefCount);
  std::sort(overlapCoefficients.begin(), overlapCoefficients.end(), [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  size_t overlapCoefficientCutoff = 1000;
  size_t curRank = 0;
  size_t curIndex = 0;
  double curCoefficient = overlapCoefficients[0].second;
  for (size_t i = 1; i < overlapCoefficients.size(); ++i) {
    if (overlapCoefficients[i].second != curCoefficient) {
      ++curRank;
    }
    if (curRank >= overlapCoefficientCutoff) break;
    curIndex = i;
  }

  return overlapCoefficients;

  // std::cout << "overlapCoefficientCutoff: " << overlapCoefficientCutoff << " last index to include: " << curIndex << std::endl;
  // exit(0);
}

void mgsr::mgsrPlacer::traverseTreeHelper(MgsrLiteNode* node) {
  // Update progress if tracker is available
  if (progressTracker) {
    progressTracker->incrementProgress(threadId);
  }

  // **** Update seeds ****
  std::unordered_set<uint64_t> affectedSeedmers;
  updateSeeds(node, affectedSeedmers);

  // **** Update gapMap ****
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapBacktracks;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapBlocksBacktracks;
  updateGapMap(node, gapMapBacktracks, gapMapBlocksBacktracks); 

  // **** Revert gapMap inversions ****
  revertGapMapInversions(gapMapBlocksBacktracks, gapMap);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapMapBlocksBacktracks); // gapMapBlocksBacktracks is no longer needed... clear memory

  auto nodeDfsIndex = curDfsIndex;
  for (MgsrLiteNode *child : node->children) {
    ++curDfsIndex;
    traverseTreeHelper(child);
  }

  // Backtrack seeds and delayedRefSeedmerStatus
  backtrackSeeds(node, nodeDfsIndex);
  
  // Backtrack gapMap
  for (auto it = gapMapBacktracks.rbegin(); it != gapMapBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }
}

void mgsr::mgsrPlacer::traverseTree() {
  gapMap.clear();
  gapMap.insert(std::make_pair(0, liteTree->blockScalarRanges.back().second));

  curDfsIndex = 0;
  traverseTreeHelper(liteTree->root);
}


void mgsr::mgsrPlacer::placeReadsHelper(MgsrLiteNode* node) {
  // Update progress if tracker is available
  if (progressTracker) {
    progressTracker->incrementProgress(threadId);
  }

  // **** Update seeds ****
  std::unordered_set<uint64_t> affectedSeedmers;
  updateSeeds(node, affectedSeedmers);

  // **** Update gapMap ****
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapBacktracks;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapBlocksBacktracks;
  updateGapMap(node, gapMapBacktracks, gapMapBlocksBacktracks);


  // **** Start placing reads ****
  // unordered_map<readIndex, pair<vector<affectedSeedmerIndexOnRead>, pair<allUniqueToNonUnique, allNonUniqueToUnique>>>
  absl::flat_hash_map<uint32_t, std::pair<std::vector<mgsr::affectedSeedmerInfo>, std::pair<bool, bool>>> readToAffectedSeedmerIndex;
  size_t binaryOverlapKminmerCountBacktract = binaryOverlapKminmerCount;
  fillReadToAffectedSeedmerIndex(readToAffectedSeedmerIndex, affectedSeedmers);
  if (threadId == 0) {
    kminmerOverlapCoefficients[node->identifier] = static_cast<double>(binaryOverlapKminmerCount) / static_cast<double>(hashToPositionMap.size());
  }

  std::vector<std::pair<size_t, int32_t>> readScoresBacktrack;
  std::vector<std::pair<size_t, std::vector<mgsr::Minichain>>> readMinichainsBacktrack;

  std::vector<readScoreDelta> currentNodeScoreDeltas;
  if (affectedSeedmers.empty()) {
    // If the node is identical to its parent, add it to the identical group
    if (node->parent != nullptr && threadId == 0) {
      const std::string& parentID = node->parent->identifier;
      auto identicalGroupIt = identicalNodeToGroup.find(parentID);
      if (identicalGroupIt == identicalNodeToGroup.end()) {
        identicalNodeToGroup[node->identifier] = parentID;
        identicalGroups[parentID] = {parentID, node->identifier};
      } else {
        identicalNodeToGroup[node->identifier] = identicalGroupIt->second;
        identicalGroups[identicalGroupIt->second].push_back(node->identifier);
      }
    }
  } else if (positionMap.empty()){
    // node for some reason has no seeds... delete all read scores and minichains
    // preallocate memory for score deltas
    size_t numReadsToReset = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
      if (reads[i].readType != mgsr::ReadType::PASS) continue;
      numReadsToReset += readScores[i] != 0;
    }
    readScoresBacktrack.reserve(numReadsToReset);
    readMinichainsBacktrack.reserve(numReadsToReset);
    currentNodeScoreDeltas.reserve(numReadsToReset);

    for (size_t i = 0; i < reads.size(); ++i) {
      if (reads[i].readType != mgsr::ReadType::PASS) continue;
      if (readScores[i] != 0) {
        readScoresBacktrack.emplace_back(i, readScores[i]);
        if (!lowMemory) {
          currentNodeScoreDeltas.emplace_back(i, 0);
        } else {
          currentNodeScoreDeltas.emplace_back(i, -readScores[i]);
        }
        setReadScore(i, 0);
        readMinichainsBacktrack.emplace_back(i, reads[i].minichains);
        reads[i].minichains.clear();
      }
    }
  } else {
    // Update read scores and minichains
    readScoresBacktrack.resize(readToAffectedSeedmerIndex.size());
    readMinichainsBacktrack.resize(readToAffectedSeedmerIndex.size());
    currentNodeScoreDeltas.resize(readToAffectedSeedmerIndex.size());

    size_t backtrackIndex = 0;
    for (auto& [readIndex, affectedSeedmerIndexInfo] : readToAffectedSeedmerIndex) {
      /*chaining debug*/if (false) std::cout << "[" << node->identifier << " " << curDfsIndex << "] processingreadIndex: " << readIndex << std::endl;

      mgsr::Read& curRead = reads[readIndex];
      if (curRead.readType != mgsr::ReadType::PASS) continue;

      auto& affectedSeedmerInfos = affectedSeedmerIndexInfo.first;
      const auto& [allUniqueToNonUnique, allNonUniqueToUnique] = affectedSeedmerIndexInfo.second;
      if (affectedSeedmerInfos.empty()) {
        std::cerr << "Error: affectedSeedmerIndices is empty" << std::endl;
        exit(1);
      }

      sortAffectedSeedmerInfos(affectedSeedmerInfos);

      readMinichainsBacktrack[backtrackIndex].first = readIndex;
      readMinichainsBacktrack[backtrackIndex].second = curRead.minichains;

      if (allUniqueToNonUnique || allNonUniqueToUnique) {
        /*chaining debug*/if (false) {
          std::cout << "\tread is allUniqueToNonUnique " << allUniqueToNonUnique << " or allNonUniqueToUnique " << allNonUniqueToUnique << "... updating minichains" << std::endl;
          for (const auto& minichain : curRead.minichains) {
            std::cout << "\t\t(" << minichain.begIndex << ", " << minichain.endIndex << ", " << minichain.rev << ") ";
          }
          std::cout << std::endl;
        }
        updateMinichains(readIndex, affectedSeedmerInfos, allUniqueToNonUnique, allNonUniqueToUnique);
      } else if (affectedSeedmerInfos.size() < curRead.seedmersList.size() / 2) {
        updateMinichainsMixed(readIndex, affectedSeedmerInfos);
      } else {
        /*chaining debug*/if (false) std::cout << "\tread is not allUniqueToNonUnique or allNonUniqueToUnique... reinitializing minichains" << std::endl;
        initializeReadMinichains(curRead);
      }
      int32_t pseudoScore = getReadPseudoScore(curRead);

      readScoresBacktrack[backtrackIndex].first         = readIndex;
      readScoresBacktrack[backtrackIndex].second        = readScores[readIndex];
      currentNodeScoreDeltas[backtrackIndex].readIndex  = readIndex;
      if (lowMemory) {
        currentNodeScoreDeltas[backtrackIndex].scoreDelta = pseudoScore - readScores[readIndex];
      } else {
        currentNodeScoreDeltas[backtrackIndex].scoreDelta = pseudoScore;
      }
      setReadScore(readIndex, pseudoScore);



      ++backtrackIndex;

      // int64_t bruteForceScore = getReadBruteForceScore(readIndex, hashCoordInfoCacheTable);
      // if (pseudoScore != bruteForceScore) {
      //   std::cerr << "[" << node->identifier << " " << curDfsIndex << "] readIndex: " << readIndex << " dynamic: " << pseudoScore << " != brute force: " << bruteForceScore << std::endl;
      //   exit(1);
      // }

      // if (curRead.duplicates.size() > excludeDuplicatesThreshold * curRead.seedmersList.size()) readTypes[readIndex] = mgsr::readType::HIGH_DUPLICATES;
    }
  }
  


  if (lowMemory && currentNodeScoreDeltas.size() > 0) {
    auto& currentNodeScoreDeltasGrouped = node->readScoreDeltasLowMemory[threadId];
    currentNodeScoreDeltasGrouped.reserve(currentNodeScoreDeltas.size());
    std::sort(currentNodeScoreDeltas.begin(), currentNodeScoreDeltas.end(), [&](const auto& a, const auto& b) {
      return a.readIndex < b.readIndex;
    });
    currentNodeScoreDeltasGrouped.emplace_back(mgsr::readScoreDeltaLowMemory{
      .readIndex = currentNodeScoreDeltas[0].readIndex,
      .scoreDelta = currentNodeScoreDeltas[0].scoreDelta,
    });
    auto currentGroup = &currentNodeScoreDeltasGrouped.back();
    for (size_t i = 1; i < currentNodeScoreDeltas.size(); ++i) {
      auto [currentScoreDeltaIndex, currentScoreDelta] = currentNodeScoreDeltas[i];
      bool startNew = false;

      int16_t scoreDeltaDiff = currentNodeScoreDeltas[i].scoreDelta - currentGroup->scoreDelta;

      if (currentScoreDeltaIndex > currentGroup->readIndex + 16) {
        startNew = true;
      } else if (currentNodeScoreDeltas[i].scoreDelta > currentGroup->scoreDelta && scoreDeltaDiff > 7) {
        startNew = true;
      } else if (currentNodeScoreDeltas[i].scoreDelta < currentGroup->scoreDelta && scoreDeltaDiff < -7) {
        startNew = true;
      }

      if (startNew) {
        currentNodeScoreDeltasGrouped.emplace_back(mgsr::readScoreDeltaLowMemory{
          .readIndex = currentScoreDeltaIndex,
          .scoreDelta = currentScoreDelta,
        });
        currentGroup = &currentNodeScoreDeltasGrouped.back();
      } else {
        currentGroup->encodeTrailingDelta(scoreDeltaDiff, currentScoreDeltaIndex);
      }
    }
    currentNodeScoreDeltasGrouped.shrink_to_fit();
    numGroupsUpdate += currentNodeScoreDeltasGrouped.size();
    numReadsUpdate += currentNodeScoreDeltas.size();
    std::vector<readScoreDelta>().swap(currentNodeScoreDeltas);
  } else if (!lowMemory) {
    node->readScoreDeltas[threadId] = std::move(currentNodeScoreDeltas);
  }




  // **** Revert gapMap inversions ****
  revertGapMapInversions(gapMapBlocksBacktracks, gapMap);

  // clear memory that are no longer needed
  decltype(gapMapBlocksBacktracks){}.swap(gapMapBlocksBacktracks);
  decltype(readToAffectedSeedmerIndex){}.swap(readToAffectedSeedmerIndex);

  auto nodeDfsIndex = curDfsIndex;
  for (MgsrLiteNode *child : node->children) {
    ++curDfsIndex;
    placeReadsHelper(child);
  }

  // Backtrack seeds and delayedRefSeedmerStatus
  backtrackSeeds(node, nodeDfsIndex);


  // Backtrack gapMap
  for (auto it = gapMapBacktracks.rbegin(); it != gapMapBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  // Backtrack kminmer count and overlap coefficient
  binaryOverlapKminmerCount = binaryOverlapKminmerCountBacktract;

  // Backtrack read scores
  for (size_t i = 0; i < readScoresBacktrack.size(); ++i) {
    const auto& [readIdx, score] = readScoresBacktrack[i];
    setReadScore(readIdx, score);
  }

  // Backtrack read minichains
  for (const auto& [readIdx, minichains] : readMinichainsBacktrack) {
    reads[readIdx].minichains = std::move(minichains);
  }

}

void mgsr::mgsrPlacer::placeReads() {
  gapMap.clear();
  gapMap.insert(std::make_pair(0, liteTree->blockScalarRanges.back().second));

  curDfsIndex = 0;
  preallocateHashCoordInfoCacheTable(0, reads.size());
  placeReadsHelper(liteTree->root);
}

void mgsr::mgsrPlacer::scoreReadsHelper(mgsr::MgsrLiteNode* node) {
  // Update progress if tracker is available
  if (progressTracker) {
    progressTracker->incrementProgress(threadId);
  }

  std::unordered_map<uint32_t, mgsr::ModifiedReadInfo> modifiedReads;
  
  // UPDATE SEEDS AND READS AT THE SAME TIME
  const auto& curSeedDeltas = node->seedDeltas;
  const auto& seedInfos = liteTree->seedInfos;
  for (const auto [seedIndex, toDelete] : curSeedDeltas) {
    if (seedmerToReads.find(seedInfos[seedIndex].hash) == seedmerToReads.end()) {
      continue;
    }
    if (toDelete) {
      size_t deletedSeedHash = seedInfos[seedIndex].hash;
      size_t deletedSeedRev  = seedInfos[seedIndex].isReverse;
      auto deletedKminmerCountIt = kminmerOnRefCount.find(deletedSeedHash);
      auto& [forwardCount, reverseCount] = deletedKminmerCountIt->second;

      if (deletedSeedRev) {
        if (reverseCount == 1) {
          if (forwardCount == 0) {
            kminmerOnRefCount.erase(deletedKminmerCountIt);
          } else {
            reverseCount--;
          }
          auto seedmerToReadsIt = seedmerToReads.find(deletedSeedHash);
          if (seedmerToReadsIt != seedmerToReads.end()) {
            for (const auto& [readIndex, seedmerIndexOnRead] : seedmerToReadsIt->second) {
              auto& read = reads[readIndex];
              if (read.readType != mgsr::ReadType::PASS) continue;
              modifiedReads.try_emplace(readIndex, mgsr::ModifiedReadInfo{read.numForwardMatching, read.numReverseMatching});
              if (read.seedmersList[seedmerIndexOnRead].rev == deletedSeedRev) {
                read.numForwardMatching--;
                if (read.numForwardMatching < 0) {
                  std::cerr << "Error! Read " << readIndex << " has negative forward matching count!" << std::endl;
                  exit(1);
                }
              } else {
                read.numReverseMatching--;
                if (read.numReverseMatching < 0) {
                  std::cerr << "Error! Read " << readIndex << " has negative reverse matching count!" << std::endl;
                  exit(1);
                }
              }
            }
          }
        } else {
          reverseCount--;
        }
      } else {
        if (forwardCount == 1) {
          if (reverseCount == 0) {
            kminmerOnRefCount.erase(deletedKminmerCountIt);
          } else {
            forwardCount--;
          }
          auto seedmerToReadsIt = seedmerToReads.find(deletedSeedHash);
          if (seedmerToReadsIt != seedmerToReads.end()) {
            for (const auto& [readIndex, seedmerIndexOnRead] : seedmerToReadsIt->second) {
              auto& read = reads[readIndex];
              if (read.readType != mgsr::ReadType::PASS) continue;
              modifiedReads.try_emplace(readIndex, mgsr::ModifiedReadInfo{read.numForwardMatching, read.numReverseMatching});
              if (read.seedmersList[seedmerIndexOnRead].rev == deletedSeedRev) {
                read.numForwardMatching--;
                if (read.numForwardMatching < 0) {
                  std::cerr << "Error! Read " << readIndex << " has negative forward matching count!" << std::endl;
                  exit(1);
                }
              } else {
                read.numReverseMatching--;
                if (read.numReverseMatching < 0) {
                  std::cerr << "Error! Read " << readIndex << " has negative reverse matching count!" << std::endl;
                  exit(1);
                }
              }
            }
          }
        } else {
          forwardCount--;
        }
      }
    } else {
      size_t insertedSeedHash = seedInfos[seedIndex].hash;
      bool insertedSeedRev = seedInfos[seedIndex].isReverse;
      auto [insertedKminmerCountIt, insertedKminmerCountInserted] = insertedSeedRev ? kminmerOnRefCount.try_emplace(insertedSeedHash, 0, 1) : kminmerOnRefCount.try_emplace(insertedSeedHash, 1, 0);
      if (!insertedKminmerCountInserted) {
        if (insertedSeedRev) {
          insertedKminmerCountIt->second.second++;
        } else {
          insertedKminmerCountIt->second.first++;
        }
      }
      if (insertedSeedRev && insertedKminmerCountIt->second.second == 1 || 
         !insertedSeedRev && insertedKminmerCountIt->second.first  == 1
      ) {
        auto seedmerToReadsIt = seedmerToReads.find(insertedSeedHash);
        if (seedmerToReadsIt != seedmerToReads.end()) {
          for (const auto& [readIndex, seedmerIndexOnRead] : seedmerToReadsIt->second) {
            auto& read = reads[readIndex];
            if (read.readType != mgsr::ReadType::PASS) continue;
            modifiedReads.try_emplace(readIndex, mgsr::ModifiedReadInfo{read.numForwardMatching, read.numReverseMatching});
            if (read.seedmersList[seedmerIndexOnRead].rev == seedInfos[seedIndex].isReverse) {
              read.numForwardMatching++;
              if (read.numForwardMatching > static_cast<int32_t>(read.seedmersList.size())) {
                std::cerr << "Error! Read " << readIndex << " has forward matching count greater than total seedmers!" << std::endl;
                exit(1);
              }
            } else {
              read.numReverseMatching++;
              if (read.numReverseMatching > static_cast<int32_t>(read.seedmersList.size())) {
                std::cerr << "Error! Read " << readIndex << " has reverse matching count greater than total seedmers!" << std::endl;
                exit(1);
              }
            }
          }
        }
      }
    }
  }



  // STORE READ SCORE INDEX
  if (modifiedReads.empty()) {
    if (node->parent != nullptr && threadId == 0) {
      const std::string& parentID = node->parent->identifier;
      auto identicalGroupIt = identicalNodeToGroup.find(parentID);
      if (identicalGroupIt == identicalNodeToGroup.end()) {
        identicalNodeToGroup[node->identifier] = parentID;
        identicalGroups[parentID] = {parentID, node->identifier};
      } else {
        identicalNodeToGroup[node->identifier] = identicalGroupIt->second;
        identicalGroups[identicalGroupIt->second].push_back(node->identifier);
      }
    }
  } else {
    if (lowMemory) {
      std::vector<uint32_t> modifiedReadIndices;
      modifiedReadIndices.reserve(modifiedReads.size());
      for (const auto& [readIndex, modifiedReadInfo] : modifiedReads) {
        modifiedReadIndices.push_back(readIndex);
      }
      std::sort(modifiedReadIndices.begin(), modifiedReadIndices.end());

      auto& currentNodeScoreDeltasGrouped = node->readScoreDeltasLowMemory[threadId];
      currentNodeScoreDeltasGrouped.reserve(modifiedReads.size());
      
      const auto readIndex = modifiedReadIndices[0];
      int32_t curScore = std::max(reads[readIndex].numForwardMatching, reads[readIndex].numReverseMatching);
      reads[readIndex].maxScore = std::max(reads[readIndex].maxScore, curScore);
      currentNodeScoreDeltasGrouped.emplace_back(mgsr::readScoreDeltaLowMemory{
        .readIndex = readIndex,
        .scoreDelta = curScore
      });
      auto currentGroup = &currentNodeScoreDeltasGrouped.back();
      for (size_t i = 1; i < modifiedReadIndices.size(); ++i) {
        const auto readIndex = modifiedReadIndices[i];
        int32_t curScore = std::max(reads[readIndex].numForwardMatching, reads[readIndex].numReverseMatching);
        reads[readIndex].maxScore = std::max(reads[readIndex].maxScore, curScore);
        int16_t scoreDeltaDiff = curScore - currentGroup->scoreDelta;

        bool startNew = false;
        if (readIndex > currentGroup->readIndex + 16) {
          startNew = true;
        } else if (curScore > currentGroup->scoreDelta && scoreDeltaDiff > 7) {
          startNew = true;
        } else if (curScore < currentGroup->scoreDelta && scoreDeltaDiff < -7) {
          startNew = true;
        }

        if (startNew) {
          currentNodeScoreDeltasGrouped.emplace_back(mgsr::readScoreDeltaLowMemory{
            .readIndex = readIndex,
            .scoreDelta = curScore
          });
          currentGroup = &currentNodeScoreDeltasGrouped.back();
        } else {
          currentGroup->encodeTrailingDelta(scoreDeltaDiff, readIndex);
        }
      }
      currentNodeScoreDeltasGrouped.shrink_to_fit();
    } else {
      auto& currentNodeScoreDeltas = node->readScoreDeltas[threadId];
      currentNodeScoreDeltas.reserve(modifiedReads.size());
      for (const auto& [readIndex, _] : modifiedReads) {
        int32_t score = std::max(reads[readIndex].numForwardMatching, reads[readIndex].numReverseMatching);
        reads[readIndex].maxScore = std::max(reads[readIndex].maxScore, score);
        currentNodeScoreDeltas.emplace_back(readIndex, score);
      }
    }
  }
  
  auto nodeDfsIndex = node->dfsIndex;
  for (MgsrLiteNode *child : node->children) {
    ++curDfsIndex;
    scoreReadsHelper(child);
  }

  // BACKTRACK
  for (const auto [seedIndex, toDelete] : curSeedDeltas) {
    if (seedmerToReads.find(seedInfos[seedIndex].hash) == seedmerToReads.end()) {
      continue;
    }
    if (!toDelete) {
      size_t deletedSeedHash = seedInfos[seedIndex].hash;
      bool   deletedSeedRev  = seedInfos[seedIndex].isReverse;
      auto deletedKminmerCountIt = kminmerOnRefCount.find(deletedSeedHash);
      auto& [forwardCount, reverseCount] = deletedKminmerCountIt->second;
      if (deletedSeedRev) {
        if (reverseCount == 1) {
          if (forwardCount == 0) {
            kminmerOnRefCount.erase(deletedKminmerCountIt);
          } else {
            reverseCount--;
          }
        } else {
          reverseCount--;
        }
      } else {
        if (forwardCount == 1) {
          if (reverseCount == 0) {
            kminmerOnRefCount.erase(deletedKminmerCountIt);
          } else {
            forwardCount--;
          }
        } else {
          forwardCount--;
        }
      }
    } else {
      size_t insertedSeedHash = seedInfos[seedIndex].hash;
      bool   insertedSeedRev  = seedInfos[seedIndex].isReverse;
      auto [insertedKminmerCountIt, insertedKminmerCountInserted] = insertedSeedRev ? kminmerOnRefCount.try_emplace(insertedSeedHash, 0, 1) : kminmerOnRefCount.try_emplace(insertedSeedHash, 1, 0);
      if (!insertedKminmerCountInserted) {
        if (insertedSeedRev) {
          insertedKminmerCountIt->second.second++;
        } else {
          insertedKminmerCountIt->second.first++;
        }
      }
    }
  }

  for (const auto& [readIndex, modifiedReadInfo] : modifiedReads) {
    reads[readIndex].numForwardMatching = modifiedReadInfo.forwardOriginalScore;
    reads[readIndex].numReverseMatching = modifiedReadInfo.reverseOriginalScore;
  }

}

void mgsr::mgsrPlacer::scoreReads() {
  scoreReadsHelper(liteTree->root);
}


void mgsr::mgsrPlacer::preallocateHashCoordInfoCacheTable(uint32_t startReadIndex, uint32_t endReadIndex) {
  std::unordered_set<size_t> uniqueHashes;
  for (const auto& read : reads) {
    for (const auto& seedmer : read.seedmersList) {
      uniqueHashes.insert(seedmer.hash);
    }
  }

  hashCoordInfoCacheTable.reserve(uniqueHashes.size());
  for (const auto& hash : uniqueHashes) {
    hashCoordInfoCacheTable.try_emplace(hash);
  }
}

void mgsr::ThreadsManager::computeKminmerCoverageHelper(
  MgsrLiteNode* node,
  std::vector<uint32_t>& readScores,
  const absl::flat_hash_map<size_t, std::vector<std::pair<uint32_t, uint32_t>>>& seedmerToReads,
  absl::flat_hash_map<size_t, std::unordered_set<uint32_t>>& coveredKminmers,
  absl::flat_hash_map<size_t, uint32_t>& refKminmers,
  size_t& dfsIndex
) {
  size_t readScoresBacktrackIndex = 0;
  std::vector<std::pair<uint32_t, int16_t>> readScoresBacktrack;

  // Apply read score deltas to readScores
  for (size_t threadId = 0; threadId < numThreads; ++threadId) {
    const auto threadReadStart = threadRanges[threadId].first;
    if (lowMemory) {
      const auto& curNodeThreadScoreDeltasLowMemory = node->readScoreDeltasLowMemory[threadId];
      for (const auto& scoreDelta : curNodeThreadScoreDeltasLowMemory) {
        const auto [trailingDelta, readIndex, numTrailing, currentScoreDelta] = scoreDelta;
        if (readScoresBacktrackIndex == readScoresBacktrack.size()) {
          readScoresBacktrack.emplace_back(threadReadStart + readIndex, readScores[threadReadStart + readIndex]);
        } else {
          readScoresBacktrack[readScoresBacktrackIndex].first = threadReadStart + readIndex;
          readScoresBacktrack[readScoresBacktrackIndex].second = readScores[threadReadStart + readIndex];
        }
        readScoresBacktrackIndex++;

        readScores[threadReadStart + readIndex] += currentScoreDelta;
        for (size_t i = 0; i < numTrailing; ++i) {
          int16_t curDecodedTrailingDelta = scoreDelta.decodeTrailingDelta(i);
          // if curDecodedTrailingDelta == -8, it means the read is not changed, so we don't need to apply the trailing delta
          if (curDecodedTrailingDelta > -8) {
            if (readScoresBacktrackIndex == readScoresBacktrack.size()) {
              readScoresBacktrack.emplace_back(threadReadStart + readIndex + i + 1, readScores[threadReadStart + readIndex + i + 1]);
            } else {
              readScoresBacktrack[readScoresBacktrackIndex].first = threadReadStart + readIndex + i + 1;
              readScoresBacktrack[readScoresBacktrackIndex].second = readScores[threadReadStart + readIndex + i + 1];
            }
            readScoresBacktrackIndex++;

            readScores[threadReadStart + readIndex + i + 1] += currentScoreDelta + curDecodedTrailingDelta;
          }
        }
      }
    } else {
      const auto& curNodeThreadScoreDeltas = node->readScoreDeltas[threadId];
      for (const auto& scoreDelta : curNodeThreadScoreDeltas) {
        if (readScoresBacktrackIndex == readScoresBacktrack.size()) {
          readScoresBacktrack.emplace_back(threadReadStart + scoreDelta.readIndex, readScores[threadReadStart + scoreDelta.readIndex]);
        } else {
          readScoresBacktrack[readScoresBacktrackIndex].first = threadReadStart + scoreDelta.readIndex;
          readScoresBacktrack[readScoresBacktrackIndex].second = readScores[threadReadStart + scoreDelta.readIndex];
        }
        readScoresBacktrackIndex++;

        readScores[threadReadStart + scoreDelta.readIndex] = scoreDelta.scoreDelta;
      }
    }
  }

  size_t coveredKminmersBacktrackIndex = 0;
  std::vector<mgsr::kminmerCoverageBacktrack> coveredKminmersBacktrack;
  const auto& curSeedDeltas = node->seedDeltas;
  const auto& seedInfos = liteTree->seedInfos;
  for (size_t i = 0; i < curSeedDeltas.size(); i++) {
    const auto [seedIndex, toDelete] = curSeedDeltas[i];
    if (i != curSeedDeltas.size() - 1 && seedInfos[seedIndex].startPos == seedInfos[curSeedDeltas[i + 1].first].startPos) {
      // substitution
      const auto& delSeed = seedInfos[seedIndex];
      auto delRefKminmerIt = refKminmers.find(delSeed.hash);
      if (delRefKminmerIt == refKminmers.end()) {
        std::cerr << "Error: seed " << seedIndex << " not found in refKminmers" << std::endl;
        exit(1);
      }
      delRefKminmerIt->second--;
      if (delRefKminmerIt->second == 0) {
        refKminmers.erase(delRefKminmerIt);
        if (coveredKminmers.find(delSeed.hash) != coveredKminmers.end()) {
          for (const auto& readIndex : coveredKminmers[delSeed.hash]) {
            if (coveredKminmersBacktrackIndex == coveredKminmersBacktrack.size()) {
              coveredKminmersBacktrack.emplace_back(delSeed.hash, readIndex, false);
            } else {
              auto& coveredKminmerBacktrack = coveredKminmersBacktrack[coveredKminmersBacktrackIndex];
              coveredKminmerBacktrack.seedmer = delSeed.hash;
              coveredKminmerBacktrack.readIndex = readIndex;
              coveredKminmerBacktrack.toDelete = false;
            }
            coveredKminmersBacktrackIndex++;
          }
        }
        coveredKminmers.erase(delSeed.hash);
      }
  
      const auto& newSeed = seedInfos[curSeedDeltas[i + 1].first];
      auto [newRefKminmerIt, newRefKminmerInserted] = refKminmers.try_emplace(newSeed.hash, 1);
      if (!newRefKminmerInserted) newRefKminmerIt->second++;
      i++;
    } else {
      if (toDelete) {
        // deletion
        const auto& seed = seedInfos[seedIndex];
        auto refKminmerIt = refKminmers.find(seed.hash);
        if (refKminmerIt == refKminmers.end()) {
          std::cerr << "Error: seed " << seedIndex << " not found in refKminmers" << std::endl;
          exit(1);
        }
        refKminmerIt->second--;
        if (refKminmerIt->second == 0) {
          refKminmers.erase(refKminmerIt);
          if (coveredKminmers.find(seed.hash) != coveredKminmers.end()) {
            for (const auto& readIndex : coveredKminmers[seed.hash]) {
              if (coveredKminmersBacktrackIndex == coveredKminmersBacktrack.size()) {
                coveredKminmersBacktrack.emplace_back(seed.hash, readIndex, false);
              } else {
                auto& coveredKminmerBacktrack = coveredKminmersBacktrack[coveredKminmersBacktrackIndex];
                coveredKminmerBacktrack.seedmer = seed.hash;
                coveredKminmerBacktrack.readIndex = readIndex;
                coveredKminmerBacktrack.toDelete = false;
              }
              coveredKminmersBacktrackIndex++;
            }
          }
          coveredKminmers.erase(seed.hash);
        }
      } else {
        // insertion
        const auto& seed = seedInfos[seedIndex];
        auto [refKminmerIt, refKminmerInserted] = refKminmers.try_emplace(seed.hash, 1);
        if (!refKminmerInserted) refKminmerIt->second++;
      }
    }
  }

  for (size_t i = 0; i < readScoresBacktrackIndex; ++i) {
    const auto [readIndex, oldScore] = readScoresBacktrack[i];
    const auto newScore = readScores[readIndex];
    const auto maxScore = reads[readIndex].maxScore;
    if (oldScore != maxScore && newScore == maxScore) {
      for (const auto& [seedmer, _] : reads[readIndex].uniqueSeedmers) {
        if (refKminmers.find(seedmer) == refKminmers.end()) continue;
        auto [coveredKminmerIt, coveredKminmerInserted] = coveredKminmers.try_emplace(seedmer, std::unordered_set<uint32_t>({readIndex}));
        if (!coveredKminmerInserted) coveredKminmerIt->second.insert(readIndex);
        if (coveredKminmersBacktrackIndex == coveredKminmersBacktrack.size()) {
          coveredKminmersBacktrack.emplace_back(seedmer, readIndex, true);
        } else {
          auto& coveredKminmerBacktrack = coveredKminmersBacktrack[coveredKminmersBacktrackIndex];
          coveredKminmerBacktrack.seedmer = seedmer;
          coveredKminmerBacktrack.readIndex = readIndex;
          coveredKminmerBacktrack.toDelete = true;
        }
        coveredKminmersBacktrackIndex++;
      }
    } else if (oldScore == maxScore && newScore != maxScore) {
      for (const auto& [seedmer, _] : reads[readIndex].uniqueSeedmers) {
        auto coveredKminmerIt = coveredKminmers.find(seedmer);
        if (refKminmers.find(seedmer) == refKminmers.end() || coveredKminmerIt == coveredKminmers.end()) continue;
        coveredKminmerIt->second.erase(readIndex);
        if (coveredKminmerIt->second.empty()) coveredKminmers.erase(coveredKminmerIt);
        if (coveredKminmersBacktrackIndex == coveredKminmersBacktrack.size()) {
          coveredKminmersBacktrack.emplace_back(seedmer, readIndex, false);
        } else {
          auto& coveredKminmerBacktrack = coveredKminmersBacktrack[coveredKminmersBacktrackIndex];
          coveredKminmerBacktrack.seedmer = seedmer;
          coveredKminmerBacktrack.readIndex = readIndex;
          coveredKminmerBacktrack.toDelete = false;
        }
        coveredKminmersBacktrackIndex++;
      }
    } else if (oldScore == maxScore && newScore == maxScore) {
      for (const auto& [seedmer, _] : reads[readIndex].uniqueSeedmers) {
        if (refKminmers.find(seedmer) == refKminmers.end()) continue;
        auto [coveredKminmerIt, coveredKminmerInserted] = coveredKminmers.try_emplace(seedmer, std::unordered_set<uint32_t>({readIndex}));
        if (coveredKminmerInserted) {
          if (coveredKminmersBacktrackIndex == coveredKminmersBacktrack.size()) {
            coveredKminmersBacktrack.emplace_back(seedmer, readIndex, true);
          } else {
            auto& coveredKminmerBacktrack = coveredKminmersBacktrack[coveredKminmersBacktrackIndex];
            coveredKminmerBacktrack.seedmer = seedmer;
            coveredKminmerBacktrack.readIndex = readIndex;
            coveredKminmerBacktrack.toDelete = true;
          }
          coveredKminmersBacktrackIndex++;
        } else {
          if (coveredKminmerIt->second.find(readIndex) == coveredKminmerIt->second.end()) {
            if (coveredKminmersBacktrackIndex == coveredKminmersBacktrack.size()) {
              coveredKminmersBacktrack.emplace_back(seedmer, readIndex, true);
            } else {
              auto& coveredKminmerBacktrack = coveredKminmersBacktrack[coveredKminmersBacktrackIndex];
              coveredKminmerBacktrack.seedmer = seedmer;
              coveredKminmerBacktrack.readIndex = readIndex;
              coveredKminmerBacktrack.toDelete = true;
            }
            coveredKminmersBacktrackIndex++;
            coveredKminmerIt->second.insert(readIndex);
          } 
        }
      }
    }
  }



  // // compute coveredKminmersBruteForce
  // std::unordered_map<size_t, std::unordered_set<uint32_t>> coveredKminmersBruteForce;
  // for (size_t i = 0; i < reads.size(); ++i) {
  //   if (readScores[i] == reads[i].maxScore) {
  //     for (const auto& seedmer : reads[i].uniqueSeedmers) {
  //       if (refKminmers.find(seedmer.first) != refKminmers.end()) {
  //         coveredKminmersBruteForce[seedmer.first].insert(i);
  //       }
  //     }
  //   }
  // }

  // if (coveredKminmersBruteForce.size() != coveredKminmers.size()) {
  //   std::cout << "Covered kminmers size mismatch: brute force " << coveredKminmersBruteForce.size() << " != dynamic " << coveredKminmers.size() << std::endl;
  //   exit(1);
  // }
  // for (const auto& [seedmer, reads] : coveredKminmersBruteForce) {
  //   auto coveredKminmerIt = coveredKminmers.find(seedmer);
  //   if (coveredKminmerIt == coveredKminmers.end()) {
  //     std::cout << "Covered kminmer " << seedmer << " not found in coveredKminmers" << std::endl;
  //     exit(1);
  //   }

  //   if (reads.size() != coveredKminmerIt->second.size()) {
  //     std::cout << "Covered kminmer " << seedmer << " reads size mismatch: brute force " << reads.size() << " != dynamic " << coveredKminmerIt->second.size() << std::endl;
  //     exit(1);
  //   }
  //   for (const auto& readIndex : reads) {
  //     if (coveredKminmerIt->second.find(readIndex) == coveredKminmerIt->second.end()) {
  //       std::cout << "Covered kminmer " << seedmer << " read index " << readIndex << " not found in coveredKminmers" << std::endl;
  //       exit(1);
  //     }
  //   }
  // }


  kminmerCoverage[node->identifier] = static_cast<double>(coveredKminmers.size()) / static_cast<double>(refKminmers.size());


  if (dfsIndex % 100 == 0) {
    std::cout << "\rComputing Kminmer Coverage " << dfsIndex << " / " << liteTree->allLiteNodes.size() << std::flush;
  }
  auto curNodeDfsIndex = dfsIndex;
  for (auto child : node->children) {
    ++dfsIndex;
    computeKminmerCoverageHelper(child, readScores, seedmerToReads, coveredKminmers, refKminmers, dfsIndex);
  }

  // backtrack seeds
  for (size_t i = 0; i < curSeedDeltas.size(); i++) {
    const auto [seedIndex, toDelete] = curSeedDeltas[i];
    if (i != curSeedDeltas.size() - 1 && seedInfos[seedIndex].startPos == seedInfos[curSeedDeltas[i + 1].first].startPos) {
      // substitute
      const auto& insSeed = seedInfos[seedIndex];
      auto [insertedRefKminmerIt, insertedRefKminmerInserted] = refKminmers.try_emplace(insSeed.hash, 1);
      if (!insertedRefKminmerInserted) insertedRefKminmerIt->second++;
      
      const auto& delSeed = seedInfos[curSeedDeltas[i + 1].first];
      auto delRefKminmerIt = refKminmers.find(delSeed.hash);
      if (delRefKminmerIt == refKminmers.end()) {
        std::cerr << "Error: backtracking seed " << curSeedDeltas[i + 1].first << " not found in refKminmers" << std::endl;
        exit(1);
      }
      delRefKminmerIt->second--;
      if (delRefKminmerIt->second == 0) refKminmers.erase(delRefKminmerIt);
      i++;
    } else {
      if (toDelete) {
        // add back deleted seed
        const auto& seed = seedInfos[seedIndex];
        auto [refKminmerIt, refKminmerInserted] = refKminmers.try_emplace(seed.hash, 1);
        if (!refKminmerInserted) refKminmerIt->second++;
      } else {
        // delete inserted seed
        const auto& seed = seedInfos[seedIndex];
        auto refKminmerIt = refKminmers.find(seed.hash);
        if (refKminmerIt == refKminmers.end()) {
          std::cerr << "Error: backtracking seed " << seedIndex << " not found in refKminmers" << std::endl;
          exit(1);
        }
        refKminmerIt->second--;
        if (refKminmerIt->second == 0) refKminmers.erase(refKminmerIt);
      }
    }
  }

  for (size_t i = 0; i < readScoresBacktrackIndex; ++i) {
    const auto& [readIdx, score] = readScoresBacktrack[i];
    readScores[readIdx] = score;
  }

  for (size_t i = 0; i < coveredKminmersBacktrackIndex; ++i) {
    const auto& [hash, readIndex, toErase] = coveredKminmersBacktrack[coveredKminmersBacktrackIndex - i - 1];
    if (toErase) {
      auto coveredKminmerIt = coveredKminmers.find(hash);
      if (coveredKminmerIt == coveredKminmers.end()) {
        std::cerr << "Error: backtracking coveredKminmer " << hash << " not found in coveredKminmers" << std::endl;
        exit(1);
      }
      coveredKminmerIt->second.erase(readIndex);
      if (coveredKminmerIt->second.empty()) coveredKminmers.erase(coveredKminmerIt);
    } else {
      auto [coveredKminmerIt, coveredKminmerInserted] = coveredKminmers.try_emplace(hash, std::unordered_set<uint32_t>({readIndex}));
      if (!coveredKminmerInserted) coveredKminmerIt->second.insert(readIndex);
    }
  }

}

void mgsr::ThreadsManager::computeKminmerCoverage() {
  absl::flat_hash_map<size_t, std::vector<std::pair<uint32_t, uint32_t>>> seedmerToReads;
  for (uint32_t i = 0; i < reads.size(); ++i) {
    for (const auto& seedmer : reads[i].uniqueSeedmers) {
      for (const auto& seedmerIndex : seedmer.second) {
        seedmerToReads[seedmer.first].emplace_back(i, seedmerIndex);
      }
    }
  }

  std::vector<uint32_t> readScores(reads.size(), 0);
  absl::flat_hash_map<size_t, std::unordered_set<uint32_t>> coveredKminmers;
  absl::flat_hash_map<size_t, uint32_t> refKminmers;
  size_t dfsIndex = 0;

  computeKminmerCoverageHelper(liteTree->root, readScores, seedmerToReads, coveredKminmers, refKminmers, dfsIndex);



}


std::vector<uint32_t> mgsr::ThreadsManager::getScoresAtNode(const std::string& nodeId) const {
  std::vector<uint32_t> curNodeScores(reads.size(), 0);
  getScoresAtNode(nodeId, curNodeScores);
  return curNodeScores;
}

void mgsr::ThreadsManager::getScoresAtNode(const std::string& nodeId, std::vector<uint32_t>& curNodeScores) const {
  if (curNodeScores.size() != reads.size()) {
    curNodeScores.resize(reads.size(), 0);
  }

  MgsrLiteNode* currentNode = liteTree->allLiteNodes[nodeId];

  // Get the path from root the current node
  std::vector<MgsrLiteNode*> nodePath;
  while (currentNode->parent != nullptr) {
    nodePath.push_back(currentNode);
    currentNode = currentNode->parent;
  }
  nodePath.push_back(currentNode);
  std::reverse(nodePath.begin(), nodePath.end());

  for (const auto& node : nodePath) {
    auto nodeDfsIndex = node->dfsIndex;
    for (size_t threadId = 0; threadId < numThreads; ++threadId) {
      const auto threadReadStart = threadRanges[threadId].first;
      if (lowMemory) {
        const auto& curNodeThreadScoreDeltasLowMemory = node->readScoreDeltasLowMemory[threadId];
        for (const auto& scoreDelta : curNodeThreadScoreDeltasLowMemory) {
          const auto [trailingDelta, readIndex, numTrailing, currentScoreDelta] = scoreDelta;
          curNodeScores[threadReadStart + readIndex] += currentScoreDelta;
          for (size_t i = 0; i < numTrailing; ++i) {
            int16_t curDecodedTrailingDelta = scoreDelta.decodeTrailingDelta(i);
            // if curDecodedTrailingDelta == -8, it means the read is not changed, so we don't need to apply the trailing delta
            if (curDecodedTrailingDelta > -8) {
              curNodeScores[threadReadStart + readIndex + i + 1] += currentScoreDelta + curDecodedTrailingDelta;
            }
          }
        }
      } else {
        const auto& curNodeThreadScoreDeltas = node->readScoreDeltas[threadId];
        for (const auto& scoreDelta : curNodeThreadScoreDeltas) {
          curNodeScores[threadReadStart + scoreDelta.readIndex] = scoreDelta.scoreDelta;
        }
      }
    }
  }
}

void mgsr::ThreadsManager::printStats() {
  std::cout << "Thread\tnumReads\tnumInitialized\tnumAdded\tnumRemoved\tnumUpdated" << std::endl;
  for (size_t i = 0; i < threadRanges.size(); ++i) {
    std::cout << i << "\t" << threadRanges[i].second - threadRanges[i].first << "\t" << readMinichainsInitialized[i] << "\t" << readMinichainsAdded[i] << "\t" << readMinichainsRemoved[i] << "\t" << readMinichainsUpdated[i] << std::endl;
  }
}

mgsr::squareEM::squareEM(
  mgsr::ThreadsManager& threadsManager,
  mgsr::MgsrLiteTree& liteTree,
  const std::string& prefix,
  uint32_t overlapCoefficientCutoff
) {
  this->prefix = prefix;
  numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

  auto& kminmerOverlapCoefficients = threadsManager.kminmerOverlapCoefficients;
  auto& readSeedmersDuplicatesIndex = threadsManager.readSeedmersDuplicatesIndex;
  auto& reads = threadsManager.reads;
  size_t numReads = reads.size();


  identicalGroups.swap(threadsManager.identicalGroups);
  identicalNodeToGroup.swap(threadsManager.identicalNodeToGroup);


  std::vector<std::pair<std::string, double>> kminmerOverlapCoefficientsVector;
  for (const auto& [nodeId, kminmerOverlapCoefficient] : kminmerOverlapCoefficients) {
    if (identicalNodeToGroup.find(nodeId) != identicalNodeToGroup.end()) continue;
    kminmerOverlapCoefficientsVector.emplace_back(nodeId, kminmerOverlapCoefficient);
  }

  std::sort(kminmerOverlapCoefficientsVector.begin(), kminmerOverlapCoefficientsVector.end(), [](const auto& a, const auto& b) {
    if (a.second == b.second) {
      return a.first < b.first;
    }
    return a.second > b.second;
  });


  std::vector<std::string> significantOverlapNodeIds{kminmerOverlapCoefficientsVector[0].first};
  size_t curRank = 0;
  for (size_t i = 1; i < kminmerOverlapCoefficientsVector.size(); ++i) {
    const auto& [nodeId, kminmerOverlapCoefficient] = kminmerOverlapCoefficientsVector[i];
    if (kminmerOverlapCoefficient != kminmerOverlapCoefficientsVector[i - 1].second) {
      ++curRank;
    }
    if (curRank >= overlapCoefficientCutoff) break;
    // if (curRank >= 9999999) break;
    significantOverlapNodeIds.push_back(nodeId);
  }


  // get score matrix to find ambiguous nodes with identical scores
  std::vector<std::vector<uint32_t>> scoreMatrix(significantOverlapNodeIds.size(), std::vector<uint32_t>(numReads, 0));
  const size_t chunkSize = (significantOverlapNodeIds.size() + numThreads - 1) / numThreads;
  std::vector<std::pair<size_t, size_t>> threadRanges(numThreads);
  for (size_t i = 0; i < numThreads; ++i) {
    size_t start = i * chunkSize;
    size_t end = (i == numThreads - 1) ? significantOverlapNodeIds.size() : (i + 1) * chunkSize;
    if (start < significantOverlapNodeIds.size()) {
      threadRanges[i].first = start;
      threadRanges[i].second = end;
    }
  }


  // tbb::parallel_for(size_t(0), numThreads, [&](size_t threadIdx) {
  //   const auto& range = threadRanges[threadIdx];
  //   for (size_t i = range.first; i < range.second; ++i) {
  //     auto significantNodeId = significantOverlapNodeIds[i];
  //     threadsManager.getScoresAtNode(significantNodeId, scoreMatrix[i]);
  //     double curTestScore = 0;
  //     for (size_t j = 0; j < scoreMatrix[i].size(); ++j) {
  //       if (scoreMatrix[i][j] == reads[j].maxScore && reads[j].maxScore > 0) {
  //         reads[j].epp++;
  //       }
  //     }
  //   }
  // });

  // std::vector<double> testScores(significantOverlapNodeIds.size(), 0);
  // std::vector<std::vector<uint32_t>> readToMaxScoreNodeIndices(numReads);
  // tbb::parallel_for(size_t(0), numThreads, [&](size_t threadIdx) {
  //   const auto& range = threadRanges[threadIdx];
  //   for (size_t i = range.first; i < range.second; ++i) {
  //     double curTestScore = 0;
  //     bool toDebug = false;
  //     if (significantOverlapNodeIds[i] == "England/MILK-344FEB3/2022|OV817379.1|2022-01-26") {
  //       toDebug = false;
  //     }
  //     for (size_t j = 0; j < scoreMatrix[i].size(); ++j) {
  //       if (scoreMatrix[i][j] == reads[j].maxScore && reads[j].maxScore > 0) {
  //         readToMaxScoreNodeIndices[j].push_back(i);
  //         double curReadScore = static_cast<double>(readSeedmersDuplicatesIndex[j].size()) / ((reads[j].seedmersList.size() - reads[j].maxScore + 1) * pow(static_cast<double>(reads[j].epp), 2));
  //         curTestScore += curReadScore;
  //         if (toDebug) {
  //           std::cout << "Read " << j << " parsimonious. Score:" << scoreMatrix[i][j] << ". Max:" << reads[j].maxScore << " epp:" << reads[j].epp << " duplicates:" << readSeedmersDuplicatesIndex[j].size() << std::endl;
  //         }
  //       } else if (reads[j].maxScore > 0 && toDebug) {
  //         std::cout << "Read " << j << " NOT parsimonious. Score:" << scoreMatrix[i][j] << ". Max:" << reads[j].maxScore << " epp:" << reads[j].epp << " duplicates:" << readSeedmersDuplicatesIndex[j].size() << std::endl;
  //       }
  //     }
  //     testScores[i] = curTestScore;
  //   }
  // });



  std::vector<std::vector<uint32_t>> readToMaxScoreNodeIndices(numReads);
  tbb::parallel_for(size_t(0), numThreads, [&](size_t threadIdx) {
    const auto& range = threadRanges[threadIdx];
    for (size_t i = range.first; i < range.second; ++i) {
      auto significantNodeId = significantOverlapNodeIds[i];
      threadsManager.getScoresAtNode(significantNodeId, scoreMatrix[i]);
    }
  });

  // std::ofstream ofs(prefix + ".testScores.txt");
  // for (size_t i = 0; i < testScores.size(); ++i) {
  //   ofs << significantOverlapNodeIds[i];
  //   if (identicalGroups.find(significantOverlapNodeIds[i]) != identicalGroups.end()) {
  //     for (const auto& member : identicalGroups[significantOverlapNodeIds[i]]) {
  //       ofs << "," << member;
  //     }
  //   }
  //   ofs << " " << std::fixed << std::setprecision(10) 
  //       << kminmerOverlapCoefficientsVector[i].second << " "
  //       << testScores[i] << std::endl;
  // }
  // ofs.close();

  // uint32_t readsRemaining;
  // for (size_t i = 0; i < numReads; ++i) {
  //   if (!readToMaxScoreNodeIndices[i].empty()) {
  //     readsRemaining++;
  //   }
  // }
  // std::unordered_map<uint32_t, double> selectedNodeIndices;
  // while (readsRemaining > 0) {
  //   double curMaxScore = -1;
  //   uint32_t curMaxScoreIndex = 0;
  //   for (size_t i = 0; i < testScores.size(); ++i) {
  //     if (testScores[i] > curMaxScore && selectedNodeIndices.find(i) == selectedNodeIndices.end()) {
  //       curMaxScore = testScores[i];
  //       curMaxScoreIndex = i;
  //     }
  //   }
  //   selectedNodeIndices[curMaxScoreIndex] = curMaxScore;
  //   std::cout << "Selected node " << significantOverlapNodeIds[curMaxScoreIndex] << " with score " << curMaxScore << std::endl;

  //   std::vector<uint32_t> selectedNodeScores(numReads, 0);
  //   threadsManager.getScoresAtNode(significantOverlapNodeIds[curMaxScoreIndex], selectedNodeScores);
  //   for (size_t i = 0; i < numReads; ++i) {
  //     if (readToMaxScoreNodeIndices[i].empty()) continue;
  //     if (selectedNodeScores[i] == reads[i].maxScore && reads[i].maxScore > 0) {
  //       double curReadScore = static_cast<double>(readSeedmersDuplicatesIndex[i].size()) / ((reads[i].seedmersList.size() - reads[i].maxScore + 1) * pow(static_cast<double>(reads[i].epp), 2));
  //       for (auto otherNodeIndex : readToMaxScoreNodeIndices[i]) {
  //         if (selectedNodeIndices.find(otherNodeIndex) != selectedNodeIndices.end()) continue;
  //         testScores[otherNodeIndex] -= curReadScore;
  //       }
  //       std::vector<uint32_t>().swap(readToMaxScoreNodeIndices[i]);
  //       --readsRemaining;
  //     }
  //   }
  // }
  
  // std::vector<std::pair<std::string, double>> selectedNodeIndicesVector;
  // for (const auto& [index, score] : selectedNodeIndices) {
  //   selectedNodeIndicesVector.emplace_back(significantOverlapNodeIds[index], score);
  // }
  // std::sort(selectedNodeIndicesVector.begin(), selectedNodeIndicesVector.end(), [](const auto& a, const auto& b) {
  //   return a.second > b.second;
  // });
  // std::ofstream ofsSelectedNodesTestScores(prefix + ".selectedNodesTestScores.txt");
  // for (const auto& [nodeId, score] : selectedNodeIndicesVector) {
  //   ofsSelectedNodesTestScores << nodeId;
  //   if (identicalGroups.find(nodeId) != identicalGroups.end()) {
  //     for (const auto& member : identicalGroups[nodeId]) {
  //       ofsSelectedNodesTestScores << "," << member;
  //     }
  //   }
  //   ofsSelectedNodesTestScores << " " << score << std::endl;
  // }
  // ofsSelectedNodesTestScores.close();



  // exit(0);
  


  // clear memory that are no longer needed
  for (auto [nodeId, node] : liteTree.allLiteNodes) {
    decltype(node->readScoreDeltas)().swap(node->readScoreDeltas);
    decltype(node->readScoreDeltasLowMemory)().swap(node->readScoreDeltasLowMemory);
  }



  std::unordered_map<std::vector<uint32_t>, std::vector<std::string_view>, mgsr::VectorHash> scoresToNodeIds;
  for (uint32_t i = 0; i < scoreMatrix.size(); ++i) {
    scoresToNodeIds[scoreMatrix[i]].push_back(significantOverlapNodeIds[i]);
  }

  for (const auto& [scores, nodeIds] : scoresToNodeIds) {
    std::string nodeIdToKeep(nodeIds[0]);
    if (nodeIds.size() > 1) {
      for (size_t i = 1; i < nodeIds.size(); ++i) {
        std::string nodeIdToMerge(nodeIds[i]);
        identicalGroups[nodeIdToKeep].push_back(nodeIdToMerge);
        identicalNodeToGroup[nodeIdToMerge] = nodeIdToKeep;
        if (identicalGroups.find(nodeIdToMerge) != identicalGroups.end()) {
          for (const auto& member : identicalGroups[nodeIdToMerge]) {
            identicalGroups[nodeIdToKeep].push_back(member);
            identicalNodeToGroup[member] = nodeIdToKeep;
          }
          identicalGroups.erase(nodeIdToMerge);
        }
      }
    }
  }

  nodes.resize(scoresToNodeIds.size());
  probs.resize(threadsManager.numPassedReads, scoresToNodeIds.size());
  std::vector<std::pair<std::vector<uint32_t>, std::vector<std::string_view>>> scoresToNodeIdsVector(scoresToNodeIds.size());
  size_t scoresToNodeIdsVecIndex = 0;
  for (auto& [scores, nodeIds] : scoresToNodeIds) {
    scoresToNodeIdsVector[scoresToNodeIdsVecIndex].first = std::move(scores);
    scoresToNodeIdsVector[scoresToNodeIdsVecIndex].second = std::move(nodeIds);
    ++scoresToNodeIdsVecIndex;
  }
  const size_t chunkSizeScoresToNodeIdsVector = (scoresToNodeIdsVector.size() + numThreads - 1) / numThreads;
  std::vector<std::pair<size_t, size_t>> threadRangesScoresToNodeIdsVector(numThreads);
  for (size_t i = 0; i < numThreads; ++i) {
    size_t start = i * chunkSizeScoresToNodeIdsVector;
    size_t end = (i == numThreads - 1) ? scoresToNodeIdsVector.size() : (i + 1) * chunkSizeScoresToNodeIdsVector;
    if (start < scoresToNodeIdsVector.size()) {
      threadRangesScoresToNodeIdsVector[i].first = start;
      threadRangesScoresToNodeIdsVector[i].second = end;
    }
  }
  tbb::parallel_for(size_t(0), numThreads, [&](size_t threadIdx) {
    const auto& range = threadRangesScoresToNodeIdsVector[threadIdx];
    for (size_t i = range.first; i < range.second; ++i) {
    const auto& [scores, nodeIds] = scoresToNodeIdsVector[i];
    nodes[i] = nodeIds[0];
    size_t passedReadIndex = 0;
    for (size_t j = 0; j < numReads; ++j) {
      if (reads[j].readType != mgsr::ReadType::PASS) continue;
      probs(passedReadIndex, i) = pow(errorRate, reads[j].seedmersList.size() - scores[j]) * pow(1 - errorRate, scores[j]);
        ++passedReadIndex;
      }
    }
  });


  
  numNodes = nodes.size();
  props = Eigen::VectorXd::Constant(numNodes, 1.0 / static_cast<double>(numNodes));
  props0 = Eigen::VectorXd::Zero(numNodes);
  props1 = Eigen::VectorXd::Zero(numNodes);
  props2 = Eigen::VectorXd::Zero(numNodes);
  propsSq = Eigen::VectorXd::Zero(numNodes);
  r = Eigen::VectorXd::Zero(numNodes);
  v = Eigen::VectorXd::Zero(numNodes);
  denoms = Eigen::VectorXd::Zero(threadsManager.numPassedReads);
  inverseDenoms = Eigen::VectorXd::Zero(threadsManager.numPassedReads);
  readDuplicates = Eigen::VectorXd::Zero(threadsManager.numPassedReads);

  size_t passedReadIndex = 0;
  for (size_t i = 0; i < readSeedmersDuplicatesIndex.size(); ++i) {
    if (reads[i].readType != mgsr::ReadType::PASS) continue;
    readDuplicates(passedReadIndex) = readSeedmersDuplicatesIndex[i].size();
    ++passedReadIndex;
  }

  const size_t propsChunkSize = (numNodes + numThreads - 1) / numThreads;
  threadsRangeByProps.resize(numThreads);
  for (size_t i = 0; i < numThreads; ++i) {
    size_t start = i * propsChunkSize;
    size_t end = (i == numThreads - 1) ? numNodes : (i + 1) * propsChunkSize;
    if (start < numNodes) {
      threadsRangeByProps[i].first = start;
      threadsRangeByProps[i].second = end;
    }
  }

  std::cout << "Probs matrix size: " << probs.rows() << "x" << probs.cols() << std::endl;
  std::cout << "Props vector size: " << props.size() << std::endl;
}
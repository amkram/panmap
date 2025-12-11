

#include "index_single_mode.hpp"
#include "panmap_utils.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include "zstd_compression.hpp"
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
#include <chrono>

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


  // check block sequence objects and coordinates
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequenceDynamic = blockSequences.sequence;
  const std::vector<char>& blockExistsDynamic = blockSequences.blockExists;
  const std::vector<char>& blockStrandDynamic = blockSequences.blockStrand;
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
            if (index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") and global scalar coord  " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].second[j] << " ?= " << sequenceBruteForce[blockId][i].second[j] << std::endl;
              std::cerr << "Degapped scalar coord: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
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
            if (index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) and global scalar coord " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].first << " ?= " << sequenceBruteForce[blockId][i].first << std::endl;
              std::cerr << "Degapped scalar coord: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
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
            if (index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", -1) and global scalar coord " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].first << " ?= " << sequenceBruteForce[blockId][i].first << std::endl;
              std::cerr << "Degapped scalar coord: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", -1)... passed: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
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
            if (index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) != localScalarCoordBruteForce) {
              std::cerr << "Degapped scalar coord mismatch at global coord (" << blockId << ", " << i << ", " << j << ") and global scalar coord " << globalScalarCoord << std::endl;
              std::cerr << "Nuc: " << sequenceDynamic[blockId][i].second[j] << " ?= " << sequenceBruteForce[blockId][i].second[j] << std::endl;
              std::cerr << "Degapped scalar coord: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " != " << localScalarCoordBruteForce << std::endl;
              std::cerr << "Local block coord: " << globalScalarCoord - globalCoords.getBlockStartScalar(blockId) << std::endl;
              std::exit(1);
            } else {
              if (printCorrect && node->identifier == nodeToDebug) std::cout << "\tIdentical scalar coord at (" << blockId << ", " << i << ", " << j << ")... passed: " << index_single_mode::degapGlobal(globalScalarCoord, degapCoordIndex) << " == " << localScalarCoordBruteForce << std::endl;
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
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersBruteForceRevcomp = seeding::rollingSyncmers(seeding::revcomp(ungappedSequence), k, s, open, t, false);
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmersDynamic;
  for (size_t i = 0; i < refOnSyncmers.size(); i++) {
    if (refOnSyncmers[i].has_value()) {
      const auto& [hash, endPos, isReverse] = refOnSyncmers[i].value();
      syncmersDynamic.emplace_back(std::make_tuple(hash, isReverse, true, index_single_mode::degapGlobal(i, degapCoordIndex)));
    }
  }

  if (syncmersBruteForce.size() != syncmersBruteForceRevcomp.size()) {
    std::cout << "Syncmer brute force and revcomp count mismatch: forward " << syncmersBruteForce.size() << " != reverse " << syncmersBruteForceRevcomp.size() << std::endl;
    std::exit(1);
  } else {
    for (size_t i = 0; i < syncmersBruteForce.size(); i++) {
      const auto& [hashFwd, isReverseFwd, isSeedFwd, startPosFwd] = syncmersBruteForce[i];
      const auto& [hashRev, isReverseRev, isSeedRev, startPosRev] = syncmersBruteForceRevcomp[syncmersBruteForceRevcomp.size() - i - 1];
      if (hashFwd != hashRev || isReverseFwd == isReverseRev || startPosFwd != ungappedSequence.size() - k - startPosRev) {
        std::cout << "Syncmer brute force and revcomp mismatch at " << i << "th syncmer: forward (" << hashFwd << ", " << startPosFwd << ", " << isReverseFwd << ") != reverse (" << hashRev << ", " << index_single_mode::degapGlobal(ungappedSequence.size() - k - startPosRev, degapCoordIndex) << ", " << isReverseRev << ")" << std::endl;
        std::exit(1);
      }
    }
    std::cout << node->identifier << " seq len: " << ungappedSequence.size() << "... num syncmers: " << syncmersBruteForce.size() << std::flush;
  }


  // check all syncmers
  if (syncmersDynamic.size() != syncmersBruteForce.size()) {
    std::cout << "Syncmer count mismatch: dynamic " << syncmersDynamic.size() << " != brute force " << syncmersBruteForce.size() << std::endl;
    std::cout << "Dynamic syncmers: ";
    for (const auto& syncmer : syncmersDynamic) {
      std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << index_single_mode::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
    }
    std::cout << std::endl;
    std::cout << "Brute force syncmers: ";
    for (const auto& syncmer : syncmersBruteForce) {
      std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << index_single_mode::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
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
      std::cout << "Dynamic syncmers: ";
      for (const auto& syncmer : syncmersDynamic) {
        std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << index_single_mode::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
      }
      std::cout << std::endl;
      std::cout << "Brute force syncmers: ";
      for (const auto& syncmer : syncmersBruteForce) {
        std::cout << "(" << std::get<0>(syncmer) << ", " << std::get<3>(syncmer) << ", " << std::get<1>(syncmer) << "," << index_single_mode::regapGlobal(std::get<3>(syncmer), regapCoordIndex) << ") ";
      }
      std::cout << std::endl;
      std::exit(1);
    }
    if (index_single_mode::regapGlobal(startPos, regapCoordIndex) != *curSyncmerOnMapIt) {
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
  //     kminmersDynamic.emplace_back(std::make_tuple(hash, index_single_mode::degapGlobal(startPos, degapCoordIndex), index_single_mode::degapGlobal(endPos, degapCoordIndex), isReverse));
  //   }
  // }

  // if (kminmersDynamic.size() != kminmersBruteForce.size()) {
  //   std::cout << "K-min-mer count mismatch: dynamic " << kminmersDynamic.size() << " != brute force " << kminmersBruteForce.size() << std::endl;
  //   std::cout << "Dynamic k-min-mers: ";
  //   for (const auto& kminmer : kminmersDynamic) {
  //     std::cout << "(" << std::get<0>(kminmer) << ", " << std::get<1>(kminmer) << ", " << index_single_mode::regapGlobal(std::get<1>(kminmer), regapCoordIndex) << ", " << std::get<2>(kminmer) << ", " << index_single_mode::regapGlobal(std::get<2>(kminmer), regapCoordIndex) << ", " << std::get<3>(kminmer) << ") ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << "Brute force k-min-mers: ";
  //   for (const auto& kminmer : kminmersBruteForce) {
  //     std::cout << "(" << std::get<0>(kminmer) << ", " << std::get<1>(kminmer) << ", " << index_single_mode::regapGlobal(std::get<1>(kminmer), regapCoordIndex) << ", " << std::get<2>(kminmer) << ", " << index_single_mode::regapGlobal(std::get<2>(kminmer), regapCoordIndex) << ", " << std::get<3>(kminmer) << ") ";
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
      index_single_mode::degapGlobal(seedInfos[infoIndex].startPos, degapCoordIndex),
      index_single_mode::degapGlobal(seedInfos[infoIndex].endPos, degapCoordIndex),
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

uint64_t index_single_mode::degapGlobal(const uint64_t& globalCoord, const std::map<uint64_t, uint64_t>& degapCoordsIndex) {
  auto coordIt = degapCoordsIndex.upper_bound(globalCoord);
  if (coordIt == degapCoordsIndex.begin()) {
      return 0;
  }
  return globalCoord - std::prev(coordIt)->second;
}

uint64_t index_single_mode::regapGlobal(const uint64_t& localCoord, const std::map<uint64_t, uint64_t>& regapCoordsIndex) {
  auto coordIt = regapCoordsIndex.upper_bound(localCoord);
  if (coordIt == regapCoordsIndex.begin()) {
      return 0;
  }
  return localCoord + std::prev(coordIt)->second;
}

void index_single_mode::updateGapMapStep(
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

void index_single_mode::updateGapMap(
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

void index_single_mode::invertGapMap(
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

void index_single_mode::revertGapMapInversions(
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

void index_single_mode::makeCoordIndex(
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

std::vector<panmapUtils::NewSyncmerRange> index_single_mode::IndexBuilder::computeNewSyncmerRangesWalk(
  panmanUtils::Node* node,
  size_t dfsIndex,
  const panmapUtils::BlockSequences& blockSequences,
  const std::vector<char>& blockExistsDelayed,
  const std::vector<char>& blockStrandDelayed,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersChangeRecord,
  const std::vector<std::optional<seeding::rsyncmer_t>>& refOnSyncmers,
  std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers
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

std::vector<panmapUtils::NewSyncmerRange> index_single_mode::IndexBuilder::computeNewSyncmerRangesJump(
  panmanUtils::Node* node,
  size_t dfsIndex,
  const panmapUtils::BlockSequences& blockSequences,
  const std::vector<char>& blockExistsDelayed,
  const std::vector<char>& blockStrandDelayed,
  const panmapUtils::GlobalCoords& globalCoords,
  const std::map<uint64_t, uint64_t>& gapMap,
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>>& localMutationRanges,
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>>& blockOnSyncmersChangeRecord,
  const std::vector<std::optional<seeding::rsyncmer_t>>& refOnSyncmers,
  std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers
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

std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> index_single_mode::IndexBuilder::computeNewKminmerRanges(
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
  const uint64_t dfsIndex
) {
  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges;
  auto l = indexBuilder.getL();
  if (refOnSyncmersMap.size() < l) {
    // no k-min-mers to compute, erase all k-min-mers
    return newKminmerRanges;
  }

  
  std::sort(refOnSyncmersChangeRecord.begin(), refOnSyncmersChangeRecord.end(), [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  if (l == 1) {
    for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
      if (changeType == panmapUtils::seedChangeType::DEL) {
        continue;
      } else {
        auto it = refOnSyncmersMap.find(syncmerPos);
        if (it == refOnSyncmersMap.end()) {
          std::cerr << "Error: syncmer position not found in refOnSyncmersMap when computing new k-min-mer ranges.\n";
          exit(1);
        }
        newKminmerRanges.emplace_back(it, it);
      }
    }
    return newKminmerRanges;
  }


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
    
    
    // expand to the left
    int offset = 1;
    if (l == 2) {
      if (curBegIt != refOnSyncmersMap.begin() && !newKminmerRanges.empty() && *curBegIt <= *(newKminmerRanges.back().second)) {
          curBegIt = newKminmerRanges.back().first;
          newKminmerRanges.pop_back();
      }
    } else {
      while (offset < l - 1 && curBegIt != refOnSyncmersMap.begin()) {
        if (!newKminmerRanges.empty() && *curBegIt <= *(newKminmerRanges.back().second)) {
          curBegIt = newKminmerRanges.back().first;
          newKminmerRanges.pop_back();
          break;
        }
        --curBegIt;
        offset++;
      }
    }

    offset = 1;
    if (l == 2) {
      while (curEndIt != refOnSyncmersMap.end()) {
        if (syncmerChangeIndex != refOnSyncmersChangeRecord.size() - 1) {
          const auto& [nextSyncmerPos, nextChangeType, nextRsyncmer] = refOnSyncmersChangeRecord[syncmerChangeIndex + 1];
          if (nextChangeType != panmapUtils::seedChangeType::DEL) {
            if (*curEndIt >= nextSyncmerPos) {
              auto nextIt = refOnSyncmersMap.lower_bound(nextSyncmerPos);
              curEndIt = nextSyncmerPos == *refOnSyncmersMap.rbegin() ? refOnSyncmersMap.end() : std::next(nextIt);
              syncmerChangeIndex++;
              continue;
            }
          } else {
            if (nextSyncmerPos < *refOnSyncmersMap.rbegin()) {
              if (*curEndIt >= nextSyncmerPos) {
                curEndIt = std::next(curEndIt);
                syncmerChangeIndex++;
                continue;
              }
            }
            //  else {
            //   curEndIt = refOnSyncmersMap.end();
            //   break;
            // }
          }
        }
        break;
      }
    } else {
      // expand to the right
      while (offset < l - 1 && curEndIt != refOnSyncmersMap.end()) {
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
    }


    newKminmerRanges.emplace_back(curBegIt, curEndIt);
    if (curEndIt == refOnSyncmersMap.end()) {
      break;
    }
    syncmerChangeIndex++;
  }
  return newKminmerRanges;
}

// Overload for parallel building that uses BuildState instead of member variables
std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> index_single_mode::IndexBuilder::computeNewKminmerRanges(
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
  BuildState& state,
  const uint64_t dfsIndex
) {
  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges;
  auto l = indexBuilder.getL();
  if (state.refOnSyncmersMap.size() < static_cast<size_t>(l)) {
    return newKminmerRanges;
  }

  std::sort(refOnSyncmersChangeRecord.begin(), refOnSyncmersChangeRecord.end(), [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  if (l == 1) {
    for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
      if (changeType == panmapUtils::seedChangeType::DEL) {
        continue;
      } else {
        auto it = state.refOnSyncmersMap.find(syncmerPos);
        if (it == state.refOnSyncmersMap.end()) {
          std::cerr << "Error: syncmer position not found in refOnSyncmersMap when computing new k-min-mer ranges.\n";
          exit(1);
        }
        newKminmerRanges.emplace_back(it, it);
      }
    }
    return newKminmerRanges;
  }

  int64_t syncmerChangeIndex = 0;
  while (syncmerChangeIndex < static_cast<int64_t>(refOnSyncmersChangeRecord.size())) {
    const auto& [syncmerPos, changeType, rsyncmer] = refOnSyncmersChangeRecord[syncmerChangeIndex];
    std::set<uint64_t>::iterator curBegIt, curEndIt;
    if (changeType == panmapUtils::seedChangeType::DEL) {
      if (syncmerPos < *state.refOnSyncmersMap.begin()) {
        syncmerChangeIndex++;
        continue;
      } else {
        curEndIt = state.refOnSyncmersMap.upper_bound(syncmerPos);
        curBegIt = std::prev(curEndIt);
      }
    } else {
      auto it = state.refOnSyncmersMap.lower_bound(syncmerPos);
      curBegIt = syncmerPos == *state.refOnSyncmersMap.begin() ? state.refOnSyncmersMap.begin() : std::prev(it);
      curEndIt = syncmerPos == *state.refOnSyncmersMap.rbegin() ? state.refOnSyncmersMap.end() : std::next(it);
    }
    
    int offset = 1;
    if (l == 2) {
      if (curBegIt != state.refOnSyncmersMap.begin() && !newKminmerRanges.empty() && *curBegIt <= *(newKminmerRanges.back().second)) {
          curBegIt = newKminmerRanges.back().first;
          newKminmerRanges.pop_back();
      }
    } else {
      while (offset < l - 1 && curBegIt != state.refOnSyncmersMap.begin()) {
        if (!newKminmerRanges.empty() && *curBegIt <= *(newKminmerRanges.back().second)) {
          curBegIt = newKminmerRanges.back().first;
          newKminmerRanges.pop_back();
          break;
        }
        --curBegIt;
        offset++;
      }
    }

    offset = 1;
    if (l == 2) {
      while (curEndIt != state.refOnSyncmersMap.end()) {
        if (syncmerChangeIndex != static_cast<int64_t>(refOnSyncmersChangeRecord.size()) - 1) {
          const auto& [nextSyncmerPos, nextChangeType, nextRsyncmer] = refOnSyncmersChangeRecord[syncmerChangeIndex + 1];
          if (nextChangeType != panmapUtils::seedChangeType::DEL) {
            if (*curEndIt >= nextSyncmerPos) {
              auto nextIt = state.refOnSyncmersMap.lower_bound(nextSyncmerPos);
              curEndIt = nextSyncmerPos == *state.refOnSyncmersMap.rbegin() ? state.refOnSyncmersMap.end() : std::next(nextIt);
              syncmerChangeIndex++;
              continue;
            }
          } else {
            if (nextSyncmerPos < *state.refOnSyncmersMap.rbegin()) {
              if (*curEndIt >= nextSyncmerPos) {
                curEndIt = std::next(curEndIt);
                syncmerChangeIndex++;
                continue;
              }
            }
          }
        }
        break;
      }
    } else {
      while (offset < l - 1 && curEndIt != state.refOnSyncmersMap.end()) {
        if (syncmerChangeIndex != static_cast<int64_t>(refOnSyncmersChangeRecord.size()) - 1) {
          const auto& [nextSyncmerPos, nextChangeType, nextRsyncmer] = refOnSyncmersChangeRecord[syncmerChangeIndex + 1];
          if (nextChangeType != panmapUtils::seedChangeType::DEL) {
            if (*curEndIt >= nextSyncmerPos) {
              auto nextIt = state.refOnSyncmersMap.lower_bound(nextSyncmerPos);
              curEndIt = nextSyncmerPos == *state.refOnSyncmersMap.rbegin() ? state.refOnSyncmersMap.end() : std::next(nextIt);
              syncmerChangeIndex++;
              offset = 1;
              continue;
            }
          } else {
            if (nextSyncmerPos < *state.refOnSyncmersMap.rbegin()) {
              if (*curEndIt >= nextSyncmerPos) {
                curEndIt = std::next(curEndIt);
                syncmerChangeIndex++;
                offset = 1;
                continue;
              }
            }
          }
        }
        ++curEndIt;
        offset++;
      }
    }

    newKminmerRanges.emplace_back(curBegIt, curEndIt);
    if (curEndIt == state.refOnSyncmersMap.end()) {
      break;
    }
    syncmerChangeIndex++;
  }
  return newKminmerRanges;
}

// Process a single node without recursing to children and without backtracking
// Used for building up state along a path in parallel processing
void index_single_mode::IndexBuilder::processSingleNodeNoBacktrack(
  panmanUtils::Node *node,
  std::unordered_set<std::string_view>& emptyNodes,
  panmapUtils::BlockSequences &blockSequences,
  std::vector<char> &blockExistsDelayed,
  std::vector<char> &blockStrandDelayed,
  panmapUtils::GlobalCoords &globalCoords,
  std::map<uint64_t, uint64_t> &gapMap,
  std::unordered_set<uint64_t> &invertedBlocks,
  uint64_t dfsIndex
) {
  // Same as buildIndexHelper but: no recursion, no backtracking
  std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> blockMutationRecord;
  std::vector<std::tuple<panmapUtils::Coordinate, char, char>> nucMutationRecord;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapUpdates;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunBlockInversionBacktracks;
  std::vector<std::pair<uint64_t, bool>> invertedBlocksBacktracks;
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> localMutationRanges;
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>> refOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>> blockOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, uint64_t>> refOnKminmersChangeRecord;
  std::vector<uint32_t> potentialSyncmerDeletions;

  if (node->blockMutation.empty() && node->nucMutation.empty()) {
    emptyNodes.insert(std::string_view(node->identifier));
  }

  applyMutations(node, dfsIndex, blockSequences, invertedBlocks, globalCoords, localMutationRanges, 
                 blockMutationRecord, nucMutationRecord, gapRunUpdates, invertedBlocksBacktracks, 
                 potentialSyncmerDeletions, blockExistsDelayed, blockStrandDelayed);

  std::sort(gapRunUpdates.begin(), gapRunUpdates.end(), 
            [](const auto& a, const auto& b) { return a.second.first < b.second.first; });
  updateGapMap(node, dfsIndex, gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  std::vector<uint64_t> invertedBlocksVec(invertedBlocks.begin(), invertedBlocks.end());
  std::sort(invertedBlocksVec.begin(), invertedBlocksVec.end());
  for (const auto& blockId : invertedBlocksVec) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
    invertGapMap(gapMap, {beg, end}, gapRunBlockInversionBacktracks, gapMapUpdates);
  }

  // Compute new syncmer ranges
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = 
      computeNewSyncmerRangesJump(node, dfsIndex, blockSequences, blockExistsDelayed, blockStrandDelayed, 
                                   globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord, refOnSyncmers, blockOnSyncmers);

  // Process syncmers
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    if (localRangeSeq.size() >= static_cast<size_t>(indexBuilder.getK())) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(localRangeSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = localRangeCoordToGlobalScalarCoords[startPos];
        auto endPosGlobal = localRangeCoordToGlobalScalarCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = localRangeCoordToBlockId[startPos];
        bool wasSeed = refOnSyncmers[startPosGlobal].has_value();
        if (!wasSeed && isSeed) {
          refOnSyncmersMap.insert(startPosGlobal);
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
          if (blockOnSyncmers[curBlockId].empty()) blockOnSyncmers.erase(curBlockId);
        } else if (wasSeed && isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::SUB, refOnSyncmers[startPosGlobal].value());
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }
    for (uint64_t pos : seedsToDelete) {
      if (refOnSyncmers[pos].has_value()) {
        refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers[pos].value());
        refOnSyncmers[pos] = std::nullopt;
        refOnSyncmersMap.erase(pos);
      }
    }
  }

  // Handle potential syncmer deletions
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

  // Handle block deletions
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      if (blockOnSyncmers.find(blockId) != blockOnSyncmers.end()) {
        for (uint64_t pos : blockOnSyncmers[blockId]) {
          refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers[pos].value());
          refOnSyncmers[pos] = std::nullopt;
          refOnSyncmersMap.erase(pos);
        }
        blockOnSyncmers.erase(blockId);
      }
    }
  }

  // Compute k-min-mers
  auto k = indexBuilder.getK();
  auto l = indexBuilder.getL();
  std::vector<uint64_t> deletedSeedHashes;
  std::vector<uint64_t> addedSeedHashes;
  std::vector<std::pair<uint64_t, uint64_t>> substitutedSeedHashes;

  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges = 
      computeNewKminmerRanges(refOnSyncmersChangeRecord, dfsIndex);

  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto [curIt, endIt] = newKminmerRanges[i];
    auto indexingIt = curIt;
    size_t forwardHash = 0;
    size_t reverseHash = 0;

    std::vector<size_t> startingSyncmerHashes;
    for (size_t j = 0; j < static_cast<size_t>(l); j++) {
      if (curIt == refOnSyncmersMap.end()) break;
      else if (endIt != refOnSyncmersMap.end() && *curIt == *endIt && j != static_cast<size_t>(l) - 1) break;
      startingSyncmerHashes.push_back(refOnSyncmers[*curIt].value().hash);
      if (j != static_cast<size_t>(l) - 1) ++curIt;
    }
    if (startingSyncmerHashes.size() < static_cast<size_t>(l)) continue;

    if (l == 1) {
      auto syncmerRev = refOnSyncmers[*curIt].value().isReverse;
      if (syncmerRev) {
        forwardHash = std::numeric_limits<size_t>::max();
        reverseHash = startingSyncmerHashes[0];
      } else {
        forwardHash = startingSyncmerHashes[0];
        reverseHash = std::numeric_limits<size_t>::max();
      }
    } else {
      for (size_t j = 0; j < static_cast<size_t>(l); j++) {
        forwardHash = seeding::rol(forwardHash, k) ^ startingSyncmerHashes[j];
        reverseHash = seeding::rol(reverseHash, k) ^ startingSyncmerHashes[l - j - 1];
      }
    }

    auto& curRefOnKminmer = refOnKminmers[*indexingIt];
    uint64_t newHash = std::min(forwardHash, reverseHash);
    if (forwardHash != reverseHash) {
      bool substitution = curRefOnKminmer.has_value();
      if (substitution) {
        uint64_t oldHash = curRefOnKminmer.value();
        substitutedSeedHashes.emplace_back(oldHash, newHash);
      } else {
        addedSeedHashes.push_back(newHash);
      }
      curRefOnKminmer = newHash;
    } else {
      if (curRefOnKminmer.has_value()) {
        deletedSeedHashes.push_back(curRefOnKminmer.value());
        curRefOnKminmer = std::nullopt;
      }
    }

    if (curIt == endIt) continue;
    
    while (curIt != endIt) {
      ++curIt;
      if (curIt == refOnSyncmersMap.end()) break;
      forwardHash = seeding::rol(forwardHash, k) ^ seeding::rol(refOnSyncmers[*indexingIt].value().hash, k * l) ^ refOnSyncmers[*curIt].value().hash;
      reverseHash = seeding::ror(reverseHash, k) ^ seeding::ror(refOnSyncmers[*indexingIt].value().hash, k) ^ seeding::rol(refOnSyncmers[*curIt].value().hash, k * (l-1));
      ++indexingIt;

      auto& curRefOnKminmer2 = refOnKminmers[*indexingIt];
      uint64_t newHash2 = std::min(forwardHash, reverseHash);
      if (forwardHash != reverseHash) {
        bool substitution = curRefOnKminmer2.has_value();
        if (substitution) {
          substitutedSeedHashes.emplace_back(curRefOnKminmer2.value(), newHash2);
        } else {
          addedSeedHashes.push_back(newHash2);
        }
        curRefOnKminmer2 = newHash2;
      } else {
        if (curRefOnKminmer2.has_value()) {
          deletedSeedHashes.push_back(curRefOnKminmer2.value());
          curRefOnKminmer2 = std::nullopt;
        }
      }
    }
  }

  // Handle end-of-range k-minmer deletions
  if (!newKminmerRanges.empty() && newKminmerRanges.back().second == refOnSyncmersMap.end()) {
    auto delIt = newKminmerRanges.back().second;
    for (size_t j = 0; j < static_cast<size_t>(l) - 1; j++) {
      --delIt;
      if (refOnKminmers[*delIt].has_value()) {
        deletedSeedHashes.push_back(refOnKminmers[*delIt].value());
        refOnKminmers[*delIt] = std::nullopt;
      }
      if (delIt == refOnSyncmersMap.begin()) break;
    }
  }

  // Handle deleted syncmers that still have k-minmers
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL && refOnKminmers[syncmerPos].has_value()) {
      deletedSeedHashes.push_back(refOnKminmers[syncmerPos].value());
      refOnKminmers[syncmerPos] = std::nullopt;
    }
  }

  // Track seed hash counts
  auto& curNodeCounts = nodeSeedCounts[dfsIndex];
  if (node->parent != nullptr && nodeToDfsIndex.count(node->parent->identifier)) {
    curNodeCounts = nodeSeedCounts[nodeToDfsIndex[node->parent->identifier]];
  }
  for (uint64_t hash : addedSeedHashes) curNodeCounts[hash]++;
  for (uint64_t hash : deletedSeedHashes) {
    auto& count = curNodeCounts[hash];
    count--;
    if (count <= 0) curNodeCounts.erase(hash);
  }
  for (const auto& [oldHash, newHash] : substitutedSeedHashes) {
    auto& oldCount = curNodeCounts[oldHash];
    oldCount--;
    if (oldCount <= 0) curNodeCounts.erase(oldHash);
    curNodeCounts[newHash]++;
  }

  // Revert gap map inversions (no full backtrack - we keep the state)
  revertGapMapInversions(gapRunBlockInversionBacktracks, gapMap);

  // Update delayed block states (keep for next node)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExistsDelayed[blockId] = newExists;
    blockStrandDelayed[blockId] = newStrand;
  }

  nodeToDfsIndex[node->identifier] = dfsIndex;
}

void index_single_mode::IndexBuilder::buildIndexHelper(
  panmanUtils::Node *node,
  std::unordered_set<std::string_view>& emptyNodes,
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

  if (node->blockMutation.empty() && node->nucMutation.empty()) {
    emptyNodes.insert(std::string_view(node->identifier));
  }

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
    newSyncmerRanges = computeNewSyncmerRangesJump(node, dfsIndex, blockSequences,  blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord, refOnSyncmers, blockOnSyncmers);
  } else {
    newSyncmerRanges = computeNewSyncmerRangesWalk(node, dfsIndex, blockSequences,  blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord, refOnSyncmers, blockOnSyncmers);
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
  // Use hash-based tracking (not indices) for LiteIndex compatibility
  std::vector<uint64_t> deletedSeedHashes;
  std::vector<uint64_t> addedSeedHashes;
  std::vector<std::pair<uint64_t, uint64_t>> substitutedSeedHashes;
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
    if (l == 1) {
      auto syncmerRev = refOnSyncmers[*curIt].value().isReverse;
      if (syncmerRev) {
        forwardHash = std::numeric_limits<size_t>::max();
        reverseHash = startingSyncmerHashes[0];
      } else {
        forwardHash = startingSyncmerHashes[0];
        reverseHash = std::numeric_limits<size_t>::max();
      }
      if (reverseHash == forwardHash) {
        std::cerr << "Error: syncmer hash collision detected in k-min-mer computation.\n";
        exit(1);
      }
    } else {
      for (size_t j = 0; j < l; j++) {
        forwardHash = seeding::rol(forwardHash, k) ^ startingSyncmerHashes[j];
        reverseHash = seeding::rol(reverseHash, k) ^ startingSyncmerHashes[l - j - 1];
      }
    }

    auto& curRefOnKminmer = refOnKminmers[*indexingIt];
    uint64_t newHash = std::min(forwardHash, reverseHash);
    if (forwardHash != reverseHash) {
      bool substitution = curRefOnKminmer.has_value();
      if (substitution) {
        uint64_t oldHash = curRefOnKminmer.value();
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, oldHash);
        substitutedSeedHashes.emplace_back(oldHash, newHash);
      } else {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
        addedSeedHashes.push_back(newHash);
      }
      // Store hash directly instead of index
      curRefOnKminmer = newHash;
    } else {
      if (curRefOnKminmer.has_value()) {
        uint64_t oldHash = curRefOnKminmer.value();
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
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
      uint64_t newHash2 = std::min(forwardHash, reverseHash);
      if (forwardHash != reverseHash) {
        bool substitution = curRefOnKminmer.has_value();
        if (substitution) {
          uint64_t oldHash = curRefOnKminmer.value();
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, oldHash);
          substitutedSeedHashes.emplace_back(oldHash, newHash2);
        } else {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
          addedSeedHashes.push_back(newHash2);
        }
        curRefOnKminmer = newHash2;
      } else {
        if (curRefOnKminmer.has_value()) {
          uint64_t oldHash = curRefOnKminmer.value();
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, oldHash);
          deletedSeedHashes.push_back(oldHash);
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
        uint64_t oldHash = refOnKminmers[*delIt].value();
        refOnKminmersChangeRecord.emplace_back(*delIt, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        refOnKminmers[*delIt] = std::nullopt;
      }
      if (delIt == refOnSyncmersMap.begin()) break;
    }
  }
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL && refOnKminmers[syncmerPos].has_value()) {
      uint64_t oldHash = refOnKminmers[syncmerPos].value();
      refOnKminmersChangeRecord.emplace_back(syncmerPos, panmapUtils::seedChangeType::DEL, oldHash);
      deletedSeedHashes.push_back(oldHash);
      refOnKminmers[syncmerPos] = std::nullopt;
    }
  }

  // No sorting needed - we're using hashes directly, not positions
 
  //  Track seed hash counts for this node (for LiteIndex format)
  //  Start with parent's counts (if not root)
  auto& curNodeCounts = nodeSeedCounts[dfsIndex];
  if (node->parent != nullptr && nodeToDfsIndex.count(node->parent->identifier)) {
    curNodeCounts = nodeSeedCounts[nodeToDfsIndex[node->parent->identifier]];
  }
  // Apply additions (using hashes directly)
  for (uint64_t hash : addedSeedHashes) {
    curNodeCounts[hash]++;
  }
  // Apply deletions
  for (uint64_t hash : deletedSeedHashes) {
    auto& count = curNodeCounts[hash];
    count--;
    if (count <= 0) curNodeCounts.erase(hash);
  }
  // Apply substitutions (old deleted, new added)
  for (const auto& [oldHash, newHash] : substitutedSeedHashes) {
    auto& oldCount = curNodeCounts[oldHash];
    oldCount--;
    if (oldCount <= 0) curNodeCounts.erase(oldHash);
    curNodeCounts[newHash]++;
  }
  
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
    buildIndexHelper(child, emptyNodes, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks, dfsIndex);
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


void index_single_mode::IndexBuilder::buildIndex() {
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
  std::unordered_set<std::string_view> emptyNodes;
  uint64_t dfsIndex = 0;
  buildIndexHelper(T->root, emptyNodes, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks, dfsIndex);
  std::cout << std::endl;

  size_t numNodes = T->allNodes.size();

  // NOTE: SeedInfo list is no longer populated - it was never used by lite placement.
  // The seed hashes are stored per-node in the struct-of-arrays format below.

  // Add block infos to index
  capnp::List<BlockRange>::Builder blockRangesBuilder = liteTreeBuilder.initBlockRanges(globalCoords.globalCoords.size());
  for (size_t i = 0; i < globalCoords.globalCoords.size(); i++) {
    blockRangesBuilder[i].setRangeBeg(globalCoords.getBlockStartScalar(i));
    blockRangesBuilder[i].setRangeEnd(globalCoords.getBlockEndScalar(i));
  }

  // Add node to index
  capnp::List<LiteNode>::Builder nodesBuilder = liteTreeBuilder.initLiteNodes(numNodes);
  for (const auto& [nodeId, node] : T->allNodes) {
    auto idx = nodeToDfsIndex[nodeId];
    nodesBuilder[idx].setId(nodeId);
    if (node->parent == nullptr) {
      nodesBuilder[idx].setParentIndex(0);
    } else {
      nodesBuilder[idx].setParentIndex(nodeToDfsIndex[node->parent->identifier]);
    }
    if (emptyNodes.find(nodeId) != emptyNodes.end()) {
      nodesBuilder[idx].setIdenticalToParent(true);
    } else {
      nodesBuilder[idx].setIdenticalToParent(false);
    }
  }

  // Build struct-of-arrays format for seed changes (LiteIndex V3)
  // First pass: collect all changes and compute totals
  std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>> nodeChanges(numNodes);
  uint64_t totalChanges = 0;
  uint32_t largestNodeChangeCount = 0;

  // Build reverse lookup: dfsIndex -> Node*
  std::vector<panmanUtils::Node*> dfsIndexToNode(numNodes, nullptr);
  for (const auto& [nodeId, node] : T->allNodes) {
    auto it = nodeToDfsIndex.find(nodeId);
    if (it != nodeToDfsIndex.end()) {
      dfsIndexToNode[it->second] = node;
    }
  }

  std::cout << "Computing seed deltas..." << std::endl;

  // Compute per-node deltas sequentially
  for (size_t nodeIdx = 0; nodeIdx < numNodes; ++nodeIdx) {
    const auto& curCounts = nodeSeedCounts[nodeIdx];

    // Find parent's counts using the reverse lookup
    const std::unordered_map<uint64_t, int64_t>* parentCounts = nullptr;
    panmanUtils::Node* node = dfsIndexToNode[nodeIdx];
    if (node && node->parent != nullptr) {
      auto parentIt = nodeToDfsIndex.find(node->parent->identifier);
      if (parentIt != nodeToDfsIndex.end()) {
        parentCounts = &nodeSeedCounts[parentIt->second];
      }
    }

    auto& changes = nodeChanges[nodeIdx];

    // Seeds in current node
    for (const auto& kv : curCounts) {
      uint64_t hash = kv.first;
      int64_t childCount = kv.second;
      int64_t parentCount = 0;
      if (parentCounts) {
        auto it = parentCounts->find(hash);
        if (it != parentCounts->end()) parentCount = it->second;
      }
      if (childCount != parentCount) changes.emplace_back(hash, parentCount, childCount);
    }

    // Seeds only in parent (deleted in child)
    if (parentCounts) {
      for (const auto& kv : *parentCounts) {
        uint64_t hash = kv.first;
        int64_t parentCount = kv.second;
        if (curCounts.find(hash) == curCounts.end()) changes.emplace_back(hash, parentCount, 0);
      }
    }

    // Aggregate totals inline
    totalChanges += changes.size();
    if (changes.size() > largestNodeChangeCount) largestNodeChangeCount = changes.size();
    if (nodeIdx % 1000 == 0) {
      std::cout << "\rProcessed node " << nodeIdx << "/" << numNodes << std::flush;
    }
  }
  std::cout << std::endl;

  // Set totals
  indexBuilder.setTotalSeedChanges(totalChanges);
  indexBuilder.setLargestNodeChangeCount(largestNodeChangeCount);

  // Allocate flat arrays
  auto seedChangeHashesBuilder = indexBuilder.initSeedChangeHashes(totalChanges);
  auto seedChangeParentCountsBuilder = indexBuilder.initSeedChangeParentCounts(totalChanges);
  auto seedChangeChildCountsBuilder = indexBuilder.initSeedChangeChildCounts(totalChanges);
  auto nodeChangeOffsetsBuilder = indexBuilder.initNodeChangeOffsets(numNodes + 1);

  // Write flat arrays
  uint32_t offset = 0;
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    nodeChangeOffsetsBuilder.set(nodeIdx, offset);
    for (const auto& [hash, parentCount, childCount] : nodeChanges[nodeIdx]) {
      seedChangeHashesBuilder.set(offset, hash);
      seedChangeParentCountsBuilder.set(offset, parentCount);
      seedChangeChildCountsBuilder.set(offset, childCount);
      offset++;
    }
  }
  nodeChangeOffsetsBuilder.set(numNodes, offset);

  // Compute and write pre-indexed genome metrics
  auto genomeMagnitudeSquaredBuilder = indexBuilder.initGenomeMagnitudeSquared(numNodes);
  auto genomeUniqueSeedCountBuilder = indexBuilder.initGenomeUniqueSeedCount(numNodes);
  auto genomeTotalSeedFrequencyBuilder = indexBuilder.initGenomeTotalSeedFrequency(numNodes);

  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    const auto& counts = nodeSeedCounts[nodeIdx];
    double magnitudeSquared = 0.0;
    int64_t totalFrequency = 0;
    
    for (const auto& [hash, count] : counts) {
      magnitudeSquared += static_cast<double>(count) * static_cast<double>(count);
      totalFrequency += count;
    }
    
    genomeMagnitudeSquaredBuilder.set(nodeIdx, magnitudeSquared);
    genomeUniqueSeedCountBuilder.set(nodeIdx, counts.size());
    genomeTotalSeedFrequencyBuilder.set(nodeIdx, totalFrequency);
  }

  std::cout << "Finished building index! Total seed changes: " << totalChanges << std::endl;
}

void index_single_mode::IndexBuilder::writeIndex(const std::string& path) {
  std::cout << "Serializing index..." << std::endl;
  
  // Serialize Cap'n Proto message to flat array
  kj::ArrayPtr<const kj::ArrayPtr<const capnp::word>> segments = outMessage.getSegmentsForOutput();
  
  // Calculate total size
  size_t totalWords = 0;
  for (auto seg : segments) {
    totalWords += seg.size();
  }
  size_t totalBytes = totalWords * sizeof(capnp::word);
  
  // Allocate buffer for serialized message using messageToFlatArray
  kj::Array<capnp::word> flatArray = capnp::messageToFlatArray(outMessage);
  const void* data = flatArray.begin();
  size_t dataSize = flatArray.size() * sizeof(capnp::word);
  
  std::cout << "Writing ZSTD-compressed index to " << path << " (" << dataSize << " bytes uncompressed)..." << std::endl;
  
  // Use ZSTD compression to write
  if (!panmap_zstd::compressToFile(data, dataSize, path, 3, 0)) {
    std::cerr << "Error: failed to write compressed index to " << path << std::endl;
    std::exit(1);
  }
  
  std::cout << "Index written to " << path << std::endl;
}

// ============================================================================
// Parallel Index Building Implementation
// ============================================================================

uint64_t index_single_mode::IndexBuilder::computeSubtreeSize(panmanUtils::Node* node) {
  uint64_t size = 1;  // Count this node
  for (auto* child : node->children) {
    size += computeSubtreeSize(child);
  }
  subtreeSizes_[node->identifier] = size;
  return size;
}

void index_single_mode::IndexBuilder::processNode(
  panmanUtils::Node *node,
  BuildState& state,
  panmapUtils::GlobalCoords &globalCoords,
  std::unordered_set<std::string_view>& localEmptyNodes,
  uint64_t dfsIndex,
  BacktrackInfo* backtrackInfo
) {
  // Local storage if backtrackInfo is null
  std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> localBlockMutationRecord;
  std::vector<std::tuple<panmapUtils::Coordinate, char, char>> localNucMutationRecord;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> localGapRunBacktracks;
  std::vector<std::pair<uint64_t, bool>> localInvertedBlocksBacktracks;
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>> localRefOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, uint64_t, panmapUtils::seedChangeType>> localBlockOnSyncmersChangeRecord;
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, uint64_t>> localRefOnKminmersChangeRecord;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> localGapRunBlockInversionBacktracks;

  // References to use
  auto& blockMutationRecord = backtrackInfo ? backtrackInfo->blockMutationRecord : localBlockMutationRecord;
  auto& nucMutationRecord = backtrackInfo ? backtrackInfo->nucMutationRecord : localNucMutationRecord;
  auto& gapRunBacktracks = backtrackInfo ? backtrackInfo->gapRunBacktracks : localGapRunBacktracks;
  auto& invertedBlocksBacktracks = backtrackInfo ? backtrackInfo->invertedBlocksBacktracks : localInvertedBlocksBacktracks;
  auto& refOnSyncmersChangeRecord = backtrackInfo ? backtrackInfo->refOnSyncmersChangeRecord : localRefOnSyncmersChangeRecord;
  auto& blockOnSyncmersChangeRecord = backtrackInfo ? backtrackInfo->blockOnSyncmersChangeRecord : localBlockOnSyncmersChangeRecord;
  auto& refOnKminmersChangeRecord = backtrackInfo ? backtrackInfo->refOnKminmersChangeRecord : localRefOnKminmersChangeRecord;
  auto& gapRunBlockInversionBacktracks = backtrackInfo ? backtrackInfo->gapRunBlockInversionBacktracks : localGapRunBlockInversionBacktracks;

  if (backtrackInfo) {
      backtrackInfo->clear();
  }

  // For building gap map
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapMapUpdates;
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>> gapRunUpdates;

  // For computing new syncmers
  std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> localMutationRanges;

  // Nuc deletions on block without block mutation
  std::vector<uint32_t> potentialSyncmerDeletions;

  if (node->blockMutation.empty() && node->nucMutation.empty()) {
    localEmptyNodes.insert(std::string_view(node->identifier));
  }

  blockMutationRecord.reserve(node->blockMutation.size());
  nucMutationRecord.reserve(node->nucMutation.size() * 6);
  potentialSyncmerDeletions.reserve(node->nucMutation.size() * 6);
  localMutationRanges.reserve(node->blockMutation.size() + node->nucMutation.size() * 6);
  gapRunUpdates.reserve(node->nucMutation.size() * 6 + node->blockMutation.size() * 10);
  
  applyMutations(node, dfsIndex, state.blockSequences, state.invertedBlocks, globalCoords, 
                 localMutationRanges, blockMutationRecord, nucMutationRecord, gapRunUpdates, 
                 invertedBlocksBacktracks, potentialSyncmerDeletions, 
                 state.blockExistsDelayed, state.blockStrandDelayed);
  
  blockMutationRecord.shrink_to_fit();
  nucMutationRecord.shrink_to_fit();
  potentialSyncmerDeletions.shrink_to_fit();
  localMutationRanges.shrink_to_fit();
  gapRunUpdates.shrink_to_fit();

  std::sort(gapRunUpdates.begin(), gapRunUpdates.end(), 
            [&](const auto& a, const auto& b) { return a.second.first < b.second.first; });
  updateGapMap(node, dfsIndex, state.gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  std::vector<uint64_t> invertedBlocksVec(state.invertedBlocks.begin(), state.invertedBlocks.end());
  std::sort(invertedBlocksVec.begin(), invertedBlocksVec.end());
  for (const auto& blockId : invertedBlocksVec) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
    invertGapMap(state.gapMap, {beg, end}, gapRunBlockInversionBacktracks, gapMapUpdates);
  }

  bool useJump = true;
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges;
  if (useJump) {
    newSyncmerRanges = computeNewSyncmerRangesJump(node, dfsIndex, state.blockSequences, 
                                                    state.blockExistsDelayed, state.blockStrandDelayed, 
                                                    globalCoords, state.gapMap, localMutationRanges, 
                                                    blockOnSyncmersChangeRecord, state.refOnSyncmers, state.blockOnSyncmers);
  } else {
    newSyncmerRanges = computeNewSyncmerRangesWalk(node, dfsIndex, state.blockSequences, 
                                                    state.blockExistsDelayed, state.blockStrandDelayed, 
                                                    globalCoords, state.gapMap, localMutationRanges, 
                                                    blockOnSyncmersChangeRecord, state.refOnSyncmers, state.blockOnSyncmers);
  }

  // Processing syncmers - matches original exactly
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    if (localRangeSeq.size() >= static_cast<size_t>(indexBuilder.getK())) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(localRangeSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = localRangeCoordToGlobalScalarCoords[startPos];
        auto endPosGlobal = localRangeCoordToGlobalScalarCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = localRangeCoordToBlockId[startPos];
        bool wasSeed = state.refOnSyncmers[startPosGlobal].has_value();
        if (!wasSeed && isSeed) {
          state.refOnSyncmersMap.insert(startPosGlobal);
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::ADD, seeding::rsyncmer_t());
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::ADD);
          state.refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
          state.blockOnSyncmers[curBlockId].insert(startPosGlobal);
        } else if (wasSeed && !isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::DEL, state.refOnSyncmers[startPosGlobal].value());
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::DEL);
          state.refOnSyncmers[startPosGlobal] = std::nullopt;
          state.refOnSyncmersMap.erase(startPosGlobal);
          state.blockOnSyncmers[curBlockId].erase(startPosGlobal);
          if (state.blockOnSyncmers[localRangeCoordToBlockId[startPos]].empty()) {
            state.blockOnSyncmers.erase(localRangeCoordToBlockId[startPos]);
          }
        } else if (wasSeed && isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::SUB, state.refOnSyncmers[startPosGlobal].value());
          state.refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }

    for (uint64_t pos : seedsToDelete) {
      if (!state.refOnSyncmers[pos].has_value()) {
        std::cerr << "Error: refOnSyncmers[" << pos << "] is null" << std::endl;
        std::exit(1);
      }
      refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, state.refOnSyncmers[pos].value());
      state.refOnSyncmers[pos] = std::nullopt;
      state.refOnSyncmersMap.erase(pos);
    }
  }

  // Handle potential syncmer deletions (matches original)
  if (useJump) {
    for (uint32_t pos : potentialSyncmerDeletions) {
      if (state.refOnSyncmers[pos].has_value()) {
        refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, state.refOnSyncmers[pos].value());
        state.refOnSyncmers[pos] = std::nullopt;
        const auto blockId = globalCoords.getBlockIdFromScalar(pos);
        state.blockOnSyncmers[blockId].erase(pos);
        if (state.blockOnSyncmers[blockId].empty()) state.blockOnSyncmers.erase(blockId);
        blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);
        state.refOnSyncmersMap.erase(pos);
      }
    }
  }

  // Handle block mutations that delete entire blocks (matches original)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      if (state.blockOnSyncmers.find(blockId) != state.blockOnSyncmers.end()) {
        for (uint64_t pos : state.blockOnSyncmers[blockId]) {
          refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, state.refOnSyncmers[pos].value());
          blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);
          state.refOnSyncmers[pos] = std::nullopt;
          state.refOnSyncmersMap.erase(pos);
        }
        state.blockOnSyncmers.erase(blockId);
      }
    }
  }

  // Processing k-min-mers - matches original exactly
  auto k = indexBuilder.getK();
  auto l = indexBuilder.getL();
  // LOCK-FREE: Store hashes directly instead of indices
  std::vector<uint64_t> deletedSeedHashes;
  std::vector<uint64_t> addedSeedHashes;
  std::vector<std::pair<uint64_t, uint64_t>> substitutedSeedHashes; // <oldHash, newHash>

  std::vector<std::pair<std::set<uint64_t>::iterator, std::set<uint64_t>::iterator>> newKminmerRanges = 
      computeNewKminmerRanges(refOnSyncmersChangeRecord, state, dfsIndex);

  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto beg = *newKminmerRanges[i].first;
    auto end = *newKminmerRanges[i].second;
    if (newKminmerRanges[i].second != state.refOnSyncmersMap.end() && beg > end) {
      std::cerr << "Error: beg (" << beg << ") > end (" << end << ") in node " << node->identifier << std::endl;
      std::exit(1);
    }
    if (i != 0 && *newKminmerRanges[i-1].second > beg) {
      std::cerr << "Error: newKminmerRanges[" << i-1 << "].second > newKminmerRanges[" << i << "].first" << std::endl;
      std::exit(1);
    }
  }

  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto [curIt, endIt] = newKminmerRanges[i];
    auto indexingIt = curIt;
    size_t forwardHash = 0;
    size_t reverseHash = 0;

    std::vector<size_t> startingSyncmerHashes;
    for (size_t j = 0; j < static_cast<size_t>(l); j++) {
      if (curIt == state.refOnSyncmersMap.end()) {
        break;
      } else if (endIt != state.refOnSyncmersMap.end() && *curIt == *endIt && j != static_cast<size_t>(l) - 1) {
        break;
      }
      startingSyncmerHashes.push_back(state.refOnSyncmers[*curIt].value().hash);
      if (j != static_cast<size_t>(l) - 1) ++curIt;
    }
    if (startingSyncmerHashes.size() < static_cast<size_t>(l)) {
      continue;
    }
    if (l == 1) {
      auto syncmerRev = state.refOnSyncmers[*curIt].value().isReverse;
      if (syncmerRev) {
        forwardHash = std::numeric_limits<size_t>::max();
        reverseHash = startingSyncmerHashes[0];
      } else {
        forwardHash = startingSyncmerHashes[0];
        reverseHash = std::numeric_limits<size_t>::max();
      }
      if (reverseHash == forwardHash) {
        std::cerr << "Error: syncmer hash collision detected in k-min-mer computation.\n";
        exit(1);
      }
    } else {
      for (size_t j = 0; j < static_cast<size_t>(l); j++) {
        forwardHash = seeding::rol(forwardHash, k) ^ startingSyncmerHashes[j];
        reverseHash = seeding::rol(reverseHash, k) ^ startingSyncmerHashes[l - j - 1];
      }
    }

    auto& curRefOnKminmer = state.refOnKminmers[*indexingIt];
    uint64_t newHash = std::min(forwardHash, reverseHash);
    if (forwardHash != reverseHash) {
      bool substitution = curRefOnKminmer.has_value();
      if (substitution) {
        // Get old hash from stored value and add to substitution list
        uint64_t oldHash = curRefOnKminmer.value(); // This is now the hash, not index
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, oldHash);
        substitutedSeedHashes.emplace_back(oldHash, newHash);
      } else {
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
        addedSeedHashes.push_back(newHash);
      }
      // Store hash directly instead of index
      curRefOnKminmer = newHash;
    } else {
      if (curRefOnKminmer.has_value()) {
        uint64_t oldHash = curRefOnKminmer.value();
        refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        curRefOnKminmer = std::nullopt;
      }
    }

    if (curIt == endIt) continue;
    
    while (curIt != endIt) {
      ++curIt;
      if (curIt == state.refOnSyncmersMap.end()) break;
      forwardHash = seeding::rol(forwardHash, k) ^ seeding::rol(state.refOnSyncmers[*indexingIt].value().hash, k * l) ^ state.refOnSyncmers[*curIt].value().hash;
      reverseHash = seeding::ror(reverseHash, k) ^ seeding::ror(state.refOnSyncmers[*indexingIt].value().hash, k)     ^ seeding::rol(state.refOnSyncmers[*curIt].value().hash, k * (l-1));
      ++indexingIt;

      auto& curRefOnKminmer2 = state.refOnKminmers[*indexingIt];
      uint64_t newHash2 = std::min(forwardHash, reverseHash);
      if (forwardHash != reverseHash) {
        bool substitution = curRefOnKminmer2.has_value();
        if (substitution) {
          uint64_t oldHash = curRefOnKminmer2.value();
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::SUB, oldHash);
          substitutedSeedHashes.emplace_back(oldHash, newHash2);
        } else {
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
          addedSeedHashes.push_back(newHash2);
        }
        // Store hash directly
        curRefOnKminmer2 = newHash2;
      } else {
        if (curRefOnKminmer2.has_value()) {
          uint64_t oldHash = curRefOnKminmer2.value();
          refOnKminmersChangeRecord.emplace_back(*indexingIt, panmapUtils::seedChangeType::DEL, oldHash);
          deletedSeedHashes.push_back(oldHash);
          curRefOnKminmer2 = std::nullopt;
        }
      }
    }
  }

  // Handle end-of-range k-minmer deletions (matches original)
  if (!newKminmerRanges.empty() && newKminmerRanges.back().second == state.refOnSyncmersMap.end()) {
    auto delIt = newKminmerRanges.back().second;
    for (size_t j = 0; j < static_cast<size_t>(l) - 1; j++) {
      --delIt;
      if (state.refOnKminmers[*delIt].has_value()) {
        uint64_t oldHash = state.refOnKminmers[*delIt].value();
        refOnKminmersChangeRecord.emplace_back(*delIt, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        state.refOnKminmers[*delIt] = std::nullopt;
      }
      if (delIt == state.refOnSyncmersMap.begin()) break;
    }
  }

  // Handle deleted syncmers that still have k-minmers (matches original)
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL && state.refOnKminmers[syncmerPos].has_value()) {
      uint64_t oldHash = state.refOnKminmers[syncmerPos].value();
      refOnKminmersChangeRecord.emplace_back(syncmerPos, panmapUtils::seedChangeType::DEL, oldHash);
      deletedSeedHashes.push_back(oldHash);
      state.refOnKminmers[syncmerPos] = std::nullopt;
    }
  }

  // No sorting needed - we're using hashes directly, not positions

  // Store seed deltas for this node (NO COPYING from parent!)
  // nodeSeedDeltas is a shared_ptr so all clones share the same storage
  auto& curDelta = (*state.nodeSeedDeltas)[dfsIndex];
  curDelta.addedHashes = std::move(addedSeedHashes);
  curDelta.deletedHashes = std::move(deletedSeedHashes);
  curDelta.substitutedHashes = std::move(substitutedSeedHashes);

  // Revert gap map inversions (for proper state for children)
  revertGapMapInversions(gapRunBlockInversionBacktracks, state.gapMap);

  // Update delayed block states (for children)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    state.blockExistsDelayed[blockId] = newExists;
    state.blockStrandDelayed[blockId] = newStrand;
  }

  // nodeToDfsIndex is pre-computed, no need to update here

  // Update progress
  uint64_t processed = ++processedNodes_;
  if (processed % 100 == 0 || processed == T->allNodes.size()) {
    std::cout << "\rProcessed " << processed << "/" << T->allNodes.size() << " nodes" << std::flush;
  }
}

// Backtrack all changes made by processNode (reverses state to pre-processNode)
void index_single_mode::IndexBuilder::backtrackNode(
  BuildState& state,
  const BacktrackInfo& backtrackInfo
) {
  // Revert block mutations (sequence exists/strand and delayed states)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : backtrackInfo.blockMutationRecord) {
    state.blockSequences.blockExists[blockId] = oldExists;
    state.blockSequences.blockStrand[blockId] = oldStrand;
    state.blockExistsDelayed[blockId] = oldExists;
    state.blockStrandDelayed[blockId] = oldStrand;
  }

  // Revert nucleotide mutations
  for (const auto& [coord, oldNuc, newNuc] : backtrackInfo.nucMutationRecord) {
    state.blockSequences.setSequenceBase(coord, oldNuc);
  }

  // Revert inverted blocks
  for (const auto& [blockId, del] : backtrackInfo.invertedBlocksBacktracks) {
    if (del) {
      state.invertedBlocks.erase(blockId);
    } else {
      state.invertedBlocks.insert(blockId);
    }
  }

  // Revert gap map changes (in reverse order)
  for (auto it = backtrackInfo.gapRunBacktracks.rbegin(); it != backtrackInfo.gapRunBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      state.gapMap.erase(range.first);
    } else {
      state.gapMap[range.first] = range.second;
    }
  }

  // Revert syncmer changes
  for (const auto& [pos, changeType, rsyncmer] : backtrackInfo.refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      // Was added... need to delete
      state.refOnSyncmers[pos] = std::nullopt;
      state.refOnSyncmersMap.erase(pos);
    } else {
      // Was deleted or replaced... need to restore
      state.refOnSyncmers[pos] = rsyncmer;
      state.refOnSyncmersMap.insert(pos);
    }
  }
  
  // Revert block-on-syncmers changes
  for (const auto& [blockId, pos, changeType] : backtrackInfo.blockOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      // Was added... need to delete
      state.blockOnSyncmers[blockId].erase(pos);
      if (state.blockOnSyncmers[blockId].empty()) {
        state.blockOnSyncmers.erase(blockId);
      }
    } else {
      state.blockOnSyncmers[blockId].insert(pos);
    }
  }

  // Revert k-minmer changes
  for (const auto& [pos, changeType, kminmerHash] : backtrackInfo.refOnKminmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      state.refOnKminmers[pos] = std::nullopt;
    } else {
      state.refOnKminmers[pos] = kminmerHash;
    }
  }
}

// Sequential helper for processing a subtree (uses backtracking like original buildIndexHelper)
// This matches the original sequential DFS pattern: process node, recurse children, then backtrack
void index_single_mode::IndexBuilder::processSubtreeSequential(
  panmanUtils::Node *node,
  BuildState& state,
  panmapUtils::GlobalCoords &globalCoords,
  std::unordered_set<std::string_view>& localEmptyNodes,
  uint64_t dfsIndex,
  BacktrackInfo* /* unused - kept for API compatibility */
) {
  // Always create local backtrack info for this node (matches original buildIndexHelper pattern)
  BacktrackInfo nodeBacktrackInfo;
  
  // Process this node, recording changes for backtracking
  processNode(node, state, globalCoords, localEmptyNodes, dfsIndex, &nodeBacktrackInfo);

  // Recursively process all children (exactly like original buildIndexHelper)
  if (!node->children.empty()) {
    uint64_t offset = dfsIndex + 1;
    for (size_t i = 0; i < node->children.size(); i++) {
      auto* child = node->children[i];
      uint64_t childIdx = offset;
      offset += subtreeSizes_[child->identifier];
      
      // Recursively process this child's subtree
      processSubtreeSequential(child, state, globalCoords, localEmptyNodes, childIdx, nullptr);
    }
  }

  // Backtrack this node's changes AFTER processing all children
  // This matches the original buildIndexHelper pattern
  backtrackNode(state, nodeBacktrackInfo);
}

// Parallel recursive helper - spawns parallel tasks at fork points with large subtrees
// This avoids the overhead of cloning state for small subtrees
void index_single_mode::IndexBuilder::processSubtreeParallel(
  panmanUtils::Node *node,
  BuildState& state,
  panmapUtils::GlobalCoords &globalCoords,
  std::unordered_set<std::string_view>& localEmptyNodes,
  uint64_t dfsIndex,
  tbb::spin_mutex& emptyNodesMutex,
  size_t parallelThreshold  // Only parallelize if combined subtree size >= this
) {
  // Always record changes for backtracking - needed to restore state for siblings
  BacktrackInfo nodeBacktrackInfo;
  
  // Process this node
  processNode(node, state, globalCoords, localEmptyNodes, dfsIndex, &nodeBacktrackInfo);

  // Process children if any
  if (!node->children.empty()) {
    // Calculate total subtree size of children to decide on parallelization
    size_t totalChildSubtreeSize = 0;
    size_t numChildren = node->children.size();
    for (auto* child : node->children) {
      totalChildSubtreeSize += subtreeSizes_[child->identifier];
    }
    
    // Parallelize if:
    // 1. More than one child (fork point)
    // 2. Combined subtree size is large enough
    // 3. At least 2 children have reasonably sized subtrees
    bool shouldParallelize = false;
    if (numChildren > 1 && totalChildSubtreeSize >= parallelThreshold) {
      size_t nonTrivialChildren = 0;
      size_t minSubtreeForParallel = std::max(size_t(50), parallelThreshold / 8);
      for (auto* child : node->children) {
        if (subtreeSizes_[child->identifier] >= minSubtreeForParallel) {
          nonTrivialChildren++;
        }
      }
      shouldParallelize = (nonTrivialChildren >= 2);
    }
    
    if (shouldParallelize) {
      // Parallel processing at this fork point
      // Clone state for each child except the last (which reuses current state)
      std::vector<std::shared_ptr<BuildState>> childStates;
      std::vector<std::unordered_set<std::string_view>> childEmptyNodes;
      
      childStates.reserve(numChildren - 1);
      childEmptyNodes.resize(numChildren);
      
      // Prepare child offsets
      std::vector<uint64_t> childOffsets;
      childOffsets.reserve(numChildren);
      uint64_t offset = dfsIndex + 1;
      for (auto* child : node->children) {
        childOffsets.push_back(offset);
        offset += subtreeSizes_[child->identifier];
      }
      
      // Clone states for all children except last
      // Sort children by subtree size descending so largest subtrees start first
      std::vector<size_t> childOrder(numChildren);
      std::iota(childOrder.begin(), childOrder.end(), 0);
      std::sort(childOrder.begin(), childOrder.end(), [&](size_t a, size_t b) {
        return subtreeSizes_[node->children[a]->identifier] > subtreeSizes_[node->children[b]->identifier];
      });
      
      for (size_t i = 0; i + 1 < numChildren; i++) {
        auto clonedState = std::make_shared<BuildState>(state);
        // Share the nodeSeedDeltas storage (it's already a shared_ptr)
        clonedState->nodeSeedDeltas = state.nodeSeedDeltas;
        childStates.push_back(clonedState);
      }
      
      // Process children in parallel (in order of decreasing subtree size)
      tbb::task_group taskGroup;
      
      for (size_t i = 0; i + 1 < numChildren; i++) {
        size_t childIdx_order = childOrder[i];
        auto* child = node->children[childIdx_order];
        uint64_t childDfsIdx = childOffsets[childIdx_order];
        auto& childState = *childStates[i];
        auto& childEmpty = childEmptyNodes[childIdx_order];
        
        taskGroup.run([this, child, &childState, &globalCoords, &childEmpty, childDfsIdx, 
                       &emptyNodesMutex, parallelThreshold]() {
          processSubtreeParallel(child, childState, globalCoords, childEmpty, childDfsIdx, 
                                 emptyNodesMutex, parallelThreshold);
        });
      }
      
      // Process the smallest subtree (last in sorted order) using current state
      size_t lastChildIdx_order = childOrder.back();
      auto* lastChild = node->children[lastChildIdx_order];
      uint64_t lastChildDfsIdx = childOffsets[lastChildIdx_order];
      processSubtreeParallel(lastChild, state, globalCoords, childEmptyNodes[lastChildIdx_order], 
                             lastChildDfsIdx, emptyNodesMutex, parallelThreshold);
      
      // Wait for all parallel tasks to complete
      taskGroup.wait();
      
      // Merge empty nodes from all children
      for (auto& childEmpty : childEmptyNodes) {
        localEmptyNodes.insert(childEmpty.begin(), childEmpty.end());
      }
    } else {
      // Sequential processing - subtrees too small to benefit from parallelism
      uint64_t offset = dfsIndex + 1;
      for (auto* child : node->children) {
        uint64_t childIdx = offset;
        offset += subtreeSizes_[child->identifier];
        processSubtreeParallel(child, state, globalCoords, localEmptyNodes, childIdx, 
                               emptyNodesMutex, parallelThreshold);
      }
    }
  }

  // Backtrack this node's changes - ALWAYS needed to restore state for siblings
  // This matches the original buildIndexHelper pattern
  backtrackNode(state, nodeBacktrackInfo);
}

void index_single_mode::IndexBuilder::buildIndexParallel(int numThreads) {
  if (numThreads <= 0) {
    numThreads = tbb::this_task_arena::max_concurrency();
  }
  
  // If single-threaded or tree is too small, use sequential
  if (numThreads == 1 || T->root->children.empty()) {
    std::cout << "Using sequential build (threads=" << numThreads << ")" << std::endl;
    buildIndex();
    return;
  }

  std::cout << "Building index in parallel with " << numThreads << " threads..." << std::endl;
  
  auto startTotal = std::chrono::high_resolution_clock::now();
  
  // Initialize state (same as sequential buildIndex)
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences);
  refOnSyncmers.resize(globalCoords.lastScalarCoord + 1);
  refOnKminmers.resize(globalCoords.lastScalarCoord + 1);

  // Step 1: Pre-compute subtree sizes
  computeSubtreeSize(T->root);

  // Step 2: Compute DFS indices for all nodes (needed for nodeSeedDeltas mapping)
  uint64_t dfsCounter = 0;
  std::function<void(panmanUtils::Node*)> computeDfsIndices = [&](panmanUtils::Node* node) {
    nodeToDfsIndex[node->identifier] = dfsCounter++;
    for (auto* child : node->children) {
      computeDfsIndices(child);
    }
  };
  computeDfsIndices(T->root);

  // Step 3: Create initial BuildState
  BuildState state;
  state.blockSequences = std::move(blockSequences);
  state.blockExistsDelayed = state.blockSequences.blockExists;
  state.blockStrandDelayed = state.blockSequences.blockStrand;
  state.gapMap = std::map<uint64_t, uint64_t>{{0, globalCoords.lastScalarCoord}};
  state.refOnSyncmers.resize(globalCoords.lastScalarCoord + 1);
  state.refOnKminmers.resize(globalCoords.lastScalarCoord + 1);
  // Initialize nodeSeedDeltas - stores only changes at each node (NO full counts!)
  state.nodeSeedDeltas = std::make_shared<std::vector<NodeSeedDelta>>(T->allNodes.size());

  // Step 4: Calculate parallelization threshold
  // Lower threshold = more parallelism but more cloning overhead
  // Higher threshold = less parallelism but less overhead
  // Target: create enough parallel tasks for good load balancing without excessive cloning
  size_t totalNodes = T->allNodes.size();
  size_t parallelThreshold = std::max(size_t(50), totalNodes / (numThreads * 8));
  std::cout << "Tree has " << totalNodes << " nodes, parallel threshold set to " << parallelThreshold << std::endl;

  // Step 5: Process tree using dynamic parallelization
  std::unordered_set<std::string_view> emptyNodes;
  tbb::spin_mutex emptyNodesMutex;
  
  auto startDfs = std::chrono::high_resolution_clock::now();
  
  // Use tbb::task_arena to ensure we're using the right number of threads
  tbb::task_arena arena(numThreads);
  arena.execute([&]() {
    processSubtreeParallel(T->root, state, globalCoords, emptyNodes, 0, emptyNodesMutex, parallelThreshold);
  });
  
  auto endDfs = std::chrono::high_resolution_clock::now();
  auto dfsDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endDfs - startDfs).count();
  std::cout << "DFS processing completed in " << dfsDurationMs << " ms" << std::endl;

  // Step 6: Compute node changes directly from deltas using a single DFS pass
  // Instead of building full counts for all nodes, we walk the tree once
  // maintaining a running count and computing changes on-the-fly
  std::cout << "Computing seed changes from deltas (single pass)..." << std::endl;
  
  size_t numNodes = T->allNodes.size();
  std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>> nodeChanges(numNodes);
  std::atomic<uint64_t> totalChanges{0};
  std::atomic<uint32_t> largestNodeChangeCount{0};
  
  // Use iterative DFS with explicit stack to avoid recursion overhead
  // Stack entries: (node, phase) where phase 0 = enter, 1 = process children done
  struct DfsEntry {
    panmanUtils::Node* node;
    size_t childIndex;  // Which child we're processing next
  };
  
  std::vector<DfsEntry> stack;
  stack.reserve(100);  // Reasonable initial depth
  
  // Running count - we modify this as we traverse
  std::unordered_map<uint64_t, int64_t> runningCounts;
  
  stack.push_back({T->root, 0});
  
  while (!stack.empty()) {
    auto& entry = stack.back();
    panmanUtils::Node* node = entry.node;
    uint64_t nodeIdx = nodeToDfsIndex[node->identifier];
    
    if (entry.childIndex == 0) {
      // First visit to this node - apply delta and compute changes
      const auto& delta = (*state.nodeSeedDeltas)[nodeIdx];
      auto& changes = nodeChanges[nodeIdx];
      
      // Track which hashes are modified at this node
      std::unordered_set<uint64_t> modifiedHashes;
      
      // Collect all modified hashes first
      for (uint64_t hash : delta.addedHashes) {
        modifiedHashes.insert(hash);
      }
      for (uint64_t hash : delta.deletedHashes) {
        modifiedHashes.insert(hash);
      }
      for (const auto& [oldHash, newHash] : delta.substitutedHashes) {
        modifiedHashes.insert(oldHash);
        modifiedHashes.insert(newHash);
      }
      
      // Record parent counts for all modified hashes BEFORE applying delta
      std::unordered_map<uint64_t, int64_t> parentCounts;
      for (uint64_t hash : modifiedHashes) {
        auto it = runningCounts.find(hash);
        parentCounts[hash] = (it != runningCounts.end()) ? it->second : 0;
      }
      
      // Apply all delta changes to running counts
      for (uint64_t hash : delta.addedHashes) {
        runningCounts[hash]++;
      }
      for (uint64_t hash : delta.deletedHashes) {
        auto it = runningCounts.find(hash);
        if (it != runningCounts.end()) {
          it->second--;
          if (it->second <= 0) runningCounts.erase(it);
        }
      }
      for (const auto& [oldHash, newHash] : delta.substitutedHashes) {
        // Delete old
        auto oldIt = runningCounts.find(oldHash);
        if (oldIt != runningCounts.end()) {
          oldIt->second--;
          if (oldIt->second <= 0) runningCounts.erase(oldIt);
        }
        // Add new
        runningCounts[newHash]++;
      }
      
      // Now compare final child counts with parent counts for each modified hash
      for (uint64_t hash : modifiedHashes) {
        int64_t parentCount = parentCounts[hash];
        int64_t childCount = 0;
        auto it = runningCounts.find(hash);
        if (it != runningCounts.end()) childCount = it->second;
        
        if (parentCount != childCount) {
          changes.emplace_back(hash, parentCount, childCount);
        }
      }
      
      totalChanges += changes.size();
      uint32_t localSize = changes.size();
      uint32_t prevMax = largestNodeChangeCount.load();
      while (prevMax < localSize && !largestNodeChangeCount.compare_exchange_weak(prevMax, localSize)) {}
    }
    
    // Process children
    if (entry.childIndex < node->children.size()) {
      panmanUtils::Node* child = node->children[entry.childIndex];
      entry.childIndex++;
      stack.push_back({child, 0});
    } else {
      // All children done - backtrack (undo delta)
      const auto& delta = (*state.nodeSeedDeltas)[nodeIdx];
      
      // Undo substitutions (in reverse: undo add new, then undo delete old)
      for (auto it = delta.substitutedHashes.rbegin(); it != delta.substitutedHashes.rend(); ++it) {
        const auto& [oldHash, newHash] = *it;
        // Undo add new
        auto newIt = runningCounts.find(newHash);
        if (newIt != runningCounts.end()) {
          newIt->second--;
          if (newIt->second <= 0) runningCounts.erase(newIt);
        }
        // Undo delete old
        runningCounts[oldHash]++;
      }
      
      // Undo deletions (re-add)
      for (auto it = delta.deletedHashes.rbegin(); it != delta.deletedHashes.rend(); ++it) {
        runningCounts[*it]++;
      }
      
      // Undo additions (re-delete)
      for (auto it = delta.addedHashes.rbegin(); it != delta.addedHashes.rend(); ++it) {
        auto countIt = runningCounts.find(*it);
        if (countIt != runningCounts.end()) {
          countIt->second--;
          if (countIt->second <= 0) runningCounts.erase(countIt);
        }
      }
      
      stack.pop_back();
    }
  }
  
  // Free delta storage - no longer needed
  state.nodeSeedDeltas.reset();
  
  std::cout << "Processed " << numNodes << " nodes" << std::endl;

  // Step 7: Build index output
  LiteTree::Builder liteTreeBuilder = indexBuilder.initLiteTree();

  capnp::List<BlockRange>::Builder blockRangesBuilder = liteTreeBuilder.initBlockRanges(globalCoords.globalCoords.size());
  for (size_t i = 0; i < globalCoords.globalCoords.size(); i++) {
    blockRangesBuilder[i].setRangeBeg(globalCoords.getBlockStartScalar(i));
    blockRangesBuilder[i].setRangeEnd(globalCoords.getBlockEndScalar(i));
  }

  capnp::List<LiteNode>::Builder nodesBuilder = liteTreeBuilder.initLiteNodes(numNodes);
  for (const auto& [nodeId, node] : T->allNodes) {
    auto idx = nodeToDfsIndex[nodeId];
    nodesBuilder[idx].setId(nodeId);
    if (node->parent == nullptr) {
      nodesBuilder[idx].setParentIndex(0);
    } else {
      nodesBuilder[idx].setParentIndex(nodeToDfsIndex[node->parent->identifier]);
    }
    nodesBuilder[idx].setIdenticalToParent(emptyNodes.find(nodeId) != emptyNodes.end());
  }

  std::cout << "\nFinished building index! Total seed changes: " << totalChanges.load() << std::endl;

  // Set totals
  indexBuilder.setTotalSeedChanges(totalChanges.load());
  indexBuilder.setLargestNodeChangeCount(largestNodeChangeCount.load());

  // Allocate and write flat arrays (sequential - Cap'n Proto builders aren't thread-safe)
  auto seedChangeHashesBuilder = indexBuilder.initSeedChangeHashes(totalChanges.load());
  auto seedChangeParentCountsBuilder = indexBuilder.initSeedChangeParentCounts(totalChanges.load());
  auto seedChangeChildCountsBuilder = indexBuilder.initSeedChangeChildCounts(totalChanges.load());
  auto nodeChangeOffsetsBuilder = indexBuilder.initNodeChangeOffsets(numNodes + 1);

  uint32_t writeOffset = 0;
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    nodeChangeOffsetsBuilder.set(nodeIdx, writeOffset);
    for (const auto& [hash, parentCount, childCount] : nodeChanges[nodeIdx]) {
      seedChangeHashesBuilder.set(writeOffset, hash);
      seedChangeParentCountsBuilder.set(writeOffset, parentCount);
      seedChangeChildCountsBuilder.set(writeOffset, childCount);
      writeOffset++;
    }
  }
  nodeChangeOffsetsBuilder.set(numNodes, writeOffset);

  // Compute genome metrics IN PARALLEL
  auto genomeMagnitudeSquaredBuilder = indexBuilder.initGenomeMagnitudeSquared(numNodes);
  auto genomeUniqueSeedCountBuilder = indexBuilder.initGenomeUniqueSeedCount(numNodes);
  auto genomeTotalSeedFrequencyBuilder = indexBuilder.initGenomeTotalSeedFrequency(numNodes);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numNodes),
    [&](const tbb::blocked_range<size_t>& range) {
      for (size_t nodeIdx = range.begin(); nodeIdx != range.end(); ++nodeIdx) {
        const auto& counts = nodeSeedCounts[nodeIdx];
        double magnitudeSquared = 0.0;
        int64_t totalFrequency = 0;
        
        for (const auto& [hash, count] : counts) {
          magnitudeSquared += static_cast<double>(count) * static_cast<double>(count);
          totalFrequency += count;
        }
        
        genomeMagnitudeSquaredBuilder.set(nodeIdx, magnitudeSquared);
        genomeUniqueSeedCountBuilder.set(nodeIdx, counts.size());
        genomeTotalSeedFrequencyBuilder.set(nodeIdx, totalFrequency);
      }
    }, tbb::auto_partitioner());

  auto endTotal = std::chrono::high_resolution_clock::now();
  auto totalDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTotal - startTotal).count();
  
  std::cout << "Total parallel build time: " << totalDurationMs << " ms" << std::endl;
}

int index_single_mode::open_file(const std::string& path) {
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd == -1) {
    std::cerr << "Error: failed to open file " << path << std::endl;
    std::exit(1);
  }
  return fd;
}



#include "index_single_mode.hpp"
#include "panmap_utils.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include "zstd_compression.hpp"
#include "logging.hpp"
#include "gap_map_utils.hpp"
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
#include <atomic>
#include <mutex>
#include <iomanip>



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
  const SyncmerMap& refOnSyncmers,
  std::unordered_map<uint32_t, std::unordered_set<uint64_t>>& blockOnSyncmers,
  uint64_t firstNonGapScalar,
  uint64_t lastNonGapScalar
) {
  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges;
  if (localMutationRanges.empty()) {
    return newSyncmerRanges;
  }

  const std::vector<char>& blockExists = blockSequences.blockExists;
  const std::vector<char>& blockStrand = blockSequences.blockStrand;
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;

  // Schwartzian transform: compute scalar coords once, sort by them
  if (localMutationRanges.size() > 1) {
    std::vector<std::pair<int64_t, size_t>> sortKeys;
    sortKeys.reserve(localMutationRanges.size());
    for (size_t i = 0; i < localMutationRanges.size(); ++i) {
      sortKeys.emplace_back(globalCoords.getScalarFromCoord(localMutationRanges[i].first, blockStrand[localMutationRanges[i].first.primaryBlockId]), i);
    }
    std::sort(sortKeys.begin(), sortKeys.end());
    std::vector<std::pair<panmapUtils::Coordinate, panmapUtils::Coordinate>> sortedRanges;
    sortedRanges.reserve(localMutationRanges.size());
    for (const auto& [key, idx] : sortKeys) {
      sortedRanges.push_back(std::move(localMutationRanges[idx]));
    }
    localMutationRanges = std::move(sortedRanges);
  }

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
      // rightGapMapIt is already at begin, no action needed
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
          // rightGapMapIt is already at begin, no action needed
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
    
    // Pre-calculate expected range size for reservations
    const size_t estimatedSize = (curEndCoordScalar > curCoordScalar) ? 
                                  (curEndCoordScalar - curCoordScalar + 1) : 256;
    
    // Clear and reserve instead of swap (avoids allocation/deallocation)
    localRangeSeq.clear();
    localRangeSeq.reserve(estimatedSize);
    localRangeCoordToGlobalScalarCoords.clear();
    localRangeCoordToGlobalScalarCoords.reserve(estimatedSize);
    seedsToDelete.clear();
    // seedsToDelete is typically small, no need to reserve
    localRangeCoordToBlockId.clear();
    localRangeCoordToBlockId.reserve(estimatedSize);
    
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
      } else if (refOnSyncmers.contains(curCoordScalar)) {
        // Only delete seeds that are INSIDE the genome extent (not in flank regions)
        // Flank regions are before firstNonGapScalar or after lastNonGapScalar
        // Seeds in flank regions are "missing data" not "true gaps"
        if (curCoordScalar >= firstNonGapScalar && curCoordScalar <= lastNonGapScalar) {
          blockOnSyncmers[curCoord.primaryBlockId].erase(curCoordScalar);
          if (blockOnSyncmers[curCoord.primaryBlockId].empty()) blockOnSyncmers.erase(curCoord.primaryBlockId);
          blockOnSyncmersChangeRecord.emplace_back(curCoord.primaryBlockId, curCoordScalar, panmapUtils::seedChangeType::DEL);
          seedsToDelete.push_back(curCoordScalar);
        }
        // If outside extent, the seed is inherited from parent (no deletion)
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
        if (refOnSyncmers.contains(localRangeCoordToGlobalScalarCoords[j])) {
          auto curGlobalScalarCoord = localRangeCoordToGlobalScalarCoords[j];
          // Only delete seeds that are INSIDE the genome extent (not in flank regions)
          if (curGlobalScalarCoord >= firstNonGapScalar && curGlobalScalarCoord <= lastNonGapScalar) {
            auto curBlockId = localRangeCoordToBlockId[j];
            blockOnSyncmers[curBlockId].erase(curGlobalScalarCoord);
            if (blockOnSyncmers[curBlockId].empty()) blockOnSyncmers.erase(curBlockId);
            blockOnSyncmersChangeRecord.emplace_back(curBlockId, curGlobalScalarCoord, panmapUtils::seedChangeType::DEL);
            seedsToDelete.push_back(localRangeCoordToGlobalScalarCoords[j]);
          }
        }
      }
    }
  }

  return newSyncmerRanges;
}

std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> index_single_mode::IndexBuilder::computeNewKminmerRanges(
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
  const uint64_t dfsIndex
) {
  std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> newKminmerRanges;
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
          output::error("syncmer position not found in refOnSyncmersMap when computing new k-min-mer ranges");
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
    index_single_mode::SyncmerSet::iterator curBegIt, curEndIt;
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
std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> index_single_mode::IndexBuilder::computeNewKminmerRanges(
  std::vector<std::tuple<uint64_t, panmapUtils::seedChangeType, seeding::rsyncmer_t>>& refOnSyncmersChangeRecord,
  BuildState& state,
  const uint64_t dfsIndex
) {
  std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> newKminmerRanges;
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
          output::error("syncmer position not found in refOnSyncmersMap when computing new k-min-mer ranges");
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
    index_single_mode::SyncmerSet::iterator curBegIt, curEndIt;
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
  panmapUtils::applyMutations(node, dfsIndex, blockSequences, invertedBlocks, globalCoords, localMutationRanges, blockMutationRecord, nucMutationRecord, gapRunUpdates, invertedBlocksBacktracks, potentialSyncmerDeletions, blockExistsDelayed, blockStrandDelayed, imputeAmb_);
  blockMutationRecord.shrink_to_fit();
  nucMutationRecord.shrink_to_fit();
  potentialSyncmerDeletions.shrink_to_fit();
  localMutationRanges.shrink_to_fit();
  gapRunUpdates.shrink_to_fit();



  std::sort(gapRunUpdates.begin(), gapRunUpdates.end(), [&](const auto& a, const auto& b) { return a.second.first < b.second.first; });
  panmapUtils::updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  std::vector<uint64_t> invertedBlocksVec(invertedBlocks.begin(), invertedBlocks.end());
  std::sort(invertedBlocksVec.begin(), invertedBlocksVec.end());
  for (const auto& blockId : invertedBlocksVec) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);

    gap_map::invertGapMap(gapMap, {beg, end}, gapRunBlockInversionBacktracks, gapMapUpdates);
  }

  // Compute genome extent from gapMap for flank masking
  uint64_t firstNonGapScalar = 0;
  uint64_t lastNonGapScalar = UINT64_MAX;
  if (extentGuard_) {
    auto [fng, lng] = panmapUtils::computeExtentFromGapMap(gapMap, globalCoords.lastScalarCoord, 0);
    firstNonGapScalar = fng;
    lastNonGapScalar = lng;
  }

  // Compute hard mask boundaries (properly skipping non-gap bases, not scalar positions)
  const int flankMaskBp = getFlankMaskBp();
  uint64_t hardMaskStart = firstNonGapScalar;
  uint64_t hardMaskEnd = lastNonGapScalar;
  if (flankMaskBp > 0) {
    auto [hms, hme] = panmapUtils::computeExtentFromGapMap(gapMap, globalCoords.lastScalarCoord, flankMaskBp);
    hardMaskStart = hms;
    hardMaskEnd = hme;
  }

  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = computeNewSyncmerRangesJump(
      node, dfsIndex, blockSequences, blockExistsDelayed, blockStrandDelayed, 
      globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord, 
      refOnSyncmers, blockOnSyncmers, firstNonGapScalar, lastNonGapScalar);

  // processing syncmers
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    
    // Apply HPC if enabled: compress localRangeSeq and remap coordinate arrays
    std::string effectiveSeq;
    std::vector<uint64_t> effectiveCoords;
    std::vector<uint64_t> effectiveBlockIds;
    if (hpc_) {
      auto [hpcSeq, hpcMapping] = seeding::hpcCompressWithMapping(localRangeSeq);
      effectiveSeq = std::move(hpcSeq);
      effectiveCoords.reserve(hpcMapping.size());
      effectiveBlockIds.reserve(hpcMapping.size());
      for (size_t idx : hpcMapping) {
        effectiveCoords.push_back(localRangeCoordToGlobalScalarCoords[idx]);
        effectiveBlockIds.push_back(localRangeCoordToBlockId[idx]);
      }
    } else {
      effectiveSeq = localRangeSeq;
      effectiveCoords = localRangeCoordToGlobalScalarCoords;
      effectiveBlockIds = localRangeCoordToBlockId;
    }
    
    if (effectiveSeq.size() >= indexBuilder.getK()) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(effectiveSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = effectiveCoords[startPos];
        auto endPosGlobal = effectiveCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = effectiveBlockIds[startPos];
        bool wasSeed = refOnSyncmers.contains(startPosGlobal);
        // Hard mask: skip all seed operations in flanked regions
        if (startPosGlobal < hardMaskStart || startPosGlobal > hardMaskEnd) continue;
        if (!wasSeed && isSeed) {
          auto it = refOnSyncmersMap.insert(startPosGlobal).first;
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::ADD, seeding::rsyncmer_t());
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::ADD);
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
          blockOnSyncmers[curBlockId].insert(startPosGlobal);
        } else if (wasSeed && !isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::DEL, refOnSyncmers.at(startPosGlobal));
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::DEL);
          refOnSyncmers.erase(startPosGlobal);
          refOnSyncmersMap.erase(startPosGlobal);
          blockOnSyncmers[curBlockId].erase(startPosGlobal);
          if (blockOnSyncmers[curBlockId].empty()) {
            blockOnSyncmers.erase(curBlockId);
          }
        } else if (wasSeed && isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::SUB, refOnSyncmers.at(startPosGlobal));
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }

    for (uint64_t pos : seedsToDelete) {
      // Hard mask: skip deletion in flanked regions
      if (pos < hardMaskStart || pos > hardMaskEnd) continue;
      if (!refOnSyncmers.contains(pos)) {
        output::error("refOnSyncmers[{}] is null", pos);
        std::exit(1);
      }
      refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers.at(pos));
      refOnSyncmers.erase(pos);
      refOnSyncmersMap.erase(pos);
      // blockOnSyncmers and blockOnSyncmersChangeRecord are updated in computeNewSyncmerRanges()
    }
  }

  // Handle potential syncmer deletions
  for (uint32_t pos : potentialSyncmerDeletions) {
    // Hard mask: skip deletion in flanked regions
    if (pos < hardMaskStart || pos > hardMaskEnd) continue;
    if (refOnSyncmers.contains(pos)) {
      refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers.at(pos));
      refOnSyncmers.erase(pos);
      const auto blockId = globalCoords.getBlockIdFromScalar(pos);
      blockOnSyncmers[blockId].erase(pos);
      if (blockOnSyncmers[blockId].empty()) blockOnSyncmers.erase(blockId);
      blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);
      refOnSyncmersMap.erase(pos);
    }
  }


  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      if (blockOnSyncmers.find(blockId) != blockOnSyncmers.end()) {
        std::vector<uint64_t> positionsToKeep;
        for (uint64_t pos : blockOnSyncmers[blockId]) {
          // Hard mask: skip deletion in flanked regions
          if (pos < hardMaskStart || pos > hardMaskEnd) {
            positionsToKeep.push_back(pos);
            continue;
          }
          refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers.at(pos));
          blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);

          refOnSyncmers.erase(pos);
          refOnSyncmersMap.erase(pos);
        }
        if (positionsToKeep.empty()) {
          blockOnSyncmers.erase(blockId);
        } else {
          blockOnSyncmers[blockId].clear();
          for (uint64_t pos : positionsToKeep) {
            blockOnSyncmers[blockId].insert(pos);
          }
        }
      }
    }
  }

  // processing k-min-mers
  std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> newKminmerRanges = computeNewKminmerRanges(refOnSyncmersChangeRecord, dfsIndex);

  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto beg = *newKminmerRanges[i].first;
    auto end = *newKminmerRanges[i].second;
    if (newKminmerRanges[i].second != refOnSyncmersMap.end() && beg > end) {
      output::error("beg ({}) > end ({}) in node {}", beg, end, node->identifier);
      std::exit(1);
    }
    if (i != 0 && *newKminmerRanges[i-1].second > beg) {
      output::error("newKminmerRanges[{}].second ({}) > newKminmerRanges[{}].first ({}) in node {}", 
                   i-1, *newKminmerRanges[i-1].second, i, *newKminmerRanges[i].first, node->identifier);
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
      startingSyncmerHashes.push_back(refOnSyncmers.at(*curIt).hash);
      if (j != l - 1) ++curIt;
    }
    if (startingSyncmerHashes.size() < l) {
      continue;
    }
    if (l == 1) {
      auto syncmerRev = refOnSyncmers.at(*curIt).isReverse;
      if (syncmerRev) {
        forwardHash = std::numeric_limits<size_t>::max();
        reverseHash = startingSyncmerHashes[0];
      } else {
        forwardHash = startingSyncmerHashes[0];
        reverseHash = std::numeric_limits<size_t>::max();
      }
      if (reverseHash == forwardHash) {
        output::error("syncmer hash collision detected in k-min-mer computation");
        exit(1);
      }
    } else {
      for (size_t j = 0; j < l; j++) {
        forwardHash = seeding::rol(forwardHash, k) ^ startingSyncmerHashes[j];
        reverseHash = seeding::rol(reverseHash, k) ^ startingSyncmerHashes[l - j - 1];
      }
    }

    uint64_t indexingPos = *indexingIt;
    uint64_t newHash = std::min(forwardHash, reverseHash);
    if (forwardHash != reverseHash) {
      bool substitution = refOnKminmers.contains(indexingPos);
      if (substitution) {
        uint64_t oldHash = refOnKminmers.at(indexingPos);
        refOnKminmersChangeRecord.emplace_back(indexingPos, panmapUtils::seedChangeType::SUB, oldHash);
        substitutedSeedHashes.emplace_back(oldHash, newHash);
      } else {
        refOnKminmersChangeRecord.emplace_back(indexingPos, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
        addedSeedHashes.push_back(newHash);
      }
      // Store hash directly
      refOnKminmers[indexingPos] = newHash;
    } else {
      if (refOnKminmers.contains(indexingPos)) {
        uint64_t oldHash = refOnKminmers.at(indexingPos);
        refOnKminmersChangeRecord.emplace_back(indexingPos, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        refOnKminmers.erase(indexingPos);
      }
    }

    if (curIt == endIt) continue;
    
    while (curIt != endIt) {
      ++curIt;
      if (curIt == refOnSyncmersMap.end()) break;
      forwardHash = seeding::rol(forwardHash, k) ^ seeding::rol(refOnSyncmers.at(*indexingIt).hash, k * l) ^ refOnSyncmers.at(*curIt).hash;
      reverseHash = seeding::ror(reverseHash, k) ^ seeding::ror(refOnSyncmers.at(*indexingIt).hash, k)     ^ seeding::rol(refOnSyncmers.at(*curIt).hash, k * (l-1));
      ++indexingIt;

      uint64_t indexingPos2 = *indexingIt;
      uint64_t newHash2 = std::min(forwardHash, reverseHash);
      if (forwardHash != reverseHash) {
        bool substitution = refOnKminmers.contains(indexingPos2);
        if (substitution) {
          uint64_t oldHash = refOnKminmers.at(indexingPos2);
          refOnKminmersChangeRecord.emplace_back(indexingPos2, panmapUtils::seedChangeType::SUB, oldHash);
          substitutedSeedHashes.emplace_back(oldHash, newHash2);
        } else {
          refOnKminmersChangeRecord.emplace_back(indexingPos2, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
          addedSeedHashes.push_back(newHash2);
        }
        refOnKminmers[indexingPos2] = newHash2;
      } else {
        if (refOnKminmers.contains(indexingPos2)) {
          uint64_t oldHash = refOnKminmers.at(indexingPos2);
          refOnKminmersChangeRecord.emplace_back(indexingPos2, panmapUtils::seedChangeType::DEL, oldHash);
          deletedSeedHashes.push_back(oldHash);
          refOnKminmers.erase(indexingPos2);
        }
      }
    }
  }
  if (!newKminmerRanges.empty() && newKminmerRanges.back().second == refOnSyncmersMap.end()) {
    auto delIt = newKminmerRanges.back().second;
    for (size_t j = 0; j < l - 1; j++) {
      --delIt;
      uint64_t delPos = *delIt;
      if (refOnKminmers.contains(delPos)) {
        uint64_t oldHash = refOnKminmers.at(delPos);
        refOnKminmersChangeRecord.emplace_back(delPos, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        refOnKminmers.erase(delPos);
      }
      if (delIt == refOnSyncmersMap.begin()) break;
    }
  }
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL && refOnKminmers.contains(syncmerPos)) {
      uint64_t oldHash = refOnKminmers.at(syncmerPos);
      refOnKminmersChangeRecord.emplace_back(syncmerPos, panmapUtils::seedChangeType::DEL, oldHash);
      deletedSeedHashes.push_back(oldHash);
      refOnKminmers.erase(syncmerPos);
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
  
  gap_map::revertGapMapChanges(gapRunBlockInversionBacktracks, gapMap);
  std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>().swap(gapRunBlockInversionBacktracks); // gapRunBlockInversionBacktracks is no longer needed... clear memory

  // update delayed block states
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExistsDelayed[blockId] = newExists;
    blockStrandDelayed[blockId] = newStrand;
  }

  nodeToDfsIndex[node->identifier] = dfsIndex;
  output::progress(fmt::format("dfsIndex: {}", dfsIndex));
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
      refOnSyncmers.erase(pos);
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

  for (const auto& [pos, changeType, kminmerHash] : refOnKminmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::ADD) {
      refOnKminmers.erase(pos);
    } else {
      refOnKminmers[pos] = kminmerHash;
    }
  }
}


void index_single_mode::IndexBuilder::buildIndex() {
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences);
  // Note: refOnSyncmers and refOnKminmers are now sparse maps, no resize needed

  std::vector<char> blockExistsDelayed = blockSequences.blockExists;
  std::vector<char> blockStrandDelayed = blockSequences.blockStrand;

  std::map<uint64_t, uint64_t> gapMap{{0, globalCoords.lastScalarCoord}};
  std::unordered_set<uint64_t> invertedBlocks;

  // add lite tree to index
  LiteTree::Builder liteTreeBuilder = indexBuilder.initLiteTree();
  std::unordered_set<std::string_view> emptyNodes;
  uint64_t dfsIndex = 0;
  buildIndexHelper(T->root, emptyNodes, blockSequences, blockExistsDelayed, blockStrandDelayed, globalCoords, gapMap, invertedBlocks, dfsIndex);
  output::progress_clear();

  size_t numNodes = T->allNodes.size();

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

  output::step("Computing seed deltas...");

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
      output::progress(fmt::format("Processed node {}/{}", nodeIdx, numNodes));
    }
  }
  output::progress_clear();

  // Cap'n Proto list element limit (~536M for UInt64). Split into primary + overflow.
  constexpr size_t CAPNP_SPLIT = 500'000'000;
  size_t numSegments = (totalChanges + CAPNP_SPLIT - 1) / CAPNP_SPLIT;

  // Allocate flat arrays
  auto seedChangeHashesBuilder = indexBuilder.initSeedChangeHashes(numSegments);
  auto seedChangeParentCountsBuilder = indexBuilder.initSeedChangeParentCounts(numSegments);
  auto seedChangeChildCountsBuilder = indexBuilder.initSeedChangeChildCounts(numSegments);
  auto nodeChangeOffsetsBuilder = indexBuilder.initNodeChangeOffsets(numNodes + 1);

  for (size_t seg = 0; seg < numSegments; seg++) {
    auto segStart = seg * CAPNP_SPLIT;
    auto segEnd = std::min((uint64_t)(segStart + CAPNP_SPLIT), totalChanges);
    auto segSize = segEnd - segStart;
    seedChangeHashesBuilder.init(seg, segSize);
    seedChangeParentCountsBuilder.init(seg, segSize);
    seedChangeChildCountsBuilder.init(seg, segSize);
  }
  
  if (numSegments > 1) {
    output::step("Large index: splitting {} seed changes into {} segments", totalChanges, numSegments);
  }

  // Sort each node's changes by hash for better compression
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    std::sort(nodeChanges[nodeIdx].begin(), nodeChanges[nodeIdx].end(),
             [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
  }

  // Write flat arrays
  uint64_t offset = 0;
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    nodeChangeOffsetsBuilder.set(nodeIdx, offset);
    for (const auto& [hash, parentCount, childCount] : nodeChanges[nodeIdx]) {
      if (parentCount < INT16_MIN || parentCount > INT16_MAX || childCount < INT16_MIN || childCount > INT16_MAX) {
        output::error("Seed count overflow at node {}: parentCount={}, childCount={} (Int16 range: [{}, {}]). "
                      "This genome has too many repeated k-min-mers for the current index format.",
                      nodeIdx, parentCount, childCount, INT16_MIN, INT16_MAX);
        std::exit(1);
      }
      auto seg = offset / CAPNP_SPLIT;
      auto segOffset = offset % CAPNP_SPLIT;
      seedChangeHashesBuilder[seg].set(segOffset, hash);
      seedChangeParentCountsBuilder[seg].set(segOffset, static_cast<int16_t>(parentCount));
      seedChangeChildCountsBuilder[seg].set(segOffset, static_cast<int16_t>(childCount));
      offset++;
    }
  }
  nodeChangeOffsetsBuilder.set(numNodes, offset);

  output::done(fmt::format("Index built ({} seed changes)", totalChanges));
}

// ============================================================================
// Substitution Spectrum Computation
// ============================================================================

static int nucToIdx(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}

void index_single_mode::IndexBuilder::computeSubstitutionSpectrum() {
  output::step("Computing substitution spectrum from tree...");

  // 4x4 counts: subCounts[from][to]
  std::vector<std::vector<int64_t>> subCounts(4, std::vector<int64_t>(4, 0));
  int64_t numBranches = 0;

  // Set up block sequences from root
  panmapUtils::BlockSequences blockSequences(T);
  auto& sequence = blockSequences.sequence;

  // Simple DFS: apply block+nuc mutations, count substitutions, recurse, undo
  std::function<void(panmanUtils::Node*)> dfs = [&](panmanUtils::Node* node) {
    std::vector<std::pair<panmapUtils::Coordinate, char>> nucUndoRecord;
    std::vector<std::pair<uint32_t, bool>> blockUndoRecord;

    // Apply block mutations
    for (const auto& blockMut : node->blockMutation) {
      uint32_t blockId = blockMut.primaryBlockId;
      bool oldExists = blockSequences.blockExists[blockId];
      bool isInsertion = blockMut.blockMutInfo;
      blockUndoRecord.emplace_back(blockId, oldExists);
      blockSequences.blockExists[blockId] = isInsertion;
    }

    if (node != T->root) {
      numBranches++;

      for (const auto& nucMutation : node->nucMutation) {
        int length = nucMutation.mutInfo >> 4;
        uint32_t type = nucMutation.mutInfo & 0x7;
        bool isSub = (type == panmanUtils::NucMutationType::NS ||
                      type == panmanUtils::NucMutationType::NSNPS);

        for (int i = 0; i < length; i++) {
          panmapUtils::Coordinate pos(nucMutation, i);
          if (pos.nucPosition >= sequence[pos.primaryBlockId].size()) continue;

          char oldNuc = blockSequences.getSequenceBase(pos);
          int newNucCode = (nucMutation.nucs >> (4 * (5 - i))) & 0xF;
          char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
          nucUndoRecord.emplace_back(pos, oldNuc);
          blockSequences.setSequenceBase(pos, newNuc);

          if (isSub && blockSequences.blockExists[pos.primaryBlockId]) {
            int oldIdx = nucToIdx(oldNuc);
            int newIdx = nucToIdx(newNuc);
            if (oldIdx >= 0 && newIdx >= 0 && oldIdx != newIdx) {
              subCounts[oldIdx][newIdx]++;
            }
          }
        }
      }
    }

    for (auto* child : node->children) {
      dfs(child);
    }

    // Backtrack
    for (auto it = nucUndoRecord.rbegin(); it != nucUndoRecord.rend(); ++it) {
      blockSequences.setSequenceBase(it->first, it->second);
    }
    for (auto it = blockUndoRecord.rbegin(); it != blockUndoRecord.rend(); ++it) {
      blockSequences.blockExists[it->first] = it->second;
    }
  };

  dfs(T->root);

  // Estimate genome length from sampled leaf nodes using existing functions
  int64_t genomeLen = 0;
  {
    std::vector<int64_t> lengths;
    // Collect leaf node identifiers
    std::vector<std::string> leafIds;
    for (const auto& [id, node] : T->allNodes) {
      if (node->children.empty()) leafIds.push_back(id);
    }
    // Sample up to 10 leaves evenly spaced
    size_t sampleCount = std::min<size_t>(10, leafIds.size());
    size_t step = std::max<size_t>(1, leafIds.size() / sampleCount);
    for (size_t i = 0; i < leafIds.size() && lengths.size() < sampleCount; i += step) {
      std::vector<std::vector<std::pair<char, std::vector<char>>>> seq;
      std::vector<char> bExists, bStrand;
      std::unordered_map<int, int> bLengths;
      panmapUtils::getSequenceFromReference(T, seq, bExists, bStrand, bLengths, leafIds[i]);
      std::string ungapped = panmapUtils::getStringFromSequence(seq, bLengths, bExists, bStrand, false);
      lengths.push_back((int64_t)ungapped.size());
    }
    if (!lengths.empty()) {
      std::sort(lengths.begin(), lengths.end());
      genomeLen = lengths[lengths.size() / 2]; // median
    }
  }

  // Count base composition from the median-length sample (approximate with uniform A/C/G/T)
  // For rate normalization we just need per-base counts; use genome/4 as fallback
  std::vector<int64_t> baseCounts(4, genomeLen / 4);

  auto substMatBuilder = indexBuilder.initSubstitutionMatrix(16);
  
  if (numBranches > 0 && genomeLen > 0) {
    int64_t totalSubs = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        if (i != j) totalSubs += subCounts[i][j];

    output::step("Substitution spectrum: {} total substitutions across {} branches ({} bp genome)",
                 totalSubs, numBranches, genomeLen);

    for (int from = 0; from < 4; from++) {
      int64_t baseCount = baseCounts[from];
      for (int to = 0; to < 4; to++) {
        double rate;
        if (from == to) {
          double offDiagSum = 0;
          for (int j = 0; j < 4; j++) {
            if (j != from && baseCount > 0) {
              offDiagSum += (double)subCounts[from][j] / (numBranches * baseCount);
            }
          }
          rate = 1.0 - offDiagSum;
        } else {
          rate = (baseCount > 0) ? (double)subCounts[from][to] / (numBranches * baseCount) : 0.0;
        }
        substMatBuilder.set(from * 4 + to, rate);
      }
    }

    // Log the matrix
    const char* bases = "ACGT";
    for (int from = 0; from < 4; from++) {
      std::string row;
      for (int to = 0; to < 4; to++) {
        if (!row.empty()) row += "  ";
        row += fmt::format("{}→{}: {:.6e}", bases[from], bases[to], substMatBuilder[from * 4 + to]);
      }
      output::debug("{}", row);
    }
  } else {
    // No data: set identity matrix 
    for (int i = 0; i < 16; i++) {
      substMatBuilder.set(i, (i / 4 == i % 4) ? 1.0 : 0.0);
    }
    output::step("No branches in tree, using identity substitution matrix");
  }

  output::done("Substitution spectrum computed");
}

void index_single_mode::IndexBuilder::writeIndex(const std::string& path, int numThreads, int zstdLevel) {
  output::step("Serializing index...");
  
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
  
  output::step("Writing index ({} bytes uncompressed)...", dataSize);
  
  // Use ZSTD compression to write
  if (!panmap_zstd::compressToFile(data, dataSize, path, zstdLevel, numThreads)) {
    output::error("failed to write compressed index to {}", path);
    std::exit(1);
  }
  
  output::done("Index written to " + path);
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
  BacktrackInfo* backtrackInfo,
  bool skipNodeChanges
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
  std::vector<std::pair<uint64_t, int64_t>> localRunningCountChanges;

  // References to use
  auto& blockMutationRecord = backtrackInfo ? backtrackInfo->blockMutationRecord : localBlockMutationRecord;
  auto& nucMutationRecord = backtrackInfo ? backtrackInfo->nucMutationRecord : localNucMutationRecord;
  auto& gapRunBacktracks = backtrackInfo ? backtrackInfo->gapRunBacktracks : localGapRunBacktracks;
  auto& invertedBlocksBacktracks = backtrackInfo ? backtrackInfo->invertedBlocksBacktracks : localInvertedBlocksBacktracks;
  auto& refOnSyncmersChangeRecord = backtrackInfo ? backtrackInfo->refOnSyncmersChangeRecord : localRefOnSyncmersChangeRecord;
  auto& blockOnSyncmersChangeRecord = backtrackInfo ? backtrackInfo->blockOnSyncmersChangeRecord : localBlockOnSyncmersChangeRecord;
  auto& refOnKminmersChangeRecord = backtrackInfo ? backtrackInfo->refOnKminmersChangeRecord : localRefOnKminmersChangeRecord;
  auto& gapRunBlockInversionBacktracks = backtrackInfo ? backtrackInfo->gapRunBlockInversionBacktracks : localGapRunBlockInversionBacktracks;
  auto& runningCountChanges = backtrackInfo ? backtrackInfo->runningCountChanges : localRunningCountChanges;

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
  localMutationRanges.reserve(std::max(size_t(32), node->blockMutation.size() + node->nucMutation.size() * 6));
  gapRunUpdates.reserve(std::max(size_t(64), node->nucMutation.size() * 6 + node->blockMutation.size() * 10));
  // Reserve for change records - these can be large
  refOnSyncmersChangeRecord.reserve(128);  // avg is 76-153
  blockOnSyncmersChangeRecord.reserve(128);  // avg is 66-112
  
  panmapUtils::applyMutations(node, dfsIndex, state.blockSequences, state.invertedBlocks, globalCoords, 
                 localMutationRanges, blockMutationRecord, nucMutationRecord, gapRunUpdates, 
                 invertedBlocksBacktracks, potentialSyncmerDeletions, 
                 state.blockExistsDelayed, state.blockStrandDelayed, imputeAmb_);
  
  blockMutationRecord.shrink_to_fit();
  nucMutationRecord.shrink_to_fit();
  potentialSyncmerDeletions.shrink_to_fit();
  localMutationRanges.shrink_to_fit();
  gapRunUpdates.shrink_to_fit();

  std::sort(gapRunUpdates.begin(), gapRunUpdates.end(), 
            [&](const auto& a, const auto& b) { return a.second.first < b.second.first; });
  panmapUtils::updateGapMap(state.gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  std::vector<uint64_t> invertedBlocksVec(state.invertedBlocks.begin(), state.invertedBlocks.end());
  std::sort(invertedBlocksVec.begin(), invertedBlocksVec.end());
  for (const auto& blockId : invertedBlocksVec) {
    uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
    uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
    gap_map::invertGapMap(state.gapMap, {beg, end}, gapRunBlockInversionBacktracks, gapMapUpdates);
  }

  // Update genome extent from gapMap - used to determine flank regions
  // Save previous extent for backtracking
  if (backtrackInfo) {
    backtrackInfo->prevFirstNonGapScalar = state.firstNonGapScalar;
    backtrackInfo->prevLastNonGapScalar = state.lastNonGapScalar;
  }

  if (extentGuard_) {
    // Pass flankSize=0 here - flank masking is applied separately via hardMaskStart/hardMaskEnd below
    auto [newFirstNonGap, newLastNonGap] = panmapUtils::computeExtentFromGapMap(state.gapMap, globalCoords.lastScalarCoord, 0);
    
    // Debug: log extent changes for nodes we care about
    static bool debugFlankMasking = (std::getenv("DEBUG_FLANK_MASKING") != nullptr);
    if (debugFlankMasking) {
      if (state.firstNonGapScalar != newFirstNonGap || state.lastNonGapScalar != newLastNonGap) {
        output::debug("[FLANK] Node {} (dfs={}): extent [{},{}] -> [{},{}]",
                     node->identifier, dfsIndex, 
                     state.firstNonGapScalar, state.lastNonGapScalar,
                     newFirstNonGap, newLastNonGap);
      }
    }
    
    state.firstNonGapScalar = newFirstNonGap;
    state.lastNonGapScalar = newLastNonGap;
  }
  
  // Compute hard mask boundaries relative to each genome's extent
  // Seeds in hard-masked regions are completely ignored (no adds, no deletes)
  const int flankMaskBp = getFlankMaskBp();
  
  // Use computeExtentFromGapMap with flankSize to correctly skip non-gap bases
  // (not scalar positions) when determining mask boundaries
  uint64_t hardMaskStart = state.firstNonGapScalar;
  uint64_t hardMaskEnd = state.lastNonGapScalar;
  if (flankMaskBp > 0) {
    auto [hms, hme] = panmapUtils::computeExtentFromGapMap(state.gapMap, globalCoords.lastScalarCoord, flankMaskBp);
    hardMaskStart = hms;
    hardMaskEnd = hme;
  }
  
  // Debug logging for specific nodes
  static bool debugMasking = (std::getenv("DEBUG_MASK") != nullptr);
  if (debugMasking && (node->identifier.find("node_9981") != std::string::npos || 
                        node->identifier.find("node_11963") != std::string::npos ||
                        node->identifier.find("node_11964") != std::string::npos)) {
    output::debug("[MASK-DEBUG] {} extent=[{},{}] hardMask=[{},{}] flankMaskBp={}",
                 node->identifier, state.firstNonGapScalar, state.lastNonGapScalar,
                 hardMaskStart, hardMaskEnd, flankMaskBp);
  }

  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = computeNewSyncmerRangesJump(
      node, dfsIndex, state.blockSequences, state.blockExistsDelayed, state.blockStrandDelayed, 
      globalCoords, state.gapMap, localMutationRanges, blockOnSyncmersChangeRecord, 
      state.refOnSyncmers, state.blockOnSyncmers, state.firstNonGapScalar, state.lastNonGapScalar);

  // Processing syncmers with dual masking:
  // 1. HARD MASK: First/last flankMaskBp positions - completely ignore all seeds
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    
    // Apply HPC if enabled: compress localRangeSeq and remap coordinate arrays
    std::string effectiveSeq;
    std::vector<uint64_t> effectiveCoords;
    std::vector<uint64_t> effectiveBlockIds;
    if (hpc_) {
      auto [hpcSeq, hpcMapping] = seeding::hpcCompressWithMapping(localRangeSeq);
      effectiveSeq = std::move(hpcSeq);
      effectiveCoords.reserve(hpcMapping.size());
      effectiveBlockIds.reserve(hpcMapping.size());
      for (size_t idx : hpcMapping) {
        effectiveCoords.push_back(localRangeCoordToGlobalScalarCoords[idx]);
        effectiveBlockIds.push_back(localRangeCoordToBlockId[idx]);
      }
    } else {
      effectiveSeq = localRangeSeq;
      effectiveCoords = localRangeCoordToGlobalScalarCoords;
      effectiveBlockIds = localRangeCoordToBlockId;
    }
    
    if (effectiveSeq.size() >= static_cast<size_t>(indexBuilder.getK())) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(effectiveSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = effectiveCoords[startPos];
        auto endPosGlobal = effectiveCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = effectiveBlockIds[startPos];
        bool wasSeed = state.refOnSyncmers.contains(startPosGlobal);
        
        // Check masking conditions:
        // 1. Hard mask: position < hardMaskStart OR position > hardMaskEnd
        bool isHardMasked = (startPosGlobal < hardMaskStart || startPosGlobal > hardMaskEnd);
        
        // Hard-masked: skip ALL seed operations (adds and deletes)
        if (isHardMasked) {
          continue;
        }
        
        if (!wasSeed && isSeed) {
          state.refOnSyncmersMap.insert(startPosGlobal);
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::ADD, seeding::rsyncmer_t());
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::ADD);
          state.refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
          state.blockOnSyncmers[curBlockId].insert(startPosGlobal);
        } else if (wasSeed && !isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::DEL, state.refOnSyncmers.at(startPosGlobal));
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::DEL);
          state.refOnSyncmers.erase(startPosGlobal);
          state.refOnSyncmersMap.erase(startPosGlobal);
          state.blockOnSyncmers[curBlockId].erase(startPosGlobal);
          if (state.blockOnSyncmers[curBlockId].empty()) {
            state.blockOnSyncmers.erase(curBlockId);
          }
        } else if (wasSeed && isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::SUB, state.refOnSyncmers.at(startPosGlobal));
          state.refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }

    for (uint64_t pos : seedsToDelete) {
      // Check masking conditions
      bool isHardMasked = (pos < hardMaskStart || pos > hardMaskEnd);
      
      // Hard-masked: skip deletion
      if (isHardMasked) {
        continue;
      }
      
      if (!state.refOnSyncmers.contains(pos)) {
        output::error("refOnSyncmers[{}] is null", pos);
        std::exit(1);
      }
      refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, state.refOnSyncmers.at(pos));
      state.refOnSyncmers.erase(pos);
      state.refOnSyncmersMap.erase(pos);
    }
  }

  // Handle potential syncmer deletions
  for (uint32_t pos : potentialSyncmerDeletions) {
    if (state.refOnSyncmers.contains(pos)) {
      // Check masking conditions
      bool isHardMasked = (pos < hardMaskStart || pos > hardMaskEnd);
      
      if (isHardMasked) {
        continue;
      }
      refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, state.refOnSyncmers.at(pos));
      state.refOnSyncmers.erase(pos);
      const auto blockId = globalCoords.getBlockIdFromScalar(pos);
      state.blockOnSyncmers[blockId].erase(pos);
      if (state.blockOnSyncmers[blockId].empty()) state.blockOnSyncmers.erase(blockId);
      blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);
      state.refOnSyncmersMap.erase(pos);
    }
  }

  // Handle block mutations that delete entire blocks (matches original)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      if (state.blockOnSyncmers.find(blockId) != state.blockOnSyncmers.end()) {
        std::vector<uint64_t> positionsToKeep;
        for (uint64_t pos : state.blockOnSyncmers[blockId]) {
          // Check masking conditions
          bool isHardMasked = (pos < hardMaskStart || pos > hardMaskEnd);
          
          if (isHardMasked) {
            positionsToKeep.push_back(pos);
            continue;
          }
          refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, state.refOnSyncmers.at(pos));
          blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);
          state.refOnSyncmers.erase(pos);
          state.refOnSyncmersMap.erase(pos);
        }
        if (positionsToKeep.empty()) {
          state.blockOnSyncmers.erase(blockId);
        } else {
          // Keep only the masked/flank positions
          state.blockOnSyncmers[blockId].clear();
          for (uint64_t pos : positionsToKeep) {
            state.blockOnSyncmers[blockId].insert(pos);
          }
        }
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

  std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> newKminmerRanges = 
      computeNewKminmerRanges(refOnSyncmersChangeRecord, state, dfsIndex);

  for (size_t i = 0; i < newKminmerRanges.size(); i++) {
    auto beg = *newKminmerRanges[i].first;
    auto end = *newKminmerRanges[i].second;
    if (newKminmerRanges[i].second != state.refOnSyncmersMap.end() && beg > end) {
      output::error("beg ({}) > end ({}) in node {}", beg, end, node->identifier);
      std::exit(1);
    }
    if (i != 0 && *newKminmerRanges[i-1].second > beg) {
      output::error("newKminmerRanges[{}].second > newKminmerRanges[{}].first", i-1, i);
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
      startingSyncmerHashes.push_back(state.refOnSyncmers.at(*curIt).hash);
      if (j != static_cast<size_t>(l) - 1) ++curIt;
    }
    if (startingSyncmerHashes.size() < static_cast<size_t>(l)) {
      continue;
    }
    if (l == 1) {
      auto syncmerRev = state.refOnSyncmers.at(*curIt).isReverse;
      if (syncmerRev) {
        forwardHash = std::numeric_limits<size_t>::max();
        reverseHash = startingSyncmerHashes[0];
      } else {
        forwardHash = startingSyncmerHashes[0];
        reverseHash = std::numeric_limits<size_t>::max();
      }
      if (reverseHash == forwardHash) {
        output::error("syncmer hash collision detected in k-min-mer computation");
        exit(1);
      }
    } else {
      for (size_t j = 0; j < static_cast<size_t>(l); j++) {
        forwardHash = seeding::rol(forwardHash, k) ^ startingSyncmerHashes[j];
        reverseHash = seeding::rol(reverseHash, k) ^ startingSyncmerHashes[l - j - 1];
      }
    }

    uint64_t indexingPos = *indexingIt;
    uint64_t newHash = std::min(forwardHash, reverseHash);
    if (forwardHash != reverseHash) {
      auto existingIt = state.refOnKminmers.find(indexingPos);
      bool substitution = (existingIt != state.refOnKminmers.end());
      if (substitution) {
        // Get old hash from stored value and add to substitution list
        uint64_t oldHash = existingIt->second;
        refOnKminmersChangeRecord.emplace_back(indexingPos, panmapUtils::seedChangeType::SUB, oldHash);
        substitutedSeedHashes.emplace_back(oldHash, newHash);
      } else {
        refOnKminmersChangeRecord.emplace_back(indexingPos, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
        addedSeedHashes.push_back(newHash);
      }
      // Store hash directly instead of index
      state.refOnKminmers[indexingPos] = newHash;
    } else {
      auto existingIt = state.refOnKminmers.find(indexingPos);
      if (existingIt != state.refOnKminmers.end()) {
        uint64_t oldHash = existingIt->second;
        refOnKminmersChangeRecord.emplace_back(indexingPos, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        state.refOnKminmers.erase(existingIt);
      }
    }

    if (curIt == endIt) continue;
    
    while (curIt != endIt) {
      ++curIt;
      if (curIt == state.refOnSyncmersMap.end()) break;
      forwardHash = seeding::rol(forwardHash, k) ^ seeding::rol(state.refOnSyncmers.at(*indexingIt).hash, k * l) ^ state.refOnSyncmers.at(*curIt).hash;
      reverseHash = seeding::ror(reverseHash, k) ^ seeding::ror(state.refOnSyncmers.at(*indexingIt).hash, k)     ^ seeding::rol(state.refOnSyncmers.at(*curIt).hash, k * (l-1));
      ++indexingIt;

      uint64_t indexingPos2 = *indexingIt;
      uint64_t newHash2 = std::min(forwardHash, reverseHash);
      if (forwardHash != reverseHash) {
        auto existingIt2 = state.refOnKminmers.find(indexingPos2);
        bool substitution = (existingIt2 != state.refOnKminmers.end());
        if (substitution) {
          uint64_t oldHash = existingIt2->second;
          refOnKminmersChangeRecord.emplace_back(indexingPos2, panmapUtils::seedChangeType::SUB, oldHash);
          substitutedSeedHashes.emplace_back(oldHash, newHash2);
        } else {
          refOnKminmersChangeRecord.emplace_back(indexingPos2, panmapUtils::seedChangeType::ADD, std::numeric_limits<uint64_t>::max());
          addedSeedHashes.push_back(newHash2);
        }
        // Store hash directly
        state.refOnKminmers[indexingPos2] = newHash2;
      } else {
        auto existingIt2 = state.refOnKminmers.find(indexingPos2);
        if (existingIt2 != state.refOnKminmers.end()) {
          uint64_t oldHash = existingIt2->second;
          refOnKminmersChangeRecord.emplace_back(indexingPos2, panmapUtils::seedChangeType::DEL, oldHash);
          deletedSeedHashes.push_back(oldHash);
          state.refOnKminmers.erase(existingIt2);
        }
      }
    }
  }

  // Handle end-of-range k-minmer deletions (matches original)
  if (!newKminmerRanges.empty() && newKminmerRanges.back().second == state.refOnSyncmersMap.end()) {
    auto delIt = newKminmerRanges.back().second;
    for (size_t j = 0; j < static_cast<size_t>(l) - 1; j++) {
      --delIt;
      auto existingIt = state.refOnKminmers.find(*delIt);
      if (existingIt != state.refOnKminmers.end()) {
        uint64_t oldHash = existingIt->second;
        refOnKminmersChangeRecord.emplace_back(*delIt, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        state.refOnKminmers.erase(existingIt);
      }
      if (delIt == state.refOnSyncmersMap.begin()) break;
    }
  }

  // Handle deleted syncmers that still have k-minmers (matches original)
  for (const auto& [syncmerPos, changeType, rsyncmer] : refOnSyncmersChangeRecord) {
    if (changeType == panmapUtils::seedChangeType::DEL) {
      auto existingIt = state.refOnKminmers.find(syncmerPos);
      if (existingIt != state.refOnKminmers.end()) {
        uint64_t oldHash = existingIt->second;
        refOnKminmersChangeRecord.emplace_back(syncmerPos, panmapUtils::seedChangeType::DEL, oldHash);
        deletedSeedHashes.push_back(oldHash);
        state.refOnKminmers.erase(existingIt);
      }
    }
  }

  // No sorting needed - we're using hashes directly, not positions

  // Compute node changes directly using running counts
  // This combines what was previously a separate post-traversal pass
  // Skip if we're just walking a path to set up state (path nodes already processed)
  if (!skipNodeChanges) {
    auto& changes = (*state.nodeChanges)[dfsIndex];
    
    // Collect all hashes that are modified at this node
    std::unordered_set<uint64_t> modifiedHashes;
    modifiedHashes.reserve(addedSeedHashes.size() + deletedSeedHashes.size() + 
                           substitutedSeedHashes.size() * 2);
    
    for (uint64_t hash : addedSeedHashes) modifiedHashes.insert(hash);
    for (uint64_t hash : deletedSeedHashes) modifiedHashes.insert(hash);
    for (const auto& [oldHash, newHash] : substitutedSeedHashes) {
      modifiedHashes.insert(oldHash);
      modifiedHashes.insert(newHash);
    }
    
    // Record parent counts BEFORE applying changes
    std::unordered_map<uint64_t, int64_t> parentCounts;
    parentCounts.reserve(modifiedHashes.size());
    for (uint64_t hash : modifiedHashes) {
      auto it = state.runningCounts.find(hash);
      parentCounts[hash] = (it != state.runningCounts.end()) ? it->second : 0;
    }
    
    // Apply changes to running counts and record for backtracking
    for (uint64_t hash : addedSeedHashes) {
      state.runningCounts[hash]++;
      runningCountChanges.emplace_back(hash, +1);
    }
    for (uint64_t hash : deletedSeedHashes) {
      auto it = state.runningCounts.find(hash);
      if (it != state.runningCounts.end()) {
        it->second--;
        if (it->second <= 0) state.runningCounts.erase(it);
      }
      runningCountChanges.emplace_back(hash, -1);
    }
    for (const auto& [oldHash, newHash] : substitutedSeedHashes) {
      // Delete old
      auto oldIt = state.runningCounts.find(oldHash);
      if (oldIt != state.runningCounts.end()) {
        oldIt->second--;
        if (oldIt->second <= 0) state.runningCounts.erase(oldIt);
      }
      runningCountChanges.emplace_back(oldHash, -1);
      // Add new
      state.runningCounts[newHash]++;
      runningCountChanges.emplace_back(newHash, +1);
    }
    
    // Compare parent vs child counts to compute changes
    for (uint64_t hash : modifiedHashes) {
      int64_t parentCount = parentCounts[hash];
      int64_t childCount = 0;
      auto it = state.runningCounts.find(hash);
      if (it != state.runningCounts.end()) childCount = it->second;
      
      if (parentCount != childCount) {
        changes.emplace_back(hash, parentCount, childCount);
      }
    }
  } else {
    // Even when skipping nodeChanges, we still need to update runningCounts
    // for correct state during subtree traversal
    for (uint64_t hash : addedSeedHashes) {
      state.runningCounts[hash]++;
      runningCountChanges.emplace_back(hash, +1);
    }
    for (uint64_t hash : deletedSeedHashes) {
      auto it = state.runningCounts.find(hash);
      if (it != state.runningCounts.end()) {
        it->second--;
        if (it->second <= 0) state.runningCounts.erase(it);
      }
      runningCountChanges.emplace_back(hash, -1);
    }
    for (const auto& [oldHash, newHash] : substitutedSeedHashes) {
      auto oldIt = state.runningCounts.find(oldHash);
      if (oldIt != state.runningCounts.end()) {
        oldIt->second--;
        if (oldIt->second <= 0) state.runningCounts.erase(oldIt);
      }
      runningCountChanges.emplace_back(oldHash, -1);
      state.runningCounts[newHash]++;
      runningCountChanges.emplace_back(newHash, +1);
    }
  }

  // Revert gap map inversions (for proper state for children)
  gap_map::revertGapMapChanges(gapRunBlockInversionBacktracks, state.gapMap);

  // Update delayed block states (for children)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    state.blockExistsDelayed[blockId] = newExists;
    state.blockStrandDelayed[blockId] = newStrand;
  }

  // nodeToDfsIndex is pre-computed, no need to update here

  // Update progress with instrumentation stats
  // Only count nodes where we actually computed nodeChanges (not path traversal duplicates)
  if (skipNodeChanges) return;
  
  uint64_t processed = ++processedNodes_;
  uint64_t totalNodes = T->allNodes.size();
  
  // Log progress every 1000 nodes or every 2 seconds, whichever comes first
  auto now = std::chrono::steady_clock::now();
  bool shouldLog = (processed % 1000 == 0) || (processed == totalNodes);
  if (!shouldLog && processed > 100) {
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastProgressLog_).count();
    shouldLog = (elapsed >= 2000);
  }
  
  if (shouldLog) {
    auto elapsedSinceStart = std::chrono::duration_cast<std::chrono::milliseconds>(now - buildStartTime_).count();
    double nodesPerSec = (elapsedSinceStart > 0) ? (processed * 1000.0 / elapsedSinceStart) : 0;
    
    lastProgressLog_ = now;
    uint64_t totalClones = totalClonesCreated_.load(std::memory_order_relaxed);
    
    double pct = 100.0 * processed / totalNodes;
    output::progress(fmt::format("[{:.1f}%] {}/{} ({:.0f} n/s) chunks:{}", 
                                 pct, processed, totalNodes, nodesPerSec, totalClones));
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
      state.refOnSyncmers.erase(pos);
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
      state.refOnKminmers.erase(pos);
    } else {
      state.refOnKminmers[pos] = kminmerHash;
    }
  }
  
  // Revert running count changes (in reverse order)
  for (auto it = backtrackInfo.runningCountChanges.rbegin(); 
       it != backtrackInfo.runningCountChanges.rend(); ++it) {
    const auto& [hash, delta] = *it;
    // Undo by applying the opposite delta
    if (delta > 0) {
      // Was added, need to decrement
      auto countIt = state.runningCounts.find(hash);
      if (countIt != state.runningCounts.end()) {
        countIt->second--;
        if (countIt->second <= 0) state.runningCounts.erase(countIt);
      }
    } else {
      // Was deleted, need to increment
      state.runningCounts[hash]++;
    }
  }
  
  // Revert genome extent
  state.firstNonGapScalar = backtrackInfo.prevFirstNonGapScalar;
  state.lastNonGapScalar = backtrackInfo.prevLastNonGapScalar;
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
  size_t parallelThreshold,  // Only parallelize if combined subtree size >= this
  size_t depth  // Current recursion depth for instrumentation
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
    // 3. At least one child has a reasonably sized subtree (worth spawning a task for)
    //    AND the remaining subtree is also worth processing in parallel
    bool shouldParallelize = false;
    if (numChildren > 1 && totalChildSubtreeSize >= parallelThreshold) {
      // Find the two largest subtrees
      size_t largest = 0, secondLargest = 0;
      for (auto* child : node->children) {
        size_t sz = subtreeSizes_[child->identifier];
        if (sz > largest) {
          secondLargest = largest;
          largest = sz;
        } else if (sz > secondLargest) {
          secondLargest = sz;
        }
      }
      // Parallelize if the second largest subtree is worth spawning a task for
      // Use a lower threshold: 1/64 of total nodes or 100, whichever is larger
      size_t minSubtreeForParallel = std::max(size_t(100), parallelThreshold / 64);
      shouldParallelize = (secondLargest >= minSubtreeForParallel);
    }
    
    if (shouldParallelize) {
      // TBB limits concurrent tasks to numThreads automatically
      // Each task clones state, but only numThreads clones exist at once
      // No need for explicit slot counting - TBB's scheduler handles it
      
      // Prepare child offsets for DFS indexing
      std::vector<uint64_t> childOffsets;
      childOffsets.reserve(numChildren);
      uint64_t offset = dfsIndex + 1;
      for (auto* child : node->children) {
        childOffsets.push_back(offset);
        offset += subtreeSizes_[child->identifier];
      }
      
      // Sort children by subtree size descending so largest subtrees start first
      std::vector<size_t> childOrder(numChildren);
      std::iota(childOrder.begin(), childOrder.end(), 0);
      std::sort(childOrder.begin(), childOrder.end(), [&](size_t a, size_t b) {
        return subtreeSizes_[node->children[a]->identifier] > subtreeSizes_[node->children[b]->identifier];
      });
      
      // Find how many children are worth parallelizing (subtree >= minSubtreeForParallel)
      size_t minSubtreeForParallel = std::max(size_t(100), parallelThreshold / 64);
      size_t numWorthParallelizing = 0;
      for (size_t i = 0; i < numChildren; i++) {
        size_t childIdx = childOrder[i];
        if (subtreeSizes_[node->children[childIdx]->identifier] >= minSubtreeForParallel) {
          numWorthParallelizing++;
        } else {
          break;  // Already sorted descending, rest are smaller
        }
      }
      
      // Need at least 2 children worth parallelizing (one runs on current thread)
      if (numWorthParallelizing >= 2) {
        // Spawn one task per large child (instead of worker pool)
        // Each task acquires its own slot, processes subtree, releases slot
        // This allows nested forks to reuse freed slots
        std::vector<std::unordered_set<std::string_view>> childEmptyNodes(numChildren);
        
        // Update instrumentation
        parallelForkPoints_.fetch_add(1, std::memory_order_relaxed);
        
        tbb::task_group taskGroup;
        
        // Spawn tasks for children that can run in parallel
        // Main thread will process one child, so spawn tasks for the others
        for (size_t i = 1; i < numWorthParallelizing; i++) {
          size_t childIdx = childOrder[i];
          auto* child = node->children[childIdx];
          uint64_t childDfsIdx = childOffsets[childIdx];
          
          taskGroup.run([this, child, childIdx, childDfsIdx, &state, 
                         &childEmptyNodes, &globalCoords, &emptyNodesMutex, 
                         parallelThreshold, depth]() {
            // Track clone for instrumentation
            activeClones_.fetch_add(1, std::memory_order_relaxed);
            totalClonesCreated_.fetch_add(1, std::memory_order_relaxed);
            
            // Update peak clones
            int newCount = activeClones_.load(std::memory_order_relaxed);
            int oldPeak = peakActiveClones_.load(std::memory_order_relaxed);
            while (newCount > oldPeak && 
                   !peakActiveClones_.compare_exchange_weak(oldPeak, newCount, std::memory_order_relaxed)) {}
            
            // Clone state for this child
            BuildState childState(state);
            childState.nodeChanges = state.nodeChanges;  // Share nodeChanges storage
            
            processSubtreeParallel(child, childState, globalCoords, 
                                   childEmptyNodes[childIdx], childDfsIdx,
                                   emptyNodesMutex, parallelThreshold, depth + 1);
            
            // Release clone tracking
            activeClones_.fetch_sub(1, std::memory_order_relaxed);
          });
        }
        
        // Main thread processes the first (largest) child with a clone
        {
          size_t childIdx = childOrder[0];
          auto* child = node->children[childIdx];
          uint64_t childDfsIdx = childOffsets[childIdx];
          
          BuildState childState(state);
          childState.nodeChanges = state.nodeChanges;  // Share nodeChanges storage
          
          processSubtreeParallel(child, childState, globalCoords,
                                 childEmptyNodes[childIdx], childDfsIdx,
                                 emptyNodesMutex, parallelThreshold, depth + 1);
        }
        
        taskGroup.wait();
        
        // Process remaining small children sequentially (no cloning, use backtracking)
        for (size_t i = numWorthParallelizing; i < numChildren; i++) {
          size_t childIdx = childOrder[i];
          auto* child = node->children[childIdx];
          uint64_t childDfsIdx = childOffsets[childIdx];
          processSubtreeParallel(child, state, globalCoords, childEmptyNodes[childIdx],
                                 childDfsIdx, emptyNodesMutex, parallelThreshold, depth + 1);
        }
        
        // Merge empty nodes from all children
        for (auto& childEmpty : childEmptyNodes) {
          localEmptyNodes.insert(childEmpty.begin(), childEmpty.end());
        }
        
        // Backtrack and return - we've handled all children
        backtrackNode(state, nodeBacktrackInfo);
        return;
      }
      // Fall through to sequential processing if not enough children to parallelize
      sequentialFallbacks_.fetch_add(1, std::memory_order_relaxed);
    }
    
    {
      // Sequential processing - subtrees too small to benefit from parallelism
      uint64_t offset = dfsIndex + 1;
      for (auto* child : node->children) {
        uint64_t childIdx = offset;
        offset += subtreeSizes_[child->identifier];
        processSubtreeParallel(child, state, globalCoords, localEmptyNodes, childIdx, 
                               emptyNodesMutex, parallelThreshold, depth + 1);
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
  
  // Handle empty tree case
  if (T->root->children.empty()) {
    output::debug("Using sequential build (empty tree)");
    buildIndex();
    return;
  }
  
  output::stage("Building index");
  output::step("Parallel build with {} threads", numThreads);
  
  // Reset instrumentation stats
  processedNodes_.store(0, std::memory_order_relaxed);
  totalClonesCreated_.store(0, std::memory_order_relaxed);
  buildStartTime_ = std::chrono::steady_clock::now();
  lastProgressLog_ = buildStartTime_;
  
  auto startTotal = std::chrono::high_resolution_clock::now();
  
  // Initialize state (same as sequential buildIndex)
  panmapUtils::BlockSequences blockSequences(T);
  panmapUtils::GlobalCoords globalCoords(blockSequences);

  // Step 1: Pre-compute subtree sizes
  computeSubtreeSize(T->root);

  size_t totalNodes = T->allNodes.size();
  output::debug("Tree has {} nodes", totalNodes);

  // Step 2: Compute DFS indices for all nodes
  uint64_t dfsCounter = 0;
  std::function<void(panmanUtils::Node*)> computeDfsIndices = [&](panmanUtils::Node* node) {
    nodeToDfsIndex[node->identifier] = dfsCounter++;
    for (auto* child : node->children) {
      computeDfsIndices(child);
    }
  };
  computeDfsIndices(T->root);

  // Step 3: Collect nodes in DFS order and divide into N equal chunks
  // Simple approach: each chunk owns a contiguous range of DFS indices
  std::vector<panmanUtils::Node*> dfsOrder(totalNodes);
  std::function<void(panmanUtils::Node*, uint64_t&)> collectDfsOrder = [&](panmanUtils::Node* node, uint64_t& idx) {
    dfsOrder[idx++] = node;
    for (auto* child : node->children) {
      collectDfsOrder(child, idx);
    }
  };
  {
    uint64_t idx = 0;
    collectDfsOrder(T->root, idx);
  }
  
  // Divide into N chunks (excluding root which we process separately)
  size_t numChunks = static_cast<size_t>(numThreads);
  size_t nodesPerChunk = (totalNodes - 1 + numChunks - 1) / numChunks;  // Ceiling division, -1 for root
  
  struct ChunkInfo {
    uint64_t startDfs;  // First DFS index this chunk owns (inclusive)
    uint64_t endDfs;    // Last DFS index this chunk owns (exclusive)
  };
  
  std::vector<ChunkInfo> chunks;
  for (size_t i = 0; i < numChunks; i++) {
    uint64_t start = 1 + i * nodesPerChunk;  // +1 to skip root
    uint64_t end = std::min(start + nodesPerChunk, static_cast<uint64_t>(totalNodes));
    if (start < totalNodes) {
      chunks.push_back({start, end});
    }
  }
  
  output::debug("Partitioned into {} chunks of ~{} nodes each", chunks.size(), nodesPerChunk);

  // Step 4: Create template BuildState 
  BuildState templateState;
  templateState.blockSequences = std::move(blockSequences);
  templateState.blockExistsDelayed = templateState.blockSequences.blockExists;
  templateState.blockStrandDelayed = templateState.blockSequences.blockStrand;
  templateState.gapMap = std::map<uint64_t, uint64_t>{{0, globalCoords.lastScalarCoord}};
  // nodeChanges will be set up next

  // Step 5: Process root node to establish initial state
  std::unordered_set<std::string_view> rootEmptyNodes;
  BacktrackInfo rootBacktrackInfo;
  // Create shared storage for nodeChanges
  templateState.nodeChanges = std::make_shared<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>>(totalNodes);
  processNode(T->root, templateState, globalCoords, rootEmptyNodes, 0, &rootBacktrackInfo);

  // Step 6: Make N copies upfront, each walks to its range start then processes its range
  std::vector<std::unordered_set<std::string_view>> chunkEmptyNodes(chunks.size());
  
  // Shared storage - each node slot is written by exactly one chunk
  auto sharedNodeChanges = std::make_shared<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>>(totalNodes);
  
  // Root was already processed
  (*sharedNodeChanges)[0] = std::move((*templateState.nodeChanges)[0]);
  
  auto startDfs = std::chrono::high_resolution_clock::now();
  
  tbb::task_arena arena(numThreads);
  arena.execute([&]() {
    tbb::parallel_for(size_t(0), chunks.size(),
      [&](size_t i) {
        const ChunkInfo& chunk = chunks[i];
        
        // Clone state from template (after root was processed)
        BuildState chunkState(templateState);
        chunkState.nodeChanges = sharedNodeChanges;  // Write to shared storage
        // Genome metrics are already shared via templateState copy
        
        totalClonesCreated_.fetch_add(1, std::memory_order_relaxed);
        
        // We need to walk from root to our first node, then DFS through our range
        // Key insight: nodes 0..startDfs-1 are ancestors or earlier siblings we must walk through
        // We do a DFS but only compute nodeChanges for nodes in our range [startDfs, endDfs)
        
        // Lambda to do DFS with range checking
        std::function<void(panmanUtils::Node*, uint64_t)> processDfsRange = 
          [&](panmanUtils::Node* node, uint64_t dfsIdx) {
            bool inMyRange = (dfsIdx >= chunk.startDfs && dfsIdx < chunk.endDfs);
            
            BacktrackInfo backtrack;
            // Process this node - compute nodeChanges only if in our range
            processNode(node, chunkState, globalCoords, chunkEmptyNodes[i], 
                        dfsIdx, &backtrack, !inMyRange);
            
            // Recurse into children if any of them might be in our range
            if (!node->children.empty()) {
              uint64_t childDfs = dfsIdx + 1;
              for (auto* child : node->children) {
                uint64_t childSubtreeEnd = childDfs + subtreeSizes_[child->identifier];
                // Only recurse if this subtree overlaps with our range
                if (childDfs < chunk.endDfs && childSubtreeEnd > chunk.startDfs) {
                  processDfsRange(child, childDfs);
                }
                childDfs = childSubtreeEnd;
              }
            }
            
            // Backtrack
            backtrackNode(chunkState, backtrack);
          };
        
        // Start DFS from root's children (root already processed)
        uint64_t childDfs = 1;
        for (auto* child : T->root->children) {
          uint64_t childSubtreeEnd = childDfs + subtreeSizes_[child->identifier];
          // Only process this subtree if it overlaps with our range
          if (childDfs < chunk.endDfs && childSubtreeEnd > chunk.startDfs) {
            processDfsRange(child, childDfs);
          }
          childDfs = childSubtreeEnd;
        }
      }
    );
  });
  
  auto endDfs = std::chrono::high_resolution_clock::now();
  auto dfsDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endDfs - startDfs).count();
  
  // Use shared nodeChanges
  auto& nodeChanges = *sharedNodeChanges;
  
  // Merge empty nodes from all chunks
  std::unordered_set<std::string_view> emptyNodes = std::move(rootEmptyNodes);
  for (auto& chunkEmpty : chunkEmptyNodes) {
    emptyNodes.insert(chunkEmpty.begin(), chunkEmpty.end());
  }
  
  // Print summary
  uint64_t finalNodes = processedNodes_.load();
  double finalNodesPerSec = (dfsDurationMs > 0) ? (finalNodes * 1000.0 / dfsDurationMs) : 0;
  
  output::progress_clear();
  output::debug("Build summary: {} nodes ({:.0f}/sec), {} ms, {} chunks", 
               finalNodes, finalNodesPerSec, dfsDurationMs, chunks.size());

  // nodeChanges were computed during processNode calls - now use them for index building
  size_t numNodes = T->allNodes.size();
  
  // Compute totals for index building
  uint64_t totalChanges = 0;
  uint32_t largestNodeChangeCount = 0;
  for (size_t i = 0; i < numNodes; i++) {
    uint32_t localSize = static_cast<uint32_t>(nodeChanges[i].size());
    totalChanges += localSize;
    if (localSize > largestNodeChangeCount) largestNodeChangeCount = localSize;
  }
  
  output::debug("Total seed changes: {}, max per node: {}", totalChanges, largestNodeChangeCount);

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

  output::done(fmt::format("Index built ({} seed changes. {} nodes)", totalChanges, numNodes));
  

  auto serializeStart = std::chrono::high_resolution_clock::now();

  // Cap'n Proto list element limit (~536M for UInt64). Split into primary + overflow.
  constexpr uint32_t CAPNP_SPLIT = 500'000'000;
  size_t numSegments = (totalChanges + CAPNP_SPLIT - 1) / CAPNP_SPLIT;

  // Allocate and write flat arrays (sequential - Cap'n Proto builders aren't thread-safe)
  auto seedChangeHashesBuilder = indexBuilder.initSeedChangeHashes(numSegments);
  auto seedChangeParentCountsBuilder = indexBuilder.initSeedChangeParentCounts(numSegments);
  auto seedChangeChildCountsBuilder = indexBuilder.initSeedChangeChildCounts(numSegments);
  auto nodeChangeOffsetsBuilder = indexBuilder.initNodeChangeOffsets(numNodes + 1);
  
  for (size_t seg = 0; seg < numSegments; seg++) {
    auto segStart = seg * CAPNP_SPLIT;
    auto segEnd = std::min((uint64_t)(segStart + CAPNP_SPLIT), totalChanges);
    auto segSize = segEnd - segStart;
    seedChangeHashesBuilder.init(seg, segSize);
    seedChangeParentCountsBuilder.init(seg, segSize);
    seedChangeChildCountsBuilder.init(seg, segSize);
  }
  if (numSegments > 1) {
    output::step("Large index: splitting {} seed changes into {} segments", totalChanges, numSegments);
  }

  // Sort each node's changes by hash for better compression
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    std::sort(nodeChanges[nodeIdx].begin(), nodeChanges[nodeIdx].end(),
             [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
  }

  output::step("Writing seed changes to index...");
  uint64_t offset = 0;
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    nodeChangeOffsetsBuilder.set(nodeIdx, offset);
    for (const auto& [hash, parentCount, childCount] : nodeChanges[nodeIdx]) {
      if (parentCount < INT16_MIN || parentCount > INT16_MAX || childCount < INT16_MIN || childCount > INT16_MAX) {
        output::error("Seed count overflow at node {}: parentCount={}, childCount={} (Int16 range: [{}, {}]). "
                      "This genome has too many repeated k-min-mers for the current index format.",
                      nodeIdx, parentCount, childCount, INT16_MIN, INT16_MAX);
        std::exit(1);
      }
      auto seg = offset / CAPNP_SPLIT;
      auto segOffset = offset % CAPNP_SPLIT;
      seedChangeHashesBuilder[seg].set(segOffset, hash);
      seedChangeParentCountsBuilder[seg].set(segOffset, static_cast<int16_t>(parentCount));
      seedChangeChildCountsBuilder[seg].set(segOffset, static_cast<int16_t>(childCount));
      offset++;
    }
  }
  nodeChangeOffsetsBuilder.set(numNodes, offset);
  
  auto serializeEnd = std::chrono::high_resolution_clock::now();
  auto serializeMs = std::chrono::duration_cast<std::chrono::milliseconds>(serializeEnd - serializeStart).count();
  output::done("Wrote seed changes", serializeMs);

  auto endTotal = std::chrono::high_resolution_clock::now();
  auto totalDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTotal - startTotal).count();
  
  output::info("Total build time: {} ms", totalDurationMs);
}

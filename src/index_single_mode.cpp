

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
#include <atomic>
#include <mutex>
#include <iomanip>

// Size statistics for optimization
namespace {
  struct SizeStats {
    std::atomic<uint64_t> localRangeSeqTotal{0};
    std::atomic<uint64_t> localRangeSeqMax{0};
    std::atomic<uint64_t> localRangeSeqCount{0};
    std::atomic<uint64_t> seedsToDeleteTotal{0};
    std::atomic<uint64_t> seedsToDeleteMax{0};
    std::atomic<uint64_t> newSyncmerRangesTotal{0};
    std::atomic<uint64_t> newSyncmerRangesMax{0};
    std::atomic<uint64_t> refOnSyncmersChangeRecordTotal{0};
    std::atomic<uint64_t> refOnSyncmersChangeRecordMax{0};
    std::atomic<uint64_t> blockOnSyncmersChangeRecordTotal{0};
    std::atomic<uint64_t> blockOnSyncmersChangeRecordMax{0};
    std::atomic<uint64_t> localMutationRangesTotal{0};
    std::atomic<uint64_t> localMutationRangesMax{0};
    std::atomic<uint64_t> gapRunUpdatesTotal{0};
    std::atomic<uint64_t> gapRunUpdatesMax{0};
    std::atomic<uint64_t> gapMapSizeTotal{0};
    std::atomic<uint64_t> gapMapSizeMax{0};
    std::atomic<uint64_t> refOnSyncmersMapSizeTotal{0};
    std::atomic<uint64_t> refOnSyncmersMapSizeMax{0};
    std::atomic<uint64_t> nodeCount{0};
    
    void updateMax(std::atomic<uint64_t>& maxVal, uint64_t newVal) {
      uint64_t current = maxVal.load();
      while (newVal > current && !maxVal.compare_exchange_weak(current, newVal)) {}
    }
    
    void print() {
      if (nodeCount == 0) return;
      std::cerr << "\n=== SIZE STATISTICS ===" << std::endl;
      std::cerr << "Nodes processed: " << nodeCount << std::endl;
      std::cerr << "localRangeSeq: avg=" << (localRangeSeqTotal/std::max(1UL, localRangeSeqCount.load())) << " max=" << localRangeSeqMax << std::endl;
      std::cerr << "seedsToDelete: avg=" << (seedsToDeleteTotal/nodeCount) << " max=" << seedsToDeleteMax << std::endl;
      std::cerr << "newSyncmerRanges: avg=" << (newSyncmerRangesTotal/nodeCount) << " max=" << newSyncmerRangesMax << std::endl;
      std::cerr << "refOnSyncmersChangeRecord: avg=" << (refOnSyncmersChangeRecordTotal/nodeCount) << " max=" << refOnSyncmersChangeRecordMax << std::endl;
      std::cerr << "blockOnSyncmersChangeRecord: avg=" << (blockOnSyncmersChangeRecordTotal/nodeCount) << " max=" << blockOnSyncmersChangeRecordMax << std::endl;
      std::cerr << "localMutationRanges: avg=" << (localMutationRangesTotal/nodeCount) << " max=" << localMutationRangesMax << std::endl;
      std::cerr << "gapRunUpdates: avg=" << (gapRunUpdatesTotal/nodeCount) << " max=" << gapRunUpdatesMax << std::endl;
      std::cerr << "gapMap size: avg=" << (gapMapSizeTotal/nodeCount) << " max=" << gapMapSizeMax << std::endl;
      std::cerr << "refOnSyncmersMap size: avg=" << (refOnSyncmersMapSizeTotal/nodeCount) << " max=" << refOnSyncmersMapSizeMax << std::endl;
      std::cerr << "======================\n" << std::endl;
    }
  };
  
  SizeStats g_stats;
  
  // Parallelization instrumentation
  struct ParallelStats {
    std::atomic<uint64_t> nodesInParallelSubtrees{0};   // Nodes processed by spawned workers
    std::atomic<uint64_t> nodesInSequentialPaths{0};    // Nodes on main thread only
    std::atomic<uint64_t> cloneTimeNs{0};               // Time spent cloning state
    std::atomic<uint64_t> processNodeTimeNs{0};         // Time in processNode()
    std::atomic<uint64_t> parallelChildrenTotal{0};     // Sum of children at fork points
    std::atomic<uint64_t> workersSpawnedTotal{0};       // Total workers across all forks
    std::atomic<uint64_t> maxForkDepth{0};              // Deepest nested fork
    std::atomic<uint64_t> subtreeSizesProcessedParallel{0}; // Sum of subtree sizes in parallel
    std::atomic<uint64_t> forksByDepth[16]{};           // Histogram of fork depths
    
    void updateMax(std::atomic<uint64_t>& maxVal, uint64_t newVal) {
      uint64_t current = maxVal.load();
      while (newVal > current && !maxVal.compare_exchange_weak(current, newVal)) {}
    }
    
    void recordFork(size_t depth, size_t numChildren, size_t numWorkers, size_t totalSubtreeSize) {
      if (depth < 16) forksByDepth[depth].fetch_add(1, std::memory_order_relaxed);
      updateMax(maxForkDepth, depth);
      parallelChildrenTotal.fetch_add(numChildren, std::memory_order_relaxed);
      workersSpawnedTotal.fetch_add(numWorkers, std::memory_order_relaxed);
      subtreeSizesProcessedParallel.fetch_add(totalSubtreeSize, std::memory_order_relaxed);
    }
    
    void print(uint64_t totalNodes, uint64_t totalForks, uint64_t wallTimeMs) {
      std::cout << "\n=== Parallelization Analysis ===" << std::endl;
      
      double parallelPct = 100.0 * subtreeSizesProcessedParallel / std::max(1UL, totalNodes);
      std::cout << "Nodes in parallel subtrees: " << subtreeSizesProcessedParallel 
                << " (" << std::fixed << std::setprecision(1) << parallelPct << "% of tree)" << std::endl;
      
      if (totalForks > 0) {
        std::cout << "Avg children per fork: " << (parallelChildrenTotal / totalForks) << std::endl;
        std::cout << "Avg workers per fork: " << (workersSpawnedTotal / totalForks) << std::endl;
      }
      std::cout << "Max fork depth: " << maxForkDepth << std::endl;
      
      std::cout << "Fork depth distribution: ";
      for (int i = 0; i < 16 && forksByDepth[i] > 0; i++) {
        std::cout << "d" << i << ":" << forksByDepth[i] << " ";
      }
      std::cout << std::endl;
      
      if (cloneTimeNs > 0) {
        double cloneTimeSec = cloneTimeNs / 1e9;
        double clonePct = 100.0 * cloneTimeSec / (wallTimeMs / 1000.0);
        std::cout << "Clone overhead: " << std::setprecision(2) << cloneTimeSec 
                  << "s (" << clonePct << "% of wall time)" << std::endl;
      }
      
      // Theoretical speedup analysis
      double sequentialFraction = 1.0 - (parallelPct / 100.0);
      std::cout << "Amdahl's Law prediction (sequential fraction=" << std::setprecision(3) << sequentialFraction << "):" << std::endl;
      for (int t : {4, 8, 16, 32, 48}) {
        double speedup = 1.0 / (sequentialFraction + (1.0 - sequentialFraction) / t);
        std::cout << "  " << t << " threads: max " << std::setprecision(2) << speedup << "x speedup" << std::endl;
      }
      
      std::cout << "================================\n" << std::endl;
    }
  };
  
  ParallelStats g_parallelStats;
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
      
      // N imputation: Skip XX->N mutations during indexing
      // When a node has 'N' (ambiguous base), we treat it as inheriting the parent's base
      // This effectively imputes the N with the ancestral value for syncmer computation
      if (newNuc == 'N') continue;
      
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

  // Compute genome extent from gapMap for flank masking
  // Pass flankSize=0 - actual flank masking is applied separately
  auto [firstNonGapScalar, lastNonGapScalar] = computeExtentFromGapMap(gapMap, globalCoords.lastScalarCoord, 0);

  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = computeNewSyncmerRangesJump(
      node, dfsIndex, blockSequences, blockExistsDelayed, blockStrandDelayed, 
      globalCoords, gapMap, localMutationRanges, blockOnSyncmersChangeRecord, 
      refOnSyncmers, blockOnSyncmers, firstNonGapScalar, lastNonGapScalar);

  // processing syncmers
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    if (localRangeSeq.size() >= indexBuilder.getK()) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(localRangeSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = localRangeCoordToGlobalScalarCoords[startPos];
        auto endPosGlobal = localRangeCoordToGlobalScalarCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = localRangeCoordToBlockId[startPos];
        bool wasSeed = refOnSyncmers.contains(startPosGlobal);
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
          if (blockOnSyncmers[localRangeCoordToBlockId[startPos]].empty()) {
            blockOnSyncmers.erase(localRangeCoordToBlockId[startPos]);
          }
        } else if (wasSeed && isSeed) {
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::SUB, refOnSyncmers.at(startPosGlobal));
          refOnSyncmers[startPosGlobal] = {hash, endPosGlobal, isReverse};
        }
      }
    }

    for (uint64_t pos : seedsToDelete) {
      if (!refOnSyncmers.contains(pos)) {
        std::cerr << "Error: refOnSyncmers[" << pos << "] is null" << std::endl;
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
        for (uint64_t pos : blockOnSyncmers[blockId]) {
          refOnSyncmersChangeRecord.emplace_back(pos, panmapUtils::seedChangeType::DEL, refOnSyncmers.at(pos));
          blockOnSyncmersChangeRecord.emplace_back(blockId, pos, panmapUtils::seedChangeType::DEL);

          refOnSyncmers.erase(pos);
          refOnSyncmersMap.erase(pos);
        }
        blockOnSyncmers.erase(blockId);
      }
    }
  }

  // processing k-min-mers
  std::vector<std::pair<index_single_mode::SyncmerSet::iterator, index_single_mode::SyncmerSet::iterator>> newKminmerRanges = computeNewKminmerRanges(refOnSyncmersChangeRecord, dfsIndex);

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
        std::cerr << "Error: syncmer hash collision detected in k-min-mer computation.\n";
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

void index_single_mode::IndexBuilder::writeIndex(const std::string& path, int numThreads) {
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
  if (!panmap_zstd::compressToFile(data, dataSize, path, 3, numThreads)) {
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

  // Update genome extent from gapMap - used to determine flank regions
  // Save previous extent for backtracking
  if (backtrackInfo) {
    backtrackInfo->prevFirstNonGapScalar = state.firstNonGapScalar;
    backtrackInfo->prevLastNonGapScalar = state.lastNonGapScalar;
  }
  // Pass flankSize=0 here - flank masking is applied separately via hardMaskStart/hardMaskEnd below
  auto [newFirstNonGap, newLastNonGap] = computeExtentFromGapMap(state.gapMap, globalCoords.lastScalarCoord, 0);
  
  // Debug: log extent changes for nodes we care about
  static bool debugFlankMasking = (std::getenv("DEBUG_FLANK_MASKING") != nullptr);
  if (debugFlankMasking) {
    if (state.firstNonGapScalar != newFirstNonGap || state.lastNonGapScalar != newLastNonGap) {
      std::cerr << "[FLANK] Node " << node->identifier << " (dfs=" << dfsIndex << "): "
                << "extent [" << state.firstNonGapScalar << "," << state.lastNonGapScalar << "] -> "
                << "[" << newFirstNonGap << "," << newLastNonGap << "]" << std::endl;
    }
  }
  
  state.firstNonGapScalar = newFirstNonGap;
  state.lastNonGapScalar = newLastNonGap;
  
  // Compute hard mask boundaries relative to each genome's extent
  // Seeds in hard-masked regions are completely ignored (no adds, no deletes)
  // Seeds in flank regions (between hard mask and extent) have deletions skipped (imputed)
  const int flankMaskBp = getFlankMaskBp();
  
  // Hard mask is relative to the genome's extent (first/last N bp of actual sequence)
  // Not relative to global scalar coordinates
  const uint64_t hardMaskStart = state.firstNonGapScalar + static_cast<uint64_t>(flankMaskBp);
  const uint64_t hardMaskEnd = state.lastNonGapScalar >= static_cast<uint64_t>(flankMaskBp) 
                                ? state.lastNonGapScalar - flankMaskBp 
                                : state.firstNonGapScalar;
  
  // Debug logging for specific nodes
  static bool debugMasking = (std::getenv("DEBUG_MASK") != nullptr);
  if (debugMasking && (node->identifier.find("node_9981") != std::string::npos || 
                        node->identifier.find("node_11963") != std::string::npos ||
                        node->identifier.find("node_11964") != std::string::npos)) {
    std::cerr << "[MASK-DEBUG] " << node->identifier << " extent=[" << state.firstNonGapScalar 
              << "," << state.lastNonGapScalar << "] hardMask=[" << hardMaskStart << "," << hardMaskEnd 
              << "] flankMaskBp=" << flankMaskBp << std::endl;
  }

  std::vector<panmapUtils::NewSyncmerRange> newSyncmerRanges = computeNewSyncmerRangesJump(
      node, dfsIndex, state.blockSequences, state.blockExistsDelayed, state.blockStrandDelayed, 
      globalCoords, state.gapMap, localMutationRanges, blockOnSyncmersChangeRecord, 
      state.refOnSyncmers, state.blockOnSyncmers, state.firstNonGapScalar, state.lastNonGapScalar);

  // Processing syncmers with dual masking:
  // 1. HARD MASK: First/last flankMaskBp positions - completely ignore all seeds
  for (const auto& syncmerRange : newSyncmerRanges) {
    const auto& [begCoord, endCoord, localRangeSeq, localRangeCoordToGlobalScalarCoords, localRangeCoordToBlockId, seedsToDelete] = syncmerRange;
    if (localRangeSeq.size() >= static_cast<size_t>(indexBuilder.getK())) {
      for (auto [hash, isReverse, isSeed, startPos] : seeding::rollingSyncmers(localRangeSeq, indexBuilder.getK(), indexBuilder.getS(), indexBuilder.getOpen(), indexBuilder.getT(), true)) {
        auto startPosGlobal = localRangeCoordToGlobalScalarCoords[startPos];
        auto endPosGlobal = localRangeCoordToGlobalScalarCoords[startPos + indexBuilder.getK() - 1];
        auto curBlockId = localRangeCoordToBlockId[startPos];
        bool wasSeed = state.refOnSyncmers.contains(startPosGlobal);
        
        // Check masking conditions:
        // 1. Hard mask: position < hardMaskStart OR position > hardMaskEnd
        bool isHardMasked = (startPosGlobal < hardMaskStart || startPosGlobal > hardMaskEnd);
        // 2. Flank imputation: position outside genome extent but not hard-masked
        bool isInFlank = !isHardMasked && (startPosGlobal < state.firstNonGapScalar || startPosGlobal > state.lastNonGapScalar);
        
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
          // if (isInFlank) {
          //   // Skip this deletion - keep the seed imputed
          //   continue;
          // }
          refOnSyncmersChangeRecord.emplace_back(startPosGlobal, panmapUtils::seedChangeType::DEL, state.refOnSyncmers.at(startPosGlobal));
          blockOnSyncmersChangeRecord.emplace_back(curBlockId, startPosGlobal, panmapUtils::seedChangeType::DEL);
          state.refOnSyncmers.erase(startPosGlobal);
          state.refOnSyncmersMap.erase(startPosGlobal);
          state.blockOnSyncmers[curBlockId].erase(startPosGlobal);
          if (state.blockOnSyncmers[localRangeCoordToBlockId[startPos]].empty()) {
            state.blockOnSyncmers.erase(localRangeCoordToBlockId[startPos]);
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
      bool isInFlank = !isHardMasked && (pos < state.firstNonGapScalar || pos > state.lastNonGapScalar);
      
      // Hard-masked: skip deletion
      if (isHardMasked) {
        continue;
      }
      // Flank imputation: skip deletion
      // if (isInFlank) {
      //   continue;
      // }
      
      if (!state.refOnSyncmers.contains(pos)) {
        std::cerr << "Error: refOnSyncmers[" << pos << "] is null" << std::endl;
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
      bool isInFlank = !isHardMasked && (pos < state.firstNonGapScalar || pos > state.lastNonGapScalar);
      
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
          bool isInFlank = !isHardMasked && (pos < state.firstNonGapScalar || pos > state.lastNonGapScalar);
          
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
        std::cerr << "Error: syncmer hash collision detected in k-min-mer computation.\n";
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

  // Collect size statistics for optimization
  {
    g_stats.nodeCount++;
    
    // Aggregate seedsToDelete across all syncmer ranges
    uint64_t totalSeedsToDelete = 0;
    uint64_t maxLocalRangeSeq = 0;
    for (const auto& sr : newSyncmerRanges) {
      totalSeedsToDelete += sr.seedsToDelete.size();
      g_stats.localRangeSeqTotal += sr.localRangeSeq.size();
      g_stats.localRangeSeqCount++;
      if (sr.localRangeSeq.size() > maxLocalRangeSeq) maxLocalRangeSeq = sr.localRangeSeq.size();
    }
    g_stats.updateMax(g_stats.localRangeSeqMax, maxLocalRangeSeq);
    g_stats.seedsToDeleteTotal += totalSeedsToDelete;
    g_stats.updateMax(g_stats.seedsToDeleteMax, totalSeedsToDelete);
    
    g_stats.newSyncmerRangesTotal += newSyncmerRanges.size();
    g_stats.updateMax(g_stats.newSyncmerRangesMax, newSyncmerRanges.size());
    
    g_stats.refOnSyncmersChangeRecordTotal += refOnSyncmersChangeRecord.size();
    g_stats.updateMax(g_stats.refOnSyncmersChangeRecordMax, refOnSyncmersChangeRecord.size());
    
    g_stats.blockOnSyncmersChangeRecordTotal += blockOnSyncmersChangeRecord.size();
    g_stats.updateMax(g_stats.blockOnSyncmersChangeRecordMax, blockOnSyncmersChangeRecord.size());
    
    g_stats.localMutationRangesTotal += localMutationRanges.size();
    g_stats.updateMax(g_stats.localMutationRangesMax, localMutationRanges.size());
    
    g_stats.gapRunUpdatesTotal += gapRunUpdates.size();
    g_stats.updateMax(g_stats.gapRunUpdatesMax, gapRunUpdates.size());
    
    g_stats.gapMapSizeTotal += state.gapMap.size();
    g_stats.updateMax(g_stats.gapMapSizeMax, state.gapMap.size());
    
    g_stats.refOnSyncmersMapSizeTotal += state.refOnSyncmersMap.size();
    g_stats.updateMax(g_stats.refOnSyncmersMapSizeMax, state.refOnSyncmersMap.size());
  }

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
  revertGapMapInversions(gapRunBlockInversionBacktracks, state.gapMap);

  // Update delayed block states (for children)
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    state.blockExistsDelayed[blockId] = newExists;
    state.blockStrandDelayed[blockId] = newStrand;
  }

  // nodeToDfsIndex is pre-computed, no need to update here

  // Compute genome metrics from runningCounts (only when we own this node)
  if (!skipNodeChanges && state.genomeMagnitudeSquared) {
    double magnitudeSquared = 0.0;
    int64_t totalFrequency = 0;
    for (const auto& [hash, count] : state.runningCounts) {
      magnitudeSquared += static_cast<double>(count) * static_cast<double>(count);
      totalFrequency += count;
    }
    (*state.genomeMagnitudeSquared)[dfsIndex] = magnitudeSquared;
    (*state.genomeUniqueSeedCount)[dfsIndex] = state.runningCounts.size();
    (*state.genomeTotalSeedFrequency)[dfsIndex] = totalFrequency;
    
    // IDF collection: For leaf nodes, record which seeds they contain
    // A node is a leaf if it has no children
    if (state.isLeafNode && node->children.empty()) {
      (*state.isLeafNode)[dfsIndex] = true;
      
      // If we have IDF tracking, collect seeds for this leaf
      if (state.leafSeedGenomeCounts) {
        for (const auto& [hash, count] : state.runningCounts) {
          if (count > 0) {
            (*state.leafSeedGenomeCounts)[hash]++;
          }
        }
      }
    }
  }

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
    std::cout << "\r[" << std::fixed << std::setprecision(1) << pct << "%] " 
              << processed << "/" << totalNodes << " (" << std::setprecision(0) << nodesPerSec << " n/s) | "
              << "chunks:" << totalClones
              << "          " << std::flush;
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
        g_parallelStats.recordFork(depth, numWorthParallelizing, 
                                   static_cast<int>(numWorthParallelizing - 1),
                                   totalChildSubtreeSize);
        
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
            auto cloneStart = std::chrono::high_resolution_clock::now();
            BuildState childState(state);
            childState.nodeChanges = state.nodeChanges;  // Share nodeChanges storage
            auto cloneEnd = std::chrono::high_resolution_clock::now();
            g_parallelStats.cloneTimeNs.fetch_add(
                std::chrono::duration_cast<std::chrono::nanoseconds>(cloneEnd - cloneStart).count(),
                std::memory_order_relaxed);
            
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
          
          auto cloneStart = std::chrono::high_resolution_clock::now();
          BuildState childState(state);
          childState.nodeChanges = state.nodeChanges;  // Share nodeChanges storage
          auto cloneEnd = std::chrono::high_resolution_clock::now();
          g_parallelStats.cloneTimeNs.fetch_add(
              std::chrono::duration_cast<std::chrono::nanoseconds>(cloneEnd - cloneStart).count(),
              std::memory_order_relaxed);
          
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
    std::cout << "Using sequential build (empty tree)" << std::endl;
    buildIndex();
    return;
  }
  
  std::cout << "Building index in parallel with " << numThreads << " threads..." << std::endl;
  
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
  std::cout << "Tree has " << totalNodes << " nodes" << std::endl;

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
  
  std::cout << "Partitioned into " << chunks.size() << " chunks of ~" 
            << nodesPerChunk << " nodes each" << std::endl;

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
  // Create shared storage for nodeChanges and genome metrics
  templateState.nodeChanges = std::make_shared<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>>(totalNodes);
  templateState.genomeMagnitudeSquared = std::make_shared<std::vector<double>>(totalNodes);
  templateState.genomeUniqueSeedCount = std::make_shared<std::vector<uint64_t>>(totalNodes);
  templateState.genomeTotalSeedFrequency = std::make_shared<std::vector<int64_t>>(totalNodes);
  templateState.isLeafNode = std::make_shared<std::vector<bool>>(totalNodes, false);
  // Note: leafSeedGenomeCounts will be set per-chunk below
  processNode(T->root, templateState, globalCoords, rootEmptyNodes, 0, &rootBacktrackInfo);

  // Step 6: Make N copies upfront, each walks to its range start then processes its range
  std::vector<std::unordered_set<std::string_view>> chunkEmptyNodes(chunks.size());
  
  // Per-chunk IDF maps - each chunk collects IDF data independently, merged at end
  std::vector<std::unordered_map<uint64_t, uint32_t>> chunkIdfMaps(chunks.size());
  
  // Shared storage - each node slot is written by exactly one chunk
  auto sharedNodeChanges = std::make_shared<std::vector<std::vector<std::tuple<uint64_t, int64_t, int64_t>>>>(totalNodes);
  auto sharedGenomeMagnitudeSquared = templateState.genomeMagnitudeSquared;  // Already created, share it
  auto sharedGenomeUniqueSeedCount = templateState.genomeUniqueSeedCount;
  auto sharedIsLeafNode = templateState.isLeafNode;
  auto sharedGenomeTotalSeedFrequency = templateState.genomeTotalSeedFrequency;
  
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
        chunkState.isLeafNode = sharedIsLeafNode;    // Share leaf tracking
        // Set per-chunk IDF map
        chunkState.leafSeedGenomeCounts = std::make_shared<std::unordered_map<uint64_t, uint32_t>>();
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
        
        // Save chunk's IDF map for merging after parallel section
        chunkIdfMaps[i] = std::move(*chunkState.leafSeedGenomeCounts);
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
  
  // Merge IDF maps from all chunks
  std::unordered_map<uint64_t, uint32_t> mergedIdfMap;
  uint32_t totalLeafGenomes = 0;
  for (size_t i = 0; i < totalNodes; i++) {
    if ((*sharedIsLeafNode)[i]) totalLeafGenomes++;
  }
  
  // Reserve space for merging
  size_t estimatedSeeds = 0;
  for (const auto& chunkMap : chunkIdfMaps) {
    estimatedSeeds += chunkMap.size();
  }
  mergedIdfMap.reserve(estimatedSeeds);
  
  // Merge counts
  for (auto& chunkMap : chunkIdfMaps) {
    for (const auto& [hash, count] : chunkMap) {
      mergedIdfMap[hash] += count;
    }
  }
  
  std::cout << "IDF: " << totalLeafGenomes << " leaf genomes, " 
            << mergedIdfMap.size() << " unique seeds" << std::endl;
  
  // Print summary
  uint64_t finalNodes = processedNodes_.load();
  double finalNodesPerSec = (dfsDurationMs > 0) ? (finalNodes * 1000.0 / dfsDurationMs) : 0;
  
  std::cout << std::endl;
  std::cout << "=== Parallel Build Summary ===" << std::endl;
  std::cout << "  Nodes processed: " << finalNodes << " (" << std::fixed << std::setprecision(0) << finalNodesPerSec << " nodes/sec)" << std::endl;
  std::cout << "  DFS time: " << dfsDurationMs << " ms" << std::endl;
  std::cout << "  Chunks: " << chunks.size() << ", clones: " << totalClonesCreated_.load() << std::endl;
  std::cout << "===============================" << std::endl;

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
  
  std::cout << "Total seed changes: " << totalChanges << ", max per node: " << largestNodeChangeCount << std::endl;

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

  std::cout << "\nFinished building index! Total seed changes: " << totalChanges << std::endl;

  auto serializeStart = std::chrono::high_resolution_clock::now();

  // Set totals
  indexBuilder.setTotalSeedChanges(totalChanges);
  indexBuilder.setLargestNodeChangeCount(largestNodeChangeCount);

  // Allocate and write flat arrays (sequential - Cap'n Proto builders aren't thread-safe)
  auto seedChangeHashesBuilder = indexBuilder.initSeedChangeHashes(totalChanges);
  auto seedChangeParentCountsBuilder = indexBuilder.initSeedChangeParentCounts(totalChanges);
  auto seedChangeChildCountsBuilder = indexBuilder.initSeedChangeChildCounts(totalChanges);
  auto nodeChangeOffsetsBuilder = indexBuilder.initNodeChangeOffsets(numNodes + 1);

  std::cout << "Writing seed changes to index..." << std::endl;
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
  
  auto serializeEnd = std::chrono::high_resolution_clock::now();
  auto serializeMs = std::chrono::duration_cast<std::chrono::milliseconds>(serializeEnd - serializeStart).count();
  std::cout << "Wrote seed changes in " << serializeMs << " ms" << std::endl;

  // Copy pre-computed genome metrics to index builders
  std::cout << "Writing genome metrics..." << std::endl;
  auto metricsStart = std::chrono::high_resolution_clock::now();
  auto genomeMagnitudeSquaredBuilder = indexBuilder.initGenomeMagnitudeSquared(numNodes);
  auto genomeUniqueSeedCountBuilder = indexBuilder.initGenomeUniqueSeedCount(numNodes);
  auto genomeTotalSeedFrequencyBuilder = indexBuilder.initGenomeTotalSeedFrequency(numNodes);

  // Metrics were computed during main DFS - just copy them
  for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
    genomeMagnitudeSquaredBuilder.set(nodeIdx, (*sharedGenomeMagnitudeSquared)[nodeIdx]);
    genomeUniqueSeedCountBuilder.set(nodeIdx, (*sharedGenomeUniqueSeedCount)[nodeIdx]);
    genomeTotalSeedFrequencyBuilder.set(nodeIdx, (*sharedGenomeTotalSeedFrequency)[nodeIdx]);
  }

  auto metricsEnd = std::chrono::high_resolution_clock::now();
  auto metricsMs = std::chrono::duration_cast<std::chrono::milliseconds>(metricsEnd - metricsStart).count();
  std::cout << "Wrote genome metrics in " << metricsMs << " ms" << std::endl;

  // Write IDF data for seed weighting
  std::cout << "Writing IDF data..." << std::endl;
  auto idfStart = std::chrono::high_resolution_clock::now();
  
  // Sort seeds by hash for binary search at load time
  std::vector<std::pair<uint64_t, uint32_t>> sortedIdf(mergedIdfMap.begin(), mergedIdfMap.end());
  std::sort(sortedIdf.begin(), sortedIdf.end(), 
            [](const auto& a, const auto& b) { return a.first < b.first; });
  
  indexBuilder.setTotalLeafGenomes(totalLeafGenomes);
  auto idfHashesBuilder = indexBuilder.initIdfSeedHashes(sortedIdf.size());
  auto idfCountsBuilder = indexBuilder.initIdfGenomeCounts(sortedIdf.size());
  
  for (size_t i = 0; i < sortedIdf.size(); i++) {
    idfHashesBuilder.set(i, sortedIdf[i].first);
    idfCountsBuilder.set(i, sortedIdf[i].second);
  }
  
  auto idfEnd = std::chrono::high_resolution_clock::now();
  auto idfMs = std::chrono::duration_cast<std::chrono::milliseconds>(idfEnd - idfStart).count();
  std::cout << "Wrote IDF data (" << sortedIdf.size() << " unique seeds) in " << idfMs << " ms" << std::endl;

  auto endTotal = std::chrono::high_resolution_clock::now();
  auto totalDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTotal - startTotal).count();
  
  std::cout << "Total parallel build time: " << totalDurationMs << " ms" << std::endl;
  
  // Print size statistics for optimization
  g_stats.print();
}

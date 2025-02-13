#ifndef __MGSR_HPP
#define __MGSR_HPP

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tbb/concurrent_vector.h>
#include <boost/icl/interval_map.hpp>
#include <boost/icl/split_interval_map.hpp>
#include "index.capnp.h"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include "haplotype_filter.hpp"
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Dense>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/combinable.h>
#include <tbb/task_group.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>


using namespace boost::icl;
typedef std::pair<std::vector<std::tuple<size_t*, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> readSeedmers_t;
typedef std::tuple<int32_t, int32_t, int32_t, int32_t, bool, int32_t> match_t;

using namespace haplotype_filter;

namespace mgsr {

  struct positionInfo {
    int64_t endPos;
    size_t fhash;
    size_t rhash;
    bool rev;
  };

  struct IteratorComparator {
      bool operator()(const std::map<int32_t, positionInfo>::iterator& lhs,
                      const std::map<int32_t, positionInfo>::iterator& rhs) const {
          return lhs->first < rhs->first;
      }
  };

  struct seedmers {
    //       beg                 end      fhash    rhash    rev
    std::map<int32_t, positionInfo> positionMap;
    //                 hash                       begs
    std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>> hashToPositionsMap;

    void addPosition(const int32_t& beg, const int32_t& end, const size_t& fhash, const size_t& rhash, const bool& rev, std::unordered_set<size_t>& affectedSeedmers) {
      auto it = positionMap.emplace(beg, positionInfo(end, fhash, rhash, rev)).first;
      if (fhash != rhash) {
        size_t minHash = std::min(fhash, rhash);
        hashToPositionsMap[minHash].insert(it);
        affectedSeedmers.insert(minHash);
      }
    }

    void subPosition(std::map<int32_t, positionInfo>::iterator& it, const int32_t& end, const size_t& fhash, const size_t& rhash, const bool& rev, std::unordered_set<size_t>& affectedSeedmers) {
      const auto& [oldEnd, oldFHash, oldRHash, oldRev] = it->second;
      if (oldFHash != oldRHash) {
        size_t minHash = std::min(oldFHash, oldRHash);
        auto hashToPositionIt = hashToPositionsMap.find(minHash);
        hashToPositionIt->second.erase(it);
        if (hashToPositionIt->second.empty()) hashToPositionsMap.erase(hashToPositionIt);
        affectedSeedmers.insert(minHash);
      }

      it->second.endPos = end;
      it->second.fhash = fhash;
      it->second.rhash = rhash;
      it->second.rev = rev;
      if (fhash != rhash) {
        size_t minHash = std::min(fhash, rhash);
        hashToPositionsMap[minHash].insert(it);
        affectedSeedmers.insert(minHash);
      }
    }

    void delPosition(std::map<int32_t, positionInfo>::iterator& it, std::unordered_set<size_t>& affectedSeedmers) {
      const auto& [end, fhash, rhash, rev] = it->second;
      if (fhash != rhash) {
        size_t minHash = std::min(fhash, rhash);
        auto hashToPositionIt = hashToPositionsMap.find(minHash);
        hashToPositionIt->second.erase(it);
        if (hashToPositionIt->second.empty()) hashToPositionsMap.erase(hashToPositionIt);
        affectedSeedmers.insert(minHash);
      }
      positionMap.erase(it);
    }

    int32_t getBegFromHash(const size_t& hash, const bool& safe=false) {
      auto hashToPositionIt = hashToPositionsMap.find(hash);
      if (safe && (hashToPositionIt == hashToPositionsMap.end() || hashToPositionIt->second.size() != 1)) return -1;
      return (*(hashToPositionIt->second.begin()))->first;
    }

    int32_t getEndFromHash(const size_t& hash, const bool& safe=false) {
      auto hashToPositionIt = hashToPositionsMap.find(hash);
      if (safe && (hashToPositionIt == hashToPositionsMap.end() || hashToPositionIt->second.size() != 1)) return -1;
      auto positionIt = *(hashToPositionIt->second.begin());
      return positionIt->second.endPos;
    }
  };

  struct readSeedmer {
    const size_t hash;
    const int64_t begPos;
    const int64_t endPos;
    const bool rev;
    const int32_t iorder;
  };

  class Read {
    public:
    std::vector<readSeedmer> seedmersList;
    std::unordered_map<size_t, std::vector<int32_t>> uniqueSeedmers;
    boost::icl::split_interval_map<int32_t, int> matches;
    std::unordered_set<int32_t> duplicates;
    size_t readIndex;
    // std::vector<bool> duplicates;
    // std::vector<bool> absentees;
    // size_t numDuplicates = 0;
    // size_t numAbsentees = 0;
  };

  int64_t degapGlobal(const int64_t& globalCoord, const std::map<int64_t, int64_t>& degapCoordsIndex) {
    auto coordIt = degapCoordsIndex.upper_bound(globalCoord);
    if (coordIt == degapCoordsIndex.begin()) {
        return 0;
    }
    return globalCoord - std::prev(coordIt)->second;
  }

  int64_t regapGlobal(const int64_t& localCoord, const std::map<int64_t, int64_t>& regapCoordsIndex) {
    auto coordIt = regapCoordsIndex.upper_bound(localCoord);
    if (coordIt == regapCoordsIndex.begin()) {
        return 0;
    }
    return localCoord + std::prev(coordIt)->second;
  }

  template<typename Iterator>
  Iterator safe_next(Iterator it, const Iterator& end, const uint8_t steps) {
    // Check if advancing will go past the end iterator
    uint8_t steps_taken = 0;
    while (steps_taken < steps && it != end) {
      ++it;
      ++steps_taken;
    }
    return it;
  }

  // Gap Mutations Processing Function
  template <typename GapMutationsType>
  void processGapMutations(
    const GapMutationsType& perNodeGapMutations_Index, 
    const int64_t& dfsIndex, 
    std::map<int64_t, int64_t>& gapMap, 
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapRunBacktracks, 
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapRunBlocksBacktracks,
    const std::unordered_set<int64_t>& inverseBlockIds, 
    const std::vector<std::pair<int64_t, int64_t>>& blockRanges, 
    mutableTreeData& data
  ) {
    ::capnp::List<MapDelta>::Reader gapMutationsList;
    if constexpr (std::is_same_v<GapMutationsType, ::capnp::List<GapMutations>::Reader>) {
      gapMutationsList = perNodeGapMutations_Index[dfsIndex].getDeltas();
    }

    for (int i = 0; i < gapMutationsList.size(); ++i) {
      const auto& gapMutation = gapMutationsList[i];
      const int32_t& pos = gapMutation.getPos();
      const auto& maybeValue = gapMutation.getMaybeValue();
      if (maybeValue.isValue()) {
        if (gapMap.find(pos) != gapMap.end()) {
          gapRunBacktracks.emplace_back(false, std::make_pair(pos, gapMap[pos]));
        } else {
          gapRunBacktracks.emplace_back(true, std::make_pair(pos, maybeValue.getValue()));
        }
        gapMap[pos] = maybeValue.getValue();
      } else {
        gapRunBacktracks.emplace_back(false, std::make_pair(pos, gapMap[pos]));
        gapMap.erase(pos);
      }
    }

    for (const auto& blockId : inverseBlockIds) {
      std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> tmpGapMapUpdates;
      if (data.blockExists[blockId].first) {
        invertGapMap(gapMap, blockRanges[blockId], gapRunBlocksBacktracks, tmpGapMapUpdates);
      }
    }
  }

  // Seed Mutations Processing Function
  template <typename SeedMutationsType>
  void processSeedMutations(
    const SeedMutationsType& perNodeSeedMutations_Index, 
    int64_t& dfsIndex, 
    std::map<uint32_t, seeding::onSeedsHash>& onSeedsHashMap, 
    std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>& seedChanges, 
    Tree* T, 
    const int& seedK, 
    const globalCoords_t& globalCoords, 
    CoordNavigator& navigator, 
    mutableTreeData& data, 
    std::map<int64_t, int64_t>& gapMap, 
    const std::vector<std::pair<int64_t, int64_t>>& blockRanges
  ) {
    auto currBasePositions = perNodeSeedMutations_Index[dfsIndex].getBasePositions();
    auto currPerPosMasks = perNodeSeedMutations_Index[dfsIndex].getPerPosMasks();
    thread_local static std::string seedBuffer;

    for (int i = 0; i < currBasePositions.size(); ++i) {
      int64_t pos = currBasePositions[i];
      uint64_t tritMask = currPerPosMasks[i];
      for (int k = 0; k < 32; ++k) {
        uint8_t ternaryNumber = (tritMask >> (k * 2)) & 0x3;
        int64_t seedPos = pos - k;
        auto seedIt = onSeedsHashMap.find(seedPos);

        if (ternaryNumber == 1) { // on -> off
          if (seedIt != onSeedsHashMap.end()) {
            const auto& [oldSeed, oldEndPos, oldIsReverse] = seedIt->second;
            seedChanges.emplace_back(std::make_tuple(seedPos, true, false, oldSeed, std::nullopt, oldIsReverse, std::nullopt, oldEndPos, std::nullopt));
          }
        } else if (ternaryNumber == 2) { // Handle insertion/change (off -> on or on -> on)
          size_t result_hash;
          int64_t result_end_pos = -1;
          bool result_is_reverse;
          size_t newSeedHash;
          bool newIsReverse;
          HotSeedIndex hotSeedIndexDummy(0); // todo use hotSeedIndex
          if(!seed_annotated_tree::getSeedAt(false, hotSeedIndexDummy, seedBuffer, result_hash, result_end_pos, result_is_reverse, seedPos, T, seedK, dfsIndex, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap, blockRanges)) {
            continue;
          }
          // getSeedAt fills result_hash and result_end_pos if cache hit
          // else fills seedBuffer and result_end_pos, need to hashSeq(seedBuffer)
          if (result_end_pos != -1) {
            newSeedHash = result_hash;
            newIsReverse = result_is_reverse;
          } else {
            auto [newSeedFHash, newSeedRHash] = hashSeq(seedBuffer);
            if (newSeedFHash < newSeedRHash) {
              newSeedHash  = newSeedFHash;
              newIsReverse = false;
            } else {
              newSeedHash  = newSeedRHash;
              newIsReverse = true;
            }
          }

          if (seedIt != onSeedsHashMap.end()) { // on -> on
            const auto& [oldSeed, oldEndPos, oldIsReverse] = seedIt->second;
            seedChanges.emplace_back(std::make_tuple(seedPos, true, true, oldSeed, newSeedHash, oldIsReverse, newIsReverse, oldEndPos, result_end_pos));
          } else { // off -> on
            seedChanges.emplace_back(std::make_tuple(seedPos, false, true, std::nullopt, newSeedHash, std::nullopt, newIsReverse, std::nullopt, result_end_pos));
          }
        }
      }
    }
  }

  template <typename SeedMutationsType, typename GapMutationsType>
  void processNodeMutations(
    const GapMutationsType& perNodeGapMutations_Index,
    const SeedMutationsType& perNodeSeedMutations_Index,
    int64_t& dfsIndex,
    std::map<int64_t, int64_t>& gapMap,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapRunBacktracks,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapRunBlocksBacktracks,
    const std::unordered_set<int64_t>& inverseBlockIds,
    const std::vector<std::pair<int64_t, int64_t>>& blockRanges,
    mutableTreeData& data,
    std::map<uint32_t, seeding::onSeedsHash>& onSeedsHashMap,
    std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>& seedChanges,
    Tree* T,
    const int& seedK,
    const globalCoords_t& globalCoords,
    CoordNavigator& navigator,
    const size_t& num_cpus
  ) {
    if (false) {
      // More than 1 CPU: Use task group to run tasks in parallel
      tbb::task_group tg;
      tg.run([&] {
        processGapMutations(
          perNodeGapMutations_Index, dfsIndex, gapMap, gapRunBacktracks, 
          gapRunBlocksBacktracks, inverseBlockIds, blockRanges, data
        );
      });
  
      tg.run([&] {
        processSeedMutations(
          perNodeSeedMutations_Index, dfsIndex, onSeedsHashMap, seedChanges, 
          T, seedK, globalCoords, navigator, data, gapMap, blockRanges
        );
      });
      try {
        tg.wait();
      } catch (const std::exception& e) {
        // Handle any exception thrown by tasks
        std::cerr << "Exception occurred: " << e.what() << std::endl;
      }

    } else {
      // Only 1 CPU: Run tasks sequentially to avoid task group overhead
      processGapMutations(
        perNodeGapMutations_Index, dfsIndex, gapMap, gapRunBacktracks, 
        gapRunBlocksBacktracks, inverseBlockIds, blockRanges, data
      );

      processSeedMutations(
        perNodeSeedMutations_Index, dfsIndex, onSeedsHashMap, seedChanges, 
        T, seedK, globalCoords, navigator, data, gapMap, blockRanges
      );
    }
  }

  void updateSeedsMapAndBlocks(const std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>& seedChanges, std::map<uint32_t, seeding::onSeedsHash>& onSeedsHashMap, const std::vector<int64_t>& scalarCoordToBlockId, std::vector<std::unordered_set<int>>& BlocksToSeeds) {
    for (const auto& change : seedChanges) {
      const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse, oldEndPos, newEndPos] = change;
      auto seedIt = onSeedsHashMap.find(pos);
      if (oldVal && newVal) { // seed at same pos changed
        seedIt->second.hash      = newSeed.value();
        seedIt->second.endPos    = newEndPos.value();
        seedIt->second.isReverse = newIsReverse.value();
      } else if (oldVal && !newVal) { // seed on to off
        onSeedsHashMap.erase(seedIt);
        int blockId = scalarCoordToBlockId[pos];
        BlocksToSeeds[blockId].erase(pos);
      } else if (!oldVal && newVal) { // seed off to on
        onSeedsHashMap[pos] = {newSeed.value(), newEndPos.value(), newIsReverse.value()};
        int blockId = scalarCoordToBlockId[pos];
        BlocksToSeeds[blockId].insert(pos);
      }
    }
  }

  void updateSeedmersIndex(const std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>& seedChanges, 
                          std::map<uint32_t, seeding::onSeedsHash>& onSeedsHashMap,
                          mgsr::seedmers& seedmersIndex,
                          std::unordered_set<size_t>& affectedSeedmers,
                          const int& seedK,
                          const int& seedL,
                          std::vector<std::tuple<int32_t, int32_t, size_t, size_t, bool>>& backTrackPositionMapChAdd,
                          std::vector<int32_t>& backTrackPositionMapErase
  ) {
    auto& positionMap = seedmersIndex.positionMap;
    auto& hashToPositionsMap = seedmersIndex.hashToPositionsMap;



    if (onSeedsHashMap.size() < seedL) {
      if (positionMap.empty()) return;
      auto positionMapIt = positionMap.begin();
      while (positionMapIt != positionMap.end()) {
        const auto& [toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev] = positionMapIt->second;
        backTrackPositionMapChAdd.emplace_back(std::make_tuple(positionMapIt->first, toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev));
        seedmersIndex.delPosition(positionMapIt, affectedSeedmers);
        ++positionMapIt;
      }
      return;
    }

    int64_t maxBegCoord = std::prev(onSeedsHashMap.end(), seedL)->first;
    std::unordered_set<int64_t> processedSeedBegs;
    
    for (const auto& change : seedChanges) {
      const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse, oldEndPos, newEndPos] = change;
      uint8_t leftStepsTaken = 0;
      std::map<uint32_t, seeding::onSeedsHash>::iterator changedSeedIt      = onSeedsHashMap.lower_bound(pos);
      std::map<uint32_t, seeding::onSeedsHash>::iterator firstKminmerSeedIt = changedSeedIt;
      std::map<uint32_t, seeding::onSeedsHash>::iterator lastKminmerSeedIt  = changedSeedIt;
      while (leftStepsTaken < seedL - 1 && firstKminmerSeedIt != onSeedsHashMap.begin()) {
        --firstKminmerSeedIt;
        ++leftStepsTaken;
        if (processedSeedBegs.find(firstKminmerSeedIt->first) != processedSeedBegs.end()) {
          ++firstKminmerSeedIt;
          --leftStepsTaken;
          break;
        }
      }
      lastKminmerSeedIt = safe_next(changedSeedIt, onSeedsHashMap.end(), seedL - leftStepsTaken - 1);
      auto curKminmerPositionIt = positionMap.lower_bound(firstKminmerSeedIt->first);
      while (lastKminmerSeedIt != onSeedsHashMap.end()) {
        if (!newVal && firstKminmerSeedIt == changedSeedIt) break;
        bool isReplacement = curKminmerPositionIt == positionMap.end() ? false : curKminmerPositionIt->first == firstKminmerSeedIt->first;

        if (!positionMap.empty() && positionMap.begin()->first <= firstKminmerSeedIt->first && firstKminmerSeedIt != onSeedsHashMap.begin()) {
          const auto& [prevEnd, prevForwardKminmerHash, prevReverseKminmerHash, prevRev] = std::prev(curKminmerPositionIt)->second;

          size_t curforwardHash = rol(prevForwardKminmerHash, seedK) ^ rol(std::prev(firstKminmerSeedIt)->second.hash, seedK * seedL) ^ lastKminmerSeedIt->second.hash;
          size_t curReverseHash = ror(prevReverseKminmerHash, seedK) ^ ror(std::prev(firstKminmerSeedIt)->second.hash, seedK)         ^ rol(lastKminmerSeedIt->second.hash, seedK * (seedL-1));        
          
          if (isReplacement) {
            const auto& [oldEnd, oldFHash, oldRHash, oldRev] = curKminmerPositionIt->second;
            backTrackPositionMapChAdd.emplace_back(std::make_tuple(curKminmerPositionIt->first, oldEnd, oldFHash, oldRHash, oldRev));
            seedmersIndex.subPosition(curKminmerPositionIt, lastKminmerSeedIt->second.endPos, curforwardHash, curReverseHash, curReverseHash < curforwardHash, affectedSeedmers);
          } else {
            backTrackPositionMapErase.emplace_back(firstKminmerSeedIt->first);
            seedmersIndex.addPosition(firstKminmerSeedIt->first, lastKminmerSeedIt->second.endPos, curforwardHash, curReverseHash, curReverseHash < curforwardHash, affectedSeedmers);
          }

          processedSeedBegs.insert(firstKminmerSeedIt->first);
        } else {
          auto curItForward = firstKminmerSeedIt;
          auto curItReverse = lastKminmerSeedIt;
          size_t forwardKminmerHash = 0;
          size_t reverseKminmerHash = 0;
          while (true) {
            forwardKminmerHash = rol(forwardKminmerHash, seedK) ^ curItForward->second.hash;
            reverseKminmerHash = rol(reverseKminmerHash, seedK) ^ curItReverse->second.hash;
            if (curItForward == lastKminmerSeedIt) break;
            ++curItForward;
            --curItReverse;
          }
          
          if (isReplacement) {
            const auto& [oldEnd, oldFHash, oldRHash, oldRev] = curKminmerPositionIt->second;
            backTrackPositionMapChAdd.emplace_back(std::make_tuple(curKminmerPositionIt->first, oldEnd, oldFHash, oldRHash, oldRev));
            seedmersIndex.subPosition(curKminmerPositionIt, lastKminmerSeedIt->second.endPos, forwardKminmerHash, reverseKminmerHash, reverseKminmerHash < forwardKminmerHash, affectedSeedmers);
          } else {
            backTrackPositionMapErase.emplace_back(firstKminmerSeedIt->first);
            seedmersIndex.addPosition(firstKminmerSeedIt->first, lastKminmerSeedIt->second.endPos, forwardKminmerHash, reverseKminmerHash, reverseKminmerHash < forwardKminmerHash, affectedSeedmers);
          }
          
          processedSeedBegs.insert(firstKminmerSeedIt->first);
        }

        if (firstKminmerSeedIt == changedSeedIt) break;
        ++firstKminmerSeedIt;
        ++lastKminmerSeedIt;
        if (curKminmerPositionIt != positionMap.end()) ++curKminmerPositionIt;
      }

      if (!newVal) {
        auto toEraseIt = positionMap.find(pos);
        if (toEraseIt != positionMap.end()) {
          const auto& [toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev] = toEraseIt->second;
          backTrackPositionMapChAdd.emplace_back(std::make_tuple(toEraseIt->first, toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev));
          processedSeedBegs.insert(toEraseIt->first);
          seedmersIndex.delPosition(toEraseIt, affectedSeedmers);
        }
      }
    }

    if (positionMap.empty()) return;
    auto lastPositionMapIt = std::prev(positionMap.end());
    while (lastPositionMapIt->first > maxBegCoord) {
      const auto& [toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev] = lastPositionMapIt->second;
      backTrackPositionMapChAdd.emplace_back(std::make_tuple(lastPositionMapIt->first, toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev));
      seedmersIndex.delPosition(lastPositionMapIt, affectedSeedmers);
      --lastPositionMapIt;
    }
  }

  int32_t extend(
    int64_t& curEnd, mgsr::Read& curRead, int rev, std::map<int32_t, mgsr::positionInfo>& positionMap,
    std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap,
    int32_t qidx, std::map<int32_t, mgsr::positionInfo>::const_iterator refPositionIt, int32_t c
  ) {
    if (qidx == curRead.seedmersList.size() - 1) return c;
    const auto& [nhash, nqbeg, nqend, nqrev, nqidx] = curRead.seedmersList[qidx+1];

    auto nextHashToPositionIt = hashToPositionsMap.find(nhash);
    if (nextHashToPositionIt != hashToPositionsMap.end()) {
      if (nextHashToPositionIt->second.size() < 2) {
        if (nextHashToPositionIt->second.size() != 1) {
          std::cerr << "Error: hash " << nhash << " has multiple positions" << std::endl;
          exit(1);
        }

        auto curRefPositionIt = *nextHashToPositionIt->second.begin();
        char nrev = nqrev == curRefPositionIt->second.rev? 1 : 2;
        if (rev == nrev) {
          if (rev == 1) {
            auto nextRefPositionIt = std::next(refPositionIt);
            while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
            if (curRefPositionIt->first == nextRefPositionIt->first) {
              ++c;
              curEnd = nqidx;
              return extend(curEnd, curRead, rev, positionMap, hashToPositionsMap, nqidx, curRefPositionIt, c);
            }
          } else if (rev == 2) {
            auto prevRefPositionIt = std::prev(refPositionIt);
            while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
            if (curRefPositionIt->first == prevRefPositionIt->first) {
              ++c;
              curEnd = nqidx;
              return extend(curEnd, curRead, rev, positionMap, hashToPositionsMap, nqidx, curRefPositionIt, c);
            }
          }
        }
      }
    }
    return c;
  }

  void initializeMatches(mgsr::Read& curRead, std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap) {
    int32_t i = 0;
    while (i < curRead.seedmersList.size()) {
      const auto& [hash, qbeg, qend, qrev, qidx] = curRead.seedmersList[i];
      int32_t c = 1;
      auto hashToPositionIt = hashToPositionsMap.find(hash);
      if (hashToPositionIt != hashToPositionsMap.end()) {
        if (hashToPositionIt->second.size() < 2) {
          if (hashToPositionIt->second.size() != 1) {
            std::cerr << "Error: hash " << hash << " has multiple positions" << std::endl;
            exit(1);
          }
          auto curRefPositionIt = *(hashToPositionIt->second.begin());
          int64_t curEnd = i;
          int rev = qrev == curRefPositionIt->second.rev? 1 : 2;
          c = extend(curEnd, curRead, rev, positionMap, hashToPositionsMap, qidx, curRefPositionIt, c);
          curRead.matches.add({discrete_interval<int32_t>::closed(i, curEnd), rev});
        } else {
          curRead.duplicates.insert(i);
        }
      }
      i += c; 
    }
  }

  bool isColinear(const std::pair<boost::icl::discrete_interval<int32_t>, int>& match1, const std::pair<boost::icl::discrete_interval<int32_t>, int>& match2, const mgsr::Read& curRead, mgsr::seedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex, const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap, const int& dfsIndex) {
    bool rev1 = match1.second == 1 ? false : true;
    bool rev2 = match2.second == 1 ? false : true;
    if (rev1 != rev2) {
      std::cerr << "Error: Invalid direction in isColinear_test()" << std::endl;
      exit(1);
    }

    const auto& first1 = curRead.seedmersList[boost::icl::first(match1.first)];
    const auto& first2 = curRead.seedmersList[boost::icl::first(match2.first)];
    const auto& last1  = curRead.seedmersList[boost::icl::last(match1.first)];
    const auto& last2  = curRead.seedmersList[boost::icl::last(match2.first)];

    if (rev1 == false) {
      // forward direction
      const auto& rglobalbeg1 = seedmersIndex.getBegFromHash(first1.hash);
      const auto& rglobalend1 = seedmersIndex.getEndFromHash(last1.hash);
      const auto& rglobalbeg2 = seedmersIndex.getBegFromHash(first2.hash);
      // const auto& rglobalend2 = (positionMap.find(*(hashToPositionsMap.find(*last2.hash)->second.begin()))->second).endPos;
      
      const auto& qbeg1 = first1.begPos;
      const auto& qend1 = last1.endPos;
      const auto& qbeg2 = first2.begPos;
      const auto& qend2 = last2.endPos;

      auto rbeg1 = degapGlobal(rglobalbeg1, degapCoordIndex);
      auto rend1 = degapGlobal(rglobalend1, degapCoordIndex);
      auto rbeg2 = degapGlobal(rglobalbeg2, degapCoordIndex);
      // auto rend2 = degapGlobal(rglobalend2, degapCoordIndex);

      int32_t qgap = abs(qbeg2 - qend1);
      int32_t rgap = abs(rbeg2 - rend1);
      if (rbeg1 < rbeg2 && abs(qgap - rgap) < maximumGap) return true;

    } else {
      // reverse direction
      auto rglobalbeg1 = seedmersIndex.getBegFromHash(last1.hash);
      // auto rglobalend1 = (positionMap.find(*(hashToPositionsMap.find(*first1.hash)->second.begin()))->second).endPos;
      auto rglobalbeg2 = seedmersIndex.getBegFromHash(last2.hash);
      auto rglobalend2 = seedmersIndex.getEndFromHash(first2.hash);

      const auto& qbeg1 = first1.begPos;
      const auto& qend1 = last1.endPos;
      const auto& qbeg2 = first2.begPos;
      const auto& qend2 = last2.endPos;

      auto rbeg1 = degapGlobal(rglobalbeg1, degapCoordIndex);
      // auto rend1 = degapGlobal(rglobalend1, degapCoordIndex);
      auto rbeg2 = degapGlobal(rglobalbeg2, degapCoordIndex);
      auto rend2 = degapGlobal(rglobalend2, degapCoordIndex);

      int32_t qgap = abs(qbeg2 - qend1);
      int32_t rgap = abs(rbeg1 - rend2);
      if (rbeg2 < rbeg1 && abs(qgap - rgap) < maximumGap) return true;
    }

    return false;
  }

  int64_t getPseudoScore(
    const mgsr::Read& curRead, mgsr::seedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex,
    const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap, const int& minimumCount, const int& minimumScore,
    const bool& rescueDuplicates, const double& rescueDuplicatesThreshold, const int& dfsIndex
  ) {
    int64_t pseudoScore = 0;
    const boost::icl::discrete_interval<int32_t>* firstMatch = nullptr;
    const boost::icl::discrete_interval<int32_t>* lastMatch  = nullptr;
    bool reversed = false;

    // simple cases
    if (curRead.matches.empty()) {
      return 0;
    } else if (curRead.matches.size() == 1) {
      pseudoScore = boost::icl::length(curRead.matches.begin()->first);
      firstMatch  = &curRead.matches.begin()->first;
      lastMatch   = &curRead.matches.begin()->first;
      reversed    = curRead.matches.begin()->second == 1 ? false : true;
    } else {
      // find longest interval
      std::pair<boost::icl::discrete_interval<int32_t>, int> longestInterval = *curRead.matches.begin();
      int longestIdx = 0;
      int curIdx = 1;
      for (auto it = std::next(curRead.matches.begin()); it != curRead.matches.end(); ++it) {
        const auto& curInterval = it;
        if (boost::icl::length(curInterval->first) > boost::icl::length(longestInterval.first)) {
          longestInterval = *curInterval;
          longestIdx = curIdx;
        }
        ++curIdx;
      }
      reversed = longestInterval.second == 1 ? false : true;

      auto longestQbeg = curRead.seedmersList[boost::icl::first(longestInterval.first)].begPos;
      curIdx = 0;
      for (const auto& curInterval : curRead.matches) {
        if (curIdx == longestIdx) {
          pseudoScore += boost::icl::length(curInterval.first);
          if (firstMatch == nullptr) firstMatch = &curInterval.first;
          lastMatch = &curInterval.first;
          ++curIdx;
          continue;
        }

        if (curInterval.second != longestInterval.second) {
          ++curIdx;
          continue;
        }

        auto curQbeg = curRead.seedmersList[boost::icl::first(curInterval.first)].begPos;
        if (longestQbeg < curQbeg) {
          if (isColinear(longestInterval, curInterval, curRead, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap, dfsIndex)) {
            pseudoScore += boost::icl::length(curInterval.first);
            if (firstMatch == nullptr) firstMatch = &curInterval.first;
            lastMatch = &curInterval.first;
          }
        } else if (longestQbeg > curQbeg) {
          if (isColinear(curInterval, longestInterval, curRead, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap, dfsIndex)) {
            pseudoScore += boost::icl::length(curInterval.first);
            if (firstMatch == nullptr) firstMatch = &curInterval.first;
            lastMatch = &curInterval.first;
          }
        } else {
          std::cerr << "Error: Invalid direction in getPseudoScore()" << std::endl;
          exit(1);
        }

        ++curIdx;
      }
    }

    const auto& curReadDuplicates = curRead.duplicates;
    if (rescueDuplicates && !curReadDuplicates.empty() && curReadDuplicates.size() <= rescueDuplicatesThreshold * curRead.seedmersList.size()) {
      const size_t& leftBoundHash  = reversed ? curRead.seedmersList[boost::icl::last(*lastMatch)].hash : curRead.seedmersList[boost::icl::first(*firstMatch)].hash;
      const size_t& rightBoundHash = reversed ? curRead.seedmersList[boost::icl::first(*firstMatch)].hash : curRead.seedmersList[boost::icl::last(*lastMatch)].hash;

      int32_t leftBoundGlobal = seedmersIndex.getBegFromHash(leftBoundHash);
      int32_t rightBoundGlobal = seedmersIndex.getEndFromHash(rightBoundHash);

      if (leftBoundGlobal >= rightBoundGlobal) {
        std::cerr << "Error: Invalid bounds in getPseudoScore()" << std::endl;
        exit(1);
      }
      int32_t leftBoundLocal = std::max(static_cast<int64_t>(0), degapGlobal(leftBoundGlobal, degapCoordIndex) - 150);
      int32_t rightBoundLocal = degapGlobal(rightBoundGlobal, degapCoordIndex) + 150;

      leftBoundGlobal = regapGlobal(leftBoundLocal, regapCoordIndex);
      rightBoundGlobal = regapGlobal(rightBoundLocal, regapCoordIndex);

      for (const auto& duplicateIndex : curReadDuplicates) {
        const auto& duplicateHash = curRead.seedmersList[duplicateIndex].hash;
        const auto& duplicateRev = curRead.seedmersList[duplicateIndex].rev;
        auto hashToPositionIt = seedmersIndex.hashToPositionsMap.find(duplicateHash);
        if (hashToPositionIt == seedmersIndex.hashToPositionsMap.end()) continue;
        for (const auto& positionMapIt : hashToPositionIt->second) {
          const auto& rbeg = positionMapIt->first;
          const auto& rend = positionMapIt->second.endPos;
          const auto& rrev = positionMapIt->second.rev;
          bool curReversedMatch = rrev != duplicateRev;
          if (rbeg >= leftBoundGlobal && rend <= rightBoundGlobal && curReversedMatch == reversed) {
            ++pseudoScore;
            break;
          }
        }
      }
    }
    return pseudoScore;
  }

  std::pair<bool, std::vector<size_t>> checkRedo(
    const std::unordered_map<size_t, std::vector<int32_t>>& a, const std::unordered_set<size_t>& b, const int& threshold
  ) {
    std::vector<size_t> affectedSeedmers;
    uint32_t numAffected = 0;
    if (a.size() < b.size()) {
      for (const auto& [h, _] : a) {
        if (b.find(h) != b.end()) {
          affectedSeedmers.push_back(h);
          ++numAffected;
          if (numAffected > threshold) return {true, {}};
        }
      }
    } else {
      for (const auto& h : b) {
        if (a.find(h) != a.end()) {
          affectedSeedmers.push_back(h);
          ++numAffected;
          if (numAffected > threshold) return {true, {}};
        }
      }
    }

    return {false, affectedSeedmers};
  }

  bool isContained(const boost::icl::split_interval_map<int32_t, int>& matches, int32_t number) {
    auto it = matches.find(number);
    if (it != matches.end()) {
      return true;
    }
    return false;
  }


  double getExp(const Eigen::MatrixXd& probs, const Eigen::VectorXd& props, const Eigen::VectorXd& numReadDuplicates) {
      assert(props.size() == probs.cols());

      Eigen::VectorXd readSums = probs * props;
      double llh = (numReadDuplicates.array() * readSums.array().log()).sum();

      return llh;
  }


  Eigen::VectorXd getMax(const Eigen::MatrixXd& probs, const Eigen::VectorXd& props, const Eigen::VectorXd& numReadDuplicates, const int32_t& totalReads) {
    size_t numNodes = probs.cols();

    // std::cerr << "calculating demons in getMax" << std::endl;
    Eigen::VectorXd denoms = probs * props;
    // std::cerr << "setting up newProps "<< std::endl;
    Eigen::VectorXd newProps(numNodes);
    newProps.setZero();
    size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

    // std::cerr << "starting parallel for in getMax" << std::endl;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numNodes, numNodes/num_cpus),
      [&](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i < range.end(); ++i) {
          Eigen::VectorXd ratios = (probs.col(i).array() * props[i]) / denoms.array();
          double newProp = (numReadDuplicates.array() * ratios.array()).sum();
          newProp /= totalReads;
          newProps(i) = newProp;
        }
      }); 


    return newProps;
  }


  void normalize(Eigen::VectorXd& props) {    
    for (int i = 0; i < props.size(); ++i) {
      if (props(i) <= 0) {
          props(i) = std::numeric_limits<double>::min();
      }
    }
    double sum = props.sum();
    props /= sum;
  }

  void updateInsigCounts(const Eigen::VectorXd& props, std::vector<int>& insigCounts, const double& insigProp, size_t totalNodes) {
    for (int i = 0; i < props.size(); ++i) {
      if (props(i) <= insigProp) {
        ++insigCounts[i];
      } else {
        insigCounts[i] = 0;
      }
    }
  }


  void squarem_test_1(
    const std::vector<std::string>& nodes, const Eigen::MatrixXd& probs,
    const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    const Eigen::VectorXd& numReadDuplicates, const int32_t& numReads, 
    Eigen::VectorXd& props, double& llh, int& curit, bool& converged,
    const size_t& iterations, const size_t& prefilterIterations,
    std::vector<int>& insigCounts, const double& insigProp, size_t totalNodes
  ) {
    assert(nodes.size() == probs.cols());
    assert(nodes.size() == props.size());
    size_t curIteration = 1;
    while (true) {
      std::cout << "\riteration " << curIteration << std::flush;
      // std::cerr << "getting max at theta1" << std::endl;
      Eigen::VectorXd theta1 = getMax(probs, props, numReadDuplicates, numReads);
      // std::cerr << "normalizing theta1" << std::endl;
      normalize(theta1);
      // std::cerr << "getting max at theta2" << std::endl;
      Eigen::VectorXd theta2 = getMax(probs, theta1, numReadDuplicates, numReads);
      // std::cerr << "normalizing theta2" << std::endl;
      normalize(theta2);

      // std::cerr << "calculating r_norm and v_norm" << std::endl;
      Eigen::VectorXd r = theta1 - props;
      Eigen::VectorXd v = theta2 - theta1 - r;
      double r_norm = r.norm();
      double v_norm = v.norm();

      double alpha;
      if (r_norm == 0 || v_norm == 0) {
        alpha = 0;
      } else {
        alpha = - r_norm / v_norm;
      }
      double newllh;

      
      Eigen::VectorXd theta_p;
      if (alpha > -1) {
        alpha = -1;
        // std::cerr << "calculating theta_p alpha > -1" << std::endl;
        theta_p = props - 2 * alpha * r + alpha * alpha * v;
        // std::cerr << "getting max at theta_p" << std::endl;
        props = getMax(probs, theta_p, numReadDuplicates, numReads);
        // std::cerr << "normalizing props" << std::endl;
        normalize(props);
        // std::cerr << "calculating newllh" << std::endl;
        newllh = getExp(probs, props, numReadDuplicates);
      } else {
        // std::cerr << "calculating theta_p alpha <= -1" << std::endl;
        theta_p = props - 2 * alpha * r + alpha * alpha * v;
        // std::cerr << "getting max at theta_p" << std::endl;
        auto newProps = getMax(probs, theta_p, numReadDuplicates, numReads);
        // std::cerr << "normalizing newProps" << std::endl;
        normalize(newProps);
        // std::cerr << "calculating newllh" << std::endl;
        newllh = getExp(probs, newProps, numReadDuplicates);
        if (newllh >= llh) {
          props = std::move(newProps);
        } else {
          while (llh - newllh > 0.0001) {
            // std::cerr << "it " << curit << "\t" << newllh << "\t" << llh << std::endl;
            alpha = (alpha - 1) / 2;
            // std::cerr << "calculating theta_p alpha <= -1 inside while loop" << std::endl;
            theta_p = props - 2 * alpha * r + alpha * alpha * v;
            // std::cerr << "getting max at theta_p" << std::endl;
            newProps = getMax(probs, theta_p, numReadDuplicates, numReads);
            // std::cerr << "normalizing newProps" << std::endl;
            normalize(newProps);
            // std::cerr << "calculating newllh" << std::endl;
            newllh   = getExp(probs, newProps, numReadDuplicates);
          }
          props = std::move(newProps);
        }
      }
      // std::cerr << "updating insigCounts" << std::endl;
      if (iterations != std::numeric_limits<size_t>::max()) {
        updateInsigCounts(props, insigCounts, insigProp, totalNodes);
      }
      // std::cerr << "it " << curit << "\t" << newllh << "\t" << llh << std::endl;
      if (newllh - llh < 0.0001) {
        llh = newllh;
        converged = true;
        break;
      } else if (curIteration == iterations) {
        llh = newllh;
        break;
      } else if (curIteration == prefilterIterations) {
        llh = newllh;
        break;
      }

      llh = newllh;
      ++curit;
      ++curIteration;
    }
    ++curit;
  }

  //squarem test 1: periodically drop nodes with very low abundance
  void squaremHelper_test_1(
    Tree *T, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    const std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex, const std::vector<bool>& lowScoreReads,
    const int32_t& numReads, const size_t& numLowScoreReads, const std::vector<bool>& excludeReads,
    const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors,
    const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets, Eigen::MatrixXd& probs,
    std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const std::string& preEMFilterMethod, const int& preEMFilterNOrder,
    const int& emFilterRound, const int& checkFrequency, const int& removeIteration, const double& insigProp,
    const int& roundsRemove, const double& removeThreshold, std::string excludeNode
  ) {
    if (excludeNode.empty()) {
      std::stringstream msg;
      msg << "starting to set up EM" << "\n";
      std::cout << msg.str();
    } else {
      std::stringstream msg;
      msg << "starting to set up EM excluding " << excludeNode << "\n";
      std::cout << msg.str();
    }

    std::cout << "pre-EM filter nodes size: " << allScores.size() - leastRecentIdenticalAncestors.size() << std::endl;
    std::cerr << "pre-EM filter nodes size: " << allScores.size() - leastRecentIdenticalAncestors.size() << "\n" << std::endl;
    if (preEMFilterMethod == "null") {  
      haplotype_filter::noFilter(nodes, probs, allScores, leastRecentIdenticalAncestors, lowScoreReads, numLowScoreReads, excludeNode, excludeReads);
    } else if (preEMFilterMethod == "uhs") {
      haplotype_filter::filter_method_1(nodes, probs, allScores, leastRecentIdenticalAncestors, identicalSets, lowScoreReads, numLowScoreReads, excludeNode, excludeReads, T, preEMFilterNOrder);
    } else if (preEMFilterMethod == "hsc1") {
      haplotype_filter::filter_method_2(nodes, probs, allScores, leastRecentIdenticalAncestors, lowScoreReads, numLowScoreReads, excludeNode, excludeReads, 1);
    } else if (preEMFilterMethod == "hsc2") {
      size_t numExcludedReads = std::count(excludeReads.begin(), excludeReads.end(), true);
      haplotype_filter::filter_method_2(nodes, probs, allScores, leastRecentIdenticalAncestors, lowScoreReads, numLowScoreReads, excludeNode, excludeReads, static_cast<size_t>(static_cast<double>(allScores.begin()->second.size() - numExcludedReads) * 0.01));
    } else {
      std::cerr << "pre-EM filter method not recognized" << std::endl;
      exit(1);
    }
    std::string filtedNodesFile = "filted_nodes.txt";
    std::ofstream filtedNodesStream(filtedNodesFile);
    for (const auto& node : nodes) {
      if (leastRecentIdenticalAncestors.find(node) != leastRecentIdenticalAncestors.end()) { 
        // sanity check
        std::cerr << "Error: Node " << node << " has a least recent identical ancestor." << std::endl;
        exit(1);
      }
      filtedNodesStream << node << std::endl;
    }
    std::cout << "post-EM filter nodes size: " << nodes.size() << std::endl;
    std::cerr << "post-EM filter nodes size: " << nodes.size() << "\n" << std::endl;
    std::cout << "post-EM filter reduced nodes size: " << allScores.size() - leastRecentIdenticalAncestors.size() - nodes.size() << std::endl;
    std::cerr << "post-EM filter reduced nodes size: " << allScores.size() - leastRecentIdenticalAncestors.size() - nodes.size() << "\n" << std::endl;

    props = Eigen::VectorXd::Constant(nodes.size(), 1.0 / static_cast<double>(nodes.size()));
    size_t totalNodes = allScores.size() - leastRecentIdenticalAncestors.size();
    size_t numExcludedReads = std::count(excludeReads.begin(), excludeReads.end(), true);
    Eigen::VectorXd readDuplicates(allScores.begin()->second.size() - numLowScoreReads - numExcludedReads);
    
    size_t indexReadDuplicates = 0;
    size_t numHighScoreReads = 0;
    for (size_t i = 0; i < readSeedmersDuplicatesIndex.size(); ++i) {
      if (excludeReads[i]) continue;
      if (!lowScoreReads[i]) {
        readDuplicates(indexReadDuplicates) = readSeedmersDuplicatesIndex[i].size();
        numHighScoreReads += readSeedmersDuplicatesIndex[i].size();
        ++indexReadDuplicates;
      }
    }

    std::cout << "readDuplicates size: " << readDuplicates.size() << std::endl;
    std::cout << "probs size: " << probs.rows() << " x " << probs.cols() << std::endl;

    if (excludeNode.empty()) {
      std::stringstream msg;
      msg << "starting EM estimation of haplotype proportions" << "\n";
      std::cout << msg.str();
      std::cerr << msg.str();
    } else {
      std::stringstream msg;
      msg << "starting EM estimation of haplotype proportions excluding " << excludeNode << "\n";
      std::cout << msg.str();
      std::cerr << msg.str();
    }


    int curit = 0;
    bool converged = false;
    int32_t filterRoundCount = 0;
    size_t prefilterIterations = 5;
    std::vector<int> insigCounts(nodes.size());
    std::cout << "\npre-filter round for " << prefilterIterations << " iterations" << std::endl;
    std::cerr << "\npre-filter round for " << prefilterIterations << " iterations" << std::endl;
    llh = getExp(probs, props, readDuplicates);
    squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, checkFrequency, prefilterIterations, insigCounts, insigProp, totalNodes);
    while (true && nodes.size() > std::max(static_cast<int>(totalNodes) / 20, 100)) {
      if (filterRoundCount >= emFilterRound) break;
      std::cout << "\nfilter round " << filterRoundCount + 1 << " out of " << emFilterRound << std::endl;
      std::cerr << "\nfilter round " << filterRoundCount + 1 << " out of " << emFilterRound << std::endl;
      ++filterRoundCount;
      llh = getExp(probs, props, readDuplicates);
      squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, checkFrequency, std::numeric_limits<size_t>::max(), insigCounts, insigProp, totalNodes);
      if (converged) {
        break;
      }
      std::cout << "\nfiltering round " << filterRoundCount << std::endl;
      std::cerr << "\nfiltering round " << filterRoundCount << std::endl;
      std::vector<size_t> significantIndices;
      std::vector<std::string> sigNodes;
      for (size_t i = 0; i < nodes.size(); ++i) {
        if (insigCounts[i] < removeIteration) {
          significantIndices.push_back(i);
        }
      }

      if (significantIndices.size() == nodes.size()) {
        continue;
      }

      for (size_t idx : significantIndices) {
        sigNodes.push_back(nodes[idx]);
      }

      Eigen::MatrixXd sigProbs(probs.rows(), sigNodes.size());
      Eigen::VectorXd sigProps(sigNodes.size());
      for (size_t i = 0; i < significantIndices.size(); ++i) {
        sigProbs.col(i) = probs.col(significantIndices[i]);
        sigProps(i) = props(i);
      }
      std::cout << "dropped " << nodes.size() - sigNodes.size() << " during EM" << std::endl;
      std::cerr << "dropped " << nodes.size() - sigNodes.size() << " during EM" << std::endl;
      std::cout << sigNodes.size() << " nodes left" << std::endl;
      std::cerr << sigNodes.size() << " nodes left" << std::endl;
      nodes = std::move(sigNodes);
      probs = std::move(sigProbs);
      props = std::move(sigProps);
      normalize(props);
      insigCounts.assign(nodes.size(), 0);
      if (nodes.size() <= std::max(static_cast<int>(totalNodes) / 20, 100) || filterRoundCount >= emFilterRound) {
        break;
      }
    }

    if (!converged) {
      std::vector<int> insigCounts(nodes.size());
      std::cout << "start full EM" << std::endl;
      std::cerr << "start full EM" << std::endl;
      llh = getExp(probs, props, readDuplicates);
      squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), insigCounts, insigProp, totalNodes);
      assert(converged == true);
    }

    for (int32_t i = 0; i < roundsRemove; ++i) {
      std::vector<size_t> significantIndices;
      std::vector<std::string> sigNodes;

      for (size_t i = 0; i < props.size(); ++i) {
          if (props(i) >= removeThreshold) {
              significantIndices.push_back(i);
          }
      }
      if (significantIndices.size() == nodes.size()) break;
      std::cout << "\nremove round " << i + 1 << std::endl;
      std::cerr << "\nremove round " << i + 1 << std::endl;
      for (size_t idx : significantIndices) {
          sigNodes.push_back(nodes[idx]);
      }

      Eigen::MatrixXd sigProbs(probs.rows(), significantIndices.size());
      sigProbs.resize(probs.rows(), significantIndices.size());
      for (size_t i = 0; i < significantIndices.size(); ++i) {
          sigProbs.col(i) = probs.col(significantIndices[i]);
      }
      Eigen::VectorXd sigProps = Eigen::VectorXd::Constant(sigNodes.size(), 1.0 / static_cast<double>(sigNodes.size()));
      llh = getExp(sigProbs, sigProps, readDuplicates);
      bool converged = false;
      size_t iterations = std::numeric_limits<size_t>::max();
      std::vector<int> insigCounts(sigNodes.size());
      squarem_test_1(sigNodes, sigProbs, identicalSets, readDuplicates, numHighScoreReads, sigProps, llh, curit, converged, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), insigCounts, insigProp, totalNodes);
      assert(converged);
      nodes = std::move(sigNodes);
      probs = std::move(sigProbs);
      props = std::move(sigProps);
    }

    if (!excludeNode.empty()) {
      Eigen::VectorXd curProbs(allScores.begin()->second.size() - numLowScoreReads);
      size_t indexCurProbs = 0;
      const auto& curNode = allScores.at(excludeNode);
      for (size_t i = 0; i < curNode.size(); ++i) {
        if (!lowScoreReads[i]) {
          const auto& score = curNode[i];
          curProbs(indexCurProbs) = score.second;
          ++indexCurProbs;
        }
      }

      probs = probs, curProbs;
      props = props, 0.0;
      llh = getExp(probs, props, readDuplicates);
    }

    if (excludeNode.empty()) {
      std::stringstream msg;
      msg << "\nFinished EM estimation of haplotype proportions. Total EM iterations: " << curit << "\n";
      std::cout << msg.str();
      std::cerr << msg.str();
    } else {
      std::stringstream msg;
      msg << "\nFinished EM estimation of haplotype proportions excluding " << excludeNode << ". Total EM iterations: " << curit<< "\n";
      std::cout << msg.str();
      std::cerr << msg.str();
    }
  }

  void assignReadsToNodes(
    const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    const std::vector<std::string>& nodes, const Eigen::MatrixXd& probs, const Eigen::VectorXd& props,
    const std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex, std::unordered_map<std::string, std::vector<size_t>>& assignedReads
  ) {
    size_t rowindex = 0;
    for (size_t i = 0; i < readSeedmersDuplicatesIndex.size(); ++i) {
      const Eigen::VectorXd& curProbs = probs.row(rowindex);
      double curmax = curProbs.maxCoeff();
      for (size_t j = 0; j < curProbs.size(); ++j) {
        if (curProbs(j) == curmax) {
          for (size_t k = 0; k < readSeedmersDuplicatesIndex[i].size(); ++k) {
            assignedReads[nodes[j]].push_back(readSeedmersDuplicatesIndex[i][k]);
          }
        }
      }
      ++rowindex;
    }
  }
}

#endif

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
#include <eigen3/Eigen/Dense>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>


using namespace boost::icl;
typedef std::pair<std::vector<std::tuple<size_t*, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> readSeedmers_t;
typedef std::tuple<int32_t, int32_t, int32_t, int32_t, bool, int32_t> match_t;

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

  template <typename SeedMutationsType, typename GapMutationsType>
  void processNodeMutations(const GapMutationsType& perNodeGapMutations_Index, const SeedMutationsType& perNodeSeedMutations_Index, int64_t dfsIndex, std::map<int64_t, int64_t>& gapMap, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapRunBacktracks, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapRunBlocksBacktracks, const std::unordered_set<int64_t>& inverseBlockIds, const std::vector<std::pair<int64_t, int64_t>>& blockRanges, mutableTreeData& data, std::map<uint32_t, seeding::onSeedsHash>& onSeedsHashMap, std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>>& seedChanges, Tree* T, int seedK, const globalCoords_t& globalCoords, CoordNavigator& navigator) {
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

    auto currBasePositions = perNodeSeedMutations_Index[dfsIndex].getBasePositions();
    auto currPerPosMasks = perNodeSeedMutations_Index[dfsIndex].getPerPosMasks();
    std::vector<std::vector<int8_t>> masks;
    for (int i = 0; i < currBasePositions.size(); ++i) {
      int64_t pos = currBasePositions[i];
      uint64_t tritMask = currPerPosMasks[i];
      for (int k = 0; k < 32; ++k) {
        uint8_t ternaryNumber = (tritMask >> (k * 2)) & 0x3;
        auto seedIt = onSeedsHashMap.find(pos - k);
        if (ternaryNumber == 1) { // on -> off
          // Handle deletion
          const auto& [oldSeed, oldEndPos, oldIsReverse] = seedIt->second;
          seedChanges.emplace_back(std::make_tuple(pos - k, true, false, oldSeed, std::nullopt, oldIsReverse, std::nullopt, oldEndPos, std::nullopt));
        } else if (ternaryNumber == 2) {
          // Handle insertion/change
          auto [newSeed, newSeedEndPos] = seed_annotated_tree::getSeedAt(pos - k, T, seedK, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap, blockRanges);
          auto [newSeedFHash, newSeedRHash] = hashSeq(newSeed);
          size_t newSeedHash;
          bool newIsReverse;
          if (newSeedFHash < newSeedRHash) {
            newSeedHash  = newSeedFHash;
            newIsReverse = false;
          } else {
            newSeedHash  = newSeedRHash;
            newIsReverse = true;
          }

          if (seedIt != onSeedsHashMap.end()) { // on -> on
            const auto& [oldSeed, oldEndPos, oldIsReverse] = seedIt->second;
            seedChanges.emplace_back(std::make_tuple(pos - k, true, true, oldSeed, newSeedHash, oldIsReverse, newIsReverse, oldEndPos, newSeedEndPos));
          } else { // off -> on
            seedChanges.emplace_back(std::make_tuple(pos - k, false, true, std::nullopt, newSeedHash, std::nullopt, newIsReverse, std::nullopt, newSeedEndPos));
          }
        }
      }
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

    auto lastPositionMapIt = std::prev(positionMap.end());
    while (lastPositionMapIt->first > maxBegCoord) {
      const auto& [toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev] = lastPositionMapIt->second;
      backTrackPositionMapChAdd.emplace_back(std::make_tuple(lastPositionMapIt->first, toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev));
      seedmersIndex.delPosition(lastPositionMapIt, affectedSeedmers);
      --lastPositionMapIt;
    }
  }

  int32_t extend(int64_t& curEnd, mgsr::Read& curRead, int rev, std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap, int32_t qidx, std::map<int32_t, mgsr::positionInfo>::const_iterator refPositionIt, int32_t c) {
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

  bool isColinear(const std::pair<boost::icl::discrete_interval<int32_t>, int>& match1, const std::pair<boost::icl::discrete_interval<int32_t>, int>& match2, const mgsr::Read& curRead, mgsr::seedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex, const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap) {
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
    const mgsr::Read& curRead, mgsr::seedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex, const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap, const int& minimumCount, const int& minimumScore
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

      // find intervals colinear with longest interval and add length to pseudoScore
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

        // check if same direction
        if (curInterval.second != longestInterval.second) {
          ++curIdx;
          continue;
        }

        auto curQbeg = curRead.seedmersList[boost::icl::first(curInterval.first)].begPos;
        // check if within maximum gap
        if (longestQbeg < curQbeg) {
          // longest query beg before current query beg
          if (isColinear(longestInterval, curInterval, curRead, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap)) {
            pseudoScore += boost::icl::length(curInterval.first);
            if (firstMatch == nullptr) firstMatch = &curInterval.first;
            lastMatch = &curInterval.first;
          }
        } else if (longestQbeg > curQbeg) {
          // longest query beg after current query beg
          if (isColinear(curInterval, longestInterval, curRead, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap)) {
            pseudoScore += boost::icl::length(curInterval.first);
            if (firstMatch == nullptr) firstMatch = &curInterval.first;
            lastMatch = &curInterval.first;
          }
        } else {
          // error
          std::cerr << "Error: Invalid direction in getPseudoScore()" << std::endl;
          exit(1);
        }

        ++curIdx;
      }
    }

    bool rescueDuplicates = false;
    const auto& curReadDuplicates = curRead.duplicates;
    if (rescueDuplicates && !curReadDuplicates.empty() && curReadDuplicates.size() < 5) {
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

  std::pair<bool, std::vector<size_t>> checkRedo(const std::unordered_map<size_t, std::vector<int32_t>>& a, const std::unordered_set<size_t>& b, int32_t threshold) {
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

    Eigen::VectorXd denoms = probs * props;

    Eigen::VectorXd newProps(numNodes);
    newProps.setZero();
    size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numNodes, numNodes/num_cpus),
      [&](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i < range.end(); ++i) {
          Eigen::VectorXd ratios = (probs.col(i).array() * props[i]) / denoms.array();
          double newProp = (numReadDuplicates.array() * ratios.array()).sum();
          newProp /= totalReads;
          newProps(i) = newProp;
        }
      }); 

    // for (size_t i = 0; i < numNodes; ++i) {
    //   Eigen::VectorXd ratios = (probs.col(i).array() * props[i]) / denoms.array();
    //   double newProp = (numReadDuplicates.array() * ratios.array()).sum();
    //   newProp /= totalReads;
    //   newProps(i) = newProp;
    // }

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

  void updateInsigCounts(const Eigen::VectorXd& props, std::vector<size_t>& insigCounts, size_t totalNodes) {
    double insigProp =  (1.0 / static_cast<double>(totalNodes)) / 10.0;
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
    Eigen::VectorXd& props, double& llh, int& curit, bool& converged, size_t iterations, std::vector<size_t>& insigCounts, size_t totalNodes
    ) {
    assert(nodes.size() == probs.cols());
    assert(nodes.size() == props.size());
    size_t curIteration = 1;
    while (true) {
      Eigen::VectorXd theta1 = getMax(probs, props, numReadDuplicates, numReads);
      normalize(theta1);
      Eigen::VectorXd theta2 = getMax(probs, theta1, numReadDuplicates, numReads);
      normalize(theta2);

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
        theta_p = props - 2 * alpha * r + alpha * alpha * v;
        props = getMax(probs, theta_p, numReadDuplicates, numReads);
        normalize(props);
        newllh = getExp(probs, props, numReadDuplicates);
      } else {
        theta_p = props - 2 * alpha * r + alpha * alpha * v;
        auto newProps = getMax(probs, theta_p, numReadDuplicates, numReads);
        normalize(newProps);
        newllh = getExp(probs, newProps, numReadDuplicates);
        if (newllh >= llh) {
          props = std::move(newProps);
        } else {
          while (llh - newllh > 0.0001) {
            // std::cerr << "it " << curit << "\t" << newllh << "\t" << llh << std::endl;
            alpha = (alpha - 1) / 2;
            theta_p = props - 2 * alpha * r + alpha * alpha * v;
            newProps = getMax(probs, theta_p, numReadDuplicates, numReads);
            normalize(newProps);
            newllh   = getExp(probs, newProps, numReadDuplicates);
          }
          props = std::move(newProps);
        }
      }

      updateInsigCounts(props, insigCounts, totalNodes);
      // std::cerr << "it " << curit << "\t" << newllh << "\t" << llh << std::endl;
      if (newllh - llh < 0.0001) {
        llh = newllh;
        converged = true;
        break;
      } else if (curIteration == iterations) {
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
    Tree *T, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, const std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex,
    const std::vector<bool>& lowScoreReads, const int32_t& numReads, const size_t& numLowScoreReads, const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    Eigen::MatrixXd& probs, std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const int32_t& roundsRemove, const double& removeThreshold, std::string exclude
  ) {
    if (exclude.empty()) {
      std::stringstream msg;
      msg << "starting to set up EM" << "\n";
      std::cerr << msg.str();
    } else {
      std::stringstream msg;
      msg << "starting to set up EM excluding " << exclude << "\n";
      std::cerr << msg.str();
    }


    if (!exclude.empty()) {
      probs.resize(allScores.begin()->second.size() - numLowScoreReads, allScores.size() - leastRecentIdenticalAncestors.size() - 1);
    } else {
      probs.resize(allScores.begin()->second.size() - numLowScoreReads, allScores.size() - leastRecentIdenticalAncestors.size());
    }
    size_t colIndex = 0;
    for (const auto& node : allScores) {
      if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
      if (!exclude.empty() && node.first == exclude) continue;
      std::vector<double> curProbs;
      size_t rowIndex = 0;
      for (size_t i = 0; i < node.second.size(); ++i) {
        if (!lowScoreReads[i]) {
          const auto& score = node.second[i];
          probs(rowIndex, colIndex) = score.second;
          ++rowIndex;
        }
      }
      nodes.push_back(node.first);
      ++colIndex;
    }

    // std::cerr << "num nodes " << nodes.size() << std::endl;
    props = Eigen::VectorXd::Constant(nodes.size(), 1.0 / static_cast<double>(nodes.size()));
    size_t totalNodes = nodes.size();
    Eigen::VectorXd readDuplicates(allScores.begin()->second.size() - numLowScoreReads);
    
    size_t indexReadDuplicates = 0;
    size_t numHighScoreReads = 0;
    for (size_t i = 0; i < readSeedmersDuplicatesIndex.size(); ++i) {
      if (!lowScoreReads[i]) {
        readDuplicates(indexReadDuplicates) = readSeedmersDuplicatesIndex[i].size();
        numHighScoreReads += readSeedmersDuplicatesIndex[i].size();
        ++indexReadDuplicates;
      }
    }

    if (exclude.empty()) {
      std::stringstream msg;
      msg << "starting EM estimation of haplotype proportions" << "\n";
      std::cerr << msg.str();
    } else {
      std::stringstream msg;
      msg << "starting EM estimation of haplotype proportions excluding " << exclude << "\n";
      std::cerr << msg.str();
    }


    int curit = 0;
    bool converged = false;
    int32_t filter_round = 0;
    size_t check_freq = 20;
    size_t remove_count = 20;
    std::vector<size_t> insigCounts(nodes.size());
    while (true) {
      std::cerr << "filter round " << filter_round + 1 << std::endl;
      ++filter_round;
      llh = getExp(probs, props, readDuplicates);
      squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, check_freq, insigCounts, totalNodes);
      if (converged) {
        break;
      }
      
      std::vector<size_t> significantIndices;
      std::vector<std::string> sigNodes;
      for (size_t i = 0; i < nodes.size(); ++i) {
        if (insigCounts[i] < remove_count) {
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
      std::cerr << "dropped " << nodes.size() - sigNodes.size() << " during EM" << std::endl;
      std::cerr << sigNodes.size() << " nodes left" << std::endl;
      nodes = std::move(sigNodes);
      probs = std::move(sigProbs);
      props = std::move(sigProps);
      normalize(props);
      insigCounts.assign(nodes.size(), 0);
      if (nodes.size() <= std::max(static_cast<int>(totalNodes) / 20, 100) || filter_round >= 10) {
        break;
      }
    }

    if (!converged) {
      std::vector<size_t> insigCounts(nodes.size());
      std::cerr << "start full EM" << std::endl;
      llh = getExp(probs, props, readDuplicates);
      squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, std::numeric_limits<size_t>::max(), insigCounts, totalNodes);
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
      std::cerr << "remove round " << i + 1 << std::endl;

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
      std::vector<size_t> insigCounts(sigNodes.size());
      squarem_test_1(sigNodes, sigProbs, identicalSets, readDuplicates, numHighScoreReads, sigProps, llh, curit, converged, iterations, insigCounts, totalNodes);
      assert(converged);
      nodes = std::move(sigNodes);
      probs = std::move(sigProbs);
      props = std::move(sigProps);
    }

    if (!exclude.empty()) {
      Eigen::VectorXd curProbs(allScores.begin()->second.size() - numLowScoreReads);
      size_t indexCurProbs = 0;
      const auto& curNode = allScores.at(exclude);
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

    if (exclude.empty()) {
      std::stringstream msg;
      msg << "Finished EM estimation of haplotype proportions. Total EM iterations: " << curit << "\n";
      std::cerr << msg.str();
    } else {
      std::stringstream msg;
      msg << "Finished EM estimation of haplotype proportions excluding " << exclude << ". Total EM iterations: " << curit<< "\n";
      std::cerr << msg.str();
    }
  }
}

#endif

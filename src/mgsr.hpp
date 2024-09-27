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

using namespace boost::icl;
typedef std::pair<std::vector<std::tuple<size_t*, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> readSeedmers_t;
typedef std::tuple<int32_t, int32_t, int32_t, int32_t, bool, int32_t> match_t;

static inline int64_t degapGlobal(const int64_t& globalCoord, const std::map<int64_t, int64_t>& coordsIndex) {
    auto coordIt = coordsIndex.upper_bound(globalCoord);
    if (coordIt == coordsIndex.begin()) {
        return 0;
    }
    return globalCoord - std::prev(coordIt)->second;
}

namespace mgsr {

  struct positionInfo {
    int64_t endPos;
    size_t fhash;
    size_t rhash;
    bool rev;
  };

  struct seedmers {
      //       beg                 end      fhash    rhash    rev
      std::map<int32_t, positionInfo> positionMap;
      //                 hash                       begs
      std::unordered_map<size_t, std::set<int32_t>> hashToPositionsMap;
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
                          std::map<int32_t, mgsr::positionInfo>& positionMap,
                          std::unordered_map<size_t, std::set<int32_t>>& hashToPositionsMap,
                          std::unordered_set<size_t>& affectedSeedmers,
                          const int& seedK,
                          const int& seedL,
                          const std::map<int64_t, int64_t>& coordIndex,
                          std::vector<std::tuple<int32_t, int32_t, size_t, size_t, bool>>& backTrackPositionMapChAdd,
                          std::vector<int32_t>& backTrackPositionMapErase) {
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
          if (isReplacement) {
            const auto& [oldEnd, oldFHash, oldRHash, oldRev] = curKminmerPositionIt->second;
            backTrackPositionMapChAdd.emplace_back(std::make_tuple(curKminmerPositionIt->first, oldEnd, oldFHash, oldRHash, oldRev));

            if (oldFHash != oldRHash) {
              size_t minHash = std::min(oldFHash, oldRHash);
              auto hashToPositionIt = hashToPositionsMap.find(minHash);
              hashToPositionIt->second.erase(curKminmerPositionIt->first);
              size_t ohashNum = hashToPositionIt->second.size();
              if (ohashNum <= 1) {
                affectedSeedmers.insert(minHash);
                if (ohashNum == 0) {
                  hashToPositionsMap.erase(hashToPositionIt);
                }
              }
            }
          } else {
            backTrackPositionMapErase.emplace_back(firstKminmerSeedIt->first);
          }

          size_t curforwardHash = rol(prevForwardKminmerHash, seedK) ^ rol(std::prev(firstKminmerSeedIt)->second.hash, seedK * seedL) ^ lastKminmerSeedIt->second.hash;
          size_t curReverseHash = ror(prevReverseKminmerHash, seedK) ^ ror(std::prev(firstKminmerSeedIt)->second.hash, seedK)         ^ rol(lastKminmerSeedIt->second.hash, seedK * (seedL-1));        
          
          if (isReplacement) {
            curKminmerPositionIt->second.endPos = lastKminmerSeedIt->second.endPos;
            curKminmerPositionIt->second.fhash = curforwardHash;
            curKminmerPositionIt->second.rhash = curReverseHash;
            curKminmerPositionIt->second.rev = curReverseHash < curforwardHash;
          } else {
            positionMap[firstKminmerSeedIt->first] = mgsr::positionInfo(lastKminmerSeedIt->second.endPos, curforwardHash, curReverseHash, curReverseHash < curforwardHash);
          }

          if (curforwardHash != curReverseHash) {
            hashToPositionsMap[std::min(curforwardHash, curReverseHash)].insert(firstKminmerSeedIt->first);
            affectedSeedmers.insert(std::min(curforwardHash, curReverseHash));
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
            
            if (oldFHash != oldRHash) {
              size_t minHash = std::min(oldFHash, oldRHash);
              auto hashToPositionIt = hashToPositionsMap.find(minHash);
              hashToPositionIt->second.erase(curKminmerPositionIt->first);
              size_t ohashNum = hashToPositionIt->second.size();
              if (ohashNum <= 1) {
                affectedSeedmers.insert(minHash);
                if (ohashNum == 0) {
                  hashToPositionsMap.erase(minHash);
                }
              }
            }

            curKminmerPositionIt->second.endPos = lastKminmerSeedIt->second.endPos;
            curKminmerPositionIt->second.fhash = forwardKminmerHash;
            curKminmerPositionIt->second.rhash = reverseKminmerHash;
            curKminmerPositionIt->second.rev = reverseKminmerHash < forwardKminmerHash;
          } else {
            backTrackPositionMapErase.emplace_back(firstKminmerSeedIt->first);
            positionMap[firstKminmerSeedIt->first] = mgsr::positionInfo(lastKminmerSeedIt->second.endPos, forwardKminmerHash, reverseKminmerHash, reverseKminmerHash < forwardKminmerHash);
          }

          if (forwardKminmerHash != reverseKminmerHash) {
            hashToPositionsMap[std::min(forwardKminmerHash, reverseKminmerHash)].insert(firstKminmerSeedIt->first);
            affectedSeedmers.insert(std::min(forwardKminmerHash, reverseKminmerHash));
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
          if (toEraseFHash != toEraseRHash) {
            size_t minHash = std::min(toEraseFHash, toEraseRHash);
            auto hashToPositionIt = hashToPositionsMap.find(minHash);
            hashToPositionIt->second.erase(pos);
            size_t ohashNum = hashToPositionIt->second.size();
            if (ohashNum <= 1) {
              affectedSeedmers.insert(minHash);
              if (ohashNum == 0) {
                hashToPositionsMap.erase(minHash);
              }
            }
          }
          positionMap.erase(toEraseIt);
        }
      }
    }

    auto lastPositionMapIt = std::prev(positionMap.end());
    while (lastPositionMapIt->first > maxBegCoord) {
      const auto& [toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev] = lastPositionMapIt->second;
      backTrackPositionMapChAdd.emplace_back(std::make_tuple(lastPositionMapIt->first, toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev));
      if (toEraseFHash != toEraseRHash) {
        size_t minHash = std::min(toEraseFHash, toEraseRHash);
        auto hashToPositionIt = hashToPositionsMap.find(minHash);
        hashToPositionIt->second.erase(lastPositionMapIt->first);
        size_t ohashNum = hashToPositionIt->second.size();
        if (ohashNum <= 1) {
          affectedSeedmers.insert(minHash);
          if (ohashNum == 0) {
            hashToPositionsMap.erase(minHash);
          }
        }
      }
      positionMap.erase(lastPositionMapIt);
      --lastPositionMapIt;
    }
  }

  int32_t extend(int64_t& curEnd, mgsr::Read& curRead, int rev, std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<int32_t>>& hashToPositionsMap, int32_t qidx, std::map<int32_t, mgsr::positionInfo>::const_iterator refPositionIt, int32_t c) {
    if (qidx == curRead.seedmersList.size() - 1) return c;
    const auto& [nhash, nqbeg, nqend, nqrev, nqidx] = curRead.seedmersList[qidx+1];

    auto nextHashToPositionIt = hashToPositionsMap.find(nhash);
    if (nextHashToPositionIt != hashToPositionsMap.end()) {
      if (nextHashToPositionIt->second.size() < 2) {
        if (nextHashToPositionIt->second.size() != 1) {
          std::cerr << "Error: hash " << nhash << " has multiple positions" << std::endl;
          exit(1);
        }

        const auto& rbeg = *(hashToPositionsMap.find(nhash)->second.begin());
        auto curRefPositionIt = positionMap.find(rbeg);
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

  void initializeMatches(mgsr::Read& curRead, std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<int32_t>>& hashToPositionsMap) {
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
          const auto& rbeg = *(hashToPositionIt->second.begin());
          auto curRefPositionIt = positionMap.find(rbeg);
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

  bool isColinear(const std::pair<boost::icl::discrete_interval<int32_t>, int>& match1, const std::pair<boost::icl::discrete_interval<int32_t>, int>& match2, const mgsr::Read& curRead, std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<int32_t>>& hashToPositionsMap, const std::map<int64_t, int64_t>& coordsIndex, const int& maximumGap) {
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
      const auto& rglobalbeg1 = *(hashToPositionsMap.find(first1.hash)->second.begin());
      const auto& rglobalend1 = (positionMap.find(*(hashToPositionsMap.find(last1.hash)->second.begin()))->second).endPos;
      const auto& rglobalbeg2 = *(hashToPositionsMap.find(first2.hash)->second.begin());
      // const auto& rglobalend2 = (positionMap.find(*(hashToPositionsMap.find(*last2.hash)->second.begin()))->second).endPos;
      
      const auto& qbeg1 = first1.begPos;
      const auto& qend1 = last1.endPos;
      const auto& qbeg2 = first2.begPos;
      const auto& qend2 = last2.endPos;
      auto rbeg1 = degapGlobal(rglobalbeg1, coordsIndex);
      auto rend1 = degapGlobal(rglobalend1, coordsIndex);
      auto rbeg2 = degapGlobal(rglobalbeg2, coordsIndex);
      // auto rend2 = degapGlobal(rglobalend2, coordsIndex);

      int32_t qgap = abs(qbeg2 - qend1);
      int32_t rgap = abs(rbeg2 - rend1);
      if (rbeg1 < rbeg2 && abs(qgap - rgap) < maximumGap) return true;

    } else {
      // reverse direction
      auto rglobalbeg1 = *(hashToPositionsMap.find(last1.hash)->second.begin());
      // auto rglobalend1 = (positionMap.find(*(hashToPositionsMap.find(*first1.hash)->second.begin()))->second).endPos;
      auto rglobalbeg2 = *(hashToPositionsMap.find(last2.hash)->second.begin());
      auto rglobalend2 = (positionMap.find(*(hashToPositionsMap.find(first2.hash)->second.begin()))->second).endPos;

      const auto& qbeg1 = first1.begPos;
      const auto& qend1 = last1.endPos;
      const auto& qbeg2 = first2.begPos;
      const auto& qend2 = last2.endPos;

      auto rbeg1 = degapGlobal(rglobalbeg1, coordsIndex);
      // auto rend1 = degapGlobal(rglobalend1, coordsIndex);
      auto rbeg2 = degapGlobal(rglobalbeg2, coordsIndex);
      auto rend2 = degapGlobal(rglobalend2, coordsIndex);

      int32_t qgap = abs(qbeg2 - qend1);
      int32_t rgap = abs(rbeg1 - rend2);
      if (rbeg2 < rbeg1 && abs(qgap - rgap) < maximumGap) return true;
    }

    return false;
  }

  int64_t getPseudoScore(
    const mgsr::Read& curRead, std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<int32_t>>& hashToPositionsMap,
    const std::map<int64_t, int64_t>& coordsIndex, const int& maximumGap, const int& minimumCount, const int& minimumScore
  ) {
    // simple cases
    if (curRead.matches.empty()) {
      return 0;
    } else if (curRead.matches.size() == 1) {
      return boost::icl::length(curRead.matches.begin()->first);
    }

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

    int64_t pseudoScore = 0;
    // find intervals colinear with longest interval and add length to pseudoScore
    const boost::icl::discrete_interval<int32_t>* firstMatch = nullptr;
    const boost::icl::discrete_interval<int32_t>* lastMatch = nullptr;
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
        if (isColinear(longestInterval, curInterval, curRead, positionMap, hashToPositionsMap, coordsIndex, maximumGap)) {
          pseudoScore += boost::icl::length(curInterval.first);
          if (firstMatch == nullptr) firstMatch = &curInterval.first;
          lastMatch = &curInterval.first;
        }
      } else if (longestQbeg > curQbeg) {
        // longest query beg after current query beg
        if (isColinear(curInterval, longestInterval, curRead, positionMap, hashToPositionsMap, coordsIndex, maximumGap)) {
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
    // const auto& first1 = curRead.seedmersList[boost::icl::first(*firstMatch)];
    // const auto& last1 = curRead.seedmersList[boost::icl::last(*lastMatch)];
    // const auto& first2 = curRead.seedmersList[boost::icl::first(*firstMatch)];
    // const auto& last2 = curRead.seedmersList[boost::icl::last(*lastMatch)];

    // auto rglobalbeg1 = *(hashToPositionsMap.find(*first1.hash)->second.begin());
    // auto rglobalend1 = (positionMap.find(*(hashToPositionsMap.find(*last1.hash)->second.begin()))->second).endPos;
    // auto rglobalbeg2 = *(hashToPositionsMap).find(*first2.hash)->second.begin();
    // auto rglobalend2 = (positionMap.find(*(hashToPositionsMap.find(*last2.hash)->second.begin()))->second).endPos;

    // get range of psuedoChain + readLen padding
    // for duplicate in duplicates:
      // for begs of duplicated hash:
        // if colinear pseudochain: add to colinear duplicates

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

}

#endif

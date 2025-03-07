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


namespace mgsr {
  enum readType : uint8_t {
    PASS,
    HIGH_DUPLICATES,
    IDENTICAL_SCORE_ACROSS_NODES
  };
  
  enum SeedmerStatus : uint8_t {
    EXIST_UNIQUE,
    EXIST_DUPLICATE,
    NOT_EXIST
  };

  enum SeedmerChangeType : uint8_t {
    EXIST_UNIQUE_TO_EXIST_UNIQUE = 0,
    EXIST_UNIQUE_TO_EXIST_DUPLICATE = 1,
    EXIST_UNIQUE_TO_NOT_EXIST = 2,
    EXIST_DUPLICATE_TO_EXIST_UNIQUE = 3,
    EXIST_DUPLICATE_TO_EXIST_DUPLICATE = 4,
    EXIST_DUPLICATE_TO_NOT_EXIST = 5,
    NOT_EXIST_TO_EXIST_UNIQUE = 6,
    NOT_EXIST_TO_EXIST_DUPLICATE = 7,
    NOT_EXIST_TO_NOT_EXIST = 8
  };

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

  struct readSeedmers {
    std::unordered_map<size_t, std::vector<std::pair<uint32_t, std::vector<uint32_t>>>> seedmerToReads;

    std::unordered_map<size_t, SeedmerStatus> seedmerStatus;

    void updateSeedmerStatus(const size_t& hash, const mgsr::SeedmerStatus& status) {
      seedmerStatus[hash] = status;
    }

    void updateSeedmerStatus(const size_t& hash, const std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap) {
      if (hashToPositionsMap.find(hash) == hashToPositionsMap.end()) {
        seedmerStatus[hash] = SeedmerStatus::NOT_EXIST;
      } else if (hashToPositionsMap.at(hash).size() == 1) {
        seedmerStatus[hash] = SeedmerStatus::EXIST_UNIQUE;
      } else if (hashToPositionsMap.at(hash).size() > 1) {
        seedmerStatus[hash] = SeedmerStatus::EXIST_DUPLICATE;
      } else {
        std::cerr << "Error: Invalid seedmer status." << std::endl;
        exit(1);
      }
    }

  };

  struct refSeedmers {
    //       beg                 end      fhash    rhash    rev
    std::map<int32_t, positionInfo> positionMap;
    //                 hash                       begs
    std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>> hashToPositionsMap;

    std::unordered_map<size_t, SeedmerStatus> seedmerStatus;

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
      const auto& beg = it->first;
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
      const auto& beg = it->first;
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

    template <typename T>
    std::pair<SeedmerStatus, bool> getCurrentSeedmerStatus(const T& hashOrIterator) {
      decltype(hashToPositionsMap)::iterator hashToPositionIt;
      
      if constexpr (std::is_same_v<T, size_t>) {
        hashToPositionIt = hashToPositionsMap.find(hashOrIterator);
      } else {
        hashToPositionIt = hashOrIterator;
      }
      
      if (hashToPositionIt == hashToPositionsMap.end()) {
        return std::make_pair(SeedmerStatus::NOT_EXIST, false);
      } else if (hashToPositionIt->second.size() == 1) {
        return std::make_pair(SeedmerStatus::EXIST_UNIQUE, (*(hashToPositionIt->second.begin()))->second.rev);
      } else if (hashToPositionIt->second.size() > 1) {
        return std::make_pair(SeedmerStatus::EXIST_DUPLICATE, false);
      } else {
        std::cerr << "Error: Invalid seedmer status." << std::endl;
        exit(1);
      }
    }

  };

  struct readSeedmer {
    const size_t hash;
    const int64_t begPos;
    const int64_t endPos;
    const bool rev;
    const int32_t iorder;
  };

  struct SeedmerState {
    bool match;
    bool rev;
    bool inChain;
  };
  
  typedef uint64_t minichain_t;


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
  bool compareMinichainByBeg(const minichain_t& a, const minichain_t& b) {
    return ((a >> 1) & 0x7FFFFFFF) < ((b >> 1) & 0x7FFFFFFF);
  }
  bool compareMinichainByEnd(const minichain_t& a, const minichain_t& b) {
    return ((a >> 32) & 0x7FFFFFFF) < ((b >> 32) & 0x7FFFFFFF);
  }
  class Read {
    public:
    std::vector<readSeedmer> seedmersList;
    std::vector<SeedmerState> seedmerStates;
    std::unordered_map<size_t, std::vector<uint32_t>> uniqueSeedmers;
    //if (rev) minichain_t |= 1ull;
    // minichain_t |= static_cast<uint64_t>(beg) << 1;
    // minichain_t |= static_cast<uint64_t>(end) << 32;

    // if minichain_t & 1 == 1 -> reversed
    // else -> forward
    // beg = (minichain_t >> 1) & 0x7FFFFFFF
    // end = (minichain_t >> 32) & 0x7FFFFFFF

    std::vector<minichain_t> minichains;
    boost::icl::split_interval_map<int32_t, int> matches;
    std::unordered_set<int32_t> duplicates;
    size_t readIndex;
    int64_t pseudoChainScores;
    double pseudoChainProb;

    // std::vector<bool> duplicates;
    // std::vector<bool> absentees;
    // size_t numDuplicates = 0;
    // size_t numAbsentees = 0;


    uint32_t extendMinichain(
      uint64_t& curEnd, bool rev, std::map<int32_t, mgsr::positionInfo>& positionMap,
      std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap,
      int32_t qidx, std::map<int32_t, mgsr::positionInfo>::const_iterator refPositionIt, uint64_t c
    ) {
      if (qidx == seedmersList.size() - 1) return c;
      const auto& [nhash, nqbeg, nqend, nqrev, nqidx] = seedmersList[qidx+1];

      auto nextHashToPositionIt = hashToPositionsMap.find(nhash);
      if (nextHashToPositionIt != hashToPositionsMap.end()) {
        if (nextHashToPositionIt->second.size() < 2) {
          if (nextHashToPositionIt->second.size() != 1) {
            std::cerr << "Error: hash " << nhash << " has multiple positions" << std::endl;
            exit(1);
          }

          auto curRefPositionIt = *nextHashToPositionIt->second.begin();
          bool nrev = nqrev != curRefPositionIt->second.rev;
          if (rev == nrev) {
            if (rev) {
              auto prevRefPositionIt = std::prev(refPositionIt);
              while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
              if (curRefPositionIt->first == prevRefPositionIt->first) {
                ++c;
                curEnd = nqidx;
                return extendMinichain(curEnd, rev, positionMap, hashToPositionsMap, nqidx, curRefPositionIt, c);
              }
            } else {
              auto nextRefPositionIt = std::next(refPositionIt);
              while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
              if (curRefPositionIt->first == nextRefPositionIt->first) {
                ++c;
                curEnd = nqidx;
                return extendMinichain(curEnd, rev, positionMap, hashToPositionsMap, nqidx, curRefPositionIt, c);
              }
            }
          }
        }
      }
      return c;
    }

    void initializeMinichains(std::map<int32_t, mgsr::positionInfo>& positionMap, std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap) {
      uint64_t i = 0;
      while (i < seedmersList.size()) {
        const auto& [hash, qbeg, qend, qrev, qidx] = seedmersList[i];
        uint64_t c = 1;
        auto hashToPositionIt = hashToPositionsMap.find(hash);
        if (hashToPositionIt != hashToPositionsMap.end()) {
          if (hashToPositionIt->second.size() < 2) {
            if (hashToPositionIt->second.size() != 1) {
              std::cerr << "Error: hash " << hash << " has multiple positions" << std::endl;
              exit(1);
            }
            auto curRefPositionIt = *(hashToPositionIt->second.begin());
            uint64_t curEnd = i;
            bool rev = qrev != curRefPositionIt->second.rev;
            c = extendMinichain(curEnd, rev, positionMap, hashToPositionsMap, qidx, curRefPositionIt, c);
            minichain_t minichain = (curEnd << 32) | (i << 1) | (rev ? 1ULL : 0ULL);
            minichains.push_back(minichain);
          } else {
            duplicates.insert(i);
          }
        }
        i += c; 
      }
    }


    void extendChainRemoval(uint64_t& c, uint64_t& curEnd, const std::vector<uint64_t>& affectedSeedmerIndexCodes) {
      if (curEnd == seedmersList.size() - 1 || c == affectedSeedmerIndexCodes.size()) return;

      uint32_t nextIndexOnSeedmerList = curEnd + 1;
      uint32_t nextIndexOnAffectedSeedmerIndexCodes = c;

      const uint64_t affectedSeedmerIndexCode = affectedSeedmerIndexCodes[nextIndexOnAffectedSeedmerIndexCodes];
      uint32_t affectedSeedmerIndex = (affectedSeedmerIndexCode >> 9) & 0xFFFFFFFF;
      mgsr::SeedmerChangeType SeedmerChangeType = static_cast<mgsr::SeedmerChangeType>((affectedSeedmerIndexCode >> 1) & 0xFF);

      if (nextIndexOnSeedmerList != affectedSeedmerIndex) return;

      if (SeedmerChangeType == mgsr::SeedmerChangeType::EXIST_UNIQUE_TO_EXIST_DUPLICATE || 
        SeedmerChangeType == mgsr::SeedmerChangeType::EXIST_UNIQUE_TO_NOT_EXIST
      ) {
        ++curEnd;
        ++c;
        extendChainRemoval(c, curEnd, affectedSeedmerIndexCodes);
      }
      return;
    }

    void extendChainAddition(uint64_t& c, uint64_t& curEnd, const std::vector<uint64_t>& affectedSeedmerIndexCodes, bool chainRev,
      std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap,
      std::map<int32_t, mgsr::positionInfo>& positionMap, std::map<int32_t, mgsr::positionInfo>::const_iterator refPositionIt
    ) {
      if (curEnd == seedmersList.size() - 1 || c == affectedSeedmerIndexCodes.size()) return;
      
      uint32_t nextIndexOnSeedmerList = curEnd + 1;
      uint32_t nextIndexOnAffectedSeedmerIndexCodes = c;

      const uint64_t affectedSeedmerIndexCode = affectedSeedmerIndexCodes[nextIndexOnAffectedSeedmerIndexCodes];
      uint32_t affectedSeedmerIndex = (affectedSeedmerIndexCode >> 9) & 0xFFFFFFFF;
      mgsr::SeedmerChangeType SeedmerChangeType = static_cast<mgsr::SeedmerChangeType>((affectedSeedmerIndexCode >> 1) & 0xFF);
      bool refRev = affectedSeedmerIndexCode & 1;
      if (nextIndexOnSeedmerList != affectedSeedmerIndex) return;

      if (SeedmerChangeType == mgsr::SeedmerChangeType::NOT_EXIST_TO_EXIST_UNIQUE ||
        SeedmerChangeType == mgsr::SeedmerChangeType::EXIST_DUPLICATE_TO_EXIST_UNIQUE
      ) {
        bool nextRev = refRev != seedmersList[affectedSeedmerIndex].rev;
        if (nextRev == chainRev) {
          auto curRefPositionIt = *(hashToPositionsMap.find(seedmersList[affectedSeedmerIndex].hash)->second.begin());
          if (chainRev) {
            auto prevRefPositionIt = std::prev(refPositionIt);
            while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
            if (curRefPositionIt->first == prevRefPositionIt->first) {
              ++c;
              ++curEnd;
              return extendChainAddition(c, curEnd, affectedSeedmerIndexCodes, chainRev, hashToPositionsMap, positionMap, curRefPositionIt);
            }
          } else {
            auto nextRefPositionIt = std::next(refPositionIt);
            while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
            if (curRefPositionIt->first == nextRefPositionIt->first) {
              ++c;
              ++curEnd;
              return extendChainAddition(c, curEnd, affectedSeedmerIndexCodes, chainRev, hashToPositionsMap, positionMap, curRefPositionIt);
            }
          }
        }
      }

      return;
    }


    void updateMinichains(size_t readIndex, const std::vector<uint64_t>& affectedSeedmerIndexCodes, std::unordered_map<size_t, std::set<std::map<int32_t, positionInfo>::iterator, IteratorComparator>>& hashToPositionsMap, std::map<int32_t, mgsr::positionInfo>& positionMap) {
      // decode affectedSeedmerIndexCode:
      // 0th bit is 1 -> reversed
      // 1st to 8th bit is the seedmerChangeType
      // 9th to 40th bit is the affectedSeedmerIndex
      std::vector<std::pair<minichain_t, bool>> updateMinichains;
      uint32_t i = 0;
      std::cout << "affectedSeedmerIndexCodes.size(): " << affectedSeedmerIndexCodes.size() << std::endl;
      while (i < affectedSeedmerIndexCodes.size()) {
        const uint64_t affectedSeedmerIndexCode = affectedSeedmerIndexCodes[i];
        uint32_t affectedSeedmerIndex = (affectedSeedmerIndexCode >> 9) & 0xFFFFFFFF;
        mgsr::SeedmerChangeType SeedmerChangeType = static_cast<mgsr::SeedmerChangeType>((affectedSeedmerIndexCode >> 1) & 0xFF);
        bool refRev = affectedSeedmerIndexCode & 1;
        uint64_t c = i + 1;
        uint64_t curEnd = affectedSeedmerIndex;

        if (SeedmerChangeType == mgsr::SeedmerChangeType::EXIST_UNIQUE_TO_EXIST_DUPLICATE || SeedmerChangeType == mgsr::SeedmerChangeType::EXIST_UNIQUE_TO_NOT_EXIST) {
          // match to no match -> remove from minichains
          extendChainRemoval(c, curEnd, affectedSeedmerIndexCodes);
          // encode minichain_t
          minichain_t minichain = (curEnd << 32) | (affectedSeedmerIndex << 1) | (0ULL);
          updateMinichains.push_back(std::make_pair(minichain, false));
        } else if (SeedmerChangeType == mgsr::SeedmerChangeType::EXIST_DUPLICATE_TO_EXIST_UNIQUE || SeedmerChangeType == mgsr::SeedmerChangeType::NOT_EXIST_TO_EXIST_UNIQUE) {
          // no match to match -> create minichains
          bool rev = refRev != seedmersList[affectedSeedmerIndex].rev;
          auto positionItFromCurrentHash = *(hashToPositionsMap.find(seedmersList[affectedSeedmerIndex].hash)->second.begin());
          extendChainAddition(c, curEnd, affectedSeedmerIndexCodes, rev, hashToPositionsMap, positionMap, positionItFromCurrentHash);
          // encode minichain_t
          minichain_t minichain = (curEnd << 32) | (affectedSeedmerIndex << 1) | (rev ? 1ULL : 0ULL);
          updateMinichains.push_back(std::make_pair(minichain, true));
        }
        i += (curEnd - affectedSeedmerIndex + 1);
      }


      // add/merge/delete minichains
      for (const auto& [minichain, toAdd] : updateMinichains) {
        uint64_t minichainBeg = (minichain >> 1) & 0x7FFFFFFF;
        uint64_t minichainEnd = (minichain >> 32) & 0x7FFFFFFF;
        bool minichainRev = minichain & 1;
        std::cout << "minichainBeg: " << minichainBeg << ", minichainEnd: " << minichainEnd << ", minichainRev: " << minichainRev << ", toAdd: " << toAdd << std::endl;
        if (toAdd) {
          bool addRangeRev = minichain & 1;
          uint64_t addRangeBeg = (minichain >> 1) & 0x7FFFFFFF;
          uint64_t addRangeEnd = (minichain >> 32) & 0x7FFFFFFF;
          if (minichains.size() == 0) {
            minichains.push_back(minichain);
          } else if (minichains.size() == 1) {
            auto& originalMinichain = minichains[0];
            uint64_t originalMinichainBeg = (originalMinichain >> 1) & 0x7FFFFFFF;
            uint64_t originalMinichainEnd = (originalMinichain >> 32) & 0x7FFFFFFF;
            uint64_t originalMinichainRev = originalMinichain & 1;

            if (originalMinichainBeg != 0 && addRangeEnd == originalMinichainBeg - 1) {
              std::cout << "Attempting to merge to the right" << std::endl;
              if (addRangeRev == originalMinichainRev) {
                auto newRangeEndKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[addRangeEnd].hash)->second.begin());
                auto originalRangeBegKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[originalMinichainBeg].hash)->second.begin());
                if (originalMinichainRev) {
                  auto prevRefPositionIt = std::prev(newRangeEndKminmerPositionIt);
                  while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
                  if (prevRefPositionIt->first == originalRangeBegKminmerPositionIt->first) {
                    // merge
                    originalMinichain = (originalMinichainEnd << 32) | (addRangeBeg << 1) | (originalMinichainRev ? 1ULL : 0ULL);
                  } else {
                    // add without merging 
                    minichains.insert(minichains.begin(), minichain);
                  }
                } else {
                  auto nextRefPositionIt = std::next(newRangeEndKminmerPositionIt);
                  while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
                  if (nextRefPositionIt->first == originalRangeBegKminmerPositionIt->first) {
                    // merge
                    originalMinichain = (originalMinichainEnd << 32) | (addRangeBeg << 1) | (originalMinichainRev ? 1ULL : 0ULL);
                  } else {
                    // add without merging
                   minichains.insert(minichains.begin(), minichain);
                  }
                }
              } else {
                // add without merging
                std::cout << "Merge failed due to opposite direction" << std::endl;
                minichains.insert(minichains.begin(), minichain);
              }
            } else if (addRangeBeg == originalMinichainEnd + 1) {
              std::cout << "Attempting to merge to the left" << std::endl;
              if (addRangeRev == originalMinichainRev) {
                auto newRangeBegKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[addRangeBeg].hash)->second.begin());
                auto originalRangeEndKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[originalMinichainEnd].hash)->second.begin());
                if (originalMinichainRev) {
                  auto prevRefPositionIt = std::prev(originalRangeEndKminmerPositionIt);
                  while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
                  if (prevRefPositionIt->first == originalRangeEndKminmerPositionIt->first) {
                    // merge
                    originalMinichain = (addRangeEnd << 32) | (originalMinichainBeg << 1) | (originalMinichainRev ? 1ULL : 0ULL);
                  } else {
                    // add without merging
                    minichains.push_back(minichain);
                  }
                } else {
                  auto nextRefPositionIt = std::next(originalRangeEndKminmerPositionIt);
                  while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
                  if (nextRefPositionIt->first == newRangeBegKminmerPositionIt->first) {
                    // merge
                    originalMinichain = (addRangeEnd << 32) | (originalMinichainBeg << 1) | (originalMinichainRev ? 1ULL : 0ULL);
                  } else {
                    // add without merging
                    minichains.push_back(minichain);
                  }
                }
              } else {
                // add without merging
                std::cout << "Merge failed due to opposite direction" << std::endl;
                minichains.push_back(minichain);
              }
            } else {
              // add without merging
              if (addRangeEnd < originalMinichainBeg) {
                minichains.insert(minichains.begin(), minichain);
              } else {
                minichains.push_back(minichain);
              }
            }
          } else {
            bool leftMinichainExists = false;
            bool rightMinichainExists = false;
            auto rightMinichainIt = std::upper_bound(minichains.begin(), minichains.end(), minichain, compareMinichainByBeg);
            auto leftMinichainIt = minichains.end();
            if (rightMinichainIt != minichains.end()) {
              rightMinichainExists = true;
            }
            if (rightMinichainIt != minichains.begin()) {
              leftMinichainIt = std::prev(rightMinichainIt);
              leftMinichainExists = true;
            }
            std::cout << "leftMinichainExists: " << leftMinichainExists << ", rightMinichainExists: " << rightMinichainExists << std::endl;

            bool mergeLeft = false;
            bool mergeRight = false;
            uint64_t leftMinichainBeg;
            uint64_t leftMinichainEnd;
            bool leftMinichainRev;
            uint64_t rightMinichainBeg;
            uint64_t rightMinichainEnd;
            bool rightMinichainRev;
            if (leftMinichainExists) {
              leftMinichainBeg = (*leftMinichainIt >> 1) & 0x7FFFFFFF;
              leftMinichainEnd = (*leftMinichainIt >> 32) & 0x7FFFFFFF;
              leftMinichainRev = *leftMinichainIt & 1;
              if (addRangeRev == leftMinichainRev) {
                if (addRangeBeg == leftMinichainEnd + 1) {
                  auto newRangeBegKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[addRangeBeg].hash)->second.begin());
                  auto leftRangeEndKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[leftMinichainEnd].hash)->second.begin());
                  if (leftMinichainRev) {
                    auto prevRefPositionIt = std::prev(leftRangeEndKminmerPositionIt);
                    while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
                    if (prevRefPositionIt->first == newRangeBegKminmerPositionIt->first) {
                      // merge
                      mergeLeft = true;
                    } else {
                      // add without merging 
                      mergeLeft = false;
                    }
                  } else {
                    auto nextRefPositionIt = std::next(newRangeBegKminmerPositionIt);
                    while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
                    if (nextRefPositionIt->first == leftRangeEndKminmerPositionIt->first) {
                      // merge
                      mergeLeft = true;
                    } else {
                      // add without merging
                      mergeLeft = false;
                    }
                  }
                } else {
                  // add without merging
                  mergeLeft = false;
                }
              } else {
                // add without merging
                mergeLeft = false;
              }
            }
            
            if (rightMinichainExists) {
              rightMinichainBeg = (*rightMinichainIt >> 1) & 0x7FFFFFFF;
              rightMinichainEnd = (*rightMinichainIt >> 32) & 0x7FFFFFFF;
              rightMinichainRev = *rightMinichainIt & 1;
              if (addRangeRev == rightMinichainRev) {
                auto newRangeEndKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[addRangeEnd].hash)->second.begin());
                auto rightRangeBegKminmerPositionIt = *(hashToPositionsMap.find(seedmersList[rightMinichainBeg].hash)->second.begin());
                if (rightMinichainRev) {
                  auto prevRefPositionIt = std::prev(newRangeEndKminmerPositionIt);
                  while (prevRefPositionIt->second.fhash == prevRefPositionIt->second.rhash) --prevRefPositionIt;
                  if (prevRefPositionIt->first == rightRangeBegKminmerPositionIt->first) {
                    // merge
                    mergeRight = true;
                  } else {
                    // add without merging 
                    mergeRight = false;
                  }
                } else {
                  auto nextRefPositionIt = std::next(newRangeEndKminmerPositionIt);
                  while (nextRefPositionIt->second.fhash == nextRefPositionIt->second.rhash) ++nextRefPositionIt;
                  if (nextRefPositionIt->first == rightRangeBegKminmerPositionIt->first) {
                    // merge
                    mergeRight = true;
                  } else {
                    // add without merging
                    mergeRight = false;
                  }
                }
              } else {
                // add without merging
                mergeRight = false;
              }
            }
            std::cout << "mergeLeft: " << mergeLeft << ", mergeRight: " << mergeRight << std::endl;
            if (mergeLeft && mergeRight) {
              *leftMinichainIt = (rightMinichainEnd << 32) | (leftMinichainBeg << 1) | (addRangeRev ? 1ULL : 0ULL);
              minichains.erase(rightMinichainIt);
            } else if (mergeLeft) {
              *leftMinichainIt = (addRangeEnd << 32) | (leftMinichainBeg << 1) | (addRangeRev ? 1ULL : 0ULL);
            } else if (mergeRight) {
              *rightMinichainIt = (rightMinichainEnd << 32) | (addRangeBeg << 1) | (addRangeRev ? 1ULL : 0ULL);
            } else {
              if (!leftMinichainExists) {
                minichains.insert(minichains.begin(), minichain);
              } else if (!rightMinichainExists) {
                minichains.push_back(minichain);
              } else {
                minichains.insert(leftMinichainIt+1, minichain);
              }
            }
          }
        } else {
          uint64_t removeRangeBeg = (minichain >> 1) & 0x7FFFFFFF;
          uint64_t removeRangeEnd = (minichain >> 32) & 0x7FFFFFFF;
          if (minichains.size() == 1) {
            auto& originalMinichain = minichains[0];
            uint64_t originalMinichainBeg = (originalMinichain >> 1) & 0x7FFFFFFF;
            uint64_t originalMinichainEnd = (originalMinichain >> 32) & 0x7FFFFFFF;
            uint64_t originalMinichainRev = originalMinichain & 1;
            if (originalMinichainBeg == removeRangeBeg) {
              if (originalMinichainEnd == removeRangeEnd) {
                // remove the minichain
                minichains.clear();
              } else {
                // edit the beg
                originalMinichain = (originalMinichainEnd << 32) | ((removeRangeEnd + 1) << 1) | (originalMinichainRev ? 1ULL : 0ULL);
              }
            } else if (originalMinichainEnd == removeRangeEnd) {
              if (originalMinichainBeg == removeRangeBeg) {
                // remove the minichain
                minichains.clear();
              } else {
                // edit the end
                originalMinichain = ((removeRangeBeg - 1) << 32) | (originalMinichainBeg << 1) | (originalMinichainRev ? 1ULL : 0ULL);
              }
            } else {
              // need to split a minichain
              originalMinichain = ((removeRangeBeg - 1) << 32) | (originalMinichainBeg << 1) | (originalMinichainRev ? 1ULL : 0ULL);
              uint64_t newMinichain = (originalMinichainEnd << 32) | ((removeRangeEnd + 1) << 1) | (originalMinichainRev ? 1ULL : 0ULL);
              minichains.push_back(newMinichain);
            }
          } else {
            auto currentOriginalMinichainIt = std::upper_bound(minichains.begin(), minichains.end(), minichain, compareMinichainByBeg);
            --currentOriginalMinichainIt;
            decltype(currentOriginalMinichainIt) endOfOverlappingOriginalMinichainIt;
            uint64_t currentOriginalMinichainBeg = (*currentOriginalMinichainIt >> 1) & 0x7FFFFFFF;
            uint64_t currentOriginalMinichainEnd = (*currentOriginalMinichainIt >> 32) & 0x7FFFFFFF;
            uint64_t currentOriginalMinichainRev = *currentOriginalMinichainIt & 1;
            std::cout << "currentOriginalMinichainBeg: " << currentOriginalMinichainBeg << ", currentOriginalMinichainEnd: " << currentOriginalMinichainEnd << ", currentOriginalMinichainRev: " << currentOriginalMinichainRev << ", removeRangeBeg: " << removeRangeBeg << ", removeRangeEnd: " << removeRangeEnd << std::endl;

            std::cout << "Line: " << __LINE__ << std::endl;
            if (removeRangeEnd > currentOriginalMinichainEnd) {
              std::cout << "Line: " << __LINE__ << std::endl;
              uint32_t numToErase = 0;
              std::cout << "Line: " << __LINE__ << std::endl;
              auto it = currentOriginalMinichainIt;
              std::cout << "Line: " << __LINE__ << std::endl;
              if (currentOriginalMinichainBeg == removeRangeBeg) {
                std::cout << "Line: " << __LINE__ << std::endl;
                // remove the minichain
                ++numToErase;
                ++currentOriginalMinichainIt;
              } else {
                std::cout << "Line: " << __LINE__ << std::endl;
                // edit the end then step right
                *currentOriginalMinichainIt = ((removeRangeBeg - 1) << 32) | (currentOriginalMinichainBeg << 1) | (currentOriginalMinichainRev ? 1ULL : 0ULL);
                std::cout << "Line: " << __LINE__ << std::endl;
                ++currentOriginalMinichainIt;
                ++it;
              }

              std::cout << "Line: " << __LINE__ << std::endl;
              currentOriginalMinichainEnd = (*currentOriginalMinichainIt >> 32) & 0x7FFFFFFF;
              std::cout << "Before while loop: " << currentOriginalMinichainEnd << std::endl;
              std::cout << "Line: " << __LINE__ << std::endl;
              while (currentOriginalMinichainIt != minichains.end() && currentOriginalMinichainEnd <= removeRangeEnd) {
                std::cout << "Line: " << __LINE__ << std::endl;
                ++numToErase;
                std::cout << "Line: " << __LINE__ << std::endl;
                ++currentOriginalMinichainIt;
                currentOriginalMinichainEnd = (*currentOriginalMinichainIt >> 32) & 0x7FFFFFFF;
              }

              std::cout << "Line: " << __LINE__ << std::endl;
              currentOriginalMinichainBeg = (*currentOriginalMinichainIt >> 1) & 0x7FFFFFFF;
              std::cout << "Line: " << __LINE__ << std::endl;
              if (currentOriginalMinichainIt != minichains.end() && currentOriginalMinichainBeg <= removeRangeEnd) {
                std::cout << "Line: " << __LINE__ << std::endl;
                *currentOriginalMinichainIt = (currentOriginalMinichainEnd << 32) | ((removeRangeEnd + 1) << 1) | (currentOriginalMinichainRev ? 1ULL : 0ULL);
              }

              std::cout << "Line: " << __LINE__ << std::endl;
              std::cout << "numToErase: " << numToErase << std::endl;
              for (size_t i = 0; i < numToErase; ++i) {
                auto curit = it+i;
                uint64_t curitMinichain = *curit;
                uint64_t curitMinichainBeg = (curitMinichain >> 1) & 0x7FFFFFFF;
                uint64_t curitMinichainEnd = (curitMinichain >> 32) & 0x7FFFFFFF;
                uint64_t curitMinichainRev = curitMinichain & 1;
                std::cout << "To erase: " << curitMinichainBeg << " " << curitMinichainEnd << " " << curitMinichainRev << std::endl;
              }
              minichains.erase(it, it+numToErase);

              // endOfOverlappingOriginalMinichainIt = std::upper_bound(minichains.begin(), minichains.end(), minichain, compareMinichainByEnd);
              // --endOfOverlappingOriginalMinichainIt;
              // uint64_t endOfOverlappingOriginalMinichainEnd = (*endOfOverlappingOriginalMinichainIt >> 32) & 0x7FFFFFFF;
              // std::cout << "endOfOverlappingOriginalMinichainEnd: " << endOfOverlappingOriginalMinichainEnd << std::endl;
              // if (endOfOverlappingOriginalMinichainEnd == removeRangeEnd) {
              //   minichains.erase(currentOriginalMinichainIt, endOfOverlappingOriginalMinichainIt+1);
              // } else {
              //   uint64_t endOfOverlappingOriginalMinichainBeg = (*endOfOverlappingOriginalMinichainIt >> 1) & 0x7FFFFFFF;
              //   uint64_t endOfOverlappingOriginalMinichainRev = *endOfOverlappingOriginalMinichainIt & 1;
              //   *endOfOverlappingOriginalMinichainIt = (endOfOverlappingOriginalMinichainEnd << 32) | ((removeRangeEnd + 1) << 1) | (endOfOverlappingOriginalMinichainRev ? 1ULL : 0ULL);
              //   minichains.erase(currentOriginalMinichainIt, endOfOverlappingOriginalMinichainIt);
              // }
            } else {
              // same as if there is one minichain but instead of clear, use erase
              if (currentOriginalMinichainBeg == removeRangeBeg) {
                if (currentOriginalMinichainEnd == removeRangeEnd) {
                  // erase the minichain
                  minichains.erase(currentOriginalMinichainIt);
                } else {
                  // edit the beg
                  *currentOriginalMinichainIt = (currentOriginalMinichainEnd << 32) | ((removeRangeEnd + 1) << 1) | (currentOriginalMinichainRev ? 1ULL : 0ULL);
                }
              } else if (currentOriginalMinichainEnd  == removeRangeEnd) {
                if (currentOriginalMinichainBeg == removeRangeBeg) {
                  // erase the minichain
                  minichains.erase(currentOriginalMinichainIt);
                } else {
                  // edit the end
                  *currentOriginalMinichainIt = ((removeRangeBeg - 1) << 32) | (currentOriginalMinichainBeg << 1) | (currentOriginalMinichainRev ? 1ULL : 0ULL);
                }
              } else {
                *currentOriginalMinichainIt = ((removeRangeBeg - 1) << 32) | (currentOriginalMinichainBeg << 1) | (currentOriginalMinichainRev ? 1ULL : 0ULL);
                uint64_t newMinichain = (currentOriginalMinichainEnd << 32) | ((removeRangeEnd + 1) << 1) | (currentOriginalMinichainRev ? 1ULL : 0ULL);
                minichains.insert(currentOriginalMinichainIt+1, newMinichain);
              }
            }
          }
        }
        for (const auto& minichain : minichains) {
          uint64_t minichainBeg = (minichain >> 1) & 0x7FFFFFFF;
          uint64_t minichainEnd = (minichain >> 32) & 0x7FFFFFFF;
          bool minichainIsReversed = minichain & 1;
          std::cout << "Current read " << readIndex << " Minichain: " << minichainBeg << " " << minichainEnd << " " << minichainIsReversed << std::endl;
        }
      }

    }

    size_t nextValidRefSeedmer(std::map<int32_t, mgsr::positionInfo>::iterator& currentPositionIt, mgsr::refSeedmers& seedmersIndex, bool rev) {
      if (rev) {
        while (currentPositionIt != seedmersIndex.positionMap.begin()) {
          --currentPositionIt;
          if (currentPositionIt->second.fhash != currentPositionIt->second.rhash) {
            return std::min(currentPositionIt->second.fhash, currentPositionIt->second.rhash);
          }
        }
      } else {
        ++currentPositionIt;
        while (currentPositionIt != seedmersIndex.positionMap.end()) {
          if (currentPositionIt->second.fhash != currentPositionIt->second.rhash) {
            return std::min(currentPositionIt->second.fhash, currentPositionIt->second.rhash);
          }
          ++currentPositionIt;
        }
      }

      return std::numeric_limits<size_t>::max();
    }

    bool isColinearFromMinichains(const bool& rev, const uint64_t beg1, const uint64_t end1, const uint64_t beg2, const uint64_t end2, mgsr::refSeedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex, const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap, const int& dfsIndex) {
      const auto& first1 = seedmersList[beg1];
      const auto& last1 = seedmersList[end1];
      const auto& first2 = seedmersList[beg2];
      const auto& last2 = seedmersList[end2];

      if (rev) {
        // forward direction
        const auto& rglobalbeg1 = seedmersIndex.getBegFromHash(first1.hash);
        const auto& rglobalend1 = seedmersIndex.getEndFromHash(last1.hash);
        const auto& rglobalbeg2 = seedmersIndex.getBegFromHash(first2.hash);
        
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

    int64_t getPsuedoChainScoreFromMinichains(const size_t& readIndex, mgsr::refSeedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex, const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap, const int& dfsIndex) {
      int64_t pseudoChainScore = 0;
      if (minichains.empty()) {
        return 0;
      } else if (minichains.size() == 1) {
        const uint64_t currentMinichain = minichains[0];
        uint64_t currentMinichainBeg = (currentMinichain >> 1) & 0x7FFFFFFF;
        uint64_t currentMinichainEnd = (currentMinichain >> 32) & 0x7FFFFFFF;
        return currentMinichainEnd - currentMinichainBeg + 1;
      } else {
        // find longest minichain
        uint64_t longestMinichainLength = 0;
        int32_t longestMinichainIndex = -1;
        uint64_t longestMinichainCode = 0;
        for (int32_t i = 0; i < minichains.size(); ++i) {
          const uint64_t currentMinichain = minichains[i];
          uint64_t currentMinichainBeg = (currentMinichain >> 1) & 0x7FFFFFFF;
          uint64_t currentMinichainEnd = (currentMinichain >> 32) & 0x7FFFFFFF;
          if (currentMinichainEnd - currentMinichainBeg + 1 > longestMinichainLength) {
            longestMinichainIndex = i;
            longestMinichainCode = currentMinichain;
            longestMinichainLength = currentMinichainEnd - currentMinichainBeg + 1;
          }
        }

        bool longestMinichainIsReversed = longestMinichainCode & 1;
        uint64_t longestMinichainBeg = (longestMinichainCode >> 1) & 0x7FFFFFFF;
        uint64_t longestMinichainEnd = (longestMinichainCode >> 32) & 0x7FFFFFFF;

        for (int32_t i = 0; i < minichains.size(); ++i) {
          const uint64_t currentMinichain = minichains[i];
          uint64_t currentMinichainBeg = (currentMinichain >> 1) & 0x7FFFFFFF;
          uint64_t currentMinichainEnd = (currentMinichain >> 32) & 0x7FFFFFFF;
          bool currentMinichainIsReversed = currentMinichain & 1;
          if (i == longestMinichainIndex) {
            pseudoChainScore += longestMinichainEnd - longestMinichainBeg + 1;
          }

          if (currentMinichainIsReversed != longestMinichainIsReversed) {
            continue;
          }

          if (longestMinichainIndex < i) {
            if (isColinearFromMinichains(longestMinichainIsReversed, longestMinichainBeg, longestMinichainEnd, currentMinichainBeg, currentMinichainEnd, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap, dfsIndex)) {
              pseudoChainScore += currentMinichainEnd - currentMinichainBeg + 1;
            }
          } else if (longestMinichainIndex > i) {
            if (isColinearFromMinichains(longestMinichainIsReversed, currentMinichainBeg, currentMinichainEnd, longestMinichainBeg, longestMinichainEnd, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap, dfsIndex)) {
              pseudoChainScore += currentMinichainEnd - currentMinichainBeg + 1;
            }
          }
        }
      }
      return pseudoChainScore;
    }
  };

  class ReadScores {
    public:
      Tree* T;
      std::vector<std::vector<std::tuple<size_t, int32_t, double>>> perNodeScoreDeltasIndex;
      std::unordered_map<std::string, int64_t> nodeToDfsIndex;
      std::vector<std::pair<int32_t, double>> scores;
      int32_t totalScore;

      ReadScores(Tree* T, size_t numReads, size_t numNodes) {
        this->T = T;
        scores.resize(numReads, std::make_pair(0, 0.0));
        perNodeScoreDeltasIndex.resize(numNodes);
        totalScore = 0;
      }

      void setScore(const size_t& readIndex, const int32_t& score, const double& prob, size_t numDuplicates) {
        totalScore += (score - scores[readIndex].first) * numDuplicates;
        scores[readIndex] = std::make_pair(score, prob);
      }

      void reserveMutationsIndex(size_t nodeIndex, size_t numChangedReads) {
        perNodeScoreDeltasIndex[nodeIndex].reserve(numChangedReads);
      }

      void assignDfsIndex(const std::string& nodeIdentifier, const int64_t& dfsIndex) {
        nodeToDfsIndex[nodeIdentifier] = dfsIndex;
      }

      void addScoreMutation(size_t nodeIndex, size_t readIndex, int32_t score, double prob) {
        perNodeScoreDeltasIndex[nodeIndex].emplace_back(std::make_tuple(readIndex, score, prob));
      }

      std::pair<int32_t, double> getScoreAtCurrentNode(const size_t& readIndex) {
        return scores[readIndex];
      }

      std::vector<std::pair<int32_t, double>> getScoresAtNode(const std::string& nodeIdentifier) {
        panmanUtils::Node* currentNode = T->allNodes[nodeIdentifier];
        std::vector<panmanUtils::Node*> nodePath;
        while (currentNode->parent != nullptr) {
          nodePath.push_back(currentNode);
          currentNode = currentNode->parent;
        }
        nodePath.push_back(currentNode);
        std::reverse(nodePath.begin(), nodePath.end());
        
        std::vector<std::pair<int32_t, double>> nodeScores(scores.size(), std::make_pair(0, 0.0));
        for (const auto& node : nodePath) {
          for (const auto& scoreDelta : perNodeScoreDeltasIndex.at(nodeToDfsIndex.at(node->identifier))) {
            nodeScores[std::get<0>(scoreDelta)].first = std::get<1>(scoreDelta);
            nodeScores[std::get<0>(scoreDelta)].second = std::get<2>(scoreDelta);
          }
        }
        return nodeScores;
      }

      bool identicalReadScores(const std::string& node1Identifier, const std::string& node2Identifier) {
        panmanUtils::Node* currentNode1 = T->allNodes[node1Identifier];
        panmanUtils::Node* currentNode2 = T->allNodes[node2Identifier];
        std::vector<panmanUtils::Node*> nodePath1;
        std::vector<panmanUtils::Node*> nodePath2;

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


        std::vector<std::pair<int32_t, double>> lcaScores = getScoresAtNode(nodePath1[lcaIndex]->identifier);
        std::vector<std::pair<int32_t, double>> currNode1Scores = lcaScores;
        std::vector<std::pair<int32_t, double>> currNode2Scores = lcaScores;

        std::vector<size_t> changedReadsIndices;
        for (size_t i = lcaIndex + 1; i < nodePath1.size(); ++i) {
          Node* currNode = nodePath1[i];
          for (const auto& scoreDelta : perNodeScoreDeltasIndex.at(nodeToDfsIndex.at(currNode->identifier))) {
            changedReadsIndices.push_back(std::get<0>(scoreDelta));
            currNode1Scores[std::get<0>(scoreDelta)].first = std::get<1>(scoreDelta);
            currNode1Scores[std::get<0>(scoreDelta)].second = std::get<2>(scoreDelta);
          }
        }

        for (size_t i = lcaIndex + 1; i < nodePath2.size(); ++i) {
          Node* currNode = nodePath2[i];
          for (const auto& scoreDelta : perNodeScoreDeltasIndex.at(nodeToDfsIndex.at(currNode->identifier))) {
            changedReadsIndices.push_back(std::get<0>(scoreDelta));
            currNode2Scores[std::get<0>(scoreDelta)].first = std::get<1>(scoreDelta);
            currNode2Scores[std::get<0>(scoreDelta)].second = std::get<2>(scoreDelta);
          }
        }

        for (const auto& readIndex : changedReadsIndices) {
          if (currNode1Scores[readIndex].first != currNode2Scores[readIndex].first) {
            return false;
          }
        }
        return true;
      }

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
    const int64_t& dfsIndex, 
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
          auto [newSeed, newSeedEndPos] = seed_annotated_tree::getSeedAt(seedPos, T, seedK, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap, blockRanges);
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
            seedChanges.emplace_back(std::make_tuple(seedPos, true, true, oldSeed, newSeedHash, oldIsReverse, newIsReverse, oldEndPos, newSeedEndPos));
          } else { // off -> on
            seedChanges.emplace_back(std::make_tuple(seedPos, false, true, std::nullopt, newSeedHash, std::nullopt, newIsReverse, std::nullopt, newSeedEndPos));
          }
        }
      }
    }
  }

  template <typename SeedMutationsType, typename GapMutationsType>
  void processNodeMutations(
    const GapMutationsType& perNodeGapMutations_Index,
    const SeedMutationsType& perNodeSeedMutations_Index,
    const int64_t& dfsIndex,
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
                          mgsr::refSeedmers& seedmersIndex,
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

  void initializeSeedmerStates(mgsr::Read& curRead, mgsr::refSeedmers& seedmersIndex) {
    auto& curSeedmerList = curRead.seedmersList;
    auto& curSeedmerStates = curRead.seedmerStates;
    for (size_t i = 0; i < curSeedmerList.size(); ++i) {
      auto hashToPositionIt = seedmersIndex.hashToPositionsMap.find(curSeedmerList[i].hash);
      auto [status, refRev] = seedmersIndex.getCurrentSeedmerStatus(hashToPositionIt);
      if (status == mgsr::SeedmerStatus::EXIST_UNIQUE) {
        curSeedmerStates[i].match = true;
        curSeedmerStates[i].rev = curSeedmerList[i].rev != refRev;
      } else if (status == mgsr::SeedmerStatus::EXIST_DUPLICATE) {
        curSeedmerStates[i].match = false;
        curRead.duplicates.insert(i);
      } else {
        curSeedmerStates[i].match = false;
      }
    }
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

  bool isColinear(const std::pair<boost::icl::discrete_interval<int32_t>, int>& match1, const std::pair<boost::icl::discrete_interval<int32_t>, int>& match2, const mgsr::Read& curRead, mgsr::refSeedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex, const std::map<int64_t, int64_t>& regapCoordIndex, const int& maximumGap, const int& dfsIndex) {
    bool rev1 = match1.second == 1 ? false : true;
    bool rev2 = match2.second == 1 ? false : true;
    if (rev1 != rev2) {
      std::cerr << "Error: Invalid direction in isColinear" << std::endl;
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
    const mgsr::Read& curRead, mgsr::refSeedmers& seedmersIndex, const std::map<int64_t, int64_t>& degapCoordIndex,
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
      int32_t leftBoundLocal = std::max(static_cast<int64_t>(0), degapGlobal(leftBoundGlobal, degapCoordIndex) - 120);
      int32_t rightBoundLocal = degapGlobal(rightBoundGlobal, degapCoordIndex) + 120;

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

    Eigen::VectorXd ratios(numNodes);
    Eigen::VectorXd inverse_denoms = denoms.array().inverse();
    for (size_t i = 0; i < numNodes; ++i) {
      const auto& col_i = probs.col(i);
      ratios.array() = col_i.array() * props[i] * inverse_denoms.array();
      double newProp = (numReadDuplicates.array() * ratios.array()).sum();
      newProp /= totalReads;
      newProps(i) = newProp;
    }


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

  void exclude_noninformative_reads(std::vector<mgsr::readType>& readTypes, const std::vector<std::pair<std::string, std::vector<std::pair<int32_t, double>>>>& ProbableNodeScores) {
    std::vector<size_t> uninformativeReads;

    for (size_t i = 0; i < ProbableNodeScores[0].second.size(); ++i) {
      int32_t curScore = -1;
      bool informative = false;
      for (const auto& [node, scores] : ProbableNodeScores) {
        if (curScore == -1) {
          curScore = scores[i].first;
        } else {
          if (scores[i].first != curScore) {
            informative = true;
            break;
          }
        }
      }
      if (!informative) {
        uninformativeReads.push_back(i);
        readTypes[i] = mgsr::readType::IDENTICAL_SCORE_ACROSS_NODES;
      }
    }

    std::cout << "There are " << uninformativeReads.size() << " reads that have identical scores across all probable nodes" << std::endl;
    std::cerr << "There are " << uninformativeReads.size() << " reads that have identical scores across all probable nodes" << std::endl;
  }

  void filter_by_mbc(
    std::vector<std::string>& nodes, Eigen::MatrixXd& probs, mgsr::ReadScores& readScores,
    const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    const std::vector<bool>& lowScoreReads, const size_t& numLowScoreReads, const std::string& excludeNode, std::vector<mgsr::readType>& readTypes,
    const std::unordered_map<std::string, double>& kminmer_binary_coverage, const int& preEMFilterMBCNum, const bool& save_kminmer_binary_coverage, const std::string& prefix
  ) {
    std::cerr << "Filter method mbc: filter out haplotypes that do not have a unique best read score" << std::endl;

    std::vector<std::pair<std::string, double>> kminmer_binary_coverage_vec;
    for (const auto& node : kminmer_binary_coverage) {
      if (leastRecentIdenticalAncestor.find(node.first) != leastRecentIdenticalAncestor.end()) continue;
      double curCoverage = node.second;
      if (identicalSets.find(node.first) != identicalSets.end()) {
        for (const auto& identicalNode : identicalSets.at(node.first)) {
          if (kminmer_binary_coverage.at(identicalNode) > curCoverage) {
            curCoverage = kminmer_binary_coverage.at(identicalNode);
          }
        }
      }
      kminmer_binary_coverage_vec.emplace_back(std::make_pair(node.first, curCoverage));
    }

    std::sort(kminmer_binary_coverage_vec.begin(), kminmer_binary_coverage_vec.end(), [](const auto& a, const auto& b) {
      return a.second > b.second;
    });

    if (save_kminmer_binary_coverage) {
      std::ofstream kminmer_binary_coverage_file(prefix + ".kminmer_binary_coverage.txt");
      for (const auto& [node, coverage] : kminmer_binary_coverage_vec) {
        kminmer_binary_coverage_file << node;
        if (identicalSets.find(node) != identicalSets.end()) {
          for (const auto& identicalNode : identicalSets.at(node)) {
            kminmer_binary_coverage_file << "," << identicalNode;
          }
        }
        kminmer_binary_coverage_file << "\t" << coverage << std::endl;
      }
    }

    std::vector<std::string> probableNodes;
    int numProbableNodes = 0;
    for (const auto& [node, coverage] : kminmer_binary_coverage_vec) {
      if (coverage == 1.0) {
        probableNodes.push_back(node);
      } else if (numProbableNodes < preEMFilterMBCNum) {
        probableNodes.push_back(node);
        ++numProbableNodes;
      }
    }

    std::vector<std::pair<std::string, std::vector<std::pair<int32_t, double>>>> probableNodeScores;
    for (const auto& node : probableNodes) {
      const auto& curNodeScores = readScores.getScoresAtNode(node);
      probableNodeScores.emplace_back(std::make_pair(node, curNodeScores));
    }

    exclude_noninformative_reads(readTypes, probableNodeScores);

    size_t numExcludedReads = readTypes.size() - std::count(readTypes.begin(), readTypes.end(), mgsr::readType::PASS);

    std::cout << "Excluding " << numExcludedReads << " reads in total" << std::endl;
    std::cerr << "Excluding " << numExcludedReads << " reads in total" << std::endl;
    
    probs.resize(readScores.scores.size() - numExcludedReads, probableNodes.size());

    size_t colIndex = 0;
    for (const auto& [node, scores] : probableNodeScores) {
      if (leastRecentIdenticalAncestor.find(node) != leastRecentIdenticalAncestor.end()) {
        std::cerr << "Error: Node " << node << " has a least recent identical ancestor." << std::endl;
        exit(1);
      }
      size_t rowIndex = 0;
      for (size_t i = 0; i < scores.size(); ++i) {
        if (readTypes[i] != mgsr::readType::PASS) continue;
        probs(rowIndex, colIndex) = scores[i].second;
        ++rowIndex;
      }
      nodes.push_back(node);
      ++colIndex;
    }

    std::cerr << "Finished mbc filter: " << nodes.size() << " nodes" << std::endl;

  }

  //squarem test 1: periodically drop nodes with very low abundance
  void squaremHelper_test_1(
    Tree *T, mgsr::ReadScores& readScores,
    const std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex, const std::vector<bool>& lowScoreReads,
    const int32_t& numReads, const size_t& numLowScoreReads, std::vector<mgsr::readType>& readTypes,
    std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors,
    std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets, Eigen::MatrixXd& probs,
    std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const std::string& preEMFilterMethod, const int& preEMFilterNOrder, const int& preEMFilterMBCNum,
    const int& emFilterRound, const int& checkFrequency, const int& removeIteration, const double& insigProp,
    const int& roundsRemove, const double& removeThreshold, const bool& leafNodesOnly, const std::unordered_map<std::string, double>& kminmer_binary_coverage,
    std::string excludeNode, const bool& save_kminmer_binary_coverage, const std::string& prefix
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

    std::cout << "pre-EM filter nodes size: " << readScores.nodeToDfsIndex.size() - leastRecentIdenticalAncestors.size() << std::endl;
    std::cerr << "pre-EM filter nodes size: " << readScores.nodeToDfsIndex.size() - leastRecentIdenticalAncestors.size() << "\n" << std::endl;
    if (preEMFilterMethod == "null") {  
      // haplotype_filter::noFilter(nodes, probs, allScores, leastRecentIdenticalAncestors, lowScoreReads, numLowScoreReads, excludeNode, excludeReads);
    } else if (preEMFilterMethod == "mbc") {
      filter_by_mbc(nodes, probs, readScores, leastRecentIdenticalAncestors, identicalSets, lowScoreReads, numLowScoreReads, excludeNode, readTypes, kminmer_binary_coverage, preEMFilterMBCNum, save_kminmer_binary_coverage, prefix);
    } else {
      std::cerr << "pre-EM filter method not recognized" << std::endl;
      exit(1);
    }
    std::string filteredNodesFile = prefix + "_filtered_nodes.txt";
    std::ofstream filtedNodesStream(filteredNodesFile);
    for (const auto& node : nodes) {
      if (leastRecentIdenticalAncestors.find(node) != leastRecentIdenticalAncestors.end()) { 
        // sanity check
        std::cerr << "Error: Node " << node << " has a least recent identical ancestor." << std::endl;
        exit(1);
      }
      filtedNodesStream << node << std::endl;
    }
    filtedNodesStream.close();

    std::cout << "post-EM filter nodes size: " << nodes.size() << std::endl;
    std::cerr << "post-EM filter nodes size: " << nodes.size() << "\n" << std::endl;
    std::cout << "post-EM filter reduced nodes size: " << readScores.nodeToDfsIndex.size() - leastRecentIdenticalAncestors.size() - nodes.size() << std::endl;
    std::cerr << "post-EM filter reduced nodes size: " << readScores.nodeToDfsIndex.size() - leastRecentIdenticalAncestors.size() - nodes.size() << "\n" << std::endl;

    props = Eigen::VectorXd::Constant(nodes.size(), 1.0 / static_cast<double>(nodes.size()));
    size_t totalNodes = readScores.nodeToDfsIndex.size() - leastRecentIdenticalAncestors.size();
    size_t numExcludedReads = readTypes.size() - std::count(readTypes.begin(), readTypes.end(), mgsr::readType::PASS);
    Eigen::VectorXd readDuplicates(readScores.scores.size() - numLowScoreReads - numExcludedReads);
    
    size_t indexReadDuplicates = 0;
    size_t numHighScoreReads = 0;
    for (size_t i = 0; i < readSeedmersDuplicatesIndex.size(); ++i) {
      if (readTypes[i] != mgsr::readType::PASS) continue;
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
    // int32_t filterRoundCount = 0;
    // size_t prefilterIterations = 5;
    // std::vector<int> insigCounts(nodes.size());
    // std::cout << "\npre-filter round for " << prefilterIterations << " iterations" << std::endl;
    // std::cerr << "\npre-filter round for " << prefilterIterations << " iterations" << std::endl;
    // llh = getExp(probs, props, readDuplicates);
    // squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, checkFrequency, prefilterIterations, insigCounts, insigProp, totalNodes);
    // while (true && nodes.size() > std::max(static_cast<int>(totalNodes) / 20, 100)) {
    //   if (filterRoundCount >= emFilterRound) break;
    //   std::cout << "\nfilter round " << filterRoundCount + 1 << " out of " << emFilterRound << std::endl;
    //   std::cerr << "\nfilter round " << filterRoundCount + 1 << " out of " << emFilterRound << std::endl;
    //   ++filterRoundCount;
    //   llh = getExp(probs, props, readDuplicates);
    //   squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, checkFrequency, std::numeric_limits<size_t>::max(), insigCounts, insigProp, totalNodes);
    //   if (converged) {
    //     break;
    //   }
    //   std::cout << "\nfiltering round " << filterRoundCount << std::endl;
    //   std::cerr << "\nfiltering round " << filterRoundCount << std::endl;
    //   std::vector<size_t> significantIndices;
    //   std::vector<std::string> sigNodes;
    //   for (size_t i = 0; i < nodes.size(); ++i) {
    //     if (insigCounts[i] < removeIteration) {
    //       significantIndices.push_back(i);
    //     }
    //   }

    //   if (significantIndices.size() == nodes.size()) {
    //     continue;
    //   }

    //   for (size_t idx : significantIndices) {
    //     sigNodes.push_back(nodes[idx]);
    //   }

    //   Eigen::MatrixXd sigProbs(probs.rows(), sigNodes.size());
    //   Eigen::VectorXd sigProps(sigNodes.size());
    //   for (size_t i = 0; i < significantIndices.size(); ++i) {
    //     sigProbs.col(i) = probs.col(significantIndices[i]);
    //     sigProps(i) = props(i);
    //   }
    //   std::cout << "dropped " << nodes.size() - sigNodes.size() << " during EM" << std::endl;
    //   std::cerr << "dropped " << nodes.size() - sigNodes.size() << " during EM" << std::endl;
    //   std::cout << sigNodes.size() << " nodes left" << std::endl;
    //   std::cerr << sigNodes.size() << " nodes left" << std::endl;
    //   nodes = sigNodes;
    //   probs = sigProbs;
    //   props = sigProps;
    //   normalize(props);
    //   insigCounts.assign(nodes.size(), 0);
    //   if (nodes.size() <= std::max(static_cast<int>(totalNodes) / 20, 100) || filterRoundCount >= emFilterRound) {
    //     break;
    //   }
    // }

    if (!converged) {
      std::vector<int> insigCounts(nodes.size());
      std::cout << "start full EM" << std::endl;
      std::cerr << "start full EM" << std::endl;
      llh = getExp(probs, props, readDuplicates);
      squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), insigCounts, insigProp, totalNodes);
      assert(converged == true);
    }
    
    std::vector<size_t> haplotypeGroupIndices;
    std::vector<double> haplotypeGroupAbundances;
    std::vector<std::pair<std::string, double>> sortedNodes(nodes.size());
    std::vector<size_t> sortedNodesIndices(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        sortedNodes.at(i) = {nodes[i], props(i)};
        sortedNodesIndices[i] = i;
    }
    std::sort(sortedNodesIndices.begin(), sortedNodesIndices.end(), [&](size_t a, size_t b) {
        return props(a) > props(b);
    });
    std::sort(sortedNodes.begin(), sortedNodes.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
        return a.second > b.second;
    });

    std::pair<size_t, std::vector<std::string>> haplotypeGroup{0, {sortedNodes[0].first}};
    haplotypeGroupAbundances.push_back(sortedNodes[0].second);
    double currGroupIndividualsAbundance = sortedNodes[0].second;
    for (size_t i = 1; i < sortedNodes.size(); ++i) {
      const auto& currNode = sortedNodes[i];
      if (currNode.second > 1e-10 && std::abs(currNode.second - currGroupIndividualsAbundance) < 1e-10) {
        haplotypeGroup.second.push_back(currNode.first);
        haplotypeGroupAbundances.back() += currNode.second;
      } else {
        std::string representativeNode = haplotypeGroup.second[0];
        for (size_t j = 1; j < haplotypeGroup.second.size(); ++j) {
          std::string currNode = haplotypeGroup.second[j];
          if (identicalSets.find(currNode) != identicalSets.end()) {
            for (const auto& identicalNode : identicalSets.at(currNode)) {
              leastRecentIdenticalAncestors[identicalNode] = representativeNode;
              identicalSets[representativeNode].insert(identicalNode);
            }
            identicalSets.erase(currNode);
          }
          leastRecentIdenticalAncestors[currNode] = representativeNode;
          identicalSets[representativeNode].insert(currNode);
        }
        haplotypeGroupIndices.push_back(haplotypeGroup.first);
        haplotypeGroup.first = i;
        haplotypeGroup.second.clear();
        haplotypeGroup.second.push_back(currNode.first);
        currGroupIndividualsAbundance = currNode.second;
        haplotypeGroupAbundances.push_back(currNode.second);
      }
    }

    if (!haplotypeGroup.second.empty()) {
      std::string representativeNode = haplotypeGroup.second[0];
      for (size_t j = 1; j < haplotypeGroup.second.size(); ++j) {
        std::string currNode = haplotypeGroup.second[j];
        if (identicalSets.find(currNode) != identicalSets.end()) {
          for (const auto& identicalNode : identicalSets.at(currNode)) {
            leastRecentIdenticalAncestors[identicalNode] = representativeNode;
            identicalSets[representativeNode].insert(identicalNode);
          }
          identicalSets.erase(currNode);
        }
        leastRecentIdenticalAncestors[currNode] = representativeNode;
        identicalSets[representativeNode].insert(currNode);
      }
      haplotypeGroupIndices.push_back(haplotypeGroup.first);
    }

    if (haplotypeGroupIndices.size() != haplotypeGroupAbundances.size()) {
      std::cerr << "Error: haplotypeGroupIndices.size() != haplotypeGroupAbundances.size()" << std::endl;
      exit(1);
    }

    std::vector<std::string> haplotypeGroupNodes(haplotypeGroupIndices.size());
    Eigen::MatrixXd haplotypeGroupProbs(probs.rows(), haplotypeGroupIndices.size());
    Eigen::VectorXd haplotypeGroupProps(haplotypeGroupIndices.size());
    for (size_t i = 0; i < haplotypeGroupIndices.size(); ++i) {
      haplotypeGroupNodes[i] = nodes[sortedNodesIndices[haplotypeGroupIndices[i]]];
      haplotypeGroupProbs.col(i) = probs.col(sortedNodesIndices[haplotypeGroupIndices[i]]);
      haplotypeGroupProps(i) = haplotypeGroupAbundances[i];
    }
    nodes = std::move(haplotypeGroupNodes);
    probs = std::move(haplotypeGroupProbs);
    props = std::move(haplotypeGroupProps);



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
      nodes = sigNodes;
      probs = sigProbs;
      props = sigProps;
    }

    if (!excludeNode.empty()) {
      Eigen::VectorXd curProbs(readScores.scores.size() - numLowScoreReads);
      size_t indexCurProbs = 0;
      const auto& curNode = readScores.getScoresAtNode(excludeNode);
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
    const std::vector<std::string>& nodes,
    const std::vector<Read>& reads,
    const std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex,
    std::unordered_map<std::string, std::vector<std::pair<size_t, int32_t>>>& assignedReads) {
      
    int32_t tolerance = 5;
    for (size_t i = 0; i < readSeedmersDuplicatesIndex.size(); ++i) {
      int32_t maxScore = std::numeric_limits<int32_t>::min();
      std::vector<std::pair<size_t, int32_t>> bestNodes;

      for (size_t j = 0; j < nodes.size(); ++j) {
        const auto& nodeScores = allScores.at(nodes[j]);
        int32_t curScore = nodeScores[i].first;
        if (curScore > maxScore) {
          maxScore = curScore;
        }
      }

      for (size_t j = 0; j < nodes.size(); ++j) {
        const auto& nodeScores = allScores.at(nodes[j]);
        int32_t curScore = nodeScores[i].first;
        if (curScore >= maxScore - tolerance) {
          bestNodes.push_back(std::make_pair(j, curScore));
        }
      }

      for (const auto& node : bestNodes) {
        const auto& nodeIdx = node.first;
        const auto& nodeScore = node.second;
        for (size_t readIdx : readSeedmersDuplicatesIndex[i]) {
          assignedReads[nodes[nodeIdx]].push_back(std::make_pair(readIdx, nodeScore));
        }
      }
    }
  }

}

#endif

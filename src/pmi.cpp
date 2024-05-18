#include "pmi.hpp"
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "tree.hpp"
#include <__config>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <sys/_types/_int32_t.h>
#include <sys/_types/_int64_t.h>
#include <unordered_set>

using namespace seeding;
using namespace pmi;
using namespace PangenomeMAT;
using namespace tree;
#include <easy/profiler.h>

/* Helpers */
void applyMutations(
    mutableTreeData &data, blockMutData_t &blockMutData,
    std::set<std::tuple<int, int, int, int>, decltype(rangeCmp)> &recompRanges,
    nucMutData_t &nucMutData, Tree *T, const Node *node,
    const globalCoords_t &globalCoords, seedIndex &index) {

  EASY_FUNCTION(profiler::colors::Magenta);

  EASY_BLOCK("applyMutations");
  blockExists_t &blockExists = data.blockExists;
  blockStrand_t &blockStrand = data.blockStrand;
  auto &coordToTuple = data.coordToTuple; // map { coordToTuple[x] => (xBlockId,
                                          // xNucPosition, xNucGapPosition) }

  
  for (const auto &mutation : node->blockMutation) {
    int32_t blockId = mutation.primaryBlockId;
    bool type = mutation.blockMutInfo;
    bool inversion = mutation.inversion;
    if (type == 1) {
      // insertion
      bool oldStrand;
      bool oldMut;
      oldStrand = blockStrand[blockId].first;
      oldMut = blockExists[blockId].first; //NICO:This should always be false???
      blockExists[blockId].first = true;
      // if insertion of inverted block takes place, the strand is backwards
      blockStrand[blockId].first = !inversion;
      blockMutData.push_back(
          std::make_tuple(blockId, oldMut, oldStrand, true, !inversion));
    } else {
      bool oldMut;
      bool oldStrand;
      if (inversion) {
        oldStrand = blockStrand[blockId].first;
        oldMut = blockExists[blockId].first;
        blockStrand[blockId].first = !oldStrand;
        blockMutData.push_back(
            std::make_tuple(blockId, oldMut, oldStrand, oldMut, !oldStrand));
      } else {
        // Actually a deletion
        oldStrand = blockStrand[blockId].first;
        oldMut = blockExists[blockId].first; //would be true
        blockExists[blockId].first = false;
        // resetting strand to true during deletion
        blockStrand[blockId].first = true;
        blockMutData.push_back(
            std::make_tuple(blockId, oldMut, oldStrand, false, true));
      }
    }
  }
  
  for (size_t i = 0; i < node->nucMutation.size(); i++) {
    int32_t blockId = node->nucMutation[i].primaryBlockId;
    int32_t nucPosition = node->nucMutation[i].nucPosition;
    int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
    uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
    char newVal = '-';
    size_t globalCoord = tree::getGlobalCoordinate(
        blockId, nucPosition, nucGapPosition, globalCoords);

    if (type < 3) { // Either S, I or D
      int len = ((node->nucMutation[i].mutInfo) >> 4);
      if (type == PangenomeMAT::NucMutationType::NS) {
        // Substitution
        if (nucGapPosition != -1) {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId]
                              .first[nucPosition]
                              .second[nucGapPosition + j];
            newVal = PangenomeMAT::getNucleotideFromCode(
                ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
            data.sequence[blockId]
                .first[nucPosition]
                .second[nucGapPosition + j] = newVal;
            nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                                 nucGapPosition + j, oldVal,
                                                 newVal, globalCoord, len));
          }
        } else {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId].first[nucPosition + j].first;
            newVal = PangenomeMAT::getNucleotideFromCode(
                ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
            data.sequence[blockId].first[nucPosition + j].first = newVal;
            nucMutData.push_back(std::make_tuple(blockId, nucPosition + j,
                                                 nucGapPosition, oldVal, newVal,
                                                 globalCoord, len));
          }
        }
      } else if (type == PangenomeMAT::NucMutationType::NI) {
        // Insertion
        if (nucGapPosition != -1) {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId]
                              .first[nucPosition]
                              .second[nucGapPosition + j];
            newVal = PangenomeMAT::getNucleotideFromCode(
                ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
            data.sequence[blockId]
                .first[nucPosition]
                .second[nucGapPosition + j] = newVal;
            nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                                 nucGapPosition + j, oldVal,
                                                 newVal, globalCoord, len));
          }
        } else {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId].first[nucPosition + j].first;
            const int nucCode =
                ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF;
            newVal = PangenomeMAT::getNucleotideFromCode(
                ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
            data.sequence[blockId].first[nucPosition + j].first = newVal;
            nucMutData.push_back(std::make_tuple(blockId, nucPosition + j,
                                                 nucGapPosition, oldVal, newVal,
                                                 globalCoord, len));
          }
        }
      } else if (type == PangenomeMAT::NucMutationType::ND) {
        // Deletion
        if (nucGapPosition != -1) {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId]
                              .first[nucPosition]
                              .second[nucGapPosition + j];
            data.sequence[blockId]
                .first[nucPosition]
                .second[nucGapPosition + j] = '-';
            nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                                 nucGapPosition + j, oldVal,
                                                 '-', globalCoord, len));
          }
        } else {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId].first[nucPosition + j].first;
            data.sequence[blockId].first[nucPosition + j].first = '-';
            nucMutData.push_back(std::make_tuple(blockId, nucPosition + j,
                                                 nucGapPosition, oldVal, '-',
                                                 globalCoord, len));
          }
        }
      }
    } else {
      int len = 0;
      if (type == PangenomeMAT::NucMutationType::NSNPS) {
        // SNP Substitution
        newVal = PangenomeMAT::getNucleotideFromCode(
            ((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (nucGapPosition != -1) {
          char oldVal =
              data.sequence[blockId].first[nucPosition].second[nucGapPosition];
          data.sequence[blockId].first[nucPosition].second[nucGapPosition] =
              newVal;
          nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                               nucGapPosition, oldVal, newVal,
                                               globalCoord, len));
        } else {
          char oldVal = data.sequence[blockId].first[nucPosition].first;
          data.sequence[blockId].first[nucPosition].first = newVal;
          nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                               nucGapPosition, oldVal, newVal,
                                               globalCoord, len));
        }
      } else if (type == PangenomeMAT::NucMutationType::NSNPI) {
        // SNP Insertion
        len = 1;
        newVal = PangenomeMAT::getNucleotideFromCode(
            ((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (nucGapPosition != -1) {
          char oldVal =
              data.sequence[blockId].first[nucPosition].second[nucGapPosition];
          data.sequence[blockId].first[nucPosition].second[nucGapPosition] =
              newVal;
          nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                               nucGapPosition, oldVal, newVal,
                                               globalCoord, len));
        } else {
          char oldVal = data.sequence[blockId].first[nucPosition].first;
          data.sequence[blockId].first[nucPosition].first = newVal;
          nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                               nucGapPosition, oldVal, newVal,
                                               globalCoord, len));
        }
      } else if (type == PangenomeMAT::NucMutationType::NSNPD) {
        // SNP Deletion
        if (nucGapPosition != -1) {
          char oldVal =
              data.sequence[blockId].first[nucPosition].second[nucGapPosition];
          data.sequence[blockId].first[nucPosition].second[nucGapPosition] =
              '-';
          nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                               nucGapPosition, oldVal, '-',
                                               globalCoord, len));
        } else {
          char oldVal = data.sequence[blockId].first[nucPosition].first;
          data.sequence[blockId].first[nucPosition].first = '-';
          nucMutData.push_back(std::make_tuple(blockId, nucPosition,
                                               nucGapPosition, oldVal, '-',
                                               globalCoord, len));
        }
      }
    }
  }
  

  // block seed mutations
  for (const auto &mutation : node->blockMutation) {
    // for now, any block mutation triggers recomputation of whole block
    int32_t blockId = mutation.primaryBlockId;
    if (!data.blockExists[blockId].first) {
      //remove seeds 
      continue;
    }
    bool type = mutation.blockMutInfo;
    bool inversion = mutation.inversion;

    int64_t mutAffectedStart = std::max(
        (int64_t)0, tree::getGlobalCoordinate(blockId, 0, 0, globalCoords)); // global start of block

    int64_t mutAffectedEnd =
        std::min(data.maxGlobalCoordinate,
                 tree::getGlobalCoordinate(
                     blockId, data.sequence[blockId].first.size() - 1, -1,
                     globalCoords));                      // global end of block
                     
    int64_t mutRecomputeEnd = mutAffectedEnd; // the substring to recompute ends

    recompRanges.insert(std::make_tuple(mutAffectedStart, mutAffectedEnd,
                                        mutAffectedStart, mutRecomputeEnd));
  
  }
  
  // Nuc seed mutations
  for (size_t i = 0; i < node->nucMutation.size(); i++) {
    int32_t blockId = node->nucMutation[i].primaryBlockId;
    int32_t nucPosition = node->nucMutation[i].nucPosition;
    int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
    uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
    size_t globalCoord = std::max(
        (int64_t)0, tree::getGlobalCoordinate(blockId, nucPosition,
                                              nucGapPosition, globalCoords));

    if (type < 3) { // Either S, I or D
      int len = ((node->nucMutation[i].mutInfo) >> 4);
      int64_t mutAffectedStart =
          globalCoord; // the mutation affects all bases from here
      int64_t mutAffectedEnd =
          globalCoord + len + 2; // the mutation affects all bases up to here
      int64_t mutRecomputeEnd =
          mutAffectedEnd; // the substring to recompute ends here (this is
                          // modified later)
      mutRecomputeEnd++;
      recompRanges.insert(std::make_tuple(mutAffectedStart, mutAffectedEnd,
                                          mutAffectedStart, mutRecomputeEnd));
    } else {
      int64_t mutAffectedStart =
          globalCoord; // the mutation affects all bases from here
      int64_t mutAffectedEnd =
          globalCoord + 1; // the mutation affects all bases up to here
      int64_t mutRecomputeEnd =
          mutAffectedEnd; // the substring to recompute ends here (this is

      recompRanges.insert(std::make_tuple(mutAffectedStart, mutAffectedEnd,
                                          mutAffectedStart, mutRecomputeEnd));
    }
  }
}

void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T,
                   const Node *node, const blockMutData_t &blockMutData,
                   const nucMutData_t &nucMutData) {
  EASY_FUNCTION(profiler::colors::Green);
  blockExists_t &blockExists = data.blockExists;
  blockStrand_t &blockStrand = data.blockStrand;

  for (auto it = blockMutData.rbegin(); it != blockMutData.rend(); it++) {
    auto mutation = *it;

    blockExists[std::get<0>(mutation)].first = std::get<1>(mutation);
    blockStrand[std::get<0>(mutation)].first = std::get<2>(mutation);
  }

  // Undo nuc mutations when current node and its subtree have been
  // processed
  for (auto it = nucMutData.rbegin(); it != nucMutData.rend(); it++) {
    auto mutation = *it;
    if (std::get<2>(mutation) != -1) {
      data.sequence[std::get<0>(mutation)]
          .first[std::get<1>(mutation)]
          .second[std::get<2>(mutation)] = std::get<3>(mutation);
    } else {
      data.sequence[std::get<0>(mutation)].first[std::get<1>(mutation)].first =
          std::get<3>(mutation);
    }
  }
}

std::set<std::tuple<int, int, int, int>, decltype(rangeCmp)>
mergeOverlappingRanges(const std::set<std::tuple<int, int, int, int>,
                                      decltype(rangeCmp)> &ranges) {
  std::set<std::tuple<int, int, int, int>, decltype(rangeCmp)> mergedRanges(
      rangeCmp);
  if (ranges.empty())
    return mergedRanges;

  auto it = ranges.begin();
  int start1 = std::get<0>(*it);
  start1 = std::max(start1, 0);
  int end1 = std::get<1>(*it);
  int start2 = std::get<2>(*it);
  start2 = std::max(start2, 0);
  int end2 = std::get<3>(*it);

  for (++it; it != ranges.end(); ++it) {
    int currentStart1 = std::get<0>(*it);
    int currentEnd1 = std::get<1>(*it);
    int currentStart2 = std::get<2>(*it);
    int currentEnd2 = std::get<3>(*it);

    if (currentStart2 <= end2) {
      end2 = std::max(end2, currentEnd2);
      end1 = std::max(end1, currentEnd1);
    } else {
      mergedRanges.emplace(start1, end1, start2, end2);
      start1 = currentStart1;
      end1 = currentEnd1;
      start2 = currentStart2;
      end2 = currentEnd2;
    }
  }
  mergedRanges.emplace(start1, end1, start2, end2);

  return mergedRanges;
}

void buildHelper(mutableTreeData &data, seedMap_t &seedMap, seedIndex &index,
                 Tree *T, const Node *node,
                 const globalCoords_t &globalCoords) {
  EASY_FUNCTION(profiler::colors::Blue);
  EASY_BLOCK("buildHelper");
  blockMutData_t blockMutData;
  nucMutData_t nucMutData;
  auto &coordToTuple = data.coordToTuple;
  std::set<std::tuple<int, int, int, int>, decltype(rangeCmp)> recompRanges;
  /* Mutate with block and nuc mutations. */
  
  applyMutations(data, blockMutData, recompRanges, nucMutData, T, node,
                 globalCoords, index);
  
  // std::cout << "node " << node->identifier << std::endl;
  // std::cout << "block mutations:\n";
  // for (const auto &blockMut : blockMutData) {
  //   std::cout << std::get<0>(blockMut) << " " <<
  //   std::get<1>(blockMut)
  //   << " "
  //             << std::get<2>(blockMut) << " " <<
  //             std::get<3>(blockMut)
  //             << " "
  //             << std::get<4>(blockMut) << std::endl;
  // }
  // std::cout << "nuc mutations:\n";
  // for (const auto &nucMut : nucMutData) {
  //   std::cout << std::get<0>(nucMut) << " " << std::get<1>(nucMut) <<
  //   "
  //   "
  //             << std::get<2>(nucMut) << " " << std::get<3>(nucMut) <<
  //             "
  //             "
  //             << std::get<4>(nucMut) << " " << std::get<5>(nucMut) <<
  //             "
  //             "
  //             << std::get<6>(nucMut) << std::endl;
  // }
  std::set<std::string> outDeletions;
  std::set<std::string> outInsertions;
  std::vector<int64_t> seedsToClear;
  std::vector<std::pair<int64_t, std::string>> backtrack;

  auto mergedRanges = mergeOverlappingRanges(recompRanges);
  for (const auto &recompRange : mergedRanges) {
  }
  std::set<std::tuple<int, int, int, int>, decltype(rangeCmp)> ranges2 =
      recompRanges;

  for (const auto &recompRange : mergedRanges | std::views::reverse) {
    int64_t substringStart = std::get<2>(recompRange);
    int64_t substringStop = std::get<3>(recompRange);

    int64_t affectedStart = substringStart;
    int64_t affectedStop = substringStop;

    auto &cta = coordToTuple[substringStart];
    auto &ctb = coordToTuple[substringStop];

    int32_t seenNongapRight = 0;
    int32_t seenNongapLeft = 0;
    while (substringStart >= 0 && seenNongapLeft < index.k * index.j) {
      //while (!data.blockExists[std::get<0>(cta)].first && substringStart > 0) {
      //  substringStart--;
      //  cta = coordToTuple[substringStart];
      //}
      if (seenNongapLeft <= index.k * index.j) {
        affectedStart--;
      }
      if (data.sequence[std::get<0>(cta)].first[std::get<1>(cta)].first !=
              'x' &&
          data.sequence[std::get<0>(cta)].first[std::get<1>(cta)].first !=
              '-'  && data.blockExists[std::get<0>(cta)].first) {
        seenNongapLeft++;
      }
      for (int k = data.sequence[std::get<0>(cta)]
                       .first[std::get<1>(cta)]
                       .second.size() -
                   1;
           k >= 0; k--) {
        if (data.sequence[std::get<0>(cta)].first[std::get<1>(cta)].second[k] !=
                'x' &&
            data.sequence[std::get<0>(cta)].first[std::get<1>(cta)].second[k] !=
                '-') {
          seenNongapLeft++;
        }
      }
      cta = coordToTuple[substringStart];
      substringStart--;
    }
    while (substringStop <= data.maxGlobalCoordinate &&
           seenNongapRight < index.k * index.j) {
      while (!data.blockExists[std::get<0>(ctb)].first &&
             substringStop < data.maxGlobalCoordinate) {
        substringStop++;
        ctb = coordToTuple[substringStop];
      }
      if (seenNongapRight <= index.k * index.j) {
        affectedStop++;
      }
      if (data.sequence[std::get<0>(ctb)].first[std::get<1>(ctb)].first !=
              'x' &&
          data.sequence[std::get<0>(ctb)].first[std::get<1>(ctb)].first !=
              '-') {
        seenNongapRight++;
      }
      for (int k = 0; k < data.sequence[std::get<0>(ctb)]
                              .first[std::get<1>(ctb)]
                              .second.size();
           k++) {
        if (data.sequence[std::get<0>(ctb)].first[std::get<1>(ctb)].second[k] !=
                'x' &&
            data.sequence[std::get<0>(ctb)].first[std::get<1>(ctb)].second[k] !=
                '-') {
          seenNongapRight++;
        }
      }
      ctb = coordToTuple[substringStop];
      substringStop++;
    }
    substringStart = std::max(substringStart, (int64_t)0);
    cta = coordToTuple[substringStart];
    ctb = coordToTuple[substringStop];
    ranges2.insert(std::make_tuple(affectedStart, affectedStop, substringStart,
                                   substringStop));
  }
  auto merged2 = mergeOverlappingRanges(ranges2);
  for (const auto &recompRange : merged2 | std::views::reverse) {
    if (std::get<2>(recompRange) > std::get<3>(recompRange)
        || std::get<0>(recompRange) > std::get<1>(recompRange)) {
      //std::cout << "skipping invalid range!" << std::endl;
      continue;
    }
    int64_t rStart = std::get<2>(recompRange);
    int64_t rStop = std::get<3>(recompRange);
    if (rStop >= data.maxGlobalCoordinate) {
      rStop = data.maxGlobalCoordinate;
    }
    if (rStart < 0) {
      rStart = 0;
    }
    auto &cta = coordToTuple[rStart];
    auto &ctb = coordToTuple[rStop];
    int64_t affectedStart = std::get<0>(recompRange);
    int64_t affectedStop = std::get<1>(recompRange);
    std::tuple<int, int, int, int> pmStart = {std::get<0>(cta), -1, 0, -1};
    std::tuple<int, int, int, int> pmStop = {std::get<0>(ctb), -1,
                                             std::get<1>(ctb), -1};
    /*std::cout << "pmStart (" << std::get<0>(pmStart) << ", "
              << std::get<1>(pmStart) << ", " << std::get<2>(pmStart) << ", "
              << std::get<3>(pmStart) << ") pmStop (" << std::get<0>(pmStop)
              << ", " << std::get<1>(pmStop) << ", " << std::get<2>(pmStop)
              << ", " << std::get<3>(pmStop) << ") " << std::endl;*/
    int64_t startCoord =
        tree::getGlobalCoordinate(std::get<0>(pmStart), std::get<2>(pmStart),
                                  data.sequence[std::get<0>(pmStart)].first[std::get<2>(pmStart)].second.size() > 0
                                    ? std::get<3>(pmStart) : -1, globalCoords);
    int64_t stopCoord =
        tree::getGlobalCoordinate(std::get<0>(pmStop), std::get<2>(pmStop),
                                  -1, globalCoords);
    if (stopCoord < startCoord) {
      stopCoord = data.maxGlobalCoordinate;
    }
    /*std::cout << "startCoord " << startCoord << " stopCoord " << stopCoord
              << std::endl;*/
    std::string recomputeSeq = tree::getNucleotideSequenceFromBlockCoordinates(
        pmStart, pmStop, data.sequence, data.blockExists, data.blockStrand, T,
        node);
    /*std::cout << "RECOMP:" << recomputeSeq << std::endl;*/
    // find the next seed position downstream of the affected range
    int64_t lastDownstreamSeedPos = affectedStop;
    auto boundItr = seedMap.upper_bound(affectedStop);
    // iterator to first element with key > affectedStop if it exists
    if (boundItr == seedMap.end()) {
      //std::cout << " no downstream seed\n" << std::endl;
      lastDownstreamSeedPos = -1;
    } else {
      lastDownstreamSeedPos = boundItr->first;
      //std::cout << " downstream seed at " << lastDownstreamSeedPos << std::endl;
    }
    //todo this 2 is a hack
    for (int64_t i = recomputeSeq.size() - 1; i >= 0; i--) {
      int64_t c = startCoord + i;
      auto tup = coordToTuple[c];
      if (!data.blockExists[std::get<0>(tup)].first) {
        continue;
      }
      /*std::cout << "c: " << c << " = " << startCoord << "+ i: " << i << "(" << std::get<0>(tup) << ", "
                << std::get<1>(tup) << ", " << std::get<2>(tup) << ") "
                << recomputeSeq[i] << std::endl;*/
      if (recomputeSeq[i] == '-' || recomputeSeq[i] == 'x') {
        if (seedMap.find(c) != seedMap.end()) {
          // no longer a seed -> delete
          //std::cout << "bt " << c << " " << seedMap[c] << std::endl;
          backtrack.push_back(std::make_pair(c, seedMap[c]));
          std::string str = "";
          str += std::to_string(c);
          str += ":@";
          str += seedMap[c];
          if (seedMap[c].size() == index.j * index.k) {
            outDeletions.insert(str);
          }
          //std::cout << " del " << str << std::endl;
          seedsToClear.push_back(c);
        }
        continue;
      }
      // get the next k non-gap bases
      std::string kmer = "";
      int64_t seen_k = 0;
      int64_t k_pos = i;
      while (seen_k < index.k && k_pos < recomputeSeq.size()) {
        if (recomputeSeq[k_pos] != '-' && recomputeSeq[k_pos] != 'x'
          && data.blockExists[std::get<0>(tup)].first) {
          kmer += recomputeSeq[k_pos];
          seen_k++;
        }
        k_pos++;
      }
      //std::cout << "kmer " << kmer << std::endl;
      if (kmer.size() < index.k) {
        //std::cout << " 2small" << std::endl;
        continue;
      }
      if (seedMap.find(c) != seedMap.end()) {
        // non gap position and kmer is already a seed.
        std::string prevseedmer =
            lastDownstreamSeedPos != -1 ? seedMap[lastDownstreamSeedPos] : "";
        //std::cout << "bt " << c << " " << seedMap[c] << std::endl;
        backtrack.push_back(std::make_pair(c, seedMap[c]));

        if (seeding::is_syncmer(kmer, index.s, false)) {
          // Is it still a seed?
          //std::cout << "still seed" << std::endl;
          seedMap[c] = kmer + prevseedmer.substr(0, (index.j - 1) * index.k);
          //std::cout << "set seedMap[" << c << "] = " << seedMap[c] << std::endl;
          lastDownstreamSeedPos = c;
          if (seedMap[c].size() == index.j * index.k) {
            std::string str = "";
            str = std::to_string(c);
            str += ":";
            str += seedMap[c];
            outInsertions.insert(str);
          }
        } else {
          // no longer a seed -> delete
          //std::cout << "no longer seed, del " << c << " " << seedMap[c] << std::endl;
          seedsToClear.push_back(c);
          std::string str = "";
          str += std::to_string(c);
          str += ":@";
          str += seedMap[c];
          if (seedMap[c].size() == index.j * index.k) {
            outDeletions.insert(str);
          }
        }
      } else {
        //std::cout << "bt " << c << " " << std::endl;
        // not in seed map, could be a seed now
        backtrack.push_back(std::make_pair(c, ""));

        if (seeding::is_syncmer(kmer, index.s, false)) {
          //std::cout << "new seed " << c << " " << kmer << std::endl;
          std::string prevseedmer =
              lastDownstreamSeedPos != -1 ? seedMap[lastDownstreamSeedPos] : "";
          seedMap[c] = kmer + prevseedmer.substr(0, (index.j - 1) * index.k);
          //std::cout << "set seedMap[" << c << "] = " << seedMap[c] << std::endl;
          if (prevseedmer.size() == index.j * index.k) {
            std::string str = "";
            str += std::to_string(c);
            str += ":";
            str += seedMap[c];
            outInsertions.insert(str);
          }
          lastDownstreamSeedPos = c;
        }
      }
    }
  }
  for (const auto &pos : seedsToClear) {
    seedMap.erase(pos);
  }
  index.outStream << node->identifier << " ";
  for (const std::string &s : outDeletions) {
    index.outStream << s << " ";
  }
  for (const std::string &s : outInsertions) {
    index.outStream << s << " ";
  }
  index.outStream << "\n";

  /* Recursive step */
  for (Node *child : node->children) {
    buildHelper(data, seedMap, index, T, child, globalCoords);
  }
  for (const auto &back : backtrack) {
    if (back.second == "") {
      seedMap.erase(back.first);
    } else {
      seedMap[back.first] = back.second;
    }
  }
  /* Undo seed and sequence mutations when backtracking */
  undoMutations(data, index, T, node, blockMutData, nucMutData);
}

/* Interface implementation */
void pmi::build(seedIndex &index, Tree *T, const size_t j, const size_t k,
                const size_t s) {
  EASY_FUNCTION(profiler::colors::Magenta);
  /* Setup for seed indexing */
  tree::mutableTreeData data;
  tree::globalCoords_t globalCoords;

  
  tree::setup(data, globalCoords, T);

  index.j = j;
  index.k = k;
  index.s = s;

  seedMap_t seedMap;
  index.outStream << k << " " << s << " " << j << "\n";
  // print block seq
  for (size_t i = 0; i < data.sequence.size(); i++) {
    //std::cout << ">block " << i << " len " << data.sequence[i].first.size()<< std::endl;
    for (size_t j = 0; j < data.sequence[i].first.size(); j++) {
      for (size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
        // std::cout << "(" << i << "," << j << "," << k
        //           << "): " << data.sequence[i].first[j].second[k] << " glob:
        //           "
        //           << tree::getGlobalCoordinate(i, j, k, globalCoords)
        //           << std::endl;
      }
      // std::cout << "(" << i << "," << j
      //           << ",-1): " << data.sequence[i].first[j].first << " glob: "
      //           << tree::getGlobalCoordinate(i, j, -1, globalCoords)
      //           << std::endl;

      // std::cout << std::endl;
    }
  }
  /* Recursive traversal of tree to build the index */
  
  buildHelper(data, seedMap, index, T, T->root, globalCoords);

}
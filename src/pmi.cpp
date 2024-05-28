#include "pmi.hpp"
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "tree.hpp"
#include <__config>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string>
#include <sys/_types/_int32_t.h>
#include <sys/_types/_int64_t.h>
#include <unordered_set>

using namespace seeding;
using namespace pmi;
using namespace PangenomeMAT;
using namespace tree;

/* Helpers */
void applyMutations(mutableTreeData &data, seedMap_t &seedMap,
                    blockMutData_t &blockMutData,
                    std::vector<tupleRange> &recompRanges,
                    nucMutData_t &nucMutData, Tree *T, const Node *node,
                    const globalCoords_t &globalCoords, seedIndex &index) {

  blockExists_t &blockExists = data.blockExists;
  blockStrand_t &blockStrand = data.blockStrand;

  for (const auto &mutation : node->blockMutation) {
    int32_t blockId = mutation.primaryBlockId;
    bool type = mutation.blockMutInfo;
    bool inversion = mutation.inversion;
    if (type == 1) {
      // insertion
      bool oldStrand;
      bool oldMut;
      oldStrand = blockStrand[blockId].first;
      oldMut = blockExists[blockId].first; // NICO:This should always be
                                           // false???
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
        oldMut = blockExists[blockId].first; // would be true
        blockExists[blockId].first = false;
        // resetting strand to true during deletion
        blockStrand[blockId].first = true;
        blockMutData.push_back(
            std::make_tuple(blockId, oldMut, oldStrand, false, true));
      }
    }
    int32_t bid = blockId == 0 ? 0 : blockId - 1;
    int32_t bid2 = blockId == data.sequence.size() - 1 ? blockId : blockId + 1;
    int32_t blen = data.sequence[bid].first.size();
    recompRanges.push_back(
        {tupleCoord_t{bid, blen - 1, -1},
         tupleCoord_t{bid2, 0,
                      data.sequence[bid2].first[0].second.empty() ? -1 : 0}});
  }
  for (size_t i = 0; i < node->nucMutation.size(); i++) {
    int32_t blockId = node->nucMutation[i].primaryBlockId;
    int32_t nucPosition = node->nucMutation[i].nucPosition;
    int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
    uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
    char newVal = '-';

    if (type < 3) { // Either S, I or D
      int len = ((node->nucMutation[i].mutInfo) >> 4);

      recompRanges.push_back(
          {{blockId, nucPosition, nucGapPosition},
           {blockId,
            std::min(nucPosition + len,
                     (int32_t)data.sequence[blockId].first.size() - 1),
            -1}});

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
                                                 nucGapPosition + j, oldVal));
          }
        } else {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId].first[nucPosition + j].first;
            newVal = PangenomeMAT::getNucleotideFromCode(
                ((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
            data.sequence[blockId].first[nucPosition + j].first = newVal;
            nucMutData.push_back(std::make_tuple(blockId, nucPosition + j,
                                                 nucGapPosition, oldVal));
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
                                                 nucGapPosition + j, oldVal));
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
                                                 nucGapPosition, oldVal));
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
                                                 nucGapPosition + j, oldVal));
          }
        } else {
          for (int j = 0; j < len; j++) {
            char oldVal = data.sequence[blockId].first[nucPosition + j].first;
            data.sequence[blockId].first[nucPosition + j].first = '-';
            nucMutData.push_back(std::make_tuple(blockId, nucPosition + j,
                                                 nucGapPosition, oldVal));
          }
        }
      }
    } else {
      int len = 0;
      recompRanges.push_back(
          {{blockId, nucPosition, nucGapPosition}, {blockId, nucPosition, -1}});

      if (type == PangenomeMAT::NucMutationType::NSNPS) {
        // SNP Substitution
        newVal = PangenomeMAT::getNucleotideFromCode(
            ((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (nucGapPosition != -1) {
          char oldVal =
              data.sequence[blockId].first[nucPosition].second[nucGapPosition];
          data.sequence[blockId].first[nucPosition].second[nucGapPosition] =
              newVal;
          nucMutData.push_back(
              std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal));
        } else {
          char oldVal = data.sequence[blockId].first[nucPosition].first;
          data.sequence[blockId].first[nucPosition].first = newVal;
          nucMutData.push_back(
              std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal));
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
          nucMutData.push_back(
              std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal));
        } else {
          char oldVal = data.sequence[blockId].first[nucPosition].first;
          data.sequence[blockId].first[nucPosition].first = newVal;
          nucMutData.push_back(
              std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal));
        }
      } else if (type == PangenomeMAT::NucMutationType::NSNPD) {
        // SNP Deletion
        if (nucGapPosition != -1) {
          char oldVal =
              data.sequence[blockId].first[nucPosition].second[nucGapPosition];
          data.sequence[blockId].first[nucPosition].second[nucGapPosition] =
              '-';
          nucMutData.push_back(
              std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal));
        } else {
          char oldVal = data.sequence[blockId].first[nucPosition].first;
          data.sequence[blockId].first[nucPosition].first = '-';
          nucMutData.push_back(
              std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal));
        }
      }
    }
  }
}

void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T,
                   const Node *node, const blockMutData_t &blockMutData,
                   const nucMutData_t &nucMutData) {
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
class CoordNavigator {
public:
  CoordNavigator(const sequence_t &sequence) : sequence(sequence) {}

  bool isGap(const tupleCoord_t &coord) const {
    char c;
    if (coord.nucGapPos == -1) {
      c = sequence[coord.blockId].first[coord.nucPos].first;
    } else {
      c = sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
    }
    return c == '-' || c == 'x';
  }

  tupleCoord_t increment(const tupleCoord_t &coord) const {
    if (coord == tupleCoord_t{-1, -1, -1}) {
      return tupleCoord_t(sequence.size() - 1, sequence.back().first.size() - 1,
                          -1);
    }
    if (coord.nucGapPos != -1) {
      if (coord.nucGapPos + 1 <
          (int32_t)sequence[coord.blockId].first[coord.nucPos].second.size()) {
        return tupleCoord_t(coord.blockId, coord.nucPos, coord.nucGapPos + 1);
      } else {
        return tupleCoord_t(coord.blockId, coord.nucPos, -1);
      }
    }
    if (coord.nucPos + 1 < (int32_t)sequence[coord.blockId].first.size()) {
      if (sequence[coord.blockId].first[coord.nucPos + 1].second.empty()) {
        return tupleCoord_t(coord.blockId, coord.nucPos + 1, -1);
      } else {
        return tupleCoord_t(coord.blockId, coord.nucPos + 1, 0);
      }
    }
    if (coord.blockId + 1 < (int32_t)sequence.size()) {
      if (sequence[coord.blockId + 1].first[0].second.empty()) {
        return tupleCoord_t(coord.blockId + 1, 0, -1);
      } else {
        return tupleCoord_t(coord.blockId + 1, 0, 0);
      }
    }
  }

  tupleCoord_t decrement(const tupleCoord_t &coord) const {
    if (coord.blockId == 0 && coord.nucPos == 0 && coord.nucGapPos == 0) {
      return tupleCoord_t{0, 0, 0};
    }
    if (coord.blockId >= 0) {
      if (coord.nucPos > 0) {
        if (coord.nucGapPos > 0) {
          return tupleCoord_t(coord.blockId, coord.nucPos, coord.nucGapPos - 1);
        } else if (coord.nucGapPos == 0) {
            return tupleCoord_t(coord.blockId, coord.nucPos-1, -1);
        } else { // coord.nucGapPos == -1
          if (sequence[coord.blockId].first[coord.nucPos].second.empty()) {
            return tupleCoord_t(coord.blockId, coord.nucPos-1, -1);
          } else {
            return tupleCoord_t(coord.blockId, coord.nucPos,
                                sequence[coord.blockId].first[coord.nucPos].second.size() - 1);
          }
          return tupleCoord_t(coord.blockId, coord.nucPos - 1, -1);
        }
      } else if (coord.nucPos == 0) {
        if (coord.nucGapPos > 0) {
          return tupleCoord_t(coord.blockId, coord.nucPos, coord.nucGapPos - 1);
        } else if (coord.nucGapPos == 0) {
          return tupleCoord_t(coord.blockId-1, sequence[coord.blockId-1].first.size()-1, -1);
        } else { // coord.nucGapPos == -1
          return tupleCoord_t(coord.blockId-1, sequence[coord.blockId-1].first.size()-1, -1);
        }
      } else { // coord.nucPos == -1, shouldn't happen
        return tupleCoord_t(coord.blockId-1, sequence[coord.blockId-1].first.size()-1, -1);
      }
    }
    return tupleCoord_t{0,0,0};
  }

private:
  const sequence_t &sequence;
};

tupleCoord_t expandLeft(const CoordNavigator &navigator, tupleCoord_t coord,
                        int neededNongap, blockExists_t &blockExists) {
  int count = 0;
  // std::cout << "expandLeft..." << std::endl;
  // std::cout << "orig coord: (" << coord.blockId << ", " << coord.nucPos << ",
  // "
  //           << coord.nucGapPos << ") to => ... " << std::endl;
  while (count < neededNongap && coord >= tupleCoord_t{0, 0, 0} &&
         coord < tupleCoord_t{-1, -1, -1}) {
    if (!blockExists[coord.blockId].first) {
      coord = navigator.decrement(coord);
      continue;
    }
    // std::cout << "isGap? " << std::endl;
    if (!navigator.isGap(coord))
      count++;
    if (count >= neededNongap || coord.blockId == -1)
      return coord;
    tupleCoord_t prev = coord;

    // std::cout << " from (" << coord.blockId << ", " << coord.nucPos << ", "
    // << coord.nucGapPos << ") to => ";

    coord = navigator.decrement(coord);
  }
  return coord;
}

tupleCoord_t expandRight(const CoordNavigator &navigator, tupleCoord_t coord,
                         int neededNongap, blockExists_t &blockExists) {
  int count = 0;
  while (count < neededNongap && coord < tupleCoord_t{-1, -1, -1} &&
         coord > tupleCoord_t{0, 0, 0}) {
    if (!blockExists[coord.blockId].first) {
      // std::cout << "++ (" << coord.blockId << ", " << coord.nucPos << ", " <<
      // coord.nucGapPos << ") ";
      coord = navigator.increment(coord);
      // std::cout << "=> (" << coord.blockId << ", " << coord.nucPos << ", " <<
      // coord.nucGapPos << ")" << std::endl;
      continue;
    }
    if (!navigator.isGap(coord))
      count++;
    if (count >= neededNongap || coord.blockId == -1 ||
        coord.blockId > blockExists.size())
      return coord;
    tupleCoord_t prev = coord;
    coord = navigator.increment(coord);
  }
  return coord;
}

std::vector<tupleRange> expandAndMergeRanges(const CoordNavigator &navigator,
                                             std::vector<tupleRange> &ranges,
                                             int neededNongap,
                                             blockExists_t &blockExists) {
  if (ranges.empty())
    return {};

  std::vector<tupleRange> merged;
  tupleRange current = ranges[0];
  tupleCoord_t tmpStop = current.stop;
  current.start =
      expandLeft(navigator, current.start, neededNongap, blockExists);
  current.stop =
      expandRight(navigator, current.stop, neededNongap, blockExists);
  merged.push_back(current);

  for (size_t i = 1; i < ranges.size(); ++i) {
    // std::cout << "Merging range " << i << " which is " <<
    // ranges[i].start.blockId << ", " << ranges[i].start.nucPos << ", " <<
    // ranges[i].start.nucGapPos << " to " << ranges[i].stop.blockId << ", " <<
    // ranges[i].stop.nucPos << ", " << ranges[i].stop.nucGapPos << std::endl;
    tupleRange expandedRange = {
        expandLeft(navigator, ranges[i].start, neededNongap, blockExists),
        expandRight(navigator, ranges[i].stop, neededNongap, blockExists),
    };
    if (expandedRange.start == tupleCoord_t{-1, -1, -1}) {
      expandedRange.start = tupleCoord_t{0, 0, 0};
    }
    if (expandedRange.start <= current.stop) {
      tmpStop = current.stop;
      current.stop = std::max(current.stop, expandedRange.stop);
    } else {
      merged.push_back(current);
      tmpStop = current.stop;
      current = expandedRange;
    }
  }
  merged.push_back(current);
  return merged;
}

int64_t tupleToScalarCoord(const tupleCoord_t &coord,
                           const globalCoords_t &globalCoords) {
  if (coord.nucGapPos >= 0) {
    return globalCoords[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
  }
  return globalCoords[coord.blockId].first[coord.nucPos].first;
}

int32_t nodesDone = 0;
void buildHelper(mutableTreeData &data, seedMap_t &seedMap, seedIndex &index,
                 Tree *T, const Node *node, const globalCoords_t &globalCoords,
                 const std::vector<tupleCoord_t> &altGlobalCoords,
                 CoordNavigator &navigator) {
  blockMutData_t blockMutData;
  nucMutData_t nucMutData;

  // std::cout << "buildHelper: " << node->identifier << std::endl;

  std::vector<tupleRange> recompRanges;

  applyMutations(data, seedMap, blockMutData, recompRanges, nucMutData, T, node,
                 globalCoords, index);
  // std::cout << "applied mutations." << std::endl;
  std::set<std::string> outDeletions;
  std::set<std::string> outInsertions;
  std::vector<tupleCoord_t> seedsToClear;
  std::vector<std::pair<tupleCoord_t, std::string>> backtrack;

  std::sort(recompRanges.begin(), recompRanges.end());
  // // handle too large ranges
  // std::cout << "Sequence: " << std::endl;
  // for (size_t i = 0; i < data.sequence.size(); i++) {
  //   std::cout << "Block " << i << ": ";
  //   for (size_t j = 0; j < data.sequence[i].first.size(); j++) {
  //     for (size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
  //       std::cout << "(" << i << ", " << j << ", " << k << "): "
  //                 << data.sequence[i].first[j].second[k] << " " << globalCoords[i].first[j].second[k] << std::endl;
  //     }
  //     //std::cout << "(" << i << ", " << j << "): " << data.sequence[i].first[j].first << " " << globalCoords[i].first[j].first << std::endl;
  //    }
  //   }
  //std::cout << std::endl;
  std::vector<tupleRange> merged;
  try {
    // std::cout << "gonna merge ranges, which initially are: " << std::endl;
    // printRanges(recompRanges, globalCoords);
    merged = expandAndMergeRanges(navigator, recompRanges, index.k * index.j,
                                  data.blockExists);
    //std::cout << "merged: " << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
  }
  //std::cout << "\nNode " << node->identifier << " has " << merged.size()
            // << " merged ranges\n";
  int32_t range_i = -1;
  for (auto &recompRange : std::ranges::reverse_view(merged)) {
    range_i++;
    // std::cout << "Processing range " << range_i << " from ("
    //           << recompRange.start.blockId << ", " << recompRange.start.nucPos
    //           << ", " << recompRange.start.nucGapPos << ") to ("
    //           << recompRange.stop.blockId << ", " << recompRange.stop.nucPos
    //           << ", " << recompRange.stop.nucGapPos << ")\n";

    std::string recomputeSeq = tree::getNucleotideSequenceFromBlockCoordinates(
        recompRange.start, recompRange.stop, data.sequence, data.blockExists,
        data.blockStrand, T, node, globalCoords);
    // std::cout << "Recompute sequence: " << recomputeSeq << std::endl;
    tupleCoord_t lastDownstreamSeedPos = recompRange.stop;
    auto boundItr = seedMap.upper_bound(recompRange.stop);

    if (boundItr == seedMap.end()) {
      lastDownstreamSeedPos = tupleCoord_t{-1, -1, -1};
    } else {
      lastDownstreamSeedPos = boundItr->first;
    }
    // std::cout << "Last downstream seed pos: (" <<
    // lastDownstreamSeedPos.blockId
    //           << ", " << lastDownstreamSeedPos.nucPos << ", "
    //           << lastDownstreamSeedPos.nucGapPos << ")\n";
    int32_t str_i = recomputeSeq.size();
    for (tupleCoord_t currCoord = recompRange.stop;
         currCoord >= recompRange.start;
         currCoord = navigator.decrement(currCoord)) {
      str_i--;
      char nt = recomputeSeq[str_i];
      // std::cout << "Processing coord (" << currCoord.blockId << ", "
      //           << currCoord.nucPos << ", " << currCoord.nucGapPos
      //           << "): " << tupleToScalarCoord(currCoord, globalCoords)
      //           << " with nt " << nt << std::endl;
      if (str_i < 0) {
        break;
      }
      if (!data.blockExists[currCoord.blockId].first) {
        // std::cout << "smushing...\n";
        recomputeSeq[str_i] = '-';
      }
      // std::cout << "str[i=" << str_i << "] = " << recomputeSeq[str_i]
      //           << std::endl;
      if (!data.blockExists[currCoord.blockId].first) {
        if (seedMap.find(currCoord) != seedMap.end()) {
          // std::cout << "(+)->(-): " << seedMap[currCoord] << std::endl;
          // was a seed, no longer a seed due to block no exist -> delete
          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);
          std::string str = "";
          str += std::to_string(tupleToScalarCoord(currCoord, globalCoords));
          str += ":@";
          str += seedMap[currCoord];
          if (seedMap[currCoord].size() == index.j * index.k) {
            outDeletions.insert(str);
          }
        }
      } else if (recomputeSeq[str_i] == '-' ||
                 recomputeSeq[str_i] ==
                     'x') { // block does exist but seq is a gap
        if (seedMap.find(currCoord) != seedMap.end()) {
          // is a gap, no longer a seed -> delete
          // std::cout << "(+)->(_): " << seedMap[currCoord] << std::endl;

          std::string str = "";
          str += std::to_string(tupleToScalarCoord(currCoord, globalCoords));
          str += ":@";
          str += seedMap[currCoord];
          if (seedMap[currCoord].size() == index.j * index.k) {
            outDeletions.insert(str);
          }
          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);

        } /* else: no seed, wasn't seed, no change */
      } else {
        // std::cout << "(+)->(+): " << seedMap[currCoord] << std::endl;
        // block exists and seq is not a gap at currCoord
        // get the next k non-gap bases
        std::string kmer = "";
        int64_t seen_k = 0;
        int64_t k_pos = str_i;
        while (seen_k < index.k * index.j && k_pos < recomputeSeq.size()) {
          if (recomputeSeq[k_pos] != '-' && recomputeSeq[k_pos] != 'x') {
            kmer += recomputeSeq[k_pos];
            seen_k++;
          }
          k_pos++;
        }
        // std::cout << "kmer: " << kmer << std::endl;
        if (seedMap.find(currCoord) != seedMap.end()) {
          // non gap position and kmer is already a seed.
          std::string prevseedmer =
              lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1}
                  ? seedMap[lastDownstreamSeedPos]
                  : "";
          // std::cout << "prevseedmer: " << prevseedmer << std::endl;
          // std::cout << "bt " << c << " " << seedMap[c] << std::endl;

          if (seeding::is_syncmer(kmer, index.s, false)) {
            // Is it still a seed?
            // std::cout << "still seed" << std::endl;
            backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
            seedMap[currCoord] =
                kmer + prevseedmer.substr(0, (index.j - 1) * index.k);

            // std::cout << "set seedMap[" << c << "] = " << seedMap[c] <<
            // std::endl;
            lastDownstreamSeedPos = currCoord;
            if (seedMap[currCoord].size() == index.j * index.k) {
              std::string str = "";
              str = std::to_string(tupleToScalarCoord(currCoord, globalCoords));
              str += ":";
              str += seedMap[currCoord];
              outInsertions.insert(str);
            }
          } else {
            // no longer a seed -> delete
            // std::cout << "no longer seed, del " << c << " " << seedMap[c] <<
            // std::endl;
            std::string str = "";
            str += std::to_string(tupleToScalarCoord(currCoord, globalCoords));
            str += ":@";
            str += seedMap[currCoord];
            if (seedMap[currCoord].size() == index.j * index.k) {
              outDeletions.insert(str);
            }
            seedsToClear.push_back(currCoord);
          }
        } else {
          //  not in seed map, could be a seed now
          if (seeding::is_syncmer(kmer, index.s, false)) {
            backtrack.push_back(std::make_pair(currCoord, ""));
            std::string prevseedmer =
                lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1}
                    ? seedMap[lastDownstreamSeedPos]
                    : "";
            if (kmer.size() == index.j * index.k) {
              std::string str = "";
              str +=
                  std::to_string(tupleToScalarCoord(currCoord, globalCoords));
              str += ":";
              str += kmer;
              outInsertions.insert(str);
            }
            seedMap[currCoord] =
                kmer + prevseedmer.substr(0, (index.j - 1) * index.k);
            lastDownstreamSeedPos = currCoord;
          }
        }
      }
    }
  }
  for (const auto &pos : seedsToClear) {
    if (seedMap.find(pos) != seedMap.end()) {
      seedMap.erase(pos);
    }
  }
  index.outStream << node->identifier << " ";
  for (const std::string &s : outDeletions) {
    index.outStream << s << " ";
  }
  for (const std::string &s : outInsertions) {
    index.outStream << s << " ";
  }
  index.outStream << "\n";

  nodesDone++;
  /* Recursive step */
  for (const Node *child : node->children) {
    buildHelper(data, seedMap, index, T, child, globalCoords, altGlobalCoords,
                navigator);
  }
  // undo seed mutations
  for (const auto &back : backtrack) {
    if (seedMap.find(back.first) != seedMap.end()) {
      seedMap.erase(back.first);
    } else {
      seedMap[back.first] = back.second;
    }
  }
  /* Undo sequence mutations when backtracking */
  undoMutations(data, index, T, node, blockMutData, nucMutData);
}

/* Interface implementation */
void pmi::build(seedIndex &index, Tree *T, const size_t j, const size_t k,
                const size_t s) {
  /* Setup for seed indexing */
  tree::mutableTreeData data;
  tree::globalCoords_t globalCoords;
  std::vector<tupleCoord_t> altGlobalCoords;

  tree::setup(data, globalCoords, altGlobalCoords, T);

  index.j = j;
  index.k = k;
  index.s = s;

  seedMap_t seedMap;
  index.outStream << k << " " << s << " " << j << "\n";
  std::cout << "Building index\n";
  /* Recursive traversal of tree to build the index */
  CoordNavigator navigator(data.sequence);
  buildHelper(data, seedMap, index, T, T->root, globalCoords, altGlobalCoords,
              navigator);
}

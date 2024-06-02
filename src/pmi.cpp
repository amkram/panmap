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
                    const globalCoords_t &globalCoords, SeedmerIndex &index) {
  blockExists_t &blockExists = data.blockExists;
  blockStrand_t &blockStrand = data.blockStrand;
  
  // Block mutations
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
    recompRanges.push_back({tupleCoord_t{blockId, 0, data.sequence[blockId].first[0].second.empty() ? -1 : 0},
        tupleCoord_t{blockId, (int32_t)data.sequence[blockId].first.size() - 1, -1}});
  }

  // Nuc mutations
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

void undoMutations(mutableTreeData &data, SeedmerIndex &index, Tree *T,
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

// Go upstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandLeft(const CoordNavigator &navigator, tupleCoord_t coord,
                        int neededNongap, blockExists_t &blockExists) {
  

  int count = 0;
  //std::cout << "expandLeft..." << std::endl;
  // std::cout << "orig coord: (" << coord.blockId << ", " << coord.nucPos << ",
  // "
  //           << coord.nucGapPos << ") to => ... " << std::endl;

  // std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

  while (count < neededNongap && coord > tupleCoord_t{0, 0, 0}) {

    //std::cout << "count: " << count << " neededNongap: " << neededNongap << std::endl;
    // std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl; 

    if (!blockExists[coord.blockId].first) {
      //TODO jump down to prev block pleaesse

      //std::cout << "coord pre decrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl; 
      coord = navigator.decrement(coord);
      //std::cout << "coord post decrememnto " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

      continue;
    }
    // std::cout << "isGap? " << std::endl;
    if (!navigator.isGap(coord)) {
      //std::cout << "is not gap\n";
      count++;
    }

    //tupleCoord_t prev = coord;

    // std::cout << " from (" << coord.blockId << ", " << coord.nucPos << ", "
    // << coord.nucGapPos << ") to => ";
    //std::cout << "coord pre mad decrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl; 
    coord = navigator.decrement(coord);
    //std::cout << "coord post mado decrememnto " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

  }
  return coord;
}

// Go downstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandRight(const CoordNavigator &navigator, tupleCoord_t coord,
                         int neededNongap, blockExists_t &blockExists) {
  

  int count = 0;

  //std::cout << "ENTERING EXPANDRIGHT\n";

  //std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;


  while (count < neededNongap && coord < tupleCoord_t{-1, -1, -1}) {

    //std::cout << "count: " << count << " neededNongap: " << neededNongap << std::endl;
    //std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl; 


    if (!blockExists[coord.blockId].first) {
      //std::cout << "block doesn't exist\n";
      //TODO jump down to next block pleaesse

      //std::cout << "coord pre incrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl; 
      coord = navigator.increment(coord);
      // std::cout << "coord post instigation " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
      // std::cout << (coord < tupleCoord_t{-1, -1, -1}) << "\n";
      // std::cout << "we\n";
      continue;
    }

    if (!navigator.isGap(coord)) {
      //std::cout << "is not gap\n";
      count++;
    }
    //tupleCoord_t prev = coord;

    //std::cout << "coord pre incrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl; 
    coord = navigator.increment(coord);
    // std::cout << "coord post incrimination " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
  }
  // std::cout << "returning coord\n";
  return coord;
}

// Merges each range with overlapping ranges after expanding left and right
// by `neededNongap` non-gap nucleotides.
std::vector<tupleRange> expandAndMergeRanges(const CoordNavigator &navigator,
                                             std::vector<tupleRange> &ranges,
                                             int neededNongap,
                                             blockExists_t &blockExists) {
  

  if (ranges.empty())
    return {};

  std::vector<tupleRange> merged;
  
  tupleRange current = {
        expandLeft(navigator, ranges[0].start, neededNongap, blockExists),
        expandRight(navigator, ranges[0].stop, neededNongap, blockExists),
    };

  for (size_t i = 1; i < ranges.size(); ++i) {
    // std::cout << "Merging range " << i << " which is " <<
    // ranges[i].start.blockId << ", " << ranges[i].start.nucPos << ", " <<
    // ranges[i].start.nucGapPos << " to " << ranges[i].stop.blockId << ", " <<
    // ranges[i].stop.nucPos << ", " << ranges[i].stop.nucGapPos << std::endl;

    // std::cout << "unexp range:" << ranges[i].start.blockId << ", "
    //           << ranges[i].start.nucPos << ", " << ranges[i].start.nucGapPos
    //           << " to " << ranges[i].stop.blockId << ", "
    //           << ranges[i].stop.nucPos << ", " << ranges[i].stop.nucGapPos
    //           << std::endl;
      
    // std::cout << "WEEEEEEEEEEE\n";
    tupleRange expandedRange = {
        expandLeft(navigator, ranges[i].start, neededNongap, blockExists),
        expandRight(navigator, ranges[i].stop, neededNongap, blockExists),
    };
    // std::cout << "HELLO????????????????\n";

    // std::cout << "Expanded range:" << expandedRange.start.blockId << ", "
    //           << expandedRange.start.nucPos << ", " << expandedRange.start.nucGapPos
    //           << " to " << expandedRange.stop.blockId << ", "
    //           << expandedRange.stop.nucPos << ", " << expandedRange.stop.nucGapPos
    //           << std::endl;


    //if (expandedRange.start == tupleCoord_t{-1, -1, -1}) {
    //  expandedRange.start = tupleCoord_t{0, 0, 0};
    //}

    if (expandedRange.start <= current.stop) {
      //tmpStop = current.stop;
      current.stop = std::max(current.stop, expandedRange.stop);
    } else {
      merged.push_back(current);
      //tmpStop = current.stop;
      current = expandedRange;
    }
  }
  merged.push_back(current);

  return merged;
}

// Get a single integer representing a position in the MSA from a tupleCoord_t = {blockId, nucPos, nucGapPos}
int64_t tupleToScalarCoord(const tupleCoord_t &coord,
                           const globalCoords_t &globalCoords) {
  if (coord == tupleCoord_t{-1, -1, -1}) {
    return globalCoords.back().first.back().first;
  }
  if (coord.nucGapPos >= 0) {
    return globalCoords[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
  }
  return globalCoords[coord.blockId].first[coord.nucPos].first;
}

// Recursive function to build the seed index
void buildHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index,
                 Tree *T, const Node *node, const globalCoords_t &globalCoords,
                 CoordNavigator &navigator) {
  blockMutData_t blockMutData;
  nucMutData_t nucMutData;

  // First, a range is made marking the start -> end
  // of each block and nuc mutation. This is done while
  // applying mutations to the sequence object.
  std::vector<tupleRange> recompRanges;
  applyMutations(data, seedMap, blockMutData, recompRanges, nucMutData, T, node,
                 globalCoords, index);
  std::sort(recompRanges.begin(), recompRanges.end());
  
  
  std::vector<tupleCoord_t> seedsToClear; // seeds to clear from seedMap
  std::vector<std::pair<tupleCoord_t, std::string>> backtrack;

  std::vector<tupleRange> merged = expandAndMergeRanges(navigator, recompRanges, index.k(), data.blockExists);
  std::cout << "merged ranges: " << std::endl;
  for (auto &range : merged) {
    std::cout << "range: " << range.start.blockId << ", " << range.start.nucPos << ", " << range.start.nucGapPos << " to " << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;
    std::cout << "ntpos: " << tupleToScalarCoord(range.start, globalCoords) << " to " << tupleToScalarCoord(range.stop, globalCoords) << std::endl;
  }
  // Protobuf message for this node's mutations
  NodeSeedmerMutations *pb_node_mutations = index.add_per_node_mutations();
  pb_node_mutations->set_node_id(node->identifier);

  // Seed re-processing
  for (auto &range : std::ranges::reverse_view(merged)) {
    std::string recomputeSeq = tree::getNucleotideSequenceFromBlockCoordinates(range.start, range.stop,
      data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords);
    
    // Track the last downstream seed to stack k-mers into seedmers
    tupleCoord_t lastDownstreamSeedPos = range.stop;
    auto boundItr = seedMap.upper_bound(range.stop);
    if (boundItr == seedMap.end()) {
      lastDownstreamSeedPos = tupleCoord_t{-1, -1, -1};
    } else {
      lastDownstreamSeedPos = boundItr->first;
    }

    int32_t str_i = recomputeSeq.size();
    int32_t seen_non_gap = 0;
    for (auto currCoord = range.stop; currCoord >= range.start; currCoord = navigator.decrement(currCoord)) {
      str_i--;
      char nt = recomputeSeq[str_i];
      std::cout << "Processing coord (" << currCoord.blockId << ", "
                << currCoord.nucPos << ", " << currCoord.nucGapPos
                << "): " << tupleToScalarCoord(currCoord, globalCoords)
                << " with nt " << nt << std::endl;
      if (str_i < 0) {
        break;
      }
      // todo jump around
      if (!data.blockExists[currCoord.blockId].first) {
        std::cout << "block doesn't exist at (" << currCoord.blockId << ", "
                  << currCoord.nucPos << ", " << currCoord.nucGapPos
                  << "): " << tupleToScalarCoord(currCoord, globalCoords)
                  << std::endl;
        // std::cout << "smushing...\n";
        recomputeSeq[str_i] = '-';
      }

      if (recomputeSeq[str_i] != '-' && recomputeSeq[str_i] != 'x') {
        seen_non_gap++;
      }

      
      // std::cout << "str[i=" << str_i << "] = " << recomputeSeq[str_i]
      //           << std::endl;

      if (!data.blockExists[currCoord.blockId].first) {
        if (seedMap.find(currCoord) != seedMap.end()) {
          // std::cout << "(+)->(-): " << seedMap[currCoord] << std::endl;
          // was a seed, no longer a seed due to block no exist -> delete
          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);
          if (seedMap[currCoord].size() == index.j() * index.k()) {
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
            pb_mut->set_seq(seedMap[currCoord]);
          }
        }
      } else if (recomputeSeq[str_i] == '-' ||
                 recomputeSeq[str_i] ==
                     'x') { // block does exist but seq is a gap
        if (seedMap.find(currCoord) != seedMap.end()) {
          // is a gap, no longer a seed -> delete
          // std::cout << "(+)->(_): " << seedMap[currCoord] << std::endl;
          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
            pb_mut->set_seq(seedMap[currCoord]);
            
          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);

        } /* else: no seed, wasn't seed, no change */
      } else {
        std::cout << "(+)->(+): " << seedMap[currCoord] << std::endl;
        // block exists and seq is not a gap at currCoord
        // get the next k non-gap bases
        std::string kmer = "";
        int64_t seen_k = 0;
        int64_t k_pos = str_i;
        while (seen_k < index.k() * index.j() && k_pos < recomputeSeq.size()) {
          if (recomputeSeq[k_pos] != '-' && recomputeSeq[k_pos] != 'x') {
            kmer += recomputeSeq[k_pos];
            seen_k++;
          }
          k_pos++;
        }

        std::cout << "kmer: " << kmer << std::endl;
        if (seedMap.find(currCoord) != seedMap.end()) {
          // non gap position and kmer is already a seed.
          std::string prevseedmer =
              lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1}
                  ? seedMap[lastDownstreamSeedPos]
                  : "";
          // std::cout << "prevseedmer: " << prevseedmer << std::endl;
          // std::cout << "bt " << c << " " << seedMap[c] << std::endl;

          if (seeding::is_syncmer(kmer, index.s(), false)) {
            // Is it still a seed?
            // std::cout << "still seed" << std::endl;
            backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
            seedMap[currCoord] =
                kmer + prevseedmer.substr(0, (index.j() - 1) * index.k());

            // std::cout << "set seedMap[" << c << "] = " << seedMap[c] <<
            // std::endl;
            lastDownstreamSeedPos = currCoord;
            if (seedMap[currCoord].size() == index.j() * index.k()) {
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
              pb_mut->set_seq(seedMap[currCoord]);
            }
          } else {
            // no longer a seed -> delete
            // std::cout << "no longer seed, del " << c << " " << seedMap[c] <<
            // std::endl;
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
            pb_mut->set_seq(seedMap[currCoord]);
      
            seedsToClear.push_back(currCoord);
          }
        } else {
          //  not in seed map, could be a seed now
          if (seeding::is_syncmer(kmer, index.s(), false)) {
            backtrack.push_back(std::make_pair(currCoord, ""));
            std::string prevseedmer =
                lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1}
                    ? seedMap[lastDownstreamSeedPos]
                    : "";
            seedMap[currCoord] = kmer + prevseedmer.substr(0, (index.j() - 1) * index.k());
            if (seedMap[currCoord].size() == index.j() * index.k()) {
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
              pb_mut->set_seq(seedMap[currCoord]);
            }
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

  /* Recursive step */
  for (const Node *child : node->children) {
    //exit(0);
    buildHelper(data, seedMap, index, T, child, globalCoords, navigator);
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

/* implementation */
void pmi::build(SeedmerIndex &index, Tree *T, int j, int k, int s) {

  // Setup for seed indexing
  tree::mutableTreeData data;
  tree::globalCoords_t globalCoords;
  
  tree::setup(data, globalCoords, T);

  index.set_j(j);
  index.set_k(k);
  index.set_s(s);

  std::cout << "Building index\n";

  // Stores seed(mer)s at positions where one exists (in tuple global coords).
  // At each node, seedMap is updated to contain seeds present in the node.
  // Each element is {tupleCoord_t, string} = {<pos>, <seedmer>}
  seedMap_t seedMap;

  // For navigating the tuple coordinate space of the sequence object.
  // Each coordinate is {blockId, nucPos, nucGapPos} where each nucPos
  // indicates a "main nuc" position within the block, and...
  //  if (nucGapPos == -1): => the coordinate reflects sequence[blockId].first[nucPos].first,
  //       which is the main nucleotide char at nucPos. There is no associated gap list (see below.)
  //  Otherwise: => the coordinate reflects sequence[blockId].first[nucPos].second[nucGapPos],
  //       which is the nucleotide char at nucGapPos within the optional "gap list" at nucPos.
  //       The gap list is a list of nucleotides that come before the main nucleotide
  //       (in tuple/int global coords) at nucPos. An example ascending ordered list of 
  //       tuple coordinates is: {0,2,0}, {0,2,1}, {0,2,2}, {0,2,-1}, {0,3,-1}, {1,0,0}, ...
  //       
  CoordNavigator navigator(data.sequence);
  
  /* Recursive traversal of tree to build the index */
  buildHelper(data, seedMap, index, T, T->root, globalCoords,navigator);
}

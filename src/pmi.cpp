#include "pmi.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include "seed_annotated_tree.hpp"
#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <memory>
#include <variant>
#include <capnp/serialize.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include "index.capnp.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>


const bool debug = false;

enum Step {
  BUILD,
  PLACE,
  SPECTRUM
};

void flipCoords(int32_t blockId, globalCoords_t &globalCoords) {

  int64_t start;
  if (globalCoords[blockId].first[0].second.size() == 0) {
    start = globalCoords[blockId].first[0].first;
  } else {
    start = globalCoords[blockId].first[0].second[0];
  }
  int64_t stop = globalCoords[blockId].first.back().first;
  if (start < stop) {
    for (int32_t i = globalCoords[blockId].first.size() - 1; i >= 0; i--) {
      globalCoords[blockId].first[i].first = start;
      start++;
      for (int32_t j = globalCoords[blockId].first[i].second.size() - 1; j >= 0; j--) {
        globalCoords[blockId].first[i].second[j] = start;
        start++;
      }
    }
  } else {
    for (int32_t i = 0; i < globalCoords[blockId].first.size(); i++) {
      for (int32_t j = 0; j < globalCoords[blockId].first[i].second.size(); j++) {
        globalCoords[blockId].first[i].second[j] = stop;
        stop++;
      }
      globalCoords[blockId].first[i].first = stop;
      stop++;

    }
  } 
}

void applyMutations(mutableTreeData &data,
                    blockMutationInfo_t &blockMutationInfo,
                    std::vector<tupleRange> &recompRanges,
                    mutationInfo_t &mutationInfo, Tree *T, Node *node,
                    globalCoords_t &globalCoords,
                    CoordNavigator &navigator, const std::vector<std::pair<int64_t, int64_t>> &blockRanges,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks,
                    blockExists_t& oldBlockExists, blockStrand_t& oldBlockStrand, const bool isPlacement,
                    std::unordered_set<int64_t>& inverseBlockIds, std::vector<std::pair<bool, int64_t>>& inverseBlockIdsBacktrack)
{

  blockExists_t &blockExists = data.blockExists;
  blockStrand_t &blockStrand = data.blockStrand;


  oldBlockExists = blockExists;
  oldBlockStrand = blockStrand;


  // here
  auto &blocks = T->blocks;
  auto &gaps = T->gaps;
  auto &blockGaps = T->blockGaps;
  auto &sequenceInverted = T->sequenceInverted;
  auto &sequence = data.sequence;
  auto &circularSequences = T->circularSequences;
  // List of blocks. Each block has a nucleotide list. Along with each
  // nucleotide is a gap list.

  // Block Mutations
  for (auto mutation : node->blockMutation)
  {

    int32_t primaryBlockId = mutation.primaryBlockId;
    int32_t secondaryBlockId = mutation.secondaryBlockId;
    bool type = mutation.blockMutInfo;
    bool inversion = mutation.inversion;

    bool startingStrand = blockStrand[primaryBlockId].first;

    if (type == 1)
    {
      // insertion
      bool oldStrand;
      bool oldMut;
      if (secondaryBlockId != -1)
      {
        oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
        oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
        blockExists[primaryBlockId].second[secondaryBlockId] = true;

        // if insertion of inverted block takes place, the strand is backwards
        blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
      }
      else
      {
        oldStrand = blockStrand[primaryBlockId].first;
        oldMut = blockExists[primaryBlockId].first;
        blockExists[primaryBlockId].first = true;

        // if insertion of inverted block takes place, the strand is backwards
        blockStrand[primaryBlockId].first = !inversion;
      }
      blockMutationInfo.emplace_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, true, !inversion));
    }
    else
    {
      bool oldMut;
      bool oldStrand;
      if (inversion)
      {
        // This means that this is not a deletion, but instead an inversion
        if (secondaryBlockId != -1)
        {
          oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
          oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
          blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
        }
        else
        {
          oldStrand = blockStrand[primaryBlockId].first;
          oldMut = blockExists[primaryBlockId].first;
          blockStrand[primaryBlockId].first = !oldStrand;
        }
        if (oldMut != true) 
        {
          std::cout << "There was a problem in PanMAT generation. Please Report." << std::endl;
        }
        blockMutationInfo.emplace_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, oldMut, !oldStrand));
      }
      else
      {
        // Actually a deletion

        if (secondaryBlockId != -1)
        {
          oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
          oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
          blockExists[primaryBlockId].second[secondaryBlockId] = false;

          // resetting strand to true during deletion
          blockStrand[primaryBlockId].second[secondaryBlockId] = true;
        }
        else
        {
          oldStrand = blockStrand[primaryBlockId].first;
          oldMut = blockExists[primaryBlockId].first;
          blockExists[primaryBlockId].first = false;

          // resetting strand to true during deletion
          blockStrand[primaryBlockId].first = true;
        }
      }
      blockMutationInfo.emplace_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true));
    }

    if (!isPlacement) {

    //Push new recomb range, order depends on inversion
    tupleRange newRange = {tupleCoord_t{primaryBlockId, 0, data.sequence[primaryBlockId].first[0].second.empty() ? -1 : 0},
                            tupleCoord_t{primaryBlockId, (int32_t)data.sequence[primaryBlockId].first.size() - 1, -1}};
    
    if(! blockStrand[primaryBlockId].first){
      auto temp = newRange.start;
      newRange.start = newRange.stop;
      newRange.stop = temp;
    }

    recompRanges.emplace_back(newRange);
    }

    if (startingStrand != blockStrand[primaryBlockId].first) {
      if (blockStrand[primaryBlockId].first) {
        inverseBlockIds.erase(primaryBlockId);
        inverseBlockIdsBacktrack.emplace_back(std::make_pair(false, primaryBlockId));
      } else {
        inverseBlockIds.insert(primaryBlockId);
        inverseBlockIdsBacktrack.emplace_back(std::make_pair(true, primaryBlockId));
      }
      flipCoords(primaryBlockId, globalCoords);
    }
  }

  // Nuc mutations
  for (int64_t i = 0; i < node->nucMutation.size(); i++)
  {
    int32_t primaryBlockId = node->nucMutation[i].primaryBlockId;
    int32_t secondaryBlockId = node->nucMutation[i].secondaryBlockId;

    int32_t nucPosition = node->nucMutation[i].nucPosition;
    int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
    uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
    char newVal = '-';
    if (type < 3)
    {
      // Either S, I or D
      int len = ((node->nucMutation[i].mutInfo) >> 4);

      if (!isPlacement){
        tupleRange newRange = {tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition + len}};
      

        if (nucGapPosition == -1){
          //recompRanges.emplace_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition+len, -1}});
          newRange = {tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition+len, -1}};
        }

        if(! blockStrand[primaryBlockId].first){
          auto temp = newRange.start;
          newRange.start = newRange.stop;
          newRange.stop = temp;
        }

        recompRanges.emplace_back(newRange);
      }



      if (type == panmanUtils::NucMutationType::NS)
      {
        // Substitution

        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            }
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j];
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            }
          }
        }
      }
      else if (type == panmanUtils::NucMutationType::NI)
      {
        // Insertion
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            }
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j];
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            }
          }
        }
      }
      else if (type == panmanUtils::NucMutationType::ND)
      {
        // Deletion
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j] = '-';
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, '-'));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
            }
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j];
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = '-';
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, '-'));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              sequence[primaryBlockId].first[nucPosition + j].first = '-';
              mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
            }
          }
        }
      }
    }
    else
    {
      int len = 0;

      if (!isPlacement){
        tupleRange newRange ={tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},
                                tupleCoord_t{primaryBlockId, nucPosition + len, nucGapPosition}};
    
        if(! blockStrand[primaryBlockId].first){
          auto temp = newRange.start;
          newRange.start = newRange.stop;
          newRange.stop = temp;
        }

        recompRanges.emplace_back(newRange);
      }
      if (type == panmanUtils::NucMutationType::NSNPS)
      {
        // SNP Substitution
        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].first;
            sequence[primaryBlockId].first[nucPosition].first = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
      }
      else if (type == panmanUtils::NucMutationType::NSNPI)
      {
        // SNP Insertion
        newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].first;
            sequence[primaryBlockId].first[nucPosition].first = newVal;
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
      }
      else if (type == panmanUtils::NucMutationType::NSNPD)
      {
        // SNP Deletion
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].first;
            sequence[primaryBlockId].first[nucPosition].first = '-';
            mutationInfo.emplace_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
        }
      }
    }
  }

  if (!isPlacement){
    for (auto &mutation : mutationInfo) {
      int blockId = std::get<0>(mutation);
      if (!(oldBlockExists[blockId].first && blockExists[blockId].first)) {
        continue;
      }

      int nucPos = std::get<2>(mutation);
      int nucGapPos = std::get<3>(mutation);
      char parChar = std::get<4>(mutation) == 'x' ? '-' : std::get<4>(mutation);
      char curChar = std::get<5>(mutation) == 'x' ? '-' : std::get<5>(mutation);
      int64_t scalar = tupleToScalarCoord(tupleCoord_t{blockId, nucPos, nucGapPos}, globalCoords);
      if (!data.blockStrand[blockId].first) {
        scalar = blockRanges[blockId].first + blockRanges[blockId].second - scalar;
      }
    
      if (parChar != '-' && curChar == '-') {
        // nuc to gap
        if (!gapRunUpdates.empty() && gapRunUpdates.back().first == true && gapRunUpdates.back().second.second + 1 == scalar) {
          ++(gapRunUpdates.back().second.second);
        }
        else {
          gapRunUpdates.emplace_back(true, std::make_pair(scalar, scalar)); 
        }
      } else if (parChar == '-' && curChar != '-') {
        // gap to nuc
        if (!gapRunUpdates.empty() && gapRunUpdates.back().first == false && gapRunUpdates.back().second.second + 1 == scalar) {
          ++(gapRunUpdates.back().second.second);
        } else {
          gapRunUpdates.emplace_back(false, std::make_pair(scalar, scalar));
        }
      }
    }
  }
}

void undoMutations(mutableTreeData &data, Tree *T,
                   const Node *node, const blockMutationInfo_t &blockMutationInfo,
                   const mutationInfo_t &mutationInfo, globalCoords_t &globalCoords)
{
  auto &sequence = data.sequence;
  auto &blockExists = data.blockExists;
  auto &blockStrand = data.blockStrand;
 // Undo block mutations when current node and its subtree have been processed
    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);

            if(blockStrand[std::get<0>(mutation)].second[std::get<1>(mutation)] != std::get<3>(mutation)){
              blockStrand[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<3>(mutation);
              flipCoords(std::get<0>(mutation),globalCoords);
            }

        } else {

            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);

            if(blockStrand[std::get<0>(mutation)].first != std::get<3>(mutation)){
              blockStrand[std::get<0>(mutation)].first = std::get<3>(mutation);
              flipCoords(std::get<0>(mutation),globalCoords);
            }
        }
    }

    // Undo nuc mutations when current node and its subtree have been processed
    for(auto it = mutationInfo.rbegin(); it != mutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            if(std::get<3>(mutation) != -1) {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].second[std::get<1>(mutation)][std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        } else {
            if(std::get<3>(mutation) != -1) {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].second[std::get<3>(mutation)] = std::get<4>(mutation);
            } else {
                sequence[std::get<0>(mutation)].first[std::get<2>(mutation)].first = std::get<4>(mutation);
            }
        }
    }


}

// Go upstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandLeft(CoordNavigator &navigator, tupleCoord_t coord,
                        int neededNongap, blockExists_t &blockExists,
                                             blockStrand_t &blockStrand, tupleCoord_t stop_coord={-1,-1,-1})
{
  int count = -1;
  
  //coord
  auto start = tupleCoord_t{0,0,0};
  if(navigator.sequence[0].first[0].second.empty()) {
    start.nucGapPos = -1;
  }

  while (count < neededNongap && coord > start)
  {

    if(stop_coord.blockId != -1 && (coord.blockId < stop_coord.blockId || coord == stop_coord)){
      break;
    }

    if (!blockExists[coord.blockId].first)
    {
      //Jump down to previous block
      if(coord.blockId == 0){

        return start;
      }else{
        
        if(blockStrand[coord.blockId - 1].first){
          //not inverted, jump to top of next block
          coord = tupleCoord_t{coord.blockId - 1, (int64_t)navigator.sequence[coord.blockId - 1].first.size() - 1, -1};
        }else{
          //inverted, jump to bottom of next block
          coord.blockId -= 1;
          coord.nucPos = 0;
          coord.nucGapPos = 0;
          if(navigator.sequence[coord.blockId].first[0].second.empty()) {
            coord.nucGapPos = -1;
          }
        }
        
      }
      
      continue;
    }
    
    if (!navigator.isGap(coord))
    {
      count++;
    }

    coord = navigator.newdecrement(coord, blockStrand);
  }

  if(coord.blockId != -1 && !blockExists[coord.blockId].first){
    coord = navigator.newincrement(coord, blockStrand);
  }

  return coord;
}

// Go downstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandRight(CoordNavigator &navigator, tupleCoord_t &coord,
                         int neededNongap, blockExists_t &blockExists,
                                             blockStrand_t &blockStrand)
{

  int count = -1;
  
  while (count < neededNongap && coord < tupleCoord_t{-1, -1, -1})
  {

    if (!blockExists[coord.blockId].first)
    {

      //Jump to next block
      if(coord.blockId == navigator.sequence.size() - 1){
        return tupleCoord_t{-1,-1,-1};
      }else{

        if( ! blockStrand[coord.blockId + 1].first){
          //inverted, jump to top of next block
          coord = tupleCoord_t{coord.blockId + 1, (int64_t)navigator.sequence[coord.blockId + 1].first.size() - 1, -1};
        }else{
          //not inverted, jump to bottom of next block
          coord.blockId += 1;
          coord.nucPos = 0;
          coord.nucGapPos = 0;
          if(navigator.sequence[coord.blockId].first[0].second.empty()) {
            coord.nucGapPos = -1;
          }
        }

      }

      continue;
    }

    if (!navigator.isGap(coord))
    {
      count++;
    }

    coord = navigator.newincrement(coord, blockStrand);
  }

  if(coord.blockId != -1 && !blockExists[coord.blockId].first){
    coord = navigator.newdecrement(coord, blockStrand);
  }

  return coord;
}


void updateGapMapStep(std::map<int64_t, int64_t>& gapMap, const std::pair<bool, std::pair<int64_t, int64_t>>& update, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates, bool recordGapMapUpdates) {
  bool toGap = update.first;
  int64_t start = update.second.first;
  int64_t end = update.second.second;

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

void updateGapMap(std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& updates, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates) {
  for (const auto& update : updates) {
    updateGapMapStep(gapMap, update, backtrack, gapMapUpdates);
  }
}

std::vector<std::pair<int64_t, int64_t>> invertRanges(const std::vector<std::pair<int64_t, int64_t>>& nucRanges, const std::pair<int64_t, int64_t>& invertRange) {
  std::vector<std::pair<int64_t, int64_t>> invertedRanges;

  auto [start, end] = invertRange;

  for (auto it = nucRanges.rbegin(); it != nucRanges.rend(); ++it) {
    const auto& [curStart, curEnd] = *it;
    invertedRanges.emplace_back(start + end - curEnd, start + end - curStart);
  }

  return invertedRanges;
}

void invertGapMap(std::map<int64_t, int64_t>& gapMap, const std::pair<int64_t, int64_t>& invertRange, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates) {
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
    updateGapMapStep(gapMap, {it->first, {curBeg, curEnd}}, backtrack, gapMapUpdates, false);
    curBeg = curEnd + 1;
  }

}

void makeCoordIndex(std::map<int64_t, int64_t>& coordIndex, const std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<int64_t, int64_t>>& blockRanges) {
  int64_t totalGapSize = 0;
  if (gapMap.empty() || gapMap.begin()->first > 0) {
    coordIndex[0] == totalGapSize;
  }
  for (auto &gap : gapMap) {
    int64_t gapStart = gap.first;
    int64_t gapEnd = gap.second;
    int64_t gapSize = gapEnd - gapStart + 1;
    if (gapEnd == blockRanges.back().second) break;
    totalGapSize += gapSize;
    coordIndex[gapEnd+1] = totalGapSize;
  }
}

int64_t degapGlobal(const int64_t& globalCoord, const std::map<int64_t, int64_t>& coordsIndex) {
    auto coordIt = coordsIndex.upper_bound(globalCoord);
    if (coordIt == coordsIndex.begin()) {
        return 0;
    }
    return globalCoord - std::prev(coordIt)->second;
}

// Merges each range with overlapping ranges after expanding left and right
// by `neededNongap` non-gap nucleotides.
std::vector<tupleRange> expandAndMergeRanges(CoordNavigator &navigator,
                                             std::vector<tupleRange> &ranges,
                                             int neededNongap,
                                             blockExists_t &blockExists,
                                             blockStrand_t &blockStrand, const globalCoords_t &globalCoords)
{

  if (ranges.empty())
    return {};


  std::vector<tupleRange> merged;

  tupleRange current = ranges[0];

  for (int64_t i = 1; i < ranges.size(); ++i) {
    
    bool replace = ranges[i].start <= current.stop;

    bool flip = ranges[i].start.blockId == current.stop.blockId && ! blockStrand[current.stop.blockId].first;
    if(flip){
      replace = ranges[i].start >= current.stop;
    }

    if (replace) //Merge ranges
    {
      flip = ranges[i].stop.blockId == current.stop.blockId && ! blockStrand[current.stop.blockId].first;

      if(flip){
        current.stop = std::min(current.stop, ranges[i].stop);
      }else{
        current.stop = std::max(current.stop, ranges[i].stop);
      }
    }
    else
    {
      merged.emplace_back(current);
      current = ranges[i];
    }
  }
  merged.emplace_back(current);

  ranges = merged;

  merged.clear();


  current = {
      expandLeft(navigator, ranges[0].start, neededNongap, blockExists, blockStrand),
      expandRight(navigator, ranges[0].stop, neededNongap, blockExists, blockStrand),
  };


  for (int64_t i = 1; i < ranges.size(); ++i) {

    tupleRange toPlace = {
        current.stop,
        current.stop,
    };

    //Search for first range that is beyond current range
    auto it = std::lower_bound(ranges.begin(), ranges.end(), toPlace, [&blockStrand](const tupleRange& A, const tupleRange& B) {
      if (A.start.blockId == B.start.blockId && !blockStrand[A.start.blockId].first) {
          return B < A; // Use B < A if blocks are inverted 
      }
      return A < B; // Default comparison
    });


    // Get the index from the range vector
    int index = std::distance(ranges.begin(), it);
    
    if(i < index - 1)
      i = index - 1;

    
    tupleRange expandedRange = {
        expandLeft(navigator, ranges[i].start, neededNongap, blockExists, blockStrand, current.stop),
        expandRight(navigator, ranges[i].stop, neededNongap, blockExists, blockStrand),
    };

    
    bool replace = expandedRange.start <= current.stop;

    //accounting for inversions TODO simplify
    bool flip = expandedRange.start.blockId == current.stop.blockId && ! blockStrand[current.stop.blockId].first;
    if(flip){
      replace = expandedRange.start >= current.stop;
    }

    if (replace) //Merge ranges
    {
      flip = expandedRange.stop.blockId == current.stop.blockId && ! blockStrand[current.stop.blockId].first;

      if(flip){
        current.stop = std::min(current.stop, expandedRange.stop);
      }else{
        current.stop = std::max(current.stop, expandedRange.stop);
      }
    }
    else
    {
      merged.emplace_back(current);
      current = expandedRange;
    }
  }
  merged.emplace_back(current);

  return merged;
}

// Get a single integer representing a position in the MSA from a tupleCoord_t = {blockId, nucPos, nucGapPos}
int64_t tupleToScalarCoord(const tupleCoord_t &coord,
                           const globalCoords_t &globalCoords)
{

  if (coord == tupleCoord_t{-1, -1, -1}) {
    return globalCoords.back().first.back().first;
  }
  if (coord.nucGapPos >= 0) {
    return globalCoords[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
  }
  
  return globalCoords[coord.blockId].first[coord.nucPos].first;
}

inline uint16_t setBit(uint16_t number, uint16_t n) {
    return number | ((uint16_t)1 << n);
}

// std::vector<std::tuple<std::string, int, int>>
// extractSeedmers(const std::string &seq, const int k, const int s,
//                 const bool open) {
//   std::vector<std::tuple<std::string, int, int>> syncmers;
//   std::unordered_map<int32_t, int32_t> degap;
//   int64_t pos = 0;
//   std::string ungapped = "";
//   for (int64_t i = 0; i < seq.size(); i++) {
//     char c = seq[i];
//     degap[pos] = i;
//     if (c != '-' && c != 'x') {
//       ungapped += c;
//       pos++;
//     }
//   }
//   if (ungapped.size() < k + 1) {
//     return syncmers;
//   }
//   for (int64_t i = 0; i < ungapped.size() - k + 1; ++i) {
//     std::string kmer = ungapped.substr(i, k);
//     if (seeding::is_syncmer(kmer, s, open)) {
//       syncmers.emplace_back(std::make_tuple(kmer, degap[i], degap[i + k - 1]));
//     }
//   }

//   return syncmers;
// }

//                     kmer,        hash  , reverse, start
std::vector<std::tuple<std::string, size_t, bool, int>>
extractSeedmers(const std::string &seq, const int k, const int s, const bool open) {
  std::vector<std::tuple<std::string, size_t, bool, int>> seedmers;
  for (int64_t i = 0; i < seq.size() - k + 1; ++i) {
    std::string kmer = seq.substr(i, k);
    auto [hash, isReverse, isSyncmer] = seeding::is_syncmer_rollingHash(kmer, s, true, 1);
    if (isSyncmer) {
      seedmers.emplace_back(std::make_tuple(kmer, hash, isReverse, i));
    }
  }
  return seedmers;
}

void bruteForceCoordIndex(const std::string& gappedSeq, std::map<int32_t, int32_t>& coordIndex){
  int32_t numGaps = 0;
  int32_t numNucs = 0;

  for (int32_t i = 0; i < gappedSeq.size(); i++){
    const char &c = gappedSeq[i];
    const char &p = gappedSeq[std::max(0, i-1)];
    if (c == '-' || c == 'x') {
      ++numGaps;
    } else {
      ++numNucs;
      if (i == 0 || p == '-' || p == 'x') {
        coordIndex[i] = numGaps;
      }
    }
  }
}

bool gappity = true;
// Recursive function to build the seed index
template <typename SeedMutationsType, typename GapMutationsType>
void buildOrPlace(Step method, mutableTreeData& data, std::vector<std::optional<std::string>>& onSeeds, std::vector<std::optional<std::pair<size_t, bool>>>& onSeedsHash, SeedMutationsType& perNodeSeedMutations_Index, GapMutationsType& perNodeGapMutations_Index, int seedK, int seedS, Tree* T, Node* node, globalCoords_t& globalCoords, CoordNavigator& navigator, std::vector<int64_t>& scalarCoordToBlockId, std::vector<std::unordered_set<int>>& BlocksToSeeds, std::vector<int>& BlockSizes, std::vector<std::pair<int64_t, int64_t>>& blockRanges, int64_t& dfsIndex, std::map<int64_t, int64_t>& gapMap, std::unordered_set<int64_t>& inverseBlockIds, int64_t jacNumer, int64_t jacDenom,  std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts) {
  // Variables needed for both build and place
  // std::vector<std::tuple<int64_t, bool, bool, std::optional<std::string>, std::optional<std::string>>> seedChanges;
  std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>>> seedChanges;
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBlocksBacktracks;
  std::vector<std::pair<bool, int64_t>> inverseBlockIdsBacktrack;
  std::map<int64_t, int64_t> coordIndex;
  std::vector<tupleRange> recompRanges;
  blockExists_t oldBlockExists = data.blockExists;
  blockStrand_t oldBlockStrand = data.blockStrand;

  applyMutations(data, blockMutationInfo, recompRanges,  mutationInfo, T, node, globalCoords, navigator, blockRanges, gapRunUpdates, gapRunBacktracks, oldBlockExists, oldBlockStrand, method == Step::PLACE, inverseBlockIds, inverseBlockIdsBacktrack);

  if (method == Step::BUILD) {
    
    //                    erase           beg      end
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapMapUpdates;
    if (gappity) {
      updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);
      std::vector<int64_t> invertedBlocks;
      for (int i = 0; i < data.blockExists.size(); i++) {
        const bool& oldExists = oldBlockExists[i].first;
        const bool& newExists = data.blockExists[i].first;
        if (newExists && !data.blockStrand[i].first) {
          invertedBlocks.emplace_back(i);
        }
        if (oldExists && !newExists) {
          // on to off -> block range to all gaps
          const auto& [start, end] = blockRanges[i];
          updateGapMapStep(gapMap, {true, {start, end}}, gapRunBacktracks, gapMapUpdates);
        } else if (!oldExists && newExists) {
          // off to on -> recompute across entire block
          tupleCoord_t coord = globalCoords[i].first[0].second.empty() ? tupleCoord_t{i, 0, -1} : tupleCoord_t{i, 0, 0};
          tupleCoord_t end = tupleCoord_t{i, globalCoords[i].first.size() - 1, -1};
          if (!data.blockStrand[i].first) std::swap(coord, end);
          // std::cout << "Recomputing block " << i << " from " << tupleToScalarCoord(coord, globalCoords) << " to " << tupleToScalarCoord(end, globalCoords) << std::endl;

          auto curIt = gapMap.end();
          std::pair<int64_t, int64_t> curNucRange = {-1, -1};
          std::vector<std::pair<int64_t, int64_t>> nucRanges;
          while (true) {
            char c = coord.nucGapPos == -1 ? data.sequence[coord.blockId].first[coord.nucPos].first : data.sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
            c = c == 'x' ? '-' : c;
            int64_t scalar = tupleToScalarCoord(coord, globalCoords);
            // std::cout << scalar << " " << c << std::endl;
            if (c != '-') {
              if (curNucRange.first != -1 && curNucRange.second + 1 == scalar) {
                ++curNucRange.second;
              } else {
                if (curNucRange.first != -1) {
                  nucRanges.emplace_back(curNucRange);
                }
                curNucRange = {scalar, scalar};
              }
            }

            if (coord == end) break;
            coord = navigator.newincrement(coord, data.blockStrand);
          }
          if (curNucRange.first != -1) {
            nucRanges.emplace_back(curNucRange);
          }

          if (data.blockStrand[i].first) {
            for (const auto& range : nucRanges) {
              updateGapMapStep(gapMap, {false, range}, gapRunBacktracks, gapMapUpdates);
            }
          } else {
            std::vector<std::pair<int64_t, int64_t>> invertedRanges = invertRanges(nucRanges, blockRanges[i]);
            for (const auto& range : invertedRanges) {
              updateGapMapStep(gapMap, {false, range}, gapRunBacktracks, gapMapUpdates);
            }
          }
        }
      }


      for (const auto& blockId : invertedBlocks) {
        invertGapMap(gapMap, blockRanges[blockId], gapRunBlocksBacktracks, gapMapUpdates);
      }

      makeCoordIndex(coordIndex, gapMap, blockRanges);
      
      for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
        const auto& [del, range] = *it;
        if (del) {
          gapMap.erase(range.first);
        } else {
          gapMap[range.first] = range.second;
        }
      }

    }

    std::sort(recompRanges.begin(), recompRanges.end(), [&data](const tupleRange& A, const tupleRange& B) {
      if (A.start.blockId == B.start.blockId && !data.blockStrand[A.start.blockId].first) {
          return B < A; // Use B < A if blocks are inverted 
      }
      return A < B; // Default comparison
    });

    std::vector<tupleRange> merged;

    merged = expandAndMergeRanges(navigator, recompRanges, seedK, data.blockExists, data.blockStrand, globalCoords);
    
    tupleCoord_t start = {0, 0, 0};
    if(data.sequence[0].first[0].second.size() == 0){
      start.nucGapPos = -1;
    }
    tupleCoord_t end = tupleCoord_t{(int64_t)data.sequence.size() - 1, (int64_t)data.sequence.back().first.size() - 1, -1};

    // Seed re-processing
    for (auto &range : std::ranges::reverse_view(merged))
    {
      bool atGlobalEnd = false;
      if (range.stop >= tupleCoord_t{(int64_t)data.sequence.size() - 1, (int64_t)data.sequence.back().first.size() - 1, -1})
      {
        atGlobalEnd = true;
        range.stop = tupleCoord_t{(int64_t)data.sequence.size() - 1, (int64_t)data.sequence.back().first.size() - 1, -1};
      }
      //Get mutated sequence to re calculate seeds 
      auto answer = seed_annotated_tree::getNucleotideSequenceFromBlockCoordinates(range.start, range.stop, data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords, navigator);
      std::string seq = std::get<0>(answer);
      std::vector<int> coords = std::get<1>(answer);
      std::vector<int> gaps  = std::get<2>(answer);
      std::vector<int> deadBlocks =  std::get<3>(answer);

      //Loop through the seeds currently in dead block and delete them
      for(int i = 0; i < deadBlocks.size(); i++){
        for (auto& pos: BlocksToSeeds[deadBlocks[i]]) {
          if (onSeedsHash[pos].has_value()) {
            // std::string oldSeed = pos < onSeeds.size() && onSeeds[pos].has_value() ? onSeeds[pos].value() : "";
            auto [oldSeed, oldIsReverse] = onSeedsHash[pos].value();
            seedChanges.emplace_back(std::make_tuple(pos,
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt));
          }
        }
      }

      //Loop through seeds that now start as gaps and delete them
      for(int i = 0; i < gaps.size(); i++){
        if (onSeedsHash[gaps[i]].has_value()){
          // std::string oldSeed = gaps[i] < onSeeds.size() && onSeeds[gaps[i]].has_value() ? onSeeds[gaps[i]].value() : "";
          auto [oldSeed, oldIsReverse] = onSeedsHash[gaps[i]].value();
          seedChanges.emplace_back(std::make_tuple(gaps[i], 
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt)); 
        }
      }

      //Loop through sequence and build up seeds as we go
      if(seq.size() >= seedK){
        // vector of (hash, isReverse, isSeed, startPos)
        std::vector<std::tuple<size_t, bool, bool, int64_t>> kmers = seeding::rollingSyncmers(seq, seedK, seedS, true, 1);

        for (int64_t i = 0; i < kmers.size(); ++i) {
          const auto& [hash, isReverse, isSeed, startPos] = kmers[i];
          bool inMap = (onSeedsHash[coords[startPos]].has_value());

          if(!inMap && isSeed){
            //Add seed
            // std::string newSeed = seq.substr(i, seedK);
            size_t newSeed = hash;
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                false, // old seed off
                                                true, // new seed on
                                                std::nullopt,
                                                newSeed,
                                                std::nullopt,
                                                isReverse));
            
          }else if(inMap && !isSeed){
            //Remove Seed
            // std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            auto [oldSeed, oldIsReverse] = onSeedsHash[coords[i]].value();
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt));

          }else if(inMap && isSeed){
            // replace seed
            // std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            auto [oldSeed, oldIsReverse] = onSeedsHash[coords[i]].value();
            // std::string newSeed = seq.substr(i, seedK);
            size_t newSeed = hash;
            if(newSeed != oldSeed || oldIsReverse != isReverse){
              seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                true, // new seed on
                                                oldSeed,
                                                newSeed,
                                                oldIsReverse,
                                                isReverse));
            }
          }
        }
      }

      //If our range reaches the end of the genome, remove seeds that aren't long enough at the end
      if(atGlobalEnd && (int)coords.size() - (int)seedK + 1 >= 0){  
        for(int i = (int)coords.size() - (int)seedK + 1; i < (int)coords.size(); i++){
          if (onSeedsHash[coords[i]].has_value()) {
            // std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            auto [oldSeed, oldIsReverse] = onSeedsHash[coords[i]].value();
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt));
          }
        }
      }
    }
    
    std::sort(seedChanges.begin(), seedChanges.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) < std::get<0>(b);
    });
    int seedChangeIndex = seedChanges.size() - 1;

    std::vector<int64_t> basePositions;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> masks_all;

    while (seedChangeIndex >= 0) {
      auto [curPos, curOldVal, curNewVal, curOldSeed, curNewSeed, curOldIsReverse, curIsReverse] = seedChanges[seedChangeIndex];
      basePositions.push_back(curPos);
      std::vector<std::pair<uint64_t, uint64_t>> masks;

      // Generate ternary numbers for the current seed position
      while (seedChangeIndex >= 0) {
        auto [maskPos, maskOldVal, maskNewVal, maskOldSeed, maskNewSeed, maskOldIsReverse, maskIsReverse] = seedChanges[seedChangeIndex];
        if (curPos - maskPos >= 32) {
          break;
        }

        int8_t ternaryNumber;
        if      (maskOldVal && maskNewVal)  ternaryNumber = 2; // changed/inserted
        else if (!maskOldVal && maskNewVal) ternaryNumber = 2; // changed/inserted
        else if (maskOldVal && !maskNewVal) ternaryNumber = 1; // deleted
        else                                ternaryNumber = 0; // same
        
        masks.emplace_back(std::make_pair(ternaryNumber, curPos - maskPos));
        --seedChangeIndex;
      }
      masks_all.push_back(masks);
    }

    size_t num_masks = 0;
    for (const auto& masks : masks_all) {
      if (masks.size() > 32) {
        throw std::runtime_error("masks.size() > 32");
      }
      num_masks += masks.size();
    }
    if (basePositions.size() != masks_all.size()) {
      std::cout << "basePositions.size(): " << basePositions.size() << " masks_all.size(): " << masks_all.size() << std::endl;
      throw std::runtime_error("basePositions.size() != masks_all.size()");
    }
    if (num_masks != seedChanges.size()) {
      std::cout << "num_masks: " << num_masks << " basePositions.size(): " << basePositions.size() << " seedChanges.size(): " << seedChanges.size() << std::endl;
      throw std::runtime_error("num_masks + basePositions.size() != seedChanges.size()");
    }

    // if (node->identifier == "node_1") {
    //   for (const auto& mask : masks_all[0]) {
    //     std::cout << "mask: " << static_cast<int>(mask.first) << " " << static_cast<unsigned int>(mask.second) << std::endl;
    //   }
    // }



    if constexpr (std::is_same_v<SeedMutationsType, ::capnp::List<SeedMutations>::Builder>) {
      auto basePositionsBuilder = perNodeSeedMutations_Index[dfsIndex].initBasePositions(basePositions.size());
      auto perPosMasksBuilder = perNodeSeedMutations_Index[dfsIndex].initPerPosMasks(masks_all.size());
      for (int i = 0; i < masks_all.size(); i++) {
        const auto &masks = masks_all[i];
        // Store the ternary numbers in the perMutMasks
        uint64_t tritMask = 0;
        for (const auto& [ternaryNumber, offset] : masks) {
          tritMask |= (ternaryNumber & 0x3) << ((offset) * 2);
        }
        perPosMasksBuilder.set(masks_all.size() - i - 1, tritMask);
        basePositionsBuilder.set(masks_all.size() - i - 1, basePositions[i]);
      }
    }

    auto nodeGapBuilder = perNodeGapMutations_Index[dfsIndex];
    if constexpr (std::is_same_v<GapMutationsType, ::capnp::List<GapMutations>::Builder>) {
      auto gapMutationsBuilder = nodeGapBuilder.initDeltas(gapMapUpdates.size());
      for (int i = 0; i < gapMapUpdates.size(); ++i) {
        const auto& [erase, range] = gapMapUpdates[i];
        if (erase) {
          gapMutationsBuilder[i].setPos(range.first);
          gapMutationsBuilder[i].initMaybeValue().setNone();
        } else {
          gapMutationsBuilder[i].setPos(range.first);
          gapMutationsBuilder[i].initMaybeValue().setValue(range.second);
        }
      }
    }
    

    merged.clear();
    recompRanges.clear();
    gapRunUpdates.clear();
    gapMapUpdates.clear();
    masks_all.clear();
    basePositions.clear();
    gapRunBlocksBacktracks.clear();
  } else if (method == Step::PLACE) {
    ::capnp::List<MapDelta>::Reader gapMutationsList;
    if constexpr (std::is_same_v<GapMutationsType, ::capnp::List<GapMutations>::Reader>) {
      gapMutationsList = perNodeGapMutations_Index[dfsIndex].getDeltas();
    }
    for (int i = 0; i < gapMutationsList.size(); ++i) {
      const auto& gapMutation = gapMutationsList[i];
      int32_t pos = gapMutation.getPos();
      auto maybeValue = gapMutation.getMaybeValue();
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

    makeCoordIndex(coordIndex, gapMap, blockRanges);

    auto currBasePositions = perNodeSeedMutations_Index[dfsIndex].getBasePositions();
    auto currPerPosMasks = perNodeSeedMutations_Index[dfsIndex].getPerPosMasks();
    std::vector<std::vector<int8_t>> masks;
    for (int i = 0; i < currBasePositions.size(); ++i) {
        int64_t pos = currBasePositions[i];
        uint64_t tritMask = currPerPosMasks[i];
        for (int k = 0; k < 32; ++k) {
          uint8_t ternaryNumber = (tritMask >> (k * 2)) & 0x3;
          if (ternaryNumber == 1) { // on -> off
            // Handle deletion
            // seedChanges.emplace_back(std::make_tuple(pos - k, true, false, onSeeds[pos - k], std::nullopt));
            auto [oldSeed, oldIsReverse] = onSeedsHash[pos - k].value();
            seedChanges.emplace_back(std::make_tuple(pos - k, true, false, oldSeed, std::nullopt, oldIsReverse, std::nullopt));
          } else if (ternaryNumber == 2) {
            // Handle insertion/change
            std::string newSeed = seed_annotated_tree::getSeedAt(pos - k, T, seedK, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap, blockRanges);
            auto [newSeedFHash, newSeedRHash] = hashSeq(newSeed);
            size_t newSeedHash;
            bool newIsReverse;
            if (newSeedFHash < newSeedRHash) {
              newSeedHash = newSeedFHash;
              newIsReverse = false;
            } else {
              newSeedHash = newSeedRHash;
              newIsReverse = true;
            }

            if (onSeedsHash[pos - k].has_value()) { // on -> on
              // seedChanges.emplace_back(std::make_tuple(pos - k, true, true, onSeeds[pos - k], newSeed));
              auto [oldSeed, oldIsReverse] = onSeedsHash[pos - k].value();
              seedChanges.emplace_back(std::make_tuple(pos - k, true, true, oldSeed, newSeedHash, oldIsReverse, newIsReverse));
            } else { // off -> on
              // seedChanges.emplace_back(std::make_tuple(pos - k, false, true, std::nullopt, newSeed));
              seedChanges.emplace_back(std::make_tuple(pos - k, false, true, std::nullopt, newSeedHash, std::nullopt, newIsReverse));
            }
          }
        }
      }

    for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
      const auto& [del, range] = *it;
      if (del) {
        gapMap.erase(range.first);
      } else {
        gapMap[range.first] = range.second;
      }
    }

    gapRunBlocksBacktracks.clear();
  }
  

  for (const auto &p : seedChanges)
  {
    const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse] = p;

    if (oldVal && newVal) { // seed at same pos changed
      // onSeeds[pos] = std::get<4>(p);
      onSeedsHash[pos].value().first = newSeed.value();
      onSeedsHash[pos].value().second = newIsReverse.value();
    } else if (oldVal && !newVal) { // seed on to off
      // if (onSeeds[std::get<0>(p)].has_value() && std::get<0>(p) < onSeeds.size()) {
      //   onSeeds[std::get<0>(p)].reset();
      // }
      if (onSeedsHash[pos].has_value() && pos < onSeedsHash.size()) {
        onSeedsHash[pos].reset();
      }
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].erase(pos);
    } else if (!oldVal && newVal) { // seed off to on
      // onSeeds[std::get<0>(p)] = std::get<4>(p);
      onSeedsHash[pos] = std::make_pair(newSeed.value(), newIsReverse.value());
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].insert(pos);
    } 

    /**
     * JACCARD: (a ∩ b) / (a ∪ b)
     * WEIGHTED_JACCARD: sum(min(f(a), f(b))) / sum(max(f(a), f(b)))
     * COSINE: sum(f(a) * f(b)) / sqrt(sum(f(a)^2) * sum(f(b)^2))
    */

    if (method == Step::PLACE) {
      if (oldVal && newVal) { // seed at same pos changed
        if (readSeedCounts.find(oldSeed) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer -= 1;
        } else {
          //case 2: seed not in reads
          jacDenom -= 1;
        }
        if (readSeedCounts.find(newSeed) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer += 1;
        } else {
          //case 2: seed not in reads
          jacDenom += 1;
        }
      } else if (oldVal && !newVal) { // seed on to off
        if (readSeedCounts.find(oldSeed) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer -= 1;
        } else {
          //case 2: seed not in reads
          jacDenom -= 1;
        }
      } else if (!oldVal && newVal) { // seed off to on
        if (readSeedCounts.find(newSeed) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer += 1;
        } else {
          //case 2: seed not in reads
          jacDenom += 1;
        }
      } 
      float jac = jacDenom == 0 ? 0 : jacNumer / (float)jacDenom;
      std::cout << node->identifier << " jac: " << jac << std::endl;
    }
   
  }
  
  if (debug) {
    // print out seeds at node
    if (method == Step::PLACE) {
      std::cout << node->identifier << " place syncmers: ";
      for (int i = 0; i < onSeedsHash.size(); i++) {
        if (onSeedsHash[i].has_value()) {
          // std::cout << degapGlobal(i, coordIndex) << ":" << onSeeds[i].value() << " ";
          std::cout << degapGlobal(i, coordIndex) << ":" << onSeedsHash[i].value().first << "|" << onSeedsHash[i].value().second << " ";
        }
      }
      std::cout << std::endl;
    } else {
      std::cout << node->identifier << " build syncmers: ";
      for (int i = 0; i < onSeedsHash.size(); i++) {
        if (onSeedsHash[i].has_value()) {
          // std::cout << degapGlobal(i, coordIndex) << ":" << onSeeds[i].value() << " ";
          std::cout << degapGlobal(i, coordIndex) << ":" << onSeedsHash[i].value().first << "|" << onSeedsHash[i].value().second << " ";
        }
      }
      std::cout << std::endl;
    }



    if (method == Step::PLACE) std::cout << node->identifier << " place coordIndex: ";
    else                       std::cout << node->identifier << " build coordIndex: ";

    for (const auto& [pos, gaps] : coordIndex) {
      std::cout << pos << ":" << gaps << " ";
    }
    std::cout << std::endl;


    // std::cout << node->identifier << " true syncmers: ";
    // auto seq = seed_annotated_tree::getStringAtNode(node, T, false);
    // auto syncmers = extractSeedmers(seq, seedK, seedS, false);
    // for (const auto &[kmer, hash, isReverse, startPos] : syncmers) {
    //   // std::cout << startPos << ":" << hash << " ";
    //   std::cout << startPos << ":" << hash << "|" << isReverse << " ";
    // }
    // std::cout << std::endl;
  }



  /* Recursive step */
  dfsIndex++;
  for (Node *child : node->children) {
    buildOrPlace(
      method, data, onSeeds, onSeedsHash, perNodeSeedMutations_Index, perNodeGapMutations_Index, seedK, seedS, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, jacNumer, jacDenom, readSeedCounts
    );
  }


  for (const auto &p : seedChanges) {
    const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse] = p;
    if (oldVal && newVal) { // UNDO seed at same pos changed
      // onSeeds[pos] = std::get<3>(p);
      onSeedsHash[pos].value().first = oldSeed.value();
      onSeedsHash[pos].value().second = oldIsReverse.value();
    } else if (oldVal && !newVal) { // seed on to off
      // onSeeds[std::get<0>(p)] = std::get<3>(p);
      onSeedsHash[pos] = std::make_pair(oldSeed.value(), oldIsReverse.value());
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].insert(pos);
    } else if (!oldVal && newVal) { // UNDO seed off to on
      // if (onSeeds[std::get<0>(p)].has_value() && std::get<0>(p) < onSeeds.size()) {
      //   onSeeds[std::get<0>(p)].reset();
      // }
      if (onSeedsHash[pos].has_value() && pos < onSeedsHash.size()) {
        onSeedsHash[pos].reset();
      }
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].erase(pos);
    } 
  }
  
  // undo gapMap updates
  for (auto it = gapRunBacktracks.rbegin(); it != gapRunBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  for (const auto& [del, blockId] : inverseBlockIdsBacktrack) {
    if (del) {
      inverseBlockIds.erase(blockId);
    } else {
      inverseBlockIds.insert(blockId);
    }
  }

  /* Undo sequence mutations when backtracking */
  undoMutations(data, T, node, blockMutationInfo, mutationInfo, globalCoords);

}


void pmi::build(Tree *T, Index::Builder &index)
{
  // Setup for seed indexing
 seed_annotated_tree::mutableTreeData data;
 seed_annotated_tree::globalCoords_t globalCoords;

 seed_annotated_tree::setup(data, globalCoords, T);
  
  CoordNavigator navigator(data.sequence);

  std::vector<int> BlockSizes(data.sequence.size(),0);
  std::vector<std::pair<int64_t, int64_t>> blockRanges(data.blockExists.size());
  std::unordered_set<int64_t> inverseBlockIds;

  int32_t k = index.getK();
  int32_t s = index.getS();
  
  std::map<int64_t, int64_t> gapMap;

  gapMap[0] = tupleToScalarCoord({blockRanges.size() - 1, globalCoords[blockRanges.size() - 1].first.size() - 1, -1}, globalCoords);
  
  std::vector<int64_t> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
  auto currCoord = tupleCoord_t{0,0,0};
  if(navigator.sequence[0].first[0].second.empty()) {
    currCoord.nucGapPos = -1;
  }

  tbb::parallel_for(tbb::blocked_range<int64_t>(0, scalarCoordToBlockId.size()), [&](const tbb::blocked_range<int64_t>& range) {
    for (int64_t i = range.begin(); i < range.end(); i++) {
      scalarCoordToBlockId[i] = currCoord.blockId;
      BlockSizes[currCoord.blockId]++;
      currCoord = navigator.newincrement(currCoord, data.blockStrand);
    }
  });

  tbb::parallel_for(tbb::blocked_range<int64_t>(0, blockRanges.size()), [&](const tbb::blocked_range<int64_t>& range) {
    for (int64_t i = range.begin(); i < range.end(); i++) {
      int64_t start = globalCoords[i].first[0].second.empty() ? tupleToScalarCoord({i, 0, -1}, globalCoords) : tupleToScalarCoord({i, 0, 0}, globalCoords);
      int64_t end = tupleToScalarCoord({i, globalCoords[i].first.size() - 1, -1}, globalCoords);
      blockRanges[i] = std::make_pair(start, end);
      if (data.blockStrand[i].first) inverseBlockIds.insert(i);
    }
  });

  std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());

  ::capnp::List<SeedMutations>::Builder perNodeSeedMutations_Builder = index.initPerNodeSeedMutations(T->allNodes.size());
  ::capnp::List<GapMutations>::Builder perNodeGapMutations_Builder = index.initPerNodeGapMutations(T->allNodes.size());

  int64_t dfsIndex = 0;
  

  // std::vector<std::optional<std::string>> onSeedsString(globalCoords.back().first.back().first + 1, std::nullopt);
  std::vector<std::optional<std::string>> onSeedsString;
  std::vector<std::optional<std::pair<size_t, bool>>> onSeedsHash(globalCoords.back().first.back().first + 1, std::nullopt);

  buildOrPlace(
    Step::BUILD, data, onSeedsString, onSeedsHash, perNodeSeedMutations_Builder, perNodeGapMutations_Builder, k, s, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, jacNumer, jacDenom, readSeedCounts
  );
}

void perfect_shuffle(std::vector<std::string>& v) {
    int n = v.size();

    std::vector<std::string> canvas(n);

    for (int i = 0; i < n / 2; i++) {
        canvas[i*2] = v[i];
        canvas[i*2+1] = v[i + n/2];
    }

    v = std::move(canvas);
}

void seedsFromFastq(const int32_t& k, const int32_t& s, std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::vector<seed> &readSeeds,  const std::string &fastqPath1, const std::string &fastqPath2) {
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath1.c_str(), "r");
    if(!fp){
        std::cerr << "Error: File " << fastqPath1 << " not found" << std::endl;
        exit(0);
    }
    seq = kseq_init(fileno(fp));
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.s);
    }
    if (fastqPath2.size() > 0) {
        fp = fopen(fastqPath2.c_str(), "r");
        if(!fp){
            std::cerr << "Error: File " << fastqPath2 << " not found" << std::endl;
            exit(0);
        }
        seq = kseq_init(fileno(fp));

        line = 0;
        int forwardReads = readSequences.size();
        while ((line = kseq_read(seq)) >= 0) {
            readSequences.push_back(seeding::reverseComplement(seq->seq.s));
            readNames.push_back(seq->name.s);
            readQuals.push_back(seq->qual.s);
        }

        if (readSequences.size() != forwardReads*2){
            std::cerr << "Error: File " << fastqPath2 << " does not contain the same number of reads as " << fastqPath1 << std::endl;
            exit(0);
        }
        
        //Shuffle reads together, so that pairs are next to eatch other
        perfect_shuffle(readSequences);
        perfect_shuffle(readNames);
        perfect_shuffle(readQuals);
    }

    for (int i = 0; i < readSequences.size(); i++) {
      for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : rollingSyncmers(readSequences[i], k, s, true, 1)) {
        if (!isSyncmer) continue;
        readSeeds.emplace_back(seed{kmerHash, startPos, -1, isReverse, startPos + k - 1});
        if (readSeedCounts.find(kmerHash) == readSeedCounts.end()) readSeedCounts[kmerHash] = std::make_pair(0, 0);
        if (isReverse) ++readSeedCounts[kmerHash].second;
        else           ++readSeedCounts[kmerHash].first;
      }
    }
}

void pmi::place(Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path)
{
    // Setup for seed indexing
    seed_annotated_tree::mutableTreeData data;
    seed_annotated_tree::globalCoords_t globalCoords;

    seed_annotated_tree::setup(data, globalCoords, T);
    
    CoordNavigator navigator(data.sequence);

    std::vector<int> BlockSizes(data.sequence.size(),0);
    std::vector<std::pair<int64_t, int64_t>> blockRanges(data.blockExists.size());
    std::unordered_set<int64_t> inverseBlockIds;


    int32_t k = index.getK();
    int32_t s = index.getS();
    
    std::map<int64_t, int64_t> gapMap;

    gapMap[0] = tupleToScalarCoord({blockRanges.size() - 1, globalCoords[blockRanges.size() - 1].first.size() - 1, -1}, globalCoords);
    
    std::vector<int64_t> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
    auto currCoord = tupleCoord_t{0,0,0};
    if(navigator.sequence[0].first[0].second.empty()) {
        currCoord.nucGapPos = -1;
    }

    for (int64_t i = 0; i < scalarCoordToBlockId.size(); i++) {
        scalarCoordToBlockId[i] = currCoord.blockId;
        BlockSizes[currCoord.blockId]++;
        currCoord = navigator.newincrement(currCoord, data.blockStrand);
    }

    for (int64_t i = 0; i < blockRanges.size(); ++i) {
        int64_t start = globalCoords[i].first[0].second.empty() ? tupleToScalarCoord({i, 0, -1}, globalCoords) : tupleToScalarCoord({i, 0, 0}, globalCoords);
        int64_t end = tupleToScalarCoord({i, globalCoords[i].first.size() - 1, -1}, globalCoords);
        blockRanges[i] = std::make_pair(start, end);
        if (!data.blockStrand[i].first) inverseBlockIds.insert(i);
    }
    
    std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());
    ::capnp::List<GapMutations>::Reader perNodeGapMutations_Reader = index.getPerNodeGapMutations();
    ::capnp::List<SeedMutations>::Reader perNodeSeedMutations_Reader= index.getPerNodeSeedMutations();
  
    int64_t dfsIndex = 0;

    
    // std::vector<std::optional<std::string>> onSeedsString(globalCoords.back().first.back().first + 1, std::nullopt);
    std::vector<std::optional<std::string>> onSeedsString;
    std::vector<std::optional<std::pair<size_t, bool>>> onSeedsHash(globalCoords.back().first.back().first + 1, std::nullopt);
    
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::vector<seed> readSeedsHash;
    std::unordered_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
    seedsFromFastq(k, s, readSeedCounts, readSequences, readQuals, readNames, readSeedsHash, reads1Path, reads2Path);

    int64_t jacNumer = 0;
    int64_t jacDenom = 0;

    // buildOrPlace<decltype(perNodeSeedMutations_Reader), decltype(perNodeGapMutations_Reader)>(
    buildOrPlace<decltype(perNodeSeedMutations_Reader), decltype(perNodeGapMutations_Reader)>(
      Step::PLACE, data, onSeedsString, onSeedsHash, perNodeSeedMutations_Reader, perNodeGapMutations_Reader, k, s, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, jacNumer, jacDenom, readSeedCounts
    );
}
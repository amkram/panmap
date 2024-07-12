#include "pmi.hpp"
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "tree.hpp"
//#include <__config>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string>
//#include <sys/_types/_int32_t.h>
//#include <sys/_types/_int64_t.h>
#include <unordered_set>

using namespace seeding;
using namespace pmi;
using namespace PangenomeMAT;
using namespace tree;

/* Helpers */

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

void applyMutations(mutableTreeData &data, seedMap_t &seedMap,
                    blockMutationInfo_t &blockMutationInfo,
                    std::vector<tupleRange> &recompRanges,
                    mutationInfo_t &mutationInfo, Tree *T, Node *node,
                    globalCoords_t &globalCoords, SeedmerIndex &index, CoordNavigator &navigator)
{
  blockExists_t &blockExists = data.blockExists;
  blockStrand_t &blockStrand = data.blockStrand;
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
      blockMutationInfo.push_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, true, !inversion));
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
        blockMutationInfo.push_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, oldMut, !oldStrand));
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
      blockMutationInfo.push_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true));
    }

    //Push new recomb range, order depends on inversion
    tupleRange newRange = {tupleCoord_t{primaryBlockId, 0, data.sequence[primaryBlockId].first[0].second.empty() ? -1 : 0},
                            tupleCoord_t{primaryBlockId, (int32_t)data.sequence[primaryBlockId].first.size() - 1, -1}};
    
    if(! blockStrand[primaryBlockId].first){
      auto temp = newRange.start;
      newRange.start = newRange.stop;
      newRange.stop = temp;
    }

    recompRanges.push_back(newRange);



    if (startingStrand != blockStrand[primaryBlockId].first) {
      flipCoords(primaryBlockId, globalCoords);
    }




  }

  // Nuc mutations
  for (size_t i = 0; i < node->nucMutation.size(); i++)
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


      tupleRange newRange = {tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition + len}};
    

      if (nucGapPosition == -1){
        //recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition+len, -1}});
        newRange = {tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition+len, -1}};
      }

      if(! blockStrand[primaryBlockId].first){
        auto temp = newRange.start;
        newRange.start = newRange.stop;
        newRange.stop = temp;
      }

      recompRanges.push_back(newRange);



      if (type == PangenomeMAT::NucMutationType::NS)
      {
        // Substitution

        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
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
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            }
          }
        }
      }
      else if (type == PangenomeMAT::NucMutationType::NI)
      {
        // Insertion
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j];
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
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
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
            }
          }
        }
      }
      else if (type == PangenomeMAT::NucMutationType::ND)
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
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, '-'));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
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
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, '-'));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              sequence[primaryBlockId].first[nucPosition + j].first = '-';
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, '-'));
            }
          }
        }
      }
    }
    else
    {
      int len = 0;
      

      //recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},
      //                        tupleCoord_t{primaryBlockId, nucPosition + len, nucGapPosition}});

      tupleRange newRange ={tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},
                              tupleCoord_t{primaryBlockId, nucPosition + len, nucGapPosition}};
    
      if(! blockStrand[primaryBlockId].first){
        auto temp = newRange.start;
        newRange.start = newRange.stop;
        newRange.stop = temp;
      }

      recompRanges.push_back(newRange);




      if (type == PangenomeMAT::NucMutationType::NSNPS)
      {
        // SNP Substitution
        newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].first;
            sequence[primaryBlockId].first[nucPosition].first = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
      }
      else if (type == PangenomeMAT::NucMutationType::NSNPI)
      {
        // SNP Insertion
        newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].first;
            sequence[primaryBlockId].first[nucPosition].first = newVal;
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, newVal));
          }
        }
      }
      else if (type == PangenomeMAT::NucMutationType::NSNPD)
      {
        // SNP Deletion
        if (secondaryBlockId != -1)
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
        }
        else
        {
          if (nucGapPosition != -1)
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].second[nucGapPosition];
            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
          else
          {
            char oldVal = sequence[primaryBlockId].first[nucPosition].first;
            sequence[primaryBlockId].first[nucPosition].first = '-';
            mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition, oldVal, '-'));
          }
        }
      }
    }
  }
}

void undoMutations(mutableTreeData &data, SeedmerIndex &index, Tree *T,
                   const Node *node, const blockMutationInfo_t &blockMutationInfo,
                   const mutationInfo_t &mutationInfo, globalCoords_t &globalCoords )
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
          coord = tupleCoord_t{coord.blockId - 1, navigator.sequence[coord.blockId - 1].first.size() - 1, -1};
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
          coord = tupleCoord_t{coord.blockId + 1, navigator.sequence[coord.blockId + 1].first.size() - 1, -1};
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

  for (size_t i = 1; i < ranges.size(); ++i) {
    
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
      merged.push_back(current);
      current = ranges[i];
    }
  }
  merged.push_back(current);

  ranges = merged;

  merged.clear();


  current = {
      expandLeft(navigator, ranges[0].start, neededNongap, blockExists, blockStrand),
      expandRight(navigator, ranges[0].stop, neededNongap, blockExists, blockStrand),
  };


  for (size_t i = 1; i < ranges.size(); ++i) {

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
      merged.push_back(current);
      current = expandedRange;
    }
  }
  merged.push_back(current);


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

// Recursive function to build the seed index
void buildHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index,
                 Tree *T, Node *node, globalCoords_t &globalCoords,
                 CoordNavigator &navigator, std::vector<int> &scalarCoordToBlockId, std::vector<std::unordered_set<int>> &BlocksToSeeds, std::vector<int> &BlockSizes)
{


  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;

  // First, a range is made marking the start -> end
  // of each block and nuc mutation. This is done while
  // applying mutations to the sequence object.
  std::vector<tupleRange> recompRanges;
  applyMutations(data, seedMap, blockMutationInfo, recompRanges, mutationInfo, T, node,
                 globalCoords, index, navigator);

  

  //std::sort(recompRanges.begin(), recompRanges.end());
  std::sort(recompRanges.begin(), recompRanges.end(), [&data](const tupleRange& A, const tupleRange& B) {
    if (A.start.blockId == B.start.blockId && !data.blockStrand[A.start.blockId].first) {
        return B < A; // Use B < A if blocks are inverted 
    }
    return A < B; // Default comparison
  });


  std::vector<int> seedsToClear; // seeds to clear from seedMap
  std::vector<std::pair<int, std::string>> addSeeds;
  std::vector<std::pair<int, std::string>> backtrack;

  std::vector<tupleRange> merged;

  


  merged = expandAndMergeRanges(navigator, recompRanges, index.k(), data.blockExists, data.blockStrand, globalCoords);
  


  
  tupleCoord_t start = {0, 0, 0};
  if(data.sequence[0].first[0].second.size() == 0){
    start.nucGapPos = -1;
  }
  tupleCoord_t end = tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1};



  // Protobuf message for this node's mutations
  NodeSeedmerMutations *pb_node_mutations = index.add_per_node_mutations();
  pb_node_mutations->set_node_id(node->identifier);

 
  //std::cout << "ALLOW\n";
  // Seed re-processing
  for (auto &range : std::ranges::reverse_view(merged))
  {

    bool atGlobalEnd = false;
    if (range.stop >= tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1})
    {
      atGlobalEnd = true;
      range.stop = tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1};
    }

    
    

    //Get mutated sequence to re calculate seeds 

    auto answer = tree::getNucleotideSequenceFromBlockCoordinates(range.start, range.stop, data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords, navigator);
    std::string seq = std::get<0>(answer);
    std::vector<int> coords = std::get<1>(answer);
    std::vector<int> gaps  = std::get<2>(answer);
    std::vector<int> deadBlocks =  std::get<3>(answer);
    std::string kmer;

    //Loop through the seeds currently in dead block and delete them
    for(int i = 0; i < deadBlocks.size(); i++){

      for (auto& pos: BlocksToSeeds[deadBlocks[i]]) {
        backtrack.push_back(std::make_pair(pos, seedMap[pos]));
        seedsToClear.push_back(pos);
        if (seedMap[pos].size() == index.k())
          {
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(pos);
            pb_mut->set_seq(seedMap[pos]);
        }
      }
    }


    //Loop through seeds that now start as gaps and delete them
    for(int i = 0; i < gaps.size(); i++){
      if (seedMap.find(gaps[i]) != seedMap.end()){
        SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
        pb_mut->set_is_deletion(true);
        pb_mut->set_pos(gaps[i]);
        pb_mut->set_seq(seedMap[gaps[i]]);

        backtrack.push_back(std::make_pair(gaps[i], seedMap[gaps[i]]));
        seedsToClear.push_back(gaps[i]);
      }
    }




    //Loop through sequence and build up seeds as we go
    
    if(seq.size() >= index.k()){

      //Check first k bp for smers
      std::string min_s = seq.substr(0,index.s());
      int first_min_coord = 0;
      int last_min_coord = 0;
      int Ns = seq[0] == 'N';
      
      for(int i = 1; i < index.k() - index.s() + 1; i++){
        std::string smer = seq.substr(i,index.s());
        if(seq[i]=='N')
          Ns++;
        
        if(smer < min_s){
          min_s = smer;
          first_min_coord = i;
          last_min_coord = i;
        }else if(smer == min_s){
          last_min_coord = i;
        }
        
      }
      for(int i = index.k() - index.s() + 1; i < index.k() ; i++){
        if(seq[i]=='N')
          Ns++;
      }

      

      //Check the rest for smers and kmer seeds
      for(int i = index.k(); i < seq.size() + 1; i++) {
        
        //Processing kmer starting at i - k in seq

        bool inMap = (seedMap.find(coords[i - index.k()]) != seedMap.end());
        bool isSeed = (first_min_coord == i - index.k() || last_min_coord == i - index.s()) && Ns <= index.k()/2;


        if(!inMap && isSeed){
          //Add seed
          std::string kmer = seq.substr(i - index.k(), index.k());

          backtrack.push_back(std::make_pair(coords[i - index.k()], ""));
          addSeeds.push_back(std::make_pair(coords[i - index.k()], kmer));
          
          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(false);
          pb_mut->set_pos(coords[i - index.k()]);
          pb_mut->set_seq(kmer);

        }else if(inMap && !isSeed){
          //Remove Seed

          backtrack.push_back(std::make_pair(coords[i - index.k()], seedMap[coords[i - index.k()]]));
          seedsToClear.push_back(coords[i - index.k()]);

          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(true);
          pb_mut->set_pos(coords[i - index.k()]);
          pb_mut->set_seq(seedMap[coords[i - index.k()]]);    

        }else if(inMap && isSeed){
          //Changed Seed
          std::string kmer = seq.substr(i - index.k(), index.k());
          if(kmer != seedMap[coords[i - index.k()]]){

          backtrack.push_back(std::make_pair(coords[i - index.k()], seedMap[coords[i - index.k()]]));
          addSeeds.push_back(std::make_pair(coords[i - index.k()], kmer));

          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(false);
          pb_mut->set_pos(coords[i - index.k()]);
          pb_mut->set_seq(kmer);
          }

        }

        
        //Updating smers for kmer starting at i - index.k() + 1
        if(i < seq.size()){
          if(seq[i] == 'N')
            Ns++;
          if(seq[i - index.k()] == 'N')
            Ns--;

          if(first_min_coord == i - index.k()){
            // Were losing the lowest smer, Re search for lowest
            min_s = seq.substr(i - index.k() + 1,index.s());
            first_min_coord = i - index.k() + 1;

            for(int j = 1; j < index.k() - index.s() + 1; j++){
              std::string smer = seq.substr(i - index.k() + 1 + j,index.s());
              if(smer < min_s){
                min_s = smer;
                first_min_coord = i - index.k() + 1 + j;
                last_min_coord = first_min_coord;
              }else if(smer == min_s){
                last_min_coord = i - index.k() + 1 + j;
              }
            }

          }else{
            //Test new smer to see if its the lowest
            std::string smer = seq.substr(i - index.s() + 1,index.s());
            if(smer < min_s){
              min_s = smer;
              first_min_coord = i - index.s() + 1;
              last_min_coord = i - index.s() + 1;
            }else if(smer == min_s){
              last_min_coord = i - index.s() + 1;
            }
          }
        }

      }

    }


    
    

    //If our range reaches the end of the genome, remove seeds that aren't long enough at the end
    if(atGlobalEnd && (int)coords.size() - (int)index.k() + 1 >= 0){
      for(int i = (int)coords.size() - (int)index.k() + 1; i < (int)coords.size(); i++){
        if (seedMap.find(coords[i]) != seedMap.end()){
          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(true);
          pb_mut->set_pos(coords[i]);
          pb_mut->set_seq(seedMap[coords[i]]);

          backtrack.push_back(std::make_pair(coords[i], seedMap[coords[i]]));
          seedsToClear.push_back(coords[i]);
        }
      }
    }



  } 

  
  for (const auto &pos : seedsToClear)
  {
    if (seedMap.find(pos) != seedMap.end())
    {

      seedMap.erase(pos);

      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].erase(pos);
    }
  }

  
  
  for (const auto &seed : addSeeds) {
    seedMap[seed.first] = seed.second;

    int blockId = scalarCoordToBlockId[seed.first];
    BlocksToSeeds[blockId].insert(seed.first);

  }

  


  /* Recursive step */
  for (Node *child : node->children) {
    
    buildHelper(data, seedMap, index, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes);
  }

  
  // undo seed mutations
  for (const auto &back : backtrack)
  {
    int blockId = scalarCoordToBlockId[back.first];

    if(back.second == ""){
      seedMap.erase(back.first);

      BlocksToSeeds[blockId].erase(back.first);
    }else{
      seedMap[back.first] = back.second;

      BlocksToSeeds[blockId].insert(back.first);
    }
  }

  /* Undo sequence mutations when backtracking */
  undoMutations(data, index, T, node, blockMutationInfo, mutationInfo, globalCoords);
}





/* implementation */
void pmi::build(SeedmerIndex &index, Tree *T, int j, int k, int s)
{
  
  // Setup for seed indexing
  tree::mutableTreeData data;
  tree::globalCoords_t globalCoords;

  tree::setup(data, globalCoords, T);
  

  

  index.set_j(j);
  index.set_k(k);
  index.set_s(s);


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




  std::vector<int> BlockSizes(data.sequence.size(),0);
  

  std::vector<int> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
  auto currCoord = tupleCoord_t{0,0,0};
  if(navigator.sequence[0].first[0].second.empty()) {
    currCoord.nucGapPos = -1;
  }


  for(int i = 0; i < scalarCoordToBlockId.size(); i++){
    
    scalarCoordToBlockId[i] = currCoord.blockId;
    BlockSizes[currCoord.blockId] ++; 
    currCoord = navigator.newincrement(currCoord, data.blockStrand);
  }
  

  std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());
  

  /* Recursive traversal of tree to build the index */
  buildHelper(data, seedMap, index, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes);
}

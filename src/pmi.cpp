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
#include <sys/_types/_int64_t.h>
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

    if (inversion) {
      flipCoords(primaryBlockId, globalCoords);
    }

    recompRanges.push_back({tupleCoord_t{primaryBlockId, 0, data.sequence[primaryBlockId].first[0].second.empty() ? -1 : 0},
                            tupleCoord_t{primaryBlockId, (int32_t)data.sequence[primaryBlockId].first.size() - 1, -1}});

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


      if (nucGapPosition != -1){
        recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition + len}});
      }else{
        recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition+len, -1}});
      }


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
      

      recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},
                              tupleCoord_t{primaryBlockId, nucPosition + len, nucGapPosition}});

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
                   const mutationInfo_t &mutationInfo)
{
  auto &sequence = data.sequence;
  auto &blockExists = data.blockExists;
  auto &blockStrand = data.blockStrand;
 // Undo block mutations when current node and its subtree have been processed
    for(auto it = blockMutationInfo.rbegin(); it != blockMutationInfo.rend(); it++) {
        auto mutation = *it;
        if(std::get<1>(mutation) != -1) {
            blockExists[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<2>(mutation);
            blockStrand[std::get<0>(mutation)].second[std::get<1>(mutation)] = std::get<3>(mutation);
        } else {
            blockExists[std::get<0>(mutation)].first = std::get<2>(mutation);
            blockStrand[std::get<0>(mutation)].first = std::get<3>(mutation);
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
                        int neededNongap, blockExists_t &blockExists)
{
  
  int count = 0;

  while (count < neededNongap && coord > tupleCoord_t{0, 0, 0})
  {

    if (!blockExists[coord.blockId].first)
    {
      //Jump down to previous block
      if(coord.blockId == 0){
        return tupleCoord_t{0,0,0};
      }else{
        coord = tupleCoord_t{coord.blockId - 1, navigator.sequence[coord.blockId - 1].first.size() - 1, -1};
      }
      
      continue;
    }
    
    if (!navigator.isGap(coord))
    {
      count++;
    }

    coord = navigator.decrement(coord);

  }
  return coord;
}

// Go downstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandRight(CoordNavigator &navigator, tupleCoord_t &coord,
                         int neededNongap, blockExists_t &blockExists)
{

  int count = 0;

  while (count < neededNongap && coord < tupleCoord_t{-1, -1, -1})
  {

    if (!blockExists[coord.blockId].first)
    {

      //Jump to next block
      if(coord.blockId == navigator.sequence.size() - 1){
        return tupleCoord_t{-1,-1,-1};
      }else{

        coord.blockId += 1;
        coord.nucPos = 0;
        coord.nucGapPos = 0;
        if(navigator.sequence[coord.blockId].first[0].second.empty()) {
          coord.nucGapPos = -1;
        }

      }

      continue;
    }

    if (!navigator.isGap(coord))
    {
      count++;
    }
    

    coord = navigator.increment(coord);
  }

  return coord;
}




// Merges each range with overlapping ranges after expanding left and right
// by `neededNongap` non-gap nucleotides.
std::vector<tupleRange> expandAndMergeRanges(CoordNavigator &navigator,
                                             std::vector<tupleRange> &ranges,
                                             int neededNongap,
                                             blockExists_t &blockExists)
{

  if (ranges.empty())
    return {};

  std::vector<tupleRange> merged;


  tupleRange current = {
      expandLeft(navigator, ranges[0].start, neededNongap, blockExists),
      expandRight(navigator, ranges[0].stop, neededNongap, blockExists),
  };


  for (size_t i = 1; i < ranges.size(); ++i) {

    
    tupleRange expandedRange = {
        expandLeft(navigator, ranges[i].start, neededNongap, blockExists),
        expandRight(navigator, ranges[i].stop, neededNongap, blockExists),
    };

    if (expandedRange.start <= current.stop) //Merge rangess
    {

      current.stop = std::max(current.stop, expandedRange.stop);
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
                 CoordNavigator &navigator, std::vector<int> &scalarCoordToBlockId, std::vector<std::unordered_set<int>> &BlocksToSeeds)
{

  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;

  // First, a range is made marking the start -> end
  // of each block and nuc mutation. This is done while
  // applying mutations to the sequence object.
  std::vector<tupleRange> recompRanges;
  applyMutations(data, seedMap, blockMutationInfo, recompRanges, mutationInfo, T, node,
                 globalCoords, index, navigator);
  

  std::sort(recompRanges.begin(), recompRanges.end());

  std::vector<int> seedsToClear; // seeds to clear from seedMap
  std::vector<std::pair<int, std::string>> addSeeds;
  std::vector<std::pair<int, std::string>> backtrack;

  std::vector<tupleRange> merged;
  merged = expandAndMergeRanges(navigator, recompRanges, index.k(), data.blockExists);

  // REMOVE
  tupleCoord_t start = {0, 0, 0};
  tupleCoord_t end = {-1, -1, -1};
  tupleRange fullseqRange = {start, end};
  merged.clear();
  merged.push_back(fullseqRange);
    // YEAH

// Protobuf message for this node's mutations
  NodeSeedmerMutations *pb_node_mutations = index.add_per_node_mutations();
  pb_node_mutations->set_node_id(node->identifier);

  std::cout << node->identifier  << std::endl;

  if (node->identifier == "MT506703.1") {
    std::cout << "MT506703.1" << std::endl;
    for (int32_t i = 0; i < globalCoords.size(); i++) {
      std::cout << "Block " << i << " ";
      std::cout << "STRAND: " << data.blockStrand[i].first << "\n";
      for (int32_t j = 0; j < globalCoords[i].first.size(); j++) {
        for (int32_t k = 0; k < globalCoords[i].first[j].second.size(); k++) {
          std::cout << globalCoords[i].first[j].second[k] << " ";
        }
        std::cout << globalCoords[i].first[j].first << " ";
      }
      std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "OLD INC\n\n";
    int i = 0;
    for(auto coord = start; coord < end; coord = navigator.increment(coord)) {
      std::cout << "OLD INC " <<  coord.blockId << " " << coord.nucPos << " " << coord.nucGapPos << " " << tupleToScalarCoord(coord, globalCoords) << " " << i << std::endl;
      i++;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "New INC\n\n";
    i = 0;
    for(auto coord = start; coord < end; coord = navigator.newincrement(coord, data.blockStrand)) {
      std::cout << "New INC " << coord.blockId << " " << coord.nucPos << " " << coord.nucGapPos << " " << tupleToScalarCoord(coord, globalCoords) << " " << i << std::endl;
      i++;
    }

  }

  // Seed re-processing
  for (auto &range : std::ranges::reverse_view(merged))
  {

    std::string recomputeSeq = tree::getNucleotideSequenceFromBlockCoordinates(range.start, range.stop, data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords, navigator);
    
    // Track the last downstream seed to stack k-mers into seedmers
    tupleCoord_t lastDownstreamSeedPos = range.stop;

    bool atGlobalEnd = false;
    if (range.stop >= tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1})
    {
      atGlobalEnd = true;
      range.stop = tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1};
    }

    
    int32_t seen_non_gap = 0;
    int32_t str_i = tupleToScalarCoord(range.stop, globalCoords) - tupleToScalarCoord(range.start, globalCoords);
    int32_t startScalar = tupleToScalarCoord(range.start, globalCoords);


    for ( ; str_i >= 0; str_i--)
    {
      if (str_i < 0)
      {
        break;
      }
      char nt = recomputeSeq[str_i];


      if (!data.blockExists[scalarCoordToBlockId[str_i + startScalar]].first) //Block doesnt exist, remove seeds
      {

        //Loop through the deleted seeds
        for (auto& pos: BlocksToSeeds[scalarCoordToBlockId[str_i + startScalar]]) {
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
        
        if(scalarCoordToBlockId[str_i + startScalar] > 0){
          str_i = tupleToScalarCoord(tupleCoord_t{scalarCoordToBlockId[str_i + startScalar] - 1, data.sequence[scalarCoordToBlockId[str_i + startScalar] - 1].first.size() - 1, -1}, globalCoords)  - startScalar;
        }else{
          break;
        }

      }

      if (seen_non_gap < index.k()) 
      {
        //Seed in map yet we dont have enough non-gaps for a seed, so we remove i
        if (atGlobalEnd && seedMap.find(str_i + startScalar) != seedMap.end())
        {
          backtrack.push_back(std::make_pair(str_i + startScalar, seedMap[str_i + startScalar]));
          seedsToClear.push_back(str_i + startScalar);
          if (seedMap[str_i + startScalar].size() == index.k())
          {
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(str_i + startScalar);
            pb_mut->set_seq(seedMap[str_i + startScalar]);
          }
        }

        if (recomputeSeq[str_i] != '-' && recomputeSeq[str_i] != 'x') {
          seen_non_gap++;
        }
        
        if (seen_non_gap < index.k() && str_i > 0) {
                                 
          continue; 
        }
      }

      
      if (data.blockExists[scalarCoordToBlockId[str_i + startScalar]].first && recomputeSeq[str_i] == '-' ||
               recomputeSeq[str_i] ==
                   'x')
      { // block does exist but seq is a gap
        
        if (seedMap.find(str_i + startScalar) != seedMap.end())
        {
          
          // is a gap, no longer a seed -> delete

          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(true);
          pb_mut->set_pos(str_i + startScalar);
          pb_mut->set_seq(seedMap[str_i + startScalar]);

          backtrack.push_back(std::make_pair(str_i + startScalar, seedMap[str_i + startScalar]));
          seedsToClear.push_back(str_i + startScalar);

        } /* else: no seed, wasn't seed, no change */

      }
      else if(data.blockExists[scalarCoordToBlockId[str_i + startScalar]].first)
      {
        // block exists and seq is not a gap at currCoord

        std::string kmer = "";
        int64_t seen_k = 0;
        int64_t k_pos = str_i;
        while (seen_k < index.k() && k_pos < recomputeSeq.size())
        {
          if (recomputeSeq[k_pos] != '-' && recomputeSeq[k_pos] != 'x')
          {
            kmer += recomputeSeq[k_pos];
            seen_k++;
          }
          k_pos++;
        }
      
        
        if (seedMap.find(str_i + startScalar) != seedMap.end())
        {
          // non gap position and kmer is already a seed.
          //std::string prevseedmer = lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1} ? seedMap[lastDownstreamSeedPos] : "";
          
          if (seeding::is_syncmer(kmer, index.s(), false))
          {
            // Is it still a seed?
            
            backtrack.push_back(std::make_pair(str_i + startScalar, seedMap[str_i + startScalar]));
            addSeeds.push_back(std::make_pair(str_i + startScalar, kmer));

            //lastDownstreamSeedPos = currCoord;
            if (kmer.size() == index.k())
            {
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(str_i + startScalar);
              pb_mut->set_seq(kmer);
            }
          }
          else
          {
            backtrack.push_back(std::make_pair(str_i + startScalar, seedMap[str_i + startScalar]));

            // no longer a seed -> delete

            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(str_i + startScalar);
            pb_mut->set_seq(seedMap[str_i + startScalar]);

            seedsToClear.push_back(str_i + startScalar);
          }
        }
        else
        {
          //  not in seed map, could be a seed now
          if (seeding::is_syncmer(kmer, index.s(), false))
          {
            backtrack.push_back(std::make_pair(str_i + startScalar, ""));

            //std::string prevseedmer = lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1} ? seedMap[lastDownstreamSeedPos] : "";

            addSeeds.push_back(std::make_pair(str_i + startScalar, kmer));
            //seedMap[currCoord] = kmer + prevseedmer.substr(0, (index.j() - 1) * index.k());
            if (kmer.size() == index.k())
            {
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(str_i + startScalar);
              pb_mut->set_seq(kmer);
            }
            //lastDownstreamSeedPos = currCoord;
          }
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
    
    buildHelper(data, seedMap, index, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds);
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
  undoMutations(data, index, T, node, blockMutationInfo, mutationInfo);
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

  std::cout << "Building index" << std::endl;

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


  std::vector<int> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
  auto currCoord = tupleCoord_t{0,0,0};

  for(int i = 0; i < scalarCoordToBlockId.size(); i++){
    
    scalarCoordToBlockId[i] = currCoord.blockId;

    currCoord = navigator.increment(currCoord);
  }


  std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());

  /* Recursive traversal of tree to build the index */
  buildHelper(data, seedMap, index, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds);
}

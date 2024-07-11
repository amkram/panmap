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
  


  // Protobuf message for this node's mutations
  NodeSeedmerMutations *pb_node_mutations = index.add_per_node_mutations();
  pb_node_mutations->set_node_id(node->identifier);


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

/* panmam stuff */

void updateCoordsInfo(
  std::vector<int32_t>& degap, std::vector<int32_t>& regap, std::map<int32_t, int32_t>& coordsIndex,
  const mutableTreeData& data, const std::vector<std::pair<size_t, size_t>>& blocksEndCoords,
  size_t& ungappedSeqSize, Tree *T, Node *node
  ) {

  std::string gappedSeq = tree::getStringAtNode(node, T, true);
  std::string ungappedSeq = "";
  int32_t numGaps = 0;

  for (int32_t i = 0; i < gappedSeq.size(); i++) {
    char &c = gappedSeq[i];
    char &p = gappedSeq[std::max(0, i-1)];
    degap.push_back(ungappedSeq.size());
    if (c != '-') {
      ungappedSeq += c;
      regap.push_back(i);
      if (i == 0 || p == '-') {
        coordsIndex[i] = numGaps;
      }
    } else {
      ++numGaps;
    }
  }

  // degap.resize(gappedSeq.size(), -1);
  // assert(data.blockExists.size() == blocksEndCoords.size());
  // for (size_t blockId = 0; blockId < blocksEndCoords.size(); ++blockId) {
  //   auto [blockStart, blockSize] = blocksEndCoords[blockId];
  //   if (blockSize == std::numeric_limits<size_t>::max()) blockSize = gappedSeq.size() - blockStart;
    
  //   if (data.blockExists[blockId].first) {
  //     char c = gappedSeq[blockStart];
  //     char p = (blockId == 0) ? gappedSeq[0] : gappedSeq[blocksEndCoords[blockId-1].first + blocksEndCoords[blockId-1].second - 1];
  //     degap[blockStart] = ungappedSeq.size();

  //     if (c != '-') {
  //       ungappedSeq += c;
  //       regap.push_back(blockStart);
  //       if (blockStart == 0 || p == '-') {
  //         coordsIndex[blockStart] = numGaps;
  //       }
  //     } else {
  //       ++numGaps;
  //     }

  //     for (size_t i = blockStart + 1; i < blockStart + blockSize; ++i) {
  //       c = gappedSeq[i];
  //       p = gappedSeq[i-1];
  //       degap[i] = ungappedSeq.size();

  //       if (c != '-') {
  //         ungappedSeq += c;
  //         regap.push_back(i);
  //         if (i == 0 || p == '-') {
  //           coordsIndex[i] = numGaps;
  //         }
  //       } else {
  //         ++numGaps;
  //       }
  //     }
  //   } else {
  //     numGaps += blockSize;
  //   }
  // }


  ungappedSeqSize = ungappedSeq.size();


}

void initializeSeedmersIndex(
  mgsr::seedmers& seedmersIndex, const auto& seedsChanges, std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds,
  std::vector<std::tuple<size_t, int32_t, int32_t, bool, bool>>& seedmersChangesBackTrack, const int& l, const int& k, std::stringstream& seedmersOutStream
  ) {
  // first seedmer
  int startIdx = 0;
  size_t cacheForwardH  = 0;
  size_t cacheReversedH = 0;
  size_t cacheMin;
  bool   rev;
  while (seedmersIndex.positionMap.empty()) {
    assert(std::get<3>(seedsChanges[startIdx]) == false);
    curSeeds[std::get<1>(seedsChanges[startIdx])] = std::make_pair(std::get<2>(seedsChanges[startIdx]),std::get<0>(seedsChanges[startIdx]));
    
    cacheForwardH  = 0;
    cacheReversedH = 0;
    for (int i = startIdx; i < startIdx + l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(seedsChanges[i]);
    for (int i = startIdx + l - 1; i > startIdx - 1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(seedsChanges[i]);

    // skip direction ambiguous
    if (cacheForwardH < cacheReversedH) {
      cacheMin = cacheForwardH;
      rev = false;
    } else if (cacheReversedH < cacheForwardH) {
      cacheMin = cacheReversedH;
      rev = true;
    } else {
      ++startIdx;
      continue;
    }

    // seedmersIndex.positionMap[std::get<1>(seedsChanges[startIdx])] = std::make_pair(cacheMin, rev);
    seedmersIndex.positionMap[std::get<1>(seedsChanges[startIdx])] = std::make_tuple(std::get<2>(seedsChanges[startIdx+l-1]), cacheMin, rev);
    seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(seedsChanges[startIdx]));
    seedmersOutStream << std::get<1>(seedsChanges[startIdx]) << ":+:"
                      << std::get<2>(seedsChanges[startIdx+l-1]) << ","
                      << cacheMin << ","
                      << rev << " ";
    
    ++startIdx;
  }

  size_t mask = 0;
  for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

  // rest of seedmers
  for (int i = startIdx; i < seedsChanges.size() - l + 1; ++i) {
    assert(std::get<3>(seedsChanges[i]) == false);

    curSeeds[std::get<1>(seedsChanges[i])] = std::make_pair(std::get<2>(seedsChanges[i]),std::get<0>(seedsChanges[i]));

    cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(seedsChanges[i+l-1]);
    cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(seedsChanges[i+l-1]) << (2 * k * (l - 1)));

    // skip direction ambiguous
    if (cacheForwardH < cacheReversedH) {
      cacheMin = cacheForwardH;
      rev = false;
    } else if (cacheReversedH < cacheForwardH) {
      cacheMin = cacheReversedH;
      rev = true;
    } else {
      continue;
    }

    // seedmersIndex.positionMap[std::get<1>(seedsChanges[i])] = std::make_pair(cacheMin, rev);
    seedmersIndex.positionMap[std::get<1>(seedsChanges[i])] = std::make_tuple(std::get<2>(seedsChanges[i+l-1]), cacheMin, rev);
    seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(seedsChanges[i]));
    seedmersOutStream << std::get<1>(seedsChanges[i]) << ":+:"
                      << std::get<2>(seedsChanges[i+l-1]) << ","
                      << cacheMin << ","
                      << rev << " ";
  }

  // rest of seeds
  for (int i = seedsChanges.size() - l + 1; i < seedsChanges.size(); i++) {
    assert(std::get<3>(seedsChanges[i]) == false);
    curSeeds[std::get<1>(seedsChanges[i])] = std::make_pair(std::get<2>(seedsChanges[i]),std::get<0>(seedsChanges[i]));
  }
}

void updateCurSeeds(std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, const auto& seedsChanges, const auto& seedGlobalEndCoorChanges) {
  for (const auto& change : seedsChanges) {
    auto [seedHash, seedBeg, seedEnd, seedDel] = change;
    if (seedDel) {
      assert(curSeeds.count(seedBeg));
      curSeeds.erase(seedBeg);
    } else {
      curSeeds[seedBeg] = std::make_pair(seedEnd, seedHash);
    }
  }

  for (const auto& endCoorChange : seedGlobalEndCoorChanges) {
    assert(curSeeds.count(endCoorChange.first));
    curSeeds[endCoorChange.first].first = endCoorChange.second;
  }
}

int32_t degapGlobalFromIndex(const int32_t& globalCoord, const std::map<int32_t, int32_t>& coordsIndex) {
    auto coordIt = coordsIndex.upper_bound(globalCoord);
    if (coordIt == coordsIndex.begin()) {
        return 0;
    }
    --coordIt;
    return globalCoord - coordIt->second;
}

void writeTrueSeeds(const std::string& seq, const int32_t& k, const int32_t& s, const std::string& outFile) {
  std::ofstream of(outFile);
  for (size_t i = 0; i < seq.size() - k + 1; ++i) {
    const std::string& kmer = seq.substr(i, k);
    if (seeding::is_syncmer(kmer, s, false)) {
      auto [seedHash, seedLegal] = getHash(kmer);
      if (seedLegal) {
        of << i << "\t"
           << i + k - 1 << "\t"
           << seedHash << std::endl;
      }
    }
  }
}

static std::vector<std::tuple<size_t, int32_t, int32_t, bool>> syncmersSketch(const std::string& seq, const int k, const int s, const bool open) {
    std::vector<std::tuple<size_t, int32_t, int32_t, bool>> syncmers;
    for (int i = 0; i < static_cast<int>(seq.size()) - k + 1; ++i) {
        std::string kmer = seq.substr(i, k);
        if (!seeding::is_syncmer(kmer, s, open)) continue;

        std::pair<size_t, bool> kmerHash = getHash(kmer);
        if (!kmerHash.second) continue;

        syncmers.emplace_back(std::make_tuple(kmerHash.first, i, i + k - 1, false));
    }
    return syncmers;
}



void writeTrueSeedmers(const std::string& seq, const int32_t& k, const int32_t& s, const int32_t& l, const std::string& outFile) {
  auto syncmers = syncmersSketch(seq, k, s, false);
  mgsr::seedmers seedmersIndex;

  // first seedmer
  int startIdx = 0;
  size_t cacheForwardH  = 0;
  size_t cacheReversedH = 0;
  size_t cacheMin;
  bool   rev;
  while (seedmersIndex.positionMap.empty() && startIdx < static_cast<int>(syncmers.size())) {
      assert(std::get<3>(syncmers[startIdx]) == false);
      cacheForwardH  = 0;
      cacheReversedH = 0;
      bool breakWhile = false;
      for (int i = startIdx; i < startIdx + l; ++i) {
        if (i >= syncmers.size()) {
          breakWhile = true;
          break;
        }
        cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(syncmers[i]);
      }
      if (breakWhile) break;
      for (int i = startIdx + l - 1; i > startIdx - 1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(syncmers[i]);
  
      // skip direction ambiguous
      if (cacheForwardH < cacheReversedH) {
          cacheMin = cacheForwardH;
          rev = false;
      } else if (cacheReversedH < cacheForwardH) {
          cacheMin = cacheReversedH;
          rev = true;
      } else {
          ++startIdx;
          continue;
      }

      // seedmersIndex.positionMap[std::get<1>(syncmers[startIdx])] = std::make_pair(cacheMin, rev);
      seedmersIndex.positionMap[std::get<1>(syncmers[startIdx])] = std::make_tuple(std::get<2>(syncmers[startIdx+l-1]), cacheMin, rev);
      seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(syncmers[startIdx]));
      ++startIdx;
  }

  size_t mask = 0;
  for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

  // rest of seedmers
  for (int i = startIdx; i < static_cast<int>(syncmers.size()) - l + 1; ++i) {
      assert(std::get<3>(syncmers[i]) == false);

      cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);
      cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(syncmers[i+l-1]) << (2 * k * (l - 1)));

      // skip direction ambiguous
      if (cacheForwardH < cacheReversedH) {
          cacheMin = cacheForwardH;
          rev = false;
      } else if (cacheReversedH < cacheForwardH) {
          cacheMin = cacheReversedH;
          rev = true;
      } else {
          continue;
      }

      // seedmersIndex.positionMap[std::get<1>(syncmers[i])] = std::make_pair(cacheMin, rev);
      seedmersIndex.positionMap[std::get<1>(syncmers[i])] = std::make_tuple(std::get<2>(syncmers[i+l-1]), cacheMin, rev);
      seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(syncmers[i]));
  }

  std::ofstream of(outFile);
  for (const auto& curPosition : seedmersIndex.positionMap) {
    of << curPosition.first << "\t"
       << std::get<0>(curPosition.second) << "\t"
       << std::get<1>(curPosition.second) << "\t"
       << std::get<2>(curPosition.second) << "\t"
       << seedmersIndex.seedmersMap.at(std::get<1>(curPosition.second)).size() << std::endl;
  }

}

void buildSeedmersHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index, mgsr::seedmers& seedmersIndex,
  std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, Tree *T, Node *node, globalCoords_t &globalCoords, CoordNavigator &navigator,
  std::vector<int> &scalarCoordToBlockId, std::vector<std::unordered_set<int>> &BlocksToSeeds, std::stringstream& seedmersOutStream,
  const std::vector<std::pair<size_t, size_t>>& blocksEndCoords, size_t& parentUngappedSeqSize, std::map<int32_t, int32_t>& parentCoordsIndex)
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

  std::unordered_set<int32_t> to_delete;
  //                     hash     beg     end      del
  std::vector<std::tuple<size_t, int32_t, int32_t, bool>> seedsChanges;
  std::vector<std::tuple<size_t, int32_t, int32_t, bool>> seedsChangesBackTrack;
  std::vector<std::pair<int32_t, int32_t>> seedGlobalEndCoorChanges;
  std::vector<std::pair<int32_t, int32_t>> seedGlobalEndCoorChangesBackTrack;

  //                     hash    beg      end      rev   del
  std::vector<std::tuple<size_t, int32_t, int32_t, bool, bool>> seedmersChangesBackTrack;
  std::vector<std::pair<int32_t, int32_t>> seedmersGlobalEndCoorChangesBackTrack;



  std::vector<tupleRange> merged;
  merged = expandAndMergeRanges(navigator, recompRanges, index.k(), data.blockExists);

  // Protobuf message for this node's mutations
  NodeSeedmerMutations *pb_node_mutations = index.add_per_node_mutations();
  pb_node_mutations->set_node_id(node->identifier);

  std::vector<int32_t> degap;
  std::vector<int32_t> regap;
  std::map<int32_t, int32_t> coordsIndex = parentCoordsIndex;
  std::map<int32_t, int32_t> parentCoordsIndexBackTrack;
  size_t ungappedSeqSize = parentUngappedSeqSize;
  size_t parentUngappedSeqSizeBackTrack;

  if (!merged.empty()) {
    coordsIndex.clear();
    ungappedSeqSize = 0;
    updateCoordsInfo(degap, regap, coordsIndex, data, blocksEndCoords, ungappedSeqSize, T, node);
    parentCoordsIndexBackTrack = std::move(parentCoordsIndex);
    parentUngappedSeqSizeBackTrack = parentUngappedSeqSize;
    parentCoordsIndex = coordsIndex;
    parentUngappedSeqSize = ungappedSeqSize;
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

            auto [seedHash, seedLegal] = getHash(seedMap[pos]);
            assert(seedLegal);
            seedsChanges.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], true));
            seedsChangesBackTrack.emplace_back(std::make_tuple(curSeeds.at(pos).second, pos, curSeeds.at(pos).first, false));
            to_delete.insert(pos);
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
            auto pos = str_i + startScalar;
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(pos);
            pb_mut->set_seq(seedMap[pos]);

            auto [seedHash, seedLegal] = getHash(seedMap[pos]);
            assert(seedLegal);
            seedsChanges.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], true));
            seedsChangesBackTrack.emplace_back(std::make_tuple(curSeeds.at(pos).second, pos, curSeeds.at(pos).first, false));
            to_delete.insert(pos);
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
          auto pos = str_i + startScalar;
          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(true);
          pb_mut->set_pos(pos);
          pb_mut->set_seq(seedMap[pos]);

          auto [seedHash, seedLegal] = getHash(seedMap[pos]);
          assert(seedLegal);
          seedsChanges.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], true));
          seedsChangesBackTrack.emplace_back(std::make_tuple(curSeeds.at(pos).second, pos, curSeeds.at(pos).first, false));
          to_delete.insert(pos);

          backtrack.push_back(std::make_pair(pos, seedMap[pos]));
          seedsToClear.push_back(pos);

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
          auto [seedHash, seedLegal] = getHash(kmer);
          if (seeding::is_syncmer(kmer, index.s(), false) && seedLegal)
          {
            // Is it still a seed?
            if (kmer != seedMap.find(str_i + startScalar)->second || to_delete.find(str_i + startScalar) != to_delete.end()) {
              // Is it a new seed?
              backtrack.push_back(std::make_pair(str_i + startScalar, seedMap[str_i + startScalar]));
              addSeeds.push_back(std::make_pair(str_i + startScalar, kmer));

              //lastDownstreamSeedPos = currCoord;
              auto pos = str_i + startScalar;
              if (kmer.size() == index.k())
              {
                SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
                pb_mut->set_is_deletion(false);
                pb_mut->set_pos(pos);
                pb_mut->set_seq(kmer);
                if (to_delete.find(pos) != to_delete.end()) {
                  assert(std::get<1>(seedsChanges.back()) == pos);
                  assert(std::get<1>(seedsChangesBackTrack.back()) == pos);
                  seedsChanges.pop_back();
                  seedsChangesBackTrack.pop_back();
                }
                seedsChanges.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], false));
                seedsChangesBackTrack.emplace_back(std::make_tuple(curSeeds.at(pos).second, pos, curSeeds.at(pos).first, false));
              }
            } else {
              auto pos = str_i + startScalar;
              if (to_delete.find(pos) == to_delete.end()) {
                seedGlobalEndCoorChanges.emplace_back(std::make_pair(pos, regap[degap[pos] + index.k() - 1]));
                seedGlobalEndCoorChangesBackTrack.emplace_back(std::make_pair(pos, curSeeds.at(pos).first));
              }
            }
          }
          else
          {
            auto pos = str_i + startScalar;
            backtrack.push_back(std::make_pair(pos, seedMap[pos]));

            // no longer a seed -> delete

            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(pos);
            pb_mut->set_seq(seedMap[pos]);


            auto [seedHash, seedLegal] = getHash(seedMap[pos]);
            assert(seedLegal);
            seedsChanges.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], true));
            seedsChangesBackTrack.emplace_back(std::make_tuple(curSeeds.at(pos).second, pos, curSeeds.at(pos).first, false));
            to_delete.insert(pos);

            seedsToClear.push_back(pos);
          }
        }
        else
        {
          //  not in seed map, could be a seed now
          auto [seedHash, seedLegal] = getHash(kmer);
          if (seeding::is_syncmer(kmer, index.s(), false) && seedLegal)
          {
            backtrack.push_back(std::make_pair(str_i + startScalar, ""));

            //std::string prevseedmer = lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1} ? seedMap[lastDownstreamSeedPos] : "";

            addSeeds.push_back(std::make_pair(str_i + startScalar, kmer));
            //seedMap[currCoord] = kmer + prevseedmer.substr(0, (index.j() - 1) * index.k());
            if (kmer.size() == index.k())
            {
              auto pos = str_i + startScalar;
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(pos);
              pb_mut->set_seq(kmer);

              seedsChanges.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], false));
              seedsChangesBackTrack.emplace_back(std::make_tuple(seedHash, pos, regap[degap[pos] + index.k() - 1], true));

            }
            //lastDownstreamSeedPos = currCoord;
          }
        }
      }
      
    } 
  } 

  // seedmers
  // std::cerr << "node " << node->identifier << std::endl;

  auto reversedSeedsChanges = std::ranges::reverse_view(std::views::all(seedsChanges));
  auto reversedSeedGlobalEndCoorChanges = std::ranges::reverse_view(std::views::all(seedGlobalEndCoorChanges));
  updateCurSeeds(curSeeds, reversedSeedsChanges, reversedSeedGlobalEndCoorChanges);

  seedmersOutStream << node->identifier << ":" << ungappedSeqSize << " ";
  
  std::unordered_set<int32_t> seenBegs;
  size_t  cacheReversedH;
  size_t  cacheForwardH;
  size_t  cacheMin;
  bool    rev;
  int     count;
  for (const auto& change : reversedSeedsChanges) {
    auto [seedHash, seedBeg, seedEnd, seedDel] = change;

    decltype(curSeeds.begin()) curSeedIt;
    if (seedDel) curSeedIt = curSeeds.lower_bound(seedBeg);
    else         curSeedIt = curSeeds.find(seedBeg);

    auto localStartSeed = curSeedIt;
    auto localEndSeed   = curSeedIt;
    int offset = 0;
    while (offset < index.j() - 1 && localStartSeed != curSeeds.begin()) {
      --localStartSeed;
      ++offset;
    }

    auto curAffectedSeed = localStartSeed;
    if (!seedDel) ++localEndSeed; 
    while (curAffectedSeed != localEndSeed) {
      if (seenBegs.find(curAffectedSeed->first) != seenBegs.end()) {
        ++curAffectedSeed;
        continue;
      }
      seenBegs.insert(curAffectedSeed->first);


      count = 0;
      cacheForwardH  = 0;
      cacheReversedH = 0;
      auto curLocalSeed = curAffectedSeed;
      while (count < index.j() && curLocalSeed != curSeeds.end()) {
        cacheForwardH  = (cacheForwardH << (2 * index.k())) + curLocalSeed->second.second;
        cacheReversedH = cacheReversedH + (curLocalSeed->second.second << (2 * index.k() * count));
        ++curLocalSeed;
        ++count;
      }
      --curLocalSeed;
      
      if (count != index.j()) {
        auto positionMapIt = seedmersIndex.positionMap.find(curAffectedSeed->first);
        if (positionMapIt != seedmersIndex.positionMap.end()) {
          auto begToErase  = curAffectedSeed->first;
          auto hashToErase = std::get<1>(positionMapIt->second);
          // std::vector<std::tuple<size_t, int32_t, int32_t, bool, bool>> seedmersChangesBackTrack;
          // std::vector<std::pair<int32_t, int32_t>> seedmersGlobalEndCoorChangesBackTrack;
          seedmersChangesBackTrack.emplace_back(std::make_tuple(
            std::get<1>(positionMapIt->second),
            positionMapIt->first, 
            std::get<0>(positionMapIt->second),
            std::get<2>(positionMapIt->second),
            false
          ));

          assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
          seedmersIndex.seedmersMap[hashToErase].erase(begToErase);
          if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
          assert(seedmersIndex.positionMap.erase(begToErase));
          seedmersIndex.positionMap.erase(begToErase);
          seedmersOutStream << begToErase << ":-:" << hashToErase << " ";
        }
        break;
      }

      if (cacheForwardH < cacheReversedH) {
        cacheMin = cacheForwardH;
        rev = false;
      } else if (cacheReversedH < cacheForwardH) {
        cacheMin = cacheReversedH;
        rev = true;
      } else {
        auto positionMapIt = seedmersIndex.positionMap.find(curAffectedSeed->first);
        if (positionMapIt != seedmersIndex.positionMap.end()) {
          auto begToErase  = curAffectedSeed->first;
          auto hashToErase = std::get<1>(positionMapIt->second);

          seedmersChangesBackTrack.emplace_back(std::make_tuple(
            std::get<1>(positionMapIt->second),
            positionMapIt->first, 
            std::get<0>(positionMapIt->second),
            std::get<2>(positionMapIt->second),
            false
          ));

          assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
          seedmersIndex.seedmersMap[hashToErase].erase(begToErase);
          if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
          assert(seedmersIndex.positionMap.erase(begToErase));
          seedmersIndex.positionMap.erase(begToErase);
          seedmersOutStream << begToErase << ":-:" << hashToErase << " ";
        }
        ++curAffectedSeed;
        continue;
      }
      
      auto potentialDel = seedmersIndex.positionMap.find(curAffectedSeed->first);
      if (potentialDel != seedmersIndex.positionMap.end()) {
        // backtrack replace
        seedmersChangesBackTrack.emplace_back(std::make_tuple(
          std::get<1>(potentialDel->second),
          potentialDel->first, 
          std::get<0>(potentialDel->second),
          std::get<2>(potentialDel->second),
          false
        ));

        seedmersIndex.seedmersMap.find(std::get<1>(potentialDel->second))->second.erase(potentialDel->first);
        if (seedmersIndex.seedmersMap[std::get<1>(potentialDel->second)].empty()) seedmersIndex.seedmersMap.erase(std::get<1>(potentialDel->second));
      } else {
        // backtrack delete
        seedmersChangesBackTrack.emplace_back(std::make_tuple(
          cacheMin,
          curAffectedSeed->first, 
          std::get<0>(potentialDel->second),
          std::get<2>(potentialDel->second),
          true
        ));
      }

      seedmersIndex.positionMap[curAffectedSeed->first] = std::make_tuple(curLocalSeed->second.first, cacheMin, rev);
      seedmersIndex.seedmersMap[cacheMin].insert(curAffectedSeed->first);



      seedmersOutStream << curAffectedSeed->first << ":+:"
                        << curLocalSeed->second.first << ","
                        << cacheMin << ","
                        << rev << " ";
      ++curAffectedSeed;
    }

    if (seedDel) {
      auto positionMapIt = seedmersIndex.positionMap.find(seedBeg);
      if (positionMapIt != seedmersIndex.positionMap.end()) {
        auto begToErase  = seedBeg;
        auto hashToErase = std::get<1>(positionMapIt->second);

        seedmersChangesBackTrack.emplace_back(std::make_tuple(
          std::get<1>(positionMapIt->second),
          positionMapIt->first, 
          std::get<0>(positionMapIt->second),
          std::get<2>(positionMapIt->second),
          false
        ));

        assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
        seedmersIndex.seedmersMap[hashToErase].erase(begToErase);
        if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
        assert(seedmersIndex.positionMap.erase(begToErase));
        seedmersIndex.positionMap.erase(begToErase);
        seedmersOutStream << begToErase << ":-:" << hashToErase << " ";
      }
    }
  }

  for (const auto& endCoorChange : reversedSeedGlobalEndCoorChanges) {
    auto curSeedIt = curSeeds.find(endCoorChange.first);
    int count = 0;
    while (count < index.j() - 1) {
      if (curSeedIt == curSeeds.begin()) break;
      ++count;
      --curSeedIt;
    }
    if (count != index.j() - 1) continue;
    auto curPositionIt = seedmersIndex.positionMap.find(curSeedIt->first);
    if (curPositionIt == seedmersIndex.positionMap.end()) continue;

    if (std::get<0>(curPositionIt->second) != endCoorChange.second) {
      seedmersGlobalEndCoorChangesBackTrack.emplace_back(std::make_pair(curPositionIt->first, std::get<0>(curPositionIt->second)));\
      curPositionIt->second = std::make_tuple(endCoorChange.second, std::get<1>(curPositionIt->second), std::get<2>(curPositionIt->second));
      seedmersOutStream << curPositionIt->first << ":e:" << endCoorChange.second << " ";
    }
  }

  seedmersOutStream << "c:";
  for (const auto& coord : coordsIndex) {
      seedmersOutStream << coord.first << "," << coord.second << ",";
  }
  seedmersOutStream << "\n";



  // std::string outFile = std::string("../dev/panmam/rsv/true_seedmers/") + node->identifier + std::string(".true.seedmers");
  // writeTrueSeedmers(tree::getStringAtNode(node, T, false), index.k(), index.s(), index.j(), outFile);

  // std::ofstream of(std::string("../dev/panmam/rsv/test_seedmers/") + node->identifier + std::string(".test.seedmers"));
  // for (const auto& curPosition : seedmersIndex.positionMap) {
  //   of << degapGlobalFromIndex(curPosition.first, coordsIndex) << "\t"
  //       << degapGlobalFromIndex(std::get<0>(curPosition.second), coordsIndex) << "\t"
  //       << std::get<1>(curPosition.second) << "\t"
  //       << std::get<2>(curPosition.second) << "\t"
  //       << seedmersIndex.seedmersMap.at(std::get<1>(curPosition.second)).size() << std::endl;
  // }


  
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
    buildSeedmersHelper(data, seedMap, index, seedmersIndex, curSeeds, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, seedmersOutStream, blocksEndCoords, parentUngappedSeqSize, parentCoordsIndex);
  }

  // undo parent coordsIndex and seqSize
  if (!parentCoordsIndexBackTrack.empty()) {
    parentCoordsIndex = std::move(parentCoordsIndexBackTrack);
    parentUngappedSeqSize = parentUngappedSeqSizeBackTrack;
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


  // undo curSeeds end coor changes
  for (const auto& back : seedGlobalEndCoorChangesBackTrack) {
    curSeeds.at(back.first).first = back.second;
  }
  // undo curSeeds changes
  for (const auto& back : seedsChangesBackTrack) {
    const auto& [hash, beg, end, del] = back;
    if (del) {
      assert(curSeeds.erase(beg));
      curSeeds.erase(beg);
    } else {
      curSeeds[beg] = std::make_pair(end, hash);
    }
  }


  // undo seedmers changes
  for (const auto& back : seedmersChangesBackTrack) {
    const auto& [hash, beg, end, rev, del] = back;
    if (del) {
      auto seedmersMapIt = seedmersIndex.seedmersMap.find(hash);
      seedmersIndex.positionMap.erase(beg);
      seedmersMapIt->second.erase(beg);
      if (seedmersMapIt->second.empty()) {
        seedmersIndex.seedmersMap.erase(hash);
      }
    } else {
      auto positionMapIt = seedmersIndex.positionMap.find(beg);
      if (positionMapIt != seedmersIndex.positionMap.end()) {
          auto [oend, ohash, orev] = positionMapIt->second;
          auto seedmersMapIt = seedmersIndex.seedmersMap.find(ohash);

          seedmersMapIt->second.erase(beg);
          if (seedmersMapIt->second.empty()) {
              seedmersIndex.seedmersMap.erase(ohash);
          }
      }
      seedmersIndex.positionMap[beg] = std::make_tuple(end, hash, rev);
      seedmersIndex.seedmersMap[hash].insert(beg);
    }
  }
  
  // undo seedmers end coor changes
  for (const auto& back : seedmersGlobalEndCoorChangesBackTrack) {
    const auto& [beg, end] = back;
    auto positionMapIt = seedmersIndex.positionMap.find(beg);
    auto [cend, chash, crev] = positionMapIt->second;
    seedmersIndex.positionMap[beg] = std::make_tuple(end, chash, crev);
  }

  /* Undo sequence mutations when backtracking */
  undoMutations(data, index, T, node, blockMutationInfo, mutationInfo);
}



void pmi::buildSeedmers(SeedmerIndex &index, Tree *T, int j, int k, int s, std::stringstream& seedmersOutStream) {
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;

    tree::setup(data, globalCoords, T);
    index.set_j(j);
    index.set_k(k);
    index.set_s(s);

    seedMap_t seedMap;
    CoordNavigator navigator(data.sequence);
    std::vector<int> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
    auto currCoord = tupleCoord_t{0,0,0};

    for(int i = 0; i < scalarCoordToBlockId.size(); i++){
        scalarCoordToBlockId[i] = currCoord.blockId;

        currCoord = navigator.increment(currCoord);
    }


    std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());

    std::vector<std::pair<size_t, size_t>> blocksEndCoords(globalCoords.size());
    for (size_t i = 0; i < globalCoords.size(); ++i) {
        size_t curBeg = globalCoords[i].first[0].second.empty() ? globalCoords[i].first[0].first : globalCoords[i].first[0].second[0];
        size_t curSize;
        if (i + 1 >= globalCoords.size()) {
            curSize = std::numeric_limits<size_t>::max();
        } else {
            curSize = globalCoords[i+1].first[0].second.empty() ? globalCoords[i+1].first[0].first - curBeg : globalCoords[i+1].first[0].second[0] - curBeg;
        }
        blocksEndCoords[i] = {curBeg, curSize};
    }

    mgsr::seedmers seedmersIndex;
    std::map<int32_t, std::pair<int32_t, size_t>> curSeeds;
    size_t parentUngappedSeqSize;
    std::map<int32_t, int32_t> parentCoordsIndex;
    seedmersOutStream << k << " " << s << " " << j << "\n";
    buildSeedmersHelper(data, seedMap, index, seedmersIndex, curSeeds, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, seedmersOutStream, blocksEndCoords, parentUngappedSeqSize, parentCoordsIndex);
}
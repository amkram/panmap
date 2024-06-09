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
  std::cout << "NODE: " << node->identifier << "\n";
  for (auto mutation : node->blockMutation)
  {

    int32_t primaryBlockId = mutation.primaryBlockId;
    int32_t secondaryBlockId = mutation.secondaryBlockId;
    bool type = mutation.blockMutInfo;
    bool inversion = mutation.inversion;

    std::cout << "issa mutante " << primaryBlockId << " " << secondaryBlockId << " " << type << " " << inversion <<"\n";

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

      tupleCoord_t endboundary = tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition};
      for(int j = 0; j <= len; j++){
        endboundary = navigator.increment(endboundary);
      }


      std::cout << "NORMal mutation " << len << "\n";
      std::cout << primaryBlockId << " "<< nucPosition << " "<< nucGapPosition << "\n";
      std::cout << endboundary.blockId << " "<< endboundary.nucPos << " "<< endboundary.nucGapPos << "\n";

      //Alex look here

      recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},    endboundary});
      
      //recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},  tupleCoord_t{primaryBlockId, nucPosition+len, -1}});


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
      
      std::cout << "Hell yeah " << len << "\n";
      std::cout << primaryBlockId << " "<< nucPosition << " "<< nucGapPosition << "\n";

      recompRanges.push_back({tupleCoord_t{primaryBlockId, nucPosition, nucGapPosition},
                              tupleCoord_t{primaryBlockId, nucPosition + len, nucGapPosition}});       //TODO inefficient and wrong?

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
  // std::cout << "expandLeft..." << std::endl;
  //  std::cout << "orig coord: (" << coord.blockId << ", " << coord.nucPos << ",
  //  "
  //            << coord.nucGapPos << ") to => ... " << std::endl;

  // std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

  while (count < neededNongap && coord > tupleCoord_t{0, 0, 0})
  {

    // std::cout << "count: " << count << " neededNongap: " << neededNongap << std::endl;
    //  std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

    if (!blockExists[coord.blockId].first)
    {
      // TODO jump down to prev block pleaesse

      // std::cout << "coord pre decrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
      coord = navigator.decrement(coord);
      // std::cout << "coord post decrememnto " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

      continue;
    }
    // std::cout << "isGap? " << std::endl;
    if (!navigator.isGap(coord))
    {
      // std::cout << "is not gap\n";
      count++;
    }

    // tupleCoord_t prev = coord;

    // std::cout << " from (" << coord.blockId << ", " << coord.nucPos << ", "
    // << coord.nucGapPos << ") to => ";
    // std::cout << "coord pre mad decrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
    coord = navigator.decrement(coord);
    // std::cout << "coord post mado decrememnto " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
  }
  return coord;
}

// Go downstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandRight(CoordNavigator &navigator, tupleCoord_t &coord,
                         int neededNongap, blockExists_t &blockExists)
{

  int count = 0;

  // std::cout << "ENTERING EXPANDRIGHT\n";

  // std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

  while (count < neededNongap && coord < tupleCoord_t{-1, -1, -1})
  {

    // std::cout << "count: " << count << " neededNongap: " << neededNongap << std::endl;
    // std::cout << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;

    if (!blockExists[coord.blockId].first) // TODO jump down to next block 
    {
      // std::cout << "block doesn't exist\n";
      

      // std::cout << "coord pre incrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
      coord = navigator.increment(coord);

      //coord.blockId += 1;
      //coord.nucPos = 0;
      //coord.nucGapPos = 0;

      // std::cout << "coord post instigation " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
      // std::cout << (coord < tupleCoord_t{-1, -1, -1}) << "\n";
      // std::cout << "we\n";
      continue;
    }

    if (!navigator.isGap(coord))
    {
      // std::cout << "is not gap\n";
      count++;
    }
    // tupleCoord_t prev = coord;

    // std::cout << "coord pre incrememnt " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
    coord = navigator.increment(coord);
    // std::cout << "coord post incrimination " << coord.blockId << ", " << coord.nucPos << ", " << coord.nucGapPos << std::endl;
  }
  // std::cout << "returning coord\n";
  //if(coord < tupleCoord_t{-1, -1, -1} ) {
    //coord = navigator.increment(coord);
  //}

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

  for (size_t i = 1; i < ranges.size(); ++i)
  {
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

    // if (expandedRange.start == tupleCoord_t{-1, -1, -1}) {
    //   expandedRange.start = tupleCoord_t{0, 0, 0};
    // }

    /*
    std::cout << "$%^*^% :" << expandedRange.start.blockId << ", "
               << expandedRange.start.nucPos << ", " << expandedRange.start.nucGapPos
               << " to " << current.stop.blockId << ", "
               << current.stop.nucPos << ", " << current.stop.nucGapPos
                << std::endl;
                */

    if (expandedRange.start <= current.stop)
    {
      // tmpStop = current.stop;
      /*
      std::cout << "comparino:" << expandedRange.stop.blockId << ", "
               << expandedRange.stop.nucPos << ", " << expandedRange.stop.nucGapPos
               << " to " << current.stop.blockId << ", "
               << current.stop.nucPos << ", " << current.stop.nucGapPos
                << std::endl;
                */

      current.stop = std::max(current.stop, expandedRange.stop);

      //std::cout << "ooosbijk:" <<  current.stop.blockId << ", "<< current.stop.nucPos << ", " << current.stop.nucGapPos  << std::endl;
                
    }
    else
    {
      merged.push_back(current);
      // tmpStop = current.stop;
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



  //std::cout << coord.blockId << " " << coord.nucPos << " " << coord.nucGapPos << "\n"; 

  if (coord == tupleCoord_t{-1, -1, -1})
  {

    //std::cout <<  globalCoords.back().first.back().first << " a\n";
    return globalCoords.back().first.back().first;
  }
  if (coord.nucGapPos >= 0)
  {
    //std::cout << globalCoords[coord.blockId].first[coord.nucPos].second[coord.nucGapPos] << " b\n";
    return globalCoords[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
  }
  //std::cout << globalCoords[coord.blockId].first[coord.nucPos].first << " c \n";
  
  return globalCoords[coord.blockId].first[coord.nucPos].first;
}

// Recursive function to build the seed index
void buildHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index,
                 Tree *T, Node *node, globalCoords_t &globalCoords,
                 CoordNavigator &navigator)
{
  
  // std::cout << "buildhelper in node: " << node->identifier << std::endl;

  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;
  tupleCoord_t start = {0, 0, 0};
  tupleCoord_t end = {-1, -1, -1};
  

  // First, a range is made marking the start -> end
  // of each block and nuc mutation. This is done while
  // applying mutations to the sequence object.
  std::vector<tupleRange> recompRanges;
  applyMutations(data, seedMap, blockMutationInfo, recompRanges, mutationInfo, T, node,
                 globalCoords, index, navigator);
  
  /*
  std::cout << "buildhelper in node: " << node->identifier << std::endl;
  std::cout << "blocks:" << std::endl;
  for (int i = 0; i < data.blockExists.size(); i++) {
    std::cout << i << ": " << data.blockExists[i].first << '\t';
  }
  std::cout << std::endl;
  */

  std::sort(recompRanges.begin(), recompRanges.end());

  bool is653 = false;


  //std::cout << "just to be sure, should be 33354 : " <<  tupleToScalarCoord(tupleCoord_t{12, 130, -1}, globalCoords) << "\n";
  //std::cout << "This is " << node->identifier << " and the 33354 seed is:" << seedMap[tupleCoord_t{-1, -1, -1}] << "\n";

   if(node->identifier == "OL776362.1")
  {
    is653 = true;
    
    //std::cout << "WE ARE AT NODE 3 OHHH YEAHHHHH....." << std::endl;
    //std::string seq = getNucleotideSequenceFromBlockCoordinates(start, end,data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords);
    //std::cout << "str: " << seq << std::endl;

    std::cout << "Here be the ranges, pre merged\n";
    for (auto &range : recompRanges) {
      std::cout << "range: " << range.start.blockId << ", " << range.start.nucPos << ", " << range.start.nucGapPos << " to " << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;
      std::cout << "ntpos: " << tupleToScalarCoord(range.start, globalCoords) << " to " << tupleToScalarCoord(range.stop, globalCoords) << std::endl;
    }
  }

  

  
  
  


   

  std::vector<tupleCoord_t> seedsToClear; // seeds to clear from seedMap
  std::vector<std::pair<tupleCoord_t, std::string>> addSeeds;
  std::vector<std::pair<tupleCoord_t, std::string>> backtrack;

  std::vector<tupleRange> merged = expandAndMergeRanges(navigator, recompRanges, index.k(), data.blockExists);

  
  if(is653){
    std::cout << "Here be the ranges, POST merged\n";
    for (auto &range : merged) {
      std::cout << "range: " << range.start.blockId << ", " << range.start.nucPos << ", " << range.start.nucGapPos << " to " << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;
      std::cout << "ntpos: " << tupleToScalarCoord(range.start, globalCoords) << " to " << tupleToScalarCoord(range.stop, globalCoords) << std::endl;
    }
  

  
    //exit(0);
  }


  //To see if downstream things are fine    //TODO delete
  //tupleRange fullseqRange = {start, end};
  //merged.clear();
  //merged.push_back(fullseqRange);




  // if(is3){
  //   merged.clear();
  //   merged.push_back(fullseqRange);
  // }

  if(is653){
    std::cout << "merged ranges: " << std::endl;
    for (auto &range : merged) {
      std::cout << "range: " << range.start.blockId << ", " << range.start.nucPos << ", " << range.start.nucGapPos << " to " << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;
      std::cout << "ntpos: " << tupleToScalarCoord(range.start, globalCoords) << " to " << tupleToScalarCoord(range.stop, globalCoords) << std::endl;
    }
  }
  // Protobuf message for this node's mutations
  NodeSeedmerMutations *pb_node_mutations = index.add_per_node_mutations();
  pb_node_mutations->set_node_id(node->identifier);

  std::cout << ">" << node->identifier << "\n";

  if (is653) { //TODO DELETE
    std::cout << "blocks:" << std::endl;
    for (int i = 0; i < data.blockExists.size(); i++) {
      std::cout << i << ": " << data.blockExists[i].first << '\t';
    }
    std::cout << std::endl;



  }




  // Seed re-processing
  for (auto &range : std::ranges::reverse_view(merged))
  {
    //if (range.stop == tupleCoord_t{-1, -1, -1}) {
    //  range.stop = tupleCoord_t{}
     //}
    // std::cout << "RANGE START: " << tupleToScalarCoord(range.start, globalCoords) << " RANGE STOP: " << tupleToScalarCoord(range.stop, globalCoords) << std::endl;
    std::string recomputeSeq = tree::getNucleotideSequenceFromBlockCoordinates(range.start, range.stop, data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords, navigator);


    if(is653){
      std::cout << " recomputeSeq " << recomputeSeq << "\n";
    }

    //std::cerr << ">" << node->identifier << "\n";
    //std::cout << "recomputeSeq " << recomputeSeq << "\n";
    //std::cout << "tuplerange: " << range.start.blockId << ", " << range.start.nucPos << ", " << range.start.nucGapPos << " to " << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;
    //std::cout << "range: " << tupleToScalarCoord(range.start, globalCoords) << " to " << tupleToScalarCoord(range.stop, globalCoords) << std::endl;

    
    // std::cout << "tuplerange: " << range.start.blockId << ", " << range.start.nucPos << ", " << range.start.nucGapPos << " to " << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;
    //std::cout << "range len: " << (tupleToScalarCoord(range.stop, globalCoords) - tupleToScalarCoord(range.start, globalCoords)) << std::endl;
    //std::cout << "recomp len: " << recomputeSeq.size() << std::endl;
    // std::cout << "first nt: " << recomputeSeq[0] << std::endl;
    // std::cout << "last nt: " << recomputeSeq[recomputeSeq.size() - 1] << std::endl;
    // Track the last downstream seed to stack k-mers into seedmers
    tupleCoord_t lastDownstreamSeedPos = range.stop;
    auto boundItr = seedMap.upper_bound(range.stop);
    if (boundItr == seedMap.end())
    {
      lastDownstreamSeedPos = tupleCoord_t{-1, -1, -1};
    }
    else
    {
      lastDownstreamSeedPos = boundItr->first;
    }


    
    bool atGlobalEnd = false;
    if (range.stop >= tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1})
    {
      atGlobalEnd = true;
      range.stop = tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1};
    }

    //std::cout << "NEWEND SHOULD PASS 12: "  << range.stop.blockId << ", " << range.stop.nucPos << ", " << range.stop.nucGapPos << std::endl;

    int32_t seen_non_gap = 0;
    int32_t str_i = tupleToScalarCoord(range.stop, globalCoords) - tupleToScalarCoord(range.start, globalCoords);


    //std::cout << "str_i should sbe 1 - this length " << str_i << " " << recomputeSeq.size() << "\n";
    //note to selves:: changed   currCoord >= range.start   to    currCoord > range.start
    for (auto currCoord = range.stop; currCoord >= range.start; currCoord = navigator.decrement(currCoord))
    {
      if (str_i < 0)
      {
        // std::cout << "break\n";
        break;
      }
      char nt = recomputeSeq[str_i];


      //std::cout << "We are at:\n";
      //std::cout << "Processing coord (" << currCoord.blockId << ", "  << currCoord.nucPos << ", " << currCoord.nucGapPos      << "): " << tupleToScalarCoord(currCoord, globalCoords)<< " with nt " << nt << std::endl;

      if(is653){
        
      }

      if (!data.blockExists[currCoord.blockId].first)
      {
        //std::cout << "block doesn't exist  " << currCoord.blockId << "\n";


        if (seedMap.find(currCoord) != seedMap.end())
        {
          //std::cout << "(+)->(-): " << seedMap[currCoord] << std::endl;
          // was a seed, no longer a seed due to block no exist -> delete
          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);
          if (seedMap[currCoord].size() == index.k())
          {
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
            pb_mut->set_seq(seedMap[currCoord]);
          }
        }else{
          //std::cout << "block dont exist but not in seedmap "  << std::endl;
        }
      }



      
      if (seen_non_gap < index.k())
      {

        //We need this for the end of the genome, but otherwise its wrong
        if (atGlobalEnd && seedMap.find(currCoord) != seedMap.end())  //Seed in map yet we dont have enough non-gaps for a seed, so we remove it
        {
          //std::cout << "(+)->(-): " << seedMap[currCoord] << std::endl;
          
          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);
          if (seedMap[currCoord].size() == index.k())
          {
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
            pb_mut->set_seq(seedMap[currCoord]);
          }
        }
        
        


        if (recomputeSeq[str_i] != '-' && recomputeSeq[str_i] != 'x')
        {
          // std::cout << "seen_non_gap++\n";
          seen_non_gap++;
        }
        // std::cout << "..skip, seen_non_gap < " << index.k() << std::endl;
        if (seen_non_gap < index.k() && str_i > 0) {
          str_i--;                           
          continue;  //
        }
      }
      //std::cout << "keep going, seen_non_gap >= " << index.k() << std::endl;
      //std::cout << "str[i=" << str_i << "] = " << recomputeSeq[str_i] << std::endl;

      if (!data.blockExists[currCoord.blockId].first) //TODO 
      {
      }
      else if (recomputeSeq[str_i] == '-' ||
               recomputeSeq[str_i] ==
                   'x')
      { // block does exist but seq is a gap
        // std::cout << "exists, is gap\n";
        if (seedMap.find(currCoord) != seedMap.end())
        {
          // std::cout << "in seedMap\n";
          // is a gap, no longer a seed -> delete
          // std::cout << "(+)->(_): " << seedMap[currCoord] << std::endl;
          SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
          pb_mut->set_is_deletion(true);
          pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
          pb_mut->set_seq(seedMap[currCoord]);

          backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
          seedsToClear.push_back(currCoord);

        } /* else: no seed, wasn't seed, no change */

      }
      else
      {
        // std::cout << "(+)->(+): " << seedMap[currCoord] << std::endl;
        // block exists and seq is not a gap at currCoord
        // std::cout << "block exists, not gap\n";
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
        // if (kmer == "AATCTTAGAACCAGA")
        // {
          if (is653) {
            std::cout << "__ KMER TIME __ \n";
            std::cout << "node: " << node->identifier << std::endl;
            std::cout << "kmer: " << kmer << std::endl;
            std::cout << "coord: " << currCoord.blockId << ", " << currCoord.nucPos << ", " << currCoord.nucGapPos << std::endl;
            std::cout << "str_i: " << str_i << std::endl;
            //std::cout << "recomputeSeq: " << recomputeSeq << std::endl;
            std::cout << "scalar coord: " << tupleToScalarCoord(currCoord, globalCoords) << std::endl;
            std::cout << "----------" << std::endl;
          }
        // }
        // std::cout << "kmer: " << kmer << std::endl;

        // std::cout << "kmer: " << kmer << std::endl;
        if (seedMap.find(currCoord) != seedMap.end())
        {
          // std::cout << "not gap + in seedMap\n";
          // non gap position and kmer is already a seed.
          std::string prevseedmer =
              lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1}
                  ? seedMap[lastDownstreamSeedPos]
                  : "";
          // std::cout << "prevseedmer: " << prevseedmer << std::endl;
          // std::cout << "bt " << c << " " << seedMap[c] << std::endl;

          if (seeding::is_syncmer(kmer, index.s(), false))
          {
            //std::cout << "IS A SEED " << kmer << "\n";
            // Is it still a seed?
            // std::cout << "still seed" << std::endl;
            backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));
            addSeeds.push_back(std::make_pair(currCoord, kmer));

            // std::cout << "set seedMap[" << c << "] = " << seedMap[c] <<
            // std::endl;
            lastDownstreamSeedPos = currCoord;
            if (kmer.size() == index.k())
            {
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
              pb_mut->set_seq(kmer);
            }
          }
          else
          {
            backtrack.push_back(std::make_pair(currCoord, seedMap[currCoord]));

            // no longer a seed -> delete
            // std::cout << "no longer seed, del " << seedMap[currCoord] << std::endl;
            // std::endl;
            SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
            pb_mut->set_is_deletion(true);
            pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
            pb_mut->set_seq(seedMap[currCoord]);

            seedsToClear.push_back(currCoord);
          }
        }
        else
        {
          //  not in seed map, could be a seed now
          if (seeding::is_syncmer(kmer, index.s(), false))
          {

            //std::cout << "IS A SEED " << kmer << "\n";

            backtrack.push_back(std::make_pair(currCoord, ""));
            std::string prevseedmer =
                lastDownstreamSeedPos != tupleCoord_t{-1, -1, -1}
                    ? seedMap[lastDownstreamSeedPos]
                    : "";
            addSeeds.push_back(std::make_pair(currCoord, kmer));
//            seedMap[currCoord] = kmer + prevseedmer.substr(0, (index.j() - 1) * index.k());
            if (kmer.size() == index.k())
            {
              SeedmerMutation *pb_mut = pb_node_mutations->add_mutations();
              pb_mut->set_is_deletion(false);
              pb_mut->set_pos(tupleToScalarCoord(currCoord, globalCoords));
              pb_mut->set_seq(kmer);
            }
            lastDownstreamSeedPos = currCoord;
          }
        }
      }
      str_i--;
    } //End of sequence loop
  } //End or ranges loop




  for (const auto &pos : seedsToClear)
  {
    if (seedMap.find(pos) != seedMap.end())
    {
      seedMap.erase(pos);
    }
  }


  for (const auto &seed : addSeeds) {
    seedMap[seed.first] = seed.second;
  }

  if (is653) {
    exit(0);
  }

  /* Recursive step */
  for (Node *child : node->children)
  {
 
    /*std::cout << "That was " << node->identifier << " and here are the seeds\n";
    for (const auto &seed : seedMap)
    {
      std::cout << seed.second << " " <<  tupleToScalarCoord(seed.first, globalCoords) << "\n";
    }*/

    
    buildHelper(data, seedMap, index, T, child, globalCoords, navigator);
  }


  // undo seed mutations
  for (const auto &back : backtrack)
  {
    if(back.second == ""){
      seedMap.erase(back.first);
    }else{
      seedMap[back.first] = back.second;
    }
    /*
    if (seedMap.find(back.first) != seedMap.end())
    {
      seedMap.erase(back.first);
    }
    else
    {
      seedMap[back.first] = back.second;
    }*/
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
  buildHelper(data, seedMap, index, T, T->root, globalCoords, navigator);
}

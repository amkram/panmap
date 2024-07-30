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

/*
if (pos == 0) {
    if (newValOneBack == '-') {
        gapRuns[pos + 1] = curIt->second - 1;
    }
    gapRuns.erase(curIt);
} else if (pos == globalCoords.back().first.back().first - 1) {
    if (oldValOneForward == '-') {
        auto prevIt = std::prev(curIt);
        --(prevIt->second);
    } else {
        gapRuns.erase(curIt);
    }
} else {
    if (newValOneBack != '-' && oldValOneForward == '-') {
        // At start of gap run
        gapRuns[pos + 1] = curIt->second - 1;
        gapRuns.erase(curIt);
    } else if (newValOneBack == '-' && oldValOneForward != '-') {
        // At end of gap run
        auto prevIt = std::prev(curIt);
        --(prevIt->second);
    } else if (newValOneBack == '-' && oldValOneForward == '-') {
        // In middle of gap run
        auto prevIt = std::prev(curIt);
        int leftGapRunLen = pos - prevIt->first;
        int rightGapRunLen = prevIt->second - leftGapRunLen - 1;
        prevIt->second = leftGapRunLen;
        gapRuns[pos + 1] = rightGapRunLen;
    } else {
        // Single base
        gapRuns.erase(curIt);
    }
}

for (const auto& change : nucChange) {
    int pos = change.pos;
    auto curit = gapRuns.lower_bound(pos);

    if (old[pos] == '-' && new[pos] != '-') {
        if (curit != gapRuns.end() && curit->first == pos) {
            // Exact match with the start of a gap run
            if (curit->second == 1) {
                // Single gap, remove the run
                gapRuns.erase(curit);
            } else {
                // Start of the run, reduce the run length
                int newStart = pos + 1;
                int newLength = curit->second - 1;
                gapRuns.erase(curit);
                gapRuns[newStart] = newLength;
            }
        } else {
            // Inside or at the end of a gap run
            auto previt = std::prev(curit);
            int start = previt->first;
            int length = previt->second;

            if (start <= pos && pos < start + length) {
                // The gap is within a run
                int leftGapRunLen = pos - start;
                int rightGapRunLen = length - leftGapRunLen - 1;

                if (leftGapRunLen == 0) {
                    // At the beginning of a gap run
                    if (rightGapRunLen > 0) {
                        // Update the start of the gap run
                        gapRuns[pos + 1] = rightGapRunLen;
                    }
                    gapRuns.erase(previt);
                } else if (rightGapRunLen == 0) {
                    // At the end of a gap run
                    previt->second = leftGapRunLen;
                } else {
                    // In the middle of a gap run, split the run
                    previt->second = leftGapRunLen;
                    gapRuns[pos + 1] = rightGapRunLen;
                }
            }
        }
    }
}
*/

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
                    globalCoords_t &globalCoords, ::capnp::List<Mutations>::Builder &indexedSeedMutations,
                    CoordNavigator &navigator, std::vector<std::pair<bool, std::pair<int, int>>> &gapRunUpdates,
                    std::vector<std::pair<bool, std::pair<int, int>>> &gapRunBacktracks)
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
  for (auto &mutation : mutationInfo) {
    int blockId = std::get<0>(mutation);
    int nucPos = std::get<2>(mutation);
    int nucGapPos = std::get<3>(mutation);
    char parChar = std::get<4>(mutation) == 'x' ? '-' : std::get<4>(mutation);
    char curChar = std::get<5>(mutation) == 'x' ? '-' : std::get<5>(mutation);
    int64_t scalar = tupleToScalarCoord(tupleCoord_t{blockId, nucPos, nucGapPos}, globalCoords);
    // std::cout << "parChar: " << parChar << " curChar: " << curChar << std::endl;
    if (parChar != '-' && curChar == '-') {
      // nuc to gap
        if (!gapRunUpdates.empty() && gapRunUpdates.back().first == true && gapRunUpdates.back().second.second + 1 == scalar) {
          gapRunUpdates.back().second.second += 1;
        }
        else {
          gapRunUpdates.emplace_back(true, std::make_pair(scalar, scalar)); 
        }
      } else if (parChar == '-' && curChar != '-') {
        // gap to nuc
        if (!gapRunUpdates.empty() && gapRunUpdates.back().first == false && gapRunUpdates.back().second.second + 1 == scalar) {
          gapRunUpdates.back().second.second += 1;
        } else {
          gapRunUpdates.emplace_back(false, std::make_pair(scalar, scalar));
        }
      }
    }
  }

void undoMutations(mutableTreeData &data, ::capnp::List<Mutations>::Builder &indexedSeedMutations, Tree *T,
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


void updateGapMap(std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<bool, std::pair<int, int>>>& updates, std::vector<std::pair<bool, std::pair<int, int>>>& backtrack) {
  for (const auto& update : updates) {
    bool toGap = update.first;
    int start = update.second.first;
    int end = update.second.second;
    auto rightIt = gapMap.upper_bound(start);
    auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

    if (toGap) {
      // add gap range
      if (gapMap.empty()) {
        gapMap[start] = end;
        backtrack.emplace_back(true, std::make_pair(start, end));
        return;
      }
      
      decltype(rightIt) curIt;

      // curIt starts outside of any range
      if (leftIt == gapMap.end() || (rightIt == gapMap.end() && start > leftIt->second) || (leftIt != gapMap.end() && start > leftIt->second && rightIt != gapMap.end() && end < rightIt->first)) {
        if (leftIt != gapMap.end() && start == leftIt->second + 1) {
          // 1 base after left range and merge with left
          curIt = leftIt;
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          curIt->second = end;
        } else {
          // insert new range
          auto tmpIt = gapMap.emplace(start, end);
          curIt = tmpIt.first;
          backtrack.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        }
      } else {
        curIt = leftIt;
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        curIt->second = end;
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
          gapMap.erase(tmpIt);
        } else if (nextIt->first <= end + 1) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          curIt->second = nextIt->second;
          backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
          gapMap.erase(nextIt);
          break;
        } else {
          break;
        }
      }
    } else {
      // remove gap range
      if (gapMap.empty() || (leftIt == gapMap.end() && end < leftIt->first) || (rightIt == gapMap.end() && start > leftIt->second)) {
        return;
      }

      decltype(rightIt) curIt;
      decltype(rightIt) nextIt;
      if (leftIt == gapMap.end() || (leftIt != gapMap.end() && start > leftIt->second && rightIt != gapMap.end() && start < rightIt->first)) {
        // curIt starts outside of any range
        curIt = rightIt;

        // ends within the curIt range
        if (end <= curIt->second) {
          if (end == curIt->second) {
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            gapMap.erase(curIt);
          } else {
            gapMap[end+1] = curIt->second;
            backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            gapMap.erase(curIt);
          }
          return;
        } else {
          nextIt = std::next(curIt);
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          gapMap.erase(curIt);
        }
        
      } else {
        // curIt starts inside of a range
        curIt = leftIt;
        
        if (end <= curIt->second) {
          // contained in the curIt range
          if (start == curIt->first && end == curIt->second) {
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            gapMap.erase(curIt);
          } else if (start == curIt->first) {
            gapMap[end + 1] = curIt->second;
            backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            gapMap.erase(curIt);
          } else if (end == curIt->second) {
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            curIt->second = start - 1;
          } else {
            gapMap[end + 1] = curIt->second;
            backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            curIt->second = start - 1;
          }
          return;
        } else {
          if (start == curIt->first) {
            nextIt = std::next(curIt);
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            gapMap.erase(curIt);
          } else {
            backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
            curIt->second = start - 1;
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
          gapMap.erase(tmpIt);
        } else {
          gapMap[end + 1] = nextIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, nextIt->second));
          backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
          gapMap.erase(nextIt);
          break;
        }
      }
    }
  }
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
std::vector<std::pair<int64_t, std::optional<int64_t>>> encodeDelta(
    const std::map<int64_t, int64_t>& map1,
    const std::map<int64_t, int64_t>& map2) 
{
    std::vector<std::pair<int64_t, std::optional<int64_t>>> delta;

    // Iterate through the first map
    for (const auto& [key, value] : map1) {
        auto it = map2.find(key);
        if (it != map2.end()) {
            // Key exists in both maps
            if (it->second != value) {
                // Values are different (update scenario)
                delta.emplace_back(key, value);
            }
        } else {
            // Key only exists in map1 (deletion scenario)
            delta.emplace_back(key, std::nullopt); // Use nullopt to indicate deletion
        }
    }

    // Iterate through the second map to find additions
    for (const auto& [key, value] : map2) {
        if (map1.find(key) == map1.end()) {
            // Key only exists in map2 (addition scenario)
            delta.emplace_back(key, value);
        }
    }

    return delta;
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

std::vector<std::tuple<std::string, int, int>>
extractSeedmers(const std::string &seq, const int k, const int s,
                const bool open) {
  std::vector<std::tuple<std::string, int, int>> syncmers;
  std::unordered_map<int32_t, int32_t> degap;
  int64_t pos = 0;
  std::string ungapped = "";
  for (int64_t i = 0; i < seq.size(); i++) {
    char c = seq[i];
    degap[pos] = i;
    if (c != '-' && c != 'x') {
      ungapped += c;
      pos++;
    }
  }
  if (ungapped.size() < k + 1) {
    return syncmers;
  }
  for (size_t i = 0; i < ungapped.size() - k + 1; ++i) {
    std::string kmer = ungapped.substr(i, k);
    if (seeding::is_syncmer(kmer, s, open)) {
      syncmers.emplace_back(std::make_tuple(kmer, degap[i], degap[i + k - 1]));
    }
  }

  return syncmers;
}

bool debug = false;
bool gappity = true;
// Recursive function to build the seed index
void buildHelper(mutableTreeData &data, seedMap_t &seedMap, ::capnp::List<Mutations>::Builder &indexedSeedMutations, ::capnp::List<Deltas>::Builder &indexedGapMutations,
                 int32_t &seedK, int32_t &seedS,
                 Tree *T, Node *node, globalCoords_t &globalCoords,
                 CoordNavigator &navigator, std::vector<int> &scalarCoordToBlockId, std::vector<std::unordered_set<int>> &BlocksToSeeds, std::vector<int> &BlockSizes,
                 int64_t &dfsIndex, posWidth &width, std::map<int64_t, int64_t> &gapMap)
{

  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;

  // First, a range is made marking the start -> end
  // of each block and nuc mutation. This is done while
  // applying mutations to the sequence object.

  //          std::pair<nucToGap, std::pair<beg, end>>
  std::vector<std::pair<bool, std::pair<int, int>>> gapRunUpdates;
  //          std::pair<del,  std::pair<beg, end>>
  std::vector<std::pair<bool, std::pair<int, int>>> gapRunBacktracks;
  std::vector<tupleRange> recompRanges;


  std::map<int64_t, int64_t> parentGapMap = gapMap;

  applyMutations(data, seedMap, blockMutationInfo, recompRanges, mutationInfo, T, node,
                 globalCoords, indexedSeedMutations, navigator, gapRunUpdates, gapRunBacktracks);
  
  
  // std::cout << "gapRunUpdates: " << gapRunUpdates.size() << std::endl;
  // for (auto &update : gapRunUpdates) {
  //   std::cout << update.first << " " << update.second.first << " " << update.second.second << std::endl;
  // }
  
  // apply gapRunUpdates to GapMap
    std::vector<std::pair<bool, std::pair<int, int>>> gapRunOffBlocksBacktracks;
    std::vector<std::pair<bool, std::pair<int, int>>> gapRunOffBlocksUpdates;
    std::map<int64_t, int64_t> coordIndex;

  if (gappity) {
    updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks);

    for (int i = 0; i < data.blockExists.size(); i++) {
      if (!data.blockExists[i].first) {
        int64_t start = globalCoords[i].first[0].first;
        int64_t end = globalCoords[i].first.back().first;
        gapRunOffBlocksUpdates.emplace_back(true, std::make_pair(start, end));
      }
    }
    updateGapMap(gapMap, gapRunOffBlocksUpdates, gapRunOffBlocksBacktracks);
  //       --AAACC----TTA
  //       (2,2) (11,6)
    //       first po of a nuc run, number of gaps before this position
    //       localCoord = globalcoord - coordIndex.some_bound->second

    int64_t totalGapSize = 0;
    if (gapMap.empty() || gapMap.begin()->first > 0) {
      coordIndex[0] == totalGapSize;
    }
    for (auto &gap : gapMap) {
      int64_t gapStart = gap.first;
      int64_t gapEnd = gap.second;
      int64_t gapSize = gapEnd - gapStart + 1;
      totalGapSize += gapSize;
      coordIndex[gapEnd+1] = totalGapSize;
    }

    // should this be here?
    for (auto it = gapRunOffBlocksBacktracks.rbegin(); it != gapRunOffBlocksBacktracks.rend(); ++it) {
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

  std::vector<int> seedsToClear; // seeds to clear from seedMap
  std::vector<std::pair<int, std::string>> addSeeds;
  std::vector<std::pair<int, std::string>> backtrack;

  std::vector<tupleRange> merged;

  merged = expandAndMergeRanges(navigator, recompRanges, seedK, data.blockExists, data.blockStrand, globalCoords);
  
  tupleCoord_t start = {0, 0, 0};
  if(data.sequence[0].first[0].second.size() == 0){
    start.nucGapPos = -1;
  }
  tupleCoord_t end = tupleCoord_t{data.sequence.size() - 1, data.sequence.back().first.size() - 1, -1};

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
      }
    }

    //Loop through seeds that now start as gaps and delete them
    for(int i = 0; i < gaps.size(); i++){
      if (seedMap.find(gaps[i]) != seedMap.end()){
        backtrack.push_back(std::make_pair(gaps[i], seedMap[gaps[i]]));
        seedsToClear.push_back(gaps[i]);
      }
    }

    //Loop through sequence and build up seeds as we go
    if(seq.size() >= seedK){

      //Check first k bp for smers
      std::string min_s = seq.substr(0,seedS);
      int first_min_coord = 0;
      int last_min_coord = 0;
      int Ns = seq[0] == 'N';
      
      for(int i = 1; i < seedK - seedS + 1; i++){
        std::string smer = seq.substr(i,seedS);
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
      for(int i = seedK - seedS + 1; i < seedK ; i++){
        if(seq[i]=='N')
          Ns++;
      }

      //Check the rest for smers and kmer seeds
      for(int i = seedK; i < seq.size() + 1; i++) {
        
        //Processing kmer starting at i - k in seq
        bool inMap = (seedMap.find(coords[i - seedK]) != seedMap.end());
        bool isSeed = (first_min_coord == i - seedK || last_min_coord == i - seedS) && Ns <= seedK/2;

        if(!inMap && isSeed){
          //Add seed
          std::string kmer = seq.substr(i - seedK, seedK);

          backtrack.push_back(std::make_pair(coords[i - seedK], ""));
          addSeeds.push_back(std::make_pair(coords[i - seedK], kmer));
          

        }else if(inMap && !isSeed){
          //Remove Seed

          backtrack.push_back(std::make_pair(coords[i - seedK], seedMap[coords[i - seedK]]));
          seedsToClear.push_back(coords[i - seedK]);

        }else if(inMap && isSeed){
          //Changed Seed
          std::string kmer = seq.substr(i - seedK, seedK);
          if(kmer != seedMap[coords[i - seedK]]){
            backtrack.push_back(std::make_pair(coords[i - seedK], seedMap[coords[i - seedK]]));
            addSeeds.push_back(std::make_pair(coords[i - seedK], kmer));
          }

        }

        //Updating smers for kmer starting at i - seedK + 1
        if(i < seq.size()){
          if(seq[i] == 'N')
            Ns++;
          if(seq[i - seedK] == 'N')
            Ns--;

          if(first_min_coord == i - seedK){
            // Were losing the lowest smer, Re search for lowest
            min_s = seq.substr(i - seedK + 1,seedS);
            first_min_coord = i - seedK + 1;

            for(int j = 1; j < seedK - seedS + 1; j++){
              std::string smer = seq.substr(i - seedK + 1 + j,seedS);
              if(smer < min_s){
                min_s = smer;
                first_min_coord = i - seedK + 1 + j;
                last_min_coord = first_min_coord;
              }else if(smer == min_s){
                last_min_coord = i - seedK + 1 + j;
              }
            }
          }else{
            //Test new smer to see if its the lowest
            std::string smer = seq.substr(i - seedS + 1,seedS);
            if(smer < min_s){
              min_s = smer;
              first_min_coord = i - seedS + 1;
              last_min_coord = i - seedS + 1;
            }else if(smer == min_s){
              last_min_coord = i - seedS + 1;
            }
          }
        }
      }
    }

    //If our range reaches the end of the genome, remove seeds that aren't long enough at the end
    if(atGlobalEnd && (int)coords.size() - (int)seedK + 1 >= 0){
      for(int i = (int)coords.size() - (int)seedK + 1; i < (int)coords.size(); i++){
        if (seedMap.find(coords[i]) != seedMap.end()){

          backtrack.push_back(std::make_pair(coords[i], seedMap[coords[i]]));
          seedsToClear.push_back(coords[i]);
        }
      }
    }
  }

  if (debug) {
    std::cout << "- real -" << std::endl;
    std::string thisNode = tree::getStringAtNode(node, T, true);
    auto syncmers = seeding::syncmerize(thisNode, seedK, seedS, false, true, 0);
    for (int i = 0; i < syncmers.size(); i++) {
      std::cout << syncmers[i].pos << " " << syncmers[i].seq << std::endl;
    }
    std::cout << "- MY seedMap -" << std::endl;
    for (auto &seed : seedMap) {
      std::cout << seed.first << " " << seed.second << std::endl;
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
  
  std::vector<int32_t> capnpDelNormal;
  std::vector<std::pair<int32_t, std::bitset<64>>> capnpDelOffset;
  size_t i = 0;
  size_t n = seedsToClear.size();

  while (i < n) {
    int pos = seedsToClear[i];
    std::bitset<64> downstream(0);
    i++; 
    if (i >= n) {
      break; 
    }

    int curr = seedsToClear[i]; 
    std::vector<int32_t> downstreamPos;

    while (curr < pos + 64 && i < n) {
      downstream[curr - pos] = 1; 
      downstreamPos.emplace_back(curr); 
      i++; 

      if (i < n) {
        curr = seedsToClear[i]; 
      }
    }

    if (downstream.none()) {
      for (int j = 0; j < downstreamPos.size(); j++) {
        capnpDelNormal.emplace_back(downstreamPos[j]); 
      }
    } else {
      capnpDelOffset.emplace_back(std::make_pair(pos, downstream)); 
    }
  }

  std::sort(addSeeds.begin(), addSeeds.end());
  std::vector<int32_t> capnpAddNormal;
  std::vector<std::pair<int32_t, std::bitset<64>>> capnpAddOffset;
  i = 0;
  n = addSeeds.size();

  while (i < n) {
    int64_t pos = addSeeds[i].first;
    std::bitset<64> downstream(0);
    i++; 
    if (i >= n) {
      break; 
    }

    int64_t curr = addSeeds[i].first; 
    std::vector<int32_t> downstreamPos;

    while (curr < pos + 64 && i < n) {
      downstream[curr - pos] = 1;
      downstreamPos.emplace_back(curr); 
      i++; 

      if (i < n) {
        curr = addSeeds[i].first; 
      }
    }

    if (downstream.none()) {
      for (int32_t j = 0; j < downstreamPos.size(); j++) {
        capnpAddNormal.emplace_back(downstreamPos[j]); 
      }
    } else {
      capnpAddOffset.emplace_back(std::make_pair(pos, downstream)); 
    }
  }
  std::vector<std::pair<int64_t, std::optional<int64_t>>> delta;
  if (gappity) {
    delta = encodeDelta(parentGapMap, gapMap);
    ::capnp::List<Delta>::Builder deltas = indexedGapMutations[dfsIndex].initChanges(delta.size());
    for (int32_t i = 0; i < delta.size(); i++) {
      if (delta[i].second.has_value()) {
        deltas[i].initOptional().setValue(delta[i].second.value());
      } else {
        deltas[i].initOptional().setVoid();
      }
    }
  }

  ::capnp::List<InsertionWithOffset>::Builder insertionsOffset = indexedSeedMutations[dfsIndex].initInsertionsWithOffset(capnpAddOffset.size());
  ::capnp::List<DeletionWithOffset>::Builder deletionsOffset = indexedSeedMutations[dfsIndex].initDeletionsWithOffset(capnpDelOffset.size());
  ::capnp::List<Insertion>::Builder insertions = indexedSeedMutations[dfsIndex].initInsertions(capnpAddNormal.size());
  ::capnp::List<Deletion>::Builder deletions = indexedSeedMutations[dfsIndex].initDeletions(capnpDelNormal.size());
  
  for (int32_t i = 0; i < capnpDelOffset.size(); i++) {
    deletionsOffset[i].setBitset(capnpDelOffset[i].second.to_ullong());
    deletionsOffset[i].initPos().setPos32(capnpDelOffset[i].first);
  }
  for (int32_t i = 0; i < capnpAddOffset.size(); i++) {
    insertionsOffset[i].setBitset(capnpAddOffset[i].second.to_ullong());
    insertionsOffset[i].initPos().setPos32(capnpAddOffset[i].first);
  }
  for (int32_t i = 0; i < capnpDelNormal.size(); i++) {
    deletions[i].initPos().setPos32(capnpDelNormal[i]);
  }
  for (int32_t i = 0; i < capnpAddNormal.size(); i++) {
    insertions[i].initPos().setPos32(capnpAddNormal[i]);
  }

  delta.clear();
  capnpAddOffset.clear();
  capnpAddNormal.clear();
  capnpDelOffset.clear();
  capnpDelNormal.clear();
  seedsToClear.clear();
  addSeeds.clear();
  merged.clear();
  recompRanges.clear();
  coordIndex.clear();
  parentGapMap.clear();
  gapRunUpdates.clear();

  dfsIndex++;
  /* Recursive step */
  for (Node *child : node->children) {
    
    buildHelper(data, seedMap, indexedSeedMutations, indexedGapMutations, seedK, seedS, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, dfsIndex, width, gapMap);
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
  if (gappity) {
    // undo gapMap updates
    for (auto it = gapRunBacktracks.rbegin(); it != gapRunBacktracks.rend(); ++it) {
      const auto& [del, range] = *it;
      if (del) {
        gapMap.erase(range.first);
      } else {
        gapMap[range.first] = range.second;
      }
    }
  }

  /* Undo sequence mutations when backtracking */
  undoMutations(data, indexedSeedMutations, T, node, blockMutationInfo, mutationInfo, globalCoords);
}

/* implementation */
void pmi::build(Tree *T, Index::Builder &index)
{
  // Setup for seed indexing
  tree::mutableTreeData data;
  tree::globalCoords_t globalCoords;

  tree::setup(data, globalCoords, T);
  
  int32_t k = index.getK();
  int32_t s = index.getS();

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
  posWidth width = globalCoords.back().first.back().first < 4294967296
      ? posWidth::pos32 
      : posWidth::pos64;

  switch(width){
    case posWidth::pos32:
      index.setWidth(32);
      break;
    case posWidth::pos64:
      index.setWidth(64);
      break;
  }

  /* Recursive traversal of tree to build the index */
  std::map<int64_t, int64_t> gapMap;

  tupleCoord_t coord = {0,0,data.sequence.front().second.size() == 0 ? -1 : 0};
  auto curIt = gapMap.end();

  while (coord < tupleCoord_t{-1, -1, -1})
  {

    char c = coord.nucGapPos == -1 ? data.sequence[coord.blockId].first[coord.nucPos].first : data.sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
    int64_t scalar = tupleToScalarCoord(coord, globalCoords);

    if (c == '-' || c == 'x') {
      if (!gapMap.empty() && curIt->second + 1 == scalar) {
        ++curIt->second;
      } else {
        auto tmpIt = gapMap.emplace(scalar, scalar);
        curIt = tmpIt.first;
      }
    }
    if (!data.blockExists[coord.blockId].first)
    {
      if(coord.blockId == navigator.sequence.size() - 1){
        break;
      } else {
        if( ! data.blockStrand[coord.blockId + 1].first){
          coord = tupleCoord_t{coord.blockId + 1, navigator.sequence[coord.blockId + 1].first.size() - 1, -1};
        }else{
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
    coord = navigator.newincrement(coord, data.blockStrand);
  }
  
  if(coord.blockId != -1 && !data.blockExists[coord.blockId].first){
    coord = navigator.newdecrement(coord, data.blockStrand);
  }

  ::capnp::List<Mutations>::Builder indexedSeedMutations = index.initPerNodeSeedMutations(T->allNodes.size());
  ::capnp::List<Deltas>::Builder indexedGapMutations = index.initPerNodeGapMutations(T->allNodes.size());
  int64_t dfsIndex = 0; 
  buildHelper(data, seedMap, indexedSeedMutations, indexedGapMutations, k, s, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, dfsIndex, width, gapMap);
}
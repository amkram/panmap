#include "pmi.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"
#include "seed_annotated_tree.hpp"
#include "conversion.hpp"

#include <iomanip>
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
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_sort.h>
#include <boost/icl/interval_map.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Dense>
// #include <omp.h>
#include <numeric>
extern "C" {
#include <bwa/bwa.h>
}

using namespace boost::icl;
using namespace mgsr;
namespace fs = boost::filesystem;

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

void makeCoordIndex(std::map<int64_t, int64_t>& degapCoordIndex, std::map<int64_t, int64_t>& regapCoordIndex, const std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<int64_t, int64_t>>& blockRanges) {
  int64_t totalGapSize = 0;
  if (gapMap.empty() || gapMap.begin()->first > 0) {
    degapCoordIndex[0] == totalGapSize;
    regapCoordIndex[0] == totalGapSize;
  }
  for (auto &gap : gapMap) {
    int64_t gapStart = gap.first;
    int64_t gapEnd = gap.second;
    int64_t gapSize = gapEnd - gapStart + 1;
    if (gapEnd == blockRanges.back().second) break;
    totalGapSize += gapSize;
    degapCoordIndex[gapEnd+1] = totalGapSize;
    regapCoordIndex[gapEnd+1-totalGapSize] = totalGapSize;
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
// extractSyncmers(const std::string &seq, const int k, const int s,
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
extractSyncmers(const std::string &seq, const int k, const int s, const int t, const bool open) {
  std::vector<std::tuple<std::string, size_t, bool, int>> seedmers;
  for (int64_t i = 0; i < seq.size() - k + 1; ++i) {
    std::string kmer = seq.substr(i, k);
    auto [hash, isReverse, isSyncmer] = seeding::is_syncmer_rollingHash(kmer, s, open, t);
    if (isSyncmer) {
      seedmers.emplace_back(std::make_tuple(kmer, hash, isReverse, i));
    }
  }
  return seedmers;
}

std::vector<std::tuple<std::string, size_t, bool, int, int>>
extractSeedmers(const std::string& seq, const int k, const int s, const int t, const int l, const bool open) {
  std::vector<std::tuple<std::string, size_t, bool, int>> syncmers;
  for (int64_t i = 0; i < seq.size() - k + 1; ++i) {
    std::string kmer = seq.substr(i, k);
    auto [hash, isReverse, isSyncmer] = seeding::is_syncmer_rollingHash(kmer, s, open, t);
    if (isSyncmer) {
      syncmers.emplace_back(std::make_tuple(kmer, hash, isReverse, i));
    }
  }
  
  std::vector<std::tuple<std::string, size_t, bool, int, int>> seedmers;
  size_t forwardRolledHash = 0;
  size_t reverseRolledHash = 0;
  for (size_t i = 0; i < syncmers.size() - l + 1; ++i) {
    std::string seedmer = "";
    std::string seedmerToHash = "";
    std::string seedmerToHashReverse = "";
    for (size_t j = 0; j < l; ++j) {
      seedmer += std::get<0>(syncmers[i+j]);
      seedmerToHash += std::get<2>(syncmers[i+j]) ? seeding::revcomp(std::get<0>(syncmers[i+j])) : std::get<0>(syncmers[i+j]);
      seedmerToHashReverse += std::get<2>(syncmers[i+l-j-1]) ? seeding::revcomp(std::get<0>(syncmers[i+l-j-1])) : std::get<0>(syncmers[i+l-j-1]);
    }
    if (i == 0) {
      for (size_t j = 0; j < l; ++j) {
        forwardRolledHash = rol(forwardRolledHash, k) ^ std::get<1>(syncmers[i+j]);
        reverseRolledHash = rol(reverseRolledHash, k) ^ std::get<1>(syncmers[i+l-j-1]);
      }
    } else {
      forwardRolledHash = rol(forwardRolledHash, k) ^ rol(std::get<1>(syncmers[i-1]), k * l) ^ std::get<1>(syncmers[i+l-1]);
      reverseRolledHash = ror(reverseRolledHash, k) ^ ror(std::get<1>(syncmers[i-1]), k)     ^ rol(std::get<1>(syncmers[i+l-1]), k * (l-1));
    }


    if (forwardRolledHash != seeding::hashSeq(seedmerToHash).first) {
      std::cout << "Forward hash mismatch" << std::endl;
      exit(1);
    }
    if (reverseRolledHash != seeding::hashSeq(seedmerToHashReverse).first) {
      std::cout << "Reverse hash mismatch" << std::endl;
      exit(1);
    }
    if (forwardRolledHash < reverseRolledHash) {
      seedmers.emplace_back(std::make_tuple(seedmer, forwardRolledHash, false, std::get<3>(syncmers[i]), std::get<3>(syncmers[i+l-1]) + k - 1));
    } else if (reverseRolledHash < forwardRolledHash) {
      seedmers.emplace_back(std::make_tuple(seedmer, reverseRolledHash, true, std::get<3>(syncmers[i]), std::get<3>(syncmers[i+l-1]) + k - 1));
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

bool debug = false;
// Recursive function to build the seed index
template <typename SeedMutationsType, typename GapMutationsType>
void buildOrPlace(Step method, mutableTreeData& data, std::vector<std::pair<std::string, float>> &placementScores, std::vector<std::optional<std::string>>& onSeeds, std::vector<std::optional<seeding::onSeedsHash>>& onSeedsHash, SeedMutationsType& perNodeSeedMutations_Index, GapMutationsType& perNodeGapMutations_Index, int seedK, int seedS, int seedT, bool open, int seedL, Tree* T, Node* node, globalCoords_t& globalCoords, CoordNavigator& navigator, std::vector<int64_t>& scalarCoordToBlockId, std::vector<std::unordered_set<int>>& BlocksToSeeds, std::vector<int>& BlockSizes, std::vector<std::pair<int64_t, int64_t>>& blockRanges, int64_t& dfsIndex, std::map<int64_t, int64_t>& gapMap, std::unordered_set<int64_t>& inverseBlockIds, int64_t jacNumer, int64_t jacDenom,  std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts, std::unordered_map<std::string, int64_t> &dfsIndexes) {
  // Variables needed for both build and place
  std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>> seedChanges;
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBlocksBacktracks;
  std::vector<std::pair<bool, int64_t>> inverseBlockIdsBacktrack;
  std::map<int64_t, int64_t> degapCoordIndex;
  std::map<int64_t, int64_t> regapCoordIndex;
  std::vector<tupleRange> recompRanges;
  blockExists_t oldBlockExists = data.blockExists;
  blockStrand_t oldBlockStrand = data.blockStrand;

  applyMutations(data, blockMutationInfo, recompRanges,  mutationInfo, T, node, globalCoords, navigator, blockRanges, gapRunUpdates, gapRunBacktracks, oldBlockExists, oldBlockStrand, method == Step::PLACE, inverseBlockIds, inverseBlockIdsBacktrack);

  if (method == Step::BUILD) {
    //                    erase           beg      end
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapMapUpdates;

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

    makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, blockRanges);
    
    for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
      const auto& [del, range] = *it;
      if (del) {
        gapMap.erase(range.first);
      } else {
        gapMap[range.first] = range.second;
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
            auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[pos].value();
            seedChanges.emplace_back(std::make_tuple(pos,
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt,
                                                oldEndPos,
                                                std::nullopt));
          }
        }
      }

      //Loop through seeds that now start as gaps and delete them
      for(int i = 0; i < gaps.size(); i++){
        if (onSeedsHash[gaps[i]].has_value()){
          // std::string oldSeed = gaps[i] < onSeeds.size() && onSeeds[gaps[i]].has_value() ? onSeeds[gaps[i]].value() : "";
          auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[gaps[i]].value();
          seedChanges.emplace_back(std::make_tuple(gaps[i], 
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt,
                                                oldEndPos,
                                                std::nullopt)); 
        }
      }

      //Loop through sequence and build up seeds as we go
      if(seq.size() >= seedK){
        // vector of (hash, isReverse, isSeed, startPos)
        std::vector<std::tuple<size_t, bool, bool, int64_t>> kmers = seeding::rollingSyncmers(seq, seedK, seedS, open, seedT, true);

        for (int64_t i = 0; i < kmers.size(); ++i) {
          const auto& [hash, isReverse, isSeed, startPos] = kmers[i];
          bool inMap = (onSeedsHash[coords[startPos]].has_value());

          if(!inMap && isSeed){
            //Add seed
            // std::string newSeed = seq.substr(i, seedK);
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                false, // old seed off
                                                true, // new seed on
                                                std::nullopt,
                                                hash,
                                                std::nullopt,
                                                isReverse,
                                                std::nullopt,
                                                coords[i+seedK-1]));
            
          }else if(inMap && !isSeed){
            //Remove Seed
            // std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[coords[i]].value();
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt,
                                                oldEndPos,
                                                std::nullopt));

          }else if(inMap && isSeed){
            // replace seed
            // std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[coords[i]].value();
            // std::string newSeed = seq.substr(i, seedK);
            size_t newSeed = hash;
            if(newSeed != oldSeed || oldIsReverse != isReverse || oldEndPos != coords[i+seedK-1]){
              seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                true, // new seed on
                                                oldSeed,
                                                newSeed,
                                                oldIsReverse,
                                                isReverse,
                                                oldEndPos,
                                                coords[i+seedK-1]));
            }
          }
        }
      }

      //If our range reaches the end of the genome, remove seeds that aren't long enough at the end
      if(atGlobalEnd && (int)coords.size() - (int)seedK + 1 >= 0){  
        for(int i = (int)coords.size() - (int)seedK + 1; i < (int)coords.size(); i++){
          if (onSeedsHash[coords[i]].has_value()) {
            // std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[coords[i]].value();
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt,
                                                oldIsReverse,
                                                std::nullopt,
                                                oldEndPos,
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
      const auto& [curPos, curOldVal, curNewVal, curOldSeed, curNewSeed, curOldIsReverse, curIsReverse, curOldEndPos, curEndPos] = seedChanges[seedChangeIndex];
      basePositions.push_back(curPos);
      std::vector<std::pair<uint64_t, uint64_t>> masks;

      // Generate ternary numbers for the current seed position
      while (seedChangeIndex >= 0) {
        const auto& [maskPos, maskOldVal, maskNewVal, maskOldSeed, maskNewSeed, maskOldIsReverse, maskIsReverse, maskOldEndPos, maskEndPos] = seedChanges[seedChangeIndex];
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

    makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, blockRanges);

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
            auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[pos - k].value();
            seedChanges.emplace_back(std::make_tuple(pos - k, true, false, oldSeed, std::nullopt, oldIsReverse, std::nullopt, oldEndPos, std::nullopt));
          } else if (ternaryNumber == 2) {
            // Handle insertion/change
            auto [newSeed, newEndPos] = seed_annotated_tree::getSeedAt(pos - k, T, seedK, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap, blockRanges);
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
              auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[pos - k].value();
              seedChanges.emplace_back(std::make_tuple(pos - k, true, true, oldSeed, newSeedHash, oldIsReverse, newIsReverse, oldEndPos, newEndPos));
            } else { // off -> on
              seedChanges.emplace_back(std::make_tuple(pos - k, false, true, std::nullopt, newSeedHash, std::nullopt, newIsReverse, std::nullopt, newEndPos));
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
    const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse, oldEndPos, newEndPos] = p;

    if (oldVal && newVal) { // seed at same pos changed
      onSeedsHash[pos].value().hash = newSeed.value();
      onSeedsHash[pos].value().endPos = newEndPos.value();
      onSeedsHash[pos].value().isReverse = newIsReverse.value();
    } else if (oldVal && !newVal) { // seed on to off
      if (onSeedsHash[pos].has_value() && pos < onSeedsHash.size()) {
        onSeedsHash[pos].reset();
      }
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].erase(pos);
    } else if (!oldVal && newVal) { // seed off to on
      onSeedsHash[pos] = {newSeed.value(), newEndPos.value(), newIsReverse.value()};
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].insert(pos);
    } 

    /**
     * JACCARD: (a ∩ b) / (a ∪ b)
     * WEIGHTED_JACCARD: sum(min(f(a), f(b))) / sum(max(f(a), f(b)))
     * COSINE: sum(f(a) * f(b)) / sqrt(sum(f(a)^2) * sum(f(b)^2))
    */

    if (method == Step::PLACE) {
      auto oldSeedVal = oldSeed.value_or(0);
      auto newSeedVal = newSeed.value_or(0);
      if (oldVal && newVal) { // seed at same pos changed
        if (oldSeedVal && readSeedCounts.find(oldSeedVal) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer -= 1;
        } else {
          //case 2: seed not in reads
          jacDenom -= 1;
        }
        if (newSeedVal && readSeedCounts.find(newSeedVal) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer += 1;
        } else {
          //case 2: seed not in reads
          jacDenom += 1;
        }
      } else if (oldSeedVal && oldVal && !newVal) { // seed on to off
        if (readSeedCounts.find(oldSeedVal) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer -= 1;
        } else {
          //case 2: seed not in reads
          jacDenom -= 1;
        }
      } else if (newSeedVal && !oldVal && newVal) { // seed off to on
        if (readSeedCounts.find(newSeedVal) != readSeedCounts.end()) {
          //case 1: seed is in reads
          jacNumer += 1;
        } else {
          //case 2: seed not in reads
          jacDenom += 1;
        }
      }
      float jac = jacDenom == 0 ? 0 : jacNumer / (float)jacDenom;
    }
  }
  if (method == Step::PLACE) {
    float jac = jacDenom == 0 ? 0 : jacNumer / (float)jacDenom;
    placementScores.emplace_back(std::make_pair(node->identifier, jac));
  }

  if (debug) {
    // print out seeds at node
    if (method == Step::PLACE) {
      std::cout << node->identifier << " place syncmers: ";
      for (int i = 0; i < onSeedsHash.size(); i++) {
        if (onSeedsHash[i].has_value()) {
          std::cout << mgsr::degapGlobal(i, degapCoordIndex) << ":" << onSeedsHash[i].value().hash << "|" << onSeedsHash[i].value().isReverse << " ";
        }
      }
      std::cout << std::endl;
    } else {
      std::cout << node->identifier << " build syncmers: ";
      for (int i = 0; i < onSeedsHash.size(); i++) {
        if (onSeedsHash[i].has_value()) {
          std::cout << mgsr::degapGlobal(i, degapCoordIndex) << ":" << onSeedsHash[i].value().hash << "|" << onSeedsHash[i].value().isReverse << " ";
        }
      }
      std::cout << std::endl;
    }



    // if (method == Step::PLACE) std::cout << node->identifier << " place coordIndex: ";
    // else                       std::cout << node->identifier << " build coordIndex: ";

    // for (const auto& [pos, gaps] : coordIndex) {
    //   std::cout << pos << ":" << gaps << " ";
    // }
    // std::cout << std::endl;

    if (method == Step::BUILD) {
      std::cout << node->identifier << " true syncmers: ";
      auto seq = seed_annotated_tree::getStringAtNode(node, T, false);
      if (seq.size() > 0) {
        auto syncmers = extractSyncmers(seq, seedK, seedS, seedT, open);
        for (const auto &[kmer, hash, isReverse, startPos] : syncmers) {
          std::cout << startPos << ":" << hash << "|" << isReverse << " ";
        }
      }
      std::cout << std::endl;
    }
  }



  /* Recursive step */
  dfsIndexes[node->identifier] = dfsIndex;
  dfsIndex++;
  for (Node *child : node->children) {
    buildOrPlace(
      method, data, placementScores, onSeeds, onSeedsHash, perNodeSeedMutations_Index, perNodeGapMutations_Index, seedK, seedS, seedT, open, seedL, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, jacNumer, jacDenom, readSeedCounts, dfsIndexes
    );
  }


  for (const auto &p : seedChanges) {
    const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse, oldEndPos, newEndPos] = p;
    if (oldVal && newVal) { // UNDO seed at same pos changed
      onSeedsHash[pos].value().hash = oldSeed.value();
      onSeedsHash[pos].value().endPos = oldEndPos.value();
      onSeedsHash[pos].value().isReverse = oldIsReverse.value();
    } else if (oldVal && !newVal) { // seed on to off
      onSeedsHash[pos] = {oldSeed.value(), oldEndPos.value(), oldIsReverse.value()};
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].insert(pos);
    } else if (!oldVal && newVal) { // UNDO seed off to on
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

void getBestNodeSeeds(std::vector<Node *> &rpath, mutableTreeData& data, std::vector<std::optional<std::string>>& onSeeds, std::vector<std::optional<seeding::onSeedsHash>>& onSeedsHash, ::capnp::List<SeedMutations>::Reader &perNodeSeedMutations_Index, ::capnp::List<GapMutations>::Reader &perNodeGapMutations_Index, int seedK, int seedS, int seedT, bool open, int seedL, Tree* T, globalCoords_t& globalCoords, CoordNavigator& navigator, std::vector<int64_t>& scalarCoordToBlockId, std::vector<std::unordered_set<int>>& BlocksToSeeds, std::vector<int>& BlockSizes, std::vector<std::pair<int64_t, int64_t>>& blockRanges, std::map<int64_t, int64_t>& gapMap, std::unordered_set<int64_t>& inverseBlockIds, std::unordered_map<std::string, int64_t> &dfsIndexes) {
  
  for (Node *node : rpath) {

    std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>> seedChanges;
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

    applyMutations(data, blockMutationInfo, recompRanges,  mutationInfo, T, node, globalCoords, navigator, blockRanges, gapRunUpdates, gapRunBacktracks, oldBlockExists, oldBlockStrand, Step::PLACE, inverseBlockIds, inverseBlockIdsBacktrack);

      ::capnp::List<MapDelta>::Reader gapMutationsList;
      gapMutationsList = perNodeGapMutations_Index[dfsIndexes[node->identifier]].getDeltas();
      
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

      // makeCoordIndex(coordIndex, gapMap, blockRanges);

      auto currBasePositions = perNodeSeedMutations_Index[dfsIndexes[node->identifier]].getBasePositions();
      auto currPerPosMasks = perNodeSeedMutations_Index[dfsIndexes[node->identifier]].getPerPosMasks();
      std::vector<std::vector<int8_t>> masks;
      for (int i = 0; i < currBasePositions.size(); ++i) {
          int64_t pos = currBasePositions[i];
          uint64_t tritMask = currPerPosMasks[i];
          for (int k = 0; k < 32; ++k) {
            uint8_t ternaryNumber = (tritMask >> (k * 2)) & 0x3;
            if (ternaryNumber == 1) { // on -> off
              // Handle deletion
              // seedChanges.emplace_back(std::make_tuple(pos - k, true, false, onSeeds[pos - k], std::nullopt));
              auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[pos - k].value();
              seedChanges.emplace_back(std::make_tuple(pos - k, true, false, oldSeed, std::nullopt, oldIsReverse, std::nullopt, oldEndPos, std::nullopt));
            } else if (ternaryNumber == 2) {
              // Handle insertion/change
              auto [newSeed, newEndPos] = seed_annotated_tree::getSeedAt(pos - k, T, seedK, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap, blockRanges);
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
                auto [oldSeed, oldEndPos, oldIsReverse] = onSeedsHash[pos - k].value();
                seedChanges.emplace_back(std::make_tuple(pos - k, true, true, oldSeed, newSeedHash, oldIsReverse, newIsReverse, oldEndPos, newEndPos));
              } else { // off -> on
                seedChanges.emplace_back(std::make_tuple(pos - k, false, true, std::nullopt, newSeedHash, std::nullopt, newIsReverse, std::nullopt, newEndPos));
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
    

    for (const auto &p : seedChanges)
    {
      const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse, oldEndPos, newEndPos] = p;

      if (oldVal && newVal) { // seed at same pos changed
        onSeedsHash[pos].value().hash = newSeed.value();
        onSeedsHash[pos].value().endPos = newEndPos.value();
        onSeedsHash[pos].value().isReverse = newIsReverse.value();
      } else if (oldVal && !newVal) { // seed on to off
        if (onSeedsHash[pos].has_value() && pos < onSeedsHash.size()) {
          onSeedsHash[pos].reset();
        }
        int blockId = scalarCoordToBlockId[pos];
        BlocksToSeeds[blockId].erase(pos);
      } else if (!oldVal && newVal) { // seed off to on
        onSeedsHash[pos] = {newSeed.value(), newEndPos.value(), newIsReverse.value()};
        int blockId = scalarCoordToBlockId[pos];
        BlocksToSeeds[blockId].insert(pos);
      }
    }
  }
}

void fillDfsIndexes(Tree *T, Node *node, int64_t &dfsIndex, std::unordered_map<std::string, int64_t> &dfsIndexes) {
  dfsIndexes[node->identifier] = dfsIndex;
  dfsIndex++;
  for (Node *child : node->children) {
    fillDfsIndexes(T, child, dfsIndex, dfsIndexes);
  }
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
  int32_t t = index.getT();
  bool open = index.getOpen();
  int32_t l = index.getL();
  
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

  for (int64_t i = 0; i < blockRanges.size(); i++) {
    int64_t start = globalCoords[i].first[0].second.empty() ? tupleToScalarCoord({i, 0, -1}, globalCoords) : tupleToScalarCoord({i, 0, 0}, globalCoords);
    int64_t end = tupleToScalarCoord({i, globalCoords[i].first.size() - 1, -1}, globalCoords);
    blockRanges[i] = std::make_pair(start, end);
    if (data.blockStrand[i].first) inverseBlockIds.insert(i);
  }

  std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());

  ::capnp::List<SeedMutations>::Builder perNodeSeedMutations_Builder = index.initPerNodeSeedMutations(T->allNodes.size());
  ::capnp::List<GapMutations>::Builder perNodeGapMutations_Builder = index.initPerNodeGapMutations(T->allNodes.size());

  int64_t dfsIndex = 0;
  

  // std::vector<std::optional<std::string>> onSeedsString(globalCoords.back().first.back().first + 1, std::nullopt);
  std::vector<std::optional<std::string>> onSeedsString;
  std::vector<std::optional<seeding::onSeedsHash>> onSeedsHash(globalCoords.back().first.back().first + 1, std::nullopt);
  int64_t jacNumer = 0;
  int64_t jacDenom = 0;
  std::unordered_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
  std::vector<std::pair<std::string, float>> placementScoresDummy;
  std::unordered_map<std::string, int64_t> dfsIndexes;
  buildOrPlace(
    Step::BUILD, data, placementScoresDummy, onSeedsString, onSeedsHash, perNodeSeedMutations_Builder, perNodeGapMutations_Builder, k, s, t, open, l, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, jacNumer, jacDenom, readSeedCounts, dfsIndexes
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

void seedsFromFastq(const int32_t& k, const int32_t& s, const int32_t& t, const bool& open, const int32_t& l, std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::vector<std::vector<seed>> &readSeeds,  const std::string &fastqPath1, const std::string &fastqPath2) {
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
            readSequences.push_back(reverseComplement(seq->seq.s));
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
      std::vector<seeding::seed> curReadSeeds;
      for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : rollingSyncmers(readSequences[i], k, s, open, t, false)) {
        if (!isSyncmer) continue;
        //curReadSeeds.emplace_back(seed{kmerHash, startPos, -1, isReverse, startPos + k - 1});
        curReadSeeds.emplace_back(seed{kmerHash, startPos + k - 1, -1, isReverse, 0});
        if (readSeedCounts.find(kmerHash) == readSeedCounts.end()) readSeedCounts[kmerHash] = std::make_pair(0, 0);
        if (isReverse) ++readSeedCounts[kmerHash].second;
        else           ++readSeedCounts[kmerHash].first;
      }
      readSeeds.push_back(std::move(curReadSeeds));
    }
}

int extractPosition(char *line) {
    int fieldCount = 0;
    char *ptr = line;
    char *fieldStart = ptr;

    while (*ptr != '\0') {
        if (*ptr == '\t') {
            fieldCount++;
            if (fieldCount == 4) {
                // We've found the end of the 4th field
                *ptr = '\0';  // Temporarily terminate the 4th field string

                // Convert the 4th field to an integer
                int position = atoi(fieldStart);

                *ptr = '\t';  // Restore the original character
                return position;
            }
            // Move to the start of the next field
            fieldStart = ptr + 1;
        }
        ptr++;
    }

    // Check if the line ends exactly after the 4th field
    if (fieldCount == 3) {
        // Line ends after the 4th field
        int position = atoi(fieldStart);
        return position;
    }

    // If we reach here, the line doesn't have at least 4 fields
    throw std::runtime_error("Line does not have at least 4 fields.");
}

void prepareAndRunBwa(
  const std::vector<std::string>& idx_args, const std::vector<std::string>& aln_args1, const std::vector<std::string>& aln_args2, const std::vector<std::string>& samaln_args, const std::string& reads1Path, const std::string& reads2Path, std::vector<std::pair<int, char*>> &samAlignmentPairs, std::vector<std::string> &samHeaders
  ) {
  std::vector<std::vector<char>> idx_argv_buf;
  std::vector<char*> idx_argv;

  std::cerr << "Constructing idx_argv..." << std::endl;
  for (const auto& arg : idx_args) {
    std::cerr << "Processing argument: " << arg << std::endl;
    idx_argv_buf.push_back(std::vector<char>(arg.begin(), arg.end()));
    idx_argv_buf.back().push_back('\0');  // Ensure null termination
    idx_argv.push_back(idx_argv_buf.back().data());  // Push the C-string pointer
  }
  idx_argv.push_back(nullptr);  // Null-terminate the argv array

  std::cerr << "Constructed idx_argv:" << std::endl;
  for (size_t i = 0; i < idx_argv.size(); ++i) {
    if (idx_argv[i] != nullptr)
      std::cerr << "idx_argv[" << i << "]: " << idx_argv[i] << std::endl;
    else
      std::cerr << "idx_argv[" << i << "]: (null)" << std::endl;
  }

  std::vector<std::vector<char>> aln_argv_buf1;
  std::vector<char*> aln_argv1;
  std::cerr << "Constructing aln_argv1..." << std::endl;
  for (const auto& arg : aln_args1) {
    std::cerr << "Processing argument: " << arg << std::endl;
    aln_argv_buf1.push_back(std::vector<char>(arg.begin(), arg.end()));
    aln_argv_buf1.back().push_back('\0');
    aln_argv1.push_back(aln_argv_buf1.back().data());
  }
  aln_argv1.push_back(nullptr);

  std::cerr << "Constructed aln_argv1:" << std::endl;
  for (size_t i = 0; i < aln_argv1.size(); ++i) {
    if (aln_argv1[i] != nullptr)
      std::cerr << "aln_argv1[" << i << "]: " << aln_argv1[i] << std::endl;
    else
      std::cerr << "aln_argv1[" << i << "]: (null)" << std::endl;
  }

  std::vector<std::vector<char>> aln_argv_buf2;
  std::vector<char*> aln_argv2;
  if (aln_args2.size() > 0) {
    std::cerr << "Constructing aln_argv2..." << std::endl;
    for (const auto& arg : aln_args2) {
      std::cerr << "Processing argument: " << arg << std::endl;
      aln_argv_buf2.push_back(std::vector<char>(arg.begin(), arg.end()));
      aln_argv_buf2.back().push_back('\0');
      aln_argv2.push_back(aln_argv_buf2.back().data());
    }
    aln_argv2.push_back(nullptr);

    std::cerr << "Constructed aln_argv2:" << std::endl;
    for (size_t i = 0; i < aln_argv2.size(); ++i) {
      if (aln_argv2[i] != nullptr)
        std::cerr << "aln_argv2[" << i << "]: " << aln_argv2[i] << std::endl;
      else
        std::cerr << "aln_argv2[" << i << "]: (null)" << std::endl;
    }
  }

  std::vector<std::vector<char>> samaln_argv_buf;
  std::vector<char*> samaln_argv;
  std::cerr << "Constructing samaln_argv..." << std::endl;
  for (const auto& arg : samaln_args) {
    std::cerr << "Processing argument: " << arg << std::endl;
    samaln_argv_buf.push_back(std::vector<char>(arg.begin(), arg.end()));
    samaln_argv_buf.back().push_back('\0');
    samaln_argv.push_back(samaln_argv_buf.back().data());
  }
  samaln_argv.push_back(nullptr);

  std::cerr << "Constructed samaln_argv:" << std::endl;
  for (size_t i = 0; i < samaln_argv.size(); ++i) {
    if (samaln_argv[i] != nullptr)
      std::cerr << "samaln_argv[" << i << "]: " << samaln_argv[i] << std::endl;
    else
      std::cerr << "samaln_argv[" << i << "]: (null)" << std::endl;
  }

  std::cerr << "About to call run_bwa with idx_argv..." << std::endl;

  try {
    std::cerr << "Running run_bwa with idx_argv..." << std::endl;
    optind = 1;
    if (run_bwa(idx_argv.size() - 1, idx_argv.data()) != 0) {
      throw std::runtime_error("BWA index failed");
    }
    std::cerr << "Finished run_bwa with idx_argv." << std::endl;

    std::cerr << "Running run_bwa with aln_argv1..." << std::endl;
    optind = 1;
    int stdoutFd = dup(fileno(stdout));
    if (run_bwa(aln_argv1.size() - 1, aln_argv1.data()) != 0) {
      throw std::runtime_error("BWA aln failed");
    }
    fflush(stdout);
    dup2(stdoutFd, fileno(stdout));
    close(stdoutFd);
    std::cerr << "Finished run_bwa with aln_argv1." << std::endl;

    if (aln_args2.size() > 0) {
      std::cerr << "Running run_bwa with aln_argv2..." << std::endl;
      optind = 1;
      int stdoutFd2 = dup(fileno(stdout));
      if (run_bwa(aln_argv2.size() - 1, aln_argv2.data()) != 0) {
        throw std::runtime_error("BWA aln failed");
      }
      fflush(stdout);
      dup2(stdoutFd2, fileno(stdout));
      close(stdoutFd2);
      std::cerr << "Finished run_bwa with aln_argv2." << std::endl;
    }

    std::cerr << "Running run_bwa with samaln_argv..." << std::endl;
    optind = 1;
    FILE *tempFile = std::tmpfile();
    if (!tempFile) {
      std::cerr << "Failed to create temporary file for capturing bwa SAM alignments." << std::endl;
      return;
    }
    int stdoutFd3 = dup(fileno(stdout));
    if (stdoutFd3 == -1) {
      std::cerr << "Failed to duplicate stdout file descriptor for capturing bwa SAM alignments." << std::endl;
      return;
    }
    if (dup2(fileno(tempFile), fileno(stdout)) == -1) {
      std::cerr << "Failed to redirect stdout to temporary file for capturing bwa SAM alignments." << std::endl;
      return;
    }
    if (run_bwa(samaln_argv.size() - 1, samaln_argv.data()) != 0) {
      throw std::runtime_error("BWA samaln failed");
    }
    fflush(stdout);
    dup2(stdoutFd3, fileno(stdout));
    close(stdoutFd3);
    rewind(tempFile); 
    char buffer[512];
    size_t i = 0;
    while (fgets(buffer, sizeof(buffer), tempFile)) {
      char *line = new char[strlen(buffer) + 1];
      std::strcpy(line, buffer);
      if (aln_args2.size() > 0) {
        if (i < 4) {
          samHeaders.push_back(line);
        } else {
          size_t len = strlen(line);
          if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';
          int pos = extractPosition(line);
          samAlignmentPairs.push_back(std::make_pair(pos, line));
        }
      } else {
        if (i < 3) {
          samHeaders.push_back(line);
        } else {
          size_t len = strlen(line);
          if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';
          int pos = extractPosition(line);
          samAlignmentPairs.push_back(std::make_pair(pos, line));
        }
      }
      ++i;
    }
    fclose(tempFile);
    std::cerr << "Finished run_bwa with samaln_argv." << std::endl;
    
  } catch (...) {
    std::cerr << "run_bwa() caused an exception!" << std::endl;
  }
}

void pmi::place(
  Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path,
  seed_annotated_tree::mutationMatrices &mutMat, std::string prefix,std::string refFileName,
  std::string samFileName, std::string bamFileName, std::string mpileupFileName, std::string vcfFileName,
  std::string aligner, const std::string& refNode, const bool& save_jaccard)
{
    // Setup for seed indexing
    seed_annotated_tree::mutableTreeData data;

    seed_annotated_tree::globalCoords_t globalCoords;

    seed_annotated_tree::setup(data, globalCoords, T);
    seed_annotated_tree::mutableTreeData bestNodeData = data;

    CoordNavigator navigator(data.sequence);
    CoordNavigator bestNodeNavigator(bestNodeData.sequence);

    std::vector<int> BlockSizes(data.sequence.size(),0);
    std::vector<std::pair<int64_t, int64_t>> blockRanges(data.blockExists.size());
    std::unordered_set<int64_t> inverseBlockIds;


    int32_t k = index.getK();
    int32_t s = index.getS();
    int32_t t = index.getT();
    bool open = index.getOpen();
    int32_t l = index.getL();
    
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
    std::vector<int> bestNodeBlockSizes = BlockSizes;
    std::vector<std::pair<int64_t, int64_t>> bestNodeBlockRanges = blockRanges;
    std::unordered_set<int64_t> bestNodeInverseBlockIds = inverseBlockIds;
    std::map<int64_t, int64_t> bestNodeGapMap = gapMap;
    std::vector<std::unordered_set<int>> bestNodeBlocksToSeeds = BlocksToSeeds;
    ::capnp::List<GapMutations>::Reader perNodeGapMutations_Reader = index.getPerNodeGapMutations();
    ::capnp::List<SeedMutations>::Reader perNodeSeedMutations_Reader= index.getPerNodeSeedMutations();
  
    
    // std::vector<std::optional<std::string>> onSeedsString(globalCoords.back().first.back().first + 1, std::nullopt);
    std::vector<std::optional<std::string>> onSeedsString;
    std::vector<std::optional<seeding::onSeedsHash>> onSeedsHash(globalCoords.back().first.back().first + 1, std::nullopt);
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::vector<std::vector<seed>> readSeeds;
    std::unordered_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
    seedsFromFastq(k, s, t, open, l, readSeedCounts, readSequences, readQuals, readNames, readSeeds, reads1Path, reads2Path);

    bool pairedEndReads = reads2Path.size();


    int64_t jacNumer = 0;
    int64_t jacDenom = 0;    
    std::vector<std::pair<std::string, float>> placementScores;
    std::string bestNodeId;
    Node *bestNode;
    Node *curr;
    int64_t dfsIndex = 0;
    std::unordered_map<std::string, int64_t> dfsIndexes;
    
    if (refNode.empty()) {
      for (const auto& count : readSeedCounts) {
        jacDenom += count.second.first + count.second.second;
        // std::cout << count.first << " " << count.second.first << " " << count.second.second << std::endl;
      }

      // buildOrPlace<decltype(perNodeSeedMutations_Reader), decltype(perNodeGapMutations_Reader)>(
      buildOrPlace<decltype(perNodeSeedMutations_Reader), decltype(perNodeGapMutations_Reader)>(
        Step::PLACE, data, placementScores, onSeedsString, onSeedsHash, perNodeSeedMutations_Reader, perNodeGapMutations_Reader, k, s, t, open, l, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, jacNumer, jacDenom, readSeedCounts, dfsIndexes
      );

      std::sort(placementScores.begin(), placementScores.end(), [](auto &left, auto &right) {
        return left.second > right.second;
      });

      if (save_jaccard) {
        std::ofstream jaccardFile(prefix + ".jaccard.txt");
        for (const auto& score : placementScores) {
          jaccardFile << score.first << "\t" << score.second << "\n";
        }
        jaccardFile.close();
      }

      bestNodeId = placementScores[0].first;
      bestNode = T->allNodes[bestNodeId];
      curr = bestNode;
    } else {
      fillDfsIndexes(T, T->root, dfsIndex, dfsIndexes);
      bestNodeId = refNode;
      bestNode = T->allNodes[refNode];
      curr = bestNode;
    }

    std::vector<Node *> rpath;
    while (curr != nullptr) {
      rpath.push_back(curr);
      curr = curr->parent;
    }
    std::reverse(rpath.begin(), rpath.end());
    std::vector<std::optional<std::string>> bestNodeOnSeedsString;
    std::vector<std::optional<seeding::onSeedsHash>> bestNodeOnSeedsHash(globalCoords.back().first.back().first + 1, std::nullopt);
    getBestNodeSeeds(rpath, bestNodeData, bestNodeOnSeedsString, bestNodeOnSeedsHash, perNodeSeedMutations_Reader, perNodeGapMutations_Reader, k, s, t, open, l, T, globalCoords, bestNodeNavigator, scalarCoordToBlockId, bestNodeBlocksToSeeds, bestNodeBlockSizes, bestNodeBlockRanges, bestNodeGapMap, bestNodeInverseBlockIds, dfsIndexes);

    if (refNode.empty()) {
      std::cout << "best nods: " << bestNodeId << std::endl;
      std::cout << "best node score: " << placementScores[0].second << std::endl;
    } else {
      std::cout << "specified reference node: " << refNode << std::endl;
    }
    // Here bestNodeOnSeedsHash contains the best node's seeds




    std::string bestMatchSequence = "";
    std::string gappedSeq = T->getStringFromReference(bestNode->identifier, true);
    std::vector<int32_t> degap;
    for (int32_t i = 0; i < gappedSeq.size(); i ++) {
        char &c = gappedSeq[i];
        degap.push_back(bestMatchSequence.size());
        if (c != '-') {
            bestMatchSequence += c;
        }
    }


    /*
    std::unordered_map<size_t, std::vector<int32_t>> seedToRefPositions;
    for(int i = 0; i < bestNodeOnSeedsHash.size(); i++){
      if(bestNodeOnSeedsHash[i].has_value()){
        size_t seed = bestNodeOnSeedsHash[i].value().hash;

        if (seedToRefPositions.find(seed) == seedToRefPositions.end()) {
            seedToRefPositions[seed] = {};
        }
        seedToRefPositions[seed].push_back(degap[i]);
      }
    }
    */

    std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>>  seedToRefPositions;
    for(int i = 0; i < bestNodeOnSeedsHash.size(); i++){
      if(bestNodeOnSeedsHash[i].has_value()){
        size_t seed = bestNodeOnSeedsHash[i].value().hash;
        bool reversed = bestNodeOnSeedsHash[i].value().isReverse;
        int pos = degap[i];

        if (seedToRefPositions.find(seed) == seedToRefPositions.end()) {
            std::vector<uint32_t> a;
            std::vector<uint32_t> b;
            seedToRefPositions[seed] = std::make_pair(a,b);
        }

        if(reversed){
          seedToRefPositions[seed].second.push_back(pos);
        }else{
          seedToRefPositions[seed].first.push_back(pos);
        }
        
      }
    }

    //std::string refFileName = "REFERENCE";
    //Print out Reference
    if (refFileName.size() == 0) {
      refFileName = "panmap.reference.fa";
    }
    if(refFileName.size() > 0){
        std::ofstream outFile{refFileName};

        if (outFile.is_open()) {
            
            outFile << ">ref\n";
            outFile << bestMatchSequence << "\n";

            std::cerr << "Wrote reference fasta to " << refFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << refFileName << std::endl;
        }
    }

    if (aligner == "minimap2") {
      //Create SAM
      std::vector<char *> samAlignments;
      std::string samHeader;

      createSam(
          readSeeds,                 
          readSequences,             
          readQuals,                 
          readNames,                 
          bestMatchSequence,         
          seedToRefPositions,        
          samFileName,               
          k,                         
          pairedEndReads,            
          
          samAlignments,             
          samHeader                  
      );
      //Convert to BAM
      sam_hdr_t *header;
      bam1_t **bamRecords;

      createBam(
          samAlignments,
          samHeader,
          bamFileName,

          header,
          bamRecords
      );

      createMplpBcf(
        prefix,
        refFileName,
        bestMatchSequence,
        bamFileName,
        mpileupFileName
      );

      createVcfWithMutationMatrices(
        prefix,
        mpileupFileName,
        mutMat,
        vcfFileName,
        0.0011
      );

      // //std::string mpileupFileName = "MPILEUP";
      // //Convert to Mplp
      // char *mplpString;

      // createMplp(
      //     bestMatchSequence,
      //     header,
      //     bamRecords,
      //     samAlignments.size(),
      //     mpileupFileName,

      //     mplpString
      // );

      // //std::string vcfFileName = "VCF";
      // //Convert to VCF
      // createVcf(
      //     mplpString,
      //     mutMat,
      //     vcfFileName,
      //     false
      // );
    } else {
      // align with bwa aln
      std::cout << "Preparing arguments for bwa..." << std::endl;

      std::vector<std::string> idx_args = {"bwa", "index", refFileName};
      std::vector<std::string> aln_args1;
      std::vector<std::string> aln_args2;
      std::vector<std::string> samaln_args;
      std::vector<std::pair<int, char*>> samAlignmentPairs;
      std::vector<std::string> samHeaders;
      if (reads2Path.empty()) {
        aln_args1 = {"bwa", "aln", "-l", "1024", "-n", "0.01", "-o", "2", "-f", reads1Path + ".tmp.sai", refFileName, reads1Path};
        samaln_args = {"bwa", "samse", refFileName, reads1Path + ".tmp.sai", reads1Path};
      } else {
        aln_args1 = {"bwa", "aln", "-l", "1024", "-n", "0.01", "-o", "2", "-f", reads1Path + ".tmp.sai", refFileName, reads1Path};
        aln_args2 = {"bwa", "aln", "-l", "1024", "-n", "0.01", "-o", "2", "-f", reads2Path + ".tmp.sai", refFileName, reads2Path};
        samaln_args = {"bwa", "sampe", refFileName, reads1Path + ".tmp.sai", reads2Path + ".tmp.sai", reads1Path, reads2Path};
      }

      prepareAndRunBwa(idx_args, aln_args1, aln_args2, samaln_args, reads1Path, reads2Path, samAlignmentPairs, samHeaders);


      std::sort(samAlignmentPairs.begin(), samAlignmentPairs.end(), [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
          return a.first < b.first;
      });

      std::vector<char*> samAlignments(samAlignmentPairs.size());
      for (size_t i = 0; i < samAlignmentPairs.size(); ++i) {
        samAlignments[i] = samAlignmentPairs[i].second;
      }


      std::ofstream samOut{samFileName};
      for (const auto& header : samHeaders) {
        samOut << header;
      }
      for (const auto& line : samAlignments) {
        samOut << line << "\n";
      }
      samOut.close();
      std::cout << "Wrote sam data to " << samFileName << std::endl;
    }
}



template <typename SeedMutationsType, typename GapMutationsType>
void place_per_read_DFS(
  mutableTreeData& data, std::map<uint32_t, seeding::onSeedsHash>& onSeedsHashMap, mgsr::refSeedmers& seedmersIndex,
  SeedMutationsType& perNodeSeedMutations_Index, GapMutationsType& perNodeGapMutations_Index,
  std::vector<mgsr::Read>& reads, const std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex,
  mgsr::readSeedmers& readSeedmersIndex,
  std::unordered_map<std::string, std::string>& identicalPairs,
  int seedK, int seedS, int seedT, int seedL, bool openSyncmers, Tree* T, Node* node, globalCoords_t& globalCoords, CoordNavigator& navigator,
  std::vector<int64_t>& scalarCoordToBlockId, std::vector<std::unordered_set<int>>& BlocksToSeeds, std::vector<int>& BlockSizes,
  std::vector<std::pair<int64_t, int64_t>>& blockRanges, int64_t& dfsIndex, std::map<int64_t, int64_t>& gapMap,
  std::unordered_set<int64_t>& inverseBlockIds, const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate,
  const int& redoReadThreshold, const bool& recalculateScore, const bool& rescueDuplicates, const double& rescueDuplicatesThreshold, const double& excludeDuplicatesThreshold, std::vector<bool>& excludeReads,
  std::unordered_map<std::string, double>& kminmer_binary_coverage, mgsr::ReadScores& readScores, std::vector<std::pair<std::string, int32_t>>& totalScores,
  double& binary_intersect_kminmer_count, double& ref_kminmer_count, std::ofstream& debugOut
) {
  size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
  std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>, std::optional<size_t>, std::optional<bool>, std::optional<bool>, std::optional<int64_t>, std::optional<int64_t>>> seedChanges;
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBlocksBacktracks;
  std::vector<std::pair<bool, int64_t>> inverseBlockIdsBacktrack;
  std::vector<std::pair<size_t, std::pair<int32_t, double>>> readScoresBacktrack;
  std::vector<std::pair<size_t, mgsr::SeedmerStatus>> seedmerStatusBacktrack;
  double ref_kminmer_count_backtrack = 0.0;
  double binary_intersect_kminmer_count_backtrack = 0.0;

  std::map<int64_t, int64_t> degapCoordIndex;
  std::map<int64_t, int64_t> regapCoordIndex;
  std::vector<tupleRange> recompRanges;
  blockExists_t oldBlockExists = data.blockExists;
  blockStrand_t oldBlockStrand = data.blockStrand;
  Step method = Step::PLACE;

  applyMutations(data, blockMutationInfo, recompRanges,  mutationInfo, T, node, globalCoords, navigator, blockRanges, gapRunUpdates, gapRunBacktracks, oldBlockExists, oldBlockStrand, method == Step::PLACE, inverseBlockIds, inverseBlockIdsBacktrack);
  recompRanges.clear();
  
  processNodeMutations(perNodeGapMutations_Index, perNodeSeedMutations_Index, dfsIndex, gapMap, gapRunBacktracks, gapRunBlocksBacktracks, inverseBlockIds, blockRanges, data, onSeedsHashMap, seedChanges, T, seedK, globalCoords, navigator, num_cpus);

  makeCoordIndex(degapCoordIndex, regapCoordIndex, gapMap, blockRanges);

  for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  gapRunBlocksBacktracks.clear();
  oldBlockExists.clear();
  oldBlockStrand.clear();
  gapRunUpdates.clear();
  std::sort(seedChanges.begin(), seedChanges.end(), [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });
  updateSeedsMapAndBlocks(seedChanges, onSeedsHashMap, scalarCoordToBlockId, BlocksToSeeds);

  //                     beg      end      fhash   rhash   rev  
  std::vector<std::tuple<int32_t, int32_t, size_t, size_t, bool>> backTrackPositionMapChAdd;
  std::vector<int32_t> backTrackPositionMapErase;

  std::unordered_set<size_t> affectedSeedmers;
  auto& positionMap = seedmersIndex.positionMap;
  auto& hashToPositionsMap = seedmersIndex.hashToPositionsMap;
  updateSeedmersIndex(seedChanges, onSeedsHashMap, seedmersIndex, affectedSeedmers, seedK, seedL, backTrackPositionMapChAdd, backTrackPositionMapErase);

  auto& seedmerToReads = readSeedmersIndex.seedmerToReads;
  
  auto& oldSeedmerStatus = seedmersIndex.seedmerStatus;
  for (const size_t& hash : affectedSeedmers) {
    mgsr::SeedmerStatus refSeemderOldStatus;
    if (oldSeedmerStatus.find(hash) == oldSeedmerStatus.end()) {
      refSeemderOldStatus = mgsr::SeedmerStatus::NOT_EXIST;
    } else {
      refSeemderOldStatus = oldSeedmerStatus.at(hash);
    }
    const mgsr::SeedmerStatus& refSeemderNewStatus = seedmersIndex.getCurrentSeedmerStatus(hash);

    seedmerStatusBacktrack.emplace_back(std::make_pair(hash, refSeemderOldStatus));

    if (refSeemderOldStatus == refSeemderNewStatus) continue;

    if (refSeemderOldStatus == mgsr::SeedmerStatus::NOT_EXIST) {
      ref_kminmer_count += 1.0;
      ref_kminmer_count_backtrack -= 1.0;
      if (seedmerToReads.find(hash) != seedmerToReads.end()) {
        binary_intersect_kminmer_count += 1.0;
        binary_intersect_kminmer_count_backtrack -= 1.0;
      }
      oldSeedmerStatus[hash] = refSeemderNewStatus;
    } else if (refSeemderOldStatus == mgsr::SeedmerStatus::EXIST_UNIQUE) {
      if (refSeemderNewStatus == mgsr::SeedmerStatus::NOT_EXIST) {
        ref_kminmer_count -= 1.0;
        ref_kminmer_count_backtrack += 1.0;
        if (seedmerToReads.find(hash) != seedmerToReads.end()) {
          binary_intersect_kminmer_count -= 1.0;
          binary_intersect_kminmer_count_backtrack += 1.0;
        }
        oldSeedmerStatus.erase(hash);
      }
    } else if (refSeemderOldStatus == mgsr::SeedmerStatus::EXIST_DUPLICATE) {
      if (refSeemderNewStatus == mgsr::SeedmerStatus::NOT_EXIST) {
        ref_kminmer_count -= 1.0;
        ref_kminmer_count_backtrack += 1.0;
        if (seedmerToReads.find(hash) != seedmerToReads.end()) {
          binary_intersect_kminmer_count -= 1.0;
          binary_intersect_kminmer_count_backtrack += 1.0;
        }
        oldSeedmerStatus.erase(hash);
      }
    } else {
      std::cout << "Error: invalid seedmer status" << std::endl;
      exit(1);
    }
  }

  kminmer_binary_coverage[node->identifier] = binary_intersect_kminmer_count / std::min(ref_kminmer_count, static_cast<double>(seedmerToReads.size()));

  if (debug) {
    // print out seeds at node
    std::cout << std::endl;
    if (method == Step::PLACE) {
      std::cout << node->identifier << " place seedmers: ";
      for (const auto& seedmer : positionMap) {
        const auto& beg = seedmer.first;
        const auto& [end, fhash, rhash, rev] = seedmer.second;
        if (fhash != rhash) {
          std::cout << mgsr::degapGlobal(beg, degapCoordIndex) << "-" << mgsr::degapGlobal(end, degapCoordIndex) << ":" << std::min(fhash, rhash) << "|" << rev << " ";
          // std::cout << beg << "|" << mgsr::degapGlobal(beg, coordIndex) << "-" << end << "|" << mgsr::degapGlobal(end, coordIndex) << ":" << std::min(fhash, rhash) << "|" << rev << " ";
        }
      }
    }

    for (const auto& [hash, positions] : hashToPositionsMap) {
      if (positions.size() == 0) {
        std::cout << "Error: hashToPositionsMap contains empty positions" << std::endl;
        exit(1);
      }
      for (const auto& position : positions) {
        if (position->second.fhash != position->second.rhash) {
          if (std::min(position->second.fhash, position->second.rhash) != hash) {
            std::cout << "Error: min(fhash, rhash) != hash" << std::endl;
            exit(1);
          }
        } else {
          std::cout << "Error: fhash == rhash" << std::endl;
          exit(1);
        }
      }
    }
    std::cout << std::endl;


    std::cout << node->identifier << " true seedmers: ";
    auto seq = seed_annotated_tree::getStringAtNode(node, T, false);
    auto seedmers = extractSeedmers(seq, seedK, seedS, seedT, seedL, openSyncmers);
    for (const auto &[seedmer, hash, isReverse, startPos, endPos] : seedmers) {
      std::cout << startPos << "-" << endPos << ":" << hash << "|" << isReverse << " ";
    }
    std::cout << std::endl;
  }


  tbb::concurrent_vector<std::pair<size_t, boost::icl::split_interval_map<int32_t, int>>> readBackTrack;
  tbb::concurrent_vector<std::pair<size_t, std::unordered_set<int32_t>>> readDuplicateSetsBackTrack;
  tbb::concurrent_vector<std::pair<size_t, std::vector<std::pair<int32_t, bool>>>> readDuplicatesBackTrack;
  // std::cout << node->identifier << " SCORE: ";
  if (node->identifier == T->root->identifier) {
    readScores.assignDfsIndex(node->identifier, dfsIndex);
    readScoresBacktrack.reserve(reads.size());
    for (size_t i = 0; i < reads.size(); ++i) {
      mgsr::Read& curRead = reads[i];
      curRead.matches.clear();
      curRead.duplicates.clear();
      initializeMatches(curRead, positionMap, hashToPositionsMap);
      int64_t pseudoScore = getPseudoScore(curRead, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap, minimumCount, minimumScore, rescueDuplicates, rescueDuplicatesThreshold, dfsIndex);
      double  pseudoProb  = pow(errorRate, curRead.seedmersList.size() - pseudoScore) * pow(1-errorRate, pseudoScore);

      readScoresBacktrack.emplace_back(std::make_pair(i, readScores.getScoreAtCurrentNode(i)));
      readScores.setScore(i, pseudoScore, pseudoProb, readSeedmersDuplicatesIndex[i].size());
      readScores.addScoreMutation(dfsIndex, i, pseudoScore, pseudoProb);

      if (curRead.duplicates.size() > excludeDuplicatesThreshold * curRead.seedmersList.size()) excludeReads[i] = true;
    }
    totalScores.emplace_back(std::make_pair(node->identifier, readScores.totalScore));
  } else {
    readScores.assignDfsIndex(node->identifier, dfsIndex);
    if (affectedSeedmers.empty()) {
      identicalPairs[node->identifier] = node->parent->identifier;
    } else if (positionMap.empty()) {
      for (size_t i = 0; i < reads.size(); ++i) {
        readBackTrack.emplace_back(std::make_pair(i, reads[i].matches));
        readDuplicateSetsBackTrack.emplace_back(std::make_pair(i, reads[i].duplicates));
        reads[i].matches.clear();
        reads[i].duplicates.clear();

        if (readScores.scores[i].first != 0) {
          readScoresBacktrack.emplace_back(std::make_pair(i, readScores.scores[i]));
          readScores.setScore(i, 0, 0, readSeedmersDuplicatesIndex[i].size());
          readScores.addScoreMutation(dfsIndex, i, 0, 0);
        }
      }
    } else {
      std::unordered_map<uint32_t, std::vector<int32_t>> readToAffectedSeedmerIndex;
      std::vector<std::pair<uint32_t, std::vector<int32_t>>> readToAffectedSeedmerIndexVec;

      if (redoReadThreshold > 0) {
        for (const size_t& affectedSeedmer : affectedSeedmers) {
          const auto& affectedSeedmerToReads = seedmerToReads.find(affectedSeedmer);
          if (affectedSeedmerToReads == seedmerToReads.end()) continue;
          for (const auto& [readIndex, affectedSeedmerIndices] : affectedSeedmerToReads->second) {
            for (const auto& affectedSeedmerIndex : affectedSeedmerIndices) {
              readToAffectedSeedmerIndex[readIndex].push_back(affectedSeedmerIndex);
            }
          }
        }
      } else {
        for (const size_t& affectedSeedmer : affectedSeedmers) {
          const auto& affectedSeedmerToReads = seedmerToReads.find(affectedSeedmer);
          if (affectedSeedmerToReads == seedmerToReads.end()) continue;
          for (const auto& [readIndex, affectedSeedmerIndices] : affectedSeedmerToReads->second) {
            const auto& readtoAffectedSeedmerIndexIt = readToAffectedSeedmerIndex.find(readIndex);
            if (readtoAffectedSeedmerIndexIt == readToAffectedSeedmerIndex.end()) {
              readToAffectedSeedmerIndex[readIndex].push_back(0);
            }
          }
        }

        // for (const size_t& affectedSeedmer : affectedSeedmers) {
        //   const auto& affectedSeedmerToReads = seedmerToReads.find(affectedSeedmer);
        //   if (affectedSeedmerToReads == seedmerToReads.end()) continue;
        //   for (const auto& [readIndex, affectedSeedmerIndices] : affectedSeedmerToReads->second) {
        //     for (const auto& affectedSeedmerIndex : affectedSeedmerIndices) {
        //       readToAffectedSeedmerIndex[readIndex].push_back(affectedSeedmerIndex);
        //     }
        //   }
        // }
      }

      readToAffectedSeedmerIndexVec.reserve(readToAffectedSeedmerIndex.size());
      for (const auto& [readIndex, affectedSeedmerIndices] : readToAffectedSeedmerIndex) {
          readToAffectedSeedmerIndexVec.emplace_back(readIndex, std::move(affectedSeedmerIndices));
      }
      readToAffectedSeedmerIndex.clear();

      for (size_t i = 0; i < readToAffectedSeedmerIndexVec.size(); ++i) {
        const auto& [readIndex, affectedSeedmerIndices] = readToAffectedSeedmerIndexVec[i];
        mgsr::Read& curRead = reads[readIndex];

        if (affectedSeedmerIndices.empty()) {
          std::cout << "Error: affectedSeedmerIndices is empty" << std::endl;
          exit(1);
        }

        readBackTrack.emplace_back(std::make_pair(readIndex, curRead.matches));
        if (affectedSeedmerIndices.size() > redoReadThreshold) {
          readDuplicateSetsBackTrack.emplace_back(std::make_pair(readIndex, curRead.duplicates));
          curRead.matches.clear();
          curRead.duplicates.clear();
          initializeMatches(curRead, positionMap, hashToPositionsMap);
        } else {
          std::vector<std::pair<int32_t, bool>> curReadDuplicatesBackTrack;

          for (const auto& index : affectedSeedmerIndices) {
            const size_t& affectedSeedmer = curRead.seedmersList[index].hash;

            const auto& affectedSeedmerHashToPositionIt = hashToPositionsMap.find(affectedSeedmer);
            
            // if affected seedmer exists and unique in ref
            if (affectedSeedmerHashToPositionIt != hashToPositionsMap.end() && affectedSeedmerHashToPositionIt->second.size() == 1) {
              if (curRead.duplicates.find(index) != curRead.duplicates.end()) {
                curRead.duplicates.erase(index);
                curReadDuplicatesBackTrack.push_back({index, false});
              }

              // if not inside a range in read matches
              if (isContained(curRead.matches, index)) {
                auto curIntervalIt = curRead.matches.find(index);
                curRead.matches.subtract({discrete_interval<int32_t>::closed(index, index), curIntervalIt->second});
              }

              const auto& affectedPositionRefIt = *(affectedSeedmerHashToPositionIt->second.begin());
              int currev = curRead.seedmersList[index].rev == affectedPositionRefIt->second.rev ? 1 : 2;
              // check left and right flanking bases and merge if possible
              auto leftFlank = curRead.matches.find(index - 1);
              auto rightFlank = curRead.matches.find(index + 1);
              bool mergedWithLeft = false;
              bool mergedWithRight = false;

              if (leftFlank != curRead.matches.end() && leftFlank->second == currev) {
                // left flank exists check if mergeable -> left adjacent kminmers matching and in the same direction
                const size_t& leftFlankSeedmer = curRead.seedmersList[index-1].hash;
                const auto& leftFlankSeedmerRefIt = hashToPositionsMap.find(leftFlankSeedmer);
                if (leftFlankSeedmerRefIt != hashToPositionsMap.end() && leftFlankSeedmerRefIt->second.size() == 1) {
                  // leftFlankSeedmer exists and is unique in ref
                  size_t prevSeedmerRef = 0;
                  if (leftFlank->second == 1) {
                    auto prevIt = std::prev(affectedPositionRefIt);
                    while (prevIt->second.fhash == prevIt->second.rhash) {
                      if (prevIt == positionMap.begin()) {
                        prevSeedmerRef = std::numeric_limits<size_t>::max();
                        break;
                      }
                      --prevIt;
                    }
                    if (prevSeedmerRef == 0) prevSeedmerRef = std::min(prevIt->second.fhash, prevIt->second.rhash);
                  } else {
                    auto nextIt = std::next(affectedPositionRefIt);
                    while (nextIt->second.fhash == nextIt->second.rhash) {
                      ++nextIt;
                      if (nextIt == positionMap.end()) {
                        prevSeedmerRef = std::numeric_limits<size_t>::max();
                        break;
                      }
                    }
                    if (prevSeedmerRef == 0) prevSeedmerRef = std::min(nextIt->second.fhash, nextIt->second.rhash);
                  }
                  if (prevSeedmerRef == leftFlankSeedmer) {
                    // merge
                    auto newBeg = first(leftFlank->first);
                    auto newEnd = last(leftFlank->first) + 1;
                    auto newRev = leftFlank->second;

                    curRead.matches.subtract({leftFlank->first, leftFlank->second});
                    curRead.matches.add({discrete_interval<int32_t>::closed(newBeg, newEnd), newRev});
                    mergedWithLeft = true;
                  }
                }
              }

              if (rightFlank != curRead.matches.end() && rightFlank->second == currev) {
                // right flank exists check if mergeable -> right adjacent kminmers matching and in the same direction
                const size_t& rightFlankSeedmer = curRead.seedmersList[index+1].hash;
                const auto& rightFlankSeedmerRefIt = hashToPositionsMap.find(rightFlankSeedmer);
                if (rightFlankSeedmerRefIt != hashToPositionsMap.end() && rightFlankSeedmerRefIt->second.size() == 1) {
                  // rightFlankSeedmer exists and is unique in ref
                  size_t nextSeedmerRef = 0;
                  if (rightFlank->second == 1) {
                    auto nextIt = std::next(affectedPositionRefIt);
                    while (nextIt->second.fhash == nextIt->second.rhash) {
                      ++nextIt;
                      if (nextIt == positionMap.end()) {
                        nextSeedmerRef = std::numeric_limits<size_t>::max();
                        break;
                      }
                    }
                    if (nextSeedmerRef == 0) nextSeedmerRef = std::min(nextIt->second.fhash, nextIt->second.rhash);
                  } else {
                    auto prevIt = std::prev(affectedPositionRefIt);
                    while (prevIt->second.fhash == prevIt->second.rhash) {
                      if (prevIt == positionMap.begin()) {
                        nextSeedmerRef = std::numeric_limits<size_t>::max();
                        break;
                      }
                      --prevIt;
                    }
                    if (nextSeedmerRef == 0) nextSeedmerRef = std::min(prevIt->second.fhash, prevIt->second.rhash);
                  }
                  if (nextSeedmerRef == rightFlankSeedmer) {
                    // merge
                    auto newBeg = mergedWithLeft ? first(curRead.matches.find(index)->first) : first(rightFlank->first) - 1;
                    auto newEnd = last(rightFlank->first);
                    auto newRev = rightFlank->second;

                    curRead.matches.subtract({rightFlank->first, rightFlank->second});
                    if (mergedWithLeft) {
                      curRead.matches.subtract({discrete_interval<int32_t>::closed(newBeg, index), newRev});
                    }
                    curRead.matches.add({discrete_interval<int32_t>::closed(newBeg, newEnd), newRev});
                    mergedWithRight = true;
                  }
                }
              }

              if (!mergedWithLeft && !mergedWithRight) {
                // start new range
                curRead.matches.add({discrete_interval<int32_t>::closed(index, index), currev});
              }
            
            } else if (affectedSeedmerHashToPositionIt == hashToPositionsMap.end() || affectedSeedmerHashToPositionIt->second.size() > 1) {
              if (affectedSeedmerHashToPositionIt == hashToPositionsMap.end()) {
                if (curRead.duplicates.find(index) != curRead.duplicates.end()) {
                  curRead.duplicates.erase(index);
                  curReadDuplicatesBackTrack.push_back({index, false});
                }
              } else {
                if (curRead.duplicates.find(index) == curRead.duplicates.end()) {
                  curRead.duplicates.insert(index);
                  curReadDuplicatesBackTrack.push_back({index, true});
                }
              }

              // if inside a range in read matches
              if (isContained(curRead.matches, index)) {
                // split range
                auto curIntervalIt = curRead.matches.find(index);
                curRead.matches.subtract({discrete_interval<int32_t>::closed(index, index), curIntervalIt->second});
              }
            }
          }
          if (!curReadDuplicatesBackTrack.empty()) {
            readDuplicatesBackTrack.emplace_back(std::make_pair(readIndex, std::move(curReadDuplicatesBackTrack)));
          }
        }
        int64_t pseudoScore = getPseudoScore(curRead, seedmersIndex, degapCoordIndex, regapCoordIndex, maximumGap, minimumCount, minimumScore, rescueDuplicates, rescueDuplicatesThreshold, dfsIndex);
        double  pseudoProb  = pow(errorRate, curRead.seedmersList.size() - pseudoScore) * pow(1 - errorRate, pseudoScore);
        if (pseudoScore != readScores.scores[readIndex].first) {
          readScoresBacktrack.emplace_back(std::make_pair(readIndex, readScores.scores[readIndex]));
          readScores.setScore(readIndex, pseudoScore, pseudoProb, readSeedmersDuplicatesIndex[readIndex].size());
          readScores.addScoreMutation(dfsIndex, readIndex, pseudoScore, pseudoProb);
        }
        if (curRead.duplicates.size() > excludeDuplicatesThreshold * curRead.seedmersList.size()) excludeReads[readIndex] = true;
      }
      totalScores.emplace_back(std::make_pair(node->identifier, readScores.totalScore));
    }
  }

  // std::cout << std::endl;

  /* Recursive step */
  dfsIndex++;
  std::cout << "\rprocessed " << dfsIndex << " / " <<  T->allNodes.size() << " haplotypes" << std::flush;
  for (Node *child : node->children) {
    place_per_read_DFS(
      data, onSeedsHashMap, seedmersIndex, perNodeSeedMutations_Index, perNodeGapMutations_Index, reads, readSeedmersDuplicatesIndex, readSeedmersIndex,
      identicalPairs, seedK, seedS, seedT, seedL, openSyncmers, T, child, globalCoords, navigator,
      scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, maximumGap,
      minimumCount, minimumScore, errorRate, redoReadThreshold, recalculateScore, rescueDuplicates, rescueDuplicatesThreshold, excludeDuplicatesThreshold, excludeReads,
      kminmer_binary_coverage, readScores, totalScores, binary_intersect_kminmer_count, ref_kminmer_count, debugOut
    );
  }

  for (const auto &p : seedChanges) {
    const auto& [pos, oldVal, newVal, oldSeed, newSeed, oldIsReverse, newIsReverse, oldEndPos, newEndPos] = p;
    auto seedIt = onSeedsHashMap.find(pos);
    if (oldVal && newVal) { // UNDO seed at same pos changed
      seedIt->second.hash      = oldSeed.value();
      seedIt->second.endPos    = oldEndPos.value();
      seedIt->second.isReverse = oldIsReverse.value();
    } else if (oldVal && !newVal) { // seed on to off
      onSeedsHashMap[pos] = {oldSeed.value(), oldEndPos.value(), oldIsReverse.value()};
      int blockId = scalarCoordToBlockId[pos];
      BlocksToSeeds[blockId].insert(pos);
    } else if (!oldVal && newVal) { // UNDO seed off to on
      onSeedsHashMap.erase(seedIt);
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

  // undo positionMap adds/changes
  for (const auto& [pos, end, fhash, rhash, rev] : backTrackPositionMapChAdd) {
    auto oldPositionMapIt = positionMap.find(pos);
    if (oldPositionMapIt != positionMap.end()) {
      const auto& [oldEnd, oldFHash, oldRHash, oldRev] = oldPositionMapIt->second;
      if (oldFHash != oldRHash) {
        size_t minHash = std::min(oldFHash, oldRHash);
        hashToPositionsMap[minHash].erase(oldPositionMapIt);
        if (hashToPositionsMap[minHash].empty()) {
          hashToPositionsMap.erase(minHash);
        }
      }

      oldPositionMapIt->second = {end, fhash, rhash, rev};
      if (fhash != rhash) {
        hashToPositionsMap[std::min(fhash, rhash)].insert(oldPositionMapIt);
      }
    } else {
      auto newPositionMapIt = positionMap.emplace(pos, positionInfo(end, fhash, rhash, rev)).first;
      if (fhash != rhash) {
        hashToPositionsMap[std::min(fhash, rhash)].insert(newPositionMapIt);
      }
    }
  }

  // undo positionMap erases
  for (const auto& pos : backTrackPositionMapErase) {
    auto toEraseIt = positionMap.find(pos);
    const auto& [toEraseEnd, toEraseFHash, toEraseRHash, toEraseRev] = toEraseIt->second;
    if (toEraseFHash != toEraseRHash) {
      size_t minHash = std::min(toEraseFHash, toEraseRHash);
      hashToPositionsMap[minHash].erase(toEraseIt);
      if (hashToPositionsMap[minHash].empty()) {
        hashToPositionsMap.erase(minHash);
      }
    }
    positionMap.erase(toEraseIt);
  }

  // undo read matches changes
  for (size_t i = 0; i < readBackTrack.size(); ++i) {
    const auto& [readIdx, matches] = readBackTrack[i];
    reads[readIdx].matches = std::move(matches);
  }

  // undo read duplicate sets changes
  for (size_t i = 0; i < readDuplicateSetsBackTrack.size(); ++i) {
    const auto& [readIdx, duplicates] = readDuplicateSetsBackTrack[i];
    reads[readIdx].duplicates = std::move(duplicates);
  }

  // undo read duplicates changes
  for (size_t i = 0; i < readDuplicatesBackTrack.size(); ++i) {
    const auto& [readIdx, changes] = readDuplicatesBackTrack[i];
    for (const auto& [hash, del] : changes) {
      if (del) {
        reads[readIdx].duplicates.erase(hash);
      } else {
        reads[readIdx].duplicates.insert(hash);
      }
    }
  }

  // undo seedmer status changes
  for (const auto& [hash, status] : seedmerStatusBacktrack) {
    if (status == mgsr::SeedmerStatus::NOT_EXIST) {
      oldSeedmerStatus.erase(hash);
    } else {
      oldSeedmerStatus[hash] = status;
    }
  }

  // undo kminmer counts
  ref_kminmer_count += ref_kminmer_count_backtrack;
  binary_intersect_kminmer_count += binary_intersect_kminmer_count_backtrack;


  // undo read scores changes
  for (size_t i = 0; i < readScoresBacktrack.size(); ++i) {
    const auto& [readIdx, score] = readScoresBacktrack[i];
    readScores.setScore(readIdx, score.first, score.second, readSeedmersDuplicatesIndex[readIdx].size());
  }

  /* Undo sequence mutations when backtracking */
  undoMutations(data, T, node, blockMutationInfo, mutationInfo, globalCoords);
}


void seedmersFromFastq(
  const std::string& fastqPath1, const std::string& fastqPath2, std::vector<mgsr::Read>& reads,
  mgsr::readSeedmers& readSeedmersIndex,
  std::vector<std::vector<size_t>>& readSeedmersDuplicatesIndex, std::vector<std::string>& readSequences,
  std::vector<std::string>& readQuals, std::vector<std::string>& readNames, std::vector<std::vector<seeding::seed>>& readSeeds,
  const int32_t& k, const int32_t& s, const int32_t& t, const int32_t& l, const bool& openSyncmers
) {
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
      readSequences.push_back(seq->seq.s);
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

  size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

  // index duplicate reads
  std::vector<size_t> sortedReadSequencesIndices(readSequences.size());
  for (size_t i = 0; i < readSequences.size(); ++i) sortedReadSequencesIndices[i] = i;
  tbb::parallel_sort(sortedReadSequencesIndices.begin(), sortedReadSequencesIndices.end(), [&readSequences](size_t i1, size_t i2) {
    return readSequences[i1] < readSequences[i2];
  });

  // index duplicate reads
  std::vector<std::pair<std::string, std::vector<size_t>>> dupReadsIndex;
  std::string prevSeq = readSequences[sortedReadSequencesIndices[0]];
  dupReadsIndex.emplace_back(std::make_pair(prevSeq, std::vector<size_t>{0}));
  for (size_t i = 1; i < sortedReadSequencesIndices.size(); ++i) {
    std::string currSeq = readSequences[sortedReadSequencesIndices[i]];
    if (currSeq == prevSeq) {
      dupReadsIndex.back().second.push_back(i);
    } else {
      dupReadsIndex.emplace_back(std::make_pair(currSeq, std::vector<size_t>{i}));
    }
    prevSeq = std::move(currSeq);
  }

  // seedmers for each unique read sequence
  std::vector<mgsr::Read> uniqueReadSeedmers(dupReadsIndex.size());
  tbb::parallel_for(tbb::blocked_range<size_t>(0, dupReadsIndex.size(), dupReadsIndex.size() / num_cpus),
    [&](const tbb::blocked_range<size_t>& range){
      for (size_t i = range.begin(); i < range.end(); ++i) {
        const auto& syncmers = seeding::rollingSyncmers(dupReadsIndex[i].first, k, s, openSyncmers, t, false);
        mgsr::Read& curRead = uniqueReadSeedmers[i];
        if (syncmers.size() < l) continue;

        size_t forwardRolledHash = 0;
        size_t reverseRolledHash = 0;
        // first kminmer
        for (size_t i = 0; i < l; ++i) {
          forwardRolledHash = rol(forwardRolledHash, k) ^ std::get<0>(syncmers[i]);
          reverseRolledHash = rol(reverseRolledHash, k) ^ std::get<0>(syncmers[l-i-1]);
        }

        int32_t iorder = 0;
        if (forwardRolledHash != reverseRolledHash) {
          size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
          curRead.uniqueSeedmers.emplace(minHash, std::vector<int32_t>{0});
          curRead.seedmersList.emplace_back(mgsr::readSeedmer{
            minHash, std::get<3>(syncmers[0]), std::get<3>(syncmers[l-1])+k-1, reverseRolledHash < forwardRolledHash, iorder});
          ++iorder;
        }

        // rest of kminmer
        for (size_t i = 1; i < syncmers.size()-l+1; ++i) {
          if (!std::get<2>(syncmers[i-1]) || !std::get<2>(syncmers[i+l-1])) {
            std::cout << "invalid syncmer" << std::endl;
            exit(0);
          }
          const size_t& prevSyncmerHash = std::get<0>(syncmers[i-1]);
          const size_t& nextSyncmerHash = std::get<0>(syncmers[i+l-1]);
          forwardRolledHash = rol(forwardRolledHash, k) ^ rol(prevSyncmerHash, k * l) ^ nextSyncmerHash;
          reverseRolledHash = ror(reverseRolledHash, k) ^ ror(prevSyncmerHash, k)     ^ rol(nextSyncmerHash, k * (l-1));

          if (forwardRolledHash != reverseRolledHash) {
            size_t minHash = std::min(forwardRolledHash, reverseRolledHash);
            auto uniqueSeedmersIt = curRead.uniqueSeedmers.find(minHash);
            if (uniqueSeedmersIt == curRead.uniqueSeedmers.end()) {
              curRead.uniqueSeedmers.emplace(minHash, std::vector<int32_t>{i});
              curRead.seedmersList.emplace_back(mgsr::readSeedmer{
                minHash, std::get<3>(syncmers[i]), std::get<3>(syncmers[i+l-1])+k-1, reverseRolledHash < forwardRolledHash, iorder});
              ++iorder;
            } else {
              uniqueSeedmersIt->second.push_back(i);
              curRead.seedmersList.emplace_back(mgsr::readSeedmer{
                uniqueSeedmersIt->first, std::get<3>(syncmers[i]), std::get<3>(syncmers[i+l-1])+k-1, reverseRolledHash < forwardRolledHash, iorder});
              ++iorder;
            }
          }
        }
      }
  });

  std::vector<size_t> sortedUniqueReadSeedmersIndices(uniqueReadSeedmers.size());
  for (size_t i = 0; i < uniqueReadSeedmers.size(); ++i) sortedUniqueReadSeedmersIndices[i] = i;
  tbb::parallel_sort(sortedUniqueReadSeedmersIndices.begin(), sortedUniqueReadSeedmersIndices.end(), [&uniqueReadSeedmers](size_t i1, size_t i2) {
    const auto& lhs = uniqueReadSeedmers[i1].seedmersList;
    const auto& rhs = uniqueReadSeedmers[i2].seedmersList;
    
    // First, compare the sizes of the seedmersList
    if (lhs.size() != rhs.size()) {
      return lhs.size() < rhs.size();
    }
    
    // If sizes are equal, compare hash values first
    for (size_t i = 0; i < lhs.size(); ++i) {
      if (lhs[i].hash != rhs[i].hash) {
        return lhs[i].hash < rhs[i].hash;
      }
    }
    
    // If all hash values are equal, compare other fields
    return std::lexicographical_compare(
      lhs.begin(), lhs.end(),
      rhs.begin(), rhs.end(),
      [](const mgsr::readSeedmer& a, const mgsr::readSeedmer& b) {
        if (a.begPos != b.begPos) return a.begPos < b.begPos;
        if (a.endPos != b.endPos) return a.endPos < b.endPos;
        if (a.rev != b.rev) return a.rev < b.rev;
        return a.iorder < b.iorder;
      }
    );
  });

  reads.emplace_back(std::move(uniqueReadSeedmers[sortedUniqueReadSeedmersIndices[0]]));
  readSeedmersDuplicatesIndex.emplace_back(std::vector<size_t>());
  for (const auto& seqSortedIndex : dupReadsIndex[sortedUniqueReadSeedmersIndices[0]].second) {
    readSeedmersDuplicatesIndex.back().push_back(sortedReadSequencesIndices[seqSortedIndex]);
  }

  for (size_t i = 1; i < sortedUniqueReadSeedmersIndices.size(); ++i) {
    const auto& currSeedmers = uniqueReadSeedmers[sortedUniqueReadSeedmersIndices[i]];

    if (!(currSeedmers.seedmersList.size() == reads.back().seedmersList.size() &&
          std::equal(currSeedmers.seedmersList.begin(), currSeedmers.seedmersList.end(), reads.back().seedmersList.begin(), reads.back().seedmersList.end(),
                     [](const mgsr::readSeedmer& a, const mgsr::readSeedmer& b) {
                         return a.hash == b.hash && a.begPos == b.begPos && a.endPos == b.endPos && a.rev == b.rev && a.iorder == b.iorder;
                     }))) {
      reads.emplace_back(std::move(uniqueReadSeedmers[sortedUniqueReadSeedmersIndices[i]]));
      readSeedmersDuplicatesIndex.emplace_back(std::vector<size_t>());
    }
    for (const auto& seqSortedIndex : dupReadsIndex[sortedUniqueReadSeedmersIndices[i]].second) {
      readSeedmersDuplicatesIndex.back().push_back(sortedReadSequencesIndices[seqSortedIndex]);
    }
  }

  auto& seedmerToReads = readSeedmersIndex.seedmerToReads;
  // auto& seedmerStatus = readSeedmersIndex.seedmerStatus;

  for (int32_t i = 0; i < reads.size(); ++i) {
    for (const auto& seedmer : reads[i].uniqueSeedmers) {
      seedmerToReads[seedmer.first].push_back(std::make_pair(i, std::move(seedmer.second)));
      // if (seedmerStatus.find(seedmer.first) == seedmerStatus.end()) {
      //   seedmerStatus[seedmer.first] = mgsr::SeedmerStatus::NOT_EXIST;
      // }
    }
    // reads[i].uniqueSeedmers.clear();
  }
}

bool identicalReadScores(const tbb::concurrent_vector<std::pair<int32_t, double>>& scores1, const tbb::concurrent_vector<std::pair<int32_t, double>>& scores2) {
    assert(scores1.size() == scores2.size());
    for (size_t i = 0; i < scores1.size(); ++i) {
        if (scores1[i].first != scores2[i].first) return false;
    }
    return true;
}

void updateIdenticalSeedmerSets(
  const std::unordered_set<std::string>& identicalGroup,
  mgsr::ReadScores& readScores,
  std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor,
  std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets
) {
  std::unordered_set<std::string> seenNodes;
  std::unordered_set<std::string> unseenNodes = identicalGroup;
  for (const std::string& currNode : identicalGroup) {
    if (seenNodes.find(currNode) != seenNodes.end()) continue;
    if (leastRecentIdenticalAncestor.find(currNode) != leastRecentIdenticalAncestor.end()) {
      std::cerr << "Error: Node " << currNode << " already has a least recent identical ancestor." << std::endl;
      exit(1);
    }
    seenNodes.insert(currNode);
    unseenNodes.erase(currNode);
    std::unordered_set<std::string> identicals;
    for (const std::string& idenNode : unseenNodes) {
      if (readScores.identicalReadScores(currNode, idenNode)) {
        identicals.insert(idenNode);
        identicalSets[currNode].insert(idenNode);
        leastRecentIdenticalAncestor[idenNode] = currNode;
        if (identicalSets.find(idenNode) != identicalSets.end()) {
          for (const auto& idenOffspring : identicalSets[idenNode]) {
            leastRecentIdenticalAncestor[idenOffspring] = currNode;
            identicalSets[currNode].insert(idenOffspring);
          }
          identicalSets.erase(idenNode);
        }
      }
    }
    for (const auto& identical : identicals) {
      seenNodes.insert(identical);
      unseenNodes.erase(identical);
    }
  }
}


void getConsensus(const char *lineStart, size_t lineLength) {
    char *line = new char[lineLength + 1];
    std::strncpy(line, lineStart, lineLength);
    line[lineLength] = '\0';

    std::vector<std::string> fields;
    stringSplit(line, '\t', fields);
    std::cout << fields[0] << "\t" << fields[1] << std::endl;

    delete[] line;
}

void getConsensusHelper(char *mplpString) {
  char *lineStart = mplpString;
  char *ptr = mplpString;

  while (*ptr != '\0') {
    if (*ptr == '\n') {
      getConsensus(lineStart, ptr - lineStart);
      lineStart = ptr + 1;
    }
    ++ptr;
  }

  if (lineStart != ptr) {
    getConsensus(lineStart, ptr - lineStart);
  }
}

void pmi::place_per_read(
  Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path,
  const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate,
  const int& redoReadThreshold, const bool& recalculateScore, const bool& rescueDuplicates,
  const double& rescueDuplicatesThreshold, const double& excludeDuplicatesThreshold,
  const std::string& preEMFilterMethod, const int& preEMFilterNOrder, const int& preEMFilterMBCNum, const int& emFilterRound,
  const int& checkFrequency, const int& removeIteration, const double& insigPropArg, const int& roundsRemove,
  const double& removeThreshold, const bool& leafNodesOnly, const bool& callSubconsensus, const std::string& prefix, const bool& save_kminmer_binary_coverage
)
{

  FILE* errorLog = freopen((prefix + ".error.log").c_str(), "w", stderr);
  if (!errorLog) {
      throw std::runtime_error("Failed to redirect stderr to error.log");
  }

  std::cout << "Wrote error log file: " << prefix + ".error.log" << std::endl;
  std::cerr << "Wrote error log file: " << prefix + ".error.log" << std::endl;

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
  int32_t t = index.getT();
  int32_t l = index.getL();
  bool openSyncmers = index.getOpen();

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
  
  std::map<uint32_t, seeding::onSeedsHash> onSeedsHashMap;    
  std::vector<std::string> readSequences;
  std::vector<std::string> readQuals;
  std::vector<std::string> readNames;
  std::vector<std::vector<seed>> readSeeds;
  std::vector<mgsr::Read> reads;
  std::vector<std::vector<size_t>> readSeedmersDuplicatesIndex;
  mgsr::readSeedmers readSeedmersIndex;

  auto read_processing_start = std::chrono::high_resolution_clock::now();
  seedmersFromFastq(reads1Path, reads2Path, reads, readSeedmersIndex, readSeedmersDuplicatesIndex, readSequences, readQuals, readNames, readSeeds, k, s, t, l, openSyncmers);   
  auto read_processing_end = std::chrono::high_resolution_clock::now();
  std::cerr << "Read processing time: " << std::chrono::duration_cast<std::chrono::milliseconds>(read_processing_end - read_processing_start).count() << " milliseconds" << std::endl;

  std::cerr << "Total reads: " << readNames.size() << std::endl;
  std::cout << "Total reads: " << readNames.size() << std::endl;
  std::cerr << "Total unique read kminmer sets: " << reads.size() << std::endl;
  std::cout << "Total unique read kminmer sets: " << reads.size() << std::endl;
  if (reads.size() != readSeedmersDuplicatesIndex.size()) {
    std::cerr << "Error: readSeedmersDuplicatesIndex size does not match reads size" << std::endl;
    exit(0);
  }

  mgsr::refSeedmers seedmersIndex;
  std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>> allScores;
  std::unordered_map<std::string, std::unordered_set<std::string>> identicalSets;
  std::unordered_map<std::string, std::string> leastRecentIdenticalAncestor;
  std::unordered_map<std::string, std::string> identicalPairs;
  mgsr::ReadScores readScores(T, reads.size(), T->allNodes.size());
  std::vector<std::pair<std::string, int32_t>> totalScores;
  totalScores.reserve(T->allNodes.size());

  // A = ref kminmers, B = read kminmers
  // kminmer_binary_coverage = A intersect B / A
  std::unordered_map<std::string, double> kminmer_binary_coverage;
  std::vector<bool> excludeReads(reads.size(), false);

  std::cerr << "start scoring DFS" << std::endl;
  std::cout << "start scoring DFS" << std::endl;
  
  auto start_time = std::chrono::high_resolution_clock::now();

  double binary_intersect_kminmer_count = 0;
  double ref_kminmer_count = 0;

  std::ofstream debugOut((prefix + ".debug.log").c_str());
  place_per_read_DFS<decltype(perNodeSeedMutations_Reader), decltype(perNodeGapMutations_Reader)>(
    data, onSeedsHashMap, seedmersIndex, perNodeSeedMutations_Reader, perNodeGapMutations_Reader, reads, readSeedmersDuplicatesIndex,
    readSeedmersIndex, identicalPairs, k, s, t, l, openSyncmers, T, T->root, globalCoords, navigator, scalarCoordToBlockId,
    BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap, inverseBlockIds, maximumGap, minimumCount, minimumScore,
    errorRate, redoReadThreshold, recalculateScore, rescueDuplicates, rescueDuplicatesThreshold, excludeDuplicatesThreshold,
    excludeReads, kminmer_binary_coverage, readScores, totalScores, binary_intersect_kminmer_count, ref_kminmer_count, debugOut
  );

  auto end_time = std::chrono::high_resolution_clock::now();

  std::cerr << "\nPseudo-chaining score execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " milliseconds" << std::endl;
  std::cout << "\nPseudo-chaining score execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " milliseconds" << std::endl;

  std::cout << "finished scoring DFS" << std::endl;
  std::cerr << "finished scoring DFS" << std::endl;


  onSeedsHashMap.clear();
  seedmersIndex.positionMap.clear();
  seedmersIndex.hashToPositionsMap.clear();
  readSeedmersIndex.seedmerToReads.clear();
  readSeedmersIndex.seedmerStatus.clear();

  for (const auto& pair : identicalPairs) {
      std::unordered_set<std::string> curIdenticals;
      std::string curNode = pair.first;
      std::string curParent = pair.second;
      curIdenticals.insert(curNode);
      while (identicalPairs.find(curParent) != identicalPairs.end()) {
          curNode = curParent;   
          curParent = identicalPairs[curParent];
          curIdenticals.insert(curNode);      
      }
      for (const auto& node : curIdenticals) {
          identicalSets[curParent].insert(node);
      }
  }

  for (const auto& set : identicalSets) {
      for (const auto& offspring : set.second) {
          leastRecentIdenticalAncestor[offspring] = set.first;
      }
  }


  std::cout << "First round of duplication removal: " << leastRecentIdenticalAncestor.size() << std::endl;
  std::cerr << "First round of duplication removal: " << leastRecentIdenticalAncestor.size() << std::endl;



  std::vector<std::pair<std::string, int32_t>> totalScoresSorted;
  totalScoresSorted.reserve(totalScores.size() - leastRecentIdenticalAncestor.size());
  for (const auto& node : totalScores) {
      if (leastRecentIdenticalAncestor.find(node.first) != leastRecentIdenticalAncestor.end()) continue;
      totalScoresSorted.emplace_back(std::make_pair(node.first, node.second));
  }
  std::sort(totalScoresSorted.begin(), totalScoresSorted.end(), [](const auto &a, const auto &b) {
      return a.second > b.second;
  });

  std::unordered_set<std::string> identicalGroup{totalScoresSorted[0].first};
  int32_t currGroupScore = totalScoresSorted[0].second;
  for (size_t i = 1; i < totalScoresSorted.size(); ++i) {
    const auto& currScore = totalScoresSorted[i];
    if (currScore.second == currGroupScore) {
      identicalGroup.insert(currScore.first);
    } else {
      if (!identicalGroup.empty()) {
        updateIdenticalSeedmerSets(identicalGroup, readScores, leastRecentIdenticalAncestor, identicalSets);
      }
      std::unordered_set<std::string>().swap(identicalGroup);
      identicalGroup.insert(currScore.first);
      currGroupScore = currScore.second;
    }
  }
  if (!identicalGroup.empty()) {
    updateIdenticalSeedmerSets(identicalGroup, readScores, leastRecentIdenticalAncestor, identicalSets);
  }


  // Sanity check
  for (const auto& node : identicalSets) {
    if (leastRecentIdenticalAncestor.find(node.first) != leastRecentIdenticalAncestor.end()) {
      std::cerr << "Error: Node " << node.first << " is in identicalSets but has a least recent identical ancestor." << std::endl;
      exit(1);
    }
    for (const auto& identicalNode : node.second) {
      if (leastRecentIdenticalAncestor.find(identicalNode) != leastRecentIdenticalAncestor.end()) {
        if (leastRecentIdenticalAncestor.at(identicalNode) != node.first) {
          std::cerr << "Error: Node " << identicalNode << " has a least recent identical ancestor of " << leastRecentIdenticalAncestor.at(identicalNode) << " but is in the identical set of " << node.first << "." << std::endl;
          exit(1);
        }
      } else {
        std::cerr << "Error: Node " << identicalNode << " is in identicalSets but does not have a least recent identical ancestor." << std::endl;
        exit(1);
      }
    }
  }

  std::cout << "Second round of duplication removal: " << leastRecentIdenticalAncestor.size() << std::endl;
  std::cerr << "Second round of duplication removal: " << leastRecentIdenticalAncestor.size() << "\n" << std::endl;


  size_t numReads = readSequences.size();
  std::atomic<size_t> numLowScoreReads = 0;
  std::vector<bool> lowScoreReads(reads.size(), false);
  std::vector<std::string> nodes;
  Eigen::MatrixXd probs;
  Eigen::VectorXd props;
  double llh;
  double insigProp =  insigPropArg <= 0 ? (1.0 / static_cast<double>(T->allNodes.size())) / 10.0 : insigPropArg;
  
  auto start = std::chrono::high_resolution_clock::now();

  size_t numcpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

  Eigen::setNbThreads(numcpus);


  mgsr::squaremHelper_test_1(
    T, readScores, readSeedmersDuplicatesIndex, lowScoreReads, numReads, numLowScoreReads, excludeReads,
    leastRecentIdenticalAncestor, identicalSets, probs, nodes, props, llh, preEMFilterMethod, preEMFilterNOrder, preEMFilterMBCNum,
    emFilterRound, checkFrequency, removeIteration, insigProp, roundsRemove, removeThreshold, leafNodesOnly, kminmer_binary_coverage, "", save_kminmer_binary_coverage, prefix);
  
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "EM execution time: " << elapsed.count() << " seconds" << std::endl;
  std::cerr << "EM execution time: " << elapsed.count() << " seconds" << std::endl;
  
  std::vector<std::pair<std::string, double>> sortedOut(nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) {
      sortedOut.at(i) = {nodes[i], props(i)};
  }

  std::sort(sortedOut.begin(), sortedOut.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
      return a.second > b.second;
  });

  std::string abundanceOutFile = prefix + ".abundance";
  std::ofstream abundanceOut(abundanceOutFile);
  for (size_t i = 0; i < sortedOut.size(); ++i) {
    const auto& node = sortedOut[i];
    abundanceOut << node.first;
    if (identicalSets.find(node.first) != identicalSets.end()) {
      for (const auto& identicalNode : identicalSets.at(node.first)) {
        abundanceOut << "," << identicalNode;
      }
    }
    abundanceOut << "\t" << std::fixed << std::setprecision(20) << node.second << "\n";
  }
  abundanceOut.close();

  std::cout << "Wrote abundance file: " << abundanceOutFile << std::endl;
  std::cerr << "Wrote abundance file: " << abundanceOutFile << std::endl;

  exit(0);
  if (callSubconsensus) {
    // calling consensus
    std::cout << "Calling consensus" << std::endl;
    std::cerr << "Calling consensus" << std::endl;

    std::unordered_map<std::string, std::vector<std::pair<size_t, int32_t>>> assignedReads;
    mgsr::assignReadsToNodes(allScores, nodes, reads, readSeedmersDuplicatesIndex, assignedReads);

    boost::filesystem::path readsDir(prefix + "_reads");
    if (!boost::filesystem::exists(readsDir)) {
      if (!boost::filesystem::create_directory(readsDir)) {
        std::cerr << "Error: Unable to create directory " << readsDir.string() << std::endl;
        exit(1);
      }
    }
    
    std::vector<std::vector<std::pair<int8_t, int32_t>>> assignedReferences(readNames.size());
    for (size_t i = 0; i < sortedOut.size(); ++i) {
      const auto& node = sortedOut[i];
      // write reads to fastq
      std::cout << "Writing assigned reads assigned to " << node.first << std::endl;
      std::cerr << "Writing assigned reads assigned to " << node.first << std::endl;
      std::string fastqPath1 = "";
      std::string fastqPath2 = "";
      if (reads2Path.size() > 0) {
        std::string sanitizedNodeFirst = node.first;
        std::replace(sanitizedNodeFirst.begin(), sanitizedNodeFirst.end(), '/', '_');
        std::replace(sanitizedNodeFirst.begin(), sanitizedNodeFirst.end(), '|', '_');
        fastqPath1 = (readsDir / (sanitizedNodeFirst + "_R1.fastq")).string();
        fastqPath2 = (readsDir / (sanitizedNodeFirst + "_R2.fastq")).string();
        std::ofstream fastqOut1(fastqPath1);
        std::ofstream fastqOut2(fastqPath2);
        for (const auto& read : assignedReads[node.first]) {
          const auto& readIdx = read.first;
          const auto& readScore = read.second;
          if (readIdx % 2 != 0) continue;
          fastqOut1 << "@" << readNames[readIdx] << "\n" << readSequences[readIdx] << "\n+\n" << readQuals[readIdx] << "\n";
          fastqOut2 << "@" << readNames[readIdx + 1] << "\n" << readSequences[readIdx + 1] << "\n+\n" << readQuals[readIdx + 1] << "\n";
          assignedReferences[readIdx].push_back(std::make_pair(i, readScore));
          assignedReferences[readIdx + 1].push_back(std::make_pair(i, readScore));
        }
        fastqOut1.close();
        fastqOut2.close();
      } else {
        std::string sanitizedNodeFirst = node.first;
        std::replace(sanitizedNodeFirst.begin(), sanitizedNodeFirst.end(), '/', '_');
        std::replace(sanitizedNodeFirst.begin(), sanitizedNodeFirst.end(), '|', '_');
        fastqPath1 = (readsDir / (sanitizedNodeFirst + ".fastq")).string();
        std::ofstream fastqOut(fastqPath1);
        for (const auto& read : assignedReads[node.first]) {
          const auto& readIdx = read.first;
          const auto& readScore = read.second;
          fastqOut << "@" << readNames[readIdx] << "\n" << readSequences[readIdx] << "\n+\n" << readQuals[readIdx] << "\n";
          assignedReferences[readIdx].push_back(std::make_pair(i, readScore));
        }
        fastqOut.close();
      }
    }
    std::ofstream assignedReferencesOut(prefix + "_assigned_references.tsv");
    for (size_t i = 0; i < assignedReferences.size(); ++i) {
      assignedReferencesOut << readNames[i] << "\t";
      for (size_t j = 0; j < assignedReferences[i].size(); ++j) {
        assignedReferencesOut << sortedOut[assignedReferences[i][j].first].first << "," << assignedReferences[i][j].second << " ";
      }
      assignedReferencesOut << "\n";
    }
    assignedReferencesOut.close();


    exit(0);
    // for (const auto& node : sortedOut) {
    //   std::cout << "Calling consensus for " << node.first << std::endl;
    //   // write reference fastas
    //   std::cerr << "Writing reference fastas for " << node.first << std::endl;
    //   std::string refPath = prefix + "." + node.first + ".fasta";
    //   std::string refSeq = T->getStringFromReference(node.first, false);
    //   std::ofstream refOut(refPath);
    //   refOut << ">" << node.first << "\n" << refSeq << "\n";
    //   refOut.close();
    //   std::cerr << "Finished writing reference fastas for " << node.first << std::endl;
    //   std::cout << "Finished writing reference fastas for " << node.first << std::endl;


    //   // write reads to fastq
    //   std::cerr << "Writing assigned reads assigned to " << node.first << std::endl;
    //   std::string fastqPath1 = "";
    //   std::string fastqPath2 = "";
    //   if (reads2Path.size() > 0) {
    //     fastqPath1 = prefix + "." + node.first + "_R1.fastq";
    //     fastqPath2 = prefix + "." + node.first + "_R2.fastq";
    //     std::unordered_set<size_t> assigned;
    //     std::ofstream fastqOut1(fastqPath1);
    //     std::ofstream fastqOut2(fastqPath2);
    //     for (size_t readIdx : assignedReads[node.first]) {
    //       readIdx = readIdx % 2 == 0 ? readIdx : readIdx - 1;
    //       if (assigned.find(readIdx) != assigned.end()) continue;
    //       fastqOut1 << "@" << readNames[readIdx] << "\n" << readSequences[readIdx] << "\n+\n" << readQuals[readIdx] << "\n";
    //       fastqOut2 << "@" << readNames[readIdx + 1] << "\n" << readSequences[readIdx + 1] << "\n+\n" << readQuals[readIdx + 1] << "\n";
    //       assigned.insert(readIdx);
    //     }
    //     fastqOut1.close();
    //     fastqOut2.close();
    //   } else {
    //     fastqPath1 = prefix + "." + node.first + ".fastq";
    //     std::ofstream fastqOut(fastqPath1);
    //     for (const auto& readIndex : assignedReads[node.first]) {
    //       fastqOut << "@" << readNames[readIndex] << "\n" << readSequences[readIndex] << "\n+\n" << readQuals[readIndex] << "\n";
    //     }
    //     fastqOut.close();
    //   }
    //   std::cerr << "Finished writing assigned reads assigned to " << node.first << std::endl;
    //   std::cout << "Finished writing assigned reads assigned to " << node.first << std::endl;
    //   std::cerr << "Running bwa aln for " << node.first << std::endl;
    //   std::string samPath = prefix + "." + node.first + ".sam";
    //   std::vector<std::string> idx_args = {"bwa", "index", refPath};
    //   std::vector<std::string> aln_args1;
    //   std::vector<std::string> aln_args2;
    //   std::vector<std::string> samaln_args;
    //   std::vector<std::pair<int, char*>> samAlignmentPairs;
    //   std::vector<std::string> samHeaders;
    //   if (fastqPath2.size() > 0) { 
    //     aln_args1 = {"bwa", "aln", "-l", "1024", "-n", "0.01", "-o", "2", "-f", fastqPath1 + ".tmp.sai", refPath, fastqPath1};
    //     aln_args2 = {"bwa", "aln", "-l", "1024", "-n", "0.01", "-o", "2", "-f", fastqPath2 + ".tmp.sai", refPath, fastqPath2};
    //     samaln_args = {"bwa", "sampe", refPath, fastqPath1 + ".tmp.sai", fastqPath2 + ".tmp.sai", fastqPath1, fastqPath2};
    //   } else {
    //     aln_args1 = {"bwa", "aln", "-l", "1024", "-n", "0.01", "-o", "2", "-f", fastqPath1 + ".tmp.sai", refPath, fastqPath1};
    //     samaln_args = {"bwa", "samse", refPath, fastqPath1 + ".tmp.sai", fastqPath1};
    //   }
      
    //   prepareAndRunBwa(idx_args, aln_args1, aln_args2, samaln_args, fastqPath1, fastqPath2, samAlignmentPairs, samHeaders);
    //   std::cerr << "Finished running bwa aln for " << node.first << std::endl;
    //   std::cout << "Finished running bwa aln for " << node.first << std::endl;


    //   std::sort(samAlignmentPairs.begin(), samAlignmentPairs.end(), [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
    //       return a.first < b.first;
    //   });

    //   std::vector<char*> samAlignments(samAlignmentPairs.size());
    //   for (size_t i = 0; i < samAlignmentPairs.size(); ++i) {
    //     samAlignments[i] = samAlignmentPairs[i].second;
    //   }

    //   std::ofstream samOut{samPath};
    //   for (const auto& header : samHeaders) {
    //     samOut << header;
    //   }
    //   for (const auto& line : samAlignments) {
    //     samOut << line << "\n";
    //   }
    //   samOut.close();
    //   std::cout << "Wrote sam data to " << samPath << std::endl;

    //   sam_hdr_t *header;
    //   bam1_t **bamRecords;
    //   std::string bamPath = prefix + "." + node.first + ".bam";
    //   std::cout << "Creating bam file for " << node.first << std::endl;
    //   std::string samHeader = samHeaders[0].substr(0, samHeaders[0].size() - 1);
    //   createBam(
    //       samAlignments,
    //       samHeader,
    //       bamPath,
    //       header,
    //       bamRecords
    //   );
    //   std::cout << "Finished creating bam file for " << node.first << std::endl;

    //   int numAlignments = samAlignments.size();
    //   samAlignments.clear();
    //   samAlignmentPairs.clear();

    //   std::cout << "Creating mpileup file for " << node.first << std::endl;

    //   char *mplpString;
    //   std::string mpileupPath = prefix + "." + node.first + ".mpileup";
    //   createMplp(
    //       refSeq,
    //       header,
    //       bamRecords,
    //       numAlignments,
    //       mpileupPath,
    //       mplpString
    //   );
    //   std::cout << "Finished creating mpileup file for " << node.first << std::endl;

    //   // getConsensusHelper(mplpString);

    //   // // delete intermediate files
    //   // std::cerr << "Cleaning up intermediate files for " << node.first << std::endl;
    //   // fs::remove(refPath);
    //   // fs::remove(refPath + ".amb");
    //   // fs::remove(refPath + ".ann");
    //   // fs::remove(refPath + ".bwt");
    //   // fs::remove(refPath + ".pac");
    //   // fs::remove(refPath + ".sa");
      

    //   // get vcf
    //   // resolve genotype conflicts
    //   // write consensus fasta
    // }    
  }

  if (fclose(errorLog) != 0) {
      throw std::runtime_error("Failed to close error.log");
  }
  freopen("/dev/tty", "w", stderr); // For Linux/Mac
}
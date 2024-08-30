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
#include <fstream>
#include <nlohmann/json.hpp>
#include <variant>
#include <curl/curl.h>  // Include libcurl for HTTP requests
#include "seed_annotated_tree.hpp"  // Ensure this header is included for blockStrand_t
#include <capnp/serialize.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include "index.capnp.h"
  #include <tbb/parallel_for.h>
  #include <tbb/blocked_range.h>

const bool DEBUG = false;


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

void sendDataToServer(const nlohmann::json& jsonData) {
    CURL* curl;
    CURLcode res;

    curl_global_init(CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init();
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, "http://localhost:5000/update-data");
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, jsonData.dump().c_str());

        res = curl_easy_perform(curl);
        if(res != CURLE_OK) {
            fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
        }
        curl_easy_cleanup(curl);
    }
    curl_global_cleanup();
}

void exportDataToJson(Tree* T) {
    nlohmann::json jsonData;
    for (const auto& node : T->allNodes) {
        jsonData["nodes"].push_back({
            {"id", node->identifier},
            {"genome", "example_genome_data"},  // Replace with actual genome data
            {"seeds", "example_seeds_data"},    // Replace with actual seeds data
            {"mutations", "example_mutations_data"}  // Replace with actual mutations data
        });
        sendDataToServer(jsonData);  // Send data to server
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
                    blockExists_t& oldBlockExists, blockStrand_t& oldBlockStrand, const bool isPlacement)
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

    if (startingStrand != blockStrand[primaryBlockId].first) {
      flipCoords(primaryBlockId, globalCoords);
      }
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
    if (!oldBlockExists[blockId].first && blockExists[blockId].first) {
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
                                             blockStrand_t &blockStrand, tupleCoord_t stop_coord)
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


void updateGapMapStep(std::map<int64_t, int64_t>& gapMap, const std::pair<bool, std::pair<int64_t, int64_t>>& update, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack) {
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
      } else {
        // insert new range
        auto tmpIt = gapMap.emplace(start, end);
        curIt = tmpIt.first;
        backtrack.emplace_back(true, std::make_pair(curIt->first, curIt->second));
      }
    } else {
      curIt = leftIt;
      if (end <= curIt->second) {
        return;
      }
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

void updateGapMap(std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& updates, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack) {
  for (const auto& update : updates) {
    updateGapMapStep(gapMap, update, backtrack);
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

void invertGapMap(std::map<int64_t, int64_t>& gapMap, const std::pair<int64_t, int64_t>& invertRange, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack) {
  const auto& [start, end] = invertRange;

  auto rightIt = gapMap.upper_bound(start);
  auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

  bool rightItExists = rightIt != gapMap.end();
  bool leftItExists = leftIt != gapMap.end();

  // completely inside or outside a gap range -> do nothing
  if (
    gapMap.empty() || // empty gap map
    (!leftItExists && end < rightIt->first) || // completely left of first gap range
    (!rightItExists && start > leftIt->second) || // completely right of last gap range
    (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first) || // completely between two gap ranges
    (leftItExists && start >= leftIt->first && end <= leftIt->second) // completely inside a gap range
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
    updateGapMapStep(gapMap, {it->first, {curBeg, curEnd}}, backtrack);
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
  for (int64_t i = 0; i < ungapped.size() - k + 1; ++i) {
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
template <typename SeedMutationsType, typename GapMutationsType>
void buildOrPlace(Step method, mutableTreeData& data, std::vector<bool>& seedVec, std::vector<std::optional<std::string>>& onSeeds, SeedMutationsType& perNodeSeedMutations_Index, GapMutationsType& perNodeGapMutations_Index, int seedK, int seedS, Tree* T, Node* node, globalCoords_t& globalCoords, CoordNavigator& navigator, std::vector<int64_t>& scalarCoordToBlockId, std::vector<std::unordered_set<int>>& BlocksToSeeds, std::vector<int>& BlockSizes, std::vector<std::pair<int64_t, int64_t>>& blockRanges, int dfsIndex, std::map<int64_t, int64_t>& gapMap) {
  // Variables needed for both build and place
  std::vector<std::tuple<int64_t, bool, bool, std::optional<std::string>, std::optional<std::string>>> seedChanges;
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::map<int64_t, int64_t> coordIndex;
  std::vector<tupleRange> recompRanges;
  blockExists_t oldBlockExists = data.blockExists;
  blockStrand_t oldBlockStrand = data.blockStrand;

  applyMutations(data, blockMutationInfo, recompRanges,  mutationInfo, T, node, globalCoords, navigator, blockRanges, gapRunUpdates, gapRunBacktracks, oldBlockExists, oldBlockStrand, method == Step::PLACE);

  if (method == Step::BUILD) {
    
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBlocksBacktracks;

    if (gappity) {
      updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks);
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
          updateGapMapStep(gapMap, {true, {start, end}}, gapRunBacktracks);
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
              updateGapMapStep(gapMap, {false, range}, gapRunBacktracks);
            }
          } else {
            std::vector<std::pair<int64_t, int64_t>> invertedRanges = invertRanges(nucRanges, blockRanges[i]);
            for (const auto& range : invertedRanges) {
              updateGapMapStep(gapMap, {false, range}, gapRunBacktracks);
            }
          }
        }
      }

      for (const auto& blockId : invertedBlocks) {
        invertGapMap(gapMap, blockRanges[blockId], gapRunBlocksBacktracks);
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
          if (seedVec[pos] == true) {
            std::string oldSeed = pos < onSeeds.size() && onSeeds[pos].has_value() ? onSeeds[pos].value() : "";
            seedChanges.emplace_back(std::make_tuple(pos,
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt));
          }
        }
      }

      //Loop through seeds that now start as gaps and delete them
      for(int i = 0; i < gaps.size(); i++){
        if (seedVec[gaps[i]] == true){
          std::string oldSeed = gaps[i] < onSeeds.size() && onSeeds[gaps[i]].has_value() ? onSeeds[gaps[i]].value() : "";
          seedChanges.emplace_back(std::make_tuple(gaps[i], 
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt)); 
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
          bool inMap = (seedVec[coords[i - seedK]] == true);
          bool isSeed = (first_min_coord == i - seedK || last_min_coord == i - seedS) && Ns <= seedK/2;

          if(!inMap && isSeed){
            //Add seed
            std::string newSeed = seq.substr(i - seedK, seedK);
            seedChanges.emplace_back(std::make_tuple(coords[i - seedK],
                                                false, // old seed off
                                                true, // new seed on
                                                std::nullopt,
                                                newSeed));
            
          }else if(inMap && !isSeed){
            //Remove Seed
            std::string oldSeed = coords[i - seedK] < onSeeds.size() && onSeeds[coords[i - seedK]].has_value() ? onSeeds[coords[i - seedK]].value() : "";
            seedChanges.emplace_back(std::make_tuple(coords[i - seedK],
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt));

          }else if(inMap && isSeed){
            std::string oldSeed = coords[i - seedK] < onSeeds.size() && onSeeds[coords[i - seedK]].has_value() ? onSeeds[coords[i - seedK]].value() : "";
            std::string newSeed = seq.substr(i - seedK, seedK);
            if(newSeed != oldSeed){
              seedChanges.emplace_back(std::make_tuple(coords[i - seedK],
                                                true, // old seed on
                                                true, // new seed on
                                                oldSeed,
                                                newSeed));
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
          if (seedVec[coords[i]] == true){
            std::string oldSeed = coords[i] < onSeeds.size() && onSeeds[coords[i]].has_value() ? onSeeds[coords[i]].value() : "";
            seedChanges.emplace_back(std::make_tuple(coords[i],
                                                true, // old seed on
                                                false, // new seed off
                                                oldSeed,
                                                std::nullopt));
          }
        }
      }
    }
    

    int seedChangeIndex = seedChanges.size() - 1;
    int maskIndex = 0;

    std::vector<int64_t> basePositions;
    std::vector<std::vector<std::vector<int8_t>>> masks_all;

    while (seedChangeIndex >= 0) {
      auto [pos, oldVal, newVal, oldSeed, newSeed] = seedChanges[seedChangeIndex];
      basePositions.push_back(pos);
      std::vector<std::vector<int8_t>> masks;
      std::vector<int8_t> mask;
      int maskIndex = 0;
      // Generate ternary numbers for the current seed position
        while (seedChangeIndex >= 0) {
          auto [pos, oldVal, newVal, oldSeed, newSeed] = seedChanges[seedChangeIndex];
          int8_t ternaryNumber = 0;

          if (oldVal && newVal) {
              ternaryNumber = 2; // changed/inserted
          } else if (oldVal && !newVal) {
              ternaryNumber = 1; // deleted
          } else if (!oldVal && newVal) {
              ternaryNumber = 2; // changed/inserted
          } else {
              ternaryNumber = 0; // same
          }

          if (maskIndex < 20) {
            mask.push_back(ternaryNumber);
          } else {
            if (std::all_of(mask.begin(), mask.end(), [](int8_t n) { return n == 0; })) {
              seedChangeIndex--;
              break;
            }
            masks.push_back(mask);
            mask.clear();
          }
          maskIndex++;
          seedChangeIndex--;
        }
        masks_all.push_back(masks);
    }



    if constexpr (std::is_same_v<SeedMutationsType, ::capnp::List<SeedMutations>::Builder>) {
        auto basePositionsBuilder = perNodeSeedMutations_Index[dfsIndex].initBasePositions(basePositions.size());
        auto perPosMasksBuilder = perNodeSeedMutations_Index[dfsIndex].initPerPosMasks(masks_all.size());
        for (int i = 0; i < masks_all.size(); i++) {
          const auto &masks = masks_all[i];
          auto maskBuilder = perPosMasksBuilder[i].initMasks(masks.size());
          // Store the ternary numbers in the perMutMasks
          for (int j = 0; j < masks.size(); j++) {
            const auto& mask = masks[j];
            uint32_t tritMask = 0;
            for (int k = 0; k < mask.size(); k++) {
              tritMask |= (mask[k] & 0x3) << (k * 2);
            }
            maskBuilder.set(j, tritMask);
          }
        }
        for (int i = 0; i < basePositions.size(); i++) {
          basePositionsBuilder.set(i, basePositions[i]);
        }
    }
    auto nodeGapBuilder = perNodeGapMutations_Index[dfsIndex];
    if constexpr (std::is_same_v<GapMutationsType, ::capnp::List<GapMutations>::Builder>) {
      auto gapMutationsBuilder = nodeGapBuilder.initDeltas(gapRunBacktracks.size());
      for (auto it = gapRunBacktracks.rbegin(); it != gapRunBacktracks.rend(); ++it) {
        const auto& [del, range] = *it;
        if (del) {
          gapMutationsBuilder[gapMutationsBuilder.size() - 1].setPos(range.first);
          gapMutationsBuilder[gapMutationsBuilder.size() - 1].initMaybeValue().setValue(range.second);
          gapMap[range.first] = range.second;
        } else {
          gapMutationsBuilder[gapMutationsBuilder.size() - 1].setPos(range.first);
          gapMutationsBuilder[gapMutationsBuilder.size() - 1].initMaybeValue().setNone();
          gapMap.erase(range.first);
        }
      }
    }
    

    merged.clear();
    recompRanges.clear();
    gapRunUpdates.clear();
    
    
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
        gapMap[pos] = maybeValue.getValue();
        gapRunBacktracks.emplace_back(true, std::make_pair(pos, maybeValue.getValue()));
      } else {
        gapMap.erase(pos);
        gapRunBacktracks.emplace_back(false, std::make_pair(pos, 0));
      }
    }

    auto currBasePositions = perNodeSeedMutations_Index[dfsIndex].getBasePositions();
    auto currPerPosMasks = perNodeSeedMutations_Index[dfsIndex].getPerPosMasks();
    std::vector<std::vector<int8_t>> masks;
    for (int i = 0; i < currBasePositions.size(); ++i) {
        int64_t pos = currBasePositions[i];
        auto masks = currPerPosMasks[i].getMasks();

        for (int j = 0; j < masks.size(); ++j) {
            uint32_t tritMask = masks[j];
            for (int k = 0; k < 32; k += 2) {
                uint8_t ternaryNumber = (tritMask >> k) & 0x3;
                if (ternaryNumber == 1) { // on -> off
                    // Handle deletion
                    seedChanges.emplace_back(std::make_tuple(pos + k / 2, true, false, onSeeds[pos + k / 2], std::nullopt));
                } else if (ternaryNumber == 2) {
                    // Handle insertion/change
                    std::string newSeed = seed_annotated_tree::getSeedAt(pos + k / 2, T, seedK, data.scalarToTupleCoord, data.sequence, data.blockExists, data.blockStrand, globalCoords, navigator, gapMap);
                    if (seedVec[pos + k / 2]) { // on -> on
                        seedChanges.emplace_back(std::make_tuple(pos + k / 2, true, true, onSeeds[pos + k / 2], newSeed));
                    } else { // off -> on
                        seedChanges.emplace_back(std::make_tuple(pos + k / 2, false, true, std::nullopt, newSeed));
                    }
                }
            }
        }
      }
  }
  

    for (const auto &p : seedChanges)
    {
      bool oldVal = std::get<1>(p);
      bool newVal = std::get<2>(p);
      int64_t pos = std::get<0>(p);
      if (oldVal && newVal) { // seed at same pos changed
        seedVec[pos] = true;
        onSeeds[pos] = std::get<4>(p);
      } else if (oldVal && !newVal) { // seed on to off
        seedVec[std::get<0>(p)] = false;
        if (onSeeds[std::get<0>(p)].has_value() && std::get<0>(p) < onSeeds.size()) {
          onSeeds[std::get<0>(p)].reset();
        }
        int blockId = scalarCoordToBlockId[std::get<0>(p)];
        BlocksToSeeds[blockId].erase(std::get<0>(p));
      } else if (!oldVal && newVal) { // seed off to on
        seedVec[std::get<0>(p)] = true;
        onSeeds[std::get<0>(p)] = std::get<4>(p);
        int blockId = scalarCoordToBlockId[std::get<0>(p)];
        BlocksToSeeds[blockId].insert(std::get<0>(p));
      } 
  }
  
  if (debug) {
    // print out seeds at node
    std::cout << node->identifier << " build syncmers: ";
    // for (const auto& [pos, numGaps] : gapMap) {
    //   std::cout << pos << ":" << numGaps << " ";
    // }
    for (int i = 0; i < seedVec.size(); i++) {
      if (seedVec[i]) {
        std::cout << degapGlobal(i, coordIndex) << ":" << onSeeds[i].value() << " ";
      }
    }

    std::cout << std::endl;
    std::cout << node->identifier << " true syncmers: ";
    tupleCoord_t startCoord= {0,0,0};
    tupleCoord_t endCoord = {data.sequence.size() - 1, data.sequence.back().first.size() - 1, data.sequence.back().first.back().second.empty() ? -1 : 0};
    auto [seq, coords, gaps, deadBlocks] = seed_annotated_tree::getNucleotideSequenceFromBlockCoordinates(
      startCoord, endCoord, data.sequence, data.blockExists, data.blockStrand, T, node, globalCoords, navigator);
    
    auto syncmers = extractSeedmers(seq, seedK, seedS, /*open=*/false);
    for (const auto &[kmer, startPos, endPos] : syncmers) {
      std::cout << startPos << ":" << kmer << " ";
    }
    std::cout << std::endl;
  }

  /* Recursive step */
  dfsIndex++;
  for (Node *child : node->children) {
    buildOrPlace(
      method, data, seedVec, onSeeds, perNodeSeedMutations_Index, perNodeGapMutations_Index, seedK, seedS, T, child, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap
    );
  }


  for (const auto &p : seedChanges) {
    bool oldVal = std::get<1>(p);
    bool newVal = std::get<2>(p);
    int64_t pos = std::get<0>(p);
    if (oldVal && newVal) { // UNDO seed at same pos changed
      seedVec[pos] = true;
      onSeeds[pos] = std::get<3>(p);
    } else if (oldVal && !newVal) { // seed on to off
      seedVec[std::get<0>(p)] = true;
      onSeeds[std::get<0>(p)] = std::get<3>(p);
      int blockId = scalarCoordToBlockId[std::get<0>(p)];
      BlocksToSeeds[blockId].insert(std::get<0>(p));
    } else if (!oldVal && newVal) { // UNDO seed off to on
      seedVec[std::get<0>(p)] = false;
      if (onSeeds[std::get<0>(p)].has_value() && std::get<0>(p) < onSeeds.size()) {
        onSeeds[std::get<0>(p)].reset();
      }
      int blockId = scalarCoordToBlockId[std::get<0>(p)];
      BlocksToSeeds[blockId].erase(std::get<0>(p));
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
    }
  });

  std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());

  ::capnp::List<SeedMutations>::Builder perNodeSeedMutations_Builder = index.initPerNodeSeedMutations(T->allNodes.size());
  ::capnp::List<GapMutations>::Builder perNodeGapMutations_Builder = index.initPerNodeGapMutations(T->allNodes.size());

  int64_t dfsIndex = 0;

  std::vector<bool> seedVec(globalCoords.back().first.back().first + 1, false);
  std::vector<std::optional<std::string>> onSeeds(globalCoords.back().first.back().first + 1, std::nullopt);

  buildOrPlace(
    Step::BUILD, data, seedVec, onSeeds, perNodeSeedMutations_Builder, perNodeGapMutations_Builder, k, s, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap
  );
}

void pmi::place(Tree *T, Index::Reader &index)
{
    // Setup for seed indexing
    seed_annotated_tree::mutableTreeData data;
    seed_annotated_tree::globalCoords_t globalCoords;

    seed_annotated_tree::setup(data, globalCoords, T);
    
    CoordNavigator navigator(data.sequence);

    std::vector<int> BlockSizes(data.sequence.size(),0);
    std::vector<std::pair<int64_t, int64_t>> blockRanges(data.blockExists.size());

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
    }

    std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());

    tupleCoord_t coord = {0, 0, globalCoords[0].first[0].second.empty() ? -1 : 0};
    auto curIt = gapMap.end();

    while (coord < tupleCoord_t{-1, -1, -1}) {
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
        coord = navigator.newincrement(coord, data.blockStrand);
    }

    if (coord.blockId != -1 && !data.blockExists[coord.blockId].first) {
        coord = navigator.newdecrement(coord, data.blockStrand);
    }

    ::capnp::List<SeedMutations>::Reader perNodeSeedMutations_Reader= index.getPerNodeSeedMutations();
    ::capnp::List<GapMutations>::Reader perNodeGapMutations_Reader = index.getPerNodeGapMutations();

    int64_t dfsIndex = 0;

    std::vector<bool> seedVec(globalCoords.back().first.back().first + 1, false);
    std::vector<std::optional<std::string>> onSeeds(globalCoords.back().first.back().first + 1, std::nullopt);
    
    buildOrPlace<decltype(perNodeSeedMutations_Reader), decltype(perNodeGapMutations_Reader)>(
      Step::PLACE, data, seedVec, onSeeds, perNodeSeedMutations_Reader, perNodeGapMutations_Reader, k, s, T, T->root, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, dfsIndex, gapMap
    );
}

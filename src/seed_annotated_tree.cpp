#include "seed_annotated_tree.hpp"
#include "pmi.hpp"
#include "seeding.hpp"
#include <cmath>
#include <iostream>
#include <numeric>
#include <boost/icl/interval_set.hpp>
//#include <sys/_types/_int64_t.h>

using namespace boost::icl;
using namespace panman;
using namespace seed_annotated_tree;

std::chrono::time_point<std::chrono::high_resolution_clock> global_timer =
    std::chrono::high_resolution_clock::now();


double time_stamp() {
  std::chrono::time_point<std::chrono::high_resolution_clock> newtime =
      std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> duration = newtime - global_timer;
  global_timer = newtime;
  return duration.count();
}


std::string seed_annotated_tree::getConsensus(Tree *T) {
  TIME_FUNCTION;
  std::string consensus = "";
  for (size_t i = 0; i < T->blocks.size(); i++) {
    for (size_t j = 0; j < T->blocks[i].consensusSeq.size(); j++) {
      uint32_t c = T->blocks[i].consensusSeq[j];
      bool endFlag = false;
      for (size_t k = 0; k < 8; k++) {
        const int nucCode = (c >> (4 * (7 - k))) & 15;
        if (nucCode == 0) {
          endFlag = true;
          break;
        }
        consensus += getNucleotideFromCode(nucCode);
      }
      if (endFlag) {
        break;
      }
    }
  }
  return consensus;
}

bool seed_annotated_tree::getSeedAt(bool useHotSeeds, HotSeedIndex& hotSeedIndex, std::string &resultString, size_t &resultHash, 
  int64_t &resultEndPos, bool &resultIsReverse, const int64_t &pos, Tree *T, const int32_t& k,
  int64_t &dfsIndex,
    std::unordered_map<int64_t, tupleCoord_t> &scalarToTupleCoord,
    const sequence_t &sequence, const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const globalCoords_t &globalCoords, CoordNavigator &navigator,
    std::map<int64_t, int64_t> &gapRuns, const std::vector<std::pair<int64_t, int64_t>>& blockRanges) {
    

    if (useHotSeeds && hotSeedIndex.getSeed(pos, dfsIndex, resultHash, resultEndPos, resultIsReverse)) {
        std::cout << "DEBUG: Found seed in hot index" << std::endl;
        return true;
    }
    tupleCoord_t currCoord = scalarToTupleCoord[pos];
    if (!blockStrand[currCoord.blockId].first) {
      currCoord = scalarToTupleCoord[blockRanges[currCoord.blockId].first + blockRanges[currCoord.blockId].second - pos];
    }
    std::string seq;
    int32_t prevBlockId = currCoord.blockId;
    while (currCoord < tupleCoord_t{-1,-1,-1}){
      if (blockExists[currCoord.blockId].first) {
        char c = '-';
        if (currCoord.nucGapPos == -1){
          c = sequence[currCoord.blockId].first[currCoord.nucPos].first;
        }else{
          c = sequence[currCoord.blockId].first[currCoord.nucPos].second[currCoord.nucGapPos];
        }
        if(c == 'x') c = '-';
        if(! blockStrand[currCoord.blockId].first){
          switch(c){
            case 'A':
              c = 'T';
              break;
            case 'T':
              c = 'A';
              break;
            case 'G':
              c = 'C';
              break;
            case 'C':
              c = 'G';
              break;
            
            case 'Y':
              c = 'R';
              break;
            case 'R':
              c = 'Y';
              break;
            case 'K':
              c = 'M';
              break;
            case 'M':
              c = 'K';
              break;
            case 'D':
              c = 'H';
              break;
            case 'H':
              c = 'D';
              break;
            case 'V':
              c = 'B';
              break;
            case 'B':
              c = 'V';
              break;
          }
        }

        if(c != '-'){
          seq.push_back(c);
          if (seq.size() == k) {
            resultString = seq;
            resultEndPos = tupleToScalarCoord(currCoord, globalCoords);
            return true;
          }
        } else {
          int scalar = tupleToScalarCoord(currCoord, globalCoords);
          if (scalarToTupleCoord.find(gapRuns.find(scalar)->second) != scalarToTupleCoord.end()) {
            currCoord = scalarToTupleCoord[gapRuns.find(scalar)->second];
            if (!blockStrand[currCoord.blockId].first) {
              currCoord = scalarToTupleCoord[blockRanges[currCoord.blockId].first + blockRanges[currCoord.blockId].second - gapRuns.find(scalar)->second];
            }
            currCoord = navigator.newincrement(currCoord, blockStrand);
            prevBlockId = currCoord.blockId;
            continue;
          } else {
            currCoord = tupleCoord_t{-1, -1, -1};
          }
        }
      } else {
        if(blockStrand[currCoord.blockId].first){
            currCoord = tupleCoord_t{currCoord.blockId, navigator.sequence()[currCoord.blockId].first.size() - 1, -1};
          }else{
            currCoord.nucPos = 0;
            currCoord.nucGapPos = 0;
            if(navigator.sequence()[currCoord.blockId].first[0].second.empty()) {
              currCoord.nucGapPos = -1;
            }
        }
        if(currCoord.blockId == navigator.sequence().size() - 1){
          break;
        }
      }
      currCoord = navigator.newincrement(currCoord, blockStrand);
      if (prevBlockId != currCoord.blockId){
        if (!blockExists[currCoord.blockId].first) {
          int scalar = tupleToScalarCoord(currCoord, globalCoords);
          if (!blockStrand[currCoord.blockId].first) scalar = blockRanges[currCoord.blockId].first + blockRanges[currCoord.blockId].second - scalar;
          currCoord = scalarToTupleCoord[gapRuns.find(scalar)->second];
          if (!blockStrand[currCoord.blockId].first) {
            currCoord = scalarToTupleCoord[blockRanges[currCoord.blockId].first + blockRanges[currCoord.blockId].second - gapRuns.find(scalar)->second];
          }
          currCoord = navigator.newincrement(currCoord, blockStrand);
          prevBlockId = currCoord.blockId;
          continue;
        }
        if (blockStrand[currCoord.blockId].first) {
          currCoord = globalCoords[currCoord.blockId].first[0].second.empty() ? tupleCoord_t{currCoord.blockId, 0, -1} : tupleCoord_t{currCoord.blockId, 0, 0};
        } else {
          currCoord = tupleCoord_t{currCoord.blockId, globalCoords[currCoord.blockId].first.size() - 1, -1};
        }
      }
      prevBlockId = currCoord.blockId;
    }

    if (seq.size() != k) {
      throw std::runtime_error("unexpected kmer length during syncmer reconstruction");
    }
    resultString = seq;
    resultEndPos = tupleToScalarCoord(currCoord, globalCoords);
    return true;
}
  
// coords (blockId, nucPosition, nucGapPosition)
// Sequence includes both boundary coordinates
//Sequence, scalarCoords of sequence, scalarCoords of Gaps, dead blocks
std::tuple<std::string, std::vector<int>, std::vector<int>, std::vector<int>> seed_annotated_tree::getNucleotideSequenceFromBlockCoordinates(
    tupleCoord_t &start, tupleCoord_t &end, const sequence_t &sequence,
    const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const Tree *T, const Node *node, const globalCoords_t &globalCoords, CoordNavigator &navigator) {
    TIME_FUNCTION;
    if (end == tupleCoord_t{-1,-1,-1})
    {
      end = tupleCoord_t{navigator.sequence().size() - 1, navigator.sequence().back().first.size() - 1, -1};
    }

    long size = tupleToScalarCoord(end, globalCoords) - tupleToScalarCoord(start, globalCoords) + 1;
    if(size < 0){
      size = -size;
      auto temp = end;
      end = start;
      start = temp;
    }
    


    std::string seq;
    std::vector<int> coords;
    int i = 0;

    std::vector<int> gaps;
    int j = 0;


    std::vector<int> deadBlocks;



    // build sequence by iterating through {blockId, nucPosition, nucGapPosition} coords
    
    for(auto currCoord = start; currCoord != end; currCoord = navigator.newincrement(currCoord, blockStrand) ){
      

      if (blockExists[currCoord.blockId].first) {
        
        char c = '-';
        if (currCoord.nucGapPos == -1){
          
          c = navigator.sequence()[currCoord.blockId].first[currCoord.nucPos].first;
        }else{
          
          c = navigator.sequence()[currCoord.blockId].first[currCoord.nucPos].second[currCoord.nucGapPos];
        }
        
        
        if(c == 'x') c = '-';
        if(! blockStrand[currCoord.blockId].first){
          c = getComplementCharacter(c);
        }
      int scalar = tupleToScalarCoord(currCoord, globalCoords);

        if(c != '-'){
          seq.push_back(c);
          coords.push_back(scalar);
          i++;
        }else{
          gaps.push_back(scalar);
          j++;
        }
      }else{
        deadBlocks.push_back(currCoord.blockId);

        //jump to start of next block
        if( blockStrand[currCoord.blockId].first){
            //not inverted, jump to top of this block
            currCoord = tupleCoord_t{currCoord.blockId, navigator.sequence()[currCoord.blockId].first.size() - 1, -1};
          }else{
            //inverted, jump to bottom of this block
            //currCoord.blockId += 1;
            currCoord.nucPos = 0;
            currCoord.nucGapPos = 0;
            if(navigator.sequence()[currCoord.blockId].first[0].second.empty()) {
              currCoord.nucGapPos = -1;
            }
        }
        
        if(currCoord.blockId == navigator.sequence().size() - 1 || currCoord.blockId == end.blockId){
          break;
        }

      }

    }
    
    //Adding end character
    if(blockExists[end.blockId].first){
        char c = '-';
        if (end.nucGapPos == -1){
          c = navigator.sequence()[end.blockId].first[end.nucPos].first;
        } else {
          c = navigator.sequence()[end.blockId].first[end.nucPos].second[end.nucGapPos];
        }
        if(c == 'x') c = '-';
        if(!blockStrand[end.blockId].first){
          c = getComplementCharacter(c);
        }
        int scalar = tupleToScalarCoord(end, globalCoords);
        if(c != '-'){
          seq.push_back(c);
          coords.push_back(scalar);
          i++;
        } else {
          gaps.push_back(scalar);
          j++;
        }
    } else {
      deadBlocks.push_back(end.blockId);
    }

    return std::make_tuple(seq, coords, gaps, deadBlocks);
}

/*
- regap not used anywhere previously -> changed how regap is constructed for the
purpose of fixing syncmer/seedemr update bug and tracking end positions.
*/

void getAllStringsHelper(std::unordered_map<std::string, std::string> &strings,
                         mutableTreeData &data, Tree *T, const Node *node,
                         globalCoords_t &globalCoords) {

  // /*  Mutate with block and nuc mutations */
  // blockMutData_t blockMutData;
  // nucMutData_t nucMutData;
  // pmi::seedIndex blank;
  // blank.k = blank.s = blank.j = 0;
  // seedMap_t seedMap;
  // std::vector<tupleRange> recompRangesUnused;

  // applyMutations(data, seedMap, blockMutData, recompRangesUnused,
  //                 nucMutData, T, node, globalCoords, blank);
  // /* Use current state of mutableTreeData to decode node's sequence */
  // std::string seq =seed_annotated_tree::getStringFromCurrData(data, T, node, true);

  // strings[node->identifier] = seq;

  // /* Recursive step */
  // for (Node *child : node->children) {
  //   getAllStringsHelper(strings, data, T, child, globalCoords);
  // }

  // /* Undo mutations when backtracking */
  // undoMutations(data, blank, T, node, blockMutData, nucMutData);
}

std::unordered_map<std::string, std::string> seed_annotated_tree::getAllNodeStrings(Tree *T) {

  std::unordered_map<std::string, std::string> strings;
  seed_annotated_tree::mutableTreeData data;
  seed_annotated_tree::globalCoords_t globalCoords;
  setup(data, globalCoords, T);

  getAllStringsHelper(strings, data, T, T->root, globalCoords);

  return strings;
}
std::string seed_annotated_tree::getStringFromCurrData(mutableTreeData &data, Tree *T,
                                        const Node *node, const bool aligned) {

  // T should be const but [] operator on T->sequenceInverted is non-const
  std::string line;
  if (node == nullptr) { // consensus sequence (all blocks on) rather than a
                         // node in the tree
    for (size_t i = 0; i < T->blocks.size(); i++) {
      if (data.blockStrand[i].first) {
        for (size_t j = 0; j < data.sequence[i].first.size(); j++) {
          for (size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
            if (data.sequence[i].first[j].second[k] != '-') {
              line += data.sequence[i].first[j].second[k];
            } else if (aligned) {
              line += '-';
            }
          }
          if (data.sequence[i].first[j].first != '-' &&
              data.sequence[i].first[j].first != 'x') {
            line += data.sequence[i].first[j].first;
          } else if (aligned) {
            line += '-';
          }
        }
      } else {
        for (size_t j = data.sequence[i].first.size() - 1; j + 1 > 0; j--) {
          if (data.sequence[i].first[j].first != '-' &&
              data.sequence[i].first[j].first != 'x') {
            line += getComplementCharacter(data.sequence[i].first[j].first);
          } else if (aligned) {
            line += '-';
          }
          for (size_t k = data.sequence[i].first[j].second.size() - 1;
               k + 1 > 0; k--) {
            if (data.sequence[i].first[j].second[k] != '-') {
              line +=
                  getComplementCharacter(data.sequence[i].first[j].second[k]);
            } else if (aligned) {
              line += '-';
            }
          }
        }
      }
    }
    return line;
  }
  for (size_t i = 0; i < data.blockExists.size(); i++) {
    if (data.blockExists[i].first) {
      if (data.blockStrand[i].first) {
        for (size_t j = 0; j < data.sequence[i].first.size(); j++) {
          for (size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
            if (data.sequence[i].first[j].second[k] != '-') {
              line += data.sequence[i].first[j].second[k];
            } else if (aligned) {
              line += '-';
            }
          }
          if (data.sequence[i].first[j].first != '-' &&
              data.sequence[i].first[j].first != 'x') {
            line += data.sequence[i].first[j].first;
          } else if (aligned) {
            line += '-';
          }
        }
      } else {
        for (size_t j = data.sequence[i].first.size() - 1; j + 1 > 0; j--) {
          if (data.sequence[i].first[j].first != '-' &&
              data.sequence[i].first[j].first != 'x') {
            line += getComplementCharacter(data.sequence[i].first[j].first);
          } else if (aligned) {
            line += '-';
          }
          for (size_t k = data.sequence[i].first[j].second.size() - 1;
               k + 1 > 0; k--) {
            if (data.sequence[i].first[j].second[k] != '-') {
              line +=
                  getComplementCharacter(data.sequence[i].first[j].second[k]);
            } else if (aligned) {
              line += '-';
            }
          }
        }
      }
    } else if (aligned) {
      for (size_t j = 0; j < data.sequence[i].first.size(); j++) {
        for (size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
          line += '-';
        }
        line += '-';
      }
    }
  }

  return line;
}
void seed_annotated_tree::setupGlobalCoordinates(int64_t &ctr, globalCoords_t &globalCoords,
                                  const BlockGapList &blockGaps,
                                  const std::vector<Block> &blocks,
                                  const std::vector<GapList> &gaps,
                                  const sequence_t &sequence,
                                  std::unordered_map<int64_t, tupleCoord_t> &scalarToTupleCoord) {
    TIME_FUNCTION;
    // Add size validation
    if (blocks.empty()) {
        throw std::runtime_error("No blocks in tree");
    }
    
    int maxBlockId = 0;
    for (const auto& block : blocks) {
        maxBlockId = std::max(maxBlockId, static_cast<int>(block.primaryBlockId));
    }
    
    // Ensure globalCoords size matches sequence
    if (maxBlockId >= sequence.size()) {
        throw std::runtime_error("Block ID exceeds sequence size");
    }
    
    globalCoords.resize(maxBlockId + 1);

    // Assigning block gaps
    for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
        globalCoords[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
    }

    // Initialize first vectors based on sequence sizes
    for (size_t i = 0; i < sequence.size(); i++) {
        globalCoords[i].first.resize(sequence[i].first.size());
    }

    // Now safe to assign nucleotide gaps
    for (size_t i = 0; i < gaps.size(); i++) {
        int32_t blockId = gaps[i].primaryBlockId;
        for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            globalCoords[blockId].first[pos].second.resize(len, 0);
        }
    }

    // Assigning coordinates
    ctr = 0;
    for (size_t i = 0; i < globalCoords.size(); i++) {
        for (size_t j = 0; j < globalCoords[i].first.size(); j++) {
            for (size_t k = 0; k < globalCoords[i].first[j].second.size(); k++) {
                globalCoords[i].first[j].second[k] = ctr;
                scalarToTupleCoord[ctr] = {i, j, k};
                ctr++;
            }
            globalCoords[i].first[j].first = ctr;
            scalarToTupleCoord[ctr] = {i, j, -1};
            ctr++;
        }
    }
}
void seed_annotated_tree::setup(mutableTreeData &data, globalCoords_t &globalCoords,
                 const Tree *T) {
    TIME_FUNCTION;
    const BlockGapList &blockGaps = T->blockGaps;
    const std::vector<GapList> &gaps = T->gaps;
    const std::vector<Block> &blocks = T->blocks;

    sequence_t sequence(blocks.size() + 1);
    blockExists_t blockExists(blocks.size() + 1, {false, {}});
    blockStrand_t blockStrand(
        blocks.size() + 1,
        {true,
         {}}); // TODO: strand is always tru for now and also not used anywhere
               // and also not set anywhere and also not used in the code and also
               // not set in the code and n  also not used in the code
    int32_t maxBlock = 0;

    for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
      sequence[blockGaps.blockPosition[i]].second.resize(
          blockGaps.blockGapLength[i]);
      blockExists[blockGaps.blockPosition[i]].second.resize(
          blockGaps.blockGapLength[i], false);
      blockStrand[blockGaps.blockPosition[i]].second.resize(
          blockGaps.blockGapLength[i], true);
    }

    for (size_t i = 0; i < blocks.size(); i++) {
      int32_t b = ((int32_t)blocks[i].primaryBlockId);
      maxBlock = std::max(maxBlock, b);
      for (size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
        bool stop = false;
        for (size_t k = 0; k < 8; k++) {
          const int nc = (((blocks[i].consensusSeq[j]) >> (4 * (7 - k))) & 15);
          if (nc == 0) {
            stop = true;
            break;
          }
          const char c = panmanUtils::getNucleotideFromCode(nc);
          sequence[b].first.push_back({c, {}});
        }
        if (stop) {
          break;
        }
      }
      sequence[b].first.push_back({'x', {}});
    }

    sequence.resize(maxBlock + 1);
    blockExists.resize(maxBlock + 1);
    blockStrand.resize(maxBlock + 1);

    // Assigning nucleotide gaps in blocks
    for (size_t i = 0; i < gaps.size(); i++) {
      int32_t id = (gaps[i].primaryBlockId);
      for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
        int len = gaps[i].nucGapLength[j];
        int pos = gaps[i].nucPosition[j];
        sequence[id].first[pos].second.resize(len, '-');
      }
    }

    data.sequence = sequence;
    data.blockExists = blockExists;
    data.blockStrand = blockStrand;
    setupGlobalCoordinates(data.maxGlobalCoordinate, globalCoords, blockGaps,
                           blocks, gaps, sequence, data.scalarToTupleCoord);
}
int64_t seed_annotated_tree::getGlobalCoordinate(const int blockId, const int nucPosition,
                                  const int nucGapPosition,
                                  const globalCoords_t &globalCoords) {
    TIME_FUNCTION;
    if (blockId > globalCoords.size() - 1 ||
        blockId == globalCoords.size() - 1 &&
            nucPosition > globalCoords.back().first.size() - 1) {
      return globalCoords[globalCoords.size() - 1]
          .first[globalCoords.back().first.size() - 1]
          .first;
    }

    if (nucGapPosition == -1) {
      return globalCoords[blockId].first[nucPosition].first;
    }

    if (globalCoords[blockId].first[nucPosition].second.size() == 0) {
      return globalCoords[blockId].first[nucPosition].first;
    }

    return globalCoords[blockId].first[nucPosition].second[nucGapPosition];
}
static int getIndexFromNucleotide(char nuc) {

  switch (nuc) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  case '*':
    return 4;
  default:
    return 5;
  }
  return 5;
}
static size_t getBeg(const std::string &s1, const std::string &s2,
                     size_t window, double threshold) {

  if (s1.empty()) {
    return 0;
  }

  size_t numAlign = 0;
  size_t numMatch = 0;
  size_t beg = 0;
  size_t idx = 0;
  std::queue<size_t> begs;
  while (idx < s1.size()) {
    if (s1[idx] == '-' && s2[idx] == '-') {
      ++idx;
      continue;
    }
    if (beg == 0) {
      beg = idx;
    } else {
      begs.push(idx);
    }
    if (s1[idx] == s2[idx]) {
      ++numMatch;
    }
    ++numAlign;
    ++idx;

    if (numAlign == window) {
      double pcid = static_cast<double>(numMatch) / static_cast<double>(window);
      if (pcid >= threshold && s1[beg] == s2[beg]) {
        return beg;
      }

      if (s1[beg] == s2[beg]) {
        --numMatch;
      }
      --numAlign;
      beg = begs.front();
      begs.pop();
    }
  }

  return s1.size() - 1;
}
static size_t getEnd(const std::string &s1, const std::string &s2,
                     size_t window, double threshold) {

  if (s1.empty()) {
    return 0;
  }

  size_t numAlign = 0;
  size_t numMatch = 0;
  size_t end = s1.size();
  size_t idx = s1.size() - 1;
  std::queue<size_t> ends;

  while (true) {
    if (s1[idx] == '-' && s2[idx] == '-') {
      if (idx == 0) {
        break;
      }
      --idx;
      continue;
    }
    if (end == s1.size()) {
      end = idx;
    } else {
      ends.push(idx);
    }
    if (s1[idx] == s2[idx]) {
      ++numMatch;
    }
    ++numAlign;

    if (idx == 0) {
      break;
    }
    --idx;

    if (numAlign == window) {
      double pcid = static_cast<double>(numMatch) / static_cast<double>(window);
      if (pcid >= threshold && s1[end] == s2[end]) {
        return end;
      }

      if (s1[end] == s2[end]) {
        --numMatch;
      }
      --numAlign;
      end = ends.front();
      ends.pop();
    }
  }

  return 0;
}
std::pair<size_t, size_t>seed_annotated_tree::getMaskCoorsForMutmat(const std::string &s1,
                                                      const std::string &s2,
                                                      size_t window,
                                                      double threshold) {

  assert(s1.size() == s2.size());
  if (window == 0 || threshold == 0.0) {
    return std::make_pair<size_t, size_t>(0, s1.size() - 1);
  }
  return std::make_pair<size_t, size_t>(getBeg(s1, s2, window, threshold),
                                        getEnd(s1, s2, window, threshold));
}

void seed_annotated_tree::writeMutationMatrices(const mutationMatrices &mutMat,
                                 std::ofstream &mmfout) {
    TIME_FUNCTION;
    for (const std::vector<double> &row : mutMat.submat) {
      for (const double &prob : row) {
        mmfout << prob << " ";
      }
      mmfout << "\n";
    }

    for (const auto& [size, count] : mutMat.insmat) {
      mmfout << size << ":" << count << " ";
    }
    mmfout << "\n";

    for (const auto& [size, count] : mutMat.delmat) {
      mmfout << size << ":" << count << " ";
    }
    mmfout << "\n";
}

// void buildMutationMatricesHelper(mutationMatrices &mutMat, Tree *T, Node* node, size_t window,
//                            double threshold) {
//   if (node->parent != nullptr) {
//     std::string curSeq = getStringAtNode(node, T, true);
//     std::vector<int64_t> curBaseCounts(4, 0);
//     for (size_t i = 0; i < curSeq.size(); i++) {
//       char c = curSeq[i];
//       if (getIndexFromNucleotide(c) < 4) {
//         curBaseCounts[getIndexFromNucleotide(c)]++;
//       }
//     }
//     std::cout << node->identifier << " true base counts: " << curBaseCounts[0] << " " << curBaseCounts[1] << " " << curBaseCounts[2] << " " << curBaseCounts[3] << std::endl;
//     std::string parSeq = getStringAtNode(node->parent, T, true);
//     size_t insLen = 0;
//     size_t delLen = 0;
//     std::pair<size_t, size_t> edgeCoor =seed_annotated_tree::getMaskCoorsForMutmat(curSeq, parSeq, window, threshold);
//     if (edgeCoor.second > edgeCoor.first) {
//       for (size_t i = edgeCoor.first; i < edgeCoor.second + 1; i++) {
//         if (parSeq[i] == '-' && curSeq[i] == '-') {
//           continue;
//         } else if (parSeq[i] != '-' && curSeq[i] == '-') {
//           delLen++;
//           if (insLen > 0) {
//             if (insLen > mutMat.insmat.size() - 1) {
//               mutMat.insmat.resize(insLen + 1);
//             }
//             mutMat.insmat[insLen]++;
//             mutMat.total_insmut++;
//             insLen = 0;
//           }
//         } else if (parSeq[i] == '-' && curSeq[i] != '-') {
//           insLen++;
//           if (delLen > 0) {
//             if (delLen > mutMat.delmat.size() - 1) {
//               mutMat.delmat.resize(delLen + 1);
//             }
//             mutMat.delmat[delLen]++;
//             mutMat.total_delmut++;
//             delLen = 0;
//           }
//         } else {
//           if (insLen > mutMat.insmat.size() - 1) {
//             mutMat.insmat.resize(insLen + 1);
//           }
//           mutMat.insmat[insLen]++;
//           mutMat.total_insmut++;
//           insLen = 0;

//           if (delLen > mutMat.delmat.size() - 1) {
//             mutMat.delmat.resize(delLen + 1);
//           }
//           mutMat.delmat[delLen]++;
//           mutMat.total_delmut++;
//           delLen = 0;

//           int parNucIdx = getIndexFromNucleotide(parSeq[i]);
//           int curNucIdx = getIndexFromNucleotide(curSeq[i]);
//           if (parNucIdx > 3 || curNucIdx > 3) {
//             continue;
//           }
//           mutMat.submat[parNucIdx][curNucIdx]++;
//           mutMat.total_submuts[parNucIdx]++;
//         }
//       }
//     }
//     curSeq.clear();
//     parSeq.clear();
//   }


//   for (Node *child : node->children) {
//     buildMutationMatricesHelper(mutMat, T, child, window, threshold);
//   }

// }

// voidseed_annotated_tree::fillMutationMatricesFromTree(mutationMatrices &mutMat, Tree *T,
//                                         size_t window, double threshold) {
//   buildMutationMatricesHelper(mutMat, T, T->root, window, threshold); 

//   for (auto i = 0; i < 4; i++) {
//     for (auto j = 0; j < 4; j++) {
//       std::cout << mutMat.submat[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }

//   for (auto i = 0; i < mutMat.insmat.size(); ++i) {
//     std::cout << mutMat.insmat[i] << " ";
//   }
//   std::cout << std::endl;
//   // deletion
//   for (auto i = 0; i < mutMat.delmat.size(); ++i) {
//     std::cout << mutMat.delmat[i] << " ";
//   }
//   std::cout << std::endl;

//   // insertion
//   for (auto i = 0; i < mutMat.insmat.size(); ++i) {
//     mutMat.insmat[i] = -10 * log10f(mutMat.insmat[i] / mutMat.total_insmut);
//   }
//   // deletion
//   for (auto i = 0; i < mutMat.delmat.size(); ++i) {
//     mutMat.delmat[i] = -10 * log10f(mutMat.delmat[i] / mutMat.total_delmut);
//   }
//   // substitution
//   for (auto i = 0; i < 4; i++) {
//     for (auto j = 0; j < 4; j++) {
//       mutMat.submat[i][j] =
//           -10 * log10f(mutMat.submat[i][j] / mutMat.total_submuts[i]);
//     }
//   }
//   mutMat.filled = true;
// }

void clearBlockBaseCounts(
  int32_t primaryBlockId, mutableTreeData &data, globalCoords_t &globalCoords, CoordNavigator &navigator,
  std::vector<int64_t> &curBaseCounts, std::vector<int64_t> &parentBaseCountsBacktrack) {

  tupleCoord_t coord = globalCoords[primaryBlockId].first[0].second.empty() ? tupleCoord_t{primaryBlockId, 0, -1} : tupleCoord_t{primaryBlockId, 0, 0};
  tupleCoord_t end = tupleCoord_t{primaryBlockId, globalCoords[primaryBlockId].first.size() - 1, -1};
  if (!data.blockStrand[primaryBlockId].first) std::swap(coord, end);
  while (true) {
    char c = coord.nucGapPos == -1 ? data.sequence[coord.blockId].first[coord.nucPos].first : data.sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
    if (!data.blockStrand[primaryBlockId].first) c = getComplementCharacter(c);
    int64_t scalar = tupleToScalarCoord(coord, globalCoords);

    if (getIndexFromNucleotide(c) <= 3) {
      --curBaseCounts[getIndexFromNucleotide(c)];
      ++parentBaseCountsBacktrack[getIndexFromNucleotide(c)];

    }
    if (coord == end) break;
    coord = navigator.newincrement(coord, data.blockStrand);
  }
}

static void applyMutations(mutableTreeData &data,
                    blockMutationInfo_t &blockMutationInfo,
                    mutationInfo_t &mutationInfo, Tree *T, Node *node,
                    globalCoords_t &globalCoords,
                    CoordNavigator &navigator, const std::vector<std::pair<int64_t, int64_t>> &blockRanges,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
                    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks,
                    blockExists_t& oldBlockExists, blockStrand_t& oldBlockStrand,
                    std::vector<int64_t> &curBaseCounts, std::vector<int64_t> &parentBaseCountsBacktrack,
                    std::vector<std::vector<int64_t>> &subCount, bool& isMutatedNuc)
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
          clearBlockBaseCounts(primaryBlockId, data, globalCoords, navigator, curBaseCounts, parentBaseCountsBacktrack);
          
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
          clearBlockBaseCounts(primaryBlockId, data, globalCoords, navigator, curBaseCounts, parentBaseCountsBacktrack);
          oldStrand = blockStrand[primaryBlockId].first;
          oldMut = blockExists[primaryBlockId].first;
          blockExists[primaryBlockId].first = false;
          // resetting strand to true during deletion
          blockStrand[primaryBlockId].first = true;
        }
      }
      blockMutationInfo.push_back(std::make_tuple(mutation.primaryBlockId, mutation.secondaryBlockId, oldMut, oldStrand, false, true));
    }

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
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
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
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
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
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
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
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition].second[nucGapPosition + j] = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition, nucGapPosition + j, oldVal, newVal));
            }
          }
          else
          {
            for (int j = 0; j < len; j++)
            {
              char oldVal = sequence[primaryBlockId].first[nucPosition + j].first;
              newVal = panmanUtils::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
              sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              mutationInfo.push_back(std::make_tuple(primaryBlockId, secondaryBlockId, nucPosition + j, nucGapPosition, oldVal, newVal));
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
      else if (type == panmanUtils::NucMutationType::NSNPD)
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
    if (oldBlockExists[blockId].first && blockExists[blockId].first && oldBlockStrand[blockId].first == blockStrand[blockId].first) {
      // on to on -> collect gap runs and nuc runs
      int nucPos = std::get<2>(mutation);
      int nucGapPos = std::get<3>(mutation);
      char parChar = std::get<4>(mutation) == 'x' ? '-' : std::get<4>(mutation);
      char curChar = std::get<5>(mutation) == 'x' ? '-' : std::get<5>(mutation);

      int64_t scalar = tupleToScalarCoord(tupleCoord_t{blockId, nucPos, nucGapPos}, globalCoords);
      if (!data.blockStrand[blockId].first) {
        scalar = blockRanges[blockId].first + blockRanges[blockId].second - scalar;
      }

      if (!oldBlockStrand[blockId].first) parChar = getComplementCharacter(parChar);
      if (!blockStrand[blockId].first) curChar =  getComplementCharacter(curChar);

      if (getIndexFromNucleotide(parChar) <= 3) {
        --curBaseCounts[getIndexFromNucleotide(parChar)];
        ++parentBaseCountsBacktrack[getIndexFromNucleotide(parChar)];
      }
      if (getIndexFromNucleotide(curChar) <= 3) {
        ++curBaseCounts[getIndexFromNucleotide(curChar)];
        --parentBaseCountsBacktrack[getIndexFromNucleotide(curChar)];
      }
      if (getIndexFromNucleotide(parChar) <= 3 && getIndexFromNucleotide(curChar) <= 3) {
        isMutatedNuc = true;
        ++subCount[getIndexFromNucleotide(parChar)][getIndexFromNucleotide(curChar)];
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

static void undoMutations(mutableTreeData &data, Tree *T,
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

interval_set<int64_t> gapMapToNucRunSet(const std::map<int64_t, int64_t> &gapMap, const std::vector<std::pair<int64_t, int64_t>>& blockRanges) {
  interval_set<int64_t> curNucRunSet;
  int64_t start = -1;
  for (const auto& [curStart, curEnd] : gapMap) {
    if (start == -1) {
      if (curStart != 0) {
        curNucRunSet.add(interval<int64_t>::closed(0, curStart - 1));
      }
    } else {
      curNucRunSet.add(interval<int64_t>::closed(start, curStart - 1));
    }
    start = curEnd + 1;
  }
  if (start <= blockRanges.back().second) {
    curNucRunSet.add(interval<int64_t>::closed(start, blockRanges.back().second));
  }
  return curNucRunSet;
}

void buildMutationMatricesHelper_test(
  mutationMatrices &mutMat, Tree *T, Node* node, mutableTreeData &data, std::map<int64_t, int64_t> &gapMap,
  globalCoords_t &globalCoords, CoordNavigator &navigator, std::vector<int> &scalarCoordToBlockId,
  std::vector<std::unordered_set<int>> &BlocksToSeeds, std::vector<int> &BlockSizes,
  const std::vector<std::pair<int64_t, int64_t>>& blockRanges,
  std::vector<int64_t> &parentBaseCounts, std::vector<int64_t> &totalBaseCounts,
  std::vector<std::vector<int64_t>> &subCount, std::unordered_map<int64_t, int64_t> &insCount, std::unordered_map<int64_t, int64_t> &delCount
) {
  blockMutationInfo_t blockMutationInfo;
  mutationInfo_t mutationInfo;

  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapMapUpdates;

  blockExists_t oldBlockExists;
  blockStrand_t oldBlockStrand;

  std::vector<int64_t> curBaseCounts = parentBaseCounts;
  std::vector<int64_t> parentBaseCountsBacktrack(4);
  bool isMutatedNuc = false;

  applyMutations(data, blockMutationInfo, mutationInfo, T, node, globalCoords,
    navigator, blockRanges, gapRunUpdates, gapRunBacktracks, oldBlockExists,
    oldBlockStrand, curBaseCounts, parentBaseCountsBacktrack, subCount, isMutatedNuc);


  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBlocksBacktracks;

  auto parentGapMap = gapMap;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> tmpGapRunBlocksBacktracks;
  for (size_t i = 0; i < oldBlockStrand.size(); ++i) {
    if (!oldBlockStrand[i].first) {
      invertGapMap(parentGapMap, blockRanges[i], tmpGapRunBlocksBacktracks, gapMapUpdates);
    }
  }
  interval_set<int64_t> parNucRunSet = gapMapToNucRunSet(parentGapMap, blockRanges);
  interval_set<int64_t> flippedSet;
  parentGapMap.clear();
  tmpGapRunBlocksBacktracks.clear();

  updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);
  std::vector<int64_t> invertedBlocks;
  for (int i = 0; i < data.blockExists.size(); i++) {
    const bool& oldExists = oldBlockExists[i].first;
    const bool& newExists = data.blockExists[i].first;
    if (newExists && !data.blockStrand[i].first) {
      invertedBlocks.push_back(i);
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
      std::pair<int64_t, int64_t> curNucRange = {-1, -1};
      std::vector<std::pair<int64_t, int64_t>> nucRanges;
      while (true) {
        char c = coord.nucGapPos == -1 ? data.sequence[coord.blockId].first[coord.nucPos].first : data.sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
        c = c == 'x' ? '-' : c;
        if (!data.blockStrand[i].first) c = getComplementCharacter(c);

        int64_t scalar = tupleToScalarCoord(coord, globalCoords);

        if (c != '-') {
          if (getIndexFromNucleotide(c) <= 3) {
            ++curBaseCounts[getIndexFromNucleotide(c)];
            --parentBaseCountsBacktrack[getIndexFromNucleotide(c)];
          }

          if (curNucRange.first != -1 && curNucRange.second + 1 == scalar) {
            ++curNucRange.second;
          } else {
            if (curNucRange.first != -1) {
              nucRanges.push_back(curNucRange);
            }
            curNucRange = {scalar, scalar};
          }
        }

        if (coord == end) break;
        coord = navigator.newincrement(coord, data.blockStrand);
      }

      if (curNucRange.first != -1) {
        nucRanges.push_back(curNucRange);
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
    } else if (oldExists && newExists && (data.blockStrand[i].first != oldBlockStrand[i].first)) {
      // on to on -> but strand flipped
      flippedSet.add(interval<int64_t>::closed(blockRanges[i].first, blockRanges[i].second));

      tupleCoord_t coord = globalCoords[i].first[0].second.empty() ? tupleCoord_t{i, 0, -1} : tupleCoord_t{i, 0, 0};
      tupleCoord_t end = tupleCoord_t{i, globalCoords[i].first.size() - 1, -1};
      if (!data.blockStrand[i].first) std::swap(coord, end);
      while (true) {
        char c = coord.nucGapPos == -1 ? data.sequence[coord.blockId].first[coord.nucPos].first : data.sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
        c = c == 'x' ? '-' : c;
        if (!data.blockStrand[i].first) c = getComplementCharacter(c);
        int64_t scalar = tupleToScalarCoord(coord, globalCoords);

        if (getIndexFromNucleotide(c) <= 3) {
          ++curBaseCounts[getIndexFromNucleotide(c)];
          --parentBaseCountsBacktrack[getIndexFromNucleotide(c)];
        }

        if (coord == end) break;
        coord = navigator.newincrement(coord, data.blockStrand);
      }
    }
  }

  for (const auto& blockId : invertedBlocks) {
    invertGapMap(gapMap, blockRanges[blockId], gapRunBlocksBacktracks, gapMapUpdates);
  }

  if (node->parent != nullptr) {
    // make nuc run interval set from gap map
    interval_set<int64_t> curNucRunSet = gapMapToNucRunSet(gapMap, blockRanges);

    // remove flipped blocks from nucRunSet so they are not counted as indels
    curNucRunSet -= flippedSet;
    parNucRunSet -= flippedSet;

    // XOR ranges between cur and par
    interval_set<int64_t> xorSet = curNucRunSet ^ parNucRunSet;

    // AND ranges bewteen xor and par -> deletion
    interval_set<int64_t> deletionSet = xorSet & parNucRunSet;

    // AND ranges between xor and cur -> insertion
    interval_set<int64_t> insertionSet = xorSet & curNucRunSet;

    // std::cout << node->identifier << " test base counts: " << curBaseCounts[0] << " " << curBaseCounts[1] << " " << curBaseCounts[2] << " " << curBaseCounts[3] << std::endl;
    
    if (isMutatedNuc || !insertionSet.empty() || !deletionSet.empty()) {
      for (const auto& interval : insertionSet) {
        int64_t size = boost::icl::last(interval) - boost::icl::first(interval) + 1;
        assert(size > 0);
        if (insCount.find(size) == insCount.end()) {
          insCount[size] = 0;
        }
        ++insCount[size];
      }

      for (const auto& interval : deletionSet) {
        int64_t size = boost::icl::last(interval) - boost::icl::first(interval) + 1;
        assert(size > 0);
        if (delCount.find(size) == delCount.end()) {
          delCount[size] = 0;
        }
        ++delCount[size];
      }

      for (size_t i = 0; i < 4; ++i) {
        totalBaseCounts[i] += curBaseCounts[i];
      }
    }
  }

  parentBaseCounts = curBaseCounts;


  for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  gapRunUpdates.clear();

  for (Node *child : node->children) {
    buildMutationMatricesHelper_test(mutMat, T, child, data, gapMap, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, parentBaseCounts, totalBaseCounts, subCount, insCount, delCount);
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

  // undo base counts
  for (int i = 0; i < 4; i++) {
    parentBaseCounts[i] += parentBaseCountsBacktrack[i];
  }

  undoMutations(data, T, node, blockMutationInfo, mutationInfo, globalCoords);
}

void seed_annotated_tree::fillMutationMatricesFromTree_test(
  mutationMatrices &mutMat, Tree *T, const std::string& path
) {
    TIME_FUNCTION;
    seed_annotated_tree::mutableTreeData data;
    seed_annotated_tree::globalCoords_t globalCoords;
    seed_annotated_tree::setup(data, globalCoords, T);

    CoordNavigator navigator(data.sequence);

    std::vector<int> BlockSizes(data.sequence.size(),0);
    std::vector<std::pair<int64_t, int64_t>> blockRanges(data.blockExists.size());
    
    std::vector<int> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
    auto currCoord = tupleCoord_t{0,0,0};
    if(navigator.sequence()[0].first[0].second.empty()) {
      currCoord.nucGapPos = -1;
    }

    for(int i = 0; i < scalarCoordToBlockId.size(); i++){
      scalarCoordToBlockId[i] = currCoord.blockId;
      BlockSizes[currCoord.blockId] ++; 
      currCoord = navigator.newincrement(currCoord, data.blockStrand);
    }

    for (size_t i = 0; i < blockRanges.size(); ++i) {
      int64_t start = globalCoords[i].first[0].second.empty() ? tupleToScalarCoord({i, 0, -1}, globalCoords) : tupleToScalarCoord({i, 0, 0}, globalCoords);
      int64_t end = tupleToScalarCoord({i, globalCoords[i].first.size() - 1, -1}, globalCoords);
      blockRanges[i] = std::make_pair(start, end);
    }

    std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());
    std::map<int64_t, int64_t> gapMap;
    gapMap[0] = tupleToScalarCoord({blockRanges.size() - 1, globalCoords[blockRanges.size() - 1].first.size() - 1, -1}, globalCoords);
    
    std::vector<int64_t> parentBaseCounts(4);
    std::vector<int64_t> totalBaseCounts(4);
    std::vector<std::vector<int64_t>> subCount(4, std::vector<int64_t>(4, 0));
    std::unordered_map<int64_t, int64_t> insCount;
    std::unordered_map<int64_t, int64_t> delCount;

    // update baseCounts
    //   update nucMut when block on -> on
    //   recompute when block off->on or on->off
    // update totalBaseCounts
    //   add baseCounts to totalbaseCounts if there is any indel or substitution
    //   do nothing if there is no indel or substitution
    buildMutationMatricesHelper_test(mutMat, T, T->root, data, gapMap, globalCoords, navigator, scalarCoordToBlockId, BlocksToSeeds, BlockSizes, blockRanges, parentBaseCounts, totalBaseCounts, subCount, insCount, delCount);

    int64_t totalNucCounts = 0;
    int64_t totalInsCounts = 0;
    int64_t totalDelCounts = 0;
    for (const auto& count : totalBaseCounts) totalNucCounts += count;
    for (const auto& [size, count] : insCount) totalInsCounts += count;
    for (const auto& [size, count] : delCount) totalDelCounts += count;
    insCount[0] = totalNucCounts - totalInsCounts;
    delCount[0] = totalNucCounts - totalDelCounts;

    for (int i = 0; i < 4; ++i) {
      mutMat.submat[i][i] = static_cast<double>(totalBaseCounts[i]);
    }
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (i != j) {
          mutMat.submat[i][j] = static_cast<double>(subCount[i][j]);
          mutMat.submat[i][i] -= mutMat.submat[i][j];
        }
      }
    }

    for (const auto& [size, count] : insCount) {
      mutMat.insmat[size] = static_cast<double>(count);
    }
    for (const auto& [size, count] : delCount) {
      mutMat.delmat[size] = static_cast<double>(count);
    }

    // insertion
    for (auto [size, count] : mutMat.insmat) {
      mutMat.insmat[size] = -10 * log10f(count / static_cast<double>(totalNucCounts));
    }
    // deletion
    for (auto [size, count] : mutMat.delmat) {
      mutMat.delmat[size] = -10 * log10f(count / static_cast<double>(totalNucCounts));
    }

    // substitution
    for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) {
        mutMat.submat[i][j] = -10 * log10f(mutMat.submat[i][j] / static_cast<double>(totalBaseCounts[i]));
      }
    }
    mutMat.filled = true;

    std::ofstream outFile(path);
    for (auto i = 0; i < 4; i++) {
      for (auto j = 0; j < 4; j++) {
        outFile << mutMat.submat[i][j] << " ";
      }
      outFile << "\n";
    }

    std::vector<std::pair<int64_t, double>> insmat_sorted(mutMat.insmat.begin(), mutMat.insmat.end());
    std::sort(insmat_sorted.begin(), insmat_sorted.end(), [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
      return a.first < b.first;
    });
    for (const auto& [size, logProb] : insmat_sorted) {
      outFile << size << ":" << logProb << " ";
    }
    outFile << "\n";

    std::vector<std::pair<int64_t, double>> delmat_sorted(mutMat.delmat.begin(), mutMat.delmat.end());
    std::sort(delmat_sorted.begin(), delmat_sorted.end(), [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
      return a.first < b.first;
    });
    for (const auto& [size, logProb] : delmat_sorted) {
      outFile << size << ":" << logProb << " ";
    }
    outFile << "\n";
    outFile.close();
}

void seed_annotated_tree::fillMutationMatricesFromFile(mutationMatrices &mutMat,
                                        std::ifstream &inf) {
    TIME_FUNCTION;
    std::string line;
    int idx = 0;
    while (getline(inf, line)) {
      std::vector<double> probs;
      std::vector<std::string> fields;
      stringSplit(line, ' ', fields);
      for (const auto &f : fields) {
        probs.push_back(std::stod(f));
      }
      if (probs.size() == 0) {
        break;
      }

      if (idx < 4) {
        if (probs.size() != 4) {
          throw std::invalid_argument(
              "Received invalid mutation matrix (.mm) file");
        }

        std::vector<double> probs;
        for (const auto &f : fields) {
          probs.push_back(std::stod(f));
        }
        mutMat.submat[idx] = std::move(probs);
      } else if (idx == 4) {
        if (probs.size() < 1) {
          throw std::invalid_argument(
              "Received invalid mutation matrix (.mm) file");
        }

        for (const auto& f : fields) {
          std::vector<std::string> subFields;
          stringSplit(f, ':', subFields);
          int64_t size = std::stoll(subFields[0]);
          double prob = std::stod(subFields[1]);
          mutMat.insmat[size] = prob;
        }
      } else if (idx == 5) {
        if (probs.size() < 1) {
          throw std::invalid_argument(
              "Received invalid mutation matrix (.mm) file");
        }

      for (const auto& f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.delmat[size] = prob;
      }
    }
    idx++;
  }

  if (idx != 6) {
    throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
  }
  mutMat.filled = true;
}

std::string seed_annotated_tree::getStringAtNode(Node *node, Tree *T, bool aligned) {

  Node *referenceNode = nullptr;

  for (auto u : T->allNodes) {
    if (u.first == node->identifier) {
      referenceNode = u.second;
      break;
    }
  }
  auto &allNodes = T->allNodes;
  auto &blocks = T->blocks;
  auto &gaps = T->gaps;
  std::string reference = node->identifier;
  auto &blockGaps = T->blockGaps;
  auto &sequenceInverted = T->sequenceInverted;
  auto &root = T->root;
  auto &circularSequences = T->circularSequences;

  if (referenceNode == nullptr) {
    return "Error: Reference sequence with matching name not found!";
  }

  std::vector<panmanUtils::Node *> path;
  Node *it = referenceNode;

  while (it != root) {
    path.push_back(it);
    it = it->parent;
  }
  path.push_back(root);

  // List of blocks. Each block has a nucleotide list. Along with each
  // nucleotide is a gap list.
  std::vector<
      std::pair<std::vector<std::pair<char, std::vector<char>>>,
                std::vector<std::vector<std::pair<char, std::vector<char>>>>>>
      sequence(blocks.size() + 1);
  std::vector<std::pair<bool, std::vector<bool>>> blockExists(blocks.size() + 1,
                                                              {false, {}});
  blockStrand_t blockStrand(blocks.size() + 1, {true, {}});

  // Assigning block gaps
  for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
    sequence[blockGaps.blockPosition[i]].second.resize(
        blockGaps.blockGapLength[i]);
    blockExists[blockGaps.blockPosition[i]].second.resize(
        blockGaps.blockGapLength[i], false);
    blockStrand[blockGaps.blockPosition[i]].second.resize(
        blockGaps.blockGapLength[i], true);
  }

  int32_t maxBlockId = 0;

  // Create block consensus sequences
  for (size_t i = 0; i < blocks.size(); i++) {
    int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
    int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);
    maxBlockId = std::max(maxBlockId, primaryBlockId);

    for (size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
      bool endFlag = false;
      for (size_t k = 0; k < 8; k++) {
        const int nucCode =
            (((blocks[i].consensusSeq[j]) >> (4 * (7 - k))) & 15);

        if (nucCode == 0) {
          endFlag = true;
          break;
        }
        const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

        if (secondaryBlockId != -1) {
          sequence[primaryBlockId].second[secondaryBlockId].push_back(
              {nucleotide, {}});
        } else {
          sequence[primaryBlockId].first.push_back({nucleotide, {}});
        }
      }
      if (endFlag) {
        break;
      }
    }

    // End character to incorporate for gaps at the end
    if (secondaryBlockId != -1) {
      sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
    } else {
      sequence[primaryBlockId].first.push_back({'x', {}});
    }
  }

  sequence.resize(maxBlockId + 1);
  blockExists.resize(maxBlockId + 1);
  blockStrand.resize(maxBlockId + 1);

  // Assigning nucleotide gaps
  for (size_t i = 0; i < gaps.size(); i++) {
    int32_t primaryBId = (gaps[i].primaryBlockId);
    int32_t secondaryBId = (gaps[i].secondaryBlockId);

    for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
      int len = gaps[i].nucGapLength[j];
      int pos = gaps[i].nucPosition[j];

      if (secondaryBId != -1) {
        sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
      } else {
        sequence[primaryBId].first[pos].second.resize(len, '-');
      }
    }
  }

  // Get all blocks on the path
  for (auto node = path.rbegin(); node != path.rend(); node++) {
    for (auto mutation : (*node)->blockMutation) {
      int primaryBlockId = mutation.primaryBlockId;
      int secondaryBlockId = mutation.secondaryBlockId;
      int type = (mutation.blockMutInfo);
      bool inversion = mutation.inversion;

      if (type == panmanUtils::BlockMutationType::BI) {
        if (secondaryBlockId != -1) {
          blockExists[primaryBlockId].second[secondaryBlockId] = true;

          // if insertion of inverted block takes place, the strand is backwards
          blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
        } else {
          blockExists[primaryBlockId].first = true;

          // if insertion of inverted block takes place, the strand is backwards
          blockStrand[primaryBlockId].first = !inversion;
        }
      } else {
        if (inversion) {
          // This is not actually a deletion but an inversion
          if (secondaryBlockId != -1) {
            blockStrand[primaryBlockId].second[secondaryBlockId] =
                !blockStrand[primaryBlockId].second[secondaryBlockId];
          } else {
            blockStrand[primaryBlockId].first =
                !blockStrand[primaryBlockId].first;
          }
        } else {
          // Actually a deletion
          if (secondaryBlockId != -1) {
            blockExists[primaryBlockId].second[secondaryBlockId] = false;
            blockStrand[primaryBlockId].second[secondaryBlockId] = true;
          } else {
            blockExists[primaryBlockId].first = false;
            blockStrand[primaryBlockId].first = true;
          }
        }
      }
    }
  }

  // Apply nucleotide mutations
  for (auto node = path.rbegin(); node != path.rend(); node++) {

    for (size_t i = 0; i < (*node)->nucMutation.size(); i++) {

      int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
      int32_t secondaryBlockId = (*node)->nucMutation[i].secondaryBlockId;

      if (secondaryBlockId != -1) {
        if (!blockExists[primaryBlockId].second[secondaryBlockId]) {
          continue;
        }
      } else {
        if (!blockExists[primaryBlockId].first) {
          continue;
        }
      }

      int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
      int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
      uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
      char newVal = '-';

      if (type < 3) {

        int len = (((*node)->nucMutation[i].mutInfo) >> 4);

        if (type == panmanUtils::NucMutationType::NS) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1)
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            }
            else
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = newVal;
              }
            }
          }
          else
          {
            if (nucGapPosition != -1)
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            }
            else
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              }
            }
          }
        }
        else if (type == panmanUtils::NucMutationType::NI) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1)
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            }
            else
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = newVal;
              }
            }
          }
          else
          {
            if (nucGapPosition != -1)
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            }
            else
            {
              for (int j = 0; j < len; j++)
              {
                newVal = panmanUtils::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              }
            }
          }
        }
        else if (type == panmanUtils::NucMutationType::ND) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1)
            {
              for (int j = 0; j < len; j++)
              {
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = '-';
              }
            }
            else
            {
              for (int j = 0; j < len; j++)
              {
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = '-';
              }
            }
          }
          else
          {
            if (nucGapPosition != -1)
            {
              for (int j = 0; j < len; j++)
              {
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = '-';
              }
            }
            else
            {
              for (int j = 0; j < len; j++)
              {
                sequence[primaryBlockId].first[nucPosition + j].first = '-';
              }
            }
          }
        }
      }
      else
      {
        if (type == panmanUtils::NucMutationType::NSNPS) {
          newVal = panmanUtils::getNucleotideFromCode(
              (((*node)->nucMutation[i].nucs) >> 20) & 0xF);
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1)
            {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = newVal;
            }
            else
            {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = newVal;
            }
          }
          else
          {
            if (nucGapPosition != -1)
            {
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = newVal;
            }
            else
            {
              sequence[primaryBlockId].first[nucPosition].first = newVal;
            }
          }
        }
        else if (type == panmanUtils::NucMutationType::NSNPI) {
          newVal = panmanUtils::getNucleotideFromCode(
              (((*node)->nucMutation[i].nucs) >> 20) & 0xF);
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1)
            {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = newVal;
            }
            else
            {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = newVal;
            }
          }
          else
          {
            if (nucGapPosition != -1)
            {
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = newVal;
            }
            else
            {
              sequence[primaryBlockId].first[nucPosition].first = newVal;
            }
          }
        }
        else if (type == panmanUtils::NucMutationType::NSNPD) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1)
            {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = '-';
            }
            else
            {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = '-';
            }
          }
          else
          {
            if (nucGapPosition != -1)
            {
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = '-';
            }
            else
            {
              sequence[primaryBlockId].first[nucPosition].first = '-';
            }
          }
        }
      }
    }
  }

  if (!aligned &&
      T->rotationIndexes.find(reference) != T->rotationIndexes.end() &&
      T->rotationIndexes[reference] != 0) {
    int ctr = -1, rotInd = 0;
    for (size_t i = 0; i < blockExists.size(); i++) {
      if (blockExists[i].first) {
        ctr++;
      }
      if (ctr == T->rotationIndexes[reference]) {
        rotInd = i;
        break;
      }
    }
    rotate(sequence.begin(), sequence.begin() + rotInd, sequence.end());
    rotate(blockExists.begin(), blockExists.begin() + rotInd,
           blockExists.end());
    rotate(blockStrand.begin(), blockStrand.begin() + rotInd,
           blockStrand.end());
  }

  if (sequenceInverted.find(reference) != sequenceInverted.end() &&
      sequenceInverted[reference]) {
    reverse(sequence.begin(), sequence.end());
    reverse(blockExists.begin(), blockExists.end());
    reverse(blockStrand.begin(), blockStrand.end());
  }

  std::string sequenceString;
  for (size_t i = 0; i < sequence.size(); i++) {
    // Iterate through gap blocks - CURRENTLY NOT BEING USED
    for (size_t j = 0; j < sequence[i].second.size(); j++) {
      if (blockExists[i].second[j]) {
        if (blockStrand[i].second[j]) {
          // If forward strand
          for (size_t k = 0; k < sequence[i].second[j].size(); k++) {
            for (size_t w = 0; w < sequence[i].second[j][k].second.size();
                 w++) {
              if (sequence[i].second[j][k].second[w] == 'x' ||
                  sequence[i].second[j][k].second[w] == '-') {
                if (aligned) {
                  sequenceString += '-';
                }
              } else {
                sequenceString += sequence[i].second[j][k].second[w];
              }
            }
            if (sequence[i].second[j][k].first == 'x' ||
                sequence[i].second[j][k].first == '-') {
              if (aligned) {
                sequenceString += '-';
              }
            } else {
              sequenceString += sequence[i].second[j][k].first;
            }
          }
        } else {
          for (size_t k = sequence[i].second[j].size() - 1; k + 1 > 0; k--) {
            // If reverse strand
            if (sequence[i].second[j][k].first == 'x' ||
                sequence[i].second[j][k].first == '-') {
              if (aligned) {
                sequenceString += '-';
              }
            } else {
              sequenceString += sequence[i].second[j][k].first;
            }
            for (size_t w = sequence[i].second[j][k].second.size() - 1;
                 w + 1 > 0; w--) {
              if (sequence[i].second[j][k].second[w] == 'x' ||
                  sequence[i].second[j][k].second[w] == '-') {
                if (aligned) {
                  sequenceString += '-';
                }
              } else {
                sequenceString += sequence[i].second[j][k].second[w];
              }
            }
          }
        }
      } else {
        if (aligned) {
          for (size_t k = 0; k < sequence[i].second[j].size(); k++) {
            for (size_t w = 0; w < sequence[i].second[j][k].second.size();
                 w++) {
              sequenceString += '-';
            }
            sequenceString += '-';
          }
        }
      }
    }

    // Main block
    if (blockExists[i].first) {
      if (blockStrand[i].first) {
        for (size_t j = 0; j < sequence[i].first.size(); j++) {
          for (size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
            if (sequence[i].first[j].second[k] == 'x' ||
                sequence[i].first[j].second[k] == '-') {
              // This shouldn't be possible but I'm still keeping it since it
              // doesn't hurt
              if (aligned) {
                sequenceString += '-';
              }
            } else {
              sequenceString += sequence[i].first[j].second[k];
            }
          }
          if (sequence[i].first[j].first == 'x' ||
              sequence[i].first[j].first == '-') {
            if (aligned) {
              sequenceString += '-';
            }
          } else {
            sequenceString += sequence[i].first[j].first;
          }
        }
      } else {
        // If reverse strand
        for (size_t j = sequence[i].first.size() - 1; j + 1 > 0; j--) {
          if (sequence[i].first[j].first == 'x' ||
              sequence[i].first[j].first == '-') {
            if (aligned) {
              sequenceString += '-';
            }
          } else {
            sequenceString +=
                getComplementCharacter(sequence[i].first[j].first);
          }
          for (size_t k = sequence[i].first[j].second.size() - 1; k + 1 > 0;
               k--) {
            if (sequence[i].first[j].second[k] == 'x' ||
                sequence[i].first[j].second[k] == '-') {
              // This shouldn't be possible but I'm still keeping it since it
              // doesn't hurt
              if (aligned) {
                sequenceString += '-';
              }
            } else {
              sequenceString +=
                  getComplementCharacter(sequence[i].first[j].second[k]);
            }
          }
        }
      }
    } else {
      if (aligned) {
        for (size_t j = 0; j < sequence[i].first.size(); j++) {
          for (size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
            sequenceString += '-';
          }
          sequenceString += '-';
        }
      }
    }
  }

  int offset = 0;
  if (!aligned &&
      circularSequences.find(reference) != circularSequences.end()) {
    offset = circularSequences[reference];
  }
  if (offset == 0) {
    return sequenceString;
  } else {
    return sequenceString.substr(offset) + sequenceString.substr(0, offset);
  }
}
#include "tree.hpp"
#include "pmi.hpp"
#include "seeding.hpp"
#include <cmath>
#include <iostream>
//#include <sys/_types/_int64_t.h>

using namespace PangenomeMAT;
using namespace tree;

std::chrono::time_point<std::chrono::high_resolution_clock> global_timer =
    std::chrono::high_resolution_clock::now();
void time_stamp() {
  std::chrono::time_point<std::chrono::high_resolution_clock> newtime =
      std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> duration = newtime - global_timer;
  std::cerr << "timing " << duration.count() << "\n\n";
  global_timer = newtime;
}

std::string tree::getConsensus(Tree *T) {
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

// coords (blockId, nucPosition, nucGapPosition)
// Sequence includes both boundary coordinates
//TODO be more efficient for non existant blocks
std::string tree::getNucleotideSequenceFromBlockCoordinates(
    const tupleCoord_t &start, tupleCoord_t &end, const sequence_t &sequence,
    const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const Tree *T, const Node *node, const globalCoords_t &globalCoords, CoordNavigator &navigator) {
  // TODO handle these
  // const auto &rotationIndexes = T->rotationIndexes;
  // const auto &sequenceInverted = T->sequenceInverted;
  // const auto &circularSequences = T->circularSequences;
    std::string seq = "";

    if (end == tupleCoord_t{-1,-1,-1})
    {
      end= tupleCoord_t{sequence.size() - 1, sequence.back().first.size() - 1, -1};
    }
  
    for(auto currCoord = start; currCoord <= end; currCoord = navigator.increment(currCoord) ){

      if(blockExists[currCoord.blockId].first){
        if (currCoord.nucGapPos == -1){
          seq += sequence[currCoord.blockId].first[currCoord.nucPos].first;
        }else{
          seq += sequence[currCoord.blockId].first[currCoord.nucPos].second[currCoord.nucGapPos];
        }
      }else{
        seq += '-';
      }
    }

    return seq;
}


/*
std::string tree::badgerNucleotideSequenceFromBlockCoordinates(
    const tupleCoord_t &start, tupleCoord_t &end, const sequence_t &sequence,
    const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const Tree *T, const Node *node, const globalCoords_t &globalCoords) {

  // TODO handle these
  // const auto &rotationIndexes = T->rotationIndexes;
  // const auto &sequenceInverted = T->sequenceInverted;
  // const auto &circularSequences = T->circularSequences;

  //std::cout << "THISNODE " << node->identifier << "\n";

  const auto &startBlockId = start.blockId;
  const auto &startNuc = start.nucPos;
  const auto &endBlockId = end.blockId;
  auto &endNuc = end.nucPos;

  std::string sequenceString;
  // todo implement inversions
  if (end == tupleCoord_t{-1, -1, -1}) {
    end.blockId = sequence.size() - 1;
    end.nucPos = sequence[end.blockId].first.size() - 1;
    end.nucGapPos = -1;
  }

  for (int32_t i = startBlockId; i <= endBlockId; i++) {


    // std::cout << "\n\ni: " << i;
    if (!blockExists[i].first) { // all gaps
                                 // std::cout << "NOT exist => ";
      if (i == startBlockId) {
        //  std::cout << "is start => ";
        if (startBlockId == endBlockId) {
          /// std::cout << "is end => ";
          for (int32_t j = std::max(0, startNuc);
               j <= std::min(endNuc, (int32_t)sequence[i].first.size() - 1);
               j++) {
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              sequenceString += '-';
            }
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
            //  std::cout << "(" << i << ", " << j << ", -1): -" << std::endl;
            sequenceString += '-';
          }
        } else {
          //  std::cout << "is not end => ";
          for (int32_t j = std::max(0, startNuc);
               j <= (int32_t)sequence[i].first.size() - 1; j++) {
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              //  std::cout << "(" << i << ", " << j << ", " << k << "): -" <<
              //  std::endl;
              sequenceString += '-';
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              
            }
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
            //  std::cout << "(" << i << ", " << j << ", -1): -" << std::endl;
            sequenceString += '-';
          }
        }
      } else {
        // std::cout << "is not start => ";
        if (i == endBlockId) {
          //  std::cout << "is end => ";
          for (int32_t j = 0;
               j <= std::min(endNuc, (int32_t)sequence[i].first.size() - 1);
               j++) {
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              // std::cout << "(" << i << ", " << j << ", " << k << "): -" <<
              // std::endl;
              sequenceString += '-';
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              
            }
            // std::cout << "(" << i << ", " << j << ", -1): -" << std::endl;
            sequenceString += '-';
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
          }
        } else {
          //  std::cout << "is not end => ";
          for (int32_t j = 0; j < sequence[i].first.size(); j++) {
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              //  std::cout << "(" << i << ", " << j << ", " << k << "): -" <<
              //  std::endl;
              sequenceString += '-';
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              
            }
            //  std::cout << "(" << i << ", " << j << ", -1): -" << std::endl;
            sequenceString += '-';
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
          }
        }
      }
    } else {
      //  std::cout << "exists => ";
      // block exists
      if (i == startBlockId) {
        //  std::cout << "is start => ";
        if (startBlockId == endBlockId) {
          //  std::cout << "is end => ";
          for (int32_t j = std::max(0, startNuc);
               j <= std::min(endNuc, (int32_t)sequence[i].first.size() - 1);
               j++) {
            char nuc;
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              nuc = sequence[i].first[j].second[k];
              // std::cout << "(" << i << ", " << j << ", " << k << "): " << nuc
              // << std::endl;
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              

              sequenceString += (nuc == 'x') ? '-' : nuc;
              // std::cout << "AAA ? " << nuc << " what? " << (nuc == 'x') <<
              // "\n"; std::cout << "now sequenceString: " << sequenceString <<
              // "\n";
            }
            nuc = sequence[i].first[j].first;
            // std::cout << "(" << i << ", " << j << ", -1): " << nuc <<
            // std::endl;
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
            sequenceString += (nuc == 'x') ? '-' : nuc;
            // std::cout << "huh ? " << nuc << " what? " << (nuc == 'x') <<
            // "\n"; std::cout << "now sequenceString: " << sequenceString <<
            // "\n";
          }
        } else {
          // std::cout << "is not end => ";
          for (int32_t j = std::max(0, startNuc); j < sequence[i].first.size();
               j++) {
            char nuc;
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              nuc = sequence[i].first[j].second[k];
              // std::cout << "(" << i << ", " << j << ", " << k << "): " << nuc
              // << std::endl;
              sequenceString += nuc == 'x' ? '-' : nuc;
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              
            }
            nuc = sequence[i].first[j].first;
            sequenceString += nuc == 'x' ? '-' : nuc;
            // std::cout << "(" << i << ", " << j << ", -1): " << nuc <<
            // std::endl;
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
          }
        }
      } else {
        //  std::cout << "is not start => ";
        if (i == endBlockId) {
          //  std::cout << "is end => ";
          for (int32_t j = 0;
               j <= std::min(endNuc, (int32_t)sequence[i].first.size() - 1);
               j++) {
            char nuc;
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              nuc = sequence[i].first[j].second[k];
              // std::cout << "(" << i << ", " << j << ", " << k << "): " << nuc
              // << std::endl;
              sequenceString += (nuc == 'x') ? '-' : nuc;
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              
              // std::cout << "NUC: " << nuc << " is x? " << (nuc == 'x') <<
              // "\n"; std::cout << "now sequenceString: " << sequenceString <<
              // "\n";
            }
            nuc = sequence[i].first[j].first;
            //  std::cout << "(" << i << ", " << j << ", -1): " << nuc <<
            //  std::endl;
            sequenceString += (nuc == 'x') ? '-' : nuc;
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
            // std::cout << "TEE: " << nuc << " is x? " << (nuc == 'x') << "\n";
            // std::cout << "now sequenceString: " << sequenceString << "\n";
          }
        } else {
          //  std::cout << "is not end => ";
          for (int32_t j = 0; j < sequence[i].first.size(); j++) {
            char nuc;
            for (int32_t k = 0; k < sequence[i].first[j].second.size(); k++) {
              nuc = sequence[i].first[j].second[k];
              //  std::cout << "(" << i << ", " << j << ", " << k << "): " <<
              //  nuc << std::endl;
              sequenceString += nuc == 'x' ? '-' : nuc;
              int64_t gc = getGlobalCoordinate(i, j, k, globalCoords);
              
            }
            nuc = sequence[i].first[j].first;
            //  std::cout << "(" << i << ", " << j << ", -1): " << nuc <<
            //  std::endl;
            sequenceString += nuc == 'x' ? '-' : nuc;
            int64_t gc = getGlobalCoordinate(i, j, -1, globalCoords);
            
          }
        }
      }
    }
  }
  return sequenceString;
}
*/








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
  // std::string seq = tree::getStringFromCurrData(data, T, node, true);

  // strings[node->identifier] = seq;

  // /* Recursive step */
  // for (Node *child : node->children) {
  //   getAllStringsHelper(strings, data, T, child, globalCoords);
  // }

  // /* Undo mutations when backtracking */
  // undoMutations(data, blank, T, node, blockMutData, nucMutData);
}

std::unordered_map<std::string, std::string> tree::getAllNodeStrings(Tree *T) {

  std::unordered_map<std::string, std::string> strings;
  tree::mutableTreeData data;
  tree::globalCoords_t globalCoords;
  setup(data, globalCoords, T);

  getAllStringsHelper(strings, data, T, T->root, globalCoords);

  return strings;
}
std::string tree::getStringFromCurrData(mutableTreeData &data, Tree *T,
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
void tree::setupGlobalCoordinates(int64_t &ctr, globalCoords_t &globalCoords,
                                  const BlockGapList &blockGaps,
                                  const std::vector<Block> &blocks,
                                  const std::vector<GapList> &gaps,
                                  const sequence_t &sequence) {

  globalCoords.resize(blocks.size() + 1);
  // Assigning block gaps
  for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
    globalCoords[blockGaps.blockPosition[i]].second.resize(
        blockGaps.blockGapLength[i]);
  }
  int32_t maxBlockId = 0;
  for (size_t i = 0; i < blocks.size(); i++) {
    int32_t blockId = ((int32_t)blocks[i].primaryBlockId);
    maxBlockId = std::max(maxBlockId, blockId);
    for (size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
      bool endFlag = false;
      for (size_t k = 0; k < 8; k++) {
        const int nucCode =
            (((blocks[i].consensusSeq[j]) >> (4 * (7 - k))) & 15);
        if (nucCode == 0) {
          endFlag = true;
          break;
        }
        globalCoords[blockId].first.push_back({0, {}});
      }
      if (endFlag) {
        break;
      }
    }
    globalCoords[blockId].first.push_back({0, {}});
  }
  globalCoords.resize(maxBlockId + 1);
  // Assigning nucleotide gaps
  for (size_t i = 0; i < gaps.size(); i++) {
    int32_t blockId = (gaps[i].primaryBlockId);
    for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
      int len = gaps[i].nucGapLength[j];
      int pos = gaps[i].nucPosition[j];
      globalCoords[blockId].first[pos].second.resize(len, 0);
    }
  }
  // Assigning coordinates
  ctr = 0;
  size_t i;
  size_t j;
  size_t k;
  for (i = 0; i < globalCoords.size(); i++) {
    for (j = 0; j < globalCoords[i].first.size(); j++) {
      for (k = 0; k < globalCoords[i].first[j].second.size(); k++) {
        globalCoords[i].first[j].second[k] = ctr;
        ctr++;
      }
      globalCoords[i].first[j].first = ctr;
      ctr++;
    }
  }
}
void tree::setup(mutableTreeData &data, globalCoords_t &globalCoords,
                 const Tree *T) {

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
    // std::cerr << b << "hello " << i << "\n";
    maxBlock = std::max(maxBlock, b);
    for (size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
      bool stop = false;
      for (size_t k = 0; k < 8; k++) {
        const int nc = (((blocks[i].consensusSeq[j]) >> (4 * (7 - k))) & 15);
        if (nc == 0) {
          stop = true;
          break;
        }
        const char c = PangenomeMAT::getNucleotideFromCode(nc);
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
                         blocks, gaps, sequence);
}
int64_t tree::getGlobalCoordinate(const int blockId, const int nucPosition,
                                  const int nucGapPosition,
                                  const globalCoords_t &globalCoords) {

  if (blockId > globalCoords.size() - 1 ||
      blockId == globalCoords.size() - 1 &&
          nucPosition > globalCoords.back().first.size() - 1) {
    return globalCoords[globalCoords.size() - 1]
        .first[globalCoords.back().first.size() - 1]
        .first;
  }

  if (nucGapPosition == -1) {
    // std::cout << "accessing globalCoords[" << blockId << "].first[" <<
    // nucPosition << "].first\n";
    return globalCoords[blockId].first[nucPosition].first;
  }

  if (globalCoords[blockId].first[nucPosition].second.size() == 0) {
    return globalCoords[blockId].first[nucPosition].first;
  }
  // std::cout << "** accessing globalCoords[" << blockId << "].first[" <<
  // nucPosition << "].second[" << nucGapPosition << "]\n";
  return globalCoords[blockId].first[nucPosition].second[nucGapPosition];
}
static int getIndexFromNucleotide(char nuc) {

  switch (nuc) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
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
std::pair<size_t, size_t> tree::getMaskCoorsForMutmat(const std::string &s1,
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
void buildMutationMatrices(mutationMatrices &mutMat, Tree *T, size_t window,
                           double threshold) {

  std::unordered_map<std::string, std::string> alignedSequences =
      getAllNodeStrings(T);
  for (const auto &sequence : alignedSequences) {
    std::string parentId;
    if (T->allNodes[sequence.first]->parent == nullptr ||
        T->allNodes[sequence.first]->nucMutation.size() == 0) {
      continue;
    } else {
      parentId = T->allNodes[sequence.first]->parent->identifier;
    }

    const std::string &curSeq = sequence.second;
    const std::string &parSeq = alignedSequences[parentId];
    size_t insLen = 0;
    size_t delLen = 0;
    std::pair<size_t, size_t> edgeCoor =
        tree::getMaskCoorsForMutmat(curSeq, parSeq, window, threshold);
    assert(edgeCoor.second >= edgeCoor.first);
    for (size_t i = edgeCoor.first; i < edgeCoor.second + 1; i++) {
      if (parSeq[i] == '-' && curSeq[i] == '-') {
        continue;
      } else if (parSeq[i] != '-' && curSeq[i] == '-') {
        delLen++;
        if (insLen > 0) {
          if (insLen > mutMat.insmat.size() - 1) {
            mutMat.insmat.resize(insLen + 1);
          }
          mutMat.insmat[insLen]++;
          mutMat.total_insmut++;
          insLen = 0;
        }
      } else if (parSeq[i] == '-' && curSeq[i] != '-') {
        insLen++;
        if (delLen > 0) {
          if (delLen > mutMat.delmat.size() - 1) {
            mutMat.delmat.resize(delLen + 1);
          }
          mutMat.delmat[delLen]++;
          mutMat.total_delmut++;
          delLen = 0;
        }
      } else {
        if (insLen > mutMat.insmat.size() - 1) {
          mutMat.insmat.resize(insLen + 1);
        }
        mutMat.insmat[insLen]++;
        mutMat.total_insmut++;
        insLen = 0;

        if (delLen > mutMat.delmat.size() - 1) {
          mutMat.delmat.resize(delLen + 1);
        }
        mutMat.delmat[delLen]++;
        mutMat.total_delmut++;
        delLen = 0;

        int parNucIdx = getIndexFromNucleotide(parSeq[i]);
        int curNucIdx = getIndexFromNucleotide(curSeq[i]);
        if (parNucIdx > 3 || curNucIdx > 3) {
          continue;
        }
        mutMat.submat[parNucIdx][curNucIdx]++;
        mutMat.total_submuts[parNucIdx]++;
      }
    }
  }

  // insertion
  for (auto i = 0; i < mutMat.insmat.size(); ++i) {
    mutMat.insmat[i] = -10 * log10f(mutMat.insmat[i] / mutMat.total_insmut);
  }
  // deletion
  for (auto i = 0; i < mutMat.delmat.size(); ++i) {
    mutMat.delmat[i] = -10 * log10f(mutMat.delmat[i] / mutMat.total_delmut);
  }
  // substitution
  for (auto i = 0; i < 4; i++) {
    for (auto j = 0; j < 4; j++) {
      mutMat.submat[i][j] =
          -10 * log10f(mutMat.submat[i][j] / mutMat.total_submuts[i]);
    }
  }
}
void tree::writeMutationMatrices(const mutationMatrices &mutMat,
                                 std::ofstream &mmfout) {
  for (const std::vector<double> &row : mutMat.submat) {
    for (const double &prob : row) {
      mmfout << prob << " ";
    }
    mmfout << "\n";
  }
  for (const double &prob : mutMat.insmat) {
    mmfout << prob << " ";
  }
  mmfout << "\n";
  for (const double &prob : mutMat.delmat) {
    mmfout << prob << " ";
  }
  mmfout << "\n";
}

void tree::fillMutationMatricesFromTree(mutationMatrices &mutMat, Tree *T,
                                        size_t window, double threshold) {
  buildMutationMatrices(mutMat, T, window, threshold);
  mutMat.filled = true;
}

void tree::fillMutationMatricesFromFile(mutationMatrices &mutMat,
                                        std::ifstream &inf) {
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
            "Received invalid mutamtion matrix (.mm) file");
      }
      mutMat.submat[idx] = std::move(probs);
    } else if (idx == 4) {
      if (probs.size() < 1) {
        throw std::invalid_argument(
            "Received invalid mutamtion matrix (.mm) file");
      }
      mutMat.insmat = std::move(probs);
    } else if (idx == 5) {
      if (probs.size() < 1) {
        throw std::invalid_argument(
            "Received invalid mutamtion matrix (.mm) file");
      }
      mutMat.delmat = std::move(probs);
    }
    idx++;
  }

  if (idx != 6) {
    throw std::invalid_argument("Received invalid mutamtion matrix (.mm) file");
  }
  mutMat.filled = true;
}

std::string tree::getStringAtNode(Node *node, Tree *T, bool aligned) {

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

  std::vector<PangenomeMAT::Node *> path;
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
        const char nucleotide = PangenomeMAT::getNucleotideFromCode(nucCode);

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

      if (type == PangenomeMAT::BlockMutationType::BI) {
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

        if (type == PangenomeMAT::NucMutationType::NS) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            } else {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = newVal;
              }
            }
          } else {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            } else {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              }
            }
          }
        } else if (type == PangenomeMAT::NucMutationType::NI) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            } else {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = newVal;
              }
            }
          } else {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = newVal;
              }
            } else {
              for (int j = 0; j < len; j++) {
                newVal = PangenomeMAT::getNucleotideFromCode(
                    (((*node)->nucMutation[i].nucs) >> (4 * (5 - j))) & 0xF);
                sequence[primaryBlockId].first[nucPosition + j].first = newVal;
              }
            }
          }
        } else if (type == PangenomeMAT::NucMutationType::ND) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition]
                    .second[nucGapPosition + j] = '-';
              }
            } else {
              for (int j = 0; j < len; j++) {
                sequence[primaryBlockId]
                    .second[secondaryBlockId][nucPosition + j]
                    .first = '-';
              }
            }
          } else {
            if (nucGapPosition != -1) {
              for (int j = 0; j < len; j++) {
                sequence[primaryBlockId]
                    .first[nucPosition]
                    .second[nucGapPosition + j] = '-';
              }
            } else {
              for (int j = 0; j < len; j++) {
                sequence[primaryBlockId].first[nucPosition + j].first = '-';
              }
            }
          }
        }
      } else {
        if (type == PangenomeMAT::NucMutationType::NSNPS) {
          newVal = PangenomeMAT::getNucleotideFromCode(
              (((*node)->nucMutation[i].nucs) >> 20) & 0xF);
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = newVal;
            } else {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = newVal;
            }
          } else {
            if (nucGapPosition != -1) {
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = newVal;
            } else {
              sequence[primaryBlockId].first[nucPosition].first = newVal;
            }
          }
        } else if (type == PangenomeMAT::NucMutationType::NSNPI) {
          newVal = PangenomeMAT::getNucleotideFromCode(
              (((*node)->nucMutation[i].nucs) >> 20) & 0xF);
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = newVal;
            } else {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = newVal;
            }
          } else {
            if (nucGapPosition != -1) {
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = newVal;
            } else {
              sequence[primaryBlockId].first[nucPosition].first = newVal;
            }
          }
        } else if (type == PangenomeMAT::NucMutationType::NSNPD) {
          if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .second[nucGapPosition] = '-';
            } else {
              sequence[primaryBlockId]
                  .second[secondaryBlockId][nucPosition]
                  .first = '-';
            }
          } else {
            if (nucGapPosition != -1) {
              sequence[primaryBlockId]
                  .first[nucPosition]
                  .second[nucGapPosition] = '-';
            } else {
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

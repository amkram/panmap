#ifndef __TREE_HPP
#define __TREE_HPP

#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>

void time_stamp();

using namespace seeding;

inline auto seed_cmp = [](const std::pair<int32_t, std::string> &a,
                          const std::pair<int32_t, std::string> &b) {
  if (a.second != b.second) {
    return a.second < b.second;
  }
  return a.first < b.first;
};

struct tupleCoord_t {
    int blockId, nucPos, nucGapPos;

    // Constructor
    tupleCoord_t(const int32_t &blockId, const int32_t &nucPos,
                 const int32_t &nucGapPos)
        : blockId(blockId), nucPos(nucPos), nucGapPos(nucGapPos) {}

    bool operator<(const tupleCoord_t &rhs) const {
        if (blockId == -1 && nucPos == -1 && nucGapPos == -1) return false;
        if (rhs.blockId == -1 && rhs.nucPos == -1 && rhs.nucGapPos == -1) return true;
        if (blockId != rhs.blockId) return blockId < rhs.blockId;
        if (nucPos != rhs.nucPos) return nucPos < rhs.nucPos;
        //if (nucGapPos != -1 && rhs.nucGapPos != -1) return nucGapPos < rhs.nucGapPos;
        if (nucGapPos == -1 && rhs.nucGapPos != -1) return false;
        if (nucGapPos != -1 && rhs.nucGapPos == -1) return true;
        return nucGapPos < rhs.nucGapPos;
    }

    bool operator<=(const tupleCoord_t &rhs) const {
        return *this < rhs || *this == rhs;
    }

    bool operator==(const tupleCoord_t &rhs) const {
        return blockId == rhs.blockId && nucPos == rhs.nucPos && nucGapPos == rhs.nucGapPos;
    }

    bool operator>=(const tupleCoord_t &rhs) const {
        return !(*this < rhs);
    }

    bool operator>(const tupleCoord_t &rhs) const {
        return !(*this < rhs || *this == rhs);
    }

};


struct TupleHash { //TODO
  std::size_t operator()(const tupleCoord_t& s) const noexcept {

    std::size_t seed = 0;

    seed ^= s.blockId + 0x9e3779b9 + (seed<<6) + (seed>>2);
    seed ^= s.nucPos + 0x9e3779b9 + (seed<<6) + (seed>>2);
    seed ^= s.nucGapPos + 0x9e3779b9 + (seed<<6) + (seed>>2);

    return seed; // or use boost::hash_combine
  }
};


struct tupleRange {
    tupleCoord_t start;
    tupleCoord_t stop;

    bool operator<(const tupleRange &rhs) const {
        return start < rhs.start;
    }
};
class CoordNavigator {
public:
  CoordNavigator(sequence_t &sequence) : sequence(sequence) {}

  bool isGap(const tupleCoord_t &coord) {
    char c;
    if (coord.nucGapPos == -1) {
      c = sequence[coord.blockId].first[coord.nucPos].first;
    } else {
      c = sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
    }
    return c == '-' || c == 'x';
  }

  
  tupleCoord_t increment(tupleCoord_t &givencoord) {
    tupleCoord_t coord = givencoord;
    //std::cout << "honda acoord.nucPos: " << coord.nucPos << std::endl;
    if (coord.nucGapPos == -1){
      coord.nucPos++;
      //std::cout << "coord.nucPos: " << coord.nucPos << std::endl;

      if(coord.nucPos >= sequence[coord.blockId].first.size()){
        //std::cout << "yeah... that just happened....\n";
        coord.blockId++;
        //std::cout << "coord.blockId: " << coord.blockId << std::endl;
        if(coord.blockId >= sequence.size()){
          //std::cout << "yeah... that also just happened....\n";
          return tupleCoord_t{-1,-1,-1};
        }
        coord.nucPos = 0;
      }

      if(!sequence[coord.blockId].first[coord.nucPos].second.empty()) {
        coord.nucGapPos = 0;
      }
      return coord;
    }else{
      coord.nucGapPos++;
      if(coord.nucGapPos >= sequence[coord.blockId].first[coord.nucPos].second.size()){
        coord.nucGapPos = -1;
      }
      return coord;
    }
  }

  tupleCoord_t decrement(tupleCoord_t &givencoord) {
    tupleCoord_t coord = givencoord;
    if (coord.nucGapPos == -1) {
      if(sequence[coord.blockId].first[coord.nucPos].second.empty()){
        coord.nucPos--;
        if(coord.nucPos < 0){
          coord.blockId--;
          if(coord.blockId < 0){
            return tupleCoord_t{0,0,0};
          }
          coord.nucPos = sequence[coord.blockId].first.size() - 1;
        }
        return coord;
      }else{
        coord.nucGapPos = sequence[coord.blockId].first[coord.nucPos].second.size() - 1;
        return coord;
      }
    }else{
      coord.nucGapPos--;
      if(coord.nucGapPos < 0){
        coord.nucGapPos = -1;
        coord.nucPos--;
        if(coord.nucPos < 0){
          coord.blockId--;
          if(coord.blockId < 0){
            return tupleCoord_t{0,0,0};
          }
          coord.nucPos = sequence[coord.blockId].first.size() - 1;
        }
        return coord;
      }
      return coord;
    }
  }


public:
  sequence_t &sequence;
};


typedef std::unordered_map<
    std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>>
    seedmerIndex_t;

/* Helpers for interacting with panmats */
namespace tree {
using namespace PangenomeMAT;

typedef std::vector<
    std::pair<std::vector<std::pair<int, std::vector<int>>>,
              std::vector<std::vector<std::pair<int, std::vector<int>>>>>>
    globalCoords_t;
typedef std::vector<std::tuple<int32_t, int32_t, bool, bool, bool, bool>> blockMutationInfo_t;
typedef std::vector<
    std::tuple<int32_t, int32_t, int32_t, int32_t, char, char>>
    mutationInfo_t;

struct mutableTreeData {
  // These fields are intended to be mutated at each node during a DFS
  sequence_t sequence; // the main object encoding the MSA
  int64_t maxGlobalCoordinate;
  std::string ungappedConsensus; // not used

  std::vector<seed> seeds; // dynamic vector of seeds in each node's sequence
  std::vector<seedmer> seedmers;
  std::unordered_map<std::string, bool> variableSeeds;         // seeds in the consensus that mutate at least once
  blockExists_t blockExists; // tracks if blocks are "on" at a node
  blockStrand_t blockStrand; // tracks strand of blocks
};

struct mutationMatrices {
  // Store mutation matrices
  std::vector<std::vector<double>> submat; // 4 x 4 substitution rate matrix
  std::vector<double> insmat = {0}; // 1 x N insertion rate by length matrix
  std::vector<double> delmat = {0}; // 1 x N deletion rate by length matrix

  // Stores total number of mutations
  bool filled = false;
  std::vector<double> total_submuts;
  double total_insmut = 0;
  double total_delmut = 0;

  mutationMatrices() {
    // initialize mutationMatrices object and intialize the correct size for
    // substitution amtrix
    total_submuts.resize(4);
    submat.resize(4);
    for (size_t i = 0; i < 4; ++i) {
      submat[i].resize(4);
    }
  }
};
/* Interface */
void removeIndices(std::vector<seed> &v, std::stack<int32_t> &rm);
std::string getConsensus(Tree *T); // ungapped!

std::unordered_map<std::string, std::string> getAllNodeStrings(Tree *T);
std::string getStringFromCurrData(mutableTreeData &data, Tree *T,
                                  const Node *node, const bool aligned);

int64_t getGlobalCoordinate(const int blockId, const int nucPosition,
                            const int nucGapPosition,
                            const globalCoords_t &globalCoords);
void setup(mutableTreeData &data, globalCoords_t &globalCoords, const Tree *T);

void setupGlobalCoordinates(
    int64_t &ctr, globalCoords_t &globalCoords,
    const BlockGapList &blockGaps, const std::vector<Block> &blocks,
    const std::vector<GapList> &gaps,
    const sequence_t &sequence);

// Fill mutation matrices from tree or file
std::pair<size_t, size_t> getMaskCoorsForMutmat(const std::string &s1,
                                                const std::string &s2,
                                                size_t window,
                                                double threshold);
void fillMutationMatricesFromTree(mutationMatrices &mutMat, Tree *T,
                                  size_t window, double threshold);
void fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf);

// Build mutation matrices by traversing through all parent-child pairs
void writeMutationMatrices(const mutationMatrices &mutMat,
                           std::ofstream &mmfout);

std::string getNucleotideSequenceFromBlockCoordinates(
    const tupleCoord_t &start, tupleCoord_t &end,
    const sequence_t &sequence,
    const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const Tree *T, const Node *node, const globalCoords_t &globalCoords, CoordNavigator &navigator);

std::string getStringAtNode(Node *node, Tree *T, bool aligned);

} // namespace tree
#endif
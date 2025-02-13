#ifndef __TREE_HPP
#define __TREE_HPP

#include "panmanUtils.hpp"
#include "seeding.hpp"
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "timing.hpp"
#include <capnp/common.h>
#include <capnp/message.h>
#include <capnp/serialize.h>
#include <capnp/serialize-packed.h>
#include "index.capnp.h"

double time_stamp();

using namespace seeding;

inline auto seed_cmp = [](const std::pair<int32_t, std::string> &a,
                          const std::pair<int32_t, std::string> &b) {
  if (a.second != b.second) {
    return a.second < b.second;
  }
  return a.first < b.first;
};

inline void fillDfsIndexes(panmanUtils::Tree *T, panmanUtils::Node *node, int64_t &dfsIndex, std::unordered_map<std::string, int64_t> &dfsIndexes) {  
  TIME_FUNCTION;
  dfsIndexes[node->identifier] = dfsIndex;
  dfsIndex++;
  for (panmanUtils::Node *child : node->children) {
    fillDfsIndexes(T, child, dfsIndex, dfsIndexes);
  }
}

struct tupleCoord_t {
    int blockId, nucPos, nucGapPos;

    // Constructor
    tupleCoord_t(const int64_t &blockId, const int64_t &nucPos,
                 const int64_t &nucGapPos)
        : blockId(blockId), nucPos(nucPos), nucGapPos(nucGapPos) {}

    // Default constructor
    tupleCoord_t() : blockId(-1), nucPos(-1), nucGapPos(-1) {}
    tupleCoord_t(const tupleCoord_t &other) : blockId(other.blockId), nucPos(other.nucPos), nucGapPos(other.nucGapPos) {}

    bool operator<(const tupleCoord_t &rhs) const {
        if (blockId == -1 && nucPos == -1 && nucGapPos == -1) return false;
        if (rhs.blockId == -1 && rhs.nucPos == -1 && rhs.nucGapPos == -1) return true;
        if (blockId != rhs.blockId) return blockId < rhs.blockId;
        if (nucPos != rhs.nucPos) return (nucPos < rhs.nucPos);
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

static inline bool compareNucMuts(const panmanUtils::NucMut &a, const panmanUtils::NucMut &b) {
  return tupleCoord_t{a.primaryBlockId, a.nucPosition, a.nucGapPosition} < tupleCoord_t{b.primaryBlockId, b.nucPosition, b.nucGapPosition};
}

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
private:
    sequence_t* sequence_ptr;  // Use pointer instead of reference

public:
    // Default constructor now allowed - sets null pointer
    CoordNavigator() : sequence_ptr(nullptr) {}
    
    // Main constructor
    explicit CoordNavigator(sequence_t& seq) : sequence_ptr(&seq) {}
    
    // Copy/move constructors
    CoordNavigator(const CoordNavigator& other) = default;
    CoordNavigator(CoordNavigator&& other) noexcept = default;
    
    // Assignment operators
    CoordNavigator& operator=(const CoordNavigator& other) = default;
    CoordNavigator& operator=(CoordNavigator&& other) noexcept = default;

    // Method to initialize/change the sequence
    void setSequence(sequence_t& seq) {
        sequence_ptr = &seq;
    }

    // Access methods now use pointer
    sequence_t& sequence() {
        if (!sequence_ptr) throw std::runtime_error("Accessing uninitialized navigator");
        return *sequence_ptr;
    }

    // All existing methods stay exactly the same, just replace sequence with sequence()
    bool isGap(const tupleCoord_t &coord) {
        char c;
        if (coord.nucGapPos == -1) {
            c = sequence()[coord.blockId].first[coord.nucPos].first;
        } else {
            c = sequence()[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
        }
        return c == '-' || c == 'x';
    }

    tupleCoord_t increment(tupleCoord_t &givencoord) {
        tupleCoord_t &coord = givencoord;
        
        if (coord.nucGapPos == -1){
            coord.nucPos++;

            if(coord.nucPos >= sequence()[coord.blockId].first.size()){
                coord.blockId++;
                
                if(coord.blockId >= sequence().size()){
                    return tupleCoord_t{-1,-1,-1};
                }
                coord.nucPos = 0;
            }

            if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                coord.nucGapPos = 0;
            }
            return coord;
        } else {
            coord.nucGapPos++;
            if(coord.nucGapPos >= sequence()[coord.blockId].first[coord.nucPos].second.size()){
                coord.nucGapPos = -1;
            }
            return coord;
        }
    }

    tupleCoord_t decrement(tupleCoord_t &givencoord) {
        tupleCoord_t &coord = givencoord;
        if (coord.nucGapPos == -1) {
            if(sequence()[coord.blockId].first[coord.nucPos].second.empty()){
                coord.nucPos--;
                if(coord.nucPos < 0){
                    coord.blockId--;
                    if(coord.blockId < 0){
                        return tupleCoord_t{0,0,0};
                    }
                    coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                }
                return coord;
            } else {
                coord.nucGapPos = sequence()[coord.blockId].first[coord.nucPos].second.size() - 1;
                return coord;
            }
        } else {
            coord.nucGapPos--;
            if(coord.nucGapPos < 0){
                coord.nucGapPos = -1;
                coord.nucPos--;
                if(coord.nucPos < 0){
                    coord.blockId--;
                    if(coord.blockId < 0){
                        return tupleCoord_t{0,0,0};
                    }
                    coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                }
                return coord;
            }
            return coord;
        }
    }

    tupleCoord_t newincrement(tupleCoord_t &givencoord, const blockStrand_t &blockStrand) {
        TIME_FUNCTION;
        if (!sequence_ptr) {
            std::cerr << "ERROR: sequence_ptr is null" << std::endl;
            throw std::runtime_error("Navigator sequence not initialized");
        }

        tupleCoord_t coord = givencoord;
        
        // std::cerr << "DEBUG: newincrement called with coord: blockId=" << coord.blockId 
        //           << ", nucPos=" << coord.nucPos 
        //           << ", nucGapPos=" << coord.nucGapPos 
        //           << ", sequence size=" << sequence().size() << std::endl;

        // Check blockId bounds
        if (coord.blockId < 0 || coord.blockId >= sequence().size()) {
            std::cerr << "ERROR: blockId " << coord.blockId << " out of bounds for sequence size " << sequence().size() << std::endl;
            return tupleCoord_t{-1,-1,-1};
        }

        if(blockStrand[coord.blockId].first) {
            if (coord.nucGapPos == -1) {
                coord.nucPos++;

                // std::cerr << "DEBUG: After increment - nucPos=" << coord.nucPos 
                //           << ", first.size()=" << sequence()[coord.blockId].first.size() << std::endl;

                // Check nucPos bounds
                if(coord.nucPos >= sequence()[coord.blockId].first.size()) {
                    // std::cerr << "DEBUG: nucPos exceeded block size, moving to next block" << std::endl;
                    coord.blockId++;
                    
                    if(coord.blockId >= sequence().size()) {
                        std::cerr << "DEBUG: Reached end of sequence" << std::endl;
                        return tupleCoord_t{-1,-1,-1};
                    }
                    if(!blockStrand[coord.blockId].first) {
                        coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                        coord.nucGapPos = -1;
                        return coord;
                    }
                    coord.nucPos = 0;
                }

                if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                    coord.nucGapPos = 0;
                }
                return coord;
            }else{
                coord.nucGapPos++;
                if(coord.nucGapPos >= sequence()[coord.blockId].first[coord.nucPos].second.size()){
                    coord.nucGapPos = -1;
                }
                return coord;
            }

        }else{

            if (coord.nucGapPos == -1) {
                if(sequence()[coord.blockId].first[coord.nucPos].second.empty()){
                    coord.nucPos--;
                    if(coord.nucPos < 0){

                        coord.blockId++;
                    
                        if(coord.blockId >= sequence().size()){
                            return tupleCoord_t{-1,-1,-1};
                        }
                        if(!blockStrand[coord.blockId].first){
                            coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                            coord.nucGapPos = -1;
                            return coord;
                        }
                        coord.nucPos = 0;
                        if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                            coord.nucGapPos = 0;
                        }
                        return coord;

                    }
                    return coord;
                }else{
                    coord.nucGapPos = sequence()[coord.blockId].first[coord.nucPos].second.size() - 1;
                    return coord;
                }
            }else{
                coord.nucGapPos--;
                if(coord.nucGapPos < 0){
                    coord.nucGapPos = -1;
                    coord.nucPos--;
                    if(coord.nucPos < 0){

                        coord.blockId++;
                    
                        if(coord.blockId >= sequence().size()){
                            return tupleCoord_t{-1,-1,-1};
                        }
                        if(!blockStrand[coord.blockId].first){
                            coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                            coord.nucGapPos = -1;
                            return coord;
                        }
                        coord.nucPos = 0;
                        if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                            coord.nucGapPos = 0;
                        }
                        return coord;

                }
                return coord;
            }
            return coord;
        }
      }
    }


    tupleCoord_t newdecrement(tupleCoord_t &givencoord, const blockStrand_t &blockStrand) {
        TIME_FUNCTION;
        tupleCoord_t coord = givencoord;

        if(blockStrand[coord.blockId].first){

        if (coord.nucGapPos == -1) {
            if(sequence()[coord.blockId].first[coord.nucPos].second.empty()){
                coord.nucPos--;
                if(coord.nucPos < 0){
                    coord.blockId--;
                    if(coord.blockId < 0){
                        return tupleCoord_t{0,0,0};
                    }

                    if(blockStrand[coord.blockId].first){
                        coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                        coord.nucGapPos = -1;
                        return coord;
                    }
                    coord.nucPos = 0;
                    if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                        coord.nucGapPos = 0;
                    }
                    return coord;

                }
                return coord;
            }else{
                coord.nucGapPos = sequence()[coord.blockId].first[coord.nucPos].second.size() - 1;
                return coord;
            }
        }else{
            coord.nucGapPos--;
            if(coord.nucGapPos < 0){
                coord.nucGapPos = -1;
                coord.nucPos--;
                if(coord.nucPos < 0){

                    //Jump to previous block
                    coord.blockId--;
                    if(coord.blockId < 0){
                        return tupleCoord_t{0,0,0};
                    }
                    if(blockStrand[coord.blockId].first){
                        coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                        coord.nucGapPos = -1;
                        return coord;
                    }
                    coord.nucPos = 0;
                    if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                        coord.nucGapPos = 0;
                    }
                    return coord;

                }
                return coord;
            }
            return coord;
        }

        }else{

            
        

        tupleCoord_t coord = givencoord;
        
        if (coord.nucGapPos == -1){
            coord.nucPos++;

            if(coord.nucPos >= sequence()[coord.blockId].first.size()){
                
                //Jump to previous block
                coord.blockId--;
                if(coord.blockId < 0){
                    return tupleCoord_t{0,0,0};
                }
                if(blockStrand[coord.blockId].first){
                    coord.nucPos = sequence()[coord.blockId].first.size() - 1;
                    coord.nucGapPos = -1;
                    return coord;
                }
                coord.nucPos = 0;
                if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                    coord.nucGapPos = 0;
                }
                return coord;
            }

            if(!sequence()[coord.blockId].first[coord.nucPos].second.empty()) {
                coord.nucGapPos = 0;
            }
            return coord;
        }else{
            coord.nucGapPos++;
            if(coord.nucGapPos >= sequence()[coord.blockId].first[coord.nucPos].second.size()){
                coord.nucGapPos = -1;
            }
            return coord;
        }

    }
    }
};


typedef std::unordered_map<
    std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>>
    seedmerIndex_t;

/* Helpers for interacting with panmats */
namespace seed_annotated_tree {
using namespace panmanUtils;

typedef std::vector<
    std::pair<std::vector<std::pair<int64_t, std::vector<int64_t>>>,
              std::vector<std::vector<std::pair<int64_t, std::vector<int64_t>>>>>>
    globalCoords_t;
int64_t tupleToScalarCoord(const tupleCoord_t &coord,
                           const globalCoords_t &globalCoords);

typedef std::vector<std::tuple<int32_t, int32_t, bool, bool, bool, bool>> blockMutationInfo_t;
typedef std::vector<
    std::tuple<int32_t, int32_t, int32_t, int32_t, char, char>>
    mutationInfo_t;

class HotSeedIndex {
public:
    typedef size_t hashedKmer_t;
    enum PositionState {
        IMMUTABLE,
        HOT,
        NORMAL
    };
    std::vector<PositionState> positionStates; 
    std::unordered_map<int64_t, std::unordered_map<int32_t, std::tuple<hashedKmer_t, int64_t, bool>>> hotSeeds; // pos -> (dfsIndex -> (kmer, endPos, isReverse))
    std::unordered_map<int64_t, std::tuple<hashedKmer_t, int64_t, bool>> immutableSeeds; // pos -> (kmer, endPos, isReverse)
    std::vector<int64_t> hotCounts; // Track how many times a position has been mutated
    // Constructor
    HotSeedIndex(int k, int hotThreshold, int64_t genomeLength) : k(k), hotThreshold(hotThreshold) {
        positionStates.resize(genomeLength, IMMUTABLE);
        hotCounts.resize(genomeLength, 0);
        immutableSeeds.reserve(genomeLength);
    }

    HotSeedIndex(int k) : k(k) {}
    
    // Record a mutation at a position
    void recordMutation(int64_t pos) {
        positionStates[pos] = NORMAL;
        hotCounts[pos]++;
        if (hotCounts[pos] >= hotThreshold) {
            positionStates[pos] = HOT;
        }
    }

    void indexSeedIfHotOrImmutable(int64_t pos, int32_t dfsIndex, hashedKmer_t kmer, int64_t endPos, bool isReverse) {
        if (positionStates[pos] == IMMUTABLE) {
            immutableSeeds[pos] = {kmer, endPos, isReverse};
        } else if (positionStates[pos] == HOT) {
            hotSeeds[pos][dfsIndex] = {kmer, endPos, isReverse};
        }
    }

    // Get seed for a position and DFS index
    bool getSeed(int64_t pos, int32_t dfsIndex, hashedKmer_t& seed, int64_t& endPos, bool& isReverse) const {
        if (pos >= hotCounts.size()) return false;
        if (immutableSeeds.find(pos) != immutableSeeds.end()) {
            // Position is immutable, return stored seed
            std::tie(seed, endPos, isReverse) = immutableSeeds.at(pos);
            return true;
        }
        
        // Position is hot, look up in hot seeds map
        auto posIt = hotSeeds.find(pos);
        if (posIt == hotSeeds.end()) return false;
        
        auto dfsIt = posIt->second.find(dfsIndex);
        if (dfsIt == posIt->second.end()) return false;
        
        std::tie(seed, endPos, isReverse) = dfsIt->second;
        return true;
    }
    
    // Serialize to Cap'n Proto
    void serialize(::HotSeedIndexSerial::Builder builder) const {
        // Pack immutable seeds
        auto immutableBuilder = builder.initImmutableSeeds(immutableSeeds.size());
        size_t i = 0;
        for (const auto& [pos, seedData] : immutableSeeds) {
            auto entry = immutableBuilder[i];
            auto [seed, endPos, isReverse] = seedData;
            entry.setSeed(seed);
            entry.setEndPosition(endPos);
            entry.setIsReverse(isReverse);
            i++;
        }
        
        // Pack hot seeds
        std::vector<std::tuple<int32_t, int64_t, hashedKmer_t, int64_t, bool>> hotEntries;
        for (const auto& [pos, dfsMap] : hotSeeds) {
            for (const auto& [dfsIndex, seedData] : dfsMap) {
                auto [seed, endPos, isReverse] = seedData;
                hotEntries.push_back({dfsIndex, pos, seed, endPos, isReverse});
            }
        }
        
        auto hotBuilder = builder.initHotSeeds(hotEntries.size());
        for (size_t i = 0; i < hotEntries.size(); i++) {
            auto entry = hotBuilder[i];
            auto [dfsIndex, pos, seed, endPos, isReverse] = hotEntries[i];
            entry.setDfsIndex(dfsIndex);
            entry.setPosition(pos);
            entry.setKmer(seed);
            entry.setEndPosition(endPos);
            entry.setIsReverse(isReverse);
        }
    }
    
    // Deserialize from Cap'n Proto
    void deserialize(::HotSeedIndexSerial::Reader reader) {
        // Read immutable seeds
        auto immutableReader = reader.getImmutableSeeds();
        immutableSeeds.reserve(immutableReader.size());
        for (size_t i = 0; i < immutableReader.size(); i++) {
            auto entry = immutableReader[i];
            immutableSeeds[i] = {entry.getSeed(), entry.getEndPosition(), entry.getIsReverse()};
        }
        
        // Read hot seeds
        auto hotReader = reader.getHotSeeds();
        hotSeeds.reserve(hotReader.size());
        for (auto entry : hotReader) {
            hotSeeds[entry.getPosition()][entry.getDfsIndex()] = {
                entry.getKmer(),
                entry.getEndPosition(),
                entry.getIsReverse()
            };
        }
    }

private:
    int k;  // k-mer length
    int hotThreshold; // times a position must be mutated to be considered hot
};

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
  std::unordered_map<int64_t, tupleCoord_t> scalarToTupleCoord;

  mutableTreeData() {
    maxGlobalCoordinate = 0;
  }
};

struct mutationMatrices {
  // Store mutation matrices
  std::vector<std::vector<double>> submat; // 4 x 4 substitution rate matrix
  std::unordered_map<int64_t, double> insmat; // 1 x N insertion rate by length matrix
  std::unordered_map<int64_t, double> delmat; // 1 x N deletion rate by length matrix

  bool filled = false;

  double maxInsLogProb = 40;
  double maxDelLogProb = 40;

  mutationMatrices() {
    // initialize mutationMatrices object and intialize the correct size for
    // substitution amtrix
    // total_submuts.resize(4);
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
    const sequence_t &sequence,
    std::unordered_map<int64_t, tupleCoord_t> &scalarToTupleCoord);



bool getSeedAt(bool useHotSeedIndex, HotSeedIndex& hotSeedIndex, std::string &seedBuffer, size_t &result_hash, int64_t &result_end_pos, bool &result_is_reverse, const int64_t &pos, Tree *T, const int32_t& k,
    int64_t &dfsIndex, std::unordered_map<int64_t, tupleCoord_t> &scalarToTupleCoord,
    const sequence_t &sequence, const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const globalCoords_t &globalCoords, CoordNavigator &navigator,
    std::map<int64_t, int64_t> &gapRuns, const std::vector<std::pair<int64_t, int64_t>>& blockRanges);

// Fill mutation matrices from tree or file
std::pair<size_t, size_t> getMaskCoorsForMutmat(const std::string &s1,
                                                const std::string &s2,
                                                size_t window,
                                                double threshold);
// void fillMutationMatricesFromTree(mutationMatrices &mutMat, Tree *T,
//                                   size_t window, double threshold);
void fillMutationMatricesFromTree_test(mutationMatrices &mutMat, Tree *T, const std::string& path);
void fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf);

// Build mutation matrices by traversing through all parent-child pairs
void writeMutationMatrices(const mutationMatrices &mutMat,
                           std::ofstream &mmfout);



std::tuple<std::string, std::vector<int>, std::vector<int>, std::vector<int>> getNucleotideSequenceFromBlockCoordinates(
    tupleCoord_t &start, tupleCoord_t &end, const sequence_t &sequence,
    const blockExists_t &blockExists, const blockStrand_t &blockStrand,
    const Tree *T, const Node *node, const globalCoords_t &globalCoords, CoordNavigator &navigator);



std::string getStringAtNode(Node *node, Tree *T, bool aligned);



} // namespace seed_annotated_tree
#endif
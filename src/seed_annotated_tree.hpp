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

    // Default constructor
    tupleCoord_t() : blockId(-1), nucPos(-1), nucGapPos(-1) {}
    
    // Copy constructor
    tupleCoord_t(const tupleCoord_t &other) : blockId(other.blockId), nucPos(other.nucPos), nucGapPos(other.nucGapPos) {}
    
    // Add constructor for direct initialization
    tupleCoord_t(int64_t b, int64_t n, int64_t g) : blockId(b), nucPos(n), nucGapPos(g) {}

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
    
    // Struct to hold seed lookup info
    struct SeedInfo {
        hashedKmer_t kmer;
        int64_t endPos;
        bool isReverse;
        
        bool operator==(const SeedInfo& other) const {
            return kmer == other.kmer && 
                   endPos == other.endPos && 
                   isReverse == other.isReverse;
        }
    };
    
    struct SeedInfoHash {
        size_t operator()(const SeedInfo& info) const {
            size_t h = std::hash<hashedKmer_t>{}(info.kmer);
            h ^= std::hash<int64_t>{}(info.endPos) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<bool>{}(info.isReverse) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    // Hash function for position-dfsIndex pair
    struct PosIndexHash {
        size_t operator()(const std::pair<int64_t, int64_t>& p) const {
            size_t h = std::hash<int64_t>{}(p.first);
            h ^= std::hash<int64_t>{}(p.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    // Maps seed info to its access count
    std::unordered_map<SeedInfo, int64_t, SeedInfoHash> seedAccessCounts;
    
    // Fast lookup from (pos, dfsIndex) to seed info
    std::unordered_map<std::pair<int64_t, int64_t>, SeedInfo, PosIndexHash> posToSeed;
    
    size_t maxCachedSeeds;

    HotSeedIndex(int k, size_t maxSeeds = 10000) : k(k), maxCachedSeeds(maxSeeds) {}

    // Record a seed access during build phase
    void recordSeedAccess(int64_t pos, int64_t dfsIndex, hashedKmer_t kmer, int64_t endPos, bool isReverse) {
        SeedInfo info{kmer, endPos, isReverse};
        auto& count = seedAccessCounts[info];
        count++;

        // Add/update position mapping
        posToSeed[{pos, dfsIndex}] = info;

        // If we've exceeded maxCachedSeeds, remove least frequently accessed
        if (seedAccessCounts.size() > maxCachedSeeds) {
            // Find least accessed seed
            auto minIt = std::min_element(seedAccessCounts.begin(), seedAccessCounts.end(),
                [](const auto& a, const auto& b) {
                    return a.second < b.second;
                });

            // Remove all position mappings for this seed
            for (auto it = posToSeed.begin(); it != posToSeed.end();) {
                if (it->second == minIt->first) {
                    it = posToSeed.erase(it);
                } else {
                    ++it;
                }
            }

            // Remove the seed itself
            seedAccessCounts.erase(minIt);
        }
    }

    // Try to get a cached seed during placement phase
    bool getSeed(int64_t pos, int64_t dfsIndex, hashedKmer_t& kmer, int64_t& endPos, bool& isReverse) {
        auto it = posToSeed.find({pos, dfsIndex});
        if (it != posToSeed.end()) {
            const auto& info = it->second;
            kmer = info.kmer;
            endPos = info.endPos;
            isReverse = info.isReverse;
            return true;
        }
        return false;
    }

    // Serialize to Cap'n Proto
    void serialize(::HotSeedIndexSerial::Builder builder) const {
        auto entriesBuilder = builder.initEntries(seedAccessCounts.size());
        size_t i = 0;
        
        // Group positions by seed info
        std::unordered_map<SeedInfo, std::vector<std::pair<int64_t, int64_t>>, SeedInfoHash> seedToPositions;
        for (const auto& [pos_idx, seed] : posToSeed) {
            seedToPositions[seed].push_back(pos_idx);
        }

        for (const auto& [info, count] : seedAccessCounts) {
            auto entry = entriesBuilder[i];
            entry.setKmer(info.kmer);
            entry.setEndPos(info.endPos);
            entry.setIsReverse(info.isReverse);
            entry.setAccessCount(count);
            
            const auto& positions = seedToPositions[info];
            auto positionsBuilder = entry.initPositions(positions.size());
            for (size_t j = 0; j < positions.size(); j++) {
                positionsBuilder[j].setPos(positions[j].first);
                positionsBuilder[j].setDfsIndex(positions[j].second);
            }
            i++;
        }
    }

    // Deserialize from Cap'n Proto
    void deserialize(::HotSeedIndexSerial::Reader reader) {
        seedAccessCounts.clear();
        posToSeed.clear();
        
        auto entries = reader.getEntries();
        for (auto entry : entries) {
            SeedInfo info{entry.getKmer(), entry.getEndPos(), entry.getIsReverse()};
            seedAccessCounts[info] = entry.getAccessCount();
            
            auto positions = entry.getPositions();
            for (auto pos : positions) {
                posToSeed[{pos.getPos(), pos.getDfsIndex()}] = info;
            }
        }
    }

private:
    int k;
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
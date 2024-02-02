#pragma once
#include "PangenomeMAT.hpp"
#include "seed.hpp"
#include <iostream>
#include <vector>
#include "minimap2_src/mmpriv.h"
#include "minimap2_src/minimap.h"
#include "minimap2_src/kseq.h"
#include "minimap2_src/kalloc.h"
#include "minimap2_src/khash.h"
#include "minimap2_src/kvec.h"

using namespace seed;

/* Helpers for interacting with panmats */
namespace tree {
    using namespace PangenomeMAT;

    typedef std::vector< std::pair< std::vector< std::pair< int, std::vector< int > > >, std::vector< std::vector< std::pair< int, std::vector< int > > > > > > globalCoords_t;
    typedef std::vector< std::tuple< int32_t, bool, bool, bool, bool > > blockMutData_t;
    typedef std::vector<std::tuple< int32_t, int32_t, int32_t, char, char, int32_t, int32_t> > nucMutData_t;
    typedef std::tuple<int32_t, int32_t, int32_t> range_t;


    struct mutableTreeData { 
        // These fields are intendet to be mutated at each node during a DFS
        sequence_t sequence; // the main object encoding the MSA
        std::string gappedConsensus;
        std::string ungappedConsensus;
        std::vector<int32_t> degap; // ungapped[degap[i]] => gapped[i]
        std::vector<int32_t> regap; // gapped[regap[i]] => ungapped[i]
        std::vector<kmer_t> seeds; // dynamic vector of seeds in each node's sequence
        std::vector<jkmer> jkmers;
        std::unordered_map<int32_t, std::pair<int32_t, std::string>> seedMap;
        std::unordered_map<std::string, bool> variableSeeds; // seeds in the consensus that mutate at least once
        blockExists_t blockExists; // tracks if blocks are "on" at a node
        blockStrand_t blockStrand; // tracks strand of blocks
    };
    
    /* Interface */

    void updateConsensus(mutableTreeData &data, Tree *T);
   
    void removeIndices(std::vector<kmer_t>& v, std::stack<int32_t>& rm);
    std::string getConsensus(Tree *T); // ungapped!

    std::unordered_map<std::string, std::string> getAllNodeStrings(Tree *T);
    std::string getStringFromCurrData(mutableTreeData data, Tree *T, const Node *node, const bool aligned);

    size_t getGlobalCoordinate(const int blockId, const int nucPosition, const int nucGapPosition, const globalCoords_t &globalCoords);
    void setup(mutableTreeData &data, globalCoords_t &globalCoords, Tree *T);
}
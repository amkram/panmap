#ifndef __TREE_HPP
#define __TREE_HPP

#pragma once
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include <iostream>
#include <vector>
#include "minimap2_src/mmpriv.h"
#include "minimap2_src/minimap.h"
#include "minimap2_src/kseq.h"
#include "minimap2_src/kalloc.h"
#include "minimap2_src/khash.h"
#include "minimap2_src/kvec.h"



void time_stamp();




using namespace seeding;

inline auto seed_cmp = [](const std::pair<int32_t, std::string> &a, const std::pair<int32_t, std::string> &b) {
    if(a.second != b.second) {
        return a.second < b.second;
    }
    return a.first < b.first;
};
typedef std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>> seedmerIndex_t;

/* Helpers for interacting with panmats */
namespace tree {
    using namespace PangenomeMAT;

    typedef std::vector< std::pair< std::vector< std::pair< int, std::vector< int > > >, std::vector< std::vector< std::pair< int, std::vector< int > > > > > > globalCoords_t;
    typedef std::vector< std::tuple< int32_t, bool, bool, bool, bool > > blockMutData_t;
    typedef std::vector<std::tuple< int32_t, int32_t, int32_t, char, char, int32_t, int32_t> > nucMutData_t;
    typedef std::tuple<int32_t, int32_t, int32_t> range_t;


    struct mutableTreeData { 
        // These fields are intended to be mutated at each node during a DFS
        sequence_t sequence; // the main object encoding the MSA
        std::string gappedConsensus;
        std::string ungappedConsensus;
        std::vector<int32_t> degap; // ungapped[degap[i]] => gapped[i]
        std::vector<int32_t> regap; // gapped[regap[i]] => ungapped[i]
        std::vector<seed> seeds; // dynamic vector of seeds in each node's sequence
        std::vector<seedmer> seedmers;
        std::unordered_map<int32_t, std::pair<int32_t, std::string>> seedMap;
        std::unordered_map<std::string, bool> variableSeeds; // seeds in the consensus that mutate at least once
        blockExists_t blockExists; // tracks if blocks are "on" at a node
        blockStrand_t blockStrand; // tracks strand of blocks
    };
    
    struct mutationMatrices {
        // Store mutation matrices
        std::vector< std::vector<double> > submat; // 4 x 4 substitution rate matrix
        std::vector<double> insmat = {0}; // 1 x N insertion rate by length matrix
        std::vector<double> delmat = {0}; // 1 x N deletion rate by length matrix
        
        // Stores total number of mutations
        bool filled = false;
        std::vector<double> total_submuts;
        double total_insmut = 0;
        double total_delmut = 0;
        
        mutationMatrices() {
            // initialize mutationMatrices object and intialize the correct size for substitution amtrix
            total_submuts.resize(4);
            submat.resize(4);
            for (size_t i = 0; i < 4; ++i) {
                submat[i].resize(4);
            }
        }  
    };

    /* Interface */
    void updateConsensus(mutableTreeData &data, Tree *T, std::map<int32_t, int32_t> *coordsIndex=nullptr);
   
    void removeIndices(std::vector<seed>& v, std::stack<int32_t>& rm);
    std::string getConsensus(Tree *T); // ungapped!

    std::unordered_map<std::string, std::string> getAllNodeStrings(Tree *T);
    std::string getStringFromCurrData(mutableTreeData &data, Tree *T, const Node *node, const bool aligned);

    size_t getGlobalCoordinate(const int blockId, const int nucPosition, const int nucGapPosition, const globalCoords_t &globalCoords);
    void setup(mutableTreeData &data, globalCoords_t &globalCoords, Tree *T);

    // Fill mutation matrices from tree or file
    std::pair<size_t, size_t> getMaskCoorsForMutmat(const std::string& s1, const std::string& s2, size_t window, double threshold);
    void fillMutationMatricesFromTree(mutationMatrices& mutMat, Tree* T, size_t window, double threshold);
    void fillMutationMatricesFromFile(mutationMatrices& mutMat, std::ifstream& inf);

    // Build mutation matrices by traversing through all parent-child pairs
    void writeMutationMatrices(const mutationMatrices& mutMat, std::ofstream& mmfout);
}

#endif
#pragma once
#include "PangenomeMAT.hpp"
#include "tree.hpp"

using namespace PangenomeMAT;
using namespace seed;
using namespace tree;

typedef std::unordered_map<int32_t, std::pair<int32_t, std::string>> seedMap_t;

namespace pmi { // functions and types for seed indexing

    struct seedIndex {
        std::stringstream outStream;
        std::vector<kmer_t> consensusSeeds;
        std::unordered_map<int32_t, std::string> consensusJkmers;
        std::unordered_map<int32_t, std::string> currJkmers;
        std::unordered_map<std::string, std::vector<kmer_t>> deletions;
        std::unordered_map<std::string, std::vector<kmer_t>> insertions;
        int32_t k;
        int32_t s;
    };
  
    /* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
    void build(seedIndex &index, Tree *T, const size_t l, const size_t k, const size_t s);

    /* Writes a seedIndex to fout */
    void write(std::ofstream &fout, Tree *T, seedIndex &index);

    /* Loads indexFile, storing it in index  */
    void load(seedIndex &index, std::ifstream &indexFile);
}

/* Expose some pmi.cpp helpers for unit testing */
using namespace pmi;

void buildHelper(mutableTreeData &data, seedMap_t seedMap, seedIndex &index, Tree *T, const Node *node, const int32_t l, const size_t k, const size_t s, const globalCoords_t &globalCoords);
void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData);
void applyMutations(mutableTreeData &data, blockMutData_t &blockMutData, nucMutData_t &nucMutData, Tree *T, const Node *node, const globalCoords_t &globalCoords);
range_t getRecomputePositions(const range_t &p, const std::string &gappedSequence, const int32_t k);
std::vector<range_t> getAffectedRanges(mutableTreeData &data, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData, std::string &seq, Tree *T, const int32_t k, const globalCoords_t &globalCoords);
void addSeeds(std::vector<kmer_t> &seeds, seedIndex &index, const std::string nodeId, const std::unordered_map<std::string, kmer_t> &newSeeds);
void recomputeSeeds(mutableTreeData &data, std::unordered_map<std::string, kmer_t> &newSeeds, const std::vector<range_t> &ranges, const std::string &sequence, int32_t k, int32_t s);
void discardSeeds(std::vector<kmer_t> &seeds, std::unordered_map<std::string, kmer_t> &newSeeds, seedIndex &index, const std::vector<range_t>& B, const std::string &seq, const std::string nid, const size_t k);
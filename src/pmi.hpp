#pragma once
#include "PangenomeMAT.hpp"
#include "tree.hpp"

using namespace PangenomeMAT;
using namespace seed;
using namespace tree;


namespace pmi { // functions and types for seed indexing

    struct seedIndex {
        std::vector<kmer_t> consensusSeeds;
        std::unordered_map<std::string, std::vector<kmer_t>> deletions;
        std::unordered_map<std::string, std::vector<kmer_t>> insertions;
    };
  
    /* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
    void build(seedIndex &index, Tree *T, const size_t k, const size_t s);
    /* Writes a seedIndex to fout */
    void write(std::ofstream &fout, Tree *T, seedIndex &index);
    /* Loads indexFile, storing it in seedIndex si */
    void load(seedIndex &index, const Node *root, const std::ifstream &indexFile);
}

/* Expose some pmi.cpp helpers for unit testing */
using namespace pmi;

void buildHelper(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const size_t k, const size_t s, const globalCoords_t &globalCoords);
void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData);
void applyMutations(mutableTreeData &data, blockMutData_t &blockMutData, nucMutData_t &nucMutData, Tree *T, const Node *node, const globalCoords_t &globalCoords);

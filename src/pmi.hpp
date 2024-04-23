#ifndef __PMI_HPP
#define __PMI_HPP

#pragma once
#include "seeding.hpp"
#include "tree.hpp"

using namespace PangenomeMAT;
using namespace seeding;
using namespace tree;

typedef std::unordered_map<int32_t, std::pair<int32_t, std::string>> seedMap_t;
namespace pmi { // functions and types for seed indexing

    struct seedIndex {
        std::stringstream outStream;
        std::vector<seed> consensusSeeds;
        std::unordered_map<int32_t, std::string> consensusseedmers;
        std::unordered_map<int32_t, std::string> currseedmers;
        std::unordered_map<std::string, std::vector<seed>> deletions;
        std::unordered_map<std::string, std::vector<seed>> insertions;
        int32_t k;
        int32_t s;
    };

    /* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
    void build(seedIndex &index, Tree *T, const size_t j, const size_t k, const size_t s);
}

/* Expose some pmi.cpp helpers for unit testing */
using namespace pmi;

void buildHelper(mutableTreeData &data, seedMap_t seedMap, seedIndex &index, Tree *T, const Node *node, const int32_t l, const size_t k, const size_t s, const globalCoords_t &globalCoords);
void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData);
void applyMutations(mutableTreeData &data, blockMutData_t &blockMutData, nucMutData_t &nucMutData, Tree *T, const Node *node, const globalCoords_t &globalCoords);
range_t getRecomputePositions(const range_t &p, const std::string &gappedSequence, const int32_t k);
std::vector<range_t> getAffectedRanges(mutableTreeData &data, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData, std::string &seq, Tree *T, const int32_t k, const globalCoords_t &globalCoords);

#endif
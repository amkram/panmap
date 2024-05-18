#ifndef __PMI_HPP
#define __PMI_HPP

#include <__config>
#pragma once
#include "seeding.hpp"
#include "tree.hpp"

using namespace PangenomeMAT;
using namespace seeding;
using namespace tree;

typedef std::map<int64_t, std::string> seedMap_t;
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
        int32_t j;
    };

    /* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
    void build(seedIndex &index, Tree *T, const size_t j, const size_t k, const size_t s);
}

/* Expose some pmi.cpp helpers for unit testing also tree.cpp uses applyMutations for now */
using namespace pmi;

void buildHelper(mutableTreeData &data, seedMap_t &seedMap, seedIndex &index, Tree *T, const Node *node, const globalCoords_t &globalCoords);
void applyMutations(mutableTreeData &data, blockMutData_t &blockMutData, std::set<std::tuple<int, int, int, int>, decltype(rangeCmp)> &recompRanges, nucMutData_t &nucMutData, Tree *T, const Node *node, const globalCoords_t &globalCoords, seedIndex &index);
void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData);

#endif
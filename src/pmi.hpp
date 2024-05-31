#ifndef __PMI_HPP
#define __PMI_HPP

#include <__config>
#pragma once
#include "seeding.hpp"
#include "tree.hpp"
#include "index.pb.h"

using namespace PangenomeMAT;
using namespace seeding;
using namespace tree;

typedef std::map<tupleCoord_t, std::string> seedMap_t;

namespace pmi { // functions and types for seed indexing

/* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
void build(SeedmerIndex &index, Tree *T, int j, int k, int s);

} // namespace pmi

/* Expose some pmi.cpp helpers for unit testing also tree.cpp uses
 * applyMutations for now */
using namespace pmi;

void buildHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index,
                 Tree *T, const Node *node, const globalCoords_t &globalCoords,
                 CoordNavigator &navigator);
void applyMutations(mutableTreeData &data, seedMap_t &seedMap,
                    blockMutData_t &blockMutData,
                    std::vector<tupleRange> &recompRanges,
                    nucMutData_t &nucMutData, Tree *T, const Node *node,
                    globalCoords_t &globalCoords, const SeedmerIndex &index);
void undoMutations(mutableTreeData &data, SeedmerIndex &index, Tree *T,
                   const Node *node, const blockMutData_t &blockMutData,
                   const nucMutData_t &nucMutData);


#endif
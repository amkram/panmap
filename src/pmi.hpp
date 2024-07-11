#ifndef __PMI_HPP
#define __PMI_HPP

//#include <__config>
#pragma once
#include "seeding.hpp"
#include "tree.hpp"
#include "index.pb.h"
#include "mgsr.hpp"
#include <unordered_map>

using namespace PangenomeMAT;
using namespace seeding;
using namespace tree;

//typedef std::unordered_map<tupleCoord_t, std::string, TupleHash> seedMap_t;
typedef std::unordered_map<int, std::string> seedMap_t;

namespace pmi { // functions and types for seed indexing

/* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
void build(SeedmerIndex &index, Tree *T, int j, int k, int s);

void buildSeedmers(SeedmerIndex &index, Tree *T, int j, int k, int s, std::stringstream& seedmersOutStream);

} // namespace pmi

/* Expose some pmi.cpp helpers for unit testing also tree.cpp uses
 * applyMutations for now */
using namespace pmi;

void buildHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index,
                 Tree *T, const Node *node, const globalCoords_t &globalCoords,
                 CoordNavigator &navigator);
void applyMutations(mutableTreeData &data, seedMap_t &seedMap,
                    blockMutationInfo_t &blockMutData,
                    std::vector<tupleRange> &recompRanges,
                    mutationInfo_t &nucMutData, Tree *T, const Node *node,
                    globalCoords_t &globalCoords, const SeedmerIndex &index);
void undoMutations(mutableTreeData &data, SeedmerIndex &index, Tree *T,
                   const Node *node, const blockMutationInfo_t &blockMutData,
                   const mutationInfo_t &nucMutData);


// Go upstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandLeft(const CoordNavigator &navigator, tupleCoord_t coord,
                        int neededNongap, blockExists_t &blockExists);

// Go downstream until neededNongap nucleotides are seen and return the new coord.
tupleCoord_t expandRight(const CoordNavigator &navigator, tupleCoord_t coord,
                         int neededNongap, blockExists_t &blockExists);

// Merges each range with overlapping ranges after expanding left and right
// by `neededNongap` non-gap nucleotides.
std::vector<tupleRange> expandAndMergeRanges(const CoordNavigator &navigator, std::vector<tupleRange> &ranges, int neededNongap, blockExists_t &blockExists);
int64_t tupleToScalarCoord(const tupleCoord_t &coord, const globalCoords_t &globalCoords);

void buildSeedmersHelper(mutableTreeData &data, seedMap_t &seedMap, SeedmerIndex &index, mgsr::seedmers& seedmersIndex,
                 std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, Tree *T, Node *node, globalCoords_t &globalCoords, CoordNavigator &navigator,
                 std::vector<int> &scalarCoordToBlockId, std::vector<std::unordered_set<int>> &BlocksToSeeds, std::stringstream& seedmersOutStream, const std::vector<std::pair<size_t, size_t>>& blocksEndCoords, size_t& parentUngappedSeqSize, std::map<int32_t, int32_t>& parentCoordsIndex);

#endif
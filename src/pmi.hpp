#ifndef __PMI_HPP
#define __PMI_HPP

//#include <__config>
#pragma once
#include "seeding.hpp"
#include "seed_annotated_tree.hpp"
#include "index.capnp.h"
#include <unordered_map>

using namespace panmanUtils;
using namespace seeding;
using namespace seed_annotated_tree;




enum posWidth {pos16, pos32, pos64};

namespace pmi { // functions and types for seed indexing

    void build(Tree *T, Index::Builder &index);
    void place(
      Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path,
      seed_annotated_tree::mutationMatrices &mutMat, std::string prefix, std::string refFileName, std::string samFileName,
      std::string bamFileName, std::string mpileupFileName, std::string vcfFileName, std::string aligner,
      const std::string& refNode, const bool& save_jaccard
    );

    void parallel_tester(Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path, const std::string &prefix);

    void place_per_read(
      Tree *T, Index::Reader &index, const std::string &reads1Path, const std::string &reads2Path,
      const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate,
      const int& redoReadThreshold, const bool& recalculateScore, const bool& rescueDuplicates,
      const double& rescueDuplicatesThreshold, const double& excludeDuplicatesThreshold,
      const std::string& preEMFilterMethod, const int& preEMFilterNOrder, const int& preEMFilterMBCNum, const int& emFilterRound, const int& checkFrequency,
      const int& removeIteration, const double& insigProb, const int& roundsRemove, const double& removeThreshold,
      const bool& leafNodesOnly, const bool& callSubconsensus, const std::string& prefix, const bool& save_kminmer_binary_coverage);

} // namespace pmi

/* Expose some pmi.cpp helpers for unit testing also tree.cpp uses
 * applyMutations for now */
using namespace pmi;

// void buildHelper(mutableTreeData &data, seedMap_t &seedMap, int32_t &k, int32_t &s, ::capnp::List<Mutations>::Builder &mutations,
//                  Tree *T, const Node *node, const globalCoords_t &globalCoords,
//                  CoordNavigator &navigator, int64_t &dfsIndex, posWidth &width, std::map<int64_t, int64_t> &gapRuns);
// void applyMutations(mutableTreeData &data, seedMap_t &seedMap,
//                     blockMutationInfo_t &blockMutData,
//                     std::vector<tupleRange> &recompRanges,
//                     mutationInfo_t &nucMutData, Tree *T, const Node *node,
//                     globalCoords_t &globalCoords, const ::capnp::List<Mutations>::Builder &mutations);
// void undoMutations(mutableTreeData &data, ::capnp::List<Mutations>::Builder &mutations, Tree *T,
//                     Node *node, const blockMutationInfo_t &blockMutData,
//                    const mutationInfo_t &nucMutData);
void flipCoords(int32_t blockId, globalCoords_t &globalCoords);
void updateGapMapStep(std::map<int64_t, int64_t>& gapMap, const std::pair<bool, std::pair<int64_t, int64_t>>& update, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates, bool recordGapMapUpdates=true);
void updateGapMap(std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& updates, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates);
std::vector<std::pair<int64_t, int64_t>> invertRanges(const std::vector<std::pair<int64_t, int64_t>>& nucRanges, const std::pair<int64_t, int64_t>& invertRange);
void invertGapMap(std::map<int64_t, int64_t>& gapMap, const std::pair<int64_t, int64_t>& invertRange, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates);
void makeCoordIndex(std::map<int64_t, int64_t>& coordIndex, const std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<int64_t, int64_t>>& blockRanges);

// // Go upstream until neededNongap nucleotides are seen and return the new coord.
// tupleCoord_t expandLeft(const CoordNavigator &navigator, tupleCoord_t coord,
//                         int neededNongap, blockExists_t &blockExists);

// // Go downstream until neededNongap nucleotides are seen and return the new coord.
// tupleCoord_t expandRight(const CoordNavigator &navigator, tupleCoord_t coord,
//                          int neededNongap, blockExists_t &blockExists, blockStrand_t &blockStrand);

// // Merges each range with overlapping ranges after expanding left and right
// // by `neededNongap` non-gap nucleotides.
// std::vector<tupleRange> expandAndMergeRanges(const CoordNavigator &navigator, std::vector<tupleRange> &ranges, int neededNongap, blockExists_t &blockExists);
int64_t tupleToScalarCoord(const tupleCoord_t &coord, const globalCoords_t &globalCoords);


#endif
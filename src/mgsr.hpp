#ifndef __MGSR_HPP
#define __MGSR_HPP

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "PangenomeMAT.hpp"
#include "tree.hpp"
#include "pmi.hpp"



typedef size_t hash_t;

namespace mgsr {
    struct seedmers {
        //       beg                hash    rev
        std::map<int32_t, std::pair<size_t, bool>> positionMap;
        //                 hash                       begs
        std::unordered_map<size_t, std::unordered_set<int32_t>> seedmersMap;
    };

    // pmi index file with k = k, s = s, l = 1
    // k
    // l
    // tree
    void accio(PangenomeMAT::Tree *T, std::ifstream& indexFile, size_t k, size_t l);
    void buildSeedmer(pmi::seedIndex &Index, PangenomeMAT::Tree *T, const size_t l, const size_t k, const size_t s, std::stringstream& seedmersOutStream);
    void buildSeedmerHelper(tree::mutableTreeData &data, seedMap_t seedMap, pmi::seedIndex &index, mgsr::seedmers seedmersIndex, std::map<int32_t, std::pair<int32_t, size_t>> curSeeds, PangenomeMAT::Tree *T, const PangenomeMAT::Node *node, const int32_t l, const int32_t k, const int32_t s, const tree::globalCoords_t &globalCoords, std::stringstream& seedmersOutStream);
}




#endif

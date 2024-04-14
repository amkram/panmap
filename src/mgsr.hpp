#ifndef __MGSR_HPP
#define __MGSR_HPP

#include <unordered_map>
#include "PangenomeMAT.hpp"
#include "tree.hpp"
#include "pmi.hpp"



typedef size_t hash_t;

namespace mgsr {
    struct seedmer;
    struct seedmer {
        hash_t hash; // hash
        int32_t beg; // starting position
        int32_t end; // ending position
        size_t num; // number of seedmers with the same hash
        bool   rev; // reversed
        seedmer* prev; // ptr to previous seedmer
        seedmer* next; // ptr to next seedmer
        
    };

    struct seedmers {
        std::unordered_map<hash_t, seedmer> seedmerMap; // seedmers
        seedmer* firstSeedmer; // starting kmm
    };

    // pmi index file with k = k, s = s, l = 1
    // k
    // l
    // tree
    void accio(PangenomeMAT::Tree *T, std::ifstream& indexFile, size_t k, size_t l);
    void buildSeedmer(pmi::seedIndex &Index, PangenomeMAT::Tree *T, const size_t l, const size_t k, const size_t s);
    void buildSeedmerHelper(tree::mutableTreeData &data, seedMap_t seedMap, pmi::seedIndex &index, mgsr::seedmers seedmersIndex, PangenomeMAT::Tree *T, const PangenomeMAT::Node *node, const int32_t l, const int32_t k, const int32_t s, const tree::globalCoords_t &globalCoords);
}




#endif
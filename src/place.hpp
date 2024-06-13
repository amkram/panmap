#pragma once
#include "tree.hpp"
#include "index.pb.h"

namespace place {
    
    enum class ScoringMethod {
        NUM_SEED_HITS, // raw number of matches per node
        JACCARD // Todo implement
    };

    void placeIsolate( SeedmerIndex &index, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, PangenomeMAT::Tree *T,  bool use_root);

}
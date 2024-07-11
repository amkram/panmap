#pragma once
#include "tree.hpp"
#include "index.pb.h"

namespace place {
    
    enum class ScoringMethod {
        NUM_SEED_HITS, // raw number of matches per node
        JACCARD // Todo implement
    };

    void placeIsolate( SeedmerIndex &index, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, PangenomeMAT::Tree *T,  bool use_root);
    void placeMetagenomics(PangenomeMAT::Tree *T, const tree::mutationMatrices& mutMat, const int32_t& accioK, const int32_t& accioS, const int32_t& accioL, const std::string& defaultKmiPath, const std::string& reads1File, const std::string& reads2File, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, const std::string& prefix, const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate, const bool& confidence, const int32_t& roundsRemove, const double& removeThreshold);

}
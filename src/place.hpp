#ifndef __PLACE_HPP
#define __PLACE_HPP

#pragma once
#include "tree.hpp"
#include "PangenomeMAT.hpp"

namespace place {

    void placeIsolate(std::ifstream &indexFile, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, const std::string& prefix, const bool& makeSam, const bool& makeBam, const bool& makeMPileup, const bool& makeVCF, const bool& makeRef, PangenomeMAT::Tree *T,  bool use_root);
    void placeAccio(PangenomeMAT::Tree *T, const tree::mutationMatrices& mutMat, const int32_t& accioK, const int32_t& accioS, const int32_t& accioL, const std::string& defaultKmiPath, const std::string& reads1File, const std::string& reads2File, const std::string& prefix, const bool& makeSam, const bool& makeBam, const bool& makeMPileup, const bool& makeVCF, const bool& makeRef, const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate, const bool& confidence, const int32_t& roundsRemove, const double& removeThreshold);
}
#endif
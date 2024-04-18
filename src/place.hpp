#ifndef __PLACE_HPP
#define __PLACE_HPP

#pragma once
#include "tree.hpp"
#include "PangenomeMAT.hpp"

namespace place {

    void placeIsolate( std::ifstream &indexFile, const tree::mutationMatrices& mutMat, const std::string &reads1Path, const std::string &reads2Path, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName, std::string &refFileName, PangenomeMAT::Tree *T,  bool use_root);

}

#endif
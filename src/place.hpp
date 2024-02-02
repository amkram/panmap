#pragma once
#include "tree.hpp"

namespace place {

    void placeIsolate( std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path, PangenomeMAT::Tree *T, int32_t k, int32_t s);

}
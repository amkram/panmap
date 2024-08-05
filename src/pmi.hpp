#ifndef __PMI_HPP
#define __PMI_HPP

//#include <__config>
#pragma once
#include "seeding.hpp"
#include "tree.hpp"
#include "index.capnp.h"
#include <unordered_map>

using namespace PangenomeMAT;
using namespace seeding;
using namespace tree;



enum posWidth {pos16, pos32, pos64};

namespace pmi { // functions and types for seed indexing

    /* Indexes T with syncmers parameterized by (k,s). Stores result in si. */
    void build(Tree *T, Index::Builder &index);

} // namespace pmi

/* Expose some pmi.cpp helpers for unit testing also tree.cpp uses
 * applyMutations for now */
using namespace pmi;

int64_t tupleToScalarCoord(const tupleCoord_t &coord, const globalCoords_t &globalCoords);


#endif
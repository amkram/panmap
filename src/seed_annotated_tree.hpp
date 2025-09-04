#pragma once

#include "coordinates.hpp"  // Include this file for full definition of CoordinateTraverser
#include "gap_map.hpp"
#include "index.capnp.h"
#include "seeding.hpp"
#include <capnp/common.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <capnp/serialize.h>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <absl/container/flat_hash_map.h> // Add include
#include <absl/container/flat_hash_set.h> // Add include

using namespace seeding;
using namespace panmanUtils;

// Forward declare types from coordinates namespace
namespace coordinates {
class CoordinateManager;
class CoordinateTraverser;
struct tupleCoord_t;
struct CoordRange;
struct GapUpdate; // Forward declare GapUpdate as well
using GapRange = std::pair<int64_t, int64_t>;
using GapMap = std::map<int64_t, int64_t>;
using sequence_t = std::vector<
    std::pair<std::vector<std::pair<char, std::vector<char>>>,
              std::vector<std::vector<std::pair<char, std::vector<char>>>>>>;
using blockExists_t = std::vector<std::pair<bool, std::vector<bool>>>;
using blockStrand_t = std::vector<std::pair<bool, std::vector<bool>>>;
} // namespace coordinates

// Now use the coordinate types
using coordinates::CoordRange;
using coordinates::GapMap;
using coordinates::GapRange;
using coordinates::GapUpdate; // Keep this line
using coordinates::sequence_t;
using coordinates::tupleCoord_t;

// Parameters for indexing and placement
struct PanmapParams {
  int32_t k; // k-mer length
  int32_t s; // syncmer length
  int32_t t; // threshold
  bool open; // open/closed syncmers
  int32_t l; // minimum length
};

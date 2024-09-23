#ifndef __MGSR_HPP
#define __MGSR_HPP

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tbb/concurrent_vector.h>
#include <boost/icl/interval_map.hpp>
#include <boost/icl/split_interval_map.hpp>
#include "panmanUtils.hpp"

typedef std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> readSeedmers_t;
typedef std::tuple<int32_t, int32_t, int32_t, int32_t, bool, int32_t> match_t;
namespace mgsr {

  struct positionInfo {
    int32_t endPos;
    size_t fhash;
    size_t rhash;
    bool rev;
  };

  struct seedmers {
      //       beg                 end      fhash    rhash    rev
      std::map<int32_t, positionInfo> positionMap;
      //                 hash                       begs
      std::unordered_map<size_t, std::set<int32_t>> hashToPositionsMap;
  };

  class Read {
    public:
    std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> seedmers;
    std::unordered_map<size_t, std::vector<int32_t>> seedmersIndex;
    boost::icl::split_interval_map<int32_t, int> matches;
    std::unordered_set<int32_t> duplicates;
    // std::vector<bool> duplicates;
    // std::vector<bool> absentees;
    // size_t numDuplicates = 0;
    // size_t numAbsentees = 0;

    Read(std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>>&& _seedmers)
      : seedmers(std::move(_seedmers)) {
        // duplicates.resize(seedmers.first.size(), false);
        // absentees.resize(seedmers.first.size(), false);
    }
  }; 
}

#endif

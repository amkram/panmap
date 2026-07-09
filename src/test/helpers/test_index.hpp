#pragma once

/**
 * @file test_index.hpp
 * @brief TestIndex RAII fixture + loadIndex(): the one place tests build/load a
 *        panmap index.
 */

#include "panmanUtils.hpp"
#include "panmap_utils.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace ts {

// Owns a loaded index. The LiteIndex seed-change struct-of-arrays (stored segmented
// as List(List(...)) in the file) are flattened into plain vectors on load, so tests
// address them with a flat index and never touch capnp readers or re-traverse the
// message. The LiteTree is self-contained after initialize() (copies structure).
class IndexData {
   public:
    std::unique_ptr<panmapUtils::LiteTree> liteTree;
    int k = 0, s = 0, t = 0, l = 0;
    bool open = false, hpc = false;
    int flankMaskBp = 0;  // not stored in the index; set by TestIndex from the build param

    std::vector<uint64_t> hashes;
    std::vector<int64_t> parentCounts;
    std::vector<int64_t> childCounts;
    std::vector<uint32_t> offsets;  // size == numNodes + 1

    size_t numNodesPlusOne() const { return offsets.size(); }

    size_t numChanges() const { return offsets.empty() ? 0 : offsets.back(); }

    uint64_t hashAt(size_t i) const { return hashes[i]; }

    int64_t parentCountAt(size_t i) const { return parentCounts[i]; }

    int64_t childCountAt(size_t i) const { return childCounts[i]; }

    uint32_t nodeChangeOffset(size_t node) const { return offsets[node]; }
};

// Decompress + parse an existing .idx, returning an owned IndexData.
IndexData loadIndex(const std::string& path);

// Build an index from a tree into a temp file, load it, and clean up on destruction.
// numThreads <= 1 uses the sequential buildIndex(); > 1 uses buildIndexParallel().
class TestIndex {
   public:
    TestIndex(panmanUtils::Tree* tree,
              int k = 15,
              int s = 8,
              int t = 0,
              int l = 3,
              bool open = false,
              bool hpc = false,
              int flankMaskBp = 0,
              int numThreads = 1,
              int zstdLevel = 1);
    ~TestIndex();

    TestIndex(const TestIndex&) = delete;
    TestIndex& operator=(const TestIndex&) = delete;

    IndexData& data() { return data_; }

    const std::string& path() const { return path_; }

   private:
    std::string path_;
    IndexData data_;
};

}  // namespace ts

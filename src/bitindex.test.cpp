#include "pmi.hpp"
#include "panmanUtils.hpp"
#include "seeding.hpp"
#include "seed_annotated_tree.hpp"
#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string>
#include <unordered_set>
#include "bm.h"
#include "bmserial.h"
#include "bmsparsevec_serial.h"
#include "bmsparsevec_util.h"
#include "bmintervals.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmstrsparsevec.h"
#include "bmaggregator.h"
#include <vector>
#include <fstream>
#include <memory>
#include <variant>

void loadIndexFromFile(const std::string& filename, bm::bvector<>& allBitsets) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    std::vector<unsigned char> compressedData((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    bm::deserialize(allBitsets, compressedData.data());
    ifs.close();
}

void decodeBitsets(const bm::bvector<>& allBitsets, std::vector<bm::bvector<>>& compressedSeedDeltas, size_t numBitsets, size_t bitsetSize) {
    compressedSeedDeltas.resize(numBitsets);
    for (size_t i = 0; i < numBitsets; ++i) {
        bm::bvector<> bitset;
        for (size_t j = 0; j < bitsetSize; ++j) {
            if (allBitsets.test(i * bitsetSize + j)) {
                bitset.set(j);
            }
        }
        compressedSeedDeltas[i] = std::move(bitset);
    }
}

void readGapMapFromFile(std::vector<std::map<int64_t, int64_t>>& gapMaps, const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    while (ifs.peek() != EOF) {
        int64_t size;
        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
        std::map<int64_t, int64_t> gapMap;
        for (int64_t i = 0; i < size; ++i) {
            int64_t key, value;
            ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
            ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
            gapMap[key] = value;
        }
        gapMaps.push_back(gapMap);
    }
    ifs.close();
}

void traverseAndScore(Tree* T, const bm::bvector<>& allBitsets, const std::vector<std::map<int64_t, int64_t>>& gapMaps, size_t numBitsets, size_t bitsetSize) {
    std::vector<bm::bvector<>> compressedSeedDeltas;
    decodeBitsets(allBitsets, compressedSeedDeltas, numBitsets, bitsetSize);

    seed_annotated_tree::mutableTreeData data;
    seed_annotated_tree::globalCoords_t globalCoords;
    seed_annotated_tree::setup(data, globalCoords, T);

    CoordNavigator navigator(data.sequence);
    std::vector<int> BlockSizes(data.sequence.size(), 0);
    std::vector<std::pair<int64_t, int64_t>> blockRanges(data.blockExists.size());

    int32_t k = 0; // Set appropriate value
    int32_t s = 0; // Set appropriate value

    std::map<int64_t, int64_t> gapMap;
    std::vector<std::optional<int64_t>> gapVec;

    gapMap[0] = tupleToScalarCoord({blockRanges.size() - 1, globalCoords[blockRanges.size() - 1].first.size() - 1, -1}, globalCoords);

    std::vector<int> scalarCoordToBlockId(globalCoords.back().first.back().first + 1);
    auto currCoord = tupleCoord_t{0, 0, 0};
    if (navigator.sequence[0].first[0].second.empty()) {
        currCoord.nucGapPos = -1;
    }

    for (int64_t i = 0; i < scalarCoordToBlockId.size(); i++) {
        scalarCoordToBlockId[i] = currCoord.blockId;
        BlockSizes[currCoord.blockId]++;
        currCoord = navigator.newincrement(currCoord, data.blockStrand);
    }

    for (int64_t i = 0; i < blockRanges.size(); ++i) {
        int64_t start = globalCoords[i].first[0].second.empty() ? tupleToScalarCoord({i, 0, -1}, globalCoords) : tupleToScalarCoord({i, 0, 0}, globalCoords);
        int64_t end = tupleToScalarCoord({i, globalCoords[i].first.size() - 1, -1}, globalCoords);
        blockRanges[i] = std::make_pair(start, end);
    }

    std::vector<std::unordered_set<int>> BlocksToSeeds(data.sequence.size());
    posWidth width = globalCoords.back().first.back().first < 4294967296 ? posWidth::pos32 : posWidth::pos64;

    switch (width) {
        case posWidth::pos32:
            // Set width to 32
            break;
        case posWidth::pos64:
            // Set width to 64
            break;
    }

    tupleCoord_t coord = {0, 0, globalCoords[0].first[0].second.empty() ? -1 : 0};
    auto curIt = gapMap.end();

    while (coord < tupleCoord_t{-1, -1, -1}) {
        char c = coord.nucGapPos == -1 ? data.sequence[coord.blockId].first[coord.nucPos].first : data.sequence[coord.blockId].first[coord.nucPos].second[coord.nucGapPos];
        int64_t scalar = tupleToScalarCoord(coord, globalCoords);
        if (c == '-' || c == 'x') {
            if (!gapMap.empty() && curIt->second + 1 == scalar) {
                ++curIt->second;
            } else {
                auto tmpIt = gapMap.emplace(scalar, scalar);
                curIt = tmpIt.first;
            }
        }
        coord = navigator.newincrement(coord, data.blockStrand);
    }

    if (coord.blockId != -1 && !data.blockExists[coord.blockId].first) {
        coord = navigator.newdecrement(coord, data.blockStrand);
    }

    // Perform DFS and score each node
    std::function<void(TreeNode*, bm::bvector<>&, bm::bvector<>&)> dfs = [&](TreeNode* node, bm::bvector<>& seedVec, bm::bvector<>& parentSeedVec) {
        // Construct seedVec and gapMap at this node using the index
        bm::bvector<> delta = compressedSeedDeltas[node->id];
        seedVec.bit_xor(delta);

        // Perform scoring (not implemented here)

        for (TreeNode* child : node->children) {
            dfs(child, seedVec, delta);
        }

        // Backtrack
        seedVec.bit_xor(delta);
    };

    bm::bvector<> seedVec;
    bm::bvector<> parentSeedVec;
    dfs(T->root, seedVec, parentSeedVec);
}

int main() {
    
    // Load the index from file
    bm::bvector<> allBitsets;
    loadIndexFromFile("compressedSeedDeltas.bin", allBitsets);

    // Load gap maps from file
    std::vector<std::map<int64_t, int64_t>> gapMaps;
    readGapMapFromFile(gapMaps, "gapMaps.bin");

    // Traverse and score the tree
    traverseAndScore(&T, allBitsets, gapMaps, 5, 1000); // Adjust numBitsets and bitsetSize as needed

    return 0;
}
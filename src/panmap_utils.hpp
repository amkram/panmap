#pragma once

#include "panmanUtils.hpp"
#include "capnp/message.h"
#include "capnp/serialize-packed.h"
#include "index_lite.capnp.h"
#include "logging.hpp"
#include "gap_map_utils.hpp"
#include <string>
#include <iostream>
#include <span>
#include <algorithm>
#include <cstdint>
#include <vector>
#include <map>
#include <array>
#include <unordered_set>
#include <tuple>
#include <climits>
#include <random>

namespace panmapUtils {

enum seedChangeType { ADD, DEL, SUB };

void getSequenceFromReference(panmanUtils::Tree* tree,
                              std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
                              std::vector<char>& blockExists,
                              std::vector<char>& blockStrand,
                              std::unordered_map<int, int>& blockLengths,
                              std::string reference);

std::string getStringFromSequence(const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
                                  const std::unordered_map<int, int>& blockLengths,
                                  const std::vector<char>& blockExists,
                                  const std::vector<char>& blockStrand,
                                  bool aligned);

std::string getStringFromReference(panmanUtils::Tree* tree, std::string reference, bool aligned);

void simulateSNPsOnSequence(
  std::string& sequence,
  std::vector<std::tuple<char,char, uint32_t>>& snpRecords,
  uint32_t numsnps,
  std::mt19937& rng
);

class LiteNode {
   public:
    std::string identifier;
    LiteNode* parent = nullptr;
    std::vector<LiteNode*> children;

    // Node index in DFS order
    uint32_t nodeIndex = 0;

    // Seed changes from parent: range into the tree's zero-copy seed-change SoA.
    uint64_t seedChangeOffset = 0;
    uint32_t seedChangeSize = 0;
};

class LiteTree {
   public:
    LiteNode* root = nullptr;
    std::unordered_map<std::string, LiteNode*> allLiteNodes;
    std::vector<std::pair<uint32_t, uint32_t>> blockScalarRanges;

    // dfsIndex -> LiteNode*
    std::vector<LiteNode*> dfsIndexToNode;

    // Zero-copy seed-change SoA: raw pointers into the mapped decompressed index, one
    // per capnp segment (counts are Int16 on disk). Nodes hold (offset,size) ranges;
    // SEED_CHANGE_SEGMENT must match the index writer's CAPNP_SPLIT.
    static constexpr uint64_t SEED_CHANGE_SEGMENT = 500'000'000ULL;
    std::vector<const uint64_t*> segSeedHash;
    std::vector<const int16_t*> segSeedParent;
    std::vector<const int16_t*> segSeedChild;

    // Iterate a node's seed changes, splitting at segment boundaries. Non-hot callers.
    template <typename F>
    void forEachSeedChange(uint64_t offset, uint64_t size, F&& f) const {
        uint64_t pos = offset, remaining = size;
        while (remaining > 0) {
            const uint64_t seg = pos / SEED_CHANGE_SEGMENT;
            const uint64_t local = pos - seg * SEED_CHANGE_SEGMENT;
            const uint64_t n = std::min<uint64_t>(remaining, SEED_CHANGE_SEGMENT - local);
            const uint64_t* H = segSeedHash[seg] + local;
            const int16_t* P = segSeedParent[seg] + local;
            const int16_t* C = segSeedChild[seg] + local;
            for (uint64_t i = 0; i < n; ++i) f(H[i], static_cast<int64_t>(P[i]), static_cast<int64_t>(C[i]));
            pos += n;
            remaining -= n;
        }
    }

    // Skip re-loading seed changes in batch mode
    bool seedChangesLoaded = false;

    ~LiteTree() {
        for (auto& pair : allLiteNodes) {
            delete pair.second;
        }
    }

    void initialize(::LiteTree::Reader liteTreeReader);

    uint32_t getBlockStartScalar(const uint32_t blockId) const;
    uint32_t getBlockEndScalar(const uint32_t blockId) const;

    std::string resolveNodeId(uint32_t nodeIndex) const {
        if (nodeIndex < dfsIndexToNode.size() && dfsIndexToNode[nodeIndex]) {
            return dfsIndexToNode[nodeIndex]->identifier;
        }
        return "";
    }
};

struct Coordinate {
    int32_t nucPosition;
    int32_t nucGapPosition;
    int32_t primaryBlockId;
    // (secondaryBlockId removed: an old PanMAN feature no longer used; always -1)

    Coordinate() {}

    Coordinate(int nucPosition, int nucGapPosition, int primaryBlockId) {
        this->nucPosition = nucPosition;
        this->nucGapPosition = nucGapPosition;
        this->primaryBlockId = primaryBlockId;
    }

    Coordinate(const panmanUtils::NucMut& nm, int offset) {
        nucPosition = nm.nucPosition;
        nucGapPosition = nm.nucGapPosition;
        primaryBlockId = nm.primaryBlockId;
        moveForward(offset);
    }

    Coordinate(const panmanUtils::NucMut& nm) : Coordinate(nm, 0) {}

    friend std::ostream& operator<<(std::ostream& os, const Coordinate& coord) {
        os << "(" << coord.primaryBlockId << ", " << coord.nucPosition << ", " << coord.nucGapPosition << ")";
        return os;
    }

    void setSequenceBase(std::vector<std::vector<std::pair<char, std::vector<char>>>>& seq, char newNuc) const {
        if (nucGapPosition != -1) {
            seq[primaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
        } else {
            seq[primaryBlockId][nucPosition].first = newNuc;
        }
    }

    void moveForward(int offset) {
        if (nucGapPosition == -1) {
            nucPosition += offset;
        } else {
            nucGapPosition += offset;
        }
    }

    bool operator==(const Coordinate& other) const {
        return nucPosition == other.nucPosition && nucGapPosition == other.nucGapPosition &&
               primaryBlockId == other.primaryBlockId;
    }

    bool operator<(const Coordinate& other) const {
        if (primaryBlockId != other.primaryBlockId) return primaryBlockId < other.primaryBlockId;
        if (nucPosition != other.nucPosition) return nucPosition < other.nucPosition;

        if (nucGapPosition != other.nucGapPosition) {
            auto adjustGap = [](int32_t gap) {
                return gap == -1 ? INT32_MAX : gap;
            };
            return adjustGap(nucGapPosition) < adjustGap(other.nucGapPosition);
        }

        return false;
    }

    bool operator!=(const Coordinate& other) const { return !(*this == other); }

    bool operator<=(const Coordinate& other) const { return *this < other || *this == other; }

    bool operator>(const Coordinate& other) const { return !(*this <= other); }

    bool operator>=(const Coordinate& other) const { return !(*this < other); }
};

struct NewSyncmerRange {
    Coordinate begCoord;
    Coordinate endCoord;
    std::string localRangeSeq;
    std::vector<uint32_t> localRangeCoordToGlobalScalarCoords;
    std::vector<uint32_t> localRangeCoordToBlockId;
    std::vector<uint64_t> seedsToDelete;
};

// Decode a block's 8-nibble-packed consensusSeq, invoking emit(nucCode) for each nucleotide
// (codes 1..15) until the first zero nibble, which terminates the block.
template <class Seq, class F>
inline void forEachConsensusNuc(const Seq& consensusSeq, F&& emit) {
    for (size_t i = 0; i < consensusSeq.size(); i++) {
        for (size_t j = 0; j < 8; j++) {
            const int nucCode = ((consensusSeq[i] >> (4 * (7 - j))) & 15);
            if (nucCode == 0) return;
            emit(nucCode);
        }
    }
}

struct BlockSequences {
    // Flat layout (equivalent to the old vector<vector<pair<char, vector<char>>>> sequence):
    //   mainSeq[b][p]   consensus nucleotide at position p of block b
    //   gapSeq[b]       block b's gap ("insertion") chars, concatenated by position (CSR values)
    //   gapStart[b][p]  offset of position p's gaps in gapSeq[b]; size mainSeq[b].size()+1 (CSR prefix sums)
    std::vector<std::vector<char>> mainSeq;
    std::vector<std::vector<char>> gapSeq;
    std::vector<std::vector<uint32_t>> gapStart;
    std::vector<char> blockExists;
    std::vector<char> blockStrand;
    uint64_t totalLength;

    BlockSequences() {}

    BlockSequences(panmanUtils::Tree* tree) {
        mainSeq.resize(tree->blocks.size() + 1);
        blockExists.resize(tree->blocks.size() + 1, false);
        blockStrand.resize(tree->blocks.size() + 1, true);
        totalLength = 0;

        int32_t maxBlockId = 0;

        for (size_t i = 0; i < tree->blocks.size(); i++) {
            const auto& curBlock = tree->blocks[i];
            int32_t primaryBlockId = ((int32_t)curBlock.primaryBlockId);
            if (i != static_cast<size_t>(primaryBlockId)) {
                output::error("Block ID mismatch: primaryBlockId={} i={}", primaryBlockId, i);
                std::exit(1);
            }
            maxBlockId = std::max(maxBlockId, primaryBlockId);
            forEachConsensusNuc(curBlock.consensusSeq, [&](int nucCode) {
                mainSeq[primaryBlockId].push_back(panmanUtils::getNucleotideFromCode(nucCode));
                totalLength++;
            });
            // end sentinel to anchor trailing gaps
            mainSeq[primaryBlockId].push_back('x');
        }

        mainSeq.resize(maxBlockId + 1);
        blockExists.resize(maxBlockId + 1);
        blockStrand.resize(maxBlockId + 1);

        // Build the per-position gap CSR: counts -> prefix sums -> '-'-filled gap chars.
        gapStart.resize(mainSeq.size());
        gapSeq.resize(mainSeq.size());
        for (size_t b = 0; b < mainSeq.size(); b++) {
            gapStart[b].assign(mainSeq[b].size() + 1, 0);
        }
        for (size_t i = 0; i < tree->gaps.size(); i++) {
            const auto& curGap = tree->gaps[i];
            int32_t primaryBId = (curGap.primaryBlockId);
            for (size_t j = 0; j < curGap.nucPosition.size(); j++) {
                int len = curGap.nucGapLength[j];
                int pos = curGap.nucPosition[j];
                gapStart[primaryBId][pos + 1] = len;   // count at pos (last wins, as the old resize did)
                totalLength += len;
            }
        }
        for (size_t b = 0; b < mainSeq.size(); b++) {
            for (size_t p = 1; p < gapStart[b].size(); p++) {
                gapStart[b][p] += gapStart[b][p - 1];
            }
            gapSeq[b].assign(gapStart[b].back(), '-');
        }
    }

    int numBlocks() const { return mainSeq.size(); }

    size_t blockLength(int blockId) const { return mainSeq[blockId].size(); }

    uint32_t gapLength(int blockId, int pos) const { return gapStart[blockId][pos + 1] - gapStart[blockId][pos]; }

    char mainBase(int blockId, int pos) const { return mainSeq[blockId][pos]; }

    char gapBase(int blockId, int pos, int gapPos) const {
        return gapSeq[blockId][gapStart[blockId][pos] + gapPos];
    }

    const bool getBlockStrand(int blockId) const { return blockStrand[blockId]; }

    const bool getBlockExists(int blockId) const { return blockExists[blockId]; }

    char getSequenceBase(const Coordinate& coord) const {
        if (coord.nucGapPosition != -1) {
            return gapSeq[coord.primaryBlockId][gapStart[coord.primaryBlockId][coord.nucPosition] +
                                                coord.nucGapPosition];
        } else {
            return mainSeq[coord.primaryBlockId][coord.nucPosition];
        }
    }

    void setSequenceBase(const Coordinate& coord, char newNuc) {
        if (coord.nucGapPosition != -1) {
            gapSeq[coord.primaryBlockId][gapStart[coord.primaryBlockId][coord.nucPosition] + coord.nucGapPosition] =
                newNuc;
        } else {
            mainSeq[coord.primaryBlockId][coord.nucPosition] = newNuc;
        }
    }
};

struct BlockEdgeCoord {
    Coordinate start;
    Coordinate end;
    uint64_t startScalar;
    uint64_t endScalar;
};

struct GlobalCoords {
    // Flat CSR layout (was vector<vector<pair<int64_t, vector<int64_t>>>> globalCoords):
    //   mainScalar[b][p]  main scalar coord (-1 for the 'x' sentinel)
    //   gapScalar[b]      gap scalar coords concatenated by position (CSR values)
    //   gapStart[b][p]    CSR offsets into gapScalar[b]; size mainScalar[b].size()+1
    std::vector<std::vector<int64_t>> mainScalar;
    std::vector<std::vector<int64_t>> gapScalar;
    std::vector<std::vector<uint32_t>> gapStart;
    std::vector<Coordinate> scalarToCoord;
    std::vector<BlockEdgeCoord> blockEdgeCoords;
    std::tuple<int64_t, int64_t, int64_t> firstTupleCoord;
    std::tuple<int64_t, int64_t, int64_t> lastTupleCoord;
    Coordinate firstCoord;
    Coordinate lastCoord;
    int64_t lastScalarCoord;

    size_t numBlocks() const { return mainScalar.size(); }
    size_t nPos(int64_t b) const { return mainScalar[b].size(); }
    uint32_t nGap(int64_t b, int64_t p) const { return gapStart[b][p + 1] - gapStart[b][p]; }
    int64_t mainScalarAt(int64_t b, int64_t p) const { return mainScalar[b][p]; }
    int64_t gapScalarAt(int64_t b, int64_t p, int64_t k) const { return gapScalar[b][gapStart[b][p] + k]; }

    GlobalCoords(const BlockSequences& blockSequences) {
        scalarToCoord.reserve(blockSequences.totalLength);
        int64_t curScalarCoord = 0;
        int nb = blockSequences.numBlocks();
        gapStart = blockSequences.gapStart;   // gap positions align 1:1 with BlockSequences
        mainScalar.resize(nb);
        gapScalar.resize(nb);
        for (int i = 0; i < nb; i++) {
            size_t blen = blockSequences.blockLength(i);
            mainScalar[i].resize(blen);
            gapScalar[i].resize(gapStart[i].back());
            for (size_t j = 0; j < blen; j++) {
                // process gap nucs first
                uint32_t gaps = blockSequences.gapLength(i, j);
                for (uint32_t k = 0; k < gaps; k++) {
                    gapScalar[i][gapStart[i][j] + k] = curScalarCoord;
                    scalarToCoord.emplace_back(j, k, i);
                    curScalarCoord++;
                }

                // process main nuc
                if (blockSequences.mainBase(i, j) == 'x') {
                    mainScalar[i][j] = -1;
                } else {
                    mainScalar[i][j] = curScalarCoord;
                    scalarToCoord.emplace_back(j, -1, i);
                    curScalarCoord++;
                }
            }
        }

        lastScalarCoord = curScalarCoord - 1;

        size_t lastBlen = blockSequences.blockLength(nb - 1);
        if (blockSequences.gapLength(nb - 1, lastBlen - 1) == 0) {
            // last block's last nuc gap is empty, take the second to last main nuc
            lastTupleCoord = std::make_tuple(nb - 1, lastBlen - 2, -1);
            lastCoord = Coordinate(lastBlen - 2, -1, nb - 1);
        } else {
            // last block's last nuc gap is not empty, take the last nuc gap
            uint32_t lastGaps = blockSequences.gapLength(nb - 1, lastBlen - 1);
            lastTupleCoord = std::make_tuple(nb - 1, lastBlen - 1, lastGaps - 1);
            lastCoord = Coordinate(lastBlen - 1, lastGaps - 1, nb - 1);
        }

        if (blockSequences.gapLength(0, 0) == 0) {
            // first block's first nuc gap is empty, take the first main nuc
            firstTupleCoord = std::make_tuple(0, 0, -1);
            firstCoord = Coordinate(0, -1, 0);
        } else {
            firstTupleCoord = std::make_tuple(0, 0, 0);
            firstCoord = Coordinate(0, 0, 0);
        }

        blockEdgeCoords.resize(nb);
        for (int i = 0; i < nb; i++) {
            blockEdgeCoords[i].start = getBlockStartCoord(i);
            blockEdgeCoords[i].end = getBlockEndCoord(i);
            blockEdgeCoords[i].startScalar = getBlockStartScalar(i);
            blockEdgeCoords[i].endScalar = getBlockEndScalar(i);
        }

        if (getScalarFromTuple(firstTupleCoord) != 0) {
            logging::err("firstScalarCoord != 0");
            std::exit(1);
        }
        if (getScalarFromTuple(lastTupleCoord) != lastScalarCoord) {
            logging::err("lastScalarCoord != getScalarFromTuple(lastTupleCoord)");
            std::exit(1);
        }
    }

    std::tuple<int64_t, int64_t, int64_t> getBlockStartTuple(int64_t blockId) const {
        if (nGap(blockId, 0) == 0) {
            return std::make_tuple(blockId, 0, -1);
        }
        return std::make_tuple(blockId, 0, 0);
    }

    Coordinate getBlockStartCoord(int64_t blockId) const {
        int32_t nucPosition = 0;
        int32_t nucGapPosition = (nGap(blockId, 0) == 0) ? -1 : 0;
        return Coordinate(nucPosition, nucGapPosition, blockId);
    }

    std::tuple<int64_t, int64_t, int64_t> getBlockEndTuple(int64_t blockId) const {
        int64_t last = nPos(blockId) - 1;
        if (nGap(blockId, last) == 0) {
            return std::make_tuple(blockId, last - 1, -1);
        }
        return std::make_tuple(blockId, last, nGap(blockId, last) - 1);
    }

    Coordinate getBlockEndCoord(int64_t blockId) const {
        int64_t last = nPos(blockId) - 1;
        int32_t nucPosition;
        int32_t nucGapPosition;
        if (nGap(blockId, last) == 0) {
            nucPosition = last - 1;
            nucGapPosition = -1;
        } else {
            nucPosition = last;
            nucGapPosition = nGap(blockId, last) - 1;
        }
        return Coordinate(nucPosition, nucGapPosition, blockId);
    }

    int64_t getBlockStartScalar(int64_t blockId) const {
        auto [primaryBlockId, nucPos, nucGapPos] = getBlockStartTuple(blockId);
        if (nucGapPos == -1) {
            return mainScalarAt(primaryBlockId, nucPos);
        } else {
            return gapScalarAt(primaryBlockId, nucPos, nucGapPos);
        }
    }

    int64_t getBlockEndScalar(int64_t blockId) const {
        auto [primaryBlockId, nucPos, nucGapPos] = getBlockEndTuple(blockId);
        if (nucGapPos == -1) {
            return mainScalarAt(primaryBlockId, nucPos);
        } else {
            return gapScalarAt(primaryBlockId, nucPos, nucGapPos);
        }
    }

    int64_t getScalarFromTuple(const std::tuple<int64_t, int64_t, int64_t>& tupleCoord, bool blockStrand = true) const {
        const auto& [blockId, nucPos, nucGapPos] = tupleCoord;
        if (nucGapPos == -1) {
            if (!blockStrand) {
                return getBlockStartScalar(blockId) + getBlockEndScalar(blockId) - mainScalarAt(blockId, nucPos);
            }
            return mainScalarAt(blockId, nucPos);
        }

        if (!blockStrand) {
            return getBlockStartScalar(blockId) + getBlockEndScalar(blockId) -
                   gapScalarAt(blockId, nucPos, nucGapPos);
        }
        return gapScalarAt(blockId, nucPos, nucGapPos);
    }

    int64_t getScalarFromCoord(const Coordinate& coord, bool blockStrand = true) const {
        if (coord.nucGapPosition == -1) {
            const int64_t mainSc = mainScalarAt(coord.primaryBlockId, coord.nucPosition);
            return blockStrand ? mainSc
                               : getBlockStartScalar(coord.primaryBlockId) + getBlockEndScalar(coord.primaryBlockId) -
                                     mainSc;
        }

        const int64_t gapPosScalar = gapScalarAt(coord.primaryBlockId, coord.nucPosition, coord.nucGapPosition);
        return blockStrand
                   ? gapPosScalar
                   : getBlockStartScalar(coord.primaryBlockId) + getBlockEndScalar(coord.primaryBlockId) - gapPosScalar;
    }

    Coordinate getCoordFromScalar(uint64_t scalar, bool blockStrand = true) const {
        if (scalar >= scalarToCoord.size()) {
            output::error("In getCoordFromScalar(), scalar {} is out of range", scalar);
            exit(1);
        }
        Coordinate forwardCoord = scalarToCoord[scalar];
        if (blockStrand) {
            return forwardCoord;
        }
        return scalarToCoord[getBlockStartScalar(forwardCoord.primaryBlockId) +
                             getBlockEndScalar(forwardCoord.primaryBlockId) - scalar];
    }

    uint32_t getBlockIdFromScalar(uint64_t scalar) const { return scalarToCoord[scalar].primaryBlockId; }

    void stepRightCoordinate(std::tuple<int64_t, int64_t, int64_t>& coord) const {
        const auto& [blockId, nucPos, nucGapPos] = coord;
        if (coord == lastTupleCoord) {
            coord = std::make_tuple(-1, -1, -1);
            return;
        }

        if (coord == getBlockEndTuple(blockId)) {
            coord = getBlockStartTuple(blockId + 1);
            return;
        }

        if (nucGapPos == -1) {
            // cur coord is a main nuc, step right to the next main nuc
            if (nGap(blockId, nucPos + 1) == 0) {
                // next main nuc has no gap nucs -> return the main nuc
                coord = std::make_tuple(blockId, nucPos + 1, -1);
            } else {
                coord = std::make_tuple(blockId, nucPos + 1, 0);
            }

        } else {
            // cur coord is at a gap nuc
            if (nucGapPos == nGap(blockId, nucPos) - 1) {
                // cur gap nuc is the last gap nuc -> return the main nuc
                coord = std::make_tuple(blockId, nucPos, -1);
            } else {
                coord = std::make_tuple(blockId, nucPos, nucGapPos + 1);
            }
        }
    }

    void stepLeftCoordinate(std::tuple<int64_t, int64_t, int64_t>& coord) const {
        const auto& [blockId, nucPos, nucGapPos] = coord;
        if (coord == firstTupleCoord) {
            coord = std::make_tuple(-1, -1, -1);
            return;
        }

        if (coord == getBlockStartTuple(blockId)) {
            coord = getBlockEndTuple(blockId - 1);
            return;
        }

        if (nucGapPos == -1) {
            // cur coord is a main nuc
            if (nGap(blockId, nucPos) == 0) {
                // cur coord has no gap nucs, return the previous main nuc
                coord = std::make_tuple(blockId, nucPos - 1, -1);
            } else {
                // cur coord has gap nucs, return the last gap nuc
                coord = std::make_tuple(blockId, nucPos, nGap(blockId, nucPos) - 1);
            }
        } else {
            if (nucGapPos == 0) {
                coord = std::make_tuple(blockId, nucPos - 1, -1);
            } else {
                coord = std::make_tuple(blockId, nucPos, nucGapPos - 1);
            }
        }
    }

    void stepRightCoordinate(Coordinate& coord) const {
        if (coord == lastCoord) {
            coord = Coordinate(-1, -1, -1);
            return;
        }

        if (coord == blockEdgeCoords[coord.primaryBlockId].end) {
            coord = blockEdgeCoords[coord.primaryBlockId + 1].start;
            return;
        }

        if (coord.nucGapPosition == -1) {
            coord.nucPosition++;
            if (nGap(coord.primaryBlockId, coord.nucPosition) != 0) {
                // next main nuc has gap nucs -> step to the first gap nuc
                coord.nucGapPosition = 0;
            }
        } else {
            // cur coord is at a gap nuc
            if (coord.nucGapPosition == nGap(coord.primaryBlockId, coord.nucPosition) - 1) {
                // cur gap nuc is the last gap nuc -> step to the main nuc
                coord.nucGapPosition = -1;
            } else {
                coord.nucGapPosition++;
            }
        }
    }

    void stepLeftCoordinate(Coordinate& coord) const {
        if (coord == firstCoord) {
            coord = Coordinate(-1, -1, -1);
            return;
        }

        if (coord == blockEdgeCoords[coord.primaryBlockId].start) {
            coord = blockEdgeCoords[coord.primaryBlockId - 1].end;
            return;
        }

        if (coord.nucGapPosition == -1) {
            // cur coord is a main nuc
            if (nGap(coord.primaryBlockId, coord.nucPosition) == 0) {
                // cur coord has no gap nucs, return the previous main nuc
                coord.nucPosition--;
            } else {
                // cur coord has gap nucs, return the last gap nuc
                coord.nucGapPosition = nGap(coord.primaryBlockId, coord.nucPosition) - 1;
            }
        } else {
            if (coord.nucGapPosition == 0) {
                coord.nucPosition--;
                coord.nucGapPosition = -1;
            } else {
                coord.nucGapPosition--;
            }
        }
    }

    void stepForwardScalar(Coordinate& coord, const std::vector<char>& blockStrand) const {
        const int64_t originalBlockId = coord.primaryBlockId;
        if (blockStrand[originalBlockId]) {
            stepRightCoordinate(coord);
            if (coord.primaryBlockId == originalBlockId) {
                return;
            } else if (coord.primaryBlockId == originalBlockId + 1) {
                if (!blockStrand[coord.primaryBlockId]) {
                    coord = blockEdgeCoords[coord.primaryBlockId].end;
                }
                return;
            } else {
                output::error(
                    "Stepping right over more than one block: {} to {}", originalBlockId, coord.primaryBlockId);
                std::exit(1);
            }
        } else {
            stepLeftCoordinate(coord);
            if (coord.primaryBlockId == originalBlockId) {
                return;
            } else if (coord.primaryBlockId == originalBlockId - 1) {
                if (blockStrand[originalBlockId + 1]) {
                    coord = blockEdgeCoords[originalBlockId + 1].start;
                } else {
                    coord = blockEdgeCoords[originalBlockId + 1].end;
                }
                return;
            } else {
                output::error(
                    "Stepping left over more than one block: {} to {}", originalBlockId, coord.primaryBlockId);
                std::exit(1);
            }
        }
    }

    Coordinate stepForwardScalar(const Coordinate& coord, const std::vector<char>& blockStrand) const {
        Coordinate nextCoord = coord;
        stepForwardScalar(nextCoord, blockStrand);
        return nextCoord;
    }

    void stepBackwardScalar(Coordinate& coord, const std::vector<char>& blockStrand) const {
        const int64_t originalBlockId = coord.primaryBlockId;
        if (blockStrand[originalBlockId]) {
            stepLeftCoordinate(coord);
            if (coord.primaryBlockId == originalBlockId) {
                return;
            } else if (coord.primaryBlockId == originalBlockId - 1) {
                if (!blockStrand[coord.primaryBlockId]) {
                    coord = blockEdgeCoords[coord.primaryBlockId].start;
                }
            } else {
                output::error(
                    "Stepping left over more than one block: {} to {}", originalBlockId, coord.primaryBlockId);
                std::exit(1);
            }
        } else {
            stepRightCoordinate(coord);
            if (coord.primaryBlockId == originalBlockId) {
                return;
            } else if (coord.primaryBlockId == originalBlockId + 1) {
                if (blockStrand[originalBlockId - 1]) {
                    coord = blockEdgeCoords[originalBlockId - 1].end;
                } else {
                    coord = blockEdgeCoords[originalBlockId - 1].start;
                }
            } else {
                output::error(
                    "Stepping right over more than one block: {} to {}", originalBlockId, coord.primaryBlockId);
                std::exit(1);
            }
        }
    }

    Coordinate stepBackwardScalar(const Coordinate& coord, const std::vector<char>& blockStrand) const {
        Coordinate nextCoord = coord;
        stepBackwardScalar(nextCoord, blockStrand);
        return nextCoord;
    }
};

inline bool isCanonical(char nuc) {
    return (nuc == 'A' || nuc == 'T' || nuc == 'C' || nuc == 'G');
}

// True if oldNuc is canonical and newNuc is an ambiguous IUPAC code (not gap '-' or sentinel 'x').
inline bool canonicalToAmb(char oldNuc, char newNuc) {
    return (newNuc != '-' && newNuc != 'x' && isCanonical(oldNuc) && !isCanonical(newNuc));
}

// Shared between lite and MGSR index building

inline void applyMutations(panmanUtils::Node* node,
                           size_t dfsIndex,
                           BlockSequences& blockSequences,
                           std::unordered_set<uint64_t>& invertedBlocks,
                           GlobalCoords& globalCoords,
                           std::vector<std::pair<Coordinate, Coordinate>>& localMutationRanges,
                           std::vector<std::tuple<uint32_t, bool, bool, bool, bool>>& blockMutationRecord,
                           std::vector<std::tuple<Coordinate, char, char>>& nucMutationRecord,
                           std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapRunUpdates,
                           std::vector<std::pair<uint64_t, bool>>& invertedBlocksBacktracks,
                           std::vector<uint32_t>& potentialSyncmerDeletions,
                           const std::vector<char>& oldBlockExists,
                           const std::vector<char>& oldBlockStrand,
                           bool imputeAmb = false) {
    std::vector<char>& blockExists = blockSequences.blockExists;
    std::vector<char>& blockStrand = blockSequences.blockStrand;

    for (const auto& blockMutation : node->blockMutation) {
        const int32_t blockId = blockMutation.primaryBlockId;
        const bool isInsertion = blockMutation.blockMutInfo;
        const bool isInversion = blockMutation.inversion;
        const bool oldExists = blockExists[blockId];
        const bool oldStrand = blockStrand[blockId];

        if (isInsertion) {
            blockExists[blockId] = true;
            blockStrand[blockId] = !isInversion;
            if (!blockStrand[blockId]) {
                invertedBlocks.insert(blockId);
                invertedBlocksBacktracks.emplace_back(blockId, true);
            }
        } else if (isInversion) {
            blockStrand[blockId] = !blockStrand[blockId];
            if (!blockStrand[blockId]) {
                invertedBlocks.insert(blockId);
                invertedBlocksBacktracks.emplace_back(blockId, true);
            } else {
                invertedBlocks.erase(blockId);
                invertedBlocksBacktracks.emplace_back(blockId, false);
            }
        } else {
            blockExists[blockId] = false;
            blockStrand[blockId] = true;
            if (!oldStrand) {
                invertedBlocks.erase(blockId);
                invertedBlocksBacktracks.emplace_back(blockId, false);
            }
        }
        blockMutationRecord.emplace_back(blockId, oldExists, oldStrand, blockExists[blockId], blockStrand[blockId]);

        const auto& curBlockEdgeCoords = globalCoords.blockEdgeCoords[blockId];
        if (blockStrand[blockId]) {
            localMutationRanges.emplace_back(curBlockEdgeCoords.start, curBlockEdgeCoords.end);
        } else {
            localMutationRanges.emplace_back(curBlockEdgeCoords.end, curBlockEdgeCoords.start);
        }
    }

    for (const auto& nucMutation : node->nucMutation) {
        int length = nucMutation.mutInfo >> 4;
        int blockId;
        int lastOffset = -1;

        for (int i = 0; i < length; i++) {
            Coordinate pos = Coordinate(nucMutation, i);
            if ((pos.nucPosition == blockSequences.blockLength(pos.primaryBlockId) - 1 && pos.nucGapPosition == -1) ||
                (pos.nucPosition >= blockSequences.blockLength(pos.primaryBlockId))) {
                continue;
            }
            lastOffset = i;
            blockId = pos.primaryBlockId;
            const char oldNuc = blockSequences.getSequenceBase(pos);
            const int newNucCode = (nucMutation.nucs >> (4 * (5 - i))) & 0xF;
            const char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);

            if (oldNuc == newNuc) continue;

            // Skip canonical->ambiguous mutations so syncmers inherit the parent's base.
            if (imputeAmb && canonicalToAmb(oldNuc, newNuc)) continue;

            blockSequences.setSequenceBase(pos, newNuc);
            nucMutationRecord.emplace_back(pos, oldNuc, newNuc);

            if (oldBlockExists[pos.primaryBlockId] && blockExists[pos.primaryBlockId]) {
                const int64_t scalarCoord = globalCoords.getScalarFromCoord(pos);
                if (newNuc == '-') {
                    if (!gapRunUpdates.empty() && gapRunUpdates.back().first == true &&
                        gapRunUpdates.back().second.second + 1 == scalarCoord) {
                        ++(gapRunUpdates.back().second.second);
                    } else {
                        gapRunUpdates.emplace_back(true, std::make_pair(scalarCoord, scalarCoord));
                    }
                    if (blockExists[blockId] && oldBlockExists[blockId] &&
                        blockStrand[blockId] == oldBlockStrand[blockId]) {
                        potentialSyncmerDeletions.push_back(
                            globalCoords.getScalarFromCoord(pos, blockStrand[pos.primaryBlockId]));
                    }
                } else if (oldNuc == '-') {
                    if (!gapRunUpdates.empty() && gapRunUpdates.back().first == false &&
                        gapRunUpdates.back().second.second + 1 == scalarCoord) {
                        ++(gapRunUpdates.back().second.second);
                    } else {
                        gapRunUpdates.emplace_back(false, std::make_pair(scalarCoord, scalarCoord));
                    }
                }
            }
        }
        if (lastOffset != -1 && blockExists[blockId] && oldBlockExists[blockId] &&
            blockStrand[blockId] == oldBlockStrand[blockId]) {
            if (blockStrand[blockId]) {
                localMutationRanges.emplace_back(Coordinate(nucMutation, 0), Coordinate(nucMutation, lastOffset));
            } else {
                localMutationRanges.emplace_back(Coordinate(nucMutation, lastOffset), Coordinate(nucMutation, 0));
            }
        }
    }

    for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
        if (oldExists && !newExists) {
            uint64_t beg = (uint64_t)globalCoords.getBlockStartScalar(blockId);
            uint64_t end = (uint64_t)globalCoords.getBlockEndScalar(blockId);
            gapRunUpdates.emplace_back(true, std::make_pair(beg, end));
        } else if (!oldExists && newExists) {
            Coordinate coord = globalCoords.blockEdgeCoords[blockId].start;
            Coordinate end = globalCoords.blockEdgeCoords[blockId].end;
            std::pair<int64_t, int64_t> curNucRange = {-1, -1};
            while (true) {
                char nuc = blockSequences.getSequenceBase(coord);
                nuc = nuc == 'x' ? '-' : nuc;
                int64_t scalar = globalCoords.getScalarFromCoord(coord);
                if (nuc != '-') {
                    if (curNucRange.first != -1 && curNucRange.second + 1 == scalar) {
                        ++curNucRange.second;
                    } else {
                        if (curNucRange.first != -1) {
                            gapRunUpdates.emplace_back(
                                false, std::make_pair((uint64_t)curNucRange.first, (uint64_t)curNucRange.second));
                        }
                        curNucRange = {scalar, scalar};
                    }
                }

                if (coord == end) break;
                globalCoords.stepRightCoordinate(coord);
            }
            if (curNucRange.first != -1) {
                gapRunUpdates.emplace_back(false,
                                           std::make_pair((uint64_t)curNucRange.first, (uint64_t)curNucRange.second));
            }
        }
    }
}

inline void updateGapMap(std::map<uint64_t, uint64_t>& gapMap,
                         const std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& updates,
                         std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& backtrack,
                         std::vector<std::pair<bool, std::pair<uint64_t, uint64_t>>>& gapMapUpdates) {
    for (const auto& update : updates) {
        gap_map::updateGapMapStep(
            gapMap, update.second.first, update.second.second, update.first, backtrack, gapMapUpdates, true);
    }
}

// Returns (firstNonGapScalar, lastNonGapScalar). flankSize non-gap bases are skipped at
// each end. Returns (UINT64_MAX, 0) if the genome is all gaps or too short.

inline std::pair<uint64_t, uint64_t>
computeExtentFromGapMap(const std::map<uint64_t, uint64_t>& gapMap, uint64_t lastScalarCoord, uint64_t flankSize = 0) {
    auto findPositionAfterNonGaps = [&](uint64_t start, uint64_t count) -> uint64_t {
        uint64_t nonGapCount = 0;
        uint64_t pos = start;

        while (pos <= lastScalarCoord && nonGapCount < count) {
            auto it = gapMap.upper_bound(pos);
            if (it != gapMap.begin()) {
                --it;
                if (it->second >= pos) {
                    pos = it->second + 1;
                    continue;
                }
            }
            nonGapCount++;
            if (nonGapCount >= count) {
                return pos;
            }
            pos++;
        }
        return UINT64_MAX;
    };

    auto findPositionBeforeNonGaps = [&](uint64_t start, uint64_t count) -> uint64_t {
        uint64_t nonGapCount = 0;
        int64_t pos = static_cast<int64_t>(start);

        while (pos >= 0 && nonGapCount < count) {
            auto it = gapMap.upper_bound(static_cast<uint64_t>(pos));
            if (it != gapMap.begin()) {
                --it;
                if (it->second >= static_cast<uint64_t>(pos)) {
                    pos = static_cast<int64_t>(it->first) - 1;
                    continue;
                }
            }
            nonGapCount++;
            if (nonGapCount >= count) {
                return static_cast<uint64_t>(pos);
            }
            pos--;
        }
        return 0;
    };

    uint64_t firstNonGap = 0;
    uint64_t lastNonGap = lastScalarCoord;

    auto it = gapMap.find(0);
    if (it != gapMap.end()) {
        if (it->second >= lastScalarCoord) {
            return {UINT64_MAX, 0};
        }
        firstNonGap = it->second + 1;
    }

    if (!gapMap.empty()) {
        auto lastIt = std::prev(gapMap.end());
        if (lastIt->second == lastScalarCoord) {
            lastNonGap = lastIt->first - 1;
        }
    }

    if (flankSize > 0) {
        uint64_t maskedFirst = findPositionAfterNonGaps(firstNonGap, flankSize);
        uint64_t maskedLast = findPositionBeforeNonGaps(lastNonGap, flankSize);

        if (maskedFirst == UINT64_MAX || maskedLast == 0 || maskedFirst > maskedLast) {
            return {UINT64_MAX, 0};
        }

        firstNonGap = maskedFirst;
        lastNonGap = maskedLast;
    }

    return {firstNonGap, lastNonGap};
}

}  // namespace panmapUtils

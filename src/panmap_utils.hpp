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
#include <map>
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

class LiteNode {
   public:
    std::string identifier;
    LiteNode* parent = nullptr;
    std::vector<LiteNode*> children;

    // Node index in DFS order (used for efficient lookups)
    uint32_t nodeIndex = 0;

    // Seed changes from parent to this node (hash, parentCount, childCount)
    std::span<const std::tuple<uint64_t, int64_t, int64_t>> seedChanges;

    // Placement scores (populated during placement)
    float logRawScore = 0.0f;
    float logCosineScore = 0.0f;
    float containmentScore = 0.0f;
    float weightedContainmentScore = 0.0f;
    float logContainmentScore = 0.0f;
};

class LiteTree {
   public:
    LiteNode* root = nullptr;
    std::unordered_map<std::string, LiteNode*> allLiteNodes;
    std::vector<std::pair<uint32_t, uint32_t>> blockScalarRanges;
    std::unordered_map<std::string, uint32_t> nodeToDfsIndex;

    // Index-based lookup (dfsIndex -> LiteNode*)
    std::vector<LiteNode*> dfsIndexToNode;

    // Storage for all seed changes (flat array, nodes reference spans into this)
    std::vector<std::tuple<uint64_t, int64_t, int64_t>> allSeedChanges;

    // True once seed changes have been loaded from the index (skip re-loading in batch mode)
    bool seedChangesLoaded = false;

    ~LiteTree() {
        for (auto& pair : allLiteNodes) {
            delete pair.second;
        }
    }

    void initialize(::LiteTree::Reader liteTreeReader);

    uint32_t getBlockStartScalar(const uint32_t blockId) const;
    uint32_t getBlockEndScalar(const uint32_t blockId) const;

    // Resolve node index to node ID string
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
    int32_t secondaryBlockId;

    Coordinate() {}

    // Create a Coordinate by position
    Coordinate(int nucPosition, int nucGapPosition, int primaryBlockId, int secondaryBlockId) {
        this->nucPosition = nucPosition;
        this->nucGapPosition = nucGapPosition;
        this->primaryBlockId = primaryBlockId;
        this->secondaryBlockId = secondaryBlockId;
    }

    // Create a Coordinate with an offset
    Coordinate(const panmanUtils::NucMut& nm, int offset) {
        nucPosition = nm.nucPosition;
        nucGapPosition = nm.nucGapPosition;
        primaryBlockId = nm.primaryBlockId;
        secondaryBlockId = nm.secondaryBlockId;
        moveForward(offset);
    }

    // Create a Coordinate copying a NucMut
    Coordinate(const panmanUtils::NucMut& nm) : Coordinate(nm, 0) {}

    friend std::ostream& operator<<(std::ostream& os, const Coordinate& coord) {
        os << "(" << coord.primaryBlockId << ", " << coord.nucPosition << ", " << coord.nucGapPosition << ")";
        return os;
    }

    // Get base corresponding to this Coordinate's position within a sequence_t
    char getSequenceBase(
        const std::vector<std::pair<std::vector<std::pair<char, std::vector<char>>>,
                                    std::vector<std::vector<std::pair<char, std::vector<char>>>>>>& seq) const {
        if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
                return seq[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            } else {
                return seq[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            }
        } else {
            if (nucGapPosition != -1) {
                return seq[primaryBlockId].first[nucPosition].second[nucGapPosition];
            } else {
                return seq[primaryBlockId].first[nucPosition].first;
            }
        }
    }

    char getSequenceBase(const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence) const {
        if (nucGapPosition != -1) {
            return sequence[primaryBlockId][nucPosition].second[nucGapPosition];
        } else {
            return sequence[primaryBlockId][nucPosition].first;
        }
    }

    char getSequenceBase(const std::vector<std::pair<char, std::vector<char>>>& blockSequence) const {
        if (nucGapPosition != -1) {
            return blockSequence[nucPosition].second[nucGapPosition];
        } else {
            return blockSequence[nucPosition].first;
        }
    }

    // Set base corresponding to this Coordinate's position within a sequence_t
    void setSequenceBase(std::vector<std::pair<std::vector<std::pair<char, std::vector<char>>>,
                                               std::vector<std::vector<std::pair<char, std::vector<char>>>>>>& seq,
                         char newNuc) const {
        if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
                seq[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId].second[secondaryBlockId][nucPosition].first = newNuc;
            }
        } else {
            if (nucGapPosition != -1) {
                seq[primaryBlockId].first[nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId].first[nucPosition].first = newNuc;
            }
        }
    }

    void setSequenceBase(std::vector<std::vector<std::pair<char, std::vector<char>>>>& seq, char newNuc) const {
        if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
                seq[primaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId][nucPosition].first = newNuc;
            }
        } else {
            if (nucGapPosition != -1) {
                seq[primaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId][nucPosition].first = newNuc;
            }
        }
    }

    void setSequenceBase(std::vector<std::pair<char, std::vector<char>>>& blockSequence, char newNuc) const {
        if (secondaryBlockId != -1) {
            if (nucGapPosition != -1) {
                blockSequence[nucPosition].second[nucGapPosition] = newNuc;
            } else {
                blockSequence[nucPosition].first = newNuc;
            }
        } else {
            if (nucGapPosition != -1) {
                blockSequence[nucPosition].second[nucGapPosition] = newNuc;
            } else {
                blockSequence[nucPosition].first = newNuc;
            }
        }
    }

    // Move "offset" steps forward
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
    std::vector<uint64_t> localRangeCoordToGlobalScalarCoords;
    std::vector<uint64_t> localRangeCoordToBlockId;
    std::vector<uint64_t> seedsToDelete;
};

struct BlockSequences {
    std::vector<std::vector<std::pair<char, std::vector<char>>>> sequence;
    std::vector<char> blockExists;
    std::vector<char> blockStrand;
    uint64_t totalLength;

    BlockSequences() {}

    BlockSequences(panmanUtils::Tree* tree) {
        sequence.resize(tree->blocks.size() + 1);
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
            for (size_t j = 0; j < curBlock.consensusSeq.size(); j++) {
                bool endFlag = false;
                for (size_t k = 0; k < 8; k++) {
                    const int nucCode = (((curBlock.consensusSeq[j]) >> (4 * (7 - k))) & 15);
                    if (nucCode == 0) {
                        endFlag = true;
                        break;
                    }
                    const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
                    sequence[primaryBlockId].push_back({nucleotide, {}});
                    totalLength++;
                }

                if (endFlag) {
                    break;
                }
            }
            // End character to incorporate for gaps at the end
            sequence[primaryBlockId].push_back({'x', {}});
        }

        sequence.resize(maxBlockId + 1);
        blockExists.resize(maxBlockId + 1);
        blockStrand.resize(maxBlockId + 1);

        for (size_t i = 0; i < tree->gaps.size(); i++) {
            const auto& curGap = tree->gaps[i];
            int32_t primaryBId = (curGap.primaryBlockId);
            int32_t secondaryBId = (curGap.secondaryBlockId);
            for (size_t j = 0; j < curGap.nucPosition.size(); j++) {
                int len = curGap.nucGapLength[j];
                int pos = curGap.nucPosition[j];
                sequence[primaryBId][pos].second.resize(len, '-');
                totalLength += len;
            }
        }
    }

    int numBlocks() const { return sequence.size(); }

    const bool getBlockStrand(int blockId) const { return blockStrand[blockId]; }

    const bool getBlockExists(int blockId) const { return blockExists[blockId]; }

    char getSequenceBase(const Coordinate& coord) const {
        if (coord.nucGapPosition != -1) {
            return sequence[coord.primaryBlockId][coord.nucPosition].second[coord.nucGapPosition];
        } else {
            return sequence[coord.primaryBlockId][coord.nucPosition].first;
        }
    }

    void setSequenceBase(const Coordinate& coord, char newNuc) {
        if (coord.nucGapPosition != -1) {
            sequence[coord.primaryBlockId][coord.nucPosition].second[coord.nucGapPosition] = newNuc;
        } else {
            sequence[coord.primaryBlockId][coord.nucPosition].first = newNuc;
        }
    }

    std::string getSequenceStringByBlockId(int blockId, bool aligned) const {
        std::string seq;
        const auto& curBlock = sequence[blockId];
        for (size_t i = 0; i < curBlock.size(); i++) {
            if (blockStrand[blockId]) {
                // Gap nucs
                for (size_t j = 0; j < curBlock[i].second.size(); j++) {
                    if (curBlock[i].second[j] != '-') {
                        seq += curBlock[i].second[j];
                    } else if (aligned) {
                        seq += '-';
                    }
                }

                // Main nuc
                if (curBlock[i].first != '-' && curBlock[i].first != 'x') {
                    seq += curBlock[i].first;
                } else if (aligned && curBlock[i].first != 'x') {
                    seq += '-';
                }
            } else {
                // Main nuc first
                if (curBlock[i].first != '-' && curBlock[i].first != 'x') {
                    seq += panmanUtils::getComplementCharacter(curBlock[i].first);
                } else if (aligned && curBlock[i].first != 'x') {
                    seq += '-';
                }

                // Gap nucs in reverse
                for (size_t j = curBlock[i].second.size() - 1; j + 1 > 0; j--) {
                    if (curBlock[i].second[j] != '-') {
                        seq += panmanUtils::getComplementCharacter(curBlock[i].second[j]);
                    } else if (aligned) {
                        seq += '-';
                    }
                }
            }
        }
        return seq;
    }

    std::string getSequenceString(bool aligned) const {
        std::string seq;
        for (size_t i = 0; i < numBlocks(); i++) {
            seq += getSequenceStringByBlockId(i, aligned);
        }
        return seq;
    }

    std::string getAlignedSequenceStringByBlockId(int blockId) const {
        return getSequenceStringByBlockId(blockId, true);
    }

    std::string getUnalignedSequenceStringByBlockId(int blockId) const {
        return getSequenceStringByBlockId(blockId, false);
    }

    std::string getAlignedSequenceString() const { return getSequenceString(true); }

    std::string getUnalignedSequenceString() const { return getSequenceString(false); }
};

struct BlockEdgeCoord {
    Coordinate start;
    Coordinate end;
    uint64_t startScalar;
    uint64_t endScalar;
};

struct GlobalCoords {
    std::vector<std::vector<std::pair<int64_t, std::vector<int64_t>>>> globalCoords;
    std::vector<Coordinate> scalarToCoord;
    std::vector<BlockEdgeCoord> blockEdgeCoords;
    std::tuple<int64_t, int64_t, int64_t> firstTupleCoord;
    std::tuple<int64_t, int64_t, int64_t> lastTupleCoord;
    Coordinate firstCoord;
    Coordinate lastCoord;
    int64_t lastScalarCoord;

    GlobalCoords(const BlockSequences& blockSequences) {
        const auto& sequence = blockSequences.sequence;
        scalarToCoord.reserve(blockSequences.totalLength);
        int64_t curScalarCoord = 0;
        globalCoords.resize(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++) {
            globalCoords[i].resize(sequence[i].size());
            for (size_t j = 0; j < sequence[i].size(); j++) {
                // process gap nucs first
                globalCoords[i][j].second.resize(sequence[i][j].second.size());
                for (size_t k = 0; k < sequence[i][j].second.size(); k++) {
                    globalCoords[i][j].second[k] = curScalarCoord;
                    scalarToCoord.emplace_back(j, k, i, -1);
                    curScalarCoord++;
                }

                // process main nuc
                if (sequence[i][j].first == 'x') {
                    // skip if main nuc is x
                    globalCoords[i][j].first = -1;
                } else {
                    globalCoords[i][j].first = curScalarCoord;
                    scalarToCoord.emplace_back(j, -1, i, -1);
                    curScalarCoord++;
                }
            }
        }

        lastScalarCoord = curScalarCoord - 1;

        const auto& lastBlock = sequence.back();
        if (lastBlock.back().second.empty()) {
            // last block's last nuc gap is empty, take the second to last main nuc
            lastTupleCoord = std::make_tuple(sequence.size() - 1, lastBlock.size() - 2, -1);
            lastCoord = Coordinate(lastBlock.size() - 2, -1, sequence.size() - 1, -1);
        } else {
            // last block's last nuc gap is not empty, take the last nuc gap
            lastTupleCoord =
                std::make_tuple(sequence.size() - 1, lastBlock.size() - 1, lastBlock.back().second.size() - 1);
            lastCoord = Coordinate(lastBlock.size() - 1, lastBlock.back().second.size() - 1, sequence.size() - 1, -1);
        }

        const auto& firstBlock = sequence[0];
        if (firstBlock[0].second.empty()) {
            // first block's first nuc gap is empty, take the first main nuc
            firstTupleCoord = std::make_tuple(0, 0, -1);
            firstCoord = Coordinate(0, -1, 0, -1);
        } else {
            // first block's first nuc gap is not empty, take the first gap nuc
            firstTupleCoord = std::make_tuple(0, 0, 0);
            firstCoord = Coordinate(0, 0, 0, -1);
        }

        blockEdgeCoords.resize(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++) {
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
        const auto& curBlock = globalCoords[blockId];
        if (curBlock[0].second.empty()) {
            return std::make_tuple(blockId, 0, -1);
        }
        return std::make_tuple(blockId, 0, 0);
    }

    Coordinate getBlockStartCoord(int64_t blockId) const {
        const auto& curBlock = globalCoords[blockId];
        int32_t nucPosition = 0;
        int32_t secondaryBlockId = -1;

        int32_t nucGapPosition;
        if (curBlock[0].second.empty()) {
            nucGapPosition = -1;
        } else {
            nucGapPosition = 0;
        }
        return Coordinate(nucPosition, nucGapPosition, blockId, secondaryBlockId);
    }

    std::tuple<int64_t, int64_t, int64_t> getBlockEndTuple(int64_t blockId) const {
        const auto& curBlock = globalCoords[blockId];
        if (curBlock.back().second.empty()) {
            return std::make_tuple(blockId, curBlock.size() - 2, -1);
        }
        return std::make_tuple(blockId, curBlock.size() - 1, curBlock.back().second.size() - 1);
    }

    Coordinate getBlockEndCoord(int64_t blockId) const {
        const auto& curBlock = globalCoords[blockId];
        int32_t nucPosition;
        int32_t nucGapPosition;
        int32_t secondaryBlockId = -1;
        if (curBlock.back().second.empty()) {
            nucPosition = curBlock.size() - 2;
            nucGapPosition = -1;
        } else {
            nucPosition = curBlock.size() - 1;
            nucGapPosition = curBlock.back().second.size() - 1;
        }
        return Coordinate(nucPosition, nucGapPosition, blockId, secondaryBlockId);
    }

    int64_t getBlockStartScalar(int64_t blockId) const {
        auto [primaryBlockId, nucPos, nucGapPos] = getBlockStartTuple(blockId);
        if (nucGapPos == -1) {
            return globalCoords[primaryBlockId][nucPos].first;
        } else {
            return globalCoords[primaryBlockId][nucPos].second[nucGapPos];
        }
    }

    int64_t getBlockEndScalar(int64_t blockId) const {
        auto [primaryBlockId, nucPos, nucGapPos] = getBlockEndTuple(blockId);
        if (nucGapPos == -1) {
            return globalCoords[primaryBlockId][nucPos].first;
        } else {
            return globalCoords[primaryBlockId][nucPos].second[nucGapPos];
        }
    }

    int64_t getScalarFromTuple(const std::tuple<int64_t, int64_t, int64_t>& tupleCoord, bool blockStrand = true) const {
        const auto& [blockId, nucPos, nucGapPos] = tupleCoord;
        if (nucGapPos == -1) {
            if (!blockStrand) {
                return getBlockStartScalar(blockId) + getBlockEndScalar(blockId) - globalCoords[blockId][nucPos].first;
            }
            return globalCoords[blockId][nucPos].first;
        }

        if (!blockStrand) {
            return getBlockStartScalar(blockId) + getBlockEndScalar(blockId) -
                   globalCoords[blockId][nucPos].second[nucGapPos];
        }
        return globalCoords[blockId][nucPos].second[nucGapPos];
    }

    int64_t getScalarFromCoord(const Coordinate& coord, bool blockStrand = true) const {
        const auto& nucPosInfo = globalCoords[coord.primaryBlockId][coord.nucPosition];
        if (coord.nucGapPosition == -1) {
            return blockStrand ? nucPosInfo.first
                               : getBlockStartScalar(coord.primaryBlockId) + getBlockEndScalar(coord.primaryBlockId) -
                                     nucPosInfo.first;
        }

        const int64_t gapPosScalar = nucPosInfo.second[coord.nucGapPosition];
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
            // go to the next block
            coord = getBlockStartTuple(blockId + 1);
            return;
        }

        // not end of block
        if (nucGapPos == -1) {
            // cur coord is a main nuc, step right to the next main nuc
            if (globalCoords[blockId][nucPos + 1].second.empty()) {
                // next main nuc has no gap nucs -> return the main nuc
                coord = std::make_tuple(blockId, nucPos + 1, -1);
            } else {
                // next main nuc has gap nucs -> return the first gap nuc
                coord = std::make_tuple(blockId, nucPos + 1, 0);
            }

        } else {
            // cur coord is at a gap nuc
            if (nucGapPos == globalCoords[blockId][nucPos].second.size() - 1) {
                // cur gap nuc is the last gap nuc -> return the main nuc
                coord = std::make_tuple(blockId, nucPos, -1);
            } else {
                // cur gap nuc is not the last gap nuc -> return the next gap nuc
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
            // go to the previous block
            coord = getBlockEndTuple(blockId - 1);
            return;
        }

        // not start of block
        if (nucGapPos == -1) {
            // cur coord is a main nuc
            if (globalCoords[blockId][nucPos].second.empty()) {
                // cur coord has no gap nucs, return the previous main nuc
                coord = std::make_tuple(blockId, nucPos - 1, -1);
            } else {
                // cur coord has gap nucs, return the last gap nuc
                coord = std::make_tuple(blockId, nucPos, globalCoords[blockId][nucPos].second.size() - 1);
            }
        } else {
            // cur coord is at a gap nuc
            if (nucGapPos == 0) {
                // cur gap coord is the first gap nuc, return the previous main nuc
                coord = std::make_tuple(blockId, nucPos - 1, -1);
            } else {
                // cur gap coord is not the first gap nuc, return the previous gap nuc
                coord = std::make_tuple(blockId, nucPos, nucGapPos - 1);
            }
        }
    }

    void stepRightCoordinate(Coordinate& coord) const {
        if (coord == lastCoord) {
            coord = Coordinate(-1, -1, -1, -1);
            return;
        }

        if (coord == blockEdgeCoords[coord.primaryBlockId].end) {
            // go to the next block
            coord = blockEdgeCoords[coord.primaryBlockId + 1].start;
            return;
        }

        // not end of block
        if (coord.nucGapPosition == -1) {
            // cur coord is a main nuc, step right to the next main nuc
            coord.nucPosition++;
            const auto& nextNucGaps = globalCoords[coord.primaryBlockId][coord.nucPosition].second;
            if (!nextNucGaps.empty()) {
                // next main nuc has gap nucs -> step to the first gap nuc
                coord.nucGapPosition = 0;
            }
        } else {
            // cur coord is at a gap nuc
            const auto& currentNucGaps = globalCoords[coord.primaryBlockId][coord.nucPosition].second;
            if (coord.nucGapPosition == currentNucGaps.size() - 1) {
                // cur gap nuc is the last gap nuc -> step to the main nuc
                coord.nucGapPosition = -1;
            } else {
                // cur gap nuc is not the last gap nuc -> step to the next gap nuc
                coord.nucGapPosition++;
            }
        }
    }

    void stepLeftCoordinate(Coordinate& coord) const {
        if (coord == firstCoord) {
            coord = Coordinate(-1, -1, -1, -1);
            return;
        }

        if (coord == blockEdgeCoords[coord.primaryBlockId].start) {
            // go to the previous block
            coord = blockEdgeCoords[coord.primaryBlockId - 1].end;
            return;
        }

        // not start of block
        if (coord.nucGapPosition == -1) {
            // cur coord is a main nuc
            const auto& currentNucGaps = globalCoords[coord.primaryBlockId][coord.nucPosition].second;
            if (currentNucGaps.empty()) {
                // cur coord has no gap nucs, return the previous main nuc
                coord.nucPosition--;
            } else {
                // cur coord has gap nucs, return the last gap nuc
                coord.nucGapPosition = globalCoords[coord.primaryBlockId][coord.nucPosition].second.size() - 1;
            }
        } else {
            // cur coord is at a gap nuc
            if (coord.nucGapPosition == 0) {
                // cur gap coord is the first gap nuc, return the previous main nuc
                coord.nucPosition--;
                coord.nucGapPosition = -1;
            } else {
                // cur gap coord is not the first gap nuc, return the previous gap nuc
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

// Returns true if the mutation is from a canonical base to an ambiguous IUPAC code
// (excluding gap '-' and sentinel 'x')
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
    std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence = blockSequences.sequence;

    // process block mutations
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

    // process nuc mutations
    for (const auto& nucMutation : node->nucMutation) {
        int length = nucMutation.mutInfo >> 4;
        int blockId;
        int lastOffset = -1;

        for (int i = 0; i < length; i++) {
            Coordinate pos = Coordinate(nucMutation, i);
            if ((pos.nucPosition == sequence[pos.primaryBlockId].size() - 1 && pos.nucGapPosition == -1) ||
                (pos.nucPosition >= sequence[pos.primaryBlockId].size())) {
                continue;
            }
            lastOffset = i;
            blockId = pos.primaryBlockId;
            const char oldNuc = blockSequences.getSequenceBase(pos);
            const int newNucCode = (nucMutation.nucs >> (4 * (5 - i))) & 0xF;
            const char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);

            if (oldNuc == newNuc) continue;

            // Imputation: when enabled, skip mutations where a canonical base
            // becomes an ambiguous IUPAC code (not gap/sentinel).
            // This effectively inherits the parent's base for syncmer computation.
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

// Returns (firstNonGapScalar, lastNonGapScalar) - the bounds of actual sequence data
// flankSize: number of non-gap bases to skip at each end (default 0)
// If genome is all gaps or too short, returns (UINT64_MAX, 0)

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

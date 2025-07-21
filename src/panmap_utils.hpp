#pragma once

#include "panmanUtils.hpp"
#include "logging.hpp"
#include <string>

namespace panmapUtils {


void getSequenceFromReference(
  panmanUtils::Tree* tree,
  std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  std::vector<bool>& blockExists,
  std::vector<bool>& blockStrand,
  std::unordered_map<int, int>& blockLengths,
  std::string reference
);

std::string getStringFromSequence(
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  const std::vector<bool>& blockExists,
  const std::vector<bool>& blockStrand,
  bool aligned
);

std::string getStringFromReference(
  panmanUtils::Tree* tree,
  std::string reference,
  bool aligned
);

struct Coordinate {
    int32_t nucPosition;
    int32_t nucGapPosition;
    int32_t primaryBlockId;
    int32_t secondaryBlockId;

    // Default constructor
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

    // Get base corresponding to this Coordinate's position within a sequence_t
    char getSequenceBase(const std::vector<std::pair<std::vector<std::pair<char,std::vector<char>>>, std::vector<std::vector<std::pair<char, std::vector<char>>>>>>& seq) const {
        if(secondaryBlockId != -1) {
            if(nucGapPosition != -1) {
                return seq[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            } else {
                return seq[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            }
        } else {
            if(nucGapPosition != -1) {
                return seq[primaryBlockId].first[nucPosition].second[nucGapPosition];
            } else {
                return seq[primaryBlockId].first[nucPosition].first;
            }
        }
    }

    char getSequenceBase(const std::vector< std::vector< std::pair< char, std::vector< char > > > >& sequence) const {
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
    void setSequenceBase( std::vector<std::pair<std::vector<std::pair<char,std::vector<char>>>, std::vector<std::vector<std::pair<char, std::vector<char>>>>>>& seq, char newNuc) const {
        if(secondaryBlockId != -1) {
            if(nucGapPosition != -1) {
                seq[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId].second[secondaryBlockId][nucPosition].first = newNuc;
            }
        } else {
            if(nucGapPosition != -1) {
                seq[primaryBlockId].first[nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId].first[nucPosition].first = newNuc;
            }
        }
    }

    void setSequenceBase(std::vector< std::vector< std::pair< char, std::vector< char > > > >& seq, char newNuc) const {
        if(secondaryBlockId != -1) {
            if(nucGapPosition != -1) {
                seq[primaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId][nucPosition].first = newNuc;
            }
        } else {
            if(nucGapPosition != -1) {
                seq[primaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId][nucPosition].first = newNuc;
            }
        }
    }

    void setSequenceBase(std::vector<std::pair<char, std::vector<char>>>& blockSequence, char newNuc) const {
      if(secondaryBlockId != -1) {
        if(nucGapPosition != -1) {
          blockSequence[nucPosition].second[nucGapPosition] = newNuc;
        } else {
          blockSequence[nucPosition].first = newNuc;
        }
      } else {
        if(nucGapPosition != -1) {
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
        return nucPosition == other.nucPosition &&
               nucGapPosition == other.nucGapPosition &&
               primaryBlockId == other.primaryBlockId &&
               secondaryBlockId == other.secondaryBlockId;
    }
};

struct BlockSequences {
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequence;
  std::vector<bool> blockExists;
  std::vector<bool> blockStrand;
  
  BlockSequences() {}

  BlockSequences(panmanUtils::Tree *tree) {
    sequence.resize(tree->blocks.size() + 1);
    blockExists.resize(tree->blocks.size() + 1, false);
    blockStrand.resize(tree->blocks.size() + 1, true);

    int32_t maxBlockId = 0;

    for(size_t i = 0; i < tree->blocks.size(); i++) {
      const auto& curBlock = tree->blocks[i];
      int32_t primaryBlockId = ((int32_t)curBlock.primaryBlockId);
      if (i != static_cast<size_t>(primaryBlockId)) {
        std::cerr << "primaryBlockId: " << primaryBlockId << " i: " << i << std::endl;
        std::exit(1);
      }
      maxBlockId = std::max(maxBlockId, primaryBlockId);
      for(size_t j = 0; j < curBlock.consensusSeq.size(); j++) {
        bool endFlag = false;
        for(size_t k = 0; k < 8; k++) {
          const int nucCode = (((curBlock.consensusSeq[j]) >> (4*(7 - k))) & 15);

          if(nucCode == 0) {
            endFlag = true;
            break;
          }
          // len++;
          const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
          sequence[primaryBlockId].push_back({nucleotide, {}});
        }

        if(endFlag) {
          break;
        }
      }
      // End character to incorporate for gaps at the end
      sequence[primaryBlockId].push_back({'x', {}});
    }

    // resize in case of early stop
    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    blockStrand.resize(maxBlockId + 1);

    // assign gaps
    for(size_t i = 0; i < tree->gaps.size(); i++) {
      const auto& curGap = tree->gaps[i];
      int32_t primaryBId = (curGap.primaryBlockId);
      int32_t secondaryBId = (curGap.secondaryBlockId);
      for(size_t j = 0; j < curGap.nucPosition.size(); j++) {
        int len = curGap.nucGapLength[j];
        int pos = curGap.nucPosition[j];
        sequence[primaryBId][pos].second.resize(len, '-');
      }
    }
  }

  int numBlocks() const {
    return sequence.size();
  }
  
  const bool getBlockStrand(int blockId) const {
    return blockStrand[blockId];
  }

  const bool getBlockExists(int blockId) const {
    return blockExists[blockId];
  }

  char getSequenceBase(const Coordinate& coord) const {
    if (coord.nucGapPosition != -1) {
      return sequence[coord.primaryBlockId][coord.nucPosition].second[coord.nucGapPosition];
    } else {
      return sequence[coord.primaryBlockId][coord.nucPosition].first;
    }
  }
  
  void setSequenceBase(const Coordinate& coord, char newNuc) {
    if(coord.nucGapPosition != -1) {
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

  std::string getAlignedSequenceString() const {
    return getSequenceString(true);
  }

  std::string getUnalignedSequenceString() const {
    return getSequenceString(false);
  }
};

struct GlobalCoords {
  std::vector<std::vector<std::pair<int64_t, std::vector<int64_t>>>> globalCoords;
  std::tuple<int64_t, int64_t, int64_t> firstTupleCoord;
  std::tuple<int64_t, int64_t, int64_t> lastTupleCoord;
  int64_t lastScalarCoord;

  GlobalCoords(const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence) {
    int64_t curScalarCoord = 0;
    globalCoords.resize(sequence.size());
    for (size_t i = 0; i < sequence.size(); i++) {
      globalCoords[i].resize(sequence[i].size());
      for (size_t j = 0; j < sequence[i].size(); j++) {
        // process gap nucs first
        globalCoords[i][j].second.resize(sequence[i][j].second.size());
        for (size_t k = 0; k < sequence[i][j].second.size(); k++) {
          globalCoords[i][j].second[k] = curScalarCoord;
          curScalarCoord++;
        }
        
        // process main nuc
        if (sequence[i][j].first == 'x') {
          // skip if main nuc is x
          globalCoords[i][j].first = -1;
        } else {
          globalCoords[i][j].first = curScalarCoord;
          curScalarCoord++;
        }
      }
    }

    lastScalarCoord = curScalarCoord - 1;
    const auto& lastBlock = sequence.back();
    if (lastBlock.back().second.empty()) {
      // last block's last nuc gap is empty, take the second to last main nuc
      lastTupleCoord = std::make_tuple(sequence.size() - 1, lastBlock.size() - 2, -1);
    } else {
      // last block's last nuc gap is not empty, take the last nuc gap
      lastTupleCoord = std::make_tuple(sequence.size() - 1, lastBlock.size() - 1, lastBlock.back().second.size() - 1);
    }
    
    const auto& firstBlock = sequence[0];
    if (firstBlock[0].second.empty()) {
      // first block's first nuc gap is empty, take the first main nuc
      firstTupleCoord = std::make_tuple(0, 0, -1);
    } else {
      // first block's first nuc gap is not empty, take the first gap nuc
      firstTupleCoord = std::make_tuple(0, 0, 0);
    }

    // sanity check
    if (getScalarFromTuple(firstTupleCoord) != 0) {
      logging::err("firstScalarCoord != 0");
      std::exit(1);
    }
    if (getScalarFromTuple(lastTupleCoord) != lastScalarCoord) {
      logging::err("lastScalarCoord != getScalarFromTuple(lastTupleCoord)");
      std::exit(1);
    }
  }

  int64_t getScalarFromTuple(const std::tuple<int64_t, int64_t, int64_t> &tupleCoord) const {
    const auto& [blockId, nucPos, nucGapPos] = tupleCoord;
    if (nucGapPos == -1) {
      return globalCoords[blockId][nucPos].first;
    }
    return globalCoords[blockId][nucPos].second[nucGapPos];
  }

  int64_t getScalarFromCoord(const Coordinate& coord) const {
    if (coord.nucGapPosition == -1) {
      return globalCoords[coord.primaryBlockId][coord.nucPosition].first;
    }
    return globalCoords[coord.primaryBlockId][coord.nucPosition].second[coord.nucGapPosition];
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
    return getScalarFromTuple(getBlockStartTuple(blockId));
  }

  int64_t getBlockEndScalar(int64_t blockId) const {
    return getScalarFromTuple(getBlockEndTuple(blockId));
  }

  std::tuple<int64_t, int64_t, int64_t> stepRight(const std::tuple<int64_t, int64_t, int64_t> &coord) const {
    const auto& [blockId, nucPos, nucGapPos] = coord;
    if (coord == lastTupleCoord) {
      return std::make_tuple(-1, -1, -1);
    }

    if (coord == getBlockEndTuple(blockId)) {
      // go to the next block
      return getBlockStartTuple(blockId + 1);
    }

    // not end of block
    if (nucGapPos == -1) {
      // cur coord is a main nuc, step right to the next main nuc
      if (globalCoords[blockId][nucPos + 1].second.empty()) {
        // next main nuc has no gap nucs -> return the main nuc
        return std::make_tuple(blockId, nucPos + 1, -1);
      } else {
        // next main nuc has gap nucs -> return the first gap nuc
        return std::make_tuple(blockId, nucPos + 1, 0);
      }

    } else {
      // cur coord is at a gap nuc
      if (nucGapPos == globalCoords[blockId][nucPos].second.size() - 1) {
        // cur gap nuc is the last gap nuc -> return the main nuc
        return std::make_tuple(blockId, nucPos, -1);
      } else {
        // cur gap nuc is not the last gap nuc -> return the next gap nuc
        return std::make_tuple(blockId, nucPos, nucGapPos + 1);
      }
    }
    return std::make_tuple(-1, -1, -1);
  }

  std::tuple<int64_t, int64_t, int64_t> stepLeft(const std::tuple<int64_t, int64_t, int64_t> &coord) const {
    const auto& [blockId, nucPos, nucGapPos] = coord;
    if (coord == firstTupleCoord) {
      return std::make_tuple(-1, -1, -1);
    }

    if (coord == getBlockStartTuple(blockId)) {
      // go to the previous block
      return getBlockEndTuple(blockId - 1);
    }

    // not start of block
    if (nucGapPos == -1) {
      // cur coord is a main nuc
      if (globalCoords[blockId][nucPos].second.empty()) {
        // cur coord has no gap nucs, return the previous main nuc
        return std::make_tuple(blockId, nucPos - 1, -1);
      } else {
        // cur coord has gap nucs, return the last gap nuc
        return std::make_tuple(blockId, nucPos, globalCoords[blockId][nucPos].second.size() - 1);
      }
    } else {
      // cur coord is at a gap nuc
      if (nucGapPos == 0) {
        // cur gap coord is the first gap nuc, return the previous main nuc
        return std::make_tuple(blockId, nucPos - 1, -1);
      } else {
        // cur gap coord is not the first gap nuc, return the previous gap nuc
        return std::make_tuple(blockId, nucPos, nucGapPos - 1);
      }
    }
    return std::make_tuple(-1, -1, -1);
  }

  Coordinate stepRight(const Coordinate &coord) const {
    std::tuple<int64_t, int64_t, int64_t> tupleCoord = std::make_tuple(coord.primaryBlockId, coord.nucPosition, coord.nucGapPosition);
    tupleCoord = stepRight(tupleCoord);
    return Coordinate(std::get<1>(tupleCoord), std::get<2>(tupleCoord), std::get<0>(tupleCoord), coord.secondaryBlockId);
  }

  Coordinate stepLeft(const Coordinate &coord) const {
    std::tuple<int64_t, int64_t, int64_t> tupleCoord = std::make_tuple(coord.primaryBlockId, coord.nucPosition, coord.nucGapPosition);
    tupleCoord = stepLeft(tupleCoord);
    return Coordinate(std::get<1>(tupleCoord), std::get<2>(tupleCoord), std::get<0>(tupleCoord), coord.secondaryBlockId);
  }

};


}
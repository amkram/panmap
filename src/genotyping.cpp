#include "genotyping.hpp"
#include "conversion.hpp"
#include "panmanUtils.hpp"
#include "panmap_utils.hpp"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <istream>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <boost/icl/interval_set.hpp>

using namespace std;
using namespace genotyping;

// Function implementations from the header file
void genotyping::stringSplit(const std::string& str, char delimiter, std::vector<std::string>& out) {
  out.clear();
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter)) {
    if (!token.empty()) {
      out.push_back(token);
    }
  }
}

void genotyping::fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf) {
  std::string line;
  int idx = 0;
  while (getline(inf, line)) {
    std::vector<double> probs;
    std::vector<std::string> fields;
    stringSplit(line, ' ', fields);
    
    if (fields.empty()) {
      break;
    }

    if (idx < 4) {
      if (fields.size() != 4) {
        throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
      }

      for (const auto &f : fields) {
        probs.push_back(std::stod(f));
      }
      mutMat.submat[idx] = std::move(probs);
    } else if (idx == 4) {
      if (fields.empty()) {
        throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
      }

      for (const auto& f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        if (subFields.size() != 2) {
          throw std::invalid_argument("Invalid format in mutation matrix file");
        }
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.insmat[size] = prob;
      }
    } else if (idx == 5) {
      if (fields.empty()) {
        throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
      }

      for (const auto& f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        if (subFields.size() != 2) {
          throw std::invalid_argument("Invalid format in mutation matrix file");
        }
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.delmat[size] = prob;
      }
    }
    idx++;
  }

  if (idx != 6) {
    throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
  }
  
  // Set maximum log probabilities for insertions and deletions
  mutMat.maxInsLogProb = 100.0; // Default high penalty
  mutMat.maxDelLogProb = 100.0; // Default high penalty
  
  // Calculate actual maximum values if matrices are not empty
  if (!mutMat.insmat.empty()) {
    mutMat.maxInsLogProb = std::max_element(
      mutMat.insmat.begin(), mutMat.insmat.end(),
      [](const auto& a, const auto& b) { return a.second < b.second; }
    )->second;
  }
  
  if (!mutMat.delmat.empty()) {
    mutMat.maxDelLogProb = std::max_element(
      mutMat.delmat.begin(), mutMat.delmat.end(),
      [](const auto& a, const auto& b) { return a.second < b.second; }
    )->second;
  }
  
  mutMat.filled = true;
}

void genotyping::buildMutationMatricesHelper(
    mutationMatrices &mutMat,
    panmanUtils::Tree *tree,
    panmanUtils::Node *node,
    std::vector<int64_t> &parentBaseCounts,
    std::vector<int64_t> &totalBaseCounts,
    std::vector<std::vector<int64_t>> &subCount,
    std::unordered_map<int64_t, int64_t> &insCount,
    std::unordered_map<int64_t, int64_t> &delCount) {
    
  if (!node) return;
  
  // Process current node mutations
  std::vector<int64_t> currentBaseCounts(4, 0);
  
  // Get sequence for this node
  std::string nodeSeq = tree->getStringFromReference(node->identifier, false, true);
  
  // Count bases in current node
  for (char c : nodeSeq) {
    if (c == 'A' || c == 'a') currentBaseCounts[0]++;
    else if (c == 'C' || c == 'c') currentBaseCounts[1]++;
    else if (c == 'G' || c == 'g') currentBaseCounts[2]++;
    else if (c == 'T' || c == 't') currentBaseCounts[3]++;
  }
  
  // If this is not the root, compare with parent
  if (node != tree->root && node->parent) {
    std::string parentSeq = tree->getStringFromReference(node->parent->identifier, false, true);
    
    // Simple sequence alignment and mutation counting
    // In a real implementation, this would need proper sequence alignment
    size_t i = 0, j = 0;
    while (i < parentSeq.length() && j < nodeSeq.length()) {
      char p = parentSeq[i];
      char n = nodeSeq[j];
      
      // Skip non-ACGT characters
      if ((p != 'A' && p != 'C' && p != 'G' && p != 'T' && 
           p != 'a' && p != 'c' && p != 'g' && p != 't')) {
        i++;
        continue;
      }
      
      if ((n != 'A' && n != 'C' && n != 'G' && n != 'T' && 
           n != 'a' && n != 'c' && n != 'g' && n != 't')) {
        j++;
        continue;
      }
      
      // Convert to upper case for comparison
      p = std::toupper(p);
      n = std::toupper(n);
      
      // Map nucleotides to indices
      int pIdx = (p == 'A') ? 0 : ((p == 'C') ? 1 : ((p == 'G') ? 2 : 3));
      int nIdx = (n == 'A') ? 0 : ((n == 'C') ? 1 : ((n == 'G') ? 2 : 3));
      
      // Count substitutions
      if (p != n) {
        subCount[pIdx][nIdx]++;
        totalBaseCounts[pIdx]++;
      }
      
      i++;
      j++;
    }
    
    // Count insertions
    if (nodeSeq.length() > parentSeq.length()) {
      int64_t insSize = nodeSeq.length() - parentSeq.length();
      insCount[insSize]++;
    }
    
    // Count deletions
    if (parentSeq.length() > nodeSeq.length()) {
      int64_t delSize = parentSeq.length() - nodeSeq.length();
      delCount[delSize]++;
    }
  }
  
  // Update parent base counts for children
  parentBaseCounts = currentBaseCounts;
  
  // Recursively process all children
  for (auto* child : node->children) {
    buildMutationMatricesHelper(mutMat, tree, child, parentBaseCounts, 
                           totalBaseCounts, subCount, insCount, delCount);
  }
}

static int getIndexFromNucleotide(char nuc) {
  switch (nuc) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  case '*':
    return 4;
  default:
    return 5;
  }
  return 5;
}

static char getNucFromTuple(const std::tuple<int64_t, int64_t, int64_t> &tupleCoord, const std::vector<std::vector<std::pair<char, std::vector<char>>>> &sequence) {
  const auto& [blockId, nucPos, nucGapPos] = tupleCoord;
  if (nucGapPos == -1) {
    return sequence[blockId][nucPos].first;
  }
  return sequence[blockId][nucPos].second[nucGapPos];
}



void clearBlockBaseCounts(
  int32_t blockId,
  bool blockStrand,
  std::vector<int64_t> &curBaseCounts,
  std::vector<int64_t> &parentBaseCountsBacktrack,
  const panmapUtils::GlobalCoords &globalCoords,
  const std::vector<std::vector<std::pair<char, std::vector<char>>>> &sequence
) {
  std::tuple<int64_t, int64_t, int64_t> coord = globalCoords.getBlockStartTuple(blockId);
  std::tuple<int64_t, int64_t, int64_t> end   = globalCoords.getBlockEndTuple(blockId);

  while (true) {
    char c = getNucFromTuple(coord, sequence);
    if (!blockStrand) {
      c = panmanUtils::getComplementCharacter(c);
    }

    if (getIndexFromNucleotide(c) <= 3) {
      --curBaseCounts[getIndexFromNucleotide(c)];
      ++parentBaseCountsBacktrack[getIndexFromNucleotide(c)];
    }

    if (coord == end) {
      break;
    }

    globalCoords.stepRightCoordinate(coord);
    if (coord == std::make_tuple(-1, -1, -1)) {
      logging::err("stepRight returned -1, -1, -1");
      std::exit(1);
    }
  }
}

static void applyMutations(
    panmanUtils::Tree *tree,
    panmanUtils::Node *node,
    panmapUtils::BlockSequences &blockSequences,
    std::vector<char> &oldBlockExists,
    std::vector<char> &oldBlockStrand,
    std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> &blockMutationRecord,
    std::vector<std::tuple<panmapUtils::Coordinate, char, char>> &nucMutationRecord,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunUpdates,
    std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> &gapRunBacktracks,
    std::unordered_set<int64_t> &invertedBlocks,
    std::vector<std::pair<int64_t, bool>> &invertedBlocksBacktracks,
    std::vector<int64_t> &curBaseCounts,
    std::vector<int64_t> &parentBaseCountsBacktrack,
    std::vector<std::vector<int64_t>> &subCount,
    bool &isMutatedNuc,
    const panmapUtils::GlobalCoords &globalCoords,
    const std::vector<std::pair<int64_t, int64_t>> &blockRanges
) {
  std::vector<char>& blockExists = blockSequences.blockExists;
  std::vector<char>& blockStrand = blockSequences.blockStrand;

  std::vector<std::vector<std::pair<char, std::vector<char>>>> &sequence = blockSequences.sequence;
  for (const auto& blockMutation : node->blockMutation) {
    const int32_t& blockId  = blockMutation.primaryBlockId;
    const bool& isInsertion = blockMutation.blockMutInfo;
    const bool& isInversion = blockMutation.inversion;
    bool oldExists = blockExists[blockId];
    bool oldStrand = blockStrand[blockId];
    if (isInsertion) {
      // insertion
      blockExists[blockId] = true;
      blockStrand[blockId] = !isInversion;
      if (!blockStrand[blockId]) {
        invertedBlocks.insert(blockId);
        if (oldStrand) {
          invertedBlocksBacktracks.emplace_back(blockId, true);
        }
      }
    } else if (isInversion) {
      // simple inversion
      blockStrand[blockId] = !blockStrand[blockId];
      if (!blockStrand[blockId]) {
        invertedBlocks.insert(blockId);
        if (oldStrand) {
          invertedBlocksBacktracks.emplace_back(blockId, true);
        }
      } else {
        invertedBlocks.erase(blockId);
        if (!oldStrand) {
          invertedBlocksBacktracks.emplace_back(blockId, false);
        }
      }
      std::vector<int64_t> oldCurBaseCounts = curBaseCounts;
      clearBlockBaseCounts(blockId, oldBlockStrand[blockId], curBaseCounts, parentBaseCountsBacktrack, globalCoords, sequence);
    } else {
      // deletion
      blockExists[blockId] = false;
      blockStrand[blockId] = true;
      if (!oldStrand) {
        invertedBlocks.erase(blockId);
        invertedBlocksBacktracks.emplace_back(blockId, false);
      }
      std::vector<int64_t> oldCurBaseCounts = curBaseCounts;
      clearBlockBaseCounts(blockId, oldBlockStrand[blockId], curBaseCounts, parentBaseCountsBacktrack, globalCoords, sequence);
    }
    blockMutationRecord.emplace_back(std::make_tuple(blockId, oldExists, oldStrand, blockExists[blockId], blockStrand[blockId]));
  }

  for (const auto& nucMutation : node->nucMutation) {
    int length = nucMutation.mutInfo >> 4;
    for (int i = 0; i < length; i++) {
      panmapUtils::Coordinate pos = panmapUtils::Coordinate(nucMutation, i);
      char originalNuc = blockSequences.getSequenceBase(pos);
      int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
      char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
      blockSequences.setSequenceBase(pos, newNuc);
      nucMutationRecord.emplace_back(std::make_tuple(pos, originalNuc, newNuc));
    }
  }

  for (auto &nucMutation : nucMutationRecord) {
    const auto& [coord, originalNuc, newNuc] = nucMutation;
    int blockId = coord.primaryBlockId;
    if (oldBlockExists[blockId] && blockExists[blockId] && oldBlockStrand[blockId] == blockStrand[blockId]) {
      // on to on -> collect gap runs and nuc runs
      char parChar = originalNuc == 'x' ? '-' : originalNuc;
      char curChar = newNuc == 'x' ? '-' : newNuc;

      int64_t scalar = globalCoords.getScalarFromCoord(coord);
      if (!blockStrand[blockId]) {
        scalar = blockRanges[blockId].first + blockRanges[blockId].second - scalar;
      }

      if (!oldBlockStrand[blockId]) parChar = panmanUtils::getComplementCharacter(parChar);
      if (!blockStrand[blockId])    curChar =  panmanUtils::getComplementCharacter(curChar);

      if (getIndexFromNucleotide(parChar) <= 3) {
        --curBaseCounts[getIndexFromNucleotide(parChar)];
        ++parentBaseCountsBacktrack[getIndexFromNucleotide(parChar)];
      }
      if (getIndexFromNucleotide(curChar) <= 3) {
        ++curBaseCounts[getIndexFromNucleotide(curChar)];
        --parentBaseCountsBacktrack[getIndexFromNucleotide(curChar)];
      }
      if (getIndexFromNucleotide(parChar) <= 3 && getIndexFromNucleotide(curChar) <= 3) {
        isMutatedNuc = true;
        ++subCount[getIndexFromNucleotide(parChar)][getIndexFromNucleotide(curChar)];
      }


      if (parChar != '-' && curChar == '-') {
        // nuc to gap
        if (!gapRunUpdates.empty() && gapRunUpdates.back().first == true && gapRunUpdates.back().second.second + 1 == scalar) {
          ++(gapRunUpdates.back().second.second);
        }
        else {
          gapRunUpdates.emplace_back(true, std::make_pair(scalar, scalar)); 
        }
      } else if (parChar == '-' && curChar != '-') {
        // gap to nuc
        if (!gapRunUpdates.empty() && gapRunUpdates.back().first == false && gapRunUpdates.back().second.second + 1 == scalar) {
          ++(gapRunUpdates.back().second.second);
        } else {
          gapRunUpdates.emplace_back(false, std::make_pair(scalar, scalar));
        }
      }
    }
  }
}

static void updateGapMapStep(std::map<int64_t, int64_t>& gapMap, const std::pair<bool, std::pair<int64_t, int64_t>>& update, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates, bool recordGapMapUpdates=true) {
  bool toGap = update.first;
  int64_t start = update.second.first;
  int64_t end = update.second.second;

  auto rightIt = gapMap.upper_bound(start);
  auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

  bool rightItExists = rightIt != gapMap.end();
  bool leftItExists = leftIt != gapMap.end();

  if (toGap) {
    // add gap range
    if (gapMap.empty()) {
      gapMap[start] = end;
      backtrack.emplace_back(true, std::make_pair(start, end));
      if (recordGapMapUpdates) {
        gapMapUpdates.emplace_back(false, std::make_pair(start, end));
      }
      return;
    }
    
    decltype(rightIt) curIt;

    // curIt starts outside of any range
    if (!leftItExists || (!rightItExists && start > leftIt->second) || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
      if (leftItExists && start == leftIt->second + 1) {
        // 1 base after left range and merge with left
        curIt = leftIt;
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        curIt->second = end;
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        }
      } else {
        // insert new range
        auto tmpIt = gapMap.emplace(start, end);
        curIt = tmpIt.first;
        backtrack.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        }
      }
    } else {
      curIt = leftIt;
      if (end <= curIt->second) {
        return;
      }
      backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
      curIt->second = end;
      if (recordGapMapUpdates) {
        gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
      }
    }

    auto nextIt = std::next(curIt);
    while (true) {
      if (nextIt == gapMap.end()) {
        break;
      }

      if (nextIt->second <= curIt->second) {
        auto tmpIt = nextIt;
        nextIt = std::next(nextIt);
        backtrack.emplace_back(false, std::make_pair(tmpIt->first, tmpIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(true, std::make_pair(tmpIt->first, tmpIt->second));
        }
        gapMap.erase(tmpIt);
      } else if (nextIt->first <= end + 1) {
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        curIt->second = nextIt->second;
        backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          gapMapUpdates.emplace_back(true, std::make_pair(nextIt->first, nextIt->second));
        }
        gapMap.erase(nextIt);
        break;
      } else {
        break;
      }
    }
  } else {
    // remove gap range
    if (gapMap.empty() || (!leftItExists && end < rightIt->first) || (!rightItExists && start > leftIt->second)) {
      return;
    }

    decltype(rightIt) curIt;
    decltype(rightIt) nextIt;
    if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
      // curIt starts outside of any range
      curIt = rightIt;

      if (end < curIt->first) {
        return;
      }

      // ends within the curIt range
      if (end <= curIt->second) {
        if (end == curIt->second) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else {
          gapMap[end+1] = curIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        }
        return;
      } else {
        nextIt = std::next(curIt);
        backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        }
        gapMap.erase(curIt);
      }
      
    } else {
      // curIt starts inside of a range
      curIt = leftIt;
      
      if (end <= curIt->second) {
        // contained in the curIt range
        if (start == curIt->first && end == curIt->second) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else if (start == curIt->first) {
          gapMap[end + 1] = curIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else if (end == curIt->second) {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          curIt->second = start - 1;
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          }
        } else {
          gapMap[end + 1] = curIt->second;
          backtrack.emplace_back(true, std::make_pair(end+1, curIt->second));
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(end+1, curIt->second));
            gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, start-1));
          }
          curIt->second = start - 1;
        }
        return;
      } else {
        if (start == curIt->first) {
          nextIt = std::next(curIt);
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          }
          gapMap.erase(curIt);
        } else {
          backtrack.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          curIt->second = start - 1;
          if (recordGapMapUpdates) {
            gapMapUpdates.emplace_back(false, std::make_pair(curIt->first, curIt->second));
          }
          nextIt = std::next(curIt);
        }
      }
    }

    
    while (true) {
      if (nextIt == gapMap.end()) {
        break;
      }

      if (nextIt->first > end) {
        break;
      } else if (nextIt->second <= end) {
        auto tmpIt = nextIt;
        nextIt = std::next(nextIt);
        backtrack.emplace_back(false, std::make_pair(tmpIt->first, tmpIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(true, std::make_pair(tmpIt->first, tmpIt->second));
        }
        gapMap.erase(tmpIt);
      } else {
        gapMap[end + 1] = nextIt->second;
        backtrack.emplace_back(true, std::make_pair(end+1, nextIt->second));
        backtrack.emplace_back(false, std::make_pair(nextIt->first, nextIt->second));
        if (recordGapMapUpdates) {
          gapMapUpdates.emplace_back(false, std::make_pair(end+1, nextIt->second));
          gapMapUpdates.emplace_back(true, std::make_pair(nextIt->first, nextIt->second));
        }
        gapMap.erase(nextIt);
        break;
      }
    }
  }
}

static void updateGapMap(std::map<int64_t, int64_t>& gapMap, const std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& updates, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates) {
  for (const auto& update : updates) {
    updateGapMapStep(gapMap, update, backtrack, gapMapUpdates);
  }
}

static void invertGapMap(std::map<int64_t, int64_t>& gapMap, const std::pair<int64_t, int64_t>& invertRange, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& backtrack, std::vector<std::pair<bool, std::pair<int64_t, int64_t>>>& gapMapUpdates) {
  const auto& [start, end] = invertRange;

  auto rightIt = gapMap.upper_bound(start);
  auto leftIt = (rightIt == gapMap.begin()) ? gapMap.end() : std::prev(rightIt);

  bool rightItExists = rightIt != gapMap.end();
  bool leftItExists = leftIt != gapMap.end();

  // completely inside or outside a gap range -> do nothing
  if (
    gapMap.empty() || // empty gap map
    (!leftItExists && end < rightIt->first) || // completely left of first gap range
    (!rightItExists && start > leftIt->second) // completely right of last gap range
    // (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first) || // completely between two gap ranges
    // (leftItExists && start >= leftIt->first && end <= leftIt->second) // completely inside a gap range
  ) {
    return;
  }
  
  //                    gaps            beg      end
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> blockRuns;
  if (!leftItExists || (leftItExists && start > leftIt->second && rightItExists && start < rightIt->first)) {
    // start outside of a range
    auto curIt = rightIt;
    blockRuns.emplace_back(false, std::make_pair(start, curIt->first - 1));
    if (end <= curIt->second) {
      blockRuns.emplace_back(true, std::make_pair(blockRuns.back().second.second + 1, end));
    } else {
      blockRuns.emplace_back(true, std::make_pair(blockRuns.back().second.second + 1, curIt->second));
      curIt = std::next(curIt);
      while (true) {
        if (curIt == gapMap.end()) {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
          break;
        } else if (end < curIt->first) {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
          break;
        } else if (end > curIt->second) {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
          blockRuns.emplace_back(true, std::make_pair(curIt->first, curIt->second));
          curIt = std::next(curIt);
        } else {
          blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
          blockRuns.emplace_back(true, std::make_pair(curIt->first, end));
          break;
        }
      }
    }
  } else {
    // start inside of a range
    auto curIt = leftIt;
    blockRuns.emplace_back(true, std::make_pair(start, curIt->second));
    curIt = std::next(curIt);
    while (true) {
      if (curIt == gapMap.end()) {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
        break;
      } else if (end < curIt->first) {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, end));
        break;
      } else if (end > curIt->second) {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
        blockRuns.emplace_back(true, std::make_pair(curIt->first, curIt->second));
        curIt = std::next(curIt);
      } else {
        blockRuns.emplace_back(false, std::make_pair(blockRuns.back().second.second + 1, curIt->first - 1));
        blockRuns.emplace_back(true, std::make_pair(curIt->first, end));
        break;
      }
    }
  }

  int64_t curBeg = blockRuns.front().second.first;
  for (auto it = blockRuns.rbegin(); it != blockRuns.rend(); ++it) {
    int64_t curEnd = curBeg + (it->second.second - it->second.first);
    updateGapMapStep(gapMap, {it->first, {curBeg, curEnd}}, backtrack, gapMapUpdates, false);
    curBeg = curEnd + 1;
  }
}

static std::vector<std::pair<int64_t, int64_t>> invertRanges(const std::vector<std::pair<int64_t, int64_t>>& nucRanges, const std::pair<int64_t, int64_t>& invertRange) {
  std::vector<std::pair<int64_t, int64_t>> invertedRanges;

  auto [start, end] = invertRange;

  for (auto it = nucRanges.rbegin(); it != nucRanges.rend(); ++it) {
    const auto& [curStart, curEnd] = *it;
    invertedRanges.emplace_back(start + end - curEnd, start + end - curStart);
  }

  return invertedRanges;
}

static boost::icl::interval_set<int64_t> gapMapToNucRunSet(const std::map<int64_t, int64_t> &gapMap, const std::vector<std::pair<int64_t, int64_t>>& blockRanges) {
  boost::icl::interval_set<int64_t> curNucRunSet;
  int64_t start = -1;
  for (const auto& [curStart, curEnd] : gapMap) {
    if (start == -1) {
      if (curStart != 0) {
        curNucRunSet.add(boost::icl::interval<int64_t>::closed(0, curStart - 1));
      }
    } else {
      curNucRunSet.add(boost::icl::interval<int64_t>::closed(start, curStart - 1));
    }
    start = curEnd + 1;
  }
  if (start <= blockRanges.back().second) {
    curNucRunSet.add(boost::icl::interval<int64_t>::closed(start, blockRanges.back().second));
  }
  return curNucRunSet;
}

void undoMutations(
  panmapUtils::BlockSequences &blockSequences,
  std::vector<char> &oldBlockExists,
  std::vector<char> &oldBlockStrand,
  std::unordered_set<int64_t> &invertedBlocks,
  const std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> &blockMutationRecord,
  const std::vector<std::tuple<panmapUtils::Coordinate, char, char>> &nucMutationRecord,
  const std::vector<std::pair<int64_t, bool>> &invertedBlocksBacktracks
) {
  std::vector<char>& blockExists = blockSequences.blockExists;
  std::vector<char>& blockStrand = blockSequences.blockStrand;

  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    blockExists[blockId] = oldExists;
    blockStrand[blockId] = oldStrand;
    oldBlockExists[blockId] = oldExists;
    oldBlockStrand[blockId] = oldStrand;
  }

  for (const auto& [coord, oldNuc, newNuc] : nucMutationRecord) {
    blockSequences.setSequenceBase(coord, oldNuc);
  }

  for (const auto& [blockId, erase] : invertedBlocksBacktracks) {
    if (erase) {
      invertedBlocks.erase(blockId);
    } else {
      invertedBlocks.insert(blockId);
    }
  }
}

void buildMutationMatricesHelper_test(
  mutationMatrices &mutMat,
  panmanUtils::Tree *tree,
  panmanUtils::Node *node,
  panmapUtils::BlockSequences &blockSequences,
  std::vector<char> &oldBlockExists,
  std::vector<char> &oldBlockStrand,
  std::map<int64_t, int64_t> &gapMap,
  std::unordered_set<int64_t> &invertedBlocks,
  std::vector<int64_t> &parentBaseCounts,
  std::vector<int64_t> &totalBaseCounts,
  std::vector<std::vector<int64_t>> &subCount,
  std::unordered_map<int64_t, int64_t> &insCount,
  std::unordered_map<int64_t, int64_t> &delCount,
  const panmapUtils::GlobalCoords &globalCoords,
  const std::vector<std::pair<int64_t, int64_t>> &blockRanges
) {
  // // for checking with brute force. DELETE THESE AFTER CONFIRMING EVERYTHING WORKS CORRECTLY!!!!!! ----------------------------
  // panmapUtils::BlockSequences oldBlockSequences = blockSequences;
  // std::vector<int64_t> oldParentBaseCounts = parentBaseCounts;
  // std::vector<std::vector<int64_t>> oldSubCount = subCount;
  // // ------------------------------------------------------------------------------------------------


  std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> blockMutationRecord;
  std::vector<std::tuple<panmapUtils::Coordinate, char, char>> nucMutationRecord;

  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunUpdates;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBacktracks;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapMapUpdates;
  std::vector<std::pair<int64_t, bool>> invertedBlocksBacktracks;

  std::vector<int64_t> curBaseCounts = parentBaseCounts;
  std::vector<int64_t> parentBaseCountsBacktrack(4);
  bool isMutatedNuc = false;

  applyMutations(
    tree, node,
    blockSequences, 
    oldBlockExists, oldBlockStrand,
    blockMutationRecord, nucMutationRecord,
    gapRunUpdates, gapRunBacktracks,
    invertedBlocks, invertedBlocksBacktracks,
    curBaseCounts, parentBaseCountsBacktrack,
    subCount, isMutatedNuc,
    globalCoords, blockRanges);

  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> gapRunBlocksBacktracks;

  auto parentGapMap = gapMap;
  std::vector<std::pair<bool, std::pair<int64_t, int64_t>>> tmpGapRunBlocksBacktracks;
  for (size_t i = 0; i < oldBlockStrand.size(); ++i) {
    if (!oldBlockStrand[i]) {
      invertGapMap(parentGapMap, blockRanges[i], tmpGapRunBlocksBacktracks, gapMapUpdates);
    }
  }

  boost::icl::interval_set<int64_t> parNucRunSet = gapMapToNucRunSet(parentGapMap, blockRanges);
  boost::icl::interval_set<int64_t> flippedSet;
  parentGapMap.clear();
  tmpGapRunBlocksBacktracks.clear();

  updateGapMap(gapMap, gapRunUpdates, gapRunBacktracks, gapMapUpdates);

  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    if (oldExists && !newExists) {
      // on to off -> block range to all gaps
      const auto& [start, end] = blockRanges[blockId];
      updateGapMapStep(gapMap, {true, {start, end}}, gapRunBacktracks, gapMapUpdates);
    } else if (!oldExists && newExists) {
      // off to on -> recompute across entire block
      panmapUtils::Coordinate coord = globalCoords.getBlockStartCoord(blockId);
      panmapUtils::Coordinate end = globalCoords.getBlockEndCoord(blockId);
      if (!blockSequences.getBlockStrand(blockId)) std::swap(coord, end);
      std::pair<int64_t, int64_t> curNucRange = {-1, -1};
      std::vector<std::pair<int64_t, int64_t>> nucRanges;
      std::vector<int64_t> oldCurBaseCounts = curBaseCounts;
      while (true) {
        char c = blockSequences.getSequenceBase(coord);
        c = c == 'x' ? '-' : c;
        if (!blockSequences.getBlockStrand(blockId)) c = panmanUtils::getComplementCharacter(c);
        int64_t scalar = globalCoords.getScalarFromCoord(coord);
        if (!blockSequences.getBlockStrand(blockId)) {
          scalar = blockRanges[blockId].first + blockRanges[blockId].second - scalar;
        }

        if (c != '-') {
          if (getIndexFromNucleotide(c) <= 3) {
            ++curBaseCounts[getIndexFromNucleotide(c)];
            --parentBaseCountsBacktrack[getIndexFromNucleotide(c)];
          }

          if (curNucRange.first != -1 && curNucRange.second + 1 == scalar) {
            ++curNucRange.second;
          } else {
            if (curNucRange.first != -1) {
              nucRanges.push_back(curNucRange);
            }
            curNucRange = {scalar, scalar};
          }
        }

        if (coord == end) break;
        if (blockSequences.getBlockStrand(blockId)) {
          globalCoords.stepRightCoordinate(coord);
        } else {
          globalCoords.stepLeftCoordinate(coord);
        }
      }

      if (curNucRange.first != -1) {
        nucRanges.push_back(curNucRange);
      }


      if (blockSequences.blockStrand[blockId]) {
        for (const auto& range : nucRanges) {
          updateGapMapStep(gapMap, {false, range}, gapRunBacktracks, gapMapUpdates);
        }
      } else {
        std::vector<std::pair<int64_t, int64_t>> invertedRanges = invertRanges(nucRanges, blockRanges[blockId]);
        for (const auto& range : invertedRanges) {
          updateGapMapStep(gapMap, {false, range}, gapRunBacktracks, gapMapUpdates);
        }
      }
    } else if (oldExists && newExists && (blockSequences.blockStrand[blockId] != oldBlockStrand[blockId])) {
      // on to on -> but strand flipped
      flippedSet.add(boost::icl::interval<int64_t>::closed(blockRanges[blockId].first, blockRanges[blockId].second));
      std::vector<int64_t> oldCurBaseCounts = curBaseCounts;
      panmapUtils::Coordinate coord = globalCoords.getBlockStartCoord(blockId);
      panmapUtils::Coordinate end = globalCoords.getBlockEndCoord(blockId);

      while (true) {
        char c = blockSequences.getSequenceBase(coord);
        c = c == 'x' ? '-' : c;
        if (!blockSequences.getBlockStrand(blockId)) c = panmanUtils::getComplementCharacter(c);

        if (getIndexFromNucleotide(c) <= 3) {
          ++curBaseCounts[getIndexFromNucleotide(c)];
          --parentBaseCountsBacktrack[getIndexFromNucleotide(c)];
        }

        if (coord == end) break;
        globalCoords.stepRightCoordinate(coord);

      }
    }
  }

  for (const auto& blockId : invertedBlocks) {
    invertGapMap(gapMap, blockRanges[blockId], gapRunBlocksBacktracks, gapMapUpdates);
  }

  if (node->parent != nullptr) {
    // make nuc run interval set from gap map
    boost::icl::interval_set<int64_t> curNucRunSet = gapMapToNucRunSet(gapMap, blockRanges);

    // remove flipped blocks from nucRunSet so they are not counted as indels
    curNucRunSet -= flippedSet;
    parNucRunSet -= flippedSet;

    // XOR ranges between cur and par
    boost::icl::interval_set<int64_t> xorSet = curNucRunSet ^ parNucRunSet;

    // AND ranges bewteen xor and par -> deletion
    boost::icl::interval_set<int64_t> deletionSet = xorSet & parNucRunSet;

    // AND ranges between xor and cur -> insertion
    boost::icl::interval_set<int64_t> insertionSet = xorSet & curNucRunSet;

    
    if (isMutatedNuc || !insertionSet.empty() || !deletionSet.empty()) {
      for (const auto& interval : insertionSet) {
        int64_t size = boost::icl::last(interval) - boost::icl::first(interval) + 1;
        if (insCount.find(size) == insCount.end()) {
          insCount[size] = 0;
        }
        ++insCount[size];
      }

      for (const auto& interval : deletionSet) {
        int64_t size = boost::icl::last(interval) - boost::icl::first(interval) + 1;
        if (delCount.find(size) == delCount.end()) {
          delCount[size] = 0;
        }
        ++delCount[size];
      }

      for (size_t i = 0; i < 4; ++i) {
        totalBaseCounts[i] += curBaseCounts[i];
      }
    }
  }

  parentBaseCounts = curBaseCounts;

  for (auto it = gapRunBlocksBacktracks.rbegin(); it != gapRunBlocksBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  gapRunUpdates.clear();

  // checking nucs with brute force
  // {
  //   const auto& oldSequence = oldBlockSequences.sequence;
  //   const auto& curSequence = blockSequences.sequence;
  //   std::vector<std::vector<int64_t>> bruteCurNodeSubCount(4, std::vector<int64_t>(4, 0));
  //   std::vector<int64_t> bruteCurNodeTotalBaseCounts(4, 0);
  //   std::vector<std::pair<int64_t, int64_t>> bruteCurNodeIns;
  //   std::vector<std::pair<int64_t, int64_t>> bruteCurNodeDel;
  //   for (size_t blockId = 0; blockId < oldSequence.size(); blockId++) {
  //     if (!(oldBlockSequences.blockExists[blockId] || blockSequences.blockExists[blockId])) {
  //       continue;
  //     }
  //     if (blockSequences.blockExists[blockId]) { 
  //       for (size_t j = 0; j < curSequence[blockId].size(); j++) {
  //         // process nuc gap positions
  //         for (size_t k = 0; k < curSequence[blockId][j].second.size(); k++) {
  //           char curNuc = curSequence[blockId][j].second[k];
  //           if (getIndexFromNucleotide(curNuc) <= 3) {
  //             if (blockSequences.blockStrand[blockId]) {
  //               ++bruteCurNodeTotalBaseCounts[getIndexFromNucleotide(curNuc)];
  //             } else {
  //               ++bruteCurNodeTotalBaseCounts[getIndexFromNucleotide(panmanUtils::getComplementCharacter(curNuc))];
  //             }
  //           }
  //         }

  //         // process main nuc position
  //         char curNuc = curSequence[blockId][j].first;
  //         if (getIndexFromNucleotide(curNuc) <= 3) {
  //           if (blockSequences.blockStrand[blockId]) {
  //             ++bruteCurNodeTotalBaseCounts[getIndexFromNucleotide(curNuc)];
  //           } else {
  //             ++bruteCurNodeTotalBaseCounts[getIndexFromNucleotide(panmanUtils::getComplementCharacter(curNuc))];
  //           }
  //         }
  //       }
  //     }
  //     if (oldBlockSequences.blockExists[blockId] && blockSequences.blockExists[blockId] && (oldBlockSequences.blockStrand[blockId] == blockSequences.blockStrand[blockId])) {
  //       for (size_t j = 0; j < oldSequence[blockId].size(); j++) {
  //         // process nuc gap positions
  //         for (size_t k = 0; k < oldSequence[blockId][j].second.size(); k++) {
  //           char oldNuc = oldBlockSequences.blockStrand[blockId] ? oldSequence[blockId][j].second[k] : panmanUtils::getComplementCharacter(oldSequence[blockId][j].second[k]);
  //           char curNuc = blockSequences.blockStrand[blockId] ? curSequence[blockId][j].second[k] : panmanUtils::getComplementCharacter(curSequence[blockId][j].second[k]);
  //           if (oldNuc != curNuc && getIndexFromNucleotide(oldNuc) <= 3 && getIndexFromNucleotide(curNuc) <= 3) {
  //             ++bruteCurNodeSubCount[getIndexFromNucleotide(oldNuc)][getIndexFromNucleotide(curNuc)];
  //           }
  //         }

  //         // process main nuc position
  //         char oldNuc = oldBlockSequences.blockStrand[blockId] ? oldSequence[blockId][j].first : panmanUtils::getComplementCharacter(oldSequence[blockId][j].first);
  //         char curNuc = blockSequences.blockStrand[blockId] ? curSequence[blockId][j].first : panmanUtils::getComplementCharacter(curSequence[blockId][j].first);
  //         if (oldNuc != curNuc && getIndexFromNucleotide(oldNuc) <= 3 && getIndexFromNucleotide(curNuc) <= 3) {
  //           ++bruteCurNodeSubCount[getIndexFromNucleotide(oldNuc)][getIndexFromNucleotide(curNuc)];
  //         }
  //       }
  //     }
  //   }

  //   std::vector<std::vector<int64_t>> dynamicCurNodeSubCount(4, std::vector<int64_t>(4, 0));
  //   for (size_t i = 0; i < 4; i++) {
  //     for (size_t j = 0; j < 4; j++) {
  //       dynamicCurNodeSubCount[i][j] = subCount[i][j] - oldSubCount[i][j];
  //     }
  //   }

  //   for (size_t i = 0; i < 4; i++) {
  //     if (bruteCurNodeTotalBaseCounts != curBaseCounts) {
  //       logging::err("At node {}, brute force total base counts do not match dynamic total base counts for index {}, brute: {}, dynamic: {}", node->identifier, i, bruteCurNodeTotalBaseCounts[i], curBaseCounts[i]);
  //       std::exit(1);
  //     } else {
  //       logging::info("At node {}, brute force total base counts match dynamic total base counts for index {}, brute: {}, dynamic: {}", node->identifier, i, bruteCurNodeTotalBaseCounts[i], curBaseCounts[i]);
  //     }
  //   }

  //   for (size_t i = 0; i < 4; i++) {
  //     for (size_t j = 0; j < 4; j++) {
  //       if (bruteCurNodeSubCount[i][j] != dynamicCurNodeSubCount[i][j]) {
  //         logging::err("At node {}, brute force sub count do not match dynamic sub count for index {}, {}, brute: {}, dynamic: {}", node->identifier, i, j, bruteCurNodeSubCount[i][j], dynamicCurNodeSubCount[i][j]);
  //         std::exit(1);
  //       } else {
  //         logging::info("At node {}, brute force sub count match dynamic sub count for index {}, {}, brute: {}, dynamic: {}", node->identifier, i, j, bruteCurNodeSubCount[i][j], dynamicCurNodeSubCount[i][j]);
  //       }
  //     }
  //   }
  //   logging::info("Brute force and dynamic nuc counts match for node {}", node->identifier);
  // }

  // update oldBlockExists and oldBlockStrand
  for (const auto& [blockId, oldExists, oldStrand, newExists, newStrand] : blockMutationRecord) {
    oldBlockExists[blockId] = newExists;
    oldBlockStrand[blockId] = newStrand;
  }

  for (panmanUtils::Node *child : node->children) {
    buildMutationMatricesHelper_test(mutMat, tree, child, blockSequences, oldBlockExists, oldBlockStrand, gapMap, invertedBlocks, parentBaseCounts, totalBaseCounts, subCount, insCount, delCount, globalCoords, blockRanges);
  }

  // undo gapMap updates
  for (auto it = gapRunBacktracks.rbegin(); it != gapRunBacktracks.rend(); ++it) {
    const auto& [del, range] = *it;
    if (del) {
      gapMap.erase(range.first);
    } else {
      gapMap[range.first] = range.second;
    }
  }

  // undo base counts
  for (int i = 0; i < 4; i++) {
    parentBaseCounts[i] += parentBaseCountsBacktrack[i];
  }

  undoMutations(blockSequences, oldBlockExists, oldBlockStrand, invertedBlocks, blockMutationRecord, nucMutationRecord, invertedBlocksBacktracks);
}

inline std::vector<std::vector<double>>
phredMatrix2ProbMatrix(const std::vector<std::vector<double>> &phredMatrix) {
  std::vector<std::vector<double>> prob_matrix = phredMatrix;
  for (int i = 0; i < prob_matrix.size(); i++) {
    for (int j = 0; j < prob_matrix[i].size(); j++) {
      prob_matrix[i][j] = pow(10, -prob_matrix[i][j] / 10);
    }
  }
  return prob_matrix;
}

inline std::vector<std::vector<double>>
probMatrix2PhredMatrix(const std::vector<std::vector<double>> &probMatrix) {
  std::vector<std::vector<double>> phred_matrix = probMatrix;
  for (int i = 0; i < phred_matrix.size(); i++) {
    for (int j = 0; j < phred_matrix[i].size(); j++) {
      phred_matrix[i][j] = -10 * log10(phred_matrix[i][j]);
    }
  }
  return phred_matrix;
}

inline double
getAverageMutationRate(const std::vector<std::vector<double>> &matrix) {
  if (matrix.empty() || matrix.size() != matrix[0].size())
    throw std::invalid_argument("Matrix must be square and non-empty.");

  double sum = 0.0;
  size_t count = 0;
  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix.size(); ++j) {
      if (i != j) {
        sum += matrix[i][j];
        ++count;
      }
    }
  }

  if (count == 0)
    throw std::runtime_error("No off-diagonal elements to calculate average.");
  return sum / count;
}

vector<std::vector<double>>
genotyping::scaleMutationSpectrum(const mutationMatrices &mutMat,
                                  double mutationRate) {
  vector<std::vector<double>> scaled_submat_phred = mutMat.submat;
  vector<std::vector<double>> scaled_submat_prob =
      phredMatrix2ProbMatrix(scaled_submat_phred);
  double avg_mutation_rate = getAverageMutationRate(scaled_submat_prob);
  double scale_factor = mutationRate / avg_mutation_rate;

  std::vector<double> scaled_submat_prob_row_sums(scaled_submat_prob.size());
  for (int i = 0; i < scaled_submat_prob.size(); i++) {
    for (int j = 0; j < scaled_submat_prob[i].size(); j++) {
      if (i == j)
        continue;
      scaled_submat_prob[i][j] *= scale_factor;
      scaled_submat_prob_row_sums[i] += scaled_submat_prob[i][j];
    }
  }

  for (int i = 0; i < scaled_submat_prob.size(); i++) {
    scaled_submat_prob[i][i] = 1 - scaled_submat_prob_row_sums[i];
  }

  return probMatrix2PhredMatrix(scaled_submat_prob);
}

std::vector<char> parse_alts(const std::string &alts_str) {
  std::vector<char> alts;
  std::istringstream alts_stream(alts_str);
  std::string alt;
  while (std::getline(alts_stream, alt, ',')) {
    if (alt.size() > 1) {
      throw std::runtime_error("Error: alt allel parsing error..");
    }
    alts.push_back(alt[0]);
  }
  return alts;
}

std::tuple<int, std::vector<double>, std::vector<int>, std::string, std::string>
parse_sample_formats(const std::string &sample_formats_str) {
  std::vector<std::string> sample_formats;
  std::istringstream sample_formats_stream(sample_formats_str);
  std::string sample_format;
  while (std::getline(sample_formats_stream, sample_format, ':')) {
    sample_formats.push_back(sample_format);
  }

  int gt = std::stoi(sample_formats[0]);
  std::vector<double> pls;
  std::istringstream pls_stream(sample_formats[1]);
  std::string pl;
  while (std::getline(pls_stream, pl, ',')) {
    pls.push_back(std::stod(pl));
  }

  std::vector<int> ads;
  std::istringstream ads_stream(sample_formats[2]);
  std::string ad;
  while (std::getline(ads_stream, ad, ',')) {
    ads.push_back(std::stoi(ad));
  }

  return std::make_tuple(gt, pls, ads, sample_formats[1], sample_formats[2]);
}

std::string genotyping::applyMutationSpectrum(
    const std::string &line,
    const std::vector<std::vector<double>> &scaled_submat) {
  std::vector<std::string> fields;
  std::istringstream line_stream(line);
  std::string field;
  while (std::getline(line_stream, field, '\t')) {
    fields.push_back(field);
  }

  if (fields.size() < 10 || fields[0] == "#CHROM") {
    return line;
  }

  if (fields.size() != 10) {
    throw std::runtime_error(
        "Couldn't parse VCF. Unrecognized number of fields.");
  }

  if (fields[4] == ".") {
    return "";
  } else if (fields[7].substr(0, 2) != "DP") {
    if (fields[9][0] == '0')
      return "";
    else
      return line;
  } else if (getIndexFromNucleotide(fields[3][0]) > 3) {
    if (fields[9][0] == '0')
      return "";
    else
      return line;
  }

  if (fields[3].size() > 1) {
    throw std::runtime_error("Error: reference allele parsing error.");
  }

  int ref_nuc_idx = getIndexFromNucleotide(fields[3][0]);
  std::vector<char> alts = parse_alts(fields[4]);

  auto [gt, pls, ads, pls_string, ads_string] = parse_sample_formats(fields[9]);

  std::vector<double> gls;
  if (alts.size() + 1 == pls.size()) {
    gls = pls;
  } else {
    for (int i = 0; i < pls.size(); i += 2) {
      gls.push_back(pls[i]);
      if (gls.size() == alts.size() + 1) {
        break;
      }
    }
  }

  gls[0] += scaled_submat[ref_nuc_idx][ref_nuc_idx];
  for (int i = 1; i < gls.size(); i++) {
    gls[i] = gls[i] +
             scaled_submat[ref_nuc_idx][getIndexFromNucleotide(alts[i - 1])];
  }

  double min_gl = *std::min_element(gls.begin(), gls.end());
  int min_gl_index = 0;
  for (int i = 0; i < gls.size(); i++) {
    gls[i] -= min_gl;
    if (gls[i] == 0)
      min_gl_index = i;
  }

  if (min_gl_index == 0)
    return "";

  gt = min_gl_index;
  double qual = gls[0];

  std::stringstream ss;
  ss << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3]
     << "\t" << fields[4] << "\t" << std::fixed << std::setprecision(4) << qual
     << "\t" << fields[6] << "\t" << fields[7] << "\t" << fields[8] << "\t"
     << gt << ":" << pls_string << ":" << ads_string;
  return ss.str();
}
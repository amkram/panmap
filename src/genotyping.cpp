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
    state::StateManager &stateManager,
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
    buildMutationMatricesHelper(mutMat, tree, child, stateManager, parentBaseCounts, 
                           totalBaseCounts, subCount, insCount, delCount);
  }
}

// void genotyping::fillMutationMatricesFromTree_test(
//     mutationMatrices &mutMat, 
//     panmanUtils::Tree *tree, 
//     const std::string& path) {
  
//   if (!tree || !tree->root) {
//     throw std::invalid_argument("Invalid tree or root node");
//   }
  
//   logging::info("Building mutation matrices from tree...");
  
//   // Initialize state manager for coordinate access
//   auto stateManager = std::make_unique<state::StateManager>();
  
//   // Initialize base count vectors and mutation counters
//   std::vector<int64_t> parentBaseCounts(4, 0);
//   std::vector<int64_t> totalBaseCounts(4, 0);
//   std::vector<std::vector<int64_t>> subCount(4, std::vector<int64_t>(4, 0));
//   std::unordered_map<int64_t, int64_t> insCount;
//   std::unordered_map<int64_t, int64_t> delCount;
  
//   // Process the tree to build mutation statistics
//   // buildMutationMatricesHelper(mutMat, tree, tree->root, *stateManager, 
//   //                        parentBaseCounts, totalBaseCounts, subCount, insCount, delCount);
  
//   // Calculate totals
//   int64_t totalNucCounts = 0;
//   int64_t totalInsCounts = 0;
//   int64_t totalDelCounts = 0;
  
//   for (const auto& count : totalBaseCounts) totalNucCounts += count;
//   for (const auto& [size, count] : insCount) totalInsCounts += count;
//   for (const auto& [size, count] : delCount) totalDelCounts += count;
  
//   // Add no-mutation counts
//   insCount[0] = totalNucCounts - totalInsCounts;
//   delCount[0] = totalNucCounts - totalDelCounts;
  
//   // Initialize substitution matrix with base counts
//   for (int i = 0; i < 4; ++i) {
//     mutMat.submat[i][i] = static_cast<double>(totalBaseCounts[i]);
//   }
  
//   // Update substitution matrix with counts
//   for (int i = 0; i < 4; ++i) {
//     for (int j = 0; j < 4; ++j) {
//       if (i != j) {
//         mutMat.submat[i][j] = static_cast<double>(subCount[i][j]);
//         mutMat.submat[i][i] -= mutMat.submat[i][j];
//       }
//     }
//   }
  
//   // Store insertion and deletion counts
//   for (const auto& [size, count] : insCount) {
//     mutMat.insmat[size] = static_cast<double>(count);
//   }
  
//   for (const auto& [size, count] : delCount) {
//     mutMat.delmat[size] = static_cast<double>(count);
//   }
  
//   // Convert counts to log probabilities
  
//   // Insertions
//   for (auto [size, count] : mutMat.insmat) {
//     mutMat.insmat[size] = -10 * log10(count / static_cast<double>(totalNucCounts));
//   }
  
//   // Deletions
//   for (auto [size, count] : mutMat.delmat) {
//     mutMat.delmat[size] = -10 * log10(count / static_cast<double>(totalNucCounts));
//   }
  
//   // Calculate maximum log probabilities for insertions and deletions
//   mutMat.maxInsLogProb = 100.0; // Default high penalty
//   mutMat.maxDelLogProb = 100.0; // Default high penalty
  
//   if (!mutMat.insmat.empty()) {
//     mutMat.maxInsLogProb = std::max_element(
//       mutMat.insmat.begin(), mutMat.insmat.end(),
//       [](const auto& a, const auto& b) { return a.second < b.second; }
//     )->second;
//   }
  
//   if (!mutMat.delmat.empty()) {
//     mutMat.maxDelLogProb = std::max_element(
//       mutMat.delmat.begin(), mutMat.delmat.end(),
//       [](const auto& a, const auto& b) { return a.second < b.second; }
//     )->second;
//   }
  
//   // Substitutions
//   for (auto i = 0; i < 4; i++) {
//     for (auto j = 0; j < 4; j++) {
//       if (totalBaseCounts[i] > 0) {
//         mutMat.submat[i][j] = -10 * log10(mutMat.submat[i][j] / static_cast<double>(totalBaseCounts[i]));
//       } else {
//         mutMat.submat[i][j] = 100.0; // High penalty for unseen bases
//       }
//     }
//   }
  
//   mutMat.filled = true;
  
//   // Write matrices to file
//   std::ofstream outFile(path);
//   if (!outFile.is_open()) {
//     std::stringstream ss;
//     ss << "Failed to open file for writing mutation matrices: " << path
//        << ". Please check file permissions and path validity.";
//     throw std::runtime_error(ss.str());
//   }
  
//   // Write substitution matrix
//   for (auto i = 0; i < 4; i++) {
//     for (auto j = 0; j < 4; j++) {
//       outFile << mutMat.submat[i][j] << " ";
//     }
//     outFile << "\n";
//   }
  
//   // Write insertion matrix
//   std::vector<std::pair<int64_t, double>> insmat_sorted(mutMat.insmat.begin(), mutMat.insmat.end());
//   std::sort(insmat_sorted.begin(), insmat_sorted.end(), 
//             [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
//               return a.first < b.first;
//             });
  
//   for (const auto& [size, logProb] : insmat_sorted) {
//     outFile << size << ":" << logProb << " ";
//   }
//   outFile << "\n";
  
//   // Write deletion matrix
//   std::vector<std::pair<int64_t, double>> delmat_sorted(mutMat.delmat.begin(), mutMat.delmat.end());
//   std::sort(delmat_sorted.begin(), delmat_sorted.end(), 
//             [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
//               return a.first < b.first;
//             });
  
//   for (const auto& [size, logProb] : delmat_sorted) {
//     outFile << size << ":" << logProb << " ";
//   }
//   outFile << "\n";
  
//   outFile.close();
//   logging::info("Mutation matrices saved to {}", path);
// }

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
static void initializeSequence(
  panmanUtils::Tree *tree,
  std::vector<std::vector<std::pair<char, std::vector<char>>>> &sequence,
  std::vector<bool> &blockExists,
  std::vector<bool> &blockStrand
) {
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
    std::vector<bool> &oldBlockExists,
    std::vector<bool> &oldBlockStrand,
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
  std::vector<bool>& blockExists = blockSequences.blockExists;
  std::vector<bool>& blockStrand = blockSequences.blockStrand;

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
  std::vector<bool> &oldBlockExists,
  std::vector<bool> &oldBlockStrand,
  std::unordered_set<int64_t> &invertedBlocks,
  const std::vector<std::tuple<uint32_t, bool, bool, bool, bool>> &blockMutationRecord,
  const std::vector<std::tuple<panmapUtils::Coordinate, char, char>> &nucMutationRecord,
  const std::vector<std::pair<int64_t, bool>> &invertedBlocksBacktracks
) {
  std::vector<bool>& blockExists = blockSequences.blockExists;
  std::vector<bool>& blockStrand = blockSequences.blockStrand;

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
  std::vector<bool> &oldBlockExists,
  std::vector<bool> &oldBlockStrand,
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

void genotyping::fillMutationMatricesFromTree_test(
  mutationMatrices &mutMat, 
  panmanUtils::Tree *tree, 
  const std::string& path
) {
  
  if (!tree || !tree->root) {
    throw std::invalid_argument("Invalid tree or root node");
  }
  
  logging::info("Building mutation matrices from tree...");

  panmapUtils::BlockSequences blockSequences(tree);
  std::vector<bool> oldBlockExists(blockSequences.numBlocks(), false);
  std::vector<bool> oldBlockStrand(blockSequences.numBlocks(), true);

  // initialize globalCoords to converted 3d to 1d coordinates -> for general coordinate conversion
  panmapUtils::GlobalCoords globalCoords(blockSequences);

  // initialize gapMap to be all gaps -> for counting indels
  std::map<int64_t, int64_t> gapMap{{0, globalCoords.lastScalarCoord}};
  std::unordered_set<int64_t> invertedBlocks;

  // block ranges -> for counting indels
  std::vector<std::pair<int64_t, int64_t>> blockRanges;
  for (int64_t i = 0; i < blockSequences.numBlocks(); i++) {
    int64_t start = globalCoords.getBlockStartScalar(i);
    int64_t end = globalCoords.getBlockEndScalar(i);
    // sanity check
    if (i > 1 && start != blockRanges.back().second + 1) {
      logging::err("blockRanges not continuous");
      std::exit(1);
    }
    blockRanges.push_back({start, end});
  }

  std::vector<int64_t> parentBaseCounts(4);
  std::vector<int64_t> totalBaseCounts(4);
  std::vector<std::vector<int64_t>> subCount(4, std::vector<int64_t>(4, 0));
  std::unordered_map<int64_t, int64_t> insCount;
  std::unordered_map<int64_t, int64_t> delCount;

  logging::info("Initializing complete... starting to fill mutation matrices");
  buildMutationMatricesHelper_test(mutMat, tree, tree->root, blockSequences, oldBlockExists, oldBlockStrand, gapMap, invertedBlocks, parentBaseCounts, totalBaseCounts, subCount, insCount, delCount, globalCoords, blockRanges);

  logging::info("Mutation matrices successfully filled from tree");

  int64_t totalNucCounts = 0;
  int64_t totalInsCounts = 0;
  int64_t totalDelCounts = 0;
  for (const auto& count : totalBaseCounts) totalNucCounts += count;
  for (const auto& [size, count] : insCount) totalInsCounts += count;
  for (const auto& [size, count] : delCount) totalDelCounts += count;
  insCount[0] = totalNucCounts - totalInsCounts;
  delCount[0] = totalNucCounts - totalDelCounts;

  for (int i = 0; i < 4; ++i) {
    mutMat.submat[i][i] = static_cast<double>(totalBaseCounts[i]);
  }
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (i != j) {
        mutMat.submat[i][j] = static_cast<double>(subCount[i][j]);
        mutMat.submat[i][i] -= mutMat.submat[i][j];
      }
    }
  }

  for (const auto& [size, count] : insCount) {
    mutMat.insmat[size] = static_cast<double>(count);
  }
  for (const auto& [size, count] : delCount) {
    mutMat.delmat[size] = static_cast<double>(count);
  }

  // insertion
  for (auto [size, count] : mutMat.insmat) {
    mutMat.insmat[size] = -10 * log10f(count / static_cast<double>(totalNucCounts));
  }
  // deletion
  for (auto [size, count] : mutMat.delmat) {
    mutMat.delmat[size] = -10 * log10f(count / static_cast<double>(totalNucCounts));
  }

  // substitution
  for (auto i = 0; i < 4; i++) {
    for (auto j = 0; j < 4; j++) {
      mutMat.submat[i][j] = -10 * log10f(mutMat.submat[i][j] / static_cast<double>(totalBaseCounts[i]));
    }
  }
  mutMat.filled = true;

  std::ofstream outFile(path);
  for (auto i = 0; i < 4; i++) {
    for (auto j = 0; j < 4; j++) {
      outFile << mutMat.submat[i][j] << " ";
    }
    outFile << "\n";
  }

  std::vector<std::pair<int64_t, double>> insmat_sorted(mutMat.insmat.begin(), mutMat.insmat.end());
  std::sort(insmat_sorted.begin(), insmat_sorted.end(), [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
    return a.first < b.first;
  });
  for (const auto& [size, logProb] : insmat_sorted) {
    outFile << size << ":" << logProb << " ";
  }
  outFile << "\n";

  std::vector<std::pair<int64_t, double>> delmat_sorted(mutMat.delmat.begin(), mutMat.delmat.end());
  std::sort(delmat_sorted.begin(), delmat_sorted.end(), [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
    return a.first < b.first;
  });
  for (const auto& [size, logProb] : delmat_sorted) {
    outFile << size << ":" << logProb << " ";
  }
  outFile << "\n";
  outFile.close();
  logging::info("Mutation matrices successfully written to file: {}", path);
}

enum variationType { SNP = 1, INS = 2, DEL = 4 };

double phred_complement(double q) {
  double p = pow(10, (-q / 10));
  return -10 * log10(1 - p);
}

double third_phred(double q) {
  double p = pow(10, (-q / 10)) / 3.0;
  return -10 * log10(p);
}

void to_upper(string &str) {
  for (char &c : str) {
    c = toupper(c);
  }
}

static char getNucleotideFromIndex(int index) {
  switch (index) {
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  case 4:
    return '*';
  default:
    return 'N';
  }
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
  int min_gl_index;
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

double likelihood(int genotype_idx, const vector<vector<double>> &read_errs,
                  const map<string, vector<double>> &deletions,
                  const map<string, vector<double>> &insertions,
                  int variation_type) {
  vector<double> genotype_probs;
  vector<double> variants_probs;

  for (int i = 0; i < read_errs.size(); i++) {
    const auto &row = read_errs[i];
    if ((variation_type & variationType::SNP) && (genotype_idx == i)) {
      for (const auto &prob : row) {
        genotype_probs.push_back(phred_complement(prob));
      }
    } else {
      variants_probs.insert(variants_probs.end(), row.begin(), row.end());
      // for (const auto& prob : row) {
      //   variants_probs.push_back(third_phred(prob));
      // }
    }
  }

  int ins_i = 0;
  for (const auto &insertion : insertions) {
    if ((variation_type & variationType::INS) && (genotype_idx == ins_i)) {
      for (const auto &prob : insertion.second) {
        genotype_probs.push_back(phred_complement(prob));
      }
    } else {
      variants_probs.insert(variants_probs.end(), insertion.second.begin(),
                            insertion.second.end());
    }
    ins_i++;
  }

  int del_i = 0;
  for (const auto &deletion : deletions) {
    if ((variation_type & variationType::DEL) && (genotype_idx == del_i)) {
      for (const auto &prob : deletion.second) {
        genotype_probs.push_back(phred_complement(prob));
      }
    } else {
      variants_probs.insert(variants_probs.end(), deletion.second.begin(),
                            deletion.second.end());
    }
    del_i++;
  }

  double genotype_prob =
      accumulate(genotype_probs.begin(), genotype_probs.end(), 0.0);
  double variant_prob =
      accumulate(variants_probs.begin(), variants_probs.end(), 0.0);

  return genotype_prob + variant_prob;
}

vector<double>
genotype_likelihoods(const vector<vector<double>> &read_errs,
                     const map<string, vector<double>> &deletions,
                     const map<string, vector<double>> &insertions,
                     const int8_t &site_info, const char &ref_nuc) {
  vector<double> likelihoods;
  likelihoods.resize(5);
  for (auto i = 0; i < 5; ++i) {
    likelihoods[i] = numeric_limits<double>::max();
  }

  auto variation_types = site_info & 7;
  auto ref_nuc_idx = site_info >> 3;

  for (int i = 0; i < read_errs.size(); ++i) {
    if (read_errs[i].empty() && i != ref_nuc_idx) {
      continue;
    }
    likelihoods[i] =
        likelihood(i, read_errs, deletions, insertions, variationType::SNP);
  }

  if (variation_types & variationType::INS) {
    for (auto i = 0; i < insertions.size(); i++) {
      likelihoods.push_back(
          likelihood(i, read_errs, deletions, insertions, variationType::INS));
    }
  }

  if (variation_types & variationType::DEL) {
    for (auto i = 0; i < deletions.size(); i++) {
      likelihoods.push_back(
          likelihood(i, read_errs, deletions, insertions, variationType::DEL));
    }
  }

  return likelihoods;
}

vector<double>
genotype_posteriors(const vector<double> &likelihoods,
                    const map<string, vector<double>> &deletions,
                    const map<string, vector<double>> &insertions,
                    const int8_t &site_info, const mutationMatrices &mutmat) {
  vector<double> posteriors;
  posteriors.resize(4);
  for (auto i = 0; i < 4; i++) {
    posteriors[i] = numeric_limits<double>::max();
  }

  auto ref_nuc = site_info >> 3;
  for (auto i = 0; i < 4; ++i) {
    if (likelihoods[i] != numeric_limits<double>::max()) {
      if (ref_nuc == 5) {
        posteriors[i] = likelihoods[i];
      } else {
        posteriors[i] = likelihoods[i] + mutmat.submat[ref_nuc][i];
      }
    }
  }

  size_t insertion_idx = 0;
  for (const auto &insertion : insertions) {
    int64_t insSize = insertion.first.size();

    if (mutmat.insmat.find(insSize) != mutmat.insmat.end()) {
      if (mutmat.insmat.at(insSize) > mutmat.maxInsLogProb) {
        posteriors.push_back(likelihoods[5 + insertion_idx] +
                             mutmat.maxInsLogProb);
      } else {
        posteriors.push_back(likelihoods[5 + insertion_idx] +
                             mutmat.insmat.at(insSize));
      }
    } else {
      posteriors.push_back(likelihoods[5 + insertion_idx] +
                           mutmat.maxInsLogProb);
    }

    insertion_idx++;
  }

  size_t deletion_idx = 0;
  for (const auto &deletion : deletions) {
    int64_t delSize = deletion.first.size();

    if (mutmat.delmat.find(delSize) != mutmat.delmat.end()) {
      if (mutmat.delmat.at(delSize) > mutmat.maxDelLogProb) {
        posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] +
                             mutmat.maxDelLogProb);
      } else {
        posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] +
                             mutmat.delmat.at(delSize));
      }
    } else {
      posteriors.push_back(likelihoods[5 + insertion_idx + deletion_idx] +
                           mutmat.maxDelLogProb);
    }

    deletion_idx++;
  }

  double min_score = *min_element(posteriors.begin(), posteriors.end());
  for (int i = 0; i < posteriors.size(); ++i) {
    posteriors[i] -= min_score;
  }
  return posteriors;
}

genotyping::VariationSite::VariationSite(size_t sid, char ref, size_t position,
                                         int variation_types,
                                         const string &nucs,
                                         const vector<string> &insertion_seqs,
                                         const vector<string> &deletion_seqs,
                                         const string &errors,
                                         const mutationMatrices &mutMat) {
  this->site_id = sid;
  this->ref_position = position;
  size_t offset = 0;
  this->site_info = (getIndexFromNucleotide(ref) << 3) + variation_types;
  this->ref_nuc = ref;
  this->read_errs.resize(5);
  assert(errors.size() ==
         nucs.size() + insertion_seqs.size() + deletion_seqs.size());

  for (auto i = 0; i < nucs.size(); ++i) {
    if (getIndexFromNucleotide(nucs[i]) == 5) {
      this->read_errs[4].push_back(double(errors[i]) - 33.0);
    } else {
      this->read_errs[getIndexFromNucleotide(nucs[i])].push_back(
          double(errors[i]) - 33.0);
    }
  }
  offset += nucs.size();

  if (variation_types & variationType::INS) {
    for (auto i = 0; i < insertion_seqs.size(); ++i) {
      this->insertions[insertion_seqs[i]].push_back(
          double(errors[i + offset] - 33.0));
    }
    offset += insertion_seqs.size();
  }

  if (variation_types & variationType::DEL) {
    for (auto i = 0; i < deletion_seqs.size(); i++) {
      this->deletions[deletion_seqs[i]].push_back(
          double(errors[i + offset] - 33.0));
    }
  }

  this->likelihoods =
      genotype_likelihoods(this->read_errs, this->deletions, this->insertions,
                           this->site_info, this->ref_nuc);
  this->posteriors =
      genotype_posteriors(this->likelihoods, this->deletions, this->insertions,
                          this->site_info, mutMat);

  for (size_t i = 0; i < 4; i++) {
    read_depth.push_back(this->read_errs[i].size());
  }
  for (const auto &ins : this->insertions) {
    read_depth.push_back(ins.second.size());
  }
  for (const auto &del : this->deletions) {
    read_depth.push_back(del.second.size());
  }

  for (size_t i = 0; i < this->posteriors.size(); i++) {
    if (this->posteriors[i] == 0.0) {
      this->most_probable_idx = i;
      break;
    }
  }
}

int parse_readbases(string readbase_string, const string &readbase_errors,
                    char ref_nuc, string &nucs, vector<string> &insertion_seqs,
                    vector<string> &deletion_seqs, string &errs) {
  int variation_types = 0;

  regex extra_regex("\\^.{1}|\\$");
  readbase_string = regex_replace(readbase_string, extra_regex, "");

  string snp_errs, ins_errs, del_errs;
  size_t cur_start = 0, cur_idx = 0;

  string basePairs = string("ATCGatcg*");

  while (cur_start < readbase_string.size()) {

    bool is_last = (cur_start == readbase_string.size() - 1);

    if (basePairs.find(readbase_string[cur_start]) != std::string::npos) {
      // SNP
      nucs += toupper(readbase_string[cur_start]);
      variation_types |= variationType::SNP;
      snp_errs += readbase_errors[cur_idx];
      cur_start += 1;
      cur_idx++;
    } else if (readbase_string[cur_start] == '.' ||
               readbase_string[cur_start] == ',') {
      if (is_last) {
        // SNP
        nucs += ref_nuc;
        snp_errs += readbase_errors[cur_idx];
        cur_start += 1;
        cur_idx++;
      } else {
        if (readbase_string[cur_start + 1] == '-') {
          // DEL
          variation_types |= variationType::DEL;
          del_errs += readbase_errors[cur_idx];
          int indel_size = 0;

          cur_start += 2;
          while (std::isdigit((int)readbase_string[cur_start])) {

            indel_size *= 10;
            indel_size += (int)readbase_string[cur_start] - 48;
            cur_start++;
          }

          string seq = readbase_string.substr(cur_start, indel_size);
          to_upper(seq);
          deletion_seqs.push_back(seq);
          cur_start += indel_size;

          cur_idx++;
        } else if (readbase_string[cur_start + 1] == '+') {

          // INS
          variation_types |= variationType::INS;
          ins_errs += readbase_errors[cur_idx];
          int indel_size = 0;

          cur_start += 2;
          while (std::isdigit((int)readbase_string[cur_start])) {
            indel_size *= 10;
            indel_size += (int)readbase_string[cur_start] - 48;
            cur_start++;
          }

          string seq = readbase_string.substr(cur_start, indel_size);
          to_upper(seq);
          insertion_seqs.push_back(seq);

          cur_start += indel_size;

          cur_idx++;
        } else {
          // SNP
          nucs += ref_nuc;
          snp_errs += readbase_errors[cur_idx];
          cur_start += 1;
          cur_idx++;
        }
      }

    } else if (readbase_string[cur_start] ==
               '-') { // This means the same bp has an insertion and a deletion,
                      // for now we ignore the deletion TODO
      // SUBSTITION

      // variation_types |= variationType::DEL;
      // del_errs += readbase_errors[cur_idx];
      int indel_size = 0;

      cur_start++;
      while (std::isdigit((int)readbase_string[cur_start])) {

        indel_size *= 10;
        indel_size += (int)readbase_string[cur_start] - 48;
        cur_start++;
      }

      // string seq = readbase_string.substr(cur_start, indel_size);
      // to_upper(seq);
      // deletion_seqs.push_back(seq);
      cur_start += indel_size;

      // cur_idx++;
    } else {
      cur_start++;
    }
  }

  errs = snp_errs + ins_errs + del_errs;
  return variation_types;
}

pair<vector<VariationSite>, pair<size_t, size_t>>
genotyping::getVariantSites(std::istream &fin, const mutationMatrices &mutMat) {
  regex variant_pattern("[ACGTacgt\\*]+");
  vector<VariationSite> candidateVariants;
  pair<size_t, size_t> maskRange(numeric_limits<size_t>::max(), 0);
  size_t site_id = 0;
  string line;

  while (getline(fin, line)) {

    vector<string> fields;
    stringSplit(line, '\t', fields);
    string readbases_string = fields[4];
    string readbases_errors = fields[5];
    size_t coverage = stoul(fields[3]);
    size_t position = stoul(fields[1]) - 1;
    if (position < maskRange.first) {
      maskRange.first = position;
    }
    if (position > maskRange.second) {
      maskRange.second = position;
    }
    if ((coverage > 0) && regex_search(readbases_string, variant_pattern)) {
      char ref_nuc = fields[2][0];
      string errs;
      string nucs;
      vector<string> insertion_seqs;
      vector<string> deletion_seqs;

      int variation_types =
          parse_readbases(readbases_string, readbases_errors, ref_nuc, nucs,
                          insertion_seqs, deletion_seqs, errs);
      candidateVariants.emplace_back(
          VariationSite(site_id, ref_nuc, position, variation_types, nucs,
                        insertion_seqs, deletion_seqs, errs, mutMat));
      site_id++;
    }
  }

  return make_pair(candidateVariants, maskRange);
}

static double get_qual_for_ambiguous_ref(const VariationSite &site) {
  double qual = 0.0;
  for (int i = 0; i < 4; i++) {
    const auto &row = site.read_errs[i];
    qual += accumulate(row.begin(), row.end(), 0.0);
  }

  for (const auto &insertion : site.insertions) {
    qual += accumulate(insertion.second.begin(), insertion.second.end(), 0.0);
  }

  for (const auto &deletion : site.deletions) {
    qual += accumulate(deletion.second.begin(), deletion.second.end(), 0.0);
  }

  vector<double> filtered_likelihoods;
  for (size_t i = 0; i < site.likelihoods.size(); ++i) {
    if (i != 3 || site.likelihoods[i] != std::numeric_limits<double>::max()) {
      filtered_likelihoods.push_back(site.likelihoods[i]);
    }
  }

  return qual -
         *min_element(filtered_likelihoods.begin(), filtered_likelihoods.end());
}

static void printVCFLine(const VariationSite &site, std::ofstream &fout) {
  size_t position = site.ref_position + 1;
  int ref_nuc_idx = site.site_info >> 3;
  size_t readDepth = 0;
  vector<string> altAlleles;
  vector<size_t> ad;
  vector<int> pl;
  string refAllele;
  int gt;

  string quality;
  if (ref_nuc_idx == 5) {
    quality = to_string(int(get_qual_for_ambiguous_ref(site)));
  } else {
    quality = to_string(int(site.posteriors[ref_nuc_idx]));
  }

  refAllele += site.ref_nuc;

  // find longest deletion
  size_t ldl = 0; // longest deletion length
  string lds;     // longest deletion string
  for (const auto &del : site.deletions) {
    if (del.first.size() > ldl) {
      ldl = del.first.size();
      lds = del.first;
    }
  }

  if (!lds.empty()) {
    refAllele += lds;
  }

  // depth and likelihood for reference
  if (ref_nuc_idx == 5) {
    ad.push_back(0);
    pl.push_back(stoi(quality));
  } else {
    ad.push_back(site.read_depth[ref_nuc_idx]);
    pl.push_back(site.posteriors[ref_nuc_idx]);
    readDepth += site.read_depth[ref_nuc_idx];
  }
  if (pl.back() == 0.0) {
    gt = 0;
  }

  // substitutions
  for (int i = 0; i < 4; i++) {
    if (i == ref_nuc_idx) {
      continue;
    }
    if (site.posteriors[i] != numeric_limits<double>::max()) {
      altAlleles.push_back(getNucleotideFromIndex(i) + lds);
      ad.push_back(site.read_depth[i]);
      pl.push_back(site.posteriors[i]);
      readDepth += site.read_depth[i];
      if (pl.back() == 0.0) {
        gt = altAlleles.size();
      }
    }
  }

  // insertions
  size_t indelIdx = 0;
  for (const auto &ins : site.insertions) {
    altAlleles.push_back(site.ref_nuc + ins.first + lds);
    ad.push_back(site.read_depth[4 + indelIdx]);
    pl.push_back(site.posteriors[4 + indelIdx]);
    readDepth += site.read_depth[4 + indelIdx];
    if (pl.back() == 0.0) {
      gt = altAlleles.size();
    }
    indelIdx += 1;
  }

  // deletions
  for (const auto &del : site.deletions) {
    size_t delSize = del.first.size();
    if (delSize == ldl) {
      altAlleles.push_back(refAllele.substr(0, 1));
    } else {
      altAlleles.push_back(site.ref_nuc + lds.substr(delSize, ldl - delSize));
    }
    ad.push_back(site.read_depth[4 + indelIdx]);
    pl.push_back(site.posteriors[4 + indelIdx]);
    readDepth += site.read_depth[4 + indelIdx];
    if (pl.back() == 0.0) {
      gt = altAlleles.size();
    }
    indelIdx += 1;
  }

  fout << "ref"
       << "\t"             // #CHROM
       << position << "\t" // POS
       << "."
       << "\t"               // ID
       << refAllele << "\t"; // REF
  for (size_t i = 0; i < altAlleles.size() - 1; i++) {
    fout << altAlleles[i] << ",";
  }
  fout << altAlleles[altAlleles.size() - 1] << "\t"; // ALT
  fout << quality << "\t"                            // QUAL
       << "."
       << "\t"                       // FILTER
       << "DP=" << readDepth << "\t" // INFO
       << "GT:AD:GP"
       << "\t"       // FORMAT
       << gt << ":"; // SAMPLE
  for (size_t i = 0; i < ad.size() - 1; i++) {
    fout << ad[i] << ",";
  }
  fout << ad[ad.size() - 1] << ":";
  for (size_t i = 0; i < pl.size() - 1; i++) {
    if (pl[i] == -1) {
      fout << ".,";
      continue;
    }
    fout << pl[i] << ",";
  }
  fout << pl[pl.size() - 1] << endl;
}

void genotyping::printSamplePlacementVCF(std::istream &fin,
                                         const mutationMatrices &mutMat,
                                         bool keep_alts, size_t maskSize,
                                         std::ofstream &fout) {
  pair<vector<VariationSite>, pair<size_t, size_t>> variantSites =
      getVariantSites(fin, mutMat);
  const vector<VariationSite> candidateVariants = variantSites.first;
  if (candidateVariants.empty()) {
    return;
  }

  fout << "##fileformat=VCFv4.3\n"
       << "##contig=<ID=ref>\n"
       << "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read Depth>\n"
       << "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
       << "##FORMAT=<ID=AD,Number=1,Type=String,Description=Read depth for "
          "each allele>\n"
       << "##FORMAT=<ID=GP,Number=1,Type=String,Description=Genotype posterior "
          "probabilities in phred scale>\n";

  fout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

  vector<VariationSite> groupSites;

  for (const auto &curSite : candidateVariants) {
    bool skip = true;
    for (int i = 0; i < curSite.posteriors.size(); i++) {
      if ((curSite.site_info >> 3) == 5 ||
          (i != (curSite.site_info >> 3) &&
           curSite.posteriors[i] != numeric_limits<double>::max())) {
        skip = false;
        break;
      }
    }
    if (skip) {
      continue;
    }

    if (!keep_alts && curSite.most_probable_idx == (curSite.site_info >> 3)) {
      continue;
    }

    if (curSite.ref_position >= variantSites.second.first + maskSize &&
        curSite.ref_position <= variantSites.second.second - maskSize) {
      printVCFLine(curSite, fout);
    }
  }
}

void genotyping::genotype(std::string prefix, std::string refFileName,
                          std::string bestMatchSequence,
                          std::string bamFileName, std::string mpileupFileName,
                          std::string vcfFileName,
                          mutationMatrices &mutMat) {
  
  createMplpBcf(prefix, refFileName, bestMatchSequence, bamFileName,
                mpileupFileName);

  createVcfWithMutationMatrices(prefix, mpileupFileName, mutMat, vcfFileName,
                                0.0011);
}
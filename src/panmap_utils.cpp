
#include "panmap_utils.hpp"

namespace panmapUtils {

std::string seedChangeTypeToString(seedChangeType changeType) {
  switch (changeType) {
    case seedChangeType::ADD:
      return "ADD";
    case seedChangeType::DEL:
      return "DEL";
    case seedChangeType::SUB:
      return "SUB";
    default:
      return "UNKNOWN";
  }
}

void getSequenceFromReference(
  panmanUtils::Tree* tree,
  std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  std::vector<char>& blockExists,
  std::vector<char>& blockStrand,
  std::unordered_map<int, int>& blockLengths,
  std::string reference
) {
  if (tree->allNodes.find(reference) == tree->allNodes.end()) {
    logging::err("Reference sequence with matching name not found: {}", reference);
    std::exit(1);
  }
  panmanUtils::Node* referenceNode = tree->allNodes[reference];

  // get path from root to reference node
  std::vector<panmanUtils::Node*> pathFromRoot;
  panmanUtils::Node* it = referenceNode;
  while (it != tree->root) {
    pathFromRoot.push_back(it);
    it = it->parent;
  }
  pathFromRoot.push_back(tree->root);
  std::reverse(pathFromRoot.begin(), pathFromRoot.end());

  // get block sequence (blockSequence[i] = true if block i is on on the reference node)
  std::vector<char> blockSequence(tree->blocks.size() + 1, false);
  for (auto node : pathFromRoot) {
    for (const auto& blockMutation : node->blockMutation) {
      int32_t blockId = blockMutation.primaryBlockId;
      bool insertion = blockMutation.blockMutInfo;
      bool inversion = blockMutation.inversion;
      if (insertion) {
        blockSequence[blockId] = true;
      } else {
        if (!inversion) {
          blockSequence[blockId] = false;
        }
      }
    }
  }

  // initialize sequence, blockExists, blockStrand, blockLengths
  sequence.clear();
  blockExists.clear();
  blockStrand.clear();
  blockLengths.clear();
  sequence.resize(tree->blocks.size() + 1);
  blockExists.resize(tree->blocks.size() + 1, false);
  blockStrand.resize(tree->blocks.size() + 1, true);
  int32_t maxBlockId = 0;

  // fill in the skeleton of the sequence object
  for (int32_t blockId = 0; blockId < tree->blocks.size(); blockId++) {
    blockLengths[blockId] = 0;
    maxBlockId = std::max(maxBlockId, blockId);
    if (blockSequence[blockId]) {
      for (size_t i = 0; i < tree->blocks[blockId].consensusSeq.size(); i++) {
        bool endFlag = false;
        for (size_t j = 0; j < 8; j++) {
          const int nucCode = (((tree->blocks[blockId].consensusSeq[i]) >> (4*(7 - j))) & 15);
          if (nucCode == 0) {
            endFlag = true;
            break;
          }
          const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);
          sequence[blockId].push_back({nucleotide, {}});
        }
        if (endFlag) {
          break;
        }
      }

      sequence[blockId].push_back({'x', {}});
    } else {
      int len = 0;
      for (size_t i = 0; i < tree->blocks[blockId].consensusSeq.size(); i++) {
        bool endFlag = false;
        for (size_t j = 0; j < 8; j++) {
          const int nucCode = (((tree->blocks[blockId].consensusSeq[i]) >> (4*(7 - j))) & 15);
          if (nucCode == 0) {
            endFlag = true;
            break;
          }
          len++;
        }
        if (endFlag) {
          break;
        }
      }
      blockLengths[blockId] += len;
    }
  }

  sequence.resize(maxBlockId + 1);
  blockExists.resize(maxBlockId + 1);
  blockStrand.resize(maxBlockId + 1);

  // Assign nuc gaps
  auto& gaps = tree->gaps;
  for(size_t i = 0; i < gaps.size(); i++) {
    int32_t primaryBId = gaps[i].primaryBlockId;
    if (blockSequence[primaryBId]){
      for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
        int len = gaps[i].nucGapLength[j];
        int pos = gaps[i].nucPosition[j];
        sequence[primaryBId][pos].second.resize(len, '-');
      }
    } else {
      int len=0;
      for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
        len += gaps[i].nucGapLength[j];
      }
      blockLengths[primaryBId] += len;
    }
  }

  // apply mutations from root to reference node
  for (auto node : pathFromRoot) {
    // apply block mutations
    for (const auto& blockMutation : node->blockMutation) {
      int32_t blockId = blockMutation.primaryBlockId;
      bool insertion = blockMutation.blockMutInfo;
      bool inversion = blockMutation.inversion;
      if (!blockSequence[blockId]) {
        continue;
      }
      if (insertion) {
        blockExists[blockId] = true;
        blockStrand[blockId] = !inversion;
      } else {
        if (inversion) {
          blockStrand[blockId] = !blockStrand[blockId];
        } else {
          blockExists[blockId] = false;
          blockStrand[blockId] = true;
        }
      }
    }

    // apply nuc mutations
    for (const auto& nucMutation : node->nucMutation) {
      int32_t blockId = nucMutation.primaryBlockId;
      if (!blockSequence[blockId]) {
        continue;
      }
      int length = nucMutation.mutInfo >> 4;
      for (int i = 0; i < length; i++) {
        panmapUtils::Coordinate pos = panmapUtils::Coordinate(nucMutation, i);
        if (pos.nucPosition == sequence[pos.primaryBlockId].size() - 1 && pos.nucGapPosition == -1) {
          continue;
        } else if (pos.nucPosition >= sequence[pos.primaryBlockId].size()) {
          continue;
        }
        int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
        char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
        pos.setSequenceBase(sequence, newNuc);
      }
    }
  }
}

std::string getStringFromSequence(
  const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  const std::unordered_map<int, int>& blockLengths,
  const std::vector<char>& blockExists,
  const std::vector<char>& blockStrand,
  bool aligned
) {
  std::string seqString;
  for (size_t i = 0; i < blockExists.size(); i++) {
    if (blockExists[i]) {
      if (blockStrand[i]) {
        for(size_t j = 0; j < sequence[i].size(); j++) {
          // Gap nucs
          for (size_t k = 0; k < sequence[i][j].second.size(); k++) {
            if(sequence[i][j].second[k] != '-') {
              seqString += sequence[i][j].second[k];
            } else if(aligned) {
              seqString += '-';
            }
          }
          // Main nuc
          if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
            seqString += sequence[i][j].first;
          } else if (aligned && sequence[i][j].first != 'x') {
            seqString += '-';
          }
        }
      } else {
        for(size_t j = sequence[i].size()-1; j+1 > 0; j--) {
            // Main nuc first since we are iterating in reverse direction
            if(sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                seqString += panmanUtils::getComplementCharacter(sequence[i][j].first);
            } else if (aligned  && sequence[i][j].first != 'x') {
                seqString += '-';
            }

            // Gap nucs
            for(size_t k = sequence[i][j].second.size()-1; k+1 > 0; k--) {
                if(sequence[i][j].second[k] != '-') {
                    seqString += panmanUtils::getComplementCharacter(sequence[i][j].second[k]);
                } else if (aligned) {
                    seqString += '-';
                }
            }
        }
      }
    } else if (aligned){
      seqString.append(blockLengths.at(i), '-');
    }
  }
  return seqString;
}

std::string getStringFromReference(panmanUtils::Tree* tree, std::string reference, bool aligned) {
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequence;
  std::unordered_map<int, int> blockLengths;
  std::vector<char> blockExists;
  std::vector<char> blockStrand;
  getSequenceFromReference(tree, sequence, blockExists, blockStrand, blockLengths, reference);
  std::string seqString = getStringFromSequence(sequence, blockLengths, blockExists, blockStrand, aligned);
  return seqString;
}

void simulateSNPsOnSequence(
  std::string& sequence,
  std::vector<std::tuple<char, char, uint32_t>>& snpRecords,
  uint32_t numsnps,
  std::mt19937& rng
) {
  if (numsnps == 0) return;
  constexpr std::array<char, 3> notA = {'C', 'G', 'T'};
  constexpr std::array<char, 3> notC = {'A', 'G', 'T'};
  constexpr std::array<char, 3> notG = {'A', 'C', 'T'};
  constexpr std::array<char, 3> notT = {'A', 'C', 'G'};
  std::uniform_int_distribution<uint32_t> distNuc(0, 2);
  snpRecords.clear();
  snpRecords.reserve(numsnps);
  std::uniform_int_distribution<uint32_t> distPos(1000, sequence.size() - 1000);
  std::unordered_set<uint32_t> visitedPositions;
  while (snpRecords.size() < numsnps) {
    uint32_t pos = distPos(rng);
    if (visitedPositions.find(pos) != visitedPositions.end()) continue;
    visitedPositions.insert(pos);
    switch (sequence[pos]) {
      case 'A':
        snpRecords.emplace_back('A', notA[distNuc(rng)], pos);
        break;
      case 'C':
        snpRecords.emplace_back('C', notC[distNuc(rng)], pos);
        break;
      case 'G':
        snpRecords.emplace_back('G', notG[distNuc(rng)], pos);
        break;
      case 'T':
        snpRecords.emplace_back('T', notT[distNuc(rng)], pos);
        break;
      default:
        continue;
    }
  }

  std::sort(snpRecords.begin(), snpRecords.end(), [](const auto& a, const auto& b) {
    return std::get<2>(a) < std::get<2>(b);
  });

  for (const auto& [ref, alt, pos] : snpRecords) {
    sequence[pos] = alt;
  }

}

void LiteTree::cleanup() {
  for (auto& pair : allLiteNodes) {
    delete pair.second;
  }
  allLiteNodes.clear();
  blockScalarRanges.clear();
  nodeToDfsIndex.clear();
  root = nullptr;
}

uint32_t LiteTree::getBlockStartScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].first;
}

uint32_t LiteTree::getBlockEndScalar(const uint32_t blockId) const {
  return blockScalarRanges[blockId].second;
}

void LiteTree::initialize(::LiteTree::Reader liteTreeReader) {
  // initialize blockScalarRanges
  auto blockScalarRangesReader = liteTreeReader.getBlockRanges();
  blockScalarRanges.resize(blockScalarRangesReader.size());
  for (size_t i = 0; i < blockScalarRangesReader.size(); i++) {
    blockScalarRanges[i] = {blockScalarRangesReader[i].getRangeBeg(), blockScalarRangesReader[i].getRangeEnd()};
  }

  // initialize allLiteNodes
  auto liteNodesReader = liteTreeReader.getLiteNodes();
  for (size_t i = 0; i < liteNodesReader.size(); i++) {
    const auto liteNodeReader = liteNodesReader[i];
    const auto& nodeIdentifier = liteNodeReader.getId();
    const auto parentIndex = liteNodeReader.getParentIndex();
    nodeToDfsIndex.emplace(nodeIdentifier, i);
    auto [it, inserted] = allLiteNodes.emplace(nodeIdentifier, new LiteNode(nodeIdentifier, nullptr, {}));
    if (i == 0) continue;
    const auto parentNodeReader = liteNodesReader[parentIndex];
    const auto& parentNodeId = parentNodeReader.getId();
    it->second->parent = allLiteNodes[parentNodeId];
    allLiteNodes[parentNodeId]->children.push_back(it->second);
  }

  root = allLiteNodes[liteNodesReader[0].getId()];
}
//end of namespace panmapUtils
}
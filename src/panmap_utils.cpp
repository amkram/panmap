
#include "panmap_utils.hpp"
#include <queue>

namespace panmapUtils {

void getSequenceFromReference(panmanUtils::Tree* tree,
                              std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
                              std::vector<char>& blockExists,
                              std::vector<char>& blockStrand,
                              std::unordered_map<int, int>& blockLengths,
                              std::string reference) {
    if (tree->allNodes.find(reference) == tree->allNodes.end()) {
        logging::err("Reference sequence with matching name not found: {}", reference);
        std::exit(1);
    }
    panmanUtils::Node* referenceNode = tree->allNodes[reference];

    std::vector<panmanUtils::Node*> pathFromRoot;
    panmanUtils::Node* it = referenceNode;
    while (it != tree->root) {
        pathFromRoot.push_back(it);
        it = it->parent;
    }
    pathFromRoot.push_back(tree->root);
    std::reverse(pathFromRoot.begin(), pathFromRoot.end());

    // blockSequence[i] = true if block i is present on the reference node
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

    sequence.clear();
    blockExists.clear();
    blockStrand.clear();
    blockLengths.clear();
    sequence.resize(tree->blocks.size() + 1);
    blockExists.resize(tree->blocks.size() + 1, false);
    blockStrand.resize(tree->blocks.size() + 1, true);
    int32_t maxBlockId = 0;

    for (int32_t blockId = 0; blockId < tree->blocks.size(); blockId++) {
        blockLengths[blockId] = 0;
        maxBlockId = std::max(maxBlockId, blockId);
        if (blockSequence[blockId]) {
            forEachConsensusNuc(tree->blocks[blockId].consensusSeq, [&](int nucCode) {
                sequence[blockId].push_back({panmanUtils::getNucleotideFromCode(nucCode), {}});
            });

            sequence[blockId].push_back({'x', {}});
        } else {
            int len = 0;
            forEachConsensusNuc(tree->blocks[blockId].consensusSeq, [&](int) { len++; });
            blockLengths[blockId] += len;
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    blockStrand.resize(maxBlockId + 1);

    auto& gaps = tree->gaps;
    for (size_t i = 0; i < gaps.size(); i++) {
        int32_t primaryBId = gaps[i].primaryBlockId;
        if (blockSequence[primaryBId]) {
            for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                int len = gaps[i].nucGapLength[j];
                int pos = gaps[i].nucPosition[j];
                sequence[primaryBId][pos].second.resize(len, '-');
            }
        } else {
            int len = 0;
            for (size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
                len += gaps[i].nucGapLength[j];
            }
            blockLengths[primaryBId] += len;
        }
    }

    for (auto node : pathFromRoot) {
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
                int newNucCode = (nucMutation.nucs >> (4 * (5 - i))) & 0xF;
                char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
                pos.setSequenceBase(sequence, newNuc);
            }
        }
    }
}

std::string getStringFromSequence(const std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
                                  const std::unordered_map<int, int>& blockLengths,
                                  const std::vector<char>& blockExists,
                                  const std::vector<char>& blockStrand,
                                  bool aligned) {
    std::string seqString;
    for (size_t i = 0; i < blockExists.size(); i++) {
        if (blockExists[i]) {
            if (blockStrand[i]) {
                for (size_t j = 0; j < sequence[i].size(); j++) {
                    for (size_t k = 0; k < sequence[i][j].second.size(); k++) {
                        if (sequence[i][j].second[k] != '-') {
                            seqString += sequence[i][j].second[k];
                        } else if (aligned) {
                            seqString += '-';
                        }
                    }
                    if (sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        seqString += sequence[i][j].first;
                    } else if (aligned && sequence[i][j].first != 'x') {
                        seqString += '-';
                    }
                }
            } else {
                for (size_t j = sequence[i].size() - 1; j + 1 > 0; j--) {
                    // Main nuc first since iterating in reverse
                    if (sequence[i][j].first != '-' && sequence[i][j].first != 'x') {
                        seqString += panmanUtils::getComplementCharacter(sequence[i][j].first);
                    } else if (aligned && sequence[i][j].first != 'x') {
                        seqString += '-';
                    }

                    for (size_t k = sequence[i][j].second.size() - 1; k + 1 > 0; k--) {
                        if (sequence[i][j].second[k] != '-') {
                            seqString += panmanUtils::getComplementCharacter(sequence[i][j].second[k]);
                        } else if (aligned) {
                            seqString += '-';
                        }
                    }
                }
            }
        } else if (aligned) {
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
  // Sampling window: skip a 1000 bp flank when the sequence is long enough, else use
  // the whole sequence (guards size_t underflow of size()-1000 and lo>hi UB on short seqs).
  uint32_t lo, hi;
  if (sequence.size() > 2u * 1000u) {
    lo = 1000;
    hi = static_cast<uint32_t>(sequence.size()) - 1000;
  } else if (!sequence.empty()) {
    lo = 0;
    hi = static_cast<uint32_t>(sequence.size()) - 1;
  } else {
    return;
  }
  std::uniform_int_distribution<uint32_t> distPos(lo, hi);
  std::unordered_set<uint32_t> visitedPositions;
  const size_t windowSize = static_cast<size_t>(hi) - lo + 1;   // stop once every position has been tried
  while (snpRecords.size() < numsnps && visitedPositions.size() < windowSize) {
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

uint32_t LiteTree::getBlockStartScalar(const uint32_t blockId) const {
    return blockScalarRanges[blockId].first;
}

uint32_t LiteTree::getBlockEndScalar(const uint32_t blockId) const {
    return blockScalarRanges[blockId].second;
}

void LiteTree::initialize(::LiteTree::Reader liteTreeReader) {
    auto blockScalarRangesReader = liteTreeReader.getBlockRanges();
    blockScalarRanges.resize(blockScalarRangesReader.size());
    for (size_t i = 0; i < blockScalarRangesReader.size(); i++) {
        blockScalarRanges[i] = {blockScalarRangesReader[i].getRangeBeg(), blockScalarRangesReader[i].getRangeEnd()};
    }

    auto liteNodesReader = liteTreeReader.getLiteNodes();
    size_t numNodes = liteNodesReader.size();

    dfsIndexToNode.resize(numNodes, nullptr);

    for (size_t i = 0; i < numNodes; i++) {
        const auto liteNodeReader = liteNodesReader[i];
        const auto& nodeIdentifier = liteNodeReader.getId();
        const auto parentIndex = liteNodeReader.getParentIndex();
        nodeToDfsIndex.emplace(nodeIdentifier, i);
        auto [it, inserted] = allLiteNodes.emplace(nodeIdentifier, new LiteNode(nodeIdentifier, nullptr, {}));

        it->second->nodeIndex = i;
        dfsIndexToNode[i] = it->second;

        if (i == 0) continue;
        const auto parentNodeReader = liteNodesReader[parentIndex];
        const auto& parentNodeId = parentNodeReader.getId();
        it->second->parent = allLiteNodes[parentNodeId];
        allLiteNodes[parentNodeId]->children.push_back(it->second);
    }

    root = allLiteNodes[liteNodesReader[0].getId()];
}

}  // namespace panmapUtils

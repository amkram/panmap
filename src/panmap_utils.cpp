
#include "panmap_utils.hpp"

namespace panmapUtils {

void getSequenceFromReference(
  panmanUtils::Tree* tree,
  std::vector<std::vector<std::pair<char, std::vector<char>>>>& sequence,
  std::vector<bool>& blockExists,
  std::vector<bool>& blockStrand,
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
  std::vector<bool> blockSequence(tree->blocks.size() + 1, false);
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
        int newNucCode = (nucMutation.nucs >> (4*(5-i))) & 0xF;
        char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
        pos.setSequenceBase(sequence, newNuc);
      }
    }
  }
}

std::string getStringFromReference(panmanUtils::Tree* tree, std::string reference, bool aligned) {
  std::vector<std::vector<std::pair<char, std::vector<char>>>> sequence;
  std::unordered_map<int, int> blockLengths;
  std::vector<bool> blockExists;
  std::vector<bool> blockStrand;
  getSequenceFromReference(tree, sequence, blockExists, blockStrand, blockLengths, reference);
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
          } else if(aligned && sequence[i][j].first != 'x') {
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
      seqString.append(blockLengths[i], '-');
    }
  }
  return seqString;
}

//end of namespace panmapUtils
}
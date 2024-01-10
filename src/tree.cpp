#include "tree.hpp"
#include "pmi.hpp"
#include "seed.hpp"

using namespace PangenomeMAT;
using namespace tree;

std::string tree::getConsensus(Tree *T) {
    std::string consensus = ""; 
    for (size_t i = 0; i < T->blocks.size(); i++) {
        for (size_t j = 0; j < T->blocks[i].consensusSeq.size(); j++) {
            uint32_t c = T->blocks[i].consensusSeq[j];
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++) {
                const int nucCode = (c >> (4*(7 - k))) & 15;
                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                consensus += getNucleotideFromCode(nucCode);
            }
            if(endFlag) {
                break;
            }
        }
    }
    return consensus;
}
// pass by value bc we need to modify sequence object (and not persist changes) before getting nt string
std::string tree::getAlignedSequence(mutableTreeData data, Tree *T, const Node *node, const bool aligned) {
    // T should be const but [] operator on T->sequenceInverted is non-const
    if(T->rotationIndexes.find(node->identifier) != T->rotationIndexes.end()) {
        int ctr = -1, rotInd = 0;
        for(size_t i = 0; i < data.blockExists.size(); i++) {
            if(data.blockExists[i].first) {
                ctr++;
            }
            if(ctr == T->rotationIndexes[node->identifier]) {
                rotInd = i;
                break;
            }
        }
        std::rotate(data.sequence.begin(), data.sequence.begin() + rotInd, data.sequence.end());
        std::rotate(data.blockExists.begin(), data.blockExists.begin() + rotInd, data.blockExists.end());
        std::rotate(data.blockStrand.begin(), data.blockStrand.begin() + rotInd, data.blockStrand.end());
    }
    if(T->sequenceInverted.find(node->identifier) != T->sequenceInverted.end() && T->sequenceInverted[node->identifier]) {
        std::reverse(data.sequence.begin(), data.sequence.end());
        std::reverse(data.blockExists.begin(), data.blockExists.end());
        std::reverse(data.blockStrand.begin(), data.blockStrand.end());
    }

    std::string line;
    for(size_t i = 0; i < data.blockExists.size(); i++) {
        if(data.blockExists[i].first) {
            if(data.blockStrand[i].first) {
                for(size_t j = 0; j < data.sequence[i].first.size(); j++) {
                    for(size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
                        if(data.sequence[i].first[j].second[k] != '-') {
                            line += data.sequence[i].first[j].second[k];
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                    if(data.sequence[i].first[j].first != '-' && data.sequence[i].first[j].first != 'x') {
                        line += data.sequence[i].first[j].first;
                    } else if(aligned) {
                        line += '-';
                    }
                }
            } else {
                for(size_t j = data.sequence[i].first.size()-1; j+1 > 0; j--) {
                    if(data.sequence[i].first[j].first != '-' && data.sequence[i].first[j].first != 'x') {
                        line += getComplementCharacter(data.sequence[i].first[j].first);
                    } else if(aligned) {
                        line += '-';
                    }
                    for(size_t k = data.sequence[i].first[j].second.size()-1; k+1 > 0; k--) {
                        if(data.sequence[i].first[j].second[k] != '-') {
                            line += getComplementCharacter(data.sequence[i].first[j].second[k]);
                        } else if(aligned) {
                            line += '-';
                        }
                    }
                }   
            }
        } else if(aligned) {
            for(size_t j = 0; j < data.sequence[i].first.size(); j++) {
                for(size_t k = 0; k < data.sequence[i].first[j].second.size(); k++) {
                    line+='-';
                }
                line+='-';
            }
        }
    }
    return line;
}
void setupGlobalCoordinates(globalCoords_t &globalCoords, const BlockGapList &blockGaps, const std::vector<Block> &blocks, const std::vector<GapList> &gaps) {
    globalCoords.resize(blocks.size()+1);    
    // Assigning block gaps
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++){
        globalCoords[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
    }
    int32_t maxBlockId = 0;
    for(size_t i = 0; i < blocks.size(); i++){
        int32_t blockId = ((int32_t)blocks[i].primaryBlockId);
        maxBlockId = std::max(maxBlockId, blockId);
        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++){
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++){
                const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                if(nucCode == 0){
                    endFlag = true;
                    break;
                }
                globalCoords[blockId].first.push_back({0, {}});
            }
            if(endFlag){
                break;
            }
        }
        globalCoords[blockId].first.push_back({0, {}});
    }
    globalCoords.resize(maxBlockId + 1);
    // Assigning nucleotide gaps
    for(size_t i = 0; i < gaps.size(); i++){
        int32_t blockId = (gaps[i].primaryBlockId);
        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++){
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            globalCoords[blockId].first[pos].second.resize(len, 0);
        }
    }
    // Assigning coordinates
    int ctr = 0;
    for(size_t i = 0; i < globalCoords.size(); i++) {
        for(size_t j = 0; j < globalCoords[i].second.size(); j++) {
            for(size_t k = 0; k < globalCoords[i].second[j].size(); k++) {
                for(size_t w = 0; w < globalCoords[i].second[j][k].second.size(); w++) {
                    globalCoords[i].second[j][k].second[w] = ctr;
                    ctr++;
                }
                globalCoords[i].second[j][k].first = ctr;
                ctr++;
            }
        }
        for(size_t j = 0; j < globalCoords[i].first.size(); j++) {
            for(size_t k = 0; k < globalCoords[i].first[j].second.size(); k++) {
                globalCoords[i].first[j].second[k] = ctr;
                ctr++;
            }
            globalCoords[i].first[j].first = ctr;
            ctr++;
        }
    }
}
void tree::removeIndices(std::vector<kmer_t> &v, std::stack<int32_t> &rm) {
    if (rm.size() < 1) {
        return;
    }
    int32_t rmVal = rm.top();
    rm.pop();
    v.erase(
        std::remove_if(std::begin(v), std::end(v), [&](kmer_t& elem)
        {
            if (rmVal == -1) {
                return false;
            }
            if (&elem - &v[0] == rmVal) {
                if (!rm.empty()) {
                    rmVal = rm.top();
                    rm.pop();
                } else {
                    rmVal = -1;
                }
                return true;
            }
            return false;

        }),
        std::end(v)
    );
}
void tree::setup(mutableTreeData &data, globalCoords_t &globalCoords, Tree *T) {
    const BlockGapList &blockGaps = T->blockGaps;
    const std::vector< GapList > &gaps = T->gaps;
    const std::vector< Block > &blocks = T->blocks;

    sequence_t sequence(blocks.size() + 1);
    blockExists_t blockExists(blocks.size() + 1, {false, {}});
    blockStrand_t blockStrand(blocks.size() + 1, {true, {}});
    int32_t maxBlock = 0;
    
    for(size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
        sequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
        blockStrand[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], true);
    }
    for(size_t i = 0; i < blocks.size(); i++) {
        int32_t b = ((int32_t)blocks[i].primaryBlockId);
        maxBlock = std::max(maxBlock, b);
        for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
            bool stop = false;
            for(size_t k = 0; k < 8; k++) {
                const int nc = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                if(nc == 0) {
                    stop = true;
                    break;
                }
                const char c = PangenomeMAT::getNucleotideFromCode(nc);
                sequence[b].first.push_back({c, {}});
            }
            if(stop) {
                break;
            }
        }
        sequence[b].first.push_back({'x', {}});
    }

    sequence.resize(maxBlock + 1);
    blockExists.resize(maxBlock + 1);
    blockStrand.resize(maxBlock + 1);

    // Assigning nucleotide gaps in blocks
    for(size_t i = 0; i < gaps.size(); i++) {
        int32_t id = (gaps[i].primaryBlockId);
        for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
            int len = gaps[i].nucGapLength[j];
            int pos = gaps[i].nucPosition[j];
            sequence[id].first[pos].second.resize(len, '-');
        }
    }

    data.sequence = sequence;
    data.blockExists = blockExists;
    data.blockStrand = blockStrand;
    setupGlobalCoordinates(globalCoords, blockGaps, blocks, gaps);
}

size_t tree::getGlobalCoordinate(const int blockId, const int nucPosition, const int nucGapPosition, const globalCoords_t &globalCoords) {
    if(nucGapPosition == -1){
        return globalCoords[blockId].first[nucPosition].first;
    }
    return globalCoords[blockId].first[nucPosition].second[nucGapPosition];
}

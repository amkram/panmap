#include "tree.hpp"
#include "pmi.hpp"
#include "seed.hpp"
#include <cmath>

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
void tree::updateConsensus(mutableTreeData &data, Tree *T) {
    std::string consensus = tree::getStringFromCurrData(data, T, T->root, true);
    data.gappedConsensus = consensus;
    std::string ungapped = "";
    data.degap.clear();
    data.regap.clear();
    int32_t ct = 0;
    for (int32_t i = 0; i < consensus.size(); i ++) {
        char &c = consensus[i];
        data.degap.push_back(ungapped.size());
        data.regap.push_back(i + ct);
        if (c != '-') {
            ungapped += c;
        } else {
            ct++;
        }
    }
    data.ungappedConsensus = ungapped;
}
void getAllStringsHelper(std::unordered_map<std::string, std::string> &strings, mutableTreeData &data, Tree *T, const Node *node, globalCoords_t &globalCoords) {
    /*  Mutate with block and nuc mutations */
    blockMutData_t blockMutData;
    nucMutData_t nucMutData;
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    /* Use current state of mutableTreeData to decode node's sequence */
    std::string seq = tree::getStringFromCurrData(data, T, node, true);

    strings[node->identifier] = seq;

    /* Recursive step */
    for(Node* child: node->children){
        getAllStringsHelper(strings, data, T, child, globalCoords);
    }

    /* Undo mutations when backtracking */
    pmi::seedIndex blank;
    undoMutations(data, blank, T, node, blockMutData, nucMutData);
}
std::unordered_map<std::string, std::string> tree::getAllNodeStrings(Tree *T) {
    std::unordered_map<std::string, std::string> strings;
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    setup(data, globalCoords, T);

    getAllStringsHelper(strings, data, T, T->root, globalCoords);

    return strings;
}
// pass by value bc we need to modify sequence object (and not persist changes) before getting nt string
std::string tree::getStringFromCurrData(mutableTreeData data, Tree *T, const Node *node, const bool aligned) {
    // T should be const but [] operator on T->sequenceInverted is non-const
    std::string line;
    if (node == nullptr) { // consensus sequence (all blocks on) rather than a node in the tree
       for (size_t i = 0; i < T->blocks.size(); i++) {
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
       }
        return line;
    }

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

static int getIndexFromNucleotide(char nuc) {
    switch(nuc) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case '*':
            return 4;
        default:
            return 5;
    }
    return 5;
}

void buildMutationMatrices(mutationMatrices& mutMat, Tree* T) {
    std::unordered_map<std::string, std::string> alignedSequences = getAllNodeStrings(T);
    for (const auto& sequence : alignedSequences) {
        std::string parentId;
        if (T->allNodes[sequence.first]->parent == nullptr || T->allNodes[sequence.first]->nucMutation.size() == 0) {
            continue;
        } else {
            parentId = T->allNodes[sequence.first]->parent->identifier;
        }

        const std::string& curSeq = sequence.second;
        const std::string& parSeq = alignedSequences[parentId];
        size_t insLen = 0;
        size_t delLen = 0;

        for (size_t i = 0; i < parSeq.size(); i++) {
            if (parSeq[i] == '-' && curSeq[i] == '-') {
                continue;
            } else if (parSeq[i] != '-' && curSeq[i] == '-') {
                delLen++;
                if (insLen > 0) {
                    if (insLen > mutMat.insmat.size() - 1) {
                        mutMat.insmat.resize(insLen + 1);
                    }
                    mutMat.insmat[insLen]++;
                    mutMat.total_insmut++;
                    insLen = 0;
                }
            } else if (parSeq[i] == '-' && curSeq[i] != '-') {
                insLen++;
                if (delLen > 0) {
                    if (delLen > mutMat.delmat.size() - 1) {
                        mutMat.delmat.resize(delLen + 1);
                    }
                    mutMat.delmat[delLen]++;
                    mutMat.total_delmut++;
                    delLen = 0;
                }
            } else {
                if (insLen > mutMat.insmat.size() - 1) {
                    mutMat.insmat.resize(insLen + 1);
                }
                mutMat.insmat[insLen]++;
                mutMat.total_insmut++;
                insLen = 0;

                if (delLen > mutMat.delmat.size() - 1) {
                    mutMat.delmat.resize(delLen + 1);
                }
                mutMat.delmat[delLen]++;
                mutMat.total_delmut++;
                delLen = 0;

                int parNucIdx = getIndexFromNucleotide(parSeq[i]);
                int curNucIdx = getIndexFromNucleotide(curSeq[i]);
                if (parNucIdx > 3 || curNucIdx > 3) { continue; }
                mutMat.submat[parNucIdx][curNucIdx]++;
                mutMat.total_submuts[parNucIdx]++;
            }
        }
    }

    // insertion
    for (auto i = 0; i < mutMat.insmat.size(); ++i) {
        mutMat.insmat[i] = -10 * log10f64x(mutMat.insmat[i] / mutMat.total_insmut);
    }
    // deletion
    for (auto i = 0; i < mutMat.delmat.size(); ++i) {
        mutMat.delmat[i] =  -10 * log10f64x(mutMat.delmat[i] / mutMat.total_delmut);
    }
    // substitution
    for (auto i = 0; i < 4; i++) {
        for (auto j = 0; j < 4; j++) {
            mutMat.submat[i][j] = -10 * log10f64x(mutMat.submat[i][j] / mutMat.total_submuts[i]);
        }
    }
}

void tree::printMutationMatrices(Tree* T, std::ofstream* ofptr) {
    mutationMatrices mutMat = mutationMatrices();

    buildMutationMatrices(mutMat, T);

    if (ofptr == nullptr) {
        for (const auto& row : mutMat.submat) {
            for (const auto& prob : row) {
                std::cout << prob << " ";
            }
            std::cout << std::endl;
        }
        for (const auto& prob : mutMat.insmat) {
            std::cout << prob << " ";
        }
        std::cout << std::endl;
        for (const auto& prob : mutMat.delmat) {
            std::cout << prob << " ";
        }
    } else {
        for (const auto& row : mutMat.submat) {
            for (const auto& prob : row) {
                *ofptr << prob << " ";
            }
            *ofptr << "\n";
        }
        for (const auto& prob : mutMat.insmat) {
            *ofptr << prob << " ";
        }
        *ofptr << "\n";
        for (const auto& prob : mutMat.delmat) {
            *ofptr << prob << " ";
        }
    }
    
}

void tree::fillMutationMatrices(mutationMatrices& mutMat, Tree* T, std::ifstream* infptr) {
    if (infptr != nullptr) {
        // read mutation matrix from file if provided
        std::string line;
        int idx = 0;
        while(getline(*infptr, line)) {
            std::vector<double> probs;
            std::vector<std::string> fields;
            stringSplit(line, ' ', fields);
            for (const auto& f : fields) {
                probs.push_back(std::stod(f));
            }
            if (idx < 4) {
                mutMat.submat[idx] = move(probs);
            } else if (idx == 4) {
                mutMat.insmat = move(probs);
            } else {
                mutMat.delmat = move(probs);
            }
            idx++;
        }
        return;
    } else {
        buildMutationMatrices(mutMat, T);
    }

}
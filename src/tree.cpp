#pragma once
#include "tree.hpp"
#include "pmi.hpp"
#include "seeding.hpp"
#include <cmath>

using namespace PangenomeMAT;
using namespace tree;

std::chrono::time_point<std::chrono::high_resolution_clock> global_timer = std::chrono::high_resolution_clock::now();
void time_stamp(){
    std::chrono::time_point<std::chrono::high_resolution_clock> newtime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = newtime - global_timer;
    std::cerr << "timing " << duration.count() << "\n\n";
    global_timer = newtime;
}


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

/*
- regap not used anywhere previously -> changed how regap is constructed for the purpose
  of fixing syncmer/seedemr update bug and tracking end positions.
*/
void tree::updateConsensus(mutableTreeData &data, Tree *T) { //TODO use mutation data

    std::string consensus = tree::getStringFromCurrData(data, T, T->root, true);

    data.gappedConsensus = consensus;
    std::string ungapped = "";
    data.degap.clear();
    data.regap.clear();

    for (int32_t i = 0; i < consensus.size(); i ++) {
        char &c = consensus[i];
        data.degap.push_back(ungapped.size());
        if (c != '-') {
            ungapped += c;
            data.regap.push_back(i);
        }
    }

    time_stamp();
    
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


// TODO make an update string function (that takes it by reference) that uses block mutations
std::string tree::getStringFromCurrData(mutableTreeData &data, Tree *T, const Node *node, const bool aligned) {
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
    

    for(size_t i = 0; i < data.blockExists.size(); i++) {         // TODO use block mutations rather than loop through all
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
void tree::removeIndices(std::vector<seed> &v, std::stack<int32_t> &rm) {
    if (rm.size() < 1) {
        return;
    }
    int32_t rmVal = rm.top();
    rm.pop();
    v.erase(
        std::remove_if(std::begin(v), std::end(v), [&](seed& elem)
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

static size_t getBeg(const std::string& s1, const std::string& s2, size_t window, double threshold) {
    if (s1.empty()) {
        return 0;
    }

    size_t numAlign = 0;
    size_t numMatch = 0;
    size_t beg = 0;
    size_t idx = 0;
    std::queue<size_t> begs;
    while (idx < s1.size()) {
        if (s1[idx] == '-' && s2[idx] == '-') {
            ++idx;
            continue;
        } 
        if (beg == 0) {
            beg = idx;
        } else {
            begs.push(idx);
        }
        if (s1[idx] == s2[idx]) {
            ++numMatch;
        }
        ++numAlign;
        ++idx;

        if (numAlign == window) {
            double pcid = static_cast<double>(numMatch) / static_cast<double>(window);
            if (pcid >= threshold && s1[beg] == s2[beg]) {
                return beg;
            }
            
            if (s1[beg] == s2[beg]) {
                --numMatch;
            }
            --numAlign;
            beg = begs.front();
            begs.pop();
        }
    }

    return s1.size() - 1;
}

static size_t getEnd(const std::string& s1, const std::string& s2, size_t window, double threshold) {
    if (s1.empty()) {
        return 0;
    }

    size_t numAlign = 0;
    size_t numMatch = 0;
    size_t end = s1.size();
    size_t idx = s1.size() - 1;
    std::queue<size_t> ends;

    while (true) {
        if (s1[idx] == '-' && s2[idx] == '-') {
            if (idx == 0) {
                break;
            }
            --idx;
            continue;
        }
        if (end == s1.size()) {
            end = idx;
        } else {
            ends.push(idx);
        }
        if (s1[idx] == s2[idx]) {
            ++numMatch;
        }
        ++numAlign;

        if (idx == 0){
            break;
        }
        --idx;

        if (numAlign == window) {
            double pcid = static_cast<double>(numMatch) / static_cast<double>(window);
            if (pcid >= threshold && s1[end] == s2[end]) {
                return end;
            }

            if (s1[end] == s2[end]) {
                --numMatch;
            }
            --numAlign;
            end = ends.front();
            ends.pop();
        }
    }

    return 0;
}

std::pair<size_t, size_t> tree::getMaskCoorsForMutmat(const std::string& s1, const std::string& s2, size_t window, double threshold) {
    assert(s1.size() == s2.size());
    if (window == 0 || threshold == 0.0) {
        return std::make_pair<size_t, size_t>(0, s1.size()-1);
    }
    return std::make_pair<size_t, size_t>(getBeg(s1, s2, window, threshold), getEnd(s1, s2, window, threshold));
}

void buildMutationMatrices(mutationMatrices& mutMat, Tree* T, size_t window, double threshold) {
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
        std::pair<size_t, size_t> edgeCoor = tree::getMaskCoorsForMutmat(curSeq, parSeq, window, threshold);
        assert(edgeCoor.second >= edgeCoor.first);
        for (size_t i = edgeCoor.first; i < edgeCoor.second + 1; i++) {
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
        mutMat.insmat[i] = -10 * log10f(mutMat.insmat[i] / mutMat.total_insmut);
    }
    // deletion
    for (auto i = 0; i < mutMat.delmat.size(); ++i) {
        mutMat.delmat[i] =  -10 * log10f(mutMat.delmat[i] / mutMat.total_delmut);
    }
    // substitution
    for (auto i = 0; i < 4; i++) {
        for (auto j = 0; j < 4; j++) {
            mutMat.submat[i][j] = -10 * log10f(mutMat.submat[i][j] / mutMat.total_submuts[i]);
        }
    }
}

void tree::writeMutationMatrices(const mutationMatrices& mutMat, std::ofstream& mmfout) {
    for (const std::vector<double>& row : mutMat.submat) {
        for (const double& prob : row) {
            mmfout << prob << " ";
        }
        mmfout << "\n";
    }
    for (const double& prob : mutMat.insmat) {
        mmfout << prob << " ";
    }
    mmfout << "\n";
    for (const double& prob : mutMat.delmat) {
        mmfout << prob << " ";
    }
    mmfout << "\n";
}

void tree::fillMutationMatricesFromTree(mutationMatrices& mutMat, Tree* T, size_t window, double threshold) {
    buildMutationMatrices(mutMat, T, window, threshold);
    mutMat.filled = true;
}

void tree::fillMutationMatricesFromFile(mutationMatrices& mutMat, std::ifstream& inf) {
    std::string line;
    int idx = 0;
    while(getline(inf, line)) {
        std::vector<double> probs;
        std::vector<std::string> fields;
        stringSplit(line, ' ', fields);
        for (const auto& f : fields) {
            probs.push_back(std::stod(f));
        }
        if (probs.size() == 0) {
            break;
        }
        if (idx < 4) {
            if (probs.size() != 4) {
                throw std::invalid_argument("Received invalid mutamtion matrix (.mm) file");
            }
            mutMat.submat[idx] = std::move(probs);
        } else if (idx == 4) {
            if (probs.size() < 1) {
                throw std::invalid_argument("Received invalid mutamtion matrix (.mm) file");
            }
            mutMat.insmat = std::move(probs);
        } else if (idx == 5) {
            if (probs.size() < 1) {
                throw std::invalid_argument("Received invalid mutamtion matrix (.mm) file");
            }
            mutMat.delmat = std::move(probs);
        }
        idx++;
    }

    if (idx != 6) {
        throw std::invalid_argument("Received invalid mutamtion matrix (.mm) file");
    }
    mutMat.filled = true;
}

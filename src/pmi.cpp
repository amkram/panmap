#include "pmi.hpp"
#include "tree.hpp"
#include "seed.hpp"

using namespace seed;
using namespace pmi;
using namespace PangenomeMAT;
using namespace tree;

/* Helpers */
bool compareTuples(const range_t &a, const range_t &b) {
    return std::get<0>(a) < std::get<0>(b);
}
std::vector<range_t> mergePositions(const mutableTreeData &data, const nucMutData_t &nucMutData, const int32_t pad) {
    auto &blockExists = data.blockExists;
    std::vector<range_t> tuples;
    for(auto it = nucMutData.begin(); it != nucMutData.end(); it++) {
        auto mutation = *it;
        int32_t blockId = std::get<0>(mutation);
        int32_t globalCoord = std::get<5>(mutation);
        int32_t mutLen = std::get<6>(mutation);

        tuples.push_back(std::make_tuple(globalCoord, mutLen, blockId));
    }
    if (tuples.size() == 0) {
        return tuples;
    }
    std::sort(tuples.begin(), tuples.end(), compareTuples);
    std::vector<range_t> mergedTuples;
    int32_t start = std::get<0>(tuples[0]) - pad + 1;
    int32_t len = std::get<1>(tuples[0]) + pad - 1;
    mergedTuples.push_back({std::max(0, start), len, std::get<2>(tuples[0])});
    
    for (int32_t i = 1; i < tuples.size(); i++) {
        int32_t currentStart = std::max(0, std::get<0>(tuples[i]) - pad + 1);
        int32_t currentEnd = currentStart + std::get<1>(tuples[i]) + pad - 1;
        int32_t prevEnd = std::get<0>(mergedTuples.back()) + std::get<1>(mergedTuples.back());
        int32_t prevStart = std::get<0>(mergedTuples.back());

        if (currentStart > prevEnd) {
            mergedTuples.push_back({std::max(0, currentStart), currentEnd - currentStart, std::get<2>(tuples[i])});
        }
        else if (currentStart > prevStart) {
            std::get<1>(mergedTuples.back()) = currentEnd - prevStart;
        }
    }
    return mergedTuples;
}
range_t getRecomputePositions(const range_t &p, const std::string &gappedSequence, const int32_t k) {
    int32_t mutPos = std::get<0>(p);
    int32_t mutLen = std::get<1>(p);
    int32_t numSeen = 0;
    int32_t start = mutPos;

    while(numSeen < k && start > 0) {
        if (gappedSequence[start] != '-') {
            numSeen++;
            if (numSeen == k) {
                break;
            }
        }
        start--;
    }
    numSeen = 0;

    int32_t stop = std::min((int32_t) gappedSequence.size()-1, mutPos + mutLen - 1);
    while(numSeen < k && stop < gappedSequence.size()) {
        if (gappedSequence[stop] != '-') {
            numSeen++;
            if (numSeen == k) {
                break;
            }
        }
        stop++;
    }
    stop = std::min(stop, (int32_t) gappedSequence.size() - 1);
    return std::make_tuple(start, stop, std::get<2>(p));
}
std::vector<range_t> getAffectedRanges(mutableTreeData &data, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData, std::string &seq, Tree *T, const int32_t k, const globalCoords_t &globalCoords) {
    std::vector<range_t> affectedRanges;
    for (const auto &tup : mergePositions(data, nucMutData, k)) {
        range_t r = getRecomputePositions(tup, seq, k);
        affectedRanges.push_back(r);
    }
    for (const auto &p : blockMutData) { // when a block turns on or off, seeds should be recomputed
        const int32_t blockId = std::get<0>(p);
        int32_t blockStart = tree::getGlobalCoordinate(blockId, 0, -1, globalCoords);
        affectedRanges.push_back(std::tuple(blockStart, std::min(blockStart + k, (int32_t) T->blocks[std::get<0>(p)].consensusSeq.size()), blockId));
    }
    return affectedRanges;
}
std::vector<kmer_t> getUniqueSeeds(const std::string &seq, const int32_t k, const int32_t s) {
    std::vector<kmer_t> seedVec = seed::syncmerize(seq, k, s, false, true, 0);
    std::set<kmer_t> seedSet(seedVec.begin(), seedVec.end());
    seedVec.assign(seedSet.begin(), seedSet.end());
    return seedVec;
}
void applyMutations(mutableTreeData &data, blockMutData_t &blockMutData, nucMutData_t &nucMutData, Tree *T, const Node *node, const globalCoords_t &globalCoords) {

    blockExists_t &blockExists = data.blockExists;
    blockStrand_t &blockStrand = data.blockStrand;

    for(const auto &mutation: node->blockMutation) {
        int32_t blockId = mutation.primaryBlockId;
        bool type = mutation.blockMutInfo;
        bool inversion = mutation.inversion;

        if(type == 1) {
            // insertion
            bool oldStrand;
            bool oldMut;
            oldStrand = blockStrand[blockId].first;
            oldMut = blockExists[blockId].first;
            blockExists[blockId].first = true;
            // if insertion of inverted block takes place, the strand is backwards
            blockStrand[blockId].first = !inversion;
            blockMutData.push_back( std::make_tuple(blockId, oldMut, oldStrand, true, !inversion) );
        } else {
            bool oldMut;
            bool oldStrand;
            if(inversion) {
                oldStrand = blockStrand[blockId].first;
                oldMut = blockExists[blockId].first;
                blockStrand[blockId].first = !oldStrand;
                blockMutData.push_back( std::make_tuple(blockId, oldMut, oldStrand, oldMut, !oldStrand) );
            } else {
                // Actually a deletion
                oldStrand = blockStrand[blockId].first;
                oldMut = blockExists[blockId].first;
                blockExists[blockId].first = false;
                // resetting strand to true during deletion
                blockStrand[blockId].first = true;
            }
            blockMutData.push_back( std::make_tuple(blockId, oldMut, oldStrand, false, true) );
        }
    }
    // Nuc mutations
    for(size_t i = 0; i < node->nucMutation.size(); i++) {
        int32_t blockId = node->nucMutation[i].primaryBlockId;
        int32_t nucPosition = node->nucMutation[i].nucPosition;
        int32_t nucGapPosition = node->nucMutation[i].nucGapPosition;
        uint32_t type = (node->nucMutation[i].mutInfo & 0x7);
        char newVal = '-';
        size_t globalCoord = tree::getGlobalCoordinate(blockId, nucPosition, nucGapPosition, globalCoords);

        if(type < 3) { // Either S, I or D
            int len = ((node->nucMutation[i].mutInfo) >> 4);
            if(type == PangenomeMAT::NucMutationType::NS) {
                // Substitution     
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        char oldVal = data.sequence[blockId].first[nucPosition].second[nucGapPosition+j];
                        newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        data.sequence[blockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                        nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition+j, oldVal, newVal, globalCoord, len));   
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        char oldVal = data.sequence[blockId].first[nucPosition+j].first;
                        newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        data.sequence[blockId].first[nucPosition+j].first = newVal;
                        nucMutData.push_back(std::make_tuple(blockId, nucPosition + j, nucGapPosition, oldVal, newVal, globalCoord, len));   
                    }
                }
            }
            else if(type == PangenomeMAT::NucMutationType::NI) {
                // Insertion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        char oldVal = data.sequence[blockId].first[nucPosition].second[nucGapPosition+j];
                        newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        data.sequence[blockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                        nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition+j, oldVal, newVal, globalCoord, len));
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        char oldVal = data.sequence[blockId].first[nucPosition+j].first;
                        const int nucCode = ((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF;
                        newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                        data.sequence[blockId].first[nucPosition+j].first = newVal;
                        nucMutData.push_back(std::make_tuple(blockId, nucPosition + j, nucGapPosition, oldVal, newVal, globalCoord, len));   
                    }
                }
            }
            else if(type == PangenomeMAT::NucMutationType::ND) {
                // Deletion
                if(nucGapPosition != -1) {
                    for(int j = 0; j < len; j++) {
                        char oldVal = data.sequence[blockId].first[nucPosition].second[nucGapPosition+j];
                        data.sequence[blockId].first[nucPosition].second[nucGapPosition+j] = '-';
                        nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition+j, oldVal, '-', globalCoord, len));
                    }
                } else {
                    for(int j = 0; j < len; j++) {
                        char oldVal = data.sequence[blockId].first[nucPosition+j].first;
                        data.sequence[blockId].first[nucPosition+j].first = '-';
                        nucMutData.push_back(std::make_tuple(blockId, nucPosition + j, nucGapPosition, oldVal, '-', globalCoord, 0));
                    }
                }
            }
        } 
        else {
            int len = 0;
            if(type == PangenomeMAT::NucMutationType::NSNPS) {
                // SNP Substitution                
                newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                if(nucGapPosition != -1) {
                    char oldVal = data.sequence[blockId].first[nucPosition].second[nucGapPosition];
                    data.sequence[blockId].first[nucPosition].second[nucGapPosition] = newVal;
                    nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));
                } else {
                    char oldVal = data.sequence[blockId].first[nucPosition].first;
                    data.sequence[blockId].first[nucPosition].first = newVal;
                    nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));   
                }
            }
            else if(type == PangenomeMAT::NucMutationType::NSNPI) {
                // SNP Insertion
                len = 1;
                newVal = PangenomeMAT::getNucleotideFromCode(((node->nucMutation[i].nucs) >> 20) & 0xF);
                if(nucGapPosition != -1) {
                    char oldVal = data.sequence[blockId].first[nucPosition].second[nucGapPosition];
                    data.sequence[blockId].first[nucPosition].second[nucGapPosition] = newVal;
                    nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));   
                } else {
                    char oldVal = data.sequence[blockId].first[nucPosition].first;
                    data.sequence[blockId].first[nucPosition].first = newVal;
                    nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal, newVal, globalCoord, len));   
                }
            }
            else if(type == PangenomeMAT::NucMutationType::NSNPD) {
                // SNP Deletion
                if(nucGapPosition != -1) {
                    char oldVal = data.sequence[blockId].first[nucPosition].second[nucGapPosition];
                    data.sequence[blockId].first[nucPosition].second[nucGapPosition] = '-';
                    nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal, '-', globalCoord, len));
                } else {
                    char oldVal = data.sequence[blockId].first[nucPosition].first;
                    data.sequence[blockId].first[nucPosition].first = '-';
                    nucMutData.push_back(std::make_tuple(blockId, nucPosition, nucGapPosition, oldVal, '-', globalCoord, len));
                 }
            }
        }
    }
}
void discardSeeds(std::vector<kmer_t> &seeds, std::unordered_map<std::string, kmer_t> &newSeeds, seedIndex &index, const std::vector<range_t>& B, const std::string &seq, const std::string nid, const size_t k) {

    auto cmp = [](kmer_t a, kmer_t b) { return a.idx > b.idx; };
    std::set<kmer_t, decltype(cmp)> dels(cmp);

    for (int32_t i = 0; i < seeds.size(); i++) {
        kmer_t &s = seeds[i];
        for (int32_t j = 0; j < B.size(); j++) {
            const auto &b = B[j];
            if (std::get<0>(b)> s.pos) {
                break;
            } else if (s.gappedEnd <= std::get<1>(b)) {
                auto it = newSeeds.find(s.seq);
                if (it == newSeeds.end()) {
                    s.idx = i;
                    dels.insert(s);
                } else {
                    newSeeds.erase(it);
                    s.pos = it->second.pos;
                    s.gappedEnd = it->second.gappedEnd;
                }
            }
        }
    }
    
    index.deletions[nid] = std::vector<kmer_t>(dels.begin(), dels.end());
    std::stack<int32_t> rmDel;
    for (const kmer_t &d : dels) {
        rmDel.push(d.idx);
    }
    removeIndices(seeds, rmDel);
}
void undoMutations(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const blockMutData_t &blockMutData, const nucMutData_t &nucMutData) {
    blockExists_t &blockExists = data.blockExists;
    blockStrand_t &blockStrand = data.blockStrand;

    for(auto it = blockMutData.rbegin(); it != blockMutData.rend(); it++) {
        auto mutation = *it;
        blockExists[std::get<0>(mutation)].first = std::get<1>(mutation);
        blockStrand[std::get<0>(mutation)].first = std::get<2>(mutation);
    }
    // Undo nuc mutations when current node and its subtree have been processed
    for(auto it = nucMutData.rbegin(); it != nucMutData.rend(); it++) {
        auto mutation = *it;
        if(std::get<2>(mutation) != -1) {
            data.sequence[std::get<0>(mutation)].first[std::get<1>(mutation)].second[std::get<2>(mutation)] = std::get<3>(mutation);
        } else {
            data.sequence[std::get<0>(mutation)].first[std::get<1>(mutation)].first = std::get<3>(mutation);
        }
    }
    if (index.insertions[node->identifier].size() > 0) {
        data.seeds.erase(data.seeds.end() - index.insertions[node->identifier].size(), data.seeds.end());
    }
    if (index.deletions[node->identifier].size() > 0) {
        for (int32_t i = index.deletions[node->identifier].size() - 1; i >= 0; i--) {
           data.variableSeeds[index.deletions[node->identifier][i].seq] = true;
           data.seeds.insert(data.seeds.begin() + index.deletions[node->identifier][i].idx, index.deletions[node->identifier][i]);
        }
    }
}
void addSeeds(std::vector<kmer_t> &seeds, seedIndex &index, const std::string nodeId, const std::unordered_map<std::string, kmer_t> &newSeeds) {
    for (auto &s : newSeeds) {
        seeds.push_back(s.second);
        index.insertions[nodeId].push_back(s.second);
    }
}
void recomputeSeeds(mutableTreeData &data, std::unordered_map<std::string, kmer_t> &newSeeds, const std::vector<range_t> &ranges, const std::string &sequence, int32_t k, int32_t s) {
    for (const range_t &range : ranges) {
        int32_t start = std::get<0>(range);
        int32_t stop = std::get<1>(range);
        if (start >= stop) {
            continue;
        }
        int32_t seqLen = sequence.size();
        std::string redo = sequence.substr(std::max(0, start), std::min(seqLen - start, 1 + stop - start)); 
        auto redone = syncmerize(redo, k, s, false, true, std::max(0, start));
        for (const kmer_t &syncmer : redone) {
            newSeeds[syncmer.seq] = syncmer;
        }
    }
}

void writeIndexDFS(std::stringstream &ss, std::vector<kmer_t> &seeds, Node *currNode, seedIndex &index) {
    ss << currNode->identifier << "\t";
    std::stack<int32_t> rmDel;
    for (const kmer_t &d : index.deletions[currNode->identifier]) {
        rmDel.push(d.idx);
        ss << d.idx << " ";
    }
    removeIndices(seeds, rmDel);
    ss << "^";
    for (const kmer_t &s : index.insertions[currNode->identifier]) {
        seeds.push_back(s);
        ss << s.seq << " ";
    }
    ss << "\n";
    for (Node *child : currNode->children) {
        writeIndexDFS(ss, seeds, child, index);
    }

    seeds.erase(seeds.end() - index.insertions[currNode->identifier].size(), seeds.end());
    
    for (int32_t i = index.deletions[currNode->identifier].size() - 1; i >= 0; i--) {
        seeds.insert(seeds.begin() + index.deletions[currNode->identifier][i].idx, index.deletions[currNode->identifier][i]);
    }
}
void buildHelper(mutableTreeData &data, seedIndex &index, Tree *T, const Node *node, const size_t k, const size_t s, const globalCoords_t &globalCoords) {

    blockMutData_t blockMutData;
    nucMutData_t nucMutData;

    /* Mutate with block and nuc mutations. */
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    /* Use the current state of mutableTreeData to decode node's sequence. */
    std::string currNodeSequence = tree::getAlignedSequence(data, T, node, true);

    /* Find ranges overlapping mutations at this node. */
    std::vector<range_t> affectedRanges = getAffectedRanges(data, blockMutData, nucMutData, currNodeSequence, T, k, globalCoords);

    /* Recompute seeds across affected ranges. If a block is off, all seeds are dropped and not recomputed. */
    std::unordered_map<std::string, kmer_t> newSeeds;
    recomputeSeeds(data, newSeeds, affectedRanges, currNodeSequence, k, s);

    /* Discard seeds in affected ranges unless they would be added back in recomputed seeds. */
    discardSeeds(data.seeds, newSeeds, index, affectedRanges, currNodeSequence, node->identifier, k); 

    /* Add the new seeds to the dynamic seed list and store the index. */
    addSeeds(data.seeds, index, node->identifier, newSeeds);

    /* Recursive step */
    for(Node* child: node->children){

        buildHelper(data, index, T, child, k, s, globalCoords);
    }

    /* Undo seed and sequence mutations when backtracking */
    undoMutations(data, index, T, node, blockMutData, nucMutData);
}

/* Interface implementation */
void pmi::build(seedIndex &index, Tree *T, const size_t k, const size_t s) {

    /* Setup for seed indexing */
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);
    
    /* Get seeds in the consensus MSA */
    std::string consensus = tree::getConsensus(T);
    data.seeds = getUniqueSeeds(consensus, k, s);

    /* Recursive traversal of tree to build the index */
    buildHelper(data, index, T, T->root, k, s, globalCoords);
    index.consensusSeeds = data.seeds;
}


void pmi::write(std::ofstream &fout, Tree *T, seedIndex &index) {
    auto consensusSeeds = index.consensusSeeds;
    for (const kmer_t &s : index.consensusSeeds) {
        fout << s.seq << " ";
    }
    fout << "\n";
    std::stringstream ss;
    writeIndexDFS(ss, consensusSeeds, T->root, index);
    fout << ss.str();
}
void pmi::load(seedIndex &index, const Node *root, const std::ifstream &indexFile) {
}
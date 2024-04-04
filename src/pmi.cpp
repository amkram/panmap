#include "pmi.hpp"
#include "seeding.hpp"
#include <iostream>
#include <sstream>


using namespace seeding;
using namespace pmi;
using namespace PangenomeMAT;
using namespace tree;


/* Helpers */

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
}

void buildHelper(mutableTreeData &data, seedMap_t seedMap, seedIndex &index, Tree *T, const Node *node, const int32_t l, const size_t k, const size_t s, const globalCoords_t &globalCoords) {

    blockMutData_t blockMutData;
    nucMutData_t nucMutData;

    /* Mutate with block and nuc mutations. */
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    tree::updateConsensus(data, T);

    std::set<std::string> outDeletions;
    std::set<std::string> outInsertions;

    nucMutData_t extended = nucMutData;
    for (const auto &blockMut : blockMutData) {
        int32_t blockId = std::get<0>(blockMut);
        int32_t blockStart = getGlobalCoordinate(blockId, 0, -1, globalCoords);
        int32_t seen = 0;
        while (seen < k && blockStart >= 1) {
            if (data.gappedConsensus[blockStart] != '-') {
                seen++;
            }
            blockStart--;
        }

        int32_t blockStop = getGlobalCoordinate(blockId, globalCoords[blockId].first.size()-1, -1, globalCoords);
        bool newStrand = std::get<4>(blockMut);

//        std::cout << "blockMut: " << blockId << "  start: " << blockStart << "  stop: " << blockStop << "  newStrand: " << newStrand << std::endl;
        if (newStrand) {
            extended.push_back(std::make_tuple(-1, -1, -1, -1, -1, blockStart, blockStop - blockStart));
        }
    }
    std::sort(extended.begin(), extended.end(), [](const auto &a, const auto &b) {
        return std::get<5>(a) < std::get<5>(b);
    });

    std::unordered_map<int32_t, bool> seen;
    for (const auto &nucMut : extended) {
        int32_t globalCoord = std::get<5>(nucMut);
        int32_t len = std::get<6>(nucMut);
        int32_t lastSeed = -1;
        for (int32_t c = globalCoord + len; c >= globalCoord; c--) {
            if (data.gappedConsensus[c] == '-') {
                continue;
            }
            std::string kmer = data.ungappedConsensus.substr(data.degap[c], k);
            if (seen.find(c) != seen.end()) {
                continue;
            }
            seen[c] = true;
            if (seedMap.find(c) != seedMap.end()) {
                // This kmer is already a seed.
                std::string prevseedmer = seedMap[c].second;
                if (seeding::is_syncmer(kmer, s, false)) {
                    // Is it still a seed?
                    seedMap[c].second = kmer + seedMap[seedMap[c].first].second.substr(0, (l-1)*k);
                    if (seedMap[c].second == prevseedmer) {
                        continue;
                    }
                    if (seedMap[c].second.size() == l*k) {
                        std::string str = "";
                        str = std::to_string(c);
                        str += ":";
                        str += seedMap[c].second;
                        outInsertions.insert(str);
                    }
                } else {
                    std::string str = "";
                    str += std::to_string(c);
                    str += ":@";
                    str += seedMap[c].second;
                    if (seedMap[c].second.size() == l*k) {
                        outDeletions.insert(str);
                    }
                    seedMap[c].first = -1;
                    seedMap[c].second = "";
                }
            } else {
                // not in seed map, could be a seed now
                if (seeding::is_syncmer(kmer, s, false)) {
                    std::string newseedmer = kmer;
                    if (lastSeed != -1) {
                        newseedmer += seedMap[lastSeed].second.substr(0, (l-1)*k);
                    }
                    if (newseedmer.size() == l*k) {
                        std::string str = "";
                        str += std::to_string(c);
                        str += ":";
                        str += newseedmer;
                        outInsertions.insert(str);
                    }
                    seedMap[c] = std::make_pair(lastSeed, newseedmer);
                    lastSeed = c;
                }
            }
        }
    }
    index.outStream << node->identifier << " ";
    for (const std::string &s : outDeletions) {
        index.outStream << s << " ";
    }
    for (const std::string &s : outInsertions) {
        index.outStream << s << " ";
    }
    index.outStream << "\n";

    /* Recursive step */
    for (Node* child: node->children){
        buildHelper(data, seedMap, index, T, child, l, k, s, globalCoords);
    }
    /* Undo seed and sequence mutations when backtracking */
    undoMutations(data, index, T, node, blockMutData, nucMutData);
}

/* Interface implementation */
void pmi::build(seedIndex &index, Tree *T, const size_t j, const size_t k, const size_t s) {

    /* Setup for seed indexing */
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);
    tree::updateConsensus(data, T);

    seedMap_t seedMap;
    index.outStream << k << " " << s << " " << j << "\n";
 
    /* Recursive traversal of tree to build the index */
    buildHelper(data, seedMap, index, T, T->root, j, k, s, globalCoords);  
}
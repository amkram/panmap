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
        auto gc1 = globalCoords[blockId].first;
        auto bs1 = gc1[gc1.size() - 1].first;
        auto bs2 = gc1[gc1.size() - 1].second;
        int32_t blockStop = 0;
        if (bs2.size() > 0) {
            blockStop = bs2[bs2.size() - 1];
        } else {
            blockStop = bs1;
        }
        affectedRanges.push_back(std::tuple(blockStart, blockStop, blockId));
    }
    return affectedRanges;
}
std::vector<kmer_t> getUniqueSeeds(const std::string &seq, const int32_t k, const int32_t s) {
    std::vector<kmer_t> seedVec = seed::syncmerize(seq, k, s, false, true, 0);
    std::set<kmer_t> seedSet(seedVec.begin(), seedVec.end());
    seedVec.assign(seedSet.begin(), seedSet.end());
    return seedVec;
}
std::vector<kmer_t> getAllSeeds(const std::string &seq, const int32_t k, const int32_t s) {    
    return seed::syncmerize(seq, k, s, false, true, 0);
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
void discardSeeds(std::vector<kmer_t> &seeds, std::vector<kmer_t> &newSeeds, seedIndex &index, const std::vector<range_t>& B, const std::string &seq, const std::string nid, const size_t k) {

    auto cmp = [](kmer_t a, kmer_t b) { return a.idx > b.idx;  };
    std::set<kmer_t, decltype(cmp)> dels(cmp);

    for (int32_t i = 0; i < seeds.size(); i++) {
        kmer_t &s = seeds[i];
        for (int32_t j = 0; j < B.size(); j++) {
            const auto &b = B[j];
            if (std::get<0>(b) > s.pos) {
                break;
            } else if (s.gappedEnd <= std::get<1>(b)) {
                auto it = std::find(newSeeds.begin(), newSeeds.end(), s);
                s.idx = i;
                if (it == newSeeds.end() || it->pos != s.pos) {
                    dels.insert(s);
                } else {
                    newSeeds.erase(it);
                    s.pos = it->pos;
                    s.gappedEnd = it->gappedEnd;
                    s.idx = i;
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

void addSeeds(std::vector<kmer_t> &seeds, seedIndex &index, const std::string nodeId, const std::vector<kmer_t> &newSeeds) {
    for (auto &s : newSeeds) {
        seeds.push_back(s);
        index.insertions[nodeId].push_back(s);
    }
}
void recomputeSeeds(mutableTreeData &data, std::vector<kmer_t> &newSeeds, const std::vector<range_t> &ranges, const std::string &sequence, int32_t k, int32_t s) {
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
            newSeeds.push_back(syncmer);
        }
    }
}

void writeIndexDFS(std::stringstream &ss, std::vector<kmer_t> &seeds, Node *currNode, seedIndex &index) {
    ss << currNode->identifier << "\t";
    std::stack<int32_t> rmDel;
    for (const kmer_t &d : index.deletions[currNode->identifier]) {
        rmDel.push(d.idx);
        ss << d.pos << " " << d.idx << " " << d.gappedEnd << " ";
    }
    removeIndices(seeds, rmDel);
    ss << "^";
    for (const kmer_t &s : index.insertions[currNode->identifier]) {
        seeds.push_back(s);
        ss << s.seq << " " << s.pos << " " << s.idx << " " << s.gappedEnd << " ";
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

void buildHelper(mutableTreeData &data, seedMap_t seedMap, seedIndex &index, Tree *T, const Node *node, const int32_t l, const size_t k, const size_t s, const globalCoords_t &globalCoords) {

    blockMutData_t blockMutData;
    nucMutData_t nucMutData;

    /* Mutate with block and nuc mutations. */
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    tree::updateConsensus(data, T);

    std::unordered_set<std::string> out;
    nucMutData_t extended = nucMutData;
    for (const auto &blockMut : blockMutData) {
        int32_t blockId = std::get<0>(blockMut);
        int32_t start = tree::getGlobalCoordinate(blockId, 0, -1, globalCoords);
        int32_t stop = start + data.sequence[blockId].first.size() - 1;
        bool oldStrand = std::get<2>(blockMut);
        bool newStrand = std::get<4>(blockMut);
        if (newStrand) {
            extended.push_back(std::make_tuple(-1, -1, -1, -1, -1, start, data.sequence[blockId].first.size()));
        }
    }

    for (const auto &nucMut : extended) {
        int32_t globalCoord = std::get<5>(nucMut);
        int32_t len = std::get<6>(nucMut);
        int32_t lastSeed = -1;
        for (int32_t c = globalCoord + len; c >= data.regap[data.degap[globalCoord] - k]; c--) {
            if (c >= data.gappedConsensus.size() || data.gappedConsensus[c] == '-') {
                continue;
            }
            std::string kmer = data.ungappedConsensus.substr(data.degap[c], k);
            if (seedMap.find(c) != seedMap.end()) {
                std::string prevJkmer = seedMap[c].second;
                if (seed::is_syncmer(kmer, k, s)) {
                    seedMap[c].second = kmer + seedMap[seedMap[c].first].second.substr(0, (l-1)*k);
                    
                    if (seedMap[c].second == prevJkmer) {
                        continue;
                    }
                    if (seedMap[c].second.size() == l*k) {
                        std::string s = "";
                        s += std::to_string(c);
                        s += ":@";
                        s += prevJkmer;
                        out.insert(s);
                        s = std::to_string(c);
                        s += ":";
                        s += seedMap[c].second;
                        out.insert(s);
                    }
                } else {
                    std::string s = "";
                    s += std::to_string(c);
                    s += ":@";
                    s += seedMap[c].second;
                    out.insert(s);
                    seedMap[c].first = -1;
                    seedMap[c].second = "";
                }
            } else {
                // not in seed map, could be a seed now
                if (seed::is_syncmer(kmer, k, s)) {
                    std::string newJkmer = kmer;
                    if (lastSeed != -1) {
                        newJkmer += seedMap[lastSeed].second.substr(0, (l-1)*k);
                    }
                    if (newJkmer.size() == l*k) {
                        std::string s = "";
                        s += std::to_string(c);
                        s += ":";
                        s += newJkmer;
                        out.insert(s);
                    }
                    seedMap[c] = std::make_pair(lastSeed, newJkmer);
                    lastSeed = c;
                }
            }
        }
    }
    index.outStream << node->identifier << " ";
    for (const std::string &s : out) {
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
void pmi::build(seedIndex &index, Tree *T, const size_t l, const size_t k, const size_t s) {

    /* Setup for seed indexing */
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);
    tree::updateConsensus(data, T);

    /*  seed and j-k-mer tracking: 
    **     global coord => {next seed, jkmer seq}
    */
    seedMap_t seedMap;
    // std::vector<kmer_t> consensusSeeds = getAllSeeds(data.gappedConsensus, k, s);
    // //jkmers
    // std::vector<jkmer> jkmers = seed::jkmerize(consensusSeeds, 3);
    // for (const jkmer &j : jkmers) {
    //     seedMap[j.positions[0]] = std::make_pair(j.positions[1], j.seq);
    // }
    //seedmap
    // for (const auto &p : seedMap) {
    //     std::cout << p.first << " " << p.second.first << " " << p.second.second << "\n";
    // }
    index.outStream << k << " " << s << "\n";
 
    /* Recursive traversal of tree to build the index */
    buildHelper(data, seedMap, index, T, T->root, l, k, s, globalCoords);
   
}
void pmi::write(std::ofstream &fout, Tree *T, seedIndex &index) {
    fout << index.k << " " << index.s << "\n";
    auto consensusSeeds = index.consensusSeeds;
    for (const kmer_t &s : index.consensusSeeds) {
        fout << s.seq << " " << s.pos << " " << s.gappedEnd << " ";
    }
    fout << "\n";
    std::stringstream ss;
    writeIndexDFS(ss, consensusSeeds, T->root, index);
    fout << ss.str();
}
void pmi::load(seedIndex &index, std::ifstream &indexFile) {
    std::string line;

    // parse syncmer parameters
    int32_t k;
    int32_t s;
    std::getline(indexFile, line); // k s 
    std::istringstream iss(line);
    iss >> k;
    iss >> s;

    // parse consensus seeds
    std::vector<kmer_t> consensusSeeds;
    std::getline(indexFile, line);
    std::vector< std::string > rootSplt;
    PangenomeMAT::stringSplit(line, ' ', rootSplt);

    for (int32_t i = 0; i < rootSplt.size() - 3 + 1; i += 3) {
        consensusSeeds.push_back(kmer_t{rootSplt[i], std::stoi(rootSplt[i+1]), i/3, -1, false, std::stoi(rootSplt[i+2])});
    }

    // parse each node's seed insertions and deletions
    while (std::getline(indexFile, line)) {
        if (line.size() < 2) {
            continue;
        }
        std::vector<std::string> spltNid;
        std::vector<std::string> spltParts;
        std::vector<std::string> spltDel = {};
        std::vector<std::string> spltIns = {};
        stringSplit(line, '\t', spltNid);
        std::string nodeId = spltNid[0];
        stringSplit(spltNid[1], '^', spltParts);
        stringSplit(spltParts[0], ' ', spltDel);
        if (spltParts.size() > 1) {
            stringSplit(spltParts[1], ' ', spltIns);
        }
        std::vector<kmer_t> deletions;
        if (spltDel.size() > 0) {
            for (int32_t i = 0; i < spltDel.size() - 3 + 1; i += 3) {
                deletions.push_back(kmer_t{"", std::stoi(spltDel[i]), std::stoi(spltDel[i+1]), -1, false, std::stoi(spltDel[i+2])});
            }
        }
        std::vector<kmer_t> insertions;
        if (spltIns.size() > 0) {
            for (int32_t i = 0; i < spltIns.size() - 4 + 1; i += 4) {
                insertions.push_back(kmer_t{spltIns[i], std::stoi(spltIns[i+1]), std::stoi(spltIns[i+2]), -1, false, std::stoi(spltIns[i+3])});
            }
        }

        // load insertions and deletions into node maps
        index.insertions[nodeId] = insertions;
        index.deletions[nodeId] = deletions;
    }
    
    index.consensusSeeds = consensusSeeds;
    index.k = k;
    index.s = s;
}
#include <algorithm>
#include <cassert>
#include <deque>
#include <queue>
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"


void mgsr::accio(PangenomeMAT::Tree *T, std::ifstream& indexFile, size_t k, size_t l) {
    std::cout << "What's my purpose\nYou pass butter" << std::endl;
}


/*testing*/
void initializeMap(mgsr::seedmers& seedmersIndex, const std::vector<std::tuple<hash_t, int32_t, int32_t, bool>>& syncmers, std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, const int32_t l, const int32_t k) {
    // first seedmer
    int startIdx = 0;
    size_t cacheForwardH  = 0;
    size_t cacheReversedH = 0;
    size_t cacheMin;
    bool   rev;
    while (seedmersIndex.positionMap.empty()) {
        assert(std::get<3>(syncmers[startIdx]) == false);
        curSeeds[std::get<1>(syncmers[startIdx])] = std::make_pair(std::get<2>(syncmers[startIdx]),std::get<0>(syncmers[startIdx]));
        
        cacheForwardH  = 0;
        cacheReversedH = 0;
        for (int i = startIdx; i < startIdx + l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(syncmers[i]);
        for (int i = startIdx + l - 1; i > startIdx - 1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(syncmers[i]);
    
        // skip direction ambiguous
        if (cacheForwardH < cacheReversedH) {
            cacheMin = cacheForwardH;
            rev = false;
        } else if (cacheReversedH < cacheForwardH) {
            cacheMin = cacheReversedH;
            rev = true;
        } else {
            ++startIdx;
            continue;
        }

        seedmersIndex.positionMap[std::get<1>(syncmers[startIdx])] = std::make_pair(cacheMin, rev);
        seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(syncmers[startIdx]));
        ++startIdx;
    }

    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    // rest of seedmers
    for (int i = startIdx; i < syncmers.size() - l + 1; ++i) {
        assert(std::get<3>(syncmers[i]) == false);

        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<2>(syncmers[i]),std::get<0>(syncmers[i]));

        cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);
        cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(syncmers[i+l-1]) << (2 * k * (l - 1)));

        // skip direction ambiguous
        if (cacheForwardH < cacheReversedH) {
            cacheMin = cacheForwardH;
            rev = false;
        } else if (cacheReversedH < cacheForwardH) {
            cacheMin = cacheReversedH;
            rev = true;
        } else {
            continue;
        }

        seedmersIndex.positionMap[std::get<1>(syncmers[i])] = std::make_pair(cacheMin, rev);
        seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(syncmers[i]));
    }

    // rest of seeds
    for (int i = syncmers.size() - l + 1; i < syncmers.size(); i++) {
        assert(std::get<3>(syncmers[i]) == false);
        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<2>(syncmers[i]),std::get<0>(syncmers[i]));
    }
}

void updateCurSeeds(std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges) {
    for (const auto& change : syncmerChanges) {
        auto [seedHash, seedBeg, seedEnd, seedDel] = change;
        if (seedDel) {
            assert(curSeeds.erase(seedBeg));
        } else {
            curSeeds[seedBeg] = std::make_pair(seedEnd, seedHash);
        }
    }
}

static size_t btn(char b) {
    size_t n;
    switch (b) {
    case 'A':
        n = 0;
        break;
    case 'C':
        n = 1;
        break;
    case 'G':
        n = 2;
        break;
    case 'T':
        n = 3;
        break;
    default:
        throw std::invalid_argument("Kmer contains non canonical base");
        break;
    }
    return n;
}

static size_t hash(const std::string& s) {
    size_t h = 0;
    if (s.empty()) {
        return h;
    } else if (s.size() == 1) {
        h = btn(s[0]);
        return h;
    }

    h = btn(s[0]);
    for (size_t i = 1; i < s.size(); ++i) {
        h = (h << 2) + btn(s[i]);
    }
    return h;
}

static char comp(char c) {
    char compC;
    switch (c) { 
    case 'A':
        compC = 'T';
        break;
    case 'C':
        compC = 'G';
        break;
    case 'G':
        compC = 'C';
        break;
    case 'T':
        compC = 'A';
        break;
    default:
        compC = 'N';
        break;
    }
    return compC;
}

static std::string revcomp(const std::string& s) {
    std::string cs = "";
    for (int i = s.size() - 1; i > -1; --i) {
        char c = s[i];
        cs += comp(c);
    }
    return cs;
}

static std::pair<size_t, bool> getHash(const std::string& s) {
    try {
        size_t u = hash(s);
        size_t v = hash(revcomp(s));
        if (u < v) return std::make_pair(u, true);
        else if (v < u) return std::make_pair(v, true);
        return std::make_pair(0, false);
    } catch (std::invalid_argument) {
        return std::make_pair(0, false); // skip strand ambiguous
    }
    return std::make_pair(0, false);
}
/*testing*/

void writeCurSeeds(const std::string& path, const std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, const tree::mutableTreeData& data) {
        std::ofstream ofSeeds(path);
        for (const auto& seed : curSeeds) {
            ofSeeds << seed.second.second << "\t"
                    << data.degap[seed.first] << "\t"
                    << data.degap[seed.second.first] << "\t"
                    << seed.first << "\t"
                    << seed.second.first << "\n";
        }
        ofSeeds.close();
}

void writeCurSeedmers(const std::string& path, const mgsr::seedmers& seedmersIndex, const tree::mutableTreeData& data) {
    std::ofstream ofSeedmers(path);
    for (const auto& position : seedmersIndex.positionMap) {
        auto beg = position.first;
        auto hash = position.second.first;
        auto rev = position.second.second;
        ofSeedmers << hash << "\t"
                   << seedmersIndex.seedmersMap.find(hash)->second.size() << "\t"
                   << data.degap[beg] << "\t"
                   << beg << "\t"
                   << rev << std::endl;

    }
    ofSeedmers.close();
}

void mgsr::buildSeedmerHelper(tree::mutableTreeData &data, seedMap_t seedMap, pmi::seedIndex &index, mgsr::seedmers seedmersIndex, std::map<int32_t, std::pair<int32_t, size_t>> curSeeds, PangenomeMAT::Tree *T, const PangenomeMAT::Node *node, const int32_t l, const int32_t k, const int32_t s, const tree::globalCoords_t &globalCoords, std::stringstream& seedmersOutStream) {
    tree::blockMutData_t blockMutData;
    tree::nucMutData_t nucMutData;

    /* Mutate with block and nuc mutations. */
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    tree::updateConsensus(data, T);

    std::set<std::string> outDeletions;
    std::set<std::string> outInsertions;

    /*testing*/
    // hash, beg, end, del
    std::vector<std::tuple<size_t, int32_t, int32_t, bool>> syncmerChanges;
    /*testing*/


    tree::nucMutData_t extended = nucMutData;
    for (const auto &blockMut : blockMutData) {
        int32_t blockId = std::get<0>(blockMut);
        int32_t blockStart = tree::getGlobalCoordinate(blockId, 0, -1, globalCoords);
        int32_t seen = 0;
        while (seen < 2 * k && blockStart >= 1) {
            if (data.gappedConsensus[blockStart] != '-') {
                seen++;
            }
            blockStart--;
        }

        int32_t blockStop = tree::getGlobalCoordinate(blockId, globalCoords[blockId].first.size()-1, -1, globalCoords);
        bool newStrand = std::get<4>(blockMut);
        std::cerr << "bloackStart? " << blockStart << std::endl;
//        std::cout << "blockMut: " << blockId << "  start: " << blockStart << "  stop: " << blockStop << "  newStrand: " << newStrand << std::endl;
        if (newStrand) {
            extended.push_back(std::make_tuple(-1, -1, -1, -1, -1, blockStart, blockStop - blockStart));
        }
    }
    std::sort(extended.begin(), extended.end(), [](const auto &a, const auto &b) {
        return std::get<5>(a) < std::get<5>(b);
    });

    std::unordered_map<int32_t, bool> seen;
    std::unordered_map<int32_t, bool> seenNucMut;
    for (const auto &nucMut : extended) {
        int32_t globalCoord = std::get<5>(nucMut);
        
        int32_t len = std::get<6>(nucMut) - 1;
        // if (std::get<3>(nucMut) != '-' && std::get<4>(nucMut) != '-') {
        //     std::cerr << "SUBS: alt len? " << data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord << std::endl;
        // }
        if (std::get<4>(nucMut) != -1 && data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord > 0) {
            // insertion
            len = data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord;
        }
        // else if (std::get<4>(nucMut) == '-') {
        //     // deletion
        //     std::cerr << "for deletion: len = " << data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord << "?" << std::endl;
        // } 
        int32_t lastSeed = -1;

        // if (seenNucMut.find(globalCoord) != seenNucMut.end()) {
        //     continue;
        // }
        // seenNucMut[globalCoord] = true;

        std::cerr << std::get<0>(nucMut) << "\t"
                  << std::get<1>(nucMut) << "\t"
                  << std::get<2>(nucMut) << "\t"
                  << std::get<3>(nucMut) << "\t"
                  << std::get<4>(nucMut) << "\t"
                  << std::get<5>(nucMut) << "\t"
                  << std::get<6>(nucMut) << std::endl;
        std::cerr << "globalCoord: " << globalCoord << "\tnucMut og: " << std::get<6>(nucMut) << "\tlen: " << len << "\taltlen: " <<  data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord << std::endl;
        std::cerr << "data.degap[globalCoord] = " << data.degap[globalCoord] << "\tdata.degap[globalCoord] - k = " << data.degap[globalCoord] - static_cast<int32_t>(k) << std::endl;
        std::cerr << "data.regap[std::max(0, data.degap[globalCoord] - static_cast<int32_t>(k))] = " << data.regap[std::max(0, data.degap[globalCoord] - static_cast<int32_t>(k))] << std::endl;
        
        int32_t leftEndCoor = std::max(0, data.regap[std::max(0, data.degap[globalCoord] - static_cast<int32_t>(k))]);
        if (leftEndCoor > globalCoord) {
            leftEndCoor = globalCoord;
            // len = data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord;

        }
        
        std::cerr << "old from " << globalCoord + len << " to " << std::max(0, data.regap[std::max(0, data.degap[globalCoord] - static_cast<int32_t>(k))]) << std::endl;
        std::cerr << "new from " << globalCoord + len << " to " << leftEndCoor << std::endl;
        for (int32_t c = globalCoord + len; c >= leftEndCoor; c--) {
            if (seen.find(c) != seen.end()) {
                continue;
            }
            seen[c] = true;

            if (data.gappedConsensus[c] == '-') {
                if (seedMap.find(c) != seedMap.end()) {
                    // no longer a seed -> delete
                    std::string str = "";
                    str += std::to_string(c);
                    str += ":@";
                    str += seedMap[c].second;
                    if (seedMap[c].second.size() == k) {
                        outDeletions.insert(str);
                        /*testing*/
                        syncmerChanges.emplace_back(std::make_tuple((getHash(seedMap[c].second)).first, c, data.regap[data.degap[c] + k - 1], true));
                        /*testing*/
                    }
                    seedMap[c].first = -1;
                    seedMap[c].second = "";
                }
                continue;
            }
            std::string kmer = data.ungappedConsensus.substr(data.degap[c], k);

            /*testing*/
            std::pair<size_t, bool> kmerHash = getHash(kmer);
            size_t curHash = kmerHash.first;
            /*testing*/

            if (seedMap.find(c) != seedMap.end()) {
                // This kmer is already a seed.
                std::string prevseedmer = seedMap[c].second;
                if (kmerHash.second && kmer.size() == k && seeding::is_syncmer(kmer, s, false)) {
                    // Is it still a seed?
                    seedMap[c].second = kmer;
                    if (seedMap[c].second == prevseedmer) {
                        continue;
                    }
                    if (seedMap[c].second.size() == k) {
                        std::string str = "";
                        str = std::to_string(c);
                        str += ":";
                        str += seedMap[c].second;
                        outInsertions.insert(str);
                        /*testing*/
                        syncmerChanges.emplace_back(std::make_tuple(curHash, c, data.regap[data.degap[c] + k - 1], false));
                        /*testing*/
                    }
                } else {
                    std::string str = "";
                    str += std::to_string(c);
                    str += ":@";
                    str += seedMap[c].second;
                    if (seedMap[c].second.size() == k) {
                        outDeletions.insert(str);
                        /*testing*/
                        syncmerChanges.emplace_back(std::make_tuple((getHash(seedMap[c].second)).first, c, data.regap[data.degap[c] + k - 1], true));
                        /*testing*/
                    }
                    seedMap[c].first = -1;
                    seedMap[c].second = "";
                }
            } else {
                // not in seed map, could be a seed now
                if (seeding::is_syncmer(kmer, s, false) && kmerHash.second) {
                    std::string newseedmer = kmer;
                    if (newseedmer.size() == k) {
                        std::string str = "";
                        str += std::to_string(c);
                        str += ":";
                        str += newseedmer;
                        outInsertions.insert(str);
                        /*testing*/
                        syncmerChanges.emplace_back(std::make_tuple(curHash, c, data.regap[data.degap[c] + k - 1], false));
                        /*testing*/
                    }
                    seedMap[c] = std::make_pair(lastSeed, newseedmer);
                    lastSeed = c;
                }
            }
        }
    }

    /*testing--*/
    std::sort(syncmerChanges.begin(), syncmerChanges.end(), [](const auto &a, const auto &b) {
        return std::get<1>(a) < std::get<1>(b);
    });
    if (seedmersIndex.seedmersMap.size() == 0) {
        /*build seedmer for root*/
        std::cerr << node->identifier << std::endl;
        initializeMap(seedmersIndex, syncmerChanges, curSeeds, l, k);

        std::string seedsPath = "../src/test/data/mgsr/test/" + node->identifier + ".smi";
        writeCurSeeds(seedsPath, curSeeds, data);
        // std::string seedmersPath = "../src/test/data/mgsr/test/" + node->identifier + ".kmi";
        // writeCurSeedmers(seedmersPath, seedmersIndex, data);
    } else {
        std::cerr << node->identifier << std::endl;
        updateCurSeeds(curSeeds, syncmerChanges);

        std::string seedsPath = "../src/test/data/mgsr/test/" + node->identifier + ".smi";
        writeCurSeeds(seedsPath, curSeeds, data);
        
        std::unordered_set<int32_t> seenBegs;

        size_t  cacheReversedH;
        size_t  cacheForwardH;
        size_t  cacheMin;
        bool    rev;
        int     count;
        for (const auto& change : syncmerChanges) {
            auto [seedHash, seedBeg, seedEnd, seedDel] = change;
            // std::cerr << "cur syncmer change: " << seedHash << "\t" << seedBeg << "\t" << seedDel << std::endl;
        
            decltype(curSeeds.begin()) curSeedIt;
            if (seedDel) curSeedIt = curSeeds.lower_bound(seedBeg);
            else         curSeedIt = curSeeds.find(seedBeg);

            auto localStartSeed = curSeedIt;
            auto localEndSeed   = curSeedIt;
            int offset = 0;
            while (offset < l - 1 && localStartSeed != curSeeds.begin()) {
                --localStartSeed;
                ++offset;
            }

            auto curAffectedSeed = localStartSeed;
            if (!seedDel) ++localEndSeed; 
            // std::cerr << "localStartSeed: " << localStartSeed->first << "\tlocalEndSeed: " << localEndSeed->first << std::endl;
            while (curAffectedSeed != localEndSeed) {
                // std::cerr << "Inside while loop curAffectedSeed: " << curAffectedSeed->first << "\tlocalEndSeed: " << localEndSeed->first << std::endl;
                if (seenBegs.find(curAffectedSeed->first) != seenBegs.end()) {
                    ++curAffectedSeed;
                    continue;
                }
                seenBegs.insert(curAffectedSeed->first);


                count = 0;
                cacheForwardH  = 0;
                cacheReversedH = 0;
                auto curLocalSeed = curAffectedSeed;
                while (count < l && curLocalSeed != curSeeds.end()) {
                    cacheForwardH  = (cacheForwardH << (2 * k)) + curLocalSeed->second.second;
                    cacheReversedH = cacheReversedH + (curLocalSeed->second.second << (2 * k * count));
                    ++curLocalSeed;
                    ++count;
                }

                if (count != l) {
                    auto positionMapIt = seedmersIndex.positionMap.find(curAffectedSeed->first);
                    if (positionMapIt != seedmersIndex.positionMap.end()) {
                        auto begToErase  = curAffectedSeed->first;
                        auto hashToErase = positionMapIt->second.first;
                        assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
                        if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
                        assert(seedmersIndex.positionMap.erase(begToErase));
                    }
                    break;
                }

                if (cacheForwardH < cacheReversedH) {
                    cacheMin = cacheForwardH;
                    rev = false;
                } else if (cacheReversedH < cacheForwardH) {
                    cacheMin = cacheReversedH;
                    rev = true;
                } else {
                    auto positionMapIt = seedmersIndex.positionMap.find(curAffectedSeed->first);
                    if (positionMapIt != seedmersIndex.positionMap.end()) {
                        auto begToErase  = curAffectedSeed->first;
                        auto hashToErase = positionMapIt->second.first;
                        assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
                        if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
                        assert(seedmersIndex.positionMap.erase(begToErase));
                    }
                    ++curAffectedSeed;
                    continue;
                }

                if (node->identifier == "node_2") {
                    std::cout << cacheMin << "\t" << curAffectedSeed->first << std::endl;
                }
                
                auto potentialDel = seedmersIndex.positionMap.find(curAffectedSeed->first);
                if (potentialDel != seedmersIndex.positionMap.end()) {
                    seedmersIndex.seedmersMap.find(potentialDel->second.first)->second.erase(potentialDel->first);
                    if (seedmersIndex.seedmersMap[potentialDel->second.first].empty()) seedmersIndex.seedmersMap.erase(potentialDel->second.first);
                }
                seedmersIndex.positionMap[curAffectedSeed->first] = std::make_pair(cacheMin, rev);
                seedmersIndex.seedmersMap[cacheMin].insert(curAffectedSeed->first);
                ++curAffectedSeed;
            }

            // std::cerr << "checking if delete" << std::endl;
            if (seedDel) {
                if (node->identifier == "node_2") std::cout << "delete: " << seedHash << "\t" << seedBeg << std::endl;
                auto positionMapIt = seedmersIndex.positionMap.find(seedBeg);
                if (positionMapIt != seedmersIndex.positionMap.end()) {
                    auto begToErase  = seedBeg;
                    auto hashToErase = positionMapIt->second.first;
                    assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
                    if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
                    assert(seedmersIndex.positionMap.erase(begToErase));
                }
            }
        }
        std::string seedmersPath = "../src/test/data/mgsr/test/" + node->identifier + ".kmi";
        writeCurSeedmers(seedmersPath, seedmersIndex, data);
    }

    /*--testing*/
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
        // if (child->identifier == "node_3" || child->identifier == "node_502" || child->identifier == "node_1002") {continue;}
        mgsr::buildSeedmerHelper(data, seedMap, index, seedmersIndex, curSeeds, T, child, l, k, s, globalCoords, seedmersOutStream);
    }
    /* Undo seed and sequence mutations when backtracking */
    undoMutations(data, index, T, node, blockMutData, nucMutData);
}

void mgsr::buildSeedmer(pmi::seedIndex &index, PangenomeMAT::Tree *T, const size_t l, const size_t k, const size_t s, std::stringstream& seedmersOutStream) {
    /* Setup for seed indexing */
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);
    tree::updateConsensus(data, T);

    seedMap_t seedMap;
    mgsr::seedmers seedmersIndex;
    std::map<int32_t, std::pair<int32_t, size_t>> curSeeds; // map[beg] = {end, hash}
    index.outStream << k << " " << s << " " << l << "\n";
    seedmersOutStream << k << " " << s << " " << l << "\n";
    /* Recursive traversal of tree to build the index */
    mgsr::buildSeedmerHelper(data, seedMap, index, seedmersIndex, curSeeds, T, T->root, l, k, s, globalCoords, seedmersOutStream);
}
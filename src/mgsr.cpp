#include <algorithm>
#include <iterator>
#include <cassert>
#include <deque>
#include <queue>
#include <math.h>
#include <atomic>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_reduce.h>
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"
#include <eigen3/Eigen/Dense>

namespace fs = boost::filesystem;

/*testing*/
void initializeMap(mgsr::seedmers& seedmersIndex, const std::vector<std::tuple<hash_t, int32_t, int32_t, bool>>& syncmers, std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, const int32_t l, const int32_t k, std::stringstream& seedmersOutStream) {
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

        // seedmersIndex.positionMap[std::get<1>(syncmers[startIdx])] = std::make_pair(cacheMin, rev);
        seedmersIndex.positionMap[std::get<1>(syncmers[startIdx])] = std::make_tuple(std::get<2>(syncmers[startIdx+l-1]), cacheMin, rev);
        seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(syncmers[startIdx]));
        seedmersOutStream << std::get<1>(syncmers[startIdx]) << ":+:"
                          << std::get<2>(syncmers[startIdx+l-1]) << ","
                          << cacheMin << ","
                          << rev << " ";
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

        // seedmersIndex.positionMap[std::get<1>(syncmers[i])] = std::make_pair(cacheMin, rev);
        seedmersIndex.positionMap[std::get<1>(syncmers[i])] = std::make_tuple(std::get<2>(syncmers[i+l-1]), cacheMin, rev);
        seedmersIndex.seedmersMap[cacheMin].insert(std::get<1>(syncmers[i]));
        seedmersOutStream << std::get<1>(syncmers[i]) << ":+:"
                          << std::get<2>(syncmers[i+l-1]) << ","
                          << cacheMin << ","
                          << rev << " ";
    }

    // rest of seeds
    for (int i = syncmers.size() - l + 1; i < syncmers.size(); i++) {
        assert(std::get<3>(syncmers[i]) == false);
        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<2>(syncmers[i]),std::get<0>(syncmers[i]));
    }
}

void updateCurSeeds(std::map<int32_t, std::pair<int32_t, size_t>>& curSeeds, const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges, const std::unordered_map<int32_t, int32_t>& syncmerGlobalEndCoorChanges) {
    for (const auto& change : syncmerChanges) {
        auto [seedHash, seedBeg, seedEnd, seedDel] = change;
        if (seedDel) {
            assert(curSeeds.erase(seedBeg));
        } else {
            curSeeds[seedBeg] = std::make_pair(seedEnd, seedHash);
        }
    }

    for (const auto& endCoorChange : syncmerGlobalEndCoorChanges) {
        curSeeds[endCoorChange.first].first = endCoorChange.second;
    }
}

static size_t btn(char b) {
    size_t n;
    switch (b) {
    case 'A': n = 0; break;
    case 'a': n = 0; break;
    case 'C': n = 1; break;
    case 'c': n = 1; break;
    case 'G': n = 2; break;
    case 'g': n = 2; break;
    case 'T': n = 3; break;
    case 't': n = 3; break;
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
    case 'A': compC = 'T'; break;
    case 'a': compC = 't'; break;
    case 'C': compC = 'G'; break;
    case 'c': compC = 'g'; break;
    case 'G': compC = 'C'; break;
    case 'g': compC = 'c'; break;
    case 'T': compC = 'A'; break;
    case 't': compC = 'a'; break;
    default:  compC = 'N'; break;
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
        auto [end, hash, rev] = position.second;

        ofSeedmers << hash << "\t"
                   << seedmersIndex.seedmersMap.find(hash)->second.size() << "\t"
                   << data.degap[beg] << "\t"
                   << data.degap[end] << "\t"
                   << beg << "\t"
                   << end << "\t"
                   << rev << std::endl;

    }
    ofSeedmers.close();
}

void writeCurSeedmers(const std::string& path, const mgsr::seedmers& seedmersIndex) {
    std::ofstream ofSeedmers(path);
    for (const auto& position : seedmersIndex.positionMap) {
        auto beg = position.first;
        auto [end, hash, rev] = position.second;

        ofSeedmers << hash << "\t"
                   << seedmersIndex.seedmersMap.find(hash)->second.size() << "\t"
                   << beg << "\t"
                   << end << "\t"
                   << rev << std::endl;

    }
    ofSeedmers.close();
}

void mgsr::buildSeedmerHelper(tree::mutableTreeData &data, seedMap_t seedMap, pmi::seedIndex &index, mgsr::seedmers seedmersIndex, std::map<int32_t, std::pair<int32_t, size_t>> curSeeds, PangenomeMAT::Tree *T, const PangenomeMAT::Node *node, const int32_t l, const int32_t k, const int32_t s, const tree::globalCoords_t &globalCoords, std::stringstream& seedmersOutStream) {
    tree::blockMutData_t blockMutData;
    tree::nucMutData_t nucMutData;

    /* Mutate with block and nuc mutations. */
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    std::map<int32_t, int32_t> coordsIndex;
    tree::updateConsensus(data, T, &coordsIndex);

    std::set<std::string> outDeletions;
    std::set<std::string> outInsertions;

    /*testing*/
    // hash, beg, end, del
    std::vector<std::tuple<size_t, int32_t, int32_t, bool>> syncmerChanges;
    std::unordered_map<int32_t, int32_t> syncmerGlobalEndCoorChanges;
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
        // std::cerr << std::get<0>(nucMut) << "\t"
        //           << std::get<1>(nucMut) << "\t"
        //           << std::get<2>(nucMut) << "\t"
        //           << std::get<3>(nucMut) << "\t"
        //           << std::get<4>(nucMut) << "\t"
        //           << std::get<5>(nucMut) << "\t"
        //           << std::get<6>(nucMut) << std::endl;

        int32_t globalCoord = std::get<5>(nucMut);

        if (std::get<4>(nucMut) != -1) {
            if (seenNucMut.find(globalCoord) != seenNucMut.end()) {
                continue;
            }
            seenNucMut[globalCoord] = true;
        }

        int32_t len;

        if (std::get<4>(nucMut) != -1 && data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord > 0) len = data.regap[data.degap[globalCoord] + std::get<6>(nucMut) - 1] - globalCoord;
        else len = std::get<6>(nucMut) - 1;

        int32_t leftEndCoor = std::max(0, data.regap[std::max(0, data.degap[globalCoord] - static_cast<int32_t>(k))]);
        if (leftEndCoor > globalCoord) leftEndCoor = globalCoord;

        int32_t lastSeed = -1;    
        
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
                        syncmerGlobalEndCoorChanges[c] = data.regap[data.degap[c] + k - 1];
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
        // std::cerr << "initializing" << node->identifier << std::endl;
        seedmersOutStream << node->identifier << ":" << data.ungappedConsensus.size() << " ";
        initializeMap(seedmersIndex, syncmerChanges, curSeeds, l, k, seedmersOutStream);
        
        seedmersOutStream << "c:";
        for (const auto& coord : coordsIndex) {
            seedmersOutStream << coord.first << "," << coord.second << ",";
        }
        seedmersOutStream << "\n";

        // std::string dirName = "../src/test/data/mgsr/test/k" + std::to_string(k) + "_s" + std::to_string(s) + "_l" + std::to_string(l);
        // fs::create_directories(dirName);
        // std::string seedsPath = dirName + "/" + node->identifier + ".smi";
        // std::string seedmersPath = dirName + "/" + node->identifier + ".kmi";
        // writeCurSeeds(seedsPath, curSeeds, data);
        // writeCurSeedmers(seedmersPath, seedmersIndex, data);
    } else {
        // std::cerr << node->identifier << std::endl;
        updateCurSeeds(curSeeds, syncmerChanges, syncmerGlobalEndCoorChanges);

        seedmersOutStream << node->identifier << ":" << data.ungappedConsensus.size() << " ";

        std::unordered_set<int32_t> seenBegs;
        size_t  cacheReversedH;
        size_t  cacheForwardH;
        size_t  cacheMin;
        bool    rev;
        int     count;
        for (const auto& change : syncmerChanges) {
            auto [seedHash, seedBeg, seedEnd, seedDel] = change;
        
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
            while (curAffectedSeed != localEndSeed) {
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
                --curLocalSeed;
                
                if (count != l) {
                    auto positionMapIt = seedmersIndex.positionMap.find(curAffectedSeed->first);
                    if (positionMapIt != seedmersIndex.positionMap.end()) {
                        auto begToErase  = curAffectedSeed->first;
                        auto hashToErase = std::get<1>(positionMapIt->second);
                        assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
                        if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
                        assert(seedmersIndex.positionMap.erase(begToErase));
                        seedmersOutStream << begToErase << ":-:" << hashToErase << " ";
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
                        auto hashToErase = std::get<1>(positionMapIt->second);
                        assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
                        if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
                        assert(seedmersIndex.positionMap.erase(begToErase));
                        seedmersOutStream << begToErase << ":-:" << hashToErase << " ";
                    }
                    ++curAffectedSeed;
                    continue;
                }
                
                auto potentialDel = seedmersIndex.positionMap.find(curAffectedSeed->first);
                if (potentialDel != seedmersIndex.positionMap.end()) {
                    seedmersIndex.seedmersMap.find(std::get<1>(potentialDel->second))->second.erase(potentialDel->first);
                    if (seedmersIndex.seedmersMap[std::get<1>(potentialDel->second)].empty()) seedmersIndex.seedmersMap.erase(std::get<1>(potentialDel->second));
                }
                
                seedmersIndex.positionMap[curAffectedSeed->first] = std::make_tuple(curLocalSeed->second.first, cacheMin, rev);
                seedmersIndex.seedmersMap[cacheMin].insert(curAffectedSeed->first);
                seedmersOutStream << curAffectedSeed->first << ":+:"
                                  << curLocalSeed->second.first << ","
                                  << cacheMin << ","
                                  << rev << " ";
                ++curAffectedSeed;
            }

            if (seedDel) {
                auto positionMapIt = seedmersIndex.positionMap.find(seedBeg);
                if (positionMapIt != seedmersIndex.positionMap.end()) {
                    auto begToErase  = seedBeg;
                    auto hashToErase = std::get<1>(positionMapIt->second);
                    assert(seedmersIndex.seedmersMap[hashToErase].erase(begToErase));
                    if (seedmersIndex.seedmersMap[hashToErase].empty()) seedmersIndex.seedmersMap.erase(hashToErase);
                    assert(seedmersIndex.positionMap.erase(begToErase));
                    seedmersOutStream << begToErase << ":-:" << hashToErase << " ";
                }
            }
        }

        for (const auto& endCoorChange : syncmerGlobalEndCoorChanges) {
            auto curSeedIt = curSeeds.find(endCoorChange.first);
            int count = 0;
            while (count < l - 1) {
                if (curSeedIt == curSeeds.begin()) break;
                ++count;
                --curSeedIt;
            }
            if (count != l - 1) continue;
            auto curPositionIt = seedmersIndex.positionMap.find(curSeedIt->first);
            if (curPositionIt == seedmersIndex.positionMap.end()) continue;

            if (std::get<0>(curPositionIt->second) != endCoorChange.second) {
                curPositionIt->second = std::make_tuple(endCoorChange.second, std::get<1>(curPositionIt->second), std::get<2>(curPositionIt->second));
                seedmersOutStream << curPositionIt->first << ":e:" << endCoorChange.second << " ";
            }
        }

        seedmersOutStream << "c:";
        for (const auto& coord : coordsIndex) {
            seedmersOutStream << coord.first << "," << coord.second << ",";
        }
        seedmersOutStream << "\n";

        
        // std::string dirName = "../src/test/data/mgsr/test/k" + std::to_string(k) + "_s" + std::to_string(s) + "_l" + std::to_string(l);
        // fs::create_directories(dirName);
        // std::string seedsPath = dirName + "/" + node->identifier + ".smi";
        // std::string seedmersPath = dirName + "/" + node->identifier + ".kmi";
        // writeCurSeeds(seedsPath, curSeeds, data);
        // writeCurSeedmers(seedmersPath, seedmersIndex, data);
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

std::vector<std::tuple<size_t, int32_t, int32_t>> syncmersSketch(const std::string& seq, const int k, const int s, const bool open) {
    std::vector<std::tuple<size_t, int32_t, int32_t>> syncmers;
    for (size_t i = 0; i < seq.size() - k + 1; ++i) {
        std::string kmer = seq.substr(i, k);
        if (!seeding::is_syncmer(kmer, s, open)) continue;

        std::pair<size_t, bool> kmerHash = getHash(kmer);
        if (!kmerHash.second) continue;

        syncmers.emplace_back(std::make_tuple(kmerHash.first, i, i + k - 1));
    }
    return syncmers;
}
readSeedmers_t extractKminmers(const std::vector<std::tuple<size_t, int32_t, int32_t>>& syncmers, const int k, const int l) {
    std::vector<std::tuple<size_t, int32_t, int32_t, bool, int>> kminmers;
    std::unordered_set<size_t> hashes;

    // first kminmer
    size_t cacheForwardH = 0;
    for (int i = 0; i < l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(syncmers[i]);

    size_t cacheReversedH = 0;
    for (int i = l - 1; i > -1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(syncmers[i]);

    int iorder = 0;
    // Skip if strand ambiguous
    if (cacheForwardH < cacheReversedH) {
        kminmers.emplace_back(cacheForwardH,  std::get<1>(syncmers[0]), std::get<2>(syncmers[l-1]), false, iorder);
        hashes.insert(cacheForwardH);
        ++iorder;
    } else if (cacheReversedH < cacheForwardH) {
        kminmers.emplace_back(cacheReversedH, std::get<1>(syncmers[0]), std::get<2>(syncmers[l-1]), true, iorder);
        hashes.insert(cacheReversedH);
        ++iorder;
    }
    
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    for (int i = 1; i < syncmers.size() - l + 1; ++i) {
        cacheForwardH = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);

        cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(syncmers[i+l-1]) << (2 * k * (l - 1)));

        // Skip if strand ambiguous
        if (cacheForwardH < cacheReversedH) {
            kminmers.emplace_back(cacheForwardH,  std::get<1>(syncmers[i]), std::get<2>(syncmers[i+l-1]), false, iorder);
            hashes.insert(cacheForwardH);
            ++iorder;
        } else if (cacheReversedH < cacheForwardH) {
            kminmers.emplace_back(cacheReversedH, std::get<1>(syncmers[i]), std::get<2>(syncmers[i+l-1]), true, iorder);
            hashes.insert(cacheReversedH);
            ++iorder;
        }
    }

    return std::make_pair(std::move(kminmers), std::move(hashes));
}

int32_t degapGlobal(const int32_t& globalCoord, const std::map<int32_t, int32_t>& coordsIndex) {
    auto coordIt = coordsIndex.upper_bound(globalCoord);
    if (coordIt == coordsIndex.begin()) {
        return 0;
    }
    --coordIt;
    return globalCoord - coordIt->second;
}

void mutateSeedmers(mgsr::seedmers& seedmers, const std::pair<int32_t, std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>>>& index, const std::map<int32_t, int32_t>& coordsIndex, std::unordered_set<size_t>& affectedSeedmers, std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>>& seedmersToRevert) {
    // const auto& length = index.first;
    for (const auto& change : index.second) {
        auto [beg, end, hash, rev, type] = change;
        if (type == 1) {
            auto positionMapIt = seedmers.positionMap.find(beg);
            auto seedmersMapIt = seedmers.seedmersMap.find(hash);
            auto [cend, chash, crev] = positionMapIt->second;
            seedmersToRevert.emplace_back(std::make_tuple(beg, cend, chash, crev, 0));

            assert(seedmers.positionMap.erase(beg));
            assert(seedmersMapIt->second.erase(beg));
            size_t curHashNum = seedmersMapIt->second.size();
            if (curHashNum <= 1) {
                affectedSeedmers.insert(hash);
                if (curHashNum == 0) {
                    seedmers.seedmersMap.erase(hash);
                }
            }
        } else if (type == 0) {
            auto positionMapIt = seedmers.positionMap.find(beg);
            if (positionMapIt != seedmers.positionMap.end()) {
                auto [oend, ohash, orev] = positionMapIt->second;
                auto seedmersMapIt = seedmers.seedmersMap.find(ohash);

                assert(seedmersMapIt->second.erase(beg));
                size_t ohashNum = seedmersMapIt->second.size();
                if (ohashNum <= 1) {
                    affectedSeedmers.insert(ohash);
                    if (ohashNum == 0) {
                        seedmers.seedmersMap.erase(ohash);
                    }
                }
                seedmersToRevert.emplace_back(std::make_tuple(beg, oend, ohash, orev, 0));
            } else {
                seedmersToRevert.emplace_back(std::make_tuple(beg, 0, hash, false, 1));
            }
            seedmers.positionMap[beg] = std::make_tuple(end, hash, rev);
            seedmers.seedmersMap[hash].insert(beg);
            affectedSeedmers.insert(hash);
        } else if (type == 2) {
            auto positionMapIt = seedmers.positionMap.find(beg);
            auto [cend, chash, crev] = positionMapIt->second;
            seedmers.positionMap[beg] = std::make_tuple(end, chash, crev);
            seedmersToRevert.emplace_back(std::make_tuple(beg, cend, 0, false, 2));
        } else {
            throw std::invalid_argument("Error reading index file. Can't determine seedmer change type.");
        }
    }
}

void revertSeedmers(mgsr::seedmers& seedmers, const std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>>& seedmersToRevert) {
    for (int i = seedmersToRevert.size() - 1; i > -1; --i) {
        auto [beg, end, hash, rev, type] = seedmersToRevert[i];
        if (type == 1) {
            auto seedmersMapIt = seedmers.seedmersMap.find(hash);
            assert(seedmers.positionMap.erase(beg));
            assert(seedmersMapIt->second.erase(beg));
            if (seedmersMapIt->second.empty()) {
                seedmers.seedmersMap.erase(hash);
            }
        } else if (type == 0) {
            auto positionMapIt = seedmers.positionMap.find(beg);
            if (positionMapIt != seedmers.positionMap.end()) {
                auto [oend, ohash, orev] = positionMapIt->second;
                auto seedmersMapIt = seedmers.seedmersMap.find(ohash);

                assert(seedmersMapIt->second.erase(beg));
                if (seedmersMapIt->second.empty()) {
                    seedmers.seedmersMap.erase(ohash);
                }
            }
            seedmers.positionMap[beg] = std::make_tuple(end, hash, rev);
            seedmers.seedmersMap[hash].insert(beg);
        } else if (type == 2) {
            auto positionMapIt = seedmers.positionMap.find(beg);
            auto [cend, chash, crev] = positionMapIt->second;
            seedmers.positionMap[beg] = std::make_tuple(end, chash, crev);
        }
    }
}

struct seedmerHasher {
    std::size_t operator()(const std::tuple<size_t, int32_t, int32_t, bool, int32_t>& t) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, std::get<0>(t));
        // boost::hash_combine(seed, std::get<1>(t));
        // boost::hash_combine(seed, std::get<2>(t));
        boost::hash_combine(seed, std::get<3>(t));
        boost::hash_combine(seed, std::get<4>(t));
        return seed;
    }
};

struct seedmerSetHasher {
    std::size_t operator()(const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& vec) const {
        std::size_t seed = 0;
        for (const auto& t : vec) {
            boost::hash_combine(seed, seedmerHasher{}(t));
        }
        return seed;
    }
};

struct seedmerSetEqual {
    bool operator()(const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& vec1,
                    const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& vec2) const {
        if (vec1.size() != vec2.size()) return false;
        for (size_t i = 0; i < vec1.size(); ++i) {
            if (std::get<0>(vec1[i]) != std::get<0>(vec2[i]) ||
                std::get<3>(vec1[i]) != std::get<3>(vec2[i]) ||
                std::get<4>(vec1[i]) != std::get<4>(vec2[i])) {
                return false;
            }
        }
        return true;
        // return vec1 == vec2;
    }
};

static std::string reverseComplement(std::string dna_sequence) {
    std::string complement = "";
    for (char c : dna_sequence) {
        switch (c) {
            case 'A': complement += 'T'; break;
            case 'T': complement += 'A'; break;
            case 'C': complement += 'G'; break;
            case 'G': complement += 'C'; break;
            default: complement += c; break;
        }
    }
    std::reverse(complement.begin(), complement.end());
    return complement;
}

std::string toUpper(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) {
        return std::toupper(c);
    });
    return result;
}

void initializeFastq(
    const std::string &fastqPath, std::vector<std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>>>& readSeedmers,
    std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::vector<std::vector<seeding::seed>>& readSeeds, std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames, const int32_t k, const int32_t s, const int32_t l
    ) {
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    if(!fp){
        std::cerr << "Error: File " << fastqPath << " not found" << std::endl;
        exit(0);
    }

    seq = kseq_init(fileno(fp));
    std::unordered_map<std::string, std::pair<int32_t, std::vector<size_t>>> dupMarkedReads;
    std::vector<std::vector<seed>> readSeedsFwd;
    std::vector<std::vector<seed>> readSeedsBwd;
    int line;
    size_t curIndex = 0;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(toUpper(seq->seq.s));
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.s);

        std::vector<seed> syncmers = seeding::syncmerize(readSequences.back(), k, s, false, false, 0);
        readSeedsFwd.push_back(syncmers);
        std::vector<seed> syncmersReverse = seeding::syncmerize(reverseComplement(readSequences.back()), k, s, false, false, 0);
        readSeedsBwd.push_back(syncmersReverse);

        ++dupMarkedReads[readSequences.back()].first;
        dupMarkedReads[readSequences.back()].second.push_back(curIndex);
        ++curIndex;
    } 

    readSeeds.reserve( readSeedsFwd.size());
    for(int i = 0; i < readSeedsFwd.size(); i++) {
        std::vector<seed> thisReadsSeeds;
        thisReadsSeeds.reserve(readSeedsFwd[i].size() + readSeedsBwd[i].size());

        for(int j = 0; j < readSeedsFwd[i].size(); j++) {
            readSeedsFwd[i][j].reversed = false;
            readSeedsFwd[i][j].pos = readSeedsFwd[i][j].pos + k - 1; // Minimap standard
            thisReadsSeeds.push_back(readSeedsFwd[i][j]);
        }
        for(int j = 0; j < readSeedsBwd[i].size(); j++) {
            readSeedsBwd[i][j].reversed = true;
            readSeedsBwd[i][j].pos = readSeedsBwd[i][j].pos + k - 1; // Minimap standard
            thisReadsSeeds.push_back(readSeedsBwd[i][j]);
        }
        readSeeds.push_back(thisReadsSeeds);
    }
    
    std::unordered_map<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, size_t, seedmerSetHasher, seedmerSetEqual> duplicateCounts;
    for (const auto& dupMarkedRead : dupMarkedReads) {
        const std::string& seq = dupMarkedRead.first;
        const int32_t&  numDup = dupMarkedRead.second.first;
        std::vector<std::tuple<size_t, int32_t, int32_t>> curSyncmers = syncmersSketch(seq, k, s, false);

        readSeedmers_t curKminmers;
        if (curSyncmers.size() < l) {
            curKminmers = {};
        } else {
            curKminmers = extractKminmers(curSyncmers, k, l);
        }
        auto it = duplicateCounts.find(curKminmers.first);
        if (it == duplicateCounts.end()) {
            duplicateCounts[curKminmers.first] = readSeedmers.size();
            readSeedmers.push_back(std::make_pair(std::move(curKminmers.first), std::move(curKminmers.second)));
            numReadDuplicates.emplace_back(std::make_pair(numDup, std::move(dupMarkedRead.second.second)));
        } else {
            std::vector<size_t> merged(dupMarkedRead.second.second.size() + numReadDuplicates[it->second].second.size());
            std::merge(numReadDuplicates[it->second].second.begin(), numReadDuplicates[it->second].second.end(), dupMarkedRead.second.second.begin(), dupMarkedRead.second.second.end(), merged.begin());
            numReadDuplicates[it->second].first += numDup;
            numReadDuplicates[it->second].second = std::move(merged);
        }
    }
}

bool redo(const std::unordered_set<size_t>& a, const std::unordered_set<size_t>& b) {
    const auto& smallerSet = (a.size() < b.size()) ? a : b;
    const auto& largerSet  = (a.size() < b.size()) ? b : a;

    for (const auto& h : smallerSet) {
        if (largerSet.find(h) != largerSet.end()) {
            return true;
        }
    }

    return false;
}

int32_t extend(match_t& curMatch, const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& querySeedmers, const mgsr::seedmers& refSeedmers, int32_t qidx, std::map<int32_t, std::tuple<int32_t, size_t, bool>>::const_iterator refPositionIt, int32_t c) {
    if (qidx == querySeedmers.size() - 1) return c;
    const auto& [qbeg, qend, rbeg, rend, rev, _] = curMatch;
    const auto& [nhash, nqbeg, nqend, nqrev, nqidx] = querySeedmers[qidx+1];
    auto prevRefPositionIt = refPositionIt;
    auto nextRefPositionIt = refPositionIt;
    --prevRefPositionIt;
    ++nextRefPositionIt;
    if (refSeedmers.seedmersMap.count(nhash) > 0 && refSeedmers.seedmersMap.find(nhash)->second.size() < 2) {
        assert(refSeedmers.seedmersMap.find(nhash)->second.size() == 1);
        const auto& rbeg = *(refSeedmers.seedmersMap.find(nhash)->second.begin());
        auto curRefPositionIt = refSeedmers.positionMap.find(rbeg);
        const auto& [rend, rhash, rrev] = curRefPositionIt->second;
        if (rev == (nqrev != rrev)) {
            if ((rev == 0 && curRefPositionIt->first == nextRefPositionIt->first) || (rev == 1 && curRefPositionIt->first == prevRefPositionIt->first)) {
                ++c;
                if (rev == 0) {
                    curMatch = std::make_tuple(std::get<0>(curMatch), nqend, std::get<2>(curMatch), rend, std::get<4>(curMatch), c);
                } else if (rev == 1) {
                    curMatch = std::make_tuple(std::get<0>(curMatch), nqend, rbeg, std::get<3>(curMatch), std::get<4>(curMatch), c);
                }
                return extend(curMatch, querySeedmers, refSeedmers, nqidx, curRefPositionIt, c);
            }
        }
    }
    return c;
}

// typedef std::tuple<int, int, int, int, bool, int> match_t;
// query start, query end, ref start, ref end, strand, count
std::vector<match_t> match(const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& querySeedmers, const mgsr::seedmers& refSeedmers, double& duplicates) {
    std::vector<match_t> matches;
    int32_t i = 0;
    while (i < querySeedmers.size()) {
        const auto& [hash, qbeg, qend, qrev, qidx] = querySeedmers[i];
        int32_t c = 1;
        if (refSeedmers.seedmersMap.count(hash) > 0) {
            if (refSeedmers.seedmersMap.find(hash)->second.size() < 2) {
                assert(refSeedmers.seedmersMap.find(hash)->second.size() == 1);
                const auto& rbeg = *(refSeedmers.seedmersMap.find(hash)->second.begin());
                auto curRefPositionIt = refSeedmers.positionMap.find(rbeg);
                const auto& [rend, rhash, rrev] = curRefPositionIt->second;
                match_t curMatch = std::make_tuple(qbeg, qend, rbeg, rend, qrev != rrev, c);
                c = extend(curMatch, querySeedmers, refSeedmers, qidx, curRefPositionIt, c);
                matches.push_back(curMatch);
            } else {
                duplicates += 1.0;
            }
        }
        i += c; 
    }
    return matches;
}

bool isColinear(const match_t& match1, const match_t& match2, const std::map<int32_t, int32_t>& coordsIndex, int maximumGap) {
    const auto& [qbeg1, qend1, rglobalbeg1, rglobalend1, rev1, count1] = match1;
    const auto& [qbeg2, qend2, rglobalbeg2, rglobalend2, rev2, count2] = match2;
    if (rev1 != rev2) return false;
    
    auto rbeg1 = degapGlobal(rglobalbeg1, coordsIndex);
    auto rend1 = degapGlobal(rglobalend1, coordsIndex);
    auto rbeg2 = degapGlobal(rglobalbeg2, coordsIndex);
    auto rend2 = degapGlobal(rglobalend2, coordsIndex);

    if (rev1 == false) {
        int32_t qgap = abs(qbeg2 - qend1);
        int32_t rgap = abs(rbeg2 - rend1);
        if (rbeg1 < rbeg2 && abs(qgap - rgap) < maximumGap) return true;
    } else {
        int32_t qgap = abs(qbeg2 - qend1);
        int32_t rgap = abs(rbeg1 - rend2);
        if (rbeg2 < rbeg1 && abs(qgap - rgap) < maximumGap) return true;
    }
    return false;

}

std::vector<match_t> chainPseudo(const std::vector<match_t>& matches, const std::map<int32_t, int32_t>& coordsIndex, int maximumGap, int minimumCount, int minimumScore) {
    std::vector<match_t> pseudoChain;

    if (matches.size() == 0) {
        return pseudoChain;
    }
    else if (matches.size() == 1) {
        pseudoChain.push_back(matches[0]);
        return pseudoChain;
    }

    size_t maxIndex = 0;
    for (size_t i = 1; i < matches.size(); ++i) {
        if (std::get<5>(matches[i]) > std::get<5>(matches[maxIndex])) maxIndex = i;
    }

    for (size_t i = 0; i < matches.size(); ++i) {
        if (i == maxIndex) {
            pseudoChain.push_back(matches[i]);
            continue;
        }

        // typedef std::tuple<int, int, int, int, bool, int> match_t;
        // query start, query end, ref start, ref end, strand, count
        if (std::get<0>(matches[maxIndex]) < std::get<0>(matches[i])) {
            if (isColinear(matches[maxIndex], matches[i], coordsIndex, maximumGap)) pseudoChain.push_back(matches[i]);
        } else {
            if (isColinear(matches[i], matches[maxIndex], coordsIndex, maximumGap)) pseudoChain.push_back(matches[i]);
        }
    }

    return pseudoChain;
}

int scorePseudoChain(const std::vector<match_t>& pseudoChain) {
    int32_t score = 0;
    for (const match_t& match : pseudoChain) {
        score += std::get<5>(match);
    }
    return score;
}

void printReadSeedmers(const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& seedmers) {
    std::cout << "seedmers:" << std::endl;
    for (const auto& seedmer : seedmers) {
        std::cout << std::get<0>(seedmer) << "\t"
                  << std::get<1>(seedmer) << "\t"
                  << std::get<2>(seedmer) << "\t"
                  << std::get<3>(seedmer) << "\t"
                  << std::get<4>(seedmer) << std::endl;
    }
    std::cout << std::endl;
}

void printMatches(const std::vector<match_t> matches) {
    std::cout << "Matches:" << std::endl;
    for (const match_t& match : matches) {
        std::cout << std::get<0>(match) << "\t"
                  << std::get<1>(match) << "\t"
                  << std::get<2>(match) << "\t"
                  << std::get<3>(match) << "\t"
                  << std::get<4>(match) << "\t"
                  << std::get<5>(match) << std::endl;
    }
    std::cout << std::endl;
}

void scoreDFS(
    mgsr::seedmers& seedmers, const std::unordered_map<std::string, std::pair<int32_t, std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>>>>& seedmersIndex,
    const std::unordered_map<std::string, std::map<int32_t, int32_t>>& coordsIndex,
    const std::vector<std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>>>& readSeedmers,
    const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    std::unordered_map<std::string, std::string>& identicalPairs, const Node *node, Tree *T,
    const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate
    ) {
    // std::cerr << "identifier " << node->identifier << std::endl;
    std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>> seedmersToRevert;
    std::unordered_set<size_t> affectedSeedmers;
    mutateSeedmers(seedmers, seedmersIndex.find(node->identifier)->second, coordsIndex.find(node->identifier)->second, affectedSeedmers, seedmersToRevert);

    if (node->identifier == T->root->identifier) {
        allScores[node->identifier].resize(readSeedmers.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, readSeedmers.size(), 1000),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const auto& curReadSeedmers = readSeedmers[i];
                    double duplicates = 0;
                    std::vector<match_t> matches = match(curReadSeedmers.first, seedmers, duplicates);
                    std::vector<match_t> pseudoChain = chainPseudo(matches, coordsIndex.find(node->identifier)->second, maximumGap, minimumCount, minimumScore);
                    int32_t pseudoScore = scorePseudoChain(pseudoChain);
                    double maxScore = static_cast<double>(curReadSeedmers.first.size()) - duplicates;
                    // Need to change how pseudoProb is calculated
                    // pseudoProb = 0 if pseudoScore < 1/2 (estimated number of kminmers on a read)
                    double pseudoProb;
                    if (pseudoScore || maxScore == 0) {
                        pseudoProb = std::numeric_limits<double>::min();
                    } else {
                        pseudoProb = pow(errorRate, maxScore - static_cast<double>(pseudoScore)) * pow(1 - errorRate, static_cast<double>(pseudoScore));
                    }
                    // double  pseudoProb  = static_cast<double>(pseudoScore) / (static_cast<double>(curReadSeedmers.first.size() - duplicates));
                    assert(pseudoProb <= 1.0);
                    allScores[node->identifier][i] = {pseudoScore, pseudoProb};
                }
        });
    } else {
        allScores[node->identifier].resize(readSeedmers.size());
        if (affectedSeedmers.size() == 0) {
            identicalPairs[node->identifier] = node->parent->identifier;
            allScores[node->identifier] = allScores[node->parent->identifier];
        } else {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, readSeedmers.size(), 1000),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto& curReadSeedmers = readSeedmers[i];
                        if (redo(curReadSeedmers.second, affectedSeedmers)) {
                            double duplicates = 0;
                            std::vector<match_t> matches = match(curReadSeedmers.first, seedmers, duplicates);
                            std::vector<match_t> pseudoChain = chainPseudo(matches, coordsIndex.find(node->identifier)->second, maximumGap, minimumCount, minimumScore);
                            int32_t pseudoScore = scorePseudoChain(pseudoChain);
                            double maxScore = static_cast<double>(curReadSeedmers.first.size()) - duplicates;
                            double pseudoProb = pow(errorRate, maxScore - static_cast<double>(pseudoScore)) * pow(1 - errorRate, static_cast<double>(pseudoScore));                            // double  pseudoProb  = static_cast<double>(pseudoScore) / (static_cast<double>(curReadSeedmers.first.size() - duplicates));
                            assert(pseudoProb <= 1.0);
                            allScores[node->identifier][i] = {pseudoScore, pseudoProb};
                        } else {
                            allScores[node->identifier][i] = allScores[node->parent->identifier][i];
                        }
                    }
            });
        }
    }


    for (Node *child : node->children) {
        scoreDFS(seedmers, seedmersIndex, coordsIndex, readSeedmers, numReadDuplicates, allScores, identicalPairs, child, T, maximumGap, minimumCount, minimumScore, errorRate);
    }
    revertSeedmers(seedmers, seedmersToRevert);
}

bool identicalReadScores(const tbb::concurrent_vector<std::pair<int32_t, double>>& scores1, const tbb::concurrent_vector<std::pair<int32_t, double>>& scores2) {
    assert(scores1.size() == scores2.size());
    for (size_t i = 0; i < scores1.size(); ++i) {
        if (scores1[i].second != scores2[i].second) return false;
    }
    return true;
}

void updateIdenticalSeedmerSets(
    const std::unordered_set<std::string>& identicalGroup,
    const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor,
    std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets
    ) {
    std::unordered_set<std::string> seenNodes;
    std::unordered_set<std::string> unseenNodes = identicalGroup;
    for (const std::string& currNode : identicalGroup) {
        if (seenNodes.find(currNode) != seenNodes.end()) continue;
        seenNodes.insert(currNode);
        unseenNodes.erase(currNode);
        std::unordered_set<std::string> identicals;
        for (const std::string& idenNode : unseenNodes) {
            if (identicalReadScores(allScores.at(currNode), allScores.at(idenNode))) {
                leastRecentIdenticalAncestor[idenNode] = currNode;
                identicalSets[currNode].insert(idenNode);
                if (identicalSets.find(idenNode) != identicalSets.end()) {
                    for (const auto& idenOffspring : identicalSets[idenNode]) {
                        leastRecentIdenticalAncestor[idenOffspring] = currNode;
                        identicalSets[currNode].insert(idenOffspring);
                    }
                    identicalSets.erase(idenNode);
                }
                identicals.insert(idenNode);
            }
        }
        for (const auto& identical : identicals) {
            seenNodes.insert(identical);
            unseenNodes.erase(identical);
        }
    }
}

void writeScores(
    const std::string& nodeId, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, const std::vector<std::string>& readSequences,
    const std::vector<std::string>& readNames, const std::string& outPath
    ) {
    std::ofstream ofs(outPath);
    const auto& scores = allScores.at(nodeId);
    assert(scores.size() == numReadDuplicates.size());
    for (size_t i = 0; i < scores.size(); ++i) {
        const auto& duplicates = numReadDuplicates[i];
        const auto& score = scores[i];
        assert(duplicates.first == duplicates.second.size());
        for (const size_t& idx : duplicates.second) {
            ofs << readNames[idx]     << "\t"
                << score.first        << "\t"
                << score.second       << "\t"
                << readSequences[idx] << "\n";
        }
    }
    ofs.close();
}

void mgsr::scorePseudo(
    std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path,
    std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, 
    std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor,
    std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets, int32_t& numReads, Tree *T,
    std::vector<std::vector<seeding::seed>>& readSeeds, std::vector<std::string>& readSequences, std::vector<std::string>& readQuals,
    std::vector<std::string>& readNames, const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate
    ) {
    // get read seeds
    std::string line;
    std::getline(indexFile, line);
    std::vector<std::string> spltTop;
    PangenomeMAT::stringSplit(line, ' ', spltTop);
    int32_t tempK = std::stoi(spltTop[0]);
    int32_t tempS = std::stoi(spltTop[1]);
    int32_t tempJ = std::stoi(spltTop[2]);

    std::vector<std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>>> readSeedmers;
    //                 node         parent
    std::unordered_map<std::string, std::string> identicalPairs;
    //                 node         children, grandchildren, etc.                  

    //                 node                                beg      end      hash    rev   type
    std::unordered_map<std::string, std::pair<int32_t, std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>>>> seedmersIndex;
    std::unordered_map<std::string, std::map<int32_t, int32_t>> coordsIndex;
    
    std::cerr << "start reading tree seedmers index" << std::endl;
    while (std::getline(indexFile, line)) {
        std::vector<std::string> split;
        PangenomeMAT::stringSplit(line, ' ', split);
        std::string nodeInfo = split[0];
        std::vector<std::string> nodeInfoSplit;
        PangenomeMAT::stringSplit(nodeInfo, ':', nodeInfoSplit);
        std::string nid = nodeInfoSplit[0];
        int32_t length = std::stoi(nodeInfoSplit[1]);
        seedmersIndex[nid] = {length, {}};
        coordsIndex[nid] = {};
        for (int32_t i = 1; i < split.size(); ++i) {
            std::vector<std::string> metaSplit;
            PangenomeMAT::stringSplit(split[i], ':', metaSplit);
            if (metaSplit[0] == "c") {
                std::vector<std::string> coorsSplit;
                PangenomeMAT::stringSplit(metaSplit[1], ',', coorsSplit);
                assert(coorsSplit.size() % 2 == 0);
                for (size_t i = 0; i < coorsSplit.size(); i+=2) {
                    int32_t globalCoord = std::stoi(coorsSplit[i]);
                    int32_t offset = std::stoi(coorsSplit[i+1]);
                    coordsIndex[nid][globalCoord] = offset;
                }
            } else {
                int32_t beg = std::stoi(metaSplit[0]);
                if (metaSplit[1] == "+") {
                    std::vector<std::string> infoSplit;
                    PangenomeMAT::stringSplit(metaSplit[2], ',', infoSplit);
                    int32_t end = std::stoi(infoSplit[0]);
                    std::stringstream sstream(infoSplit[1]);
                    size_t hash;
                    sstream >> hash;
                    bool rev = (infoSplit[2] == "1") ? true : false;
                    seedmersIndex[nid].second.emplace_back(std::make_tuple(beg, end, hash, rev, 0)); 
                }
                else if (metaSplit[1] == "-") {
                    std::stringstream sstream(metaSplit[2]);
                    size_t hash;
                    sstream >> hash;
                    seedmersIndex[nid].second.emplace_back(std::make_tuple(beg, 0, hash, false, 1));
                } else if (metaSplit[1] == "e") {
                    int32_t newEnd = std::stoi(metaSplit[2]);
                    seedmersIndex[nid].second.emplace_back(std::make_tuple(beg, newEnd, 0, false, 2));
                } else {
                    throw std::invalid_argument("Error reading index file. Can't determine insertion/subsitution or deletion");
                }
            }
        }
    }
    std::cerr << "finished reading tree seedmers index\n" << std::endl;

    std::cerr << "start initializing read seedmers" << std::endl; 
    initializeFastq(reads1Path, readSeedmers, numReadDuplicates, readSeeds, readSequences, readQuals, readNames, tempK, tempS, tempJ);
    assert(readSeedmers.size() == numReadDuplicates.size());
    numReads = readSequences.size();
    std::cerr << "finished initializing read seedmers... total number of reads " << readSequences.size() << "\n" << std::endl;

    std::cerr << "start scoring DFS" << std::endl;
    mgsr::seedmers seedmers;
    size_t totalCount = 0;
    size_t redoCount = 0;
    
    scoreDFS(seedmers, seedmersIndex, coordsIndex, readSeedmers, numReadDuplicates, allScores, identicalPairs, T->root, T, maximumGap, minimumCount, minimumScore, errorRate);
    std::cerr << "finished scoring DFS\n" << std::endl;

    // writeScores("OP573714.1", allScores, numReadDuplicates, readSequences, readNames, "OP573714.1_scores.txt");
    // writeScores("MW671606.1", allScores, numReadDuplicates, readSequences, readNames, "MW671606.1_scores.txt");
    // writeScores("OM410575.1", allScores, numReadDuplicates, readSequences, readNames, "OM410575.1_scores.txt");


    // //                 node         parent
    // std::unordered_map<std::string, std::string> identicalPairs;
    // //                 node         children, grandchildren, etc.                  
    // std::unordered_map<std::string, std::unordered_set<std::string>> identicalSets
    for (const auto& pair : identicalPairs) {
        std::unordered_set<std::string> curIdenticals;
        std::string curNode = pair.first;
        std::string curParent = pair.second;
        curIdenticals.insert(curNode);
        while (identicalPairs.find(curParent) != identicalPairs.end()) {
            curNode = curParent;   
            curParent = identicalPairs[curParent];
            curIdenticals.insert(curNode);      
        }
        for (const auto& node : curIdenticals) {
            identicalSets[curParent].insert(node);
        }
    }

    for (const auto& set : identicalSets) {
        for (const auto& offspring : set.second) {
            leastRecentIdenticalAncestor[offspring] = set.first;
        }
    }

    std::cerr << "First round of duplication removal: " << leastRecentIdenticalAncestor.size() << std::endl;

    std::vector<std::pair<std::string, int32_t>> scores(allScores.size() - leastRecentIdenticalAncestor.size());
    for (const auto& node : allScores) {
        if (leastRecentIdenticalAncestor.find(node.first) != leastRecentIdenticalAncestor.end()) continue;
        int32_t score = 0;
        for (size_t i = 0; i < node.second.size(); ++i) {
            score += node.second[i].first * numReadDuplicates[i].first;
        }
        scores.emplace_back(std::make_pair(node.first, score));
    }

    std::sort(scores.begin(), scores.end(), [](const auto &a, const auto &b) {
        return a.second > b.second;
    });

    std::unordered_set<std::string> identicalGroup;
    for (size_t i = 0; i < scores.size() - 1; ++i) {
        const auto& currScore = scores[i];
        const auto& nextScore = scores[i+1];
        if (currScore.second == nextScore.second) {
            identicalGroup.insert(currScore.first);
            identicalGroup.insert(nextScore.first);
        } else {
            if (!identicalGroup.empty()) {
                updateIdenticalSeedmerSets(identicalGroup, allScores, leastRecentIdenticalAncestor, identicalSets);
                std::unordered_set<std::string>().swap(identicalGroup);
            }
        }
    }
    if (!identicalGroup.empty()) {
        updateIdenticalSeedmerSets(identicalGroup, allScores, leastRecentIdenticalAncestor, identicalSets);
    }
    std::cerr << "Second round of duplication removal: " << leastRecentIdenticalAncestor.size() << "\n" << std::endl;
}

double getExp(const Eigen::MatrixXd& probs, const Eigen::VectorXd& props, const Eigen::VectorXd& numReadDuplicates) {
    assert(props.size() == probs.cols());

    Eigen::VectorXd readSums = probs * props;
    double llh = (numReadDuplicates.array() * readSums.array().log()).sum();

    return llh;
}

Eigen::VectorXd getMax(const Eigen::MatrixXd& probs, const Eigen::VectorXd& props, const Eigen::VectorXd& numReadDuplicates, const int32_t& totalReads) {
    size_t numNodes = probs.cols();

    Eigen::VectorXd denoms = probs * props;

    Eigen::VectorXd newProps(numNodes);
    newProps.setZero();

    for (size_t i = 0; i < numNodes; ++i) {
        Eigen::VectorXd ratios = (probs.col(i).array() * props[i]) / denoms.array();
        double newProp = (numReadDuplicates.array() * ratios.array()).sum();
        newProp /= totalReads;
        newProps(i) = newProp;
    }

    return newProps;

}

void normalize(Eigen::VectorXd& props) {    
    for (int i = 0; i < props.size(); ++i) {
        if (props(i) <= 0) {
            props(i) = std::numeric_limits<double>::min();
        }
    }
    double sum = props.sum();
    props /= sum;
}

void squarem(
    const std::vector<std::string>& nodes, const Eigen::MatrixXd& probs,
    const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    const Eigen::VectorXd& numReadDuplicates, const int32_t& numReads, 
    Eigen::VectorXd& props, double& llh, int& curit
    ) {
    assert(nodes.size() == probs.cols());
    assert(nodes.size() == props.size());
    while (true) {
        // std::cerr << "it " << curit << std::endl;
        Eigen::VectorXd theta1 = getMax(probs, props, numReadDuplicates, numReads);
        normalize(theta1);
        Eigen::VectorXd theta2 = getMax(probs, theta1, numReadDuplicates, numReads);
        normalize(theta2);

        Eigen::VectorXd r = theta1 - props;
        Eigen::VectorXd v = theta2 - theta1 - r;
        double r_norm = r.norm();
        double v_norm = v.norm();

        double alpha;
        if (r_norm == 0 || v_norm == 0) {
            alpha = 0;
        } else {
            alpha = - r_norm / v_norm;
        }
        double newllh;

        
        Eigen::VectorXd theta_p;
        if (alpha > -1) {
            alpha = -1;
            theta_p = props - 2 * alpha * r + alpha * alpha * v;
            props = getMax(probs, theta_p, numReadDuplicates, numReads);
            normalize(props);
            newllh = getExp(probs, props, numReadDuplicates);
        } else {
            theta_p = props - 2 * alpha * r + alpha * alpha * v;
            auto newProps = getMax(probs, theta_p, numReadDuplicates, numReads);
            normalize(newProps);
            newllh = getExp(probs, newProps, numReadDuplicates);
            if (newllh >= llh) {
                props = std::move(newProps);
            } else {
                while (llh - newllh > 0.00001) {
                    // std::cerr << "alpha " << alpha << std::endl;
                    alpha = (alpha - 1) / 2;
                    theta_p = props - 2 * alpha * r + alpha * alpha * v;
                    newProps = getMax(probs, theta_p, numReadDuplicates, numReads);
                    normalize(newProps);
                    newllh   = getExp(probs, newProps, numReadDuplicates);
                }
                props = std::move(newProps);
            }
        }

        if (newllh - llh < 0.00001) {
            llh = newllh;
            break;
        }        
        llh = newllh;
        ++curit;
    }
    ++curit;
}

void mgsr::squaremHelper(
    PangenomeMAT::Tree *T, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates,
    const int32_t& numReads, const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    Eigen::MatrixXd& probs, std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const int32_t& roundsRemove, const double& removeThreshold, std::string exclude
    ) {
    if (exclude.empty()) {
        std::stringstream msg;
        msg << "starting EM estimation of haplotype proportions" << "\n";
        std::cerr << msg.str();
    } else {
        std::stringstream msg;
        msg << "starting EM estimation of haplotype proportions excluding " << exclude << "\n";
        std::cerr << msg.str();
    }


    if (!exclude.empty()) {
        probs.resize(allScores.begin()->second.size(), allScores.size() - leastRecentIdenticalAncestors.size() - 1);
    } else {
        probs.resize(allScores.begin()->second.size(), allScores.size() - leastRecentIdenticalAncestors.size());
    }
    size_t colIndex = 0;
    for (const auto& node : allScores) {
        if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
        if (!exclude.empty() && node.first == exclude) continue;
        std::vector<double> curProbs;
        size_t rowIndex = 0;
        for (const auto& score : node.second) {
            probs(rowIndex, colIndex) = score.second;
            ++rowIndex;
        }
        nodes.push_back(node.first);
        ++colIndex;
    }

    // std::cerr << "num nodes " << nodes.size() << std::endl;
    props = Eigen::VectorXd::Constant(nodes.size(), 1.0 / static_cast<double>(nodes.size()));
    Eigen::VectorXd readDuplicates(allScores.begin()->second.size());
    for (size_t i = 0; i < numReadDuplicates.size(); ++i) {
        readDuplicates(i) = numReadDuplicates[i].first;
    }

    
    int curit = 0;
    llh = getExp(probs, props, readDuplicates);
    // std::cerr << "iteration " << curit << ": " << llh << std::endl;
    squarem(nodes, probs, identicalSets, readDuplicates, numReads, props, llh, curit);

    for (int32_t i = 0; i < roundsRemove; ++i) {
        std::vector<size_t> significantIndices;
        std::vector<std::string> sigNodes;

        for (size_t i = 0; i < props.size(); ++i) {
            if (props(i) >= removeThreshold) {
                significantIndices.push_back(i);
            }
        }
        if (significantIndices.size() == nodes.size()) break;
        std::cerr << "remove round " << i + 1 << std::endl;

        for (size_t idx : significantIndices) {
            sigNodes.push_back(nodes[idx]);
        }

        Eigen::MatrixXd sigProbs(probs.rows(), significantIndices.size());
        sigProbs.resize(probs.rows(), significantIndices.size());
        for (size_t i = 0; i < significantIndices.size(); ++i) {
            sigProbs.col(i) = probs.col(significantIndices[i]);
        }
        Eigen::VectorXd sigProps = Eigen::VectorXd::Constant(sigNodes.size(), 1.0 / static_cast<double>(sigNodes.size()));
        llh = getExp(sigProbs, sigProps, readDuplicates);
        squarem(sigNodes, sigProbs, identicalSets, readDuplicates, numReads, sigProps, llh, curit);
        nodes = std::move(sigNodes);
        probs = std::move(sigProbs);
        props = std::move(sigProps);
    }

    if (!exclude.empty()) {
        Eigen::VectorXd curProbs(numReads);
        for (size_t i = 0; i < numReads; ++i) {
            curProbs(i) = allScores.at(exclude)[i].second;
        }

        probs = probs, curProbs;
        props = props, 0.0;
        llh = getExp(probs, props, readDuplicates);
    }

    if (exclude.empty()) {
        std::stringstream msg;
        msg << "Finished EM estimation of haplotype proportions. Total EM iterations: " << curit << "\n";
        std::cerr << msg.str();
    } else {
        std::stringstream msg;
        msg << "Finished EM estimation of haplotype proportions excluding " << exclude << ". Total EM iterations: " << curit<< "\n";
        std::cerr << msg.str();
    }
}


void mgsr::accio(
    const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    const std::vector<std::string>& nodes, const Eigen::MatrixXd& probs, const Eigen::VectorXd& props,
    const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates,
    std::unordered_map<std::string, std::unordered_set<size_t>>& assignedReads
    ) {
        for (size_t i = 0; i < probs.rows(); ++i) {
            const Eigen::VectorXd& curprobs = probs.row(i);
            double curmax = curprobs.maxCoeff();
            // size_t numMax = 0;
            // size_t maxIdx = std::numeric_limits<size_t>::max();
            for (size_t j = 0; j < curprobs.size(); ++j) {
                if (curmax == curprobs(j)) {
                    // ++numMax;
                    // maxIdx = j;
                    for (size_t readIndex : numReadDuplicates[i].second) {
                        assignedReads[nodes[j]].insert(readIndex);
                    }
                }
            }
            // assert(numMax >= 1);
            // if (numMax == 1) {
            //     for (size_t readIndex : numReadDuplicates[i].second) {
            //         assignedReads[nodes[maxIdx]].second.insert(readIndex);
            //     }
            // }

        }

}

std::unordered_map<std::string, std::pair<double, double>> mgsr::getReadAssignmentAccuracy(
    const std::unordered_map<std::string, std::unordered_set<size_t>>& assignedReads,
    const std::vector<std::string>& nodes, const std::vector<std::string>& readNames,
    const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor
    ) {
    //                                        correct total
    std::unordered_map<std::string, std::pair<double, double>> readAssignmentAccuracy;
    for (size_t i = 0; i < readNames.size(); ++i) {
        std::vector<std::string> split;
        PangenomeMAT::stringSplit(readNames[i], '_', split);
        std::string trueNode = split[0];
        if (leastRecentIdenticalAncestor.find(trueNode) != leastRecentIdenticalAncestor.end()) {
            const std::string& identicalAncestor = leastRecentIdenticalAncestor.at(trueNode);
            readAssignmentAccuracy[identicalAncestor].second += 1;
            if (assignedReads.find(identicalAncestor) != assignedReads.end() && assignedReads.at(identicalAncestor).find(i) != assignedReads.at(identicalAncestor).end()) {
                readAssignmentAccuracy[identicalAncestor].first += 1;
            }
        } else {
            readAssignmentAccuracy[trueNode].second += 1;
            if (assignedReads.find(trueNode) != assignedReads.end() && assignedReads.at(trueNode).find(i) != assignedReads.at(trueNode).end()) {
                readAssignmentAccuracy[trueNode].first += 1;
            }
        }
    }

    return readAssignmentAccuracy;
}
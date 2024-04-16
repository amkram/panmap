#include <cassert>
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"


void mgsr::accio(PangenomeMAT::Tree *T, std::ifstream& indexFile, size_t k, size_t l) {
    std::cout << "What's my purpose\nYou pass butter" << std::endl;
}

void resolveSeedmerIndexConflict(mgsr::seedmers& seedmerIndex, hash_t h, mgsr::seedmer** previous) {
    if (seedmerIndex.seedmerMap[h].num > 1) {
        /* hash already has 2 or more counts, aka already resolved */
        seedmerIndex.seedmerMap[h].num++;
    } else if (seedmerIndex.firstSeedmer == *previous) {
        /* previous == start seedmer */
        assert(seedmerIndex.firstSeedmer->hash == h);
        seedmerIndex.seedmerMap[h].num++;
        seedmerIndex.firstSeedmer = nullptr;
        *previous = nullptr;
    } else {
        /* hash has only one count and needs to resolved */
        if (h == seedmerIndex.firstSeedmer->hash) {
            // conflicts with start seedmer
            seedmerIndex.firstSeedmer->num++;
            seedmerIndex.firstSeedmer = seedmerIndex.firstSeedmer->next;
            seedmerIndex.firstSeedmer->prev->next = nullptr;
            seedmerIndex.firstSeedmer->prev = nullptr;
        } else if (h == (*previous)->hash) {
            // conflicts with previous seedmer
            mgsr::seedmer* tmpptr = *previous;
            tmpptr->num++;
            *previous = tmpptr->prev;
            (*previous)->next = nullptr;
            tmpptr->prev = nullptr;
        } else {
            // conflicts with intermediate seedmer
            mgsr::seedmer* tmpptr = &seedmerIndex.seedmerMap[h];
            tmpptr->num++;
            tmpptr->prev->next = tmpptr->next;
            tmpptr->next->prev = tmpptr->prev;
            tmpptr->prev = nullptr;
            tmpptr->next = nullptr;
        }
    }
}

void buildLocalSeedmers(mgsr::seedmers& localSeedmerIndex, const std::pair<std::pair<size_t, size_t>,std::vector<size_t>>& group, const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, const int32_t l, const int32_t k) {
    // build chain from left unaffected to right unaffected...  
    for (const size_t idx : group) {
        const std::tuple<size_t, int32_t, int32_t, bool>& syncmer = syncmerChanges[idx];
        if (std::get<3>(syncmer)) continue;

    }
}

std::vector<std::pair<std::pair<size_t, size_t>,std::vector<size_t>>> getGroups(const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, const int32_t l) {
    std::vector<std::pair<std::pair<size_t, size_t>,std::vector<size_t>>> groupsIndex;
    if (syncmerChanges.empty()) return groupsIndex;

    std::pair<std::pair<size_t, size_t>,std::vector<size_t>> curGroup;
    auto prevlb = curSeeds.begin();
    auto curlb  = curSeeds.begin();

    for (size_t i = 0; i < syncmerChanges.size(); ++i) {
        const auto& syncmer = syncmerChanges[i];
        if (curGroup.second.empty()) {
            curGroup.second.push_back(i);
            prevlb = curSeeds.lower_bound(std::get<1>(syncmer));
            auto leftUnaffected = prevlb;
            for (int i = 0; i < l && leftUnaffected != curSeeds.begin(); ++i) --leftUnaffected;
            curGroup.first.first = leftUnaffected->first;
            continue;
        }

        curlb = curSeeds.lower_bound(std::get<1>(syncmer));
        if (prevlb == curlb) {
            curGroup.second.push_back(i);
        } else {
            size_t numSteps = 0;
            auto prevlbtmp = prevlb;
            while (prevlbtmp != curlb) {
                if (numSteps > l) break;
                ++prevlbtmp;
                ++numSteps;
            }
            if (numSteps <= l) {
                curGroup.second.push_back(i);
            } else {
                auto rightUnaffected = prevlb;
                if (prevlb->first == std::get<1>(syncmerChanges[curGroup.second.back()])) ++rightUnaffected;
                curGroup.first.second = rightUnaffected->first;
                
                groupsIndex.push_back(curGroup);
                curGroup.second.clear();
                curGroup.second.push_back(i);
                prevlb = curSeeds.lower_bound(std::get<1>(syncmer));

                auto leftUnaffected = prevlb;
                for (int i = 0; i < l && leftUnaffected != curSeeds.begin(); ++i) --leftUnaffected;
                curGroup.first.first = leftUnaffected->first;
            }
        }
        prevlb = curlb;
    }
    groupsIndex.push_back(curGroup);
    return groupsIndex;
}

/*testing*/
void initializeMap(mgsr::seedmers& seedmersIndex, const std::vector<std::tuple<hash_t, int32_t, int32_t, bool>>& syncmers, std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, const int32_t l, const int32_t k) {
    // first seedmer
    assert(std::get<3>(syncmers[0]) == false);
    curSeeds[std::get<1>(syncmers[0])] = std::make_pair(std::get<0>(syncmers[0]),std::get<2>(syncmers[0]));

    size_t cacheForwardH = 0;
    for (int i = 0; i < l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(syncmers[i]);

    size_t cacheReversedH = 0;
    for (int i = l - 1; i > -1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(syncmers[i]);

    // Skip if strand ambiguous
    if (cacheForwardH < cacheReversedH) {
        seedmersIndex.seedmerMap[cacheForwardH] = {
            cacheForwardH,              // hash
            std::get<1>(syncmers[0]),   // start position
            std::get<2>(syncmers[l-1]), // end position
            1,                          // num of identical hash
            false,                      // reversed
            nullptr,                    // ptr to previous
            nullptr                     // ptr to next
        };
        seedmersIndex.firstSeedmer = &seedmersIndex.seedmerMap[cacheForwardH];
    } else if (cacheReversedH < cacheForwardH) {
        seedmersIndex.seedmerMap[cacheReversedH] = {
            cacheReversedH,
            std::get<1>(syncmers[0]),
            std::get<2>(syncmers[l-1]),
            1,
            true,
            nullptr,
            nullptr
        };
        seedmersIndex.firstSeedmer = &seedmersIndex.seedmerMap[cacheReversedH];
    }
    
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    mgsr::seedmer* previousSeedmer = seedmersIndex.firstSeedmer;

    for (int i = 1; i < syncmers.size() - l + 1; ++i) {
        assert(std::get<3>(syncmers[i]) == false);

        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<0>(syncmers[i]),std::get<2>(syncmers[i]));

        cacheForwardH = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);
        cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(syncmers[i+l-1]) << (2 * k * (l - 1)));
        
        hash_t curHash;
        bool rev;
        if (cacheForwardH < cacheReversedH) {
            curHash = cacheForwardH;
            rev = false;
        } else if (cacheReversedH < cacheForwardH) {
            curHash = cacheReversedH;
            rev = true;
        } else {
            continue; // Skip if strand ambiguous
        }

        if (seedmersIndex.seedmerMap.count(curHash) > 0) {
            resolveSeedmerIndexConflict(seedmersIndex, curHash, &previousSeedmer);
        } else {
            seedmersIndex.seedmerMap[curHash] = {
                curHash,
                std::get<1>(syncmers[i]),
                std::get<2>(syncmers[i+l-1]),
                1,
                rev,
                nullptr,
                nullptr
            };

            if (seedmersIndex.firstSeedmer != nullptr) {
                assert(previousSeedmer != nullptr);
                seedmersIndex.seedmerMap[curHash].prev = previousSeedmer;
                previousSeedmer->next = &seedmersIndex.seedmerMap[curHash];
            } else {
                seedmersIndex.firstSeedmer = &seedmersIndex.seedmerMap[curHash];
            }
            previousSeedmer = &seedmersIndex.seedmerMap[curHash];
        }
    }

    /*rest of syncmers*/
    for (int i = syncmers.size() - l + 1; i < syncmers.size(); i++) {
        assert(std::get<3>(syncmers[i]) == false);
        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<0>(syncmers[i]),std::get<2>(syncmers[i]));
    }

}

size_t btn(char b) {
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

size_t hash(const std::string& s) {
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

char comp(char c) {
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

std::string revcomp(const std::string& s) {
    std::string cs = "";
    for (int i = s.size() - 1; i > -1; --i) {
        char c = s[i];
        cs += comp(c);
    }
    return cs;
}

std::pair<size_t, bool> getHash(const std::string& s) {
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

void mgsr::buildSeedmerHelper(tree::mutableTreeData &data, seedMap_t seedMap, pmi::seedIndex &index, mgsr::seedmers seedmersIndex, std::map<int32_t, std::pair<size_t, int32_t>> curSeeds, PangenomeMAT::Tree *T, const PangenomeMAT::Node *node, const int32_t l, const int32_t k, const int32_t s, const tree::globalCoords_t &globalCoords, std::stringstream& seedmersOutStream) {
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
        while (seen < k && blockStart >= 1) {
            if (data.gappedConsensus[blockStart] != '-') {
                seen++;
            }
            blockStart--;
        }

        int32_t blockStop = tree::getGlobalCoordinate(blockId, globalCoords[blockId].first.size()-1, -1, globalCoords);
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

            /*testing*/
            std::pair<size_t, bool> kmerHash = getHash(kmer);
            if (!kmerHash.second) continue;
            size_t curHash = kmerHash.first;
            /*testing*/

            if (seen.find(c) != seen.end()) {
                continue;
            }
            seen[c] = true;
            if (seedMap.find(c) != seedMap.end()) {
                // This kmer is already a seed.
                std::string prevseedmer = seedMap[c].second;
                if (seeding::is_syncmer(kmer, s, false)) {
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
                if (seeding::is_syncmer(kmer, s, false)) {
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

    if (seedmersIndex.seedmerMap.size() == 0) {
        /*build seedmer for root*/
        std::sort(syncmerChanges.begin(), syncmerChanges.end(), [](const auto &a, const auto &b) {
            return std::get<1>(a) < std::get<1>(b);
        });

        initializeMap(seedmersIndex, syncmerChanges, curSeeds, l, k);

        seedmersOutStream << node->identifier << " ";
        mgsr::seedmer* curSeedmer = seedmersIndex.firstSeedmer;
        while (curSeedmer != nullptr) {
            std::string s = std::to_string(curSeedmer->hash) + ","
                          + std::to_string(curSeedmer->beg) + ","
                          + std::to_string(curSeedmer->end) + ","
                          + std::to_string(curSeedmer->rev);

            seedmersOutStream << s << " ";
            curSeedmer = curSeedmer->next;
        }
        seedmersOutStream << "\n";
    } 
    else {
        if (node->identifier == "node_2") {
        /*update seedmer and curSeeds*/
    
        std::sort(syncmerChanges.begin(), syncmerChanges.end(), [](const auto &a, const auto &b) {
            return std::get<1>(a) < std::get<1>(b);
        });

        // group changes in syncmers
            // changes within each group has syncmer distance < l-1
            // each group separated by l-1 or more syncmers

        //vector< pair< pair<leftUnaffect, rightUnaffcted>, vector<index> > >
        std::cout << syncmerChanges.size() << std::endl;         
        std::vector<std::pair<std::pair<size_t, size_t>,std::vector<size_t>>> groupsIndex = getGroups(syncmerChanges, curSeeds, l);
        std::cout << node->identifier << std::endl;
        for (const auto& group : groupsIndex) {
            std::cout << group.first.first << ":" << group.first.second << "\t";
            for (const auto& idx : group.second) {
                std::cout << idx << " ";
            }
            std::cout << std::endl;
        }

        // for each group
        for (const std::pair<std::pair<size_t, size_t>,std::vector<size_t>>& group : groupsIndex) {
            // apply changes to seedmerIndex
                // make new local seedmer chain in group
                    // resolve conflicts if ensued
            mgsr::seedmers localSeedmerIndex;
            buildLocalSeedmers(localSeedmerIndex, group, syncermerChanges, curSeeds, l, k); // To speed things up: maybe instead of making localSeedmerIndex just update the curSeeds and record the range that's changed 
                // identify first unaffected seedmer on the left
                // identify first unaffected ssedmer on the right
                // remove old local seedmer chain that will be replaced
                    // if a removed seedmer has same hash elsewhere and becomes unique after removal
                        // restore hash that's elsewhere and record change in stringstream
                    // else
                        // erase
                // merge local seedmer chain to the main chain and record change in stringstream
                    // resolve conflicts if ensued and record change in stringstream
            // apply changes to  curSeeds
            // outstream changes to seedmerOutStream
        }


        seedmersOutStream << node->identifier << " ";
        }
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
    std::map<int32_t, std::pair<size_t, int32_t>> curSeeds; // map[beg] = {hash, end}
    index.outStream << k << " " << s << " " << l << "\n";
    seedmersOutStream << k << " " << s << " " << l << "\n";
    /* Recursive traversal of tree to build the index */
    mgsr::buildSeedmerHelper(data, seedMap, index, seedmersIndex, curSeeds, T, T->root, l, k, s, globalCoords, seedmersOutStream);
}
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

void resolveSeedmersIndexConflict(mgsr::seedmers& seedmersIndex, hash_t h, std::pair<size_t, bool>& previous, std::pair<size_t, bool>& endSeedmer) {
    std::cout << "h: " << h << "\th_num: " << seedmersIndex.seedmerMap[h].num << "\tprevious: " << previous.first << std::endl; 
    if (seedmersIndex.seedmerMap[h].num > 1) {
        /* hash already has 2 or more counts, aka already resolved */
        std::cout << "--rsi-A--" << std::endl;
        seedmersIndex.seedmerMap[h].num++;
        if (previous.second == true) seedmersIndex.seedmerMap[h].prev.insert(previous.first);
        else seedmersIndex.seedmerMap[h].prev.insert(std::numeric_limits<size_t>::max());
    } else if (seedmersIndex.firstSeedmer.second != false && seedmersIndex.firstSeedmer.first == h && seedmersIndex.firstSeedmer.first == previous.first) {
        /* previous == start seedmer */
        std::cout << "--rsi-B--" << std::endl;
        assert(seedmersIndex.firstSeedmer.first == h);
        seedmersIndex.seedmerMap[h].num++;
        seedmersIndex.seedmerMap[h].prev.insert(previous.first);
        seedmersIndex.firstSeedmer.second = false;
    } else {
        /* hash has only one count and needs to resolved */
        if (seedmersIndex.firstSeedmer.second != false && h == seedmersIndex.firstSeedmer.first) {
            // conflicts with start seedmer -> start seedmer loses prev and next.. Practically disconnected until rescued
            std::cout << "--rsi-C--" << std::endl;
            auto firstSeedmerIt = seedmersIndex.seedmerMap.find(seedmersIndex.firstSeedmer.first);
            firstSeedmerIt->second.num++;
            firstSeedmerIt->second.prev.insert(previous.first);
            seedmersIndex.firstSeedmer.first = firstSeedmerIt->second.next.first;

            firstSeedmerIt = seedmersIndex.seedmerMap.find(seedmersIndex.firstSeedmer.first);
            assert(firstSeedmerIt->second.num == 1);
            seedmersIndex.seedmerMap[*(firstSeedmerIt->second.prev.begin())].next = {0, false};
            firstSeedmerIt->second.prev.clear();
            firstSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
        } else if (previous.second != false && h == previous.first) {
            // conflicts with previous seedmer
            std::cout << "--rsi-D--" << std::endl;
            auto previousSeedmerIt = seedmersIndex.seedmerMap.find(previous.first);
            assert(previousSeedmerIt->second.num == 1);
            previous.first = *(previousSeedmerIt->second.prev.begin());

            previousSeedmerIt = seedmersIndex.seedmerMap.find(previous.first);
            seedmersIndex.seedmerMap[previousSeedmerIt->second.next.first].prev.insert(previousSeedmerIt->second.next.first);
            seedmersIndex.seedmerMap[previousSeedmerIt->second.next.first].num++;
            previousSeedmerIt->second.next = {0, false};
        } else if (seedmersIndex.lastSeedmer.second != false && h == seedmersIndex.lastSeedmer.first) {
            // conflict with last seedmer
            std::cout << "--rsi-E--" << std::endl;
            auto lastSeedmerIt = seedmersIndex.seedmerMap.find(seedmersIndex.lastSeedmer.first);
            std::cout << "prev size: " << lastSeedmerIt->second.prev.size() << "\t"
                      << "beg size: "  << lastSeedmerIt->second.begs.size() << "\t"
                      << "num: "       << lastSeedmerIt->second.num << std::endl;
            std::cout << "prevs: ";
            for (const auto& p : lastSeedmerIt->second.prev)  std::cout << p << "\t";
            std::cout << std::endl;
            std::cout << "begs: ";
            for (const auto& b : lastSeedmerIt->second.begs)  std::cout << b << "\t";
            std::cout << std::endl;
            assert(lastSeedmerIt->second.num == 1);
            lastSeedmerIt->second.num++;
            // assert(lastSeedmerIt->second.prev.size() == 1);

            if (endSeedmer.second == true && endSeedmer.first == seedmersIndex.lastSeedmer.first) {
                std::cout << "endSeedmer.second == true && endSeedmer.first == seedmersIndex.lastSeedmer.first" << std::endl;
                if (previous.second == true) lastSeedmerIt->second.prev.insert(previous.first);
                else lastSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
                endSeedmer.second = false;
            } else {
                std::cout << "NOT endSeedmer.second == true && endSeedmer.first == seedmersIndex.lastSeedmer.first" << std::endl;
                std::cout << "prev size: " << lastSeedmerIt->second.prev.size() << "\t"
                          << "beg size: "  << lastSeedmerIt->second.begs.size() << "\t"
                          << "num: "       << lastSeedmerIt->second.num << std::endl;
                seedmersIndex.lastSeedmer.first = *(lastSeedmerIt->second.prev.begin());
                if (previous.second == true) lastSeedmerIt->second.prev.insert(previous.first);
                else lastSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
                lastSeedmerIt = seedmersIndex.seedmerMap.find(seedmersIndex.lastSeedmer.first);
                lastSeedmerIt->second.next = {0, false};
                std::cout << "new lastSeedmerIt: " << lastSeedmerIt->second.hash << "\tnext:" << lastSeedmerIt->second.next.first << std::endl;
            }
        } else if (endSeedmer.second != false && h == endSeedmer.first) {
            // conflicts with end seedmer
            std::cout << "--rsi-F--" << std::endl;
            auto endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmer.first);
            endSeedmerIt->second.num++;
            if (previous.second == true) endSeedmerIt->second.prev.insert(previous.first);
            else endSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
            endSeedmer.first = endSeedmerIt->second.next.first;
            endSeedmerIt->second.next = {0, false};
            endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmer.first);
            endSeedmerIt->second.prev.clear();
        } else {
            std::cout << "--rsi-G--" << std::endl;
            // conflicts with intermediate seedmer
            auto curSeedmerIt  = seedmersIndex.seedmerMap.find(h);
            auto prevSeedmerIt = seedmersIndex.seedmerMap.find(*(curSeedmerIt->second.prev.begin()));
            auto nextSeedmerIt = seedmersIndex.seedmerMap.find(curSeedmerIt->second.next.first);
            prevSeedmerIt->second.next.first = nextSeedmerIt->second.hash;
            nextSeedmerIt->second.prev.clear();
            nextSeedmerIt->second.prev.insert(prevSeedmerIt->second.hash);
            curSeedmerIt->second.num++;
            if (previous.second == true) curSeedmerIt->second.prev.insert(previous.first);
            else curSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
            curSeedmerIt->second.next = {0, false};
        }
    }
}

void removeSeedmersBetweenLocators(std::pair<size_t, bool>& begSeedmer, std::pair<size_t, bool>& endSeedmer, const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>& locators, mgsr::seedmers& seedmersIndex, std::stringstream& seedmersOutStream, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, const int32_t l, const int32_t k) {
    std::cout << "start to remove seedmers between locators" << std::endl;
    auto endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmer.first);
    // auto begSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmer.first);

    decltype(seedmersIndex.positionMap.begin()) curSeedmerIt;
    if (endSeedmer.second == false) {
        curSeedmerIt = seedmersIndex.positionMap.end();
        std::cout << "endSeedmer.second == false" << std::endl;
    } else {
        curSeedmerIt = seedmersIndex.positionMap.find(*(endSeedmerIt->second.begs.begin()));
        if (curSeedmerIt == seedmersIndex.positionMap.begin()) return;
    }
    
    if (begSeedmer.second) {
        std::cout << "beSeedmer begs size: " <<  seedmersIndex.seedmerMap[begSeedmer.first].begs.size() << std::endl;
        std::cout << "begSeedmer\t" << *(seedmersIndex.seedmerMap[begSeedmer.first].begs.begin()) << std::endl;
    }
    std::cout << "begSeedmerHash\t" << begSeedmer.first << "\t" << begSeedmer.second << std::endl;
    std::cout << "endSeedmerHash\t" << endSeedmer.first << "\t" << endSeedmer.second << std::endl;

    std::cout << "curSeedmerIt: " << curSeedmerIt->first << "\t" << curSeedmerIt->second.first << "\t" << curSeedmerIt->second.second << std::endl;
    std::cout << "foo" << std::endl;
    --curSeedmerIt;
    std::cout << "boo" << std::endl;
    std::cout << "curSeedmerIt: " << curSeedmerIt->first << "\t" << curSeedmerIt->second.first << "\t" << curSeedmerIt->second.second << std::endl;

    while (true) {
        std::cout << "--2--" << std::endl;
        if (begSeedmer.second != false && seedmersIndex.seedmerMap[curSeedmerIt->second.second].hash == begSeedmer.first) break;
        auto curSeedmer = &seedmersIndex.seedmerMap[curSeedmerIt->second.second];
        std::cout << "inside while loop, curSeedmer: " << curSeedmer->hash << "\tnum: "<< curSeedmer->num << std::endl;
        std::cout << "--3--" << std::endl;
        if (curSeedmer->num == 1) {
            std::cout << "--CCCCC--" << std::endl;
            std::cout << "positionMap begin: " << seedmersIndex.positionMap.begin()->first << std::endl;
            if (curSeedmerIt == seedmersIndex.positionMap.begin() || curSeedmer->hash == seedmersIndex.firstSeedmer.first) {
                std::cout << "1\t" << curSeedmerIt->first << "\t" << curSeedmerIt->second.second << std::endl;
                assert(*(curSeedmer->prev.begin()) == std::numeric_limits<size_t>::max());
                if (curSeedmer->next.second == true) {
                    auto nextSeedmerIt = seedmersIndex.seedmerMap.find(curSeedmer->next.first);
                    nextSeedmerIt->second.prev.clear();
                    nextSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
                    std::cout << "2\t" << nextSeedmerIt->first << "\t" << nextSeedmerIt->second.prev.size() << std::endl;
                }

            } else {
                std::cout << curSeedmer->prev.size() << std::endl;
                for (const auto& p : curSeedmer->prev) std::cout << p << "\t";
                std::cout << std::endl;
                assert(curSeedmer->prev.size() == 1);
                std::cout << curSeedmerIt->first << std::endl;
                auto prevSeedmerIt = seedmersIndex.seedmerMap.find(*(curSeedmer->prev.begin()));
                if (curSeedmer->next.second == true) {
                    auto nextSeedmerIt = seedmersIndex.seedmerMap.find(curSeedmer->next.first);
                    prevSeedmerIt->second.next.first = nextSeedmerIt->second.hash;
                    nextSeedmerIt->second.prev.clear();
                    nextSeedmerIt->second.prev.insert(prevSeedmerIt->second.hash);
                } else {
                    prevSeedmerIt->second.next = {0, false};
                }
            }
            assert(seedmersIndex.seedmerMap.erase(curSeedmer->hash));
        } else if (curSeedmer->num == 2) {
            std::cout << "--DDDDD--" << std::endl;
            for (const auto& it : curSeedmer->begs) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
            for (const auto& it : curSeedmer->prev) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
            assert(curSeedmer->begs.size() == 2);

            // erase closest previous
            int32_t toEraseBeg = curSeedmerIt->first;
            auto toErase = seedmersIndex.positionMap.find(toEraseBeg);
            assert(toErase != seedmersIndex.positionMap.end());
            
            size_t prevToErase;
            if (toErase == seedmersIndex.positionMap.begin()) {
                prevToErase = std::numeric_limits<size_t>::max();
            } else {
                auto prevIt = toErase;
                --prevIt;
                while (true) {
                    std::cout << "cur prevIt: " << prevIt->first << "\t" << prevIt->second.second << std::endl; 
                    if (curSeedmer->prev.find(prevIt->second.second) != curSeedmer->prev.end()) {
                        prevToErase = prevIt->second.second;
                        break;
                    }
                    --prevIt;
                }
            }
            
            if (curSeedmer->prev.size() < curSeedmer->begs.size()) {
                assert(prevToErase == *(curSeedmer->prev.begin()));
            } else if (curSeedmer->prev.size() == curSeedmer->begs.size()) {
                assert(curSeedmer->prev.erase(prevToErase));
            } else {
                throw std::runtime_error("seedmer.prev.size() > seedmer.begs.size()... something went wrong");
            }
            assert(curSeedmer->begs.erase(toEraseBeg));

            if (*(curSeedmer->prev.begin()) == std::numeric_limits<size_t>::max() || *(curSeedmer->prev.begin()) == 0) {
                // rescuing start seedmer
                std::cout << "rescuring start seedmer" << std::endl;
                assert(seedmersIndex.firstSeedmer.second == true);
                auto newFirstSeedmerIt = seedmersIndex.positionMap.find(*(curSeedmer->begs.begin()));
                auto newFirstSeedmer = seedmersIndex.seedmerMap.find(newFirstSeedmerIt->second.second);
                auto oldFirstSeedmer = seedmersIndex.seedmerMap.find(seedmersIndex.firstSeedmer.first);
                newFirstSeedmer->second.next = {oldFirstSeedmer->second.hash, true};
                oldFirstSeedmer->second.prev.clear();
                oldFirstSeedmer->second.prev.insert(newFirstSeedmer->second.hash);
                seedmersIndex.firstSeedmer.first = newFirstSeedmer->second.hash;
            } else {
                auto prevSeedmerIt = seedmersIndex.positionMap.find(*(curSeedmer->begs.begin()));
                --prevSeedmerIt;
                auto leftSeedmer = seedmersIndex.seedmerMap.find(prevSeedmerIt->second.second);
                while (leftSeedmer->second.next.second == false) {
                    --prevSeedmerIt;
                    leftSeedmer = seedmersIndex.seedmerMap.find(prevSeedmerIt->second.second);
                }
                auto rightSeedmer = seedmersIndex.seedmerMap.find(leftSeedmer->second.next.first);
                leftSeedmer->second.next = {curSeedmer->hash, true};
                curSeedmer->next = {rightSeedmer->second.hash, true};
                curSeedmer->prev.clear();
                curSeedmer->prev.insert(leftSeedmer->second.hash);
                rightSeedmer->second.prev.clear();
                rightSeedmer->second.prev.insert(curSeedmer->hash);
            }
            auto curSeed = curSeeds.find(curSeedmerIt->first);
            int32_t count = 0;
            while (count < l-1) {
                ++curSeed;
                ++count;
            }
            curSeedmer->end = curSeed->second.second;
            curSeedmer->num--;

            for (const auto& it : curSeedmer->begs) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
            for (const auto& it : curSeedmer->prev) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "--EEEEE--" << std::endl;
            // erase closest previous
            std::cout << curSeedmerIt->first << "\t" << curSeedmerIt->second.second << std::endl;
            for (const auto& it : curSeedmer->begs) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
            for (const auto& it : curSeedmer->prev) {
                std::cout << it << " ";
            }
            std::cout << std::endl;

            int32_t toEraseBeg = curSeedmerIt->first;

            auto toErase = seedmersIndex.positionMap.find(toEraseBeg);
            assert(toErase != seedmersIndex.positionMap.end());
            std::cout << "toErase: " << toErase->first << "\t" << toErase->second.second << std::endl; 
            
            if (toErase == seedmersIndex.positionMap.begin()) std::cout << "Got you1!" << std::endl;
            std::cout << "first seedmer hash: " << seedmersIndex.firstSeedmer.first << std::endl;
            
            size_t prevToErase;
            if (toErase == seedmersIndex.positionMap.begin()) {
                prevToErase = std::numeric_limits<size_t>::max();
            } else {
                auto prevIt = toErase;
                --prevIt;
                while (true) {
                    std::cout << "cur prevIt: " << prevIt->first << "\t" << prevIt->second.second << std::endl; 
                    if (curSeedmer->prev.find(prevIt->second.second) != curSeedmer->prev.end()) {
                        prevToErase = prevIt->second.second;
                        break;
                    }
                    --prevIt;
                }
            }

            
            std::cout << "foo\t1" << std::endl;
            if (curSeedmer->prev.size() == curSeedmer->begs.size()) {
                assert(curSeedmer->prev.erase(prevToErase));
            } else {
                assert(curSeedmer->prev.size() < curSeedmer->begs.size());
                std::vector<size_t> prevHashes;
                prevHashes.reserve(curSeedmer->begs.size());
                auto itBeg = curSeedmer->begs.begin();
                while (itBeg != curSeedmer->begs.end()) {
                    bool prevFound = false;
                    int32_t curBeg = *itBeg;
                    std::cout << "curBeg: " << curBeg << std::endl;
                    auto curPrevIt = seedmersIndex.positionMap.find(curBeg);
                    std::cout << "cur prevIt: " << curPrevIt->first << "\t" << curPrevIt->second.second << std::endl; 
                    if (curPrevIt != seedmersIndex.positionMap.begin()) {
                        --curPrevIt;
                        while (true) {
                        std::cout << "cur prevIt: " << curPrevIt->first << "\t" << curPrevIt->second.second << std::endl; 
                            if (curSeedmer->prev.find(curPrevIt->second.second) != curSeedmer->prev.end()) {
                                prevHashes.push_back(curPrevIt->second.second);
                                prevFound = true;
                                break;
                            }
                            if (curPrevIt == seedmersIndex.positionMap.begin()) break;
                            --curPrevIt;
                        }
                    }
                    if (!prevFound) prevHashes.push_back(std::numeric_limits<size_t>::max());
                    ++itBeg;
                }
                std::cout << "foo\t2" << std::endl;
                for (const auto ph : prevHashes) std::cout << ph << "\t";
                std::cout << std::endl;
                auto numHash = std::count(prevHashes.begin(), prevHashes.end(), prevToErase);
                if (numHash == 1) assert(curSeedmer->prev.erase(prevToErase));
            }
            assert(curSeedmer->begs.erase(toEraseBeg));
            curSeedmer->num--;

            for (const auto& it : curSeedmer->begs) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
            for (const auto& it : curSeedmer->prev) {
                std::cout << it << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "--FFFFF--" << std::endl;
        auto begToErase = curSeedmerIt->first;
        std::cout << "seedmer at " << begToErase << ": " << curSeedmerIt->second.second << " erased" << std::endl;
        std::cout << std::endl;
        if (begSeedmer.second == false && curSeedmerIt == seedmersIndex.positionMap.begin()) {
            assert(seedmersIndex.positionMap.erase(begToErase));
            std::cout << "Break!" << std::endl;
            break;
        }
        --curSeedmerIt;
        assert(seedmersIndex.positionMap.erase(begToErase));
    }
}


void mergeSeedmerIndex(const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& newSeedmers, const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>& locators, mgsr::seedmers& seedmersIndex, std::stringstream& seedmersOutStream, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, const int32_t l, const int32_t k) {
    if (locators.first.second && locators.second.second) {
        std::cout << "--AAAAA--" << std::endl;
        /* Remove seedmers between locators on seedmerIndex */
        auto [begSeedmerHash, endSeedmerHash] = locators;
        assert(seedmersIndex.seedmerMap.find(begSeedmerHash.first)->second.num == 1 && seedmersIndex.seedmerMap.find(endSeedmerHash.first)->second.num == 1);

        removeSeedmersBetweenLocators(begSeedmerHash, endSeedmerHash, locators, seedmersIndex, seedmersOutStream, curSeeds, l, k);
        std::cout << "Start merging now!" << std::endl;
        /* Add new seedmers */
        auto previousSeedmerHash = begSeedmerHash;
        auto endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmerHash.first);
        assert(endSeedmerIt->second.num == 1 && endSeedmerIt->second.prev.size() == 1);
        std::cout << "before merging new seedmers, endSeedmerIt: " << endSeedmerIt->second.hash << " " << endSeedmerIt->second.prev.size() << " " << endSeedmerIt->second.num << std::endl;
        endSeedmerIt->second.prev.clear();
        for (const auto& newSeedmer : newSeedmers) {
            auto [curHash, curBeg, curEnd, curRev] = newSeedmer;
            std::cout << "current new seedmer: " << curBeg << "\t" << curHash << std::endl;
            if (seedmersIndex.seedmerMap.count(curHash) > 0) {
                std::cout << "conflict when merging seedmers!\t" << curBeg << "\t" << curHash << std::endl;
                resolveSeedmersIndexConflict(seedmersIndex, curHash, previousSeedmerHash, endSeedmerHash);
                seedmersIndex.seedmerMap[curHash].begs.insert(curBeg);
            } else {
                seedmersIndex.seedmerMap[curHash] = {
                    curHash,
                    {curBeg},
                    curEnd,
                    1,
                    curRev,
                    {},
                    std::make_pair(0, false)
                };
                
                if (seedmersIndex.firstSeedmer.second != false) {
                    assert(previousSeedmerHash.second != false && seedmersIndex.seedmerMap[previousSeedmerHash.first].num == 1);
                    seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {curHash, true};
                } else {
                    seedmersIndex.firstSeedmer = std::make_pair(curHash, true);
                    seedmersIndex.seedmerMap[curHash].prev.insert(std::numeric_limits<size_t>::max());
                }
                if (previousSeedmerHash.second != false) seedmersIndex.seedmerMap[curHash].prev.insert(previousSeedmerHash.first);
                previousSeedmerHash = std::make_pair(curHash, true);
            }

            seedmersIndex.positionMap[curBeg] = std::make_pair(curEnd, curHash);
        }
        // what if newSeedmers.size() == 0? Need to address
        endSeedmerIt->second.prev.insert(previousSeedmerHash.first);
        if (endSeedmerHash.second == true) {
            if (endSeedmerIt->second.hash == endSeedmerHash.first) {
                std::cout << "endSeemderHash did NOT change" << std::endl;
                seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {endSeedmerHash.first, true};
            } else {
                std::cout << "You've changed :(" << std::endl;
                assert(endSeedmerIt->second.num > 1);
                endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmerHash.first);
                seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {endSeedmerHash.first, true};
                endSeedmerIt->second.prev.insert(previousSeedmerHash.first);
            }
            std::cout << endSeedmerIt->second.prev.size() << "\t" << endSeedmerIt->second.num << std::endl;
            assert(endSeedmerIt->second.prev.size() == 1 && endSeedmerIt->second.num == 1);
        } else {
            // endSeedmer == lastSeedmer
            auto newLastSeedmerIt = seedmersIndex.seedmerMap.find(previousSeedmerHash.first);
            seedmersIndex.lastSeedmer.first = newLastSeedmerIt->second.hash;
            newLastSeedmerIt->second.next = {0, false};
        }

    } else if (locators.first.second) {
        std::cout << "--HHHHH--" << std::endl;
        // new seedmers at end
        auto begSeedmerHash = locators.first;
        std::pair<size_t, bool> endSeedmerHash = {0, false};
        assert(seedmersIndex.seedmerMap.find(begSeedmerHash.first)->second.num == 1);
        removeSeedmersBetweenLocators(begSeedmerHash, endSeedmerHash, locators, seedmersIndex, seedmersOutStream, curSeeds, l, k);
        std::cout << "--GGGGG--" << std::endl;
        seedmersIndex.lastSeedmer = {0, false};

        
        /* Add new seedmers */
        auto previousSeedmerHash = begSeedmerHash;
        for (const auto& newSeedmer : newSeedmers) {
            auto [curHash, curBeg, curEnd, curRev] = newSeedmer;
            if (seedmersIndex.seedmerMap.count(curHash) > 0) {
                resolveSeedmersIndexConflict(seedmersIndex, curHash, previousSeedmerHash, endSeedmerHash);
                seedmersIndex.seedmerMap[curHash].begs.insert(curBeg);
            } else {
                seedmersIndex.seedmerMap[curHash] = {
                    curHash,
                    {curBeg},
                    curEnd,
                    1,
                    curRev,
                    {},
                    std::make_pair(0, false)
                };
                
                if (seedmersIndex.firstSeedmer.second != false) {
                    assert(previousSeedmerHash.second != false && seedmersIndex.seedmerMap[previousSeedmerHash.first].num == 1);
                    seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {curHash, true};
                } else {
                    seedmersIndex.firstSeedmer = std::make_pair(curHash, true);
                    seedmersIndex.seedmerMap[curHash].prev.insert(std::numeric_limits<size_t>::max());
                }
                if (previousSeedmerHash.second != false) seedmersIndex.seedmerMap[curHash].prev.insert(previousSeedmerHash.first);
                previousSeedmerHash = std::make_pair(curHash, true);
            }

            seedmersIndex.positionMap[curBeg] = std::make_pair(curEnd, curHash);
        }
        seedmersIndex.lastSeedmer = previousSeedmerHash;
    } else if (locators.second.second) {
        // new seedmers at beginning
        std::cout << "--IIIII--" << std::endl;
        auto endSeedmerHash = locators.second;
        std::pair<size_t, bool> begSeedmerHash = {0, false};
        assert(seedmersIndex.seedmerMap.find(endSeedmerHash.first)->second.num == 1);
        removeSeedmersBetweenLocators(begSeedmerHash, endSeedmerHash, locators, seedmersIndex, seedmersOutStream, curSeeds, l, k);
        std::cout << "--GGGGG--" << std::endl;
        seedmersIndex.firstSeedmer = {0, false};
        
        /* Add new seedmers */
        std::pair<size_t, bool> previousSeedmerHash = {0, false};
        auto endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmerHash.first);
        assert(endSeedmerIt->second.num == 1);
        endSeedmerIt->second.prev.clear();
        for (const auto& newSeedmer : newSeedmers) {
            auto [curHash, curBeg, curEnd, curRev] = newSeedmer;
            std::cout << "new seed: " << curHash << "\t" << curBeg << std::endl;
            if (seedmersIndex.seedmerMap.count(curHash) > 0) {
                resolveSeedmersIndexConflict(seedmersIndex, curHash, previousSeedmerHash, endSeedmerHash);
                seedmersIndex.seedmerMap[curHash].begs.insert(curBeg);
            } else {
                seedmersIndex.seedmerMap[curHash] = {
                    curHash,
                    {curBeg},
                    curEnd,
                    1,
                    curRev,
                    {},
                    std::make_pair(0, false)
                };
                
                if (seedmersIndex.firstSeedmer.second != false) {
                    assert(previousSeedmerHash.second != false && seedmersIndex.seedmerMap[previousSeedmerHash.first].num == 1);
                    seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {curHash, true};
                } else {
                    seedmersIndex.firstSeedmer = std::make_pair(curHash, true);
                    seedmersIndex.seedmerMap[curHash].prev.insert(std::numeric_limits<size_t>::max());
                }
                if (previousSeedmerHash.second != false) seedmersIndex.seedmerMap[curHash].prev.insert(previousSeedmerHash.first);
                previousSeedmerHash = std::make_pair(curHash, true);
            }

            seedmersIndex.positionMap[curBeg] = std::make_pair(curEnd, curHash);
        }
        std::cout << "Chain end" << std::endl;
        if (previousSeedmerHash.second == true) {
            std::cout << "newSeedmers.size() > 0" << std::endl;
            endSeedmerIt->second.prev.insert(previousSeedmerHash.first);
            if (endSeedmerIt->second.hash == endSeedmerHash.first) {
                seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {endSeedmerHash.first, true};
            } else {
                assert(endSeedmerIt->second.num > 1);
                endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmerHash.first);
                seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {endSeedmerHash.first, true};
                endSeedmerIt->second.prev.insert(previousSeedmerHash.first);
            }
            std::cout << endSeedmerIt->second.prev.size() << std::endl;
            std::cout << endSeedmerIt->second.num << std::endl;
            assert(endSeedmerIt->second.prev.size() == 1 && endSeedmerIt->second.num == 1);
        } else {
            endSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
            if (endSeedmerIt->second.hash == endSeedmerHash.first) {
                seedmersIndex.firstSeedmer = endSeedmerHash;
                assert(endSeedmerIt->second.prev.size() == 1 && endSeedmerIt->second.num == 1);
            } else {
                assert(endSeedmerIt->second.num > 1);
                endSeedmerIt = seedmersIndex.seedmerMap.find(endSeedmerHash.first);
                assert(endSeedmerIt->second.num == 1);
                endSeedmerIt->second.prev.clear();
                endSeedmerIt->second.prev.insert(std::numeric_limits<size_t>::max());
                seedmersIndex.firstSeedmer = endSeedmerHash;
            }
        }
    } else {
        throw std::runtime_error("Can't find locators to insert new seedmers");
    }

    // std::cout << "Great Success! - Borat" << std::endl;
    // seedmersOutStream << " ";
}

std::map<int32_t, std::pair<size_t, int32_t>> getLocalSeeds(const std::pair<std::pair<std::pair<size_t, size_t>, size_t>, std::vector<size_t>>& group, const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds) {
    std::map<int32_t, std::pair<size_t, int32_t>> localSeeds;
    auto curlb = curSeeds.lower_bound(std::get<1>(syncmerChanges[group.second.front()]));
    auto endlb = curSeeds.lower_bound(std::get<1>(syncmerChanges[group.second.back() ]));

    if (curlb != curSeeds.end() && endlb == curSeeds.end()) {
        std::cout << "--JJJJJJ--" << std::endl;
        // partially end
        while (curlb != curSeeds.end()) {
            localSeeds[curlb->first] = curlb->second;
            ++curlb;
        }
    } else if (endlb != curSeeds.end()) {
        std::cout << "--KKKKK--" << std::endl;
        // partially beginning or intermediate
        if (endlb == curSeeds.begin()) {
            std::cout << "--LLLLL--" << std::endl;
            assert(curlb == endlb);
            if (endlb->first == std::get<1>(syncmerChanges[group.second.back()])) {
                localSeeds[curlb->first] = curlb->second;
            }
        } else {
            std::cout << "--MMMMM--" << std::endl;
            std::cout << endlb->first << std::endl;
            std::cout << curlb->first << std::endl;
    
            if (curlb == endlb) {
                std::cout << "coocoo" << std::endl;
                if (curlb->first == std::get<1>(syncmerChanges[group.second.back()])) {
                    localSeeds[curlb->first] = curlb->second;
                }
            } else {
                if (endlb->first != std::get<1>(syncmerChanges[group.second.back()])) --endlb;
                while (true) {
                    localSeeds[curlb->first] = curlb->second;
                    if (curlb == endlb) break;
                    ++curlb;
                }   
            }
        }
    }



    for (auto idx : group.second) {
        auto [seedHash, seedBeg, seedEnd, seedDel] = syncmerChanges[idx];
        if (seedDel) {
            assert(localSeeds.erase(seedBeg));
        } else {
            localSeeds[seedBeg] = {seedHash, seedEnd};
        }
    }
    return localSeeds;
}

std::vector<std::tuple<size_t, int32_t, int32_t, bool>> buildNewSeedmers(const mgsr::seedmers& seedmersIndex, const std::pair<std::pair<std::pair<size_t, size_t>, size_t>,std::vector<size_t>>& group, const std::map<int32_t, std::pair<size_t, int32_t>>& localSeeds, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>& locators, const int32_t l, const int32_t k) {
    std::deque<decltype(curSeeds.begin())> leftSeeds;
    std::deque<decltype(curSeeds.begin())> rightSeeds;
    std::cout << group.first.first.first << "\t" << group.first.first.second << std::endl;
    std::cout << group.first.second << std::endl;
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    // collect leftSeeds
    if (group.first.first.first != std::numeric_limits<size_t>::max()) {
        auto curSeed = curSeeds.find(group.first.first.first);
        if (curSeed == curSeeds.begin()) {
            // starting from leftend
            for (int i = 0; i < group.first.first.second; i++) {
                leftSeeds.push_back(curSeed);
                ++curSeed;
            }
        } else {
            // starting from non-leftend
            size_t cacheReversedH = 0;
            size_t cacheForwardH  = 0;
            size_t cacheMin = 0;
            size_t offset = 0;
            auto leftUnaffectedSeed = curSeed;
            while (true) {
                std::cout << "cur left unaffected seed: " << curSeed->first << std::endl;
                if (offset == 0) {
                    int32_t counter = 0;
                    while (counter < l) {
                        cacheForwardH = (cacheForwardH << (2 * k)) + curSeed->second.first;
                        leftSeeds.push_back(curSeed);
                        ++curSeed;
                        ++counter;
                    }
                    --curSeed;
                    while (counter > 0) {
                        cacheReversedH = (cacheReversedH << (2 * k)) + curSeed->second.first;
                        --curSeed;
                        --counter;
                    }
                    curSeed = leftUnaffectedSeed;
                } else {
                    cacheForwardH  = (cacheForwardH >> (2 * k)) + (curSeed->second.first << (2 * k * (l - 1)));
                    cacheReversedH = ((cacheReversedH & mask) << (k * 2)) + curSeed->second.first;
                }

                if (cacheForwardH < cacheReversedH) {
                    cacheMin = cacheForwardH;
                } else if (cacheReversedH < cacheForwardH) {
                    cacheMin = cacheReversedH;
                } else {
                    // strand ambiguous
                    std::cout << "strand ambiguous" << std::endl;
                    ++offset;
                    --curSeed;
                    leftSeeds.push_front(curSeed);
                    if (curSeed == curSeeds.begin()) break;
                    continue;
                }
                
                assert(seedmersIndex.seedmerMap.count(cacheMin) > 0);
                if (seedmersIndex.seedmerMap.at(cacheMin).num > 1) {
                    std::cout << cacheMin << " not unique!" << std::endl;                  
                    ++offset;
                    --curSeed;
                    leftSeeds.push_front(curSeed);
                    if (curSeed == curSeeds.begin()) break;
                    continue;
                }
                locators.first.first  = cacheMin;
                locators.first.second = true;
                break;
            }
        }
    }

    std::cout << "left seeds size: " << leftSeeds.size() << std::endl;
    for (const auto& seed : leftSeeds) std::cout << seed->first << " ";
    std::cout << std::endl;


    // collect rightSeeds
    if (group.first.second != std::numeric_limits<size_t>::max()) {
        auto curSeed = curSeeds.find(group.first.second);
        auto rightUnaffectedSeed = curSeed;
        size_t cacheReversedH = 0;
        size_t cacheForwardH  = 0;
        size_t cacheMin = 0;
        size_t offset = 0;
        size_t counter = 0;
        std::cout << "cur right unaffected seed: " << rightUnaffectedSeed->first << std::endl;
        while (counter < l && curSeed != curSeeds.end()) {
            ++curSeed;
            ++counter;
        }

        curSeed = rightUnaffectedSeed;
        if (counter < l) {
            while (curSeed != curSeeds.end()) {
                rightSeeds.push_back(curSeed);
                ++curSeed;
            }
        } else {
            while (true) {
                auto nextSeed = curSeed;
                ++nextSeed;
                if (offset == 0) {
                    int32_t counter = 0;
                    while (counter < l) {
                        cacheForwardH = (cacheForwardH << (2 * k)) + curSeed->second.first;
                        rightSeeds.push_back(curSeed);
                        std::cout << "in offset == 0: -----------------pushed curSeed to rightSeeds: " << curSeed->first << std::endl;
                        ++curSeed;
                        ++counter;
                    }
                    --curSeed;
                    while (counter > 0) {
                        cacheReversedH = (cacheReversedH << (2 * k)) + curSeed->second.first;
                        --curSeed;
                        --counter;
                    }
                    curSeed = rightUnaffectedSeed;
                } else {
                    std::cout << "previous cacheForward: " << cacheForwardH << std::endl;
                    std::cout << "previous cacheReversed: " << cacheReversedH << std::endl;
                    std::cout << "curSeed hash: " << curSeed->second.first << std::endl;
                    cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + nextSeed->second.first;
                    cacheReversedH = (cacheReversedH >> (2 * k)) + (nextSeed->second.first << (2 * k * (l - 1)));
                    std::cout << curSeed->first << "\t" << cacheForwardH << "\t" << cacheReversedH << std::endl;
                }
                std::cout << curSeed->first << "\t" << cacheForwardH << "\t" << cacheReversedH << std::endl;
                if (cacheForwardH < cacheReversedH) {
                    cacheMin = cacheForwardH;
                } else if (cacheReversedH < cacheForwardH) {
                    cacheMin = cacheReversedH;
                } else {
                    // strand ambiguous
                    std::cout << curSeed->first << " strand ambiguous" << std::endl;
                    ++offset;
                    ++curSeed;
                    ++nextSeed;
                    if (nextSeed == curSeeds.end()) break;
                    rightSeeds.push_back(nextSeed);
                    std::cout << "in strand ambiguous: -----------------pushed curSeed to rightSeeds: " << nextSeed->first << std::endl;
                    continue;
                }

                assert(seedmersIndex.seedmerMap.count(cacheMin) > 0);
                if (seedmersIndex.seedmerMap.at(cacheMin).num > 1) {
                    std::cout << curSeed->first << ": " << cacheMin << " not unique" << std::endl;
                    ++offset;
                    ++curSeed;
                    ++nextSeed;
                    if (nextSeed == curSeeds.end()) break;
                    rightSeeds.push_back(nextSeed);
                    std::cout << "in not unique: -----------------pushed curSeed to rightSeeds: " << nextSeed->first << std::endl;
                    continue;
                }
                locators.second.first  = cacheMin;
                locators.second.second = true;
                break;
            }
        }
    }

    std::cout << "right seeds size: " << rightSeeds.size() << std::endl;
    for (const auto& seed : rightSeeds) std::cout << seed->first << " ";
    std::cout << std::endl;

    std::vector<std::tuple<size_t, int32_t, int32_t, bool>> newSeedmers;
    std::queue<int32_t> begs;
    size_t  cacheReversedH = 0;
    size_t  cacheForwardH  = 0;
    size_t  cacheMin = 0;
    int32_t rendoff = 0;
    int32_t lstart = 0;
    int32_t count = 0;
    bool    rev = false;

    std::cout << "begin making new seeedmers" << std::endl;
    //leftSeeds
    if (locators.first.second) lstart = 1; 
    for (size_t i = lstart; i < leftSeeds.size(); ++i) {
        auto curSeed = leftSeeds[i];
        if (count < l - 1) {
            cacheForwardH  = (cacheForwardH << (2 * k)) + curSeed->second.first;
            cacheReversedH = cacheReversedH + (curSeed->second.first << (2 * k * count));
            begs.push(curSeed->first);
            ++count;
            continue;
        }

        if (count == l - 1) {
            cacheForwardH  = (cacheForwardH << (2 * k)) + curSeed->second.first;
            cacheReversedH = cacheReversedH + (curSeed->second.first << (2 * k * count));
        } else {
            cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + curSeed->second.first;
            cacheReversedH = (cacheReversedH >> (2 * k)) + (curSeed->second.first << (2 * k * (l - 1)));
        }
        
        begs.push(curSeed->first);
        ++count;

        if (cacheForwardH < cacheReversedH) {
            cacheMin = cacheForwardH;
            rev = false;
        } else if (cacheReversedH < cacheForwardH) {
            cacheMin = cacheReversedH;
            rev = true;
        } else {
            continue; // Skip if strand ambiguous
        }
        
        newSeedmers.push_back(std::make_tuple(cacheMin, begs.front(), curSeed->second.second, rev));
        begs.pop();
    }
    std::cout << "left seeds finished" << std::endl;

    //group
    auto localSeedIt = localSeeds.begin();
    while (localSeedIt != localSeeds.end()) {
        auto seedHash = localSeedIt->second.first;
        auto seedBeg  = localSeedIt->first;
        auto seedEnd  = localSeedIt->second.second;
        if (count < l - 1) {
            cacheForwardH  = (cacheForwardH << (2 * k)) + seedHash;
            cacheReversedH = cacheReversedH + (seedHash << (2 * k * count));
            begs.push(seedBeg);
            ++count;
            ++localSeedIt;
            continue;
        }

        if (count == l - 1) {
            cacheForwardH  = (cacheForwardH << (2 * k)) + seedHash;
            cacheReversedH = cacheReversedH + (seedHash << (2 * k * count));
        } else {
            cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + seedHash;
            cacheReversedH = (cacheReversedH >> (2 * k)) + (seedHash << (2 * k * (l - 1)));
        }
        begs.push(seedBeg);
        ++count;

        if (cacheForwardH < cacheReversedH) {
            cacheMin = cacheForwardH;
            rev = false;
        } else if (cacheReversedH < cacheForwardH) {
            cacheMin = cacheReversedH;
            rev = true;
        } else {
            begs.pop();
            ++localSeedIt;
            continue; // Skip if strand ambiguous
        }

        newSeedmers.push_back(std::make_tuple(cacheMin, begs.front(), seedEnd, rev));
        begs.pop();

        ++localSeedIt;
    }

    std::cout << "group finished" << std::endl;
    //rightSeeds
    if (locators.second.second) rendoff = 1;
    for (size_t i = 0; i < rightSeeds.size() - rendoff; ++i) {
        auto curSeed = rightSeeds[i];
        std::cout << "cur right seed\t" << curSeed->second.first << "\t" << curSeed->first << std::endl;
        if (count < l - 1) {
            cacheForwardH  = (cacheForwardH << (2 * k)) + curSeed->second.first;
            cacheReversedH = cacheReversedH + (curSeed->second.first << (2 * k * count));
            begs.push(curSeed->first);
            ++count;
            continue;
        }

        if (count == l - 1) {
            cacheForwardH  = (cacheForwardH << (2 * k)) + curSeed->second.first;
            cacheReversedH = cacheReversedH + (curSeed->second.first << (2 * k * count));
        } else {
            cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + curSeed->second.first;
            cacheReversedH = (cacheReversedH >> (2 * k)) + (curSeed->second.first << (2 * k * (l - 1)));
        }
        begs.push(curSeed->first);
        ++count;
        if (cacheForwardH < cacheReversedH) {
            cacheMin = cacheForwardH;
            rev = false;
        } else if (cacheReversedH < cacheForwardH) {
            cacheMin = cacheReversedH;
            rev = true;
        } else {
            continue; // Skip if strand ambiguous
        }

        newSeedmers.push_back(std::make_tuple(cacheMin, begs.front(), curSeed->second.second, rev));
        begs.pop();
    }
    return newSeedmers;
}


std::vector< std::pair< std::pair< std::pair<size_t, size_t>, size_t >, std::vector<size_t> > > getGroups(const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges, const std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds, const int32_t l) {
    std::vector<std::pair<std::pair<std::pair<size_t, size_t>, size_t>,std::vector<size_t>>> groupsIndex;
    if (syncmerChanges.empty()) return groupsIndex;

    //                 < < <leftUnaffected, offset>, right unaffected >, group_indices >
    std::pair< std::pair< std::pair<size_t, size_t>, size_t >, std::vector<size_t> > curGroup;

    decltype(curSeeds.begin()) prevlb;
    decltype(curSeeds.begin()) curlb;

    if (syncmerChanges.size() == 1) {
        curlb = curSeeds.lower_bound(std::get<1>(syncmerChanges[0]));
        curGroup.second.push_back(0);

        if (curlb == curSeeds.begin()) {
            curGroup.first.first.first = std::numeric_limits<size_t>::max();
            if (curlb->first == std::get<1>(syncmerChanges[0])) ++curlb;
            curGroup.first.second = curlb->first;
        } else {
            auto leftUnaffected = curlb;
            size_t offset = 0;
            for (int i = 0; i < l && leftUnaffected != curSeeds.begin(); ++i) {
                --leftUnaffected;
                ++offset;
            }
            curGroup.first.first.first = leftUnaffected->first;
            curGroup.first.first.second = offset;
            
            auto rightUnaffected = curlb;
            if (curlb == curSeeds.end()) {
                curGroup.first.second = std::numeric_limits<size_t>::max();
            } else  {
                if (curlb->first == std::get<1>(syncmerChanges[0])) ++rightUnaffected;
                if (rightUnaffected == curSeeds.end()) curGroup.first.second = std::numeric_limits<size_t>::max();
                else curGroup.first.second = rightUnaffected->first;
            }
        }

        groupsIndex.push_back(curGroup);
        return groupsIndex;
    }

    for (size_t i = 0; i < syncmerChanges.size(); ++i) {
        const auto& syncmer = syncmerChanges[i];
        if (curGroup.second.empty()) {
            curGroup.second.push_back(i);
            prevlb = curSeeds.lower_bound(std::get<1>(syncmer));
            auto leftUnaffected = prevlb;
            if (leftUnaffected == curSeeds.begin()) {
                curGroup.first.first.first = std::numeric_limits<size_t>::max();
            } else {
                size_t offset = 0;
                for (int i = 0; i < l && leftUnaffected != curSeeds.begin(); ++i) {
                    --leftUnaffected;
                    ++offset;
                }
                curGroup.first.first.first = leftUnaffected->first;
                curGroup.first.first.second = offset;
            }
            continue;
        }

        curlb = curSeeds.lower_bound(std::get<1>(syncmer));
        if (prevlb == curlb) {
            curGroup.second.push_back(i);
        } else {
            size_t numSteps = 0;
            auto prevlbtmp = prevlb;
            while (prevlbtmp != curlb) {
                if (numSteps > l+1) break;
                ++prevlbtmp;
                ++numSteps;
            }
            if (numSteps <= l+1) {
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
                curGroup.first.first.first = leftUnaffected->first;
            }
        }
        prevlb = curlb;
    }

    auto rightUnaffected = curlb;
    if (rightUnaffected == curSeeds.end()) {
        curGroup.first.second = std::numeric_limits<size_t>::max();
    } else {
        if (curlb->first == std::get<1>(syncmerChanges[curGroup.second.back()])) ++rightUnaffected;
        curGroup.first.second = rightUnaffected->first;
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
            {std::get<1>(syncmers[0])},   // start position
            std::get<2>(syncmers[l-1]), // end position
            1,                          // num of identical hash
            false,                      // reversed
            {std::numeric_limits<size_t>::max()},                         // ptr to previous
            std::make_pair(0, false)                     // ptr to next
        };
        seedmersIndex.firstSeedmer = std::make_pair(cacheForwardH, true);
        seedmersIndex.positionMap[std::get<1>(syncmers[0])] = std::make_pair(std::get<2>(syncmers[l-1]), cacheForwardH);
    } else if (cacheReversedH < cacheForwardH) {
        seedmersIndex.seedmerMap[cacheReversedH] = {
            cacheReversedH,
            {std::get<1>(syncmers[0])},
            std::get<2>(syncmers[l-1]),
            1,
            true,
            {std::numeric_limits<size_t>::max()},
            std::make_pair(0, false)
        };
        seedmersIndex.firstSeedmer = std::make_pair(cacheReversedH, true);
        seedmersIndex.positionMap[std::get<1>(syncmers[0])] = std::make_pair(std::get<2>(syncmers[l-1]), cacheReversedH);
    }
    
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    std::pair<size_t, bool> previousSeedmerHash = std::make_pair(seedmersIndex.firstSeedmer.first, true);
    std::pair<size_t, bool> emptySeedmer = {0, false};
    for (int i = 1; i < syncmers.size() - l + 1; ++i) {
        assert(std::get<3>(syncmers[i]) == false);

        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<0>(syncmers[i]),std::get<2>(syncmers[i]));

        cacheForwardH  = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);
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
            resolveSeedmersIndexConflict(seedmersIndex, curHash, previousSeedmerHash, emptySeedmer);
            seedmersIndex.seedmerMap[curHash].begs.insert(std::get<1>(syncmers[i]));
        } else {
            seedmersIndex.seedmerMap[curHash] = {
                curHash,
                {std::get<1>(syncmers[i])},
                std::get<2>(syncmers[i+l-1]),
                1,
                rev,
                {},
                std::make_pair(0, false)
            };

            if (seedmersIndex.firstSeedmer.second != false) {
                assert(previousSeedmerHash.second != false && seedmersIndex.seedmerMap[previousSeedmerHash.first].num == 1);
                seedmersIndex.seedmerMap[previousSeedmerHash.first].next = {curHash, true};
            } else {
                seedmersIndex.firstSeedmer = std::make_pair(curHash, true);
                seedmersIndex.seedmerMap[curHash].prev.insert(std::numeric_limits<size_t>::max());
            }
            if (previousSeedmerHash.second != false) seedmersIndex.seedmerMap[curHash].prev.insert(previousSeedmerHash.first);
            previousSeedmerHash = std::make_pair(curHash, true);
        }

        seedmersIndex.positionMap[std::get<1>(syncmers[i])] = std::make_pair(std::get<2>(syncmers[i+l-1]), curHash);
    }
    seedmersIndex.lastSeedmer = previousSeedmerHash;
    /*rest of syncmers*/
    for (int i = syncmers.size() - l + 1; i < syncmers.size(); i++) {
        assert(std::get<3>(syncmers[i]) == false);
        curSeeds[std::get<1>(syncmers[i])] = std::make_pair(std::get<0>(syncmers[i]),std::get<2>(syncmers[i]));
    }

}

void updateCurSeeds(const std::pair<std::pair<std::pair<size_t, size_t>, size_t>, std::vector<size_t>>& group, const std::vector<std::tuple<size_t, int32_t, int32_t, bool>>& syncmerChanges, std::map<int32_t, std::pair<size_t, int32_t>>& curSeeds) {
    for (const auto& idx : group.second) {
        const std::tuple<size_t, int32_t, int32_t, bool>& change = syncmerChanges[idx];
        auto [seedHash, seedBeg, seedEnd, seedDel] = change;
        if (seedDel) {
            assert(curSeeds.erase(seedBeg));
        } else {
            curSeeds[seedBeg] = std::make_pair(seedHash, seedEnd);
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
        for (int32_t c = globalCoord + len; c >= std::max(0, data.regap[std::max(0, data.degap[globalCoord] - static_cast<int32_t>(k * l))]); c--) {
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
        std::cout << node->identifier << " All" << std::endl;
        for (const auto& seedmer : seedmersIndex.seedmerMap) {
            std::cout << seedmer.second.hash << "\t";
            for (auto b : seedmer.second.begs) {
                std::cout << b << ",";
            }
            std::cout << "\t"
                      << seedmer.second.end << "\t"
                      << seedmer.second.num << "\t"
                      << seedmer.second.rev << "\t";
            for (auto p : seedmer.second.prev) {
                std::cout << p << ",";
            }
            std::cout << "\t" << seedmer.second.next.first << ":" << seedmer.second.next.second <<std::endl;  
        }
        std::cout << "\n---------\n" << std::endl;

        auto curSeedmerHash = seedmersIndex.firstSeedmer;
        auto curSeedmerIt   = seedmersIndex.seedmerMap.find(curSeedmerHash.first);
        std::cout << node->identifier << " Valid only" << std::endl;
        while (true) {
            assert(curSeedmerIt->second.begs.size() == 1);
            std::string s = std::to_string(curSeedmerIt->second.hash) + ","
                          + std::to_string(*(curSeedmerIt->second.begs.begin())) + ","
                          + std::to_string(curSeedmerIt->second.end) + ","
                          + std::to_string(curSeedmerIt->second.rev);
            seedmersOutStream << s << " ";
            

            std::cout << curSeedmerIt->second.hash << "\t"
                      << *(curSeedmerIt->second.begs.begin()) << "\t"
                      << curSeedmerIt->second.end << "\t"
                      << curSeedmerIt->second.num << "\t"
                      << curSeedmerIt->second.rev << "\t";
            if (curSeedmerIt->second.prev.size() > 0) std::cout << *(curSeedmerIt->second.prev.begin()) << "\t";
            else std::cout << "\\\t";
            std::cout << curSeedmerIt->second.next.first << ":" << curSeedmerIt->second.next.second << std::endl;

            if (curSeedmerIt->second.next.second == false) break;
            curSeedmerIt = seedmersIndex.seedmerMap.find(curSeedmerIt->second.next.first);
        }
        std::cout << "---------" << std::endl;
        seedmersOutStream << "\n";
    } else {
        // if (node->identifier == "node_2") {
        std::cerr << node->identifier << std::endl;
        if (true) {
            /*update seedmer and curSeeds*/
            std::cout << node->identifier << std::endl;
            std::sort(syncmerChanges.begin(), syncmerChanges.end(), [](const auto &a, const auto &b) {
                return std::get<1>(a) < std::get<1>(b);
            });

            /*
            group changes in syncmers
                changes within each group has syncmer distance < l-1
                each group separated by l-1 or more syncmers
            */

            
            std::cout << "Total syncmer changes " << syncmerChanges.size() << std::endl;
            //                          < < <leftunaffeected, offset>, rightunaffected > group_indices>
            std::vector<std::pair<std::pair<std::pair<size_t, size_t>, size_t>, std::vector<size_t>>> groupsIndex = getGroups(syncmerChanges, curSeeds, l);

            std::cout << node->identifier << " groups" << std::endl;
            for (const auto& group : groupsIndex) {
                std::cout << group.first.first.first << ":" << group.first.second << "\t";
                for (const auto& idx : group.second) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;
            }
            seedmersOutStream << node->identifier << " ";
            for (const auto& group : groupsIndex) {
                /*
                apply changes to seedmerIndex
                    make new local seedmer chain in group
                        resolve conflicts if ensued
                */
                std::cout << "processing group: ";
                for (const auto& idx : group.second) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;
                std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> locators = {{0, false}, {0, false}};
                std::cout << "start building localSeeds" << std::endl;
                std::map<int32_t, std::pair<size_t, int32_t>> localSeeds = getLocalSeeds(group, syncmerChanges, curSeeds);
                for (const auto& localSeed : localSeeds) {
                    std::cout << localSeed.first << "\t" << localSeed.second.second << "\t" << localSeed.second.first << std::endl;
                }
                std::cout << "start building newSeedmers" << std::endl;
                std::vector<std::tuple<size_t, int32_t, int32_t, bool>> newSeedmers = buildNewSeedmers(seedmersIndex, group, localSeeds, curSeeds, locators, l, k);
                for (const auto& newSeedmer : newSeedmers) {
                    std::cout << std::get<0>(newSeedmer) << "\t"
                                << std::get<1>(newSeedmer) << "\t"
                                << std::get<2>(newSeedmer) << "\t"
                                << std::get<3>(newSeedmer) << std::endl;
                }
                std::cout << "\n-----" << std::endl;

                // std::cout << "ALL:" << std::endl;
                // for (const auto& seedmer : localSeedmersIndex.seedmerMap) {
                //     std::cout << seedmer.second.hash << "\t"
                //             << seedmer.second.beg << "\t"
                //             << seedmer.second.end << "\t"
                //             << seedmer.second.num << "\t"
                //             << seedmer.second.rev << "\t";
                //     for (auto p : seedmer.second.prev) {
                //         std::cout << p << ",";
                //     }
                //     std::cout << "\t" << seedmer.second.next << std::endl;
                // }

                // std::cout << "\n--------\n" << std::endl;

                // std::cout << "VALID:" << std::endl;
                // std::cout << "left locator: " << locators.first.first << "\t" << locators.first.second << std::endl;
                // auto curSeedmer = localSeedmersIndex.firstSeedmer;
                // while (curSeedmer != nullptr) {
                //     assert(curSeedmer->num == 1);
                //     std::cout << curSeedmer->hash << "\t"
                //             << curSeedmer->beg << "\t"
                //             << curSeedmer->end << "\t"
                //             << curSeedmer->num << "\t"
                //             << curSeedmer->rev << "\t";
                //     if (curSeedmer->prev.size() > 0) std::cout << *(curSeedmer->prev.begin()) << "\t";
                //     else std::cout << "\\\t";
                //     std::cout << curSeedmer->next << std::endl;
                //     curSeedmer = curSeedmer->next;
                // }
                // std::cout << "right locator: " << locators.second.first << "\t" << locators.second.second << std::endl;
                // std::cout << "\n--------\n" << std::endl;
                std::cout << "start merging newSeedmers into seedmerIndex" << std::endl;
                mergeSeedmerIndex(newSeedmers, locators, seedmersIndex, seedmersOutStream, curSeeds, l, k);
                std::cout << "after merging, new start seedmer is: " << seedmersIndex.firstSeedmer.first << "\t" << *(seedmersIndex.seedmerMap[seedmersIndex.firstSeedmer.first].begs.begin()) << std::endl;
                /*apply changes to curSeeds*/
                updateCurSeeds(group, syncmerChanges, curSeeds);
            }

            std::cout << node->identifier << " All" << std::endl;
            for (const auto& seedmer : seedmersIndex.seedmerMap) {
                std::cout << seedmer.second.hash << "\t";
                for (auto b : seedmer.second.begs) {
                    std::cout << b << ",";
                }
                std::cout << "\t"
                        << seedmer.second.end << "\t"
                        << seedmer.second.num << "\t"
                        << seedmer.second.rev << "\t";
                for (auto p : seedmer.second.prev) {
                    std::cout << p << ",";
                }
                std::cout << "\t" << seedmer.second.next.first << ":" << seedmer.second.next.second <<std::endl;  
            }
            std::cout << "\n---------\n" << std::endl;

            auto curSeedmerHash = seedmersIndex.firstSeedmer;
            auto curSeedmerIt   = seedmersIndex.seedmerMap.find(curSeedmerHash.first);
            std::cout << node->identifier << " Valid only" << std::endl;
            while (true) {
                std::cout << curSeedmerIt->second.hash << "\t"
                        << *(curSeedmerIt->second.begs.begin()) << "\t"
                        << curSeedmerIt->second.end << "\t"
                        << curSeedmerIt->second.num << "\t"
                        << curSeedmerIt->second.rev << "\t";
                assert(curSeedmerIt->second.begs.size() == 1 && curSeedmerIt->second.num == 1 && curSeedmerIt->second.prev.size() <= 1);
                if (curSeedmerIt->second.prev.size() > 0) std::cout << *(curSeedmerIt->second.prev.begin()) << "\t";
                else std::cout << "\\\t";
                std::cout << curSeedmerIt->second.next.first << ":" << curSeedmerIt->second.next.second << std::endl;

                if (curSeedmerIt->second.next.second == false) break;
                curSeedmerIt = seedmersIndex.seedmerMap.find(curSeedmerIt->second.next.first);
            }
            std::cout << "---------" << std::endl;
            seedmersOutStream << "\n";
        }

        std::cout << node->identifier << " Seeds" << std::endl;
        for (const auto& curSeed : curSeeds) {
            std::cout << curSeed.first << "\t"
                      << curSeed.second.second << "\t"
                      << curSeed.second.first << std::endl;
        }
        std::cout << "---------" << std::endl;
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
    std::map<int32_t, std::pair<size_t, int32_t>> curSeeds; // map[beg] = {hash, end}
    index.outStream << k << " " << s << " " << l << "\n";
    seedmersOutStream << k << " " << s << " " << l << "\n";
    /* Recursive traversal of tree to build the index */
    mgsr::buildSeedmerHelper(data, seedMap, index, seedmersIndex, curSeeds, T, T->root, l, k, s, globalCoords, seedmersOutStream);
}
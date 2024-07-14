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
#include <tbb/global_control.h>
#include "PangenomeMAT.hpp"
#include "seeding.hpp"
#include "mgsr.hpp"
#include <eigen3/Eigen/Dense>

namespace fs = boost::filesystem;
using namespace PangenomeMAT;

void mgsr::purpose() {
    std::cerr << "butter" << std::endl;
}

std::vector<std::tuple<size_t, int32_t, int32_t>> syncmersSketch(const std::string& seq, const int k, const int s, const bool open) {
    std::vector<std::tuple<size_t, int32_t, int32_t>> syncmers;
    for (size_t i = 0; i < seq.size() - k + 1; ++i) {
        std::string kmer = seq.substr(i, k);
        if (!seeding::is_syncmer(kmer, s, open)) continue;

        std::pair<size_t, bool> kmerHash = seeding::getHash(kmer);
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
            seedmers.positionMap.erase(beg);
            assert(seedmersMapIt->second.erase(beg));
            seedmersMapIt->second.erase(beg);
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
                seedmersMapIt->second.erase(beg);
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
            seedmers.positionMap.erase(beg);
            assert(seedmersMapIt->second.erase(beg));
            seedmersMapIt->second.erase(beg);
            if (seedmersMapIt->second.empty()) {
                seedmers.seedmersMap.erase(hash);
            }
        } else if (type == 0) {
            auto positionMapIt = seedmers.positionMap.find(beg);
            if (positionMapIt != seedmers.positionMap.end()) {
                auto [oend, ohash, orev] = positionMapIt->second;
                auto seedmersMapIt = seedmers.seedmersMap.find(ohash);

                assert(seedmersMapIt->second.erase(beg));
                seedmersMapIt->second.erase(beg);
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


void scoreDFS(
    mgsr::seedmers& seedmers, const std::unordered_map<std::string, std::pair<int32_t, std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>>>>& seedmersIndex,
    const std::unordered_map<std::string, std::map<int32_t, int32_t>>& coordsIndex,
    const std::vector<std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>>>& readSeedmers,
    const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::vector<bool>& lowScoreReads,
    std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores,
    std::unordered_map<std::string, std::string>& identicalPairs, const Node *node, Tree *T, std::atomic<size_t>& numLowScoreReads,
    const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate
    ) {
    // std::cerr << "identifier " << node->identifier << std::endl;
    std::vector<std::tuple<int32_t, int32_t, size_t, bool, int16_t>> seedmersToRevert;
    std::unordered_set<size_t> affectedSeedmers;
    mutateSeedmers(seedmers, seedmersIndex.find(node->identifier)->second, coordsIndex.find(node->identifier)->second, affectedSeedmers, seedmersToRevert);    
    size_t num_cpus = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    if (node->identifier == T->root->identifier) {
        allScores[node->identifier].resize(readSeedmers.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, readSeedmers.size(), readSeedmers.size() / num_cpus),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const auto& curReadSeedmers = readSeedmers[i];
                    double duplicates = 0;
                    std::vector<match_t> matches = match(curReadSeedmers.first, seedmers, duplicates);
                    std::vector<match_t> pseudoChain = chainPseudo(matches, coordsIndex.find(node->identifier)->second, maximumGap, minimumCount, minimumScore);
                    int32_t pseudoScore = scorePseudoChain(pseudoChain);
                    double maxScore = static_cast<double>(curReadSeedmers.first.size()) - duplicates;
                    // expected number seedmer per read is (150 - k + 1) / ((k - s + 1) / 2) - l + 1
                    // Need to change how pseudoProb is calculated
                    // pseudoProb = 0 if pseudoScore < 1/2 (estimated number of kminmers on a read)
                    double pseudoProb;
                    if (pseudoScore < 45 / 2) {
                        pseudoProb = std::numeric_limits<double>::min();
                    } else {
                        pseudoProb = pow(errorRate, maxScore - static_cast<double>(pseudoScore)) * pow(1 - errorRate, static_cast<double>(pseudoScore));
                        if (lowScoreReads[i]) {
                            lowScoreReads[i] = false;
                            --numLowScoreReads;
                        }
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
            tbb::parallel_for(tbb::blocked_range<size_t>(0, readSeedmers.size(), readSeedmers.size() / num_cpus),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto& curReadSeedmers = readSeedmers[i];
                        if (redo(curReadSeedmers.second, affectedSeedmers)) {
                            double duplicates = 0;
                            std::vector<match_t> matches = match(curReadSeedmers.first, seedmers, duplicates);
                            std::vector<match_t> pseudoChain = chainPseudo(matches, coordsIndex.find(node->identifier)->second, maximumGap, minimumCount, minimumScore);
                            int32_t pseudoScore = scorePseudoChain(pseudoChain);
                            double maxScore = static_cast<double>(curReadSeedmers.first.size()) - duplicates;
                            double pseudoProb;
                            if (pseudoScore < 45 / 2) {
                                pseudoProb = std::numeric_limits<double>::min();
                            } else {
                                pseudoProb = pow(errorRate, maxScore - static_cast<double>(pseudoScore)) * pow(1 - errorRate, static_cast<double>(pseudoScore));
                                if (lowScoreReads[i]) {
                                    lowScoreReads[i] = false;
                                    --numLowScoreReads;
                                }
                            }
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
        scoreDFS(seedmers, seedmersIndex, coordsIndex, readSeedmers, numReadDuplicates, lowScoreReads, allScores, identicalPairs, child, T, numLowScoreReads, maximumGap, minimumCount, minimumScore, errorRate);
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
    std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates, std::vector<bool>& lowScoreReads,
    std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestor,
    std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets, Tree *T,
    std::vector<std::vector<seeding::seed>>& readSeeds, std::vector<std::string>& readSequences, std::vector<std::string>& readQuals,
    std::vector<std::string>& readNames, std::atomic<size_t>& numLowScoreReads, const int& maximumGap, const int& minimumCount, const int& minimumScore, const double& errorRate
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
                if (metaSplit.size() < 2) continue;
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
    std::cerr << "readNames.size(): " << readNames.size() << std::endl;
    lowScoreReads.resize(readSeedmers.size(), true);
    std::cerr << "finished initializing read seedmers... total number of reads " << readSequences.size() << "\n" << std::endl;

    std::cerr << "start scoring DFS" << std::endl;
    mgsr::seedmers seedmers;
    size_t totalCount = 0;
    size_t redoCount = 0;
    numLowScoreReads = readSeedmers.size();
    scoreDFS(seedmers, seedmersIndex, coordsIndex, readSeedmers, numReadDuplicates, lowScoreReads, allScores, identicalPairs, T->root, T, numLowScoreReads, maximumGap, minimumCount, minimumScore, errorRate);
    std::cerr << "finished scoring DFS\n" << std::endl;



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

    std::vector<std::pair<std::string, int32_t>> scores;
    scores.reserve(allScores.size() - leastRecentIdenticalAncestor.size());
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
    const std::vector<bool>& lowScoreReads, const int32_t& numReads, const size_t& numLowScoreReads, const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    Eigen::MatrixXd& probs, std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const int32_t& roundsRemove, const double& removeThreshold, std::string exclude
    ) {
    if (exclude.empty()) {
        std::stringstream msg;
        msg << "starting to set up EM" << "\n";
        std::cerr << msg.str();
    } else {
        std::stringstream msg;
        msg << "starting to set up EM excluding " << exclude << "\n";
        std::cerr << msg.str();
    }


    if (!exclude.empty()) {
        probs.resize(allScores.begin()->second.size() - numLowScoreReads, allScores.size() - leastRecentIdenticalAncestors.size() - 1);
    } else {
        probs.resize(allScores.begin()->second.size() - numLowScoreReads, allScores.size() - leastRecentIdenticalAncestors.size());
    }
    size_t colIndex = 0;
    for (const auto& node : allScores) {
        if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
        if (!exclude.empty() && node.first == exclude) continue;
        std::vector<double> curProbs;
        size_t rowIndex = 0;
        for (size_t i = 0; i < node.second.size(); ++i) {
            if (!lowScoreReads[i]) {
                const auto& score = node.second[i];
                probs(rowIndex, colIndex) = score.second;
                ++rowIndex;
            }
        }
        nodes.push_back(node.first);
        ++colIndex;
    }

    // std::cerr << "num nodes " << nodes.size() << std::endl;
    props = Eigen::VectorXd::Constant(nodes.size(), 1.0 / static_cast<double>(nodes.size()));
    Eigen::VectorXd readDuplicates(allScores.begin()->second.size() - numLowScoreReads);
    
    size_t indexReadDuplicates = 0;
    int32_t numHighScoreReads = 0;
    for (size_t i = 0; i < numReadDuplicates.size(); ++i) {
        if (!lowScoreReads[i]) {
            readDuplicates(indexReadDuplicates) = numReadDuplicates[i].first;
            numHighScoreReads += numReadDuplicates[i].first;
            ++indexReadDuplicates;
        }
    }

    if (exclude.empty()) {
        std::stringstream msg;
        msg << "starting EM estimation of haplotype proportions" << "\n";
        std::cerr << msg.str();
    } else {
        std::stringstream msg;
        msg << "starting EM estimation of haplotype proportions excluding " << exclude << "\n";
        std::cerr << msg.str();
    }

    int curit = 0;
    llh = getExp(probs, props, readDuplicates);
    // std::cerr << "iteration " << curit << ": " << llh << std::endl;
    squarem(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit);

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
        squarem(sigNodes, sigProbs, identicalSets, readDuplicates, numHighScoreReads, sigProps, llh, curit);
        nodes = std::move(sigNodes);
        probs = std::move(sigProbs);
        props = std::move(sigProps);
    }

    if (!exclude.empty()) {
        Eigen::VectorXd curProbs(allScores.begin()->second.size() - numLowScoreReads);
        size_t indexCurProbs = 0;
        const auto& curNode = allScores.at(exclude);
        for (size_t i = 0; i < curNode.size(); ++i) {
            if (!lowScoreReads[i]) {
                const auto& score = curNode[i];
                curProbs(indexCurProbs) = score.second;
                ++indexCurProbs;
            }
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
    const std::vector<bool>& lowScoreReads, std::unordered_map<std::string, std::unordered_set<size_t>>& assignedReads
    ) {
    size_t rowindex = 0;
    for (size_t i = 0; i < numReadDuplicates.size(); ++i) {
        if (lowScoreReads[i]) continue;
        const Eigen::VectorXd& curprobs = probs.row(rowindex);
        ++rowindex;
        double curmax = curprobs.maxCoeff();
        for (size_t j = 0; j < curprobs.size(); ++j) {
            if (curmax == curprobs(j)) {
                for (size_t readIndex : numReadDuplicates[i].second) {
                    assignedReads[nodes[j]].insert(readIndex);
                }
            }
        }
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



// squaremHelper_test_1: remove haplotype whose abundance is estimated to be at std::numeric_limits<double>::min() for more than m iterations in n iterations
// default to 5 iterations rn



void updateInsigCounts(const Eigen::VectorXd& props, std::vector<size_t>& insigCounts, size_t totalNodes) {
    double insigProp =  (1.0 / static_cast<double>(totalNodes)) / 10;
    for (int i = 0; i < props.size(); ++i) {
        if (props(i) <= insigProp) {
            ++insigCounts[i];
        } else {
            insigCounts[i] = 0;
        }
    }
}

void squarem_test_1(
    const std::vector<std::string>& nodes, const Eigen::MatrixXd& probs,
    const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    const Eigen::VectorXd& numReadDuplicates, const int32_t& numReads, 
    Eigen::VectorXd& props, double& llh, int& curit, bool& converged, size_t iterations, std::vector<size_t>& insigCounts, size_t totalNodes
    ) {
    assert(nodes.size() == probs.cols());
    assert(nodes.size() == props.size());
    size_t curIteration = 1;
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

        updateInsigCounts(props, insigCounts, totalNodes);
        if (newllh - llh < 0.00001) {
            llh = newllh;
            converged = true;
            break;
        } else if (curIteration == iterations) {
            llh = newllh;
            break;
        }
        llh = newllh;
        ++curit;
        ++curIteration;
    }
    ++curit;
}

void mgsr::squaremHelper_test_1(
    PangenomeMAT::Tree *T, const std::unordered_map<std::string, tbb::concurrent_vector<std::pair<int32_t, double>>>& allScores, const std::vector<std::pair<int32_t, std::vector<size_t>>>& numReadDuplicates,
    const std::vector<bool>& lowScoreReads, const int32_t& numReads, const size_t& numLowScoreReads, const std::unordered_map<std::string, std::string>& leastRecentIdenticalAncestors, const std::unordered_map<std::string, std::unordered_set<std::string>>& identicalSets,
    Eigen::MatrixXd& probs, std::vector<std::string>& nodes, Eigen::VectorXd& props, double& llh, const int32_t& roundsRemove, const double& removeThreshold, std::string exclude
    ) {
    if (exclude.empty()) {
        std::stringstream msg;
        msg << "starting to set up EM" << "\n";
        std::cerr << msg.str();
    } else {
        std::stringstream msg;
        msg << "starting to set up EM excluding " << exclude << "\n";
        std::cerr << msg.str();
    }


    if (!exclude.empty()) {
        probs.resize(allScores.begin()->second.size() - numLowScoreReads, allScores.size() - leastRecentIdenticalAncestors.size() - 1);
    } else {
        probs.resize(allScores.begin()->second.size() - numLowScoreReads, allScores.size() - leastRecentIdenticalAncestors.size());
    }
    size_t colIndex = 0;
    for (const auto& node : allScores) {
        if (leastRecentIdenticalAncestors.find(node.first) != leastRecentIdenticalAncestors.end()) continue;
        if (!exclude.empty() && node.first == exclude) continue;
        std::vector<double> curProbs;
        size_t rowIndex = 0;
        for (size_t i = 0; i < node.second.size(); ++i) {
            if (!lowScoreReads[i]) {
                const auto& score = node.second[i];
                probs(rowIndex, colIndex) = score.second;
                ++rowIndex;
            }
        }
        nodes.push_back(node.first);
        ++colIndex;
    }

    // std::cerr << "num nodes " << nodes.size() << std::endl;
    props = Eigen::VectorXd::Constant(nodes.size(), 1.0 / static_cast<double>(nodes.size()));
    size_t totalNodes = nodes.size();
    Eigen::VectorXd readDuplicates(allScores.begin()->second.size() - numLowScoreReads);
    
    size_t indexReadDuplicates = 0;
    int32_t numHighScoreReads = 0;
    for (size_t i = 0; i < numReadDuplicates.size(); ++i) {
        if (!lowScoreReads[i]) {
            readDuplicates(indexReadDuplicates) = numReadDuplicates[i].first;
            numHighScoreReads += numReadDuplicates[i].first;
            ++indexReadDuplicates;
        }
    }

    if (exclude.empty()) {
        std::stringstream msg;
        msg << "starting EM estimation of haplotype proportions" << "\n";
        std::cerr << msg.str();
    } else {
        std::stringstream msg;
        msg << "starting EM estimation of haplotype proportions excluding " << exclude << "\n";
        std::cerr << msg.str();
    }

    int curit = 0;
    llh = getExp(probs, props, readDuplicates);
    // std::cerr << "iteration " << curit << ": " << llh << std::endl;
    // bool& converged, size_t iterations, std::vector<size_t>& insigCounts
    bool converged = false;
    size_t iterations = 20;
    size_t remove_count = 20;
    for (size_t i = 0; i < 3; ++i) {
        std::cerr << "Here" << std::endl;
        std::vector<size_t> insigCounts(nodes.size());
        squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, iterations, insigCounts, totalNodes);
        if (converged) {
            break;
        }
        
        std::vector<size_t> significantIndices;
        std::vector<std::string> sigNodes;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (insigCounts[i] < remove_count) {
                significantIndices.push_back(i);
            }
        }

        if (significantIndices.size() == nodes.size()) continue;
        for (size_t idx : significantIndices) {
            sigNodes.push_back(nodes[idx]);
        }

        Eigen::MatrixXd sigProbs(probs.rows(), significantIndices.size());
        for (size_t i = 0; i < significantIndices.size(); ++i) {
            sigProbs.col(i) = probs.col(significantIndices[i]);
        }
        Eigen::VectorXd sigProps = Eigen::VectorXd::Constant(sigNodes.size(), 1.0 / static_cast<double>(sigNodes.size()));
        std::cerr << "dropped " << nodes.size() - sigNodes.size() << " during EM" << std::endl;
        nodes = std::move(sigNodes);
        probs = std::move(sigProbs);
        props = std::move(sigProps);
    }

    if (!converged) {
        std::vector<size_t> insigCounts(nodes.size());
        squarem_test_1(nodes, probs, identicalSets, readDuplicates, numHighScoreReads, props, llh, curit, converged, std::numeric_limits<size_t>::max(), insigCounts, totalNodes);
        assert(converged == true);
    }

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
        bool converged = false;
        size_t iterations = std::numeric_limits<size_t>::max();
        std::vector<size_t> insigCounts(sigNodes.size());
        squarem_test_1(sigNodes, sigProbs, identicalSets, readDuplicates, numHighScoreReads, sigProps, llh, curit, converged, iterations, insigCounts, totalNodes);
        assert(converged);
        nodes = std::move(sigNodes);
        probs = std::move(sigProbs);
        props = std::move(sigProps);
    }

    if (!exclude.empty()) {
        Eigen::VectorXd curProbs(allScores.begin()->second.size() - numLowScoreReads);
        size_t indexCurProbs = 0;
        const auto& curNode = allScores.at(exclude);
        for (size_t i = 0; i < curNode.size(); ++i) {
            if (!lowScoreReads[i]) {
                const auto& score = curNode[i];
                curProbs(indexCurProbs) = score.second;
                ++indexCurProbs;
            }
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

/*
normalize:
    if prop = 0:
        low_sig_freq += 1
        ...
    else:
        if low_sig_freq > 0
        low_sig_freq = 0
        ...
    ...

squarem:
    ...
    if curit == 50:
        return;
    else if meet  convergence criteria:
        converged = true
        return;
    ...

squaremHelper:
    bool converged = false
    int n = 50
    while not converged:
        low_sig_freq = []  //keep track of how many times a node has been at low freq
        squarem(params, converged, n, low_sig_freq)
        remove low_sig_freq nodes from props and probs

*/
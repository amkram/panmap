#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "util.hpp"
#include <cmath>

using namespace PangenomeMAT;
using namespace tree;
using namespace seeding;
using namespace util;

void mutateSeedmerMap(std::unordered_map<int32_t, std::string> &seedmers, std::string nid, seedmerIndex_t &index) {
    for (const auto &op : index[nid]) {
        std::string seedmer = op.second;
        if (seedmer == "" || seedmer == "@") {
            continue;
        }
        int32_t pos = op.first;
        if (seedmer[0] == '@') { // deletion
            if (seedmers.find(pos) != seedmers.end()) {
                seedmers.erase(seedmers.find(pos));
            }
        } else {
            seedmers[pos] = seedmer;
        }
    }
}
void revertSeedmerMap(std::unordered_map<int32_t, std::string> &seedmers, std::string nid, seedmerIndex_t &index) {
    for (auto it = index[nid].rbegin(); it != index[nid].rend(); it++) {
        std::string seedmer = it->second;
        int32_t pos = it->first;
        if (seedmer[0] == '@') { // deletion
            seedmers[pos] = seedmer.substr(1);
        } else {
            if (seedmers.find(pos) != seedmers.end()) {
                seedmers.erase(seedmers.find(pos));
            }
        }
    }
}
void getPhyloCounts(seedmerIndex_t &index, std::unordered_map<std::string, int32_t> &counts, std::unordered_map<int32_t, std::string> &seedmers, Tree *T, Node *node) {
    mutateSeedmerMap(seedmers, node->identifier, index);
    std::unordered_map<std::string, bool> seen;
    for (const auto &seed : seedmers) {
        if (seed.second == "") {
            continue;
        }
        if (seen.find(seed.second) == seen.end()) {
            seen[seed.second] = true;
        } else {
            continue;
        }
        counts[seed.second] += 1;
    }
    for (Node *child : node->children) {
        getPhyloCounts(index, counts, seedmers, T, child);
    }
    revertSeedmerMap(seedmers, node->identifier, index);
}

void placeDFS(std::ofstream *out, std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>> &index, std::unordered_map<int32_t, std::string> &seedmers, std::map<std::string, float> &scores, std::unordered_map<std::string, int32_t> &seedmerCounts, std::unordered_map<std::string, int32_t> &phyloCounts, const Node *node, Tree *T, const std::string optionalTarget, std::unordered_map<int32_t, std::string> *optionalOutputSeedmers) {
    mutateSeedmerMap(seedmers, node->identifier, index);

    if (optionalTarget != "" && optionalTarget == node->identifier) {
        *optionalOutputSeedmers = seedmers;
        return;
    }
    float score = 0;
    std::unordered_map<std::string, int32_t> nodeSeedmerCounts;
    for (const auto &seed : seedmers) {
        if (nodeSeedmerCounts.find(seed.second) == nodeSeedmerCounts.end()) {
            nodeSeedmerCounts[seed.second] = 1;
        } else {
            nodeSeedmerCounts[seed.second] += 1;
        }
    }
    for (const auto &seed : seedmers) {
        if (seedmerCounts.find(seed.second) != seedmerCounts.end()) {
            int32_t den = phyloCounts[seed.second] * nodeSeedmerCounts[seed.second];
            score += std::min(5.0, (double) seedmerCounts[seed.second]) / (double) den;
        }
    }

    scores[node->identifier] = score;
    if (out != nullptr) {
        *out << node->identifier << " " << score << "\n";
    }
    for (Node *child : node->children) {
        placeDFS(out, index, seedmers, scores, seedmerCounts, phyloCounts, child, T, optionalTarget, optionalOutputSeedmers);
    }
    revertSeedmerMap(seedmers, node->identifier, index);
}

std::string reverseComplement(std::string dna_sequence) {
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
void placeHelper(std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(seed_cmp)>> &index, std::unordered_map<int32_t, std::string> &seedmers, std::map<std::string, float> &scores, std::unordered_map<std::string, int32_t> &seedmerCounts, std::unordered_map<std::string, int32_t> &phyloCounts, const Node *node, Tree *T, const std::string optionalTarget) {
    util::scopedTimer();
    std::ofstream out("placement.out");

    placeDFS(&out, index, seedmers, scores, seedmerCounts, phyloCounts, T->root, T, "", nullptr);
}

seedmerIndex_t seedsFromFastq(std::ifstream &indexFile, int32_t *k, int32_t *s, int32_t *j, std::unordered_map<std::string, int32_t> &seedmerCounts, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, const std::string &fastqPath) {
    util::scopedTimer();
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    seq = kseq_init(fileno(fp));
    seedmerIndex_t seedmerIndex;
    std::string line0;
    std::getline(indexFile, line0);
    std::vector<std::string> spltTop;
    stringSplit(line0, ' ', spltTop);
    int32_t tempK = std::stoi(spltTop[0]);
    int32_t tempS = std::stoi(spltTop[1]);
    int32_t tempJ = std::stoi(spltTop[2]);
    
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.s);
    }
    for (int i = 0; i < readSequences.size(); i++) {        
        std::string seq = readSequences[i];
        std::string name = readNames[i];
        std::string rc = reverseComplement(seq);
        std::vector<seed> syncmers = syncmerize(seq, tempK, tempS, false, true, 0);
        std::vector<seed> syncmersReverse = syncmerize(rc, tempK, tempS, false, true, 0);
        std::vector<seedmer> seedmers = seedmerize(syncmers, tempJ);
        std::vector<seedmer> seedmersReverse = seedmerize(syncmersReverse, tempJ);  
        std::unordered_map<std::string, bool> seen;
        for (const auto &m : seedmers) {
            if (seen.find(m.seq) == seen.end()) {
                seen[m.seq] = true;
            } else {
                continue;
            }
            if (seedmerCounts.find(m.seq) == seedmerCounts.end()) {
                seedmerCounts[m.seq] = 1;
            } else {
                seedmerCounts[m.seq] += 1;
            }
        }
        for (const auto &m : seedmersReverse) {
            if (seen.find(m.seq) == seen.end()) {
                seen[m.seq] = true;
            } else {
                continue;
            }
            if (seedmerCounts.find(m.seq) == seedmerCounts.end()) {
                seedmerCounts[m.seq] = 1;
            } else {
                seedmerCounts[m.seq] += 1;
            }
        }
    }

    while (std::getline(indexFile, line0)) {
        std::vector<std::string> splt;
        stringSplit(line0, ' ', splt);
        std::string nid = splt[0];
        seedmerIndex[nid] = {};
        for (int32_t i = 1; i < splt.size(); i++) {
            std::vector<std::string> splt2;
            stringSplit(splt[i], ':', splt2);
            int32_t pos = std::stoi(splt2[0]);
            std::string seedmer = splt2[1];
            seedmerIndex[nid].insert(std::make_pair(pos, seedmer));
        }
    }
    *k = tempK;
    *s = tempS;
    *j = tempJ;
    return seedmerIndex;
}

void place::placeIsolate(std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path, Tree *T) {
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);

    /* Read processing */
    std::cout << "\n◠◡ Processing reads ... " << std::flush;
    std::vector<std::string> readSequences;
    std::vector<std::string> readQuals;
    std::vector<std::string> readNames;
    std::unordered_map<std::string, int32_t> seedmerCounts;
    int32_t k, s, j;

    seedmerIndex_t seedmerIndex = seedsFromFastq(indexFile, &k, &s, &j, seedmerCounts, readSequences, readQuals, readNames, reads1Path);

    /* Sample placement */
    std::cout << "⋌⋋ Placing sample ... " << std::flush;
    std::map<std::string, float> scores;
    std::unordered_map<std::string, int32_t> phyloCounts;
    std::unordered_map<int32_t, std::string> dynamicSeedmersPhylo;
    getPhyloCounts(seedmerIndex, phyloCounts, dynamicSeedmersPhylo, T, T->root);

    std::unordered_map<int32_t, std::string> dynamicSeedmersPlace;
    placeHelper(seedmerIndex, dynamicSeedmersPlace, scores, seedmerCounts, phyloCounts, T->root, T, "");

    /* Setup target sequence and seeds */
    std::vector<std::pair<std::string, float>> targetNodes;
    std::copy(scores.begin(), scores.end(), back_inserter<std::vector<std::pair<std::string, float>>>(targetNodes));
    std::sort(targetNodes.begin(), targetNodes.end(), [] (auto &left, auto &right) { return left.second > right.second; });
    // redo DFS with target node -> returns early with dynamicSeedmers in target node's state
    std::string bestMatch = targetNodes[0].first;
    std::string bestMatchSequence = "";
    std::string gappedSeq = T->getStringFromReference(bestMatch, true);
    std::vector<int32_t> degap;
    for (int32_t i = 0; i < gappedSeq.size(); i ++) {
        char &c = gappedSeq[i];
        degap.push_back(bestMatchSequence.size());
        if (c != '-') {
            bestMatchSequence += c;
        }
    }
    std::unordered_map<std::string, std::vector<int32_t>> seedToRefPositions;
    std::unordered_map<int32_t, std::string> targetSeedmers;
    std::unordered_map<int32_t, std::string> dynamicSeedmersTarget;
    placeDFS(nullptr, seedmerIndex, dynamicSeedmersTarget, scores, seedmerCounts, phyloCounts, T->root, T, bestMatch, &targetSeedmers);
    for (const auto &seedmer : targetSeedmers) {
        if (seedmer.second == "" ) {
            continue;
        }
        int32_t refPos = seedmer.first;
        std::string seed = seedmer.second.substr(0, k);
        if (seedToRefPositions.find(seed) == seedToRefPositions.end()) {
            seedToRefPositions[seed] = {};
        }
        seedToRefPositions[seed].push_back(degap[refPos]);
    }
    std::cout << "\n\nBest match: " << bestMatch << " (" << targetNodes[0].second << ")\n";

    /* Alignment to target */

    /*  @nico: at this point,
    **    - bestMatch should contain the target node id
    **    - bestMatchSequence has the target node's sequence without gaps
    **    - readSequences, readQuals, readNames contain the relevant info
    **    - seedToRefPositions maps a seed sequence (just seeds not seed-mers) to all matching positions in bestMatchSequence
    */
}

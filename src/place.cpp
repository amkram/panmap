#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "tree.hpp"
#include "util.hpp"

using namespace PangenomeMAT;
using namespace tree;
using namespace seed;
using namespace util;

struct dynamicJaccard {
    size_t intersectionSize;
    size_t unionSize;
    float jaccardIndex;
};

void updateJaccard(dynamicJaccard &dj, const std::vector<std::string> &readSeeds, const std::vector<std::string> &nodeSeeds, const std::vector<kmer_t> &deletedSeeds, const std::vector<kmer_t> &insertedSeeds) {
    std::unordered_map<std::string, int32_t> readSeedCounts;
    std::unordered_map<std::string, int32_t> nodeSeedCounts;
    int32_t max = 0;
    for (const std::string &s : readSeeds) {
        if (readSeedCounts.find(s) == readSeedCounts.end()) {
            readSeedCounts[s] = 1;
        } else {
            readSeedCounts[s] += 1;
        }
        if (readSeedCounts[s] > max) {
            max = readSeedCounts[s];
        }
    }
    for (const std::string &s : nodeSeeds) {
        if (nodeSeedCounts.find(s) == nodeSeedCounts.end()) {
            nodeSeedCounts[s] = 1;
        } else {
            nodeSeedCounts[s] += 1;
        }
        if (readSeedCounts[s] > max) {
            max = readSeedCounts[s];
        }
    }
    for (const kmer_t &syncmer : deletedSeeds) {
        
        if (std::find(readSeeds.begin(), readSeeds.end(), syncmer.seq) != readSeeds.end()) {
            dj.intersectionSize -= readSeedCounts[syncmer.seq] / nodeSeedCounts[syncmer.seq];
            dj.unionSize -= readSeedCounts[syncmer.seq] / nodeSeedCounts[syncmer.seq];
        } else {
            dj.unionSize -= readSeedCounts[syncmer.seq] / nodeSeedCounts[syncmer.seq];
        }
    }
    for (const kmer_t &syncmer : insertedSeeds) {
        if (std::find(readSeeds.begin(), readSeeds.end(), syncmer.seq) != readSeeds.end()) {
            dj.intersectionSize += readSeedCounts[syncmer.seq] / nodeSeedCounts[syncmer.seq];;
            dj.unionSize += readSeedCounts[syncmer.seq] / nodeSeedCounts[syncmer.seq];;
        } else {
            dj.unionSize += readSeedCounts[syncmer.seq] / nodeSeedCounts[syncmer.seq];;
        }
    }
    dj.jaccardIndex = (float) dj.intersectionSize / dj.unionSize;
}
void initializeJaccard(dynamicJaccard &dj, const std::vector<kmer_t> &readSeeds, std::vector<kmer_t> nodeSeeds) {
    dj.intersectionSize = 0;
    dj.unionSize = 0;
    std::vector<std::string> v;
    std::sort(nodeSeeds.begin(), nodeSeeds.end());
    int32_t max = 0;
    std::unordered_map<std::string, int32_t> readSeedCounts;
    std::unordered_map<std::string, int32_t> nodeSeedCounts;

    for (const kmer_t &s : readSeeds) {
        if (readSeedCounts.find(s.seq) == readSeedCounts.end()) {
            readSeedCounts[s.seq] = 1;
        } else {
            readSeedCounts[s.seq] += 1;
        }
        if (readSeedCounts[s.seq] > max) {
            max = readSeedCounts[s.seq];
        }
    }
    for (const kmer_t &s : nodeSeeds) {
        if (nodeSeedCounts.find(s.seq) == nodeSeedCounts.end()) {
            nodeSeedCounts[s.seq] = 1;
        } else {
            nodeSeedCounts[s.seq] += 1;
        }
    }

    std::sort(nodeSeeds.begin(), nodeSeeds.end(), [] (auto &left, auto &right) {
        return left.seq < right.seq;
    });
    for (int32_t i = 0; i < nodeSeeds.size(); i++) {
        

    }
    // for (const std::string &s : nodeSeeds) {
    //     dj.intersectionSize += readSeedCounts[s] / nodeSeedCounts[s];
    //     std::cout << s << " " << readSeedCounts[s] << " " << nodeSeedCounts[s] << "\n";
    //     // score increments by (# in reads) / (# in node)
    //     // e.g. if a seed appears once in the node, 4 / 1 = 4
    //     // or appearing 3 times: 4 / 3 = 1 (integer division) TODO: float?
    // }
    dj.unionSize = nodeSeeds.size() * max + readSeeds.size() - dj.intersectionSize;
    dj.jaccardIndex = (float) dj.intersectionSize / dj.unionSize;

    std::cout << "Jaccard init: " << dj.intersectionSize << " / " << dj.unionSize << " = " << dj.jaccardIndex << "\n";
}
void placeDFS(std::vector<kmer_t> &nodeSeeds, dynamicJaccard &dj, std::unordered_map<std::string, float> &scores, seedIndex &index, const std::vector<kmer_t> &readSeeds, const Node *node, Tree *T) {
    std::stack<int32_t> rmDel;
    for (const kmer_t &d : index.deletions[node->identifier]) {
        rmDel.push(d.idx);
    }
    removeIndices(nodeSeeds, rmDel);
    for (const kmer_t &s : index.insertions[node->identifier]) {
        nodeSeeds.push_back(s);
    }
    std::vector<std::string> nodeSeedStrings;
    for (const kmer_t &s : nodeSeeds) {
        nodeSeedStrings.push_back(s.seq);
    }
    std::vector<std::string> readSeedStrings;
    for (const auto &p : readSeeds) {
        readSeedStrings.push_back(p.seq);
    }
    if (node == T->root) {
        initializeJaccard(dj, readSeeds, nodeSeeds);
    }
    updateJaccard(dj, readSeedStrings, nodeSeedStrings, index.deletions[node->identifier], index.insertions[node->identifier]);
    scores[node->identifier] = dj.jaccardIndex;

    for (Node *child : node->children) {
        placeDFS(nodeSeeds, dj, scores, index, readSeeds, child, T);
    }

    nodeSeeds.erase(nodeSeeds.end() - index.insertions[node->identifier].size(), nodeSeeds.end());
    for (int32_t i = index.deletions[node->identifier].size() - 1; i >= 0; i--) {
        nodeSeeds.insert(nodeSeeds.begin() + index.deletions[node->identifier][i].idx, index.deletions[node->identifier][i]);
    }
}
auto cmp = [](const std::pair<int32_t, std::string> &a, const std::pair<int32_t, std::string> &b) {
    if(a.first != b.first) {
        return a.first < b.first;
    }
    return a.second < b.second;
};

void getPhyloCounts(std::unordered_map<std::string, int32_t> &counts, std::unordered_map<int32_t, std::string> jkmers, Tree *T, Node *node, std::unordered_map<std::string, std::string> &debugStrings, int32_t debug_k, int32_t debug_s) {
    
    std::vector<kmer_t> debug_syncmers = syncmerize(debugStrings[node->identifier], debug_k, debug_s, false, true, 0);
    std::vector<jkmer> debug_jkmers = jkmerize(debug_syncmers, 3);

    int32_t max = 0;
    std::unordered_map<std::string, bool> seen;
    for (const auto &jk : debug_jkmers) {
        if (jk.seq == "" ) {
            continue;
        }
        if (seen.find(jk.seq) == seen.end()) {
            counts[jk.seq] += 1;
            seen[jk.seq] = true;
        }
    }
    for (Node *child : node->children) {
        getPhyloCounts(counts, jkmers, T, child, debugStrings, debug_k, debug_s);
    }
}
void placeDFS2(std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(cmp)>> &index, std::unordered_map<int32_t, std::string> jkmers, std::map<std::string, float> &scores, std::unordered_map<std::string, int32_t> &jkmerCounts, std::unordered_map<std::string, int32_t> &phyloCounts, const Node *node, Tree *T, std::unordered_map<std::string, std::string> &debugStrings, int32_t debug_k, int32_t debug_s) {
    for (const auto &op : index[node->identifier]) {
        std::string jkmer = op.second;
        int32_t pos = op.first;
        if (jkmer[0] == '@') { // deletion
            jkmers[pos] = "";
        } else {
            jkmers[pos] = jkmer;
        }
    }
    // std::vector<kmer_t> debugNodeSyncmersAln = seed::syncmerize(debugStrings[node->identifier], debug_k, debug_s, false, true, 0);
    // std::vector<kmer_t> debugNodeSyncmers = seed::syncmerize(debugStrings[node->identifier], debug_k, debug_s, false, true, 0);
    // std::vector<jkmer> debugNodeJkmers = seed::jkmerize(debugNodeSyncmers, 3);
    float jkmerSum = 0;
    int32_t numOn = 0;
    for (const auto &jk : jkmers) {
        if (jk.second == "" ) {
            continue;
        }
        numOn++;
        if (jkmerCounts.find(jk.second) != jkmerCounts.end()) {
            jkmerSum += jkmerCounts[jk.second];
        }
    }
    for (const auto &jk : jkmers) {
        if (jk.second == "" ) {
            continue;
        }
    }
    float score = 0;
    int count = 0;
    int32_t sumReadHits = 0;
    int32_t sumNodeHits = 0;
    float nodeJkmerCoverage = 0;
    std::unordered_map<std::string, int32_t> nodeJkmerCounts;

    std::vector<kmer_t> debug_syncmers = syncmerize(debugStrings[node->identifier], debug_k, debug_s, false, true, 0);
    std::vector<jkmer> debug_jkmers = jkmerize(debug_syncmers, 3);

    int32_t max = 0;
    
    for (const auto &jk : debug_jkmers) {
        if (jk.seq == "" ) {
            continue;
        }
        nodeJkmerCounts[jk.seq] += 1;
    }
    for (const auto &jk : debug_jkmers) {
        if (jk.seq == "" ) {
            continue;
        }
        if (jkmerCounts.find(jk.seq) != jkmerCounts.end()) {
            score += std::min(1.0, std::log2(jkmerCounts[jk.seq])/2) / nodeJkmerCounts[jk.seq];
            sumNodeHits += jkmerCounts[jk.seq];
            if (jkmerCounts[jk.seq] > max) {
                max = jkmerCounts[jk.seq];
            }
            int32_t den = phyloCounts[jk.seq];
            sumReadHits += std::min(1.0, std::log2(jkmerCounts[jk.seq])) / den;
        }
    }
    nodeJkmerCoverage = ((float) sumReadHits) / (debug_jkmers.size());
    // # of jkmers in node, # node jkmers hit reads, # read hits to node, unscaled score, score, coverage
    scores[node->identifier] = nodeJkmerCoverage;
    for (Node *child : node->children) {
        placeDFS2(index, jkmers, scores, jkmerCounts, phyloCounts, child, T, debugStrings, debug_k, debug_s);
    }
}
std::string reverse_complement(std::string dna_sequence) {
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


std::vector<kmer_t> seedsFromFastq(std::vector<read_t> &reads, std::unordered_map<std::string, int32_t> &counts, const std::string &fastqPath, const size_t k, const size_t s) {
    util::scopedTimer();
    std::cout << "\n◠◡  Processing reads ... ";
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    seq = kseq_init(fileno(fp));
    std::vector<std::string> input;
    std::vector<std::string> input_quality;
    std::vector<std::string> input_names;
    
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        std::string this_seq  = seq->seq.s;
        std::string this_qual = seq->qual.s;
        std::string this_name = seq->name.s;

        input.push_back(this_seq);
        input_quality.push_back(this_qual);
        input_names.push_back(this_name);
    }
    float est_coverage = 0;
    bool open = false;
    
    std::vector<kmer_t> syncmers;

    reads.resize(input.size());

    for (int i = 0; i < input.size(); i++) {        
        read_t this_read;
        std::string seq = input[i];
        std::string name = input_names[i];
        
        this_read.seq = seq;
        this_read.qseq = input_quality[i];
        this_read.name = name;

        std::string rc = reverse_complement(seq);
        std::vector<kmer_t> these = syncmerize(seq, k, s, false, true, 0);
        std::vector<kmer_t> these_rc = syncmerize(rc, k, s, false, true, 0);
        
        for (const auto &m : these) {
            if (counts.find(m.seq) == counts.end()) {
                counts[m.seq] = 1;
            } else {
                counts[m.seq] += 1;
            }
            if (counts[m.seq] >= est_coverage) {
                syncmers.push_back(m);
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, false, -1});
            }
        }
        for (const auto &m : these_rc) {
            if (counts.find(m.seq) == counts.end()) {
                counts[m.seq] = 1;
            } else {
                counts[m.seq] += 1;
            }
            if (counts[m.seq] >= est_coverage) {
                syncmers.push_back(m);
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, true, -1});
            }
        }
        reads[i] = this_read;
    }
    return syncmers;
}


void place::placeIsolate( std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path, Tree *T, int32_t k, int32_t s) {


    seedIndex index;
//    pmi::load(index, indexFile);

    std::vector<Block> &blocks = T->blocks;
    util::scopedTimer();

    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);

    FILE *fp;
    kseq_t *seq;
    fp = fopen(reads1Path.c_str(), "r");
    seq = kseq_init(fileno(fp));
    std::vector<std::string> input;
    std::vector<std::string> input_names;

    std::unordered_map<std::string, int32_t> jkmerCounts;
    
    std::unordered_map<std::string, std::set<std::pair<int32_t, std::string>, decltype(cmp)>> jkIndex;
    std::unordered_map<int32_t, std::string> dynamicJkmers;
    // read each line in file
    std::string line0;
    std::getline(indexFile, line0);
    std::vector<std::string> spltTop;
    stringSplit(line0, ' ', spltTop);
    k = std::stoi(spltTop[0]);
    s = std::stoi(spltTop[1]);
    std::cout << "index parameters (k,s) = (" << k << "," << s << ")\n";
    while (std::getline(indexFile, line0)) {
        std::vector<std::string> splt;
        stringSplit(line0, ' ', splt);
        std::string nid = splt[0];
        jkIndex[nid] = {};
        for (int32_t i = 1; i < splt.size(); i++) {
            std::vector<std::string> splt2;
            stringSplit(splt[i], ':', splt2);
            int32_t pos = std::stoi(splt2[0]);
            std::string jkmer = splt2[1];
            jkIndex[nid].insert(std::make_pair(pos, jkmer));
        }
    }
    std::cout << "\n⋌⋋  Placing sample ... ";
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        std::string this_seq  = seq->seq.s;
        std::string this_name = seq->name.s;

        input.push_back(this_seq);
        input_names.push_back(this_name);
    }

    for (int i = 0; i < input.size(); i++) {        
        std::string seq = input[i];
        std::string name = input_names[i];
        
        std::string rc = reverse_complement(seq);

        std::vector<kmer_t> syncmers = syncmerize(seq, k, s, false, true, 0);
        std::vector<kmer_t> syncmers_rc = syncmerize(rc, k, s, false, true, 0);

        std::vector<jkmer> jkmers = jkmerize(syncmers, 3);
        std::vector<jkmer> jkmers_rc = jkmerize(syncmers_rc, 3);

        for (const auto &m : jkmers) {
            if (jkmerCounts.find(m.seq) == jkmerCounts.end()) {
                jkmerCounts[m.seq] = 1;
            } else {
                jkmerCounts[m.seq] += 1;
            }
        }
        for (const auto &m : jkmers_rc) {
            if (jkmerCounts.find(m.seq) == jkmerCounts.end()) {
                jkmerCounts[m.seq] = 1;
            } else {
                jkmerCounts[m.seq] += 1;
            }
        }
    }
    std::map<std::string, float> scores;
    std::unordered_map<std::string, std::string> debugStrings = tree::getAllNodeStrings(T);
    std::unordered_map<std::string, int32_t> phyloCounts;
    getPhyloCounts(phyloCounts, dynamicJkmers, T, T->root, debugStrings, k, s);

    placeDFS2(jkIndex, dynamicJkmers, scores, jkmerCounts, phyloCounts, T->root, T, debugStrings, k, s);
    std::vector<std::pair<std::string, float>> v;
    std::copy(scores.begin(), scores.end(),
       back_inserter<std::vector<std::pair<std::string, float> > >(v)
    );
    std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
        return left.second < right.second;
    });
    
    for (size_t i = 0; i < v.size(); ++i) {
       std::cout << v[i].first << ": " << v[i].second << "\n";
    }
}



#include <algorithm>
#include "place.hpp"
#include "pmi.hpp"
#include "tree.hpp"

using namespace PangenomeMAT;
using namespace tree;
using namespace seed;

#define PRINT_TIME(start, end) std::cout << "✔︎ done (took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms)" << std::endl;

struct dynamicJaccard {
    size_t intersectionSize;
    size_t unionSize;
    float jaccardIndex;
};

void updateJaccard(dynamicJaccard &dj, const std::unordered_map<std::string, bool> &readSeeds, const std::vector<kmer_t> &deletedSeeds, const std::vector<kmer_t> &insertedSeeds) {
    for (const kmer_t &syncmer : deletedSeeds) {
        if (readSeeds.find(syncmer.seq) != readSeeds.end()) {
            dj.intersectionSize -= 1;
            dj.unionSize -= 1;
        } else {
            dj.unionSize -= 1;
        }
    }
    for (const kmer_t &syncmer : insertedSeeds) {
        if (readSeeds.find(syncmer.seq) != readSeeds.end()) {
            dj.intersectionSize += 1;
            dj.unionSize += 1;
        } else {
            dj.unionSize += 1;
        }
    }
    dj.jaccardIndex = (float) dj.intersectionSize / dj.unionSize;
}
void initializeJaccard(dynamicJaccard &dj, std::vector<std::string> &nodeSeeds, std::vector<std::string> &readSeeds) {
    std::vector<std::string> v;
    std::sort(nodeSeeds.begin(), nodeSeeds.end());
    std::sort(readSeeds.begin(), readSeeds.end());
    std::set_intersection(nodeSeeds.begin(), nodeSeeds.end(), readSeeds.begin(), readSeeds.end(), std::back_inserter(v));
    std::cout << "Jaccard init:  " << v.size() << "\n";
    dj.intersectionSize = v.size();
    dj.unionSize = nodeSeeds.size() + readSeeds.size() - dj.intersectionSize;
    dj.jaccardIndex = (float) dj.intersectionSize / dj.unionSize;
}
void placeDFS(std::vector<kmer_t> &nodeSeeds, dynamicJaccard &dj, std::unordered_map<std::string, float> &scores, seedIndex &index, const std::unordered_map<std::string, bool> &readSeeds, const Node *node, Tree *T) {
    std::stack<int32_t> rmDel;
    for (const kmer_t &d : index.deletions[node->identifier]) {
        rmDel.push(d.idx);
    }
    removeIndices(nodeSeeds, rmDel);
    for (const kmer_t &s : index.insertions[node->identifier]) {
        nodeSeeds.push_back(s);
    }
    if (node == T->root) {
        std::vector<std::string> nodeSeedStrings;
        for (const kmer_t &s : nodeSeeds) {
            nodeSeedStrings.push_back(s.seq);
        }
        std::vector<std::string> readSeedStrings;
        for (const auto &p : readSeeds) {
            readSeedStrings.push_back(p.first);
        }
        initializeJaccard(dj, readSeedStrings, nodeSeedStrings);
    }
    updateJaccard(dj, readSeeds, index.deletions[node->identifier], index.insertions[node->identifier]);
    scores[node->identifier] = dj.jaccardIndex;

    for (Node *child : node->children) {
        placeDFS(nodeSeeds, dj, scores, index, readSeeds, child, T);
    }

    nodeSeeds.erase(nodeSeeds.end() - index.insertions[node->identifier].size(), nodeSeeds.end());
    for (int32_t i = index.deletions[node->identifier].size() - 1; i >= 0; i--) {
        nodeSeeds.insert(nodeSeeds.begin() + index.deletions[node->identifier][i].idx, index.deletions[node->identifier][i]);
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
std::set<kmer_t> seedsFromFastq(std::vector<read_t> &reads, const std::string &fastqPath, const size_t k, const size_t s) {
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    seq = kseq_init(fileno(fp));
    std::vector<std::string> input;
    std::vector<std::string> input_names;
    
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        std::string this_seq  = seq->seq.s;
        std::string this_name = seq->name.s;

        input.push_back(this_seq);
        input_names.push_back(this_name);
    }
    float est_coverage = 0;
    bool open = false;
    
    std::set<kmer_t> syncmers;
    std::unordered_map<std::string, int> counts;
    std::unordered_map<std::string, int> counts_rc;

    reads.resize(input.size());

    for (int i = 0; i < input.size(); i++) {        
        read_t this_read;
        std::string seq = input[i];
        std::string name = input_names[i];
        
        this_read.seq = seq;
        this_read.name = name;

        std::string rc = reverse_complement(seq);
        std::vector<kmer_t> these = syncmerize(seq, k, s, false, false, 0);
        std::vector<kmer_t> these_rc = syncmerize(rc, k, s, false, false, 0);
        
        for (const auto &m : these) {
            if (counts.find(m.seq) == counts.end()) {
                counts[m.seq] = 1;
            } else {
                counts[m.seq] += 1;
            }
            if (counts[m.seq] > est_coverage) {
                syncmers.insert(m);
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, false, -1});
            }
        }
        
        for (const auto &m : these_rc) {
            if (counts_rc.find(m.seq) == counts_rc.end()) {
                counts_rc[m.seq] = 1;
            } else {
                counts_rc[m.seq] += 1;
            }
            if (counts_rc[m.seq] > est_coverage) {
                syncmers.insert(m);
                this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, true, -1});
            }
        }
        reads[i] = this_read;
    }
    return syncmers;
}

void place::placeIsolate(const std::ifstream &indexFile, const std::string &reads1Path, const std::string &reads2Path, Tree *T, int32_t k, int32_t s) {

    seedIndex index;
    pmi::load(index, T->root, indexFile);
    
    std::cout << "⟗  Processing reads..." << std::endl;
    std::vector<read_t> reads;
    std::vector<Block> &blocks = T->blocks;
    auto fastqStart = std::chrono::high_resolution_clock::now();
    std::set<kmer_t> readSeeds = seedsFromFastq(reads, reads1Path, k, s);
    auto fastqStop = std::chrono::high_resolution_clock::now();
    std::set<kmer_t> consensusSeeds = std::set<kmer_t>(index.consensusSeeds.begin(), index.consensusSeeds.end());

    PRINT_TIME(fastqStart, fastqStop);
    
    auto placeStart = std::chrono::high_resolution_clock::now();
    std::cout << "⋌⋋⋋ Placing sample..." << std::endl;
    std::unordered_map<std::string, bool> off;
    
    // Don't start jaccard computation until we're at the tree root.
    // Here we are still computing on the consensus MSA which != root necessarily
    struct dynamicJaccard dj;
    
    std::unordered_map<std::string, float> scores;
    std::unordered_map<std::string, bool> readSeedsMap;
    for (const auto &k : readSeeds) {
        readSeedsMap[k.seq] = true;
    }

    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    tree::setup(data, globalCoords, T);

    placeDFS(index.consensusSeeds, dj, scores, index, readSeedsMap, T->root, T);

    auto placeStop = std::chrono::high_resolution_clock::now();

    PRINT_TIME(placeStart, placeStop);

    std::vector<std::pair<std::string, float>> v;
    for ( const auto &p : scores ) {
        v.push_back(std::make_pair(p.first, p.second));
    } 
    std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
        return left.second > right.second;
    });

    std::string best_match = v[0].first;

    for (auto &val : v) {
        std::cout << val.first << "\t" << val.second << "\n";
    }
}
#ifndef __SEED_HPP
#define __SEED_HPP

#pragma once
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <mutex>
#include "kseq.h"

KSEQ_INIT(int, read)


/*
minimap2:

anchor: (x, y, w)
    x: ref end position
    y: read end position
    w: k-mer size

a(j, i):  min {
            min ( yi - yj, xi - xj ),
            w[i]
        }
b(j, i): gc( (yi-yj) - (xi-xj) )

gap cost gc(n) =  0.01 * w * n + 0.5 * log2(n)

list of anchors sorted by x:
    chaining score
    f(i) = max of {
            -> max for 1 <= j <= i of f(j) + a(j, i) - b(j, i)
            -> and w[i]
    }



*/

namespace seed {

struct kmer_t {

    std::string seq;
    int32_t pos;  // start position in consensus MSA
    int32_t rpos; // pos on reference (placement dependent)
    int32_t idx;  // index of kmer in vector used in indexing DFS
    int32_t qpos; // query pos
    bool reversed;
    int32_t gappedEnd;

    bool operator<(const kmer_t &rhs) const {
        return pos < rhs.pos;
    };
    bool operator==(const kmer_t &rhs) const {
        return pos == rhs.pos;
    };
};

struct jkmer {
    int32_t j;
    int32_t k;
    std::string seq; // concatenated string of s kmers
    std::vector<int32_t> positions; // positions of s kmers
};

class KHash {
public:
    size_t operator()(const kmer_t& t) const
    {
        const std::hash<std::string> h;
        return h(t.seq);
    }
    
};

struct read_t {
    std::string seq;     //genetic sequence
    std::string qseq;    //quality string
    std::vector<kmer_t> kmers;
    std::string name;
    //std::vector<int> read_coord;
    //std::vector<int> ref_coord;
    //std::vector<bool> reversed;
};


inline bool is_syncmer(const std::string &seq, const int s, const bool open) {
    if (seq.size() < s) {
        return false;
    }
    std::string min(s, 'Z');
    for (size_t i = 0; i < seq.size() -  s + 1; i++) {
        std::string submer = seq.substr(i, s);
        if (submer < min) {
            min = submer;
        }
    }
    if (open) {
        if (min == seq.substr(0, s)) {
            return true;
        }
    } else {
        if (min == seq.substr(0, s) || min == seq.substr(seq.length()-s, s)) {
            return true;
        }
    }
    return false;
}

inline std::vector<jkmer> jkmerize(const std::vector<kmer_t> &kmers, const int32_t j) {
    std::vector <jkmer> ret;
    if (kmers.size() < j) {
        return ret;
    }
    for (size_t i = 0; i < kmers.size() - j + 1; i++) {
        jkmer jk = {j, static_cast<int32_t>(kmers[i].seq.length()), ""};
        for (size_t k = 0; k < j; k++) {
            jk.positions.push_back(kmers[i+k].pos);
            jk.seq += kmers[i+k].seq;
        }
        ret.push_back(jk);
    }
    return ret; 
}

inline std::string getNextSyncmer(std::string &seq, const int32_t currPos, const int32_t k, const int32_t s) {

    for (int32_t i = currPos; i < seq.size() - k + 1; i++) {
        std::string kmer = seq.substr(i, k);
        if (is_syncmer(kmer, s, false)) {
            return kmer;
        }
    }
    return "";
}
inline std::vector<kmer_t> syncmerize(const std::string &seq, const int32_t k, const int32_t s, const bool open, const bool aligned, const int32_t pad) {
    std::mutex mtx;
    std::vector<kmer_t> ret;
    mtx.lock();
    int32_t seqLen = seq.size();
    mtx.unlock();
    if (seqLen < k) {
        return ret;
    }
    if (aligned) {
        std::unordered_map<int32_t, int32_t> degap;
        int32_t pos = 0;
        std::string ungapped = "";
        for (int32_t i = 0; i < seqLen; i++) {
            char c = seq[i];
            degap[pos] = i;
            if (c != '-') {
                ungapped += c;
                pos++;
            }
        }

        if (ungapped.size() < k + 1) {
            return ret;
        }
        
        for(int32_t i = 0; i < ungapped.size() - k + 1; i++) {
            std::string kmer = ungapped.substr(i, k);
            if (is_syncmer(kmer, s, open)) {
                ret.push_back(kmer_t{kmer, degap[i]+pad, -1, -1, false, degap[i+k-1]+pad});
            }
        }
    } else {
        for (int32_t i = 0; i < seqLen - k + 1; i++) {
            mtx.lock();
            std::string kmer = seq.substr(i, k);
            mtx.unlock();
            if (is_syncmer(kmer, s, open)) {
                ret.push_back(kmer_t{kmer, i+pad, -1, -1, false, i+k+pad});
            }
        }
    }
    return ret;
}
}

// std::set<kmer_t> PangenomeMAT::syncmersFromFastq(std::string fastqPath,  std::vector<read_t> &reads, size_t k, size_t s) {
//     FILE *fp;
//     kseq_t *seq;
//     fp = fopen(fastqPath.c_str(), "r");
//     seq = kseq_init(fileno(fp));
//     std::vector<std::string> input;
//     std::vector<std::string> input_names;
    
//     int line;
//     while ((line = kseq_read(seq)) >= 0) {
//         std::string this_seq  = seq->seq.s;
//         std::string this_name = seq->name.s;

//         input.push_back(this_seq);
//         input_names.push_back(this_name);
//     }
//     float est_coverage = 1; //TODO change this to 1
//     bool open = false;
    
//     std::set<kmer_t> syncmers;
//     std::unordered_map<std::string, int> counts;
//     std::unordered_map<std::string, int> counts_rc;

//     std::cerr << "length: " << input.size() << "\n";
//     reads.resize(input.size());

//     for (int i = 0; i < input.size(); i++) {        
//         read_t this_read;
//         std::string seq = input[i];
//         std::string name = input_names[i];
        
//         this_read.seq = seq;
//         this_read.name = name;


//         std::string rc = reverse_complement(seq);
//         std::vector<kmer_t> these = syncmerize(seq, k, s, false, false, 0);
//         std::vector<kmer_t> these_rc = syncmerize(rc, k, s, false, false, 0);
        
        
//         for (const auto &m : these) {
//             if (counts.find(m.seq) == counts.end()) {
//                 counts[m.seq] = 1;
//             } else {
//                 counts[m.seq] += 1;
//             }
//             if (counts[m.seq] > est_coverage) {
//                 syncmers.insert(m);

//  //               m.pos = m.pos + k - 1;
// //                m.reversed = false;
//                 this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, false, -1});
//                 //this_read.read_coord.push_back(m.pos + k - 1);
//                 //this_read.reversed.push_back(false);
//             }
//         }
        
//         for (const auto &m : these_rc) {
//             if (counts_rc.find(m.seq) == counts_rc.end()) {
//                 counts_rc[m.seq] = 1;
//             } else {
//                 counts_rc[m.seq] += 1;
//             }
//             if (counts_rc[m.seq] > est_coverage) {
//                 syncmers.insert(m);
//                 this_read.kmers.push_back(kmer_t{m.seq, m.pos + (int32_t) k - 1, -1, 0, true, -1});
//             }
//         }
//         reads[i] = this_read;
//     }
   
//     return syncmers;
// }

#endif
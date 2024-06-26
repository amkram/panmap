#ifndef __SEEDING_HPP
#define __SEEDING_HPP

#pragma once
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <mutex>
#include <zlib.h>  
#include <stdio.h>  
#include <minimap2/kseq.h>
#include <unordered_map>


KSEQ_INIT(int, read)

namespace seeding {

struct seed {

    std::string seq;
    int32_t pos;  // start position in sequence
    int32_t idx;  // index of kmer in vector used in indexing DFS
    bool reversed;
    int32_t gappedEnd;
    int32_t rpos; // position on reference (only used during alignment)

    bool operator<(const seed &rhs) const {
        return pos < rhs.pos;
    };
    bool operator==(const seed &rhs) const {
        return pos == rhs.pos;
    };
};

struct seedmer {
    int32_t j;
    int32_t k;
    std::string seq; // concatenated string of s kmers
    std::vector<int32_t> positions; // positions of s kmers
};

class KHash {   
public:
    size_t operator()(const seed& t) const  
    {
        const std::hash<std::string> h;
        return h(t.seq);
    }
    
};

struct read_t {
    std::string seq; // read sequence
    std::string qual; // quality string
    std::vector<int32_t> seedPositions; // positions of seed starts in read
    std::string name; // read id
};

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

inline std::pair<size_t, bool> getHash(const std::string& s) {
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

inline std::vector<seedmer> seedmerize(const std::vector<seed> &kmers, const int32_t j) {
    std::vector <seedmer> ret;
    if (kmers.size() < j) {
        return ret;
    }
    for (size_t i = 0; i < kmers.size() - j + 1; i++) {
        seedmer jk = {j, static_cast<int32_t>(kmers[i].seq.length()), ""};
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
inline std::vector<seed> syncmerize(const std::string &seq, const int32_t k, const int32_t s, const bool open, const bool aligned, const int32_t pad) {
    std::mutex mtx;
    std::vector<seed> ret;
    int32_t seqLen = seq.size();
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
                ret.push_back(seed{kmer, degap[i]+pad, -1, false, degap[i+k-1]+pad});
            }
        }
    } else {
        for (int32_t i = 0; i < seqLen - k + 1; i++) {
            mtx.lock();
            std::string kmer = seq.substr(i, k);
            mtx.unlock();
            if (is_syncmer(kmer, s, open)) {
                ret.push_back(seed{kmer, i+pad, -1, false, i+k+pad});
            }
        }
    }
    return ret;
}
}

#endif
#pragma once

#include <cstdint>
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#include <immintrin.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <string_view>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <deque>
#include <stdio.h>
#include <tuple>
#include <zlib.h>
#include <cctype>

#include <htslib/kseq.h>
KSEQ_INIT(int, read)

#include <absl/container/flat_hash_map.h>

namespace seeding {

// Interleave paired-end reads: R1_0, R2_0, R1_1, R2_1, ...
template <typename T>
inline void perfect_shuffle(std::vector<T>& v) {
    const size_t n = v.size();
    if (n < 2) return;

    std::vector<T> canvas(n);
    for (size_t i = 0; i < n / 2; i++) {
        canvas[i * 2] = std::move(v[i]);
        canvas[i * 2 + 1] = std::move(v[i + n / 2]);
    }
    v = std::move(canvas);
}

struct rsyncmer_t {
    size_t hash;
    uint64_t endPos;
    bool isReverse;

    rsyncmer_t(size_t hash, uint64_t endPos, bool isReverse) : hash(hash), endPos(endPos), isReverse(isReverse) {}

    rsyncmer_t() : hash(0), endPos(0), isReverse(false) {}
};

struct uniqueKminmer_t {
    uint64_t startPos;
    uint64_t endPos;
    size_t hash;
    bool isReverse;

    uniqueKminmer_t(uint64_t startPos, uint64_t endPos, size_t hash, bool isReverse)
        : startPos(startPos), endPos(endPos), hash(hash), isReverse(isReverse) {}

    uniqueKminmer_t() : startPos(0), endPos(0), hash(0), isReverse(false) {}

    bool operator==(const uniqueKminmer_t& other) const {
        return startPos == other.startPos && endPos == other.endPos && hash == other.hash &&
               isReverse == other.isReverse;
    }
};

struct seed_t {
    size_t hash;
    int64_t pos;
    int32_t idx;    // set during indexing
    bool reversed;
    uint32_t rpos;  // set during alignment

    int64_t endPos;  // counting gaps

    bool operator<(const seed_t& rhs) const { return pos < rhs.pos; }

    bool operator==(const seed_t& rhs) const { return pos == rhs.pos; }
};

constexpr char comp(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'a': return 't';
        case 'C': return 'G';
        case 'c': return 'g';
        case 'G': return 'C';
        case 'g': return 'c';
        case 'T': return 'A';
        case 't': return 'a';
        default: return 'N';
    }
}

constexpr size_t chash(char c) {
    switch (c) {
        case 'a':
        case 'A': return 0x3c8bfbb395c60474;
        case 'c':
        case 'C': return 0x3193c18562a02b4c;
        case 'g':
        case 'G': return 0x20323ed082572324;
        case 't':
        case 'T': return 0x295549f54be24456;
        default: return 0;
    }
}

constexpr size_t rol(size_t h, size_t r) {
    return (h << r) | (h >> (64 - r));
}

constexpr size_t ror(size_t h, size_t r) {
    return (h >> r) | (h << (64 - r));
}

// Hash a sequence, returns (forward hash, reverse complement hash)
std::pair<size_t, size_t> hashSeq(const std::string& s);

// returns vector of (hash, isReverse, isSyncmer, startPos)
std::vector<std::tuple<size_t, bool, bool, int64_t>>
rollingSyncmers(std::string_view seq, int k, int s, bool open, int t = 0, bool returnAll = true);

std::string reverseComplement(std::string dna_sequence);

// Paired-end (fastqPath2 set): reverse-complements R2 and interleaves pairs.
void readFastqPaired(std::vector<std::string>& readSequences,
                     std::vector<std::string>& readQuals,
                     std::vector<std::string>& readNames,
                     const std::string& fastqPath1,
                     const std::string& fastqPath2);

// Homopolymer compression: collapse consecutive identical bases (e.g., AAACGG -> ACG)
std::string hpcCompress(const std::string& seq);

// mapping[i] = index in original seq of the i-th base in the compressed seq
std::pair<std::string, std::vector<size_t>> hpcCompressWithMapping(const std::string& seq);

}  // namespace seeding

namespace std {
template <>
struct hash<seeding::uniqueKminmer_t> {
    size_t operator()(const seeding::uniqueKminmer_t& kminmer) const {
        size_t h1 = std::hash<uint64_t>{}(kminmer.startPos);
        size_t h2 = std::hash<uint64_t>{}(kminmer.endPos);
        size_t h3 = std::hash<size_t>{}(kminmer.hash);
        size_t h4 = std::hash<bool>{}(kminmer.isReverse);

        static constexpr size_t HASH_MIX = 0x9e3779b9;
        size_t seed = 0;
        seed ^= h1 + HASH_MIX + (seed << 6) + (seed >> 2);
        seed ^= h2 + HASH_MIX + (seed << 6) + (seed >> 2);
        seed ^= h3 + HASH_MIX + (seed << 6) + (seed >> 2);
        seed ^= h4 + HASH_MIX + (seed << 6) + (seed >> 2);
        return seed;
    }
};
}  // namespace std

#include "seeding.hpp"
#include "mgsr.hpp"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <cctype>
#include <absl/container/flat_hash_map.h>

namespace seeding {

std::pair<size_t, size_t> hashSeq(const std::string& s) {
    size_t fHash = 0;
    size_t rHash = 0;
    int k = s.size();
    for (int i = 0; i < k; i++) {
        if (chash(s[i]) == 0) throw std::invalid_argument("Kmer contains non canonical base");
        fHash ^= rol(chash(s[i]), k - i - 1);
        rHash ^= rol(chash(comp(s[k - i - 1])), k - i - 1);
    }
    return std::make_pair(fHash, rHash);
}

std::vector<std::tuple<size_t, bool, bool, int64_t>>
rollingSyncmers(std::string_view seq, int k, int s, bool open, int t, bool returnAll) {
    std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmers;
    if (seq.size() < k) return syncmers;
    syncmers.reserve(seq.size() - k + 1);

    constexpr size_t max_size_t = std::numeric_limits<size_t>::max();
    size_t forwardKmerHash = 0, reverseKmerHash = 0;
    size_t forwardSmerHash = 0, reverseSmerHash = 0;

    int recentAmbiguousBaseIndex = std::numeric_limits<int>::min();

    size_t curMinSmerHashForward = max_size_t;
    size_t curMinSmerHashReverse = max_size_t;
    int curMinSmerHashIndexForward = -1;
    int curMinSmerHashIndexReverse = -1;
    std::deque<size_t> curSmerHashesForward;
    std::deque<size_t> curSmerHashesReverse;

    // first smer
    for (int i = 0; i < s; ++i) {
        size_t fwd_hash = chash(seq[i]);
        size_t rev_hash = chash(comp(seq[k - i - 1]));
        size_t smer_rev_hash = chash(comp(seq[s - i - 1]));

        if (fwd_hash == 0) recentAmbiguousBaseIndex = i;

        forwardKmerHash ^= rol(fwd_hash, k - i - 1);
        reverseKmerHash ^= rol(rev_hash, k - i - 1);
        forwardSmerHash ^= rol(fwd_hash, s - i - 1);
        reverseSmerHash ^= rol(smer_rev_hash, s - i - 1);
    }

    curSmerHashesForward.push_back(forwardSmerHash);
    curMinSmerHashForward = forwardSmerHash;
    curMinSmerHashIndexForward = 0;
    curSmerHashesReverse.push_back(reverseSmerHash);
    curMinSmerHashReverse = reverseSmerHash;
    curMinSmerHashIndexReverse = k - s;

    // first kmer
    for (int i = s; i < k; ++i) {
        size_t fwd_hash = chash(seq[i]);
        size_t rev_hash = chash(comp(seq[i]));
        size_t old_fwd_hash = chash(seq[i - s]);
        size_t old_rev_hash = chash(comp(seq[i - s]));
        size_t rev_kmer_hash = chash(comp(seq[k - i - 1]));

        if (fwd_hash == 0) recentAmbiguousBaseIndex = i;

        forwardKmerHash ^= rol(fwd_hash, k - i - 1);
        reverseKmerHash ^= rol(rev_kmer_hash, k - i - 1);
        forwardSmerHash = rol(forwardSmerHash, 1) ^ rol(old_fwd_hash, s) ^ fwd_hash;
        reverseSmerHash = ror(reverseSmerHash, 1) ^ ror(old_rev_hash, 1) ^ rol(rev_hash, s - 1);

        curSmerHashesForward.push_back(forwardSmerHash);
        if (forwardSmerHash < curMinSmerHashForward) {
            curMinSmerHashForward = forwardSmerHash;
            curMinSmerHashIndexForward = i - s + 1;
        }
        curSmerHashesReverse.push_front(reverseSmerHash);
        if (reverseSmerHash < curMinSmerHashReverse) {
            curMinSmerHashReverse = reverseSmerHash;
            curMinSmerHashIndexReverse = k - i - 1;
        }
    }

    if (recentAmbiguousBaseIndex >= 0) {
        if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
    } else {
        bool forwardIsSyncmer = false;
        bool reverseIsSyncmer = false;

        if (open) {
            if (curSmerHashesForward[t] == curMinSmerHashForward) forwardIsSyncmer = true;
            if (curSmerHashesReverse[t] == curMinSmerHashReverse) reverseIsSyncmer = true;
        } else {
            if (curSmerHashesForward[t] == curMinSmerHashForward ||
                curSmerHashesForward[k - s - t] == curMinSmerHashForward)
                forwardIsSyncmer = true;
            if (curSmerHashesReverse[t] == curMinSmerHashReverse ||
                curSmerHashesReverse[k - s - t] == curMinSmerHashReverse)
                reverseIsSyncmer = true;
        }

        if (forwardIsSyncmer || reverseIsSyncmer) {
            if (forwardKmerHash < reverseKmerHash) {
                syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, 0));
            } else if (reverseKmerHash < forwardKmerHash) {
                syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, 0));
            } else {
                if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
            }
        } else {
            if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
        }
    }

    // start rolling
    for (size_t i = k; i < seq.size(); ++i) {
        if (chash(seq[i]) == 0) recentAmbiguousBaseIndex = i;
        curSmerHashesForward.pop_front();
        curSmerHashesReverse.pop_back();
        --curMinSmerHashIndexForward;
        ++curMinSmerHashIndexReverse;
        if (curMinSmerHashIndexForward < 0) {
            curMinSmerHashForward = max_size_t;
            for (size_t j = 0; j < curSmerHashesForward.size(); ++j) {
                if (curSmerHashesForward[j] < curMinSmerHashForward) {
                    curMinSmerHashForward = curSmerHashesForward[j];
                    curMinSmerHashIndexForward = j;
                }
            }
        }
        if (curMinSmerHashIndexReverse > k - s) {
            curMinSmerHashReverse = max_size_t;
            for (size_t j = 0; j < curSmerHashesReverse.size(); ++j) {
                if (curSmerHashesReverse[j] < curMinSmerHashReverse) {
                    curMinSmerHashReverse = curSmerHashesReverse[j];
                    curMinSmerHashIndexReverse = j;
                }
            }
            curMinSmerHashIndexReverse++;
        }

        size_t fwd_hash = chash(seq[i]);
        size_t rev_hash = chash(comp(seq[i]));
        size_t old_kmer_fwd = chash(seq[i - k]);
        size_t old_kmer_rev = chash(comp(seq[i - k]));
        size_t old_smer_fwd = chash(seq[i - s]);
        size_t old_smer_rev = chash(comp(seq[i - s]));

        forwardKmerHash = rol(forwardKmerHash, 1) ^ rol(old_kmer_fwd, k) ^ fwd_hash;
        reverseKmerHash = ror(reverseKmerHash, 1) ^ ror(old_kmer_rev, 1) ^ rol(rev_hash, k - 1);
        forwardSmerHash = rol(forwardSmerHash, 1) ^ rol(old_smer_fwd, s) ^ fwd_hash;
        reverseSmerHash = ror(reverseSmerHash, 1) ^ ror(old_smer_rev, 1) ^ rol(rev_hash, s - 1);

        curSmerHashesForward.push_back(forwardSmerHash);
        if (forwardSmerHash < curMinSmerHashForward) {
            curMinSmerHashForward = forwardSmerHash;
            curMinSmerHashIndexForward = k - s;
        }
        curSmerHashesReverse.push_front(reverseSmerHash);
        if (reverseSmerHash < curMinSmerHashReverse) {
            curMinSmerHashReverse = reverseSmerHash;
            curMinSmerHashIndexReverse = 0;
        }

        if (recentAmbiguousBaseIndex >= 0 && i < recentAmbiguousBaseIndex + k) {
            if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
        } else {
            bool forwardIsSyncmer = false;
            bool reverseIsSyncmer = false;

            if (open) {
                if (curSmerHashesForward[t] == curMinSmerHashForward) forwardIsSyncmer = true;
                if (curSmerHashesReverse[t] == curMinSmerHashReverse) reverseIsSyncmer = true;
            } else {
                if (curSmerHashesForward[t] == curMinSmerHashForward ||
                    curSmerHashesForward[k - s - t] == curMinSmerHashForward)
                    forwardIsSyncmer = true;
                if (curSmerHashesReverse[t] == curMinSmerHashReverse ||
                    curSmerHashesReverse[k - s - t] == curMinSmerHashReverse)
                    reverseIsSyncmer = true;
            }

            if (forwardIsSyncmer || reverseIsSyncmer) {
                if (forwardKmerHash < reverseKmerHash) {
                    syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, i - k + 1));
                } else if (reverseKmerHash < forwardKmerHash) {
                    syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, i - k + 1));
                } else {
                    if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
                }
            } else {
                if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
            }
        }
    }
    syncmers.shrink_to_fit();
    return syncmers;
}

void readFastqPaired(std::vector<std::string>& readSequences,
                     std::vector<std::string>& readQuals,
                     std::vector<std::string>& readNames,
                     const std::string& fastqPath1,
                     const std::string& fastqPath2) {
    mgsr::FastqFile fq1(fastqPath1);
    kseq_t* seq = kseq_init(fileno(fq1.fp));
    while (kseq_read(seq) >= 0) {
        readSequences.push_back(seq->seq.s);
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.l > 0 ? seq->qual.s : std::string(seq->seq.l, 'I'));
    }
    kseq_destroy(seq);

    if (!fastqPath2.empty()) {
        mgsr::FastqFile fq2(fastqPath2);
        seq = kseq_init(fileno(fq2.fp));
        int forwardReads = readSequences.size();
        while (kseq_read(seq) >= 0) {
            // R2 is reverse-complemented; reverse its per-base quals to stay aligned with the bases.
            readSequences.push_back(reverseComplement(seq->seq.s));
            readNames.push_back(seq->name.s);
            std::string r2qual = seq->qual.l > 0 ? std::string(seq->qual.s) : std::string(seq->seq.l, 'I');
            std::reverse(r2qual.begin(), r2qual.end());
            readQuals.push_back(std::move(r2qual));
        }
        kseq_destroy(seq);

        if (readSequences.size() != static_cast<size_t>(forwardReads) * 2) {
            std::cerr << "Error: " << fastqPath2 << " does not contain the same number of reads as " << fastqPath1
                      << std::endl;
            exit(1);
        }

        perfect_shuffle(readSequences);
        perfect_shuffle(readNames);
        perfect_shuffle(readQuals);
    }
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

std::string hpcCompress(const std::string& seq) {
    // Same homopolymer-collapse rule as hpcCompressWithMapping, without the position map.
    return hpcCompressWithMapping(seq).first;
}

std::pair<std::string, std::vector<size_t>> hpcCompressWithMapping(const std::string& seq) {
    if (seq.empty()) return {{}, {}};
    std::string compressed;
    std::vector<size_t> mapping;
    compressed.reserve(seq.size());
    mapping.reserve(seq.size());
    compressed.push_back(seq[0]);
    mapping.push_back(0);
    for (size_t i = 1; i < seq.size(); ++i) {
        if (std::toupper(seq[i]) != std::toupper(seq[i - 1])) {
            compressed.push_back(seq[i]);
            mapping.push_back(i);
        }
    }
    return {compressed, mapping};
}

}  // namespace seeding

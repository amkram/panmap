#ifndef __SEEDING_HPP
#define __SEEDING_HPP

#pragma once
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <minimap2/kseq.h>
#include <mutex>
#include <optional>
#include <stdio.h>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <zlib.h>

KSEQ_INIT(int, read)

namespace seeding {

struct seed {

  // std::string seq;
  size_t hash;
  int32_t pos; // start position in sequence
  int32_t idx; // index of kmer in vector used in indexing DFS
  bool reversed;
  int32_t gappedEnd;
  uint32_t rpos; // position on reference (only used during alignment)

  bool operator<(const seed &rhs) const { return pos < rhs.pos; };
  bool operator==(const seed &rhs) const { return pos == rhs.pos; };
};

struct onSeedsHash {
  size_t hash;
  int64_t endPos;
  bool isReverse;
};

struct seedmer {
  int32_t j;
  int32_t k;
  std::string seq;                // concatenated string of s kmers
  std::vector<int32_t> positions; // positions of s kmers
};

// class KHash
// {
// public:
//     size_t operator()(const seed &t) const
//     {
//         const std::hash<std::string> h;
//         return h(t.seq);
//     }
// };

struct read_t {
  std::string seq;                    // read sequence
  std::string qual;                   // quality string
  std::vector<int32_t> seedPositions; // positions of seed starts in read
  std::string name;                   // read id
};

static size_t btn(char b) {
  size_t n;
  switch (b) {
  case 'A':
    n = 0;
    break;
  case 'a':
    n = 0;
    break;
  case 'C':
    n = 1;
    break;
  case 'c':
    n = 1;
    break;
  case 'G':
    n = 2;
    break;
  case 'g':
    n = 2;
    break;
  case 'T':
    n = 3;
    break;
  case 't':
    n = 3;
    break;
  default:
    throw std::invalid_argument("Kmer contains non canonical base");
    break;
  }
  return n;
}

static size_t hash(const std::string &s) {
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
  case 'a':
    compC = 't';
    break;
  case 'C':
    compC = 'G';
    break;
  case 'c':
    compC = 'g';
    break;
  case 'G':
    compC = 'C';
    break;
  case 'g':
    compC = 'c';
    break;
  case 'T':
    compC = 'A';
    break;
  case 't':
    compC = 'a';
    break;
  default:
    compC = 'N';
    break;
  }

  return compC;
}

static std::string revcomp(const std::string &s) {
  std::string cs = "";
  for (int i = s.size() - 1; i > -1; --i) {
    char c = s[i];
    cs += comp(c);
  }
  return cs;
}

inline std::pair<size_t, bool> getHash(const std::string &s) {
  try {
    size_t u = hash(s);
    size_t v = hash(revcomp(s));
    if (u < v)
      return std::make_pair(u, true);
    else if (v < u)
      return std::make_pair(v, true);
    return std::make_pair(0, false);
  } catch (std::invalid_argument) {
    return std::make_pair(0, false); // skip strand ambiguous
  }
  return std::make_pair(0, false);
}

inline bool is_syncmer(const std::string &seq, const int s, const bool open) {
  int NsCount = 0;
  for (size_t i = 0; i < seq.size(); i++) {
    if (seq[i] == 'N') {
      NsCount++;
    }
  }
  if (NsCount > seq.size() / 2) {
    return false;
  }
  if (seq.size() < s) {
    return false;
  }
  std::string min(s, 'Z');
  for (size_t i = 0; i < seq.size() - s + 1; i++) {
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
    if (min == seq.substr(0, s) || min == seq.substr(seq.length() - s, s)) {
      return true;
    }
  }
  return false;
}

inline size_t chash(const char &c) {
  switch (c) {
  case 'a':
  case 'A':
    return 0x3c8bfbb395c60474;
  case 'c':
  case 'C':
    return 0x3193c18562a02b4c;
  case 'g':
  case 'G':
    return 0x20323ed082572324;
  case 't':
  case 'T':
    return 0x295549f54be24456;
  default:
    // throw std::invalid_argument("Kmer contains non canonical base");
    return 0;
  }
  return 0;
}

inline size_t rol(const size_t &h, const size_t &r) {
  return (h << r) | (h >> (64 - r));
}

inline size_t ror(const size_t &h, const size_t &r) {
  return (h >> r) | (h << (64 - r));
}

// returns a tuple of (hash, isReverse, isSyncmer)
inline std::tuple<size_t, bool, bool>
is_syncmer_rollingHash(const std::string &seq, const int s, const bool open,
                       int t = 0) {
  if (seq.size() < s)
    return std::make_tuple(0, false, false);

  int k = seq.size();
  size_t forwardKmerHash = 0, reverseKmerHash = 0;
  size_t forwardSmerHash = 0, reverseSmerHash = 0;

  size_t minSmerHash = std::numeric_limits<size_t>::max();
  size_t tthSmerHash = std::numeric_limits<size_t>::max();
  size_t lastthSmerHash = std::numeric_limits<size_t>::min();

  for (int i = 0; i < s; ++i) {
    if (chash(seq[i]) == 0)
      return std::make_tuple(0, false, false);

    forwardKmerHash ^= rol(chash(seq[i]), k - i - 1);
    reverseKmerHash ^= rol(chash(comp(seq[k - i - 1])), k - i - 1);
    forwardSmerHash ^= rol(chash(seq[i]), s - i - 1);
    reverseSmerHash ^= rol(chash(comp(seq[s - i - 1])), s - i - 1);
  }

  if (forwardSmerHash < reverseSmerHash) {
    minSmerHash = forwardSmerHash;
    if (t == 0)
      tthSmerHash = forwardSmerHash;
  } else if (reverseSmerHash < forwardSmerHash) {
    minSmerHash = reverseSmerHash;
    if (t == 0)
      tthSmerHash = reverseSmerHash;
  }

  for (int i = 1; i < k - s + 1; ++i) {
    if (chash(seq[i + s - 1]) == 0)
      return std::make_tuple(0, false, false);

    forwardKmerHash ^= rol(chash(seq[i + s - 1]), k - s - i);
    reverseKmerHash ^= rol(chash(comp(seq[k - s - i])), k - s - i);
    forwardSmerHash = rol(forwardSmerHash, 1) ^ rol(chash(seq[i - 1]), s) ^
                      chash(seq[i + s - 1]);
    reverseSmerHash = ror(reverseSmerHash, 1) ^
                      ror(chash(comp(seq[i - 1])), 1) ^
                      rol(chash(comp(seq[i + s - 1])), s - 1);

    if (i == t) {
      if (forwardSmerHash < reverseSmerHash)
        tthSmerHash = forwardSmerHash;
      else if (reverseSmerHash < forwardSmerHash)
        tthSmerHash = reverseSmerHash;
    } else if (i == k - s - t) {
      if (forwardSmerHash < reverseSmerHash)
        lastthSmerHash = forwardSmerHash;
      else if (reverseSmerHash < forwardSmerHash)
        lastthSmerHash = reverseSmerHash;
    }

    if (forwardSmerHash < reverseSmerHash && forwardSmerHash < minSmerHash)
      minSmerHash = forwardSmerHash;
    else if (reverseSmerHash < forwardSmerHash && reverseSmerHash < minSmerHash)
      minSmerHash = reverseSmerHash;
  }

  if (minSmerHash == std::numeric_limits<size_t>::max())
    return std::make_tuple(0, false, false);

  if (forwardKmerHash < reverseKmerHash) {
    if (open) {
      if (minSmerHash == tthSmerHash)
        return std::make_tuple(forwardKmerHash, false, true);
    } else {
      if (minSmerHash == tthSmerHash || minSmerHash == lastthSmerHash)
        return std::make_tuple(forwardKmerHash, false, true);
    }
  } else if (reverseKmerHash < forwardKmerHash) {
    if (open) {
      if (minSmerHash == tthSmerHash)
        return std::make_tuple(reverseKmerHash, true, true);
    } else {
      if (minSmerHash == tthSmerHash || minSmerHash == lastthSmerHash)
        return std::make_tuple(reverseKmerHash, true, true);
    }
  }

  return std::make_tuple(0, false, false);
}

// Overload to support string_view
inline std::pair<size_t, size_t> hashSeq(std::string_view s) {
  size_t fHash = 0;
  size_t rHash = 0;
  int k = s.size();
  for (int i = 0; i < k; i++) {
    if (chash(s[i]) == 0)
      throw std::invalid_argument("Kmer contains non canonical base");
    fHash ^= rol(chash(s[i]), k - i - 1);
    rHash ^= rol(chash(comp(s[k - i - 1])), k - i - 1);
  }
  return std::make_pair(fHash, rHash);
}

// returns vector of (hash, isReverse, isSyncmer, startPos)
inline std::vector<std::tuple<size_t, bool, bool, int64_t>>
rollingSyncmers(const std::string_view &seq, int k, int s, bool open, int t = 0,
                bool returnAll = true) {
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmers;

  const size_t max_size_t = std::numeric_limits<size_t>::max();
  size_t forwardKmerHash = 0, reverseKmerHash = 0;
  size_t forwardSmerHash = 0, reverseSmerHash = 0;

  int recentAmbiguousBaseIndex = std::numeric_limits<int>::min();

  size_t curMinSmerHash = max_size_t;
  int curMinSmerHashIndex = -1;
  std::deque<size_t> curSmerHashes;

  // Debug flag - set to true to print k-mer strings
  bool printKmers = true;

  // first smer
  for (int i = 0; i < s; ++i) {
    if (chash(seq[i]) == 0)
      recentAmbiguousBaseIndex = i;

    forwardKmerHash ^= rol(chash(seq[i]), k - i - 1);
    reverseKmerHash ^= rol(chash(comp(seq[k - i - 1])), k - i - 1);
    forwardSmerHash ^= rol(chash(seq[i]), s - i - 1);
    reverseSmerHash ^= rol(chash(comp(seq[s - i - 1])), s - i - 1);
  }

  if (forwardSmerHash < reverseSmerHash) {
    curSmerHashes.push_back(forwardSmerHash);
    if (forwardSmerHash < curMinSmerHash) {
      curMinSmerHash = forwardSmerHash;
      curMinSmerHashIndex = 0;
    }
  } else if (reverseSmerHash < forwardSmerHash) {
    curSmerHashes.push_back(reverseSmerHash);
    if (reverseSmerHash < curMinSmerHash) {
      curMinSmerHash = reverseSmerHash;
      curMinSmerHashIndex = 0;
    }
  } else {
    curSmerHashes.push_back(max_size_t);
  }

  // first kmer
  for (int i = s; i < k; ++i) {
    if (chash(seq[i]) == 0)
      recentAmbiguousBaseIndex = i;

    forwardKmerHash ^= rol(chash(seq[i]), k - i - 1);
    reverseKmerHash ^= rol(chash(comp(seq[k - i - 1])), k - i - 1);
    forwardSmerHash =
        rol(forwardSmerHash, 1) ^ rol(chash(seq[i - s]), s) ^ chash(seq[i]);
    reverseSmerHash = ror(reverseSmerHash, 1) ^
                      ror(chash(comp(seq[i - s])), 1) ^
                      rol(chash(comp(seq[i])), s - 1);

    if (forwardSmerHash < reverseSmerHash) {
      curSmerHashes.push_back(forwardSmerHash);
      if (forwardSmerHash < curMinSmerHash) {
        curMinSmerHash = forwardSmerHash;
        curMinSmerHashIndex = i - s + 1;
      }
    } else if (reverseSmerHash < forwardSmerHash) {
      curSmerHashes.push_back(reverseSmerHash);
      if (reverseSmerHash < curMinSmerHash) {
        curMinSmerHash = reverseSmerHash;
        curMinSmerHashIndex = i - s + 1;
      }
    } else {
      curSmerHashes.push_back(max_size_t);
    }
  }

  if (recentAmbiguousBaseIndex >= 0) {
    if (returnAll)
      syncmers.emplace_back(std::make_tuple(0, false, false, 0));
  } else {
    if (forwardKmerHash < reverseKmerHash) {
      if (open) {
        if (curMinSmerHash == curSmerHashes[t]) {
          syncmers.emplace_back(
              std::make_tuple(forwardKmerHash, false, true, 0));
          if (printKmers) {
            std::string kmer = std::string(seq.substr(0, k));
            // std::cout << "SYNCMER: pos=" << 0 << " kmer=" << kmer << " hash="
            // << forwardKmerHash << " orientation=forward" << std::endl;
          }
        } else {
          if (returnAll)
            syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
        }
      } else {
        if (curMinSmerHash == curSmerHashes[t] ||
            curMinSmerHash == curSmerHashes[k - s - t]) {
          syncmers.emplace_back(
              std::make_tuple(forwardKmerHash, false, true, 0));
          if (printKmers) {
            std::string kmer = std::string(seq.substr(0, k));
            // std::cout << "SYNCMER: pos=" << 0 << " kmer=" << kmer << " hash="
            // << forwardKmerHash << " orientation=forward" << std::endl;
          }
        } else {
          if (returnAll)
            syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
        }
      }
    } else if (reverseKmerHash < forwardKmerHash) {
      if (open) {
        if (curMinSmerHash == curSmerHashes[t]) {
          syncmers.emplace_back(
              std::make_tuple(reverseKmerHash, true, true, 0));
          if (printKmers) {
            std::string kmer = std::string(seq.substr(0, k));
            // std::cout << "SYNCMER: pos=" << 0 << " kmer=" << kmer << " hash="
            // << reverseKmerHash << " orientation=reverse" << std::endl;
          }
        } else {
          if (returnAll)
            syncmers.emplace_back(std::make_tuple(max_size_t, true, false, 0));
        }
      } else {
        if (curMinSmerHash == curSmerHashes[t] ||
            curMinSmerHash == curSmerHashes[k - s - t]) {
          syncmers.emplace_back(
              std::make_tuple(reverseKmerHash, true, true, 0));
          if (printKmers) {
            std::string kmer = std::string(seq.substr(0, k));
            // std::cout << "SYNCMER: pos=" << 0 << " kmer=" << kmer << " hash="
            // << reverseKmerHash << " orientation=reverse" << std::endl;
          }
        } else {
          if (returnAll)
            syncmers.emplace_back(std::make_tuple(max_size_t, true, false, 0));
        }
      }
    } else {
      if (returnAll)
        syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
    }
  }

  // start rolling
  for (size_t i = k; i < seq.size(); ++i) {
    if (chash(seq[i]) == 0)
      recentAmbiguousBaseIndex = i;
    curSmerHashes.pop_front();
    --curMinSmerHashIndex;
    if (curMinSmerHashIndex < 0) {
      curMinSmerHash = max_size_t;
      for (size_t j = 0; j < curSmerHashes.size(); ++j) {
        if (curSmerHashes[j] < curMinSmerHash) {
          curMinSmerHash = curSmerHashes[j];
          curMinSmerHashIndex = j;
        }
      }
    }
    forwardKmerHash =
        rol(forwardKmerHash, 1) ^ rol(chash(seq[i - k]), k) ^ chash(seq[i]);
    reverseKmerHash = ror(reverseKmerHash, 1) ^
                      ror(chash(comp(seq[i - k])), 1) ^
                      rol(chash(comp(seq[i])), k - 1);
    forwardSmerHash =
        rol(forwardSmerHash, 1) ^ rol(chash(seq[i - s]), s) ^ chash(seq[i]);
    reverseSmerHash = ror(reverseSmerHash, 1) ^
                      ror(chash(comp(seq[i - s])), 1) ^
                      rol(chash(comp(seq[i])), s - 1);

    if (forwardSmerHash < reverseSmerHash) {
      curSmerHashes.push_back(forwardSmerHash);
      if (forwardSmerHash < curMinSmerHash) {
        curMinSmerHash = forwardSmerHash;
        curMinSmerHashIndex = k - s;
      }
    } else if (reverseSmerHash < forwardSmerHash) {
      curSmerHashes.push_back(reverseSmerHash);
      if (reverseSmerHash < curMinSmerHash) {
        curMinSmerHash = reverseSmerHash;
        curMinSmerHashIndex = k - s;
      }
    } else {
      curSmerHashes.push_back(max_size_t);
    }

    if (recentAmbiguousBaseIndex >= 0 && i < recentAmbiguousBaseIndex + k) {
      if (returnAll)
        syncmers.emplace_back(
            std::make_tuple(max_size_t, false, false, i - k + 1));
    } else {
      if (forwardKmerHash < reverseKmerHash) {
        if (open) {
          if (curMinSmerHash == curSmerHashes[t]) {
            syncmers.emplace_back(
                std::make_tuple(forwardKmerHash, false, true, i - k + 1));
            if (printKmers) {
              std::string kmer = std::string(seq.substr(i - k + 1, k));
              // std::cout << "SYNCMER: pos=" << i-k+1 << " kmer=" << kmer << "
              // hash=" << forwardKmerHash << " orientation=forward" <<
              // std::endl;
            }
          } else {
            if (returnAll)
              syncmers.emplace_back(
                  std::make_tuple(max_size_t, false, false, i - k + 1));
          }
        } else {
          if (curMinSmerHash == curSmerHashes[t] ||
              curMinSmerHash == curSmerHashes[k - s - t]) {
            syncmers.emplace_back(
                std::make_tuple(forwardKmerHash, false, true, i - k + 1));
            if (printKmers) {
              std::string kmer = std::string(seq.substr(i - k + 1, k));
              // std::cout << "SYNCMER: pos=" << i-k+1 << " kmer=" << kmer << "
              // hash=" << forwardKmerHash << " orientation=forward" <<
              // std::endl;
            }
          } else {
            if (returnAll)
              syncmers.emplace_back(
                  std::make_tuple(max_size_t, false, false, i - k + 1));
          }
        }
      } else if (reverseKmerHash < forwardKmerHash) {
        if (open) {
          if (curMinSmerHash == curSmerHashes[t]) {
            syncmers.emplace_back(
                std::make_tuple(reverseKmerHash, true, true, i - k + 1));
            if (printKmers) {
              std::string kmer = std::string(seq.substr(i - k + 1, k));
              // std::cout << "SYNCMER: pos=" << i-k+1 << " kmer=" << kmer << "
              // hash=" << reverseKmerHash << " orientation=reverse" <<
              // std::endl;
            }
          } else {
            if (returnAll)
              syncmers.emplace_back(
                  std::make_tuple(max_size_t, true, false, i - k + 1));
          }
        } else {
          if (curMinSmerHash == curSmerHashes[t] ||
              curMinSmerHash == curSmerHashes[k - s - t]) {
            syncmers.emplace_back(
                std::make_tuple(reverseKmerHash, true, true, i - k + 1));
            if (printKmers) {
              std::string kmer = std::string(seq.substr(i - k + 1, k));
              // std::cout << "SYNCMER: pos=" << i-k+1 << " kmer=" << kmer << "
              // hash=" << reverseKmerHash << " orientation=reverse" <<
              // std::endl;
            }
          } else {
            if (returnAll)
              syncmers.emplace_back(
                  std::make_tuple(max_size_t, true, false, i - k + 1));
          }
        }
      } else {
        if (returnAll)
          syncmers.emplace_back(
              std::make_tuple(max_size_t, false, false, i - k + 1));
      }
    }
  }

  return syncmers;
}

void seedsFromFastq(
    const int32_t &k, const int32_t &s, const int32_t &t, const bool &open,
    const int32_t &l,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::vector<std::vector<seed>> &readSeeds, const std::string &fastqPath1,
    const std::string &fastqPath2);

inline std::string getNextSyncmer(std::string &seq, const int32_t currPos,
                                  const int32_t k, const int32_t s) {

  for (int32_t i = currPos; i < seq.size() - k + 1; i++) {
    std::string kmer = seq.substr(i, k);
    if (is_syncmer(kmer, s, false)) {
      return kmer;
    }
  }
  return "";
}

inline std::string reverseComplement(std::string dna_sequence) {
  std::string complement = "";
  for (char c : dna_sequence) {
    switch (c) {
    case 'A':
      complement += 'T';
      break;
    case 'T':
      complement += 'A';
      break;
    case 'C':
      complement += 'G';
      break;
    case 'G':
      complement += 'C';
      break;
    default:
      complement += c;
      break;
    }
  }
  std::reverse(complement.begin(), complement.end());
  return complement;
}
// inline std::vector<seed> syncmerize(const std::string &seq, const int32_t k,
// const int32_t s, const bool open, const bool aligned, const int32_t pad)
// {
//     std::mutex mtx;
//     std::vector<seed> ret;
//     int32_t seqLen = seq.size();
//     if (seqLen < k)
//     {
//         return ret;
//     }
//     if (aligned)
//     {
//         std::unordered_map<int32_t, int32_t> degap;
//         int32_t pos = 0;
//         std::string ungapped = "";
//         for (int32_t i = 0; i < seqLen; i++)
//         {
//             char c = seq[i];
//             degap[pos] = i;
//             if (c != '-')
//             {
//                 ungapped += c;
//                 pos++;
//             }
//         }

//         if (ungapped.size() < k + 1)
//         {
//             return ret;
//         }
//         for (int32_t i = 0; i < ungapped.size() - k + 1; i++)
//         {
//             std::string kmer = ungapped.substr(i, k);
//             if (is_syncmer(kmer, s, open))
//             {
//                 ret.push_back(seed{kmer, degap[i] + pad, -1, false, degap[i +
//                 k - 1] + pad});
//             }
//         }
//     }
//     else
//     {

//         // Loop through sequence and build up seeds as we go
//         if (seq.size() >= k)
//         {

//             // Check first k bp for smers
//             std::string min_s = seq.substr(0, s);
//             int first_min_coord = 0;
//             int last_min_coord = 0;
//             int Ns = seq[0] == 'N';

//             for (int i = 1; i < k - s + 1; i++)
//             {
//                 std::string smer = seq.substr(i, s);
//                 if (seq[i] == 'N')
//                     Ns++;

//                 if (smer < min_s)
//                 {
//                     min_s = smer;
//                     first_min_coord = i;
//                     last_min_coord = i;
//                 }
//                 else if (smer == min_s)
//                 {
//                     last_min_coord = i;
//                 }
//             }
//             for (int i = k - s + 1; i < k; i++)
//             {
//                 if (seq[i] == 'N')
//                     Ns++;
//             }

//             // Check the rest for smers and kmer seeds
//             for (int i = k; i < seq.size() + 1; i++)
//             {

//                 // Processing kmer starting at i - k in seq
//                 bool isSeed = (first_min_coord == i - k || last_min_coord ==
//                 i - s) && Ns <= k / 2;

//                 if (isSeed)
//                 {
//                     std::string kmer = seq.substr(i - k, k);
//                     ret.push_back(seed{kmer, i - k, -1, false, i});
//                 }

//                 // Updating smers for kmer starting at i - k + 1
//                 if (i < seq.size())
//                 {
//                     if (seq[i] == 'N')
//                         Ns++;
//                     if (seq[i - k] == 'N')
//                         Ns--;

//                     if (first_min_coord == i - k)
//                     {
//                         // Were losing the lowest smer, Re search for lowest
//                         min_s = seq.substr(i - k + 1, s);
//                         first_min_coord = i - k + 1;

//                         for (int j = 1; j < k - s + 1; j++)
//                         {
//                             std::string smer = seq.substr(i - k + 1 + j, s);
//                             if (smer < min_s)
//                             {
//                                 min_s = smer;
//                                 first_min_coord = i - k + 1 + j;
//                                 last_min_coord = first_min_coord;
//                             }
//                             else if (smer == min_s)
//                             {
//                                 last_min_coord = i - k + 1 + j;
//                             }
//                         }
//                     }
//                     else
//                     {
//                         // Test new smer to see if its the lowest
//                         std::string smer = seq.substr(i - s + 1, s);
//                         if (smer < min_s)
//                         {
//                             min_s = smer;
//                             first_min_coord = i - s + 1;
//                             last_min_coord = i - s + 1;
//                         }
//                         else if (smer == min_s)
//                         {
//                             last_min_coord = i - s + 1;
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return ret;
// }

// SIMD-accelerated hash calculation for multiple k-mers at once
#ifdef __AVX2__
#include <immintrin.h>
#endif

/**
 * @brief Calculate hashes for multiple k-mers in a batch
 * @param buffer Array of packed k-mer sequences
 * @param k K-mer length
 * @param batch_size Number of k-mers to process
 * @param fwdHashes Output array for forward strand hashes
 * @param revHashes Output array for reverse strand hashes
 * @param stride Distance between consecutive k-mers in buffer
 */
inline void batchHashKmers(
    const char *buffer, // Input: array of k-mers (k*batch_size chars)
    size_t k,           // Input: k-mer length
    size_t batch_size,  // Input: number of k-mers to process
    size_t *fwdHashes,  // Output: forward hashes
    size_t *revHashes,  // Output: reverse hashes
    size_t stride = 0   // Stride between k-mers (default: k)
) {
  if (stride == 0)
    stride = k;

  // Process in smaller chunks for better cache behavior
  constexpr size_t CHUNK_SIZE = 16; // Process 16 k-mers at a time

#pragma omp parallel for
  for (size_t chunk = 0; chunk < batch_size; chunk += CHUNK_SIZE) {
    size_t end = std::min(chunk + CHUNK_SIZE, batch_size);

    for (size_t i = chunk; i < end; i++) {
      const char *kmer = buffer + i * stride;
      size_t fHash = 0;
      size_t rHash = 0;

      // Prefetch next k-mer
      if (i + 1 < end) {
        __builtin_prefetch(buffer + (i + 1) * stride, 0, 3);
      }

#ifdef __AVX2__
      // Use AVX2 intrinsics when available for processing 4 bases at once
      for (size_t j = 0; j < k; j += 4) {
        // Process 4 nucleotides at a time when possible
        if (j + 4 <= k) {
          // Load 4 bases
          __m128i bases = _mm_loadu_si32(kmer + j); // Load 4 bytes

          // Process each base
          for (size_t m = 0; m < 4; m++) {
            if (j + m >= k)
              break;

            // Compute hash for current base
            char base = kmer[j + m];
            size_t cval;
            switch (base) {
            case 'A':
            case 'a':
              cval = 0x3c8bfbb395c60474;
              break;
            case 'C':
            case 'c':
              cval = 0x3193c18562a02b4c;
              break;
            case 'G':
            case 'g':
              cval = 0x20323ed082572324;
              break;
            case 'T':
            case 't':
              cval = 0x295549f54be24456;
              break;
            default:
              cval = 0;
              throw std::invalid_argument("Kmer contains non canonical base");
            }

            fHash ^= rol(cval, k - (j + m) - 1);

            // Reverse complement hash
            char revBase = comp(kmer[k - (j + m) - 1]);
            size_t rval;
            switch (revBase) {
            case 'A':
            case 'a':
              rval = 0x3c8bfbb395c60474;
              break;
            case 'C':
            case 'c':
              rval = 0x3193c18562a02b4c;
              break;
            case 'G':
            case 'g':
              rval = 0x20323ed082572324;
              break;
            case 'T':
            case 't':
              rval = 0x295549f54be24456;
              break;
            default:
              rval = 0;
              throw std::invalid_argument("Kmer contains non canonical base");
            }

            rHash ^= rol(rval, k - (j + m) - 1);
          }
        } else {
          // Handle remaining bases (< 4)
          for (size_t m = 0; j + m < k; m++) {
            char base = kmer[j + m];
            size_t cval;
            switch (base) {
            case 'A':
            case 'a':
              cval = 0x3c8bfbb395c60474;
              break;
            case 'C':
            case 'c':
              cval = 0x3193c18562a02b4c;
              break;
            case 'G':
            case 'g':
              cval = 0x20323ed082572324;
              break;
            case 'T':
            case 't':
              cval = 0x295549f54be24456;
              break;
            default:
              cval = 0;
              throw std::invalid_argument("Kmer contains non canonical base");
            }

            fHash ^= rol(cval, k - (j + m) - 1);

            // Reverse complement hash
            char revBase = comp(kmer[k - (j + m) - 1]);
            size_t rval;
            switch (revBase) {
            case 'A':
            case 'a':
              rval = 0x3c8bfbb395c60474;
              break;
            case 'C':
            case 'c':
              rval = 0x3193c18562a02b4c;
              break;
            case 'G':
            case 'g':
              rval = 0x20323ed082572324;
              break;
            case 'T':
            case 't':
              rval = 0x295549f54be24456;
              break;
            default:
              rval = 0;
              throw std::invalid_argument("Kmer contains non canonical base");
            }

            rHash ^= rol(rval, k - (j + m) - 1);
          }
        }
      }
#else
      // Scalar implementation as fallback
      for (size_t j = 0; j < k; j++) {
        char base = kmer[j];
        size_t cval;
        switch (base) {
        case 'A':
        case 'a':
          cval = 0x3c8bfbb395c60474;
          break;
        case 'C':
        case 'c':
          cval = 0x3193c18562a02b4c;
          break;
        case 'G':
        case 'g':
          cval = 0x20323ed082572324;
          break;
        case 'T':
        case 't':
          cval = 0x295549f54be24456;
          break;
        default:
          cval = 0;
          throw std::invalid_argument("Kmer contains non canonical base");
        }

        fHash ^= rol(cval, k - j - 1);

        // Calculate reverse complement hash in parallel
        char revBase = comp(kmer[k - j - 1]);
        size_t rval;
        switch (revBase) {
        case 'A':
        case 'a':
          rval = 0x3c8bfbb395c60474;
          break;
        case 'C':
        case 'c':
          rval = 0x3193c18562a02b4c;
          break;
        case 'G':
        case 'g':
          rval = 0x20323ed082572324;
          break;
        case 'T':
        case 't':
          rval = 0x295549f54be24456;
          break;
        default:
          rval = 0;
          throw std::invalid_argument("Kmer contains non canonical base");
        }

        rHash ^= rol(rval, k - j - 1);
      }
#endif

      fwdHashes[i] = fHash;
      revHashes[i] = rHash;
    }
  }
}

/**
 * @brief Calculate canonical hash (min of forward/reverse) for multiple k-mers
 * @param buffer Array of k-mers
 * @param k K-mer length
 * @param batch_size Number of k-mers to process
 * @param hashes Output buffer for canonical hashes
 * @param orientations Output buffer for orientations (true=reverse)
 * @param stride Distance between consecutive k-mers
 */
inline void
batchHashAndSelect(const char *buffer, // Input: array of k-mers
                   size_t k,           // Input: k-mer length
                   size_t batch_size,  // Input: number of k-mers to process
                   size_t *hashes,     // Output: canonical hashes
                   bool *orientations, // Output: orientations
                   size_t stride = 0   // Stride between k-mers
) {
  // Allocate temporary arrays for both hash types
  alignas(64) static thread_local std::vector<size_t> fwdHashes;
  alignas(64) static thread_local std::vector<size_t> revHashes;

  if (fwdHashes.size() < batch_size) {
    fwdHashes.resize(batch_size);
    revHashes.resize(batch_size);
  }

  // Calculate both forward and reverse hashes
  batchHashKmers(buffer, k, batch_size, fwdHashes.data(), revHashes.data(),
                 stride);

  // Select canonical hash (minimum of forward and reverse)
  for (size_t i = 0; i < batch_size; i++) {
    if (fwdHashes[i] <= revHashes[i]) {
      hashes[i] = fwdHashes[i];
      orientations[i] = false;
    } else {
      hashes[i] = revHashes[i];
      orientations[i] = true;
    }
  }
}

// Optimized k-mer batch hashing using AVX2
#ifdef __AVX2__
inline void batchHashKmersAVX2(const char *kmers,    // Input array of k-mers
                               size_t k,             // K-mer length
                               size_t batchSize,     // Number of k-mers
                               size_t *resultHashes, // Output hashes
                               bool *resultIsReverse // Output orientations
) {
  // Base nucleotide hash constants
  const size_t A_HASH = 0x3c8bfbb395c60474;
  const size_t C_HASH = 0x3193c18562a02b4c;
  const size_t G_HASH = 0x20323ed082572324;
  const size_t T_HASH = 0x295549f54be24456;

  // Process 4 k-mers at a time
  for (size_t i = 0; i < batchSize; i += 4) {
    size_t currentBatchSize = std::min(size_t(4), batchSize - i);
    alignas(32) char revComps[4][64]; // Buffer for reverse complements

    // Step 1: Create reverse complements in SIMD using AVX2
    for (size_t j = 0; j < currentBatchSize; j++) {
      const char *kmer = kmers + (i + j) * k;

      // Load a k-mer
      for (size_t n = 0; n < k; n++) {
        // Reverse complement using AVX2
        __m256i base = _mm256_set1_epi8(kmer[k - n - 1]);

        // Create masks for each base
        __m256i isA = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('A'));
        __m256i isa = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('a'));
        __m256i isC = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('C'));
        __m256i isc = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('c'));
        __m256i isG = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('G'));
        __m256i isg = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('g'));
        __m256i isT = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('T'));
        __m256i ist = _mm256_cmpeq_epi8(base, _mm256_set1_epi8('t'));

        // Combine A and a, C and c, etc.
        __m256i isAny_A = _mm256_or_si256(isA, isa);
        __m256i isAny_C = _mm256_or_si256(isC, isc);
        __m256i isAny_G = _mm256_or_si256(isG, isg);
        __m256i isAny_T = _mm256_or_si256(isT, ist);

        // Select the complementary bases
        __m256i result = _mm256_setzero_si256();

        // A -> T
        result = _mm256_blendv_epi8(result, _mm256_set1_epi8('T'), isAny_A);

        // C -> G
        result = _mm256_blendv_epi8(result, _mm256_set1_epi8('G'), isAny_C);

        // G -> C
        result = _mm256_blendv_epi8(result, _mm256_set1_epi8('C'), isAny_G);

        // T -> A
        result = _mm256_blendv_epi8(result, _mm256_set1_epi8('A'), isAny_T);

        // Store the result (just 1 byte)
        revComps[j][n] = _mm256_extract_epi8(result, 0);
      }
    }

    // Step 2: Compute hashes for both forward and reverse
    // The hashing is complex with rolling shifts, using scalar code for
    // correctness
    for (size_t j = 0; j < currentBatchSize; j++) {
      // Get the current k-mer and its reverse complement
      const char *kmer = kmers + (i + j) * k;
      const char *revComp = revComps[j];

      // Compute forward hash
      size_t fwdHash = 0;
      for (size_t n = 0; n < k; n++) {
        size_t baseHash;
        switch (kmer[n]) {
        case 'A':
        case 'a':
          baseHash = A_HASH;
          break;
        case 'C':
        case 'c':
          baseHash = C_HASH;
          break;
        case 'G':
        case 'g':
          baseHash = G_HASH;
          break;
        case 'T':
        case 't':
          baseHash = T_HASH;
          break;
        default:
          baseHash = 0;
          break;
        }
        fwdHash ^= rol(baseHash, k - n - 1);
      }

      // Compute reverse hash
      size_t revHash = 0;
      for (size_t n = 0; n < k; n++) {
        size_t baseHash;
        switch (revComp[n]) {
        case 'A':
        case 'a':
          baseHash = A_HASH;
          break;
        case 'C':
        case 'c':
          baseHash = C_HASH;
          break;
        case 'G':
        case 'g':
          baseHash = G_HASH;
          break;
        case 'T':
        case 't':
          baseHash = T_HASH;
          break;
        default:
          baseHash = 0;
          break;
        }
        revHash ^= rol(baseHash, k - n - 1);
      }

      // Store canonical hash and orientation
      if (fwdHash <= revHash) {
        resultHashes[i + j] = fwdHash;
        resultIsReverse[i + j] = false;
      } else {
        resultHashes[i + j] = revHash;
        resultIsReverse[i + j] = true;
      }
    }
  }
}
#endif

} // namespace seeding

#endif
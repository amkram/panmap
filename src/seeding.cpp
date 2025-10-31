#include "seeding.hpp"
#include "performance.hpp"  // For memory::Buffer and memory::BufferPool
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <cctype>
#include <absl/container/flat_hash_map.h>

namespace seeding {

// Updated buffer management support
class SeedingBufferGuard {
private:
  memory::Buffer *buffer;
  memory::BufferPool &pool;

public:
  SeedingBufferGuard(memory::BufferPool &bufferPool, size_t size)
      : pool(bufferPool) {
    buffer = pool.acquire(size);
  }

  ~SeedingBufferGuard() {
    if (buffer) {
      pool.release(buffer);
    }
  }

  std::vector<char> &getSeqBuffer() { return buffer->seqBuffer; }
  std::vector<int> &getGapsBuffer() { return buffer->gapsBuffer; }
  std::vector<int> &getCoordsBuffer() { return buffer->coordsBuffer; }
  std::vector<int> &getDeadBlocksBuffer() { return buffer->deadBlocksBuffer; }
};

void perfect_shuffle(std::vector<std::string> &v) {
  
  int n = v.size();

  std::vector<std::string> canvas(n);

  for (int i = 0; i < n / 2; i++) {
    canvas[i * 2] = v[i];
    canvas[i * 2 + 1] = v[i + n / 2];
  }

  v = std::move(canvas);
}

size_t chash(const char& c) {
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

bool is_syncmer(const std::string &seq, const int s, const bool open) {
  int NsCount = 0;
  for (size_t i = 0; i < seq.size(); i++) {
    if (seq[i] == 'N') {
      NsCount++;
    }
  }
  if (NsCount > 1) {
    return false;  
  }
  if (static_cast<size_t>(s) > seq.size()) {
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

std::tuple<size_t, bool, bool> is_syncmer_rollingHash(const std::string &seq, const int s, const bool open, int t) {
  if (seq.size() < s) return std::make_tuple(0, false, false);

  int k = seq.size();
  size_t forwardKmerHash = 0, reverseKmerHash = 0;
  size_t forwardSmerHash = 0, reverseSmerHash = 0;

  size_t minSmerHash    = std::numeric_limits<size_t>::max();
  size_t tthSmerHash    = std::numeric_limits<size_t>::max();
  size_t lastthSmerHash = std::numeric_limits<size_t>::min();

  for (int i = 0; i < s; ++i) {
    if (chash(seq[i]) == 0) return std::make_tuple(0, false, false);

    forwardKmerHash ^= rol(chash(seq[i]), k-i-1);
    reverseKmerHash ^= rol(chash(comp(seq[k-i-1])), k-i-1);
    forwardSmerHash ^= rol(chash(seq[i]), s-i-1);
    reverseSmerHash ^= rol(chash(comp(seq[s-i-1])), s-i-1);
  }

  if (forwardSmerHash < reverseSmerHash) {
    minSmerHash = forwardSmerHash;
    if (t == 0) tthSmerHash = forwardSmerHash;
  } else if (reverseSmerHash < forwardSmerHash) {
    minSmerHash = reverseSmerHash;
    if (t == 0) tthSmerHash = reverseSmerHash;
  }

  for (int i = 1; i < k-s+1; ++i) {
    if (chash(seq[i+s-1]) == 0) return std::make_tuple(0, false, false);

    forwardKmerHash ^= rol(chash(seq[i+s-1]), k-s-i);
    reverseKmerHash ^= rol(chash(comp(seq[k-s-i])), k-s-i);
    forwardSmerHash = rol(forwardSmerHash, 1) ^ rol(chash(seq[i-1]), s) ^ chash(seq[i+s-1]);
    reverseSmerHash = ror(reverseSmerHash, 1) ^ ror(chash(comp(seq[i-1])), 1) ^ rol(chash(comp(seq[i+s-1])), s-1);
    
    if (i == t) {
      if      (forwardSmerHash < reverseSmerHash) tthSmerHash = forwardSmerHash;
      else if (reverseSmerHash < forwardSmerHash) tthSmerHash = reverseSmerHash;
    } else if (i == k-s-t) {
      if      (forwardSmerHash < reverseSmerHash) lastthSmerHash = forwardSmerHash;
      else if (reverseSmerHash < forwardSmerHash) lastthSmerHash = reverseSmerHash;
    }
    
    if      (forwardSmerHash < reverseSmerHash && forwardSmerHash < minSmerHash) minSmerHash = forwardSmerHash;
    else if (reverseSmerHash < forwardSmerHash && reverseSmerHash < minSmerHash) minSmerHash = reverseSmerHash;
  }

  if (minSmerHash == std::numeric_limits<size_t>::max()) return std::make_tuple(0, false, false);
  
  if (forwardKmerHash < reverseKmerHash) {
    if (open) {
      if (minSmerHash == tthSmerHash) return std::make_tuple(forwardKmerHash, false, true);
    } else {
      if (minSmerHash == tthSmerHash || minSmerHash == lastthSmerHash) return std::make_tuple(forwardKmerHash, false, true);
    }
  } else if (reverseKmerHash < forwardKmerHash) {
    if (open) {
      if (minSmerHash == tthSmerHash) return std::make_tuple(reverseKmerHash, true, true);
    } else {
      if (minSmerHash == tthSmerHash || minSmerHash == lastthSmerHash) return std::make_tuple(reverseKmerHash, true, true);
    }
  }

  return std::make_tuple(0, false, false);
}

std::pair<size_t, size_t> hashSeq(const std::string& s) {
  size_t fHash = 0;
  size_t rHash = 0;
  int k = s.size();
  for (int i = 0; i < k; i++) {
    if (chash(s[i]) == 0) throw std::invalid_argument("Kmer contains non canonical base");
    fHash ^= rol(chash(s[i]), k-i-1);
    rHash ^= rol(chash(comp(s[k-i-1])), k-i-1);
  }
  return std::make_pair(fHash, rHash);
}

std::vector<std::tuple<size_t, bool, bool, int64_t>> rollingSyncmers(std::string_view seq, int k, int s, bool open, int t, bool returnAll) {
  std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmers;
  if (seq.size() < k) return syncmers;
  syncmers.reserve(seq.size() - k + 1);

  const size_t max_size_t = std::numeric_limits<size_t>::max();
  size_t forwardKmerHash = 0, reverseKmerHash = 0;
  size_t forwardSmerHash = 0, reverseSmerHash = 0;

  int recentAmbiguousBaseIndex = std::numeric_limits<int>::min();

  size_t curMinSmerHash = max_size_t;
  int curMinSmerHashIndex = -1;
  std::deque<size_t> curSmerHashes;

  // first smer
  for (int i = 0; i < s; ++i) {
    if (chash(seq[i]) == 0) recentAmbiguousBaseIndex = i;

    forwardKmerHash ^= rol(chash(seq[i]), k-i-1);
    reverseKmerHash ^= rol(chash(comp(seq[k-i-1])), k-i-1);
    forwardSmerHash ^= rol(chash(seq[i]), s-i-1);
    reverseSmerHash ^= rol(chash(comp(seq[s-i-1])), s-i-1);
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
    if (chash(seq[i]) == 0) recentAmbiguousBaseIndex = i;

    forwardKmerHash ^= rol(chash(seq[i]), k-i-1);
    reverseKmerHash ^= rol(chash(comp(seq[k-i-1])), k-i-1);
    forwardSmerHash  = rol(forwardSmerHash, 1) ^ rol(chash(seq[i-s]), s)       ^ chash(seq[i]);
    reverseSmerHash  = ror(reverseSmerHash, 1) ^ ror(chash(comp(seq[i-s])), 1) ^ rol(chash(comp(seq[i])), s-1);

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
    if (returnAll) syncmers.emplace_back(std::make_tuple(0, false, false, 0));
  } else {
    if (forwardKmerHash < reverseKmerHash) {
      if (open) {
        if (curMinSmerHash == curSmerHashes[t]) {
          syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, 0));
        } else {
          if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
        }
      } else {
        if (curMinSmerHash == curSmerHashes[t] || curMinSmerHash == curSmerHashes[k-s-t]) {
          syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, 0));
        } else {
          if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
        }
      }
    } else if (reverseKmerHash < forwardKmerHash) {
      if (open) {
        if (curMinSmerHash == curSmerHashes[t]) {
          syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, 0));
        } else {
          if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, true, false, 0));
        }
      } else {
        if (curMinSmerHash == curSmerHashes[t] || curMinSmerHash == curSmerHashes[k-s-t]) {
          syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, 0));
        } else {
          if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, true, false, 0));
        }
      }
    } else {
      if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
    }
  }

  // start rolling
  for (size_t i = k; i < seq.size(); ++i) {
    if (chash(seq[i]) == 0) recentAmbiguousBaseIndex = i;
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
    forwardKmerHash = rol(forwardKmerHash, 1) ^ rol(chash(seq[i-k]), k) ^ chash(seq[i]);
    reverseKmerHash = ror(reverseKmerHash, 1) ^ ror(chash(comp(seq[i-k])), 1) ^ rol(chash(comp(seq[i])), k-1);
    forwardSmerHash = rol(forwardSmerHash, 1) ^ rol(chash(seq[i-s]), s) ^ chash(seq[i]);
    reverseSmerHash = ror(reverseSmerHash, 1) ^ ror(chash(comp(seq[i-s])), 1) ^ rol(chash(comp(seq[i])), s-1);

    if (forwardSmerHash < reverseSmerHash) {
      curSmerHashes.push_back(forwardSmerHash);
      if (forwardSmerHash < curMinSmerHash) {
        curMinSmerHash = forwardSmerHash;
        curMinSmerHashIndex = k-s;
      }
    } else if (reverseSmerHash < forwardSmerHash) {
      curSmerHashes.push_back(reverseSmerHash);
      if (reverseSmerHash < curMinSmerHash) {
        curMinSmerHash = reverseSmerHash;
        curMinSmerHashIndex = k-s;
      }
    } else {
      curSmerHashes.push_back(max_size_t);
    }

    if (recentAmbiguousBaseIndex >= 0 && i < recentAmbiguousBaseIndex + k) {
      if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
    } else {
      if (forwardKmerHash < reverseKmerHash) {
        if (open) {
          if (curMinSmerHash == curSmerHashes[t]) {
            syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, i - k + 1));
          } else {
            if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
          }
        } else {
          if (curMinSmerHash == curSmerHashes[t] || curMinSmerHash == curSmerHashes[k-s-t]) {
            syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, i - k + 1));
          } else {
            if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
          }
        }
      } else if (reverseKmerHash < forwardKmerHash) {
        if (open) {
          if (curMinSmerHash == curSmerHashes[t]) {
            syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, i - k + 1));
          } else {
            if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, true, false, i - k + 1));
          }
        } else {
          if (curMinSmerHash == curSmerHashes[t] || curMinSmerHash == curSmerHashes[k-s-t]) {
            syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, i - k + 1));
          } else {
            if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, true, false, i - k + 1));
          }
        }
      } else {
        if (returnAll) syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
      }
    }
  }
  syncmers.shrink_to_fit();
  return syncmers;
}



void seedsFromFastq(
    const int32_t &k, const int32_t &s, const int32_t &t, const bool &open,
    const int32_t &l,
    absl::flat_hash_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::vector<std::vector<seeding::seed_t>> &readSeeds,
    std::vector<std::vector<std::string>> &readSeedSeqs,
    const std::string &fastqPath1,
    const std::string &fastqPath2) {
  FILE *fp;
  kseq_t *seq;
  fp = fopen(fastqPath1.c_str(), "r");
  if (!fp) {
    std::cerr << "Error: File " << fastqPath1 << " not found" << std::endl;
    exit(0);
  }
  seq = kseq_init(fileno(fp));
  int line;
  while ((line = kseq_read(seq)) >= 0) {
    readSequences.push_back(seq->seq.s);
    readNames.push_back(seq->name.s);
    readQuals.push_back(seq->qual.s);
  }
  if (fastqPath2.size() > 0) {
    fp = fopen(fastqPath2.c_str(), "r");
    if (!fp) {
      std::cerr << "Error: File " << fastqPath2 << " not found" << std::endl;
      exit(0);
    }
    seq = kseq_init(fileno(fp));

    line = 0;
    int forwardReads = readSequences.size();
    while ((line = kseq_read(seq)) >= 0) {
      readSequences.push_back(seq->seq.s);
      readNames.push_back(seq->name.s);
      readQuals.push_back(seq->qual.s);
    }

    // Shuffle reads together, so that pairs are next to each other
    perfect_shuffle(readSequences);
    perfect_shuffle(readNames);
    perfect_shuffle(readQuals);
  }

  for (int i = 0; i < readSequences.size(); i++) {
    std::vector<seeding::seed_t> curReadSeeds;
    std::vector<std::string> curReadSeedSeqs;
    std::string readSeq = readSequences[i];
    for (const auto &[kmerHash, isReverse, isSyncmer, startPos] :
         rollingSyncmers(readSeq, k, s, false, t, false)) {
      if (!isSyncmer)
        continue;
      std::string kmer = readSeq.substr(startPos, k);
      // Create seed with correct field order and add endPos
      // Order: hash, pos, idx, reversed, rpos, endPos
      curReadSeeds.emplace_back(
          seed_t{
            kmerHash,                       // hash
            static_cast<int64_t>(startPos), // pos 
            -1,                             // idx
            isReverse,                      // reversed
            0,                              // rpos
            static_cast<int64_t>(startPos + k - 1) // endPos
          });
      curReadSeedSeqs.push_back(kmer);
      if (readSeedCounts.find(kmerHash) == readSeedCounts.end())
        readSeedCounts[kmerHash] = std::make_pair(0, 0);
      if (isReverse)
        ++readSeedCounts[kmerHash].second;
      else
        ++readSeedCounts[kmerHash].first;
    }
    readSeeds.push_back(std::move(curReadSeeds));
    readSeedSeqs.push_back(std::move(curReadSeedSeqs));
  }
}

void recalculateReadSeeds(
  int32_t k, int32_t s, bool open, int32_t t,
  const std::vector<std::string> &readSequences,
  std::vector<std::vector<seed_t>> &readSeeds
){
  readSeeds.clear();
  readSeeds.resize(readSequences.size());
  for (int i = 0; i < readSequences.size(); i++) {
    std::vector<seeding::seed_t> curReadSeeds;
    const std::string &readSeq = readSequences[i];
    for (const auto &[kmerHash, isReverse, isSyncmer, startPos] : rollingSyncmers(readSeq, k, s, open, t, false)) {
        if (!isSyncmer) continue;
        curReadSeeds.emplace_back(
          seeding::seed_t{
            kmerHash,                       // hash
            static_cast<int64_t>(startPos), // pos 
            -1,                             // idx
            isReverse,                      // reversed
            0,                              // rpos
            static_cast<int64_t>(startPos + k - 1) // endPos
          });
      }
    readSeeds[i] = std::move(curReadSeeds);
  }
}

std::string getNextSyncmer(std::string &seq, const int32_t currPos, const int32_t k, const int32_t s) {
  for (int32_t i = currPos; i < static_cast<int32_t>(seq.size()) - k + 1; i++) {
  std::string kmer = seq.substr(i, k);
  if (is_syncmer(kmer, s, false)) {
  return kmer;
  }
  }
  return "";
}

std::string reverseComplement(std::string dna_sequence) {
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

} // namespace seeding
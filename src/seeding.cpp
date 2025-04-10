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
#include <unordered_map>
#include <utility>
#include <vector>

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

void seedsFromFastq(
    const int32_t &k, const int32_t &s, const int32_t &t, const bool &open,
    const int32_t &l,
    std::unordered_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::vector<std::vector<seeding::seed_t>> &readSeeds, const std::string &fastqPath1,
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
      readSequences.push_back(reverseComplement(seq->seq.s));
      readNames.push_back(seq->name.s);
      readQuals.push_back(seq->qual.s);
    }

    if (readSequences.size() != forwardReads * 2) {
      std::cerr << "Error: File " << fastqPath2
                << " does not contain the same number of reads as "
                << fastqPath1 << std::endl;
      exit(0);
    }

    // Shuffle reads together, so that pairs are next to each other
    perfect_shuffle(readSequences);
    perfect_shuffle(readNames);
    perfect_shuffle(readQuals);
  }

  for (int i = 0; i < readSequences.size(); i++) {
    std::vector<seeding::seed_t> curReadSeeds;
    for (const auto &[kmerHash, isReverse, isSyncmer, startPos] :
         rollingSyncmers(readSequences[i], k, s, open, t, false)) {
      if (!isSyncmer)
        continue;
      
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
        
      if (readSeedCounts.find(kmerHash) == readSeedCounts.end())
        readSeedCounts[kmerHash] = std::make_pair(0, 0);
      if (isReverse)
        ++readSeedCounts[kmerHash].second;
      else
        ++readSeedCounts[kmerHash].first;
    }
    readSeeds.push_back(std::move(curReadSeeds));
  }
}
} // namespace seeding
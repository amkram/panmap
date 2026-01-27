#pragma once

// Include SIMD headers in global namespace
#include <cstdint>
#include <immintrin.h> // For SIMD intrinsics

// Other includes
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
#include <cctype> // Added for std::tolower

#include <htslib/kseq.h>
KSEQ_INIT(int, read)

#include <absl/container/flat_hash_map.h>



namespace seeding
{

// Perfect shuffle for interleaving paired-end reads (R1_0, R2_0, R1_1, R2_1, ...)
template <typename T>
inline void perfect_shuffle(std::vector<T>& v) {
  const size_t n = v.size();
  if (n < 2) return;
  
  std::vector<T> canvas(n);
  for (size_t i = 0; i < n / 2; i++) {
    canvas[i*2] = std::move(v[i]);
    canvas[i*2+1] = std::move(v[i + n/2]);
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

  uniqueKminmer_t(uint64_t startPos, uint64_t endPos, size_t hash, bool isReverse) : startPos(startPos), endPos(endPos), hash(hash), isReverse(isReverse) {}
  uniqueKminmer_t() : startPos(0), endPos(0), hash(0), isReverse(false) {}

  bool operator==(const uniqueKminmer_t& other) const {
    return startPos == other.startPos && endPos == other.endPos && hash == other.hash && isReverse == other.isReverse;
  }
};

// A syncmer seed, defined within a single sequence
struct seed_t {
  // Core properties
  size_t hash;        // Hash of the k-mer
  int64_t pos;        // Position in sequence (keep as pos, not startPos)
  int32_t idx;        // Index in vector (used during indexing)
  bool reversed;      // Orientation flag
  uint32_t rpos;      // Position on reference (used during alignment)
  
  // Add the end position field
  int64_t endPos;     // End position with gaps
  
  // Comparison operators
  bool operator<(const seed_t &rhs) const { return pos < rhs.pos; }
  bool operator==(const seed_t &rhs) const { return pos == rhs.pos; }
};

// Represents a seed match between self and a reference sequence
struct anchor_t : seed_t {
  int64_t referencePos; // position in reference sequence
  int64_t referenceEndPos; // end position in reference sequence
};

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
    // Instead of throwing an exception, return a special value
    // that will be handled by the calling function
    return std::numeric_limits<size_t>::max();
  }
  return n;
}

[[maybe_unused]] static size_t hash(const std::string &s) {
  size_t h = 0;
  if (s.empty()) {
    return h;
  } else if (s.size() == 1) {
    h = btn(s[0]);
    // Check for non-canonical base
    if (h == std::numeric_limits<size_t>::max()) {
      return 0; // Return 0 for non-canonical bases
    }
    return h;
  }

  h = btn(s[0]);
  // Check for non-canonical base in first position
  if (h == std::numeric_limits<size_t>::max()) {
    return 0; // Return 0 for non-canonical bases
  }
  
  for (size_t i = 1; i < s.size(); ++i) {
    size_t base = btn(s[i]);
    // Check for non-canonical base
    if (base == std::numeric_limits<size_t>::max()) {
      return 0; // Return 0 for non-canonical bases
    }
    h = (h << 2) + base;
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


bool is_syncmer(const std::string &seq, const int s, const bool open);

size_t chash(const char& c);

inline size_t rol(const size_t& h, const size_t& r) { return (h << r) | (h >> (64-r)); }

inline size_t ror(const size_t& h, const size_t& r) { return (h >> r) | (h << (64-r)); }

// Hash a sequence, returns (forward hash, reverse complement hash)
std::pair<size_t, size_t> hashSeq(const std::string& s);

// returns vector of (hash, isReverse, isSyncmer, startPos)  
std::vector<std::tuple<size_t, bool, bool, int64_t>> rollingSyncmers(std::string_view seq, int k, int s, bool open, int t = 0, bool returnAll = true);

void seedsFromFastq(
    const int32_t &k, const int32_t &s, const int32_t &t, const bool &open,
    const int32_t &l,
    absl::flat_hash_map<size_t, std::pair<size_t, size_t>> &readSeedCounts,
    std::vector<std::string> &readSequences,
    std::vector<std::string> &readQuals, std::vector<std::string> &readNames,
    std::vector<std::vector<seed_t>> &readSeeds,
    std::vector<std::vector<std::string>> &readSeedSeqs,
    const std::string &fastqPath1,
    const std::string &fastqPath2);

std::string reverseComplement(std::string dna_sequence);

// For tracking seed mutations at each node
struct SeedChange {
    int64_t pos;       // Starting position of the seed, in scalar global genomic coordinates (includes gaps)
    seed_t seed;  // Seed hash information including end position and orientation
    int64_t tritMask;  // Ternary mask for efficient serialization
    
    explicit SeedChange(int64_t position = 0, const seed_t& seedHash = {}, int64_t mask = 0)
        : pos(position), seed(seedHash), tritMask(mask) {}
    
    // Convenience accessor for end position
    int64_t endPos() const { return seed.endPos; }
    
    bool operator==(const SeedChange& other) const {
        return pos == other.pos && 
               seed == other.seed &&
               tritMask == other.tritMask;
    }
};

} // namespace seeding

namespace std {
  template <>
  struct hash<seeding::uniqueKminmer_t> {
    size_t operator()(const seeding::uniqueKminmer_t& kminmer) const {
      // Combine hashes of all member variables, including the new 'startPos'.
      size_t h1 = std::hash<uint64_t>{}(kminmer.startPos); // Hash for startPos
      size_t h2 = std::hash<uint64_t>{}(kminmer.endPos);
      size_t h3 = std::hash<size_t>{}(kminmer.hash);
      size_t h4 = std::hash<bool>{}(kminmer.isReverse);
  
      // A common pattern for combining hashes.
      // Ensure all member hashes are combined.
      size_t seed = 0;
      seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2); // Include h4
      return seed;
    }
  };
}
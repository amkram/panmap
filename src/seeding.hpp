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



namespace {

// Helper to get the complement of a base (case-insensitive)
inline char complementBase(char base) {
    char lowerBase = std::tolower(base);
    switch (lowerBase) {
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default:  return lowerBase; // Return original if not ACGT
    }
}

// Helper to check if a k-mer is a homopolymer (all same bases)
inline bool isHomopolymer(std::string_view kmer) {
    if (kmer.length() < 2) {
        return true; // Single base or empty is trivially a homopolymer
    }
    char firstBaseLower = std::tolower(kmer[0]);
    // Use std::all_of for concise check
    return std::all_of(kmer.begin() + 1, kmer.end(), 
                       [firstBaseLower](char c){ return std::tolower(c) == firstBaseLower; });
}

// Helper to check if a k-mer is its own reverse complement
inline bool isReverseComplementPalindrome(std::string_view kmer) {
    if (kmer.empty()) {
        return true; // Empty string is a palindrome
    }
    // Compare character pairs from beginning and end, checking complement
    auto front = kmer.begin();
    auto back = kmer.rbegin();
    for (size_t i = 0; i < kmer.length() / 2; ++i, ++front, ++back) {
        if (std::tolower(*front) != complementBase(*back)) {
            return false;
        }
    }
    // For odd lengths, the middle character must be its own complement (or non-ACGT)
    if (kmer.length() % 2 != 0) {
        char middle = kmer[kmer.length() / 2];
        if (std::tolower(middle) != complementBase(middle)) {
            // return false; // Optionally uncomment for strict ACGT palindromes
        }
    }
    return true;
}

} // --- End Anonymous Namespace ---

namespace seeding
{

struct rsyncmer_t {
  size_t hash;
  uint64_t endPos;
  bool isReverse;

  rsyncmer_t(size_t hash, uint64_t endPos, bool isReverse) : hash(hash), endPos(endPos), isReverse(isReverse) {}
  rsyncmer_t() : hash(0), endPos(0), isReverse(false) {}
};


struct rkminmer_t {
  size_t hash;
  uint64_t endPos;
  bool isReverse;

  rkminmer_t(size_t hash, uint64_t endPos, bool isReverse) : hash(hash), endPos(endPos), isReverse(isReverse) {}
  rkminmer_t() : hash(0), endPos(0), isReverse(false) {}
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

[[maybe_unused]] static std::string revcomp(const std::string &s) {
  std::string cs = "";
  for (int i = s.size() - 1; i > -1; --i) {
    char c = s[i];
    cs += comp(c);
  }
  return cs;
}


bool is_syncmer(const std::string &seq, const int s, const bool open);

size_t chash(const char& c);

inline size_t rol(const size_t& h, const size_t& r) { return (h << r) | (h >> (64-r)); }

inline size_t ror(const size_t& h, const size_t& r) { return (h >> r) | (h << (64-r)); }


// returns a tuple of (hash, isReverse, isSyncmer)    
std::tuple<size_t, bool, bool> is_syncmer_rollingHash(const std::string &seq, const int s, const bool open, int t = 0);


std::pair<size_t, size_t> hashSeq(const std::string& s);

// returns vector of (hash, isReverse, isSyncmer, startPos)  
std::vector<std::tuple<size_t, bool, bool, int64_t>> rollingSyncmers(const std::string& seq, int k, int s, bool open, int t = 0, bool returnAll = true);

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

void recalculateReadSeeds(
  int32_t k, int32_t s, bool open, int32_t t,
  const std::vector<std::string> &readSequences,
  std::vector<std::vector<seed_t>> &readSeeds
);

std::string getNextSyncmer(std::string &seq, const int32_t currPos, const int32_t k, const int32_t s);

std::string reverseComplement(std::string dna_sequence);

// SIMD-accelerated hash calculation for multiple k-mers at once
#ifdef __AVX2__
#endif

/**
 * @brief Calculate hashes for multiple k-mers in a batch
 * @param buffer Array of packed k-mer sequences
 * @param k K-mer length
 * @param batch_size Number of k-mers to process
 * @param fwdHashes Output array for forward strand hashes
 * @param revHashes Output array for reverse strand hashes
 * @param stride Distance between consecutive k-mers in buffer
 * @param validResults Optional array for validity flags
 */
inline void batchHashKmers(
    const char *buffer, // Input: array of k-mers (k*batch_size chars)
    size_t k,           // Input: k-mer length
    size_t batch_size,  // Input: number of k-mers to process
    size_t *fwdHashes,  // Output: forward hashes
    size_t *revHashes,  // Output: reverse hashes
    size_t stride = 0,   // Stride between k-mers (default: k)
    bool *validResults = nullptr // Optional: array to store validity flags
) {
  if (stride == 0)
    stride = k;

  // Process in smaller chunks for better cache behavior
  constexpr size_t CHUNK_SIZE = 16; // Process 16 k-mers at a time

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t chunk = 0; chunk < batch_size; chunk += CHUNK_SIZE) {
    size_t end = std::min(chunk + CHUNK_SIZE, batch_size);

    for (size_t i = chunk; i < end; i++) {
      const char *kmer = buffer + i * stride;
      size_t fHash = 0;
      size_t rHash = 0;
      bool valid = true;
      
      // Initialize result as valid (if array is provided)
      if (validResults != nullptr) {
        validResults[i] = true;
      }

      // Prefetch next k-mer
      if (i + 1 < end) {
        __builtin_prefetch(buffer + (i + 1) * stride, 0, 3);
      }

#ifdef __AVX2__
      // Use AVX2 intrinsics when available for processing 4 bases at once
      for (size_t j = 0; j < k; j += 4) {
        // Process 4 nucleotides at a time when possible
        if (j + 4 <= k) {
          // Load 4 bases but use directly without storing in variable
          // Comment out the unused variable declaration
          // __m128i bases = _mm_loadu_si32(kmer + j); // Load 4 bytes
          
          // Or use _mm_loadu_si32 directly where needed:
          __m128i bases = _mm_loadu_si32(kmer + j); // Load 4 bytes
          (void)bases; // Mark as used to prevent warning

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
              // Found an ambiguous nucleotide - update the valid flag if provided
              if (validResults != nullptr) {
                validResults[i] = false;
              }
              
              // Skip the rest of this k-mer
              j = k;
              break;
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
              valid = false;
              break;  // Skip this base but try to process others (will be marked invalid)
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
              // Found an ambiguous nucleotide - update the valid flag if provided
              if (validResults != nullptr) {
                validResults[i] = false;
              }
              
              // Skip the rest of this k-mer
              j = k;
              break;
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
              valid = false;
              break;  // Skip this base but try to process others (will be marked invalid)
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
          // Found an ambiguous nucleotide - update the valid flag if provided
          if (validResults != nullptr) {
            validResults[i] = false;
          }
          
          // Skip the rest of this k-mer
          j = k;
          break;
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
          valid = false;
          break;  // Skip this base but try to process others (will be marked invalid)
        }

        rHash ^= rol(rval, k - j - 1);
      }
#endif

      // Store results
      fwdHashes[i] = valid ? fHash : 0;
      revHashes[i] = valid ? rHash : 0;
      
      // Update validity flag if array is provided
      if (validResults != nullptr) {
        validResults[i] = valid;
      }
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
 * @param endPositions Output buffer for end positions (optional)
 * @param startPosition Input start scalar position
 * @param stride Distance between consecutive k-mers
 */
inline void
batchHashAndSelect(const char *buffer, // Input: array of k-mers
                   size_t k,           // Input: k-mer length
                   size_t batch_size,  // Input: number of k-mers to process
                   size_t *hashes,     // Output: canonical hashes
                   bool *orientations, // Output: orientations
                   int64_t *endPositions = nullptr, // Output: end positions (optional)
                   int64_t startPosition = 0,      // Input: start scalar position
                   size_t stride = 0,   // Stride between k-mers
                   bool *validResults = nullptr // Output: validity flags (optional)
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
    // Check if this k-mer is valid by looking for zeros in the hashes
    // (zeros indicate ambiguous nucleotides)
    bool isValid = (fwdHashes[i] != 0 && revHashes[i] != 0);
    
    if (isValid) {
      if (fwdHashes[i] <= revHashes[i]) {
        hashes[i] = fwdHashes[i];
        orientations[i] = false;
      } else {
        hashes[i] = revHashes[i];
        orientations[i] = true;
      }
    } else {
      // For invalid k-mers, set hash to 0 and orientation to false
      hashes[i] = 0;
      orientations[i] = false;
    }
    
    // If validResults is provided, update it
    if (validResults != nullptr) {
      validResults[i] = isValid;
    }
  }

  // If endPositions is provided, calculate them
  if (endPositions != nullptr) {
    for (size_t b = 0; b < batch_size; b++) {
      const char* kmer = buffer + b * (stride ? stride : k);
      
      // Count how many non-gap characters we see
      int nonGapCount = 0;
      int64_t currPos = startPosition;
      
      // Iterate through the k-mer
      for (size_t i = 0; i < k; i++) {
        if (kmer[i] != '-') {
          nonGapCount++;
        }
        currPos++; // Always increment position
      }
      
      // Store the end position
      endPositions[b] = currPos - 1; // -1 because we want the position of the last character
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

// Class to store syncmer seeds and their orientation
class SyncmerSeeds {
private:
    std::vector<size_t> seeds;
    std::vector<bool> isReverse;
    
public:
    // Default constructor
    SyncmerSeeds() = default;
    
    // Initialize from a sequence using rollingSyncmers
    void init(const char* seq, size_t length, int k, int s, bool open) {
        auto syncmers = rollingSyncmers(std::string(seq, length), k, s, open);
        seeds.clear();
        isReverse.clear();
        
        for (const auto& syncmer : syncmers) {
            seeds.push_back(std::get<0>(syncmer));
            isReverse.push_back(std::get<1>(syncmer));
        }
    }
    
    // Get number of seeds
    size_t numSeeds() const {
        return seeds.size();
    }
    
    // Get seed at position i
    size_t getSeedAt(size_t i) const {
        return seeds.at(i);
    }
    
    // Get isReverse flag at position i
    bool getIsReverseAt(size_t i) const {
        return isReverse.at(i);
    }
};

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
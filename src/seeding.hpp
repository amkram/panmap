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
#include <optional>
#include <tuple>
#include <deque>


KSEQ_INIT(int, read)

namespace seeding
{

    struct seed
    {

        std::string seq;
        int32_t pos; // start position in sequence
        int32_t idx; // index of kmer in vector used in indexing DFS
        bool reversed;
        int32_t gappedEnd;
        int32_t rpos; // position on reference (only used during alignment)

        bool operator<(const seed &rhs) const
        {
            return pos < rhs.pos;
        };
        bool operator==(const seed &rhs) const
        {
            return pos == rhs.pos;
        };
    };

    struct seedmer
    {
        int32_t j;
        int32_t k;
        std::string seq;                // concatenated string of s kmers
        std::vector<int32_t> positions; // positions of s kmers
    };

    class KHash
    {
    public:
        size_t operator()(const seed &t) const
        {
            const std::hash<std::string> h;
            return h(t.seq);
        }
    };

    struct read_t
    {
        std::string seq;                    // read sequence
        std::string qual;                   // quality string
        std::vector<int32_t> seedPositions; // positions of seed starts in read
        std::string name;                   // read id
    };

    static size_t btn(char b)
    {
        size_t n;
        switch (b)
        {
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

    static size_t hash(const std::string &s)
    {
        size_t h = 0;
        if (s.empty())
        {
            return h;
        }
        else if (s.size() == 1)
        {
            h = btn(s[0]);
            return h;
        }

        h = btn(s[0]);
        for (size_t i = 1; i < s.size(); ++i)
        {
            h = (h << 2) + btn(s[i]);
        }
        return h;
    }

    static char comp(char c)
    {
        char compC;
        switch (c)
        {
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

    static std::string revcomp(const std::string &s)
    {
        std::string cs = "";
        for (int i = s.size() - 1; i > -1; --i)
        {
            char c = s[i];
            cs += comp(c);
        }
        return cs;
    }

    inline std::pair<size_t, bool> getHash(const std::string &s)
    {
        try
        {
            size_t u = hash(s);
            size_t v = hash(revcomp(s));
            if (u < v)
                return std::make_pair(u, true);
            else if (v < u)
                return std::make_pair(v, true);
            return std::make_pair(0, false);
        }
        catch (std::invalid_argument)
        {
            return std::make_pair(0, false); // skip strand ambiguous
        }
        return std::make_pair(0, false);
    }

    inline bool is_syncmer(const std::string &seq, const int s, const bool open)
    {
        int NsCount = 0;
        for (size_t i = 0; i < seq.size(); i++)
        {
            if (seq[i] == 'N')
            {
                NsCount++;
            }
        }
        if (NsCount > seq.size() / 2)
        {
            return false;
        }
        if (seq.size() < s)
        {
            return false;
        }
        std::string min(s, 'Z');
        for (size_t i = 0; i < seq.size() - s + 1; i++)
        {
            std::string submer = seq.substr(i, s);
            if (submer < min)
            {
                min = submer;
            }
        }
        if (open)
        {
            if (min == seq.substr(0, s))
            {
                return true;
            }
        }
        else
        {
            if (min == seq.substr(0, s) || min == seq.substr(seq.length() - s, s))
            {
                return true;
            }
        }
        return false;
    }

    inline size_t chash(const char& c) {
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

    inline size_t rol(const size_t& h, const size_t& r) { return (h << r) | (h >> (64-r)); }

    inline size_t ror(const size_t& h, const size_t& r) { return (h >> r) | (h << (64-r)); }

    
    // returns a tuple of (hash, isReverse, isSyncmer)    
    inline std::tuple<size_t, bool, bool> is_syncmer_rollingHash(const std::string &seq, const int s, const bool open, int t = 0) {
      if (seq.size() < s) return std::make_tuple(0, false, false);
      if (!open) t = 0;

      int k = seq.size();
      size_t forwardKmerHash = 0, reverseKmerHash = 0;
      size_t forwardSmerHash = 0, reverseSmerHash = 0;

      size_t minSmerHash   = std::numeric_limits<size_t>::max();
      size_t tthSmerHash  = std::numeric_limits<size_t>::max();
      size_t lastSmerHash  = std::numeric_limits<size_t>::min();

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
        }
        
        if      (forwardSmerHash < reverseSmerHash && forwardSmerHash < minSmerHash) minSmerHash = forwardSmerHash;
        else if (reverseSmerHash < forwardSmerHash && reverseSmerHash < minSmerHash) minSmerHash = reverseSmerHash;
      }

      if (forwardSmerHash < reverseSmerHash) lastSmerHash = forwardSmerHash;
      else                                   lastSmerHash = reverseSmerHash;

      if (minSmerHash == std::numeric_limits<size_t>::max()) return std::make_tuple(0, false, false);
      
      if (forwardKmerHash < reverseKmerHash) {
        if (open) {
          if (minSmerHash == tthSmerHash) return std::make_tuple(forwardKmerHash, false, true);
        } else {
          if (minSmerHash == tthSmerHash || minSmerHash == lastSmerHash) return std::make_tuple(forwardKmerHash, false, true);
        }
      } else if (reverseKmerHash < forwardKmerHash) {
        if (open) {
          if (minSmerHash == tthSmerHash) return std::make_tuple(reverseKmerHash, true, true);
        } else {
          if (minSmerHash == tthSmerHash || minSmerHash == lastSmerHash) return std::make_tuple(reverseKmerHash, true, true);
        }
      }

      return std::make_tuple(0, false, false);
    }


    inline std::pair<size_t, size_t> hashSeq(const std::string& s) {
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

    // returns vector of (hash, isReverse, isSyncmer, startPos)  
    inline std::vector<std::tuple<size_t, bool, bool, int64_t>> rollingSyncmers(const std::string& seq, int k, int s, bool open) {
      std::vector<std::tuple<size_t, bool, bool, int64_t>> syncmers;

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
        forwardSmerHash  = rol(forwardSmerHash, 1) ^ rol(chash(seq[i-s]), s) ^ chash(seq[i]);
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
        syncmers.emplace_back(std::make_tuple(0, false, false, 0));
      } else {
        if (forwardKmerHash < reverseKmerHash) {
          if (open) {
            if (curMinSmerHash == curSmerHashes.front()) {
              syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, 0));
            } else {
              syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
            }
          } else {
            if (curMinSmerHash == curSmerHashes.front() || curMinSmerHash == curSmerHashes.back()) {
              syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, 0));
            } else {
              syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
            }
          }
        } else if (reverseKmerHash < forwardKmerHash) {
          if (open) {
            if (curMinSmerHash == curSmerHashes.front()) {
              syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, 0));
            } else {
              syncmers.emplace_back(std::make_tuple(max_size_t, true, false, 0));
            }
          } else {
            if (curMinSmerHash == curSmerHashes.front() || curMinSmerHash == curSmerHashes.back()) {
              syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, 0));
            } else {
              syncmers.emplace_back(std::make_tuple(max_size_t, true, false, 0));
            }
          }
        } else {
          syncmers.emplace_back(std::make_tuple(max_size_t, false, false, 0));
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
          syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
        } else {
          if (forwardKmerHash < reverseKmerHash) {
            if (open) {
              if (curMinSmerHash == curSmerHashes.front()) {
                syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, i - k + 1));
              } else {
                syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
              }
            } else {
              if (curMinSmerHash == curSmerHashes.front() || curMinSmerHash == curSmerHashes.back()) {
                syncmers.emplace_back(std::make_tuple(forwardKmerHash, false, true, i - k + 1));
              } else {
                syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
              }
            }
          } else if (reverseKmerHash < forwardKmerHash) {
            if (open) {
              if (curMinSmerHash == curSmerHashes.front()) {
                syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, i - k + 1));
              } else {
                syncmers.emplace_back(std::make_tuple(max_size_t, true, false, i - k + 1));
              }
            } else {
              if (curMinSmerHash == curSmerHashes.front() || curMinSmerHash == curSmerHashes.back()) {
                syncmers.emplace_back(std::make_tuple(reverseKmerHash, true, true, i - k + 1));
              } else {
                syncmers.emplace_back(std::make_tuple(max_size_t, true, false, i - k + 1));
              }
            }
          } else {
            syncmers.emplace_back(std::make_tuple(max_size_t, false, false, i - k + 1));
          }
        }
      }

      return syncmers;
    }

    inline std::vector<seedmer> seedmerize(const std::vector<seed> &kmers, const int32_t j)
    {
        std::vector<seedmer> ret;
        if (kmers.size() < j)
        {
            return ret;
        }
        for (size_t i = 0; i < kmers.size() - j + 1; i++)
        {
            seedmer jk = {j, static_cast<int32_t>(kmers[i].seq.length()), ""};
            for (size_t k = 0; k < j; k++)
            {
                jk.positions.push_back(kmers[i + k].pos);
                jk.seq += kmers[i + k].seq;
            }
            ret.push_back(jk);
        }
        return ret;
    }

    inline std::string getNextSyncmer(std::string &seq, const int32_t currPos, const int32_t k, const int32_t s)
    {

        for (int32_t i = currPos; i < seq.size() - k + 1; i++)
        {
            std::string kmer = seq.substr(i, k);
            if (is_syncmer(kmer, s, false))
            {
                return kmer;
            }
        }
        return "";
    }
    inline std::vector<seed> syncmerize(const std::string &seq, const int32_t k, const int32_t s, const bool open, const bool aligned, const int32_t pad)
    {
        std::mutex mtx;
        std::vector<seed> ret;
        int32_t seqLen = seq.size();
        if (seqLen < k)
        {
            return ret;
        }
        if (aligned)
        {
            std::unordered_map<int32_t, int32_t> degap;
            int32_t pos = 0;
            std::string ungapped = "";
            for (int32_t i = 0; i < seqLen; i++)
            {
                char c = seq[i];
                degap[pos] = i;
                if (c != '-')
                {
                    ungapped += c;
                    pos++;
                }
            }

            if (ungapped.size() < k + 1)
            {
                return ret;
            }
            for (int32_t i = 0; i < ungapped.size() - k + 1; i++)
            {
                std::string kmer = ungapped.substr(i, k);
                if (is_syncmer(kmer, s, open))
                {
                    ret.push_back(seed{kmer, degap[i] + pad, -1, false, degap[i + k - 1] + pad});
                }
            }
        }
        else
        {

            // Loop through sequence and build up seeds as we go
            if (seq.size() >= k)
            {

                // Check first k bp for smers
                std::string min_s = seq.substr(0, s);
                int first_min_coord = 0;
                int last_min_coord = 0;
                int Ns = seq[0] == 'N';

                for (int i = 1; i < k - s + 1; i++)
                {
                    std::string smer = seq.substr(i, s);
                    if (seq[i] == 'N')
                        Ns++;

                    if (smer < min_s)
                    {
                        min_s = smer;
                        first_min_coord = i;
                        last_min_coord = i;
                    }
                    else if (smer == min_s)
                    {
                        last_min_coord = i;
                    }
                }
                for (int i = k - s + 1; i < k; i++)
                {
                    if (seq[i] == 'N')
                        Ns++;
                }

                // Check the rest for smers and kmer seeds
                for (int i = k; i < seq.size() + 1; i++)
                {

                    // Processing kmer starting at i - k in seq
                    bool isSeed = (first_min_coord == i - k || last_min_coord == i - s) && Ns <= k / 2;

                    if (isSeed)
                    {
                        std::string kmer = seq.substr(i - k, k);
                        ret.push_back(seed{kmer, i - k, -1, false, i});
                    }

                    // Updating smers for kmer starting at i - k + 1
                    if (i < seq.size())
                    {
                        if (seq[i] == 'N')
                            Ns++;
                        if (seq[i - k] == 'N')
                            Ns--;

                        if (first_min_coord == i - k)
                        {
                            // Were losing the lowest smer, Re search for lowest
                            min_s = seq.substr(i - k + 1, s);
                            first_min_coord = i - k + 1;

                            for (int j = 1; j < k - s + 1; j++)
                            {
                                std::string smer = seq.substr(i - k + 1 + j, s);
                                if (smer < min_s)
                                {
                                    min_s = smer;
                                    first_min_coord = i - k + 1 + j;
                                    last_min_coord = first_min_coord;
                                }
                                else if (smer == min_s)
                                {
                                    last_min_coord = i - k + 1 + j;
                                }
                            }
                        }
                        else
                        {
                            // Test new smer to see if its the lowest
                            std::string smer = seq.substr(i - s + 1, s);
                            if (smer < min_s)
                            {
                                min_s = smer;
                                first_min_coord = i - s + 1;
                                last_min_coord = i - s + 1;
                            }
                            else if (smer == min_s)
                            {
                                last_min_coord = i - s + 1;
                            }
                        }
                    }
                }
            }
        }
        return ret;
    }

}

#endif
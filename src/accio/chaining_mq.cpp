#include <iostream>
#include <vector>
#include <unordered_map>


typedef std::tuple<size_t, int, bool> minimizer_t;
typedef std::tuple<size_t, int, int, bool, int> kminmer_t;
typedef std::tuple<int, int, int, int, bool, int> match_t;

std::vector<minimizer_t> minimizerSketch(const std::string s, int w, int k);
std::vector<kminmer_t> extractKminmers(const std::vector<minimizer_t>& minimizers, int k, int l);
std::vector<match_t> match(const std::vector<kminmer_t>& query_kminmers, const std::unordered_map<size_t, std::tuple<int, int, bool, int>>& ref_unique_kminmers);
std::vector<match_t> chainPseudo(const std::vector<match_t>& matches, const std::pair<int, int>& lens, int maximumGap, int minimumCount, int minimumScore);
std::tuple<std::vector<minimizer_t>, std::vector<minimizer_t>, std::vector<kminmer_t>, std::vector<kminmer_t>, std::unordered_map<size_t, std::tuple<int, int, bool, int>>, std::vector<match_t>, std::vector<match_t>, int>
runAll(const std::string& read, const std::string& ref, int k, int w, int l, int maximumGap, int minimumCount, int minimumScore);
void printMinimizers(const std::vector<minimizer_t>& minimizers);
void printKminmers(const std::vector<kminmer_t>& kminmers);
void printRefIndex(const std::unordered_map<size_t, std::tuple<int, int, bool, int>>& ref_unique_kminmers);
void printMatches(const std::vector<match_t> matches);

int main() {
    std::cout << "What's my purpose?\nYou pass butter\n" << std::endl;

    // Taken from @11915_0_0/1 in panmap/src/test/data/sim_variants_random_nodes/reads/node_250_variant_reads_R1.fastq 
    std::string ref1 = "ATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTG"
                       "CGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAA"
                       "TTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAA"
                       "ACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGC"
                       "CTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATG"
                       "ACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCC";
    
    std::string ref2 = "ATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTG"
                       "CGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTACCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAA"
                       "TTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAA"
                       "ACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGC"
                       "CTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATG"
                       "ACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCC";

    
    std::string r1 = "TGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTT";
              //      |||||||||||||||||||||||||||||||||||||||||| 
              // ref1 TGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTT
              //      ||||||||||||||||||||||:|||||||||||||||||||
              // ref2 TGAGAATTTGACTAAAGAAGGTACCACTACTTGTGGTTACTT
    std::string r2 = "AAGTAACCACAAGTAGTGGCACCTTCTTTAGTCAAATTCTCA"; //revcomp of r1

    std::vector<std::string> refs{ref1, ref2};
    int k = 3;
    int w = 5;
    int l = 3;
    int minimumScore  = 0;
    int minimumCount  = 0;
    int maximumGap    = 50;
    
    {
    auto r1_to_ref1 = runAll(r1, ref1, k, w, l, maximumGap, minimumCount, minimumScore);
    const auto& readMinimizers  = std::get<0>(r1_to_ref1);
    const auto& refMinimizers   = std::get<1>(r1_to_ref1);
    const auto& readKminmers    = std::get<2>(r1_to_ref1);
    const auto& refKminmers     = std::get<3>(r1_to_ref1);
    const auto& refUniqKminmers = std::get<4>(r1_to_ref1);
    const auto& matches         = std::get<5>(r1_to_ref1);
    const auto& pseudoChain     = std::get<6>(r1_to_ref1);
    
    std::cout << "read minimizers" << std::endl;
    printMinimizers(readMinimizers);
    std::cout << "ref minimizers" << std::endl;
    printMinimizers(refMinimizers);
    std::cout << "read kminmers" << std::endl;
    printKminmers(readKminmers);
    std::cout << "ref kminmers" << std::endl;
    printKminmers(refKminmers);
    std::cout << "ref unique kminmers" << std::endl;
    printRefIndex(refUniqKminmers);
    std::cout << "matches" << std::endl;
    printMatches(matches);
    std::cout << "pseudoChain" << std::endl;
    printMatches(pseudoChain);
    std::cout << "score:" << std::get<7>(r1_to_ref1) << std::endl;
    }

    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
    
    {
    auto r1_to_ref1 = runAll(r2, ref1, k, w, l, maximumGap, minimumCount, minimumScore);
    const auto& readMinimizers  = std::get<0>(r1_to_ref1);
    const auto& refMinimizers   = std::get<1>(r1_to_ref1);
    const auto& readKminmers    = std::get<2>(r1_to_ref1);
    const auto& refKminmers     = std::get<3>(r1_to_ref1);
    const auto& refUniqKminmers = std::get<4>(r1_to_ref1);
    const auto& matches         = std::get<5>(r1_to_ref1);
    const auto& pseudoChain     = std::get<6>(r1_to_ref1);
    
    std::cout << "read minimizers" << std::endl;
    printMinimizers(readMinimizers);
    std::cout << "ref minimizers" << std::endl;
    printMinimizers(refMinimizers);
    std::cout << "read kminmers" << std::endl;
    printKminmers(readKminmers);
    std::cout << "ref kminmers" << std::endl;
    printKminmers(refKminmers);
    std::cout << "ref unique kminmers" << std::endl;
    printRefIndex(refUniqKminmers);
    std::cout << "matches" << std::endl;
    printMatches(matches);
    std::cout << "pseudoChain" << std::endl;
    printMatches(pseudoChain);
    std::cout << "score:" << std::get<7>(r1_to_ref1) << std::endl;
    }
}

// Using simple hash function: A, C, G, T = 0, 1, 2, 3
char comp(char base);
std::string revcomp(const std::string& s);
size_t hash(const std::string& s);
size_t btn(char b);
bool isColinear(const match_t& match1, const match_t& match2, int maximumGap);
int scorePseudoChain(const std::vector<match_t>& pseudoChain);

std::tuple<std::vector<minimizer_t>, std::vector<minimizer_t>, std::vector<kminmer_t>, std::vector<kminmer_t>, std::unordered_map<size_t, std::tuple<int, int, bool, int>>, std::vector<match_t>, std::vector<match_t>, int>
runAll(const std::string& read, const std::string& ref, int k, int w, int l, int maximumGap, int minimumCount, int minimumScore) {
    std::vector<minimizer_t> ref_minimizers  = minimizerSketch(ref, w, k);
    std::vector<minimizer_t> read_minimizers = minimizerSketch(read, w, k);
    
    std::vector<kminmer_t> read_kminmers = extractKminmers(read_minimizers, k, l);
    std::vector<kminmer_t> ref_kminmers  = extractKminmers(ref_minimizers, k, l);
    
    std::unordered_map<size_t, std::tuple<int, int, bool, int>> ref_unique_kminmers;
    for (const kminmer_t& kminmer : ref_kminmers) {
        if (ref_unique_kminmers.count(std::get<0>(kminmer)) > 0) {
            ref_unique_kminmers[std::get<0>(kminmer)] = std::make_tuple(-1, -1, false, -1);
        } else {
            ref_unique_kminmers[std::get<0>(kminmer)] = std::make_tuple(std::get<1>(kminmer), std::get<2>(kminmer), std::get<3>(kminmer), std::get<4>(kminmer));
        }
    }

    std::vector<match_t> matches = match(read_kminmers, ref_unique_kminmers);
    std::vector<match_t> pseudoChain = chainPseudo(matches, std::make_pair(read.size(), ref.size()), maximumGap, minimumCount, minimumScore);

    int score = scorePseudoChain(pseudoChain);

    return std::make_tuple(read_minimizers, ref_minimizers, read_kminmers, ref_kminmers, ref_unique_kminmers, matches, pseudoChain, score);
}

int scorePseudoChain(const std::vector<match_t>& pseudoChain) {
    int score = 0;
    for (const match_t& match : pseudoChain) {
        score += std::get<5>(match);
    }
    return score;
}

std::vector<match_t> chainPseudo(const std::vector<match_t>& matches, const std::pair<int, int>& lens, int maximumGap, int minimumCount, int minimumScore) {
    std::vector<match_t> pseudoChain;

    if (matches.size() == 0) {
        return pseudoChain;
    }
    else if (matches.size() == 1) {
        pseudoChain.push_back(matches[0]);
        return pseudoChain;
    }

    size_t maxIndex = 0;
    for (size_t i = 1; i < matches.size(); ++i) {
        if (std::get<5>(matches[i]) > std::get<5>(matches[maxIndex])) maxIndex = i;
    }

    for (size_t i = 0; i < matches.size(); ++i) {
        if (i == maxIndex) {
            pseudoChain.push_back(matches[i]);
            continue;
        }

        if (std::get<0>(matches[maxIndex]) < std::get<0>(matches[i])) {
            if (isColinear(matches[maxIndex], matches[i], maximumGap)) pseudoChain.push_back(matches[i]);
        } else {
            if (isColinear(matches[i], matches[maxIndex], maximumGap)) pseudoChain.push_back(matches[i]);
        }
        
    }

    return pseudoChain;
}

bool isColinear(const match_t& match1, const match_t& match2, int maximumGap) {
    if (std::get<4>(match1) != std::get<4>(match2)) return false;
    
    if (std::get<4>(match1) == false) {
        int qgap = abs(std::get<0>(match2) - std::get<1>(match1));
        int rgap = abs(std::get<2>(match2) - std::get<3>(match1));
        if (std::get<2>(match1) < std::get<2>(match2) && abs(qgap - rgap) < maximumGap) return true;
    } else {
        int qgap = abs(std::get<0>(match2) - std::get<1>(match1));
        int rgap = abs(std::get<2>(match1) - std::get<3>(match2));
        if (std::get<2>(match2) < std::get<2>(match1) && abs(qgap - rgap) < maximumGap) return true;
    }
    return false;
}

// std::tuple<size_t, int, int, bool, int> kminmer_t;
// hash, start, end, reverse, i
std::vector<kminmer_t> extractKminmers(const std::vector<minimizer_t>& minimizers, int k, int l) {
    std::vector<kminmer_t> kminmers;
    
    // first kminmer
    size_t cacheForwardH = 0;
    for (int i = 0; i < l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(minimizers[i]);

    size_t cacheReversedH = 0;
    for (int i = l - 1; i > -1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(minimizers[i]);

    int iorder = 0;
    // Skip if strand ambiguous
    if (cacheForwardH < cacheReversedH) {
        kminmers.emplace_back(cacheForwardH,  std::get<1>(minimizers[0]), std::get<1>(minimizers[l-1]) + k - 1, false, iorder);
        ++iorder;
    } else if (cacheReversedH < cacheForwardH) {
        kminmers.emplace_back(cacheReversedH, std::get<1>(minimizers[0]), std::get<1>(minimizers[l-1]) + k - 1, true, iorder);
        ++iorder;
    }
    
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    for (int i = 1; i < minimizers.size() - l + 1; ++i) {
        cacheForwardH = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(minimizers[i+l-1]);

        cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(minimizers[i+l-1]) << (2 * k * (l - 1)));

        // Skip if strand ambiguous
        if (cacheForwardH < cacheReversedH) {
            kminmers.emplace_back(cacheForwardH,  std::get<1>(minimizers[i]), std::get<1>(minimizers[i+l-1]) + k - 1, false, iorder);
            ++iorder;
        } else if (cacheReversedH < cacheForwardH) {
            kminmers.emplace_back(cacheReversedH, std::get<1>(minimizers[i]), std::get<1>(minimizers[i+l-1]) + k - 1, true, iorder);
            ++iorder;
        }
    }

    return kminmers;
}

int extend(match_t& match, const std::vector<kminmer_t>& query_kminmers, const std::unordered_map<size_t, std::tuple<int, int, bool, int>>& ref_unique_kminmers, int iq, int ir, int c) {
    if (iq == query_kminmers.size() - 1) return c;
    kminmer_t nextqKminmer = query_kminmers[iq+1];
    if (ref_unique_kminmers.count(std::get<0>(nextqKminmer)) > 0 && std::get<3>(ref_unique_kminmers.at(std::get<0>(nextqKminmer))) != -1) {
        const std::tuple<int, int, bool, int>& rKminmer = ref_unique_kminmers.at(std::get<0>(nextqKminmer));
        if (std::get<4>(match) == (std::get<3>(nextqKminmer) != std::get<2>(rKminmer))) {
            if ((std::get<4>(match) == 0 && std::get<3>(rKminmer) == ir + 1) || (std::get<4>(match) == 1 && std::get<3>(rKminmer) == ir - 1)) {
                c+=1;
                if (std::get<4>(match) == 0) {
                    // sq, e'q, sr, e'r, pi, c+1
                    match = std::make_tuple(std::get<0>(match), std::get<2>(nextqKminmer), std::get<2>(match), std::get<1>(rKminmer), std::get<4>(match), c);
                } else if (std::get<4>(match) == 1) {
                    // sq, e'q, s'r, er, pi, c+1
                    match = std::make_tuple(std::get<0>(match), std::get<2>(nextqKminmer), std::get<0>(rKminmer), std::get<3>(match), std::get<4>(match), c);
                }
                return extend(match, query_kminmers, ref_unique_kminmers, std::get<4>(nextqKminmer), std::get<3>(rKminmer), c);
            }
        }
    }
    return c;
}

// typedef std::tuple<int, int, int, int, bool, int> match_t;
// query start, query end, ref start, ref end, strand, count
std::vector<match_t> match(const std::vector<kminmer_t>& query_kminmers, const std::unordered_map<size_t, std::tuple<int, int, bool, int>>& ref_unique_kminmers) {
    std::vector<match_t> matches;
    int i = 0;
    while (i < query_kminmers.size()) {
        const kminmer_t& qKminmer = query_kminmers[i];
        int c = 1;
        if (ref_unique_kminmers.count(std::get<0>(qKminmer)) > 0 && std::get<3>(ref_unique_kminmers.at(std::get<0>(qKminmer))) != -1) {
            const std::tuple<int, int, bool, int>& rKminmer = ref_unique_kminmers.at(std::get<0>(qKminmer));
            match_t curMatch = std::make_tuple(
                std::get<1>(qKminmer), // query start
                std::get<2>(qKminmer), // query end
                std::get<0>(rKminmer), // ref start
                std::get<1>(rKminmer), // ref end
                std::get<3>(qKminmer) != std::get<2>(rKminmer), // 0 if same strand; 1 if opposite
                c
            );
            c = extend(curMatch, query_kminmers, ref_unique_kminmers, std::get<4>(qKminmer), std::get<3>(rKminmer), c);
            matches.push_back(curMatch);
        }
        i+=c;
    }
    return matches;
}

/*
(h, i, r)
    h - hash value (string here)
    i - position
    r - reverse
*/
std::vector<minimizer_t> minimizerSketch(const std::string s, int w, int k) {
    std::vector<minimizer_t> minimizers;
    for (size_t i = 0; i < s.size() - w - k + 2; ++i) {
        size_t m = std::numeric_limits<size_t>::max();

        // Find min value
        for (int j = 0; j < w; ++j) {
            std::string kmer = s.substr(i+j, k);
            size_t u = hash(kmer);
            size_t v = hash(revcomp(kmer));
            if (u != v) {
                m = std::min(m, std::min(u, v));
            }
        }
        // Collect minimizers
        for (int j = 0; j < w; ++j) {
            std::string kmer = s.substr(i+j, k);
            size_t u = hash(kmer);
            size_t v = hash(revcomp(kmer));
            if (u == m && u < v) {
                std::tuple minCand(u , i+j, false);
                if (minimizers.size() == 0) minimizers.push_back(minCand);
                // else if (minimizers.back() != minCand) minimizers.push_back(minCand);
                else if (std::find(minimizers.begin(), minimizers.end(), minCand) == minimizers.end()) minimizers.push_back(minCand);
            } else if (v == m && v < u) {
                std::tuple minCand(v , i+j, true);
                if (minimizers.size() == 0) minimizers.push_back(minCand);
                // if (minimizers.back() != minCand) minimizers.push_back(minCand);
                else if (std::find(minimizers.begin(), minimizers.end(), minCand) == minimizers.end()) minimizers.push_back(minCand);
            }
        }
    }

    return minimizers;
}

size_t pow(size_t b, size_t p) {
    size_t r = 1;
    if (p == 0) return r;
    for (size_t i = 0; i < p; ++i) r *= b; 
    return r;
}

size_t btn(char b) {
    size_t n;
    switch (b) {
    case 'A':
        n = 0;
        break;
    case 'C':
        n = 1;
        break;
    case 'G':
        n = 2;
        break;
    case 'T':
        n = 3;
        break;
    }
    return n;
}

size_t hash(const std::string& s) {
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

std::string revcomp(const std::string& s) {
    std::string cs = "";
    for (int i = s.size() - 1; i > -1; --i) {
        char c = s[i];
        cs += comp(c);
    }
    return cs;
}

char comp(char c) {
    char compC;
    switch (c) { 
    case 'A':
        compC = 'T';
        break;
    case 'C':
        compC = 'G';
        break;
    case 'G':
        compC = 'C';
        break;
    case 'T':
        compC = 'A';
        break;
    default:
        compC = 'N';
        break;
    }
    return compC;
}

void printMinimizers(const std::vector<minimizer_t>& minimizers) {
    std::cout << "minimizers:" << std::endl;
    for (const auto& mini : minimizers) {
        std::cout << std::get<0>(mini) << "\t"
                  << std::get<1>(mini) << "\t"
                  << std::get<2>(mini) << std::endl;
    }
    std::cout << std::endl;
}

void printKminmers(const std::vector<kminmer_t>& kminmers) {
    std::cout << "kminmers:" << std::endl;
    for (const auto& minmer : kminmers) {
        std::cout << std::get<0>(minmer) << "\t"
                  << std::get<1>(minmer) << "\t"
                  << std::get<2>(minmer) << "\t"
                  << std::get<3>(minmer) << "\t"
                  << std::get<4>(minmer) << std::endl;
    }
    std::cout << std::endl;
}

void printRefIndex(const std::unordered_map<size_t, std::tuple<int, int, bool, int>>& ref_unique_kminmers) {
    std::cout << "ref unique kminmers:" << std::endl;
    for (const auto& minmer : ref_unique_kminmers) {
        std::cout << minmer.first << "\t"
                  << std::get<0>(minmer.second) << "\t"
                  << std::get<1>(minmer.second) << "\t"
                  << std::get<2>(minmer.second) << "\t"
                  << std::get<3>(minmer.second) << std::endl;
    }
    std::cout << std::endl;
}

void printMatches(const std::vector<match_t> matches) {
    for (const match_t& match : matches) {
        std::cout << std::get<0>(match) << "\t"
                  << std::get<1>(match) << "\t"
                  << std::get<2>(match) << "\t"
                  << std::get<3>(match) << "\t"
                  << std::get<4>(match) << "\t"
                  << std::get<5>(match) << std::endl;
    }
    std::cout << std::endl;
}
#include <boost/functional/hash.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <iostream>
#include "../mgsr.hpp"
#include "../seeding.hpp"

using namespace std;
namespace fs = boost::filesystem;

typedef std::tuple<size_t, int32_t, int32_t> syncmer_t;
typedef std::tuple<size_t, int32_t, int32_t, bool, int32_t> kminmer_t;

void makeFasta(const std::string& name, const std::string& seq, const std::string& path);
static std::vector<syncmer_t> syncmersSketch(const std::string& seq, const int k, const int s, const bool open);
static std::vector<kminmer_t> extractKminmers(const std::vector<syncmer_t>& syncmers, const int k, const int l);
static std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> extractReadKminmers(const std::vector<std::tuple<size_t, int32_t, int32_t>>& syncmers, const int k, const int l);


// make true kminmers for each seq
BOOST_AUTO_TEST_CASE(kminmers) {
    // ifstream ifs("../dev/examples/sars2k.pmat");
    // boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
    // b.push(boost::iostreams::gzip_decompressor());
    // b.push(ifs);
    // istream is(&b);
    // auto T = new PangenomeMAT::Tree(is);

    // int k = 15;
    // int s = 7;
    // int l = 2;
    // struct seedmer {
    //     size_t hash;
    //     size_t num;
    //     int32_t beg;
    //     int32_t end;
    //     bool rev;
    // };

    // for (const auto& n : T->allNodes) {
    //     std::string node_idn = n.first;
    //     std::string node_seq = T->getStringFromReference(node_idn, false);
    //     std::vector<syncmer_t> syncmers = syncmersSketch(node_seq, k, s, false);
    //     std::vector<kminmer_t> kminmers = extractKminmers(syncmers, k, l);

    //     std::unordered_map<size_t, std::vector<int32_t>> positionMap;
    //     std::map<int32_t, seedmer> kminmersMap;
    //     for (const auto& kminmer : kminmers) {
    //         auto [curHash, curBeg, curEnd, curRev, curI] = kminmer;
    //         size_t curNum = 1;
    //         if (positionMap.find(curHash) != positionMap.end()) {
    //             for (const auto& beg : positionMap[curHash]) {
    //                 kminmersMap[beg].num++;
    //             }
    //             curNum = kminmersMap[positionMap[curHash].front()].num;
    //             positionMap[curHash].push_back(curBeg);
    //         } else {
    //             positionMap[curHash] = {curBeg};
    //         }
    //         kminmersMap[curBeg] = {curHash, curNum, curBeg, curEnd, curRev};
    //     }
        
    //     std::string dirName = "../src/test/data/mgsr/true/k" + std::to_string(k) + "_s" + std::to_string(s) + "_l" + to_string(l);
    //     fs::create_directories(dirName);
    //     std::string outKmmPath = dirName + "/" + node_idn + ".kmi";
    //     std::ofstream okmmfs(outKmmPath);
    //     for (const auto& kminmer : kminmersMap) {
    //         okmmfs << kminmer.second.hash << "\t"
    //                << kminmer.second.num  << "\t"
    //                << kminmer.second.beg  << "\t"
    //                << kminmer.second.end  << "\t"
    //                << kminmer.second.rev  << "\n";
    //     }
    //     okmmfs.close();

    //     fs::create_directories(dirName);
    //     std::string outSyncPath = dirName + "/" + node_idn + ".smi";
    //     std::ofstream osyncfs(outSyncPath);
    //     for (const auto& syncmer : syncmers) {
    //         osyncfs << std::get<0>(syncmer) << "\t"
    //                 << std::get<1>(syncmer) << "\t"
    //                 << std::get<2>(syncmer) << std::endl;
    //     }
    //     osyncfs.close();
    // }
}

BOOST_AUTO_TEST_CASE(_index) {
    // ifstream ifs("../dev/examples/sars2k.pmat");
    // boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
    // b.push(boost::iostreams::gzip_decompressor());
    // b.push(ifs);
    // istream is(&b);
    // auto T = new PangenomeMAT::Tree(is);


    // size_t k = 5;
    // size_t s = 2;
    // size_t l = 6;
    // string indexFilePath = "../dev/examples/sars2k.pmat.spmi";
    // string seedmersIndexPath = "../dev/examples/sars2k_" + std::to_string(k) + "_" + std::to_string(s) + "_" + to_string(l) + ".pmat.kmi";
    // pmi::seedIndex index;
    // std::stringstream seedmersOutStream;

    // mgsr::buildSeedmer(index, T, l, k, s, seedmersOutStream);
    // std::cout << "Writing to " << indexFilePath << "..." << std::endl;
    // std::ofstream fout(indexFilePath);
    // std::ofstream smfout(seedmersIndexPath);
    // fout << index.outStream.str();
    // smfout << seedmersOutStream.str();
    // std::cout << "Done." << std::endl;

    // ifstream indexFile(indexFilePath);

    // mgsr::accio(T, indexFile, k, l);
}
void _initializeFastq(const std::string &fastqPath, std::vector<readSeedmers_t>& readSeedmers, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, const int32_t k, const int32_t s, const int32_t l);
BOOST_AUTO_TEST_CASE(_reads) {
    // int32_t numsample = 1;
    // int32_t numread   = 2;
    // int32_t k = 10;
    // int32_t s = 5;
    // int32_t l = 2;
    // std::vector<std::string> readSequences;
    // std::vector<std::string> readQuals;
    // std::vector<std::string> readNames;
    // std::vector<readSeedmers_t> readSeedmers;
    // std::string reads1Path = "../src/test/data/mgsr/simulated_mgsr_reads/" + std::to_string(numsample) + "_genomes_metagenomics_" + std::to_string(numread) + "k_reads_R1.fastq";
    // _initializeFastq(reads1Path, readSeedmers, readSequences, readQuals, readNames, k, s, l);
}

BOOST_AUTO_TEST_CASE(_score) {
    ifstream ifs("../dev/examples/sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
    b.push(boost::iostreams::gzip_decompressor());
    b.push(ifs);
    istream is(&b);
    auto T = new PangenomeMAT::Tree(is);


    size_t k = 10;
    size_t s = 5;
    size_t l = 2;
    std::string seedmersIndexPath = "../dev/examples/sars2k_" + std::to_string(k) + "_" + std::to_string(s) + "_" + std::to_string(l) + ".pmat.kmi";
    std::ifstream seedmersIndex(seedmersIndexPath);
    // std::vector<int32_t> numSamples = {1, 3, 5, 10};
    std::vector<int32_t> numSamples = {10};
    std::vector<int32_t> numReads = {100};

    for (auto numread : numReads) {
        for (auto numsample : numSamples) {
            std::string r1 = "../src/test/data/mgsr/simulated_mgsr_reads/" + std::to_string(numsample) + "_genomes_metagenomics_" + std::to_string(numread) + "k_reads_R1.fastq";
            // std::string r1 = "../src/test/data/mgsr/tmp/3_genomes_metagenomics_20k_reads_R1.fastq";
            std::string r2 = "../src/test/data/mgsr/simulated_mgsr_reads/" + std::to_string(numsample) + "_genomes_metagenomics_" + std::to_string(numread) + "k_reads_R2.fastq";
            int maximumGap = 10;
            int minimumCount = 0;
            int minimumScore = 0;
            double errorRate = 0.025;
            int32_t ignoreEnds = 50;
            int32_t numReads;
            std::unordered_map<std::string, std::vector<std::pair<int32_t, double>>> allScores;
            std::vector<std::pair<int32_t, std::vector<size_t>>> numReadDuplicates;
            std::unordered_map<std::string, std::string> leastRecentIdenticalAncestor;
            std::unordered_map<std::string, std::unordered_set<std::string>> identicalSets;
            mgsr::scorePseudo(seedmersIndex, r1, r2, allScores, numReadDuplicates, leastRecentIdenticalAncestor, identicalSets, numReads, T, maximumGap, minimumCount, minimumScore, errorRate, ignoreEnds);

            std::string outScoresPath = "../src/test/data/mgsr/pseudoScores/k" + std::to_string(k) + "_s" + std::to_string(s) + "_l" + std::to_string(l) + "/" + std::to_string(numsample) + "_" + std::to_string(numread) + "k.txt";
            ofstream outScores(outScoresPath);
            std::vector<std::pair<std::string, int32_t>> scores;
            for (const auto& node : allScores) {
                int32_t score = 0;
                for (size_t i = 0; i < node.second.size(); ++i) {
                    score += node.second[i].first * numReadDuplicates[i].first;
                }
                scores.emplace_back(std::make_pair(node.first, score));
            }
            std::sort(scores.begin(), scores.end(), [](const auto &a, const auto &b) {
                return a.second > b.second;
            });
            for (const auto& score : scores) {
                outScores << score.first << "\t" << score.second << "\n";
            }
            outScores.close();

            int32_t totalReads = 0;
            for (const auto& n : numReadDuplicates) totalReads += n.first;
            assert(totalReads == numReads);
            
       
            mgsr::squarem(T, allScores, numReadDuplicates, numReads, leastRecentIdenticalAncestor, identicalSets);
            // mgsr::em(T, allScores, numReadDuplicates, numReads, leastRecentIdenticalAncestor, identicalSets);

        }
    }
}

struct kminmerHasher {
    std::size_t operator()(const std::tuple<size_t, int32_t, int32_t, bool, int32_t>& t) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, std::get<0>(t));
        // boost::hash_combine(seed, std::get<1>(t));
        // boost::hash_combine(seed, std::get<2>(t));
        boost::hash_combine(seed, std::get<3>(t));
        boost::hash_combine(seed, std::get<4>(t));
        return seed;
    }
};

struct kminmerSetHasher {
    std::size_t operator()(const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& vec) const {
        std::size_t seed = 0;
        for (const auto& t : vec) {
            boost::hash_combine(seed, kminmerHasher{}(t));
        }
        return seed;
    }
};

struct kminmerSetEqual {
    bool operator()(const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& vec1,
                    const std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>& vec2) const {
        if (vec1.size() != vec2.size()) return false;
        for (size_t i = 0; i < vec1.size(); ++i) {
            if (std::get<0>(vec1[i]) != std::get<0>(vec2[i]) ||
                std::get<3>(vec1[i]) != std::get<3>(vec2[i]) ||
                std::get<4>(vec1[i]) != std::get<4>(vec2[i])) {
                return false;
            }
        }
        return true;
    }
};

void _initializeFastq(const std::string &fastqPath, std::vector<readSeedmers_t>& readSeedmers, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, const int32_t k, const int32_t s, const int32_t l) {
    
    FILE *fp;
    kseq_t *seq;
    fp = fopen(fastqPath.c_str(), "r");
    if(!fp){
        std::cerr << "Error: File " << fastqPath << " not found" << std::endl;
        exit(0);
    }

    seq = kseq_init(fileno(fp));
    std::unordered_map<std::string, int32_t> dupMarkedReads;
    int line;
    while ((line = kseq_read(seq)) >= 0) {
        readSequences.push_back(seq->seq.s);
        if (dupMarkedReads.find(seq->seq.s) == dupMarkedReads.end()) {
            dupMarkedReads[seq->seq.s] = 1;
        } else {
            ++dupMarkedReads[seq->seq.s];
        }
        readNames.push_back(seq->name.s);
        readQuals.push_back(seq->qual.s);
    }

    std::unordered_map<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::pair<std::unordered_set<size_t>, int32_t>, kminmerSetHasher, kminmerSetEqual> dupMarkedKminmers;

    readSeedmers.reserve(readSequences.size());
    for (const auto& dupMarkedRead : dupMarkedReads) {
        const std::string& seq = dupMarkedRead.first;
        const int32_t&  numDup = dupMarkedRead.second;
        std::vector<std::tuple<size_t, int32_t, int32_t>> curSyncmers = syncmersSketch(seq, k, s, false);
        std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> curKminmers = extractReadKminmers(curSyncmers, k, l);
        dupMarkedKminmers[curKminmers.first] = std::make_pair(curKminmers.second, dupMarkedKminmers[curKminmers.first].second + numDup);
        readSeedmers.push_back(std::move(curKminmers));
    }

    std::cerr << "num reads " << readSequences.size() << std::endl;
    std::cerr << "num unique reads " << dupMarkedReads.size() << std::endl;
    std::cerr << "num unique kminmer sets " << dupMarkedKminmers.size() << std::endl;
}

void makeFasta(const std::string& name, const std::string& seq, const std::string& path) {
    if (!fs::exists(path)) {
        std::ofstream faos(path);
        faos << '>' << name << '\n';
        size_t linesize = 80;
        for (size_t i = 0; i < seq.size(); i += linesize) {
            faos << seq.substr(i, std::min(linesize, seq.size() - i)) << '\n';
        }
        faos.close();
    }
}

char comp(char base);
std::string revcomp(const std::string& s);
size_t hashString(const std::string& s);
size_t btn(char b);
std::pair<size_t, bool> getHash(const std::string& s);

static std::vector<syncmer_t> syncmersSketch(const std::string& seq, const int k, const int s, const bool open) {
    std::vector<syncmer_t> syncmers;
    for (size_t i = 0; i < seq.size() - k + 1; ++i) {
        std::string kmer = seq.substr(i, k);
        if (!seeding::is_syncmer(kmer, s, open)) continue;

        std::pair<size_t, bool> kmerHash = getHash(kmer);
        if (!kmerHash.second) continue;

        syncmers.emplace_back(std::make_tuple(kmerHash.first, i, i + k - 1));
    }
    return syncmers;
}

static std::vector<kminmer_t> extractKminmers(const std::vector<syncmer_t>& syncmers, const int k, const int l) {
    std::vector<kminmer_t> kminmers;
    
    // first kminmer
    size_t cacheForwardH = 0;
    for (int i = 0; i < l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(syncmers[i]);

    size_t cacheReversedH = 0;
    for (int i = l - 1; i > -1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(syncmers[i]);

    int iorder = 0;
    // Skip if strand ambiguous
    if (cacheForwardH < cacheReversedH) {
        kminmers.emplace_back(cacheForwardH,  std::get<1>(syncmers[0]), std::get<2>(syncmers[l-1]), false, iorder);
        ++iorder;
    } else if (cacheReversedH < cacheForwardH) {
        kminmers.emplace_back(cacheReversedH, std::get<1>(syncmers[0]), std::get<2>(syncmers[l-1]), true, iorder);
        ++iorder;
    }
    
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    for (int i = 1; i < syncmers.size() - l + 1; ++i) {
        cacheForwardH = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);

        cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(syncmers[i+l-1]) << (2 * k * (l - 1)));

        // Skip if strand ambiguous
        if (cacheForwardH < cacheReversedH) {
            kminmers.emplace_back(cacheForwardH,  std::get<1>(syncmers[i]), std::get<2>(syncmers[i+l-1]), false, iorder);
            ++iorder;
        } else if (cacheReversedH < cacheForwardH) {
            kminmers.emplace_back(cacheReversedH, std::get<1>(syncmers[i]), std::get<2>(syncmers[i+l-1]), true, iorder);
            ++iorder;
        }
    }

    return kminmers;
}

static std::pair<std::vector<std::tuple<size_t, int32_t, int32_t, bool, int32_t>>, std::unordered_set<size_t>> extractReadKminmers(const std::vector<std::tuple<size_t, int32_t, int32_t>>& syncmers, const int k, const int l) {
    std::vector<std::tuple<size_t, int32_t, int32_t, bool, int>> kminmers;
    std::unordered_set<size_t> hashes;

    // first kminmer
    size_t cacheForwardH = 0;
    for (int i = 0; i < l; ++i) cacheForwardH = (cacheForwardH << (2 * k)) + std::get<0>(syncmers[i]);

    size_t cacheReversedH = 0;
    for (int i = l - 1; i > -1; --i) cacheReversedH = (cacheReversedH << (2 * k)) + std::get<0>(syncmers[i]);

    int iorder = 0;
    // Skip if strand ambiguous
    if (cacheForwardH < cacheReversedH) {
        kminmers.emplace_back(cacheForwardH,  std::get<1>(syncmers[0]), std::get<2>(syncmers[l-1]), false, iorder);
        hashes.insert(cacheForwardH);
        ++iorder;
    } else if (cacheReversedH < cacheForwardH) {
        kminmers.emplace_back(cacheReversedH, std::get<1>(syncmers[0]), std::get<2>(syncmers[l-1]), true, iorder);
        hashes.insert(cacheReversedH);
        ++iorder;
    }
    
    size_t mask = 0;
    for (int i = 0; i < 2 * k * (l - 1); i++) mask = (mask << 1) + 1;

    for (int i = 1; i < syncmers.size() - l + 1; ++i) {
        cacheForwardH = ((cacheForwardH & mask) << (k * 2)) + std::get<0>(syncmers[i+l-1]);

        cacheReversedH = (cacheReversedH >> (2 * k)) + (std::get<0>(syncmers[i+l-1]) << (2 * k * (l - 1)));

        // Skip if strand ambiguous
        if (cacheForwardH < cacheReversedH) {
            kminmers.emplace_back(cacheForwardH,  std::get<1>(syncmers[i]), std::get<2>(syncmers[i+l-1]), false, iorder);
            hashes.insert(cacheForwardH);
            ++iorder;
        } else if (cacheReversedH < cacheForwardH) {
            kminmers.emplace_back(cacheReversedH, std::get<1>(syncmers[i]), std::get<2>(syncmers[i+l-1]), true, iorder);
            hashes.insert(cacheReversedH);
            ++iorder;
        }
    }

    return std::make_pair(std::move(kminmers), std::move(hashes));
}

std::pair<size_t, bool> getHash(const std::string& s) {
    try {
        size_t u = hashString(s);
        size_t v = hashString(revcomp(s));
        if (u < v) return std::make_pair(u, true);
        else if (v < u) return std::make_pair(v, true);
        return std::make_pair(0, false);
    } catch (std::invalid_argument) {
        return std::make_pair(0, false); // skip strand ambiguous
    }
    return std::make_pair(0, false);
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
    default:
        throw std::invalid_argument("Kmer contains non canonical base");
        break;
    }
    return n;
}

size_t hashString(const std::string& s) {
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


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
std::vector<syncmer_t> syncmersSketch(const std::string& seq, const int k, const int s, const bool open);
std::vector<kminmer_t> extractKminmers(const std::vector<syncmer_t>& syncmers, const int k, const int l);

// make true kminmers for each seq
// BOOST_AUTO_TEST_CASE(kminmers) {
//     ifstream ifs("../dev/examples/sars2k.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
//     b.push(boost::iostreams::gzip_decompressor());
//     b.push(ifs);
//     istream is(&b);
//     auto T = new PangenomeMAT::Tree(is);

//     int k = 10;
//     int s = 5;
//     int l = 2;
//     struct seedmer {
//         size_t hash;
//         size_t num;
//         int32_t beg;
//         int32_t end;
//         bool rev;
//     };

//     for (const auto& n : T->allNodes) {
//         std::string node_idn = n.first;
//         std::string node_seq = T->getStringFromReference(node_idn, false);
//         std::vector<syncmer_t> syncmers = syncmersSketch(node_seq, k, s, false);
//         std::vector<kminmer_t> kminmers = extractKminmers(syncmers, k, l);

//         std::unordered_map<size_t, std::vector<int32_t>> positionMap;
//         std::map<int32_t, seedmer> kminmersMap;
//         for (const auto& kminmer : kminmers) {
//             auto [curHash, curBeg, curEnd, curRev, curI] = kminmer;
//             size_t curNum = 1;
//             if (positionMap.find(curHash) != positionMap.end()) {
//                 for (const auto& beg : positionMap[curHash]) {
//                     kminmersMap[beg].num++;
//                 }
//                 curNum = kminmersMap[positionMap[curHash].front()].num;
//                 positionMap[curHash].push_back(curBeg);
//             } else {
//                 positionMap[curHash] = {curBeg};
//             }
//             kminmersMap[curBeg] = {curHash, curNum, curBeg, curEnd, curRev};
//         }
        
//         std::string outKmmPath = "../src/test/data/kminmers/" + node_idn + ".kmi";
//         std::ofstream okmmfs(outKmmPath);
//         for (const auto& kminmer : kminmersMap) {
//             okmmfs << kminmer.second.hash << "\t"
//                    << kminmer.second.num  << "\t"
//                    << kminmer.second.beg  << "\t"
//                    << kminmer.second.end  << "\t"
//                    << kminmer.second.rev  << "\n";
//         }
//         okmmfs.close();

//         std::string outSyncPath = "../src/test/data/kminmers/" + node_idn + ".smi";
//         std::ofstream osyncfs(outSyncPath);
//         for (const auto& syncmer : syncmers) {
//             osyncfs << std::get<0>(syncmer) << "\t"
//                     << std::get<1>(syncmer) << "\t"
//                     << std::get<2>(syncmer) << std::endl;
//         }
//         osyncfs.close();
//     }
// }

BOOST_AUTO_TEST_CASE(tmp) {
    ifstream ifs("../dev/examples/sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
    b.push(boost::iostreams::gzip_decompressor());
    b.push(ifs);
    istream is(&b);
    auto T = new PangenomeMAT::Tree(is);

    string indexFilePath = "../dev/examples/sars2k.pmat.spmi";
    string seedmersIndexPath = "../dev/examples/sars2k.pmat.kmi";
    size_t k = 10;
    size_t s = 5;
    size_t l = 2;

    pmi::seedIndex index;
    std::stringstream seedmersOutStream;

    mgsr::buildSeedmer(index, T, l, k, s, seedmersOutStream);
    std::cout << "Writing to " << indexFilePath << "..." << std::endl;
    std::ofstream fout(indexFilePath);
    std::ofstream smfout(seedmersIndexPath);
    fout << index.outStream.str();
    smfout << seedmersOutStream.str();
    std::cout << "Done." << std::endl;

    ifstream indexFile(indexFilePath);

    mgsr::accio(T, indexFile, k, l);
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

std::vector<syncmer_t> syncmersSketch(const std::string& seq, const int k, const int s, const bool open) {
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

std::vector<kminmer_t> extractKminmers(const std::vector<syncmer_t>& syncmers, const int k, const int l) {
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
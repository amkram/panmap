#include "../pmi.hpp"
#include "../seeding.hpp"
#include "../tree.hpp"
#include "PangenomeMAT.hpp"
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_CASE(_removeIndices) {
  /* should erase elements at the given positions */
  using namespace seeding;
  std::vector<seed> a = {seed{"zero"},  seed{"one"},  seed{"two"},
                         seed{"three"}, seed{"four"}, seed{"five"}};

  // smaller values at the top
  std::stack<int32_t> indices;

  indices.push(4);
  indices.push(2);
  indices.push(0);

  std::vector<seed> expected = {seed{"one"}, seed{"three"}, seed{"five"}};

  removeIndices(a, indices);

  BOOST_TEST(a == expected);
}

std::vector<std::tuple<std::string, int, int>>
extractSeedmers(const std::string &seq, const int k, const int s, const int l,
                const bool open);
void writeSeedmersFromIndexFile(
    seedmerIndex_t &index, std::unordered_map<int32_t, std::string> &seedmers,
    Tree *T, const Node *node, const std::string &testDirName);

BOOST_AUTO_TEST_CASE(seedmersIndex) {
  // construct true seeds
  std::ifstream ifs("../dev/examples/sars2k.pmat");
  boost::iostreams::filtering_streambuf<boost::iostreams::input> b;
  b.push(boost::iostreams::gzip_decompressor());
  b.push(ifs);
  std::istream is(&b);
  auto T = new PangenomeMAT::Tree(is);

  // std::vector<std::tuple<int, int, int>> parameters = {{10, 5, 2}, {10, 5,
  // 1}, {15, 7, 3}, {15, 7, 1}};
  std::vector<std::tuple<int, int, int>> parameters = {{10, 5, 1}};

  for (const auto &param : parameters) {
    auto k = std::get<0>(param);
    auto s = std::get<1>(param);
    auto l = std::get<2>(param);

    std::string dirName =
        std::string("../src/test/data/pmi_test/true_seedmers") + "_k" +
        std::to_string(k) + "_s" + std::to_string(s) + "_l" +
        std::to_string(l) + "/";
    fs::create_directories(dirName);
    std::cerr << "Printing seedmers from sequence, k: " << k << ", s: " << s
              << ", l: " << l << " in " << dirName << std::endl;
    for (const auto &n : T->allNodes) {
      std::string node_idn = n.first;
      std::string outSeedmersPath = dirName + node_idn + ".pmi";
      if (fs::exists(outSeedmersPath))
        continue;

      std::string node_seq = T->getStringFromReference(node_idn, true);
      std::vector<std::tuple<std::string, int, int>> seedmers =
          extractSeedmers(node_seq, k, s, l, false);

      std::ofstream osdmfs(outSeedmersPath);
      osdmfs << node_seq << std::endl;

      for (const auto &seedmer : seedmers) {
        osdmfs << std::get<0>(seedmer) << "\t" << std::get<1>(seedmer) << "\t"
               << std::get<2>(seedmer) << std::endl;
      }
      osdmfs.close();
    }

    std::string testDirName =
        std::string("../src/test/data/pmi_test/test_seedmers") + "_k" +
        std::to_string(k) + "_s" + std::to_string(s) + "_l" +
        std::to_string(l) + "/";
    std::cerr << "Building seedmers from index, k: " << k << ", s: " << s
              << ", l: " << l << " in " << dirName << std::endl;

    pmi::seedIndex index;
    pmi::build(index, T, l, k, s);
    std::string indexFile = std::string("../src/test/data/pmi_test/") + "k" +
                            std::to_string(k) + "_s" + std::to_string(s) +
                            "_l" + std::to_string(l) + ".pmi";
    std::ofstream fout(indexFile);
    fout << index.outStream.str();
    fout.close();
    fs::create_directories(testDirName);

    std::ifstream indexIf(indexFile);
    std::string line;
    std::getline(indexIf, line);
    std::vector<std::string> spltTop;
    PangenomeMAT::stringSplit(line, ' ', spltTop);
    int32_t tempK = std::stoi(spltTop[0]);
    int32_t tempS = std::stoi(spltTop[1]);
    int32_t tempL = std::stoi(spltTop[2]);

    seedmerIndex_t seedmerIndex;
    while (std::getline(indexIf, line)) {
      std::vector<std::string> splt;
      stringSplit(line, ' ', splt);
      std::string nid = splt[0];
      seedmerIndex[nid] = {};
      for (int32_t i = 1; i < splt.size(); i++) {
        std::vector<std::string> splt2;
        stringSplit(splt[i], ':', splt2);
        int32_t pos = std::stoi(splt2[0]);
        std::string seedmer = splt2[1];
        seedmerIndex[nid].insert(std::make_pair(pos, seedmer));
      }
    }

    std::cerr << "Printing seedmers from index, k: " << k << ", s: " << s
              << ", l: " << l << " in " << dirName << std::endl;

    std::unordered_map<int32_t, std::string> dynamicSeedmersTarget;
    writeSeedmersFromIndexFile(seedmerIndex, dynamicSeedmersTarget, T, T->root,
                               testDirName);
  }
}

static void mutateSeedmerMap(
    std::unordered_map<int32_t, std::string> &seedmers, const std::string &nid,
    seedmerIndex_t &index,
    std::vector<std::pair<int32_t, std::string>> &seedmersToRevert) {
  // std::cout << "Forward => " << nid << "\n";
  for (const auto &op : index[nid]) {
    std::string seedmer = op.second;
    int32_t pos = op.first;
    if (seedmer[0] == '@') { // deletion
      if (seedmers.find(pos) != seedmers.end()) {
        seedmers.erase(seedmers.find(pos));
      }
    } else { // insertion
      seedmersToRevert.push_back({pos, seedmers[pos]});
      seedmers[pos] = seedmer;
    }
  }
}

static void revertSeedmerMap(
    std::unordered_map<int32_t, std::string> &seedmers, const std::string &nid,
    seedmerIndex_t &index,
    std::vector<std::pair<int32_t, std::string>> &seedmersToRevert) {
  for (auto it = index[nid].rbegin(); it != index[nid].rend(); it++) {
    std::string seedmer = it->second;
    int32_t pos = it->first;
    if (seedmer[0] == '@') { // deletion
      seedmers[pos] = seedmer.substr(1);
    } else {
      if (seedmers.find(pos) != seedmers.end()) {
        seedmers.erase(seedmers.find(pos));
      }
    }
  }
  for (const auto &p : seedmersToRevert) {
    if (p.second.size() > 1) {
      seedmers[p.first] = p.second;
    }
  }
}

void writeSeedmersFromIndexFile(
    seedmerIndex_t &index, std::unordered_map<int32_t, std::string> &seedmers,
    Tree *T, const Node *node, const std::string &testDirName) {
  std::vector<std::pair<int32_t, std::string>> seedmersToRevert;
  mutateSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);

  std::vector<std::pair<int32_t, std::string>> seedmersToWrite;
  seedmersToWrite.reserve(seedmers.size());
  for (const auto &s : seedmers)
    seedmersToWrite.emplace_back(std::make_pair(s.first, s.second));

  std::sort(seedmersToWrite.begin(), seedmersToWrite.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  std::ofstream osdmfs(testDirName + node->identifier + ".pmi");
  for (const auto &seedmer : seedmersToWrite) {
    osdmfs << seedmer.second << "\t" << seedmer.first << std::endl;
  }
  osdmfs.close();

  for (Node *child : node->children) {
    writeSeedmersFromIndexFile(index, seedmers, T, child, testDirName);
  }

  revertSeedmerMap(seedmers, node->identifier, index, seedmersToRevert);
}
std::vector<std::tuple<std::string, int, int>> extractSeedmers(const std::string &seq, const int k, const int s, const int l, const bool open) {
  std::vector<std::tuple<std::string, int, int>> syncmers;
  std::unordered_map<int32_t, int32_t> degap;
  int64_t pos = 0;
  std::string ungapped = "";
  for (int64_t i = 0; i < seq.size(); i++) {
    char c = seq[i];
    degap[pos] = i;
    if (c != '-' && c != 'x') {
      ungapped += c;
      pos++;
    }
  }
  if (ungapped.size() < k + 1) {
    return syncmers;
  }
  for (size_t i = 0; i < ungapped.size() - k + 1; ++i) {
    std::string kmer = ungapped.substr(i, k);
    if (seeding::is_syncmer(kmer, s, open)) {
      syncmers.emplace_back(std::make_tuple(kmer, degap[i], degap[i + k - 1]));
    }
  }

  std::vector<std::tuple<std::string, int, int>> seedmers;
  for (size_t i = 0; i < syncmers.size() - l + 1; ++i) {
    std::string seedmerSeq = "";
    for (size_t j = 0; j < l; ++j)
      seedmerSeq += std::get<0>(syncmers[i + j]);
    seedmers.emplace_back(std::make_tuple(seedmerSeq, std::get<1>(syncmers[i]),
                                          std::get<2>(syncmers[i + l - 1])));
  }

  return seedmers;
}
// BOOST_AUTO_TEST_CASE(_getRecomputePositions)
// {
//     using namespace tree;
//     int32_t k = 3;

//     range_t a = {6, 1, 0};
//     range_t b = {22, 3, 0};
//     range_t c = {31, 2, 0};

//     //                 length:  0               3        2
//     //                          *               ***      **
//     //      recompute range: xxxxxxx       xxxxxxxxxx  xxxxx
//     std::string gapped = "ABCD-EF-GH-IJKLM-N-O--P-QRSTUVWXYZ";
//     //                    0123456789111111111122222222223333
//     //                              012345678901234567890123

//     range_t expected_a = {3, 9, 0};
//     range_t expected_b = {17, 26, 0};
//     range_t expected_c = {29, 33, 0};
//     getRecomputePositions(expected_a, gapped, k);

//     BOOST_TEST(
//         getRecomputePositions(a, gapped, k) == expected_a
//     );
//     BOOST_TEST(
//         getRecomputePositions(b, gapped, k) == expected_b
//     );
//     BOOST_TEST(
//         getRecomputePositions(c, gapped, k) == expected_c
//     );
// }

// // BOOST_AUTO_TEST_CASE(_indexSeeds) {
// //     size_t k = 13;
// //     size_t s = 6;
// //     std::ifstream is("../sars2k.pmat");
// //     boost::iostreams::filtering_streambuf< boost::iostreams::input>
// inPMATBuffer;
// //     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
// //     inPMATBuffer.push(is);
// //     std::istream inputStream(&inPMATBuffer);
// //     Tree *T = new Tree(inputStream);
// //     std::ofstream os("./test.out");
// //     indexSeeds(T, os, k, s);
// // }

// BOOST_AUTO_TEST_CASE(accuracy) {
//     std::vector<std::string> files = {
//         "../src/test/test.fastq"
//         // "../fastq_2k/BS001339.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/BS001694.1.1kreads.fastq_R1.fastq"
//         // "../fastq_2k/BS002022.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/BS002229.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/BS003222.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/BS003231.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON085554.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON092190.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON094090.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON100161.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON104875.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON105009.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON105170.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON105456.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON106383.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON111925.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON112672.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON123083.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON131307.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON132382.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON134820.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/ON137843.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OQ193850.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OQ194908.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OQ196547.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OQ196749.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OX589733.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OX590546.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OX591210.1.1kreads.fastq_R1.fastq",
//         // "../fastq_2k/OX591273.1.1kreads.fastq_R1.fastq"
//         };

//     size_t k = 13;
//     size_t s = 7;
//     std::ifstream is("../src/test/test.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input>
//     inPMATBuffer; inPMATBuffer.push(boost::iostreams::gzip_decompressor());
//     inPMATBuffer.push(is);
//     std::istream inputStream(&inPMATBuffer);
//     Tree *T = new Tree(inputStream);
//     std::ofstream os("./test.out");
//     indexSeeds(T, os, k, s);

//     is.close();
//     os.close();
//     PangenomeMAT::Node *root = T->root;
//     struct seedIndex index;
//     std::ifstream indexFile("./test.out");
//     PangenomeMAT::loadIndex(T->root, indexFile, index);

//     for (std::string f : files) {
//         PangenomeMAT::placeSample(T, f, index, k, s);
//     }
// }

// // // BOOST_AUTO_TEST_CASE(_indexSeeds) {
// // //     size_t k = 25;
// // //     size_t s = 2;
// // //     std::ifstream is("../mammal_mito_refseq.pmat");
// // // //    std::string fastqPath = "";
// // //     boost::iostreams::filtering_streambuf< boost::iostreams::input>
// inPMATBuffer;
// // //     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
// // //     inPMATBuffer.push(is);
// // //     std::istream inputStream(&inPMATBuffer);
// // //     Tree *T = new Tree(inputStream);
// // //     std::ofstream os("./test.out");
// // //     std::cout << "indexing...\n";
// // //     indexSeeds(T, os, k, s);

// // //     std::cout << "done.\n";
// // //     is.close();
// // //     os.close();

// // //     PangenomeMAT::Node *root = T->root;
// // //     std::vector<read_t> reads;

// // //     struct seedIndex index;
// // //     std::ifstream indexFile("./test.out");

// // //     PangenomeMAT::loadIndex(T->root, indexFile, index);

// // //     auto fastq_start = std::chrono::high_resolution_clock::now();
// // //     std::set<seed> readSyncmers = syncmersFromFastq("../r1.fastq",
// reads, k, s);
// // //     auto fastq_end = std::chrono::high_resolution_clock::now();

// // //     std::cout << "fastq time: " <<
// std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end -
// fastq_start).count() << "\n";

// // //     auto place_start = std::chrono::high_resolution_clock::now();

// // //     std::set<seed> rootSyncmers =
// std::set<seed>(index.rootSeeds.begin(), index.rootSeeds.end());

// // //     std::cerr << "\n";
// // //     std::cerr << "Placing sample...\n";

// // //     struct dynamicJaccard dj;

// // //     dj.intersectionSize = intersection_size(rootSyncmers,
// readSyncmers);
// // //     dj.unionSize = rootSyncmers.size() + readSyncmers.size() -
// dj.intersectionSize;
// // //     dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;

// // //     std::cout << "root seeds: " << rootSyncmers.size() << "\n";
// // //     std::cout << "read seeds: " << readSyncmers.size() << "\n";
// // //     for (const auto &k : readSyncmers) {
// // //         std::cout << k.seq << "\n";
// // //     }
// // //     std::cout << "initial jaccard: " << dj.jaccardIndex << "\n";

// // //     std::unordered_map<std::string, float> scores;
// // //     std::unordered_map<std::string, bool> readSyncmersMap;
// // //     for (const auto &k : readSyncmers) {
// // //         readSyncmersMap[k.seq] = true;
// // //     }

// // //     placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj,
// scores);

// // //     auto place_end = std::chrono::high_resolution_clock::now();

// // //     std::cout << "place time: " <<
// std::chrono::duration_cast<std::chrono::milliseconds>(place_end -
// place_start).count() << "\n";

// // //     std::vector<std::pair<std::string, float>> v;
// // //     for ( const auto &p : scores ) {
// // //         v.push_back(std::make_pair(p.first, p.second));
// // //     }
// // //     std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
// // //         return left.second > right.second;
// // //     });

// // //     std::string best_match = v[0].first;
// // //     for (const auto &s : v) {
// // //         std::cerr << s.first << ": " << s.second << "\n";
// // //     }

// // // }

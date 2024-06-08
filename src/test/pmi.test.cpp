#include "../pmi.hpp"
#include "../seeding.hpp"
#include "../tree.hpp"
#include "PangenomeMAT.hpp"
#include "index.pb.h"
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <string>
//#include <sys/_types/_int64_t.h>

namespace fs = boost::filesystem;

std::vector<std::tuple<std::string, int, int>>
extractSeedmers(const std::string &seq, const int k, const int s, const int l,
                const bool open);
void writeSeedmersFromIndexFile(
    seedmerIndex_t &index, std::unordered_map<int32_t, std::string> &seedmers,
    Tree *T, const Node *node, const std::string &testDirName);

void findSyncmers(
    const std::string &sequence, int k, int s,
    std::unordered_map<std::string, std::vector<int32_t>> &refSeeds) {
  std::unordered_set<std::string> syncmers;
  int len = sequence.length();

  for (int i = 0; i <= len - k; ++i) {
    std::string kmer = sequence.substr(i, k);
    std::string start_smer = kmer.substr(0, s);
    std::string end_smer = kmer.substr(k - s, s);

    bool start_is_min = true, end_is_min = true;

    for (int j = 1; j < k - s; ++j) {
      std::string smer = kmer.substr(j, s);
      if (start_smer > smer) {
        start_is_min = false;
      }
      if (end_smer > smer) {
        end_is_min = false;
      }
    }

    if (start_is_min || end_is_min) {

      if (refSeeds.find(kmer) == refSeeds.end()) {
        refSeeds[kmer] = {};
      }
      refSeeds[kmer].push_back(i);
    }
  }
}

// Recursive function to build the seed index
void buildHelper2(SeedmerIndex &index, Tree *T, Node *node,
                 int32_t &pb_i, std::map<int32_t, std::string> &seedmersAlex) {
  // std::cout << "BH2 " << node->identifier << "\n";





  std::cout << "pbi " << pb_i << "\n";
  // std::cout << "size " << index.per_node_mutations_size() << "\n";



  NodeSeedmerMutations pb_node_mutations = index.per_node_mutations(pb_i);  //somthing bad

  // std::cout << "other " << pb_node_mutations.mutations_size() << "\n";

  std::string node_idn = node->identifier;
  bool is22 = (node->identifier == "node_22");


  std::vector<std::pair<int32_t, std::string>> backtrack;
  std::vector<int64_t> delSeeds;
  std::vector<std::pair<int64_t, std::string>> addSeeds;
  
  for (int mut_i = 0; mut_i < pb_node_mutations.mutations_size(); mut_i++) {
    const auto &mut = pb_node_mutations.mutations(mut_i);
    int32_t pos = mut.pos();
    if (mut.is_deletion()) {
      backtrack.push_back({pos, seedmersAlex[pos]});
      delSeeds.push_back(pos);
    } else {
      if (seedmersAlex.find(pos) != seedmersAlex.end()) {
        backtrack.push_back({pos, seedmersAlex[pos]});
      } else {
        backtrack.push_back({pos, ""});
      }
      addSeeds.push_back({pos, mut.seq()});
      
    }
  }
  for (const auto &p : delSeeds) {
    seedmersAlex.erase(p);
  }
  for (const auto &p : addSeeds) {
    seedmersAlex[p.first] = p.second;
  }

  std::string node_seq = tree::getStringAtNode(node, T, true);
  std::string node_seq_nogap = tree::getStringAtNode(node, T, false);

  std::unordered_map<int32_t, int32_t> degap;
  int64_t pos = 0;
  std::string ungapped = "";
  for (int64_t i = 0; i < node_seq.size(); i++) {
    char c = node_seq[i];
    degap[i] = pos;
    if (c != '-' && c != 'x') {
      ungapped += c;
      pos++;
    }
  }
  std::string dirName = std::string("../dev/eval-performance/");

  std::ofstream osdmfsAlex(dirName + node_idn + ".true.alex.pmi");
  osdmfsAlex << node_seq << std::endl;
  osdmfsAlex << node_seq_nogap << std::endl;
  for (const auto &pair : seedmersAlex) {
    osdmfsAlex << pair.second << "\t" << degap[pair.first] << "\t" << pair.first
               << std::endl;
  }
  osdmfsAlex.close();

  /* Recursive step */
  for (Node *child : node->children) {
    if(is22){
      //exit(0);
    }
    pb_i++;
    buildHelper2(index, T, child, pb_i, seedmersAlex);
  }
  // undo seed mutations
  for (const auto &p : backtrack) {
    if (p.second.size() > 0) {
      seedmersAlex[p.first] = p.second;
    } else {
      seedmersAlex.erase(p.first);
    }
  }
}

BOOST_AUTO_TEST_CASE(performance) {
  std::ifstream ifs("../dev/examples/sars2k.pmat");
  boost::iostreams::filtering_streambuf<boost::iostreams::input> b;
  b.push(boost::iostreams::gzip_decompressor());
  b.push(ifs);
  std::istream is(&b);
  PangenomeMAT::Tree *T = new PangenomeMAT::Tree(is);

  std::vector<std::tuple<int, int, int>> parameters = {{15, 8, 1}};

  for (const auto &param : parameters) {
    auto k = std::get<0>(param);
    auto s = std::get<1>(param);
    auto j = std::get<2>(param);

    std::string dirName = std::string("../dev/eval-performance/");
    fs::create_directories(dirName);

    SeedmerIndex index;
    pmi::build(index, T, j, k, s);
    std::ofstream fout("../dev/eval-performance/sars2k.pmi");
    index.SerializeToOstream(&fout);
    std::map<int32_t, std::string> seedmersAlex;
    int32_t pb_i = 0;

    buildHelper2(index, T, T->root, pb_i, seedmersAlex);

    exit(0);

    for (auto &n : T->allNodes) {
      std::string node_idn = n.first;
      std::string outSeedmersPath = dirName + node_idn + ".true.alan.pmi";
      std::string outSeedmersPathNico = dirName + node_idn + ".true.nico.pmi";
      Node *nod= n.second;
      
      std::string node_seq = tree::getStringAtNode(nod, T, true);
       std::string node_seq_nogap = tree::getStringAtNode(nod, T, false);

      std::vector<std::tuple<std::string, int, int>> seedmers =
          extractSeedmers(node_seq, k, s, j, false);
      std::vector<std::tuple<std::string, int, int>> seedmers_nogap =
          extractSeedmers(node_seq_nogap, k, s, j, false);
      std::ofstream osdmfs(outSeedmersPath);
      osdmfs << node_seq << std::endl;
      osdmfs << node_seq_nogap << std::endl;

      for (int i = 0; i < seedmers.size(); i++) {
        osdmfs << std::get<0>(seedmers[i]) << "\t"
               << std::get<1>(seedmers_nogap[i]) << "\t" << std::get<1>(seedmers[i])
                << std::endl;
               
      }
      osdmfs.close();

      std::unordered_map<std::string, std::vector<int32_t>> NseedToRefPositions;
      NseedToRefPositions = {};
      findSyncmers(node_seq_nogap, k, s, NseedToRefPositions);

      // Collecting reference Seeds
      std::vector<seed> refSeeds;
      for (auto kv : NseedToRefPositions) {
        for (int32_t pos : kv.second) {
          seed thisSeed;
          thisSeed.seq = kv.first;
          thisSeed.pos = pos;
          refSeeds.push_back(thisSeed);
        }
      }

      // Sorting ref seeds
      sort(refSeeds.begin(), refSeeds.end(),
           [](const seed &lhs, const seed &rhs) { return lhs.pos < rhs.pos; });

      std::ofstream osdmfsNico(outSeedmersPathNico);
      osdmfsNico << node_seq << std::endl;
      osdmfsNico << node_seq_nogap << std::endl;

      for (const auto &s : refSeeds) {
        osdmfsNico << s.seq << "\t" << s.pos << std::endl;
      }
      osdmfsNico.close();
    }
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
std::vector<std::tuple<std::string, int, int>>
extractSeedmers(const std::string &seq, const int k, const int s, const int l,
                const bool open) {
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


/*
BOOST_AUTO_TEST_CASE(getNucSeq)
{
    using namespace tree;

    size_t k = 5;
    size_t s = 3;
    tree::mutableTreeData data;
    std::string seq = "----A-A-CCG---TLEMONC--CC---------";
    const sequence_t sequence = {
      {
        { // block 0
          {'-', {}}
        },{}
      },
        //block 1
      { {
          {'A', {'-', '-', '-', 'A', '-'}},
  //(0,0,-1)^     ^ (0,0,0)
          {'-', {'-', 'C', 'C', 'G'}},
        }, {}//unused
      },
      {
        // block 2 - OFF
       {
        {'-', {}},
        {'G', {}},
       }, {}//unused
      },
      { // block 3
       {
        {'T', {}},
        {'C', {'L', 'E', 'M', 'O', 'N'}},
        {'C', {'-', '-'}},
        {'C', {}},
        {'x', {'-'}}
       }, {}//unused
      },
      {
        // block 4
       {
        {'-', {'P','L','S','-','N','O'}},
       }, {}//unused
      },
      
    };
    
    globalCoords_t globalCoords = {
      {//block 0
        {
          {0, {}},
        }, {}//unused
      },
      {//block 1
        {
          {6, {1, 2, 3, 4, 5}},
          {11, {7, 8, 9, 10}},
        }, {}//unused
      },
      {//block 2
        
        {
          {12, {}},
          {13, {}},
        }, {}//unused
      },
      { // block 3
       {
        {14, {}},
        {20, {15, 16, 17, 18, 19}},
        {23, {21, 22}},
        {24, {}},
        {26, {25}}
       }, {}//unused
      },
      {//block 4
        
        {
          {33, {27, 28, 29, 30, 31, 32}},
        }, {}//unused
      },
    };
    // block starts: 0:^       1:^       2:^
    //
    //                 ---A-A-CCG- | TLEMONC--CC--
    //                 012345...     ^11...    
    tupleCoord_t start = {0,0,0};
    tupleCoord_t end = {-1,-1,-1};
    const blockExists_t blockExists = {{{false},{}}, {{true},{}}, {{false},{}}, {{true},{}}, {{false},{}}};
    const blockStrand_t blockStrand = {{{true},{}}, {{true},{}}, {{true},{}}, {{true},{}}, {{true},{}}};
    const Tree *T;
    const Node *node;

    CoordNavigator navigator(sequence);

    std::string wholeSeq =
      tree::getNucleotideSequenceFromBlockCoordinates(start, end, sequence,
      blockExists, blockStrand, T, node, globalCoords, navigator);
    std::cout << wholeSeq << "\n";
    BOOST_TEST(wholeSeq == seq);
    //BOOST_TEST(false);
    std::cout << "sup homie, you really wanna knowme?\n";


}
*/


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

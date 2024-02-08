#include <boost/test/unit_test.hpp>
#include <stack>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <filesystem>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include "PangenomeMAT.hpp"
#include "../seeding.hpp"

using namespace PangenomeMAT;
using namespace seeding;


BOOST_AUTO_TEST_CASE(_isSyncmer) {
    
/* should return true only for syncmers */
    std::string a = "ABCDEF";

    BOOST_TEST(
        is_syncmer(a, 3, true) == true
    );
    BOOST_TEST(
		is_syncmer(a, 3, false) == true
    );
    
    std::string b = "ZZZAAAZZZ";

    BOOST_TEST(
		is_syncmer(b, 6, true) == false
    );
    BOOST_TEST(
		is_syncmer(b, 6, false) == true
    );

    std::string c = "Z";
 
    BOOST_TEST(
		is_syncmer(c, 1, true) == true
    );
}

BOOST_AUTO_TEST_CASE(_syncmerize) {
/* should pick syncmers correctly */
    std::string a = "AAAAAABB"; // k=3 s=2 open
    /*             * AAA
                   *  AAA
                   *   AAA
                   *    AAA
                   *     AAB
                   *      ABB
    */
    std::vector<seed> ret_a = {
        seed{"AAA", 0, -1},
        seed{"AAA", 1, -1},
        seed{"AAA", 2, -1},
        seed{"AAA", 3, -1},
        seed{"AAB", 4, -1},
        seed{"ABB", 5, -1}                
    };
    std::vector<seed> ret_a_padded = {
        seed{"AAA", 42, -1},
        seed{"AAA", 43, -1},
        seed{"AAA", 44, -1},
        seed{"AAA", 45, -1},
        seed{"AAB", 46, -1},
        seed{"ABB", 47, -1}                
    };
    BOOST_TEST(
		syncmerize(a, 3, 2, true, false, 0) == ret_a
    );
    BOOST_TEST(
		syncmerize(a, 3, 2, true, false, 42) == ret_a_padded
    );
    std::string b = "AABCBAAACC"; // k=6 s=3 open
  /*               * AABCBA 
                   *  ABCBAA
                       BCBAAA
                        CBAAAC
                         BAAACC
  */
    std::vector<seed> ret_b = {
        seed{"AABCBA", 0, -1},
        seed{"ABCBAA", 1, -1}
    };
    std::string c = "AABCBAAACC"; // k=6 s=3 closed
  /*               * AABCBA 
                   *  ABCBAA
                   *   BCBAAA
                        CBAAAC
                         BAAACC
  */
    std::vector<seed> ret_c = {
        seed{"AABCBA", 0, -1},
        seed{"ABCBAA", 1, -1},
        seed{"BCBAAA", 2, -1}
    };
    BOOST_TEST(
		syncmerize(b, 6, 3, true, false, 0) == ret_b
    );
    BOOST_TEST(
		syncmerize(c, 6, 3, false, false, 0) == ret_c
    );
    std::string d = "--A-A-B--CBAAAC-C-"; // k=6 s=3 closed aligned
    /*             * AABCBA 
                   *  ABCBAA
                   *   BCBAAA
                        CBAAAC
                         BAAACC
    */
     std::vector<seed> ret_d = {
        seed{"AABCBA", 2, -1},
        seed{"ABCBAA", 4, -1},
        seed{"BCBAAA", 6, -1}
    };
    BOOST_TEST(
		syncmerize(d, 6, 3, false, true, 0) == ret_d
    );
}

BOOST_AUTO_TEST_CASE(_seedmerize) {

/* should pick syncmers correctly */
    std::string a = "AAAAAABCGCGCGCGCGCGGCDSJHFDJHDJFHJJJEEEOOOIOIIOAISOIFAPSNFAJKSBNFKJWRBKJVBSAKLJDECNKDBCJBDEFBBAAAAAABB";
    
    
		std::vector<seed> syncmers = syncmerize(a, 6, 3, true, false, 0);
    for (const auto &s : syncmers) {
        std::cout << s.seq << "@" << s.pos << "\n";
    }
    std::cout << "jk:\n";
    std::vector<seedmer> seedmers = seedmerize(syncmers, 3);
    for (const auto &j : seedmers) {
        std::cout << j.seq << "@";
        for (const auto &p : j.positions) {
            std::cout << p << " ";
        }
        std::cout << "\n";
    }
}

// template<>
// struct boost::test_tools::tt_detail::print_log_value<std::pair<int,int>> {
//     void operator()(std::ostream& os, std::pair<int,int> const& value) {
//         os << '{' << value.first << ',' << value.second << '}';
//     }
// };

// BOOST_AUTO_TEST_CASE(_getRecomputePositions)
// {
//     int32_t k = 3;

//     std::pair<int32_t, int32_t> a = {6, 1};
//     std::pair<int32_t, int32_t> b = {22, 3};
//     std::pair<int32_t, int32_t> c = {31, 2};

//     //                 length:  0               3        2
//     //                          *               ***      **  
//     //      recompute range: xxxxxxx       xxxxxxxxxx  xxxxx
//     std::string gapped = "ABCD-EF-GH-IJKLM-N-O--P-QRSTUVWXYZ";
//     //                    0123456789111111111122222222223333
//     //                              012345678901234567890123

//     std::pair<int32_t, int32_t> expected_a = {3, 9};
//     std::pair<int32_t, int32_t> expected_b = {17, 26};
//     std::pair<int32_t, int32_t> expected_c = {29, 33};

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

// bool eq(tbb::concurrent_unordered_map<std::string, seed> a, tbb::concurrent_unordered_map<std::string, seed> b) {
//     if (a.size() != b.size()) {
//         return false;
//     }
//     for (auto &p : a) {
//         if (b.find(p.first) == b.end()) {
//             return false;
//         }
//         if (b[p.first].seq != p.second.seq || b[p.first].idx != p.second.idx) {
//             return false;
//         }
//     }
//     return true;
// }


// BOOST_AUTO_TEST_CASE(_discardSyncmers) {
// /*  should discard syncmers contained in any range in B 
//     unless they are going to be inserted. Then remove from
//     to_insert (based on sequence for now). Also, track variable syncmers.
// */
//     size_t k = 4;

//     std::vector<seed> seeds = {
//         seed{"AAAA", 0, -1}, // x
//         seed{"ACCA", 2, -1}, //
//         seed{"AABA", 7, -1}, // x
//         seed{"AACA", 8, -1}, // x
//         seed{"ABDA", 9, -1}, //
//         seed{"ABBA", 15,-1}, // gapped length exceeds bound  
//         seed{"CDEF", 23,-1}, 
//     };
//     std::string seq = "xxxxxxxxxxxxxxxABB---Axxxxxx"; // len 28

//     tbb::concurrent_vector<std::pair<int32_t, int32_t>> B;
//     B.push_back({0, 4});
//     B.push_back({7, 11});
//     B.push_back({14, 20});
//     B.push_back({22, 35});

//     std::vector<seed> expected_seeds = {
//         seed{"ACCA", 2, -1},
//         seed{"ABDA", 9, -1},
//         seed{"ABBA", 15,-1},
//         seed{"CDEF", 23,-1}, 
//     };

//     tbb::concurrent_unordered_map<std::string, seed> to_insert = {
//         {"XYZW", seed{"XYZW", 8, -1}},
//         {"CDEF", seed{"CDEF", 23, -1}}
//     };
//     tbb::concurrent_unordered_map<std::string, seed> expected_to_insert = {
//         {"XYZW", seed{"XYZW", 8, -1}}
//     };
 

//     seedIndex index;
    
//     std::vector<seed> expected_deletions = {
//                 seed{"AACA", 8, 3},
//                 seed{"AABA", 7, 2},
//                 seed{"AAAA", 0, 0},
//     };
//     std::string nid = "my_node_id";

  

//     std::vector<int32_t> positions;
//     for (seed s : seeds) {
//         positions.push_back(s.pos);
//     }
//     tbb::concurrent_unordered_map<int32_t, int32_t> gappedEnds = getGappedEnds(seq, k, positions);
//     discardSyncmers(seeds, B, seq, gappedEnds, to_insert, index, nid, k);

//     BOOST_TEST(
//         seeds == expected_seeds
//     );
//     BOOST_TEST(
//         eq(to_insert, expected_to_insert)
//     );

//     BOOST_TEST(
//         index.deletions[nid] == expected_deletions
//     );
//     BOOST_TEST(index.deletions[nid][0].pos == 8);
//     BOOST_TEST(index.deletions[nid][0].idx == 3);
//     BOOST_TEST(index.deletions[nid][1].pos == 7);
//     BOOST_TEST(index.deletions[nid][1].idx == 2);
//     BOOST_TEST(index.deletions[nid][2].pos == 0);
//     BOOST_TEST(index.deletions[nid][2].idx == 0);

// }

// BOOST_AUTO_TEST_CASE(_indexSeeds) {
//     size_t k = 13;
//     size_t s = 6;
//     std::ifstream is("../sars2k.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
//     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
//     inPMATBuffer.push(is);
//     std::istream inputStream(&inPMATBuffer);
//     Tree *T = new Tree(inputStream);
//     std::ofstream os("./test.out");
//     indexSeeds(T, os, k, s);
// }

// BOOST_AUTO_TEST_CASE(accuracy) {
//     std::vector<std::string> files = {
//         "../src/test/test.fastq"
        // "../fastq_2k/BS001339.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/BS001694.1.1kreads.fastq_R1.fastq"
        // "../fastq_2k/BS002022.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/BS002229.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/BS003222.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/BS003231.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON085554.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON092190.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON094090.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON100161.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON104875.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON105009.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON105170.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON105456.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON106383.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON111925.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON112672.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON123083.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON131307.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON132382.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON134820.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/ON137843.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OQ193850.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OQ194908.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OQ196547.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OQ196749.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OX589733.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OX590546.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OX591210.1.1kreads.fastq_R1.fastq",
        // "../fastq_2k/OX591273.1.1kreads.fastq_R1.fastq"
//         };s

// // BOOST_AUTO_TEST_CASE(_indexSeeds) {
// //     size_t k = 25;
// //     size_t s = 2;   
// //     std::ifstream is("../mammal_mito_refseq.pmat");
// // //    std::string fastqPath = "";
// //     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
// //     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
// //     inPMATBuffer.push(is);
// //     std::istream inputStream(&inPMATBuffer);
// //     Tree *T = new Tree(inputStream);
// //     std::ofstream os("./test.out");
// //     std::cout << "indexing...\n";
// //     indexSeeds(T, os, k, s);
    
// //     std::cout << "done.\n";
// //     is.close();
// //     os.close();

// //     PangenomeMAT::Node *root = T->root;
// //     std::vector<read_t> reads;

// //     struct seedIndex index;
// //     std::ifstream indexFile("./test.out");

// //     PangenomeMAT::loadIndex(T->root, indexFile, index);

// //     auto fastq_start = std::chrono::high_resolution_clock::now();
// //     std::set<seed> readSyncmers = syncmersFromFastq("../r1.fastq", reads, k, s);
// //     auto fastq_end = std::chrono::high_resolution_clock::now();

// //     std::cout << "fastq time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end - fastq_start).count() << "\n";

// //     auto place_start = std::chrono::high_resolution_clock::now();

// //     std::set<seed> rootSyncmers = std::set<seed>(index.rootSeeds.begin(), index.rootSeeds.end());

// //     std::cerr << "\n";
// //     std::cerr << "Placing sample...\n";


// //     struct dynamicJaccard dj;
 
// //     dj.intersectionSize = intersection_size(rootSyncmers, readSyncmers);
// //     dj.unionSize = rootSyncmers.size() + readSyncmers.size() - dj.intersectionSize;
// //     dj.jaccardIndex = (float)dj.intersectionSize / dj.unionSize;
    
    
// //     std::cout << "root seeds: " << rootSyncmers.size() << "\n";
// //     std::cout << "read seeds: " << readSyncmers.size() << "\n";
// //     for (const auto &k : readSyncmers) {
// //         std::cout << k.seq << "\n";
// //     }
// //     std::cout << "initial jaccard: " << dj.jaccardIndex << "\n";

// //     std::unordered_map<std::string, float> scores;
// //     std::unordered_map<std::string, bool> readSyncmersMap;
// //     for (const auto &k : readSyncmers) {
// //         readSyncmersMap[k.seq] = true;
// //     }

// //     placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj, scores);

// //     auto place_end = std::chrono::high_resolution_clock::now();

// //     std::cout << "place time: " << std::chrono::duration_cast<std::chrono::milliseconds>(place_end - place_start).count() << "\n";

// //     std::vector<std::pair<std::string, float>> v;
// //     for ( const auto &p : scores ) {
// //         v.push_back(std::make_pair(p.first, p.second));
// //     } 
// //     std::sort(v.begin(), v.end(), [] (auto &left, auto &right) {
// //         return left.second > right.second;
// //     });

// //     std::string best_match = v[0].first;
// //     for (const auto &s : v) {
// //         std::cerr << s.first << ": " << s.second << "\n";
// //     }

// // }


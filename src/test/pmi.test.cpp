#include <boost/test/unit_test.hpp>
#include "../pmi.hpp"


BOOST_AUTO_TEST_CASE(_removeIndices) {
/* should erase elements at the given positions */
    using namespace seeding;
    std::vector<seed> a = {
        seed{"zero"},
        seed{"one"},
        seed{"two"},
        seed{"three"},
        seed{"four"},
        seed{"five"}      
    };
    
    // smaller values at the top
    std::stack<int32_t> indices;

    indices.push(4);
    indices.push(2);
    indices.push(0);

    std::vector<seed> expected = {
        seed{"one"},
        seed{"three"},
        seed{"five"}
    };

    removeIndices(a, indices);

    BOOST_TEST(
        a == expected
    );

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
// //     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
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
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
//     inPMATBuffer.push(boost::iostreams::gzip_decompressor());
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
// // //     boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
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
// // //     std::set<seed> readSyncmers = syncmersFromFastq("../r1.fastq", reads, k, s);
// // //     auto fastq_end = std::chrono::high_resolution_clock::now();

// // //     std::cout << "fastq time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end - fastq_start).count() << "\n";

// // //     auto place_start = std::chrono::high_resolution_clock::now();

// // //     std::set<seed> rootSyncmers = std::set<seed>(index.rootSeeds.begin(), index.rootSeeds.end());

// // //     std::cerr << "\n";
// // //     std::cerr << "Placing sample...\n";


// // //     struct dynamicJaccard dj;
 
// // //     dj.intersectionSize = intersection_size(rootSyncmers, readSyncmers);
// // //     dj.unionSize = rootSyncmers.size() + readSyncmers.size() - dj.intersectionSize;
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

// // //     placeDFS(root, index.rootSeeds, readSyncmersMap, index, dj, scores);

// // //     auto place_end = std::chrono::high_resolution_clock::now();

// // //     std::cout << "place time: " << std::chrono::duration_cast<std::chrono::milliseconds>(place_end - place_start).count() << "\n";

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


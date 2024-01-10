

// BOOST_AUTO_TEST_CASE(_removeIndices) {
// /* should erase elements at the given positions */
//     std::vector<kmer_t> a = {
//         kmer_t{"zero"},
//         kmer_t{"one"},
//         kmer_t{"two"},
//         kmer_t{"three"},
//         kmer_t{"four"},
//         kmer_t{"five"}      
//     };
    
//     // smaller values at the top
//     std::stack<int32_t> indices;

//     indices.push(4);
//     indices.push(2);
//     indices.push(0);

//     std::vector<kmer_t> expected = {
//         kmer_t{"one"},
//         kmer_t{"three"},
//         kmer_t{"five"}
//     };

//     removeIndices(a, indices);

//     BOOST_TEST(
//         a == expected
//     );

// }

// // BOOST_AUTO_TEST_CASE(_getRecomputePositions)
// // {
// //     int32_t k = 3;

// //     std::pair<int32_t, int32_t> a = {6, 1};
// //     std::pair<int32_t, int32_t> b = {22, 3};
// //     std::pair<int32_t, int32_t> c = {31, 2};

// //     //                 length:  0               3        2
// //     //                          *               ***      **  
// //     //      recompute range: xxxxxxx       xxxxxxxxxx  xxxxx
// //     std::string gapped = "ABCD-EF-GH-IJKLM-N-O--P-QRSTUVWXYZ";
// //     //                    0123456789111111111122222222223333
// //     //                              012345678901234567890123

// //     std::pair<int32_t, int32_t> expected_a = {3, 9};
// //     std::pair<int32_t, int32_t> expected_b = {17, 26};
// //     std::pair<int32_t, int32_t> expected_c = {29, 33};

// //     BOOST_TEST(
// //         getRecomputePositions(a, gapped, k) == expected_a
// //     );
// //     BOOST_TEST(
// //         getRecomputePositions(b, gapped, k) == expected_b
// //     );
// //     BOOST_TEST(
// //         getRecomputePositions(c, gapped, k) == expected_c
// //     );
// // }

// // bool eq(tbb::concurrent_unordered_map<std::string, kmer_t> a, tbb::concurrent_unordered_map<std::string, kmer_t> b) {
// //     if (a.size() != b.size()) {
// //         return false;
// //     }
// //     for (auto &p : a) {
// //         if (b.find(p.first) == b.end()) {
// //             return false;
// //         }
// //         if (b[p.first].seq != p.second.seq || b[p.first].idx != p.second.idx) {
// //             return false;
// //         }
// //     }
// //     return true;
// // }


// // BOOST_AUTO_TEST_CASE(_discardSyncmers) {
// // /*  should discard syncmers contained in any range in B 
// //     unless they are going to be inserted. Then remove from
// //     to_insert (based on sequence for now). Also, track variable syncmers.
// // */
// //     size_t k = 4;

// //     std::vector<kmer_t> seeds = {
// //         kmer_t{"AAAA", 0, -1}, // x
// //         kmer_t{"ACCA", 2, -1}, //
// //         kmer_t{"AABA", 7, -1}, // x
// //         kmer_t{"AACA", 8, -1}, // x
// //         kmer_t{"ABDA", 9, -1}, //
// //         kmer_t{"ABBA", 15,-1}, // gapped length exceeds bound  
// //         kmer_t{"CDEF", 23,-1}, 
// //     };
// //     std::string seq = "xxxxxxxxxxxxxxxABB---Axxxxxx"; // len 28

// //     tbb::concurrent_vector<std::pair<int32_t, int32_t>> B;
// //     B.push_back({0, 4});
// //     B.push_back({7, 11});
// //     B.push_back({14, 20});
// //     B.push_back({22, 35});

// //     std::vector<kmer_t> expected_seeds = {
// //         kmer_t{"ACCA", 2, -1},
// //         kmer_t{"ABDA", 9, -1},
// //         kmer_t{"ABBA", 15,-1},
// //         kmer_t{"CDEF", 23,-1}, 
// //     };

// //     tbb::concurrent_unordered_map<std::string, kmer_t> to_insert = {
// //         {"XYZW", kmer_t{"XYZW", 8, -1}},
// //         {"CDEF", kmer_t{"CDEF", 23, -1}}
// //     };
// //     tbb::concurrent_unordered_map<std::string, kmer_t> expected_to_insert = {
// //         {"XYZW", kmer_t{"XYZW", 8, -1}}
// //     };
 

// //     seedIndex index;
    
// //     std::vector<kmer_t> expected_deletions = {
// //                 kmer_t{"AACA", 8, 3},
// //                 kmer_t{"AABA", 7, 2},
// //                 kmer_t{"AAAA", 0, 0},
// //     };
// //     std::string nid = "my_node_id";

  

// //     std::vector<int32_t> positions;
// //     for (kmer_t s : seeds) {
// //         positions.push_back(s.pos);
// //     }
// //     tbb::concurrent_unordered_map<int32_t, int32_t> gappedEnds = getGappedEnds(seq, k, positions);
// //     discardSyncmers(seeds, B, seq, gappedEnds, to_insert, index, nid, k);

// //     BOOST_TEST(
// //         seeds == expected_seeds
// //     );
// //     BOOST_TEST(
// //         eq(to_insert, expected_to_insert)
// //     );

// //     BOOST_TEST(
// //         index.deletions[nid] == expected_deletions
// //     );
// //     BOOST_TEST(index.deletions[nid][0].pos == 8);
// //     BOOST_TEST(index.deletions[nid][0].idx == 3);
// //     BOOST_TEST(index.deletions[nid][1].pos == 7);
// //     BOOST_TEST(index.deletions[nid][1].idx == 2);
// //     BOOST_TEST(index.deletions[nid][2].pos == 0);
// //     BOOST_TEST(index.deletions[nid][2].idx == 0);

// // }

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
// // //     std::set<kmer_t> readSyncmers = syncmersFromFastq("../r1.fastq", reads, k, s);
// // //     auto fastq_end = std::chrono::high_resolution_clock::now();

// // //     std::cout << "fastq time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fastq_end - fastq_start).count() << "\n";

// // //     auto place_start = std::chrono::high_resolution_clock::now();

// // //     std::set<kmer_t> rootSyncmers = std::set<kmer_t>(index.rootSeeds.begin(), index.rootSeeds.end());

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


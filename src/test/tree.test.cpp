#include <boost/test/unit_test.hpp>
#include <stack>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <filesystem>
#include <tbb/concurrent_unordered_map.h>
#include "panmanUtils.hpp"
#include "../pmi.hpp"

using namespace panmanUtils;
using namespace seeding;
using namespace seed_annotated_tree;
using namespace pmi;


BOOST_AUTO_TEST_CASE(sequenceReconstruction) {

    std::ifstream file("../dev/examples/sars2k.fa");
    std::unordered_map<std::string, std::string> expectedSequences;

    std::string id;
    while (std::getline(file, id)) {
        std::string seq;
        std::getline(file, seq);
        expectedSequences[id.substr(1)] = seq;
    }

    /* Panmat loading */
    std::ifstream is("../dev/examples/sars2k.pmat");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::gzip_decompressor());
    inPMATBuffer.push(is);
    std::istream inputStream(&inPMATBuffer);
    panmanUtils::Tree *T = new panmanUtils::Tree(inputStream);

    /* Get all sequences */
    std::unordered_map<std::string, std::string> resultSequences = seed_annotated_tree::getAllNodeStrings(T);
  
    /* Do they match expected sequences from aligned fasta? */
    for (auto r : expectedSequences) {
        BOOST_TEST(resultSequences[r.first] == expectedSequences[r.first]);
    }
}
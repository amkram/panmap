#include <boost/test/unit_test.hpp>
#include <stack>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <json/json.h>
#include <filesystem>
#include <tbb/concurrent_unordered_map.h>
#include "PangenomeMAT.hpp"
#include "../pmi.hpp"

using namespace PangenomeMAT;
using namespace seed;
using namespace tree;
using namespace pmi;


void getSequencesDFS(std::unordered_map<std::string, std::string> &result, mutableTreeData &data, Tree *T, const Node *node, const globalCoords_t &globalCoords) {

    blockMutData_t blockMutData;
    nucMutData_t nucMutData;

    /*  Mutate with block and nuc mutations */
    applyMutations(data, blockMutData, nucMutData, T, node, globalCoords);

    /* Use current state of mutableTreeData to decode node's sequence */
    std::string seq = tree::getAlignedSequence(data, T, node, true);

    result[node->identifier] = seq;

    /* Recursive step */
    for(Node* child: node->children){
        getSequencesDFS(result, data, T, child, globalCoords);
    }

    /* Undo mutations when backtracking */
    seedIndex blank;
    undoMutations(data, blank, T, node, blockMutData, nucMutData);
}

BOOST_AUTO_TEST_CASE(sequenceReconstruction) {

    std::ifstream file("test.aligned.fa");
    std::unordered_map<std::string, std::string> expectedSequences;

    std::string id;
    while (std::getline(file, id)) {
        std::string seq;
        std::getline(file, seq);
        expectedSequences[id.substr(1)] = seq;
    }

    /* Panmat loading */
    std::ifstream is("test.pmat");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::gzip_decompressor());
    inPMATBuffer.push(is);
    std::istream inputStream(&inPMATBuffer);
    PangenomeMAT::Tree *T = new PangenomeMAT::Tree(inputStream);

    /* Setup */
    tree::mutableTreeData data;
    tree::globalCoords_t globalCoords;
    setup(data, globalCoords, T);

    std::unordered_map<std::string, std::string> resultSequences;
    getSequencesDFS(resultSequences, data, T, T->root, globalCoords);

    /* Verify that each sequence in a test tree matches both panmat lib's
    ** getStringFromReference() and our DFS method's output
    */
    for (auto r : expectedSequences) {
        BOOST_TEST(resultSequences[r.first] == expectedSequences[r.first]);
        BOOST_TEST(resultSequences[r.first] == T->getStringFromReference(r.first));
    }
}
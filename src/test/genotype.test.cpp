#include <boost/test/unit_test.hpp>
#include <stack>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <json/json.h>
#include <filesystem>
#include <tbb/concurrent_unordered_map.h>
#include "PangenomeMAT.hpp"
#include "../pmi.hpp"
#include "../genotype.hpp"

using namespace std;
using namespace tree;
using namespace PangenomeMAT;

BOOST_AUTO_TEST_CASE(genotyping) {

    /* Panmat loading */
    std::ifstream is("../dev/examples/sars_20000.pmat");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inPMATBuffer;
    inPMATBuffer.push(boost::iostreams::gzip_decompressor());
    inPMATBuffer.push(is);
    std::istream inputStream(&inPMATBuffer);
    PangenomeMAT::Tree *T = new PangenomeMAT::Tree(inputStream);

    // make tmp pileup data
    ofstream pileup("../src/test/data/tmp.pileup");
    pileup << "ref\t1\tT\t4\t.C.^3,\tE>FF\n"
            << "ref\t2\tC\t3\tTt,\tFFF\n"
            << "ref\t3\tG\t4\t..,,\tFFFF\n"
            << "ref\t4\tA\t7\tGGgggGG\t:::::FF\n"
            << "ref\t5\tA\t4\t.+2CC.+2CC,+2CC.+2CC\tCCFE\n"
            << "ref\t6\tT\t6\t.,,.,,\tFFFFFF\n"
            << "ref\t7\tG\t4\t.-3TTA.-3TTA,-3TTA.\t:::,\n";
    pileup.close();
    ifstream inPileup("../src/test/data/tmp.pileup");


    // genotype likelihoods alone
    ofstream ofEmptyMutMat("../src/test/data/empty.tmp.mutmat");
    ofEmptyMutMat << "0.0 0.0 0.0 0.0\n"
                    << "0.0 0.0 0.0 0.0\n"
                    << "0.0 0.0 0.0 0.0\n"
                    << "0.0 0.0 0.0 0.0\n"
                    << "0.0 0.0 0.0 0.0 0.0 0.0\n"
                    << "0.0 0.0 0.0 0.0 0.0 0.0\n";
    ofEmptyMutMat.close();

    ifstream inEmptyMutMat("../src/test/data/empty.tmp.mutmat");
    mutationMatrices emptyMutMat = mutationMatrices();
    fillMutationMatrices(emptyMutMat, T, &inEmptyMutMat);
    inEmptyMutMat.close();

    streambuf* originalCoutBuffer = std::cout.rdbuf();
    stringstream buffer;
    cout.rdbuf(buffer.rdbuf());
    genotype::printSamplePlacementVCF(inPileup, emptyMutMat, false, 0);
    cout.rdbuf(originalCoutBuffer);
    string capturedVCF = buffer.str();

    unordered_map<int, vector<string> > vcfOut;
    vector<string> lines;
    stringSplit(capturedVCF, '\n', lines);
    for (const auto& line : lines) {
        if (line.substr(0, 3) != "ref") {continue; }
        vector<string> fields;
        stringSplit(line, '\t', fields);
        vcfOut[stoi(fields[1])] = {fields[4], fields.back()};
    }

    unordered_map<int, vector<string> > expOut;
    expOut[1] = {"C", "0:3,1:0,81"};
    expOut[2] = {"T", "1:1,2:36,0"};
    expOut[4] = {"G", "1:0,7:198,0"};
    expOut[5] = {"ACC", "1:0,4:140,0"};
    expOut[7] = {"G", "1:1,3:64,0"};

    inPileup.clear();
    inPileup.seekg(0, ios::beg);

    for (const auto& var : expOut) {
        int pos = var.first;
        BOOST_TEST(expOut[pos][0] == vcfOut[pos][0]);
        BOOST_TEST(expOut[pos][1] == vcfOut[pos][1]);
    }

    genotype::printSamplePlacementVCF(inPileup, emptyMutMat, true, 0);

    // // genotype with genotype prior
    // ofstream fullmMutMat("../src/test/data/full.tmp.mutmat");
    // fullMutMat << "1   20  25  23\n"
    //        << "21  1   24  12\n"
    //        << "20  22  1   23\n"
    //        << "20  21  19  1\n"

    inPileup.close();
    filesystem::remove("../src/test/data/tmp.pileup");
    filesystem::remove("../src/test/data/empty.tmp.mutmat");

    // ofstream vcfof("../talk/root.vcf");
    // printVCF(T, T->root->identifier, vcfof);


}
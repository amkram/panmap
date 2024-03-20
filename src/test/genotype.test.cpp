#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <stack>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include "PangenomeMAT.hpp"
#include "../genotype.hpp"

using namespace std;
using namespace tree;
namespace fs = boost::filesystem;

BOOST_AUTO_TEST_CASE(genotyping) {

    ifstream puif("../src/test/data/test.pileup");

    /*############################
    # genotype likelihoods alone #
    ############################*/
    ifstream emif("../src/test/data/test_empty.mm");
    mutationMatrices emptyMutMat = mutationMatrices();
    fillMutationMatricesFromFile(emptyMutMat, emif);
    emif.close();

    ofstream vcfof("../src/test/data/vcf.out.tmp");
    genotype::printSamplePlacementVCF(puif, emptyMutMat, false, 0, vcfof);
    vcfof.close();
    puif.clear();
    puif.seekg(0, ios::beg);

    unordered_map<int, vector<string> > vcfOut;
    ifstream vcfif("../src/test/data/vcf.out.tmp");
    string line;
    while (getline(vcfif, line)) {
        if (line.substr(0, 3) != "ref") {continue;}
        vector<string> fields;
        stringSplit(line, '\t', fields);
        vcfOut[stoi(fields[1])] = {fields[4], fields.back()};
    }
    vcfif.close();
    fs::remove("../src/test/data/vcf.out.tmp");

    unordered_map<int, vector<string> > expOut;
    expOut[1] = {"C", "0:3,1:0,81"};
    expOut[2] = {"T", "1:1,2:36,0"};
    expOut[4] = {"G", "1:0,7:198,0"};
    expOut[5] = {"ACC", "1:0,4:140,0"};
    expOut[7] = {"G", "1:1,3:64,0"};

    for (const auto& var : expOut) {
        int pos = var.first;
        BOOST_TEST(expOut[pos][0] == vcfOut[pos][0]);
        BOOST_TEST(expOut[pos][1] == vcfOut[pos][1]);
    }


    /*##############################
    # genotype likelihoods + Prior #
    ##############################*/
    emif.open("../src/test/data/test.mm");
    mutationMatrices mutmat = mutationMatrices();
    fillMutationMatricesFromFile(mutmat, emif);
    emif.close();

    vcfof.open("../src/test/data/vcf.out.tmp");
    genotype::printSamplePlacementVCF(puif, mutmat, false, 0, vcfof);
    vcfof.close();
    puif.close();
    
    vcfOut.clear();
    vcfif.open("../src/test/data/vcf.out.tmp");
    line.clear();
    while (getline(vcfif, line)) {
        if (line.substr(0, 3) != "ref") {continue;}
        vector<string> fields;
        stringSplit(line, '\t', fields);
        vcfOut[stoi(fields[1])] = {fields[4], fields.back()};
    }
    vcfif.close();
    fs::remove("../src/test/data/vcf.out.tmp");

    expOut.clear();
    expOut[1] = {"C", "0:3,1:0,101"};
    expOut[2] = {"T", "1:1,2:25,0"};
    expOut[4] = {"G", "1:0,7:174,0"};
    expOut[5] = {"ACC", "1:0,4:91,0"};
    expOut[7] = {"G", "1:1,3:5,0"};


    for (const auto& var : expOut) {
        int pos = var.first;
        BOOST_TEST(expOut[pos][0] == vcfOut[pos][0]);
        BOOST_TEST(expOut[pos][1] == vcfOut[pos][1]);
    }
}
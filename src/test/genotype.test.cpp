#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <stack>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include "PangenomeMAT.hpp"
#include "../genotype.hpp"
#include "../tree.hpp"

using namespace std;
using namespace tree;
namespace fs = boost::filesystem;


BOOST_AUTO_TEST_CASE(readMutationMatrices) {
    string testDataDir = "../src/test/data/genotype_test_data/";
    
    // case 1: full mutation matrix + new lines
    for (const string& mmPath : {testDataDir + "test.mm", testDataDir + "test2.mm"}) {
        mutationMatrices mutMat = mutationMatrices();
        ifstream mminf(mmPath);
        fillMutationMatricesFromFile(mutMat, mminf);
        
        BOOST_TEST(mutMat.submat.size() == 4);
        for (const vector<double>& row : mutMat.submat) {
            BOOST_TEST(row.size() == 4);
        }
        BOOST_TEST(mutMat.insmat.size() == 10);
        BOOST_TEST(mutMat.delmat.size() == 10);
    }

    // case 2: missing elements
    for (const string& mmPath : {testDataDir + "test3.mm", testDataDir + "test4.mm"}) {
        mutationMatrices mutMat = mutationMatrices();
        ifstream mminf(mmPath);
        BOOST_CHECK_THROW(fillMutationMatricesFromFile(mutMat, mminf), invalid_argument);
    }
    
}

BOOST_AUTO_TEST_CASE(genotyping) {
    cout << "\n\nCreate vcf tests (using ../src/test/data/genotype_test_data/test.pileup)" << endl;
    ifstream puif("../src/test/data/genotype_test_data/test.pileup");
    /*############################
    # genotype likelihoods alone #
    ############################*/
    cout << "Testing pileup to vcf without prior (using ../src/test/data/genotype_test_data/test_empty.mm)" << endl;
    ifstream emif("../src/test/data/genotype_test_data/test_empty.mm");
    mutationMatrices emptyMutMat = mutationMatrices();
    fillMutationMatricesFromFile(emptyMutMat, emif);
    emif.close();

    ofstream vcfof("../src/test/data/genotype_test_data/vcf_np.out.tmp");
    genotype::printSamplePlacementVCF(puif, emptyMutMat, false, 0, vcfof);
    vcfof.close();
    puif.clear();
    puif.seekg(0, ios::beg);

    cout << "Wrote vcf file to ../src/test/data/genotype_test_data/vcf_np.out.tmp" << endl;

    unordered_map<int, vector<string> > vcfOut;
    ifstream vcfif("../src/test/data/genotype_test_data/vcf_np.out.tmp");
    string line;
    while (getline(vcfif, line)) {
        if (line.substr(0, 3) != "ref") {continue;}
        vector<string> fields;
        stringSplit(line, '\t', fields);
        vcfOut[stoi(fields[1])] = {fields[4], fields.back()};
    }
    vcfif.close();

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

    cout << "Finished subtest.. Deleting ../src/test/data/genotype_test_data/vcf_np.out.tmp\n" << endl;
    fs::remove("../src/test/data/genotype_test_data/vcf_np.out.tmp");

    /*##############################
    # genotype likelihoods + Prior #
    ##############################*/
    cout << "Testing pileup to vcf with prior (using ../src/test/data/genotype_test_data/test.mm)" << endl;
    emif.open("../src/test/data/genotype_test_data/test.mm");
    mutationMatrices mutmat = mutationMatrices();
    fillMutationMatricesFromFile(mutmat, emif);
    emif.close();

    vcfof.open("../src/test/data/genotype_test_data/vcf_wp.out.tmp");
    genotype::printSamplePlacementVCF(puif, mutmat, false, 0, vcfof);
    vcfof.close();
    puif.close();

    cout << "Wrote vcf file to ../src/test/data/genotype_test_data/vcf_wp.out.tmp" << endl;
    
    vcfOut.clear();
    vcfif.open("../src/test/data/genotype_test_data/vcf_wp.out.tmp");
    line.clear();
    while (getline(vcfif, line)) {
        if (line.substr(0, 3) != "ref") {continue;}
        vector<string> fields;
        stringSplit(line, '\t', fields);
        vcfOut[stoi(fields[1])] = {fields[4], fields.back()};
    }
    vcfif.close();

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

    cout << "Finished subtest.. Deleting ../src/test/data/genotype_test_data/vcf_wp.out.tmp\n" << endl;
    fs::remove("../src/test/data/genotype_test_data/vcf_wp.out.tmp");
}
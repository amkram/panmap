#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
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

size_t getStart_brute(const std::string& s1, const std::string& s2, size_t window, double threshold);
size_t getEnd_brute(const std::string& s1, const std::string& s2, size_t window, double threshold);
size_t getBeg_perfectMatch(const std::string& s1, const std::string& s2, size_t matchLen);
size_t getEnd_perfectMatch(const std::string& s1, const std::string& s2, size_t matchLen);

BOOST_AUTO_TEST_CASE(edgeMaskingForBuildingMutMat) {
    std::ifstream ifs("../dev/examples/sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
    b.push(boost::iostreams::gzip_decompressor());
    b.push(ifs);
    std::istream is(&b);

    auto T = new PangenomeMAT::Tree(is);

    std::unordered_map<std::string, std::string> alignedSequences = getAllNodeStrings(T);
    for (const auto& sequence : alignedSequences) {
        std::string parentId;
        if (T->allNodes[sequence.first]->parent == nullptr) {
            continue;
        } else {
            parentId = T->allNodes[sequence.first]->parent->identifier;
        }
        const std::string& curSeq = sequence.second;
        const std::string& parSeq = alignedSequences[parentId];
        size_t bruteBeg = getStart_brute(curSeq, parSeq, 20, 0.8);
        size_t bruteEnd = getEnd_brute(curSeq, parSeq, 20, 0.8);
        std::pair<size_t, size_t> maskCoors = tree::getMaskCoorsForMutmat(curSeq, parSeq, 20, 0.8);
        BOOST_TEST(bruteBeg == maskCoors.first);
        BOOST_TEST(bruteEnd == maskCoors.second);
    }
}

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

size_t getStart_brute(const std::string& s1, const std::string& s2, size_t window, double threshold) {
    assert(s1.size() == s2.size());
    if (s1 == "") {
        return 0;
    }

    for (size_t i = 0; i < s1.size() - window + 1; i++) {
        std::string sub1 = "";
        std::string sub2 = "";
        size_t off = 0;
        while (sub1.size() < window) {
            size_t idx = i + off;
            if (idx >= s1.size()) {
                break;
            }
            if (s1[idx] == '-' && s2[idx] == '-') {
                off++;
                continue;
            } else {
                sub1 += s1[idx];
                sub2 += s2[idx];
            }
            off++;
        }
        if (sub1.size() < window) {
            continue;
        }

        double numMatch = 0.0;
        for (size_t i = 0; i < sub1.size(); i++) {
            if (sub1[i] == sub2[i]) {
                numMatch += 1.0;
            }
        }
        double pcid = numMatch / static_cast<double>(window);
        
        if (pcid >= threshold && sub1[0] == sub2[0]) {
            size_t nucOff = 0;
            for (size_t j = i + nucOff; j < s1.size(); j++) {
                if (s1[j] != '-' || s2[j] != '-') {
                    break;
                }
                ++nucOff;
            }
            return i + nucOff;
        }
    }

    return s1.size() - 1;
}



size_t getEnd_brute(const std::string& s1, const std::string& s2, size_t window, double threshold) {
    assert(s1.size() == s2.size());
    if (s1.empty()) {
        return 0;
    }

    for (size_t i = s1.size(); i > window - 1; i--) {
        std::string sub1 = "";
        std::string sub2 = "";
        size_t off = 0;
        while (sub1.size() < window) {
            size_t idx = i - 1 - off;
            if (idx >= s1.size()) {
                break;
            }
            if (s1[idx] == '-' && s2[idx] == '-') {
                off++;
                continue;
            } else {
                sub1 = s1[idx] + sub1;
                sub2 = s2[idx] + sub2;
            }
            off++;
        }
        if (sub1.size() < window) {
            continue;
        }

        double numMatch = 0.0;
        for (size_t j = 0; j < sub1.size(); j++) {
            if (sub1[j] == sub2[j]) {
                numMatch += 1.0;
            }
        }
        double pcid = numMatch / static_cast<double>(window);
        
        if (pcid >= threshold && sub1[window-1] == sub2[window-1]) {
            size_t nucOff = 0;
            for (size_t j = i - 1 - nucOff; j > 0; j--) {
                if (s1[j] != '-' || s2[j] != '-') {
                    break;
                }
                ++nucOff;
            }
            return i - 1 - nucOff;
        }
    }

    return 0;
}

size_t getBeg_perfectMatch(const std::string& s1, const std::string& s2, size_t matchLen) {
    if (s1 == "" || matchLen == 0) {
        return 0;
    }
    size_t start = 0;
    size_t numMatch = 0;
    for (size_t i = 0; i < s1.size(); i++) {
        if (s1[i] == '-' && s2[i] == '-') {
            continue;
        } else if (s1[i] == s2[i]) {
            numMatch++;
            if (numMatch == matchLen) {
                break;
            }
        } else {
            numMatch = 0;
            start = i + 1;
        }
    }

    if (numMatch < matchLen) {
        return s1.size() - 1;
    }
    return start;
}

size_t getEnd_perfectMatch(const std::string& s1, const std::string& s2, size_t matchLen) {
    if (s1 == "") {
        return 0;
    } else if (matchLen == 0) {
        return s1.size() - 1;
    }
    size_t end = s1.size() - 1;
    size_t numMatch = 0;
    size_t i;
    for (size_t j = 0; j < s1.size(); j++) {
        i = s1.size() - 1 - j;
        if (s1[i] == '-' && s2[i] == '-') {
            continue;
        } else if (s1[i] == s2[i]) {
            numMatch++;
            if (numMatch == matchLen) {
                break;
            }
        } else {
            numMatch = 0;
            end = i - 1;
        }
    }

    if (numMatch < matchLen) {
        return 0;
    }
    return end;
}




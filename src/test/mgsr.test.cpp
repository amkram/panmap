#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <iostream>
#include "../mgsr.hpp"

using namespace std;
namespace fs = boost::filesystem;

void makeFasta(const std::string& name, const std::string& seq, const std::string& path);

// BOOST_AUTO_TEST_CASE(seqs) {
//     ifstream ifs("../dev/examples/sars2k.pmat");
//     boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
//     b.push(boost::iostreams::gzip_decompressor());
//     b.push(ifs);
//     istream is(&b);
//     auto T = new PangenomeMAT::Tree(is);

//     std::unordered_map<std::string, std::string> alignedSequences = getAllNodeStrings(T);

//     for (const auto& sequence : alignedSequences) {
//         std::string path = "../src/test/data/aligned_sequences/";
//         makeFasta(sequence.first, sequence.second, path + sequence.first + ".fa");
//     }

// }

BOOST_AUTO_TEST_CASE(tmp) {
    ifstream ifs("../dev/examples/sars2k.pmat");
    boost::iostreams::filtering_streambuf< boost::iostreams::input> b;
    b.push(boost::iostreams::gzip_decompressor());
    b.push(ifs);
    istream is(&b);
    auto T = new PangenomeMAT::Tree(is);

    string indexFilePath = "../dev/examples/sars2k.pmat.spmi";
    string seedmersIndexPath = "../dev/examples/sars2k.pmat.kmi";
    size_t k = 10;
    size_t s = 5;
    size_t l = 2;

    pmi::seedIndex index;
    std::stringstream seedmersOutStream;

    mgsr::buildSeedmer(index, T, l, k, s, seedmersOutStream);
    std::cout << "Writing to " << indexFilePath << "..." << std::endl;
    std::ofstream fout(indexFilePath);
    std::ofstream smfout(seedmersIndexPath);
    fout << index.outStream.str();
    smfout << seedmersOutStream.str();
    std::cout << "Done." << std::endl;

    ifstream indexFile(indexFilePath);

    mgsr::accio(T, indexFile, k, l);
}

void makeFasta(const std::string& name, const std::string& seq, const std::string& path) {
    if (!fs::exists(path)) {
        std::ofstream faos(path);
        faos << '>' << name << '\n';
        size_t linesize = 80;
        for (size_t i = 0; i < seq.size(); i += linesize) {
            faos << seq.substr(i, std::min(linesize, seq.size() - i)) << '\n';
        }
        faos.close();
    }
}
// merge_indexes: merge N single-tree panmap MGSR lite indexes (.idx) into ONE combined
// LiteIndex spanning all clades, so `panmap --meta -i merged.idx reads.fq` seeds reads ONCE
// and places across all clades. Node ids are namespaced "<clade>::<id>" (clade = filename stem)
// to avoid collisions. Core merge logic lives in merge_lite_index.hpp (shared with buildMgsrIndex).
//
// Usage: merge_indexes OUT.idx IN1.idx [IN2.idx ...]
#include "merge_lite_index.hpp"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

static std::string cladeStem(const std::string& p) {
    std::string s = std::filesystem::path(p).filename().string();
    for (const std::string suf : {".panman.idx", ".idx", ".panman"}) {
        if (s.size() > suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0) {
            return s.substr(0, s.size() - suf.size());
        }
    }
    return s;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: merge_indexes OUT.idx IN1.idx [IN2.idx ...]\n";
        return 1;
    }
    const std::string outPath = argv[1];
    std::vector<std::string> inputs(argv + 2, argv + argc);
    std::vector<std::string> prefixes;
    prefixes.reserve(inputs.size());
    for (const auto& in : inputs) prefixes.push_back(cladeStem(in) + "::");

    try {
        mergeLiteIndexes(inputs, prefixes, outPath);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    std::cerr << "merged " << inputs.size() << " indexes -> " << outPath << "\n";
    return 0;
}

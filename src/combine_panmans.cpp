// combine_panmans: merge N single-tree .panman files into ONE multi-PanMAT PanMAN
// (a TreeGroup whose `trees` vector holds N disconnected PanMATs, zero complex mutations).
//
// Node identifiers are namespaced "<clade>::<id>" (clade = input filename stem) so that
// internal node ids (node_1, node_2, ...), which collide across clades, become globally
// unique — required because panmap's MgsrLiteTree::initialize exits on duplicate ids.
//
// Mirrors panmap's loadPanMAN (lzma-decompress -> TreeGroup(istream)) and panmanUtils'
// writePanMAN (lzma level-9 -> TreeGroup::writeToFile).
//
// Usage: combine_panmans OUT.panman IN1.panman [IN2.panman ...]
#include "panman.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/copy.hpp>
#include <kj/std/iostream.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

namespace bio = boost::iostreams;

static std::string cladeStem(const std::string& path) {
    return std::filesystem::path(path).stem().string();  // e.g. "clade_001"
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: combine_panmans OUT.panman IN1.panman [IN2.panman ...]\n";
        return 1;
    }
    const std::string outPath = argv[1];
    std::vector<std::string> inputs(argv + 2, argv + argc);

    // Per-file TreeGroups must outlive the Tree* collection (TreeGroup(vector<Tree*>)
    // copies Tree structs but the Node* they hold are owned by these holders).
    std::vector<std::unique_ptr<panmanUtils::TreeGroup>> holders;
    holders.reserve(inputs.size());
    std::vector<panmanUtils::Tree*> treePtrs;
    treePtrs.reserve(inputs.size());

    size_t totalNodes = 0;
    for (const auto& in : inputs) {
        const std::string clade = cladeStem(in);
        std::ifstream file(in, std::ios::binary);
        if (!file) { std::cerr << "ERROR: cannot open " << in << "\n"; return 1; }

        auto buf = std::make_unique<bio::filtering_streambuf<bio::input>>();
        buf->push(bio::lzma_decompressor());
        buf->push(file);
        std::istream stream(buf.get());

        auto tg = std::make_unique<panmanUtils::TreeGroup>(stream);
        if (tg->trees.empty()) { std::cerr << "ERROR: no trees in " << in << "\n"; return 1; }
        panmanUtils::Tree& T = tg->trees[0];

        // Namespace every node id: "<clade>::<oldid>". Rewrite Node::identifier (same Node*
        // is referenced by children pointers, so structure is preserved) and rebuild allNodes.
        std::unordered_map<std::string, panmanUtils::Node*> newMap;
        newMap.reserve(T.allNodes.size());
        for (auto& kv : T.allNodes) {
            panmanUtils::Node* n = kv.second;
            n->identifier = clade + "::" + n->identifier;
            newMap[n->identifier] = n;
        }
        T.allNodes = std::move(newMap);

        totalNodes += T.allNodes.size();
        treePtrs.push_back(&T);
        holders.push_back(std::move(tg));
    }
    std::cerr << "loaded " << treePtrs.size() << " trees, " << totalNodes << " total nodes\n";

    panmanUtils::TreeGroup combined(treePtrs);
    std::cerr << "combined TreeGroup: " << combined.trees.size() << " trees, "
              << combined.complexMutations.size() << " complex mutations\n";

    std::ofstream outputFile(outPath, std::ios::binary);
    if (!outputFile) { std::cerr << "ERROR: cannot write " << outPath << "\n"; return 1; }
    bio::filtering_streambuf<bio::output> outBuf;
    bio::lzma_params params;
    params.level = 9;
    outBuf.push(bio::lzma_compressor(params));
    outBuf.push(outputFile);
    std::ostream outstream(&outBuf);
    kj::std::StdOutputStream outputStream(outstream);
    combined.writeToFile(outputStream);
    bio::close(outBuf);
    outputFile.close();

    std::cerr << "wrote " << outPath << "\n";
    return 0;
}

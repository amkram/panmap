// merge_indexes: merge N single-tree panmap MGSR lite indexes (.idx) into ONE combined
// LiteIndex spanning all clades, so `panmap --meta -i merged.idx reads.fq` seeds reads ONCE
// and places across all clades (eliminating the N-fold read-seeding cost of N separate runs).
//
// Why this works (verified against MgsrLiteTree::initialize / mgsr.cpp):
//   * lite index uses seedInfos (seedHashes + perNodeChanges); seedChange*/nodeChangeOffsets
//     are empty (initialize exit(1)s if both seedInfos and refSeedChanges are present).
//   * perNodeChanges[i] aligns positionally with liteNodes[i]; parentIndex is a positional
//     index; liteNodes[0] is the root; nodes are preorder (parent before child).
//   * seedDeltaIndices index into the global seedHashes array.
// Merge layout: liteNodes[0] = synthetic super-root; then each clade's nodes appended in order
// with parentIndex offset by the clade's node base (clade root -> super-root 0), seedDeltaIndices
// offset by the clade's seed base, and node ids namespaced "<clade>::<id>" (clade = filename stem)
// because initialize exit(1)s on duplicate ids and internal ids (node_1...) collide across clades.
//
// Usage: merge_indexes OUT.idx IN1.idx [IN2.idx ...]
#include "index_lite.capnp.h"

#include <capnp/message.h>
#include <capnp/serialize.h>

#include <fcntl.h>
#include <unistd.h>

#include <filesystem>
#include <iostream>
#include <memory>
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

    // Open all inputs; keep fds + readers alive while copying (capnp Readers reference them).
    std::vector<int> fds;
    std::vector<std::unique_ptr<capnp::StreamFdMessageReader>> readers;
    std::vector<LiteIndex::Reader> idxs;
    std::vector<std::string> clades;

    size_t totalSeeds = 0, totalNodes = 1;  // +1 synthetic super-root
    capnp::ReaderOptions opt;
    opt.traversalLimitInWords = (uint64_t)1 << 40;  // these per-clade indexes are small

    for (const auto& in : inputs) {
        int fd = ::open(in.c_str(), O_RDONLY);
        if (fd < 0) { std::cerr << "ERROR: cannot open " << in << "\n"; return 1; }
        auto r = std::make_unique<capnp::StreamFdMessageReader>(fd, opt);
        LiteIndex::Reader idx = r->getRoot<LiteIndex>();
        if (idx.getLiteTree().getLiteNodes().size() != idx.getPerNodeChanges().size()) {
            std::cerr << "ERROR: liteNodes/perNodeChanges size mismatch in " << in << "\n";
            return 1;
        }
        totalSeeds += idx.getSeedHashes().size();
        totalNodes += idx.getLiteTree().getLiteNodes().size();
        clades.push_back(cladeStem(in));
        idxs.push_back(idx);
        readers.push_back(std::move(r));
        fds.push_back(fd);
    }

    capnp::MallocMessageBuilder msg;
    LiteIndex::Builder out = msg.initRoot<LiteIndex>();
    out.setK(idxs[0].getK());
    out.setS(idxs[0].getS());
    out.setT(idxs[0].getT());
    out.setL(idxs[0].getL());
    out.setOpen(idxs[0].getOpen());
    out.setHpc(idxs[0].getHpc());

    auto seedHashesB = out.initSeedHashes(totalSeeds);
    auto seedRevB = out.initSeedIsReverse(totalSeeds);
    auto liteTreeB = out.initLiteTree();
    auto liteNodesB = liteTreeB.initLiteNodes(totalNodes);
    auto perNodeB = out.initPerNodeChanges(totalNodes);

    // Synthetic super-root at index 0: no seeds, parent = itself (0), root of the forest.
    liteNodesB[0].setId("__superroot__");
    liteNodesB[0].setParentIndex(0);
    liteNodesB[0].setIdenticalToParent(false);
    perNodeB[0].setNodeIndex(0);
    perNodeB[0].initSeedDeltaIndices(0);
    perNodeB[0].initSeedDeltaIsDeleted(0);
    perNodeB[0].initGapRunDeltas(0);
    perNodeB[0].initInvertedBlocks(0);

    uint32_t seedOff = 0, nodeOff = 1;
    for (size_t c = 0; c < idxs.size(); c++) {
        auto idx = idxs[c];
        auto sh = idx.getSeedHashes();
        auto sr = idx.getSeedIsReverse();
        for (uint32_t i = 0; i < sh.size(); i++) {
            seedHashesB.set(seedOff + i, sh[i]);
            seedRevB.set(seedOff + i, sr[i]);
        }
        auto ln = idx.getLiteTree().getLiteNodes();
        auto pnc = idx.getPerNodeChanges();
        for (uint32_t j = 0; j < ln.size(); j++) {
            uint32_t g = nodeOff + j;
            liteNodesB[g].setId(clades[c] + "::" + std::string(ln[j].getId().cStr()));
            liteNodesB[g].setIdenticalToParent(ln[j].getIdenticalToParent());
            // clade root (local index 0) -> super-root 0; else offset local parentIndex.
            liteNodesB[g].setParentIndex(j == 0 ? 0u : nodeOff + ln[j].getParentIndex());

            auto src = pnc[j];
            auto dst = perNodeB[g];
            dst.setNodeIndex(g);
            auto sdi = src.getSeedDeltaIndices();
            auto sdd = src.getSeedDeltaIsDeleted();
            auto sdiB = dst.initSeedDeltaIndices(sdi.size());
            auto sddB = dst.initSeedDeltaIsDeleted(sdd.size());
            for (uint32_t x = 0; x < sdi.size(); x++) {
                sdiB.set(x, seedOff + sdi[x]);  // shift into global seedHashes
                sddB.set(x, sdd[x]);
            }
            auto grd = src.getGapRunDeltas();
            auto grdB = dst.initGapRunDeltas(grd.size());
            for (uint32_t x = 0; x < grd.size(); x++) {
                grdB[x].setStartPos(grd[x].getStartPos());
                grdB[x].setEndPos(grd[x].getEndPos());
                grdB[x].setToGap(grd[x].getToGap());
            }
            auto ib = src.getInvertedBlocks();
            auto ibB = dst.initInvertedBlocks(ib.size());
            for (uint32_t x = 0; x < ib.size(); x++) ibB.set(x, ib[x]);
        }
        seedOff += (uint32_t)sh.size();
        nodeOff += (uint32_t)ln.size();
    }
    // seedChangeHashes / nodeChangeOffsets / seedStartPos / seedEndPos / substitutionMatrix:
    // intentionally left empty (lite mode uses seedInfos path).

    int ofd = ::open(outPath.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (ofd < 0) { std::cerr << "ERROR: cannot write " << outPath << "\n"; return 1; }
    capnp::writeMessageToFd(ofd, msg);
    ::close(ofd);
    for (int fd : fds) ::close(fd);

    std::cerr << "merged " << idxs.size() << " indexes -> " << outPath
              << "  (" << totalNodes << " nodes incl super-root, " << totalSeeds << " seeds)\n";
    return 0;
}

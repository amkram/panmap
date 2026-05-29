// Merge N single-tree MGSR lite indexes into one combined LiteIndex (a forest under a
// synthetic super-root at node 0). Shared by the standalone `merge_indexes` tool and by
// `buildMgsrIndex` (multi-tree panmans). See MULTIPANMAT.md.
//
// Invariants honored (from MgsrLiteTree::initialize): perNodeChanges[i] aligns positionally
// with liteNodes[i]; parentIndex is a positional index; node 0 is the root; nodes are preorder;
// seedDeltaIndices index into the global seedHashes array; node ids must be unique (initialize
// exit(1)s on duplicates). Lite mode leaves seedChange*/nodeChangeOffsets empty.
#pragma once

#include "index_lite.capnp.h"

#include <capnp/message.h>
#include <capnp/serialize.h>

#include <fcntl.h>
#include <unistd.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// inPaths[i]   : path to a single-tree lite .idx (capnp, non-packed)
// prefixes[i]  : string prepended to every node id of input i ("" = leave ids unchanged).
//                Use to namespace colliding internal ids (node_1, ...) across clades.
inline void mergeLiteIndexes(const std::vector<std::string>& inPaths,
                             const std::vector<std::string>& prefixes,
                             const std::string& outPath) {
    if (inPaths.empty()) throw std::runtime_error("mergeLiteIndexes: no inputs");
    if (prefixes.size() != inPaths.size())
        throw std::runtime_error("mergeLiteIndexes: prefixes/inPaths size mismatch");

    std::vector<int> fds;
    std::vector<std::unique_ptr<capnp::StreamFdMessageReader>> readers;
    std::vector<LiteIndex::Reader> idxs;

    size_t totalSeeds = 0, totalNodes = 1;  // +1 synthetic super-root
    capnp::ReaderOptions opt;
    opt.traversalLimitInWords = (uint64_t)1 << 42;

    for (const auto& in : inPaths) {
        int fd = ::open(in.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("mergeLiteIndexes: cannot open " + in);
        auto r = std::make_unique<capnp::StreamFdMessageReader>(fd, opt);
        LiteIndex::Reader idx = r->getRoot<LiteIndex>();
        if (idx.getLiteTree().getLiteNodes().size() != idx.getPerNodeChanges().size())
            throw std::runtime_error("mergeLiteIndexes: liteNodes/perNodeChanges mismatch in " + in);
        totalSeeds += idx.getSeedHashes().size();
        totalNodes += idx.getLiteTree().getLiteNodes().size();
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

    // Synthetic super-root at index 0: no seeds, parent = itself.
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
        const std::string& pre = prefixes[c];
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
            liteNodesB[g].setId(pre + std::string(ln[j].getId().cStr()));
            liteNodesB[g].setIdenticalToParent(ln[j].getIdenticalToParent());
            liteNodesB[g].setParentIndex(j == 0 ? 0u : nodeOff + ln[j].getParentIndex());

            auto src = pnc[j];
            auto dst = perNodeB[g];
            dst.setNodeIndex(g);
            auto sdi = src.getSeedDeltaIndices();
            auto sdd = src.getSeedDeltaIsDeleted();
            auto sdiB = dst.initSeedDeltaIndices(sdi.size());
            auto sddB = dst.initSeedDeltaIsDeleted(sdd.size());
            for (uint32_t x = 0; x < sdi.size(); x++) {
                sdiB.set(x, seedOff + sdi[x]);
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

    int ofd = ::open(outPath.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (ofd < 0) throw std::runtime_error("mergeLiteIndexes: cannot write " + outPath);
    capnp::writeMessageToFd(ofd, msg);
    ::close(ofd);
    for (int fd : fds) ::close(fd);
}

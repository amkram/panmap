#include "test_index.hpp"

#include "index_single_mode.hpp"
#include "zstd_compression.hpp"

#include "capnp/message.h"
#include "capnp/serialize.h"
#include "index_lite.capnp.h"

#include <unistd.h>

#include <atomic>
#include <filesystem>
#include <stdexcept>

namespace ts {

namespace {
std::string makeTempIndexPath() {
    static std::atomic<uint64_t> counter{0};
    auto name = "panmap_testidx_" + std::to_string(::getpid()) + "_" + std::to_string(counter.fetch_add(1)) + ".idx";
    return (std::filesystem::temp_directory_path() / name).string();
}
}

IndexData loadIndex(const std::string& path) {
    std::vector<uint8_t> buffer;
    // Skip the uncompressed param header prepended by writeIndex (see readIndexHeader).
    index_single_mode::IndexParamsHeader ph;
    const size_t dataOffset =
        index_single_mode::readIndexHeader(path, ph) ? index_single_mode::kIndexHeaderSize : 0;
    if (!panmap_zstd::decompressFromFile(path, buffer, 0, dataOffset)) {
        throw std::runtime_error("Failed to decompress index: " + path);
    }

    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 1ULL << 50;  // trusted test index; effectively unlimited
    opts.nestingLimit = 1024;

    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(reinterpret_cast<const ::capnp::word*>(buffer.data()),
                                          buffer.size() / sizeof(::capnp::word)),
        opts);

    auto indexRoot = reader.getRoot<LiteIndex>();
    IndexData d;
    d.k = indexRoot.getK();
    d.s = indexRoot.getS();
    d.t = indexRoot.getT();
    d.l = indexRoot.getL();
    d.open = indexRoot.getOpen();
    d.hpc = indexRoot.getHpc();

    // Flatten segmented seed-change struct-of-arrays into vectors.
    auto offsetsReader = indexRoot.getNodeChangeOffsets();
    d.offsets.resize(offsetsReader.size());
    for (size_t i = 0; i < offsetsReader.size(); ++i) d.offsets[i] = offsetsReader[i];

    const size_t total = d.offsets.empty() ? 0 : d.offsets.back();
    d.hashes.reserve(total);
    d.parentCounts.reserve(total);
    d.childCounts.reserve(total);

    auto hashSegs = indexRoot.getSeedChangeHashes();
    auto parentSegs = indexRoot.getSeedChangeParentCounts();
    auto childSegs = indexRoot.getSeedChangeChildCounts();
    for (uint32_t seg = 0; seg < hashSegs.size(); ++seg) {
        auto hs = hashSegs[seg];
        auto ps = parentSegs[seg];
        auto cs = childSegs[seg];
        for (uint32_t j = 0; j < hs.size(); ++j) {
            d.hashes.push_back(hs[j]);
            d.parentCounts.push_back(ps[j]);
            d.childCounts.push_back(cs[j]);
        }
    }

    d.liteTree = std::make_unique<panmapUtils::LiteTree>();
    d.liteTree->initialize(indexRoot.getLiteTree());
    return d;
}

TestIndex::TestIndex(panmanUtils::Tree* tree,
                     int k,
                     int s,
                     int t,
                     int l,
                     bool open,
                     bool hpc,
                     int flankMaskBp,
                     int numThreads,
                     int zstdLevel)
    : path_(makeTempIndexPath()) {
    index_single_mode::IndexBuilder builder(tree, k, s, t, l, open, flankMaskBp, hpc);
    if (numThreads > 1) {
        builder.buildIndexParallel(numThreads);
    } else {
        builder.buildIndex();
    }
    builder.writeIndex(path_, /*numThreads=*/numThreads > 0 ? numThreads : 0, zstdLevel);

    data_ = loadIndex(path_);
    data_.flankMaskBp = flankMaskBp;
}

TestIndex::~TestIndex() {
    std::error_code ec;
    std::filesystem::remove(path_, ec);
}

const RSVPanmanFixture& sharedRSVFixture() {
    static const RSVPanmanFixture fx;   // rsv_4K loaded once for the whole binary
    return fx;
}

const IndexData& sharedRSVIndex() {
    // Built once; index build is read-only on the tree, so the shared fixture stays reusable.
    static TestIndex idx(sharedRSVFixture().tree(), /*k=*/15, /*s=*/8, /*t=*/0, /*l=*/1);
    return idx.data();
}

}  // namespace ts

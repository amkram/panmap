// Unit tests for core data-structure utilities: zstd round-trip and LiteTree.
#include <boost/test/unit_test.hpp>

#include "helpers/test_index.hpp"
#include "helpers/tree_helpers.hpp"
#include "zstd_compression.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <random>
#include <unistd.h>
#include <vector>

namespace {
std::string tmpFile(const char* tag) {
    static std::atomic<uint64_t> ctr{0};
    auto name = std::string("panmap_test_") + tag + "_" + std::to_string(::getpid()) + "_" +
                std::to_string(ctr.fetch_add(1)) + ".zst";
    return (std::filesystem::temp_directory_path() / name).string();
}

std::vector<uint8_t> randomBytes(size_t n, uint64_t seed) {
    std::mt19937_64 gen(seed);
    std::vector<uint8_t> v(n);
    for (auto& b : v) b = static_cast<uint8_t>(gen() & 0xff);
    return v;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(core_utils_tests)

BOOST_AUTO_TEST_CASE(zstd_roundtrip_identity) {
    for (size_t n : {size_t{1024}, size_t{5u * 1024 * 1024}}) {
        auto data = randomBytes(n, 0xC0FFEE + n);
        std::string path = tmpFile("zstd");

        BOOST_REQUIRE(panmap_zstd::compressToFile(data.data(), data.size(), path, /*level=*/3));

        std::vector<uint8_t> out1, out4;
        BOOST_REQUIRE(panmap_zstd::decompressFromFile(path, out1, /*threads=*/1));
        BOOST_REQUIRE(panmap_zstd::decompressFromFile(path, out4, /*threads=*/4));

        BOOST_TEST(out1.size() == data.size());
        BOOST_TEST((out1 == data));
        BOOST_TEST((out4 == out1));  // thread count must not affect decompressed bytes

        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
}

// LiteTree structural invariants, validated on a real index built from rsv_4K.
BOOST_AUTO_TEST_CASE(litetree_structure_invariants) {
    ts::RSVPanmanFixture fixture;
    ts::TestIndex idx(fixture.tree(), /*k=*/15, /*s=*/8, /*t=*/0, /*l=*/1);
    const auto& tree = *idx.data().liteTree;

    BOOST_REQUIRE(tree.root != nullptr);
    BOOST_TEST(tree.root->parent == nullptr);
    BOOST_TEST(!tree.root->children.empty());
    BOOST_TEST(tree.dfsIndexToNode.size() == tree.allLiteNodes.size());

    size_t withParent = 0;
    for (size_t i = 0; i < tree.dfsIndexToNode.size(); ++i) {
        auto* node = tree.dfsIndexToNode[i];
        BOOST_REQUIRE(node != nullptr);  // no holes in the dfs index
        if (node->parent != nullptr) {
            withParent++;
            // parent<->child links agree
            auto& sibs = node->parent->children;
            BOOST_TEST((std::find(sibs.begin(), sibs.end(), node) != sibs.end()));
        }
        BOOST_TEST(tree.resolveNodeId(static_cast<uint32_t>(i)) == node->identifier);
    }
    // Exactly one node (the root) has no parent.
    BOOST_TEST(withParent == tree.allLiteNodes.size() - 1);

    // Out-of-range index resolves to empty string.
    BOOST_TEST(tree.resolveNodeId(static_cast<uint32_t>(tree.dfsIndexToNode.size())) == "");
}

BOOST_AUTO_TEST_SUITE_END()

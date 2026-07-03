// Index round-trip tests: structural integrity + delta reconstruction == direct
// extraction. This file is the SOLE owner of the "reconstruction equals direct
// genome seed extraction" contract; test_placement assumes it.
#include <boost/test/unit_test.hpp>

#include "helpers/seed_helpers.hpp"
#include "helpers/test_index.hpp"
#include "helpers/traversal.hpp"
#include "helpers/tree_helpers.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <vector>

namespace {
// Build params used wherever we compare against raw rollingSyncmers extraction:
// l=1 (kminmers == syncmers) and flankMaskBp=0 (no masking), matching extractSeeds.
constexpr int K = 15, S = 8, L = 1;

std::vector<std::string> sampleNodeIds(const ts::RSVPanmanFixture& fix, int n, uint64_t seed) {
    auto all = fix.nodeIds();
    std::mt19937_64 gen(seed);
    std::shuffle(all.begin(), all.end(), gen);
    if (all.size() > static_cast<size_t>(n)) all.resize(n);
    return all;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(index_tests)

BOOST_AUTO_TEST_CASE(index_parameters_persisted) {
    ts::RSVPanmanFixture fix;
    ts::TestIndex idx(fix.tree(), /*k=*/19, /*s=*/6, /*t=*/0, /*l=*/3, /*open=*/false, /*hpc=*/false);
    auto& d = idx.data();
    BOOST_TEST(d.k == 19);
    BOOST_TEST(d.s == 6);
    BOOST_TEST(d.t == 0);
    BOOST_TEST(d.l == 3);
    BOOST_TEST(d.open == false);
    BOOST_TEST(d.hpc == false);
}

BOOST_AUTO_TEST_CASE(offset_and_parent_child_consistency) {
    ts::RSVPanmanFixture fix;
    ts::TestIndex idx(fix.tree(), K, S, 0, L);
    auto& d = idx.data();

    const size_t numOffsets = d.numNodesPlusOne();
    BOOST_REQUIRE(numOffsets >= 2);
    BOOST_TEST(d.nodeChangeOffset(0) == 0u);
    BOOST_TEST(d.nodeChangeOffset(numOffsets - 1) == d.numChanges());

    // Monotonic non-decreasing offsets.
    for (size_t i = 1; i < numOffsets; ++i) {
        BOOST_REQUIRE(d.nodeChangeOffset(i) >= d.nodeChangeOffset(i - 1));
    }

    // Along several root->node paths, each change's parentCount must equal the
    // running (child) count established by the previous change for that seed.
    const auto& tree = *d.liteTree;
    int inconsistencies = 0;
    for (const auto& nodeId : sampleNodeIds(fix, 8, 13)) {
        int32_t ni = ts::findNodeIndex(tree, nodeId);
        if (ni < 0) continue;
        absl::flat_hash_map<uint64_t, int64_t> running;
        for (auto* node : ts::pathToRoot(tree, ni)) {
            uint64_t start = d.nodeChangeOffset(node->nodeIndex);
            uint64_t end = d.nodeChangeOffset(node->nodeIndex + 1);
            for (uint64_t j = start; j < end; ++j) {
                uint64_t h = d.hashAt(j);
                int64_t parent = d.parentCountAt(j);
                int64_t child = d.childCountAt(j);
                int64_t expected = running.count(h) ? running[h] : 0;
                if (parent != expected) inconsistencies++;
                if (child > 0)
                    running[h] = child;
                else
                    running.erase(h);
            }
        }
    }
    BOOST_TEST(inconsistencies == 0);
}

BOOST_AUTO_TEST_CASE(delta_reconstruction_equals_direct) {
    ts::RSVPanmanFixture fix;
    ts::TestIndex idx(fix.tree(), K, S, 0, L);
    auto& d = idx.data();
    const auto& tree = *d.liteTree;

    std::vector<std::string> targets = {tree.root->identifier};
    for (const auto& id : sampleNodeIds(fix, 5, 42)) targets.push_back(id);

    for (const auto& nodeId : targets) {
        int32_t ni = ts::findNodeIndex(tree, nodeId);
        BOOST_REQUIRE(ni >= 0);
        auto reconstructed = ts::reconstructGenomeSeeds(ts::pathToRoot(tree, ni), d);

        std::string genome = fix.genomeOf(nodeId);
        auto direct = ts::extractSeeds(genome, K, S, /*countDuplicates=*/true);

        // N-imputation can only ADD seeds (N breaks syncmers in direct extraction but
        // imputed bases preserve them). So every direct seed must appear in the
        // reconstruction with the SAME count; the reconstruction may have extras.
        const bool genomeHasN = genome.find('N') != std::string::npos;
        for (const auto& [h, count] : direct) {
            auto it = reconstructed.find(h);
            BOOST_REQUIRE_MESSAGE(it != reconstructed.end(), "missing seed for node " << nodeId);
            BOOST_TEST(it->second == count);
        }
        if (!genomeHasN) {
            BOOST_TEST(reconstructed.size() == direct.size());
        } else {
            BOOST_TEST(reconstructed.size() >= direct.size());
        }
    }
}

// Exercises the production parallel builder (buildIndexParallel, the -t path) and
// checks it yields the SAME reconstructed genome seeds as the sequential builder.
BOOST_AUTO_TEST_CASE(parallel_build_matches_single_thread) {
    ts::RSVPanmanFixture fix;
    ts::TestIndex seq(fix.tree(), K, S, 0, L, false, false, 0, /*numThreads=*/1);
    ts::TestIndex par(fix.tree(), K, S, 0, L, false, false, 0, /*numThreads=*/4);

    auto& ds = seq.data();
    auto& dp = par.data();
    BOOST_TEST(ds.numChanges() == dp.numChanges());

    const auto& ts_tree = *ds.liteTree;
    const auto& tp_tree = *dp.liteTree;
    for (const auto& nodeId : sampleNodeIds(fix, 6, 5)) {
        int32_t si = ts::findNodeIndex(ts_tree, nodeId);
        int32_t pi = ts::findNodeIndex(tp_tree, nodeId);
        BOOST_REQUIRE(si >= 0 && pi >= 0);
        auto a = ts::reconstructGenomeSeeds(ts::pathToRoot(ts_tree, si), ds);
        auto b = ts::reconstructGenomeSeeds(ts::pathToRoot(tp_tree, pi), dp);
        BOOST_TEST(a.size() == b.size());
        bool equal = true;
        for (const auto& [h, c] : a) {
            auto it = b.find(h);
            if (it == b.end() || it->second != c) {
                equal = false;
                break;
            }
        }
        BOOST_TEST(equal);
    }
}

BOOST_AUTO_TEST_SUITE_END()

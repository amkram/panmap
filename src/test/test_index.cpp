#include <boost/test/unit_test.hpp>

#include "helpers/seed_helpers.hpp"
#include "helpers/test_index.hpp"
#include "helpers/traversal.hpp"
#include "helpers/tree_helpers.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
// l=1 (kminmers == syncmers), flankMaskBp=0 (no masking): matches extractSeeds.
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
    // Non-default params (vs sharedRSVIndex's 15/8/0/1) prove the values aren't hardcoded.
    ts::TestIndex idx(ts::sharedRSVFixture().tree(), /*k=*/19, /*s=*/6, /*t=*/0, /*l=*/3);
    const auto& d = idx.data();
    BOOST_TEST(d.k == 19);
    BOOST_TEST(d.s == 6);
    BOOST_TEST(d.t == 0);
    BOOST_TEST(d.l == 3);
    BOOST_TEST(d.open == false);
    BOOST_TEST(d.hpc == false);
}

BOOST_AUTO_TEST_CASE(offset_and_parent_child_consistency) {
    const auto& fix = ts::sharedRSVFixture();
    const auto& d = ts::sharedRSVIndex();

    const size_t numOffsets = d.numNodesPlusOne();
    BOOST_REQUIRE(numOffsets >= 2);
    BOOST_TEST(d.nodeChangeOffset(0) == 0u);
    BOOST_TEST(d.nodeChangeOffset(numOffsets - 1) == d.numChanges());

    for (size_t i = 1; i < numOffsets; ++i) {
        BOOST_REQUIRE(d.nodeChangeOffset(i) >= d.nodeChangeOffset(i - 1));
    }

    // Each change's parentCount must equal the running (child) count from the previous change for that seed.
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
    const auto& fix = ts::sharedRSVFixture();
    const auto& d = ts::sharedRSVIndex();
    const auto& tree = *d.liteTree;

    std::vector<std::string> targets = {tree.root->identifier};
    for (const auto& id : sampleNodeIds(fix, 5, 42)) targets.push_back(id);

    for (const auto& nodeId : targets) {
        int32_t ni = ts::findNodeIndex(tree, nodeId);
        BOOST_REQUIRE(ni >= 0);
        auto reconstructed = ts::reconstructGenomeSeeds(ts::pathToRoot(tree, ni), d);

        std::string genome = fix.genomeOf(nodeId);
        auto direct = ts::extractSeeds(genome, K, S, /*countDuplicates=*/true);

        // N-imputation only adds seeds (N breaks syncmers in direct extraction, imputed bases
        // preserve them), so every direct seed appears in the reconstruction; extras are allowed.
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

BOOST_AUTO_TEST_CASE(parallel_build_matches_single_thread) {
    const auto& fix = ts::sharedRSVFixture();
    const auto& ds = ts::sharedRSVIndex();   // the shared index IS the sequential (1-thread) build
    ts::TestIndex par(fix.tree(), K, S, 0, L, false, false, 0, /*numThreads=*/4);

    const auto& dp = par.data();
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

// rsv_4K has no block inversions, so the strand-flip seed recompute (computeNewSyncmerRangesJump,
// index_single_mode.cpp ~line 307) is otherwise untested. Inject synthetic pure inversions at
// leaves, assert incremental seeds == from-scratch. Also validated against 89 real inversions in
// tb_400.panman.
BOOST_AUTO_TEST_CASE(block_inversion_reconstruction_equals_direct) {
    ts::RSVPanmanFixture fix;
    panmanUtils::Tree* T = fix.tree();

    std::vector<std::string> leaves;
    for (const auto& id : fix.nodeIds()) {
        if (T->allNodes.at(id)->children.empty()) leaves.push_back(id);
    }
    std::mt19937_64 g(7);
    std::shuffle(leaves.begin(), leaves.end(), g);

    // Inject an inversion of the largest forward block (>= k, so the flip changes seeds) into the
    // first N leaves; rest are controls. Record pre-injection genomes.
    std::vector<std::string> injected, controls;
    std::unordered_map<std::string, std::string> genomeBefore;
    const size_t kInject = 25;
    for (const auto& id : leaves) {
        if (injected.size() >= kInject) {
            controls.push_back(id);
            continue;
        }
        std::vector<std::vector<std::pair<char, std::vector<char>>>> seq;
        std::vector<char> blockExists, blockStrand;
        std::unordered_map<int, int> blockLengths;
        panmapUtils::getSequenceFromReference(T, seq, blockExists, blockStrand, blockLengths, id);
        // getSequenceFromReference fills seq[b] with an existing block's nucleotides (blockLengths
        // tracks only absent blocks). Pick the largest forward block (>= k bases).
        int best = -1, bestLen = 0;
        for (size_t b = 0; b < blockExists.size() && b < seq.size(); ++b) {
            if (!blockExists[b] || !blockStrand[b]) continue;
            int len = static_cast<int>(seq[b].size());
            if (len >= K && len > bestLen) {
                bestLen = len;
                best = static_cast<int>(b);
            }
        }
        if (best < 0) {
            controls.push_back(id);
            continue;
        }
        genomeBefore[id] = T->getStringFromReference(id, false);
        panmanUtils::BlockMut bm;
        bm.primaryBlockId = best;
        bm.secondaryBlockId = -1;
        bm.blockMutInfo = false;  // false + inversion=true == pure block inversion of an existing block
        bm.inversion = true;
        T->allNodes.at(id)->blockMutation.push_back(bm);
        injected.push_back(id);
    }
    BOOST_TEST_MESSAGE("injected pure inversions into " << injected.size() << " leaves; controls=" << controls.size());
    BOOST_REQUIRE(!injected.empty());

    // Build the index from the mutated tree.
    ts::TestIndex idx(T, K, S, 0, L, /*open=*/false, /*hpc=*/false, /*flankMaskBp=*/0);
    auto& d = idx.data();
    const auto& tree = *d.liteTree;

    int meaningful = 0;
    for (const auto& id : injected) {
        if (T->getStringFromReference(id, false) != genomeBefore[id]) meaningful++;
    }
    BOOST_TEST_MESSAGE("meaningful (genome-changing) injections: " << meaningful << "/" << injected.size());
    BOOST_TEST(meaningful == static_cast<int>(injected.size()));

    // Multiset diff: incremental (reconstructed) vs from-scratch (direct).
    auto compareNode = [&](const std::string& nodeId, long& missing, long& extra, bool& hasN) -> bool {
        int32_t ni = ts::findNodeIndex(tree, nodeId);
        if (ni < 0) return false;
        auto reconstructed = ts::reconstructGenomeSeeds(ts::pathToRoot(tree, ni), d);
        std::string genome = T->getStringFromReference(nodeId, false);
        hasN = genome.find('N') != std::string::npos || genome.find('n') != std::string::npos;
        auto direct = ts::extractSeeds(genome, K, S, /*countDuplicates=*/true);
        missing = 0;  // seeds in the true genome that the index failed to record
        extra = 0;    // stale seeds the index kept that the true genome lacks
        for (const auto& [h, c] : direct) {
            auto it = reconstructed.find(h);
            long r = (it == reconstructed.end()) ? 0 : it->second;
            if (r < c) missing += (c - r);
        }
        for (const auto& [h, c] : reconstructed) {
            auto it = direct.find(h);
            long dd = (it == direct.end()) ? 0 : it->second;
            if (c > dd) extra += (c - dd);
        }
        return true;
    };

    auto runSet = [&](const char* label, const std::vector<std::string>& ids, size_t cap, long& mism,
                      long& missTot, long& extraTot) {
        std::vector<std::string> sample = ids;
        std::mt19937_64 g(99);
        std::shuffle(sample.begin(), sample.end(), g);
        if (sample.size() > cap) sample.resize(cap);
        mism = 0; missTot = 0; extraTot = 0;
        int reported = 0;
        for (const auto& id : sample) {
            long missing = 0, extra = 0; bool hasN = false;
            if (!compareNode(id, missing, extra, hasN)) continue;
            if (missing || extra) {
                mism++; missTot += missing; extraTot += extra;
                if (reported++ < 6)
                    BOOST_TEST_MESSAGE("  [" << label << "] MISMATCH node=" << id
                                       << " missing=" << missing << " extra=" << extra << " hasN=" << hasN);
            }
        }
        BOOST_TEST_MESSAGE(label << ": tested=" << sample.size() << " mismatched=" << mism
                           << " missingTotal=" << missTot << " extraTotal=" << extraTot);
    };

    long ctrlMism = 0, ctrlMiss = 0, ctrlExtra = 0;
    runSet("CONTROL", controls, 40, ctrlMism, ctrlMiss, ctrlExtra);
    long invMism = 0, invMiss = 0, invExtra = 0;
    runSet("INJECTED-INVERSION", injected, kInject, invMism, invMiss, invExtra);

    // Inversion nodes must satisfy the same contract as controls (baseline 0); any excess is a seed-recompute defect.
    BOOST_TEST(ctrlMiss == 0);
    BOOST_TEST(ctrlExtra == 0);
    BOOST_TEST(invMiss == 0);
    BOOST_TEST(invExtra == 0);
}

BOOST_AUTO_TEST_SUITE_END()

// Placement metric tests. The smoking-gun test asserts the live NodeMetrics getters
// EQUAL the independent GroundTruthMetrics oracle (within float epsilon) — replacing
// the old suite's misleading score>0 assertions. Genome-equality is owned by
// test_index; this file assumes it and focuses on the read-interaction metrics.
#include <boost/test/unit_test.hpp>

#include "placement.hpp"
#include "helpers/metrics_oracle.hpp"
#include "helpers/seed_helpers.hpp"
#include "helpers/test_index.hpp"
#include "helpers/traversal.hpp"
#include "helpers/tree_helpers.hpp"

#include <cmath>
#include <tuple>
#include <vector>

namespace {
constexpr int K = 15, S = 8, L = 1;
const std::string kTruthNode = "MZ515733.1";

using Change = std::tuple<uint64_t, int64_t, int64_t>;  // (hash, parentCount, childCount)

// Build a PlacementGlobalState from a read seed set, computing the same denominators
// the live getters and the oracle both divide by. Both sides use THIS state, so the
// comparison isolates the numerator logic (computeChildMetrics vs oracle::compute).
placement::PlacementGlobalState makeState(const indexUtils::SeedCountMap& readSeeds,
                                          const indexUtils::SeedCountMap& rootGenome) {
    placement::PlacementGlobalState state;
    double logMagSq = 0.0, logSum = 0.0, wcDenom = 0.0;
    for (const auto& [hash, count] : readSeeds) {
        const double lc = std::log1p(static_cast<double>(count));
        state.seedFreqInReads[hash] = count;
        state.logReadCounts[hash] = lc;
        logMagSq += lc * lc;
        logSum += lc;
        auto it = rootGenome.find(hash);
        if (it != rootGenome.end() && it->second > 0) {
            const double inv = 1.0 / static_cast<double>(it->second);
            state.seedInverseGenomeCounts[hash] = inv;
            wcDenom += inv;
        }
    }
    state.readUniqueSeedCount = readSeeds.size();
    state.logReadMagnitude = std::sqrt(logMagSq);
    state.logContainmentDenominator = logSum;
    state.weightedContainmentDenominator = wcDenom > 0.0 ? wcDenom : 1.0;
    return state;
}

placement::NodeMetrics accumulateAlongPath(const ts::IndexData& d,
                                           const std::vector<panmapUtils::LiteNode*>& path,
                                           const placement::PlacementGlobalState& state) {
    placement::NodeMetrics metrics;
    for (auto* node : path) {
        std::vector<Change> changes;
        uint64_t start = d.nodeChangeOffset(node->nodeIndex);
        uint64_t end = d.nodeChangeOffset(node->nodeIndex + 1);
        changes.reserve(end - start);
        for (uint64_t j = start; j < end; ++j) {
            changes.emplace_back(d.hashAt(j), d.parentCountAt(j), d.childCountAt(j));
        }
        placement::NodeMetrics::computeChildMetrics(metrics, changes, state);
    }
    return metrics;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(placement_tests)

// Hand-built seed-change spans exercised directly through the real computeChildMetrics,
// checked against the oracle. Pins the corners where delta bugs hide.
BOOST_AUTO_TEST_CASE(child_metrics_unit_cases) {
    const uint64_t H = 0xAAAA;  // a read seed
    const uint64_t X = 0xBBBB;  // NOT a read seed

    placement::PlacementGlobalState state;
    state.logReadCounts[H] = std::log1p(3.0);  // read count 3
    state.seedFreqInReads[H] = 3;
    state.readUniqueSeedCount = 1;
    state.logReadMagnitude = std::log1p(3.0);  // single seed
    state.logContainmentDenominator = std::log1p(3.0);
    state.weightedContainmentDenominator = 1.0;

    // (a) H appears in genome (0 -> 2): all read-interaction numerators update.
    {
        std::vector<Change> ch = {{H, 0, 2}};
        placement::NodeMetrics m;
        placement::NodeMetrics::computeChildMetrics(m, ch, state);
        auto oracle = indexUtils::GroundTruthMetrics::compute({{H, 2}}, state);
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), oracle.logRawScore(state), 1e-3);
        BOOST_CHECK_CLOSE(m.getLogCosineScore(state.logReadMagnitude), oracle.logCosineScore(state), 1e-3);
        BOOST_CHECK_CLOSE(m.getWeightedContainmentScore(state.weightedContainmentDenominator),
                          oracle.weightedContainmentScore(state),
                          1e-3);
        BOOST_CHECK_CLOSE(
            m.getLogContainmentScore(state.logContainmentDenominator), oracle.logContainmentScore(state), 1e-3);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 1u);

        // Hand-derived ground truth for r=3, g=2 (independently verified, not via the
        // oracle): logRaw = log(4)/2 / log(4) = 0.5; logCosine = log(4)log(3) /
        // (log(4)·log(3)) = 1; containment = 1/1 = 1; weightedContainment = (1/2)/1 =
        // 0.5; logContainment = log(4)/log(4) = 1.
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), 0.5, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogCosineScore(state.logReadMagnitude), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getContainmentScore(state.readUniqueSeedCount), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getWeightedContainmentScore(state.weightedContainmentDenominator), 0.5, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogContainmentScore(state.logContainmentDenominator), 1.0, 1e-3);
    }

    // (b) Only a non-read seed X changes: read-interaction numerators stay 0,
    //     but genome magnitude still updates.
    {
        std::vector<Change> ch = {{X, 0, 5}};
        placement::NodeMetrics m;
        placement::NodeMetrics::computeChildMetrics(m, ch, state);
        BOOST_CHECK_SMALL(m.getLogRawScore(state.logReadMagnitude), 1e-12);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 0u);
        BOOST_CHECK_GT(m.genomeMagnitudeSquared, 0.0);  // genome-only metric updated
    }

    // (c) freqDelta == 0 but seed present (2 -> 2): no-op for every metric.
    {
        std::vector<Change> ch = {{H, 2, 2}};
        placement::NodeMetrics m;
        placement::NodeMetrics::computeChildMetrics(m, ch, state);
        BOOST_CHECK_SMALL(m.logRawNumerator, 1e-12);
        BOOST_CHECK_SMALL(m.genomeMagnitudeSquared, 1e-12);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 0u);
    }
}

// THE smoking-gun test: for nodes along the path to a known leaf, the five live
// getters must equal the oracle computed from the reconstructed genome + same state.
BOOST_AUTO_TEST_CASE(all_metrics_equal_ground_truth_at_nodes) {
    ts::RSVPanmanFixture fix;
    ts::TestIndex idx(fix.tree(), K, S, 0, L);
    auto& d = idx.data();
    const auto& tree = *d.liteTree;

    // Reads from the truth genome; seeds extracted the same way as the genome.
    std::string truthGenome = fix.genomeOf(kTruthNode);
    auto reads = ts::generateReads(truthGenome, 150, 200, 777);
    BOOST_REQUIRE(!reads.empty());
    indexUtils::SeedCountMap readSeeds;
    for (const auto& r : reads) {
        for (const auto& [h, c] : ts::extractSeeds(r, K, S, true)) readSeeds[h] += c;
    }
    BOOST_REQUIRE(!readSeeds.empty());

    int32_t targetIdx = ts::findNodeIndex(tree, kTruthNode);
    BOOST_REQUIRE(targetIdx >= 0);
    auto fullPath = ts::pathToRoot(tree, targetIdx);
    BOOST_REQUIRE(fullPath.size() >= 2);

    // Root genome supplies the weighted-containment denominator.
    auto rootGenome = ts::reconstructGenomeSeeds({tree.root}, d);
    auto state = makeState(readSeeds, rootGenome);

    // Check several nodes along the path (root, midpoint, leaf).
    std::vector<size_t> checkpoints = {1, fullPath.size() / 2, fullPath.size()};
    int comparisons = 0;
    for (size_t prefixLen : checkpoints) {
        std::vector<panmapUtils::LiteNode*> prefix(fullPath.begin(), fullPath.begin() + prefixLen);
        auto live = accumulateAlongPath(d, prefix, state);
        auto genome = ts::reconstructGenomeSeeds(prefix, d);
        auto oracle = indexUtils::GroundTruthMetrics::compute(genome, state);

        // Genome-only sanity: oracle and live agree on the genome magnitude too.
        BOOST_CHECK_CLOSE(live.genomeMagnitudeSquared, oracle.genomeMagnitudeSquared, 1e-3);
        BOOST_CHECK_EQUAL(live.presenceIntersectionCount, oracle.presenceIntersectionCount);

        auto closeOrZero = [](double a, double b) {
            if (std::abs(b) < 1e-12)
                BOOST_CHECK_SMALL(a, 1e-9);
            else
                BOOST_CHECK_CLOSE(a, b, 1e-3);
        };
        closeOrZero(live.getLogRawScore(state.logReadMagnitude), oracle.logRawScore(state));
        closeOrZero(live.getLogCosineScore(state.logReadMagnitude), oracle.logCosineScore(state));
        closeOrZero(live.getContainmentScore(state.readUniqueSeedCount), oracle.containmentScore(state));
        closeOrZero(live.getWeightedContainmentScore(state.weightedContainmentDenominator),
                    oracle.weightedContainmentScore(state));
        closeOrZero(live.getLogContainmentScore(state.logContainmentDenominator), oracle.logContainmentScore(state));
        comparisons++;
    }
    BOOST_TEST(comparisons == 3);
}

BOOST_AUTO_TEST_SUITE_END()

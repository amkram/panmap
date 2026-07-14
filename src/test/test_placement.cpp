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

// Denominators match what the live getters and oracle both divide by, isolating numerators.
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

void computeChildMetricsFor(placement::NodeMetrics& m,
                            const std::vector<Change>& changes,
                            placement::PlacementGlobalState& state) {
    std::vector<uint64_t> hashes;
    std::vector<int16_t> parents, children;
    hashes.reserve(changes.size());
    parents.reserve(changes.size());
    children.reserve(changes.size());
    for (const auto& [h, p, c] : changes) {
        hashes.push_back(h);
        parents.push_back(static_cast<int16_t>(p));
        children.push_back(static_cast<int16_t>(c));
    }
    panmapUtils::LiteTree tree;
    tree.segSeedHash = {hashes.data()};
    tree.segSeedParent = {parents.data()};
    tree.segSeedChild = {children.data()};
    const panmapUtils::LiteTree* prev = state.liteTree;
    state.liteTree = &tree;
    placement::NodeMetrics::computeChildMetrics(m, 0, static_cast<uint32_t>(changes.size()), state);
    state.liteTree = prev;
}

placement::NodeMetrics accumulateAlongPath(const ts::IndexData& d,
                                           const std::vector<panmapUtils::LiteNode*>& path,
                                           placement::PlacementGlobalState& state) {
    placement::NodeMetrics metrics;
    for (auto* node : path) {
        std::vector<Change> changes;
        uint64_t start = d.nodeChangeOffset(node->nodeIndex);
        uint64_t end = d.nodeChangeOffset(node->nodeIndex + 1);
        changes.reserve(end - start);
        for (uint64_t j = start; j < end; ++j) {
            changes.emplace_back(d.hashAt(j), d.parentCountAt(j), d.childCountAt(j));
        }
        computeChildMetricsFor(metrics, changes, state);
    }
    return metrics;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(placement_tests)

BOOST_AUTO_TEST_CASE(child_metrics_unit_cases) {
    const uint64_t H = 0xAAAA;  // a read seed
    const uint64_t X = 0xBBBB;  // not a read seed

    placement::PlacementGlobalState state;
    state.logReadCounts[H] = std::log1p(3.0);
    state.seedFreqInReads[H] = 3;
    state.readUniqueSeedCount = 1;
    state.logReadMagnitude = std::log1p(3.0);
    state.logContainmentDenominator = std::log1p(3.0);
    state.weightedContainmentDenominator = 1.0;

    // (a) H in genome (0 -> 2).
    {
        std::vector<Change> ch = {{H, 0, 2}};
        placement::NodeMetrics m;
        computeChildMetricsFor(m, ch, state);
        auto oracle = indexUtils::GroundTruthMetrics::compute({{H, 2}}, state);
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), oracle.logRawScore(state), 1e-3);
        BOOST_CHECK_CLOSE(m.getLogCosineScore(state.logReadMagnitude), oracle.logCosineScore(state), 1e-3);
        BOOST_CHECK_CLOSE(m.getWeightedContainmentScore(state.weightedContainmentDenominator),
                          oracle.weightedContainmentScore(state),
                          1e-3);
        BOOST_CHECK_CLOSE(
            m.getLogContainmentScore(state.logContainmentDenominator), oracle.logContainmentScore(state), 1e-3);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 1u);

        // Hand-derived for r=3, g=2: logRaw=0.5, logCosine=1, containment=1, weightedContainment=0.5, logContainment=1.
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), 0.5, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogCosineScore(state.logReadMagnitude), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getContainmentScore(state.readUniqueSeedCount), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getWeightedContainmentScore(state.weightedContainmentDenominator), 0.5, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogContainmentScore(state.logContainmentDenominator), 1.0, 1e-3);
    }

    // (b) Only non-read seed X changes: read numerators stay 0, genome magnitude still updates.
    {
        std::vector<Change> ch = {{X, 0, 5}};
        placement::NodeMetrics m;
        computeChildMetricsFor(m, ch, state);
        BOOST_CHECK_SMALL(m.getLogRawScore(state.logReadMagnitude), 1e-12);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 0u);
        BOOST_CHECK_GT(m.genomeMagnitudeSquared, 0.0);
    }

    // (c) freqDelta == 0, seed present (2 -> 2): no-op for every metric.
    {
        std::vector<Change> ch = {{H, 2, 2}};
        placement::NodeMetrics m;
        computeChildMetricsFor(m, ch, state);
        BOOST_CHECK_SMALL(m.logRawNumerator, 1e-12);
        BOOST_CHECK_SMALL(m.genomeMagnitudeSquared, 1e-12);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 0u);
    }
}

BOOST_AUTO_TEST_CASE(child_metrics_multiseed_hand_derived) {
    const uint64_t H1 = 0x1111, H2 = 0x2222;
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

    // Reads: both seeds at count 3 (log(1+3)=log4); root genome: both at count 2.
    indexUtils::SeedCountMap readSeeds{{H1, 3}, {H2, 3}};
    indexUtils::SeedCountMap rootGenome{{H1, 2}, {H2, 2}};
    auto state = makeState(readSeeds, rootGenome);

    // (a) Both seeds in genome at count 2 (symmetric r=3,g=2): logRaw=1/sqrt2, other four scores=1.
    {
        std::vector<Change> ch = {{H1, 0, 2}, {H2, 0, 2}};
        placement::NodeMetrics m;
        computeChildMetricsFor(m, ch, state);
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), inv_sqrt2, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogCosineScore(state.logReadMagnitude), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getContainmentScore(state.readUniqueSeedCount), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getWeightedContainmentScore(state.weightedContainmentDenominator), 1.0, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogContainmentScore(state.logContainmentDenominator), 1.0, 1e-3);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 2u);
        auto oracle = indexUtils::GroundTruthMetrics::compute({{H1, 2}, {H2, 2}}, state);
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), oracle.logRawScore(state), 1e-3);
    }

    // (b) Genome has only H1; H2 read seed missing -> containment=1/2, read scores halve (logRaw=1/(2·sqrt2), logCosine=1/sqrt2).
    {
        std::vector<Change> ch = {{H1, 0, 2}};
        placement::NodeMetrics m;
        computeChildMetricsFor(m, ch, state);
        BOOST_CHECK_CLOSE(m.getLogRawScore(state.logReadMagnitude), 0.5 * inv_sqrt2, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogCosineScore(state.logReadMagnitude), inv_sqrt2, 1e-3);
        BOOST_CHECK_CLOSE(m.getContainmentScore(state.readUniqueSeedCount), 0.5, 1e-3);
        BOOST_CHECK_CLOSE(m.getWeightedContainmentScore(state.weightedContainmentDenominator), 0.5, 1e-3);
        BOOST_CHECK_CLOSE(m.getLogContainmentScore(state.logContainmentDenominator), 0.5, 1e-3);
        BOOST_CHECK_EQUAL(m.presenceIntersectionCount, 1u);
    }
}

BOOST_AUTO_TEST_CASE(all_metrics_equal_ground_truth_at_nodes) {
    const auto& fix = ts::sharedRSVFixture();
    const auto& d = ts::sharedRSVIndex();
    const auto& tree = *d.liteTree;

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

    auto rootGenome = ts::reconstructGenomeSeeds({tree.root}, d);
    auto state = makeState(readSeeds, rootGenome);

    std::vector<size_t> checkpoints = {1, fullPath.size() / 2, fullPath.size()};
    int comparisons = 0, nonTrivial = 0;
    for (size_t prefixLen : checkpoints) {
        std::vector<panmapUtils::LiteNode*> prefix(fullPath.begin(), fullPath.begin() + prefixLen);
        auto live = accumulateAlongPath(d, prefix, state);
        auto genome = ts::reconstructGenomeSeeds(prefix, d);
        auto oracle = indexUtils::GroundTruthMetrics::compute(genome, state);

        BOOST_CHECK_CLOSE(live.genomeMagnitudeSquared, oracle.genomeMagnitudeSquared, 1e-3);
        BOOST_CHECK_EQUAL(live.presenceIntersectionCount, oracle.presenceIntersectionCount);

        // A ~0 oracle value only pins live ~0; count nonzero comparisons so an all-zero
        // degenerate run (no read/genome overlap) can't pass with trivial checks.
        auto closeOrZero = [&](double a, double b) {
            if (std::abs(b) < 1e-12) {
                BOOST_CHECK_SMALL(a, 1e-9);
            } else {
                BOOST_CHECK_CLOSE(a, b, 1e-3);
                nonTrivial++;
            }
        };
        closeOrZero(live.getLogRawScore(state.logReadMagnitude), oracle.logRawScore(state));
        closeOrZero(live.getLogCosineScore(state.logReadMagnitude), oracle.logCosineScore(state));
        closeOrZero(live.getContainmentScore(state.readUniqueSeedCount), oracle.containmentScore(state));
        closeOrZero(live.getWeightedContainmentScore(state.weightedContainmentDenominator),
                    oracle.weightedContainmentScore(state));
        closeOrZero(live.getLogContainmentScore(state.logContainmentDenominator), oracle.logContainmentScore(state));

        // At the leaf, reads came from this genome, so containment must be substantial.
        if (prefixLen == fullPath.size()) {
            BOOST_TEST(oracle.containmentScore(state) > 0.5);
            BOOST_TEST(live.getContainmentScore(state.readUniqueSeedCount) > 0.5);
        }
        comparisons++;
    }
    BOOST_TEST(comparisons == 3);
    BOOST_TEST(nonTrivial >= 5);
}

BOOST_AUTO_TEST_CASE(read_seed_state_and_min_support) {
    // resolveMinReadSupport auto-mode (configured = -1).
    {  // High coverage: three seeds seen in >=2 reads, mean count 4 > 3 -> require >=2.
        placement::PlacementGlobalState st;
        st.seedFreqInReads[0xA1] = 5;
        st.seedFreqInReads[0xB2] = 4;
        st.seedFreqInReads[0xC3] = 3;
        BOOST_TEST(placement::resolveMinReadSupport(st.seedFreqInReads, -1) == 2);
    }
    {  // Low coverage: one multi-read seed (count 2), mean 2 <= 3 -> keep all (1).
        placement::PlacementGlobalState st;
        st.seedFreqInReads[0xA1] = 2;
        st.seedFreqInReads[0xB2] = 1;
        st.seedFreqInReads[0xC3] = 1;
        BOOST_TEST(placement::resolveMinReadSupport(st.seedFreqInReads, -1) == 1);
    }
    {  // No multi-read seeds -> estCov 0 -> keep all (1).
        placement::PlacementGlobalState st;
        st.seedFreqInReads[0xA1] = 1;
        BOOST_TEST(placement::resolveMinReadSupport(st.seedFreqInReads, -1) == 1);
    }
    {  // An explicit value is used verbatim (no auto).
        placement::PlacementGlobalState st;
        st.seedFreqInReads[0xA1] = 5;
        BOOST_TEST(placement::resolveMinReadSupport(st.seedFreqInReads, 7) == 7);
    }

    // minSupport=2 drops the singleton; denominators over survivors only. Hand-derived from a:5, b:3, c:1.
    {
        placement::PlacementGlobalState st;
        st.seedFreqInReads[0xA1] = 5;
        st.seedFreqInReads[0xB2] = 3;
        st.seedFreqInReads[0xC3] = 1;

        const size_t filtered = placement::computeReadSeedMagnitudes(st, 2);
        const double la = std::log1p(5.0), lb = std::log1p(3.0);

        BOOST_TEST(filtered == 1u);
        BOOST_TEST(st.readUniqueSeedCount == 2u);
        BOOST_TEST(st.totalReadSeedFrequency == 9);  // 5+3+1, counted before the filter
        BOOST_TEST(st.logReadCounts.size() == 2u);
        BOOST_TEST(st.logReadCounts.count(0xC3) == 0u);
        BOOST_CHECK_CLOSE(st.logReadCounts.at(0xA1), la, 1e-9);
        BOOST_CHECK_CLOSE(st.logReadCounts.at(0xB2), lb, 1e-9);
        BOOST_CHECK_CLOSE(st.logContainmentDenominator, la + lb, 1e-9);
        BOOST_CHECK_CLOSE(st.logReadMagnitude, std::sqrt(la * la + lb * lb), 1e-9);

        // minSupport=1 keeps every seed.
        placement::PlacementGlobalState all;
        all.seedFreqInReads = st.seedFreqInReads;
        BOOST_TEST(placement::computeReadSeedMagnitudes(all, 1) == 0u);
        BOOST_TEST(all.readUniqueSeedCount == 3u);
    }
}

BOOST_AUTO_TEST_SUITE_END()

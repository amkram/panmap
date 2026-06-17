#include "metrics_oracle.hpp"

#include <algorithm>
#include <cmath>

namespace indexUtils {

GroundTruthMetrics GroundTruthMetrics::compute(const SeedCountMap& nodeGenome,
                                               const placement::PlacementGlobalState& state) {
    GroundTruthMetrics m;
    for (const auto& [hash, genomeCount] : nodeGenome) {
        if (genomeCount <= 0) continue;  // only seeds present in this node's genome
        const double logGenome = std::log1p(static_cast<double>(genomeCount));

        // Genome-only metric (placement.cpp:215) — accumulates over ALL present seeds.
        m.genomeMagnitudeSquared += logGenome * logGenome;
        m.genomeUniqueSeedCount++;

        // Read-interaction metrics only for seeds also present in the (filtered) reads.
        auto it = state.logReadCounts.find(hash);
        if (it == state.logReadCounts.end()) continue;
        const double logReadCount = it->second;

        m.presenceIntersectionCount++;
        m.logRawNumerator += logReadCount / static_cast<double>(genomeCount);   // placement.cpp:240-242
        m.logCosineNumerator += logReadCount * logGenome;                       // placement.cpp:246
        m.weightedContainmentNumerator += 1.0 / static_cast<double>(genomeCount);  // placement.cpp:253-255
        m.logContainmentNumerator += logReadCount;                             // placement.cpp:260
    }
    return m;
}

double GroundTruthMetrics::logRawScore(const placement::PlacementGlobalState& state) const {
    if (state.logReadMagnitude <= 0.0) return 0.0;
    return logRawNumerator / state.logReadMagnitude;
}

double GroundTruthMetrics::logCosineScore(const placement::PlacementGlobalState& state) const {
    const double genomeMagnitude = std::sqrt(genomeMagnitudeSquared);
    if (state.logReadMagnitude <= 0.0 || genomeMagnitude <= 0.0) return 0.0;
    const double score = logCosineNumerator / (state.logReadMagnitude * genomeMagnitude);
    return std::clamp(score, 0.0, 1.0);
}

double GroundTruthMetrics::containmentScore(const placement::PlacementGlobalState& state) const {
    return (state.readUniqueSeedCount > 0)
               ? static_cast<double>(presenceIntersectionCount) / state.readUniqueSeedCount
               : 0.0;
}

double GroundTruthMetrics::weightedContainmentScore(const placement::PlacementGlobalState& state) const {
    return (state.weightedContainmentDenominator > 0.0)
               ? weightedContainmentNumerator / state.weightedContainmentDenominator
               : 0.0;
}

double GroundTruthMetrics::logContainmentScore(const placement::PlacementGlobalState& state) const {
    return (state.logContainmentDenominator > 0.0)
               ? logContainmentNumerator / state.logContainmentDenominator
               : 0.0;
}

}  // namespace indexUtils

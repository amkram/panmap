#pragma once

/**
 * @brief Independent ground-truth oracle for placement metrics.
 *
 * Reproduces the five placement scores exactly:
 *   (a) per-seed numerators use the same formulas as
 *       placement::NodeMetrics::computeChildMetrics (placement.cpp:206-260);
 *   (b) denominators come from the live PlacementGlobalState, which already encodes
 *       the minReadSupport read filter and the root-derived inverse genome counts.
 * Denominators are read from state rather than recomputed, so (live getter == oracle)
 * holds within float epsilon.
 *
 * Formulas (verbatim from placement.cpp computeChildMetrics):
 *   genomeMagnitudeSquared        = Σ_{seed∈genome}            log(1+gc)²            (placement.cpp:215)
 *   logRawNumerator               = Σ_{seed∈reads∩genome}      log(1+rc) / gc        (placement.cpp:240-242, RAW gc)
 *   logCosineNumerator            = Σ_{seed∈reads∩genome}      log(1+rc) · log(1+gc) (placement.cpp:246)
 *   weightedContainmentNumerator  = Σ_{seed∈reads∩genome}      1 / gc                (placement.cpp:253-255, node gc)
 *   logContainmentNumerator       = Σ_{seed∈reads∩genome}      log(1+rc)             (placement.cpp:260)
 *   presenceIntersectionCount     = |reads ∩ genome|
 * where rc = read count (from state), gc = this node's genome count for that seed.
 */

#include "panmap_utils.hpp"
#include "placement.hpp"

#include <absl/container/flat_hash_map.h>
#include <cstdint>
#include <string>

namespace indexUtils {

using SeedCountMap = absl::flat_hash_map<uint64_t, int64_t>;

struct GroundTruthMetrics {
    // Numerators only; denominators come from PlacementGlobalState at score time.
    double logRawNumerator = 0.0;
    double logCosineNumerator = 0.0;
    double genomeMagnitudeSquared = 0.0;
    double weightedContainmentNumerator = 0.0;
    double logContainmentNumerator = 0.0;
    size_t presenceIntersectionCount = 0;

    /**
     * @brief Compute numerators from a node's reconstructed genome seed set and the
     *        live placement state, which supplies the filtered read seeds.
     *
     * @param nodeGenome  hash -> genome count for this node (from reconstructGenomeSeeds)
     * @param state       the running PlacementGlobalState (filtered reads + log counts)
     */
    static GroundTruthMetrics compute(const SeedCountMap& nodeGenome, const placement::PlacementGlobalState& state);

    // The five scores, each divided by the live denominator in state.
    // Mirrors placement::NodeMetrics::get*Score.
    double logRawScore(const placement::PlacementGlobalState& state) const;
    double logCosineScore(const placement::PlacementGlobalState& state) const;
    double containmentScore(const placement::PlacementGlobalState& state) const;
    double weightedContainmentScore(const placement::PlacementGlobalState& state) const;
    double logContainmentScore(const placement::PlacementGlobalState& state) const;
};

}  // namespace indexUtils

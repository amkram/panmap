#include "index_utils.hpp"

namespace indexUtils {

// ============================================================================
// Ground Truth Metrics Implementation
// ============================================================================

GroundTruthMetrics GroundTruthMetrics::compute(
    const SeedCountMap& readSeeds,
    const SeedCountMap& genomeSeeds) {
    
    GroundTruthMetrics m;
    
    // Read stats
    m.readUniqueSeedCount = readSeeds.size();
    for (const auto& [hash, count] : readSeeds) {
        m.totalReadSeedFrequency += count;
        m.readMagnitude += count * count;
    }
    m.readMagnitude = std::sqrt(m.readMagnitude);
    
    // Genome stats
    m.genomeUniqueSeedCount = genomeSeeds.size();
    for (const auto& [hash, count] : genomeSeeds) {
        m.genomeTotalSeedFrequency += count;
        m.genomeMagnitudeSquared += count * count;
    }
    
    // Intersection metrics - iterate over reads and check genome
    for (const auto& [hash, readCount] : readSeeds) {
        auto it = genomeSeeds.find(hash);
        if (it != genomeSeeds.end()) {
            int64_t genomeCount = it->second;
            m.presenceIntersectionCount++;
            m.jaccardNumerator++;  // Presence-based: 1 per shared seed
            m.weightedJaccardNumerator += std::min(readCount, genomeCount);
            m.cosineNumerator += static_cast<double>(readCount) * genomeCount;
        }
    }
    
    // Compute final scores
    // Jaccard: |intersection| / |union|
    size_t genomeOnlyCount = m.genomeUniqueSeedCount - m.presenceIntersectionCount;
    size_t totalUnion = m.readUniqueSeedCount + genomeOnlyCount;
    m.jaccardScore = (totalUnion > 0) ? 
        static_cast<double>(m.jaccardNumerator) / totalUnion : 0.0;
    
    // Presence Jaccard (same as Jaccard for presence-based)
    m.presenceJaccardScore = (totalUnion > 0) ?
        static_cast<double>(m.presenceIntersectionCount) / totalUnion : 0.0;
    
    // Weighted Jaccard: Σmin / Σmax
    int64_t denominator = m.totalReadSeedFrequency + m.genomeTotalSeedFrequency - m.weightedJaccardNumerator;
    m.weightedJaccardScore = (denominator > 0) ?
        static_cast<double>(m.weightedJaccardNumerator) / denominator : 0.0;
    
    // Cosine: dot / (||a|| * ||b||)
    double genomeMagnitude = std::sqrt(m.genomeMagnitudeSquared);
    m.cosineScore = (m.readMagnitude > 0 && genomeMagnitude > 0) ?
        m.cosineNumerator / (m.readMagnitude * genomeMagnitude) : 0.0;
    m.cosineScore = std::clamp(m.cosineScore, 0.0, 1.0);
    
    return m;
}

// ============================================================================
// Node Lookup Implementation
// ============================================================================

int32_t findNodeIndex(
    const panmapUtils::LiteTree& liteTree,
    const std::string& nodeId) {
    
    for (size_t i = 0; i < liteTree.allLiteNodes.size(); i++) {
        if (liteTree.resolveNodeId(i) == nodeId) {
            return static_cast<int32_t>(i);
        }
    }
    return -1;
}

} // namespace indexUtils

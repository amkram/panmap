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
// Path Building Implementation
// ============================================================================

std::vector<panmapUtils::LiteNode*> buildPathToNode(
    const panmapUtils::LiteTree& liteTree,
    int32_t targetNodeIdx) {
    
    std::vector<panmapUtils::LiteNode*> path;
    
    if (targetNodeIdx < 0 || 
        static_cast<size_t>(targetNodeIdx) >= liteTree.dfsIndexToNode.size()) {
        return path;
    }
    
    panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
    while (current != nullptr) {
        path.push_back(current);
        current = current->parent;
    }
    std::reverse(path.begin(), path.end());
    
    return path;
}

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

// ============================================================================
// Seed Reconstruction Implementation
// ============================================================================

SeedCountMap reconstructSeedsViaDeltas(
    const std::vector<panmapUtils::LiteNode*>& pathToTarget,
    const ::capnp::List<uint64_t>::Reader& seedChangeHashes,
    const ::capnp::List<int64_t>::Reader& seedChangeChildCounts,
    const ::capnp::List<uint32_t>::Reader& nodeChangeOffsets) {
    
    SeedCountMap reconstructedSeeds;
    
    for (const auto* node : pathToTarget) {
        uint32_t nodeIdx = node->nodeIndex;
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        for (uint64_t j = startOffset; j < endOffset; j++) {
            uint64_t hash = seedChangeHashes[j];
            int64_t childCount = seedChangeChildCounts[j];
            
            if (childCount > 0) {
                reconstructedSeeds[hash] = childCount;
            } else {
                reconstructedSeeds.erase(hash);
            }
        }
    }
    
    return reconstructedSeeds;
}

SeedCountMap reconstructSeedsViaDeltas(
    const std::vector<panmapUtils::LiteNode*>& pathToTarget,
    const std::vector<std::tuple<uint64_t, int64_t, int64_t>>& allSeedChanges,
    const std::vector<uint32_t>& nodeChangeOffsets) {
    
    SeedCountMap reconstructedSeeds;
    
    for (const auto* node : pathToTarget) {
        uint32_t nodeIdx = node->nodeIndex;
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : allSeedChanges.size();
        
        for (uint64_t j = startOffset; j < endOffset; j++) {
            const auto& [hash, parentCount, childCount] = allSeedChanges[j];
            
            if (childCount > 0) {
                reconstructedSeeds[hash] = childCount;
            } else {
                reconstructedSeeds.erase(hash);
            }
        }
    }
    
    return reconstructedSeeds;
}

// ============================================================================
// Direct Seed Extraction Implementation
// ============================================================================

SeedCountMap extractSeedsFromGenome(
    const std::string& genome,
    int k, int s, bool open, int t) {
    
    SeedCountMap seeds;
    auto syncmers = seeding::rollingSyncmers(genome, k, s, open, t, false);
    
    for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
        if (isSyncmer) {
            seeds[hash]++;
        }
    }
    
    return seeds;
}

SeedCountMap extractSeedsFromReads(
    const std::vector<std::string>& reads,
    int k, int s, bool open, int t) {
    
    SeedCountMap seeds;
    
    for (const auto& read : reads) {
        auto syncmers = seeding::rollingSyncmers(read, k, s, open, t, false);
        for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
            if (isSyncmer) {
                seeds[hash]++;
            }
        }
    }
    
    return seeds;
}

// ============================================================================
// Metric Computation Implementation
// ============================================================================

void computeGenomeStats(
    const SeedCountMap& seedCounts,
    size_t& uniqueCount,
    int64_t& totalFrequency,
    double& magnitudeSquared) {
    
    uniqueCount = seedCounts.size();
    totalFrequency = 0;
    magnitudeSquared = 0.0;
    
    for (const auto& [hash, count] : seedCounts) {
        totalFrequency += count;
        magnitudeSquared += count * count;
    }
}

void computeReadStats(
    const SeedCountMap& seedCounts,
    size_t& uniqueCount,
    int64_t& totalFrequency,
    double& magnitude) {
    
    uniqueCount = seedCounts.size();
    totalFrequency = 0;
    double magnitudeSquared = 0.0;
    
    for (const auto& [hash, count] : seedCounts) {
        totalFrequency += count;
        magnitudeSquared += count * count;
    }
    magnitude = std::sqrt(magnitudeSquared);
}

// ============================================================================
// Validation Utilities Implementation
// ============================================================================

bool compareSeedMaps(
    const SeedCountMap& a,
    const SeedCountMap& b,
    int& missing,
    int& extra,
    int& countMismatch) {
    
    missing = 0;
    extra = 0;
    countMismatch = 0;
    
    // Check for seeds in b but not in a (missing from a)
    for (const auto& [hash, count] : b) {
        auto it = a.find(hash);
        if (it == a.end()) {
            missing++;
        } else if (it->second != count) {
            countMismatch++;
        }
    }
    
    // Check for seeds in a but not in b (extra in a)
    for (const auto& [hash, count] : a) {
        if (b.find(hash) == b.end()) {
            extra++;
        }
    }
    
    return (missing == 0 && extra == 0 && countMismatch == 0);
}

int verifyParentChildConsistency(
    const std::vector<panmapUtils::LiteNode*>& pathToTarget,
    const ::capnp::List<uint64_t>::Reader& seedChangeHashes,
    const ::capnp::List<int64_t>::Reader& seedChangeParentCounts,
    const ::capnp::List<int64_t>::Reader& seedChangeChildCounts,
    const ::capnp::List<uint32_t>::Reader& nodeChangeOffsets) {
    
    SeedCountMap currentCounts;
    int inconsistencies = 0;
    
    for (const auto* node : pathToTarget) {
        uint32_t nodeIdx = node->nodeIndex;
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        for (uint64_t j = startOffset; j < endOffset; j++) {
            uint64_t hash = seedChangeHashes[j];
            int64_t parentCount = seedChangeParentCounts[j];
            int64_t childCount = seedChangeChildCounts[j];
            
            // Verify: parentCount should match what we currently have
            int64_t expectedParent = currentCounts.count(hash) ? currentCounts[hash] : 0;
            if (parentCount != expectedParent) {
                inconsistencies++;
            }
            
            // Update current counts
            if (childCount > 0) {
                currentCounts[hash] = childCount;
            } else {
                currentCounts.erase(hash);
            }
        }
    }
    
    return inconsistencies;
}

} // namespace indexUtils

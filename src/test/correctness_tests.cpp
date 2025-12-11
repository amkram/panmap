/**
 * @file correctness_tests.cpp
 * @brief Comprehensive correctness tests for panmap indexing and placement workflow
 * 
 * Tests cover:
 * 1. Genome reconstruction from panman
 * 2. Placement metric computations (Jaccard, weighted Jaccard, cosine)
 * 3. Seed computation (syncmers, canonical hashes)
 * 4. Index round-trip (build → write → read)
 * 5. Delta-based metric updates vs ground truth direct computation
 * 
 * Uses rsv_4K.panman as the test dataset.
 * Uses MZ515733.1.fa as ground truth genome for metric validation.
 */

#define BOOST_TEST_MODULE PanmapCorrectnessTests
#include <boost/test/unit_test.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <tbb/global_control.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <filesystem>
#include <sstream>
#include <queue>

#include "panman.hpp"
#include "panmanUtils.hpp"
#include "../seeding.hpp"
#include "../placement.hpp"
#include "../index_single_mode.hpp"
#include "../panmap_utils.hpp"
#include "../zstd_compression.hpp"
#include "../index_utils.hpp"
#include "capnp/message.h"
#include "capnp/serialize.h"
#include "index_lite.capnp.h"

namespace fs = std::filesystem;

// Use shared types from index_utils
using indexUtils::SeedCountMap;

// ============================================================================
// Helper Functions
// ============================================================================

// Read a FASTA file and return the sequence (concatenated if multiple records)
std::string readFasta(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open FASTA file: " + path);
    }
    
    std::string sequence;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        // Remove any whitespace
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        sequence += line;
    }
    return sequence;
}

// Generate simulated reads from a genome sequence
std::vector<std::string> generateReads(const std::string& genome, int readLength, int numReads, int seed = 42) {
    std::vector<std::string> reads;
    std::mt19937 gen(seed);
    
    if (genome.size() < static_cast<size_t>(readLength)) {
        return reads;
    }
    
    std::uniform_int_distribution<size_t> posDist(0, genome.size() - readLength);
    
    for (int i = 0; i < numReads; i++) {
        size_t pos = posDist(gen);
        std::string read = genome.substr(pos, readLength);
        // Skip reads with too many Ns
        size_t nCount = std::count(read.begin(), read.end(), 'N');
        if (nCount < read.size() / 4) {
            reads.push_back(read);
        }
    }
    return reads;
}

// Compute seeds from a sequence using the same parameters as placement
std::vector<std::tuple<size_t, int64_t>> extractSeedsWithCounts(
    const std::string& seq, int k, int s, int l = 1) {
    
    // Get syncmers
    auto syncmers = seeding::rollingSyncmers(seq, k, s, false, 0, false);
    
    // Count seed frequencies
    absl::flat_hash_map<size_t, int64_t> seedCounts;
    for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
        if (isSyncmer) {
            seedCounts[hash]++;
        }
    }
    
    // Convert to vector
    std::vector<std::tuple<size_t, int64_t>> result;
    for (const auto& [hash, count] : seedCounts) {
        result.emplace_back(hash, count);
    }
    return result;
}

// Use GroundTruthMetrics from index_utils
using GroundTruthMetrics = indexUtils::GroundTruthMetrics;

// ============================================================================
// Test Fixture - loads RSV panman once for all tests
// ============================================================================

struct RSVPanmanFixture {
    std::unique_ptr<panmanUtils::TreeGroup> TG;
    panmanUtils::Tree* T;
    std::string panmanPath;
    std::string truthGenomePath;
    std::string truthGenome;
    std::string truthNodeId;
    
    RSVPanmanFixture() {
        panmanPath = "data/rsv_4K.panman";
        truthGenomePath = "data/MZ515733.1.fa";
        truthNodeId = "MZ515733.1";
        
        // Load the panman
        std::ifstream inputFile(panmanPath, std::ios::binary);
        BOOST_REQUIRE_MESSAGE(inputFile.is_open(), "Failed to open test panman: " + panmanPath);
        
        auto inBuffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
        inBuffer->push(boost::iostreams::lzma_decompressor());
        inBuffer->push(inputFile);
        
        std::istream inputStream(inBuffer.get());
        TG = std::make_unique<panmanUtils::TreeGroup>(inputStream);
        BOOST_REQUIRE_MESSAGE(!TG->trees.empty(), "Panman has no trees");
        
        T = &TG->trees[0];
        BOOST_TEST_MESSAGE("Loaded RSV panman with " << countNodes(T->root) << " nodes");
        
        // Load truth genome
        truthGenome = readFasta(truthGenomePath);
        BOOST_TEST_MESSAGE("Loaded truth genome MZ515733.1 with " << truthGenome.size() << " bases");
    }
    
    size_t countNodes(panmanUtils::Node* node) {
        if (!node) return 0;
        size_t count = 1;
        for (auto* child : node->children) {
            count += countNodes(child);
        }
        return count;
    }
    
    // Helper to get a random node from the tree
    std::string getRandomNodeId(int seed = 42) {
        std::vector<std::string> allIds;
        collectNodeIds(T->root, allIds);
        BOOST_REQUIRE(!allIds.empty());
        
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> dis(0, allIds.size() - 1);
        return allIds[dis(gen)];
    }
    
    void collectNodeIds(panmanUtils::Node* node, std::vector<std::string>& ids) {
        if (!node) return;
        ids.push_back(node->identifier);
        for (auto* child : node->children) {
            collectNodeIds(child, ids);
        }
    }
    
    // Find node index by ID in LiteTree - delegate to index_utils
    int32_t findNodeIndex(const panmapUtils::LiteTree& liteTree, const std::string& nodeId) {
        return indexUtils::findNodeIndex(liteTree, nodeId);
    }
};


// ============================================================================
// TEST SUITE 1: Genome Reconstruction
// ============================================================================

BOOST_FIXTURE_TEST_SUITE(GenomeReconstructionTests, RSVPanmanFixture)

BOOST_AUTO_TEST_CASE(test_truth_genome_matches_panman) {
    // Verify that MZ515733.1.fa matches what we get from the panman
    std::string panmanSequence = T->getStringFromReference(truthNodeId, false);
    
    BOOST_TEST_MESSAGE("Truth genome length: " << truthGenome.size());
    BOOST_TEST_MESSAGE("Panman sequence length: " << panmanSequence.size());
    
    // They should be the same (ignoring case)
    std::string truthUpper = truthGenome;
    std::string panmanUpper = panmanSequence;
    std::transform(truthUpper.begin(), truthUpper.end(), truthUpper.begin(), ::toupper);
    std::transform(panmanUpper.begin(), panmanUpper.end(), panmanUpper.begin(), ::toupper);
    
    BOOST_TEST(truthUpper == panmanUpper, 
        "Truth genome should match panman extraction for node " << truthNodeId);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================================
// TEST SUITE 2: Seed Computation (Syncmers)
// ============================================================================

BOOST_AUTO_TEST_SUITE(SeedComputationTests)

BOOST_AUTO_TEST_CASE(test_syncmer_basic) {
    // Test basic syncmer detection with known sequence
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    int k = 15;
    int s = 8;
    bool open = false;
    int t = 0;
    
    auto syncmers = seeding::rollingSyncmers(seq, k, s, open, t, false);
    
    BOOST_TEST_MESSAGE("Sequence length: " << seq.size());
    BOOST_TEST_MESSAGE("Found " << syncmers.size() << " syncmers");
    
    // Should find some syncmers
    BOOST_TEST(syncmers.size() > 0);
    
    // All syncmers should be flagged as such
    for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
        BOOST_TEST(isSyncmer == true);
    }
}

BOOST_AUTO_TEST_CASE(test_syncmer_canonical_hash) {
    // Canonical hash should be min of forward and reverse
    std::string seq = "ACGTACGTACGTACGT";
    std::string revcomp = seeding::reverseComplement(seq);
    
    auto [fwdHash, revHash] = seeding::hashSeq(seq);
    auto [rcFwdHash, rcRevHash] = seeding::hashSeq(revcomp);
    
    BOOST_TEST_MESSAGE("Sequence: " << seq);
    BOOST_TEST_MESSAGE("RevComp: " << revcomp);
    BOOST_TEST_MESSAGE("Forward hash: " << fwdHash);
    BOOST_TEST_MESSAGE("Reverse hash: " << revHash);
    
    // Canonical should be the same for both orientations
    size_t canonical1 = std::min(fwdHash, revHash);
    size_t canonical2 = std::min(rcFwdHash, rcRevHash);
    
    BOOST_TEST(canonical1 == canonical2);
}

BOOST_AUTO_TEST_CASE(test_reverse_complement) {
    std::string seq = "ACGT";
    std::string rc = seeding::reverseComplement(seq);
    
    BOOST_TEST(rc == "ACGT");  // ACGT is a palindrome
    
    std::string seq2 = "AAAA";
    std::string rc2 = seeding::reverseComplement(seq2);
    BOOST_TEST(rc2 == "TTTT");
    
    std::string seq3 = "GCGC";
    std::string rc3 = seeding::reverseComplement(seq3);
    BOOST_TEST(rc3 == "GCGC");  // Palindrome
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================================
// TEST SUITE 3: Placement Metrics - Ground Truth Validation with Real Data
// ============================================================================

BOOST_FIXTURE_TEST_SUITE(PlacementMetricTests, RSVPanmanFixture)

BOOST_AUTO_TEST_CASE(test_index_genome_seeds_match_direct_extraction) {
    // Verify that seeds reconstructed via delta traversal for MZ515733.1 match
    // what we get by directly extracting syncmers from the genome
    
    const int k = 15;
    const int s = 8;
    const int l = 1;
    std::string indexPath = "data/test_seed_validation.pmi";
    
    // Build index
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    // Read index
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    
    // Initialize LiteTree and find MZ515733.1 node
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    int32_t targetNodeIdx = findNodeIndex(liteTree, truthNodeId);
    BOOST_REQUIRE_MESSAGE(targetNodeIdx >= 0, "Could not find node " << truthNodeId << " in index");
    
    BOOST_TEST_MESSAGE("Found target node " << truthNodeId << " at index " << targetNodeIdx);
    
    // Get index arrays for delta reconstruction
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Build path from root to target node
    std::vector<panmapUtils::LiteNode*> pathToTarget;
    panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
    while (current != nullptr) {
        pathToTarget.push_back(current);
        current = current->parent;
    }
    std::reverse(pathToTarget.begin(), pathToTarget.end());
    
    BOOST_TEST_MESSAGE("Path from root to " << truthNodeId << " has " << pathToTarget.size() << " nodes");
    
    // Reconstruct seeds via delta traversal
    absl::flat_hash_map<uint64_t, int64_t> reconstructedSeeds;
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
    
    // Compute metrics from reconstructed seeds
    size_t indexUniqueSeedCount = reconstructedSeeds.size();
    int64_t indexTotalFrequency = 0;
    double indexMagnitudeSquared = 0.0;
    for (const auto& [hash, count] : reconstructedSeeds) {
        indexTotalFrequency += count;
        indexMagnitudeSquared += count * count;
    }
    
    BOOST_TEST_MESSAGE("Index metrics (via delta reconstruction) for " << truthNodeId << ":");
    BOOST_TEST_MESSAGE("  Unique seed count: " << indexUniqueSeedCount);
    BOOST_TEST_MESSAGE("  Total seed frequency: " << indexTotalFrequency);
    BOOST_TEST_MESSAGE("  Magnitude squared: " << indexMagnitudeSquared);
    
    // Now compute ground truth by extracting syncmers directly from genome
    absl::flat_hash_map<size_t, int64_t> genomeSeeds;
    auto genomeSyncmers = seeding::rollingSyncmers(truthGenome, k, s, false, 0, false);
    for (const auto& [hash, isReverse, isSyncmer, pos] : genomeSyncmers) {
        if (isSyncmer) {
            genomeSeeds[hash]++;
        }
    }
    
    // Compute ground truth metrics
    size_t truthUniqueSeedCount = genomeSeeds.size();
    int64_t truthTotalFrequency = 0;
    double truthMagnitudeSquared = 0.0;
    for (const auto& [hash, count] : genomeSeeds) {
        truthTotalFrequency += count;
        truthMagnitudeSquared += count * count;
    }
    
    BOOST_TEST_MESSAGE("Ground truth metrics for " << truthNodeId << ":");
    BOOST_TEST_MESSAGE("  Unique seed count: " << truthUniqueSeedCount);
    BOOST_TEST_MESSAGE("  Total seed frequency: " << truthTotalFrequency);
    BOOST_TEST_MESSAGE("  Magnitude squared: " << truthMagnitudeSquared);
    
    // Compare - they should match exactly
    BOOST_TEST(indexUniqueSeedCount == truthUniqueSeedCount,
        "Unique seed count mismatch: index=" << indexUniqueSeedCount << " truth=" << truthUniqueSeedCount);
    BOOST_TEST(indexTotalFrequency == truthTotalFrequency,
        "Total frequency mismatch: index=" << indexTotalFrequency << " truth=" << truthTotalFrequency);
    BOOST_TEST(std::abs(indexMagnitudeSquared - truthMagnitudeSquared) < 1e-6,
        "Magnitude squared mismatch: index=" << indexMagnitudeSquared << " truth=" << truthMagnitudeSquared);
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_CASE(test_placement_metrics_at_truth_node) {
    // CORE TEST: Verify delta-based metric accumulation matches direct computation
    //
    // Ground Truth: Directly compute metrics from:
    //   - readSeeds: syncmers extracted from simulated reads
    //   - genomeSeeds: syncmers extracted directly from MZ515733.1 genome
    //
    // What We Test: Delta-based traversal through index produces same metrics
    
    const int k = 15;
    const int s = 8;
    const int l = 1;
    const int readLength = 150;
    const int numReads = 200;
    std::string indexPath = "data/test_placement_validation.pmi";
    
    // ========================================================================
    // STEP 1: Generate reads and extract seeds
    // ========================================================================
    auto reads = generateReads(truthGenome, readLength, numReads, 42);
    BOOST_REQUIRE(!reads.empty());
    BOOST_TEST_MESSAGE("Generated " << reads.size() << " reads from " << truthNodeId);
    
    absl::flat_hash_map<size_t, int64_t> readSeeds;
    for (const auto& read : reads) {
        auto readSyncmers = seeding::rollingSyncmers(read, k, s, false, 0, false);
        for (const auto& [hash, isReverse, isSyncmer, pos] : readSyncmers) {
            if (isSyncmer) {
                readSeeds[hash]++;
            }
        }
    }
    
    size_t readUniqueSeedCount = readSeeds.size();
    int64_t totalReadSeedFrequency = 0;
    double readMagnitudeSquared = 0.0;
    for (const auto& [hash, count] : readSeeds) {
        totalReadSeedFrequency += count;
        readMagnitudeSquared += count * count;
    }
    double readMagnitude = std::sqrt(readMagnitudeSquared);
    
    BOOST_TEST_MESSAGE("Read stats: " << readUniqueSeedCount << " unique seeds, " 
        << totalReadSeedFrequency << " total, magnitude=" << readMagnitude);
    
    // ========================================================================
    // STEP 2: Extract seeds from truth genome directly
    // ========================================================================
    absl::flat_hash_map<size_t, int64_t> genomeSeeds;
    auto genomeSyncmers = seeding::rollingSyncmers(truthGenome, k, s, false, 0, false);
    
    BOOST_TEST_MESSAGE("DEBUG: truthGenome.size() = " << truthGenome.size());
    BOOST_TEST_MESSAGE("DEBUG: genomeSyncmers.size() = " << genomeSyncmers.size());
    
    // Also extract genome from panman and compare
    std::string panmanGenome = T->getStringFromReference(truthNodeId, false);
    BOOST_TEST_MESSAGE("DEBUG: panmanGenome.size() = " << panmanGenome.size());
    BOOST_TEST_MESSAGE("DEBUG: genomes match = " << (truthGenome == panmanGenome ? "YES" : "NO"));
    if (truthGenome != panmanGenome) {
        size_t diffCount = 0;
        for (size_t i = 0; i < std::min(truthGenome.size(), panmanGenome.size()); ++i) {
            if (truthGenome[i] != panmanGenome[i]) diffCount++;
        }
        diffCount += std::abs(static_cast<int>(truthGenome.size()) - static_cast<int>(panmanGenome.size()));
        BOOST_TEST_MESSAGE("DEBUG: diff count = " << diffCount);
    }
    
    size_t syncmerTrueCount = 0;
    for (const auto& [hash, isReverse, isSyncmer, pos] : genomeSyncmers) {
        if (isSyncmer) {
            genomeSeeds[hash]++;
            syncmerTrueCount++;
        }
    }
    
    BOOST_TEST_MESSAGE("DEBUG: syncmerTrueCount = " << syncmerTrueCount);
    BOOST_TEST_MESSAGE("DEBUG: genomeSeeds.size() (unique) = " << genomeSeeds.size());
    
    // ========================================================================
    // STEP 3: Compute GROUND TRUTH metrics directly from seed sets
    // ========================================================================
    auto groundTruth = GroundTruthMetrics::compute(readSeeds, genomeSeeds);
    
    BOOST_TEST_MESSAGE("=== GROUND TRUTH (direct computation) ===");
    BOOST_TEST_MESSAGE("  Genome unique seeds: " << groundTruth.genomeUniqueSeedCount);
    BOOST_TEST_MESSAGE("  Genome total freq: " << groundTruth.genomeTotalSeedFrequency);
    BOOST_TEST_MESSAGE("  Intersection count: " << groundTruth.presenceIntersectionCount);
    BOOST_TEST_MESSAGE("  Jaccard: " << groundTruth.jaccardScore);
    BOOST_TEST_MESSAGE("  Weighted Jaccard: " << groundTruth.weightedJaccardScore);
    BOOST_TEST_MESSAGE("  Cosine: " << groundTruth.cosineScore);
    
    // ========================================================================
    // STEP 4: Build index and traverse using delta accumulation
    // ========================================================================
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    int32_t targetNodeIdx = findNodeIndex(liteTree, truthNodeId);
    BOOST_REQUIRE(targetNodeIdx >= 0);
    
    // Set up placement state
    placement::PlacementGlobalState state;
    for (const auto& [hash, count] : readSeeds) {
        state.seedFreqInReads[hash] = count;
    }
    state.readUniqueSeedCount = readUniqueSeedCount;
    state.totalReadSeedFrequency = totalReadSeedFrequency;
    state.readMagnitude = readMagnitude;
    
    // Get seed change arrays from index
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeParentCounts = indexRoot.getSeedChangeParentCounts();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Find path from root to target
    std::vector<panmapUtils::LiteNode*> pathToTarget;
    panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
    while (current != nullptr) {
        pathToTarget.push_back(current);
        current = current->parent;
    }
    std::reverse(pathToTarget.begin(), pathToTarget.end());
    
    BOOST_TEST_MESSAGE("Path from root to " << truthNodeId << " has " << pathToTarget.size() << " nodes");
    
    // ========================================================================
    // STEP 5: Accumulate metrics via delta traversal
    // ========================================================================
    placement::NodeMetrics metrics;  // Start empty - represents parent of root
    
    // Arrays already loaded above
    BOOST_TEST_MESSAGE("Index arrays: " << seedChangeHashes.size() << " seed changes, " 
        << nodeChangeOffsets.size() << " offset entries");
    
    // Debug: dump first few offsets
    BOOST_TEST_MESSAGE("First 10 offsets: ");
    for (uint32_t i = 0; i < 10 && i < nodeChangeOffsets.size(); i++) {
        BOOST_TEST_MESSAGE("  nodeChangeOffsets[" << i << "] = " << nodeChangeOffsets[i]);
    }
    
    for (size_t i = 0; i < pathToTarget.size(); i++) {
        panmapUtils::LiteNode* node = pathToTarget[i];
        uint32_t nodeIdx = node->nodeIndex;
        
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        size_t numChanges = endOffset - startOffset;
        BOOST_TEST_MESSAGE("Node " << i << " (idx=" << nodeIdx << "): offset [" 
            << startOffset << ", " << endOffset << ") = " << numChanges << " changes");
        
        std::vector<std::tuple<uint64_t, int64_t, int64_t>> changes;
        changes.reserve(numChanges);
        for (uint64_t j = startOffset; j < endOffset; j++) {
            changes.emplace_back(
                seedChangeHashes[j],
                seedChangeParentCounts[j],
                seedChangeChildCounts[j]
            );
        }
        
        placement::NodeMetrics::computeChildMetrics(metrics, changes, state);
        
        // Debug: show running totals after each node
        BOOST_TEST_MESSAGE("  After node " << i << ": unique=" << metrics.genomeUniqueSeedCount 
            << " total=" << metrics.genomeTotalSeedFrequency
            << " intersect=" << metrics.presenceIntersectionCount);
    }
    
    // VERIFICATION: Manually track seed counts to verify unique count
    // Note: seedChangeChildCounts contains ABSOLUTE counts, not deltas
    absl::flat_hash_map<uint64_t, int64_t> manualSeedCounts;
    for (size_t i = 0; i < pathToTarget.size(); i++) {
        panmapUtils::LiteNode* node = pathToTarget[i];
        uint32_t nodeIdx = node->nodeIndex;
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        for (uint64_t j = startOffset; j < endOffset; j++) {
            uint64_t hash = seedChangeHashes[j];
            int64_t parentCount = seedChangeParentCounts[j];
            int64_t childCount = seedChangeChildCounts[j];
            // Update to child's absolute count
            if (childCount > 0) {
                manualSeedCounts[hash] = childCount;
            } else {
                manualSeedCounts.erase(hash);
            }
        }
    }
    
    int64_t manualUniqueCount = manualSeedCounts.size();
    int64_t manualTotalFreq = 0;
    for (const auto& [h, c] : manualSeedCounts) {
        manualTotalFreq += c;
    }
    BOOST_TEST_MESSAGE("Manual verification: unique=" << manualUniqueCount << " total=" << manualTotalFreq);
    
    // Compute final scores from accumulated metrics
    double deltaJaccard = metrics.getJaccardScore(readUniqueSeedCount);
    double deltaWeightedJaccard = metrics.getWeightedJaccardScore(totalReadSeedFrequency);
    double deltaCosine = metrics.getCosineScore(readMagnitude);
    
    BOOST_TEST_MESSAGE("=== DELTA-BASED (index traversal) ===");
    BOOST_TEST_MESSAGE("  Genome unique seeds: " << metrics.genomeUniqueSeedCount);
    BOOST_TEST_MESSAGE("  Genome total freq: " << metrics.genomeTotalSeedFrequency);
    BOOST_TEST_MESSAGE("  Intersection count: " << metrics.presenceIntersectionCount);
    BOOST_TEST_MESSAGE("  Jaccard: " << deltaJaccard);
    BOOST_TEST_MESSAGE("  Weighted Jaccard: " << deltaWeightedJaccard);
    BOOST_TEST_MESSAGE("  Cosine: " << deltaCosine);
    
    // ========================================================================
    // STEP 6: VERIFY delta-based matches ground truth
    // ========================================================================
    BOOST_TEST(metrics.genomeUniqueSeedCount == groundTruth.genomeUniqueSeedCount,
        "Genome unique: delta=" << metrics.genomeUniqueSeedCount 
        << " truth=" << groundTruth.genomeUniqueSeedCount);
    
    BOOST_TEST(metrics.genomeTotalSeedFrequency == groundTruth.genomeTotalSeedFrequency,
        "Genome total freq: delta=" << metrics.genomeTotalSeedFrequency 
        << " truth=" << groundTruth.genomeTotalSeedFrequency);
    
    BOOST_TEST(metrics.presenceIntersectionCount == groundTruth.presenceIntersectionCount,
        "Intersection: delta=" << metrics.presenceIntersectionCount 
        << " truth=" << groundTruth.presenceIntersectionCount);
    
    BOOST_TEST(std::abs(deltaJaccard - groundTruth.jaccardScore) < 1e-9,
        "Jaccard: delta=" << deltaJaccard << " truth=" << groundTruth.jaccardScore);
    
    BOOST_TEST(std::abs(deltaWeightedJaccard - groundTruth.weightedJaccardScore) < 1e-9,
        "Weighted Jaccard: delta=" << deltaWeightedJaccard << " truth=" << groundTruth.weightedJaccardScore);
    
    BOOST_TEST(std::abs(deltaCosine - groundTruth.cosineScore) < 1e-9,
        "Cosine: delta=" << deltaCosine << " truth=" << groundTruth.cosineScore);
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_CASE(test_reads_from_truth_genome_place_to_truth_node) {
    // Verify that reads from MZ515733.1 get highest scores at that node
    // This tests the full placement logic end-to-end
    
    const int k = 15;
    const int s = 8;
    const int l = 1;
    const int readLength = 150;
    const int numReads = 100;
    std::string indexPath = "data/test_placement_accuracy.pmi";
    
    // Generate reads from truth genome
    auto reads = generateReads(truthGenome, readLength, numReads, 999);
    BOOST_REQUIRE(!reads.empty());
    
    // Build index
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    // Extract seeds from reads
    absl::flat_hash_map<size_t, int64_t> readSeeds;
    for (const auto& read : reads) {
        auto readSyncmers = seeding::rollingSyncmers(read, k, s, false, 0, false);
        for (const auto& [hash, isReverse, isSyncmer, pos] : readSyncmers) {
            if (isSyncmer) {
                readSeeds[hash]++;
            }
        }
    }
    
    size_t readUniqueSeedCount = readSeeds.size();
    int64_t totalReadSeedFrequency = 0;
    double readMagnitudeSquared = 0.0;
    for (const auto& [hash, count] : readSeeds) {
        totalReadSeedFrequency += count;
        readMagnitudeSquared += count * count;
    }
    double readMagnitude = std::sqrt(readMagnitudeSquared);
    
    // Load index
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    int32_t targetNodeIdx = findNodeIndex(liteTree, truthNodeId);
    BOOST_REQUIRE(targetNodeIdx >= 0);
    
    // Set up state
    placement::PlacementGlobalState state;
    for (const auto& [hash, count] : readSeeds) {
        state.seedFreqInReads[hash] = count;
    }
    state.readUniqueSeedCount = readUniqueSeedCount;
    state.totalReadSeedFrequency = totalReadSeedFrequency;
    state.readMagnitude = readMagnitude;
    
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeParentCounts = indexRoot.getSeedChangeParentCounts();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Find best scoring node by traversing all nodes
    double bestJaccard = -1.0;
    int32_t bestJaccardNodeIdx = -1;
    double targetNodeJaccard = 0.0;
    
    // BFS traversal to compute metrics at each node
    // Pass parent metrics to children, apply deltas at each node
    std::queue<std::pair<panmapUtils::LiteNode*, placement::NodeMetrics>> queue;
    
    // Initialize with empty metrics - root's seed changes will build from empty
    placement::NodeMetrics emptyMetrics;
    queue.push({liteTree.root, emptyMetrics});
    
    while (!queue.empty()) {
        auto [node, parentMetrics] = queue.front();
        queue.pop();
        
        placement::NodeMetrics nodeMetrics = parentMetrics;
        uint32_t nodeIdx = node->nodeIndex;
        
        // Get seed changes for this node (including root)
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        std::vector<std::tuple<uint64_t, int64_t, int64_t>> changes;
        for (uint64_t j = startOffset; j < endOffset; j++) {
            changes.emplace_back(
                seedChangeHashes[j],
                seedChangeParentCounts[j],
                seedChangeChildCounts[j]
            );
        }
        
        // Apply delta updates - ALL metrics computed incrementally
        placement::NodeMetrics::computeChildMetrics(nodeMetrics, changes, state);
        
        // Compute Jaccard score
        double jaccard = nodeMetrics.getJaccardScore(readUniqueSeedCount);
        
        if (jaccard > bestJaccard) {
            bestJaccard = jaccard;
            bestJaccardNodeIdx = node->nodeIndex;
        }
        
        if (static_cast<int32_t>(node->nodeIndex) == targetNodeIdx) {
            targetNodeJaccard = jaccard;
        }
        
        // Add children to queue
        for (auto* child : node->children) {
            queue.push({child, nodeMetrics});
        }
    }
    
    std::string bestNodeId = liteTree.resolveNodeId(bestJaccardNodeIdx);
    BOOST_TEST_MESSAGE("Best Jaccard node: " << bestNodeId << " (score: " << bestJaccard << ")");
    BOOST_TEST_MESSAGE("Target node " << truthNodeId << " Jaccard: " << targetNodeJaccard);
    
    // The target node should have one of the highest scores
    // (it may not be THE highest if there are very similar genomes)
    BOOST_TEST(targetNodeJaccard > 0.5, 
        "Reads from " << truthNodeId << " should have high Jaccard at that node");
    
    // Target should be within top scores
    BOOST_TEST(targetNodeJaccard >= bestJaccard * 0.8,
        "Target node should be close to best scoring node");
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================================
// TEST SUITE 4: Delta-Based Metric Updates - Real Data Validation
// ============================================================================

BOOST_FIXTURE_TEST_SUITE(DeltaMetricTests, RSVPanmanFixture)

BOOST_AUTO_TEST_CASE(test_incremental_traversal_matches_direct_at_every_node) {
    // For each node along the path to MZ515733.1, verify that:
    // 1. The incrementally computed metrics match direct computation
    // 2. Genome-only metrics from index match direct extraction
    
    const int k = 15;
    const int s = 8;
    const int l = 1;
    const int readLength = 150;
    const int numReads = 100;
    std::string indexPath = "data/test_delta_validation.pmi";
    
    // Generate reads from truth genome
    auto reads = generateReads(truthGenome, readLength, numReads, 777);
    BOOST_REQUIRE(!reads.empty());
    
    // Extract seeds from reads
    absl::flat_hash_map<size_t, int64_t> readSeeds;
    for (const auto& read : reads) {
        auto readSyncmers = seeding::rollingSyncmers(read, k, s, false, 0, false);
        for (const auto& [hash, isReverse, isSyncmer, pos] : readSyncmers) {
            if (isSyncmer) {
                readSeeds[hash]++;
            }
        }
    }
    
    size_t readUniqueSeedCount = readSeeds.size();
    int64_t totalReadSeedFrequency = 0;
    double readMagnitudeSquared = 0.0;
    for (const auto& [hash, count] : readSeeds) {
        totalReadSeedFrequency += count;
        readMagnitudeSquared += count * count;
    }
    double readMagnitude = std::sqrt(readMagnitudeSquared);
    
    // Build index
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    int32_t targetNodeIdx = findNodeIndex(liteTree, truthNodeId);
    BOOST_REQUIRE(targetNodeIdx >= 0);
    
    // Set up state
    placement::PlacementGlobalState state;
    for (const auto& [hash, count] : readSeeds) {
        state.seedFreqInReads[hash] = count;
    }
    state.readUniqueSeedCount = readUniqueSeedCount;
    state.totalReadSeedFrequency = totalReadSeedFrequency;
    state.readMagnitude = readMagnitude;
    
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeParentCounts = indexRoot.getSeedChangeParentCounts();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Find path from root to target
    std::vector<panmapUtils::LiteNode*> pathToTarget;
    panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
    while (current != nullptr) {
        pathToTarget.push_back(current);
        current = current->parent;
    }
    std::reverse(pathToTarget.begin(), pathToTarget.end());
    
    BOOST_TEST_MESSAGE("Validating " << pathToTarget.size() << " nodes along path to " << truthNodeId);
    
    // Track cumulative genome seeds by applying deltas
    absl::flat_hash_map<size_t, int64_t> accumulatedGenomeSeeds;
    
    // For each node along the path, verify genome metrics match
    // Initialize metrics from empty - will accumulate via deltas
    placement::NodeMetrics metrics;
    
    for (size_t i = 0; i < pathToTarget.size(); i++) {
        panmapUtils::LiteNode* node = pathToTarget[i];
        uint32_t nodeIdx = node->nodeIndex;
        std::string nodeId = liteTree.resolveNodeId(nodeIdx);
        
        // Get genome sequence from panman
        std::string genomeSeq = T->getStringFromReference(nodeId, false);
        
        // Extract seeds directly
        absl::flat_hash_map<size_t, int64_t> directGenomeSeeds;
        auto genomeSyncmers = seeding::rollingSyncmers(genomeSeq, k, s, false, 0, false);
        for (const auto& [hash, isReverse, isSyncmer, pos] : genomeSyncmers) {
            if (isSyncmer) {
                directGenomeSeeds[hash]++;
            }
        }
        
        // Compute direct genome stats
        size_t directUnique = directGenomeSeeds.size();
        int64_t directTotal = 0;
        double directMagSq = 0.0;
        for (const auto& [hash, count] : directGenomeSeeds) {
            directTotal += count;
            directMagSq += count * count;
        }
        
        // Apply seed change deltas to update ALL metrics (genome and intersection)
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        std::vector<std::tuple<uint64_t, int64_t, int64_t>> changes;
        for (uint64_t j = startOffset; j < endOffset; j++) {
            changes.emplace_back(
                seedChangeHashes[j],
                seedChangeParentCounts[j],
                seedChangeChildCounts[j]
            );
        }
        
        placement::NodeMetrics::computeChildMetrics(metrics, changes, state);
        
        // Verify incremental genome metrics match direct extraction
        BOOST_TEST(metrics.genomeUniqueSeedCount == directUnique,
            "Node " << nodeId << " incremental unique: " << metrics.genomeUniqueSeedCount 
            << " vs direct: " << directUnique);
        BOOST_TEST(metrics.genomeTotalSeedFrequency == directTotal,
            "Node " << nodeId << " incremental total: " << metrics.genomeTotalSeedFrequency 
            << " vs direct: " << directTotal);
        BOOST_TEST(std::abs(metrics.genomeMagnitudeSquared - directMagSq) < 1e-6,
            "Node " << nodeId << " incremental mag^2: " << metrics.genomeMagnitudeSquared 
            << " vs direct: " << directMagSq);
        
        // Compute direct metrics for comparison
        auto directTruth = GroundTruthMetrics::compute(readSeeds, directGenomeSeeds);
        
        // Compute scores from delta-based metrics
        double deltaJaccard = metrics.getJaccardScore(readUniqueSeedCount);
        double deltaWeighted = metrics.getWeightedJaccardScore(totalReadSeedFrequency);
        double deltaCosine = metrics.getCosineScore(readMagnitude);
        
        // Verify intersection metrics match at final target node
        if (nodeId == truthNodeId) {
            BOOST_TEST_MESSAGE("Final validation at " << truthNodeId << ":");
            BOOST_TEST_MESSAGE("  Direct intersection: " << directTruth.presenceIntersectionCount);
            BOOST_TEST_MESSAGE("  Delta intersection: " << metrics.presenceIntersectionCount);
            BOOST_TEST_MESSAGE("  Direct Jaccard: " << directTruth.jaccardScore);
            BOOST_TEST_MESSAGE("  Delta Jaccard: " << deltaJaccard);
            
            BOOST_TEST(metrics.presenceIntersectionCount == directTruth.presenceIntersectionCount);
            BOOST_TEST(std::abs(deltaJaccard - directTruth.jaccardScore) < 1e-9);
            BOOST_TEST(std::abs(deltaWeighted - directTruth.weightedJaccardScore) < 1e-9);
            BOOST_TEST(std::abs(deltaCosine - directTruth.cosineScore) < 1e-9);
        }
    }
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================================
// TEST SUITE 5: Index Round-Trip
// ============================================================================

BOOST_FIXTURE_TEST_SUITE(IndexRoundTripTests, RSVPanmanFixture)

BOOST_AUTO_TEST_CASE(test_index_build_and_parameters) {
    // Build an index and verify parameters are stored correctly
    int k = 15;
    int s = 8;
    int l = 1;
    
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    
    BOOST_TEST(builder.getK() == k);
    BOOST_TEST(builder.getS() == s);
    BOOST_TEST(builder.getL() == l);
}

BOOST_AUTO_TEST_CASE(test_index_write_and_read) {
    // Build, write, and read back an index
    std::string indexPath = "data/test_rsv_index.pmi";
    
    int k = 15;
    int s = 8;
    int l = 1;
    
    // Build and write
    {
        index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
        builder.buildIndex();
        builder.writeIndex(indexPath);
    }
    
    BOOST_REQUIRE(fs::exists(indexPath));
    
    // Read back
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    
    // Verify parameters match
    BOOST_TEST(indexRoot.getK() == k);
    BOOST_TEST(indexRoot.getS() == s);
    BOOST_TEST(indexRoot.getL() == l);
    BOOST_TEST(indexRoot.getVersion() == 3);
    
    // Verify V3 arrays are present
    BOOST_TEST(indexRoot.hasSeedChangeHashes());
    BOOST_TEST(indexRoot.hasSeedChangeParentCounts());
    BOOST_TEST(indexRoot.hasSeedChangeChildCounts());
    BOOST_TEST(indexRoot.hasNodeChangeOffsets());
    
    // Verify tree structure
    auto liteTree = indexRoot.getLiteTree();
    auto nodes = liteTree.getLiteNodes();
    BOOST_TEST(nodes.size() > 0);
    
    BOOST_TEST_MESSAGE("Index has " << nodes.size() << " nodes");
    BOOST_TEST_MESSAGE("Total seed changes: " << indexRoot.getTotalSeedChanges());
    
    // Clean up
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================================
// TEST SUITE 6: End-to-End Placement (Integration)
// ============================================================================

BOOST_FIXTURE_TEST_SUITE(PlacementIntegrationTests, RSVPanmanFixture)

BOOST_AUTO_TEST_CASE(test_lite_tree_initialization) {
    std::string indexPath = "data/test_rsv_lite.pmi";
    
    // Build and write index
    {
        index_single_mode::IndexBuilder builder(T, 15, 8, 0, 1, false);
        builder.buildIndex();
        builder.writeIndex(indexPath);
    }
    
    // Read and initialize LiteTree
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    BOOST_TEST(liteTree.allLiteNodes.size() > 0);
    BOOST_TEST(liteTree.root != nullptr);
    
    // Verify node structure
    size_t nodesWithChildren = 0;
    size_t nodesWithParent = 0;
    for (const auto& [nodeId, node] : liteTree.allLiteNodes) {
        if (!node->children.empty()) nodesWithChildren++;
        if (node->parent != nullptr) nodesWithParent++;
    }
    
    BOOST_TEST_MESSAGE("Nodes with children: " << nodesWithChildren);
    BOOST_TEST_MESSAGE("Nodes with parent: " << nodesWithParent);
    
    // Root should have children but no parent
    BOOST_TEST(liteTree.root->parent == nullptr);
    BOOST_TEST(!liteTree.root->children.empty());
    
    // All non-root nodes should have a parent
    BOOST_TEST(nodesWithParent == liteTree.allLiteNodes.size() - 1);
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_CASE(test_resolve_node_id) {
    std::string indexPath = "data/test_rsv_resolve.pmi";
    
    {
        index_single_mode::IndexBuilder builder(T, 15, 8, 0, 1, false);
        builder.buildIndex();
        builder.writeIndex(indexPath);
    }
    
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    // Resolve root node ID
    std::string resolvedRootId = liteTree.resolveNodeId(0);
    
    BOOST_TEST(!resolvedRootId.empty());
    BOOST_TEST_MESSAGE("Resolved root ID: " << resolvedRootId);
    
    // Should match the original tree's root
    BOOST_TEST(resolvedRootId == T->root->identifier);
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================================
// Indexing Correctness Tests
// ============================================================================

BOOST_FIXTURE_TEST_SUITE(IndexingCorrectnessTests, RSVPanmanFixture)

/**
 * CORE TEST: Verify index DFS produces correct seed state at each node.
 * 
 * For multiple target nodes, check that:
 * 1. Reconstructing seeds by traversing root→target through index deltas
 * 2. Matches directly computing seeds from that node's genome sequence
 */
BOOST_AUTO_TEST_CASE(test_index_dfs_state_matches_genome) {
    const int k = 15;
    const int s = 8;
    const int l = 1;
    std::string indexPath = "data/test_index_dfs_state.pmi";
    
    // Build index
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    // Load index
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Test multiple nodes: root, truth node, and a few random ones
    std::vector<std::string> testNodeIds = {
        T->root->identifier,
        truthNodeId,
        getRandomNodeId(1),
        getRandomNodeId(2),
        getRandomNodeId(3)
    };
    
    int totalFailures = 0;
    
    for (const auto& nodeId : testNodeIds) {
        int32_t targetNodeIdx = findNodeIndex(liteTree, nodeId);
        if (targetNodeIdx < 0) {
            BOOST_TEST_MESSAGE("Skipping node " << nodeId << " - not found in index");
            continue;
        }
        
        // Build path from root to target
        std::vector<panmapUtils::LiteNode*> pathToTarget;
        panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
        while (current != nullptr) {
            pathToTarget.push_back(current);
            current = current->parent;
        }
        std::reverse(pathToTarget.begin(), pathToTarget.end());
        
        // Reconstruct seed state by applying deltas along path
        absl::flat_hash_map<uint64_t, int64_t> reconstructedSeeds;
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
        
        // Compute ground truth seeds directly from genome
        std::string genome = T->getStringFromReference(nodeId, false);
        absl::flat_hash_map<uint64_t, int64_t> truthSeeds;
        for (const auto& [hash, isReverse, isSyncmer, pos] : seeding::rollingSyncmers(genome, k, s, false, 0, false)) {
            if (isSyncmer) {
                truthSeeds[hash]++;
            }
        }
        
        // Compare
        bool match = (reconstructedSeeds.size() == truthSeeds.size());
        if (match) {
            for (const auto& [hash, count] : truthSeeds) {
                auto it = reconstructedSeeds.find(hash);
                if (it == reconstructedSeeds.end() || it->second != count) {
                    match = false;
                    break;
                }
            }
        }
        
        if (!match) {
            totalFailures++;
            BOOST_TEST_MESSAGE("FAIL: Node " << nodeId << " (idx=" << targetNodeIdx << ", path=" << pathToTarget.size() << " nodes)");
            BOOST_TEST_MESSAGE("  Reconstructed: " << reconstructedSeeds.size() << " unique seeds");
            BOOST_TEST_MESSAGE("  Truth: " << truthSeeds.size() << " unique seeds");
            
            // Count differences
            int missing = 0, extra = 0, mismatch = 0;
            for (const auto& [h, c] : truthSeeds) {
                auto it = reconstructedSeeds.find(h);
                if (it == reconstructedSeeds.end()) missing++;
                else if (it->second != c) mismatch++;
            }
            for (const auto& [h, c] : reconstructedSeeds) {
                if (truthSeeds.find(h) == truthSeeds.end()) extra++;
            }
            BOOST_TEST_MESSAGE("  Missing=" << missing << " Extra=" << extra << " CountMismatch=" << mismatch);
        } else {
            BOOST_TEST_MESSAGE("PASS: Node " << nodeId << " (" << truthSeeds.size() << " seeds)");
        }
    }
    
    BOOST_TEST(totalFailures == 0, "Index DFS state check failed for " << totalFailures << " nodes");
    
    fs::remove(indexPath);
}

/**
 * CORE TEST: Verify placement BFS produces correct metrics at target node.
 * 
 * Simulates actual placement traversal and checks that accumulated metrics
 * match directly computed metrics from read seeds + target genome seeds.
 */
BOOST_AUTO_TEST_CASE(test_placement_bfs_metrics_match_truth) {
    const int k = 15;
    const int s = 8;
    const int l = 1;
    const int readLength = 150;
    const int numReads = 100;
    std::string indexPath = "data/test_placement_bfs.pmi";
    
    // Build index
    index_single_mode::IndexBuilder builder(T, k, s, 0, l, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    // Load index
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeParentCounts = indexRoot.getSeedChangeParentCounts();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Generate reads from truth node
    auto reads = generateReads(truthGenome, readLength, numReads, 42);
    BOOST_REQUIRE(!reads.empty());
    
    // Extract read seeds
    absl::flat_hash_map<size_t, int64_t> readSeeds;
    for (const auto& read : reads) {
        for (const auto& [hash, isReverse, isSyncmer, pos] : seeding::rollingSyncmers(read, k, s, false, 0, false)) {
            if (isSyncmer) {
                readSeeds[hash]++;
            }
        }
    }
    
    size_t readUniqueSeedCount = readSeeds.size();
    int64_t totalReadSeedFrequency = 0;
    double readMagnitudeSquared = 0.0;
    for (const auto& [hash, count] : readSeeds) {
        totalReadSeedFrequency += count;
        readMagnitudeSquared += count * count;
    }
    double readMagnitude = std::sqrt(readMagnitudeSquared);
    
    // Set up placement state
    placement::PlacementGlobalState state;
    for (const auto& [hash, count] : readSeeds) {
        state.seedFreqInReads[hash] = count;
    }
    state.readUniqueSeedCount = readUniqueSeedCount;
    state.totalReadSeedFrequency = totalReadSeedFrequency;
    state.readMagnitude = readMagnitude;
    
    // Find target node and build path
    int32_t targetNodeIdx = findNodeIndex(liteTree, truthNodeId);
    BOOST_REQUIRE(targetNodeIdx >= 0);
    
    std::vector<panmapUtils::LiteNode*> pathToTarget;
    panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
    while (current != nullptr) {
        pathToTarget.push_back(current);
        current = current->parent;
    }
    std::reverse(pathToTarget.begin(), pathToTarget.end());
    
    // Simulate BFS/DFS traversal: accumulate metrics along path
    placement::NodeMetrics metrics;
    for (const auto* node : pathToTarget) {
        uint32_t nodeIdx = node->nodeIndex;
        uint64_t startOffset = nodeChangeOffsets[nodeIdx];
        uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
            nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
        
        std::vector<std::tuple<uint64_t, int64_t, int64_t>> changes;
        for (uint64_t j = startOffset; j < endOffset; j++) {
            changes.emplace_back(seedChangeHashes[j], seedChangeParentCounts[j], seedChangeChildCounts[j]);
        }
        
        placement::NodeMetrics::computeChildMetrics(metrics, changes, state);
    }
    
    // Compute ground truth directly
    std::string targetGenome = T->getStringFromReference(truthNodeId, false);
    absl::flat_hash_map<size_t, int64_t> genomeSeeds;
    for (const auto& [hash, isReverse, isSyncmer, pos] : seeding::rollingSyncmers(targetGenome, k, s, false, 0, false)) {
        if (isSyncmer) {
            genomeSeeds[hash]++;
        }
    }
    
    auto groundTruth = GroundTruthMetrics::compute(readSeeds, genomeSeeds);
    
    // Compare BFS-accumulated metrics vs ground truth
    BOOST_TEST_MESSAGE("BFS metrics: unique=" << metrics.genomeUniqueSeedCount 
        << " total=" << metrics.genomeTotalSeedFrequency
        << " intersect=" << metrics.presenceIntersectionCount);
    BOOST_TEST_MESSAGE("Truth: unique=" << groundTruth.genomeUniqueSeedCount 
        << " total=" << groundTruth.genomeTotalSeedFrequency
        << " intersect=" << groundTruth.presenceIntersectionCount);
    
    BOOST_TEST(metrics.genomeUniqueSeedCount == static_cast<int64_t>(groundTruth.genomeUniqueSeedCount),
        "Unique seed count: BFS=" << metrics.genomeUniqueSeedCount << " truth=" << groundTruth.genomeUniqueSeedCount);
    BOOST_TEST(metrics.genomeTotalSeedFrequency == groundTruth.genomeTotalSeedFrequency,
        "Total frequency: BFS=" << metrics.genomeTotalSeedFrequency << " truth=" << groundTruth.genomeTotalSeedFrequency);
    BOOST_TEST(metrics.presenceIntersectionCount == static_cast<int64_t>(groundTruth.presenceIntersectionCount),
        "Intersection: BFS=" << metrics.presenceIntersectionCount << " truth=" << groundTruth.presenceIntersectionCount);
    
    double bfsJaccard = metrics.getJaccardScore(readUniqueSeedCount);
    BOOST_TEST(std::abs(bfsJaccard - groundTruth.jaccardScore) < 1e-9,
        "Jaccard: BFS=" << bfsJaccard << " truth=" << groundTruth.jaccardScore);
    
    fs::remove(indexPath);
}

/**
 * Test offset structure integrity
 */
BOOST_AUTO_TEST_CASE(test_node_change_offsets_structure) {
    std::string indexPath = "data/test_offsets_structure.pmi";
    
    index_single_mode::IndexBuilder builder(T, 15, 8, 0, 1, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto liteTree = indexRoot.getLiteTree();
    auto nodes = liteTree.getLiteNodes();
    
    size_t numNodes = nodes.size();
    size_t numOffsets = nodeChangeOffsets.size();
    size_t totalChanges = seedChangeHashes.size();
    
    BOOST_TEST_MESSAGE("Nodes: " << numNodes << ", Offsets: " << numOffsets << ", Total changes: " << totalChanges);
    
    // Offsets should be monotonically non-decreasing
    bool monotonic = true;
    for (size_t i = 1; i < numOffsets; i++) {
        if (nodeChangeOffsets[i] < nodeChangeOffsets[i-1]) {
            BOOST_TEST_MESSAGE("Non-monotonic at " << i);
            monotonic = false;
            break;
        }
    }
    BOOST_TEST(monotonic, "Offsets must be monotonically non-decreasing");
    
    fs::remove(indexPath);
}

/**
 * Test parent/child count consistency along any path
 */
BOOST_AUTO_TEST_CASE(test_seed_change_parent_child_consistency) {
    const int k = 15;
    const int s = 8;
    std::string indexPath = "data/test_parent_child.pmi";
    
    index_single_mode::IndexBuilder builder(T, k, s, 0, 1, false);
    builder.buildIndex();
    builder.writeIndex(indexPath);
    
    std::vector<uint8_t> data;
    BOOST_REQUIRE(panmap_zstd::decompressFromFile(indexPath, data));
    
    ::capnp::ReaderOptions opts;
    opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
    opts.nestingLimit = 1024;
    
    ::capnp::FlatArrayMessageReader reader(
        kj::ArrayPtr<const ::capnp::word>(
            reinterpret_cast<const ::capnp::word*>(data.data()),
            data.size() / sizeof(::capnp::word)),
        opts);
    
    auto indexRoot = reader.getRoot<LiteIndex>();
    panmapUtils::LiteTree liteTree;
    liteTree.initialize(reader.getRoot<LiteIndex>().getLiteTree());
    
    auto seedChangeHashes = indexRoot.getSeedChangeHashes();
    auto seedChangeParentCounts = indexRoot.getSeedChangeParentCounts();
    auto seedChangeChildCounts = indexRoot.getSeedChangeChildCounts();
    auto nodeChangeOffsets = indexRoot.getNodeChangeOffsets();
    
    // Test with multiple paths
    std::vector<std::string> testNodeIds = {truthNodeId, getRandomNodeId(10), getRandomNodeId(20)};
    
    int totalInconsistencies = 0;
    
    for (const auto& nodeId : testNodeIds) {
        int32_t targetNodeIdx = findNodeIndex(liteTree, nodeId);
        if (targetNodeIdx < 0) continue;
        
        std::vector<panmapUtils::LiteNode*> pathToTarget;
        panmapUtils::LiteNode* current = liteTree.dfsIndexToNode[targetNodeIdx];
        while (current != nullptr) {
            pathToTarget.push_back(current);
            current = current->parent;
        }
        std::reverse(pathToTarget.begin(), pathToTarget.end());
        
        // Track current seed counts as we traverse
        absl::flat_hash_map<uint64_t, int64_t> currentCounts;
        int pathInconsistencies = 0;
        
        for (const auto* node : pathToTarget) {
            uint32_t nodeIdx = node->nodeIndex;
            uint64_t startOffset = nodeChangeOffsets[nodeIdx];
            uint64_t endOffset = (nodeIdx + 1 < nodeChangeOffsets.size()) ? 
                nodeChangeOffsets[nodeIdx + 1] : seedChangeHashes.size();
            
            for (uint64_t j = startOffset; j < endOffset; j++) {
                uint64_t hash = seedChangeHashes[j];
                int64_t parentCount = seedChangeParentCounts[j];
                int64_t childCount = seedChangeChildCounts[j];
                
                int64_t expectedParent = currentCounts.count(hash) ? currentCounts[hash] : 0;
                if (parentCount != expectedParent) {
                    pathInconsistencies++;
                }
                
                if (childCount > 0) {
                    currentCounts[hash] = childCount;
                } else {
                    currentCounts.erase(hash);
                }
            }
        }
        
        if (pathInconsistencies > 0) {
            BOOST_TEST_MESSAGE("Node " << nodeId << ": " << pathInconsistencies << " parent/child inconsistencies");
        }
        totalInconsistencies += pathInconsistencies;
    }
    
    BOOST_TEST(totalInconsistencies == 0, "Found " << totalInconsistencies << " total inconsistencies");
    
    fs::remove(indexPath);
}

BOOST_AUTO_TEST_SUITE_END()

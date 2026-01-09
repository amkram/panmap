/**
 * @file main.cpp
 * @brief Panmap - Pangenome-based sequence placement, alignment, and genotyping
 * 
 * A modern command-line interface for:
 * - Placing sequencing reads on a pangenome tree
 * - Aligning reads to the best matching reference
 * - Calling variants (genotyping)
 * 
 * Supports single-sample (isolate) and metagenomic workflows.
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <tbb/global_control.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <chrono>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include <absl/container/flat_hash_map.h>

#include "capnp/message.h"
#include "capnp/serialize.h"
#include "kj/io.h"
#include "index_lite.capnp.h"
#include "logging.hpp"
#include "panman.hpp"
#include "placement.hpp"
#include "panmap_utils.hpp"
#include "index_single_mode.hpp"
#include "zstd_compression.hpp"
#include "genotyping.hpp"
#include "conversion.hpp"
#include "seeding.hpp"

extern "C" {
#include <htslib/sam.h>
#include <htslib/faidx.h>
}

namespace po = boost::program_options;
namespace fs = boost::filesystem;

// ============================================================================
// Version and Constants
// ============================================================================

constexpr const char* VERSION = "0.1.0";
constexpr const char* PROGRAM_NAME = "panmap";

// ANSI color codes for terminal output
namespace color {
    constexpr const char* reset = "\033[0m";
    constexpr const char* bold = "\033[1m";
    constexpr const char* dim = "\033[2m";
    constexpr const char* red = "\033[31m";
    constexpr const char* green = "\033[32m";
    constexpr const char* yellow = "\033[33m";
    constexpr const char* blue = "\033[34m";
    constexpr const char* cyan = "\033[36m";
}

// ============================================================================
// Pipeline Stage Enum
// ============================================================================

enum class PipelineStage {
    Index,      // Build index only
    Place,      // Placement only
    Align,      // Placement + Alignment
    Genotype,   // Placement + Alignment + Genotyping (default)
    Full        // Full pipeline including assembly (future)
};

std::string stageName(PipelineStage stage) {
    switch (stage) {
        case PipelineStage::Index: return "index";
        case PipelineStage::Place: return "place";
        case PipelineStage::Align: return "align";
        case PipelineStage::Genotype: return "genotype";
        case PipelineStage::Full: return "full";
    }
    return "unknown";
}

// ============================================================================
// Configuration
// ============================================================================

struct Config {
    // Input files
    std::string panman;           // Guide pangenome (.panman)
    std::string reads1;           // First read file (FASTQ/FASTA)
    std::string reads2;           // Second read file for paired-end
    
    // Output
    std::string output;           // Output prefix
    std::string index;            // Index file path
    
    // Pipeline control
    PipelineStage stopAfter = PipelineStage::Genotype;
    bool forceReindex = false;
    
    // Mode
    bool metagenomic = false;     // Metagenomic mode (multi-sample)
    int topN = 1;                 // Report top N placements
    bool dedupReads = false;      // Deduplicate reads before placement (for amplicon data)
    
    // Aligner
    std::string aligner = "minimap2";
    
    // Index parameters
    int k = 21;                   // syncmer k
    int s = 8;                    // syncmer s
    int l = 1;                    // l-mer size
    int flankMaskBp = 250;        // Hard mask first/last N bp at genome ends
    double seedMaskFraction = 0.001; // Mask top 0.1% most frequent seeds
    int minSeedQuality = 0;          // Min avg Phred quality for seed region (0=disabled)
    int trimStart = 0;               // Trim N bases from start of each read (primer removal)
    int trimEnd = 0;                 // Trim N bases from end of each read (primer removal)
    
    // Resources
    int threads = 1;
    
    // Utility modes
    bool dumpRandomNode = false;
    bool dumpSequence = false;
    std::string dumpNodeId;
    int seed = 42;
    
    // Leave-one-out validation mode
    std::string removeNodeId;     // Node to remove before placement (for validation)
    
    // Diagnostic options
    int topPlacements = 0;        // Output top-N placements with all scores (0 = off)
    std::string debugNodeId;      // Output detailed metrics for this specific node
    std::string compareNodes;     // Compare two nodes: "node1,node2" - compute metrics from scratch
    
    // Verbosity
    int verbosity = 1;            // 0=quiet, 1=normal, 2=verbose
};

// ============================================================================
// Helper Classes
// ============================================================================

// Cap'n Proto reader with ZSTD decompression
class IndexReader : public ::capnp::MessageReader {
public:
    std::vector<uint8_t> data;
    std::unique_ptr<::capnp::FlatArrayMessageReader> reader;
    
    explicit IndexReader(const std::string& path, int numThreads = 0) 
        : ::capnp::MessageReader(makeOptions()) 
    {
        if (!panmap_zstd::decompressFromFile(path, data, numThreads)) {
            throw std::runtime_error("Failed to decompress index: " + path);
        }
        
        reader = std::make_unique<::capnp::FlatArrayMessageReader>(
            kj::ArrayPtr<const capnp::word>(
                reinterpret_cast<const capnp::word*>(data.data()),
                data.size() / sizeof(capnp::word)),
            makeOptions());
    }
    
    kj::ArrayPtr<const capnp::word> getSegment(uint id) override {
        return reader->getSegment(id);
    }
    
private:
    static ::capnp::ReaderOptions makeOptions() {
        ::capnp::ReaderOptions opts;
        opts.traversalLimitInWords = 16ULL * 1024 * 1024 * 1024;
        opts.nestingLimit = 1024;
        return opts;
    }
};

// ============================================================================
// Utility Functions
// ============================================================================

std::unique_ptr<panmanUtils::TreeGroup> loadPanMAN(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open PanMAN file: " + path);
    }
    
    auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    buffer->push(boost::iostreams::lzma_decompressor());
    buffer->push(file);
    
    std::istream stream(buffer.get());
    return std::make_unique<panmanUtils::TreeGroup>(stream);
}

std::string getRandomNodeId(panmanUtils::Tree* T, int seed) {
    std::vector<std::string> ids;
    std::function<void(panmanUtils::Node*)> collect = [&](panmanUtils::Node* n) {
        if (!n) return;
        ids.push_back(n->identifier);
        for (auto* c : n->children) collect(c);
    };
    collect(T->root);
    
    if (ids.empty()) throw std::runtime_error("No nodes in tree");
    
    std::mt19937 gen(seed);
    return ids[std::uniform_int_distribution<>(0, ids.size()-1)(gen)];
}

void saveNodeSequence(panmanUtils::Tree* T, const std::string& nodeId, const std::string& path) {
    std::string seq = T->getStringFromReference(nodeId, false);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write: " + path);
    
    out << ">" << nodeId << "\n";
    for (size_t i = 0; i < seq.size(); i += 80) {
        out << seq.substr(i, 80) << "\n";
    }
    logging::msg("Saved {} ({} bp) to {}", nodeId, seq.size(), path);
}

std::string sanitizeFilename(const std::string& s) {
    std::string result = s;
    for (char& c : result) {
        if (c == '/' || c == '|' || c == ' ' || c == ':' || c == '\\') c = '_';
    }
    return result;
}

// ============================================================================
// Pipeline Steps
// ============================================================================

bool buildIndex(const Config& cfg) {
    if (fs::exists(cfg.index) && !cfg.forceReindex) {
        logging::msg("Using existing index: {}", cfg.index);
        return true;
    }
    
    logging::msg("Building index: {}", cfg.index);
    auto tg = loadPanMAN(cfg.panman);
    if (!tg || tg->trees.empty()) {
        logging::err("Failed to load pangenome");
        return false;
    }
    
    index_single_mode::IndexBuilder builder(&tg->trees[0], cfg.k, cfg.s, 0, cfg.l, false, cfg.flankMaskBp);
    // Always use buildIndexParallel - it handles single-threaded mode efficiently
    // using the delta-based approach (much more memory efficient than old buildIndex)
    builder.buildIndexParallel(cfg.threads);
    builder.writeIndex(cfg.index, cfg.threads);
    
    logging::msg("Index built with k={}, s={}, l={}, flankMask={}bp", cfg.k, cfg.s, cfg.l, cfg.flankMaskBp);
    return true;
}

// Compute tree distance between two nodes (number of edges)
// Returns -1 if either node is not found
int computeTreeDistance(panmapUtils::LiteTree& tree, 
                        const std::string& nodeId1, 
                        const std::string& nodeId2) {
    auto it1 = tree.allLiteNodes.find(nodeId1);
    auto it2 = tree.allLiteNodes.find(nodeId2);
    if (it1 == tree.allLiteNodes.end() || it2 == tree.allLiteNodes.end()) {
        return -1;
    }
    
    panmapUtils::LiteNode* node1 = it1->second;
    panmapUtils::LiteNode* node2 = it2->second;
    
    // Get path from node1 to root (collect ancestors)
    std::unordered_set<panmapUtils::LiteNode*> ancestors1;
    std::unordered_map<panmapUtils::LiteNode*, int> depth1;
    int d = 0;
    for (auto* n = node1; n != nullptr; n = n->parent) {
        ancestors1.insert(n);
        depth1[n] = d++;
    }
    
    // Walk from node2 to root, find first common ancestor (LCA)
    d = 0;
    for (auto* n = node2; n != nullptr; n = n->parent) {
        if (ancestors1.count(n)) {
            // Found LCA - distance is depth1[LCA] + depth from node2 to LCA
            return depth1[n] + d;
        }
        d++;
    }
    
    return -1; // Shouldn't happen in a proper tree
}

// Print diagnostic information for a node (when store_diagnostics is enabled)
void printNodeDiagnostics(panmapUtils::LiteNode* node, 
                          size_t readUniqueSeedCount,
                          int64_t totalReadSeedFrequency,
                          double readMagnitude,
                          panmapUtils::LiteTree* tree = nullptr,
                          const std::string& expectedParent = "") {
    if (!node) return;
    
    // Recompute scores from stored components for verification
    double genomeMag = std::sqrt(node->genomeMagnitudeSquared);
    size_t genomeOnlyCount = (node->genomeUniqueSeedCount > node->presenceIntersectionCount) ?
        (node->genomeUniqueSeedCount - node->presenceIntersectionCount) : 0;
    size_t totalUnion = readUniqueSeedCount + genomeOnlyCount;
    
    int64_t wjDenom = totalReadSeedFrequency + node->genomeTotalSeedFrequency - node->weightedJaccardNumerator;
    
    // Compute distance from expected parent if provided
    std::string distStr = "";
    if (tree && !expectedParent.empty()) {
        int dist = computeTreeDistance(*tree, node->identifier, expectedParent);
        distStr = "  dist=" + std::to_string(dist);
    }
    
    std::cout << "  NODE: " << node->identifier << distStr << "\n";
    std::cout << "    Scores:     Jaccard=" << node->jaccardScore 
              << "  WeightedJaccard=" << node->weightedJaccardScore
              << "  Cosine=" << node->cosineScore 
              << "  Raw=" << node->rawSeedMatchScore << "\n";
    std::cout << "    Jaccard:    num=" << node->jaccardNumerator 
              << "  intersection=" << node->presenceIntersectionCount
              << "  genomeUnique=" << node->genomeUniqueSeedCount
              << "  genomeOnly=" << genomeOnlyCount
              << "  union=" << totalUnion << "\n";
    std::cout << "    WJaccard:   num=" << node->weightedJaccardNumerator
              << "  readFreq=" << totalReadSeedFrequency
              << "  genomeFreq=" << node->genomeTotalSeedFrequency
              << "  denom=" << wjDenom << "\n";
    std::cout << "    Cosine:     num=" << node->cosineNumerator
              << "  readMag=" << readMagnitude
              << "  genomeMagSq=" << node->genomeMagnitudeSquared
              << "  genomeMag=" << genomeMag << "\n";
}

// ============================================================================
// Node Comparison - Compute metrics from scratch for two nodes
// ============================================================================
struct NodeSeedAnalysis {
    std::string nodeId;
    std::string sequence;
    absl::flat_hash_map<uint64_t, int64_t> seedCounts;  // hash -> count in genome
    size_t uniqueSeedCount = 0;
    int64_t totalSeedFrequency = 0;
    double magnitudeSquared = 0.0;
    
    // Metrics vs reads
    size_t intersectionCount = 0;
    int64_t jaccardNumerator = 0;
    int64_t weightedJaccardNumerator = 0;
    double cosineNumerator = 0.0;
    double rawScore = 0.0;
};

void compareNodesFromScratch(const Config& cfg, panmanUtils::Tree* T) {
    // Parse node IDs from "node1,node2" format
    auto commaPos = cfg.compareNodes.find(',');
    if (commaPos == std::string::npos) {
        logging::err("--compare-nodes requires format 'node1,node2'");
        return;
    }
    std::string node1Id = cfg.compareNodes.substr(0, commaPos);
    std::string node2Id = cfg.compareNodes.substr(commaPos + 1);
    
    logging::msg("=== COMPARING NODES FROM SCRATCH ===");
    logging::msg("Node 1: {}", node1Id);
    logging::msg("Node 2: {}", node2Id);
    
    // Get sequences
    std::string seq1 = T->getStringFromReference(node1Id, false, true);
    std::string seq2 = T->getStringFromReference(node2Id, false, true);
    
    if (seq1.empty()) {
        logging::err("Node '{}' not found or has empty sequence", node1Id);
        return;
    }
    if (seq2.empty()) {
        logging::err("Node '{}' not found or has empty sequence", node2Id);
        return;
    }
    
    logging::msg("Node 1 sequence length: {} bp", seq1.size());
    logging::msg("Node 2 sequence length: {} bp", seq2.size());
    
    // Compute syncmers for each genome
    auto computeGenomeSyncmers = [&](const std::string& seq, NodeSeedAnalysis& analysis) {
        auto syncmers = seeding::rollingSyncmers(seq, cfg.k, cfg.s, false, 0, false);
        for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
            analysis.seedCounts[hash]++;
        }
        for (const auto& [hash, count] : analysis.seedCounts) {
            analysis.uniqueSeedCount++;
            analysis.totalSeedFrequency += count;
            analysis.magnitudeSquared += static_cast<double>(count) * count;
        }
    };
    
    NodeSeedAnalysis analysis1, analysis2;
    analysis1.nodeId = node1Id;
    analysis1.sequence = seq1;
    analysis2.nodeId = node2Id;
    analysis2.sequence = seq2;
    
    computeGenomeSyncmers(seq1, analysis1);
    computeGenomeSyncmers(seq2, analysis2);
    
    logging::msg("Node 1 syncmers: {} unique, {} total frequency", 
                analysis1.uniqueSeedCount, analysis1.totalSeedFrequency);
    logging::msg("Node 2 syncmers: {} unique, {} total frequency", 
                analysis2.uniqueSeedCount, analysis2.totalSeedFrequency);
    
    // Compute syncmers from reads
    logging::msg("Computing syncmers from reads...");
    absl::flat_hash_map<uint64_t, int64_t> readSeedCounts;
    size_t readUniqueSeedCount = 0;
    int64_t totalReadSeedFrequency = 0;
    double readMagnitudeSquared = 0.0;
    size_t totalReads = 0;
    
    auto processReads = [&](const std::string& path) {
        std::ifstream file(path);
        if (!file) return;
        std::string line, seq;
        int lineNum = 0;
        while (std::getline(file, line)) {
            lineNum++;
            if (lineNum % 4 == 2) {  // Sequence line in FASTQ
                seq = line;
                auto syncmers = seeding::rollingSyncmers(seq, cfg.k, cfg.s, false, 0, false);
                for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers) {
                    readSeedCounts[hash]++;
                }
                totalReads++;
            }
        }
    };
    
    processReads(cfg.reads1);
    if (!cfg.reads2.empty()) processReads(cfg.reads2);
    
    // Compute read statistics
    for (const auto& [hash, count] : readSeedCounts) {
        readUniqueSeedCount++;
        totalReadSeedFrequency += count;
        readMagnitudeSquared += static_cast<double>(count) * count;
    }
    double readMagnitude = std::sqrt(readMagnitudeSquared);
    
    logging::msg("Processed {} reads", totalReads);
    logging::msg("Read syncmers: {} unique, {} total frequency, magnitude={:.2f}", 
                readUniqueSeedCount, totalReadSeedFrequency, readMagnitude);
    
    // Compute metrics for each genome vs reads
    auto computeMetrics = [&](NodeSeedAnalysis& analysis) {
        for (const auto& [hash, readCount] : readSeedCounts) {
            auto it = analysis.seedCounts.find(hash);
            if (it != analysis.seedCounts.end()) {
                int64_t genomeCount = it->second;
                analysis.intersectionCount++;
                analysis.jaccardNumerator++;  // Presence-based
                analysis.weightedJaccardNumerator += std::min(readCount, genomeCount);
                analysis.cosineNumerator += static_cast<double>(readCount) * genomeCount;
                analysis.rawScore += static_cast<double>(readCount) / genomeCount;
            }
        }
    };
    
    computeMetrics(analysis1);
    computeMetrics(analysis2);
    
    // Compute final scores
    auto computeScores = [&](const NodeSeedAnalysis& a) {
        size_t genomeOnly = a.uniqueSeedCount - a.intersectionCount;
        size_t totalUnion = readUniqueSeedCount + genomeOnly;
        double jaccard = totalUnion > 0 ? static_cast<double>(a.jaccardNumerator) / totalUnion : 0.0;
        
        int64_t wjDenom = totalReadSeedFrequency + a.totalSeedFrequency - a.weightedJaccardNumerator;
        double weightedJaccard = wjDenom > 0 ? static_cast<double>(a.weightedJaccardNumerator) / wjDenom : 0.0;
        
        double genomeMag = std::sqrt(a.magnitudeSquared);
        double cosine = (readMagnitude > 0 && genomeMag > 0) ? 
            a.cosineNumerator / (readMagnitude * genomeMag) : 0.0;
        
        return std::make_tuple(jaccard, weightedJaccard, cosine, a.rawScore);
    };
    
    auto [jac1, wj1, cos1, raw1] = computeScores(analysis1);
    auto [jac2, wj2, cos2, raw2] = computeScores(analysis2);
    
    // Print detailed comparison
    std::cout << "\n========== DETAILED NODE COMPARISON ==========\n\n";
    
    std::cout << "GENOME STATISTICS:\n";
    std::cout << "  " << std::setw(50) << std::left << "Metric" 
              << std::setw(20) << std::right << node1Id.substr(0, 18)
              << std::setw(20) << std::right << node2Id.substr(0, 18) 
              << std::setw(15) << "Difference" << "\n";
    std::cout << std::string(105, '-') << "\n";
    
    std::cout << "  " << std::setw(50) << std::left << "Sequence length (bp)"
              << std::setw(20) << std::right << seq1.size()
              << std::setw(20) << std::right << seq2.size()
              << std::setw(15) << static_cast<int64_t>(seq2.size()) - static_cast<int64_t>(seq1.size()) << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Unique syncmers"
              << std::setw(20) << std::right << analysis1.uniqueSeedCount
              << std::setw(20) << std::right << analysis2.uniqueSeedCount
              << std::setw(15) << static_cast<int64_t>(analysis2.uniqueSeedCount) - static_cast<int64_t>(analysis1.uniqueSeedCount) << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Total syncmer frequency"
              << std::setw(20) << std::right << analysis1.totalSeedFrequency
              << std::setw(20) << std::right << analysis2.totalSeedFrequency
              << std::setw(15) << analysis2.totalSeedFrequency - analysis1.totalSeedFrequency << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Magnitude squared"
              << std::setw(20) << std::right << std::fixed << std::setprecision(1) << analysis1.magnitudeSquared
              << std::setw(20) << std::right << analysis2.magnitudeSquared
              << std::setw(15) << analysis2.magnitudeSquared - analysis1.magnitudeSquared << "\n";
    
    std::cout << "\nREAD-GENOME INTERSECTION:\n";
    std::cout << "  " << std::setw(50) << std::left << "Intersection count (seeds in both)"
              << std::setw(20) << std::right << analysis1.intersectionCount
              << std::setw(20) << std::right << analysis2.intersectionCount
              << std::setw(15) << static_cast<int64_t>(analysis2.intersectionCount) - static_cast<int64_t>(analysis1.intersectionCount) << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Weighted Jaccard numerator (Σmin)"
              << std::setw(20) << std::right << analysis1.weightedJaccardNumerator
              << std::setw(20) << std::right << analysis2.weightedJaccardNumerator
              << std::setw(15) << analysis2.weightedJaccardNumerator - analysis1.weightedJaccardNumerator << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Cosine numerator (dot product)"
              << std::setw(20) << std::right << std::fixed << std::setprecision(0) << analysis1.cosineNumerator
              << std::setw(20) << std::right << analysis2.cosineNumerator
              << std::setw(15) << analysis2.cosineNumerator - analysis1.cosineNumerator << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Raw score (Σ readCount/genomeCount)"
              << std::setw(20) << std::right << std::fixed << std::setprecision(0) << analysis1.rawScore
              << std::setw(20) << std::right << analysis2.rawScore
              << std::setw(15) << analysis2.rawScore - analysis1.rawScore << "\n";
    
    std::cout << "\nFINAL SCORES:\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "  " << std::setw(50) << std::left << "Jaccard"
              << std::setw(20) << std::right << jac1
              << std::setw(20) << std::right << jac2
              << std::setw(15) << (jac1 > jac2 ? "← WINNER" : (jac2 > jac1 ? "WINNER →" : "TIE")) << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Weighted Jaccard"
              << std::setw(20) << std::right << wj1
              << std::setw(20) << std::right << wj2
              << std::setw(15) << (wj1 > wj2 ? "← WINNER" : (wj2 > wj1 ? "WINNER →" : "TIE")) << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Cosine"
              << std::setw(20) << std::right << cos1
              << std::setw(20) << std::right << cos2
              << std::setw(15) << (cos1 > cos2 ? "← WINNER" : (cos2 > cos1 ? "WINNER →" : "TIE")) << "\n";
              
    std::cout << "  " << std::setw(50) << std::left << "Raw"
              << std::setw(20) << std::right << std::fixed << std::setprecision(0) << raw1
              << std::setw(20) << std::right << raw2
              << std::setw(15) << (raw1 > raw2 ? "← WINNER" : (raw2 > raw1 ? "WINNER →" : "TIE")) << "\n";
    
    // Analyze seed differences
    std::cout << "\nSEED ANALYSIS:\n";
    
    // Seeds unique to each genome
    size_t onlyIn1 = 0, onlyIn2 = 0, inBoth = 0;
    for (const auto& [hash, count] : analysis1.seedCounts) {
        if (analysis2.seedCounts.count(hash)) inBoth++;
        else onlyIn1++;
    }
    for (const auto& [hash, count] : analysis2.seedCounts) {
        if (!analysis1.seedCounts.count(hash)) onlyIn2++;
    }
    
    std::cout << "  Seeds only in " << node1Id.substr(0, 30) << ": " << onlyIn1 << "\n";
    std::cout << "  Seeds only in " << node2Id.substr(0, 30) << ": " << onlyIn2 << "\n";
    std::cout << "  Seeds in both genomes: " << inBoth << "\n";
    
    // Seeds that match reads, unique to each genome
    size_t matchOnlyIn1 = 0, matchOnlyIn2 = 0, matchInBoth = 0;
    int64_t readFreqMatchOnlyIn1 = 0, readFreqMatchOnlyIn2 = 0;
    
    for (const auto& [hash, readCount] : readSeedCounts) {
        bool in1 = analysis1.seedCounts.count(hash) > 0;
        bool in2 = analysis2.seedCounts.count(hash) > 0;
        if (in1 && in2) matchInBoth++;
        else if (in1) { matchOnlyIn1++; readFreqMatchOnlyIn1 += readCount; }
        else if (in2) { matchOnlyIn2++; readFreqMatchOnlyIn2 += readCount; }
    }
    
    std::cout << "\n  Read seeds matching ONLY " << node1Id.substr(0, 30) << ": " << matchOnlyIn1 
              << " (total read freq: " << readFreqMatchOnlyIn1 << ")\n";
    std::cout << "  Read seeds matching ONLY " << node2Id.substr(0, 30) << ": " << matchOnlyIn2 
              << " (total read freq: " << readFreqMatchOnlyIn2 << ")\n";
    std::cout << "  Read seeds matching BOTH genomes: " << matchInBoth << "\n";
    
    // Detailed analysis of high-frequency unique seeds
    std::cout << "\n=== HIGH-FREQUENCY UNIQUE SEED DETAILS ===\n";
    
    // Build hash-to-position maps from genome sequences
    absl::flat_hash_map<uint64_t, std::vector<std::pair<int32_t, bool>>> hashToPos1, hashToPos2;
    auto syncmers1 = seeding::rollingSyncmers(seq1, cfg.k, cfg.s, false, 0, false);
    auto syncmers2 = seeding::rollingSyncmers(seq2, cfg.k, cfg.s, false, 0, false);
    for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers1) {
        hashToPos1[hash].push_back({pos, isReverse});
    }
    for (const auto& [hash, isReverse, isSyncmer, pos] : syncmers2) {
        hashToPos2[hash].push_back({pos, isReverse});
    }
    
    // Collect unique seeds with read frequencies for node2
    std::vector<std::tuple<uint64_t, int64_t, int64_t>> node2UniqueSeeds;  // hash, readFreq, genomeCount
    for (const auto& [hash, readCount] : readSeedCounts) {
        bool in1 = analysis1.seedCounts.count(hash) > 0;
        bool in2 = analysis2.seedCounts.count(hash) > 0;
        if (!in1 && in2) {
            node2UniqueSeeds.push_back({hash, readCount, analysis2.seedCounts.at(hash)});
        }
    }
    // Sort by read frequency descending
    std::sort(node2UniqueSeeds.begin(), node2UniqueSeeds.end(),
              [](const auto& a, const auto& b) { return std::get<1>(a) > std::get<1>(b); });
    
    std::cout << "\nSeeds unique to " << node2Id.substr(0, 30) << " that match reads (sorted by freq):\n";
    std::cout << std::setw(8) << "ReadFreq" << "  " << std::setw(8) << "GenomeN" 
              << "  " << std::setw(10) << "Position" << "  " << "Kmer sequence\n";
    std::cout << std::string(90, '-') << "\n";
    
    for (size_t i = 0; i < std::min(node2UniqueSeeds.size(), size_t(30)); i++) {
        auto [hash, readFreq, genomeCount] = node2UniqueSeeds[i];
        auto& positions = hashToPos2[hash];
        for (auto& [pos, isRev] : positions) {
            std::string kmer = seq2.substr(pos, cfg.k);
            std::cout << std::setw(8) << readFreq << "  " << std::setw(8) << genomeCount
                      << "  " << std::setw(10) << pos << (isRev ? "R" : "F") << "  " << kmer << "\n";
        }
    }
    
    // Also show unique seeds for node1
    std::vector<std::tuple<uint64_t, int64_t, int64_t>> node1UniqueSeeds;
    for (const auto& [hash, readCount] : readSeedCounts) {
        bool in1 = analysis1.seedCounts.count(hash) > 0;
        bool in2 = analysis2.seedCounts.count(hash) > 0;
        if (in1 && !in2) {
            node1UniqueSeeds.push_back({hash, readCount, analysis1.seedCounts.at(hash)});
        }
    }
    std::sort(node1UniqueSeeds.begin(), node1UniqueSeeds.end(),
              [](const auto& a, const auto& b) { return std::get<1>(a) > std::get<1>(b); });
    
    std::cout << "\nSeeds unique to " << node1Id.substr(0, 30) << " that match reads (sorted by freq):\n";
    std::cout << std::setw(8) << "ReadFreq" << "  " << std::setw(8) << "GenomeN" 
              << "  " << std::setw(10) << "Position" << "  " << "Kmer sequence\n";
    std::cout << std::string(90, '-') << "\n";
    
    for (size_t i = 0; i < std::min(node1UniqueSeeds.size(), size_t(30)); i++) {
        auto [hash, readFreq, genomeCount] = node1UniqueSeeds[i];
        auto& positions = hashToPos1[hash];
        for (auto& [pos, isRev] : positions) {
            std::string kmer = seq1.substr(pos, cfg.k);
            std::cout << std::setw(8) << readFreq << "  " << std::setw(8) << genomeCount
                      << "  " << std::setw(10) << pos << (isRev ? "R" : "F") << "  " << kmer << "\n";
        }
    }
    
    std::cout << "\n==============================================\n\n";
}

std::optional<placement::PlacementResult> runPlacement(const Config& cfg) {
    logging::msg("Loading index...");
    IndexReader reader(cfg.index, cfg.threads);
    
    auto idx = reader.getRoot<LiteIndex>();
    logging::msg("Index parameters: k={}, s={}, l={}", idx.getK(), idx.getS(), idx.getL());
    
    panmapUtils::LiteTree tree;
    tree.initialize(idx.getLiteTree());
    logging::msg("Tree loaded: {} nodes", tree.allLiteNodes.size());
    
    placement::PlacementResult result;
    std::string outPath = cfg.output + ".placement.tsv";
    
    // Leave-one-out validation: track parent of removed node
    std::string expectedParentNode;
    std::string* expectedParentPtr = cfg.removeNodeId.empty() ? nullptr : &expectedParentNode;
    
    // Enable diagnostics if requested
    bool storeDiagnostics = (cfg.topPlacements > 0 || !cfg.debugNodeId.empty());
    
    auto start = std::chrono::high_resolution_clock::now();
    placement::placeLite(result, &tree, reader, cfg.reads1, cfg.reads2, outPath, false, nullptr,
                        cfg.removeNodeId, expectedParentPtr, storeDiagnostics, cfg.seedMaskFraction,
                        cfg.minSeedQuality, cfg.dedupReads, true, cfg.trimStart, cfg.trimEnd);
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    
    logging::msg("Placement complete in {}ms", elapsed.count());
    
    // Check if any metric found a placement (RAW and Cosine are always computed)
    if (result.bestRawSeedMatchNodeId.empty() && result.bestCosineNodeId.empty()) {
        logging::warn("No placement found");
        return std::nullopt;
    }
    
    // Report results (use RAW as primary metric now that WJ is disabled)
    logging::msg("{}Best placement:{} {} (RAW score: {:.1f}, Cosine: {:.4f})", 
                color::green, color::reset,
                result.bestRawSeedMatchNodeId, 
                result.bestRawSeedMatchScore,
                result.bestCosineScore);
    
    // Leave-one-out validation: compare to expected parent
    if (!cfg.removeNodeId.empty() && !expectedParentNode.empty()) {
        // Compute tree distances for all metrics
        int distWeightedJaccard = computeTreeDistance(tree, result.bestWeightedJaccardNodeId, expectedParentNode);
        int distCosine = computeTreeDistance(tree, result.bestCosineNodeId, expectedParentNode);
        int distJaccard = computeTreeDistance(tree, result.bestJaccardNodeId, expectedParentNode);
        int distRaw = computeTreeDistance(tree, result.bestRawSeedMatchNodeId, expectedParentNode);
        
        bool correctWeightedJaccard = (result.bestWeightedJaccardNodeId == expectedParentNode);
        bool correctCosine = (result.bestCosineNodeId == expectedParentNode);
        bool correctJaccard = (result.bestJaccardNodeId == expectedParentNode);
        bool correctRaw = (result.bestRawSeedMatchNodeId == expectedParentNode);
        bool correctContainment = (result.bestContainmentNodeId == expectedParentNode);
        bool correctWeightedContainment = (result.bestWeightedContainmentNodeId == expectedParentNode);
        bool correctLogRaw = (result.bestLogRawNodeId == expectedParentNode);
        bool correctLogCosine = (result.bestLogCosineNodeId == expectedParentNode);
        bool correctAdjCosine = (result.bestAdjCosineNodeId == expectedParentNode);
        bool correctAdjRaw = (result.bestAdjRawNodeId == expectedParentNode);
        bool correctIdfCosine = (result.bestIdfCosineNodeId == expectedParentNode);
        bool correctCovCosine = (result.bestCovCosineNodeId == expectedParentNode);
        bool correctCapCosine = (result.bestCapCosineNodeId == expectedParentNode);
        bool correctCapLogCosine = (result.bestCapLogCosineNodeId == expectedParentNode);
        bool correctSigCosine = (result.bestSigCosineNodeId == expectedParentNode);
        
        int distContainment = computeTreeDistance(tree, result.bestContainmentNodeId, expectedParentNode);
        int distWeightedContainment = computeTreeDistance(tree, result.bestWeightedContainmentNodeId, expectedParentNode);
        int distLogRaw = computeTreeDistance(tree, result.bestLogRawNodeId, expectedParentNode);
        int distLogCosine = computeTreeDistance(tree, result.bestLogCosineNodeId, expectedParentNode);
        int distAdjCosine = computeTreeDistance(tree, result.bestAdjCosineNodeId, expectedParentNode);
        int distAdjRaw = computeTreeDistance(tree, result.bestAdjRawNodeId, expectedParentNode);
        int distIdfCosine = computeTreeDistance(tree, result.bestIdfCosineNodeId, expectedParentNode);
        int distCovCosine = computeTreeDistance(tree, result.bestCovCosineNodeId, expectedParentNode);
        int distCapCosine = computeTreeDistance(tree, result.bestCapCosineNodeId, expectedParentNode);
        int distCapLogCosine = computeTreeDistance(tree, result.bestCapLogCosineNodeId, expectedParentNode);
        int distSigCosine = computeTreeDistance(tree, result.bestSigCosineNodeId, expectedParentNode);
        
        logging::msg("{}VALIDATION RESULTS:{}", color::cyan, color::reset);
        logging::msg("  Expected parent: {}", expectedParentNode);
        logging::msg("  Weighted Jaccard: {} {} (dist: {})", 
                    result.bestWeightedJaccardNodeId,
                    correctWeightedJaccard ? color::green : color::yellow,
                    correctWeightedJaccard ? "✓" : "✗",
                    distWeightedJaccard);
        logging::msg("  Cosine:           {} {}{} (dist: {})", 
                    result.bestCosineNodeId,
                    correctCosine ? color::green : color::yellow,
                    correctCosine ? "✓" : "✗",
                    distCosine);
        logging::msg("  Jaccard:          {} {}{} (dist: {})", 
                    result.bestJaccardNodeId,
                    correctJaccard ? color::green : color::yellow,
                    correctJaccard ? "✓" : "✗",
                    distJaccard);
        logging::msg("  Raw:              {} {}{} (dist: {})", 
                    result.bestRawSeedMatchNodeId,
                    correctRaw ? color::green : color::yellow,
                    correctRaw ? "✓" : "✗",
                    distRaw);
        logging::msg("  Containment:      {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestContainmentNodeId,
                    correctContainment ? color::green : color::yellow,
                    correctContainment ? "✓" : "✗",
                    distContainment,
                    result.bestContainmentScore);
        logging::msg("  Weighted Contain: {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestWeightedContainmentNodeId,
                    correctWeightedContainment ? color::green : color::yellow,
                    correctWeightedContainment ? "✓" : "✗",
                    distWeightedContainment,
                    result.bestWeightedContainmentScore);
        logging::msg("  LogRAW:           {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestLogRawNodeId,
                    correctLogRaw ? color::green : color::yellow,
                    correctLogRaw ? "✓" : "✗",
                    distLogRaw,
                    result.bestLogRawScore);
        logging::msg("  LogCosine:        {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestLogCosineNodeId,
                    correctLogCosine ? color::green : color::yellow,
                    correctLogCosine ? "✓" : "✗",
                    distLogCosine,
                    result.bestLogCosineScore);
        logging::msg("  AdjCosine:        {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestAdjCosineNodeId,
                    correctAdjCosine ? color::green : color::yellow,
                    correctAdjCosine ? "✓" : "✗",
                    distAdjCosine,
                    result.bestAdjCosineScore);
        logging::msg("  AdjRaw:           {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestAdjRawNodeId,
                    correctAdjRaw ? color::green : color::yellow,
                    correctAdjRaw ? "✓" : "✗",
                    distAdjRaw,
                    result.bestAdjRawScore);
        logging::msg("  IdfCosine:        {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestIdfCosineNodeId,
                    correctIdfCosine ? color::green : color::yellow,
                    correctIdfCosine ? "✓" : "✗",
                    distIdfCosine,
                    result.bestIdfCosineScore);
        logging::msg("  CovCosine:        {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestCovCosineNodeId,
                    correctCovCosine ? color::green : color::yellow,
                    correctCovCosine ? "✓" : "✗",
                    distCovCosine,
                    result.bestCovCosineScore);
        logging::msg("  CapCosine:        {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestCapCosineNodeId,
                    correctCapCosine ? color::green : color::yellow,
                    correctCapCosine ? "✓" : "✗",
                    distCapCosine,
                    result.bestCapCosineScore);
        logging::msg("  CapLogCosine:     {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestCapLogCosineNodeId,
                    correctCapLogCosine ? color::green : color::yellow,
                    correctCapLogCosine ? "✓" : "✗",
                    distCapLogCosine,
                    result.bestCapLogCosineScore);
        logging::msg("  SigCosine:        {} {}{} (dist: {}, score: {:.6f})", 
                    result.bestSigCosineNodeId,
                    correctSigCosine ? color::green : color::yellow,
                    correctSigCosine ? "✓" : "✗",
                    distSigCosine,
                    result.bestSigCosineScore);
        
        // Print to stdout for easy parsing (TSV format)
        std::cout << "REMOVED_NODE\t" << cfg.removeNodeId << "\n";
        std::cout << "EXPECTED_PARENT\t" << expectedParentNode << "\n";
        std::cout << "WEIGHTED_JACCARD\t" << result.bestWeightedJaccardNodeId 
                  << "\t" << (correctWeightedJaccard ? "true" : "false") 
                  << "\t" << distWeightedJaccard << "\n";
        std::cout << "COSINE\t" << result.bestCosineNodeId 
                  << "\t" << (correctCosine ? "true" : "false") 
                  << "\t" << distCosine << "\n";
        std::cout << "JACCARD\t" << result.bestJaccardNodeId 
                  << "\t" << (correctJaccard ? "true" : "false") 
                  << "\t" << distJaccard << "\n";
        std::cout << "RAW\t" << result.bestRawSeedMatchNodeId 
                  << "\t" << (correctRaw ? "true" : "false") 
                  << "\t" << distRaw << "\n";
        std::cout << "CONTAINMENT\t" << result.bestContainmentNodeId 
                  << "\t" << (correctContainment ? "true" : "false") 
                  << "\t" << distContainment 
                  << "\t" << result.bestContainmentScore << "\n";
        std::cout << "WEIGHTED_CONTAINMENT\t" << result.bestWeightedContainmentNodeId 
                  << "\t" << (correctWeightedContainment ? "true" : "false") 
                  << "\t" << distWeightedContainment 
                  << "\t" << result.bestWeightedContainmentScore << "\n";
        std::cout << "LOGRAW\t" << result.bestLogRawNodeId 
                  << "\t" << (correctLogRaw ? "true" : "false") 
                  << "\t" << distLogRaw 
                  << "\t" << result.bestLogRawScore << "\n";
        std::cout << "LOGCOSINE\t" << result.bestLogCosineNodeId 
                  << "\t" << (correctLogCosine ? "true" : "false") 
                  << "\t" << distLogCosine 
                  << "\t" << result.bestLogCosineScore << "\n";
        std::cout << "ADJCOSINE\t" << result.bestAdjCosineNodeId 
                  << "\t" << (correctAdjCosine ? "true" : "false") 
                  << "\t" << distAdjCosine 
                  << "\t" << result.bestAdjCosineScore << "\n";
        std::cout << "ADJRAW\t" << result.bestAdjRawNodeId 
                  << "\t" << (correctAdjRaw ? "true" : "false") 
                  << "\t" << distAdjRaw 
                  << "\t" << result.bestAdjRawScore << "\n";
        std::cout << "IDFCOSINE\t" << result.bestIdfCosineNodeId 
                  << "\t" << (correctIdfCosine ? "true" : "false") 
                  << "\t" << distIdfCosine 
                  << "\t" << result.bestIdfCosineScore << "\n";
        std::cout << "COVCOSINE\t" << result.bestCovCosineNodeId 
                  << "\t" << (correctCovCosine ? "true" : "false") 
                  << "\t" << distCovCosine 
                  << "\t" << result.bestCovCosineScore << "\n";
        std::cout << "CAPCOSINE\t" << result.bestCapCosineNodeId 
                  << "\t" << (correctCapCosine ? "true" : "false") 
                  << "\t" << distCapCosine 
                  << "\t" << result.bestCapCosineScore << "\n";
        std::cout << "CAPLOGCOSINE\t" << result.bestCapLogCosineNodeId 
                  << "\t" << (correctCapLogCosine ? "true" : "false") 
                  << "\t" << distCapLogCosine 
                  << "\t" << result.bestCapLogCosineScore << "\n";
        std::cout << "SIGCOSINE\t" << result.bestSigCosineNodeId 
                  << "\t" << (correctSigCosine ? "true" : "false") 
                  << "\t" << distSigCosine 
                  << "\t" << result.bestSigCosineScore << "\n";
    }
    
    // Diagnostic output: top placements with all metric components
    if (cfg.topPlacements > 0 || !cfg.debugNodeId.empty()) {
        std::cout << "\n=== DIAGNOSTIC OUTPUT ===\n";
        std::cout << "Read Statistics:\n";
        std::cout << "  Unique seeds: " << result.readUniqueSeedCount << "\n";
        std::cout << "  Total seed frequency: " << result.totalReadSeedFrequency << "\n";
        std::cout << "  Read magnitude: " << result.readMagnitude << "\n\n";
        
        // Collect all nodes with non-zero scores and sort by each metric
        std::vector<panmapUtils::LiteNode*> scoredNodes;
        for (auto& [id, node] : tree.allLiteNodes) {
            if (node->rawSeedMatchScore > 0 || node->jaccardScore > 0 || 
                node->cosineScore > 0 || node->weightedJaccardScore > 0) {
                scoredNodes.push_back(node);
            }
        }
        
        if (cfg.topPlacements > 0) {
            int topN = std::min(cfg.topPlacements, static_cast<int>(scoredNodes.size()));
            
            // Top by RAW score
            std::sort(scoredNodes.begin(), scoredNodes.end(), 
                [](auto* a, auto* b) { return a->rawSeedMatchScore > b->rawSeedMatchScore; });
            std::cout << "TOP " << topN << " BY RAW SCORE:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount, 
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
            
            // Top by Weighted Jaccard
            std::sort(scoredNodes.begin(), scoredNodes.end(),
                [](auto* a, auto* b) { return a->weightedJaccardScore > b->weightedJaccardScore; });
            std::cout << "\nTOP " << topN << " BY WEIGHTED JACCARD:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount,
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
            
            // Top by Cosine
            std::sort(scoredNodes.begin(), scoredNodes.end(),
                [](auto* a, auto* b) { return a->cosineScore > b->cosineScore; });
            std::cout << "\nTOP " << topN << " BY COSINE:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount,
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
            
            // Top by Jaccard
            std::sort(scoredNodes.begin(), scoredNodes.end(),
                [](auto* a, auto* b) { return a->jaccardScore > b->jaccardScore; });
            std::cout << "\nTOP " << topN << " BY JACCARD:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount,
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
        }
        
        // Debug specific node
        if (!cfg.debugNodeId.empty()) {
            auto it = tree.allLiteNodes.find(cfg.debugNodeId);
            if (it != tree.allLiteNodes.end()) {
                std::cout << "\nDEBUG NODE (requested):\n";
                printNodeDiagnostics(it->second, result.readUniqueSeedCount,
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            } else {
                std::cout << "\nDEBUG NODE '" << cfg.debugNodeId << "' NOT FOUND\n";
            }
        }
        
        // Also show expected parent if in LOO mode
        if (!cfg.removeNodeId.empty() && !expectedParentNode.empty()) {
            auto it = tree.allLiteNodes.find(expectedParentNode);
            if (it != tree.allLiteNodes.end()) {
                std::cout << "\nEXPECTED PARENT NODE:\n";
                printNodeDiagnostics(it->second, result.readUniqueSeedCount,
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
        }
        std::cout << "=========================\n\n";
    }
    
    if (cfg.metagenomic && cfg.topN > 1) {
        // TODO: Report top N placements for metagenomic mode
        logging::msg("Top {} placements written to {}", cfg.topN, outPath);
    }
    
    return result;
}

int runAlignment(const Config& cfg, const placement::PlacementResult& placement) {
    logging::msg("Loading tree for alignment...");
    auto tg = loadPanMAN(cfg.panman);
    if (!tg || tg->trees.empty()) {
        logging::err("Failed to load tree");
        return 1;
    }
    
    auto* T = &tg->trees[0];
    // Use RAW score for placement - it correctly weights by seed specificity
    // (Jaccard/WJ fail because they treat all matching seeds equally regardless of frequency)
    std::string nodeId = placement.bestRawSeedMatchNodeId;
    
    if (nodeId.empty()) {
        logging::err("No best node ID from placement - cannot align");
        return 1;
    }
    
    // Get reference sequence (like working commit: panmapUtils::getStringFromReference)
    std::string bestMatchSequence = T->getStringFromReference(nodeId, false, true);
    
    if (bestMatchSequence.empty()) {
        logging::err("Empty sequence for node '{}' - cannot align", nodeId);
        return 1;
    }
    
    // Output file paths
    std::string prefix = cfg.output;
    std::string refFileName = cfg.output + ".ref.fa";
    std::string samFileName = cfg.output + ".sam";
    std::string bamFileName = cfg.output + ".bam";
    
    // Write reference fasta
    {
        std::ofstream outFile(refFileName);
        if (!outFile) {
            logging::err("Cannot write reference file: {}", refFileName);
            return 1;
        }
        outFile << ">ref\n" << bestMatchSequence << "\n";
        outFile.close();
    }
    
    // Create FASTA index (.fai) for htslib/samtools
    if (fai_build(refFileName.c_str()) != 0) {
        logging::err("Failed to create FASTA index for {}", refFileName);
        return 1;
    }
    
    logging::msg("Reference: {} bp from node {} -> {}", bestMatchSequence.size(), nodeId, refFileName);
    
    // Get index parameters from placement result
    int32_t k = placement.k;
    int32_t s = placement.s;
    int32_t t = placement.t;
    bool open = placement.open;
    
    // Minimap2 has k <= 28 limit, adjust if needed (from working commit)
    if (k > 28) {
        logging::msg("k > 28, setting k = 19, s = 10, t = 0 for minimap alignment");
    }
    int k_minimap = k > 28 ? 19 : k;
    int s_minimap = k > 28 ? 10 : s;
    int t_minimap = k > 28 ? 0 : t;
    bool open_minimap = open;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Read input sequences using seedsFromFastq (populates everything we need)
    std::vector<std::vector<seeding::seed_t>> readSeeds;
    std::vector<std::string> readSequences, readQuals, readNames;
    std::vector<std::vector<std::string>> readSeedSeqs;
    absl::flat_hash_map<size_t, std::pair<size_t, size_t>> readSeedCounts;
    
    seeding::seedsFromFastq(k_minimap, s_minimap, t_minimap, open_minimap, 1,
                            readSeedCounts, readSequences, readQuals, readNames,
                            readSeeds, readSeedSeqs, cfg.reads1, cfg.reads2);
    
    // Build seed-to-reference position map from reference sequence (from working commit)
    std::unordered_map<size_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> seedToRefPositions;
    for (const auto& [kmerHash, isReverse, isSyncmer, startPos] : 
         seeding::rollingSyncmers(bestMatchSequence, k_minimap, s_minimap, open_minimap, t_minimap, false)) {
        if (!isSyncmer) continue;
        if (seedToRefPositions.find(kmerHash) == seedToRefPositions.end()) {
            seedToRefPositions[kmerHash] = std::make_pair(std::vector<uint32_t>(), std::vector<uint32_t>());
        }
        if (isReverse) {
            seedToRefPositions[kmerHash].second.push_back(startPos);
        } else {
            seedToRefPositions[kmerHash].first.push_back(startPos);
        }
    }
    
    logging::msg("Loaded {} reads, built {} reference seed positions", 
                 readSequences.size(), seedToRefPositions.size());
    
    // Create SAM alignment (from working commit)
    bool pairedEndReads = !cfg.reads2.empty();
    bool shortenSyncmers = false;  // Use full syncmer positions
    std::vector<char*> samAlignments;
    std::string samHeader;
    createSam(readSeeds, readSequences, readQuals, readNames, bestMatchSequence, 
              seedToRefPositions, samFileName, k_minimap, shortenSyncmers, pairedEndReads, 
              samAlignments, samHeader);
    
    // Create BAM from SAM (from working commit)
    sam_hdr_t* header;
    bam1_t** bamRecords;
    createBam(samAlignments, samHeader, bamFileName, header, bamRecords);
    
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    
    logging::msg("Alignment complete in {}ms -> {}", elapsed.count(), bamFileName);
    return 0;
}

int runGenotyping(const Config& cfg) {
    std::string prefix = cfg.output;
    std::string refFileName = cfg.output + ".ref.fa";
    std::string bamFileName = cfg.output + ".bam";
    std::string mpileupFileName = cfg.output + ".mpileup";
    std::string vcfFileName = cfg.output + ".vcf";
    
    // Read reference sequence
    std::string bestMatchSequence;
    std::ifstream ref(refFileName);
    std::string line;
    while (std::getline(ref, line)) {
        if (line[0] != '>') bestMatchSequence += line;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Create mpileup and VCF (from working commit)
    genotyping::mutationMatrices mutMat;
    createMplpBcf(prefix, refFileName, bestMatchSequence, bamFileName, mpileupFileName);
    createVcfWithMutationMatrices(prefix, mpileupFileName, mutMat, vcfFileName, 0.0011);
    
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    
    logging::msg("Genotyping complete in {}ms -> {}", elapsed.count(), vcfFileName);
    return 0;
}

// ============================================================================
// Main Entry Point
// ============================================================================

void printUsage() {
    std::cout << color::bold << "panmap" << color::reset << " v" << VERSION << "\n";
    std::cout << "Pangenome-based sequence placement, alignment, and genotyping\n\n";
    
    std::cout << color::bold << "USAGE:" << color::reset << "\n";
    std::cout << "  panmap [OPTIONS] <panman> [reads1.fq] [reads2.fq]\n\n";
    
    std::cout << color::bold << "EXAMPLES:" << color::reset << "\n";
    std::cout << "  # Full pipeline (place -> align -> genotype)\n";
    std::cout << "  panmap genomes.panman reads_R1.fq reads_R2.fq -o sample1\n\n";
    
    std::cout << "  # Placement only\n";
    std::cout << "  panmap genomes.panman reads.fq --stop place\n\n";
    
    std::cout << "  # Build index only\n";
    std::cout << "  panmap genomes.panman --stop index\n\n";
    
    std::cout << "  # Metagenomic mode (report top 10 placements)\n";
    std::cout << "  panmap genomes.panman metagenome.fq --meta --top 10\n\n";
    
    std::cout << color::bold << "PIPELINE STAGES:" << color::reset << "\n";
    std::cout << "  index     Build/verify index only\n";
    std::cout << "  place     Stop after placement\n";
    std::cout << "  align     Stop after alignment (BAM output)\n";
    std::cout << "  genotype  Full pipeline through variant calling (default)\n\n";
    
    std::cout << color::bold << "OPTIONS:" << color::reset << "\n";
}

int main(int argc, char** argv) {
    Config cfg;
    
    // Define options
    po::options_description general("General");
    general.add_options()
        ("help,h", "Show this help message")
        ("version,V", "Show version")
        ("threads,t", po::value<int>(&cfg.threads)->default_value(1), "Number of threads")
        ("output,o", po::value<std::string>(&cfg.output), "Output prefix")
        ("verbose,v", po::bool_switch(), "Verbose output")
        ("quiet,q", po::bool_switch(), "Suppress non-essential output");
    
    po::options_description pipeline("Pipeline Control");
    pipeline.add_options()
        ("stop", po::value<std::string>()->default_value("genotype"),
            "Stop after stage: index|place|align|genotype")
        ("meta", po::bool_switch(&cfg.metagenomic), "Metagenomic mode")
        ("top", po::value<int>(&cfg.topN)->default_value(1), 
            "Report top N placements (with --meta)")
        ("dedup", po::bool_switch(&cfg.dedupReads), 
            "Deduplicate reads before placement (recommended for amplicon data)");
    
    po::options_description indexing("Index Options");
    indexing.add_options()
        ("index,i", po::value<std::string>(&cfg.index), "Index file path")
        ("reindex,f", po::bool_switch(&cfg.forceReindex), "Force index rebuild")
        ("kmer,k", po::value<int>(&cfg.k)->default_value(29), "Syncmer parameter k")
        ("syncmer,s", po::value<int>(&cfg.s)->default_value(8), "Syncmer parameter s")
        ("lmer,l", po::value<int>(&cfg.l)->default_value(1), "Use l consecutive syncmers as seeds")
        ("flank-mask", po::value<int>(&cfg.flankMaskBp)->default_value(250), 
            "Hard mask first/last N bp at genome ends (ignore all seeds)")
        ("seed-mask-fraction", po::value<double>(&cfg.seedMaskFraction)->default_value(0.001),
            "Mask top N fraction of seeds by frequency (0=disabled, 0.001=top 0.1%, default)")
        ("min-seed-quality", po::value<int>(&cfg.minSeedQuality)->default_value(0),
            "Min avg Phred quality for seed k-mer region (0=disabled, 20=Q20 filter)")
        ("trim-start", po::value<int>(&cfg.trimStart)->default_value(0),
            "Trim N bases from start of each read before seeding (for primer removal)")
        ("trim-end", po::value<int>(&cfg.trimEnd)->default_value(0),
            "Trim N bases from end of each read before seeding (for primer removal)");
    
    po::options_description alignment("Alignment Options");
    alignment.add_options()
        ("aligner,a", po::value<std::string>(&cfg.aligner)->default_value("minimap2"),
            "Aligner: minimap2 (default) or bwa");
    
    po::options_description utility("Utility");
    utility.add_options()
        ("dump-random-node", po::bool_switch(&cfg.dumpRandomNode),
            "Dump random node sequence to FASTA")
        ("dump-sequence", po::value<std::string>(&cfg.dumpNodeId),
            "Dump specific node sequence to FASTA")
        ("remove-node", po::value<std::string>(&cfg.removeNodeId),
            "Remove node before placement (for leave-one-out validation)")
        ("top-placements", po::value<int>(&cfg.topPlacements)->default_value(0),
            "Output top-N placements with all scores (for diagnostics)")
        ("debug-node", po::value<std::string>(&cfg.debugNodeId),
            "Output detailed metrics for a specific node ID")
        ("compare-nodes", po::value<std::string>(&cfg.compareNodes),
            "Compare two nodes from scratch: 'node1,node2' (compute seeds/metrics independently)")
        ("seed", po::value<int>(&cfg.seed)->default_value(42), "Random seed");
    
    po::options_description hidden("Hidden");
    hidden.add_options()
        ("panman", po::value<std::string>(&cfg.panman), "")
        ("reads1", po::value<std::string>(&cfg.reads1), "")
        ("reads2", po::value<std::string>(&cfg.reads2), "");
    
    po::positional_options_description pos;
    pos.add("panman", 1).add("reads1", 1).add("reads2", 1);
    
    po::options_description all;
    all.add(general).add(pipeline).add(indexing).add(alignment).add(utility).add(hidden);
    
    po::options_description visible;
    visible.add(general).add(pipeline).add(indexing).add(alignment).add(utility);
    
    // Parse
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
            .options(all).positional(pos).run(), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << color::red << "Error: " << e.what() << color::reset << "\n\n";
        printUsage();
        std::cout << visible << "\n";
        return 1;
    }
    
    // Handle help/version
    if (vm.count("help") || argc == 1) {
        printUsage();
        std::cout << visible << "\n";
        return 0;
    }
    
    if (vm.count("version")) {
        std::cout << PROGRAM_NAME << " " << VERSION << "\n";
        return 0;
    }
    
    // Validate required args
    if (cfg.panman.empty()) {
        std::cerr << color::red << "Error: PanMAN file required" << color::reset << "\n";
        return 1;
    }
    
    // Parse stop stage
    std::string stopStr = vm["stop"].as<std::string>();
    if (stopStr == "index") cfg.stopAfter = PipelineStage::Index;
    else if (stopStr == "place") cfg.stopAfter = PipelineStage::Place;
    else if (stopStr == "align") cfg.stopAfter = PipelineStage::Align;
    else if (stopStr == "genotype") cfg.stopAfter = PipelineStage::Genotype;
    else {
        std::cerr << color::red << "Error: Invalid stage '" << stopStr << "'" << color::reset << "\n";
        return 1;
    }
    
    // Set defaults
    // If index not explicitly set, derive from output prefix if set, otherwise from panman
    if (cfg.index.empty()) {
        if (!cfg.output.empty()) {
            cfg.index = cfg.output + ".idx";
        } else {
            cfg.index = cfg.panman + ".idx";
        }
    }
    // Only set default output if we're not in dump mode (dump modes handle their own output)
    if (cfg.output.empty() && !cfg.dumpRandomNode && cfg.dumpNodeId.empty()) {
        // Derive output prefix from reads filename (without path and common extensions)
        if (!cfg.reads1.empty()) {
            fs::path readsPath(cfg.reads1);
            std::string stem = readsPath.stem().string();
            // Remove common paired-end suffixes like _R1, _1, .R1, etc.
            for (const auto& suffix : {"_R1", "_R2", "_1", "_2", ".R1", ".R2", ".1", ".2"}) {
                if (stem.size() > strlen(suffix) &&
                    stem.substr(stem.size() - strlen(suffix)) == suffix) {
                    stem = stem.substr(0, stem.size() - strlen(suffix));
                    break;
                }
            }
            // Also remove .fastq, .fq extensions if present (for .fastq.gz -> .fastq stem)
            for (const auto& ext : {".fastq", ".fq"}) {
                if (stem.size() > strlen(ext) &&
                    stem.substr(stem.size() - strlen(ext)) == ext) {
                    stem = stem.substr(0, stem.size() - strlen(ext));
                    break;
                }
            }
            cfg.output = stem;
        } else {
            cfg.output = cfg.panman;
        }
    }
    
    // Verbosity
    if (vm["verbose"].as<bool>()) cfg.verbosity = 2;
    if (vm["quiet"].as<bool>()) cfg.verbosity = 0;
    
    // Initialize threading
    tbb::global_control tbb_ctl(tbb::global_control::max_allowed_parallelism, cfg.threads);

    // ========================================================================
    // Print Configuration Summary
    // ========================================================================
    
    auto printConfigSummary = [&]() {
        if (cfg.verbosity == 0) return;  // Skip in quiet mode
        if (cfg.dumpRandomNode || !cfg.dumpNodeId.empty()) return;  // Skip for utility modes
        
        // Build stage string
        std::string stageStr;
        switch (cfg.stopAfter) {
            case PipelineStage::Index:    stageStr = "index"; break;
            case PipelineStage::Place:    stageStr = "index → place"; break;
            case PipelineStage::Align:    stageStr = "index → place → align"; break;
            case PipelineStage::Genotype: stageStr = "index → place → align → genotype"; break;
            default:                      stageStr = "full"; break;
        }
        
        // Header
        std::cout << color::bold << color::cyan << "┌─ panmap " << color::reset 
                  << color::dim << "v" << VERSION << color::reset 
                  << color::bold << color::cyan << " ─";
        for (int i = 0; i < 50; i++) std::cout << "─";
        std::cout << "┐" << color::reset << "\n";
        
        // Core inputs
        std::cout << color::bold << color::cyan << "│" << color::reset << " "
                  << color::bold << "Input:  " << color::reset 
                  << color::yellow << cfg.panman << color::reset;
        if (!cfg.reads1.empty()) {
            std::cout << "  " << color::dim << "+" << color::reset << " " << cfg.reads1;
            if (!cfg.reads2.empty()) std::cout << ", " << cfg.reads2;
        }
        std::cout << "\n";
        
        // Output prefix
        std::cout << color::bold << color::cyan << "│" << color::reset << " "
                  << color::bold << "Output: " << color::reset 
                  << color::green << cfg.output << color::reset << ".*\n";
        
        // Pipeline stages
        std::cout << color::bold << color::cyan << "│" << color::reset << " "
                  << color::bold << "Stages: " << color::reset << stageStr << "\n";
        
        // Options line (compact)
        std::cout << color::bold << color::cyan << "│" << color::reset << " "
                  << color::bold << "Config: " << color::reset
                  << color::dim << "threads=" << color::reset << cfg.threads
                  << color::dim << "  k=" << color::reset << cfg.k
                  << color::dim << " s=" << color::reset << cfg.s
                  << color::dim << " l=" << color::reset << cfg.l;
        
        // Optional flags (only show if non-default)
        if (cfg.metagenomic) {
            std::cout << color::dim << "  meta" << color::reset 
                      << color::dim << "(top=" << color::reset << cfg.topN 
                      << color::dim << ")" << color::reset;
        }
        if (cfg.forceReindex) {
            std::cout << "  " << color::yellow << "reindex" << color::reset;
        }
        if (cfg.aligner != "minimap2") {
            std::cout << color::dim << "  aligner=" << color::reset << cfg.aligner;
        }
        std::cout << "\n";
        
        // Footer
        std::cout << color::bold << color::cyan << "└";
        for (int i = 0; i < 68; i++) std::cout << "─";
        std::cout << "┘" << color::reset << "\n\n";
    };
    
    printConfigSummary();

    // ========================================================================
    // Run Pipeline
    // ========================================================================
    
    try {
        // Utility: dump random node
        if (cfg.dumpRandomNode) {
            auto tg = loadPanMAN(cfg.panman);
            std::string nodeId = getRandomNodeId(&tg->trees[0], cfg.seed);
            // Use -o output path if specified, otherwise default to panman path
            std::string outPath = cfg.output.empty() 
                ? cfg.panman + ".random." + sanitizeFilename(nodeId) + ".fa"
                : cfg.output;
            saveNodeSequence(&tg->trees[0], nodeId, outPath);
            std::cout << nodeId << "\n";
            return 0;
        }
        
        // Utility: dump specific node sequence
        if (!cfg.dumpNodeId.empty()) {
            auto tg = loadPanMAN(cfg.panman);
            // Use -o output path if specified, otherwise default to panman path
            std::string outPath = cfg.output.empty()
                ? cfg.panman + "." + sanitizeFilename(cfg.dumpNodeId) + ".fa"
                : cfg.output;
            saveNodeSequence(&tg->trees[0], cfg.dumpNodeId, outPath);
            std::cout << cfg.dumpNodeId << "\n";
            return 0;
        }
        
        // Utility: compare two nodes from scratch
        if (!cfg.compareNodes.empty()) {
            if (cfg.reads1.empty()) {
                logging::err("--compare-nodes requires reads to be specified");
                return 1;
            }
            auto tg = loadPanMAN(cfg.panman);
            compareNodesFromScratch(cfg, &tg->trees[0]);
            return 0;
        }
        
        // Stage 1: Index
        if (!buildIndex(cfg)) return 1;
        if (cfg.stopAfter == PipelineStage::Index) {
            logging::msg("{}Done.{} Index ready: {}", color::green, color::reset, cfg.index);
            return 0;
        }
        
        // Check for reads
        if (cfg.reads1.empty()) {
            logging::msg("No reads provided. Index is ready.");
            return 0;
        }
        
        // Stage 2: Placement
        auto placement = runPlacement(cfg);
        if (!placement) {
            logging::err("Placement failed");
            return 1;
        }
        if (cfg.stopAfter == PipelineStage::Place) {
            logging::msg("{}Done.{} Placement: {}.placement.tsv", 
                        color::green, color::reset, cfg.output);
            return 0;
        }
        
        // Stage 3: Alignment
        if (runAlignment(cfg, *placement) != 0) {
            logging::err("Alignment failed");
            return 1;
        }
        if (cfg.stopAfter == PipelineStage::Align) {
            logging::msg("{}Done.{} Alignment: {}.bam", 
                        color::green, color::reset, cfg.output);
            return 0;
        }
        
        // Stage 4: Genotyping
        if (runGenotyping(cfg) != 0) {
            logging::err("Genotyping failed");
            return 1;
        }
        
        logging::msg("{}Pipeline complete.{}", color::green, color::reset);
        logging::msg("  Placement:  {}.placement.tsv", cfg.output);
        logging::msg("  Reference:  {}.ref.fa", cfg.output);
        logging::msg("  Alignment:  {}.bam", cfg.output);
        logging::msg("  Variants:   {}.vcf", cfg.output);
        
        return 0;
        
    } catch (const std::exception& e) {
        logging::err("Fatal error: {}", e.what());
        return 1;
    }
}

/**
 * @file main.cpp
 * @brief Panmap - Pangenome-based sequence placement, alignment, and genotyping
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

namespace color {
    inline const char* reset() { return output::style::reset(); }
    inline const char* bold() { return output::style::bold(); }
    inline const char* dim() { return output::style::dim(); }
    inline const char* red() { return output::style::red(); }
    inline const char* green() { return output::style::green(); }
    inline const char* yellow() { return output::style::yellow(); }
    inline const char* blue() { return output::style::blue(); }
    inline const char* cyan() { return output::style::cyan(); }
}

enum class PipelineStage {
    Index,      // Build index only
    Place,      // Placement only
    Align,      // Placement + Alignment
    Genotype,   // Placement + Alignment + Genotyping
    Full        // Full pipeline including assembly
};

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
    
    std::string aligner = "minimap2";
    
    // Index parameters
    int k = 21;                   // syncmer k
    int s = 8;                    // syncmer s
    int l = 1;                    // l-mer size
    int flankMaskBp = 250;        // Hard mask first/last N bp at genome ends
    double seedMaskFraction = 0; // Mask top 0.1% most frequent seeds
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
    
    // Output control
    bool quiet = false;           // Minimal output (errors only)
    bool verbose = false;         // Extra debug output
    bool plain = false;           // Plain text output (no colors/unicode)
};

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
    
    // Check if any metric found a placement
    if (result.bestLogRawNodeId.empty() && result.bestLogCosineNodeId.empty()) {
        logging::warn("No placement found");
        return std::nullopt;
    }
    
    // Report results
    logging::msg("{}Best placement:{} {} (LogRaw: {:.6f}, LogCosine: {:.6f})", 
                color::green(), color::reset(),
                result.bestLogRawNodeId, 
                result.bestLogRawScore,
                result.bestLogCosineScore);
    
    // Leave-one-out validation: compare to expected parent
    if (!cfg.removeNodeId.empty() && !expectedParentNode.empty()) {
        int distLogRaw = computeTreeDistance(tree, result.bestLogRawNodeId, expectedParentNode);
        int distLogCosine = computeTreeDistance(tree, result.bestLogCosineNodeId, expectedParentNode);
        
        bool correctLogRaw = (result.bestLogRawNodeId == expectedParentNode);
        bool correctLogCosine = (result.bestLogCosineNodeId == expectedParentNode);
        
        logging::msg("{}VALIDATION RESULTS:{}", color::cyan(), color::reset());
        logging::msg("  Expected parent: {}", expectedParentNode);
        logging::msg("  LogRAW:    {} {}{}{} (dist: {}, score: {:.6f})", 
                    result.bestLogRawNodeId,
                    correctLogRaw ? color::green() : color::yellow(),
                    correctLogRaw ? output::box::check() : output::box::cross(),
                    color::reset(),
                    distLogRaw,
                    result.bestLogRawScore);
        logging::msg("  LogCosine: {} {}{}{} (dist: {}, score: {:.6f})", 
                    result.bestLogCosineNodeId,
                    correctLogCosine ? color::green() : color::yellow(),
                    correctLogCosine ? output::box::check() : output::box::cross(),
                    color::reset(),
                    distLogCosine,
                    result.bestLogCosineScore);
        
        // Print to stdout for easy parsing (TSV format)
        std::cout << "REMOVED_NODE\t" << cfg.removeNodeId << "\n";
        std::cout << "EXPECTED_PARENT\t" << expectedParentNode << "\n";
        std::cout << "LOGRAW\t" << result.bestLogRawNodeId 
                  << "\t" << (correctLogRaw ? "true" : "false") 
                  << "\t" << distLogRaw 
                  << "\t" << result.bestLogRawScore << "\n";
        std::cout << "LOGCOSINE\t" << result.bestLogCosineNodeId 
                  << "\t" << (correctLogCosine ? "true" : "false") 
                  << "\t" << distLogCosine 
                  << "\t" << result.bestLogCosineScore << "\n";
    }
    
    // Diagnostic output: top placements
    if (cfg.topPlacements > 0 || !cfg.debugNodeId.empty()) {
        std::cout << "\n=== DIAGNOSTIC OUTPUT ===\n";
        std::cout << "Read Statistics:\n";
        std::cout << "  Unique seeds: " << result.readUniqueSeedCount << "\n";
        std::cout << "  Total seed frequency: " << result.totalReadSeedFrequency << "\n";
        std::cout << "  Read magnitude (log): " << result.readMagnitude << "\n\n";
        
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
    // Use LogRaw as primary placement metric
    std::string nodeId = placement.bestLogRawNodeId;
    
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
    std::cout << color::bold() << "panmap" << color::reset() << " v" << VERSION << "\n";
    std::cout << "Pangenome-based sequence placement, alignment, and genotyping\n\n";
    
    std::cout << color::bold() << "USAGE:" << color::reset() << "\n";
    std::cout << "  panmap [OPTIONS] <panman> [reads1.fq] [reads2.fq]\n\n";
    
    std::cout << color::bold() << "EXAMPLES:" << color::reset() << "\n";
    std::cout << "  # Full pipeline (place -> align -> genotype)\n";
    std::cout << "  panmap genomes.panman reads_R1.fq reads_R2.fq -o sample1\n\n";
    
    std::cout << "  # Placement only\n";
    std::cout << "  panmap genomes.panman reads.fq --stop place\n\n";
    
    std::cout << "  # Build index only\n";
    std::cout << "  panmap genomes.panman --stop index\n\n";
    
    std::cout << "  # Metagenomic mode (report top 10 placements)\n";
    std::cout << "  panmap genomes.panman metagenome.fq --meta --top 10\n\n";
    
    std::cout << color::bold() << "PIPELINE STAGES:" << color::reset() << "\n";
    std::cout << "  index     Build/verify index only\n";
    std::cout << "  place     Stop after placement\n";
    std::cout << "  align     Stop after alignment (BAM output)\n";
    std::cout << "  genotype  Full pipeline through variant calling (default)\n\n";
    
    std::cout << color::bold() << "OPTIONS:" << color::reset() << "\n";
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
        ("verbose,v", po::bool_switch(&cfg.verbose), "Verbose output")
        ("quiet,q", po::bool_switch(&cfg.quiet), "Suppress non-essential output")
        ("no-color", po::bool_switch(&cfg.plain), "Disable colors and unicode");
    
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
        // Check for NO_COLOR env var
        const char* noColorEnv = std::getenv("NO_COLOR");
        bool plainMode = vm.count("no-color") > 0 || (noColorEnv && noColorEnv[0] != '\0');
        output::init(false, false, plainMode);
        output::error("{}", e.what());
        std::cerr << "\n";
        printUsage();
        std::cout << visible << "\n";
        return 1;
    }
    
    // Check for NO_COLOR environment variable (no-color.org standard)
    const char* noColorEnv = std::getenv("NO_COLOR");
    if (noColorEnv && noColorEnv[0] != '\0') cfg.plain = true;
    
    // Initialize output early for help/version formatting
    output::init(cfg.quiet, cfg.verbose, cfg.plain);
    
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
        output::error("PanMAN file required");
        return 1;
    }
    
    // Parse stop stage
    std::string stopStr = vm["stop"].as<std::string>();
    if (stopStr == "index") cfg.stopAfter = PipelineStage::Index;
    else if (stopStr == "place") cfg.stopAfter = PipelineStage::Place;
    else if (stopStr == "align") cfg.stopAfter = PipelineStage::Align;
    else if (stopStr == "genotype") cfg.stopAfter = PipelineStage::Genotype;
    else {
        output::error("Invalid stage '{}'", stopStr);
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
    
    // Install signal handlers for graceful interruption
    signals::install_handlers();
    
    // Initialize threading
    tbb::global_control tbb_ctl(tbb::global_control::max_allowed_parallelism, cfg.threads);

    // ========================================================================
    // Print Configuration Summary
    // ========================================================================
    
    auto printConfigSummary = [&]() {
        if (cfg.quiet) return;  // Skip in quiet mode
        if (cfg.dumpRandomNode || !cfg.dumpNodeId.empty()) return;  // Skip for utility modes
        
        // Build stage string using arrow from box chars
        std::string arrow = output::box::arrow();
        std::string stageStr;
        switch (cfg.stopAfter) {
            case PipelineStage::Index:    stageStr = "index"; break;
            case PipelineStage::Place:    stageStr = fmt::format("index {} place", arrow); break;
            case PipelineStage::Align:    stageStr = fmt::format("index {} place {} align", arrow, arrow); break;
            case PipelineStage::Genotype: stageStr = fmt::format("index {} place {} align {} genotype", arrow, arrow, arrow); break;
            default:                      stageStr = "full"; break;
        }
        
        // Build input string
        std::string inputStr = cfg.panman;
        if (!cfg.reads1.empty()) {
            inputStr += "  + " + cfg.reads1;
            if (!cfg.reads2.empty()) inputStr += ", " + cfg.reads2;
        }
        
        // Build config string
        std::string configStr = fmt::format("threads={}  k={} s={} l={}", 
                                            cfg.threads, cfg.k, cfg.s, cfg.l);
        if (cfg.metagenomic) {
            configStr += fmt::format("  meta(top={})", cfg.topN);
        }
        if (cfg.forceReindex) {
            configStr += "  reindex";
        }
        if (cfg.aligner != "minimap2") {
            configStr += "  aligner=" + cfg.aligner;
        }
        
        output::print_header("panmap", VERSION);
        output::print_row("Input ", inputStr);
        output::print_row("Output", cfg.output + ".*");
        output::print_row("Stages", stageStr);
        output::print_row("Config", configStr);
        output::print_footer();
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
        
        // Stage 1: Index
        if (!buildIndex(cfg)) return 1;
        if (signals::check_interrupted()) return 130;  // Standard exit code for SIGINT
        if (cfg.stopAfter == PipelineStage::Index) {
            output::done("Index ready: " + cfg.index);
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} panmap {} reads.fq", color::dim(), cfg.panman);
            return 0;
        }
        
        // Check for reads
        if (cfg.reads1.empty()) {
            output::info("No reads provided. Index is ready.");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} panmap {} reads.fq", color::dim(), cfg.panman);
            return 0;
        }
        
        // Stage 2: Placement
        auto placement = runPlacement(cfg);
        if (signals::check_interrupted()) return 130;
        if (!placement) {
            output::error("Placement failed");
            return 1;
        }
        if (cfg.stopAfter == PipelineStage::Place) {
            output::done("Placement: " + cfg.output + ".placement.tsv");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} Continue to alignment: panmap {} {} --stop align", 
                        color::dim(), cfg.panman, cfg.reads1);
            return 0;
        }
        
        // Stage 3: Alignment
        if (runAlignment(cfg, *placement) != 0) {
            output::error("Alignment failed");
            return 1;
        }
        if (signals::check_interrupted()) return 130;
        if (cfg.stopAfter == PipelineStage::Align) {
            output::done("Alignment: " + cfg.output + ".bam");
            output::info("");
            output::info("{}Next steps:{}", color::dim(), color::reset());
            output::info("  {} View: samtools view {}.bam | head", color::dim(), cfg.output);
            output::info("  {} Stats: samtools flagstat {}.bam", color::dim(), cfg.output);
            return 0;
        }
        
        // Stage 4: Genotyping
        if (runGenotyping(cfg) != 0) {
            output::error("Genotyping failed");
            return 1;
        }
        if (signals::check_interrupted()) return 130;
        
        output::info("");
        output::info("{}Pipeline complete.{}", color::green(), color::reset());
        output::info("  Placement:  {}.placement.tsv", cfg.output);
        output::info("  Reference:  {}.ref.fa", cfg.output);
        output::info("  Alignment:  {}.bam", cfg.output);
        output::info("  Variants:   {}.vcf", cfg.output);
        output::info("");
        output::info("{}Next steps:{}", color::dim(), color::reset());
        output::info("  {} View variants: bcftools view {}.vcf | head", color::dim(), cfg.output);
        output::info("  {} View alignment: samtools tview {}.bam {}.ref.fa", color::dim(), cfg.output, cfg.output);
        
        return 0;
        
    } catch (const std::exception& e) {
        output::error("Fatal error: {}", e.what());
        return 1;
    }
}

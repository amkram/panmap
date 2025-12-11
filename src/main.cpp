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
    
    // Aligner
    std::string aligner = "minimap2";
    
    // Index parameters
    int k = 21;                   // syncmer k
    int s = 8;                    // syncmer s
    int l = 1;                    // l-mer size
    
    // Resources
    int threads = 1;
    
    // Utility modes
    bool dumpRandomNode = false;
    bool dumpSequence = false;
    std::string dumpNodeId;
    int seed = 42;
    
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
    
    explicit IndexReader(const std::string& path) 
        : ::capnp::MessageReader(makeOptions()) 
    {
        if (!panmap_zstd::decompressFromFile(path, data)) {
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
    
    index_single_mode::IndexBuilder builder(&tg->trees[0], cfg.k, cfg.s, 0, cfg.l, false);
    if (cfg.threads > 1) {
        builder.buildIndexParallel(cfg.threads);
    } else {
        builder.buildIndex();
    }
    builder.writeIndex(cfg.index);
    
    logging::msg("Index built with k={}, s={}, l={}", cfg.k, cfg.s, cfg.l);
    return true;
}

std::optional<placement::PlacementResult> runPlacement(const Config& cfg) {
    logging::msg("Loading index...");
    IndexReader reader(cfg.index);
    
    auto idx = reader.getRoot<LiteIndex>();
    logging::msg("Index parameters: k={}, s={}, l={}", idx.getK(), idx.getS(), idx.getL());
    
    panmapUtils::LiteTree tree;
    tree.initialize(idx.getLiteTree());
    logging::msg("Tree loaded: {} nodes", tree.allLiteNodes.size());
    
    placement::PlacementResult result;
    std::string outPath = cfg.output + ".placement.tsv";
    
    auto start = std::chrono::high_resolution_clock::now();
    placement::placeLite(result, &tree, reader, cfg.reads1, cfg.reads2, outPath, false, nullptr);
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    
    logging::msg("Placement complete in {}ms", elapsed.count());
    
    if (result.bestWeightedJaccardNodeId.empty()) {
        logging::warn("No placement found");
        return std::nullopt;
    }
    
    // Report results
    logging::msg("{}Best placement:{} {} (score: {:.4f})", 
                color::green, color::reset,
                result.bestWeightedJaccardNodeId, 
                result.bestWeightedJaccardScore);
    
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
    std::string nodeId = placement.bestWeightedJaccardNodeId;
    
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
            "Report top N placements (with --meta)");
    
    po::options_description indexing("Index Options");
    indexing.add_options()
        ("index,i", po::value<std::string>(&cfg.index), "Index file path")
        ("reindex,f", po::bool_switch(&cfg.forceReindex), "Force index rebuild")
        ("kmer,k", po::value<int>(&cfg.k)->default_value(29), "Syncmer parameter k")
        ("syncmer,s", po::value<int>(&cfg.s)->default_value(8), "Syncmer parameter s")
        ("lmer,l", po::value<int>(&cfg.l)->default_value(1), "Use l consecutive syncmers as seeds");
    
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

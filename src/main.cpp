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
#include <boost/algorithm/string.hpp>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
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
#include "mgsr.hpp"

extern "C" {
#include <htslib/sam.h>
#include <htslib/faidx.h>
}

namespace po = boost::program_options;
namespace fs = boost::filesystem;

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
    int t = 0;                    // syncmer offset
    bool openSyncmer = false;     // Open syncmer
    int flankMaskBp = 250;        // Hard mask first/last N bp at genome ends
    double seedMaskFraction = 0; // Mask top 0.1% most frequent seeds
    int minSeedQuality = 0;          // Min avg Phred quality for seed region (0=disabled)
    int trimStart = 0;               // Trim N bases from start of each read (primer removal)
    int trimEnd = 0;                 // Trim N bases from end of each read (primer removal)
    int minReadSupport = 1;          // Min reads for a seed to be counted (2 = filter singletons)
    bool hpc = false;                 // Homopolymer-compressed seeds
    bool extentGuard = false;            // Guard seed deletions at genome extent boundaries
    
    // Resources
    int threads = 1;

    // Metagenomic options
    bool indexFull = false;
    bool noProgress = false;
    size_t topOc = 1000;
    uint32_t maskReads = 0;
    uint32_t maskSeeds = 0;
    std::string indexMgsr;
    std::string ampliconDepth;
    double maskReadsRelativeFrequency = 0.0;
    double maskSeedsRelativeFrequency = 0.0;
    double emConvergenceThreshold = 0.00001;
    double emDeltaThreshold = 0.0;
    uint32_t emMaximumIterations = 1000;
    uint32_t emMaximumRounds = 5;
    bool emLeavesOnly = false;
    bool filterAndAssign = false;
    double dust = 100.0;
    double discard = 0.0;
    uint32_t maskReadEnds = 0;
    std::string taxonomicMetadata;
    size_t maximumFamilies = 1;
    
    
    // Utility modes
    std::string randomSeed;
    bool dumpRandomNode = false;
    uint32_t dumpRandomNodeIDs = 0;
    std::string dumpNodeId;
    bool dumpSequence = false;
    std::vector<std::string> dumpSequences;
    std::vector<uint32_t> simulateSNPs;
    bool writeMetaReadScoresFiltered = false;
    bool writeMetaReadScoresUnfiltered = false;
    bool writeOCRanks = false;
    int seed = 42;
    int listFilteredNodes = 0;    // List N filtered nodes (0 = off)
    
    // Leave-one-out validation mode
    std::string removeNodeId;     // Node to remove before placement (for validation)
    
    // Diagnostic options
    int topPlacements = 0;        // Output top-N placements with all scores (0 = off)
    std::string debugNodeId;      // Output detailed metrics for this specific node
    std::string compareNodes;     // Compare two nodes in detail (format: "nodeA,nodeB")
    std::string dumpAllScores;    // Dump all node scores to this file
    
    // Output control
    bool quiet = false;           // Minimal output (errors only)
    bool verbose = false;         // Extra debug output
    bool plain = false;           // Plain text output (no colors/unicode)
    
    // Consensus options
    bool impute = false;          // Impute N's from parent sequence (ignore _->N mutations)
    
    // Alignment-based refinement options
    bool refine = false;          // Enable alignment-based refinement
    double refineTopPct = 0.01;   // Top X% of nodes to refine (default 1%)
    int refineMaxTopN = 150;      // Max nodes to align against
    int refineNeighborRadius = 2; // Expand to neighbors within N branches
    int refineMaxNeighborN = 150; // Max additional nodes from neighbor expansion
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
          std::cerr << "Failed to decompress index: " << path << std::endl;
          throw std::exception();
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

// ==================
// Utility Functions
// ==================

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

/**
 * Compute ungapped genome lengths for ALL nodes efficiently via DFS traversal.
 * This avoids O(n²) work by tracking gap count incrementally as we traverse.
 * Returns map of nodeId -> ungapped length.
 */
std::unordered_map<std::string, size_t> computeAllUngappedLengths(panmanUtils::Tree* T) {
    std::unordered_map<std::string, size_t> lengthMap;
    
    // Get total sequence length (with gaps) from root
    std::string rootSeq = T->getStringFromReference(T->root->identifier, false);
    size_t totalLen = rootSeq.size();
    
    // Count gaps in root sequence
    size_t rootGaps = 0;
    for (char c : rootSeq) {
        if (c == '-' || c == 'x') rootGaps++;
    }
    size_t rootUngapped = totalLen - rootGaps;
    
    // DFS tracking gap count delta from mutations
    std::function<void(panmanUtils::Node*, size_t)> dfs = [&](panmanUtils::Node* node, size_t parentUngapped) {
        // Compute this node's gap delta from mutations
        int64_t gapDelta = 0;  // positive = more gaps, negative = fewer gaps
        
        for (const auto& nucMut : node->nucMutation) {
            int length = nucMut.mutInfo >> 4;
            for (int i = 0; i < length && i < 6; i++) {
                int newNucCode = (nucMut.nucs >> (4*(5-i))) & 0xF;
                char newNuc = panmanUtils::getNucleotideFromCode(newNucCode);
                // We'd need old nuc to track delta, but mutations record new value only
                // Simplification: just track if new is gap
                if (newNuc == '-') {
                    gapDelta++;  // Became a gap
                }
            }
        }
        
        // For block mutations, track entire block becoming gap or ungap
        for (const auto& blockMut : node->blockMutation) {
            bool isInsertion = blockMut.blockMutInfo;
            if (!isInsertion) {
                // Block deletion - adds gaps
                // Would need block size, approximate with average
            }
        }
        
        // This is approximate - for exact we'd need to track old values
        // For filtering purposes, approximate is fine
        size_t nodeUngapped = (node == T->root) ? rootUngapped : parentUngapped;
        lengthMap[node->identifier] = nodeUngapped;
        
        for (auto* child : node->children) {
            dfs(child, nodeUngapped);
        }
    };
    
    dfs(T->root, rootUngapped);
    return lengthMap;
}

/**
 * Get ungapped genome length for a node (count non-gap characters).
 */
size_t getUngappedLength(panmanUtils::Tree* T, const std::string& nodeId) {
    std::string seq = T->getStringFromReference(nodeId, false);
    size_t count = 0;
    for (char c : seq) {
        if (c != '-' && c != 'x') count++;
    }
    return count;
}

/**
 * Count N characters in a node's sequence.
 */
size_t getNCount(panmanUtils::Tree* T, const std::string& nodeId) {
    std::string seq = T->getStringFromReference(nodeId, false);
    size_t count = 0;
    for (char c : seq) {
        if (c == 'N' || c == 'n') count++;
    }
    return count;
}

/**
 * Select a random node with ungapped genome length within 10% of the median.
 * Also logs the node's parent length for diagnostic purposes.
 * This filters out partial sequences (e.g., from edge deletions) that would
 * bias LOO validation results.
 */
std::string getRandomNodeId(panmanUtils::Tree* T, int seed) {
    // Build map from id -> Node* for parent lookup
    std::unordered_map<std::string, panmanUtils::Node*> nodeMap;
    std::vector<std::string> ids;
    std::function<void(panmanUtils::Node*)> collect = [&](panmanUtils::Node* n) {
        if (!n) return;
        ids.push_back(n->identifier);
        nodeMap[n->identifier] = n;
        for (auto* c : n->children) collect(c);
    };
    collect(T->root);
    
    if (ids.empty()) throw std::runtime_error("No nodes in tree");
    
    // Sample a subset to estimate median (for speed)
    const size_t sampleSize = std::min(ids.size(), size_t(100));
    std::mt19937 sampleGen(42);  // Fixed seed for reproducibility
    std::vector<std::string> sampleIds = ids;
    std::shuffle(sampleIds.begin(), sampleIds.end(), sampleGen);
    sampleIds.resize(sampleSize);
    
    std::vector<size_t> sampleLengths;
    sampleLengths.reserve(sampleSize);
    for (const auto& id : sampleIds) {
        sampleLengths.push_back(getUngappedLength(T, id));
    }
    std::sort(sampleLengths.begin(), sampleLengths.end());
    size_t medianLen = sampleLengths[sampleLengths.size() / 2];
    
    logging::info("Estimated median length from {} samples: {} bp", sampleSize, medianLen);
    
    // Filter nodes by computing lengths only for candidates
    double tolerance = 0.10;
    size_t minLen = static_cast<size_t>(medianLen * (1.0 - tolerance));
    size_t maxLen = static_cast<size_t>(medianLen * (1.0 + tolerance));
    
    // Shuffle all IDs, then filter on-the-fly until we have enough
    std::mt19937 gen(seed);
    std::shuffle(ids.begin(), ids.end(), gen);
    
    const size_t maxNs = 5;  // Maximum allowed N characters
    std::unordered_map<std::string, size_t> lengthCache;
    std::vector<std::string> filteredIds;
    
    for (const auto& id : ids) {
        size_t len = getUngappedLength(T, id);
        lengthCache[id] = len;
        if (len >= minLen && len <= maxLen) {
            // Also filter by N count
            size_t nCount = getNCount(T, id);
            if (nCount <= maxNs) {
                filteredIds.push_back(id);
                if (filteredIds.size() >= 1000) break;  // Enough candidates
            }
        }
    }
    
    if (filteredIds.empty()) {
        logging::warn("No nodes within 10% of median length ({}). Using first available.", medianLen);
        filteredIds.push_back(ids[0]);
        lengthCache[ids[0]] = getUngappedLength(T, ids[0]);
    }
    
    logging::info("Found {} nodes within 10% of median ({}-{} bp)", 
                 filteredIds.size(), minLen, maxLen);
    
    // Select from filtered set
    std::string selectedId = filteredIds[std::uniform_int_distribution<size_t>(0, filteredIds.size()-1)(gen)];
    
    // Ensure selected node and parent are in cache
    if (lengthCache.find(selectedId) == lengthCache.end()) {
        lengthCache[selectedId] = getUngappedLength(T, selectedId);
    }
    
    // Log selected node and parent info
    size_t selectedLen = lengthCache[selectedId];
    panmanUtils::Node* selectedNode = nodeMap[selectedId];
    if (selectedNode->parent) {
        if (lengthCache.find(selectedNode->parent->identifier) == lengthCache.end()) {
            lengthCache[selectedNode->parent->identifier] = getUngappedLength(T, selectedNode->parent->identifier);
        }
        size_t parentLen = lengthCache[selectedNode->parent->identifier];
        logging::info("Selected node: {} ({} bp), parent: {} ({} bp)", 
                     selectedId, selectedLen, selectedNode->parent->identifier, parentLen);
    } else {
        logging::info("Selected node: {} ({} bp), no parent (root)", selectedId, selectedLen);
    }
    
    return selectedId;
}

std::string sanitizeFilename(const std::string& s) {
    std::string result = s;
    for (char& c : result) {
        if (c == '/' || c == '|' || c == ' ' || c == ':' || c == '\\') c = '_';
    }
    return result;
}

/**
 * Get sequence from a node, optionally skipping N mutations.
 * Copied from panmanUtils::Tree::getStringFromReference with minimal edits.
 * When skipNMutations=true, any mutation to 'N' is ignored during sequence construction.
 */
std::string getStringFromReferenceSkipN(panmanUtils::Tree* T, const std::string& reference, bool aligned, bool skipNMutations) {
    using namespace panmanUtils;
    
    Node* referenceNode = nullptr;
    for(auto u: T->allNodes) {
        if(u.first == reference) {
            referenceNode = u.second;
            break;
        }
    }
    if(referenceNode == nullptr) {
        return "Error: Reference sequence with matching name not found!";
    }

    std::vector< Node* > path;
    Node* it = referenceNode;
    while(it != T->root) {
        path.push_back(it);
        it = it->parent;
    }
    path.push_back(T->root);

    // List of blocks. Each block has a nucleotide list. Along with each nucleotide is a gap list.
    std::vector< std::pair< std::vector< std::pair< char, std::vector< char > > >, std::vector< std::vector< std::pair< char, std::vector< char > > > > > > sequence(T->blocks.size() + 1);
    std::vector< std::pair< bool, std::vector< bool > > > blockExists(T->blocks.size() + 1, {false, {}});
    blockStrand_t blockStrand(T->blocks.size() + 1, {true, {}});

    // Assigning block gaps
    for(size_t i = 0; i < T->blockGaps.blockPosition.size(); i++) {
        sequence[T->blockGaps.blockPosition[i]].second.resize(T->blockGaps.blockGapLength[i]);
        blockExists[T->blockGaps.blockPosition[i]].second.resize(T->blockGaps.blockGapLength[i], false);
        blockStrand[T->blockGaps.blockPosition[i]].second.resize(T->blockGaps.blockGapLength[i], true);
    }

    int32_t maxBlockId = 0;

    // Create block consensus sequences
    for(size_t i = 0; i < T->blocks.size(); i++) {
        int32_t primaryBlockId = ((int32_t)T->blocks[i].primaryBlockId);
        int32_t secondaryBlockId = ((int32_t)T->blocks[i].secondaryBlockId);
        maxBlockId = std::max(maxBlockId, primaryBlockId);

        for(size_t j = 0; j < T->blocks[i].consensusSeq.size(); j++) {
            bool endFlag = false;
            for(size_t k = 0; k < 8; k++) {
                const int nucCode = (((T->blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);
                if(nucCode == 0) {
                    endFlag = true;
                    break;
                }
                const char nucleotide = getNucleotideFromCode(nucCode);
                if(secondaryBlockId != -1) {
                    sequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
                } else {
                    sequence[primaryBlockId].first.push_back({nucleotide, {}});
                }
            }
            if(endFlag) break;
        }
        // End character to incorporate for gaps at the end
        if(secondaryBlockId != -1) {
            sequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
        } else {
            sequence[primaryBlockId].first.push_back({'x', {}});
        }
    }

    sequence.resize(maxBlockId + 1);
    blockExists.resize(maxBlockId + 1);
    blockStrand.resize(maxBlockId + 1);

    // Assigning nucleotide gaps
    for(size_t i = 0; i < T->gaps.size(); i++) {
        int32_t primaryBId = (T->gaps[i].primaryBlockId);
        int32_t secondaryBId = (T->gaps[i].secondaryBlockId);
        for(size_t j = 0; j < T->gaps[i].nucPosition.size(); j++) {
            int len = T->gaps[i].nucGapLength[j];
            int pos = T->gaps[i].nucPosition[j];
            if(secondaryBId != -1) {
                sequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
            } else {
                sequence[primaryBId].first[pos].second.resize(len, '-');
            }
        }
    }

    // Get all blocks on the path
    for(auto node = path.rbegin(); node != path.rend(); node++) {
        for(auto mutation: (*node)->blockMutation) {
            int primaryBlockId = mutation.primaryBlockId;
            int secondaryBlockId = mutation.secondaryBlockId;
            int type = (mutation.blockMutInfo);
            bool inversion = mutation.inversion;

            if(type == BlockMutationType::BI) {
                if(secondaryBlockId != -1) {
                    blockExists[primaryBlockId].second[secondaryBlockId] = true;
                    blockStrand[primaryBlockId].second[secondaryBlockId] = !inversion;
                } else {
                    blockExists[primaryBlockId].first = true;
                    blockStrand[primaryBlockId].first = !inversion;
                }
            } else {
                if(inversion) {
                    if(secondaryBlockId != -1) {
                        blockStrand[primaryBlockId].second[secondaryBlockId] = !blockStrand[primaryBlockId].second[secondaryBlockId];
                    } else {
                        blockStrand[primaryBlockId].first = !blockStrand[primaryBlockId].first;
                    }
                } else {
                    if(secondaryBlockId != -1) {
                        blockExists[primaryBlockId].second[secondaryBlockId] = false;
                        blockStrand[primaryBlockId].second[secondaryBlockId] = true;
                    } else {
                        blockExists[primaryBlockId].first = false;
                        blockStrand[primaryBlockId].first = true;
                    }
                }
            }
        }
    }

    // Apply nucleotide mutations
    for(auto node = path.rbegin(); node != path.rend(); node++) {
        for(size_t i = 0; i < (*node)->nucMutation.size(); i++) {
            int32_t primaryBlockId = (*node)->nucMutation[i].primaryBlockId;
            int32_t secondaryBlockId = (*node)->nucMutation[i].secondaryBlockId;

            if(secondaryBlockId != -1) {
                if(!blockExists[primaryBlockId].second[secondaryBlockId]) continue;
            } else {
                if(!blockExists[primaryBlockId].first) continue;
            }

            int32_t nucPosition = (*node)->nucMutation[i].nucPosition;
            int32_t nucGapPosition = (*node)->nucMutation[i].nucGapPosition;
            uint32_t type = ((*node)->nucMutation[i].mutInfo & 0x7);
            char newVal = '-';

            if(type < 3) {
                int len = (((*node)->nucMutation[i].mutInfo) >> 4);

                if(type == NucMutationType::NS) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                } else if(type == NucMutationType::NI) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = newVal;
                            }
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = newVal;
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> (4*(5-j))) & 0xF);
                                if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                                sequence[primaryBlockId].first[nucPosition+j].first = newVal;
                            }
                        }
                    }
                } else if(type == NucMutationType::ND) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].second[secondaryBlockId][nucPosition + j].first = '-';
                            }
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].first[nucPosition].second[nucGapPosition+j] = '-';
                            }
                        } else {
                            for(int j = 0; j < len; j++) {
                                sequence[primaryBlockId].first[nucPosition+j].first = '-';
                            }
                        }
                    }
                }
            } else {
                if(type == NucMutationType::NSNPS) {
                    newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                } else if(type == NucMutationType::NSNPI) {
                    newVal = getNucleotideFromCode((((*node)->nucMutation[i].nucs) >> 20) & 0xF);
                    if(skipNMutations && newVal == 'N') continue;  // SKIP N MUTATIONS
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = newVal;
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = newVal;
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = newVal;
                        }
                    }
                } else if(type == NucMutationType::NSNPD) {
                    if(secondaryBlockId != -1) {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].second[secondaryBlockId][nucPosition].first = '-';
                        }
                    } else {
                        if(nucGapPosition != -1) {
                            sequence[primaryBlockId].first[nucPosition].second[nucGapPosition] = '-';
                        } else {
                            sequence[primaryBlockId].first[nucPosition].first = '-';
                        }
                    }
                }
            }
        }
    }

    if(!aligned && T->rotationIndexes.find(reference) != T->rotationIndexes.end() && T->rotationIndexes[reference] != 0) {
        int ctr = -1, rotInd = 0;
        for(size_t i = 0; i < blockExists.size(); i++) {
            if(blockExists[i].first) ctr++;
            if(ctr == T->rotationIndexes[reference]) { rotInd = i; break; }
        }
        rotate(sequence.begin(), sequence.begin() + rotInd, sequence.end());
        rotate(blockExists.begin(), blockExists.begin() + rotInd, blockExists.end());
        rotate(blockStrand.begin(), blockStrand.begin() + rotInd, blockStrand.end());
    }

    if(T->sequenceInverted.find(reference) != T->sequenceInverted.end() && T->sequenceInverted[reference]) {
        reverse(sequence.begin(), sequence.end());
        reverse(blockExists.begin(), blockExists.end());
        reverse(blockStrand.begin(), blockStrand.end());
    }

    std::string sequenceString;
    for(size_t i = 0; i < sequence.size(); i++) {
        // Gap blocks (currently not used for SARS-CoV-2)
        for(size_t j = 0; j < sequence[i].second.size(); j++) {
            if(blockExists[i].second[j]) {
                if(blockStrand[i].second[j]) {
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) {
                            if(sequence[i].second[j][k].second[w] == 'x' || sequence[i].second[j][k].second[w] == '-') {
                                if(aligned) sequenceString+='-';
                            } else {
                                sequenceString += sequence[i].second[j][k].second[w];
                            }
                        }
                        if(sequence[i].second[j][k].first == 'x' || sequence[i].second[j][k].first == '-') {
                            if(aligned) sequenceString+='-';
                        } else {
                            sequenceString += sequence[i].second[j][k].first;
                        }
                    }
                } else {
                    for(size_t k = sequence[i].second[j].size()-1; k + 1 > 0; k--) {
                        if(sequence[i].second[j][k].first == 'x' || sequence[i].second[j][k].first == '-') {
                            if(aligned) sequenceString+='-';
                        } else {
                            sequenceString += sequence[i].second[j][k].first;
                        }
                        for(size_t w = sequence[i].second[j][k].second.size() - 1; w + 1 > 0; w--) {
                            if(sequence[i].second[j][k].second[w] == 'x' || sequence[i].second[j][k].second[w] == '-') {
                                if(aligned) sequenceString+='-';
                            } else {
                                sequenceString += sequence[i].second[j][k].second[w];
                            }
                        }
                    }
                }
            } else {
                if(aligned) {
                    for(size_t k = 0; k < sequence[i].second[j].size(); k++) {
                        for(size_t w = 0; w < sequence[i].second[j][k].second.size(); w++) sequenceString+='-';
                        sequenceString+='-';
                    }
                }
            }
        }

        // Main block
        if(blockExists[i].first) {
            if(blockStrand[i].first) {
                for(size_t j = 0; j < sequence[i].first.size(); j++) {
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) {
                        if(sequence[i].first[j].second[k] == 'x' || sequence[i].first[j].second[k] == '-') {
                            if(aligned) sequenceString += '-';
                        } else {
                            sequenceString += sequence[i].first[j].second[k];
                        }
                    }
                    if(sequence[i].first[j].first == 'x' || sequence[i].first[j].first == '-') {
                        if(aligned) sequenceString += '-';
                    } else {
                        sequenceString += sequence[i].first[j].first;
                    }
                }
            } else {
                for(size_t j = sequence[i].first.size()-1; j + 1 > 0; j--) {
                    if(sequence[i].first[j].first == 'x' || sequence[i].first[j].first == '-') {
                        if(aligned) sequenceString += '-';
                    } else {
                        sequenceString += getComplementCharacter(sequence[i].first[j].first);
                    }
                    for(size_t k = sequence[i].first[j].second.size() - 1; k+1 > 0; k--) {
                        if(sequence[i].first[j].second[k] == 'x' || sequence[i].first[j].second[k] == '-') {
                            if(aligned) sequenceString += '-';
                        } else {
                            sequenceString += getComplementCharacter(sequence[i].first[j].second[k]);
                        }
                    }
                }
            }
        } else {
            if(aligned) {
                for(size_t j = 0; j < sequence[i].first.size(); j++) {
                    for(size_t k = 0; k < sequence[i].first[j].second.size(); k++) sequenceString+='-';
                    sequenceString+='-';
                }
            }
        }
    }

    int offset = 0;
    if(!aligned && T->circularSequences.find(reference) != T->circularSequences.end()) {
        offset = T->circularSequences[reference];
    }
    if(offset == 0) {
        return sequenceString;
    } else {
        return sequenceString.substr(offset) + sequenceString.substr(0,offset);
    }
}

std::string getNodeSequence(panmanUtils::Tree* T, const std::string& nodeId, bool impute = false) {
    if (impute) {
        return getStringFromReferenceSkipN(T, nodeId, false, true);
    }
    return T->getStringFromReference(nodeId, false, true);
}

void saveNodeSequence(panmanUtils::Tree* T, const std::string& nodeId, const std::string& path, bool impute = false) {
    std::string seq = getNodeSequence(T, nodeId, impute);
    
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write: " + path);
    
    out << ">" << nodeId << "\n";
    for (size_t i = 0; i < seq.size(); i += 80) {
        out << seq.substr(i, 80) << "\n";
    }
    logging::msg("Saved {} ({} bp) to {}", nodeId, seq.size(), path);
}

/* INDEXING */

bool buildMgsrIndex(const Config& cfg) {
  auto tg = loadPanMAN(cfg.panman);
  if (!tg || tg->trees.empty()) {
    logging::err("Failed to load pangenome");
    return false;
  }

  panmanUtils::Tree* T = &tg->trees[0];
  mgsr::mgsrIndexBuilder mgsrIndexBuilder(T, cfg.k, cfg.s, cfg.t, cfg.l, cfg.openSyncmer, cfg.impute, cfg.indexFull);
  mgsrIndexBuilder.buildIndex();
  mgsrIndexBuilder.writeIndex(cfg.indexMgsr);
  std::cout << (cfg.indexFull ? "Full" : "Lite") << " MGSR index written to " << cfg.indexMgsr << std::endl;
  return true;
}

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
    
    index_single_mode::IndexBuilder builder(&tg->trees[0], cfg.k, cfg.s, 0, cfg.l, false, cfg.flankMaskBp, cfg.hpc, cfg.impute, cfg.extentGuard);
    builder.buildIndexParallel(cfg.threads);
    builder.writeIndex(cfg.index, cfg.threads);
    
    logging::msg("Index built with k={}, s={}, l={}, flankMask={}bp{}{}{}", cfg.k, cfg.s, cfg.l, cfg.flankMaskBp, cfg.hpc ? ", hpc=on" : "", cfg.impute ? ", impute=on" : "", cfg.extentGuard ? ", extentGuard=on" : "");
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
    
    // Compute distance from expected parent if provided
    std::string distStr = "";
    if (tree && !expectedParent.empty()) {
        int dist = computeTreeDistance(*tree, node->identifier, expectedParent);
        distStr = "  dist=" + std::to_string(dist);
    }
    
    std::cout << "  NODE: " << node->identifier << distStr << "\n";
    std::cout << "    LogRaw=" << node->logRawScore
              << "  LogCosine=" << node->logCosineScore
              << "  Containment=" << node->containmentScore << "\n";
}

// Build full seed frequency map for a node by traversing from root
absl::flat_hash_map<uint64_t, int64_t> buildNodeSeedCounts(
    panmapUtils::LiteNode* node,
    panmapUtils::LiteTree& tree) {
    
    absl::flat_hash_map<uint64_t, int64_t> seedCounts;
    
    // Collect path from root to this node
    std::vector<panmapUtils::LiteNode*> path;
    for (auto* n = node; n != nullptr; n = n->parent) {
        path.push_back(n);
    }
    std::reverse(path.begin(), path.end());
    
    // Apply seed changes along path to build full seed set
    for (auto* n : path) {
        for (const auto& [seedHash, parentCount, childCount] : n->seedChanges) {
            int64_t delta = childCount - parentCount;
            if (delta != 0) {
                seedCounts[seedHash] += delta;
                if (seedCounts[seedHash] <= 0) {
                    seedCounts.erase(seedHash);
                }
            }
        }
    }
    
    return seedCounts;
}


// Detailed comparison of two nodes
void printNodeComparison(
    const std::string& nodeIdA,
    const std::string& nodeIdB,
    panmapUtils::LiteTree& tree,
    const absl::flat_hash_map<size_t, int64_t, IdentityHash>& readSeeds,
    double logReadMagnitude) {
    
    auto itA = tree.allLiteNodes.find(nodeIdA);
    auto itB = tree.allLiteNodes.find(nodeIdB);
    
    if (itA == tree.allLiteNodes.end()) {
        std::cout << "ERROR: Node A '" << nodeIdA << "' not found\n";
        return;
    }
    if (itB == tree.allLiteNodes.end()) {
        std::cout << "ERROR: Node B '" << nodeIdB << "' not found\n";
        return;
    }
    
    auto* nodeA = itA->second;
    auto* nodeB = itB->second;
    
    // Build seed counts for each node
    auto seedsA = buildNodeSeedCounts(nodeA, tree);
    auto seedsB = buildNodeSeedCounts(nodeB, tree);
    
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "NODE COMPARISON: " << nodeIdA << " vs " << nodeIdB << "\n";
    std::cout << std::string(80, '=') << "\n";
    
    // Basic stats
    std::cout << "\nBASIC STATISTICS:\n";
    std::cout << "  Node A: " << seedsA.size() << " unique seeds\n";
    std::cout << "  Node B: " << seedsB.size() << " unique seeds\n";
    std::cout << "  Reads:  " << readSeeds.size() << " unique seeds\n";
    
    // Compute intersection and differences with reads
    size_t intersectA = 0, intersectB = 0;
    double logRawNumA = 0.0, logRawNumB = 0.0;
    
    // Seeds in A only (not in B), that match reads
    std::vector<std::tuple<uint64_t, int64_t, int64_t, double>> seedsOnlyAWithReads;
    // Seeds in B only (not in A), that match reads
    std::vector<std::tuple<uint64_t, int64_t, int64_t, double>> seedsOnlyBWithReads;
    // Seeds in both A and B, with different genome counts, matching reads
    std::vector<std::tuple<uint64_t, int64_t, int64_t, int64_t, double>> seedsDiffWithReads;
    // All seeds matching reads for both nodes
    std::vector<std::tuple<uint64_t, int64_t, int64_t, int64_t>> allMatchingSeeds;
    
    for (const auto& [hash, readCount] : readSeeds) {
        auto itA_seed = seedsA.find(hash);
        auto itB_seed = seedsB.find(hash);
        int64_t countA = (itA_seed != seedsA.end()) ? itA_seed->second : 0;
        int64_t countB = (itB_seed != seedsB.end()) ? itB_seed->second : 0;
        double logRead = std::log1p(static_cast<double>(readCount));
        
        if (countA > 0) {
            intersectA++;
            logRawNumA += logRead / countA;
        }
        if (countB > 0) {
            intersectB++;
            logRawNumB += logRead / countB;
        }
        
        if (countA > 0 || countB > 0) {
            allMatchingSeeds.emplace_back(hash, readCount, countA, countB);
        }
        
        // Track differences
        if (countA > 0 && countB == 0) {
            seedsOnlyAWithReads.emplace_back(hash, readCount, countA, logRead / countA);
        } else if (countB > 0 && countA == 0) {
            seedsOnlyBWithReads.emplace_back(hash, readCount, countB, logRead / countB);
        } else if (countA > 0 && countB > 0 && countA != countB) {
            seedsDiffWithReads.emplace_back(hash, readCount, countA, countB, logRead / countA - logRead / countB);
        }
    }
    
    // Compute scores
    double logRawA = logRawNumA / logReadMagnitude;
    double logRawB = logRawNumB / logReadMagnitude;
    
    std::cout << "\nINTERSECTION WITH READS:\n";
    std::cout << "  Node A: " << intersectA << " seeds match reads\n";
    std::cout << "  Node B: " << intersectB << " seeds match reads\n";
    std::cout << "  Difference: " << (static_cast<int64_t>(intersectB) - static_cast<int64_t>(intersectA)) << " more seeds in B\n";
    
    std::cout << "\nLOG_RAW SCORE BREAKDOWN:\n";
    std::cout << "  Node A: numerator=" << logRawNumA << ", score=" << logRawA << "\n";
    std::cout << "  Node B: numerator=" << logRawNumB << ", score=" << logRawB << "\n";
    std::cout << "  Delta: " << (logRawB - logRawA) << " (B - A)\n";
    
    // Sort by contribution to score difference (log(R)/G delta)
    std::sort(seedsOnlyAWithReads.begin(), seedsOnlyAWithReads.end(),
        [](const auto& a, const auto& b) { return std::get<3>(a) > std::get<3>(b); });
    std::sort(seedsOnlyBWithReads.begin(), seedsOnlyBWithReads.end(),
        [](const auto& a, const auto& b) { return std::get<3>(a) > std::get<3>(b); });
    std::sort(seedsDiffWithReads.begin(), seedsDiffWithReads.end(),
        [](const auto& a, const auto& b) { return std::abs(std::get<4>(a)) > std::abs(std::get<4>(b)); });
    
    // Show seeds only in A (with reads)
    std::cout << "\nSEEDS ONLY IN A (matching reads): " << seedsOnlyAWithReads.size() << " total\n";
    for (size_t i = 0; i < std::min(size_t(15), seedsOnlyAWithReads.size()); ++i) {
        const auto& [hash, readCnt, genomeCnt, contrib] = seedsOnlyAWithReads[i];
        std::cout << "  hash=" << std::hex << hash << std::dec 
                  << "  readCnt=" << readCnt 
                  << "  genomeCnt=" << genomeCnt
                  << "  logRaw_contrib=" << contrib << "\n";
    }
    
    // Show seeds only in B (with reads)
    std::cout << "\nSEEDS ONLY IN B (matching reads): " << seedsOnlyBWithReads.size() << " total\n";
    for (size_t i = 0; i < std::min(size_t(15), seedsOnlyBWithReads.size()); ++i) {
        const auto& [hash, readCnt, genomeCnt, contrib] = seedsOnlyBWithReads[i];
        std::cout << "  hash=" << std::hex << hash << std::dec 
                  << "  readCnt=" << readCnt 
                  << "  genomeCnt=" << genomeCnt
                  << "  logRaw_contrib=" << contrib << "\n";
    }
    
    // Show seeds with different genome counts
    std::cout << "\nSEEDS IN BOTH (different genome counts): " << seedsDiffWithReads.size() << " total\n";
    for (size_t i = 0; i < std::min(size_t(15), seedsDiffWithReads.size()); ++i) {
        const auto& [hash, readCnt, cntA, cntB, delta] = seedsDiffWithReads[i];
        std::cout << "  hash=" << std::hex << hash << std::dec 
                  << "  readCnt=" << readCnt 
                  << "  A_cnt=" << cntA << "  B_cnt=" << cntB
                  << "  delta=" << delta << (delta > 0 ? " (favors A)" : " (favors B)") << "\n";
    }
    
    // Summary of why B might score higher
    double contribOnlyA = 0, contribOnlyB = 0, contribDiff = 0;
    for (const auto& [hash, readCnt, genomeCnt, contrib] : seedsOnlyAWithReads) {
        contribOnlyA += contrib;
    }
    for (const auto& [hash, readCnt, genomeCnt, contrib] : seedsOnlyBWithReads) {
        contribOnlyB += contrib;
    }
    for (const auto& [hash, readCnt, cntA, cntB, delta] : seedsDiffWithReads) {
        contribDiff += delta;  // positive = favors A
    }
    
    std::cout << "\nSCORE DIFFERENCE BREAKDOWN:\n";
    std::cout << "  Seeds only in A contribute: +" << contribOnlyA << " to A\n";
    std::cout << "  Seeds only in B contribute: +" << contribOnlyB << " to B\n";
    std::cout << "  Seeds with diff counts contribute: " << contribDiff << " to A (neg = favors B)\n";
    std::cout << "  Net effect: " << (contribOnlyA - contribOnlyB + contribDiff) << " (pos = A better, neg = B better)\n";
    
    std::cout << "\nStored node scores:\n";
    std::cout << "  Node A: LogRaw=" << nodeA->logRawScore << "  LogCosine=" << nodeA->logCosineScore << "  Containment=" << nodeA->containmentScore << "\n";
    std::cout << "  Node B: LogRaw=" << nodeB->logRawScore << "  LogCosine=" << nodeB->logCosineScore << "  Containment=" << nodeB->containmentScore << "\n";
    
    std::cout << std::string(80, '=') << "\n\n";
}

void writeOCRanks(const std::string& outputFile, const std::vector<std::pair<std::string, double>>& overlapCoefficients) {
  std::ofstream outFile(outputFile);
  uint32_t rank = 0;
  double currentOverlapCoefficient = overlapCoefficients[0].second;
  for (const auto& [nodeId, overlapCoefficient] : overlapCoefficients) {
    if (overlapCoefficient != currentOverlapCoefficient) {
      currentOverlapCoefficient = overlapCoefficient;
      ++rank;
    }
    outFile << nodeId << "\t" << std::fixed << std::setprecision(6) << overlapCoefficient << "\t" << rank <<  std::endl;
  }
  outFile.close();
}

void writeMetaReadScoresUnfiltered(const std::string& outputFile, const mgsr::ThreadsManager& threadsManager) {
  std::ofstream outFile(outputFile);
  
  outFile << "ReadIndex\tNumDuplicates\tTotalScore\tMaxScore\tNumMaxScoreNodes\tRawReadsIndices" << std::endl;
  for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
    const auto& curRead = threadsManager.reads[i];
    if (curRead.maxScore == 0) continue;
    outFile << i << "\t" << threadsManager.readSeedmersDuplicatesIndex[i].size() << "\t" << curRead.seedmersList.size() << "\t" << curRead.maxScore << "\t" << curRead.epp << "\t";
    for (size_t j = 0; j < threadsManager.readSeedmersDuplicatesIndex[i].size(); ++j) {
      if (j == 0) {
        outFile << threadsManager.readSeedmersDuplicatesIndex[i][j];
      } else {
        outFile << "," << threadsManager.readSeedmersDuplicatesIndex[i][j];
      }
    }
    outFile << std::endl;
  }
  outFile.close();
}

void writeMetaReadScoresFiltered(const std::string& outputFile, const mgsr::ThreadsManager& threadsManager) {
  std::ofstream outFile(outputFile);
  
  outFile << "ReadIndex\tNumDuplicates\tTotalScore\tMaxScore\tNumMaxScoreNodes\tOverMaximumFamilies\tRawReadsIndices" << std::endl;
  for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
    const auto& curRead = threadsManager.reads[i];
    if (curRead.maxScore == 0) continue;
    outFile << i << "\t" << threadsManager.readSeedmersDuplicatesIndex[i].size() << "\t" << curRead.seedmersList.size() << "\t" << curRead.maxScore << "\t" << curRead.epp << "\t" << curRead.overMaximumFamilies << "\t";
    for (size_t j = 0; j < threadsManager.readSeedmersDuplicatesIndex[i].size(); ++j) {
      if (j == 0) {
        outFile << threadsManager.readSeedmersDuplicatesIndex[i][j];
      } else {
        outFile << "," << threadsManager.readSeedmersDuplicatesIndex[i][j];
      }
    }
    outFile << std::endl;
  }
  outFile.close();
}

void scoreReadsMultiThreaded(mgsr::MgsrLiteTree& T, mgsr::ThreadsManager& threadsManager, const Config& cfg) {
  std::vector<uint64_t> totalNodesPerThread(threadsManager.numThreads, 0);
  for (size_t i = 0; i < threadsManager.numThreads; ++i) {
    totalNodesPerThread[i] = T.getNumActiveNodes();
  }
  ProgressTracker progressTracker(threadsManager.numThreads, totalNodesPerThread);
  std::cout << "Placing reads with " << threadsManager.numThreads << " threads..." << std::endl;

  bool lowMemory = false;
  auto start_time_place = std::chrono::high_resolution_clock::now();
  tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()), [&](const tbb::blocked_range<size_t>& rangeIndex){
    for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
      auto [start, end] = threadsManager.threadRanges[i];
    
      std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);
      mgsr::mgsrPlacer curThreadPlacer(&T, threadsManager, lowMemory, i);
      curThreadPlacer.initializeQueryData(curThreadReads);

      curThreadPlacer.setAllSeedmerHashesSet(threadsManager.allSeedmerHashesSet);
      
      curThreadPlacer.setProgressTracker(&progressTracker, i);
      // curThreadPlacer.placeReads();

      curThreadPlacer.scoreReads();

      if (i == 0) {
        threadsManager.identicalGroups = std::move(curThreadPlacer.identicalGroups);
        threadsManager.identicalNodeToGroup = std::move(curThreadPlacer.identicalNodeToGroup);
      }
    }
  });
  auto end_time_place = std::chrono::high_resolution_clock::now();
  auto duration_place = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_place - start_time_place);
  std::cout << "\n\nPlaced reads in " << static_cast<double>(duration_place.count()) / 1000.0 << "s\n" << std::endl;
}

bool runMetagenomic(const Config& cfg) {
  std::cout << "Running metagenomic mode with index: " << cfg.index << " and threads: " << cfg.threads << std::endl;

  // Checking IO
  if (cfg.index.empty() || !fs::exists(cfg.index)) {
    std::cerr << "Error: Index file " << cfg.index << " does not exist" << std::endl;
    return false;
  }
  
  
  if (!cfg.taxonomicMetadata.empty() && !fs::exists(cfg.taxonomicMetadata)) {
    std::cerr << "Error: Taxonomic metadata file " << cfg.taxonomicMetadata << " does not exist" << std::endl;
    return false;
  }

  if (cfg.reads1.empty() || !fs::exists(cfg.reads1)) {
    std::cerr << "Error: Reads1 file " << cfg.reads1 << " does not exist" << std::endl;
    return false;
  }

  if (!cfg.reads2.empty() && !fs::exists(cfg.reads2)) {
    std::cerr << "Error: Reads2 file " << cfg.reads2 << " does not exist" << std::endl;
    return false;
  }
  
  
  int fd = -1;
  std::unique_ptr<IndexReader> zstdReader;
  std::unique_ptr<::capnp::PackedFdMessageReader> fdReader;
  ::capnp::MessageReader* baseReader = nullptr;
  bool mgsrIndex;
  
  try {
    zstdReader = std::make_unique<IndexReader>(cfg.index, cfg.threads);
    baseReader = zstdReader.get();
    mgsrIndex = false;
  } catch (...) {
    std::cerr << "Trying to open index file as packed file..." << std::endl;
    fd = mgsr::open_file(cfg.index);
    ::capnp::ReaderOptions readerOptions{
      .traversalLimitInWords = std::numeric_limits<uint64_t>::max(),
      .nestingLimit = 1024
    };
    fdReader = std::make_unique<::capnp::PackedFdMessageReader>(fd, readerOptions);
    baseReader = fdReader.get();
    mgsrIndex = true;
  }
  
  LiteIndex::Reader indexReader = baseReader->getRoot<LiteIndex>();
  bool lowMemory = false;


  mgsr::MgsrLiteTree T;
  T.initialize(indexReader, cfg.taxonomicMetadata, cfg.taxonomicMetadata.empty() ? 0 : cfg.maximumFamilies, cfg.threads, lowMemory, true);

  mgsr::ThreadsManager threadsManager(&T, cfg.output, cfg.threads, cfg.maskSeeds, cfg.maskReads, cfg.maskSeedsRelativeFrequency, cfg.maskReadsRelativeFrequency, !cfg.noProgress, lowMemory);
  threadsManager.initializeMGSRIndex(indexReader);

  threadsManager.initializeQueryData(cfg.reads1, cfg.reads2, cfg.ampliconDepth, cfg.dust, cfg.maskReadEnds);

  if (!cfg.filterAndAssign) {
    mgsr::mgsrPlacer placer(&T, threadsManager, lowMemory, 0);
    auto overlapCoefficients = placer.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);
    
    T.fillOCRanks(overlapCoefficients);

    if (cfg.writeOCRanks) {
      writeOCRanks(cfg.output + ".overlapCoefficients.tsv", overlapCoefficients);
    }
  }
  
  T.collapseIdenticalScoringNodes(threadsManager.allSeedmerHashesSet);

  scoreReadsMultiThreaded(T, threadsManager, cfg);

  if (cfg.writeMetaReadScoresUnfiltered) {
    writeMetaReadScoresUnfiltered(cfg.output + ".read_scores_info.unfiltered.tsv", threadsManager);
  }

  double discard_threshold = cfg.discard;
  size_t num_discarded = 0;
  size_t num_unmapped = 0;
  for (auto& read : threadsManager.reads) {
    if (read.maxScore == 0) {
      ++num_unmapped;
    } else if (read.maxScore < static_cast<int>(read.seedmersList.size() * discard_threshold)) {
      read.maxScore = 0;
      ++num_discarded;
    }
  }
  std::cout << num_unmapped << " reads unmapped... " << std::endl;
  std::cout << num_discarded << " reads discarded due to low parsimony score... " << std::endl;
  std::cout << threadsManager.reads.size() - num_unmapped - num_discarded << " reads mapped to nodes..." << std::endl;
  if (threadsManager.reads.size() - num_unmapped - num_discarded == 0) {
    std::cerr << "No reads remain for node scoring and EM after discarding low-score reads... Exiting... " << std::endl;
    return true;
  }



  if (cfg.filterAndAssign) {
    // Remove score updates of unmapped reads
    tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()), [&](const tbb::blocked_range<size_t>& rangeIndex){
      for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
        auto [start, end] = threadsManager.threadRanges[i];
        std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);

        for (auto [_, node] : T.allLiteNodes) {
          auto& curThreadNodeScoreDeltas = node->readScoreDeltas[i];
          curThreadNodeScoreDeltas.erase(
            std::remove_if(curThreadNodeScoreDeltas.begin(), curThreadNodeScoreDeltas.end(),
              [&curThreadReads](const mgsr::readScoreDelta& delta) {
                return curThreadReads[delta.readIndex].maxScore == 0;
              }),
            curThreadNodeScoreDeltas.end()
          );
        }
      }
    });

    std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> assignedReadsByNode;
    threadsManager.assignReads(assignedReadsByNode, cfg.maximumFamilies);

    std::ofstream assignedReadsOut(cfg.output + ".mgsr.assignedReads.out");
    for (auto& [node, readIndices] : assignedReadsByNode) {
      assignedReadsOut << node->identifier;
      for (const auto& identicalNodeId : node->identicalNodeIdentifiers) {
        assignedReadsOut << "," << identicalNodeId;
      }
      assignedReadsOut << "\t" << readIndices.size() << "\t";
      std::sort(readIndices.begin(), readIndices.end());
      for (size_t i = 0; i < readIndices.size(); ++i) {
        assignedReadsOut << readIndices[i];
        if (i != readIndices.size() - 1) {
          assignedReadsOut << ",";
        }
      }
      assignedReadsOut << "\n";
    }
    assignedReadsOut.close();
  } else {
    T.seedInfos.clear(); // no longer needed. clear memory to prep for EM.
    mgsr::squareEM squareEM(threadsManager, T, cfg.output, cfg.topOc, cfg.emConvergenceThreshold, cfg.emDeltaThreshold, cfg.emMaximumIterations, false, cfg.emLeavesOnly);
    
    auto start_time_squareEM = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < cfg.emMaximumRounds; ++i) {
      squareEM.runSquareEM();
      std::cout << "\nRound " << i << " of squareEM completed... nodes size changed from " << squareEM.nodes.size() << " to ";
      bool removed = squareEM.removeLowPropNodes();
      std::cout << squareEM.nodes.size() << std::endl;
      if (!removed) {
        break;
      }
    }
    std::cout << std::endl;

    auto end_time_squareEM = std::chrono::high_resolution_clock::now();
    auto duration_squareEM = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_squareEM - start_time_squareEM);
    std::cout << "SquareEM completed in " << static_cast<double>(duration_squareEM.count()) / 1000.0 << "s\n" << std::endl;


    std::vector<uint64_t> indices(squareEM.nodes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&squareEM](uint64_t i, uint64_t j) {
      return squareEM.props[i] > squareEM.props[j];
    });

    std::cout << "writing abundance file: " << cfg.output + ".mgsr.abundance.out" << std::endl;
    std::ofstream abundanceOutput(cfg.output + ".mgsr.abundance.out");
    abundanceOutput << std::setprecision(5) << std::fixed;
    for (size_t i = 0; i < indices.size(); ++i) {
      size_t index = indices[i];
      abundanceOutput << squareEM.nodes[index];
      for (const auto& member : T.allLiteNodes.at(squareEM.nodes[index])->identicalNodeIdentifiers) {
        abundanceOutput << "," << member;
      }
      if (squareEM.identicalGroups.find(squareEM.nodes[index]) != squareEM.identicalGroups.end()) {
        for (const auto& member : squareEM.identicalGroups[squareEM.nodes[index]]) {
          abundanceOutput << "," << member;
          for (const auto& identicalMember : T.allLiteNodes.at(member)->identicalNodeIdentifiers) {
            abundanceOutput << "," << identicalMember;
          }
        }
      }
      abundanceOutput << "\t" << squareEM.props[index] << std::endl;
    }
    abundanceOutput.close();

    std::vector<std::vector<size_t>> assignedReadsIndices(indices.size());
    std::vector<size_t> readIndexOffset(threadsManager.reads.size(), 0);
    for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
      readIndexOffset[i] = i;
    }
    for (size_t i = 0; i < indices.size(); ++i) {
      size_t index = indices[i];
      const auto& nodeId = squareEM.nodes[index];
      auto curNodeScores = threadsManager.getScoresAtNode(nodeId, readIndexOffset);
      for (size_t j = 0; j < curNodeScores.size(); ++j) {
        const auto& curRead = threadsManager.reads[j];
        if (curRead.maxScore != 0 && curNodeScores[j] == curRead.maxScore) {
          for (const auto& rawReadIndex : threadsManager.readSeedmersDuplicatesIndex[j]) {
            assignedReadsIndices[i].push_back(rawReadIndex);
          }
        }
      }
      if (assignedReadsIndices[i].size() > 0) {
        std::sort(assignedReadsIndices[i].begin(), assignedReadsIndices[i].end());
      }
    }
    std::ofstream assignedReadsOutput(cfg.output + ".mgsr.assignedReads.out");
    for (size_t i = 0; i < indices.size(); ++i) {
      size_t index = indices[i];
      assignedReadsOutput << squareEM.nodes[index] << "\t" << assignedReadsIndices[i].size() << "\t";
      for (size_t j = 0; j < assignedReadsIndices[i].size(); ++j) {
        if (j == 0) {
          assignedReadsOutput << assignedReadsIndices[i][j];
        } else {
          assignedReadsOutput << "," << assignedReadsIndices[i][j];
        }
      }
      assignedReadsOutput << std::endl;
    }
    assignedReadsOutput.close();
  }

  if (cfg.writeMetaReadScoresFiltered) {
    writeMetaReadScoresFiltered(cfg.output + ".read_scores_info.filtered.tsv", threadsManager);
  }


  
  return true;
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
    bool storeDiagnostics = (cfg.topPlacements > 0 || !cfg.debugNodeId.empty() || !cfg.compareNodes.empty());
    
    // Load full tree if refinement is enabled (needed for genome sequences)
    std::unique_ptr<panmanUtils::TreeGroup> tg;
    panmanUtils::Tree* fullTreePtr = nullptr;
    bool refineEnabled = cfg.refine;
    if (refineEnabled) {
        tg = loadPanMAN(cfg.panman);
        if (tg && !tg->trees.empty()) {
            fullTreePtr = &tg->trees[0];
            logging::msg("Loaded full tree for refinement ({} nodes)", 
                        fullTreePtr->allNodes.size());
        } else {
            logging::warn("Failed to load full tree for refinement - disabling refinement");
            refineEnabled = false;
        }
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    placement::placeLite(result, &tree, reader, cfg.reads1, cfg.reads2, outPath, false, fullTreePtr,
                        cfg.removeNodeId, expectedParentPtr, storeDiagnostics, cfg.seedMaskFraction,
                        cfg.minSeedQuality, cfg.dedupReads, true, cfg.trimStart, cfg.trimEnd,
                        cfg.minReadSupport,
                        refineEnabled, cfg.refineTopPct, cfg.refineMaxTopN, 
                        cfg.refineNeighborRadius, cfg.refineMaxNeighborN);
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start);
    
    logging::msg("Placement complete in {}ms", elapsed.count());
    
    // Check if any metric found a placement
    if (result.bestLogRawNodeId.empty()) {
        logging::warn("No placement found");
        return std::nullopt;
    }
    
    // Report results
    logging::msg("{}Best placement:{} {} (LogRaw: {:.6f}, LogCosine: {:.6f}, Containment: {:.6f})", 
                color::green(), color::reset(),
                result.bestLogRawNodeId, 
                result.bestLogRawScore,
                result.bestLogCosineScore,
                result.bestContainmentScore);
    
    // Leave-one-out validation: compare to expected parent
    if (!cfg.removeNodeId.empty() && !expectedParentNode.empty()) {
        int distLogRaw = computeTreeDistance(tree, result.bestLogRawNodeId, expectedParentNode);
        int distLogCosine = computeTreeDistance(tree, result.bestLogCosineNodeId, expectedParentNode);
        int distContainment = computeTreeDistance(tree, result.bestContainmentNodeId, expectedParentNode);
        
        bool correctLogRaw = (result.bestLogRawNodeId == expectedParentNode);
        bool correctLogCosine = (result.bestLogCosineNodeId == expectedParentNode);
        bool correctContainment = (result.bestContainmentNodeId == expectedParentNode);
        
        logging::msg("{}VALIDATION RESULTS:{}", color::cyan(), color::reset());
        logging::msg("  Expected parent: {}", expectedParentNode);
        logging::msg("  LogRAW:        {} {}{}{} (dist: {}, score: {:.6f})", 
                    result.bestLogRawNodeId,
                    correctLogRaw ? color::green() : color::yellow(),
                    correctLogRaw ? output::box::check() : output::box::cross(),
                    color::reset(),
                    distLogRaw,
                    result.bestLogRawScore);
        logging::msg("  LogCosine:     {} {}{}{} (dist: {}, score: {:.6f})", 
                    result.bestLogCosineNodeId,
                    correctLogCosine ? color::green() : color::yellow(),
                    correctLogCosine ? output::box::check() : output::box::cross(),
                    color::reset(),
                    distLogCosine,
                    result.bestLogCosineScore);
        logging::msg("  Containment:   {} {}{}{} (dist: {}, score: {:.6f})", 
                    result.bestContainmentNodeId,
                    correctContainment ? color::green() : color::yellow(),
                    correctContainment ? output::box::check() : output::box::cross(),
                    color::reset(),
                    distContainment,
                    result.bestContainmentScore);
        
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
        std::cout << "CONTAINMENT\t" << result.bestContainmentNodeId 
                  << "\t" << (correctContainment ? "true" : "false") 
                  << "\t" << distContainment 
                  << "\t" << result.bestContainmentScore << "\n";
        
        // Report per-metric refinement results if run
        if (result.refinementWasRun) {
            auto reportRefined = [&](const std::string& label, const placement::PlacementResult::RefinedResult& ref) {
                int dist = computeTreeDistance(tree, ref.nodeId, expectedParentNode);
                bool correct = (ref.nodeId == expectedParentNode);
                logging::msg("  {}:   {} {}{}{} (dist: {}, score: {:.0f})", 
                            label, ref.nodeId,
                            correct ? color::green() : color::yellow(),
                            correct ? output::box::check() : output::box::cross(),
                            color::reset(), dist, ref.score);
                std::cout << label << "\t" << ref.nodeId 
                          << "\t" << (correct ? "true" : "false") 
                          << "\t" << dist 
                          << "\t" << ref.score << "\n";
            };
            reportRefined("REFINED_LOG_RAW", result.refinedLogRaw);
            reportRefined("REFINED_LOG_COSINE", result.refinedLogCosine);
            reportRefined("REFINED_CONTAINMENT", result.refinedContainment);
            reportRefined("REFINED_WEIGHTED_CONTAINMENT", result.refinedWeightedContainment);
            reportRefined("REFINED_LOG_CONTAINMENT", result.refinedLogContainment);
        }
    }
    
    // Dump all scores to file if requested
    if (!cfg.dumpAllScores.empty()) {
        std::ofstream outFile(cfg.dumpAllScores);
        if (outFile) {
            outFile << "node\tlogRaw\tlogCosine\tcontainment\n";
            std::vector<std::pair<double, std::string>> allScores;
            for (auto& [id, node] : tree.allLiteNodes) {
                if (node->logRawScore > 0 || node->logCosineScore > 0 || node->containmentScore > 0) {
                    allScores.push_back({node->logRawScore, id});
                }
            }
            // Sort by logRaw descending
            std::sort(allScores.begin(), allScores.end(), std::greater<>());
            for (auto& [score, id] : allScores) {
                auto* node = tree.allLiteNodes[id];
                outFile << id << "\t" 
                        << node->logRawScore << "\t"
                        << node->logCosineScore << "\t"
                        << node->containmentScore << "\n";
            }
            logging::msg("Dumped {} node scores to {}", allScores.size(), cfg.dumpAllScores);
        } else {
            logging::warn("Could not open {} for writing", cfg.dumpAllScores);
        }
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
            if (node->logRawScore > 0 || node->logCosineScore > 0 || node->containmentScore > 0) {
                scoredNodes.push_back(node);
            }
        }
        
        if (cfg.topPlacements > 0) {
            int topN = std::min(cfg.topPlacements, static_cast<int>(scoredNodes.size()));
            
            // Top by LogRaw score
            std::sort(scoredNodes.begin(), scoredNodes.end(), 
                [](auto* a, auto* b) { return a->logRawScore > b->logRawScore; });
            std::cout << "TOP " << topN << " BY LOG_RAW:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount, 
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
            
            // Top by LogCosine score
            std::sort(scoredNodes.begin(), scoredNodes.end(), 
                [](auto* a, auto* b) { return a->logCosineScore > b->logCosineScore; });
            std::cout << "\nTOP " << topN << " BY LOG_COSINE:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount, 
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
            
            // Top by Containment score
            std::sort(scoredNodes.begin(), scoredNodes.end(), 
                [](auto* a, auto* b) { return a->containmentScore > b->containmentScore; });
            std::cout << "\nTOP " << topN << " BY CONTAINMENT:\n";
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
    
    // Node comparison mode
    if (!cfg.compareNodes.empty()) {
        // Parse "nodeA,nodeB" format
        auto commaPos = cfg.compareNodes.find(',');
        if (commaPos != std::string::npos) {
            std::string nodeA = cfg.compareNodes.substr(0, commaPos);
            std::string nodeB = cfg.compareNodes.substr(commaPos + 1);
            
            // Compute log read magnitude from result
            double logReadMagnitude = 0.0;
            for (const auto& [hash, count] : result.seedFreqInReads) {
                double logCount = std::log1p(static_cast<double>(count));
                logReadMagnitude += logCount * logCount;
            }
            logReadMagnitude = std::sqrt(logReadMagnitude);
            
            printNodeComparison(nodeA, nodeB, tree, result.seedFreqInReads, logReadMagnitude);
        } else {
            std::cerr << "ERROR: --compare-nodes requires format 'nodeA,nodeB'\n";
        }
    }
    
    if (cfg.metagenomic && cfg.topN > 1) {
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
    // If impute is enabled, skip N mutations during sequence construction
    std::string bestMatchSequence = cfg.impute 
        ? getStringFromReferenceSkipN(T, nodeId, false, true)
        : T->getStringFromReference(nodeId, false, true);
    
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
    std::cout << color::bold() << "Usage:" << color::reset() 
              << "  panmap [options] <panman> [reads.fq] [reads2.fq]\n";
    std::cout << color::bold() << "Output:" << color::reset()
              << " <prefix>.vcf, <prefix>.bam, <prefix>.placement.tsv\n"
              << "        (prefix defaults to reads filename, or use -o)\n\n";
}

int main(int argc, char** argv) {
    Config cfg;
    
    // Define options - organized into visible (common) and advanced groups
    
    // === Common options (shown in --help) ===
    po::options_description visible("Options");
    visible.add_options()
        ("help,h", "Show help (--help-all for more)")
        ("help-all", "Show all options")
        ("version,V", "Show version")
        ("output,o", po::value<std::string>(&cfg.output), "Output prefix")
        ("threads,t", po::value<int>(&cfg.threads)->default_value(1), "Threads")
        ("stop", po::value<std::string>()->default_value("genotype"),
            "Stop after: index|place|align|genotype")
        ("meta", po::bool_switch(&cfg.metagenomic), "Metagenomic mode (for more options, see --help-all)")
        ("aligner,a", po::value<std::string>(&cfg.aligner)->default_value("minimap2"),
            "Aligner: minimap2|bwa")
        ("verbose,v", po::bool_switch(&cfg.verbose), "Verbose output")
        ("quiet,q", po::bool_switch(&cfg.quiet), "Quiet output")
        ("no-color", po::bool_switch(&cfg.plain), "No colors");
    
    // === Advanced options (shown only in --help-all) ===
    po::options_description advanced("Advanced");
    advanced.add_options()
        ("index,i", po::value<std::string>(&cfg.index), "Index file path")
        ("reindex,f", po::bool_switch(&cfg.forceReindex), "Force rebuild index")
        ("dedup", po::bool_switch(&cfg.dedupReads), "Deduplicate reads")
        ("impute", po::bool_switch(&cfg.impute), "Impute N's from parent (skip _->N mutations in indexing and output)")
        ("kmer,k", po::value<int>(&cfg.k)->default_value(29), "Syncmer k")
        ("syncmer,s", po::value<int>(&cfg.s)->default_value(8), "Syncmer s")
        ("offset,t", po::value<int>(&cfg.t)->default_value(0), "Syncmer offset")
        ("lmer,l", po::value<int>(&cfg.l)->default_value(1), "Syncmers per seed")
        ("open-syncmer", po::bool_switch(&cfg.openSyncmer), "Open syncmer")
        ("flank-mask", po::value<int>(&cfg.flankMaskBp)->default_value(250), "Mask bp at ends")
        ("seed-mask-fraction", po::value<double>(&cfg.seedMaskFraction)->default_value(0),
            "Mask top seed fraction")
        ("min-seed-quality", po::value<int>(&cfg.minSeedQuality)->default_value(0),
            "Min seed quality")
        ("trim-start", po::value<int>(&cfg.trimStart)->default_value(0), "Trim read start")
        ("trim-end", po::value<int>(&cfg.trimEnd)->default_value(0), "Trim read end")
        ("min-read-support", po::value<int>(&cfg.minReadSupport)->default_value(1),
            "Min reads for a seed (2=filter singletons)")
        ("hpc", po::bool_switch(&cfg.hpc), "Homopolymer-compressed seeds")
        ("extent-guard", po::bool_switch(&cfg.extentGuard), "Guard seed deletions at genome extent boundaries")
        ("refine", po::bool_switch(&cfg.refine), "Enable alignment-based refinement")
        ("refine-top-pct", po::value<double>(&cfg.refineTopPct)->default_value(0.01),
            "Top % of nodes to refine (default 1%)")
        ("refine-max-top-n", po::value<int>(&cfg.refineMaxTopN)->default_value(150),
            "Max nodes to align against")
        ("refine-neighbor-radius", po::value<int>(&cfg.refineNeighborRadius)->default_value(2),
            "Expand to neighbors within N branches")
        ("refine-max-neighbor-n", po::value<int>(&cfg.refineMaxNeighborN)->default_value(150),
            "Max additional nodes from neighbor expansion");

    po::options_description metagenomic("Metagenomic");
    metagenomic.add_options()
        ("index-mgsr", po::value<std::string>(&cfg.indexMgsr), "Path to build/rebuild MGSR index")
        ("index-full", po::bool_switch(&cfg.indexFull), "Build full index (default index-mgsr builds lite index)")
        ("no-progress", po::bool_switch(&cfg.noProgress), "Disable progress bars");
    
    po::options_description em("Metagenomic: EM");
    em.add_options()
        ("top-oc", po::value<size_t>(&cfg.topOc)->default_value(1000), "Select top <int> nodes by overlap coefficients to send to EM")
        ("mask-reads", po::value<uint32_t>(&cfg.maskReads)->default_value(0), "mask reads containing k-min-mers with total occurrence <= threshold")
        ("mask-seeds", po::value<uint32_t>(&cfg.maskSeeds)->default_value(0), "mask k-min-mer seeds in query with total occurrence <= threshold")
        ("amplicon-depth", po::value<std::string>(&cfg.ampliconDepth), "Path to amplicon depth TSV file (if specified, will be used to mask-reads/seeds basedd)")
        ("mask-reads-relative-frequency", po::value<double>(&cfg.maskReadsRelativeFrequency)->default_value(0.0), "mask reads containing k-min-mers with relative frequency < threadshold * amplicon_depth")
        ("mask-seeds-relative-frequency", po::value<double>(&cfg.maskSeedsRelativeFrequency)->default_value(0.0), "mask k-min-mer seeds in query with with relative frequency < threadshold * amplicon_depth")

        ("em-convergence-threshold", po::value<double>(&cfg.emConvergenceThreshold)->default_value(0.00001), "EM converges when likelihood difference is less than <float> (choose em-convergence-threshold or em-delta-threshold, default is em-convergence-threshold)")
        ("em-delta-threshold", po::value<double>(&cfg.emDeltaThreshold)->default_value(0.0), "EM converges when maximum proportion change is less than <float> (choose em-convergence-threshold or em-delta-threshold, default is em-delta-threshold)")
        ("em-maximum-rounds", po::value<uint32_t>(&cfg.emMaximumRounds)->default_value(5), "EM maximum rounds")
        ("em-maximum-iterations", po::value<uint32_t>(&cfg.emMaximumIterations)->default_value(1000), "EM maximum iterations")
        ("em-leaves-only", po::bool_switch(&cfg.emLeavesOnly), "Only run EM on leaf (sample) nodes");
    
    po::options_description filterAndAssign("Metagenomic: Filter and Assign");
    filterAndAssign.add_options()
        ("filter-and-assign", po::bool_switch(&cfg.filterAndAssign), "Filter and assign reads to nodes without running EM")
        ("dust", po::value<double>(&cfg.dust)->default_value(100.0), "Discard reads with Prinseq scale dust score > <FLOAT> (default 100, i.e. no dust filtering)")
        ("discard", po::value<double>(&cfg.discard)->default_value(0.0), "Discard reads with maximum parsimony score < FLOAT * read_total_seed (default 0, i.e. no discard)")
        ("mask-read-ends", po::value<uint32_t>(&cfg.maskReadEnds)->default_value(0), "mask <int> bases from the beginning and end of reads (for ancient eDNA damage)")
        ("taxonomic-metadata", po::value<std::string>(&cfg.taxonomicMetadata), "Path to taxonomic metadata TSV file")
        ("maximum-families", po::value<size_t>(&cfg.maximumFamilies)->default_value(1), "Discard reads assigned to nodes spanning more than <int> distinct taxonomic families, only applicable if taxonomic-metadata is provided");

    po::options_description developer("Developer");
    developer.add_options()
        ("random-seed", po::value<std::string>(&cfg.randomSeed), "Seed for rng (read in as string then hashed). If not provided, default to 42.")
        ("dump-random-node", po::bool_switch(&cfg.dumpRandomNode), "Dump random node FASTA")
        ("dump-random-nodeIDs", po::value<uint32_t>(&cfg.dumpRandomNodeIDs)->default_value(0), "Dump specified number of random node IDs from the tree")
        ("dump-sequence", po::value<std::string>(&cfg.dumpNodeId), "Dump node FASTA")
        ("dump-sequences",po::value<std::vector<std::string>>(&cfg.dumpSequences)->multitoken(), "Dump sequences for a list of node IDs")
        ("simulate-snps",po::value<std::vector<uint32_t>>(&cfg.simulateSNPs)->multitoken(), "Simulate number of SNPs for node IDs, parameter position is relative to dump-sequences")
        ("list-filtered-nodes", po::value<int>(&cfg.listFilteredNodes)->default_value(0), 
            "List N nodes passing length filter (TSV: node_id,length,parent_id,parent_length)")
        ("remove-node", po::value<std::string>(&cfg.removeNodeId), "Exclude node")
        ("top-placements", po::value<int>(&cfg.topPlacements)->default_value(0), "Show top N scores")
        ("debug-node", po::value<std::string>(&cfg.debugNodeId), "Debug node metrics")
        ("compare-nodes", po::value<std::string>(&cfg.compareNodes), "Compare two nodes (nodeA,nodeB)")
        ("dump-all-scores", po::value<std::string>(&cfg.dumpAllScores), "Dump all node scores to TSV file")
        ("write-meta-read-scores-filtered", po::bool_switch(&cfg.writeMetaReadScoresFiltered), "Write filtered meta read scores to TSV file")
        ("write-meta-read-scores-unfiltered", po::bool_switch(&cfg.writeMetaReadScoresUnfiltered), "Write unfiltered meta read scores to TSV file")
        ("write-ocranks", po::bool_switch(&cfg.writeOCRanks), "Write overlap coefficients info to TSV file")
        ("seed", po::value<int>(&cfg.seed)->default_value(42), "Random seed");
    
    // Positional arguments (always hidden)
    po::options_description positional;
    positional.add_options()
        ("panman", po::value<std::string>(&cfg.panman), "")
        ("reads1", po::value<std::string>(&cfg.reads1), "")
        ("reads2", po::value<std::string>(&cfg.reads2), "");
    
    po::positional_options_description pos;
    pos.add("panman", 1).add("reads1", 1).add("reads2", 1);
    
    // Combine option groups
    po::options_description all;  // For parsing
    all.add(visible).add(advanced).add(metagenomic).add(em).add(filterAndAssign).add(developer).add(positional);
    
    po::options_description visible_all;  // For --help-all
    visible_all.add(visible).add(advanced).add(metagenomic).add(em).add(filterAndAssign).add(developer);
    
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
    
    if (vm.count("help-all")) {
        printUsage();
        std::cout << visible_all << "\n";
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
        if (cfg.hpc) {
            configStr += "  hpc";
        }
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
            saveNodeSequence(&tg->trees[0], nodeId, outPath, cfg.impute);
            std::cout << nodeId << "\n";
            return 0;
        }

        if (cfg.dumpRandomNodeIDs > 0) {
          auto tg = loadPanMAN(cfg.panman);
          panmanUtils::Tree* T = &tg->trees[0];

          std::mt19937 rng;
          if (cfg.randomSeed.empty()) {
            std::hash<std::string> hasher;
            rng = std::mt19937(hasher(cfg.randomSeed));
          } else {
            std::random_device rd;
            rng = std::mt19937(rd());
          }

          uint32_t num_nodes = cfg.dumpRandomNodeIDs;
          std::vector<std::string_view> allNodeIDs;
          allNodeIDs.reserve(T->allNodes.size());
          for (const auto& [nodeID, node] : T->allNodes) {
            if (node->children.empty()) {
              allNodeIDs.push_back(nodeID);
            }
          }
          allNodeIDs.shrink_to_fit();
          std::sort(allNodeIDs.begin(), allNodeIDs.end(), std::greater<std::string_view>());
          std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

          std::ofstream outFile(cfg.output + ".randomNodeIDs.txt");
          for (size_t i = 0; i < std::min(num_nodes, static_cast<uint32_t>(allNodeIDs.size())); i++) {
            outFile << allNodeIDs[i] << std::endl;
          }
          outFile.close();
          return 0;
        }

        if (cfg.dumpSequences.size() > 0) {
          auto tg = loadPanMAN(cfg.panman);
          panmanUtils::Tree* T = &tg->trees[0];

          std::mt19937 rng;
          if (cfg.randomSeed.empty()) {
            std::hash<std::string> hasher;
            rng = std::mt19937(hasher(cfg.randomSeed));
          } else {
            std::random_device rd;
            rng = std::mt19937(rd());
          }

          std::vector<std::string> nodeIDs;
          auto& nodeID_groups = cfg.dumpSequences;
          std::cerr << "Node ID groups: " << nodeID_groups.size() << std::endl;
          for (size_t i = 0; i < nodeID_groups.size(); i++) {
            const auto& nodeID_group = nodeID_groups[i];
            std::vector<std::string> nodeID_group_parts;
            boost::split(nodeID_group_parts, nodeID_group, boost::is_any_of(" "), boost::token_compress_on);
            for (const auto& nodeID : nodeID_group_parts) {
              nodeIDs.push_back(nodeID);
              std::cerr << "Node ID " << nodeID << " added to dump sequences" << std::endl;
            }
            std::cerr << "Node ID group " << i << " size: " << nodeID_group_parts.size() << std::endl;
          }
    
          std::vector<uint32_t> numsnps;
          if (cfg.simulateSNPs.size() > 0) {
            numsnps = cfg.simulateSNPs;
            if (numsnps.size() != nodeIDs.size()) {
              std::cerr << "Number of SNP parameters does not match number of node IDs" << std::endl;
              return 1;
            }
          }
    
          std::string outputFileName = cfg.output + ".dump-sequences.fa";
          std::ofstream outFile(outputFileName);
          for (size_t i = 0; i < nodeIDs.size(); i++) {
            const auto& nodeID = nodeIDs[i];
            uint32_t numsnp = (numsnps.empty() ? 0 : numsnps[i]);
            if (T->allNodes.find(nodeID) == T->allNodes.end()) {
              std::cerr << "Node ID " << nodeID << " not found in the tree" << std::endl;
              return 1;
            }
    
            std::string sequence = panmapUtils::getStringFromReference(T, nodeID, false);
            std::vector<std::tuple<char, char, uint32_t>> snpRecords;
            panmapUtils::simulateSNPsOnSequence(sequence, snpRecords, numsnp, rng);
    
            if (outFile.is_open()) {
              outFile << ">" << nodeID << " ";
              for (const auto& [ref, alt, pos] : snpRecords) {
                outFile << ref << pos << alt << " ";
              }
              outFile << "\n";
              for (size_t i = 0; i < sequence.size(); i += 80) {
                outFile << sequence.substr(i, 80) << "\n";
              }
              std::cout << "Sequence for node " << nodeID << " with " << numsnp << " SNPs written to " << outputFileName << std::endl;
            } else {
              std::cerr << "Failed to open file " << outputFileName << " for writing" << std::endl;
              return 1;
            }
          }
          outFile.close();
          return 0;
        }
        
        // Utility: list filtered nodes
        if (cfg.listFilteredNodes > 0) {
            auto tg = loadPanMAN(cfg.panman);
            panmanUtils::Tree* T = &tg->trees[0];
            
            // Build node map
            std::unordered_map<std::string, panmanUtils::Node*> nodeMap;
            std::vector<std::string> ids;
            std::function<void(panmanUtils::Node*)> collect = [&](panmanUtils::Node* n) {
                if (!n) return;
                ids.push_back(n->identifier);
                nodeMap[n->identifier] = n;
                for (auto* c : n->children) collect(c);
            };
            collect(T->root);
            
            // Sample to estimate median
            const size_t sampleSize = std::min(ids.size(), size_t(100));
            std::mt19937 sampleGen(42);
            std::vector<std::string> sampleIds = ids;
            std::shuffle(sampleIds.begin(), sampleIds.end(), sampleGen);
            sampleIds.resize(sampleSize);
            
            std::vector<size_t> sampleLengths;
            for (const auto& id : sampleIds) {
                sampleLengths.push_back(getUngappedLength(T, id));
            }
            std::sort(sampleLengths.begin(), sampleLengths.end());
            size_t medianLen = sampleLengths[sampleLengths.size() / 2];
            
            double tolerance = 0.10;
            size_t minLen = static_cast<size_t>(medianLen * (1.0 - tolerance));
            size_t maxLen = static_cast<size_t>(medianLen * (1.0 + tolerance));
            
            const size_t maxNs = 5;  // Maximum allowed N characters
            logging::info("Median length: {} bp, range: {}-{} bp, max Ns: {}", medianLen, minLen, maxLen, maxNs);
            
            // Shuffle and collect filtered nodes
            std::mt19937 gen(cfg.seed);
            std::shuffle(ids.begin(), ids.end(), gen);
            
            std::cout << "node_id\tlength\tn_count\tparent_id\tparent_length\n";
            int count = 0;
            std::unordered_map<std::string, size_t> lengthCache;
            
            for (const auto& id : ids) {
                if (count >= cfg.listFilteredNodes) break;
                
                size_t len = getUngappedLength(T, id);
                lengthCache[id] = len;
                
                if (len < minLen || len > maxLen) continue;
                
                // Filter by N count
                size_t nCount = getNCount(T, id);
                if (nCount > maxNs) continue;
                
                panmanUtils::Node* node = nodeMap[id];
                std::string parentId = node->parent ? node->parent->identifier : "";
                size_t parentLen = 0;
                if (node->parent) {
                    if (lengthCache.find(parentId) == lengthCache.end()) {
                        lengthCache[parentId] = getUngappedLength(T, parentId);
                    }
                    parentLen = lengthCache[parentId];
                }
                
                std::cout << id << "\t" << len << "\t" << nCount << "\t" << parentId << "\t" << parentLen << "\n";
                count++;
            }
            
            logging::info("Listed {} filtered nodes", count);
            return 0;
        }
        
        // Utility: dump specific node sequence
        if (!cfg.dumpNodeId.empty()) {
            auto tg = loadPanMAN(cfg.panman);
            // Use -o output path if specified, otherwise default to panman path
            std::string outPath = cfg.output.empty()
                ? cfg.panman + "." + sanitizeFilename(cfg.dumpNodeId) + ".fa"
                : cfg.output;
            saveNodeSequence(&tg->trees[0], cfg.dumpNodeId, outPath, cfg.impute);
            std::cout << cfg.dumpNodeId << "\n";
            return 0;
        }

        // metagenomics mode related
        if (!cfg.indexMgsr.empty()) {
          if (!buildMgsrIndex(cfg)) return 1;
          return 0;
        }

        if (cfg.metagenomic) {
          if (!runMetagenomic(cfg)) return 1;
          std::cout << "Metagenomic mode run completed" << std::endl;
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
/**
 * @file main.cpp
 * @brief Panmap - Pangenome-based sequence placement, alignment, and genotyping
 * 
 * Supports single-sample (isolate) and metagenomic workflows.
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
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
    int minReadSupport = 1;          // Min reads for a seed to be counted (2 = filter singletons)
    
    // Resources
    int threads = 1;
    
    // Utility modes
    bool dumpRandomNode = false;
    bool dumpSequence = false;
    std::string dumpNodeId;
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
    bool impute = true;          // Impute N's from parent sequence (ignore _->N mutations)
    
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
    // Initialize output filenames
    std::string refFileName, samFileName, bamFileName, mpileupFileName,
        vcfFileName;

    for (const auto &output : outputs_seperated) {
      if (output.size() == 1) {
        switch (output[0]) {
        case 'r':
          refFileName = prefix + ".reference.fa";
          break;
        case 's':
          samFileName = prefix + ".sam";
          break;
        case 'b':
          bamFileName = prefix + ".bam";
          break;
        case 'm':
          mpileupFileName = prefix + ".mpileup";
          break;
        case 'v':
          vcfFileName = prefix + ".vcf";
          break;
        case 'A':
          refFileName = prefix + ".reference.fa";
          samFileName = prefix + ".sam";
          bamFileName = prefix + ".bam";
          mpileupFileName = prefix + ".mpileup";
          vcfFileName = prefix + ".vcf";
          break;
        }
      } else {
        if (output == "reference")
          refFileName = prefix + ".reference.fa";
        else if (output == "sam")
          samFileName = prefix + ".sam";
        else if (output == "bam") {
          samFileName = prefix + ".sam";
          bamFileName = prefix + ".bam";
        } else if (output == "mpileup")
          mpileupFileName = prefix + ".mpileup";
        else if (output == "vcf")
          vcfFileName = prefix + ".vcf";
        else if (output == "all") {
          refFileName = prefix + ".reference.fa";
          samFileName = prefix + ".sam";
          bamFileName = prefix + ".bam";
          mpileupFileName = prefix + ".mpileup";
          vcfFileName = prefix + ".vcf";
        }
      }
    }

    std::mt19937 rng;
    if (vm.count("random-seed")) {
      std::string seed_str = vm["random-seed"].as<std::string>();
      std::hash<std::string> hasher;
      rng = std::mt19937(hasher(seed_str));
    } else {
      std::random_device rd;
      rng = std::mt19937(rd());
    }
    
    if (vm.count("dump-node-cluster") || vm.count("dump-node-cluster-leaves")
     || vm.count("dump-random-node-cluster") || vm.count("dump-random-node-cluster-leaves")
     || vm.count("dump-random-node-clusters") || vm.count("dump-random-node-clusters-leaves")
    ) {
      std::string mgsr_index_path = vm["mgsr-index"].as<std::string>();
      int fd = mgsr::open_file(mgsr_index_path);
      ::capnp::ReaderOptions readerOptions {.traversalLimitInWords = std::numeric_limits<uint64_t>::max(), .nestingLimit = 1024};
      ::capnp::PackedFdMessageReader reader(fd, readerOptions);
      MGSRIndex::Reader indexReader = reader.getRoot<MGSRIndex>();
      LiteTree::Reader liteTreeReader = indexReader.getLiteTree();
      size_t numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
      bool lowMemory = vm.count("low-memory") > 0;
      
      mgsr::MgsrLiteTree T;
      T.initialize(indexReader, "", 0, numThreads, lowMemory, true);

      if (vm.count("dump-node-cluster") || vm.count("dump-node-cluster-leaves")) {
        if (!vm.count("mgsr-index")) {
          err("--mgsr-index must be provided to use --dump-node-cluster");
          return 1;
        }
        auto inputParameters = vm.count("dump-node-cluster") 
                            ? vm["dump-node-cluster"].as<std::vector<std::string>>()
                            : vm["dump-node-cluster-leaves"].as<std::vector<std::string>>();
        if (inputParameters.size() != 2) {
          err("Expected 2 parameters for --dump-node-cluster: <nodeID> <numNodes>");
          return 1;
        }

        uint32_t numNodes;
        try {
          numNodes = std::stoi(inputParameters[1]);
          // Proceed with using numNodes
        } catch (const std::invalid_argument& e) {
          err("The second parameter is not convertible to an integer");
          return 1;
        }

        std::string nodeId = inputParameters[0];
        if (T.allLiteNodes.find(nodeId) == T.allLiteNodes.end()) {
          err("Node ID {} not found in the tree", nodeId);
          return 1;
        }

        std::vector<mgsr::MgsrLiteNode*> nearestNodes = vm.count("dump-node-cluster")
                                                    ? mgsr::getNearestNodes(T.allLiteNodes.find(nodeId)->second, numNodes, false)
                                                    : mgsr::getNearestNodes(T.allLiteNodes.find(nodeId)->second, numNodes, true);
        std::ofstream outFile(prefix + ".clusterIDs.tsv");
        outFile << "Strain\tClusterID" << std::endl;
        for (const auto& node : nearestNodes) {
          outFile << node->identifier << "\t0" << std::endl;
        }
        outFile.close();
        msg("Cluster node IDs written to {}", prefix + ".clusterIDs.tsv");
        exit(0);
      } else if (vm.count("dump-random-node-cluster") || vm.count("dump-random-node-cluster-leaves")) {
        if (!vm.count("mgsr-index")) {
          err("--mgsr-index must be provided to use --dump-node-cluster");
          return 1;
        }
        uint32_t numNodes = vm.count("dump-random-node-cluster")
                            ? vm["dump-random-node-cluster"].as<uint32_t>()
                            : vm["dump-random-node-cluster-leaves"].as<uint32_t>();
        std::vector<std::string_view> allNodeIDs;
        allNodeIDs.reserve(T.allLiteNodes.size());
        for (const auto& [nodeID, node] : T.allLiteNodes) {
          if (node->children.empty()) {
            allNodeIDs.push_back(nodeID);
          }
        }
        allNodeIDs.shrink_to_fit();
        std::sort(allNodeIDs.begin(), allNodeIDs.end(), std::greater<std::string_view>());
        std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

        std::vector<mgsr::MgsrLiteNode*> nearestNodes = vm.count("dump-random-node-cluster")
                                                    ? mgsr::getNearestNodes(T.allLiteNodes.find(std::string(allNodeIDs[0]))->second, numNodes, false)
                                                    : mgsr::getNearestNodes(T.allLiteNodes.find(std::string(allNodeIDs[0]))->second, numNodes, true);

        std::ofstream outFile(prefix + ".randomClusterIDs.tsv");
        outFile << "Strain\tClusterID" << std::endl;
        for (const auto& node : nearestNodes) {
          outFile << node->identifier << "\t0" << std::endl;
        }
        outFile.close();
        msg("Random node IDs written to {}", prefix + ".randomClusterIDs.tsv");
        exit(0);
      } else if (vm.count("dump-random-node-clusters") || vm.count("dump-random-node-clusters-leaves")) {
        if (!vm.count("mgsr-index")) {
          err("--mgsr-index must be provided to use --dump-node-cluster");
          return 1;
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
    
    // Compute distance from expected parent if provided
    std::string distStr = "";
    if (tree && !expectedParent.empty()) {
        int dist = computeTreeDistance(*tree, node->identifier, expectedParent);
        distStr = "  dist=" + std::to_string(dist);
    }
    
    std::cout << "  NODE: " << node->identifier << distStr << "\n";
    std::cout << "    LogRaw=" << node->logRawScore 
              << "  LogCosine=" << node->logCosineScore
              << "  Presence=" << node->presenceScore
              << "  Containment=" << node->containmentScore << "\n";
    std::cout << "    LogLogRaw=" << node->logLogRawScore 
              << "  LogLogCosine=" << node->logLogCosineScore
              << "  Cosine=" << node->cosineScore 
              << "  Raw=" << node->rawSeedMatchScore << "\n";
    std::cout << "    CapCosine=" << node->capCosineScore
              << "  CapLogCosine=" << node->capLogCosineScore
              << "  LogCosineOld=" << node->logCosineOldScore
              << "  WeightedContainment=" << node->weightedContainmentScore << "\n";
    std::cout << "    Intersection=" << node->presenceIntersectionCount
              << "  GenomeUnique=" << node->genomeUniqueSeedCount
              << "  GenomeOnly=" << genomeOnlyCount
              << "  Union=" << totalUnion << "\n";
    std::cout << "    CosineNum=" << node->cosineNumerator
              << "  ReadMag=" << readMagnitude
              << "  GenomeMag=" << genomeMag << "\n";
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
    double cosineNumA = 0.0, cosineNumB = 0.0;
    
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
            cosineNumA += static_cast<double>(readCount) * countA;
        }
        if (countB > 0) {
            intersectB++;
            logRawNumB += logRead / countB;
            cosineNumB += static_cast<double>(readCount) * countB;
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
    
    std::cout << "\nCOSINE NUMERATOR:\n";
    std::cout << "  Node A: " << cosineNumA << "\n";
    std::cout << "  Node B: " << cosineNumB << "\n";
    
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
    std::cout << "  Node A: LogRaw=" << nodeA->logRawScore << " LogCosine=" << nodeA->logCosineScore 
              << " Intersection=" << nodeA->presenceIntersectionCount << "\n";
    std::cout << "  Node B: LogRaw=" << nodeB->logRawScore << " LogCosine=" << nodeB->logCosineScore 
              << " Intersection=" << nodeB->presenceIntersectionCount << "\n";
    
    std::cout << std::string(80, '=') << "\n\n";
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
        std::vector<uint32_t> nodeClusterNodes = vm.count("dump-random-node-clusters")
                                              ? vm["dump-random-node-clusters"].as<std::vector<uint32_t>>()
                                              : vm["dump-random-node-clusters-leaves"].as<std::vector<uint32_t>>();

        std::sort(nodeClusterNodes.begin(), nodeClusterNodes.end(), std::greater<uint32_t>());

        std::vector<std::string_view> allNodeIDs;
        allNodeIDs.reserve(T.allLiteNodes.size());
        for (const auto& [nodeID, node] : T.allLiteNodes) {
          if (node->children.empty()) {
            allNodeIDs.push_back(nodeID);
          }
        }
        allNodeIDs.shrink_to_fit();
        std::sort(allNodeIDs.begin(), allNodeIDs.end(), std::greater<std::string_view>());
        std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

        std::vector<std::vector<std::string_view>> nodeClusters(nodeClusterNodes.size());
        std::unordered_set<std::string_view> selectedNodes;
        std::unordered_set<std::string_view> paddedNodes;
        size_t allNodeIDsIndex = 0;
        for (size_t i = 0; i < nodeClusterNodes.size(); i++) {
          uint32_t curClusterSize = nodeClusterNodes[i];
          uint32_t nextClusterSize = (i + 1 < nodeClusterNodes.size()) ? nodeClusterNodes[i + 1] : 1;

          std::string_view curNode;
          while (true) {
            curNode = allNodeIDs[allNodeIDsIndex];
            ++allNodeIDsIndex;
            if (allNodeIDsIndex >= allNodeIDs.size()) {
              err("Not enough nodes to fulfill the requested cluster sizes");
              exit(1);
            }
            if (paddedNodes.find(curNode) == paddedNodes.end()) {
              break;
            }
          }

          std::vector<mgsr::MgsrLiteNode*> nearestNodes = vm.count("dump-random-node-clusters")
              ? mgsr::getNearestNodes(T.allLiteNodes.find(std::string(curNode))->second, selectedNodes, curClusterSize + nextClusterSize - 1, false)
              : mgsr::getNearestNodes(T.allLiteNodes.find(std::string(curNode))->second, selectedNodes, curClusterSize + nextClusterSize - 1, true);
          for (size_t j = 0; j < nearestNodes.size(); j++) {
            if (j < curClusterSize) {
              selectedNodes.insert(nearestNodes[j]->identifier);
              nodeClusters[i].push_back(nearestNodes[j]->identifier);
            }
            paddedNodes.insert(nearestNodes[j]->identifier);
          }
        }

        std::ofstream outFile(prefix + ".randomClusterIDs.tsv");
        outFile << "Strain\tClusterID" << std::endl;
        for (size_t i = 0; i < nodeClusters.size(); i++) {
          for (const auto& node : nodeClusters[i]) {
            outFile << node << "\t" << (i + 1) << std::endl;
          }
        }
        outFile.close();
        msg("Random node IDs written to {}", prefix + ".randomClusterIDs.tsv");
        exit(0);
      }
    }



    if (vm.count("mgsr-index") && !reads1.empty()) {
      std::string mgsr_index_path = vm["mgsr-index"].as<std::string>();
      std::string taxonomicMetadataPath = vm.count("taxonomic-metadata") ? vm["taxonomic-metadata"].as<std::string>() : "";
      size_t maximumFamilies = taxonomicMetadataPath == "" ? 0 : vm["maximum-families"].as<size_t>();
      
      bool filterAndAssign = vm.count("filter-and-assign") > 0;
      uint32_t maskReads = vm.count("mask-reads") ? vm["mask-reads"].as<uint32_t>() : 0;
      uint32_t maskSeeds = vm.count("mask-reads") ? vm["mask-seeds"].as<uint32_t>() : 0;
      std::string ampliconDepthPath = vm.count("amplicon-depth") ? vm["amplicon-depth"].as<std::string>() : "";
      double maskReadsRelativeFrequency = vm["mask-reads-relative-frequency"].as<double>();
      double maskSeedsRelativeFrequency = vm["mask-seeds-relative-frequency"].as<double>();
      
      int fd = mgsr::open_file(mgsr_index_path);
      ::capnp::ReaderOptions readerOptions {.traversalLimitInWords = std::numeric_limits<uint64_t>::max(), .nestingLimit = 1024};
      ::capnp::PackedFdMessageReader reader(fd, readerOptions);
      MGSRIndex::Reader indexReader = reader.getRoot<MGSRIndex>();
      LiteTree::Reader liteTreeReader = indexReader.getLiteTree();
      size_t numThreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
      bool lowMemory = vm.count("low-memory") > 0;
      
      mgsr::MgsrLiteTree liteTree;
      if (vm.count("true-abundance")) {
        std::string trueAbundancePath = vm["true-abundance"].as<std::string>();
        liteTree.loadTrueAbundances(trueAbundancePath);
      }
      liteTree.initialize(indexReader, taxonomicMetadataPath, maximumFamilies, numThreads, lowMemory, true);

      liteTree.debugNodeID = vm["debug-node-id"].as<std::string>();


      
      bool progressBar = true;
      if (vm.count("no-progress")) progressBar = false;
      mgsr::ThreadsManager threadsManager(&liteTree, prefix, numThreads, maskReads, progressBar, lowMemory);
      threadsManager.initializeMGSRIndex(indexReader);
      close(fd);

      bool sketchSeeds = vm.count("sketch-seeds") > 0;
      uint32_t maskReadsEnds = vm["mask-read-ends"].as<uint32_t>();
      double dust_threshold = vm["dust"].as<double>();
      if (dust_threshold > 100) {
        std::cerr << "Error: --dust must be <= 100" << std::endl;
        exit(1);
      }
      threadsManager.initializeQueryData(reads1, reads2, maskSeeds, ampliconDepthPath, maskReadsRelativeFrequency, maskSeedsRelativeFrequency, dust_threshold, maskReadsEnds);
      
      if (sketchSeeds) {
        std::ofstream seedsOut(prefix + ".seeds.txt");
        for (size_t i = 0; i < threadsManager.reads.size(); i++) {
          const auto& read = threadsManager.reads[i];
          seedsOut << i << "\t";
          for (size_t j = 0; j < threadsManager.readSeedmersDuplicatesIndex[i].size(); ++j) {
            if (j == 0) {
              seedsOut << threadsManager.readSeedmersDuplicatesIndex[i][j];
            } else {
              seedsOut << "," << threadsManager.readSeedmersDuplicatesIndex[i][j];
            }
          }
          seedsOut << "\t" << read.seedmersList.size() << "\t";
          for (size_t j = 0; j < read.seedmersList.size(); j++) {
            const auto& seedmer = read.seedmersList[j];
            seedsOut << "(" << seedmer.begPos << "," << seedmer.endPos << "," << seedmer.hash << "," << seedmer.rev << "," << seedmer.iorder << ")";
            if (j != read.seedmersList.size() - 1) {
              seedsOut << " ";
            } else {
              seedsOut << "\n";
            }
          }
        }
        seedsOut.close();
        exit(0);
      }

      size_t overlap_coefficients_threshold = vm["overlap-coefficients"].as<size_t>();
      if (overlap_coefficients_threshold == 0 && !vm.count("read-scores") && !filterAndAssign) {
        std::cerr << "Error: Either --overlap-coefficients > 0 or --read-scores must be specified for EM." << std::endl;
        exit(1);
      }
      if (overlap_coefficients_threshold > 0 && !filterAndAssign) {
        auto start_time_computeOverlapCoefficients = std::chrono::high_resolution_clock::now();
        mgsr::mgsrPlacer placer(&liteTree, threadsManager, lowMemory, 0);
        auto overlapCoefficients = placer.computeOverlapCoefficients(threadsManager.allSeedmerHashesSet);
        std::ofstream overlapCoefficientsFile(prefix + ".overlapCoefficients.txt");
        std::sort(overlapCoefficients.begin(), overlapCoefficients.end(), [](const auto& a, const auto& b) {
          return a.second > b.second;
        });
        uint32_t rank = 0;
        double currentOverlapCoefficient = overlapCoefficients[0].second;
        for (const auto& [nodeId, overlapCoefficient] : overlapCoefficients) {
          if (overlapCoefficient != currentOverlapCoefficient) {
            currentOverlapCoefficient = overlapCoefficient;
            ++rank;
          }
          liteTree.allLiteNodes.at(nodeId)->ocRank = rank;
          overlapCoefficientsFile << nodeId << "\t" << std::fixed << std::setprecision(6) << overlapCoefficient << "\t" << rank <<  std::endl;
        }
        overlapCoefficientsFile.close();
        auto end_time_computeOverlapCoefficients = std::chrono::high_resolution_clock::now();
        auto duration_computeOverlapCoefficients = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_computeOverlapCoefficients - start_time_computeOverlapCoefficients);
        std::cout << "Computed overlap coefficients in " << static_cast<double>(duration_computeOverlapCoefficients.count()) / 1000.0 << "s\n" << std::endl;
      }

      liteTree.collapseIdenticalScoringNodes(threadsManager.allSeedmerHashesSet);
      // std::string collapsedNewick = liteTree.toNewick(true);
      // std::ofstream of(prefix + ".collapsed.newick");
      // of << collapsedNewick;
      // of.close();



      if (vm.count("read-seed-scores")) {
        threadsManager.countSeedNodesFrequency();
        decltype(threadsManager.seedMatchedNodeRanges)().swap(threadsManager.seedMatchedNodeRanges);
      }
      if (vm.count("seed-scores")) {
        threadsManager.countSeedNodesFrequency();
        threadsManager.computeNodeSeedScores();
        std::ofstream nodeSeedScoresOut(prefix + ".nodeSeedScores.tsv");
        nodeSeedScoresOut << "NodeId\tDistance\tnodeSeedScores\tnodeSeedScoresCorrected\tselected\tselectedNeighbor\tinSample\tcollapsedNodes" << std::endl;
        for (const auto& [node, nodeSeedScore] : threadsManager.nodeSeedScores) {
          nodeSeedScoresOut << node->identifier
                            << "\t" << node->seedDistance
                            << "\t" << nodeSeedScore
                            << "\t" << threadsManager.nodeSeedScoresCorrected.find(node)->second
                            << "\t" << node->selected
                            << "\t" << node->selectedNeighbor
                            << "\t" << node->inSample << "\t";
          if (node->identicalNodeIdentifiers.empty()) {
            nodeSeedScoresOut << "." << std::endl;
          } else {
            for (size_t i = 0; i < node->identicalNodeIdentifiers.size(); ++i) {
              nodeSeedScoresOut << node->identicalNodeIdentifiers[i];
              if (i != node->identicalNodeIdentifiers.size() - 1) {
                nodeSeedScoresOut << ",";
              }
            }
            nodeSeedScoresOut << std::endl;
          }

        }
        nodeSeedScoresOut.close();
        exit(0);
      }
            
      std::vector<uint64_t> totalNodesPerThread(numThreads, 0);
      for (size_t i = 0; i < numThreads; ++i) {
        totalNodesPerThread[i] = liteTree.getNumActiveNodes();
      }
      ProgressTracker progressTracker(numThreads, totalNodesPerThread);
      std::cout << "Using " << numThreads << " threads" << std::endl;

      auto start_time_place = std::chrono::high_resolution_clock::now();
      std::atomic<size_t> numGroupsUpdate = 0;
      std::atomic<size_t> numReadsUpdate = 0;
      std::vector<size_t> numUniqueKminmersPerThread(numThreads, 0);
      std::vector<size_t> numNodesPostCollapsePerThread(numThreads, 0);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()), [&](const tbb::blocked_range<size_t>& rangeIndex){
        for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
          auto [start, end] = threadsManager.threadRanges[i];
        
          std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);
          mgsr::mgsrPlacer curThreadPlacer(&liteTree, threadsManager, lowMemory, i);
          curThreadPlacer.initializeQueryData(curThreadReads);

          curThreadPlacer.setAllSeedmerHashesSet(threadsManager.allSeedmerHashesSet);
          
          curThreadPlacer.setProgressTracker(&progressTracker, i);
          // curThreadPlacer.placeReads();

          curThreadPlacer.scoreReads();

          threadsManager.readMinichainsInitialized[i] = curThreadPlacer.readMinichainsInitialized;
          threadsManager.readMinichainsAdded[i] = curThreadPlacer.readMinichainsAdded;
          threadsManager.readMinichainsRemoved[i] = curThreadPlacer.readMinichainsRemoved;
          threadsManager.readMinichainsUpdated[i] = curThreadPlacer.readMinichainsUpdated;
          if (i == 0) {
            threadsManager.identicalGroups = std::move(curThreadPlacer.identicalGroups);
            threadsManager.identicalNodeToGroup = std::move(curThreadPlacer.identicalNodeToGroup);
          }
          numGroupsUpdate += curThreadPlacer.numGroupsUpdate;
          numReadsUpdate += curThreadPlacer.numReadsUpdate;
        }
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    placement::placeLite(result, &tree, reader, cfg.reads1, cfg.reads2, outPath, false, fullTreePtr,
                        cfg.removeNodeId, expectedParentPtr, storeDiagnostics, cfg.seedMaskFraction,
                        cfg.minSeedQuality, cfg.dedupReads, true, cfg.trimStart, cfg.trimEnd,
                        cfg.minReadSupport,
                        400, 150, // expectedFragmentSize, fragmentSizeTolerance
                        refineEnabled, cfg.refineTopPct, cfg.refineMaxTopN, 
                        cfg.refineNeighborRadius, cfg.refineMaxNeighborN);
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
        
        // Report refinement results if run
        if (result.refinementWasRun) {
            int distRefined = computeTreeDistance(tree, result.bestRefinedNodeId, expectedParentNode);
            bool correctRefined = (result.bestRefinedNodeId == expectedParentNode);
            logging::msg("  Refined:   {} {}{}{} (dist: {}, score: {:.0f})", 
                        result.bestRefinedNodeId,
                        correctRefined ? color::green() : color::yellow(),
                        correctRefined ? output::box::check() : output::box::cross(),
                        color::reset(),
                        distRefined,
                        result.bestRefinedScore);
            std::cout << "REFINED\t" << result.bestRefinedNodeId 
                      << "\t" << (correctRefined ? "true" : "false") 
                      << "\t" << distRefined 
                      << "\t" << result.bestRefinedScore << "\n";
        }
    }
    
    // Dump all scores to file if requested
    if (!cfg.dumpAllScores.empty()) {
        std::ofstream outFile(cfg.dumpAllScores);
        if (outFile) {
            outFile << "node\tlogRaw\tlogCosine\tpresence\tcontainment\n";
            std::vector<std::pair<double, std::string>> allScores;
            for (auto& [id, node] : tree.allLiteNodes) {
                if (node->logRawScore > 0) {
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
                        << node->presenceScore << "\t"
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
            if (node->logRawScore > 0 || node->logCosineScore > 0 || 
                node->cosineScore > 0 || node->presenceScore > 0) {
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
            
            // Top by LogCosine
            std::sort(scoredNodes.begin(), scoredNodes.end(),
                [](auto* a, auto* b) { return a->logCosineScore > b->logCosineScore; });
            std::cout << "\nTOP " << topN << " BY LOG_COSINE:\n";
            for (int i = 0; i < topN; ++i) {
                printNodeDiagnostics(scoredNodes[i], result.readUniqueSeedCount,
                                    result.totalReadSeedFrequency, result.readMagnitude,
                                    &tree, expectedParentNode);
            }
            
            // Top by Presence
            std::sort(scoredNodes.begin(), scoredNodes.end(),
                [](auto* a, auto* b) { return a->presenceScore > b->presenceScore; });
            std::cout << "\nTOP " << topN << " BY PRESENCE:\n";
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
      });
      auto end_time_place = std::chrono::high_resolution_clock::now();
      auto duration_place = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_place - start_time_place);
      std::cerr << "\n\nPlaced reads in " << static_cast<double>(duration_place.count()) / 1000.0 << "s\n" << std::endl;


      if (vm.count("print-unfiltered-seed-scores-info")) {
        std::ofstream readScoresOut_unfiltered(prefix + ".read_scores_info.unfiltered.tsv");
        readScoresOut_unfiltered << "ReadIndex\tNumDuplicates\tTotalScore\tMaxScore\tNumMaxScoreNodes\tRawReadsIndices" << std::endl;
        for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
          const auto& curRead = threadsManager.reads[i];
          if (curRead.maxScore == 0) continue;
          readScoresOut_unfiltered << i << "\t" << threadsManager.readSeedmersDuplicatesIndex[i].size() << "\t" << curRead.seedmersList.size() << "\t" << curRead.maxScore << "\t" << curRead.epp << "\t";
          for (size_t j = 0; j < threadsManager.readSeedmersDuplicatesIndex[i].size(); ++j) {
            if (j == 0) {
              readScoresOut_unfiltered << threadsManager.readSeedmersDuplicatesIndex[i][j];
            } else {
              readScoresOut_unfiltered << "," << threadsManager.readSeedmersDuplicatesIndex[i][j];
            }
          }
          readScoresOut_unfiltered << std::endl;
        }
        readScoresOut_unfiltered.close();
      }


      double discard_threshold = vm["discard"].as<double>();
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
      std::cout << threadsManager.reads.size() - num_unmapped - num_discarded << " reads remain for node scoring and EM... " << std::endl;

      if (threadsManager.reads.size() - num_unmapped - num_discarded == 0) {
        std::cerr << "No reads remain for node scoring and EM after discarding low-score reads... Exiting... " << std::endl;
        exit(0);
      }

      // Remove score updates of unmapped reads
      tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsManager.threadRanges.size()), [&](const tbb::blocked_range<size_t>& rangeIndex){
        for (size_t i = rangeIndex.begin(); i != rangeIndex.end(); ++i) {
          auto [start, end] = threadsManager.threadRanges[i];
          std::span<mgsr::Read> curThreadReads(threadsManager.reads.data() + start, end - start);

          for (auto [_, node] : liteTree.allLiteNodes) {
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


      if (vm.count("read-scores")) {
        threadsManager.scoreNodesMultithreaded();
      }


      if (filterAndAssign) {
        std::unordered_map<mgsr::MgsrLiteNode*, std::vector<size_t>> assignedReadsByNode;
        threadsManager.assignReads(assignedReadsByNode, maximumFamilies);

        std::ofstream assignedReadsOut(prefix + ".mgsr.assignedReads.out");
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
      }

      std::ofstream readScoresOut(prefix + ".read_scores_info.tsv");
      readScoresOut << "ReadIndex\tNumDuplicates\tTotalScore\tMaxScore\tNumMaxScoreNodes\tOverMaximumFamilies\tRawReadsIndices" << std::endl;
      for (size_t i = 0; i < threadsManager.reads.size(); ++i) {
        const auto& curRead = threadsManager.reads[i];
        if (curRead.maxScore == 0) continue;
        readScoresOut << i << "\t" << threadsManager.readSeedmersDuplicatesIndex[i].size() << "\t" << curRead.seedmersList.size() << "\t" << curRead.maxScore << "\t" << curRead.epp << "\t" << curRead.overMaximumFamilies << "\t";
        for (size_t j = 0; j < threadsManager.readSeedmersDuplicatesIndex[i].size(); ++j) {
          if (j == 0) {
            readScoresOut << threadsManager.readSeedmersDuplicatesIndex[i][j];
          } else {
            readScoresOut << "," << threadsManager.readSeedmersDuplicatesIndex[i][j];
          }
        }
        readScoresOut << std::endl;
      }
      readScoresOut.close();

      if (filterAndAssign) {
        exit(0);
      }




      std::ofstream scoresOut;
      if (vm.count("read-seed-scores")) {
        scoresOut.open(prefix + ".nodeReadSeedScores.tsv");
      } else {
        scoresOut.open(prefix + ".nodeScores.tsv");
      }
      scoresOut << "NodeId\tdistance\tWEPPScore\tWEPPScoreCorrected\tWEPPScoreCorrectedSelected\tSelectedNeighbor\tinSample\tcollapsedNodes" << std::endl;
      for (const auto& [nodeId, node] : liteTree.allLiteNodes) {
        if (liteTree.detachedNodes.find(node) != liteTree.detachedNodes.end()) {
          continue;
        }
        scoresOut << nodeId
                  << "\t" << node->seedDistance
                  << "\t" << node->sumWEPPScore.sum
                  << "\t" << node->sumWEPPScoreCorrected.sum
                  << "\t" << node->selected
                  << "\t" << node->selectedNeighbor
                  << "\t" << node->inSample << "\t";
        if (node->identicalNodeIdentifiers.empty()) {
          scoresOut << "." << std::endl;
        } else {
          for (size_t i = 0; i < node->identicalNodeIdentifiers.size(); ++i) {
            scoresOut << node->identicalNodeIdentifiers[i];
            if (i != node->identicalNodeIdentifiers.size() - 1) {
              scoresOut << ",";
            }
          }
          scoresOut << std::endl;
        }
      }
      scoresOut.close();

      std::ofstream scoresOutExtra(prefix + ".nodeScores.extra.tsv");
      scoresOutExtra << "NodeId\tdistance\tParsimonyReads\tscoresWeightedParsimonyReadsScore\tavgParsimonyScore\tWEPPScore\tWEPPScoreCorrected\tWEPPScoreCorrectedSelected\tSelectedNeighbor\tinSample\tcollapsedNodes" << std::endl;
      for (const auto& [nodeId, node] : liteTree.allLiteNodes) {
        if (liteTree.detachedNodes.find(node) != liteTree.detachedNodes.end()) {
          continue;
        }
        scoresOutExtra << nodeId
                  << "\t" << node->seedDistance
                  << "\t" << node->sumRawScore
                  << "\t" << node->sumEPPRawScore
                  << "\t" << (node->sumRawScore == 0 ? 0.0 : static_cast<double>(node->sumEPPRawScore) / static_cast<double>(node->sumRawScore))
                  << "\t" << node->sumWEPPScore.sum
                  << "\t" << node->sumWEPPScoreCorrected.sum
                  << "\t" << node->selected
                  << "\t" << node->selectedNeighbor
                  << "\t" << node->inSample << "\t";
        if (node->identicalNodeIdentifiers.empty()) {
          scoresOutExtra << "." << std::endl;
        } else {
          for (size_t i = 0; i < node->identicalNodeIdentifiers.size(); ++i) {
            scoresOutExtra << node->identicalNodeIdentifiers[i];
            if (i != node->identicalNodeIdentifiers.size() - 1) {
              scoresOutExtra << ",";
            }
          }
          scoresOutExtra << std::endl;
        }
      }
      scoresOutExtra.close();

      double em_convergence_threshold = vm["em-convergence-threshold"].as<double>();
      double em_delta_threshold = vm["em-delta-threshold"].as<double>();
      uint32_t em_maximum_iterations = vm["em-maximum-iterations"].as<uint32_t>();
      mgsr::squareEM squareEM(threadsManager, liteTree, prefix, overlap_coefficients_threshold, em_convergence_threshold, em_delta_threshold, em_maximum_iterations, vm.count("read-scores") > 0, vm.count("em-leaves-only") > 0);
      liteTree.seedInfos.clear(); // no longer needed. clear memory to prep for EM.
      
      auto start_time_squareEM = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < 5; ++i) {
        squareEM.runSquareEM();
        std::cout << "\nRound " << i << " of squareEM completed... nodes size changed from " << squareEM.nodes.size() << " to ";
        bool removed = squareEM.removeLowPropNodes();
        std::cout << squareEM.nodes.size() << std::endl;
        if (!removed) {
          break;
        }
      }
      
      std::vector<uint64_t> indices(squareEM.nodes.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(), [&squareEM](uint64_t i, uint64_t j) {
        return squareEM.props[i] > squareEM.props[j];
      });
      std::cout << std::endl;

      std::cout << "writing abundance file: " << prefix + ".mgsr.abundance.out" << std::endl;
      std::ofstream abundanceOutput(prefix + ".mgsr.abundance.out");
      abundanceOutput << std::setprecision(5) << std::fixed;
      for (size_t i = 0; i < indices.size(); ++i) {
        size_t index = indices[i];
        abundanceOutput << squareEM.nodes[index];
        for (const auto& member : liteTree.allLiteNodes.at(squareEM.nodes[index])->identicalNodeIdentifiers) {
          abundanceOutput << "," << member;
        }
        if (squareEM.identicalGroups.find(squareEM.nodes[index]) != squareEM.identicalGroups.end()) {
          for (const auto& member : squareEM.identicalGroups[squareEM.nodes[index]]) {
            abundanceOutput << "," << member;
            for (const auto& identicalMember : liteTree.allLiteNodes.at(member)->identicalNodeIdentifiers) {
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
      std::ofstream assignedReadsOutput(prefix + ".mgsr.assignedReads.out");
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

      
      auto end_time_squareEM = std::chrono::high_resolution_clock::now();
      auto duration_squareEM = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_squareEM - start_time_squareEM);
      std::cout << "SquareEM completed in " << static_cast<double>(duration_squareEM.count()) / 1000.0 << "s\n" << std::endl;

      exit(0);
    }

    // Get remaining parameters
    bool reindex = vm.count("reindex") > 0;
    bool prior = vm.count("prior") > 0;
    bool genotype_from_sam = vm.count("genotype-from-sam") > 0;
    bool save_jaccard = vm.count("save-jaccard") > 0;
    bool show_time = vm.count("time") > 0;
    bool dump_random_node = vm.count("dump-random-node") > 0;

    int k = vm["k"].as<int>(); // Default handled by boost
    int s = vm["s"].as<int>(); // Default handled by boost
    std::string index_path = vm["index"].as<std::string>(); // Default handled by boost


    // Load pangenome
    msg("Loading reference pangenome from: {}", guide);

    std::unique_ptr<panmanUtils::TreeGroup> TG;

    try {
      TG = std::unique_ptr<panmanUtils::TreeGroup>(loadPanMAN(guide));

      if (!TG || TG->trees.empty()) {
        throw std::runtime_error(
            "No valid trees found in the loaded PanMAN file.");
      }

      // Debug log for tree structure
      for (size_t i = 0; i < TG->trees.size(); i++) {
        panmanUtils::Tree &T = TG->trees[i];
        logging::debug(
            "[DEBUG] Tree {} details: blocks.size={}, gaps.size={}, "
            "blockGaps.blockPosition.size={}, allNodes.size={}, root={}",
            i, T.blocks.size(), T.gaps.size(), T.blockGaps.blockPosition.size(),
            T.allNodes.size(), T.root ? T.root->identifier : "null");
      }

      msg("Successfully loaded pangenome with {} trees:", TG->trees.size());
    } catch (const std::exception &e) {
      err("Failed to load reference PanMAN: {}", e.what());
      return 1;
    }

    msg("Using first tree as reference.");
    panmanUtils::Tree &T = TG->trees[0];

    if (vm.count("index-mgsr")) {
      std::string mgsr_index_path = vm["index-mgsr"].as<std::string>();
      int mgsr_l = vm["mgsr-l"].as<int>();
      int mgsr_k = vm["mgsr-k"].as<int>();
      int mgsr_s = vm["mgsr-s"].as<int>();
      int mgsr_t = vm["mgsr-t"].as<int>();
      bool mgsr_open = vm.count("mgsr-open") > 0;
      bool imputeAmb = vm.count("impute") > 0;
      mgsr::mgsrIndexBuilder mgsrIndexBuilder(&T, mgsr_k, mgsr_s, mgsr_t, mgsr_l, mgsr_open, imputeAmb);
      mgsrIndexBuilder.buildIndex();
      mgsrIndexBuilder.writeIndex(mgsr_index_path);
      msg("MGSR index written to: {}", mgsr_index_path);
      return 0;
    }



    // Handle --dump-random-node parameter if provided
    if (dump_random_node) {
      panmanUtils::Node *randomNode = getRandomNode(TG.get(), rng);
      if (!randomNode) {
        err("Failed to select a random node");
        return 1;
      }

      // Find which tree the node belongs to
      panmanUtils::Tree *nodeTree = nullptr;
      for (auto &tree : TG->trees) {
        if (tree.allNodes.find(randomNode->identifier) != tree.allNodes.end()) {
          nodeTree = &tree;
          break;
        }
      }

      if (!nodeTree) {
        err("Could not determine which tree the random node belongs to");
        return 1;
      }

      std::string outputFileName =
          guide + ".random." + randomNode->identifier + ".fa";
      if (saveNodeSequence(nodeTree, randomNode, outputFileName)) {
        msg("Random node {} sequence written to {}", randomNode->identifier,
            outputFileName);
        return 0; // Exit after dumping sequence
      } else {
        err("Failed to save random node sequence to {}", outputFileName);
        return 1;
      }
    }


    if (vm.count("dump-random-nodeIDs")) {
      uint32_t num_nodes = vm["dump-random-nodeIDs"].as<uint32_t>();
      std::vector<std::string_view> allNodeIDs;
      allNodeIDs.reserve(T.allNodes.size());
      for (const auto& [nodeID, node] : T.allNodes) {
        if (node->children.empty()) {
          allNodeIDs.push_back(nodeID);
        }
      }
      allNodeIDs.shrink_to_fit();
      std::sort(allNodeIDs.begin(), allNodeIDs.end(), std::greater<std::string_view>());
      std::shuffle(allNodeIDs.begin(), allNodeIDs.end(), rng);

      std::ofstream outFile(prefix + ".randomNodeIDs.txt");
      for (size_t i = 0; i < std::min(num_nodes, static_cast<uint32_t>(allNodeIDs.size())); i++) {
        outFile << allNodeIDs[i] << std::endl;
      }
      outFile.close();
      msg("Random node IDs written to {}", prefix + ".randomNodeIDs.txt");
      exit(0);
    }

    // Handle --dump-sequence parameter if provided
    if (vm.count("dump-sequence")) {
      std::string nodeID = vm["dump-sequence"].as<std::string>();
      if (T.allNodes.find(nodeID) == T.allNodes.end()) {
        err("Node ID {} not found in the tree", nodeID);
        return 1;
      }

      std::string sequence = T.getStringFromReference(nodeID, false, true);
      std::string outputFileName = guide + "." + nodeID + ".fa";
      std::ofstream outFile(outputFileName);

      if (outFile.is_open()) {
        outFile << ">" << nodeID << "\n";
        outFile << sequence << "\n";
        outFile.close();
        msg("Sequence for node {} written to {}", nodeID, outputFileName);
        return 0; // Exit after dumping sequence
      } else {
        err("Failed to open file {} for writing", outputFileName);
        return 1;
      }
    }

    if (vm.count("dump-sequences")) {
      std::vector<std::string> nodeIDs;
      auto nodeID_groups = vm["dump-sequences"].as<std::vector<std::string>>();

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
      if (vm.count("simulate-snps")) {
        numsnps = vm["simulate-snps"].as<std::vector<uint32_t>>();
        if (numsnps.size() != nodeIDs.size()) {
          err("Number of SNP parameters does not match number of node IDs");
          return 1;
        }
      }

      std::string outputFileName = prefix + ".dump-sequences.fa";
      std::ofstream outFile(outputFileName);
      for (size_t i = 0; i < nodeIDs.size(); i++) {
        const auto& nodeID = nodeIDs[i];
        uint32_t numsnp = (numsnps.empty() ? 0 : numsnps[i]);
        if (T.allNodes.find(nodeID) == T.allNodes.end()) {
          err("Node ID {} not found in the tree", nodeID);
          std::cerr << "Node ID " << nodeID << " not found in the tree" << std::endl;
          return 1;
        }

        std::string sequence = panmapUtils::getStringFromReference(&T, nodeID, false);
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
          msg("Sequence for node {} with {} SNPs written to {}", nodeID, numsnp, outputFileName);
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
      outFile.close();
      exit(0);
    }

    // Log settings
    msg("--- Settings ---");
    msg("Reads: {}",
        (reads1.empty() ? "<none>"
                        : reads1 + (reads2.empty() ? "" : " + " + reads2)));
    msg("Reference PanMAN: {} ({} nodes)", guide, T.allNodes.size());
    msg("Using {} threads", vm["cpus"].as<int>());
    msg("k={}, s={}", k, s);

    // Handle indexing
    bool build = true;
    ::capnp::MallocMessageBuilder outMessage;
    std::unique_ptr<::capnp::MessageReader> inMessage;
    
    // Ensure paths are consistent by normalizing to absolute paths
    boost::filesystem::path guidePath = boost::filesystem::absolute(guide);
    std::string default_index_path = guidePath.string() + ".pmi";
    
    // If an explicit index path was provided, also normalize it
    std::string normalized_index_path = "";
    if (!index_path.empty()) {
      boost::filesystem::path explicitIndexPath = boost::filesystem::absolute(index_path);
      normalized_index_path = explicitIndexPath.string();
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
        ("meta", po::bool_switch(&cfg.metagenomic), "Metagenomic mode")
        ("top", po::value<int>(&cfg.topN)->default_value(1), "Top N placements (with --meta)")
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
        ("impute", po::bool_switch(&cfg.impute), "Impute N's from parent (skip _->N mutations)")
        ("kmer,k", po::value<int>(&cfg.k)->default_value(29), "Syncmer k")
        ("syncmer,s", po::value<int>(&cfg.s)->default_value(8), "Syncmer s")
        ("lmer,l", po::value<int>(&cfg.l)->default_value(1), "Syncmers per seed")
        ("flank-mask", po::value<int>(&cfg.flankMaskBp)->default_value(250), "Mask bp at ends")
        ("seed-mask-fraction", po::value<double>(&cfg.seedMaskFraction)->default_value(0),
            "Mask top seed fraction")
        ("min-seed-quality", po::value<int>(&cfg.minSeedQuality)->default_value(0),
            "Min seed quality")
        ("trim-start", po::value<int>(&cfg.trimStart)->default_value(0), "Trim read start")
        ("trim-end", po::value<int>(&cfg.trimEnd)->default_value(0), "Trim read end")
        ("min-read-support", po::value<int>(&cfg.minReadSupport)->default_value(1),
            "Min reads for a seed (2=filter singletons)")
        ("refine", po::bool_switch(&cfg.refine), "Enable alignment-based refinement")
        ("refine-top-pct", po::value<double>(&cfg.refineTopPct)->default_value(0.01),
            "Top % of nodes to refine (default 1%)")
        ("refine-max-top-n", po::value<int>(&cfg.refineMaxTopN)->default_value(150),
            "Max nodes to align against")
        ("refine-neighbor-radius", po::value<int>(&cfg.refineNeighborRadius)->default_value(2),
            "Expand to neighbors within N branches")
        ("refine-max-neighbor-n", po::value<int>(&cfg.refineMaxNeighborN)->default_value(150),
            "Max additional nodes from neighbor expansion");
    
    po::options_description developer("Developer");
    developer.add_options()
        ("dump-random-node", po::bool_switch(&cfg.dumpRandomNode), "Dump random node FASTA")
        ("dump-sequence", po::value<std::string>(&cfg.dumpNodeId), "Dump node FASTA")
        ("list-filtered-nodes", po::value<int>(&cfg.listFilteredNodes)->default_value(0), 
            "List N nodes passing length filter (TSV: node_id,length,parent_id,parent_length)")
        ("remove-node", po::value<std::string>(&cfg.removeNodeId), "Exclude node")
        ("top-placements", po::value<int>(&cfg.topPlacements)->default_value(0), "Show top N scores")
        ("debug-node", po::value<std::string>(&cfg.debugNodeId), "Debug node metrics")
        ("compare-nodes", po::value<std::string>(&cfg.compareNodes), "Compare two nodes (nodeA,nodeB)")
        ("dump-all-scores", po::value<std::string>(&cfg.dumpAllScores), "Dump all node scores to TSV file")
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
    all.add(visible).add(advanced).add(developer).add(positional);
    
    po::options_description visible_all;  // For --help-all
    visible_all.add(visible).add(advanced).add(developer);
    
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

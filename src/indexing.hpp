#pragma once

#include "capnp/list.h"
#include "coordinates.hpp"
#include "panman.hpp"
#include <capnp/serialize.h>
#include <capnp/message.h>
#include "index.capnp.h"
#include "state.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <optional> // Needed for kmer_utils
#include <string_view> // Needed for kmer_utils
#include <tuple>
#include <unordered_map>
#include <vector>
#include <functional>
#include <memory>
#include <limits> // Needed for kmer_utils implementation
#include "placement.hpp"
#include "logging.hpp" // Needed for kmer_utils implementation
#include <absl/container/flat_hash_map.h>

/**
 * @brief Macro for simplified error handling with consistent logging
 * 
 * @param operation The operation to perform
 * @param errorMsg The error message to display if an exception occurs
 * @param returnVal The value to return if an exception occurs
 */
#define SAFE_OPERATION(operation, errorMsg, returnVal) \
    try { \
        operation; \
    } catch (const std::exception& e) { \
        std::cerr << errorMsg << e.what() << std::endl; \
        return returnVal; \
    }

/**
 * @brief Macro for consistent logging format
 * 
 * @param level The log level (INFO, WARN, ERROR, DEBUG)
 * @param format The format string (must use {} for placeholders)
 * @param ... Additional arguments to be formatted into the message
 */
#define LOG(level, format, ...) \
    logging::level(format, __VA_ARGS__)

// Forward declarations
namespace placement {
class PlacementEngine;
struct PlacementResult;
}

// ---> START kmer_utils namespace declarations <---
/**
 * @namespace kmer_utils
 * @brief Utilities for k-mer extraction, dictionary lookup, and processing
 */
namespace kmer_utils {

/**
 * @brief Lookup a k-mer in the dictionary and return its ID and end offset
 * 
 * @param kmerSeq The k-mer sequence to look up
 * @param pos The starting position of the k-mer
 * @param endPos The end position of the k-mer
 * @param dictionary The k-mer dictionary mapping sequences to IDs
 * @return Optional pair of dictionary ID and end position offset
 */
std::optional<std::pair<uint32_t, uint16_t>> lookupKmerInDictionary(
    const std::string& kmerSeq,
    int64_t pos,
    int64_t endPos,
    const absl::flat_hash_map<std::string, uint32_t>& dictionary);

/**
 * @brief Safely extract a k-mer from a position in a node
 *
 * @param stateManager StateManager pointer for sequence extraction
 * @param nodeId ID of the node to extract from
 * @param pos Starting position to extract from
 * @param k K-mer length
 * @return Optional k-mer sequence if extraction succeeds
 */
std::optional<std::string> extractKmer(
    const state::StateManager &stateManager, std::string_view nodeId, int64_t pos, int k);

/**
 * @brief More efficient version using string_view to avoid copies.
 * 
 * @param stateManager StateManager reference.
 * @param nodeId Node ID.
 * @param pos Starting position.
 * @param k K-mer length.
 * @param buffer Buffer provided by caller to store extracted k-mer.
 * @return Optional string_view of the k-mer in the buffer.
 */
std::optional<std::string_view> extractKmerView(
    const state::StateManager& stateManager,
    const std::string& nodeId,
    int64_t pos,
    int k,
    std::string& buffer);

/**
 * @brief Calculate end position offset between start and end positions
 * 
 * @param startPos Starting position
 * @param endPos Ending position
 * @return End position offset, or k-1 if invalid
 */
uint16_t calculateEndOffset(int64_t startPos, int64_t endPos, int k);

/**
 * @brief Get existing seed at a position if available
 * 
 * @param stateManager StateManager pointer for seed access
 * @param pos Position to check for seed
 * @return Optional seed information
 */
std::optional<seeding::seed_t> getSeedAtPosition(
    state::StateManager& stateManager,
    int64_t pos);

/**
 * @brief Get existing seed at a position for a specific node if available
 * 
 * @param stateManager StateManager pointer for seed access
 * @param nodeId ID of the node to check
 * @param pos Position to check for seed
 * @return Optional seed information
 */
std::optional<seeding::seed_t> getSeedAtPosition(
    state::StateManager& stateManager,
    const std::string& nodeId,
    int64_t pos);

} // namespace kmer_utils
// ---> END kmer_utils namespace declarations <---


// ---> START kmer_utils implementations <---
inline std::optional<std::pair<uint32_t, uint16_t>> kmer_utils::lookupKmerInDictionary(
    const std::string& kmerSeq,
    int64_t pos,
    int64_t endPos,
    const absl::flat_hash_map<std::string, uint32_t>& dictionary) {
    
    // Validate input parameters
    if (kmerSeq.empty() || pos < 0 || endPos < pos) {
        return std::nullopt;
    }
    
    // Look up in dictionary
    auto it = dictionary.find(kmerSeq);
    if (it == dictionary.end()) {
        return std::nullopt;
    }
    
    // Calculate end offset
    uint16_t endOffset = static_cast<uint16_t>(endPos - pos);
    
    // Return dictionary ID and end offset wrapped in optional
    // Ensure make_pair and make_optional are available (requires <utility> and <optional>)
    return std::make_optional(std::make_pair(it->second, endOffset)); 
}

inline std::optional<std::string> kmer_utils::extractKmer(
    const state::StateManager &stateManager, std::string_view nodeId, int64_t pos, int k) {

    try {
        auto result = stateManager.extractKmer(nodeId, pos, k);
        
        if (result.first.empty()) {
            return std::nullopt;
        }
        
        // No need to strip gaps anymore - the StateManager now returns an ungapped k-mer
        return result.first;
    } catch (const std::exception& e) {
        logging::warn("Error in kmer_utils::extractKmer: {}", e.what());
        return std::nullopt;
    }
}

inline std::optional<std::string_view> kmer_utils::extractKmerView(
    const state::StateManager& stateManager,
    const std::string& nodeId,
    int64_t pos,
    int k,
    std::string& buffer) {  // Buffer provided by caller to avoid allocation
    
    // Handle the case where k is invalid
    if (k <= 0) {
        logging::err("Invalid k value: {}", k);
        return std::nullopt;
    }

    try {
        // Extract sequence into the provided buffer
        // The extractKmer method now handles node-specific caching internally
        // and returns ungapped k-mers
        auto result = stateManager.extractKmer(nodeId, pos, k);
        
        if (result.first.empty()) {
            return std::nullopt;
        }
        
        // Move result into buffer and create view
        buffer = std::move(result.first);
        std::string_view view(buffer.data(), std::min(k, static_cast<int>(buffer.length())));
        
        // Return the view wrapped in optional using std::in_place
        return std::optional<std::string_view>(std::in_place, view);
    } catch (const std::exception& e) {
        logging::err("Error extracting k-mer view at position {}: {}", pos, e.what());
        return std::nullopt;
    }
}

inline uint16_t kmer_utils::calculateEndOffset(int64_t startPos, int64_t endPos, int k) {
    // Validate positions
    if (endPos < startPos) {
        // Default to k-1 if positions are invalid
        return static_cast<uint16_t>(k - 1);
    }
    
    // Calculate offset (ensuring we don't exceed uint16_t range)
    int64_t offset = endPos - startPos;
    if (offset > std::numeric_limits<uint16_t>::max()) {
        return std::numeric_limits<uint16_t>::max();
    }
    
    return static_cast<uint16_t>(offset);
}

inline std::optional<seeding::seed_t> kmer_utils::getSeedAtPosition(
    state::StateManager& stateManager,
    int64_t pos) {
    
    if (pos < 0) {
        // It's generally better practice to throw exceptions for truly invalid arguments
        // rather than returning nullopt silently, especially if a negative position is never expected.
        // Logging could also be added here.
        throw std::invalid_argument("Invalid position " + std::to_string(pos) + 
                                  " in getSeedAtPosition");
    }
    
    // Find the root node ID
    std::string rootNodeId;
    for (const auto& [nodeId, hierarchy] : stateManager.getNodeHierarchy()) {
        if (hierarchy.parentId.empty()) {
            rootNodeId = nodeId;
            break;
        }
    }
    
    if (rootNodeId.empty()) {
        logging::warn("No root node found for getSeedAtPosition({})", pos);
        return std::nullopt;
    }
    
    // Delegate to node-specific StateManager method with the root node ID
    return stateManager.getSeedAtPosition(rootNodeId, pos);
}

inline std::optional<seeding::seed_t> kmer_utils::getSeedAtPosition(
    state::StateManager& stateManager,
    const std::string& nodeId,
    int64_t pos) {
    
    if (pos < 0) {
        throw std::invalid_argument("Invalid position " + std::to_string(pos) + 
                                  " in getSeedAtPosition");
    }
    // Delegate to node-specific StateManager method
    return stateManager.getSeedAtPosition(nodeId, pos);
}


namespace indexing {

/**
 * @struct PanmapParams
 * @brief Parameters for the indexing process
 */
struct PanmapParams {
    int k = 32; // k-mer length (must be 0 < k <= 2048)
    int s = 8; // s-mer length for syncmer (must be 0 < s <= k)
    bool open = false; // Use open syncmers (true) or closed syncmers (false)
    bool debugOutput = true; // Enable debug file output (default: true)
    
    bool validate() const {
        return k > 0 && s > 0 && s < k;
    }
};


/**
 * @struct SeedMutations
 * @brief Cap'n Proto generated type for storing seed mutations between nodes
 * 
 * Stores information about k-mer seeds (syncmers) that are added or removed
 * when traversing from a parent node to a child node in the tree.
 */

/**
 * @struct GapMutations
 * @brief Cap'n Proto generated type for storing "gap map" insertions/deletions between nodes
 * 
 * Stores information about gaps that are added or removed when
 * traversing from a parent node to a child node in the tree.
 */


/**
 * @brief Perform main indexing traversal of a PanMAN
 * 
 * Traverse the tree to index seed and gap mutations for each node
 * and write them to the output builders.
 * 
 * @param tree Pointer to PanMAN tree
 * @param rootNode Starting node for traversal
 * @param perNodeSeedMutations Output Cap'n Proto builder for seed mutations
 * @param perNodeGapMutations Output Cap'n Proto builder for gap mutations
 * @param params Indexing parameters (k, s, etc.)
 * @param performBacktracking Whether to perform node backtracking (default: false)
 * @throw std::invalid_argument If tree or rootNode is null, or if params are invalid
 */
void indexingTraversal(
    panmanUtils::Tree *tree, 
    panmanUtils::Node *rootNode,
    ::capnp::List<::SeedMutations>::Builder perNodeSeedMutations,
    ::capnp::List<::GapMutations>::Builder perNodeGapMutations,
    const PanmapParams &params,
    bool performBacktracking = false,
    state::StateManager *stateManager = nullptr);

/**
 * @brief Build panmap index for a PanMAN pangenome tree
 * 
 * Creates a panmap index for a PanMAN pangenome tree,
 * including seed mutations and related data.
 * 
 * @param tree Pointer to PanMAN tree
 * @param index Cap'n Proto builder for the index
 * @param k Length of syncmer seeds
 * @param s Syncmer s parameter
 * @param message Cap'n Proto message builder for serialization
 * @param indexPath Path to write the index file
 * @throw std::invalid_argument If tree is null or k,s parameters are invalid
 */
void index(panmanUtils::Tree *tree, Index::Builder &index, int k, int s,
           ::capnp::MallocMessageBuilder& message, const std::string& indexPath);

/**
 * @brief Build additional data structures to accelerate placement
 * 
 * Creates additional index structures that speed up the placement process,
 * including:
 * - Node path information (ancestry, levels)
 * - Block information (active blocks in each node)
 * - Ancestor-descendant relationship matrix
 * 
 * @param tree Pointer to PanMAN tree
 * @param index Cap'n Proto builder for the index
 * @param stateManager State manager with precomputed block and node information
 * @throw std::invalid_argument If tree is null
 */
void buildPlacementAccelerationData(
    panmanUtils::Tree *tree,
    Index::Builder &index,
    const state::StateManager& stateManager);

/**
 * @brief Process seeds between specified coordinate ranges
 * 
 * Computes k-mer seeds (syncmers) for specific regions of the sequence.
 * Used internally by the indexing process to update seed information
 * after mutations.
 * 
 * @param stateManager State manager containing sequence information
 * @param nodeId ID of the node to process
 * @param initialRecompRanges Coordinate ranges to recompute seeds for
 * @param forceBlockId Block ID to force for ternary mask (-1 for none)
 */
void seedBetweenRanges(
    state::StateManager& stateManager,
    const std::string& nodeId,
    const std::vector<coordinates::CoordRange>& initialRecompRanges, 
    int32_t forceBlockId = -1);

/**
 * Extract seed changes from recomputation ranges
 * 
 * @param stateManager State manager for sequence access
 * @param nodeId ID of the node being processed
 * @param ranges Vector of coordinate ranges to process
 * @param k K-mer size
 * @param s S-mer size for syncmer
 * @return Vector of seed change tuples (position, oldHash, oldExists, oldEndPos, newHash, newIsReverse, newEndPos)
 */
std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t, bool, int64_t>> 
extractSeedChanges(state::StateManager& stateManager, 
                   const std::string& nodeId,
                   const std::vector<coordinates::CoordRange>& ranges,
                   int16_t k,
                   int8_t s);

/**
 * Store seed changes in the index
 * 
 * Helper function to serialize seed changes into the CapnProto index format
 * 
 * @param seedChanges Vector of seed change tuples
 * @param seedMutation CapnProto builder for the seed mutation
 */
void storeSeedChangesInIndex(
    const std::vector<std::tuple<int64_t, size_t, bool, int64_t, int64_t, bool, int64_t>>& seedChanges,
    ::SeedMutations::Builder& seedMutation);

/**
 * Process mutations for a specific node
 * 
 * Helper function to apply block and nucleotide mutations to a node
 * from either common ancestor or root node
 * 
 * @param stateManager StateManager to apply mutations to
 * @param tree Pointer to the tree
 * @param node Node to process mutations for
 * @param commonAncestor Common ancestor node (can be nullptr if none)
 * @param nodePaths Map of paths from root to each node
 * @param ab_before_output_str Pointer to store AB_BEFORE logging
 */
void processMutationsForNode(
    state::StateManager &stateManager, 
    panmanUtils::Tree *tree,
    panmanUtils::Node *node,
    panmanUtils::Node *commonAncestor,
    const std::unordered_map<std::string, std::vector<panmanUtils::Node *>> &nodePaths,
    std::string* ab_before_output_str
    );


/**
 * Process nodes by level with parallelism strategy based on level size
 * 
 * @param levelNodes Vector of nodes at the current level
 * @param processNode Function to process a single node
 */
void processNodesByLevel(
    const std::vector<panmanUtils::Node*>& levelNodes,
    const std::function<void(panmanUtils::Node*)>& processNode);

/**
 * Initialize state manager with sequence data from the reference
 * 
 * @param tree The tree containing the pangenome
 * @param rootNode The root node of the tree
 * @param kmerSize The k-mer size to set in state manager
 * @param smerSize The s-mer size to set in state manager
 * @return Initialized state manager
 */
std::unique_ptr<state::StateManager> initializeStateManager(
    panmanUtils::Tree* tree,
    panmanUtils::Node* rootNode,
    int kmerSize,
    int smerSize);

/**
 * Initialize state manager with minimal computation for placement 
 * 
 * This is an optimized version of initializeStateManager that skips some 
 * unnecessary steps during placement, since we're only interested in using
 * the pre-computed seed information from the index.
 * 
 * @param tree The tree containing the pangenome
 * @param rootNode The root node of the tree
 * @param kmerSize The k-mer size to set in state manager
 * @param smerSize The s-mer size to set in state manager
 * @return Initialized state manager with minimal computation
 */
std::unique_ptr<state::StateManager>
initializeStateManagerLight(panmanUtils::Tree* tree, panmanUtils::Node* rootNode,
                          int kmerSize, int smerSize);

/**
 * Debug function to dump detailed node and index information to files
 * 
 * @param tree The tree being processed
 * @param node The node being processed
 * @param stateManager StateManager containing node state
 * @param phase Current processing phase ("indexing" or "placement")
 * @param nodeIndex Index of the node in the current traversal (for first N nodes)
 * @param dumpFilename Base filename for dumping data (defaults to "panmap_debug" if empty)
 */
void dumpNodeDetails(
    const panmanUtils::Tree* tree,
    const panmanUtils::Node* node,
    const state::StateManager& stateManager,
    const std::string& phase,
    size_t nodeIndex,
    const std::string& dumpFilename = "panmap_debug");

/**
 * Debug function to dump complete index information to a file
 * 
 * @param index Index data to dump
 * @param dumpFilename Filename for dumping data (defaults to "panmap_debug" if empty)
 */
void dumpIndexData(
    const Index::Reader& index,
    const std::string& dumpFilename = "panmap_debug");

/**
 * Debug function to create a detailed diagnostic dump for the recomputation process
 * 
 * @param nodeId ID of the node being processed
 * @param origMutations Original mutations that caused recomputation
 * @param initialRanges Initial recomputation ranges derived from mutations
 * @param mergedRanges Merged recomputation ranges after optimization
 * @param stateManager State manager with precomputed block and node information
 * @param k K-mer size for seed generation
 * @param s S-mer size for syncmer
 */
void createRecompDumpFile(
    const std::string& nodeId,
    const std::vector<std::pair<int64_t, uint8_t>>& origMutations,
    const std::vector<coordinates::CoordRange>& initialRanges,
    const std::vector<coordinates::CoordRange>& mergedRanges,
    state::StateManager& stateManager,
    int k,
    int s);

/**
 * @brief Process seed changes and encode them in quaternary format (4 values)
 * 
 * Each seed change is encoded using 2 bits:
 * - 0: seed unchanged (off → off)
 * - 1: seed deleted (on → off)
 * - 2: seed added (off → on)
 * - 3: seed modified (on → on but changed)
 * 
 * This quaternary encoding is more efficient than the previous ternary encoding,
 * allowing more precise tracking of seed changes and more efficient storage.
 * 
 * @param seedChanges Vector of seed changes tuples
 * @param basePositions Output vector of base positions
 * @param bitMasks Output vector of quaternary-encoded bit masks
 */
void processSeedChanges(
    const std::vector<std::tuple<int64_t, bool, bool, std::optional<size_t>,
                                 std::optional<size_t>, std::optional<bool>,
                                 std::optional<bool>, std::optional<int64_t>,
                                 std::optional<int64_t>>> &seedChanges,
    std::vector<int64_t> &basePositions, std::vector<uint64_t> &bitMasks);

/**
 * @brief Decode quaternary-encoded seed changes back to original form
 * 
 * Reverses the encoding performed by processSeedChanges to recover
 * the original seed change information.
 * 
 * @param basePositions Vector of base positions
 * @param bitMasks Vector of quaternary-encoded bit masks
 * @return Vector of (pos, wasOn, isOn) tuples
 */
std::vector<std::tuple<int64_t, bool, bool>>
decodeSeedChanges(const std::vector<int64_t> &basePositions,
                  const std::vector<uint64_t> &tritMasks);

/**
 * @brief Recompute seeds for a node and optionally write to Cap'n Proto format
 *
 * This function handles:
 * 1. Expanding recomputation ranges to ensure complete k-mers
 * 2. Extracting sequences from those ranges
 * 3. Computing new seeds
 * 4. Tracking seed changes (added/deleted/modified)
 * 5. Encoding seed changes in quaternary format
 * 6. Writing to Cap'n Proto if builder is provided, or node state otherwise
 *
 * @param stateManager The state manager containing coordinate system
 * @param engine The placement engine for seed storage
 * @param node The node to recompute seeds for
 * @param k K-mer size
 * @param s S-mer size
 * @param perNodeSeedMutations Optional Cap'n Proto builder for direct writing
 * @param kmerDictionary Optional pointer to k-mer dictionary for direct writing
 * @param uniqueKmersCollector Collector for unique k-mers
 */
void recomputeSeeds(
    state::StateManager &stateManager, 
    placement::PlacementEngine &engine,
    panmanUtils::Node *node, int k, int s,
    ::capnp::List<SeedMutations>::Builder* perNodeSeedMutations,
    const std::unordered_map<std::string, uint32_t>* kmerDictionary,
    std::unordered_set<std::string>& uniqueKmersCollector,
    std::ofstream& debugSeedFile);

/**
 * Initialize block sequences in the state manager.
 *
 * @param stateManager The state manager to update
 * @param blockSequences Map of block sequences
 * @param blockRanges Map of block coordinate ranges
 */
void initializeBlockSequences(
    state::StateManager &stateManager,
    const absl::flat_hash_map<int32_t, std::string> &blockSequences,
    const absl::flat_hash_map<int32_t, coordinates::CoordRange> &blockRanges);

/**
 * Process nodes at a given level using the provided function.
 *
 * @param nodes Vector of nodes to process
 * @param processFunction Function to call for each node
 */
void processNodesByLevel(
    const std::vector<panmanUtils::Node*>& nodes,
    std::function<void(panmanUtils::Node*)> processFunction);

/**
 * Fill a vector with nodes in depth-first traversal order.
 *
 * @param node The root node to start traversal from
 * @param nodes Vector to fill with nodes
 */
void fillNodesDepthFirst(panmanUtils::Node* node, std::vector<panmanUtils::Node*>& nodes);

/**
 * Index a PanMAT tree with parallel processing.
 *
 * @param tree The tree to index
 * @param stateManager The state manager to use
 * @param commonAncestor The common ancestor node
 * @param seeding Whether to generate seeds
 * @param k K-mer size
 * @param s S-mer size
 * @param threads Number of threads to use
 * @param perNodeSeedMutations Cap'n Proto builder for seed mutations
 * @param outMessage Cap'n Proto message builder for output
 * @param indexPath Path to write the index file
 */
void parallelIndexPan(
    panmanUtils::Tree *tree, state::StateManager *stateManager, 
    panmanUtils::Node *commonAncestor,
    bool seeding, int k, int s, int threads,
    ::capnp::List<SeedMutations>::Builder *perNodeSeedMutations,
    ::capnp::MallocMessageBuilder *outMessage,
    const std::string &indexPath);

// Forward declaration of the writeCapnp function from main.cpp
extern void writeCapnp(::capnp::MallocMessageBuilder &message, const std::string &path);

void logNodeDetailsToFile(
    const panmanUtils::Tree* tree,
    const panmanUtils::Node* node,
    const state::StateManager& stateManager,
    const std::string& phase,
    size_t nodeIndex,
    const std::string& dumpFilename = "panmap_debug");

} // namespace indexing
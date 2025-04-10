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
#include <tuple>
#include <unordered_map>
#include <vector>
#include <functional>
#include <memory>
#include "placement.hpp"

// Forward declarations
namespace placement {
class PlacementEngine;
struct PlacementResult;
}

namespace indexing {

/**
 * @struct PanmapParams
 * @brief Parameters for the indexing process
 */
struct PanmapParams {
    int k = 31; // k-mer length (must be 0 < k <= 2048)
    int s = 8; // s-mer length for syncmer (must be 0 < s <= k)
    bool open = true; // Use open syncmers (true) or closed syncmers (false)
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
 * @param node Node to process mutations for
 * @param commonAncestor Common ancestor node (can be nullptr if none)
 * @param nodePaths Map of paths from root to each node
 */
void processMutationsForNode(
    state::StateManager& stateManager,
    panmanUtils::Node* node,
    panmanUtils::Node* commonAncestor,
    const std::unordered_map<std::string, std::vector<panmanUtils::Node*>>& nodePaths);


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
    const std::vector<std::tuple<int64_t, bool, bool, 
                               std::optional<size_t>, std::optional<size_t>,
                               std::optional<bool>, std::optional<bool>, 
                               std::optional<int64_t>, std::optional<int64_t>>>& seedChanges,
    std::vector<int64_t>& basePositions,
    std::vector<uint64_t>& bitMasks);

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
decodeSeedChanges(const std::vector<int64_t>& basePositions, 
                 const std::vector<uint64_t>& bitMasks);

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
 */
void recomputeSeeds(state::StateManager &stateManager, 
                   placement::PlacementEngine &engine,
                   panmanUtils::Node *node, int k, int s,
                   ::capnp::List<SeedMutations>::Builder* perNodeSeedMutations = nullptr);

/**
 * Initialize block sequences in the state manager.
 *
 * @param stateManager The state manager to update
 * @param blockSequences Map of block sequences
 * @param blockRanges Map of block coordinate ranges
 */
void initializeBlockSequences(
    state::StateManager &stateManager,
    const std::unordered_map<int32_t, std::string> &blockSequences,
    const std::unordered_map<int32_t, coordinates::CoordRange> &blockRanges);

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
 */
void parallelIndexPan(panmanUtils::Tree *tree, state::StateManager *stateManager, 
                     panmanUtils::Node *commonAncestor,
                     bool seeding, int k, int s, int threads);

// Forward declaration of the periodicallyFlushCapnp function from main.cpp
void periodicallyFlushCapnp(::capnp::MallocMessageBuilder &message,
                            const std::string &path,
                            bool forceFlush = false,
                            size_t opCount = 0,
                            size_t flushInterval = 1000);

} // namespace indexing
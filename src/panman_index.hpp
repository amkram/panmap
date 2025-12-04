#pragma once

/**
 * @file panman_index.hpp
 * @brief Fast single-node genome reconstruction from PanMAN files
 * 
 * This module provides utilities to:
 * 1. Create an index file (.pmx) alongside a panman that enables fast single-node access
 * 2. Reconstruct a single genome by loading only necessary data from disk
 * 
 * The key insight is that Cap'n Proto supports lazy deserialization, but XZ compression
 * prevents random access. We solve this by:
 * - Converting XZ to ZSTD with seekable frames for parallel decompression
 * - Creating a secondary index with node-to-offset mappings
 * - Using mmap + Cap'n Proto's lazy access to only touch required data
 */

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <cstdint>

namespace panman_index {

// ============================================================================
// Index File Format (.pmx - PanMan indeX)
// ============================================================================

/**
 * File format:
 * 
 * [Header]
 *   - Magic: "PMX1" (4 bytes)
 *   - Version: uint32_t
 *   - Flags: uint32_t (reserved)
 *   - NumNodes: uint32_t
 *   - NewickOffset: uint64_t (offset to newick string in decompressed data)
 *   - NewickLength: uint64_t
 *   - ConsensusOffset: uint64_t (offset to consensus sequence data)
 *   - GapsOffset: uint64_t (offset to gap data)
 * 
 * [Node Index Table] - sorted by node name for binary search
 *   For each node:
 *   - NameOffset: uint32_t (offset into name string table)
 *   - NameLength: uint16_t
 *   - ParentIndex: uint32_t (index in this table, UINT32_MAX for root)
 *   - PreorderIndex: uint32_t (index in the Cap'n Proto nodes list)
 *   - MutationCount: uint32_t (number of mutations on this node)
 * 
 * [Name String Table]
 *   - Concatenated null-terminated node names
 * 
 * [Sorted Name Index]
 *   - Array of uint32_t indices into Node Index Table, sorted by name
 */

constexpr uint32_t PMX_MAGIC = 0x31584D50;  // "PMX1" in little-endian
constexpr uint32_t PMX_VERSION = 1;

#pragma pack(push, 1)
struct PmxHeader {
    uint32_t magic;
    uint32_t version;
    uint32_t flags;
    uint32_t numNodes;
    uint64_t newickOffset;
    uint64_t newickLength;
    uint64_t consensusOffset;
    uint64_t gapsOffset;
    uint64_t nodeTableOffset;
    uint64_t nameTableOffset;
    uint64_t sortedIndexOffset;
};

struct PmxNodeEntry {
    uint32_t nameOffset;      // Offset into name string table
    uint16_t nameLength;      // Length of name (excluding null terminator)
    uint16_t reserved;        // Padding for alignment
    uint32_t parentIndex;     // Index of parent in node table (UINT32_MAX = root)
    uint32_t preorderIndex;   // Index in Cap'n Proto nodes list
    uint32_t mutationCount;   // Number of mutations
    uint64_t mutationOffset;  // Byte offset to mutations in Cap'n Proto data
};
#pragma pack(pop)

// ============================================================================
// Index Builder
// ============================================================================

/**
 * @brief Builds a .pmx index file for a panman
 * 
 * This also optionally recompresses the panman from XZ to ZSTD for faster loading.
 */
class IndexBuilder {
public:
    /**
     * @brief Build index for a panman file
     * 
     * @param panmanPath Path to .panman file (XZ compressed)
     * @param outputPath Path for .pmx index file (default: panmanPath + ".pmx")
     * @param recompressToZstd If true, also create a .panman.zst file with ZSTD compression
     * @return true on success
     */
    static bool build(
        const std::string& panmanPath,
        const std::string& outputPath = "",
        bool recompressToZstd = false
    );
};

// ============================================================================
// Fast Genome Extractor
// ============================================================================

/**
 * @brief Extracts a single genome from a panman using the index
 * 
 * This class is optimized for extracting one or a few genomes without
 * loading the entire panman into memory.
 */
class GenomeExtractor {
public:
    /**
     * @brief Open a panman for genome extraction
     * 
     * @param panmanPath Path to .panman file
     * @param indexPath Path to .pmx index file (default: panmanPath + ".pmx")
     */
    GenomeExtractor(
        const std::string& panmanPath,
        const std::string& indexPath = ""
    );
    
    ~GenomeExtractor();
    
    /**
     * @brief Check if the extractor is ready
     */
    bool isOpen() const { return isOpen_; }
    
    /**
     * @brief Get list of all node IDs
     */
    std::vector<std::string> getNodeIds() const;
    
    /**
     * @brief Check if a node exists
     */
    bool hasNode(const std::string& nodeId) const;
    
    /**
     * @brief Extract genome sequence for a node
     * 
     * @param nodeId Node identifier
     * @param aligned If true, include gap characters for alignment
     * @return Genome sequence, or empty string on error
     */
    std::string extractGenome(const std::string& nodeId, bool aligned = false);
    
    /**
     * @brief Get the path from root to a node
     * 
     * @param nodeId Target node identifier
     * @return Vector of node IDs from root to target (inclusive)
     */
    std::vector<std::string> getPathToNode(const std::string& nodeId) const;
    
    /**
     * @brief Get statistics about extraction
     */
    struct Stats {
        size_t bytesRead;           // Total bytes read from disk
        size_t nodesTraversed;      // Number of nodes in path
        size_t mutationsApplied;    // Total mutations applied
        double extractionTimeMs;    // Time for last extraction
    };
    Stats getStats() const { return stats_; }

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
    bool isOpen_ = false;
    Stats stats_ = {};
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert XZ-compressed panman to ZSTD with parallel decompression support
 * 
 * @param inputPath Path to XZ-compressed .panman file
 * @param outputPath Output path (default: inputPath + ".zst")
 * @param compressionLevel ZSTD compression level (1-22, default 3)
 * @param frameSize Size of each ZSTD frame for parallel decompression (default 4MB)
 * @return true on success
 */
bool convertXzToZstd(
    const std::string& inputPath,
    const std::string& outputPath = "",
    int compressionLevel = 3,
    size_t frameSize = 4 * 1024 * 1024
);

/**
 * @brief Check if a panman has an associated index
 */
bool hasIndex(const std::string& panmanPath);

/**
 * @brief Get the default index path for a panman
 */
std::string getIndexPath(const std::string& panmanPath);

} // namespace panman_index

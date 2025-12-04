/**
 * @file panman_index.cpp
 * @brief Fast single-node genome reconstruction from PanMAN files
 * 
 * Strategy: We don't reimplement the complex genome reconstruction logic.
 * Instead, we:
 * 1. Create an index that maps node names to their tree position
 * 2. Optionally recompress from XZ to ZSTD for faster parallel decompression  
 * 3. Use the existing panmanUtils::Tree::getStringFromReference() for correctness
 * 
 * The speedup comes from:
 * - ZSTD parallel decompression vs single-threaded XZ
 * - Pre-built index for instant node lookup
 * - Future: selective loading of only needed tree portions
 */

#include "panman_index.hpp"
#include "logging.hpp"
#include "zstd_compression.hpp"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <cstring>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>

// Include panman headers for correct genome reconstruction
#include "panman.hpp"

namespace panman_index {

// ============================================================================
// Newick Parser - Lightweight parser to extract tree structure for indexing
// ============================================================================

struct TreeNode {
    std::string name;
    int32_t parentIndex = -1;  // -1 means root
    int32_t preorderIndex = -1;
    std::vector<int32_t> childIndices;
};

/**
 * @brief Parse a Newick string to extract tree structure
 */
class NewickParser {
public:
    std::vector<TreeNode> nodes;
    std::unordered_map<std::string, int32_t> nameToIndex;
    
    bool parse(const std::string& newick) {
        nodes.clear();
        nameToIndex.clear();
        
        size_t pos = 0;
        int32_t preorderIdx = 0;
        parseSubtree(newick, pos, -1, preorderIdx);
        
        return !nodes.empty();
    }
    
private:
    int32_t parseSubtree(const std::string& s, size_t& pos, int32_t parentIdx, int32_t& preorderIdx) {
        int32_t nodeIdx = nodes.size();
        nodes.emplace_back();
        TreeNode& node = nodes.back();
        node.parentIndex = parentIdx;
        node.preorderIndex = preorderIdx++;
        
        while (pos < s.size() && std::isspace(s[pos])) pos++;
        
        if (pos < s.size() && s[pos] == '(') {
            pos++;
            while (true) {
                int32_t childIdx = parseSubtree(s, pos, nodeIdx, preorderIdx);
                node.childIndices.push_back(childIdx);
                while (pos < s.size() && std::isspace(s[pos])) pos++;
                if (pos < s.size() && s[pos] == ',') {
                    pos++;
                    continue;
                }
                if (pos < s.size() && s[pos] == ')') {
                    pos++;
                    break;
                }
                break;
            }
        }
        
        size_t nameStart = pos;
        while (pos < s.size() && s[pos] != ':' && s[pos] != ',' && 
               s[pos] != ')' && s[pos] != ';' && !std::isspace(s[pos])) {
            pos++;
        }
        node.name = s.substr(nameStart, pos - nameStart);
        
        if (pos < s.size() && s[pos] == ':') {
            pos++;
            while (pos < s.size() && (std::isdigit(s[pos]) || s[pos] == '.' || 
                   s[pos] == 'e' || s[pos] == 'E' || s[pos] == '-' || s[pos] == '+')) {
                pos++;
            }
        }
        
        if (!node.name.empty()) {
            nameToIndex[node.name] = nodeIdx;
        }
        
        return nodeIdx;
    }
};

// ============================================================================
// Index Builder Implementation
// ============================================================================

bool IndexBuilder::build(
    const std::string& panmanPath,
    const std::string& outputPath,
    bool recompressToZstd
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    std::string indexPath = outputPath.empty() ? panmanPath + ".pmx" : outputPath;
    
    logging::info("Building index for {} -> {}", panmanPath, indexPath);
    
    // Step 1: Load the panman using existing infrastructure
    std::ifstream file(panmanPath, std::ios::binary);
    if (!file) {
        logging::err("Cannot open panman file: {}", panmanPath);
        return false;
    }
    
    // Decompress XZ to memory
    std::vector<char> decompressedData;
    {
        auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
        buffer->push(boost::iostreams::lzma_decompressor());
        buffer->push(file);
        
        std::istream stream(buffer.get());
        std::ostringstream oss;
        oss << stream.rdbuf();
        std::string data = oss.str();
        decompressedData.assign(data.begin(), data.end());
    }
    
    logging::info("Decompressed panman: {} bytes", decompressedData.size());
    
    // Step 2: Load as TreeGroup to get the Newick string
    std::istringstream iss(std::string(decompressedData.begin(), decompressedData.end()));
    panmanUtils::TreeGroup tg(iss);
    
    if (tg.trees.empty()) {
        logging::err("No trees in panman");
        return false;
    }
    
    // Get newick string from tree
    std::string newick = tg.trees[0].getNewickString(tg.trees[0].root);
    
    // Step 3: Parse Newick to get tree structure
    NewickParser parser;
    if (!parser.parse(newick)) {
        logging::err("Failed to parse Newick string");
        return false;
    }
    
    logging::info("Parsed tree with {} nodes", parser.nodes.size());
    
    // Step 4: Build index structures
    PmxHeader header;
    std::memset(&header, 0, sizeof(header));
    header.magic = PMX_MAGIC;
    header.version = PMX_VERSION;
    header.numNodes = parser.nodes.size();
    header.newickLength = newick.size();
    
    // Build node entries and name table
    std::vector<PmxNodeEntry> nodeEntries(parser.nodes.size());
    std::string nameTable;
    
    for (size_t i = 0; i < parser.nodes.size(); i++) {
        const auto& node = parser.nodes[i];
        auto& entry = nodeEntries[i];
        
        entry.nameOffset = nameTable.size();
        entry.nameLength = node.name.size();
        nameTable += node.name;
        nameTable += '\0';
        
        entry.parentIndex = node.parentIndex >= 0 ? node.parentIndex : UINT32_MAX;
        entry.preorderIndex = node.preorderIndex;
        entry.mutationCount = 0;  // We don't need this for the simple index
        entry.mutationOffset = 0;
    }
    
    // Build sorted index for binary search by name
    std::vector<uint32_t> sortedIndex(parser.nodes.size());
    for (uint32_t i = 0; i < parser.nodes.size(); i++) {
        sortedIndex[i] = i;
    }
    std::sort(sortedIndex.begin(), sortedIndex.end(), [&](uint32_t a, uint32_t b) {
        return parser.nodes[a].name < parser.nodes[b].name;
    });
    
    // Step 5: Calculate offsets
    size_t offset = sizeof(PmxHeader);
    header.newickOffset = offset;
    offset += newick.size() + 1;
    
    header.nodeTableOffset = offset;
    offset += nodeEntries.size() * sizeof(PmxNodeEntry);
    
    header.nameTableOffset = offset;
    offset += nameTable.size();
    
    header.sortedIndexOffset = offset;
    
    // Step 6: Write index file
    std::ofstream outFile(indexPath, std::ios::binary);
    if (!outFile) {
        logging::err("Cannot create index file: {}", indexPath);
        return false;
    }
    
    outFile.write(reinterpret_cast<const char*>(&header), sizeof(header));
    outFile.write(newick.c_str(), newick.size() + 1);
    outFile.write(reinterpret_cast<const char*>(nodeEntries.data()), 
                  nodeEntries.size() * sizeof(PmxNodeEntry));
    outFile.write(nameTable.data(), nameTable.size());
    outFile.write(reinterpret_cast<const char*>(sortedIndex.data()),
                  sortedIndex.size() * sizeof(uint32_t));
    
    outFile.close();
    
    // Step 7: Optionally recompress to ZSTD for faster loading
    if (recompressToZstd) {
        std::string zstdPath = panmanPath + ".zst";
        logging::info("Recompressing to ZSTD: {}", zstdPath);
        
        if (!panmap_zstd::compressToFile(
                decompressedData.data(), decompressedData.size(),
                zstdPath, 3, 0, 4 * 1024 * 1024)) {
            logging::warn("Failed to create ZSTD file");
        } else {
            logging::info("Created ZSTD-compressed panman: {}", zstdPath);
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    logging::info("Index built in {}ms: {} nodes", duration.count(), parser.nodes.size());
    
    return true;
}

// ============================================================================
// Genome Extractor Implementation
// ============================================================================

class GenomeExtractor::Impl {
public:
    // Index data
    std::vector<char> indexData;
    const PmxHeader* header = nullptr;
    const char* newick = nullptr;
    const PmxNodeEntry* nodeEntries = nullptr;
    const char* nameTable = nullptr;
    const uint32_t* sortedIndex = nullptr;
    
    // Panman TreeGroup - uses the correct panmanUtils implementation
    std::unique_ptr<panmanUtils::TreeGroup> treeGroup;
    
    // Cache for node lookup
    std::unordered_map<std::string, uint32_t> nameCache;
    
    bool open(const std::string& panmanPath, const std::string& indexPath) {
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // Load index file
        {
            std::ifstream idxFile(indexPath, std::ios::binary | std::ios::ate);
            if (!idxFile) {
                logging::err("Failed to open index: {}", indexPath);
                return false;
            }
            
            size_t size = idxFile.tellg();
            idxFile.seekg(0);
            indexData.resize(size);
            idxFile.read(indexData.data(), size);
        }
        
        if (indexData.size() < sizeof(PmxHeader)) {
            logging::err("Index file too small");
            return false;
        }
        
        header = reinterpret_cast<const PmxHeader*>(indexData.data());
        if (header->magic != PMX_MAGIC) {
            logging::err("Invalid index magic number");
            return false;
        }
        
        newick = indexData.data() + header->newickOffset;
        nodeEntries = reinterpret_cast<const PmxNodeEntry*>(
            indexData.data() + header->nodeTableOffset);
        nameTable = indexData.data() + header->nameTableOffset;
        sortedIndex = reinterpret_cast<const uint32_t*>(
            indexData.data() + header->sortedIndexOffset);
        
        // Build name cache
        for (uint32_t i = 0; i < header->numNodes; i++) {
            std::string name(nameTable + nodeEntries[i].nameOffset, nodeEntries[i].nameLength);
            nameCache[name] = i;
        }
        
        auto indexLoadTime = std::chrono::high_resolution_clock::now();
        auto indexDuration = std::chrono::duration_cast<std::chrono::milliseconds>(indexLoadTime - startTime);
        logging::info("Index loaded in {}ms: {} nodes", indexDuration.count(), header->numNodes);
        
        // Check if ZSTD version exists for faster loading
        std::string zstdPath = panmanPath + ".zst";
        std::ifstream zstdFile(zstdPath, std::ios::binary);
        
        if (zstdFile) {
            zstdFile.close();
            logging::info("Using ZSTD-compressed panman: {}", zstdPath);
            
            std::vector<uint8_t> data;
            if (!panmap_zstd::decompressFromFile(zstdPath, data)) {
                logging::err("Failed to decompress ZSTD panman");
                return false;
            }
            
            std::string dataStr(data.begin(), data.end());
            std::istringstream iss(dataStr);
            treeGroup = std::make_unique<panmanUtils::TreeGroup>(iss);
            
        } else {
            // Fall back to XZ decompression
            logging::info("Using XZ-compressed panman: {}", panmanPath);
            
            std::ifstream file(panmanPath, std::ios::binary);
            if (!file) {
                logging::err("Cannot open panman: {}", panmanPath);
                return false;
            }
            
            auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
            buffer->push(boost::iostreams::lzma_decompressor());
            buffer->push(file);
            
            std::istream stream(buffer.get());
            treeGroup = std::make_unique<panmanUtils::TreeGroup>(stream);
        }
        
        if (treeGroup->trees.empty()) {
            logging::err("No trees in panman");
            return false;
        }
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        logging::info("Total load time: {}ms", totalDuration.count());
        
        return true;
    }
    
    int32_t findNodeIndex(const std::string& nodeId) const {
        auto it = nameCache.find(nodeId);
        if (it != nameCache.end()) {
            return it->second;
        }
        return -1;
    }
    
    std::vector<uint32_t> getPathIndices(const std::string& nodeId) const {
        std::vector<uint32_t> path;
        
        int32_t idx = findNodeIndex(nodeId);
        if (idx < 0) return path;
        
        while (idx >= 0) {
            path.push_back(idx);
            uint32_t parentIdx = nodeEntries[idx].parentIndex;
            if (parentIdx == UINT32_MAX) break;
            idx = parentIdx;
        }
        
        std::reverse(path.begin(), path.end());
        return path;
    }
};

GenomeExtractor::GenomeExtractor(
    const std::string& panmanPath,
    const std::string& indexPath
) : impl_(std::make_unique<Impl>()) {
    std::string idxPath = indexPath.empty() ? panmanPath + ".pmx" : indexPath;
    isOpen_ = impl_->open(panmanPath, idxPath);
}

GenomeExtractor::~GenomeExtractor() = default;

std::vector<std::string> GenomeExtractor::getNodeIds() const {
    std::vector<std::string> ids;
    if (!isOpen_) return ids;
    
    ids.reserve(impl_->header->numNodes);
    for (uint32_t i = 0; i < impl_->header->numNodes; i++) {
        ids.emplace_back(impl_->nameTable + impl_->nodeEntries[i].nameOffset,
                        impl_->nodeEntries[i].nameLength);
    }
    return ids;
}

bool GenomeExtractor::hasNode(const std::string& nodeId) const {
    if (!isOpen_) return false;
    return impl_->findNodeIndex(nodeId) >= 0;
}

std::vector<std::string> GenomeExtractor::getPathToNode(const std::string& nodeId) const {
    std::vector<std::string> path;
    if (!isOpen_) return path;
    
    auto indices = impl_->getPathIndices(nodeId);
    path.reserve(indices.size());
    for (uint32_t idx : indices) {
        path.emplace_back(impl_->nameTable + impl_->nodeEntries[idx].nameOffset,
                         impl_->nodeEntries[idx].nameLength);
    }
    return path;
}

std::string GenomeExtractor::extractGenome(const std::string& nodeId, bool aligned) {
    if (!isOpen_) return "";
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Use the correct panmanUtils implementation for genome reconstruction
    // This ensures 100% correctness with all edge cases handled:
    // - Block structure (primary/secondary blocks)
    // - Consensus sequence decoding (packed 4-bit nucleotides)
    // - All mutation types (NS, NI, ND, NSNPS, NSNPI, NSNPD)
    // - Gap handling (nucleotide gaps and block gaps)
    // - Strand inversions and reverse complement
    // - Rotation indexes for circular sequences
    // - Sequence inversion
    std::string result = impl_->treeGroup->trees[0].getStringFromReference(nodeId, aligned);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    stats_.extractionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    // Get path for stats
    auto pathIndices = impl_->getPathIndices(nodeId);
    stats_.nodesTraversed = pathIndices.size();
    
    // Count mutations along path (approximate)
    stats_.mutationsApplied = 0;
    for (uint32_t idx : pathIndices) {
        stats_.mutationsApplied += impl_->nodeEntries[idx].mutationCount;
    }
    
    if (result.find("Error:") == 0) {
        logging::err("Genome extraction failed: {}", result);
        return "";
    }
    
    logging::info("Extracted genome for {}: {} bp in {:.1f}ms",
                 nodeId, result.size(), stats_.extractionTimeMs);
    
    return result;
}

// ============================================================================
// Utility Functions
// ============================================================================

bool convertXzToZstd(
    const std::string& inputPath,
    const std::string& outputPath,
    int compressionLevel,
    size_t frameSize
) {
    std::string outPath = outputPath.empty() ? inputPath + ".zst" : outputPath;
    
    logging::info("Converting {} to ZSTD: {}", inputPath, outPath);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Decompress XZ
    std::ifstream file(inputPath, std::ios::binary);
    if (!file) {
        logging::err("Cannot open input file: {}", inputPath);
        return false;
    }
    
    auto buffer = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    buffer->push(boost::iostreams::lzma_decompressor());
    buffer->push(file);
    
    std::vector<char> data;
    {
        std::istream stream(buffer.get());
        std::ostringstream oss;
        oss << stream.rdbuf();
        std::string str = oss.str();
        data.assign(str.begin(), str.end());
    }
    
    auto decompressTime = std::chrono::high_resolution_clock::now();
    auto decompressDuration = std::chrono::duration_cast<std::chrono::milliseconds>(decompressTime - startTime);
    logging::info("XZ decompression: {} bytes in {}ms", data.size(), decompressDuration.count());
    
    // Compress with ZSTD
    bool success = panmap_zstd::compressToFile(data.data(), data.size(), outPath, 
                                               compressionLevel, 0, frameSize);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    if (success) {
        logging::info("Conversion complete in {}ms", totalDuration.count());
    }
    
    return success;
}

bool hasIndex(const std::string& panmanPath) {
    std::string indexPath = panmanPath + ".pmx";
    std::ifstream file(indexPath, std::ios::binary);
    if (!file) return false;
    
    uint32_t magic;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    return file.good() && magic == PMX_MAGIC;
}

std::string getIndexPath(const std::string& panmanPath) {
    return panmanPath + ".pmx";
}

} // namespace panman_index

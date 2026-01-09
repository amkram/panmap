#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <memory>

namespace panmap_zstd {

// Segment information for parallel compression/decompression
struct SegmentInfo {
    size_t offset;              // Byte offset in file
    size_t compressedSize;      // Compressed size in bytes
    size_t uncompressedSize;    // Uncompressed size in bytes
};

// Magic number for segmented ZSTD format: "PZST" (Panmap ZSTD)
constexpr uint32_t SEGMENTED_MAGIC = 0x54535A50;  // "PZST" in little-endian

/**
 * @brief Compress multiple data segments to file with parallel ZSTD compression
 * 
 * Each segment is compressed independently, allowing for parallel decompression.
 * File format: [Header][SegmentTable][Segment0][Segment1]...[SegmentN]
 * 
 * @param segments Vector of (data pointer, size) pairs
 * @param outputPath Output file path (.idx extension)
 * @param compressionLevel ZSTD compression level (1-22, default 3 for speed)
 * @param numThreads Number of compression threads (0 = auto-detect)
 * @return true on success
 */
bool compressSegmentedToFile(
    const std::vector<std::pair<const void*, size_t>>& segments,
    const std::string& outputPath,
    int compressionLevel = 3,
    int numThreads = 0
);

/**
 * @brief Decompress segmented ZSTD file with parallel decompression
 * 
 * Decompresses all segments in parallel using TBB.
 * 
 * @param inputPath Input file path
 * @param outputData Output buffer (will be resized to total uncompressed size)
 * @param numThreads Number of decompression threads (0 = auto-detect)
 * @return true on success
 */
bool decompressSegmentedFromFile(
    const std::string& inputPath,
    std::vector<uint8_t>& outputData,
    int numThreads = 0
);

/**
 * @brief Decompress segmented ZSTD file with parallel decompression into raw buffer
 * 
 * Allocates buffer using new[] to avoid zero-initialization overhead.
 * 
 * @param inputPath Input file path
 * @param outputData Reference to unique_ptr that will hold the allocated data
 * @param outputSize Reference to size_t that will hold the size
 * @param numThreads Number of decompression threads (0 = auto-detect)
 * @return true on success
 */
bool decompressSegmentedFromFile(
    const std::string& inputPath,
    std::unique_ptr<uint8_t[]>& outputData,
    size_t& outputSize,
    int numThreads = 0
);

/**
 * @brief Check if file is segmented ZSTD format
 * 
 * @param inputPath Input file path
 * @return true if file starts with segmented magic number
 */
bool isSegmentedZstd(const std::string& inputPath);

/**
 * @brief Compress Cap'n Proto message to file with parallel ZSTD compression
 * 
 * DEPRECATED: Use compressSegmentedToFile for better performance.
 * Kept for single-segment use cases.
 * 
 * @param inputData Raw message data
 * @param inputSize Size of input data
 * @param outputPath Output file path (.idx extension)
 * @param compressionLevel ZSTD compression level (1-22, default 3 for speed)
 * @param numThreads Number of compression threads (0 = auto-detect)
 * @param frameSize Size of each seekable frame in bytes (default 4MB for good parallelism)
 * @return true on success
 */
bool compressToFile(
    const void* inputData,
    size_t inputSize,
    const std::string& outputPath,
    int compressionLevel = 3,
    int numThreads = 0,
    size_t frameSize = 4 * 1024 * 1024  // 4MB frames
);

/**
 * @brief Decompress ZSTD file with parallel decompression
 * 
 * DEPRECATED: Use decompressSegmentedFromFile for better performance.
 * 
 * @param inputPath Input file path
 * @param outputData Output buffer (will be resized)
 * @param numThreads Number of decompression threads (0 = auto-detect)
 * @return true on success
 */
bool decompressFromFile(
    const std::string& inputPath,
    std::vector<uint8_t>& outputData,
    int numThreads = 0
);

/**
 * @brief Get uncompressed size without full decompression
 * 
 * @param inputPath Input file path
 * @return Uncompressed size in bytes, or 0 on error
 */
size_t getUncompressedSize(const std::string& inputPath);

} // namespace panmap_zstd

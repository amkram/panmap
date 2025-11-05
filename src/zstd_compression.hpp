#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <memory>

namespace panmap_zstd {

/**
 * @brief Compress Cap'n Proto message to file with parallel ZSTD compression
 * 
 * Uses ZSTD's multi-threaded compression with seekable format for fast parallel decompression.
 * The seekable format splits data into frames that can be decompressed independently.
 * 
 * @param inputData Raw message data
 * @param inputSize Size of input data
 * @param outputPath Output file path (.pmi extension)
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
 * Uses seekable format for parallel frame decompression.
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
 * Reads frame headers to determine original size.
 * 
 * @param inputPath Input file path
 * @return Uncompressed size in bytes, or 0 on error
 */
size_t getUncompressedSize(const std::string& inputPath);

/**
 * @brief Check if file is ZSTD compressed
 * 
 * @param inputPath Input file path
 * @return true if file starts with ZSTD magic number
 */
bool isZstdCompressed(const std::string& inputPath);

} // namespace panmap_zstd

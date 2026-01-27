#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace panmap_zstd {

/**
 * @brief Compress data to file with parallel ZSTD compression
 * 
 * @param inputData Raw data
 * @param inputSize Size of input data
 * @param outputPath Output file path
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

} // namespace panmap_zstd

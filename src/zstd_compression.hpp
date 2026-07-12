#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace panmap_zstd {

/**
 * Compress inputData to outputPath with parallel ZSTD.
 * compressionLevel: 1-22. numThreads: 0 = auto-detect. frameSize: seekable
 * frame size in bytes. Returns true on success.
 */
[[nodiscard]] bool compressToFile(const void* inputData,
                                  size_t inputSize,
                                  const std::string& outputPath,
                                  int compressionLevel = 3,
                                  int numThreads = 0,
                                  size_t frameSize = 4 * 1024 * 1024,  // 4MB frames
                                  const void* headerData = nullptr,    // optional uncompressed prefix,
                                  size_t headerSize = 0                // written verbatim before the frames
);

/**
 * Decompress inputPath into outputData (resized to fit) with parallel ZSTD.
 * numThreads: 0 = auto-detect. dataOffset: skip this many leading bytes (an
 * uncompressed header) before the ZSTD frames. Returns true on success.
 */
[[nodiscard]] bool
decompressFromFile(const std::string& inputPath, std::vector<uint8_t>& outputData, int numThreads = 0,
                   size_t dataOffset = 0);

}  // namespace panmap_zstd

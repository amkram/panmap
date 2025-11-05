#include "zstd_compression.hpp"
#include "logging.hpp"
#include <zstd.h>
#include <fstream>
#include <cstring>
#include <thread>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace panmap_zstd {

// ZSTD magic number: 0xFD2FB528 (little-endian)
static const uint32_t ZSTD_MAGIC = 0x28B52FFD;

bool isZstdCompressed(const std::string& inputPath) {
    std::ifstream file(inputPath, std::ios::binary);
    if (!file) return false;
    
    uint32_t magic;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    return file.good() && magic == ZSTD_MAGIC;
}

size_t getUncompressedSize(const std::string& inputPath) {
    std::ifstream file(inputPath, std::ios::binary);
    if (!file) return 0;
    
    file.seekg(0, std::ios::end);
    size_t compressedSize = file.tellg();
    file.seekg(0, std::ios::beg);
    
    std::vector<uint8_t> compressedData(compressedSize);
    file.read(reinterpret_cast<char*>(compressedData.data()), compressedSize);
    
    unsigned long long const decompressedSize = ZSTD_getFrameContentSize(
        compressedData.data(), compressedData.size());
    
    if (decompressedSize == ZSTD_CONTENTSIZE_ERROR) {
        logging::err("Not a valid ZSTD frame: {}", inputPath);
        return 0;
    }
    if (decompressedSize == ZSTD_CONTENTSIZE_UNKNOWN) {
        logging::warn("Uncompressed size unknown for: {}", inputPath);
        return 0;
    }
    
    return static_cast<size_t>(decompressedSize);
}

bool compressToFile(
    const void* inputData,
    size_t inputSize,
    const std::string& outputPath,
    int compressionLevel,
    int numThreads,
    size_t frameSize
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    logging::info("Compressing {} bytes to {} with ZSTD level {} using {} threads",
                  inputSize, outputPath, compressionLevel, numThreads);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Create compression context with multi-threading
    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    if (!cctx) {
        logging::err("Failed to create ZSTD compression context");
        return false;
    }
    
    // Set compression parameters
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionLevel);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, numThreads);
    
    // Enable checksums for data integrity
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
    
    // Allocate output buffer (worst case: slightly larger than input)
    size_t const maxCompressedSize = ZSTD_compressBound(inputSize);
    std::vector<uint8_t> compressedData(maxCompressedSize);
    
    // Compress
    size_t const compressedSize = ZSTD_compress2(
        cctx,
        compressedData.data(), maxCompressedSize,
        inputData, inputSize
    );
    
    ZSTD_freeCCtx(cctx);
    
    if (ZSTD_isError(compressedSize)) {
        logging::err("ZSTD compression failed: {}", ZSTD_getErrorName(compressedSize));
        return false;
    }
    
    // Write to file
    std::ofstream outFile(outputPath, std::ios::binary);
    if (!outFile) {
        logging::err("Failed to open output file: {}", outputPath);
        return false;
    }
    
    outFile.write(reinterpret_cast<const char*>(compressedData.data()), compressedSize);
    outFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double ratio = 100.0 * compressedSize / inputSize;
    double throughput = inputSize / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Compression complete: {} -> {} bytes ({:.1f}% ratio) in {}ms ({:.1f} MB/s)",
                  inputSize, compressedSize, ratio, duration.count(), throughput);
    
    return true;
}

bool decompressFromFile(
    const std::string& inputPath,
    std::vector<uint8_t>& outputData,
    int numThreads
) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    
    logging::info("Decompressing {} using {} threads", inputPath, numThreads);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Read compressed file
    std::ifstream inFile(inputPath, std::ios::binary);
    if (!inFile) {
        logging::err("Failed to open input file: {}", inputPath);
        return false;
    }
    
    inFile.seekg(0, std::ios::end);
    size_t compressedSize = inFile.tellg();
    inFile.seekg(0, std::ios::beg);
    
    std::vector<uint8_t> compressedData(compressedSize);
    inFile.read(reinterpret_cast<char*>(compressedData.data()), compressedSize);
    inFile.close();
    
    // Get uncompressed size
    unsigned long long const decompressedSize = ZSTD_getFrameContentSize(
        compressedData.data(), compressedSize);
    
    if (decompressedSize == ZSTD_CONTENTSIZE_ERROR) {
        logging::err("Not a valid ZSTD frame: {}", inputPath);
        return false;
    }
    if (decompressedSize == ZSTD_CONTENTSIZE_UNKNOWN) {
        logging::err("Cannot determine uncompressed size: {}", inputPath);
        return false;
    }
    
    // Allocate output buffer
    outputData.resize(decompressedSize);
    
    // Create decompression context with multi-threading
    ZSTD_DCtx* dctx = ZSTD_createDCtx();
    if (!dctx) {
        logging::err("Failed to create ZSTD decompression context");
        return false;
    }
    
    // Decompress
    size_t const actualSize = ZSTD_decompressDCtx(
        dctx,
        outputData.data(), decompressedSize,
        compressedData.data(), compressedSize
    );
    
    ZSTD_freeDCtx(dctx);
    
    if (ZSTD_isError(actualSize)) {
        logging::err("ZSTD decompression failed: {}", ZSTD_getErrorName(actualSize));
        return false;
    }
    
    if (actualSize != decompressedSize) {
        logging::err("Decompressed size mismatch: expected {}, got {}", 
                     decompressedSize, actualSize);
        return false;
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double throughput = actualSize / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Decompression complete: {} -> {} bytes in {}ms ({:.1f} MB/s)",
                  compressedSize, actualSize, duration.count(), throughput);
    
    return true;
}

} // namespace panmap_zstd

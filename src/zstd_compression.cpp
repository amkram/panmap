#include "zstd_compression.hpp"
#include "logging.hpp"
#include <zstd.h>
#include <fstream>
#include <thread>
#include <tbb/parallel_for.h>
#include <boost/iostreams/device/mapped_file.hpp>

namespace panmap_zstd {

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
    
    logging::info("Compressing {} bytes to {} with ZSTD level {} using {} threads (frame size: {} MB)",
                  inputSize, outputPath, compressionLevel, numThreads, frameSize / (1024 * 1024));
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Open output file
    std::ofstream outFile(outputPath, std::ios::binary);
    if (!outFile) {
        logging::err("Failed to open output file: {}", outputPath);
        return false;
    }
    
    // Create compression context with multi-threading
    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    if (!cctx) {
        logging::err("Failed to create ZSTD compression context");
        return false;
    }
    
    // Set compression parameters for parallel compression
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionLevel);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, numThreads);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_jobSize, frameSize);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
    
    // Allocate output buffer
    size_t const maxCompressedSize = ZSTD_compressBound(inputSize);
    std::vector<uint8_t> compressedData(maxCompressedSize);
    
    // Compress with parallel frames
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
    outFile.write(reinterpret_cast<const char*>(compressedData.data()), compressedSize);
    if (!outFile) {
        logging::err("Failed to write compressed data to file: {}", outputPath);
        return false;
    }
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
    
    // Use mmap for input
    boost::iostreams::mapped_file_source mappedFile;
    try {
        mappedFile.open(inputPath);
    } catch (const std::exception& e) {
        logging::err("Failed to memory map input file: {} ({})", inputPath, e.what());
        return false;
    }
    
    if (!mappedFile.is_open()) {
        logging::err("Failed to open input file: {}", inputPath);
        return false;
    }
    
    const char* compressedData = mappedFile.data();
    size_t compressedSize = mappedFile.size();
    
    // Get uncompressed size
    unsigned long long const decompressedSize = ZSTD_getFrameContentSize(
        compressedData, compressedSize);
    
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
    
    // Check if the compressed data has multiple frames
    size_t numFrames = 0;
    size_t pos = 0;
    while (pos < compressedSize) {
        unsigned long long frameSize = ZSTD_findFrameCompressedSize(compressedData + pos, compressedSize - pos);
        if (ZSTD_isError(frameSize)) {
            break;
        }
        numFrames++;
        pos += frameSize;
    }
    
    if (numFrames > 1 && numThreads > 1) {
        // Parallel decompression of multiple frames
        logging::debug("Detected {} ZSTD frames, decompressing in parallel", numFrames);
        
        std::vector<std::pair<size_t, size_t>> frameRanges;
        std::vector<std::pair<size_t, size_t>> outputRanges;
        
        pos = 0;
        size_t outPos = 0;
        while (pos < compressedSize) {
            unsigned long long frameCompSize = ZSTD_findFrameCompressedSize(compressedData + pos, compressedSize - pos);
            if (ZSTD_isError(frameCompSize)) break;
            
            unsigned long long frameDecompSize = ZSTD_getFrameContentSize(compressedData + pos, frameCompSize);
            if (ZSTD_isError(frameDecompSize) || frameDecompSize == ZSTD_CONTENTSIZE_UNKNOWN) break;
            
            frameRanges.push_back({pos, frameCompSize});
            outputRanges.push_back({outPos, frameDecompSize});
            
            pos += frameCompSize;
            outPos += frameDecompSize;
        }
        
        // Decompress frames in parallel
        std::atomic<bool> error{false};
        tbb::parallel_for(size_t(0), frameRanges.size(), [&](size_t i) {
            if (error.load()) return;
            
            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            if (!dctx) {
                error.store(true);
                return;
            }
            
            size_t result = ZSTD_decompressDCtx(
                dctx,
                outputData.data() + outputRanges[i].first, outputRanges[i].second,
                compressedData + frameRanges[i].first, frameRanges[i].second
            );
            
            ZSTD_freeDCtx(dctx);
            
            if (ZSTD_isError(result) || result != outputRanges[i].second) {
                error.store(true);
            }
        });
        
        if (error.load()) {
            logging::err("Parallel ZSTD decompression failed");
            mappedFile.close();
            return false;
        }
        
    } else {
        // Single-threaded decompression
        logging::debug("Using single-threaded decompression");
        
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        if (!dctx) {
            logging::err("Failed to create ZSTD decompression context");
            return false;
        }
        
        size_t const actualSize = ZSTD_decompressDCtx(
            dctx,
            outputData.data(), decompressedSize,
            compressedData, compressedSize
        );
        
        ZSTD_freeDCtx(dctx);
        
        if (ZSTD_isError(actualSize)) {
            logging::err("ZSTD decompression failed: {}", ZSTD_getErrorName(actualSize));
            mappedFile.close();
            return false;
        }
        
        if (actualSize != decompressedSize) {
            logging::err("Decompressed size mismatch: expected {}, got {}", 
                         decompressedSize, actualSize);
            mappedFile.close();
            return false;
        }
    }
    
    mappedFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    double throughput = decompressedSize / (1024.0 * 1024.0 * duration.count() / 1000.0);
    
    logging::info("Decompression complete: {} -> {} bytes in {}ms ({:.1f} MB/s)",
                  compressedSize, decompressedSize, duration.count(), throughput);
    
    return true;
}

} // namespace panmap_zstd

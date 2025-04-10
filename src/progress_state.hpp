#pragma once

#include <string>
#include <mutex>
#include <chrono>

namespace placement {

/**
 * Structure used to track and report progress during placement
 */
struct PlacementProgressState {
    // Mutex to protect concurrent access
    std::mutex mtx;
    
    // Status flags
    bool running = false;
    bool complete = false;
    
    // Progress tracking
    size_t nodesProcessed = 0;
    size_t totalNodes = 0;
    double progress = 0.0;
    
    // Current state
    std::string currentNodeId;
    
    // Best scores found so far
    int64_t bestHitCount = 0;
    std::string bestHitNode;
    
    double bestJaccardScore = 0.0;
    std::string bestJaccardNode;
    
    double bestCosineScore = 0.0;
    std::string bestCosineNode;
    
    // Timing information
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> lastUpdateTime;
    double elapsedTime = 0.0; // Total elapsed time in seconds
    
    // Result message
    std::string resultMessage;
};

} // namespace placement 
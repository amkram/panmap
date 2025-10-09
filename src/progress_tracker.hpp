#pragma once
#include <atomic>
#include <vector>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <mutex>

class ProgressTracker {
private:
  std::vector<std::atomic<uint64_t>> threadProgress;
  std::vector<uint64_t> threadTotals;
  mutable std::mutex outputMutex;
  size_t numThreads;
  std::chrono::steady_clock::time_point startTime;
  mutable std::chrono::steady_clock::time_point lastGlobalUpdate;
  
public:
  ProgressTracker(size_t numThreads, const std::vector<uint64_t>& totalNodesPerThread) 
    : numThreads(numThreads), threadProgress(numThreads), threadTotals(totalNodesPerThread) {
    
    for (size_t i = 0; i < numThreads; ++i) {
      threadProgress[i].store(0, std::memory_order_relaxed);
    }
    startTime = std::chrono::steady_clock::now();
    lastGlobalUpdate = startTime;
    
    // Reserve space for all thread lines
    for (size_t i = 0; i < numThreads; ++i) {
      std::cout << "T" << std::setw(2) << i << ": [  0.00%]        0/" 
                << std::setw(8) << threadTotals[i] << " (   0.0 n/s)\n";
    }
  }
  
  void updateProgress(size_t threadId, uint64_t currentNode) {
    threadProgress[threadId].store(currentNode, std::memory_order_relaxed);
    
    // Throttle updates per thread
    static thread_local auto lastUpdate = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    
    if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count() >= 100) {
      tryDisplayUpdate();
      lastUpdate = now;
    }
  }
  
  void incrementProgress(size_t threadId) {
    threadProgress[threadId].fetch_add(1, std::memory_order_relaxed);
    
    // Less frequent updates for increment calls
    static thread_local auto lastUpdate = std::chrono::steady_clock::now();
    static thread_local int updateCounter = 0;
    
    if (++updateCounter % 50 == 0) {  // Update every 50 increments
      auto now = std::chrono::steady_clock::now();
      if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count() >= 200) {
        tryDisplayUpdate();
        lastUpdate = now;
      }
    }
  }
  
private:
  void tryDisplayUpdate() const {
    std::unique_lock lock(outputMutex, std::try_to_lock);
    if (!lock.owns_lock()) return;
    
    auto now = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastGlobalUpdate).count() < 1000) {
      return;
    }
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
    
    // Save cursor position and clear from cursor to end of screen
    std::cout << "\033[s\033[J";
    
    for (size_t i = 0; i < numThreads; ++i) {
      uint64_t threadProg = threadProgress[i].load(std::memory_order_relaxed);
      uint64_t threadTotal = threadTotals[i];
      
      double percentComplete = threadTotal > 0 ? (100.0 * threadProg) / threadTotal : 0.0;
      double nodesPerSecond = elapsed > 0 ? static_cast<double>(threadProg) / elapsed : 0.0;
      
      std::cout << "T" << std::setw(2) << i << ": ["
                << std::setw(6) << std::fixed << std::setprecision(2) << percentComplete << "%] "
                << std::setw(8) << threadProg << "/" << std::setw(8) << threadTotal
                << " (" << std::setw(6) << std::setprecision(1) << nodesPerSecond << " n/s)\n";
    }
    
    // Restore cursor position
    std::cout << "\033[u" << std::flush;
    
    lastGlobalUpdate = now;
  }

public:
  void finalDisplay() {
    std::lock_guard<std::mutex> lock(outputMutex);
    
    // Move cursor down past all thread lines
    for (size_t i = 0; i < numThreads; ++i) {
      std::cout << "\n";
    }
    
    // Calculate totals
    uint64_t totalProgress = 0;
    uint64_t totalNodes = 0;
    
    for (size_t i = 0; i < numThreads; ++i) {
      totalProgress += threadProgress[i].load(std::memory_order_relaxed);
      totalNodes += threadTotals[i];
    }
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - startTime).count();
    double nodesPerSecond = elapsed > 0 ? static_cast<double>(totalProgress) / elapsed : 0.0;
    
    std::cout << "Completed: " << totalProgress << "/" << totalNodes 
              << " nodes in " << elapsed << "s (" 
              << std::setprecision(1) << nodesPerSecond << " nodes/s)" << std::endl;
  }
};

// Simpler alternative that logs to separate streams
class LoggingProgressTracker {
private:
  std::vector<std::atomic<uint64_t>> threadProgress;
  std::vector<uint64_t> threadTotals;
  std::mutex logMutex;
  size_t numThreads;
  std::chrono::steady_clock::time_point startTime;
  
public:
  LoggingProgressTracker(size_t numThreads, const std::vector<uint64_t>& totalNodesPerThread) 
    : numThreads(numThreads), threadProgress(numThreads), threadTotals(totalNodesPerThread) {
    
    for (size_t i = 0; i < numThreads; ++i) {
      threadProgress[i].store(0, std::memory_order_relaxed);
    }
    startTime = std::chrono::steady_clock::now();
  }
  
  void incrementProgress(size_t threadId) {
    uint64_t newProgress = threadProgress[threadId].fetch_add(1, std::memory_order_relaxed) + 1;
    
    // Log every 1000 nodes to reduce output volume
    if (newProgress % 1000 == 0) {
      std::lock_guard<std::mutex> lock(logMutex);
      
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - startTime).count();
      double percent = threadTotals[threadId] > 0 ? (100.0 * newProgress) / threadTotals[threadId] : 0.0;
      double rate = elapsed > 0 ? static_cast<double>(newProgress) / elapsed : 0.0;
      
      std::cout << "Thread " << threadId << ": " << std::setprecision(2) << std::fixed 
                << percent << "% (" << newProgress << "/" << threadTotals[threadId] 
                << ") " << std::setprecision(1) << rate << " nodes/s" << std::endl;
    }
  }
  
  void finalDisplay() {
    std::lock_guard<std::mutex> lock(logMutex);
    
    uint64_t totalProgress = 0;
    uint64_t totalNodes = 0;
    
    for (size_t i = 0; i < numThreads; ++i) {
      uint64_t progress = threadProgress[i].load(std::memory_order_relaxed);
      totalProgress += progress;
      totalNodes += threadTotals[i];
      
      double percent = threadTotals[i] > 0 ? (100.0 * progress) / threadTotals[i] : 0.0;
      std::cout << "Thread " << i << " final: " << std::setprecision(2) << std::fixed 
                << percent << "% (" << progress << "/" << threadTotals[i] << ")" << std::endl;
    }
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - startTime).count();
    double nodesPerSecond = elapsed > 0 ? static_cast<double>(totalProgress) / elapsed : 0.0;
    
    std::cout << "Total completed: " << totalProgress << "/" << totalNodes 
              << " nodes in " << elapsed << "s (" 
              << std::setprecision(1) << nodesPerSecond << " nodes/s)" << std::endl;
  }
};
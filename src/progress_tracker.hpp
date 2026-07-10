#pragma once
#include <atomic>
#include <vector>
#include <chrono>
#include <cstdint>
#include <memory>
#include <numeric>
#include "logging.hpp"

// Aggregates per-thread progress into a single progress bar.
class ProgressTracker {
   public:
    ProgressTracker(size_t numThreads, const std::vector<uint64_t>& totalNodesPerThread)
        : threadProgress(numThreads), numThreads_(numThreads) {
        for (size_t i = 0; i < numThreads; ++i) {
            threadProgress[i].store(0, std::memory_order_relaxed);
        }
        uint64_t total = 0;
        for (auto t : totalNodesPerThread) total += t;
        bar_ = std::make_unique<output::ProgressBar>("place", total);
    }

    void incrementProgress(size_t threadId) {
        threadProgress[threadId].fetch_add(1, std::memory_order_relaxed);
        publish();
    }

    // Flush final progress, then clear the bar.
    void finalDisplay() {
        if (!bar_) return;
        uint64_t totalProgress = 0;
        for (size_t i = 0; i < numThreads_; ++i) {
            totalProgress += threadProgress[i].load(std::memory_order_relaxed);
        }
        bar_->set(totalProgress);
        bar_->clear();
        // Don't emit an action line here; the caller emits the "Placed" line.
        bar_.reset();
    }

   private:
    void publish() {
        if (!bar_) return;
        uint64_t total = 0;
        for (size_t i = 0; i < numThreads_; ++i) {
            total += threadProgress[i].load(std::memory_order_relaxed);
        }
        bar_->set(total);
    }

    std::vector<std::atomic<uint64_t>> threadProgress;
    size_t numThreads_;
    std::unique_ptr<output::ProgressBar> bar_;
};

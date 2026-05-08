#pragma once
#include <atomic>
#include <vector>
#include <chrono>
#include <cstdint>
#include <memory>
#include <numeric>
#include "logging.hpp"

// Aggregates per-thread progress into a single sleek progress bar.
// Public API matches the previous per-thread tracker so callers don't change.
class ProgressTracker {
   public:
    ProgressTracker(size_t numThreads, const std::vector<uint64_t>& totalNodesPerThread)
        : threadProgress(numThreads),
          threadTotals(totalNodesPerThread),
          numThreads_(numThreads) {
        for (size_t i = 0; i < numThreads; ++i) {
            threadProgress[i].store(0, std::memory_order_relaxed);
        }
        uint64_t total = 0;
        for (auto t : totalNodesPerThread) total += t;
        bar_ = std::make_unique<output::ProgressBar>("place", total);
    }

    void updateProgress(size_t threadId, uint64_t currentNode) {
        threadProgress[threadId].store(currentNode, std::memory_order_relaxed);
        publish();
    }

    void incrementProgress(size_t threadId) {
        threadProgress[threadId].fetch_add(1, std::memory_order_relaxed);
        publish();
    }

    // Replace the bar with a final action line.
    void finalDisplay() {
        if (!bar_) return;
        uint64_t totalProgress = 0;
        uint64_t totalNodes = 0;
        for (size_t i = 0; i < numThreads_; ++i) {
            totalProgress += threadProgress[i].load(std::memory_order_relaxed);
            totalNodes += threadTotals[i];
        }
        bar_->set(totalProgress);
        bar_->clear();
        // Don't emit an action line here; the caller emits the canonical "Placed" line.
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
    std::vector<uint64_t> threadTotals;
    size_t numThreads_;
    std::unique_ptr<output::ProgressBar> bar_;
};

#pragma once
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace timing {

class Timer {
private:
  struct TimingStats {
    double total_time{0.0};
    double min_time{std::numeric_limits<double>::max()};
    double max_time{0.0};
    int calls{0};
    std::chrono::high_resolution_clock::time_point last_start;
  };

  // Thread-local timing data to reduce mutex contention
  static inline thread_local std::stack<std::pair<std::string,
                           std::chrono::high_resolution_clock::time_point>>
      local_timing_stack{};
  static inline thread_local std::unordered_map<std::string, TimingStats> 
      local_function_stats{};
  
  // Global aggregation data (protected by mutex)
  static inline std::mutex timing_mutex{};
  static inline std::unordered_map<std::string, TimingStats> global_function_stats{};
  static inline double min_time_to_report{0.0001}; // 100 microseconds
  
  // Aggregate thread-local stats to global (call periodically)
  static void aggregateLocalStats() {
    if (local_function_stats.empty()) return;
    
    std::lock_guard<std::mutex> lock(timing_mutex);
    for (const auto& [func_name, local_stats] : local_function_stats) {
      auto& global_stats = global_function_stats[func_name];
      global_stats.total_time += local_stats.total_time;
      global_stats.calls += local_stats.calls;
      global_stats.min_time = std::min(global_stats.min_time, local_stats.min_time);
      global_stats.max_time = std::max(global_stats.max_time, local_stats.max_time);
    }
    local_function_stats.clear(); // Clear after aggregation
  }

public:
  static void start(const std::string &function_name) {
    // Use thread-local stack - no mutex needed
    local_timing_stack.push({function_name, std::chrono::high_resolution_clock::now()});
  }

  static void startBlock(const std::string &block_name,
                         const std::string &function_name) {
    std::string full_name = function_name + "::" + block_name;
    start(full_name);
  }

  static void stop() {
    // Use thread-local stack - no mutex needed
    if (local_timing_stack.empty()) return;
    
    auto [function_name, start_time] = local_timing_stack.top();
    local_timing_stack.pop();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        end_time - start_time)
                        .count() /
                    1000000.0;

    // Update thread-local stats
    auto& stats = local_function_stats[function_name];
    stats.total_time += duration;
    stats.calls++;
    stats.min_time = std::min(stats.min_time, duration);
    stats.max_time = std::max(stats.max_time, duration);
    
    // Periodically aggregate to reduce global mutex contention
    static thread_local int aggregation_counter = 0;
    if (++aggregation_counter >= 100) { // Aggregate every 100 timings
      aggregateLocalStats();
      aggregation_counter = 0;
    }
  }

  static void setMinTimeToReport(double seconds) {
    std::lock_guard<std::mutex> lock(timing_mutex);
    min_time_to_report = seconds;
  }

  static std::string formatTime(double seconds) {
    std::stringstream ss;
    if (seconds < 0.000001) {
      ss << std::fixed << std::setprecision(2) << seconds * 1e9 << " ns";
    } else if (seconds < 0.001) {
      ss << std::fixed << std::setprecision(2) << seconds * 1e6 << " µs";
    } else if (seconds < 1.0) {
      ss << std::fixed << std::setprecision(2) << seconds * 1e3 << " ms";
    } else {
      ss << std::fixed << std::setprecision(2) << seconds << " s";
    }
    return ss.str();
  }

  static void report() {
    // First aggregate any remaining thread-local stats
    aggregateLocalStats();
    
    std::lock_guard<std::mutex> lock(timing_mutex);

    std::cout << "\n=== Timing Report ===\n\n";
    std::cout << std::left << std::setw(40) << "Function" << std::right
              << std::setw(10) << "Calls" << std::right << std::setw(15)
              << "Total" << std::right << std::setw(15) << "Average"
              << std::right << std::setw(15) << "Min" << std::right
              << std::setw(15) << "Max"
              << "\n";
    std::cout << std::string(110, '-') << "\n";

    // Sort by total time
    std::vector<std::pair<std::string, TimingStats>> sorted_stats(
        global_function_stats.begin(), global_function_stats.end());
    std::sort(sorted_stats.begin(), sorted_stats.end(),
              [](const auto &a, const auto &b) {
                return a.second.total_time > b.second.total_time;
              });

    for (const auto &[func, stats] : sorted_stats) {
      if (stats.total_time < min_time_to_report)
        continue;

      double avg_time = stats.total_time / stats.calls;

      std::cout << std::left << std::setw(40) << func << std::right
                << std::setw(10) << stats.calls << std::right << std::setw(15)
                << formatTime(stats.total_time) << std::right << std::setw(15)
                << formatTime(avg_time) << std::right << std::setw(15)
                << formatTime(stats.min_time) << std::right << std::setw(15)
                << formatTime(stats.max_time) << "\n";
    }

    std::cout << "\nNote: Only showing functions that took longer than "
              << formatTime(min_time_to_report) << " total time\n";
  }

  static void reset() {
    // First aggregate any remaining thread-local stats
    aggregateLocalStats();
    
    std::lock_guard<std::mutex> lock(timing_mutex);
    global_function_stats.clear();
    
    // Also clear thread-local data
    local_timing_stack = {};
    local_function_stats.clear();
  }
};

} // namespace timing

// Macro for easy function timing
#define TIME_FUNCTION                                                          \
  timing::Timer::start(__func__);                                              \
  auto timer_guard = std::unique_ptr<void, std::function<void(void *)>>(       \
      (void *)1, [](void *) { timing::Timer::stop(); });

// Macro for timing blocks within functions
#define TIME_BLOCK(block_name)                                                 \
  timing::Timer::startBlock(block_name, __func__);                             \
  auto timer_block_guard = std::unique_ptr<void, std::function<void(void *)>>( \
      (void *)1, [](void *) { timing::Timer::stop(); });
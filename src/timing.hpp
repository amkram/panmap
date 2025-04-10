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

  static inline std::mutex timing_mutex{};
  static inline std::unordered_map<
      std::thread::id,
      std::stack<std::pair<std::string,
                           std::chrono::high_resolution_clock::time_point>>>
      timing_stacks{};
  static inline std::unordered_map<std::string, TimingStats> function_stats{};
  static inline double min_time_to_report{0.0001}; // 100 microseconds

public:
  static void start(const std::string &function_name) {
    std::lock_guard<std::mutex> lock(timing_mutex);
    auto &stack = timing_stacks[std::this_thread::get_id()];
    stack.push({function_name, std::chrono::high_resolution_clock::now()});
  }

  static void startBlock(const std::string &block_name,
                         const std::string &function_name) {
    std::string full_name = function_name + "::" + block_name;
    start(full_name);
  }

  static void stop() {
    std::lock_guard<std::mutex> lock(timing_mutex);
    auto &stack = timing_stacks[std::this_thread::get_id()];
    if (stack.empty())
      return;

    auto [function_name, start_time] = stack.top();
    stack.pop();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                        end_time - start_time)
                        .count() /
                    1000000.0;

    auto &stats = function_stats[function_name];
    stats.total_time += duration;
    stats.min_time = std::min(stats.min_time, duration);
    stats.max_time = std::max(stats.max_time, duration);
    stats.calls++;
  }

  static void setMinTimeToReport(double seconds) {
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
    std::lock_guard<std::mutex> lock(timing_mutex);

    // Check for any unclosed timers
    for (const auto &[thread_id, stack] : timing_stacks) {
      if (!stack.empty()) {
        std::cerr << "Warning: Found " << stack.size()
                  << " unclosed timer(s) in thread " << thread_id << "\n";
      }
    }

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
        function_stats.begin(), function_stats.end());
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
    std::lock_guard<std::mutex> lock(timing_mutex);
    timing_stacks.clear();
    function_stats.clear();
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
#pragma once

// Force-enable spdlog since we've added it as a dependency
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#define SPDLOG_AVAILABLE 1

#include <chrono>
#include <iostream>
#include <string>
#include <fstream>
#include <mutex>
#include <vector>
#include <sstream>

namespace logging {

// Logging verbosity levels
enum class LogLevel {
  TRACE,
  DEBUG,
  INFO,
  WARN,
  ERROR,
  CRITICAL,
  OFF,
  QUIET,     // Restored
  NORMAL,    // Restored
  VERBOSE    // Restored
};

// Global verbosity level
static LogLevel loggingLevel = LogLevel::INFO;
static bool spdlogInitialized = false;

// Node tracking file for node_1 and node_2
static std::ofstream nodeTrackingFile;
static std::mutex nodeTrackingMutex;
static bool nodeTrackingInitialized = false;

// --- Update helper to map LogLevel to spdlog::level::level_enum ---
inline spdlog::level::level_enum toSpdlogLevel(LogLevel customLevel) {
    switch (customLevel) {
        case LogLevel::TRACE:   return spdlog::level::trace;
        case LogLevel::DEBUG:   return spdlog::level::debug;
        case LogLevel::INFO:    return spdlog::level::info;
        case LogLevel::WARN:    return spdlog::level::warn;
        case LogLevel::ERROR:   return spdlog::level::err;
        case LogLevel::CRITICAL:return spdlog::level::critical;
        case LogLevel::OFF:     return spdlog::level::off;
        // Handle QUIET, NORMAL, VERBOSE 
        case LogLevel::QUIET:   return spdlog::level::critical; // Map QUIET to critical
        case LogLevel::NORMAL:  return spdlog::level::info;     // Map NORMAL to info
        case LogLevel::VERBOSE: return spdlog::level::trace;    // Map VERBOSE to trace
        default:                return spdlog::level::info; // Default safety
    }
}
// --- End helper --- 

// Initialize spdlog
inline void initSpdlog() {
    if (!spdlogInitialized) {
        // Set pattern with timestamps, level, and message
        spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e [%^%l%$] %v");
        
        // Create default logger if none exists
        if (!spdlog::default_logger()) {
            auto logger = spdlog::stdout_color_mt("console");
            spdlog::set_default_logger(logger);
        }
        
        // ---> Set spdlog's level based on our variable <--- 
        spdlog::set_level(toSpdlogLevel(loggingLevel));
        spdlogInitialized = true;
    }
}

// Initialize node tracking log file
inline void initNodeTracking() {
    std::lock_guard<std::mutex> lock(nodeTrackingMutex);
    if (!nodeTrackingInitialized) {
        nodeTrackingFile.open("node_tracking.log", std::ios::out | std::ios::trunc);
        if (nodeTrackingFile.is_open()) {
            // Write header
            nodeTrackingFile << "TIMESTAMP | NODE | OPERATION | BLOCK | DETAILS\n";
            nodeTrackingFile << "--------------------------------------------------------------------------\n";
            nodeTrackingInitialized = true;
        } else {
            std::cerr << "WARNING: Failed to open node_tracking.log for writing" << std::endl;
        }
    }
}

// Flush all logs
inline void flushNodeTracking() {
    std::lock_guard<std::mutex> lock(nodeTrackingMutex);
    
    if (nodeTrackingFile.is_open()) {
        nodeTrackingFile.flush();
    }
}

// New function to dump complete sequence and mapping data
inline void dumpNodeSequenceAndMappings(const std::string& nodeId, const std::string& sequence, 
                                      const std::vector<std::pair<int32_t, int64_t>>& blockMappings) {
    // Only track node_1 and node_2
    if (nodeId != "node_1" && nodeId != "node_2") {
        return;
    }
    
    std::lock_guard<std::mutex> lock(nodeTrackingMutex);
    
    // Initialize if not already done
    if (!nodeTrackingInitialized) {
        initNodeTracking();
    }
    
    if (nodeTrackingFile.is_open()) {
        // Get current timestamp
        auto now = std::chrono::system_clock::now();
        auto nowMs = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
        auto epoch = nowMs.time_since_epoch();
        auto value = std::chrono::duration_cast<std::chrono::milliseconds>(epoch).count();
        
        // Write sequence information
        nodeTrackingFile << value << " | " 
                       << nodeId << " | " 
                       << "SEQUENCE_DUMP" << " | "
                       << "-1" << " | "
                       << "Full sequence length=" << sequence.length() << std::endl;
        
        // Write the full sequence (in chunks if needed)
        const size_t chunkSize = 80;
        for (size_t i = 0; i < sequence.length(); i += chunkSize) {
            size_t len = std::min(chunkSize, sequence.length() - i);
            nodeTrackingFile << value << " | " 
                           << nodeId << " | " 
                           << "SEQUENCE_CHUNK" << " | "
                           << i << " | "
                           << sequence.substr(i, len) << std::endl;
        }
        
        // Write mapping information
        nodeTrackingFile << value << " | " 
                       << nodeId << " | " 
                       << "MAPPING_DUMP" << " | "
                       << "-1" << " | "
                       << "Block mappings count=" << blockMappings.size() << std::endl;
        
        // Consolidate mappings into fewer lines
        const size_t mappingsPerLine = 10;
        for (size_t i = 0; i < blockMappings.size(); i += mappingsPerLine) {
            std::stringstream ss;
            ss << "Mappings: ";
            
            for (size_t j = i; j < std::min(i + mappingsPerLine, blockMappings.size()); j++) {
                if (j > i) ss << ", ";
                ss << "Block " << blockMappings[j].first << " → " << blockMappings[j].second;
            }
            
            nodeTrackingFile << value << " | " 
                           << nodeId << " | " 
                           << "MAPPING_CHUNK" << " | "
                           << i << " | "
                           << ss.str() << std::endl;
        }
        
        // Flush to ensure data is written
        nodeTrackingFile.flush();
    }
}

// Get current logging level
inline LogLevel getLoggingLevel() {
  return loggingLevel;
}

// Simple timer with concise output
class Timer {
public:
  // Update default level to INFO
  explicit Timer(const std::string &name, LogLevel level = LogLevel::INFO) 
      : name_(name), level_(level), start_(std::chrono::high_resolution_clock::now()) {}

  ~Timer() {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start_);
    // Use spdlog::info directly, level check is handled by spdlog
    spdlog::info(FMT_STRING("{}: {} ms"), name_, duration.count());
  }

private:
  std::string name_;
  LogLevel level_; // Keep level if needed for other logic, but not for spdlog filtering
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

// Define a better macro that doesn't conflict with std::chrono
#define TIME_OPERATION(name) \
  logging::Timer _timer_##__LINE__ { name }

// Debug log messages
template <typename... Args>
inline void debug(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::debug(fmt, std::forward<Args>(args)...);
}

inline void debug(const std::string& msg) {
  spdlog::debug(msg);
}

// Info log messages
template <typename... Args>
inline void info(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::info(fmt, std::forward<Args>(args)...);
}

inline void info(const std::string& msg) {
  spdlog::info(msg);
}

// Warning log messages
template <typename... Args>
inline void warn(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::warn(fmt, std::forward<Args>(args)...);
}

inline void warn(const std::string& msg) {
  spdlog::warn(msg);
}

// Error log messages
template <typename... Args>
inline void err(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::error(fmt, std::forward<Args>(args)...);
}

inline void err(const std::string& msg) {
  spdlog::error(msg);
}

// Critical log messages
template <typename... Args>
inline void critical(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::critical(fmt, std::forward<Args>(args)...);
}

inline void critical(const std::string& msg) {
  spdlog::critical(msg);
}

// Regular messages (normal level)
template <typename... Args>
inline void msg(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::info(fmt, std::forward<Args>(args)...);
}

inline void msg(const std::string& message) {
  spdlog::info(message);
}

// Verbose messages
template <typename... Args>
inline void verbose(fmt::format_string<Args...> fmt, Args&&... args) {
  spdlog::trace(fmt, std::forward<Args>(args)...);
}

inline void verbose(const std::string& message) {
  spdlog::trace(message);
}

} // namespace logging

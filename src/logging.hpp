#pragma once

// Try to include spdlog if available
#ifdef __has_include
  #if __has_include(<spdlog/spdlog.h>)
    #include <spdlog/spdlog.h>
    #include <spdlog/fmt/bundled/format.h>
    #define SPDLOG_AVAILABLE 1
  #else
    #define SPDLOG_AVAILABLE 0
  #endif
#else
  // Older compilers without __has_include
  #define SPDLOG_AVAILABLE 0
#endif

#include <chrono>
#include <iostream>
#include <string>
// utility header removed - not directly used
#include <sstream>
#include <fstream>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace logging {

// Logging verbosity levels
enum class LogLevel {
  TRACE,
  DEBUG,     // Debug level explicitly called out
  INFO,
  WARN,
  ERROR,
  CRITICAL,
  OFF,
  QUIET,     // Show only critical messages
  NORMAL,    // Default level, show important logs but not per-node details
  VERBOSE    // Show all logs, including detailed per-node processing
};

// Global verbosity level
static LogLevel loggingLevel = LogLevel::NORMAL;

// Node tracking file for node_1 and node_2
static std::ofstream nodeTrackingFile;
static std::mutex nodeTrackingMutex;
static bool nodeTrackingInitialized = false;

// Buffers for consolidating related log entries - removed in favor of immediate logging
// static std::unordered_map<std::string, std::string> consolidatedLogs;

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

// Forward declaration of flushNodeTracking
inline void flushNodeTracking();


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

// Function to set logging verbosity
inline void setLoggingLevel(LogLevel level) {
  loggingLevel = level;
}

// Get current logging level
inline LogLevel getLoggingLevel() {
  return loggingLevel;
}

// Helper to check if a certain type of message should be logged
inline bool shouldLog(LogLevel messageLevel) {
  return static_cast<int>(messageLevel) >= static_cast<int>(loggingLevel);
}

// Simple timer with concise output
class Timer {
public:
  explicit Timer(const std::string &name, LogLevel level = LogLevel::NORMAL)
      : name_(name), level_(level), start_(std::chrono::high_resolution_clock::now()) {}

  ~Timer() {
    if (!shouldLog(level_)) return;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start_);
    
    #if SPDLOG_AVAILABLE
    spdlog::info("{}: {} ms", name_, duration.count());
    #else
    std::cout << name_ << ": " << duration.count() << " ms" << std::endl;
    #endif
  }

private:
  std::string name_;
  LogLevel level_;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

// Define a better macro that doesn't conflict with std::chrono
#define TIME_OPERATION(name) \
  logging::Timer _timer_##__LINE__ { name }

// Logging helper functions that wrap spdlog or fall back to cout

// Regular info message (NORMAL level by default)
template <typename FormatString, typename... Args>
inline void msg(const FormatString& fmt, Args&&... args) {
  if (!shouldLog(LogLevel::NORMAL)) return;
  
  #if SPDLOG_AVAILABLE
  spdlog::info(fmt::runtime(fmt), std::forward<Args>(args)...);
  #else
  // Simple fallback implementation
  std::cout << "INFO: ";
  std::cout << fmt;
  ((std::cout << " " << args), ...);
  std::cout << std::endl;
  std::cout.flush();
  #endif
}

// Verbose logging for detailed operations (VERBOSE level)
template <typename FormatString, typename... Args>
inline void verbose(const FormatString& fmt, Args&&... args) {
  if (!shouldLog(LogLevel::VERBOSE)) return;
  
  #if SPDLOG_AVAILABLE
  spdlog::info(fmt::runtime(fmt), std::forward<Args>(args)...);
  #else
  // Simple fallback implementation
  std::cout << "VERBOSE: " << fmt << std::endl;
  std::cout.flush();
  #endif
}

// Critical info message (always shown)
template <typename FormatString, typename... Args>
inline void critical(const FormatString& fmt, Args&&... args) {
  // Critical messages show in all modes
  #if SPDLOG_AVAILABLE
  spdlog::info(fmt::runtime(fmt), std::forward<Args>(args)...);
  #else
  // Simple fallback implementation
  std::cout << "CRITICAL: ";
  std::cout << fmt;
  ((std::cout << " " << args), ...);
  std::cout << std::endl;
  std::cout.flush();
  #endif
}

// Warning message (shown in NORMAL and VERBOSE)
template <typename FormatString, typename... Args>
inline void warn(const FormatString& fmt, Args&&... args) {
  if (!shouldLog(LogLevel::NORMAL)) return;
  
  #if SPDLOG_AVAILABLE
  spdlog::warn(fmt::runtime(fmt), std::forward<Args>(args)...);
  #else
  // Simple fallback implementation
  std::cerr << "WARN: " << fmt << std::endl;
  #endif
}

// Error message (always shown)
template <typename FormatString, typename... Args>
inline void err(const FormatString& fmt, Args&&... args) {
  // Error messages show in all modes
  #if SPDLOG_AVAILABLE
  spdlog::error(fmt::runtime(fmt), std::forward<Args>(args)...);
  #else
  // Simple fallback implementation
  std::cerr << "ERROR: ";
  std::cerr << fmt;
  ((std::cerr << " " << args), ...);
  std::cerr << std::endl;
  std::cerr.flush();
  #endif
}

// Debug function (VERBOSE level only)
template <typename FormatString, typename... Args>
inline void debug(const FormatString& fmt, Args&&... args) {
  if (!shouldLog(LogLevel::VERBOSE)) return;
  
  #if SPDLOG_AVAILABLE
  spdlog::debug(fmt::runtime(fmt), std::forward<Args>(args)...);
  #else
  // Simple fallback implementation
  // std::cout << "DEBUG: " << fmt << std::endl;
  #endif
}

// Forward declarations for all logging functions
template <typename... Args>
void trace(const std::string& fmt, Args&&... args);

template <typename... Args>
void info(const std::string& fmt, Args&&... args);

// Implementation for all logging functions
template <typename... Args>
void trace(const std::string& fmt, Args&&... args) {
    // Implementation will use spdlog in the cpp file
}

template <typename... Args>
void info(const std::string& fmt, Args&&... args) {
    // Implementation will use spdlog in the cpp file
}

// Simple logging functions that don't depend on fmt
inline void debug(const std::string& message) {
    std::cout << "[DEBUG] " << message << std::endl;
    std::cout.flush();
}

inline void info(const std::string& message) {
    std::cout << "[INFO] " << message << std::endl;
    std::cout.flush();
}

inline void warn(const std::string& message) {
    std::cout << "[WARN] " << message << std::endl;
    std::cout.flush();
}

inline void err(const std::string& message) {
    std::cerr << "[ERROR] " << message << std::endl;
    std::cerr.flush();
}

inline void verbose(const std::string& message) {
    std::cout << "[VERBOSE] " << message << std::endl;
    std::cout.flush();
}

// Template versions for formatting with one argument
template <typename T>
void debug(const std::string& fmt, const T& arg) {
    std::ostringstream oss;
    size_t pos = fmt.find("{}");
    if (pos != std::string::npos) {
        oss << fmt.substr(0, pos) << arg << fmt.substr(pos + 2);
    } else {
        oss << fmt;
    }
    debug(oss.str());
}

template <typename T>
void info(const std::string& fmt, const T& arg) {
    std::ostringstream oss;
    size_t pos = fmt.find("{}");
    if (pos != std::string::npos) {
        oss << fmt.substr(0, pos) << arg << fmt.substr(pos + 2);
    } else {
        oss << fmt;
    }
    info(oss.str());
}

template <typename T>
void warn(const std::string& fmt, const T& arg) {
    std::ostringstream oss;
    size_t pos = fmt.find("{}");
    if (pos != std::string::npos) {
        oss << fmt.substr(0, pos) << arg << fmt.substr(pos + 2);
    } else {
        oss << fmt;
    }
    warn(oss.str());
}

template <typename T>
void err(const std::string& fmt, const T& arg) {
    std::ostringstream oss;
    size_t pos = fmt.find("{}");
    if (pos != std::string::npos) {
        oss << fmt.substr(0, pos) << arg << fmt.substr(pos + 2);
    } else {
        oss << fmt;
    }
    err(oss.str());
}

template <typename T>
void verbose(const std::string& fmt, const T& arg) {
    std::ostringstream oss;
    size_t pos = fmt.find("{}");
    if (pos != std::string::npos) {
        oss << fmt.substr(0, pos) << arg << fmt.substr(pos + 2);
    } else {
        oss << fmt;
    }
    verbose(oss.str());
}

// Template version for two arguments
template <typename T1, typename T2>
void debug(const std::string& fmt, const T1& arg1, const T2& arg2) {
    std::ostringstream oss;
    size_t pos1 = fmt.find("{}");
    if (pos1 != std::string::npos) {
        size_t pos2 = fmt.find("{}", pos1 + 2);
        if (pos2 != std::string::npos) {
            oss << fmt.substr(0, pos1) << arg1 << fmt.substr(pos1 + 2, pos2 - pos1 - 2) << arg2 << fmt.substr(pos2 + 2);
        } else {
            oss << fmt.substr(0, pos1) << arg1 << fmt.substr(pos1 + 2);
        }
    } else {
        oss << fmt;
    }
    debug(oss.str());
}

template <typename T1, typename T2>
void info(const std::string& fmt, const T1& arg1, const T2& arg2) {
    std::ostringstream oss;
    size_t pos1 = fmt.find("{}");
    if (pos1 != std::string::npos) {
        size_t pos2 = fmt.find("{}", pos1 + 2);
        if (pos2 != std::string::npos) {
            oss << fmt.substr(0, pos1) << arg1 << fmt.substr(pos1 + 2, pos2 - pos1 - 2) << arg2 << fmt.substr(pos2 + 2);
        } else {
            oss << fmt.substr(0, pos1) << arg1 << fmt.substr(pos1 + 2);
        }
    } else {
        oss << fmt;
    }
    info(oss.str());
}

// Template version for three arguments
template <typename T1, typename T2, typename T3>
void info(const std::string& fmt, const T1& arg1, const T2& arg2, const T3& arg3) {
    std::ostringstream oss;
    size_t pos1 = fmt.find("{}");
    if (pos1 != std::string::npos) {
        size_t pos2 = fmt.find("{}", pos1 + 2);
        if (pos2 != std::string::npos) {
            size_t pos3 = fmt.find("{}", pos2 + 2);
            if (pos3 != std::string::npos) {
                oss << fmt.substr(0, pos1) << arg1
                    << fmt.substr(pos1 + 2, pos2 - pos1 - 2) << arg2
                    << fmt.substr(pos2 + 2, pos3 - pos2 - 2) << arg3
                    << fmt.substr(pos3 + 2);
            } else {
                oss << fmt.substr(0, pos1) << arg1
                    << fmt.substr(pos1 + 2, pos2 - pos1 - 2) << arg2
                    << fmt.substr(pos2 + 2);
            }
        } else {
            oss << fmt.substr(0, pos1) << arg1 << fmt.substr(pos1 + 2);
        }
    } else {
        oss << fmt;
    }
    info(oss.str());
}

} // namespace logging
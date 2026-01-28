#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <chrono>

#ifdef _WIN32
#include <io.h>
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>
#endif

namespace output {

// ============================================================================
// Global Configuration  
// ============================================================================

struct Config {
    std::atomic<bool> quiet{false};     // Minimal output (errors only)
    std::atomic<bool> verbose{false};   // Extra debug output
    std::atomic<bool> plain{false};     // No colors or unicode
    std::atomic<bool> isTTY{true};      // Auto-detected
    
    void detectTTY() {
        isTTY = isatty(fileno(stderr)) && isatty(fileno(stdout));
    }
};

inline Config& config() {
    static Config cfg;
    return cfg;
}

// Initialize output system. Call once at startup.
inline void init(bool quiet = false, bool verbose = false, bool plain = false) {
    config().quiet = quiet;
    config().verbose = verbose;
    config().plain = plain;
    config().detectTTY();
    
    // Auto-enable plain mode if not a TTY
    if (!config().isTTY) {
        config().plain = true;
    }
    
    // Configure spdlog
    spdlog::set_pattern("%v");
    if (quiet) {
        spdlog::set_level(spdlog::level::err);
    } else if (verbose) {
        spdlog::set_level(spdlog::level::debug);
    } else {
        spdlog::set_level(spdlog::level::info);
    }
}

// ============================================================================
// Color/Style Codes
// ============================================================================

namespace style {
    inline const char* reset()   { return config().plain ? "" : "\033[0m"; }
    inline const char* bold()    { return config().plain ? "" : "\033[1m"; }
    inline const char* dim()     { return config().plain ? "" : "\033[2m"; }
    inline const char* red()     { return config().plain ? "" : "\033[31m"; }
    inline const char* green()   { return config().plain ? "" : "\033[32m"; }
    inline const char* yellow()  { return config().plain ? "" : "\033[33m"; }
    inline const char* blue()    { return config().plain ? "" : "\033[34m"; }
    inline const char* cyan()    { return config().plain ? "" : "\033[36m"; }
    inline const char* magenta() { return config().plain ? "" : "\033[35m"; }
}

// ============================================================================
// Box Drawing Characters (Unicode or ASCII fallback)
// ============================================================================

namespace box {
    inline const char* topLeft()     { return config().plain ? "+" : "┌"; }
    inline const char* topRight()    { return config().plain ? "+" : "┐"; }
    inline const char* bottomLeft()  { return config().plain ? "+" : "└"; }
    inline const char* bottomRight() { return config().plain ? "+" : "┘"; }
    inline const char* horizontal()  { return config().plain ? "-" : "─"; }
    inline const char* vertical()    { return config().plain ? "|" : "│"; }
    inline const char* arrow()       { return config().plain ? "->" : "→"; }
    inline const char* bullet()      { return config().plain ? "*" : "•"; }
    inline const char* check()       { return config().plain ? "[OK]" : "✓"; }
    inline const char* cross()       { return config().plain ? "[X]" : "✗"; }
}

// ============================================================================
// Status Indicators
// ============================================================================

inline std::string status_ok() {
    return fmt::format("{}{}{}", style::green(), box::check(), style::reset());
}

inline std::string status_fail() {
    return fmt::format("{}{}{}", style::red(), box::cross(), style::reset());
}

// ============================================================================
// Core Output Functions
// ============================================================================

// Error messages - always shown
template <typename... Args>
inline void error(fmt::format_string<Args...> fmt, Args&&... args) {
    std::cerr << style::red() << "error: " << style::reset();
    std::cerr << fmt::format(fmt, std::forward<Args>(args)...) << "\n";
}

inline void error(const std::string& msg) {
    std::cerr << style::red() << "error: " << style::reset() << msg << "\n";
}

// Warning messages - shown unless quiet
template <typename... Args>
inline void warn(fmt::format_string<Args...> fmt, Args&&... args) {
    if (config().quiet) return;
    std::cerr << style::yellow() << "warn: " << style::reset();
    std::cerr << fmt::format(fmt, std::forward<Args>(args)...) << "\n";
}

inline void warn(const std::string& msg) {
    if (config().quiet) return;
    std::cerr << style::yellow() << "warn: " << style::reset() << msg << "\n";
}

// Info messages - shown by default, hidden in quiet mode
template <typename... Args>
inline void info(fmt::format_string<Args...> fmt, Args&&... args) {
    if (config().quiet) return;
    std::cerr << fmt::format(fmt, std::forward<Args>(args)...) << "\n";
}

inline void info(const std::string& msg) {
    if (config().quiet) return;
    std::cerr << msg << "\n";
}

// Debug messages - shown only in verbose mode
template <typename... Args>
inline void debug(fmt::format_string<Args...> fmt, Args&&... args) {
    if (!config().verbose) return;
    std::cerr << style::dim() << fmt::format(fmt, std::forward<Args>(args)...) 
              << style::reset() << "\n";
}

inline void debug(const std::string& msg) {
    if (!config().verbose) return;
    std::cerr << style::dim() << msg << style::reset() << "\n";
}

// ============================================================================
// Structured Output
// ============================================================================

// Stage header - marks major pipeline steps
inline void stage(const std::string& name) {
    if (config().quiet) return;
    std::cerr << style::bold() << style::cyan() << box::arrow() << " " 
              << name << style::reset() << "\n";
}

// Step within a stage
template <typename... Args>
inline void step(fmt::format_string<Args...> fmt, Args&&... args) {
    if (config().quiet) return;
    std::cerr << "  " << fmt::format(fmt, std::forward<Args>(args)...) << "\n";
}

inline void step(const std::string& msg) {
    if (config().quiet) return;
    std::cerr << "  " << msg << "\n";
}

// Done message with optional timing
inline void done(const std::string& what, int64_t ms) {
    if (config().quiet) return;
    std::cerr << "  " << status_ok() << " " << what;
    if (ms >= 0) {
        // Format time nicely: ms for < 1s, s for < 60s, m:ss for >= 60s
        if (ms < 1000) {
            std::cerr << style::dim() << " (" << ms << "ms)" << style::reset();
        } else if (ms < 60000) {
            std::cerr << style::dim() << " (" << std::fixed << std::setprecision(1) 
                      << (ms / 1000.0) << "s)" << style::reset();
        } else {
            int64_t secs = ms / 1000;
            int64_t mins = secs / 60;
            secs = secs % 60;
            std::cerr << style::dim() << " (" << mins << "m" << std::setw(2) 
                      << std::setfill('0') << secs << "s)" << style::reset();
        }
    }
    std::cerr << "\n";
}

inline void done(const std::string& what) {
    done(what, -1);
}

// Progress indicator (overwrites line on TTY)
inline void progress(const std::string& msg) {
    if (config().quiet) return;
    if (config().isTTY && !config().plain) {
        std::cerr << "\r\033[K  " << style::dim() << msg << style::reset() << std::flush;
    }
}

// Progress with percentage and elapsed time
inline void progress_pct(const std::string& task, size_t current, size_t total, int64_t elapsed_ms = -1) {
    if (config().quiet) return;
    if (config().isTTY && !config().plain) {
        int pct = total > 0 ? static_cast<int>((current * 100) / total) : 0;
        std::cerr << "\r\033[K  " << style::dim() << task << " " << pct << "%";
        if (elapsed_ms >= 0) {
            if (elapsed_ms < 1000) {
                std::cerr << " (" << elapsed_ms << "ms)";
            } else {
                std::cerr << " (" << std::fixed << std::setprecision(1) 
                          << (elapsed_ms / 1000.0) << "s)";
            }
        }
        std::cerr << style::reset() << std::flush;
    }
}

inline void progress_clear() {
    if (config().isTTY && !config().plain) {
        std::cerr << "\r\033[K" << std::flush;
    }
}

// ============================================================================
// Summary Box
// ============================================================================

inline void print_header(const std::string& title, const std::string& version = "") {
    if (config().quiet) return;
    
    std::cerr << style::bold() << style::cyan() << box::topLeft() 
              << box::horizontal() << " " << title;
    if (!version.empty()) {
        std::cerr << style::reset() << style::dim() << " v" << version;
    }
    std::cerr << style::reset() << style::bold() << style::cyan() << " ";
    for (int i = 0; i < 50; i++) std::cerr << box::horizontal();
    std::cerr << box::topRight() << style::reset() << "\n";
}

inline void print_row(const std::string& label, const std::string& value) {
    if (config().quiet) return;
    std::cerr << style::bold() << style::cyan() << box::vertical() << style::reset()
              << " " << style::bold() << label << ": " << style::reset() << value << "\n";
}

inline void print_footer() {
    if (config().quiet) return;
    std::cerr << style::bold() << style::cyan() << box::bottomLeft();
    for (int i = 0; i < 68; i++) std::cerr << box::horizontal();
    std::cerr << box::bottomRight() << style::reset() << "\n\n";
}

// ============================================================================
// Timer Utility
// ============================================================================

class Timer {
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}
    
    int64_t elapsed_ms() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count();
    }
    
    void reset() { start_ = std::chrono::high_resolution_clock::now(); }

private:
    std::chrono::high_resolution_clock::time_point start_;
};

} // namespace output

// ============================================================================
// Signal Handling
// ============================================================================

#include <csignal>

namespace signals {

// Global flag for graceful interruption
inline std::atomic<bool>& interrupted() {
    static std::atomic<bool> flag{false};
    return flag;
}

// Check if interrupted and optionally print message
inline bool check_interrupted(bool print_message = true) {
    if (interrupted()) {
        if (print_message) {
            output::warn("Interrupted by user");
        }
        return true;
    }
    return false;
}

// Signal handler
inline void handler(int signum) {
    (void)signum;  // unused
    interrupted() = true;
    // Print newline to clean up any progress output
    std::cerr << "\n";
}

// Install signal handlers for graceful interruption
inline void install_handlers() {
    std::signal(SIGINT, handler);
    std::signal(SIGTERM, handler);
#ifndef _WIN32
    std::signal(SIGHUP, handler);
#endif
}

} // namespace signals

// ============================================================================
// Legacy Logging Namespace (for compatibility)
// ============================================================================

namespace logging {

template <typename... Args>
inline void debug(fmt::format_string<Args...> fmt, Args&&... args) {
    output::debug(fmt, std::forward<Args>(args)...);
}
inline void debug(const std::string& msg) { output::debug(msg); }

template <typename... Args>
inline void info(fmt::format_string<Args...> fmt, Args&&... args) {
    output::info(fmt, std::forward<Args>(args)...);
}
inline void info(const std::string& msg) { output::info(msg); }

template <typename... Args>
inline void warn(fmt::format_string<Args...> fmt, Args&&... args) {
    output::warn(fmt, std::forward<Args>(args)...);
}
inline void warn(const std::string& msg) { output::warn(msg); }

template <typename... Args>
inline void err(fmt::format_string<Args...> fmt, Args&&... args) {
    output::error(fmt, std::forward<Args>(args)...);
}
inline void err(const std::string& msg) { output::error(msg); }

template <typename... Args>
inline void msg(fmt::format_string<Args...> fmt, Args&&... args) {
    output::info(fmt, std::forward<Args>(args)...);
}
inline void msg(const std::string& message) { output::info(message); }

} // namespace logging

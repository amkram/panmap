#pragma once

// Force-enable spdlog since we've added it as a dependency
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#define SPDLOG_AVAILABLE 1

#include <string>

namespace logging {

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

} // namespace logging

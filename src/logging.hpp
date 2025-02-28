#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <memory>
#include <spdlog/sinks/daily_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>

static const bool DEBUG = true;

// Use the global namespace explicitly for spdlog
#define msg(...) spdlog::info(__VA_ARGS__)
#define debug(...)                                                             \
  if (DEBUG)                                                                   \
  spdlog::debug(__VA_ARGS__)
#define warn(...) spdlog::warn(__VA_ARGS__)
#define err(...) spdlog::error(__VA_ARGS__)

namespace logging {

// Simple timer with concise output
class Timer {
public:
  explicit Timer(const std::string &name)
      : name_(name), start_(std::chrono::high_resolution_clock::now()) {}

  ~Timer() {
    auto end = std::chrono::high_resolution_clock::now();
    auto ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start_)
            .count();
    msg("{}: {}ms", name_, ms);
  }

private:
  std::string name_;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

#define time(name)                                                             \
  logging::Timer _t { name }

} // namespace logging

#endif // LOGGING_HPP
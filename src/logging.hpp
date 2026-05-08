#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <mutex>

#ifdef _WIN32
#include <io.h>
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>
#include <sys/ioctl.h>
#endif

namespace output {

struct Config {
    std::atomic<bool> quiet{false};
    std::atomic<bool> verbose{false};
    std::atomic<bool> plain{false};
    std::atomic<bool> isTTY{true};

    void detectTTY() { isTTY = isatty(fileno(stderr)); }
};

inline Config& config() {
    static Config cfg;
    return cfg;
}

inline int term_width() {
#ifndef _WIN32
    struct winsize w;
    if (ioctl(STDERR_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0) return w.ws_col;
#endif
    return 80;
}

inline void init(bool quiet = false, bool verbose = false, bool plain = false) {
    config().quiet = quiet;
    config().verbose = verbose;
    config().plain = plain;
    config().detectTTY();

    if (!config().isTTY) config().plain = true;

    spdlog::set_pattern("%v");
    if (quiet)         spdlog::set_level(spdlog::level::err);
    else if (verbose)  spdlog::set_level(spdlog::level::debug);
    else               spdlog::set_level(spdlog::level::info);
}

namespace style {
inline const char* reset()   { return config().plain ? "" : "\033[0m";  }
inline const char* bold()    { return config().plain ? "" : "\033[1m";  }
inline const char* dim()     { return config().plain ? "" : "\033[2m";  }
inline const char* red()     { return config().plain ? "" : "\033[31m"; }
inline const char* green()   { return config().plain ? "" : "\033[32m"; }
inline const char* yellow()  { return config().plain ? "" : "\033[33m"; }
inline const char* blue()    { return config().plain ? "" : "\033[34m"; }
inline const char* magenta() { return config().plain ? "" : "\033[35m"; }
inline const char* cyan()    { return config().plain ? "" : "\033[36m"; }
}  // namespace style

namespace box {
inline const char* arrow() { return config().plain ? "->"   : "→"; }
inline const char* check() { return config().plain ? "[ok]" : "✓"; }
inline const char* cross() { return config().plain ? "[x]"  : "✗"; }
inline const char* dot()   { return config().plain ? "*"    : "·"; }
inline const char* bullet(){ return config().plain ? "*"    : "•"; }
}  // namespace box

// Spinner frames (Braille pattern, 10 frames at ~80ms = ~12fps)
inline const char* spinner_frame(int i) {
    static const char* frames[] = {"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"};
    static const char* plain_frames[] = {"|", "/", "-", "\\", "|", "/", "-", "\\", "|", "/"};
    return config().plain ? plain_frames[i % 10] : frames[i % 10];
}

inline std::string format_duration(int64_t ms) {
    if (ms < 0) return "";
    if (ms < 1000) return fmt::format("{}ms", ms);
    if (ms < 60000) return fmt::format("{:.1f}s", ms / 1000.0);
    int64_t s = ms / 1000;
    return fmt::format("{}m{:02d}s", s / 60, static_cast<int>(s % 60));
}

template <typename... Args>
inline void error(fmt::format_string<Args...> fmt, Args&&... args) {
    std::cerr << style::red() << "error" << style::reset() << style::dim() << ": " << style::reset()
              << fmt::format(fmt, std::forward<Args>(args)...) << "\n";
}
inline void error(const std::string& msg) {
    std::cerr << style::red() << "error" << style::reset() << style::dim() << ": " << style::reset() << msg << "\n";
}

template <typename... Args>
inline void warn(fmt::format_string<Args...> fmt, Args&&... args) {
    if (config().quiet) return;
    std::cerr << style::yellow() << "warn" << style::reset() << style::dim() << ": " << style::reset()
              << fmt::format(fmt, std::forward<Args>(args)...) << "\n";
}
inline void warn(const std::string& msg) {
    if (config().quiet) return;
    std::cerr << style::yellow() << "warn" << style::reset() << style::dim() << ": " << style::reset() << msg << "\n";
}

// Verbose-only chatter
template <typename... Args>
inline void info(fmt::format_string<Args...> fmt, Args&&... args) {
    if (!config().verbose) return;
    std::cerr << style::dim() << fmt::format(fmt, std::forward<Args>(args)...) << style::reset() << "\n";
}
inline void info(const std::string& msg) {
    if (!config().verbose) return;
    std::cerr << style::dim() << msg << style::reset() << "\n";
}

template <typename... Args>
inline void debug(fmt::format_string<Args...> fmt, Args&&... args) {
    if (!config().verbose) return;
    std::cerr << style::dim() << fmt::format(fmt, std::forward<Args>(args)...) << style::reset() << "\n";
}
inline void debug(const std::string& msg) {
    if (!config().verbose) return;
    std::cerr << style::dim() << msg << style::reset() << "\n";
}

// One-line tool+version banner. No subtitle in default mode.
inline void banner(const std::string& version, const std::string& /*subtitle*/ = "") {
    if (config().quiet) return;
    std::cerr << style::bold() << "panmap" << style::reset()
              << style::dim() << " " << version << style::reset() << "\n\n";
}

// Stage line columns:
//   "  ICN  LABEL.  SUBJECT.................  STAT.........    TIME"
//
// Subject and stat are left-padded to a min width so columns align across
// stages; time is right-padded so right edges line up. Values longer than
// the min just push the next column over (no truncation). All input strings
// must be plain text (no ANSI codes), or width math will be wrong — embed
// styling outside via the helpers below.
constexpr int kLabelWidth   = 6;
constexpr int kSubjectWidth = 22;
constexpr int kStatWidth    = 14;
constexpr int kTimeWidth    = 6;

inline std::string status_ok()   { return fmt::format("{}{}{}", style::green(), box::check(), style::reset()); }
inline std::string status_fail() { return fmt::format("{}{}{}", style::red(),   box::cross(), style::reset()); }

inline std::string pad_right(const std::string& s, int w) {
    if ((int)s.size() >= w) return s;
    return s + std::string(w - s.size(), ' ');
}
inline std::string pad_left(const std::string& s, int w) {
    if ((int)s.size() >= w) return s;
    return std::string(w - s.size(), ' ') + s;
}

// Format an integer with thousands separators (commas) for terminal display.
inline std::string fmt_count(uint64_t n) {
    std::string s = std::to_string(n);
    int insert = (int)s.size() - 3;
    while (insert > 0) {
        s.insert(insert, ",");
        insert -= 3;
    }
    return s;
}

inline void write_status_line(const std::string& icon,
                              const char* icon_color,
                              const std::string& label,
                              const std::string& subject,
                              const std::string& stat,
                              int64_t ms) {
    if (config().quiet) return;
    std::cerr << "  " << icon_color << icon << style::reset() << "  "
              << style::bold() << pad_right(label, kLabelWidth) << style::reset() << "  "
              << pad_right(subject, kSubjectWidth) << "  "
              << style::dim() << pad_right(stat, kStatWidth) << style::reset() << "  ";
    if (ms >= 0) {
        std::cerr << style::dim() << pad_left(format_duration(ms), kTimeWidth) << style::reset();
    }
    std::cerr << "\n";
}

inline void done(const std::string& label,
                 const std::string& subject,
                 const std::string& stat = "",
                 int64_t ms = -1) {
    write_status_line(box::check(), style::green(), label, subject, stat, ms);
}
inline void fail(const std::string& label, const std::string& subject, const std::string& stat = "") {
    write_status_line(box::cross(), style::red(), label, subject, stat, -1);
}

// Legacy 1- and 2-arg done() are now debug-only (kept for old callsites)
inline void done(const std::string& what) {
    if (!config().verbose) return;
    std::cerr << style::dim() << "  " << status_ok() << " " << what << style::reset() << "\n";
}
inline void done(const std::string& what, int64_t ms) {
    if (!config().verbose) return;
    std::cerr << style::dim() << "  " << status_ok() << " " << what
              << " (" << format_duration(ms) << ")" << style::reset() << "\n";
}

// Trailing blank line so the prompt has breathing room.
inline void summary(int64_t /*total_ms*/) {
    if (config().quiet) return;
    std::cerr << "\n";
}

// --- Legacy stage/step/done helpers (now debug-only) ---
inline void stage(const std::string& name) {
    if (!config().verbose) return;
    std::cerr << style::dim() << box::arrow() << " " << name << style::reset() << "\n";
}
template <typename... Args>
inline void step(fmt::format_string<Args...> fmt, Args&&... args) {
    if (!config().verbose) return;
    std::cerr << style::dim() << "  " << fmt::format(fmt, std::forward<Args>(args)...) << style::reset() << "\n";
}
inline void step(const std::string& msg) {
    if (!config().verbose) return;
    std::cerr << style::dim() << "  " << msg << style::reset() << "\n";
}
// --- Indeterminate progress message (one-line, overwritten) ---
inline void progress(const std::string& msg) {
    if (config().quiet || !config().isTTY || config().plain) return;
    std::cerr << "\r\033[K" << style::dim() << "  " << msg << style::reset() << std::flush;
}
inline void progress_clear() {
    if (config().isTTY && !config().plain) std::cerr << "\r\033[K" << std::flush;
}
inline void progress_pct(const std::string& task, size_t current, size_t total, int64_t elapsed_ms = -1) {
    if (config().quiet || !config().isTTY || config().plain) return;
    int pct = total > 0 ? static_cast<int>((current * 100) / total) : 0;
    std::cerr << "\r\033[K" << style::dim() << "  " << task << " " << pct << "%";
    if (elapsed_ms >= 0) std::cerr << " (" << format_duration(elapsed_ms) << ")";
    std::cerr << style::reset() << std::flush;
}

// === Progress bar ===
//
// In-place line:
//     <spinner> <label>  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━╸    42%  · 1.2s
//
// On completion, clear() leaves the line empty so caller can emit a done() line.
class ProgressBar {
   public:
    ProgressBar(std::string label, uint64_t total)
        : label_(std::move(label)), total_(total) {
        start_ = std::chrono::steady_clock::now();
        last_render_ = start_ - std::chrono::seconds(1);
    }

    ~ProgressBar() { clear(); }

    void set(uint64_t current) {
        current_.store(current, std::memory_order_relaxed);
        maybe_render();
    }
    void add(uint64_t delta = 1) {
        current_.fetch_add(delta, std::memory_order_relaxed);
        maybe_render();
    }
    void set_label(const std::string& l) {
        std::lock_guard<std::mutex> g(mu_);
        label_ = l;
    }

    void clear() {
        if (!config().isTTY || config().plain || config().quiet) return;
        std::lock_guard<std::mutex> g(mu_);
        std::cerr << "\r\033[K" << std::flush;
    }

    int64_t elapsed_ms() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::steady_clock::now() - start_).count();
    }

   private:
    void maybe_render() {
        if (config().quiet || !config().isTTY || config().plain) return;
        auto now = std::chrono::steady_clock::now();
        std::unique_lock<std::mutex> lk(mu_, std::try_to_lock);
        if (!lk.owns_lock()) return;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_render_).count() < 80) return;
        last_render_ = now;
        render(now);
    }

    void render(std::chrono::steady_clock::time_point now) {
        uint64_t cur = current_.load(std::memory_order_relaxed);
        if (total_ > 0 && cur > total_) cur = total_;
        double frac = total_ > 0 ? static_cast<double>(cur) / total_ : 0.0;
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count();

        int frame = static_cast<int>(elapsed / 80);
        std::string spin  = spinner_frame(frame);
        std::string pct   = fmt::format("{:>3}%", static_cast<int>(frac * 100));
        std::string ela   = format_duration(elapsed);

        // Layout: "  <spin> <label-padded> <bar>  <pct>  · <ela>"
        std::string label_padded = label_;
        if ((int)label_padded.size() < kLabelWidth) label_padded.append(kLabelWidth - label_padded.size(), ' ');

        // Reserved visible chars (no ANSI): "  X " + label + " " + bar + "  " + "NNN%" + "  · " + ela
        int reserved = 2 + 1 + 1 + (int)label_padded.size() + 1 + 2 + 4 + 4 + (int)ela.size();
        int width = term_width();
        int bar_w = width - reserved - 2;
        if (bar_w > 50) bar_w = 50;
        if (bar_w < 6) bar_w = 6;

        int filled = static_cast<int>(frac * bar_w);
        if (filled > bar_w) filled = bar_w;
        std::string bar;
        bar.reserve(bar_w * 3 + 16);
        bar += style::cyan();
        for (int i = 0; i < filled; ++i) bar += "━";
        if (filled < bar_w) {
            bar += style::reset();
            bar += style::dim();
            if (filled < bar_w) bar += "╸";
            for (int i = filled + 1; i < bar_w; ++i) bar += "─";
        }
        bar += style::reset();

        std::cerr << "\r\033[K"
                  << "  " << style::cyan() << spin << style::reset() << "  "
                  << style::bold() << label_padded << style::reset() << "  "
                  << bar << "  "
                  << style::bold() << pct << style::reset() << "  "
                  << style::dim() << pad_left(ela, kTimeWidth) << style::reset()
                  << std::flush;
    }

    std::string label_;
    uint64_t total_;
    std::atomic<uint64_t> current_{0};
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point last_render_;
    std::mutex mu_;
};

}  // namespace output

#include <csignal>

namespace signals {

inline std::atomic<bool>& interrupted() {
    static std::atomic<bool> flag{false};
    return flag;
}

inline bool check_interrupted(bool print_message = true) {
    if (interrupted()) {
        if (print_message) output::warn("Interrupted by user");
        return true;
    }
    return false;
}

inline void handler(int signum) {
    (void)signum;
    interrupted() = true;
    std::cerr << "\n";
}

inline void install_handlers() {
    std::signal(SIGINT, handler);
    std::signal(SIGTERM, handler);
#ifndef _WIN32
    std::signal(SIGHUP, handler);
#endif
}

}  // namespace signals

// Legacy logging namespace (default-visible info now requires --verbose)
namespace logging {
template <typename... Args> inline void debug(fmt::format_string<Args...> f, Args&&... a) { output::debug(f, std::forward<Args>(a)...); }
template <typename... Args> inline void info (fmt::format_string<Args...> f, Args&&... a) { output::info (f, std::forward<Args>(a)...); }
template <typename... Args> inline void warn (fmt::format_string<Args...> f, Args&&... a) { output::warn (f, std::forward<Args>(a)...); }
template <typename... Args> inline void err  (fmt::format_string<Args...> f, Args&&... a) { output::error(f, std::forward<Args>(a)...); }
template <typename... Args> inline void msg  (fmt::format_string<Args...> f, Args&&... a) { output::info (f, std::forward<Args>(a)...); }

inline void debug(const std::string& s) { output::debug(s); }
inline void info (const std::string& s) { output::info (s); }
inline void warn (const std::string& s) { output::warn (s); }
inline void err  (const std::string& s) { output::error(s); }
inline void msg  (const std::string& s) { output::info (s); }
}  // namespace logging

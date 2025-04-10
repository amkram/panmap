#pragma once

// This is a dedicated precompiled header for the panmap project
// It contains headers that rarely change during development

// Standard C++ headers
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>
#include <algorithm>
#include <functional>
#include <atomic>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <exception>
#include <array>
#include <utility>
#include <thread>

// Standard C headers
#include <csignal>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>

// Boost headers
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/program_options.hpp>

// Third-party headers
#include <tbb/global_control.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h> 
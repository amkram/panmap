#pragma once

// Compatibility header for panman source compilation
// This header provides missing includes and using declarations

// Standard library includes needed by panman
#include <cmath>
#include <utility>
#include <filesystem>

// TBB headers that panman source files need but don't include
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

// Fix for panman using 'pair' instead of 'std::pair' in template declarations
using std::pair;

# Contributing to panmap

## Building from source

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Dependencies

- CMake 3.14+
- C++20 compiler (GCC 11+, Clang 14+, AppleClang 15+)
- Boost (program_options, iostreams, filesystem, date_time)
- Protocol Buffers
- Cap'n Proto
- htslib
- zlib, Eigen3

Other dependencies (TBB, jsoncpp, spdlog, abseil, zstd, libdeflate) are fetched automatically.

## Testing

Integration tests run against test data in `src/test/data/`:

```bash
ctest --test-dir build -R integration --output-on-failure
```

Unit tests (requires Boost.Test):
```bash
cmake .. -DOPTION_BUILD_TESTS=ON
make -j$(nproc)
ctest --test-dir build --output-on-failure
```

## Code style

- Use `.clang-format` (Google-based, 4-space indent, 120 col)
- Only format `src/` files, not `src/3rdparty/`
- Prefer `static_cast<>` over C-style casts
- Use `const` references for read-only parameters
- Add `[[nodiscard]]` to functions returning error codes

## Pull requests

1. Create a feature branch from `main`
2. Ensure `ctest -R integration` passes
3. Keep commits focused and messages concise

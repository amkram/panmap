# Contributing to panmap

## Building from source

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
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

Tests are off by default. Enable them with `-DOPTION_BUILD_TESTS=ON`. There are two
CTest targets: `unit` (a Boost.Test binary, `panmap_tests`) and `e2e` (an
end-to-end pipeline script). Both run against test data in `src/test/data/`.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOPTION_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure        # unit + e2e
ctest --test-dir build -R unit --output-on-failure   # unit only
ctest --test-dir build -R e2e  --output-on-failure   # e2e only
```

See `TESTING.md` for details.

## Code style

- Use `.clang-format` (Google-based, 4-space indent, 120 col)
- Only format `src/` files, not `src/3rdparty/`
- Prefer `static_cast<>` over C-style casts
- Use `const` references for read-only parameters
- Add `[[nodiscard]]` to functions returning error codes

Formatting is checked in CI (the `format` job runs `clang-format --dry-run
--Werror` over first-party sources). Run `clang-format -i` on your changes
before opening a PR.

## Pull requests

1. Create a feature branch from `main`
2. Ensure `ctest --test-dir build --output-on-failure` passes (unit + e2e)
3. Ensure the `format` check passes (`clang-format --dry-run --Werror`)
4. Keep commits focused and messages concise

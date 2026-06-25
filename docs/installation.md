# Installation

## Bioconda (recommended)

```bash
conda install -c bioconda panmap
```

This installs `panmap` and `panmanUtils`.

---

## Docker

```bash
docker pull alanalohaucsc/panmap:latest
docker run --rm alanalohaucsc/panmap:latest --help
```

To run on local files, mount a volume:

```bash
docker run --rm -v $(pwd):/data -w /data alanalohaucsc/panmap:latest \
  ref.panman reads.fq -o sample
```

---

## Building from source

### With conda (recommended)

The repository ships an `environment.yml` with all build dependencies, so nothing
needs to be installed system-wide:

```bash
conda env create -f environment.yml
conda activate panmap
export CPATH="$CONDA_PREFIX/include" LIBRARY_PATH="$CONDA_PREFIX/lib"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=ON
cmake --build build -j
```

`-DUSE_SYSTEM_LIBS=ON` builds against the conda-provided Cap'n Proto, htslib,
Protobuf, and Abseil; the `CPATH`/`LIBRARY_PATH` exports let the bundled aligners
find the conda headers and libraries. The binary is at `build/bin/panmap`.

### With system packages

| Package | Ubuntu/Debian |
|---------|---------------|
| CMake >= 3.14 | `cmake` |
| C++17 compiler | `g++` or `clang++` |
| Protobuf | `protobuf-compiler`, `libprotobuf-dev` |
| Boost | `libboost-program-options-dev`, `libboost-iostreams-dev`, `libboost-filesystem-dev`, `libboost-system-dev`, `libboost-date-time-dev` |
| zlib | `zlib1g-dev` |
| Eigen3 | `libeigen3-dev` |

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The binary is at `build/bin/panmap`. Cap'n Proto and htslib are built from the
bundled sources automatically.

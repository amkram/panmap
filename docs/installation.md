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

### Dependencies

| Package | Ubuntu/Debian |
|---------|---------------|
| CMake >= 3.14 | `cmake` |
| C++17 compiler | `g++` or `clang++` |
| Protobuf | `protobuf-compiler`, `libprotobuf-dev` |
| Boost | `libboost-program-options-dev`, `libboost-iostreams-dev`, `libboost-filesystem-dev`, `libboost-system-dev`, `libboost-date-time-dev` |
| zlib | `zlib1g-dev` |
| htslib | `libhts-dev` |
| Cap'n Proto | `capnproto`, `libcapnp-dev` |
| Eigen3 | `libeigen3-dev` |

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The binary is at `build/bin/panmap`.

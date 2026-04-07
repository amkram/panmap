# Installation

## Bioconda (recommended)

```bash
conda install -c bioconda panmap
```

This installs `panmap` and `panmanUtils`.

---

## Docker

```bash
docker build -t panmap .
docker run --rm panmap panmap -h
```

---

## Building from source

### Dependencies

| Package | Ubuntu/Debian |
|---------|---------------|
| CMake >= 3.12 | `cmake` |
| C++17 compiler | `g++` or `clang++` |
| Protobuf | `protobuf-compiler`, `libprotobuf-dev` |
| Boost | `libboost-program-options-dev`, `libboost-iostreams-dev`, `libboost-filesystem-dev`, `libboost-system-dev`, `libboost-date-time-dev` |
| zlib | `zlib1g-dev` |
| htslib | `libhts-dev` |
| Cap'n Proto | `capnproto`, `libcapnp-dev` |
| Eigen3 | `libeigen3-dev` |

### Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The binary is at `build/bin/panmap`.

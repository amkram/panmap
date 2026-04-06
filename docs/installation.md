# Installation

## Docker (recommended)

```bash
docker build -t panmap .
```

Verify:

```bash
docker run panmap panmap -h
```

Docker bundles all dependencies -- no manual setup required.

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
| Cap'n Proto | `capnproto`, `libcapnp-dev` |
| Eigen3 | `libeigen3-dev` |
| htslib 1.20 | (bundled) |

### Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The binary is at `build/bin/panmap`.

!!! tip
    If you hit dependency issues, the Docker build is the fastest path forward.

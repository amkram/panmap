# Installation

panmap runs on any modern Linux or macOS machine (x86-64 or Apple Silicon) with
about 8 GB of RAM.

## Bioconda (recommended)

```bash
conda install -c conda-forge -c bioconda panmap
panmap -h
```

This installs `panmap` and `panmanUtils`. [Mamba](https://mamba.readthedocs.io/)
(`mamba install ...`) resolves the environment faster.

---

## Docker

[BioContainers](https://biocontainers.pro/) builds a container image for each
Bioconda release:

```bash
docker pull quay.io/biocontainers/panmap:0.1.3--0
docker run --rm quay.io/biocontainers/panmap:0.1.3--0 panmap -h
```

Pick a concrete tag from the
[tags page](https://quay.io/repository/biocontainers/panmap?tab=tags); BioContainers
images are not tagged `latest`. To run on local files, mount the working directory:

```bash
docker run --rm -v "$(pwd):/data" -w /data \
  quay.io/biocontainers/panmap:0.1.3--0 \
  panmap ref.panman reads.fq -o sample
```

You can also build the image from the repository (`~6 min`):

```bash
docker build -t panmap .
docker run --rm panmap panmap -h
```

---

## Building from source

### With conda (recommended)

The repository ships an `environment.yml` with every build and runtime
dependency:

```bash
conda env create -f environment.yml && conda activate panmap
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
build/bin/panmap -h
```

With a conda environment active, the build links against the conda-provided
Cap'n Proto, htslib, Protobuf, Abseil, TBB, and zlib, and passes the
environment's include and library paths to the bundled aligners. No extra flags
or environment variables are needed. The binary is at `build/bin/panmap`.

### Faster builds

Install `ccache` and `ninja` for faster rebuilds:

```bash
conda install -c conda-forge ccache ninja   # or: brew install ccache ninja
```

`ccache` is picked up automatically; unchanged sources are served from cache.
The bundled presets wire up Ninja plus ccache:

```bash
cmake --preset dev && cmake --build --preset dev     # Release + tests
cmake --preset debug && cmake --build --preset debug # -O0, fastest compile
```

### Without conda (system packages)

Install the dependencies below, then build. Cap'n Proto, htslib, and the other
heavier libraries are compiled from the bundled sources, so the system list is
short:

| Dependency | Ubuntu / Debian package |
|------------|-------------------------|
| CMake ≥ 3.14 | `cmake` |
| C++20 compiler | `g++` or `clang++` |
| Protobuf | `protobuf-compiler`, `libprotobuf-dev` |
| Boost | `libboost-program-options-dev`, `libboost-iostreams-dev`, `libboost-filesystem-dev`, `libboost-system-dev`, `libboost-date-time-dev` |
| zlib | `zlib1g-dev` |
| Eigen3 | `libeigen3-dev` |

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
build/bin/panmap -h
```

!!! tip "Build options"
    - `-DUSE_SYSTEM_LIBS=ON` links against system or conda copies of Cap'n Proto,
      htslib, Protobuf, and Abseil instead of building them from source. It is
      on by default inside a conda environment.
    - `-DOPTION_BUILD_TESTS=ON` builds the unit and end-to-end test suites
      (run with `ctest --test-dir build`).

# syntax=docker/dockerfile:1
# Stage 1: Build
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential cmake git wget ca-certificates ccache \
  protobuf-compiler libprotobuf-dev \
  libboost-program-options-dev libboost-iostreams-dev \
  libboost-filesystem-dev libboost-system-dev libboost-date-time-dev \
  zlib1g-dev libhts-dev capnproto libcapnp-dev libeigen3-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /build/panmap
COPY . .

# ccache lives on a BuildKit cache mount kept warm across CI runs by
# buildkit-cache-dance in the workflow, so the from-source deps (abseil, tbb,
# spdlog, jsoncpp, panman) and unchanged objects come from cache -- only changed
# sources recompile. COMPILERCHECK=content hashes the compiler bytes, not mtime,
# so a re-pulled toolchain still hits. OPTION_PORTABLE pins the ISA to x86-64-v3
# so objects cached on one runner never SIGILL on another's CPU (and the shipped
# image runs on any Haswell-or-newer host).
ENV CCACHE_DIR=/ccache CCACHE_COMPILERCHECK=content
RUN --mount=type=cache,target=/ccache \
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOPTION_PORTABLE=ON \
    -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  && cmake --build build -j"$(nproc)" \
  && cmake --install build --prefix /usr/local \
  && ccache -s

# Collect all shared library dependencies for the runtime stage
RUN mkdir -p /runtime-libs && \
  (ldd /usr/local/bin/panmap 2>/dev/null || true) | grep '=>' | awk '{print $3}' | \
  while read lib; do [ -f "$lib" ] && cp -L "$lib" /runtime-libs/; done && \
  find /build/panmap/build -name '*.so*' -exec cp -L {} /runtime-libs/ \;

# Stage 2: Runtime
FROM ubuntu:22.04

COPY --from=builder /usr/local/bin/panmap /usr/local/bin/
COPY --from=builder /runtime-libs/ /usr/lib/

RUN ldconfig && useradd -m panmap
USER panmap
ENTRYPOINT ["panmap"]

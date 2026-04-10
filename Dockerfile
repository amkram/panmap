# Stage 1: Build
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential cmake git wget ca-certificates \
  protobuf-compiler libprotobuf-dev \
  libboost-program-options-dev libboost-iostreams-dev \
  libboost-filesystem-dev libboost-system-dev libboost-date-time-dev \
  zlib1g-dev libhts-dev capnproto libcapnp-dev libeigen3-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /build/panmap
COPY . .
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  && cmake --build build -j$(nproc) \
  && cmake --install build --prefix /usr/local

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

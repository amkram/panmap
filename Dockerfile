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

# Stage 2: Runtime
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
  libboost-program-options1.74.0 libboost-iostreams1.74.0 \
  libboost-filesystem1.74.0 libboost-date-time1.74.0 \
  libprotobuf23 libcapnp-0.8.0 libhts3 zlib1g \
  && rm -rf /var/lib/apt/lists/* \
  && useradd -m panmap

COPY --from=builder /usr/local/bin/panmap /usr/local/bin/

USER panmap
ENTRYPOINT ["panmap"]

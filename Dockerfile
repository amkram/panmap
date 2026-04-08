FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  git \
  wget \
  protobuf-compiler \
  libprotobuf-dev \
  libboost-program-options-dev \
  libboost-iostreams-dev \
  libboost-filesystem-dev \
  libboost-system-dev \
  libboost-date-time-dev \
  zlib1g-dev \
  libhts-dev \
  capnproto \
  libcapnp-dev \
  libeigen3-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /home/panmap

COPY . .

RUN mkdir -p build && cd build && cmake .. && make -j$(nproc)

ENV PATH="/home/panmap/build/bin:${PATH}"

CMD ["panmap"]

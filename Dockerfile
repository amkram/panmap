FROM ubuntu:22.04
WORKDIR /panmap
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update --allow-insecure-repositories

RUN apt-get install -y \
    build-essential \
    wget \
    curl \
    ninja-build \
    libhtscodecs-dev \
    libncurses5-dev \
    libbz2-dev \
    zlib1g-dev \
    libncursesw5-dev \
    liblzma-dev \
    libboost-all-dev \
    cmake \
    libssl-dev \
    unzip \
    zip \
    make \
    git \
    pkg-config \
    coreutils \
    autotools-dev \
    lzma \
    libtool \
    flex \
    bison \
    automake \
    autoconf \
    libprotobuf-dev \
    protobuf-compiler 

COPY . .

ENV CMAKE_BUILD_PARALLEL_LEVEL=4
RUN mkdir build && cd build && cmake .. && make -j4 && make install

ENV LD_LIBRARY_PATH=/usr/local/lib

RUN chmod +x /usr/local/bin/*

CMD ["/usr/local/bin/panmap", "-h"]
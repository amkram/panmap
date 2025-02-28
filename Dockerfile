FROM ubuntu:22.04
WORKDIR /panmap
ARG DEBIAN_FRONTEND=noninteractive
ARG CPUS=8
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
    protobuf-compiler \
    libeigen3-dev \
    libomp-dev \
    libdeflate-dev

COPY ./*.* .
COPY ./src ./src
COPY ./cmake ./cmake
COPY ./dev/examples ./dev/examples
COPY ./dev/simulation_testing ./dev/simulation_testing

# ENV LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:$LD_LIBRARY_PATH

ENV CMAKE_BUILD_PARALLEL_LEVEL=${CPUS}
RUN mkdir build && cd build && cmake -DOPTION_BUILD_SIMULATE=ON .. && cmake --build . --parallel && cmake --install .

RUN chmod +x /usr/local/bin/*

CMD ["/usr/local/bin/panmap", "-h"]

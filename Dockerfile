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
  capnproto \
  libcapnp-dev \
  autoconf \
  automake \
  libtool \
  libeigen3-dev \
  && rm -rf /var/lib/apt/lists/*

# dirty fix for samtools build...
RUN cd /tmp && \
  wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
  tar -xjf htslib-1.20.tar.bz2 && \
  cd htslib-1.20 && \
  ./configure --prefix=/usr/local --disable-lzma --disable-bz2 --disable-libcurl && \
  make -j$(nproc) && \
  make install && \
  ldconfig && \
  cd / && \
  rm -rf /tmp/htslib-1.20*

WORKDIR /panmap

RUN echo 'PS1="\[\033[01;32m\]tiger\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "' >> ~/.bashrc

CMD ["/bin/bash"]
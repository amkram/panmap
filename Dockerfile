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


WORKDIR /home

RUN cd /home && \
  git clone https://github.com/amkram/panmap.git && \
  mkdir -p /home/panmap/build && \
  cd /home/panmap/build && \
  cmake .. && make -j 4 

ENV PATH="/home/panmap/build/bin:${PATH}"

CMD ["panmap"]


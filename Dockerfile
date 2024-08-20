
FROM ubuntu:22.04

WORKDIR /panmap
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y libncurses5-dev libboost-all-dev libssl-dev cmake g++ make git build-essential pkg-config coreutils autotools-dev libtool flex bison automake autoconf

COPY . .

RUN rm -rf build && mkdir build && cd build && cmake -DOPTION_DEBUG=OFF -DCMAKE_INSTALL_PREFIX=/usr/local .. && make -j && make install 
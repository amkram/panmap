
FROM ubuntu:20.04

WORKDIR /performance

COPY . .

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y python3 cmake g++ make git build-essential pkg-config libjsoncpp-dev coreutils autotools-dev libtool flex bison automake autoconf python3-dev python3-pip sudo

RUN git clone --recursive https://github.com/vgteam/vg.git
# RUN cd vg && make -j get-deps
# RUN cd vg && . ./source_me.sh && make -j
# RUN pip install snakemake

ENV PATH="$PATH:/performance/vg/bin"


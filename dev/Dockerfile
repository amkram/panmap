
FROM ubuntu:20.04

RUN mkdir /eval

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y cmake g++ make git build-essential pkg-config coreutils autotools-dev libtool flex bison automake autoconf python3-dev python3-pip sudo

RUN git clone --recursive https://github.com/vgteam/vg.git
RUN cd vg && make -j get-deps
RUN cd vg && . ./source_me.sh && make -j

RUN pip install apples

RUN git clone https://github.com/lpipes/tronko.git
RUN cd tronko/tronko-build && make -j
RUN cd tronko/tronko-assign && make -j

RUN pip install snakemake

ENV PATH="$PATH:/eval/tronko/tronko-build"
ENV PATH="$PATH:/eval/tronko/tronko-assign"
ENV PATH="$PATH:/eval/tronko/bin"
ENV PATH="$PATH:/eval/tronko/scripts"
ENV PATH="$PATH:/eval/vg/bin"


FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /performance

# COPY . .

RUN apt-get update && apt-get install -y python3 jellyfish cmake g++ make git build-essential pkg-config libjsoncpp-dev coreutils autotools-dev libtool flex bison automake autoconf python3-dev python3-pip sudo

RUN git clone --recursive https://github.com/vgteam/vg.git
RUN cd vg && . ./source_me.sh && make -j get-deps
RUN cd vg && . ./source_me.sh && make -j

RUN git clone https://github.com/turakhialab/panman
RUN cd panman && install/installationUbuntu.sh

RUN pip3 install pulp==2.7.0
RUN pip3 install snakemake
RUN pip3 install biopython

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz
RUN tar xvzf seqkit_linux_amd64.tar.gz
RUN mv ./seqkit /usr/local/bin

ENV PATH="$PATH:/performance/vg/bin"
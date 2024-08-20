
FROM ubuntu:20.04

RUN mkdir /eval

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y && apt-get install -y cmake g++ make git build-essential pkg-config coreutils autotools-dev libtool flex bison automake autoconf python3-dev python3-pip
RUN apt-get install -y libboost-all-dev make git build-essential sudo protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake gettext autopoint 
RUN apt-get install -y cereal libtool jq bsdmainutils bc rs parallel npm samtools curl unzip redland-utils librdf-dev cmake pkg-config wget gtk-doc-tools raptor2-utils rasqal-utils bison flex gawk libgoogle-perftools-dev liblz4-dev liblzma-dev libcairo2-dev libpixman-1-dev libffi-dev libcairo-dev libprotobuf-dev libboost-all-dev tabix bcftools libzstd-dev pybind11-dev python3-pybind11 gcc

RUN apt-get upgrade -y

RUN wget https://www.python.org/ftp/python/3.12.4/Python-3.12.4.tgz
RUN tar -xvzf Python-3.12.4.tgz
RUN cd Python-3.12.4 && ./configure --prefix=/usr/local && make && make install

RUN git clone https://github.com/USCiLab/cereal

RUN cd cereal && mkdir build && cd build && \
    cmake .. && \
    make -j && \
    make -j install

RUN git clone --recursive https://github.com/vgteam/vg.git
RUN cd vg && make -j get-deps
RUN cd vg && . ./source_me.sh && make -j

RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz
RUN tar -xvzf jellyfish-2.3.1.tar.gz
RUN cd jellyfish-2.3.1 && ./configure && make -j && make install

RUN git clone https://github.com/eblerjana/pangenie.git  
RUN cd pangenie && mkdir build && cd build && cmake .. && make -j

RUN git clone https://github.com/at-cg/minichain
RUN cd minichain && make -j

RUN git clone https://github.com/matthiasblanke/App-SpaM
RUN cd App-SpaM && mkdir build && cd build && cmake .. && make -j

RUN pip install apples

RUN git clone https://github.com/lpipes/tronko.git
RUN cd tronko/tronko-build && make -j
RUN cd tronko/tronko-assign && make -j

RUN apt-get install -y camlp4-extra \
gawk \
libgsl0-dev \
libsqlite3-dev \
libz-dev \
m4 \
make \
ocaml \
patch



RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.3.5/iqtree-2.3.5-Linux-intel.tar.gz
RUN tar -xvzf iqtree-2.3.5-Linux-intel.tar.gz

RUN git clone https://github.com/yatisht/usher.git
RUN cd usher && install/installUbuntu.sh

RUN pip install snakemake

ENV PATH="$PATH:/eval/iqtree-2.3.5-Linux-intel/bin"
ENV PATH="$PATH:/eval/App-SpaM/build"
ENV PATH="$PATH:/eval/jellyfish-2.3.1/bin"
ENV PATH="$PATH:/eval/minichain"
ENV PATH="$PATH:/eval/pangenie/build/src"
ENV PATH="$PATH:/eval/tronko/tronko-build"
ENV PATH="$PATH:/eval/tronko/tronko-assign"
ENV PATH="$PATH:/eval/tronko/bin"
ENV PATH="$PATH:/eval/tronko/scripts"
ENV PATH="$PATH:/eval/vg/bin"
ENV PATH="$PATH:/eval/usher/build"

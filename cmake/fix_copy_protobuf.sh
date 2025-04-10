#!/bin/bash

# Set up paths
BUILD_DIR="/private/groups/corbettlab/alex/tomG/build"
PANMAN_SRC_DIR="${BUILD_DIR}/_deps/panman-src"
PANMAN_BUILD_DIR="${BUILD_DIR}/_deps/panman-build"

# Create necessary directories
mkdir -p "${PANMAN_SRC_DIR}/src"
mkdir -p "${PANMAN_BUILD_DIR}/src"

# Copy proto files from source to build directory
for proto_file in "${PANMAN_SRC_DIR}/src/"*.proto; do
    if [ -f "$proto_file" ]; then
        cp "$proto_file" "${PANMAN_BUILD_DIR}/src/"
    fi
done

# Copy generated pb files from source to build directory
for pb_file in "${PANMAN_SRC_DIR}/src/"*.pb.*; do
    if [ -f "$pb_file" ]; then
        cp "$pb_file" "${PANMAN_BUILD_DIR}/src/"
    fi
done

# Copy Cap'n Proto files if they exist
for capnp_file in "${PANMAN_SRC_DIR}/src/"*.capnp*; do
    if [ -f "$capnp_file" ]; then
        cp "$capnp_file" "${PANMAN_BUILD_DIR}/src/"
    fi
done

echo "Protocol Buffer files copied successfully" 
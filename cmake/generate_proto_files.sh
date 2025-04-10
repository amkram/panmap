#!/bin/bash

# Set up paths
BUILD_DIR="/private/groups/corbettlab/alex/tomG/build"
PANMAN_SRC_DIR="${BUILD_DIR}/_deps/panman-src"
PANMAN_BUILD_DIR="${BUILD_DIR}/_deps/panman-build"

# Create necessary directories
mkdir -p "${PANMAN_SRC_DIR}/src"
mkdir -p "${PANMAN_BUILD_DIR}/src"

# Function to find protoc
find_protoc() {
    if command -v protoc >/dev/null 2>&1; then
        echo "Found system protoc: $(which protoc)"
        return 0
    else
        echo "Error: protoc not found in system path"
        return 1
    fi
}

# Check for protoc
find_protoc

# Copy proto files if they exist in the root of panman source
if [ -f "${PANMAN_SRC_DIR}/panman.proto" ]; then
    cp "${PANMAN_SRC_DIR}/panman.proto" "${PANMAN_SRC_DIR}/src/"
fi

if [ -f "${PANMAN_SRC_DIR}/usher.proto" ]; then
    cp "${PANMAN_SRC_DIR}/usher.proto" "${PANMAN_SRC_DIR}/src/"
fi

# Generate Protocol Buffer files
if [ -f "${PANMAN_SRC_DIR}/src/panman.proto" ]; then
    protoc --cpp_out="${PANMAN_SRC_DIR}/src" --proto_path="${PANMAN_SRC_DIR}/src" "${PANMAN_SRC_DIR}/src/panman.proto"
    # Copy generated files to build directory
    cp "${PANMAN_SRC_DIR}/src/panman.pb."* "${PANMAN_BUILD_DIR}/src/"
fi

if [ -f "${PANMAN_SRC_DIR}/src/usher.proto" ]; then
    protoc --cpp_out="${PANMAN_SRC_DIR}/src" --proto_path="${PANMAN_SRC_DIR}/src" "${PANMAN_SRC_DIR}/src/usher.proto"
    # Copy generated files to build directory
    cp "${PANMAN_SRC_DIR}/src/usher.pb."* "${PANMAN_BUILD_DIR}/src/"
fi

# Generate Cap'n Proto files if they exist
if [ -f "${PANMAN_SRC_DIR}/src/panman.capnp" ]; then
    if [ -f "${BUILD_DIR}/bin/capnp" ] && [ -f "${BUILD_DIR}/bin/capnpc-c++" ]; then
        cd "${PANMAN_SRC_DIR}/src"
        "${BUILD_DIR}/bin/capnp" compile -o"${BUILD_DIR}/bin/capnpc-c++" --src-prefix="${PANMAN_SRC_DIR}/src" "panman.capnp"
        cd "${BUILD_DIR}"
        # Copy generated files to build directory
        cp "${PANMAN_SRC_DIR}/src/panman.capnp."* "${PANMAN_BUILD_DIR}/src/"
    fi
fi

echo "Protocol Buffer file generation completed successfully" 
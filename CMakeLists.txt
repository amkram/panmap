cmake_minimum_required(VERSION 3.12)

project(panmap VERSION 0.1.0)

# Export compile_commands.json for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project options
option(OPTION_DEBUG "Compile with debug symbols and no optimizations." OFF)
option(OPTION_BUILD_TESTS "Build tests." OFF)
option(OPTION_BUILD_SIMULATE "Build simulate." ON)
option(OPTION_PROFILE "Enable profiling with Valgrind/Callgrind." ON)
option(OPTION_ENABLE_FIXED_KMER "Enable optimized k-mer size implementations (8, 16, 32)." ON)
option(OPTION_ENABLE_LTO "Enable Link Time Optimization" OFF)
option(OPTION_ENABLE_PGO "Enable Profile-Guided Optimization" OFF)
option(OPTION_PORTABLE "Build with maximum compatibility across x86-64 systems" ON)
option(USE_SYSTEM_LIBS "Use system libraries when available instead of building from source" OFF)

# Set consistent compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function -Wno-deprecated-declarations -w")

include(FetchContent)
include(ExternalProject)
include(GNUInstallDirs)

# Optimize FetchContent behavior
set(FETCHCONTENT_QUIET OFF)
set(POSITION_INDEPENDENT_CODE ON)
set(ABSL_PROPAGATE_CXX_STD ON)
# Cache third-party builds to speed up subsequent builds
set(FETCHCONTENT_BASE_DIR "${CMAKE_SOURCE_DIR}/external_cache" CACHE PATH "Cache external dependencies")

# Add the custom cmake directories to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set TBB_DIR to match old build
set(TBB_DIR ${CMAKE_CURRENT_BINARY_DIR}/oneTBB-2019_U9)

# Download and configure TBB
FetchContent_Declare(
    tbb
    URL https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
)
FetchContent_MakeAvailable(tbb)

# Build TBB from source
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libtbb.so
    COMMAND cd ${tbb_SOURCE_DIR} && make -j
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib
    COMMAND cp ${tbb_SOURCE_DIR}/build/*_release/libtbb.so* ${CMAKE_CURRENT_BINARY_DIR}/lib/
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(tbb_build DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib/libtbb.so)

# Set up TBB library path
set(TBB_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libtbb.so")
# Set up TBB include directories
include_directories(${tbb_SOURCE_DIR}/include)

# Download and configure jsoncpp
FetchContent_Declare(
    jsoncpp
    GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
    GIT_TAG 1.9.5
)
FetchContent_MakeAvailable(jsoncpp)

# Pull in spdlog so its target is available before we define panmap
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.11.0
)
FetchContent_MakeAvailable(spdlog)

# Pull in Abseil for flat_hash_map
FetchContent_Declare(
    absl
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
    GIT_TAG 20230802.1 # Use a specific tag for stability
)
FetchContent_MakeAvailable(absl)

# Find Protobuf
find_package(Protobuf REQUIRED)
if(NOT Protobuf_FOUND)
    message(FATAL_ERROR "Protobuf not found")
endif()

# Force the correct protobuf library path
set(Protobuf_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprotobuf.so")
set(Protobuf_INCLUDE_DIRS "/usr/include/google/protobuf")

# Include TBB and jsoncpp headers
include_directories(${tbb_SOURCE_DIR}/include)
include_directories(${jsoncpp_SOURCE_DIR}/include)
include_directories(${Protobuf_INCLUDE_DIRS})

# Setup for CapnProto using system-installed tools
find_program(CAPNP_EXECUTABLE NAMES capnp)
find_program(CAPNPC_CXX_EXECUTABLE NAMES capnpc-c++)

if(NOT CAPNP_EXECUTABLE OR NOT CAPNPC_CXX_EXECUTABLE)
    message(FATAL_ERROR "CapnProto executables not found. Please install CapnProto.")
endif()

message(STATUS "Using CapnProto tools: ${CAPNP_EXECUTABLE}, ${CAPNPC_CXX_EXECUTABLE}")

# Create a dummy target for consistency with the rest of the build
add_custom_target(capnp_tools ALL)

# Set up CapnProto include directories
include_directories(/usr/local/include)

# Generate CapnProto files from index.capnp
set(CAPNP_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/index.capnp")
set(CAPNP_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(CAPNP_GEN_CPP "${CAPNP_GEN_DIR}/index.capnp.c++")
set(CAPNP_GEN_H "${CAPNP_GEN_DIR}/index.capnp.h")

# Custom command to generate CapnProto files
add_custom_command(
    OUTPUT ${CAPNP_GEN_CPP} ${CAPNP_GEN_H}
    COMMAND ${CAPNP_EXECUTABLE}
        compile
        -o${CAPNPC_CXX_EXECUTABLE}
        --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
        ${CAPNP_SRC_FILE}
    DEPENDS ${CAPNP_SRC_FILE} capnp_tools
    COMMENT "Generating C++ files from index.capnp"
    VERBATIM
)


set(MGSR_CAPNP_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/mgsr_index.capnp")
set(MGSR_CAPNP_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(MGSR_CAPNP_GEN_CPP "${MGSR_CAPNP_GEN_DIR}/mgsr_index.capnp.c++")
set(MGSR_CAPNP_GEN_H "${MGSR_CAPNP_GEN_DIR}/mgsr_index.capnp.h")
# Custom command to generate CapnProto files
add_custom_command(
    OUTPUT ${MGSR_CAPNP_GEN_CPP} ${MGSR_CAPNP_GEN_H}
    COMMAND ${CAPNP_EXECUTABLE}
        compile
        -o${CAPNPC_CXX_EXECUTABLE}
        --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
        ${MGSR_CAPNP_SRC_FILE}
    DEPENDS ${MGSR_CAPNP_SRC_FILE} capnp_tools
    COMMENT "Generating C++ files from mgsr_index.capnp"
    VERBATIM
)

# Set variables to point to the generated files
set(CAPNP_SRCS ${CAPNP_GEN_CPP})
set(CAPNP_HDRS ${CAPNP_GEN_H})
set(MGSR_CAPNP_SRCS ${MGSR_CAPNP_GEN_CPP})
set(MGSR_CAPNP_HDRS ${MGSR_CAPNP_GEN_H})

# Create a target that depends on the generated files
add_custom_target(index_capnp_generator ALL
    DEPENDS ${CAPNP_GEN_CPP} ${CAPNP_GEN_H}
    DEPENDS ${MGSR_CAPNP_GEN_CPP} ${MGSR_CAPNP_GEN_H}
)

# Create a simplified version of the capnp_generate_cpp macro
function(capnp_generate_cpp SOURCES HEADERS)
  if(NOT ARGN)
    message(SEND_ERROR "capnp_generate_cpp: missing schema files")
    return()
endif()

  set(output_sources "")
  set(output_headers "")

  foreach(schema_file ${ARGN})
    get_filename_component(file_name ${schema_file} NAME_WE)
    
    set(schema_cpp "${CMAKE_CURRENT_BINARY_DIR}/${file_name}.capnp.c++")
    set(schema_h "${CMAKE_CURRENT_BINARY_DIR}/${file_name}.capnp.h")
    
    list(APPEND output_sources "${schema_cpp}")
    list(APPEND output_headers "${schema_h}")

add_custom_command(
      OUTPUT "${schema_cpp}" "${schema_h}"
      COMMAND ${CAPNP_EXECUTABLE}
          compile 
          -o${CAPNPC_CXX_EXECUTABLE}
          --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
          ${schema_file}
      DEPENDS "${schema_file}" capnp_tools
      COMMENT "Compiling Cap'n Proto schema ${schema_file}"
      VERBATIM
    )
  endforeach()
  
  set(${SOURCES} ${output_sources} PARENT_SCOPE)
  set(${HEADERS} ${output_headers} PARENT_SCOPE)
endfunction()

# ===== Setup for Panman =====
# First, create the necessary config files for panman
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmake)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmake/TBBConfig.cmake "
# TBB Config File
set(TBB_FOUND TRUE)
set(TBB_INCLUDE_DIRS \"${TBB_DIR}/include\")
set(TBB_LIBRARIES \"${TBB_LIBRARY}\")
set(TBB_LIBRARY_DIRS \"/usr/lib/x86_64-linux-gnu\")

if(NOT TARGET TBB::tbb)
    add_library(TBB::tbb SHARED IMPORTED)
    set_target_properties(TBB::tbb PROPERTIES
        IMPORTED_LOCATION \"${TBB_LIBRARY}\"
        INTERFACE_INCLUDE_DIRECTORIES \"${TBB_DIR}/include\"
    )
endif()

")

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmake/jsoncppConfig.cmake "
# JsonCpp Config File
set(JSONCPP_FOUND TRUE)
set(JSONCPP_INCLUDE_DIRS \"${CMAKE_CURRENT_BINARY_DIR}/include\")
set(JSONCPP_LIBRARIES \"jsoncpp\")

if(NOT TARGET jsoncpp_lib)
    add_library(jsoncpp_lib INTERFACE IMPORTED)
endif()

if(NOT TARGET JsonCpp::JsonCpp)
    add_library(JsonCpp::JsonCpp INTERFACE IMPORTED)
endif()

")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake/TBBBuild.cmake "
function(tbb_build)
    # Empty implementation
endfunction()
")

# Create a dummy JsonCpp library
add_library(jsoncpp INTERFACE)

# Fetch the panman source code (clone only, disable update)
FetchContent_Declare(
    panman
    GIT_REPOSITORY https://github.com/TurakhiaLab/panman
    GIT_TAG main
    GIT_SHALLOW TRUE
    UPDATE_COMMAND ""
)

# Force the population of the content
FetchContent_GetProperties(panman)
if(NOT panman_POPULATED)
    FetchContent_Populate(panman)
    
    # Get the source directory
    FetchContent_GetProperties(panman SOURCE_DIR panman_SOURCE_DIR BINARY_DIR panman_BINARY_DIR)
    
    # Generate version.hpp file for panman
    file(WRITE "${panman_SOURCE_DIR}/src/version.hpp" 
        "#ifndef VERSION_HPP\n"
        "#define VERSION_HPP\n"
        "\n"
        "#define PROJECT_VERSION \"1.0.0\"\n"
        "#define PMAT_VERSION \"1.0.0\"\n"
        "\n"
        "#endif // VERSION_HPP\n"
    )
    
    # Create a patch file for panman's CMakeLists.txt
    file(WRITE "${panman_SOURCE_DIR}/CMakeLists.txt.orig" "")
    file(READ "${panman_SOURCE_DIR}/CMakeLists.txt" PANMAN_CMAKE_CONTENT)
    file(WRITE "${panman_SOURCE_DIR}/CMakeLists.txt.orig" "${PANMAN_CMAKE_CONTENT}")
    
    # Modify the panman CMakeLists.txt
    file(WRITE "${panman_SOURCE_DIR}/CMakeLists.txt" [=[
cmake_minimum_required (VERSION 3.8)
project(panmanUtils)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function -Wno-deprecated-declarations -w")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Protobuf
find_package(Protobuf REQUIRED)
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
  endif()

# Force the correct protobuf library path
set(Protobuf_LIBRARIES "/usr/lib/x86_64-linux-gnu/libprotobuf.so")
set(Protobuf_INCLUDE_DIRS "/usr/include/google/protobuf")

include_directories(${Protobuf_INCLUDE_DIRS})
include_directories("/usr/include")

# Generate Protocol Buffers
set(PROTO_FILES 
    "${CMAKE_CURRENT_SOURCE_DIR}/panman.proto"
    "${CMAKE_CURRENT_SOURCE_DIR}/usher.proto"
)

protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${PROTO_FILES})

# Add proto files to include path
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Boost
find_package(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Define source files
set(PANMAN_SOURCES
    src/panman.cpp
    ${PROTO_SRCS}
)

# Build both a shared library and an executable
add_library(panman_lib STATIC ${PANMAN_SOURCES})
set_target_properties(panman_lib PROPERTIES OUTPUT_NAME "panman")
target_compile_options(panman_lib PRIVATE -fPIC)

# Define the panmanUtils target with actual sources
add_executable(panmanUtils
    src/panmanUtils.cpp
    ${PROTO_SRCS}
)

# Generate C++ from panman.capnp
set(PANMAN_CAPNP_SRC "${CMAKE_CURRENT_SOURCE_DIR}/panman.capnp")
set(PANMAN_CAPNP_H "${CMAKE_CURRENT_BINARY_DIR}/panman.capnp.h")
set(PANMAN_CAPNP_CPP "${CMAKE_CURRENT_BINARY_DIR}/panman.capnp.c++")

add_custom_command(
    OUTPUT "${PANMAN_CAPNP_H}" "${PANMAN_CAPNP_CPP}"
    COMMAND ${CAPNP_EXECUTABLE} compile -o${CAPNPC_CXX_EXECUTABLE} 
            --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
            ${PANMAN_CAPNP_SRC}
            -I${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS "${PANMAN_CAPNP_SRC}"
    COMMENT "Generating C++ from panman.capnp"
)

add_custom_target(generate_panman_capnp_extern ALL
    DEPENDS "${PANMAN_CAPNP_H}" "${PANMAN_CAPNP_CPP}"
)

# Add the generated files to both targets
target_sources(panman_lib PRIVATE "${PANMAN_CAPNP_CPP}")
target_sources(panmanUtils PRIVATE "${PANMAN_CAPNP_CPP}")

add_dependencies(panman_lib generate_panman_capnp_extern)
add_dependencies(panmanUtils generate_panman_capnp_extern)

# Link dependencies for the library
target_link_libraries(panman_lib
    capnp
    kj
    jsoncpp_lib
    ${Protobuf_LIBRARIES}
    ${Boost_LIBRARIES}
    ${TBB_LIBRARY}
    pthread
)

# Link dependencies for the executable (including the library)
target_link_libraries(panmanUtils
    panman_lib
    capnp
    kj
    jsoncpp_lib
    ${Protobuf_LIBRARIES}
    ${Boost_LIBRARIES}
    ${TBB_LIBRARY}
    pthread
)

# Install targets
install(TARGETS panmanUtils panman_lib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)
]=])

    # Add the subdirectory to the build
    add_subdirectory(${panman_SOURCE_DIR} ${panman_BINARY_DIR})
    
    # Add an explicit dependency to ensure the CapnProto tools are built first
    add_dependencies(generate_panman_capnp_extern capnp_tools)
    
    # Set up the include paths
    target_include_directories(panmanUtils PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/include
        ${CAPNP_BUILD_DIR}/include
        ${Boost_INCLUDE_DIRS}
    )
    
    # Create a dummy library to link against for our main project
    add_library(panman INTERFACE)
    target_include_directories(panman INTERFACE
        ${panman_SOURCE_DIR}/src
        ${panman_BINARY_DIR}
        ${CAPNP_BUILD_DIR}/include
    )
    
    # Create panman_capnp_obj target
    add_library(panman_capnp_obj OBJECT IMPORTED)
    
    message(STATUS "Panman source directory: ${panman_SOURCE_DIR}")
    message(STATUS "Panman binary directory: ${panman_BINARY_DIR}")
endif()

# Add dependency on panman_lib in _deps/panman-build/
add_dependencies(panman panman_lib)

# Define the Minimap2 build target
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 && make clean && make -j
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/bin
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/minimap2 ${CMAKE_CURRENT_BINARY_DIR}/bin/minimap2
    COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/*.o
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create a minimap2 target that depends on the library build
add_custom_target(minimap2_target 
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMENT "Building minimap2 library from src/3rdparty/minimap2"
)

# Build htslib using ExternalProject
set(HTSLIB_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib_install)
ExternalProject_Add(htslib_ext
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${HTSLIB_INSTALL_DIR} --disable-lzma --disable-bz2 --disable-libcurl
    BUILD_COMMAND make -j
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1 # htslib needs to be built in source
    UPDATE_COMMAND "" # Don't try to update the source
)
# Create an imported target for htslib
add_library(htslib STATIC IMPORTED GLOBAL)
set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION ${HTSLIB_INSTALL_DIR}/lib/libhts.a
)
add_dependencies(htslib htslib_ext)

# Build samtools using ExternalProject
set(SAMTOOLS_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/samtools_install)
ExternalProject_Add(samtools_ext
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${SAMTOOLS_INSTALL_DIR} --with-htslib=${HTSLIB_INSTALL_DIR} --disable-lzma --disable-bz2 --disable-libcurl --disable-shared --with-pic
    # Explicitly build only the libst.a target
    BUILD_COMMAND     ${CMAKE_COMMAND} -E env SHELL=/bin/bash make libst.a
    # Explicitly set INSTALL_COMMAND to nothing to prevent default make install
    INSTALL_COMMAND   ""
    BUILD_IN_SOURCE 1 # samtools needs to be built in source
    UPDATE_COMMAND "" # Don't try to update the source
    DEPENDS htslib_ext # Ensure htslib is built first
)
# Create an imported target for samtools
add_library(samtools STATIC IMPORTED GLOBAL)
# Build directory structure for in-source builds with ExternalProject
set(SAMTOOLS_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/libst.a") 
set_target_properties(samtools PROPERTIES
    IMPORTED_LOCATION ${SAMTOOLS_LIB_PATH}
)
add_dependencies(samtools samtools_ext)

# Define BWA sources
set(BWA_SOURCES
    src/3rdparty/bwa/bwa.c
    src/3rdparty/bwa/run.c
    src/3rdparty/bwa/bwase.c
    src/3rdparty/bwa/bwaseqio.c
    src/3rdparty/bwa/bwt.c
    src/3rdparty/bwa/bwtaln.c
    src/3rdparty/bwa/bwtindex.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/bwtsw2_core.c
    src/3rdparty/bwa/bwtsw2_main.c
    src/3rdparty/bwa/bwtsw2_aux.c
    src/3rdparty/bwa/bwt_lite.c
    src/3rdparty/bwa/bwtsw2_chain.c
    src/3rdparty/bwa/fastmap.c
    src/3rdparty/bwa/bwtsw2_pair.c
    src/3rdparty/bwa/utils.c
    src/3rdparty/bwa/bwape.c
    src/3rdparty/bwa/kopen.c
    src/3rdparty/bwa/pemerge.c
    src/3rdparty/bwa/maxk.c
    src/3rdparty/bwa/bwashm.c
    src/3rdparty/bwa/bntseq.c
    src/3rdparty/bwa/is.c
    src/3rdparty/bwa/bwamem.c
    src/3rdparty/bwa/bwamem_pair.c
    src/3rdparty/bwa/bwamem_extra.c
    src/3rdparty/bwa/malloc_wrap.c
    src/3rdparty/bwa/QSufSort.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/rope.c
    src/3rdparty/bwa/kstring.c
    src/3rdparty/bwa/ksw.c
    src/3rdparty/bwa/bwt.c
    src/3rdparty/bwa/bamlite.c
    src/3rdparty/bwa/bntseq.c
    src/3rdparty/bwa/bwashm.c
    src/3rdparty/bwa/bwase.c
    src/3rdparty/bwa/bwaseqio.c
    src/3rdparty/bwa/bwt_lite.c
    src/3rdparty/bwa/bwtaln.c
    src/3rdparty/bwa/bwtindex.c
    src/3rdparty/bwa/is.c
    src/3rdparty/bwa/kopen.c
    src/3rdparty/bwa/kstring.c
    src/3rdparty/bwa/ksw.c
    src/3rdparty/bwa/utils.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/rope.c
    src/3rdparty/bwa/rle.c
    src/3rdparty/bwa/bwtsw2_aux.c
    src/3rdparty/bwa/bwtsw2_chain.c
    src/3rdparty/bwa/bwtsw2_main.c
    src/3rdparty/bwa/bwtsw2_pair.c
    src/3rdparty/bwa/bwtgap.c
    src/mm_align.c
)

# Aggregate all the dependencies into a single target
add_custom_target(build-deps DEPENDS
    htslib_ext
    samtools_ext
    minimap2_target
    capnp_tools
    tbb_build
)

# Create imported targets for the dependencies
# kj and capnp_tools are already provided by FetchContent_MakeAvailable(capnp_tools)
# There's no need to create imported targets for them

add_library(minimap2 STATIC IMPORTED)
set_target_properties(minimap2 PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a")
add_dependencies(minimap2 minimap2_target)

add_library(bwa STATIC ${BWA_SOURCES})
target_include_directories(bwa PUBLIC src/3rdparty/bwa)
target_compile_options(bwa PRIVATE -fPIC)

file(GLOB BCFTOOLS_SOURCES src/3rdparty/bcftools/*.c)
add_library(bcftools STATIC ${BCFTOOLS_SOURCES})
target_include_directories(bcftools PUBLIC src/3rdparty/bcftools)
target_compile_options(bcftools PRIVATE -fPIC)

# Define panmap executable sources using the correct list
set(PANMAP_SOURCES_MAIN
    src/main.cpp
    src/mgsr.cpp
    src/conversion.cpp
    src/panmap_utils.cpp
    src/genotyping.cpp
    src/gap_map.cpp
    src/indexing.cpp
    src/alignment.cpp
    src/placement.cpp
    src/seeding.cpp
    src/mm_align.c
    src/pileup.c
    src/performance.cpp
    src/state.cpp
    src/state_helpers.cpp
)

# Find required packages before using them
find_package(Boost COMPONENTS program_options iostreams filesystem date_time system REQUIRED)
find_package(ZLIB REQUIRED)

add_executable(panmap ${PANMAP_SOURCES_MAIN})
add_dependencies(panmap build-deps index_capnp_generator tbb_build)

# Build simulate tool if enabled
if(OPTION_BUILD_SIMULATE)
    add_executable(simulate 
        src/test/simulate.cpp
        src/simulate_utils.cpp
    )
    
    add_dependencies(simulate panman_lib)
    
    target_include_directories(simulate PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/include
        ${panman_SOURCE_DIR}/src
        ${panman_BINARY_DIR}
    )
    
    target_link_libraries(simulate PRIVATE
        panman_lib
        ${Boost_LIBRARIES}
        ${PROTOBUF_LIBRARIES}
        spdlog::spdlog
        absl::flat_hash_map
        absl::flat_hash_set
    )
    
    # Set output directory to build/bin
    set_target_properties(simulate PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
endif()

# Add generated sources to panmap target
  target_sources(panmap PRIVATE ${CAPNP_SRCS} ${MGSR_CAPNP_SRCS})

# Add include directories
target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(panmap PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(panmap PUBLIC ${HTSLIB_INSTALL_DIR}/include)

# Link spdlog (provides its own includes)
target_link_libraries(panmap PRIVATE spdlog::spdlog)

# Add other FetchContent dependencies
FetchContent_Declare(
    docopt
    GIT_REPOSITORY https://github.com/docopt/docopt.cpp
    GIT_TAG master
)
FetchContent_MakeAvailable(docopt)
if(TARGET docopt)
    target_compile_options(docopt PRIVATE -w -m64 -mtune=generic)
endif()

# Pull in libdeflate
FetchContent_Declare(deflate
    GIT_REPOSITORY https://github.com/ebiggers/libdeflate
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(deflate)

# Link libraries
target_link_libraries(panmap PRIVATE 
    panman
    bwa
    bcftools
    minimap2 
    samtools 
    htslib 
    docopt
    # Link system capnp libs directly
    /usr/local/lib/libcapnp.a 
    /usr/local/lib/libkj.a 
    jsoncpp_lib
    ${Protobuf_LIBRARIES}
    Boost::program_options
    Boost::iostreams
    Boost::filesystem
    Boost::system
    ZLIB::ZLIB
    ${TBB_LIBRARY}
    spdlog::spdlog  # Add spdlog
    # Link panman_lib directly to access the implementation of panmanUtils functions
    panman_lib
    # Fix libdeflate linking
    libdeflate_shared
    # Link Abseil components
    absl::flat_hash_map
    m
    pthread
)

# Set compiler options based on build type
if(OPTION_DEBUG)
    target_compile_options(panmap PUBLIC -g -O0 -ggdb)
elseif(OPTION_PROFILE)
    # Add -pg flag for gprof profiling
    target_compile_options(panmap PUBLIC -g -O2  -fno-omit-frame-pointer)
    target_link_options(panmap PUBLIC  -fno-omit-frame-pointer)
    
    # Also enable profiling for panman_lib
    if(TARGET panman_lib)
        target_compile_options(panman_lib PUBLIC -g -O2 -fno-omit-frame-pointer)
        target_link_options(panman_lib PUBLIC -fno-omit-frame-pointer)
    endif()
else()
    target_compile_options(panmap PUBLIC -g -O3)
    
    # Enable Link Time Optimization if requested
    if(OPTION_ENABLE_LTO)
        target_compile_options(panmap PUBLIC -flto)
        target_link_options(panmap PUBLIC -flto)
        message(STATUS "Link Time Optimization (LTO): Enabled")
    endif()
    
    # Enable vectorization explicitly
    target_compile_options(panmap PUBLIC -ftree-vectorize)
endif()

# Unit testing data
configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
configure_file(src/test/data/test.nwk test.nwk COPYONLY)
configure_file(src/test/data/test.json test.json COPYONLY)
configure_file(src/test/data/test.pmat test.pmat COPYONLY)
configure_file(src/test/data/test.fastq test.fastq COPYONLY)

install(TARGETS panmap RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/ DESTINATION ${CMAKE_INSTALL_LIBDIR} FILES_MATCHING PATTERN "*.a" PATTERN "*.so" PATTERN "*.dylib" PATTERN "*.la" PATTERN "*.so.*")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin/ DESTINATION ${CMAKE_INSTALL_BINDIR} FILES_MATCHING PATTERN "*")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Create a compile_commands.json symlink for clangd
add_custom_command(
    TARGET panmap POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_SOURCE_DIR}/compile_commands.json
    COMMENT "Creating compile_commands.json symlink for clangd"
)

# Print build configuration summary
message(STATUS "")
message(STATUS "Build Configuration Summary:")
message(STATUS "----------------------------")
message(STATUS "CMake Version: ${CMAKE_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "")
message(STATUS "Optimizations:")
message(STATUS "  Fixed k-mer optimization: ${OPTION_ENABLE_FIXED_KMER}")
message(STATUS "  Debug mode: ${OPTION_DEBUG}")
message(STATUS "  Profiling: ${OPTION_PROFILE}")
message(STATUS "  Link Time Optimization: ${OPTION_ENABLE_LTO}")
message(STATUS "  Profile-Guided Optimization: ${OPTION_ENABLE_PGO}")
message(STATUS "")
message(STATUS "Components:")
message(STATUS "  Build tests: ${OPTION_BUILD_TESTS}")
message(STATUS "  Build simulate: ${OPTION_BUILD_SIMULATE}")
message(STATUS "----------------------------")
message(STATUS "")

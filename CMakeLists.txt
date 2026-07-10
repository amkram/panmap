cmake_minimum_required(VERSION 3.14)

project(panmap VERSION 0.1.2)

# Use ccache as compiler launcher when installed: unchanged translation units
# (notably mgsr.cpp, ~25s at -O3) are served from cache. No-op if ccache is
# absent; override with -DCMAKE_CXX_COMPILER_LAUNCHER=.
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM AND NOT DEFINED CMAKE_CXX_COMPILER_LAUNCHER)
    message(STATUS "ccache found: ${CCACHE_PROGRAM} (used as compiler launcher)")
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/version.h"
    @ONLY
)

# Suppress deprecation warnings from third-party dependencies
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)  # Allow FetchContent_Populate (needed for panman patching)
endif()
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)  # Use config-mode for Protobuf/Boost
    set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
endif()

# Export compile_commands.json for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

option(OPTION_DEBUG "Compile with debug symbols and no optimizations." OFF)
option(OPTION_BUILD_TESTS "Build tests." OFF)
option(OPTION_BUILD_SIMULATE "Build simulate." ON)
option(OPTION_PORTABLE "Build portable binary (x86-64-v3: AVX2/FMA, Haswell 2013+)" OFF)
# Default to system libs inside an active conda env so a plain cmake there builds
# without extra flags; otherwise build the bundled sources.
if(DEFINED ENV{CONDA_PREFIX})
    set(_use_system_libs_default ON)
else()
    set(_use_system_libs_default OFF)
endif()
option(USE_SYSTEM_LIBS "Use system libraries when available instead of building from source" ${_use_system_libs_default})

# Point the build at an active conda env ONLY when using its libraries. Injecting
# its include/lib dirs unconditionally leaks headers (e.g. an incompatible protobuf)
# into a bundled (USE_SYSTEM_LIBS=OFF) build and breaks it.
if(USE_SYSTEM_LIBS AND DEFINED ENV{CONDA_PREFIX})
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
    include_directories(SYSTEM "$ENV{CONDA_PREFIX}/include")
    link_directories("$ENV{CONDA_PREFIX}/lib")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pipe -g -ggdb3 -Wno-unused-function -Wno-deprecated-declarations -w")

include(FetchContent)
include(ExternalProject)
include(GNUInstallDirs)

set(FETCHCONTENT_QUIET OFF)
set(POSITION_INDEPENDENT_CODE ON)
set(ABSL_PROPAGATE_CXX_STD ON)
set(FETCHCONTENT_BASE_DIR "${CMAKE_SOURCE_DIR}/external_cache" CACHE PATH "Cache external dependencies")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# TBB: use the system library under USE_SYSTEM_LIBS (conda) or on macOS
# (oneTBB 2019_U9 does not support ARM); otherwise build it from source.
if(USE_SYSTEM_LIBS OR APPLE)
    find_package(TBB REQUIRED)
    set(TBB_LIBRARY TBB::tbb)
    add_custom_target(tbb_build) # no-op target for dependency compatibility
else()
    set(TBB_DIR ${CMAKE_CURRENT_BINARY_DIR}/oneTBB-2019_U9)

    FetchContent_Declare(
        tbb
        URL https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
    )
    FetchContent_MakeAvailable(tbb)

    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libtbb.so
        COMMAND ${CMAKE_COMMAND} -E env "CXXFLAGS=-Wno-changes-meaning" make -C ${tbb_SOURCE_DIR} -j compiler=gcc
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib
        COMMAND cp ${tbb_SOURCE_DIR}/build/*_release/libtbb.so* ${CMAKE_CURRENT_BINARY_DIR}/lib/
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_custom_target(tbb_build DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib/libtbb.so)

    set(TBB_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libtbb.so")
    include_directories(${tbb_SOURCE_DIR}/include)
endif()

if(USE_SYSTEM_LIBS)
    find_package(jsoncpp REQUIRED)  # provides imported target jsoncpp_lib
else()
    set(JSONCPP_WITH_TESTS OFF CACHE BOOL "" FORCE)
    set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        jsoncpp
        GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
        GIT_TAG 1.9.5
    )
    FetchContent_MakeAvailable(jsoncpp)
endif()

if(USE_SYSTEM_LIBS)
    find_package(spdlog REQUIRED)
else()
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.11.0
    )
    FetchContent_MakeAvailable(spdlog)
endif()

if(USE_SYSTEM_LIBS OR APPLE)
    find_package(absl REQUIRED)
else()
    FetchContent_Declare(
        absl
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 20230802.1
    )
    FetchContent_MakeAvailable(absl)
endif()

find_package(Protobuf REQUIRED)
if(NOT Protobuf_FOUND)
    message(FATAL_ERROR "Protobuf not found")
endif()
# Use find_library to get a concrete path (avoids IMPORTED target scope issues)
set(PROTOBUF_LIB_PATH ${PROTOBUF_LIBRARY})

# TBB and jsoncpp headers only when built from source; imported targets
# carry their own usage-requirement includes.
if(NOT USE_SYSTEM_LIBS AND NOT APPLE)
    include_directories(${tbb_SOURCE_DIR}/include)
endif()
if(NOT USE_SYSTEM_LIBS)
    include_directories(${jsoncpp_SOURCE_DIR}/include)
endif()
include_directories(SYSTEM ${Protobuf_INCLUDE_DIRS})

# Both branches define the same vars/target (CAPNP_LIBRARY, KJ_LIBRARY, CAPNP_EXECUTABLE,
# CAPNPC_CXX_EXECUTABLE, CAPNP_SCHEMA_DIR, and a 'capnp_tools' target) consumed below.
if(USE_SYSTEM_LIBS)
    find_library(CAPNP_LIBRARY NAMES capnp REQUIRED)
    find_library(KJ_LIBRARY NAMES kj REQUIRED)
    find_program(CAPNP_EXECUTABLE NAMES capnp REQUIRED)
    find_program(CAPNPC_CXX_EXECUTABLE NAMES capnpc-c++ REQUIRED)
    find_path(CAPNP_SCHEMA_DIR capnp/c++.capnp REQUIRED)
    add_custom_target(capnp_tools)  # no-op: tools come from the system/conda prefix
    message(STATUS "Using system Cap'n Proto: ${CAPNP_LIBRARY}")
    include_directories(${CAPNP_SCHEMA_DIR})
else()
    set(CAPNP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/capnproto)
    set(CAPNP_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/capnp_install)
    ExternalProject_Add(capnp_tools
        SOURCE_DIR ${CAPNP_SOURCE_DIR}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${CAPNP_INSTALL_DIR}
            -DCMAKE_INSTALL_LIBDIR=lib  # force lib/ (not lib64) so the paths below resolve on all distros
            -DCMAKE_BUILD_TYPE=Release
            -DBUILD_TESTING=OFF
        BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release -j
        INSTALL_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
    )
    set(CAPNP_EXECUTABLE ${CAPNP_INSTALL_DIR}/bin/capnp)
    set(CAPNPC_CXX_EXECUTABLE ${CAPNP_INSTALL_DIR}/bin/capnpc-c++)
    set(CAPNP_SCHEMA_DIR ${CAPNP_INSTALL_DIR}/include)
    set(CAPNP_LIBRARY ${CAPNP_INSTALL_DIR}/lib/libcapnp.a)
    set(KJ_LIBRARY ${CAPNP_INSTALL_DIR}/lib/libkj.a)
    message(STATUS "Building CapnProto from bundled source")
    include_directories(${CAPNP_INSTALL_DIR}/include)
endif()

# index_lite.capnp is the only capnp schema used.
set(LITE_CAPNP_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/index_lite.capnp")
set(LITE_CAPNP_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(LITE_CAPNP_GEN_CPP "${LITE_CAPNP_GEN_DIR}/index_lite.capnp.c++")
set(LITE_CAPNP_GEN_H "${LITE_CAPNP_GEN_DIR}/index_lite.capnp.h")
file(MAKE_DIRECTORY ${LITE_CAPNP_GEN_DIR})
add_custom_command(
    OUTPUT ${LITE_CAPNP_GEN_CPP} ${LITE_CAPNP_GEN_H}
    COMMAND ${CAPNP_EXECUTABLE}
        compile
        -o${CAPNPC_CXX_EXECUTABLE}
        -I${CAPNP_SCHEMA_DIR}
        --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}/src
        ${LITE_CAPNP_SRC_FILE}
    WORKING_DIRECTORY ${LITE_CAPNP_GEN_DIR}
    DEPENDS ${LITE_CAPNP_SRC_FILE} capnp_tools
    COMMENT "Generating C++ files from index_lite.capnp"
    VERBATIM
)

set(LITE_CAPNP_SRCS ${LITE_CAPNP_GEN_CPP})
set(LITE_CAPNP_HDRS ${LITE_CAPNP_GEN_H})

add_custom_target(index_capnp_generator ALL
    DEPENDS ${LITE_CAPNP_GEN_CPP} ${LITE_CAPNP_GEN_H}
)

# Simplified reimplementation of capnp_generate_cpp.
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

# ===== Panman =====
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmake)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmake/TBBConfig.cmake "
# TBB Config File
set(TBB_FOUND TRUE)
set(TBB_INCLUDE_DIRS \"${TBB_DIR}/include\")
set(TBB_LIBRARIES \"${TBB_LIBRARY}\")
set(TBB_LIBRARY_DIRS \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")

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

# Dummy JsonCpp library
add_library(jsoncpp INTERFACE)

# panman is built from source either way (no conda package).
if(USE_SYSTEM_LIBS OR APPLE)
    set(PANMAN_TBB_COMPAT "true")
else()
    set(PANMAN_TBB_COMPAT "false")
endif()

if(PANMAN_SOURCE_DIR OR EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/panman")
    # Source pre-placed (e.g. conda recipe's second tarball) or provided via
    # -DPANMAN_SOURCE_DIR; no network fetch.
    if(NOT PANMAN_SOURCE_DIR)
        set(PANMAN_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/panman")
    endif()
    set(panman_SOURCE_DIR "${PANMAN_SOURCE_DIR}")
    set(panman_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/panman-build")
    execute_process(COMMAND sed -i.bak "15d" "${panman_SOURCE_DIR}/src/panmanUtils.cpp")
    execute_process(COMMAND sed -i.bak "s/PROJECT_VERSION/\"1.0\"/g" "${panman_SOURCE_DIR}/src/panmanUtils.cpp")
    execute_process(COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_panman.sh "${panman_SOURCE_DIR}" ${PANMAN_TBB_COMPAT})
else()
FetchContent_Declare(
    panman
    GIT_REPOSITORY https://github.com/TurakhiaLab/panman
    GIT_TAG f6e44c9f43b6804f466665742d5a979a2675982e # using the old version for now until panman implements backward compatibility with published trees
    GIT_SHALLOW False
    UPDATE_COMMAND ""
    PATCH_COMMAND sed -i.bak "15d" <SOURCE_DIR>/src/panmanUtils.cpp &&
                  sed -i.bak "s/PROJECT_VERSION/\\\"1.0\\\"/g" <SOURCE_DIR>/src/panmanUtils.cpp &&
                  sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_panman.sh <SOURCE_DIR> ${PANMAN_TBB_COMPAT}
)
    FetchContent_GetProperties(panman)
    if(NOT panman_POPULATED)
        FetchContent_Populate(panman)
        FetchContent_GetProperties(panman SOURCE_DIR panman_SOURCE_DIR BINARY_DIR panman_BINARY_DIR)
    endif()
endif()
    
    file(WRITE "${panman_SOURCE_DIR}/src/version.hpp"
        "#ifndef VERSION_HPP\n"
        "#define VERSION_HPP\n"
        "\n"
        "#define PROJECT_VERSION \"1.0.0\"\n"
        "#define PMAT_VERSION \"1.0.0\"\n"
        "\n"
        "#endif // VERSION_HPP\n"
    )
    
    file(WRITE "${panman_SOURCE_DIR}/CMakeLists.txt.orig" "")
    file(READ "${panman_SOURCE_DIR}/CMakeLists.txt" PANMAN_CMAKE_CONTENT)
    file(WRITE "${panman_SOURCE_DIR}/CMakeLists.txt.orig" "${PANMAN_CMAKE_CONTENT}")

    file(WRITE "${panman_SOURCE_DIR}/CMakeLists.txt" [=[
cmake_minimum_required (VERSION 3.8)
project(panmanUtils)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pipe -g -ggdb3 -Wno-unused-function -Wno-deprecated-declarations -w")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Protobuf
find_package(Protobuf REQUIRED)
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()

# Find library paths explicitly to avoid linker search path issues
set(PANMAN_PROTOBUF_LIB ${PROTOBUF_LIBRARY})
set(PANMAN_CAPNP_LIB ${CAPNP_LIBRARY})
set(PANMAN_KJ_LIB ${KJ_LIBRARY})

include_directories(SYSTEM ${Protobuf_INCLUDE_DIRS})

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
            -I${CAPNP_SCHEMA_DIR}
            --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
            ${PANMAN_CAPNP_SRC}
            -I${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS "${PANMAN_CAPNP_SRC}" capnp_tools
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
    ${PANMAN_CAPNP_LIB}
    ${PANMAN_KJ_LIB}
    jsoncpp_lib
    ${PANMAN_PROTOBUF_LIB}
    ${Boost_LIBRARIES}
    ${TBB_LIBRARY}
    pthread
)

# Link dependencies for the executable (including the library)
target_link_libraries(panmanUtils
    panman_lib
    ${PANMAN_CAPNP_LIB}
    ${PANMAN_KJ_LIB}
    jsoncpp_lib
    ${PANMAN_PROTOBUF_LIB}
    ${Boost_LIBRARIES}
    ${TBB_LIBRARY}
    ${PANMAN_CAPNP_LIB}
    ${PANMAN_KJ_LIB}
    pthread
)

# Install targets
install(TARGETS panmanUtils panman_lib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)
]=])

    add_subdirectory(${panman_SOURCE_DIR} ${panman_BINARY_DIR})

    # Ensure the CapnProto tools build first.
    add_dependencies(generate_panman_capnp_extern capnp_tools)

    target_include_directories(panmanUtils PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/include
        ${Boost_INCLUDE_DIRS}
    )

    # Dummy interface library our main project links against.
    add_library(panman INTERFACE)
    target_include_directories(panman INTERFACE
        ${panman_SOURCE_DIR}/src
        ${panman_BINARY_DIR}
    )
    
    message(STATUS "Panman source directory: ${panman_SOURCE_DIR}")
    message(STATUS "Panman binary directory: ${panman_BINARY_DIR}")

add_dependencies(panman panman_lib)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
    set(MINIMAP2_MAKE_ARGS "arm_neon=1" "aarch64=1")
else()
    set(MINIMAP2_MAKE_ARGS "")
endif()
# Minimap2 uses its own Makefile, so it does not inherit CMake's include/lib dirs.
# Under USE_SYSTEM_LIBS (conda) zlib lives only in the env, so point the sub-make's
# compiler at it via CPATH/LIBRARY_PATH (gcc/clang honor these). No-op otherwise.
if(USE_SYSTEM_LIBS AND DEFINED ENV{CONDA_PREFIX})
    set(_mm2_env "CPATH=$ENV{CONDA_PREFIX}/include:$ENV{CPATH}"
                 "LIBRARY_PATH=$ENV{CONDA_PREFIX}/lib:$ENV{LIBRARY_PATH}")
else()
    set(_mm2_env "")
endif()
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND ${CMAKE_COMMAND} -E env ${_mm2_env} make -C ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 clean
    COMMAND ${CMAKE_COMMAND} -E env ${_mm2_env} make -C ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2 -j ${MINIMAP2_MAKE_ARGS}
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/bin
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/libminimap2.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/minimap2 ${CMAKE_CURRENT_BINARY_DIR}/bin/minimap2
    COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2/*.o
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(minimap2_target
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a
    COMMENT "Building minimap2 library from src/3rdparty/minimap2"
)

# Both branches define the imported target `htslib`, a no-op-able `htslib_ext`
# dependency, and HTSLIB_INCLUDE_DIR (consumed below).
if(USE_SYSTEM_LIBS)
    find_library(HTS_LIBRARY NAMES hts REQUIRED)
    find_path(HTS_INCLUDE_DIR htslib/hts.h REQUIRED)
    add_library(htslib SHARED IMPORTED GLOBAL)
    set_target_properties(htslib PROPERTIES
        IMPORTED_LOCATION "${HTS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${HTS_INCLUDE_DIR}"
    )
    add_custom_target(htslib_ext)  # no-op so 'DEPENDS/add_dependencies ... htslib_ext' resolves
    set(HTSLIB_INCLUDE_DIR "${HTS_INCLUDE_DIR}")
else()
    set(HTSLIB_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib_install)
    ExternalProject_Add(htslib_ext
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools/htslib-1.20
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${HTSLIB_INSTALL_DIR} --disable-lzma --disable-bz2 --disable-libcurl --disable-libdeflate
        BUILD_COMMAND make -j lib-static
        INSTALL_COMMAND make install
        BUILD_IN_SOURCE 1 # htslib needs to be built in source
        UPDATE_COMMAND ""
    )
    add_library(htslib STATIC IMPORTED GLOBAL)
    set_target_properties(htslib PROPERTIES
        IMPORTED_LOCATION ${HTSLIB_INSTALL_DIR}/lib/libhts.a
    )
    add_dependencies(htslib htslib_ext)
    set(HTSLIB_INCLUDE_DIR "${HTSLIB_INSTALL_DIR}/include")
endif()


set(BWA_SOURCES
    src/3rdparty/bwa/bamlite.c
    src/3rdparty/bwa/bntseq.c
    src/3rdparty/bwa/bwa.c
    src/3rdparty/bwa/bwamem.c
    src/3rdparty/bwa/bwamem_extra.c
    src/3rdparty/bwa/bwamem_pair.c
    src/3rdparty/bwa/bwape.c
    src/3rdparty/bwa/bwase.c
    src/3rdparty/bwa/bwaseqio.c
    src/3rdparty/bwa/bwashm.c
    src/3rdparty/bwa/bwt.c
    src/3rdparty/bwa/bwt_gen.c
    src/3rdparty/bwa/bwt_lite.c
    src/3rdparty/bwa/bwtaln.c
    src/3rdparty/bwa/bwtgap.c
    src/3rdparty/bwa/bwtindex.c
    src/3rdparty/bwa/bwtsw2_aux.c
    src/3rdparty/bwa/bwtsw2_chain.c
    src/3rdparty/bwa/bwtsw2_core.c
    src/3rdparty/bwa/bwtsw2_main.c
    src/3rdparty/bwa/bwtsw2_pair.c
    src/3rdparty/bwa/fastmap.c
    src/3rdparty/bwa/is.c
    src/3rdparty/bwa/kopen.c
    src/3rdparty/bwa/kstring.c
    src/3rdparty/bwa/ksw.c
    src/3rdparty/bwa/malloc_wrap.c
    src/3rdparty/bwa/maxk.c
    src/3rdparty/bwa/pemerge.c
    src/3rdparty/bwa/QSufSort.c
    src/3rdparty/bwa/rle.c
    src/3rdparty/bwa/rope.c
    src/3rdparty/bwa/run.c
    src/3rdparty/bwa/utils.c
    src/mm_align.c
)

add_custom_target(build-deps DEPENDS
    htslib_ext
    minimap2_target
    capnp_tools
    tbb_build
)

add_library(minimap2 STATIC IMPORTED)
set_target_properties(minimap2 PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/lib/libminimap2.a")
add_dependencies(minimap2 minimap2_target)

add_library(bwa STATIC ${BWA_SOURCES})
target_include_directories(bwa PUBLIC src/3rdparty/bwa)
target_compile_options(bwa PRIVATE -fPIC)

file(GLOB BCFTOOLS_SOURCES src/3rdparty/bcftools/*.c)
add_library(bcftools STATIC ${BCFTOOLS_SOURCES})
target_include_directories(bcftools PUBLIC src/3rdparty/bcftools)
if(USE_SYSTEM_LIBS)
    # build.sh rewrites bcftools' vendored-relative htslib includes to <htslib/...>;
    # point them at the conda htslib headers.
    target_include_directories(bcftools PRIVATE ${HTS_INCLUDE_DIR})
endif()
target_compile_options(bcftools PRIVATE -fPIC)

# Product sources minus main.cpp, compiled once into an OBJECT library (panmap_lib)
# so `panmap` and the test binary link the exact same objects. Single owner of the
# source list; no target re-lists these files. (Distinct from panman_lib, the PanMAN submodule.)
set(PANMAP_LIB_SOURCES
    src/conversion.cpp
    src/panmap_utils.cpp
    src/genotyping.cpp
    src/index_single_mode.cpp
    src/placement.cpp
    src/seeding.cpp
    src/mgsr.cpp
    src/mm_align.c
    src/bwa_align.c
    src/pileup.c
    src/zstd_compression.cpp
    ${LITE_CAPNP_SRCS}
)

find_package(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
find_package(ZLIB REQUIRED)

# Consumers link the target name `libzstd_static`, preserved in both branches.
if(USE_SYSTEM_LIBS)
    find_library(ZSTD_LIBRARY NAMES zstd REQUIRED)
    add_library(libzstd_static INTERFACE IMPORTED)
    set_target_properties(libzstd_static PROPERTIES INTERFACE_LINK_LIBRARIES "${ZSTD_LIBRARY}")
else()
    FetchContent_Declare(
        zstd
        GIT_REPOSITORY https://github.com/facebook/zstd.git
        GIT_TAG v1.5.7
        SOURCE_SUBDIR build/cmake
    )

    set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
    set(ZSTD_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(ZSTD_BUILD_STATIC ON CACHE BOOL "" FORCE)
    set(ZSTD_MULTITHREAD_SUPPORT ON CACHE BOOL "" FORCE)
    set(ZSTD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS_SAVED "${BUILD_SHARED_LIBS}")
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_MakeAvailable(zstd)
    set(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS_SAVED}")
endif()

add_library(panmap_lib OBJECT ${PANMAP_LIB_SOURCES})
add_dependencies(panmap_lib panman_lib minimap2_target build-deps index_capnp_generator tbb_build)

# panmap links panmap_lib (OBJECT library): pulls in its objects AND its PUBLIC usage
# requirements (include dirs, compile options, link libraries), so main.cpp compiles
# in the identical environment. See target_link_libraries below.
add_executable(panmap src/main.cpp)
add_dependencies(panmap panmap_lib build-deps index_capnp_generator tbb_build)
# build/bin matches the docs and the simulate/test targets.
set_target_properties(panmap PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(OPTION_BUILD_SIMULATE)
    add_executable(simulate
        src/test/simulate.cpp
        src/genotyping.cpp
        src/conversion.cpp
        src/mm_align.c
        src/bwa_align.c
        src/pileup.c
    )
    
    add_dependencies(simulate panman_lib minimap2_target htslib_ext)
    
    target_include_directories(simulate PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/src
        ${CMAKE_CURRENT_BINARY_DIR}/include
        ${panman_SOURCE_DIR}/src
        ${panman_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty
        ${HTSLIB_INCLUDE_DIR}
    )
    
    target_link_libraries(simulate PRIVATE
        panman_lib
        minimap2
        bwa
        bcftools
        htslib
        libdeflate_shared
        ${Boost_LIBRARIES}
        ${PROTOBUF_LIBRARIES}
        spdlog::spdlog
        absl::flat_hash_map
        absl::flat_hash_set
        ZLIB::ZLIB
        pthread
        m
    )
    
    set_target_properties(simulate PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
endif()

# Unified test suite (one Boost.Test binary + the e2e shell script).
if(OPTION_BUILD_TESTS)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)
    enable_testing()

    # Shared test-support library: the ONE place tests build/load indices, extract
    # seeds, walk paths, and compute the ground-truth metric oracle. Needs only headers
    # at compile time; symbols resolve when panmap_tests links panmap_lib.
    add_library(testsupport STATIC
        src/test/helpers/metrics_oracle.cpp
        src/test/helpers/seed_helpers.cpp
        src/test/helpers/tree_helpers.cpp
        src/test/helpers/test_index.cpp
        src/test/helpers/traversal.cpp
    )
    add_dependencies(testsupport panmap_lib panman_lib minimap2_target build-deps index_capnp_generator tbb_build)
    # Link panmap_lib PUBLIC so its headers (spdlog, panman) reach testsupport.
    target_link_libraries(testsupport PUBLIC panmap_lib)
    target_compile_definitions(testsupport PUBLIC
        TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/test/data")
    target_include_directories(testsupport PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/test
        ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty
        ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/minimap2
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/include
        ${HTSLIB_INCLUDE_DIR}
        ${panman_SOURCE_DIR}/src
        ${panman_BINARY_DIR}
    )

    # The single test executable. Links panmap_lib (objects + transitive usage
    # requirements) so it names ZERO product source files.
    add_executable(panmap_tests
        src/test/test_main.cpp
        src/test/test_seeding.cpp
        src/test/test_core_utils.cpp
        src/test/test_index.cpp
        src/test/test_placement.cpp
        src/test/test_mgsr.cpp
        src/test/test_genotyping.cpp
        src/test/test_conversion.cpp
    )
    add_dependencies(panmap_tests panmap_lib panman_lib minimap2_target build-deps index_capnp_generator tbb_build)
    target_compile_definitions(panmap_tests PRIVATE BOOST_TEST_DYN_LINK)
    target_link_libraries(panmap_tests PRIVATE
        testsupport
        panmap_lib
        Boost::unit_test_framework
    )
    set_target_properties(panmap_tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

    add_test(NAME unit COMMAND panmap_tests --log_level=test_suite WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# End-to-end tests (always available, run with: ctest -R e2e)
enable_testing()
add_test(NAME e2e
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/src/test/e2e/run_e2e.sh $<TARGET_FILE:panmap>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Include dirs live on panmap_lib as PUBLIC usage requirements so the object library's
# own sources and consumers (panmap, panmap_tests) inherit them.
target_include_directories(panmap_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(panmap_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_include_directories(panmap_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(panmap_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
target_include_directories(panmap_lib PUBLIC ${HTSLIB_INCLUDE_DIR})

# Consumers link the target name `libdeflate_shared`, preserved in both branches.
if(USE_SYSTEM_LIBS)
    find_library(DEFLATE_LIBRARY NAMES deflate REQUIRED)
    add_library(libdeflate_shared INTERFACE IMPORTED)
    set_target_properties(libdeflate_shared PROPERTIES INTERFACE_LINK_LIBRARIES "${DEFLATE_LIBRARY}")
else()
    FetchContent_Declare(deflate
        GIT_REPOSITORY https://github.com/ebiggers/libdeflate
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(deflate)
endif()

# Link libraries live on panmap_lib as PUBLIC so its sources compile with the right
# usage requirements (include dirs / compile defs from panman, htslib, etc.) AND so
# consumers (panmap, panmap_tests) inherit the full link set transitively.
target_link_libraries(panmap_lib PUBLIC
    panman
    panman_lib
    bwa
    bcftools
    minimap2
    htslib
    ${CAPNP_LIBRARY}
    ${KJ_LIBRARY}
    jsoncpp_lib
    ${PROTOBUF_LIB_PATH}
    Boost::program_options
    Boost::iostreams
    Boost::filesystem
    ZLIB::ZLIB
    libzstd_static
    ${TBB_LIBRARY}
    spdlog::spdlog
    libdeflate_shared
    absl::flat_hash_map
    m
    pthread
)

target_link_libraries(panmap PRIVATE panmap_lib)

# Compiler options by build type, on panmap_lib PUBLIC so its sources compile with
# them AND consumers inherit them.
if(OPTION_DEBUG)
    target_compile_options(panmap_lib PUBLIC -g -O0 -ggdb)
elseif(OPTION_PORTABLE)
    target_compile_options(panmap_lib PUBLIC
        -O3
        -march=x86-64-v3
        -mtune=generic
        -fno-omit-frame-pointer
        -ftree-vectorize
        -funroll-loops
        -fprefetch-loop-arrays
    )
else()
    target_compile_options(panmap_lib PUBLIC
        -g -O3
        -march=native
        -fno-omit-frame-pointer
        -ftree-vectorize
        -funroll-loops
        -fprefetch-loop-arrays
    )
endif()

# Test data, only when building tests; lets a tests-off build (e.g. conda) skip
# src/test/data entirely so it can be export-ignored.
if(OPTION_BUILD_TESTS)
    configure_file(src/test/data/test.aligned.fa test.aligned.fa COPYONLY)
    configure_file(src/test/data/test.nwk test.nwk COPYONLY)
    configure_file(src/test/data/test.json test.json COPYONLY)
    configure_file(src/test/data/test.pmat test.pmat COPYONLY)
    configure_file(src/test/data/test.fastq test.fastq COPYONLY)
endif()

install(TARGETS panmap RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_custom_command(
    TARGET panmap POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_SOURCE_DIR}/compile_commands.json
    COMMENT "Creating compile_commands.json symlink for clangd"
)

message(STATUS "")
message(STATUS "Build Configuration Summary:")
message(STATUS "----------------------------")
message(STATUS "CMake Version: ${CMAKE_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "")
message(STATUS "Optimizations:")
message(STATUS "  Debug mode: ${OPTION_DEBUG}")
message(STATUS "")
message(STATUS "Components:")
message(STATUS "  Build tests: ${OPTION_BUILD_TESTS}")
message(STATUS "  Build simulate: ${OPTION_BUILD_SIMULATE}")
message(STATUS "----------------------------")
message(STATUS "")
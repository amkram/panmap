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
include_directories(${CMAKE_SOURCE_DIR}/src)

# Download and configure jsoncpp
FetchContent_Declare(
    jsoncpp
    GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
    GIT_TAG 1.9.5
)

# Configure JsonCpp options before making it available
set(JSONCPP_WITH_TESTS OFF CACHE BOOL "Disable JsonCpp tests")
set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE BOOL "Disable JsonCpp post build tests")
set(JSONCPP_WITH_PKGCONFIG_SUPPORT OFF CACHE BOOL "Disable JsonCpp pkgconfig")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries")

FetchContent_MakeAvailable(jsoncpp)

# Ensure JsonCpp headers are available where panman expects them
file(MAKE_DIRECTORY ${jsoncpp_BINARY_DIR}/include/json)
file(COPY ${jsoncpp_SOURCE_DIR}/include/json/ DESTINATION ${jsoncpp_BINARY_DIR}/include/json/)

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

# Setup for CapnProto - try system first, then fetch if not found
find_program(CAPNP_EXECUTABLE NAMES capnp)
find_program(CAPNPC_CXX_EXECUTABLE NAMES capnpc-c++)
find_package(CapnProto QUIET)

if(CAPNP_EXECUTABLE AND CAPNPC_CXX_EXECUTABLE AND CapnProto_FOUND)
    # Use system CapnProto
    set(CAPNP_LIBRARIES CapnProto::capnp CapnProto::kj)
    set(USING_SYSTEM_CAPNP TRUE)
    message(STATUS "Using system CapnProto tools: ${CAPNP_EXECUTABLE}, ${CAPNPC_CXX_EXECUTABLE}")
    
    # Create a dummy target for consistency
    add_custom_target(capnp_tools ALL)
else()
    message(STATUS "CapnProto not found on system, fetching from source...")
    
    FetchContent_Declare(
        capnproto
        GIT_REPOSITORY https://github.com/capnproto/capnproto.git
        GIT_TAG v1.2.0 # Latest as of Sep 4 2025
    )
    
    set(BUILD_TESTING OFF CACHE BOOL "Disable testing for CapnProto" FORCE)
    FetchContent_MakeAvailable(capnproto)
    
    # Set variables for fetched CapnProto (using generator expressions)
    set(CAPNP_EXECUTABLE $<TARGET_FILE:capnp_tool>)
    set(CAPNPC_CXX_EXECUTABLE $<TARGET_FILE:capnpc_cpp>)
    set(CAPNP_LIBRARIES CapnProto::capnp CapnProto::kj)
    set(USING_SYSTEM_CAPNP FALSE)
    message(STATUS "Using fetched CapnProto")
    
    # Create a target dependency for the tools
    add_custom_target(capnp_tools ALL DEPENDS capnp_tool capnpc_cpp)
endif()

# Set up CapnProto include directories
include_directories(/usr/local/include)

# Generate CapnProto files from index.capnp
set(CAPNP_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/index.capnp")
set(CAPNP_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(CAPNP_GEN_CPP "${CAPNP_GEN_DIR}/index.capnp.c++")
set(CAPNP_GEN_H "${CAPNP_GEN_DIR}/index.capnp.h")

# Custom command to generate CapnProto files
if(USING_SYSTEM_CAPNP)
    # Use system CapnProto
    add_custom_command(
        OUTPUT ${CAPNP_GEN_CPP} ${CAPNP_GEN_H}
        COMMAND ${CAPNP_EXECUTABLE}
            compile
            -o${CAPNPC_CXX_EXECUTABLE}
            --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
            ${CAPNP_SRC_FILE}
        DEPENDS ${CAPNP_SRC_FILE}
        COMMENT "Generating C++ files from index.capnp using system CapnProto"
        VERBATIM
    )
else()
    # Use fetched CapnProto
    add_custom_command(
        OUTPUT ${CAPNP_GEN_CPP} ${CAPNP_GEN_H}
        COMMAND $<TARGET_FILE:capnp_tool>
            compile
            -o$<TARGET_FILE:capnpc_cpp>
            --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
            ${CAPNP_SRC_FILE}
        DEPENDS ${CAPNP_SRC_FILE} capnp_tools
        COMMENT "Generating C++ files from index.capnp using fetched CapnProto"
        VERBATIM
    )
endif()


set(MGSR_CAPNP_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/mgsr_index.capnp")
set(MGSR_CAPNP_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(MGSR_CAPNP_GEN_CPP "${MGSR_CAPNP_GEN_DIR}/mgsr_index.capnp.c++")
set(MGSR_CAPNP_GEN_H "${MGSR_CAPNP_GEN_DIR}/mgsr_index.capnp.h")
# Custom command to generate CapnProto files for mgsr_index.capnp
if(USING_SYSTEM_CAPNP)
    # Use system CapnProto
    add_custom_command(
        OUTPUT ${MGSR_CAPNP_GEN_CPP} ${MGSR_CAPNP_GEN_H}
        COMMAND ${CAPNP_EXECUTABLE}
            compile
            -o${CAPNPC_CXX_EXECUTABLE}
            --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
            ${MGSR_CAPNP_SRC_FILE}
        DEPENDS ${MGSR_CAPNP_SRC_FILE}
        COMMENT "Generating C++ files from mgsr_index.capnp using system CapnProto"
        VERBATIM
    )
else()
    # Use fetched CapnProto
    add_custom_command(
        OUTPUT ${MGSR_CAPNP_GEN_CPP} ${MGSR_CAPNP_GEN_H}
        COMMAND $<TARGET_FILE:capnp_tool>
            compile
            -o$<TARGET_FILE:capnpc_cpp>
            --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
            ${MGSR_CAPNP_SRC_FILE}
        DEPENDS ${MGSR_CAPNP_SRC_FILE} capnp_tools
        COMMENT "Generating C++ files from mgsr_index.capnp using fetched CapnProto"
        VERBATIM
    )
endif()

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
# Create the version.txt file that panman expects
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/version.txt "v1.0.0")

# Set up variables that panman's CMakeLists.txt expects
set(TBB_DIR ${tbb_BINARY_DIR})

# For fetched CapnProto, set up the environment so panman can find the tools
if(NOT USING_SYSTEM_CAPNP)
    # Add the CapnProto binary directory to CMAKE_PROGRAM_PATH so find_program() can locate them
    list(APPEND CMAKE_PROGRAM_PATH "${CMAKE_CURRENT_BINARY_DIR}/bin")
    
    # Make sure CapnProto CMake functions are available to panman
    # Include the CapnProto CMake modules so capnp_generate_cpp function is available
    list(APPEND CMAKE_MODULE_PATH "${capnproto_SOURCE_DIR}/c++/cmake")
    
    # IMPORTANT: Explicitly unset CapnProto_PATH to force panman to use capnp_generate_cpp
    # This ensures panman goes to the else() branch in its CMakeLists.txt
    unset(CapnProto_PATH)
    unset(CapnProto_PATH CACHE)
    
    # Also unset the executable variables for panman so it uses capnp_generate_cpp
    set(PANMAN_CAPNP_EXECUTABLE "")
    set(PANMAN_CAPNPC_CXX_EXECUTABLE "")
endif()

# Set up JsonCpp for panman - ensure the include directory exists
# The issue is that panman expects JsonCpp headers in the build directory
# but they're actually in the source directory. Create them in the build directory.

# Ensure JsonCpp headers are available where panman expects them
# This needs to happen before panman's CMakeLists.txt runs
add_custom_command(
    OUTPUT ${jsoncpp_BINARY_DIR}/include/json/json.h
    COMMAND ${CMAKE_COMMAND} -E make_directory ${jsoncpp_BINARY_DIR}/include/json
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${jsoncpp_SOURCE_DIR}/include/json ${jsoncpp_BINARY_DIR}/include/json
    DEPENDS ${jsoncpp_SOURCE_DIR}/include/json/json.h
    COMMENT "Copying JsonCpp headers to build directory"
)

add_custom_target(jsoncpp_headers ALL DEPENDS ${jsoncpp_BINARY_DIR}/include/json/json.h)

# Create a minimal jsoncppConfig.cmake file that panman can find
file(WRITE ${jsoncpp_BINARY_DIR}/jsoncppConfig.cmake "
# Minimal jsoncpp config for panman
if(NOT TARGET JsonCpp::JsonCpp)
    add_library(JsonCpp::JsonCpp ALIAS jsoncpp_static)
endif()
set(jsoncpp_FOUND TRUE)
set(jsoncpp_VERSION 1.9.5)
")

# Set paths so panman finds our config instead of system JsonCpp
set(jsoncpp_DIR ${jsoncpp_BINARY_DIR})
list(PREPEND CMAKE_PREFIX_PATH ${jsoncpp_BINARY_DIR})

# Ensure CMAKE_TOOLCHAIN_FILE is set to avoid include() errors
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE OR CMAKE_TOOLCHAIN_FILE STREQUAL "")
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/cmake/dummy_toolchain.cmake")
    # Create a dummy toolchain file
    file(WRITE "${CMAKE_TOOLCHAIN_FILE}" "# Dummy toolchain file\n")
endif()

# Create TBB config files that panman expects
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmake)
file(MAKE_DIRECTORY ${tbb_BINARY_DIR}/cmake)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cmake/TBBBuild.cmake "
function(tbb_build)
    cmake_parse_arguments(TBB \"\" \"TBB_ROOT;CONFIG_DIR\" \"MAKE_ARGS\" \${ARGN})
    # Empty implementation - TBB is already built by parent
    set(TBB_IMPORTED_TARGETS \${TBB_LIBRARY} PARENT_SCOPE)
endfunction()
")
# Also create the file in the TBB directory where panman expects it
file(WRITE ${tbb_BINARY_DIR}/cmake/TBBBuild.cmake "
function(tbb_build)
    cmake_parse_arguments(TBB \"\" \"TBB_ROOT;CONFIG_DIR\" \"MAKE_ARGS\" \${ARGN})
    # Empty implementation - TBB is already built by parent
    set(TBB_IMPORTED_TARGETS \${TBB_LIBRARY} PARENT_SCOPE)
endfunction()
")

# Set TBB variables for panman
set(TBB_IMPORTED_TARGETS ${TBB_LIBRARY})

# Create CapnProto imported targets for panman if using system CapnProto
if(CapnProto_FOUND AND CAPNP_EXECUTABLE AND CAPNPC_CXX_EXECUTABLE)
    # Create imported targets for system CapnProto libraries using find_package variables
    if(NOT TARGET CapnProto::capnp)
        add_library(CapnProto::capnp SHARED IMPORTED)
        find_library(CAPNP_LIB NAMES capnp PATHS ${CapnProto_LIBRARY_DIRS} /usr/local/lib /usr/lib)
        set_target_properties(CapnProto::capnp PROPERTIES
            IMPORTED_LOCATION "${CAPNP_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${CapnProto_INCLUDE_DIRS}"
        )
    endif()
    
    if(NOT TARGET CapnProto::kj)
        add_library(CapnProto::kj SHARED IMPORTED)
        find_library(KJ_LIB NAMES kj PATHS ${CapnProto_LIBRARY_DIRS} /usr/local/lib /usr/lib)
        set_target_properties(CapnProto::kj PROPERTIES
            IMPORTED_LOCATION "${KJ_LIB}"
            INTERFACE_INCLUDE_DIRECTORIES "${CapnProto_INCLUDE_DIRS}"
        )
    endif()
    
    set(CAPNP_LIBRARIES CapnProto::capnp CapnProto::kj)
    
    # For system CapnProto, don't set CapnProto_PATH to force panman to use capnp_generate_cpp
    # This ensures compatibility since system CapnProto may not have the modern CMake functions
    unset(CapnProto_PATH)
    unset(CapnProto_PATH CACHE)
else()
    # For fetched CapnProto, the targets are already available from FetchContent
    set(CAPNP_LIBRARIES CapnProto::capnp CapnProto::kj)
endif()

# Fetch the panman source code
FetchContent_Declare(
    panman
    GIT_REPOSITORY https://github.com/TurakhiaLab/panman
    GIT_TAG main
    GIT_SHALLOW TRUE
    UPDATE_COMMAND ""
)

# Make panman available but don't automatically add it to the build
FetchContent_MakeAvailable(panman)

# IMPORTANT: For FetchContent CapnProto, ensure panman uses capnp_generate_cpp
# We need to prevent panman from finding the CapnProto tools so it falls back to capnp_generate_cpp
if(NOT USING_SYSTEM_CAPNP)
    # Temporarily unset these variables so panman doesn't find the tools
    set(SAVED_CAPNP_EXECUTABLE ${CAPNP_EXECUTABLE})
    set(SAVED_CAPNPC_CXX_EXECUTABLE ${CAPNPC_CXX_EXECUTABLE})
    unset(CAPNP_EXECUTABLE)
    unset(CAPNPC_CXX_EXECUTABLE)
    unset(CAPNP_EXECUTABLE CACHE)
    unset(CAPNPC_CXX_EXECUTABLE CACHE)
    message(STATUS "Temporarily hiding CapnProto tools from panman to force capnp_generate_cpp usage")
else()
    # Re-find system CapnProto tools to override any panman-set variables
    find_program(CAPNP_EXECUTABLE NAMES capnp PATHS /usr/local/bin /usr/bin NO_DEFAULT_PATH)
    find_program(CAPNPC_CXX_EXECUTABLE NAMES capnpc-c++ PATHS /usr/local/bin /usr/bin NO_DEFAULT_PATH)
    message(STATUS "Restored system CapnProto tools after panman: ${CAPNP_EXECUTABLE}, ${CAPNPC_CXX_EXECUTABLE}")
endif()

# Copy panman.capnp to the build directory so capnp_generate_cpp can find it
configure_file(
    ${panman_SOURCE_DIR}/panman.capnp
    ${panman_BINARY_DIR}/panman.capnp
    COPYONLY
)

# Also copy the proto files for consistency
configure_file(
    ${panman_SOURCE_DIR}/usher.proto
    ${panman_BINARY_DIR}/usher.proto
    COPYONLY
)
configure_file(
    ${panman_SOURCE_DIR}/panman.proto
    ${panman_BINARY_DIR}/panman.proto
    COPYONLY
)

# Set CapnProto variables that panman's CMakeLists.txt expects
# These variables control how panman generates its CapnProto files
if(NOT USING_SYSTEM_CAPNP)
    # When using fetched CapnProto, we want panman to use capnp_generate_cpp, not capnp_generate
    # So we DON'T set CapnProto_PATH, which forces panman to use the else() branch
    # The tools are already available via cached CAPNP_EXECUTABLE variables
    
    # Make sure panman can find the CapnProto CMake files
    list(APPEND CMAKE_PREFIX_PATH ${capnproto_BINARY_DIR})
    
    message(STATUS "Using fetched CapnProto for panman (via capnp_generate_cpp)")
endif()

# Restore CapnProto tools for our main project after panman configuration
if(NOT USING_SYSTEM_CAPNP)
    # Restore the saved CapnProto tools for our main project's use
    set(CAPNP_EXECUTABLE ${SAVED_CAPNP_EXECUTABLE})
    set(CAPNPC_CXX_EXECUTABLE ${SAVED_CAPNPC_CXX_EXECUTABLE})
    message(STATUS "Restored fetched CapnProto tools for main project: ${CAPNP_EXECUTABLE}, ${CAPNPC_CXX_EXECUTABLE}")
endif()

# Create our own panman library from the panman source files
# Since panman only creates an executable, we need to create a library for linking
add_library(panman_interface STATIC
    ${panman_SOURCE_DIR}/src/panman.cpp
    ${panman_SOURCE_DIR}/src/annotate.cpp
    ${panman_SOURCE_DIR}/src/chaining.cpp
    ${panman_SOURCE_DIR}/src/fasta.cpp
    ${panman_SOURCE_DIR}/src/fitchSankoff.cpp
    ${panman_SOURCE_DIR}/src/gfa.cpp
    ${panman_SOURCE_DIR}/src/maf.cpp
    ${panman_SOURCE_DIR}/src/reroot.cpp
    ${panman_SOURCE_DIR}/src/rotation.cpp
    ${panman_SOURCE_DIR}/src/subnet.cpp
    ${panman_SOURCE_DIR}/src/summary.cpp
    ${panman_SOURCE_DIR}/src/vcf.cpp
    ${panman_SOURCE_DIR}/src/aaTrans.cpp
)

# Set include directories for the panman library
target_include_directories(panman_interface PUBLIC
    ${panman_SOURCE_DIR}/src
    ${panman_BINARY_DIR}
    ${tbb_SOURCE_DIR}/include  # Add TBB include directory
)

# Add TBB preprocessor definitions for panman
target_compile_definitions(panman_interface PRIVATE
    -DTBB_SUPPRESS_DEPRECATED_MESSAGES
    -DTBB_PREVIEW_GLOBAL_CONTROL  # Enable TBB global control features
)

# Add compatibility header to fix missing includes and using declarations
target_compile_options(panman_interface PRIVATE
    -include${CMAKE_CURRENT_SOURCE_DIR}/panman_compat.h
)

# Add forced include of TBB headers for panman source files
target_compile_options(panman_interface PRIVATE
    -include${tbb_SOURCE_DIR}/include/tbb/parallel_for.h
    -include${tbb_SOURCE_DIR}/include/tbb/parallel_for_each.h
    -include${tbb_SOURCE_DIR}/include/tbb/parallel_reduce.h
    -include${tbb_SOURCE_DIR}/include/tbb/parallel_sort.h
    -include${tbb_SOURCE_DIR}/include/tbb/blocked_range.h
    -include${tbb_SOURCE_DIR}/include/tbb/concurrent_vector.h
    -include${tbb_SOURCE_DIR}/include/tbb/concurrent_unordered_set.h
    -includecmath  # Add missing math functions
    -includefilesystem  # Add missing filesystem functions
    -includeutility  # Add std::pair and other utilities
)

# Link panman_interface with its dependencies
target_link_libraries(panman_interface PUBLIC
    jsoncpp_static
    ${CAPNP_LIBRARIES}
    ${Protobuf_LIBRARIES}
    ${Boost_LIBRARIES}
    ${TBB_IMPORTED_TARGETS}
)

# Generate protobuf and capnproto files for panman_interface
protobuf_generate_cpp(PANMAN_PROTO_SRCS PANMAN_PROTO_HDRS 
    ${panman_SOURCE_DIR}/usher.proto 
    ${panman_SOURCE_DIR}/panman.proto
)

# Note: panman generates its own CapnProto files via its CMakeLists.txt
# We set the CapnProto_PATH variable above to ensure it uses our fetched version

# Add generated files to the panman library
target_sources(panman_interface PRIVATE 
    ${PANMAN_PROTO_SRCS} 
    ${PANMAN_PROTO_HDRS}
    # Note: panman.capnp.c++ is handled by panman's own CMakeLists.txt
)

message(STATUS "Panman source directory: ${panman_SOURCE_DIR}")
message(STATUS "Panman binary directory: ${panman_BINARY_DIR}")

# Add dependency on panmanUtils in _deps/panman-build/
add_dependencies(panman_interface panmanUtils)

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
target_include_directories(panmap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools)
add_dependencies(panmap build-deps index_capnp_generator tbb_build)

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
    panman_interface
    bwa
    bcftools
    minimap2 
    samtools 
    htslib 
    docopt
    # Use consistent CapnProto libraries (either system or fetched)
    ${CAPNP_LIBRARIES}
    jsoncpp_static
    ${Protobuf_LIBRARIES}
    Boost::program_options
    Boost::iostreams
    Boost::filesystem
    Boost::system
    ZLIB::ZLIB
    ${TBB_LIBRARY}
    spdlog::spdlog  # Add spdlog
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
    target_compile_options(panmap PUBLIC -g -O2 -pg -fno-omit-frame-pointer)
    target_link_options(panmap PUBLIC -pg -fno-omit-frame-pointer)
    
    # Also enable profiling for panmanUtils
    if(TARGET panmanUtils)
        target_compile_options(panmanUtils PUBLIC -g -O2 -pg -fno-omit-frame-pointer)
        target_link_options(panmanUtils PUBLIC -pg -fno-omit-frame-pointer)
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

# Build simulate if requested
if(OPTION_BUILD_SIMULATE)
    set(SIMULATE_SOURCES
        src/test/simulate.cpp
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
    
    add_executable(simulate ${SIMULATE_SOURCES})
    target_include_directories(simulate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty/samtools)
    add_dependencies(simulate build-deps index_capnp_generator tbb_build)
    
    # Add generated sources to simulate target
    target_sources(simulate PRIVATE ${CAPNP_SRCS} ${MGSR_CAPNP_SRCS})
    
    # Add include directories
    target_include_directories(simulate PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_include_directories(simulate PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
    target_include_directories(simulate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/3rdparty)
    target_include_directories(simulate PUBLIC ${HTSLIB_INSTALL_DIR}/include)
    
    # Link libraries for simulate
    target_link_libraries(simulate PRIVATE 
    panman_interface
    bwa
    bcftools
    minimap2 
    samtools 
    htslib 
    docopt
    # Use consistent CapnProto libraries (either system or fetched)
    ${CAPNP_LIBRARIES}
    jsoncpp_static
    ${Protobuf_LIBRARIES}
    Boost::program_options
    Boost::iostreams
    Boost::filesystem
    Boost::system
    ZLIB::ZLIB
    ${TBB_LIBRARY}
    spdlog::spdlog
    libdeflate_shared
    absl::flat_hash_map
    m
    pthread
    )
    
    # Set compiler options for simulate
    if(OPTION_DEBUG)
        target_compile_options(simulate PUBLIC -g -O0 -ggdb)
    elseif(OPTION_PROFILE)
        target_compile_options(simulate PUBLIC -g -O2 -pg -fno-omit-frame-pointer)
        target_link_options(simulate PUBLIC -pg -fno-omit-frame-pointer)
    else()
        target_compile_options(simulate PUBLIC -g -O3)
        if(OPTION_ENABLE_LTO)
            target_compile_options(simulate PUBLIC -flto)
            target_link_options(simulate PUBLIC -flto)
        endif()
        target_compile_options(simulate PUBLIC -ftree-vectorize)
    endif()
    
    install(TARGETS simulate RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

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

#
# Find the Cap'n Proto libraries and tools for use with CMake
#
# This module defines the following variables:
#  CapnProto_FOUND           - True if Cap'n Proto is found
#  CAPNP_INCLUDE_DIRS        - Include directories for Cap'n Proto headers
#  CAPNP_LIBRARIES           - Link to libcapnp.a
#  KJ_LIBRARIES              - Link to libkj.a
#  CAPNP_EXECUTABLE          - Path to the capnp tool
#  CAPNPC_CXX_EXECUTABLE     - Path to the capnpc-c++ tool
#
# The module also defines the following function:
#  CAPNP_GENERATE_CPP - Generate C++ code from Cap'n Proto schema files
#

# Check if variables are already defined by the parent CMake script
if(CAPNP_INCLUDE_DIR AND CAPNP_LIBRARY AND KJ_LIBRARY AND CAPNP_EXECUTABLE AND CAPNPC_CXX_EXECUTABLE)
  set(CapnProto_FOUND TRUE)
  set(CAPNP_INCLUDE_DIRS "${CAPNP_INCLUDE_DIR}")
  message(STATUS "Using pre-defined Cap'n Proto paths: ${CAPNP_INCLUDE_DIRS}")
  
  # Skip the rest of the file when we've found Cap'n Proto explicitly
  return()
endif()

# We want to check if we have a built version of Cap'n Proto from the build directory
# first, because that's the most reliable.
if(EXISTS "${CMAKE_BINARY_DIR}/include/capnp/common.h")
  set(CapnProto_FOUND TRUE)
  set(CAPNP_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/include")
  set(CAPNP_LIBRARIES "${CMAKE_BINARY_DIR}/lib/libcapnp.a")
  set(KJ_LIBRARIES "${CMAKE_BINARY_DIR}/lib/libkj.a")
  set(CAPNP_EXECUTABLE "${CMAKE_BINARY_DIR}/bin/capnp")
  set(CAPNPC_CXX_EXECUTABLE "${CMAKE_BINARY_DIR}/bin/capnpc-c++")
  message(STATUS "Using Cap'n Proto from local build: ${CAPNP_INCLUDE_DIRS}")
  
  # Skip the rest of the file when we've found Cap'n Proto locally
  return()
endif()

# Find the include directory
find_path(CAPNP_INCLUDE_DIR capnp/common.h
  PATHS "${CAPNP_ROOT_DIR}/include" "$ENV{CAPNP_ROOT_DIR}/include" "${CAPNP_INCLUDE_DIRS}"
  HINTS
  ${CAPNP_INCLUDE_DIRS}
  PATH_SUFFIXES include
)

# Find the Cap'n Proto libraries
find_library(CAPNP_LIBRARY NAMES capnp
  PATHS "${CAPNP_ROOT_DIR}/lib" "$ENV{CAPNP_ROOT_DIR}/lib" "${CMAKE_CURRENT_BINARY_DIR}/lib"
  HINTS
  ${CAPNP_LIBRARIES}
  PATH_SUFFIXES lib
)

find_library(KJ_LIBRARY NAMES kj
  PATHS "${CAPNP_ROOT_DIR}/lib" "$ENV{CAPNP_ROOT_DIR}/lib" "${CMAKE_CURRENT_BINARY_DIR}/lib"
  HINTS
  ${CAPNP_LIBRARIES}
  PATH_SUFFIXES lib
)

# Find capnp executable
find_program(CAPNP_EXECUTABLE capnp
  PATHS "${CAPNP_ROOT_DIR}/bin" "$ENV{CAPNP_ROOT_DIR}/bin" "${CMAKE_CURRENT_BINARY_DIR}/bin"
  HINTS
  ${CAPNP_EXECUTABLE}
  PATH_SUFFIXES bin
)

# Find capnpc-c++ executable
find_program(CAPNPC_CXX_EXECUTABLE capnpc-c++
  PATHS "${CAPNP_ROOT_DIR}/bin" "$ENV{CAPNP_ROOT_DIR}/bin" "${CMAKE_CURRENT_BINARY_DIR}/bin"
  HINTS
  ${CAPNPC_CXX_EXECUTABLE}
  PATH_SUFFIXES bin
)

# Version detection
if(CAPNP_EXECUTABLE)
  execute_process(
    COMMAND ${CAPNP_EXECUTABLE} --version
    OUTPUT_VARIABLE CAPNP_VERSION_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" CAPNP_VERSION "${CAPNP_VERSION_STRING}")
endif()

# Handle standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CapnProto 
  REQUIRED_VARS CAPNP_INCLUDE_DIR CAPNP_LIBRARY KJ_LIBRARY CAPNP_EXECUTABLE CAPNPC_CXX_EXECUTABLE
  VERSION_VAR CAPNP_VERSION
)

# Set the variables
if(CapnProto_FOUND)
  set(CAPNP_INCLUDE_DIRS ${CAPNP_INCLUDE_DIR})
  set(CAPNP_LIBRARIES ${CAPNP_LIBRARY} ${KJ_LIBRARY})
  
  # Backward compatibility variables
  set(CAPNP_LIBRARY ${CAPNP_LIBRARY})
  set(KJ_LIBRARY ${KJ_LIBRARY})
  
  # Add definitions
  set(CAPNP_DEFINITIONS)
endif()

# Function to generate C++ code from Cap'n Proto schema files
function(CAPNP_GENERATE_CPP SOURCES HEADERS)
  if(NOT ARGN)
    message(SEND_ERROR "CAPNP_GENERATE_CPP() called without any schema files")
    return()
  endif()

  # Parse optional arguments
  set(options)
  set(oneValueArgs DIRECTORY OUTPUT_DIR)
  set(multiValueArgs IMPORT_DIRS)
  cmake_parse_arguments(CAPNP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  # Default output directory is the current build directory
  if(NOT CAPNP_OUTPUT_DIR)
    set(CAPNP_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
  
  # Create the output directory if it doesn't exist
  if(NOT IS_DIRECTORY "${CAPNP_OUTPUT_DIR}")
    file(MAKE_DIRECTORY "${CAPNP_OUTPUT_DIR}")
  endif()

  # Set up the include directories for capnp command
  set(CAPNP_COMMAND_ARGS)
  if(CAPNP_IMPORT_DIRS)
    foreach(DIR ${CAPNP_IMPORT_DIRS})
      list(APPEND CAPNP_COMMAND_ARGS "-I${DIR}")
    endforeach()
  endif()

  # Default to the schema file directory
  if(CAPNP_DIRECTORY)
    list(APPEND CAPNP_COMMAND_ARGS "-I${CAPNP_DIRECTORY}")
  endif()

  # Process each schema file
  set(${SOURCES})
  set(${HEADERS})
  foreach(SCHEMA_FILE ${CAPNP_UNPARSED_ARGUMENTS})
    # Get the file basename
    get_filename_component(FILE_NAME "${SCHEMA_FILE}" NAME_WE)
    
    # Set the output file paths
    set(SCHEMA_CPP "${CAPNP_OUTPUT_DIR}/${FILE_NAME}.c++")
    set(SCHEMA_H "${CAPNP_OUTPUT_DIR}/${FILE_NAME}.h")
    
    # Add to the list of generated files
    list(APPEND ${SOURCES} "${SCHEMA_CPP}")
    list(APPEND ${HEADERS} "${SCHEMA_H}")

    # Create the command to generate the files
    add_custom_command(
      OUTPUT "${SCHEMA_CPP}" "${SCHEMA_H}"
      COMMAND "${CAPNP_EXECUTABLE}"
      ARGS compile -o "${CAPNPC_CXX_EXECUTABLE}:${CAPNP_OUTPUT_DIR}" ${CAPNP_COMMAND_ARGS} "${SCHEMA_FILE}"
      DEPENDS "${SCHEMA_FILE}"
      COMMENT "Generating C++ from ${SCHEMA_FILE}"
      VERBATIM
    )
  endforeach()

  # Set parent scope variables
  set(${SOURCES} ${${SOURCES}} PARENT_SCOPE)
  set(${HEADERS} ${${HEADERS}} PARENT_SCOPE)
endfunction()

# Define find_package compatibility mode
set(CapnProto_FOUND ${CAPNP_FOUND})
set(CapnProto_VERSION "1.0.0") # Placeholder version
set(CAPNP_FOUND TRUE)

# Ensure we have the required Cap'n Proto components
if(NOT DEFINED CAPNP_EXECUTABLE)
  message(FATAL_ERROR "CAPNP_EXECUTABLE must be defined")
endif()

if(NOT DEFINED CAPNPC_CXX_EXECUTABLE)
  message(FATAL_ERROR "CAPNPC_CXX_EXECUTABLE must be defined")
endif()

if(NOT DEFINED CAPNP_INCLUDE_DIRS)
  set(CAPNP_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}/include")
  message(STATUS "Setting default CAPNP_INCLUDE_DIRS to ${CAPNP_INCLUDE_DIRS}")
endif()

if(NOT DEFINED CAPNP_LIBRARIES)
  set(CAPNP_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/lib/libcapnp.a" "${CMAKE_CURRENT_BINARY_DIR}/lib/libkj.a")
  message(STATUS "Setting default CAPNP_LIBRARIES to ${CAPNP_LIBRARIES}")
endif()

# Define the capnp_generate_cpp function
function(capnp_generate_cpp SRCS HDRS)
  if(NOT ARGN)
    message(SEND_ERROR "capnp_generate_cpp: missing capnp input files")
    return()
  endif()
  
  set(${SRCS})
  set(${HDRS})

  foreach(CAPNP_FILE ${ARGN})
    get_filename_component(CAPNP_FILE_ABS "${CAPNP_FILE}" ABSOLUTE)
    get_filename_component(CAPNP_NAME "${CAPNP_FILE}" NAME_WE)
    
    # Set output file paths
    set(CAPNP_OUTPUT_CPP "${CMAKE_CURRENT_BINARY_DIR}/${CAPNP_NAME}.capnp.c++")
    set(CAPNP_OUTPUT_H "${CMAKE_CURRENT_BINARY_DIR}/${CAPNP_NAME}.capnp.h")
    
    # Add custom command to generate the Cap'n Proto files
    add_custom_command(
      OUTPUT "${CAPNP_OUTPUT_H}" "${CAPNP_OUTPUT_CPP}"
      COMMAND "${CAPNP_EXECUTABLE}"
      ARGS compile -o "${CAPNPC_CXX_EXECUTABLE}" 
           --src-prefix "${CMAKE_CURRENT_SOURCE_DIR}" 
           "${CAPNP_FILE_ABS}"
      DEPENDS "${CAPNP_FILE_ABS}"
      COMMENT "Compiling Cap'n Proto file ${CAPNP_FILE}"
      VERBATIM
    )
    
    list(APPEND ${SRCS} "${CAPNP_OUTPUT_CPP}")
    list(APPEND ${HDRS} "${CAPNP_OUTPUT_H}")
  endforeach()

  set(${SRCS} ${${SRCS}} PARENT_SCOPE)
  set(${HDRS} ${${HDRS}} PARENT_SCOPE)
endfunction()

# Define the newer-style capnp_generate function for compatibility
if(NOT COMMAND capnp_generate)
  function(capnp_generate)
    # Parse arguments
    set(options "")
    set(oneValueArgs LANGUAGE TARGET)
    set(multiValueArgs PROTOS)
    cmake_parse_arguments(CAPNP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Validate arguments
    if(NOT CAPNP_LANGUAGE)
      message(FATAL_ERROR "capnp_generate requires LANGUAGE argument")
    endif()
    if(NOT CAPNP_PROTOS)
      message(FATAL_ERROR "capnp_generate requires PROTOS argument")
    endif()
    
    # Use the core capnp_generate_cpp function for C++ generation
    if(CAPNP_LANGUAGE STREQUAL "cpp")
      set(SRCS)
      set(HDRS)
      capnp_generate_cpp(SRCS HDRS ${CAPNP_PROTOS})
      
      # Add generated files to target
      if(CAPNP_TARGET)
        target_sources(${CAPNP_TARGET} PRIVATE ${SRCS} ${HDRS})
      endif()
    else()
      message(WARNING "Unsupported language ${CAPNP_LANGUAGE}")
    endif()
  endfunction()
endif()

mark_as_advanced(
  CAPNP_INCLUDE_DIR
  CAPNP_LIBRARY
  KJ_LIBRARY
  CAPNP_EXECUTABLE
  CAPNPC_CXX_EXECUTABLE
)

# Include macro for generating C++ sources from Cap'n Proto files
include(${CMAKE_CURRENT_LIST_DIR}/CapnProtoMacros.cmake) 
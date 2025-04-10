# FindProtobuf.cmake
# Find and configure Protocol Buffers for use in this project
# Sets the following variables:
#  Protobuf_FOUND - whether Protocol Buffers was found
#  Protobuf_INCLUDE_DIRS - the Protocol Buffers include directories
#  Protobuf_LIBRARIES - the Protocol Buffers libraries
#  Protobuf_PROTOC_EXECUTABLE - the Protocol Buffers compiler

# Define compatibility variables if they don't exist
if(NOT DEFINED Protobuf_FOUND)
  set(Protobuf_FOUND TRUE)
endif()

if(NOT DEFINED Protobuf_VERSION AND NOT DEFINED Protobuf_VERSION_STRING)
  set(Protobuf_VERSION "3.12.4")
  set(Protobuf_VERSION_STRING "3.12.4")
endif()

if(NOT DEFINED Protobuf_PROTOC_EXECUTABLE)
  find_program(Protobuf_PROTOC_EXECUTABLE protoc)
  if(NOT Protobuf_PROTOC_EXECUTABLE)
    message(WARNING "protoc executable not found, Protocol Buffer generation may fail")
  endif()
endif()

if(NOT DEFINED Protobuf_INCLUDE_DIRS)
  set(Protobuf_INCLUDE_DIRS "/usr/include")
endif()

if(NOT DEFINED Protobuf_LIBRARIES)
  set(Protobuf_LIBRARIES "/usr/lib/libprotobuf.so")
endif()

# Define the protobuf_generate_cpp function
function(protobuf_generate_cpp SRCS HDRS)
  set(options "")
  set(oneValueArgs EXPORT_MACRO PROTOC_OUT_DIR)
  set(multiValueArgs PROTOS)
  cmake_parse_arguments(protobuf_generate_cpp "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  if(NOT protobuf_generate_cpp_PROTOS)
    message(SEND_ERROR "protobuf_generate_cpp: missing PROTO input files")
    return()
  endif()
  
  set(outdir "${CMAKE_CURRENT_BINARY_DIR}")
  if(protobuf_generate_cpp_PROTOC_OUT_DIR)
    set(outdir "${protobuf_generate_cpp_PROTOC_OUT_DIR}")
  endif()
  
  set(_generated_srcs)
  set(_generated_hdrs)
  
  foreach(_proto ${protobuf_generate_cpp_PROTOS})
    get_filename_component(_abs_file "${_proto}" ABSOLUTE)
    get_filename_component(_basename "${_proto}" NAME_WE)
    
    set(_generated_src "${outdir}/${_basename}.pb.cc")
    set(_generated_hdr "${outdir}/${_basename}.pb.h")
    
    # Add the custom command to generate the Protocol Buffer files
    add_custom_command(
      OUTPUT "${_generated_src}" "${_generated_hdr}"
      COMMAND ${Protobuf_PROTOC_EXECUTABLE}
      ARGS --cpp_out=${outdir} --proto_path=${CMAKE_CURRENT_SOURCE_DIR} ${_abs_file}
      DEPENDS ${_abs_file}
      COMMENT "Running C++ protocol buffer compiler on ${_proto}"
      VERBATIM
    )
    
    list(APPEND _generated_srcs "${_generated_src}")
    list(APPEND _generated_hdrs "${_generated_hdr}")
  endforeach()
  
  set(${SRCS} ${_generated_srcs} PARENT_SCOPE)
  set(${HDRS} ${_generated_hdrs} PARENT_SCOPE)
endfunction()

# Define the newer-style protobuf_generate function for compatibility
if(NOT COMMAND protobuf_generate)
  function(protobuf_generate)
    # Parse arguments
    set(options "")
    set(oneValueArgs LANGUAGE TARGET)
    set(multiValueArgs PROTOS)
    cmake_parse_arguments(PROTOBUF "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Validate arguments
    if(NOT PROTOBUF_LANGUAGE)
      message(FATAL_ERROR "protobuf_generate requires LANGUAGE argument")
    endif()
    if(NOT PROTOBUF_PROTOS)
      message(FATAL_ERROR "protobuf_generate requires PROTOS argument")
    endif()
    
    # Use the core protobuf_generate_cpp function for C++ generation
    if(PROTOBUF_LANGUAGE STREQUAL "cpp")
      set(SRCS)
      set(HDRS)
      protobuf_generate_cpp(SRCS HDRS PROTOS ${PROTOBUF_PROTOS})
      
      # Add generated files to target
      if(PROTOBUF_TARGET)
        target_sources(${PROTOBUF_TARGET} PRIVATE ${SRCS} ${HDRS})
      endif()
    else()
      message(WARNING "Unsupported language ${PROTOBUF_LANGUAGE}")
    endif()
  endfunction()
endif()

# Custom FindProtobuf module for panmap

# Try to use the system module first
include(FindPackageHandleStandardArgs)

# Find the Protobuf compiler
find_program(PROTOBUF_PROTOC_EXECUTABLE
    NAMES protoc
    DOC "The Google Protocol Buffers Compiler"
)

# Find the Protobuf include directory
find_path(PROTOBUF_INCLUDE_DIR
    NAMES google/protobuf/service.h
    DOC "The Google Protocol Buffers Headers"
)

# Find the Protobuf libraries
find_library(PROTOBUF_LIBRARY
    NAMES protobuf
    DOC "The Google Protocol Buffers Library"
)

find_library(PROTOBUF_LITE_LIBRARY
    NAMES protobuf-lite
    DOC "The Google Protocol Buffers Lite Library"
)

find_library(PROTOBUF_PROTOC_LIBRARY
    NAMES protoc
    DOC "The Google Protocol Buffers Compiler Library"
)

# Set PROTOBUF_FOUND based on the availability of the components
mark_as_advanced(PROTOBUF_INCLUDE_DIR PROTOBUF_LIBRARY PROTOBUF_LITE_LIBRARY PROTOBUF_PROTOC_LIBRARY PROTOBUF_PROTOC_EXECUTABLE)

# Improved protobuf_generate_cpp function with better error handling
function(protobuf_generate_cpp SRCS HDRS)
    # Parse arguments
    set(_options)
    set(_singleargs EXPORT_MACRO)
    set(_multiargs PROTOS PROTO_PATH)
    cmake_parse_arguments(protobuf_generate_cpp "${_options}" "${_singleargs}" "${_multiargs}" ${ARGN})

    # Check if PROTOS is specified
    if(NOT protobuf_generate_cpp_PROTOS)
        # Try to use positional arguments instead
        set(protobuf_generate_cpp_PROTOS ${ARGN})
    endif()

    # Final check for PROTOS files
    if(NOT protobuf_generate_cpp_PROTOS)
        message(FATAL_ERROR "protobuf_generate_cpp: missing PROTO input files. Please specify .proto files via PROTOS parameter or as direct arguments.")
        return()
    endif()

    # Debug information about proto files
    message(STATUS "Protocol Buffer files to compile: ${protobuf_generate_cpp_PROTOS}")

    # Check if all proto files exist and are accessible
    foreach(proto ${protobuf_generate_cpp_PROTOS})
        # Try different directories to find the proto file
        set(proto_file "${proto}")
        if(NOT EXISTS "${proto_file}")
            set(proto_file "${CMAKE_CURRENT_SOURCE_DIR}/${proto}")
        endif()
        if(NOT EXISTS "${proto_file}")
            set(proto_file "${CMAKE_CURRENT_BINARY_DIR}/${proto}")
        endif()
        
        # If still not found, look in common protocol buffer directories
        if(NOT EXISTS "${proto_file}")
            message(WARNING "Protocol Buffer file '${proto}' not found at expected locations. Checking in additional directories...")
            file(GLOB_RECURSE proto_search_result 
                 "${CMAKE_SOURCE_DIR}/*/${proto}"
                 "${CMAKE_BINARY_DIR}/*/${proto}")
            if(proto_search_result)
                # Use the first found file
                list(GET proto_search_result 0 proto_file)
                message(STATUS "Found Protocol Buffer file at: ${proto_file}")
            else()
                message(FATAL_ERROR "Protocol Buffer file '${proto}' not found. Please make sure it exists and is accessible.")
                return()
            endif()
        endif()
        
        # Generate the C++ code from the proto file
        get_filename_component(proto_name "${proto}" NAME_WE)
        set(proto_src "${CMAKE_CURRENT_BINARY_DIR}/${proto_name}.pb.cc")
        set(proto_hdr "${CMAKE_CURRENT_BINARY_DIR}/${proto_name}.pb.h")
        
        get_filename_component(proto_abs_path "${proto_file}" ABSOLUTE)
        get_filename_component(proto_dir "${proto_abs_path}" DIRECTORY)
                
        add_custom_command(
            OUTPUT "${proto_src}" "${proto_hdr}"
            COMMAND ${PROTOBUF_PROTOC_EXECUTABLE}
            ARGS --cpp_out=${CMAKE_CURRENT_BINARY_DIR} --proto_path=${proto_dir} ${proto_abs_path}
            DEPENDS ${proto_abs_path} ${PROTOBUF_PROTOC_EXECUTABLE}
            COMMENT "Running C++ protocol buffer compiler on ${proto}"
            VERBATIM
        )
        
        list(APPEND ${SRCS} "${proto_src}")
        list(APPEND ${HDRS} "${proto_hdr}")
    endforeach()
    
    # Set output variables in parent scope
    set(${SRCS} ${${SRCS}} PARENT_SCOPE)
    set(${HDRS} ${${HDRS}} PARENT_SCOPE)
endfunction() 
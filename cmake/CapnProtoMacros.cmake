# CapnProto CMake Macros
# This function generates C++ source files from Cap'n Proto schema files.
function(capnp_generate_cpp SOURCES HEADERS)
  if(NOT ARGN)
    message(SEND_ERROR "capnp_generate_cpp: missing schema files")
    return()
  endif()

  if(NOT CAPNP_EXECUTABLE)
    message(SEND_ERROR "Could not locate Cap'n Proto executable (CAPNP_EXECUTABLE).")
    return()
  endif()
  
  if(NOT CAPNPC_CXX_EXECUTABLE)
    message(SEND_ERROR "Could not locate Cap'n Proto C++ compiler (CAPNPC_CXX_EXECUTABLE).")
    return()
  endif()

  # Process each schema file
  set( "" PARENT_SCOPE)
  set( "" PARENT_SCOPE)
  
  foreach(schema_file )
    # Get the file name without extension
    get_filename_component(file_name ${schema_file} NAME_WE)
    
    # Set output files
    set(schema_cpp "${CMAKE_CURRENT_BINARY_DIR}/${file_name}.capnp.c++")
    set(schema_h "${CMAKE_CURRENT_BINARY_DIR}/${file_name}.capnp.h")
    
    # Add to output lists
    list(APPEND output_sources "${schema_cpp}")
    list(APPEND output_headers "${schema_h}")

    # Add command to generate the source
    add_custom_command(
      OUTPUT "${schema_cpp}" "${schema_h}"
      COMMAND ${CAPNP_EXECUTABLE}
          compile 
          -o${CAPNPC_CXX_EXECUTABLE}
          --src-prefix=${CMAKE_CURRENT_SOURCE_DIR}
          ${schema_file}
      DEPENDS "${schema_file}"
      COMMENT "Compiling Cap'n Proto schema ${schema_file}"
      VERBATIM
    )
  endforeach()
  
  # Set output variables in the parent scope
  set( ${output_sources} PARENT_SCOPE)
  set( ${output_headers} PARENT_SCOPE)
endfunction()

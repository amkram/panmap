# Architecture detection script for CMake
# This file provides functions to detect CPU architecture

include(CheckCXXCompilerFlag)

# Detect the host architecture
function(detect_host_arch OUT_VAR)
    if(APPLE)
        # Check for Apple ARM (M1/M2)
        execute_process(
            COMMAND uname -m
            OUTPUT_VARIABLE ARCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        if(ARCH MATCHES "arm64")
            set(HAS_NEON TRUE PARENT_SCOPE)
            set(${OUT_VAR} "arm64" PARENT_SCOPE)
        elseif(ARCH MATCHES "x86_64")
            set(HAS_NEON FALSE PARENT_SCOPE)
            set(${OUT_VAR} "x86_64" PARENT_SCOPE)
        else()
            set(HAS_NEON FALSE PARENT_SCOPE)
            set(${OUT_VAR} ${ARCH} PARENT_SCOPE)
        endif()
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
        set(HAS_NEON TRUE PARENT_SCOPE)
        set(${OUT_VAR} "arm64" PARENT_SCOPE)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
        set(HAS_NEON TRUE PARENT_SCOPE)
        set(${OUT_VAR} "arm" PARENT_SCOPE)
    else()
        set(HAS_NEON FALSE PARENT_SCOPE)
        set(${OUT_VAR} ${CMAKE_SYSTEM_PROCESSOR} PARENT_SCOPE)
    endif()
endfunction() 
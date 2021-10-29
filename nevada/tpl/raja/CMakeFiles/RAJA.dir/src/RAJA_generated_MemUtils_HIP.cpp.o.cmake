# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other BLT Project Developers. See the top-level LICENSE file for details
#
# SPDX-License-Identifier: (BSD-3-Clause)

###############################################################################
# Runs commands using HIPCC
###############################################################################

###############################################################################
# This file runs the hipcc commands to produce the desired output file
# along with the dependency file needed by CMake to compute dependencies.
#
# Input variables:
#
# verbose:BOOL=<>               OFF: Be as quiet as possible (default)
#                               ON : Describe each step
# build_configuration:STRING=<> Build configuration. Defaults to Debug.
# generated_file:STRING=<>      File to generate. Mandatory argument.

if(NOT build_configuration)
    set(build_configuration Debug)
endif()
if(NOT generated_file)
    message(FATAL_ERROR "You must specify generated_file on the command line")
endif()

# Set these up as variables to make reading the generated file easier
set(HIP_HIPCC_EXECUTABLE "/opt/rocm-4.2.0/hip/bin/hipcc") # path
set(HIP_HIPCONFIG_EXECUTABLE "/opt/rocm-4.2.0/hip/bin/hipconfig") #path
set(HIP_HOST_COMPILER "/opt/rocm-4.2.0/llvm/bin/clang++") # path
set(CMAKE_COMMAND "/usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.2.1/cmake-3.14.5-2o4l6ukye2v2vfcj3hrazdptk4vrejsn/bin/cmake") # path
set(HIP_run_make2cmake "/usr/workspace/rana2/Kripke/blt/cmake/thirdparty/FindHIP/run_make2cmake.cmake") # path
set(HCC_HOME "") #path
set(HIP_CLANG_PATH "/opt/rocm-4.2.0/llvm/bin") #path
set(HIP_CLANG_PARALLEL_BUILD_COMPILE_OPTIONS "")


set(HIP_HIPCC_FLAGS -D__HIP_PLATFORM_HCC__ -D__HIP_ROCclr__ -D__HIP_PLATFORM_AMD__ -D__HIP_PLATFORM_HCC__ -D__HIP_ROCclr__ -D__HIP_PLATFORM_AMD__ -D__HIP_PLATFORM_HCC__ -D__HIP_ROCclr__ -D__HIP_PLATFORM_AMD__ -D__HIP_ARCH_GFX900__=1 -D__HIP_PLATFORM_HCC__ -D__HIP_ROCclr__ -D__HIP_PLATFORM_AMD__ -D__HIP_ARCH_GFX900__=1)
set(HIP_HIPCC_FLAGS_RELEASE -O2)
set(HIP_HIPCC_FLAGS_DEBUG -g;-O0)
set(HIP_HIPCC_FLAGS_MINSIZEREL -Os)
set(HIP_HIPCC_FLAGS_RELWITHDEBINFO -g;-O2)
set(HIP_HCC_FLAGS )
set(HIP_HCC_FLAGS_RELEASE )
set(HIP_HCC_FLAGS_DEBUG )
set(HIP_HCC_FLAGS_MINSIZEREL )
set(HIP_HCC_FLAGS_RELWITHDEBINFO )
set(HIP_CLANG_FLAGS )
set(HIP_CLANG_FLAGS_RELEASE )
set(HIP_CLANG_FLAGS_DEBUG )
set(HIP_CLANG_FLAGS_MINSIZEREL )
set(HIP_CLANG_FLAGS_RELWITHDEBINFO )
set(HIP_NVCC_FLAGS )
set(HIP_NVCC_FLAGS_RELEASE )
set(HIP_NVCC_FLAGS_DEBUG )
set(HIP_NVCC_FLAGS_MINSIZEREL )
set(HIP_NVCC_FLAGS_RELWITHDEBINFO )
#Needed to bring the HIP_HIPCC_INCLUDE_ARGS variable in scope
set(HIP_HIPCC_INCLUDE_ARGS -I/usr/workspace/rana2/Kripke/tpl/raja/include -I/usr/workspace/rana2/Kripke/nevada/tpl/raja/include -I/opt/rocm-4.2.0/hip/include -I/opt/rocm-4.2.0/hip/include -I/opt/rocm-4.2.0/include -I/opt/rocm-4.2.0/include -I/usr/workspace/rana2/Kripke/tpl/raja/tpl/camp/include) # list

set(cmake_dependency_file "/usr/workspace/rana2/Kripke/nevada/tpl/raja/CMakeFiles/RAJA.dir/src/RAJA_generated_MemUtils_HIP.cpp.o.depend") # path
set(source_file "/usr/workspace/rana2/Kripke/tpl/raja/src/MemUtils_HIP.cpp") # path
set(host_flag "FALSE") # bool

# Determine compiler and compiler flags
execute_process(COMMAND ${HIP_HIPCONFIG_EXECUTABLE} --platform OUTPUT_VARIABLE HIP_PLATFORM OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${HIP_HIPCONFIG_EXECUTABLE} --compiler OUTPUT_VARIABLE HIP_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${HIP_HIPCONFIG_EXECUTABLE} --runtime OUTPUT_VARIABLE HIP_RUNTIME OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT host_flag)
    set(__CC ${HIP_HIPCC_EXECUTABLE})
    # ROCm version 4.0 and higher will use amd for the platform, but older
    # versions use hcc. Keep both for backwards compatibility.
    if("${HIP_PLATFORM}" STREQUAL "hcc" OR "${HIP_PLATFORM}" STREQUAL "amd")
        if("${HIP_COMPILER}" STREQUAL "hcc")
            if(NOT "x${HCC_HOME}" STREQUAL "x")
                set(ENV{HCC_HOME} ${HCC_HOME})
            endif()
            set(__CC_FLAGS ${HIP_HIPCC_FLAGS} ${HIP_HCC_FLAGS} ${HIP_HIPCC_FLAGS_${build_configuration}} ${HIP_HCC_FLAGS_${build_configuration}})
        elseif("${HIP_COMPILER}" STREQUAL "clang")
            if(NOT "x${HIP_CLANG_PATH}" STREQUAL "x")
                set(ENV{HIP_CLANG_PATH} ${HIP_CLANG_PATH})
            endif()
            # Temporarily include HIP_HCC_FLAGS for HIP-Clang for PyTorch builds
            set(__CC_FLAGS ${HIP_CLANG_PARALLEL_BUILD_COMPILE_OPTIONS} ${HIP_HIPCC_FLAGS} ${HIP_HCC_FLAGS} ${HIP_CLANG_FLAGS} ${HIP_HIPCC_FLAGS_${build_configuration}} ${HIP_HCC_FLAGS_${build_configuration}} ${HIP_CLANG_FLAGS_${build_configuration}})
        endif()
    else()
        set(__CC_FLAGS ${HIP_HIPCC_FLAGS} ${HIP_NVCC_FLAGS} ${HIP_HIPCC_FLAGS_${build_configuration}} ${HIP_NVCC_FLAGS_${build_configuration}})
    endif()
else()
    set(__CC ${HIP_HOST_COMPILER})
    set(__CC_FLAGS ${CMAKE_HOST_FLAGS} ${CMAKE_HOST_FLAGS_${build_configuration}})
endif()
set(__CC_INCLUDES ${HIP_HIPCC_INCLUDE_ARGS})

# hip_execute_process - Executes a command with optional command echo and status message.
#   status     - Status message to print if verbose is true
#   command    - COMMAND argument from the usual execute_process argument structure
#   ARGN       - Remaining arguments are the command with arguments
#   HIP_result - Return value from running the command
macro(hip_execute_process status command)
    set(_command ${command})
    if(NOT "x${_command}" STREQUAL "xCOMMAND")
        message(FATAL_ERROR "Malformed call to hip_execute_process.  Missing COMMAND as second argument. (command = ${command})")
    endif()
    if(verbose)
        execute_process(COMMAND "${CMAKE_COMMAND}" -E echo -- ${status})
        # Build command string to print
        set(hip_execute_process_string)
        foreach(arg ${ARGN})
            # Escape quotes if any
            string(REPLACE "\"" "\\\"" arg ${arg})
            # Surround args with spaces with quotes
            if(arg MATCHES " ")
                list(APPEND hip_execute_process_string "\"${arg}\"")
            else()
                list(APPEND hip_execute_process_string ${arg})
            endif()
        endforeach()
        # Echo the command
        execute_process(COMMAND ${CMAKE_COMMAND} -E echo ${hip_execute_process_string})
    endif()
    # Run the command
    execute_process(COMMAND ${ARGN} RESULT_VARIABLE HIP_result)
endmacro()

# Delete the target file
hip_execute_process(
    "Removing ${generated_file}"
    COMMAND "${CMAKE_COMMAND}" -E remove "${generated_file}"
    )

# Generate the dependency file
hip_execute_process(
    "Generating dependency file: ${cmake_dependency_file}.pre"
    COMMAND "${__CC}"
    -M
    "${source_file}"
    -o "${cmake_dependency_file}.pre"
    ${__CC_FLAGS}
    ${__CC_INCLUDES}
    )

if(HIP_result)
    message(FATAL_ERROR "Error generating ${generated_file}")
endif()

# Generate the cmake readable dependency file to a temp file
hip_execute_process(
    "Generating temporary cmake readable file: ${cmake_dependency_file}.tmp"
    COMMAND "${CMAKE_COMMAND}"
    -D "input_file:FILEPATH=${cmake_dependency_file}.pre"
    -D "output_file:FILEPATH=${cmake_dependency_file}.tmp"
    -D "verbose=${verbose}"
    -P "${HIP_run_make2cmake}"
    )

if(HIP_result)
    message(FATAL_ERROR "Error generating ${generated_file}")
endif()

# Copy the file if it is different
hip_execute_process(
    "Copy if different ${cmake_dependency_file}.tmp to ${cmake_dependency_file}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${cmake_dependency_file}.tmp" "${cmake_dependency_file}"
    )

if(HIP_result)
    message(FATAL_ERROR "Error generating ${generated_file}")
endif()

# Delete the temporary file
hip_execute_process(
    "Removing ${cmake_dependency_file}.tmp and ${cmake_dependency_file}.pre"
    COMMAND "${CMAKE_COMMAND}" -E remove "${cmake_dependency_file}.tmp" "${cmake_dependency_file}.pre"
    )

if(HIP_result)
    message(FATAL_ERROR "Error generating ${generated_file}")
endif()

# Generate the output file
hip_execute_process(
    "Generating ${generated_file}"
    COMMAND "${__CC}"
    -c
    "${source_file}"
    -o "${generated_file}"
    ${__CC_FLAGS}
    ${__CC_INCLUDES}
    )

if(HIP_result)
    # Make sure that we delete the output file
    hip_execute_process(
        "Removing ${generated_file}"
        COMMAND "${CMAKE_COMMAND}" -E remove "${generated_file}"
        )
    message(FATAL_ERROR "Error generating file ${generated_file}")
else()
    if(verbose)
        message("Generated ${generated_file} successfully.")
    endif()
endif()
# vim: ts=4:sw=4:expandtab:smartindent

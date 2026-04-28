#
# Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

cmake_minimum_required(VERSION 3.23)

project(KRIPKE_SUPERBUILD LANGUAGES CXX)

include(ExternalProject)

option(ENABLE_CALIPER "Enable Caliper support" FALSE)

set(adiak_DIR "${adiak_DIR}" CACHE PATH "Path to an existing adiak CMake config directory")

set(KRIPKE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/cmake/kripke")
set(ADIAK_SOURCE_DIR "${CMAKE_SOURCE_DIR}/tpl/adiak")
set(ADIAK_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/adiak-install")

function(collect_forwarded_cmake_cache_args out_var)
  get_cmake_property(cache_vars CACHE_VARIABLES)
  set(args)

  set(skip_vars
    CMAKE_HOME_DIRECTORY
    CMAKE_CACHEFILE_DIR
    CMAKE_COMMAND
    CMAKE_CTEST_COMMAND
    CMAKE_ROOT
    CMAKE_PROJECT_NAME
    PROJECT_NAME
    CMAKE_SOURCE_DIR
    CMAKE_BINARY_DIR
    CMAKE_CURRENT_SOURCE_DIR
    CMAKE_CURRENT_BINARY_DIR
    KRIPKE_SOURCE_DIR
    ADIAK_SOURCE_DIR
    ADIAK_INSTALL_PREFIX
  )

  foreach(var IN LISTS cache_vars)
    if(var IN_LIST skip_vars)
      continue()
    endif()

    get_property(is_advanced CACHE "${var}" PROPERTY ADVANCED)
    get_property(type CACHE "${var}" PROPERTY TYPE)

    if(NOT type)
      set(type STRING)
    endif()

    if(var STREQUAL "adiak_DIR")
      continue()
    endif()

    set(value "${${var}}")
    list(APPEND args "-D${var}:${type}=${value}")
  endforeach()

  set(${out_var} "${args}" PARENT_SCOPE)
endfunction()

set(ENABLE_GTEST Off CACHE BOOL "")
set(ENABLE_TESTS Off CACHE BOOL "")
set(ENABLE_TBB Off CACHE BOOL "")
set(ENABLE_EXAMPLES Off CACHE BOOL "")
set(ENABLE_DOCUMENTATION Off CACHE BOOL "")

# Use C++17 standard
set(BLT_CXX_STD "c++17" CACHE STRING "")

if (NOT BLT_LOADED)
  if (DEFINED BLT_SOURCE_DIR)
    message(${BLT_SOURCE_DIR})
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
      message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
  else ()
    set (BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/blt" CACHE PATH "")

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
      message(FATAL_ERROR "The BLT submodule is not present.")
    endif ()
 endif ()

  set(BLT_CXX_STD "c++17" CACHE STRING "BLT requires C++17")
  set(ENABLE_ALL_WARNINGS OFF CACHE BOOL "")
  set(ENABLE_ASTYLE OFF CACHE BOOL "")
  set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")
  set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
  set(ENABLE_CPPCHECK OFF CACHE BOOL "")
  set(ENABLE_CLANGQUERY OFF CACHE BOOL "")
  set(ENABLE_FRUIT OFF CACHE BOOL "")
  set(ENABLE_FIND_MPI ON CACHE BOOL "")
  set(BLT_DOCS_TARGET_NAME "docsBLT" CACHE STRING "")

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR} CACHE PATH "")
  include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

endif()

if(ENABLE_CALIPER)
  if(adiak_DIR)
    # Use existing build if adiak_DIR is specified
    find_package(adiak REQUIRED PATHS "${adiak_DIR}" NO_DEFAULT_PATH)
    message(STATUS "Using user-provided adiak_DIR: ${adiak_DIR}")
    set(KRIPKE_ADIAK_DIR "${adiak_DIR}")
  else()
    # else check if another adiak installation is available
    find_package(adiak QUIET)
    if(NOT adiak_FOUND)
      # otherwise build the adiak submodule
      # and set adiak_DIR to point to that build
      message(STATUS "No adiak_DIR provided, building bundled adiak from ${ADIAK_SOURCE_DIR}")
      ExternalProject_Add(adiak_ep
        SOURCE_DIR  "${ADIAK_SOURCE_DIR}"
        BINARY_DIR  "${CMAKE_BINARY_DIR}/adiak-build"
        INSTALL_DIR "${ADIAK_INSTALL_PREFIX}"
        CMAKE_ARGS
          -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
          -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
          -DCMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}
          -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
          -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
          -DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}
          -DCMAKE_SHARED_LINKER_FLAGS=${CMAKE_SHARED_LINKER_FLAGS}
          -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
          -DENABLE_DOCS=OFF
          -DENABLE_TESTS=OFF
          -DENABLE_GTEST=OFF
      )
      set(KRIPKE_ADIAK_DIR "${ADIAK_INSTALL_PREFIX}/lib/cmake/adiak")
    endif()
  endif()
endif()

collect_forwarded_cmake_cache_args(FORWARDED_ARGS)

set(KRIPKE_CMAKE_ARGS
  ${FORWARDED_ARGS}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  -DENABLE_CALIPER:BOOL=${ENABLE_CALIPER}
  -DKRIPKE_SOURCE_ROOT:PATH=${CMAKE_SOURCE_DIR}
)

if(CMAKE_BUILD_TYPE)
  list(APPEND KRIPKE_CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  )
endif()

if(ENABLE_CALIPER AND KRIPKE_ADIAK_DIR)
  list(APPEND KRIPKE_CMAKE_ARGS
    -Dadiak_DIR:PATH=${KRIPKE_ADIAK_DIR}
  )
endif()

if(ENABLE_CALIPER AND NOT adiak_FOUND AND NOT adiak_DIR)
  ExternalProject_Add(kripke_ep
    SOURCE_DIR "${KRIPKE_SOURCE_DIR}"
    BINARY_DIR "${CMAKE_BINARY_DIR}/kripke-build"
    INSTALL_COMMAND ""
    DEPENDS adiak_ep
    CMAKE_ARGS ${KRIPKE_CMAKE_ARGS}
  )
else()
  ExternalProject_Add(kripke_ep
    SOURCE_DIR "${KRIPKE_SOURCE_DIR}"
    BINARY_DIR "${CMAKE_BINARY_DIR}/kripke-build"
    INSTALL_COMMAND ""
    CMAKE_ARGS ${KRIPKE_CMAKE_ARGS}
  )
endif()

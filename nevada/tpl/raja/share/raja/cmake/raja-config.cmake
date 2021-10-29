#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
# Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
#
#=== Usage ===================================================================
# This file allows RAJA to be automatically detected by other libraries
# using CMake.  To build with RAJA, you can do one of two things:
#
#   1. Set the RAJA_DIR environment variable to the root directory of the RAJA
#      installation.  If you loaded RAJA through a dotkit, this may already
#      be set, and RAJA will be autodetected by CMake.
#
#   2. Configure your project with this option:
#      -DRAJA_DIR=<RAJA install prefix>/share/
#
# If you have done either of these things, then CMake should automatically find
# and include this file when you call find_package(RAJA) from your
# CMakeLists.txt file.
#
#=== Components ==============================================================
#
# To link against these, just do, for example:
#
#   find_package(RAJA REQUIRED)
#   add_executable(foo foo.c)
#   target_link_libraries(foo RAJA)
#
# That's all!
#
if (NOT RAJA_CONFIG_LOADED)
  set(RAJA_CONFIG_LOADED TRUE)

  if (NOT DEFINED camp_DIR)
    if (EXISTS )
      set(camp_DIR )
    else ()
      set(camp_DIR /usr/local/lib/cmake/camp)
    endif ()
  endif ()
  find_package(camp REQUIRED PATHS ${camp_DIR} NO_DEFAULT_PATH)
  include(/usr/local/share/raja/cmake/RAJA.cmake)
endif()

# Export version number
set(RAJA_VERSION_MAJOR 0)
set(RAJA_VERSION_MINOR 14)
set(RAJA_VERSION_PATCHLEVEL 0)

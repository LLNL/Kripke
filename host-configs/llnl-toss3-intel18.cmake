#
# Copyright (c) 2014-23, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

set(RAJA_COMPILER "RAJA_COMPILER_INTEL" CACHE STRING "")

set(CMAKE_C_COMPILER   "/usr/tce/bin/icc-18.0.2" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/bin/icpc-18.0.2" CACHE PATH "")
set(CMAKE_LINKER       "/usr/tce/bin/icpc-18.0.2" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

set(ENABLE_OPENMP On CACHE BOOL "")
set(ENABLE_MPI On CACHE BOOL "")



##
## Copyright (c) 2016-22, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## All rights reserved.
##
##

set(RAJA_COMPILER "RAJA_COMPILER_CLANG" CACHE STRING "")

set(CMAKE_C_COMPILER   "/usr/apps/gnu/clang/2017.06.06/llnl/bin/mpiclang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/clang/2017.06.06/llnl/bin/mpiclang++" CACHE PATH "")
set(CMAKE_LINKER       "/usr/apps/gnu/clang/2017.06.06/llnl/bin/mpiclang++" CACHE PATH "")

set(CMAKE_CXX_FLAGS "-stdlib=libc++" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ffast-math" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

set(ENABLE_OPENMP On CACHE BOOL "")
set(ENABLE_MPI_WRAPPER On CACHE BOOL "")



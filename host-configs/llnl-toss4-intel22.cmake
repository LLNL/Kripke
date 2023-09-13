##
## Copyright (c) 2016-22, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## All rights reserved.
##
##

# module load intel/2022.1.0-magic

set(RAJA_COMPILER "RAJA_COMPILER_INTEL" CACHE STRING "")

set(CMAKE_C_COMPILER   "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-2022.1.0-magic/bin/mpicc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-2022.1.0-magic/bin/mpicxx" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

set(ENABLE_OPENMP On CACHE BOOL "")
set(ENABLE_MPI On CACHE BOOL "")



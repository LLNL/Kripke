#
# Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

# module load clang/14.0.6-magic

set(RAJA_COMPILER "RAJA_COMPILER_CLANG" CACHE STRING "")

set(CMAKE_C_COMPILER   "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpiclang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpiclang++" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math --gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ffast-math --gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g --gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")

set(ENABLE_CHAI On CACHE BOOL "")
set(ENABLE_CUDA On CACHE BOOL "")
set(ENABLE_OPENMP Off CACHE BOOL "")
set(ENABLE_MPI On CACHE BOOL "")

set(CMAKE_CUDA_ARCHITECTURES "90" CACHE STRING "")
set(CMAKE_CUDA_FLAGS "-restrict -gencode=arch=compute_90,code=sm_90 --expt-relaxed-constexpr -Xcompiler=--gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler=--gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O3 -lineinfo --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler=--gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_DEBUG "-O0 -g -G --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler=--gcc-toolchain=/usr/tce/packages/gcc/gcc-10.3.1" CACHE STRING "")
set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE STRING "")



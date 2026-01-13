#
# Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

# module load rocmcc/6.4.3-cce-20.0.2-magic

set(RAJA_COMPILER "RAJA_COMPILER_CLANG" CACHE STRING "")

set(BLT_MPI_INCLUDES "/opt/cray/pe/mpich/9.0.1/ofi/cray/20.0/include" CACHE STRING "")

set(CMAKE_C_COMPILER   "/usr/tce/packages/cray-mpich/cray-mpich-9.0.1-rocmcc-6.4.3-cce-20.0.2-magic/bin/mpiamdclang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/cray-mpich/cray-mpich-9.0.1-rocmcc-6.4.3-cce-20.0.2-magic/bin/mpiamdclang++" CACHE PATH "")

set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ffast-math" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,/opt/cray/pe/cce/20.0.2/cce/x86_64/lib -L/opt/cray/pe/cce/20.0.2/cce/x86_64/lib" CACHE STRING "")

set(ROCM_ROOT_DIR "/opt/rocm-6.4.3" CACHE PATH "")
set(HIP_ROOT_DIR "/opt/rocm-6.4.3/hip" CACHE PATH "")
set(HIP_PATH "/opt/rocm-6.4.3/llvm/bin" CACHE PATH "")
set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")
set(GPU_TARGETS "gfx942" CACHE STRING "")
set(AMD_GPU_TARGETS "gfx942" CACHE STRING "")

set(ENABLE_CHAI On CACHE BOOL "")
set(ENABLE_HIP On CACHE BOOL "")
set(ENABLE_OPENMP Off CACHE BOOL "")
set(ENABLE_MPI On CACHE BOOL "")

#set(CMAKE_HIPCC_FLAGS_RELEASE "-O3 --expt-extended-lambda" CACHE STRING "")
#set(CMAKE_HIPCC_FLAGS_RELWITHDEBINFO "-O3 -lineinfo --expt-extended-lambda" CACHE STRING "")
#set(CMAKE_HIPCC_FLAGS_DEBUG "-O0 -g -G --expt-extended-lambda" CACHE STRING "")
#set(CMAKE_HIPCC_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE STRING "")

# For LLNL TCE packages
set(ENABLE_MPI_WRAPPER On CACHE BOOL "")



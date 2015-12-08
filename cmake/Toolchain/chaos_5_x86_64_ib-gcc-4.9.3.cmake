set(GCC_PATH /usr/apps/gnu/4.9.3 CACHE STRING foobar FORCE)

set(CMAKE_CXX_COMPILER "${GCC_PATH}/bin/mpig++" CACHE STRING barfoo FORCE)
set(CMAKE_LINKER "${GCC_PATH}/bin/mpig++" CACHE STRING barfoo FORCE)

set(CMAKE_CXX_FLAGS "-mpi=mvapich2-gnu-2.1 -g -O3 -std=c++11 -mtune=native -mavx -DRAJA_PLATFORM_X86_SSE -DRAJA_COMPILER_GNU" CACHE STRING barfoo FORCE)

set(CMAKE_LINKER_FLAGS "-mpi=mvapich2-gnu-2.1" CACHE STRING barfoo FORCE)

# CUDA Compiler Setup
set(CUDA_NVCC_FLAGS "-arch compute_35 -std=c++11 --expt-extended-lambda" CACHE STRING foobar FORCE)
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING barfoo FORCE)


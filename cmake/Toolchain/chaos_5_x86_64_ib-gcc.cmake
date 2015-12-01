set(CMAKE_CXX_COMPILER "mpig++" CACHE STRING barfoo FORCE)
set(CMAKE_LINKER "mpig++" CACHE STRING barfoo FORCE)

set(CMAKE_CXX_FLAGS "-mpi=mvapich2-gnu-2.1 -g -O3 -mtune=native" CACHE STRING barfoo FORCE)

set(CMAKE_LINKER_FLAGS "-mpi=mvapich2-gnu-2.1" CACHE STRING barfoo FORCE)

# CUDA Compiler Setup
set(CUDA_NVCC_FLAGS "-arch compute_35" CACHE STRING foobar FORCE)
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING barfoo FORCE)


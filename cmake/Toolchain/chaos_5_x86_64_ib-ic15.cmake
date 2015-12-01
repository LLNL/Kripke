set(ICC_VER 15.0.223 CACHE STRING foobar FORCE)

set(CMAKE_C_COMPILER "mpiicc-${ICC_VER}" CACHE STRING barfoo FORCE)
set(CMAKE_CXX_COMPILER "mpiicpc-${ICC_VER}" CACHE STRING barfoo FORCE)
set(CMAKE_LINKER "mpiicpc-${ICC_VER}" CACHE STRING barfoo FORCE)

set(CMAKE_CXX_FLAGS "-mpi=mvapich2-intel-2.1 -gxx-name=g++-4.6.1 -gcc-name=gcc-4.6.1 -g -O3 -unroll-aggressive -finline-functions -axAVX -msse4.2" CACHE STRING barfoo FORCE)

set(CMAKE_LINKER_FLAGS "-mpi=mvapich2-intel-2.1 -gxx-name=g++-4.6.1 -gcc-name=gcc-4.6.1" CACHE STRING barfoo FORCE)

# CUDA Compiler Setup
set(CUDA_NVCC_FLAGS "-arch compute_35" CACHE STRING foobar FORCE)
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING barfoo FORCE)


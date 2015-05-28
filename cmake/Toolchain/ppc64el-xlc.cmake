
set(CMAKE_C_COMPILER "mpicc")

set(CMAKE_CXX_COMPILER "mpicxx") 

#set(CMAKE_NVCC_COMPILER "/usr/local/cuda-7.0/bin/nvcc")

set(CMAKE_LINKER mpicxx)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -v -O3  -std=c++0x -qsmp=omp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -O3 -std=c++0x -qsmp=omp")
#set(CMAKE_CUDA_NVCC_FLAGS "${CMAKE_CUDA_NVCC_FLAGS} -gencode arch=compute_35,code=sm_35")



set(PKG_PATH "/usr/lib/")


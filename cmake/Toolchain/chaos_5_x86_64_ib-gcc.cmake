set(CMAKE_C_COMPILER "mpigcc")
set(CMAKE_CXX_COMPILER "mpig++")
set(CMAKE_LINKER "mpig++")

set(CMAKE_C_FLAGS "-O3 -g -fopenmp -DKRIPKE_USE_OPENMP -DKRIPKE_USE_CUDA -g -mtune=native ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-O3 -g -fopenmp -DKRIPKE_USE_OPENMP -DKRIPKE_USE_CUDA -std=c++11 -mtune=native ${CMAKE_CXX_FLAGS}")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/gnu-4.9-opt")




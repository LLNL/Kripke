
set(CMAKE_C_COMPILER mpicc)

set(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_LINKER mpicxx)

set(CMAKE_C_FLAGS -O3) 

set(CMAKE_CXX_FLAGS -O3)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mtune=power8 -mcpu=power8 -ftree-vectorizer-verbose=1 -finline-functions -finline-limit=20000 -std=c++0x -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=power8 -mcpu=power8 -ftree-vectorizer-verbose=1 -finline-functions -finline-limit=20000 -std=c++0x -fopenmp")

set(PKG_PATH "/usr/lib/")


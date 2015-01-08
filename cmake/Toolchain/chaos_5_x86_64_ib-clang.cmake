
set(CMAKE_C_COMPILER "mpiclang")

set(CMAKE_CXX_COMPILER "mpiclang++")

set(CMAKE_LINKER "mpiclang++ -pie")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++11")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")


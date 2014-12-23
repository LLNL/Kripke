
set(CMAKE_C_COMPILER "mpiclang")

set(CMAKE_CXX_COMPILER "mpiclang++11")

set(CMAKE_LINKER "mpiclang++11")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -static -Wswitch -g ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -static -Wswitch -g")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")


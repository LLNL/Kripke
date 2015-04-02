set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_LINKER "mpic++")

set(CMAKE_C_FLAGS "-g -mtune=native ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-g -O0 ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-g -std=c++11 -mtune=native ${CMAKE_CXX_FLAGS}")

set(PKG_PATH "/usr/local")



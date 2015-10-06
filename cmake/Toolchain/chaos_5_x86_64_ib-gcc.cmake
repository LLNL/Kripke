set(CMAKE_C_COMPILER "mpigcc")
set(CMAKE_CXX_COMPILER "mpig++")
set(CMAKE_LINKER "mpig++")

set(CMAKE_C_FLAGS "-g -mtune=native ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "-g -std=c++11 -mtune=native ${CMAKE_CXX_FLAGS}")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/gnu-4.9-opt")



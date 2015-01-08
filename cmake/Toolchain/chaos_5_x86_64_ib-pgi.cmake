
set(CMAKE_C_COMPILER "mpipgcc")

set(CMAKE_CXX_COMPILER "mpipgCC")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g --c++11 -fast -Mipa=fast,inline")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")



set(CMAKE_C_COMPILER mpixlcxx_r)

set(CMAKE_CXX_COMPILER mpixlcxx_r)

set(CMAKE_LINKER mpixlcxx_r)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -qarch=auto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qarch=auto")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")


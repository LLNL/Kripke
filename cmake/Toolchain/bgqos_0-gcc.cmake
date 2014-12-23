
set(CMAKE_C_COMPILER mpigcc-4.7.2-fastmpi)

set(CMAKE_CXX_COMPILER mpig++-4.7.2-fastmpi)

set(CMAKE_LINKER mpig++-4.7.2-fastmpi)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcpu=a2 -mtune=a2 -finline-functions -finline-limit=20000 -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcpu=a2 -mtune=a2 -finline-functions -finline-limit=20000 -std=c++11 -ftree-vectorizer-verbose=6")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")


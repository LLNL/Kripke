cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_COMPILER mpiicc-12.1.339)
set(CMAKE_CXX_COMPILER mpiicpc-12.1.339)
set(CMAKE_LINKER mpiicpc-12.1.339)

set(SYS_TYPE $ENV{SYS_TYPE})
set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")

project(KRIPKE)

enable_language(CXX)

set(INCLUDE_PATHS ${INCLUDE_PATHS} /usr/gapps/bdiv/chaos_5_x86_64_ib/opt/google-perftools/1.9/include)
link_directories(/usr/gapps/bdiv/chaos_5_x86_64_ib/opt/google-perftools/1.9/lib)

set(KRIPKE_BASE_LIBS kripke lapack profiler unwind) 

set(KRIPKE_LIBS ${KRIPKE_BASE_LIBS})

include_directories(${INCLUDE_PATHS})

add_subdirectory(src)


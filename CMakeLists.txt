cmake_minimum_required(VERSION 2.8)

set(ICC_VER 14.0.144)

set(CMAKE_C_COMPILER mpiicc-${ICC_VER})
set(CMAKE_CXX_COMPILER mpiicpc-${ICC_VER})
set(CMAKE_LINKER mpiicpc-${ICC_VER})

set(SYS_TYPE $ENV{SYS_TYPE})
set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")

project(KRIPKE)

enable_language(CXX)

set(INCLUDE_PATHS ${INCLUDE_PATHS} 
	/usr/gapps/bdiv/chaos_5_x86_64_ib/opt/google-perftools/1.9/include
	/usr/gapps/bdiv/chaos_5_x86_64_ib/opt/papi/4.2.0/static/include
)

link_directories(
	/usr/gapps/bdiv/chaos_5_x86_64_ib/opt/google-perftools/1.9/lib
	/usr/gapps/bdiv/chaos_5_x86_64_ib/opt/papi/4.2.0/static/lib
)



set(KRIPKE_BASE_LIBS kripke lapack profiler unwind papi) 

set(KRIPKE_LIBS ${KRIPKE_BASE_LIBS})

include_directories(${INCLUDE_PATHS})

add_subdirectory(src)


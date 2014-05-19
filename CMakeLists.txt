cmake_minimum_required(VERSION 2.8)

set(ICC_VER 14.0.174)

set(CMAKE_C_COMPILER mpiicc-${ICC_VER})
set(CMAKE_CXX_COMPILER mpiicpc-${ICC_VER})
set(CMAKE_LINKER mpiicpc-${ICC_VER})

set(SYS_TYPE $ENV{SYS_TYPE})
set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/opt")


project(KRIPKE)

enable_language(CXX)

set(INCLUDE_PATHS ${INCLUDE_PATHS} 
	/usr/gapps/bdiv/${SYS_TYPE}/opt/google-perftools/1.9/include
	/usr/gapps/bdiv/${SYS_TYPE}/opt/papi/4.2.0/static/include
)

link_directories(
	/usr/gapps/bdiv/${SYS_TYPE}/opt/google-perftools/1.9/lib
	/usr/gapps/bdiv/${SYS_TYPE}/opt/papi/4.2.0/static/lib
)



set(KRIPKE_BASE_LIBS kripke lapack profiler unwind papi) 

set(KRIPKE_LIBS ${KRIPKE_BASE_LIBS})

include_directories(${INCLUDE_PATHS})



#
# OpenMP
# 
option(ENABLE_OPENMP "Turn on compiler support for OpenMP" OFF)

if(${ENABLE_OPENMP})
  find_package(OpenMP)
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
      add_definitions (-DKRIPKE_USE_OPENMP=1)
  endif()
endif()


add_subdirectory(src)


cmake_minimum_required(VERSION 3.0)

project(KRIPKE)


find_package(MPI)
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set(KRIPKE_LIBS ${KRIPKE_LIBS} ${MPI_LIBRARIES})
  
  # Dsplay details
  message(STATUS "MPI C++ Compiler Flags: ${MPI_CXX_COMPILE_FLAGS}")
  message(STATUS "MPI C++ Include Paths: ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "MPI C++ Linker Flags: ${MPI_CXX_LINK_FLAGS}")
  message(STATUS "MPI C++ Libraries: ${MPI_CXX_LIBRARIES}")  
endif()

find_package(OpenMP)

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "tarball.py;kripke-tarball;.git;.cproject;.project;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

message(STATUS "Selected Build Type ${CMAKE_BUILD_TYPE}")
message(STATUS "Selected C++ Compiler ${CMAKE_CXX_COMPILER}")
message(STATUS "Selected C++ Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Selected BDIV Package Path ${PKG_PATH}")


#
# External libraries 
#  

# Typically used (turned on by default)
#bdiv_opt_pkg(tcmalloc 0 OFF)

# Optional 
#bdiv_opt_pkg(memP 0.0 OFF)
  #bdiv_opt_pkg(papi 4.2.0 OFF)
#bdiv_opt_pkg(perftools 0 OFF)
#bdiv_opt_pkg(SILO 4.9 OFF)
#bdiv_opt_pkg(HDF5 gnu-serial-1.8.8 OFF)

# BG/Q Specific
#bdiv_opt_pkg(BGPM 20030410 OFF)


#
# Other Build Options
# 

# OpenMP build option
#include(cmake/OpenMP.cmake)

# Traverse subdirectories
include_directories(src)
add_subdirectory(src)

# Setup the kripke target
add_executable(kripke "src/kripke.cpp")
target_link_libraries(kripke ${KRIPKE_LIBS} ${KRIPKE_LIBS})




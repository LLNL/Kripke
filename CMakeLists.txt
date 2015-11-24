cmake_minimum_required(VERSION 3.1)

project(KRIPKE LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})


message(STATUS "Selected toolchain file ${CMAKE_TOOLCHAIN_FILE}")
message(STATUS "Selected Build Type ${CMAKE_BUILD_TYPE}")
message(STATUS "Selected C++ Compiler ${CMAKE_CXX_COMPILER}")
message(STATUS "Selected C++ Flags: ${CMAKE_CXX_FLAGS}")


#
# Cuda
# 
option(ENABLE_CUDA "Turn on compiler support for CUDA" OFF)
message(STATUS "CUDA Support is ${ENABLE_CUDA}")

if(${ENABLE_CUDA})
  find_package(CUDA QUIET REQUIRED)
  if (CUDA_FOUND)    
      add_definitions (-DKRIPKE_USE_CUDA -DLPlusTimes_sweep_combined_ -DLPlusTimes_sweep_LTimes_combined__  -DKRIPKE_USE_CUBLAS)
  endif()
  message(STATUS "Selected NVCC Compiler ${CUDA_NVCC_EXECUTABLE}")
  message(STATUS "Selected CUDA NVCC Flags: ${CUDA_NVCC_FLAGS}")
endif()



#
# Invoke the project, choose language, display some debug info
#

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "tarball.py;kripke-tarball;.git;.cproject;.project;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)


#
# Libraries to include
#  

# Include macros to aid in BDIV packages
#set(BDIV_LIBRARY_PATH "${KRIPKE_SOURCE_DIR}/cmake/Library")
#include(cmake/bdiv_pkg.cmake)

# Typically used (turned on by default)
#bdiv_opt_pkg(tcmalloc 0 OFF)

# Optional 
#bdiv_opt_pkg(papi 4.2.0 OFF)
#bdiv_opt_pkg(perftools 0 OFF)
#bdiv_opt_pkg(SILO 4.9 OFF)
#bdiv_opt_pkg(HDF5 gnu-serial-1.8.8 OFF)


# BG/Q Specific
#bdiv_opt_pkg(BGPM 20030410 OFF)

# BLAS 
#bdiv_opt_pkg(esslblas 0 OFF)

#HWLOC
#bdiv_opt_pkg(hwloc 0 OFF)

#HPC
#bdiv_opt_pkg(hpclib 0 OFF)



#
# Other Build Options
# 

# OpenMP build option
include(cmake/OpenMP.cmake)
message(STATUS "Selected C++ Compiler ${CMAKE_CXX_COMPILER}")
message(STATUS "Selected C++ Flags: ${CMAKE_CXX_FLAGS}")

#
# Process build files for Kripke
#
list(APPEND KRIPKE_LIBS kripkelib)
list(APPEND KRIPKE_LIBS ${KRIPKE_LIBS} ${KRIPKE_LIBS})



add_subdirectory(src)

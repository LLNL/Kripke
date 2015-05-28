project(KRIPKE)

cmake_minimum_required(VERSION 2.8)


#project(KRIPKE)

#
# Figure out toolchain files to use in this build
#

# CUDA build option
include(cmake/CUDA.cmake)


set(SYS_TYPE $ENV{SYS_TYPE})



set(CMAKE_TOOLCHAIN_PATH "./cmake/Toolchain")

set(CMAKE_C_FLAGS -O3)
set(CMAKE_CXX_FLAGS -O3)


if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)

  # If the user didn't specify a toolchain, then default the toolchain file
  if(NOT DEFINED TOOLCHAIN)  
    if(SYS_TYPE MATCHES bgqos_0)
      set(CMAKE_TOOLCHAIN_FILE "${CMAKE_TOOLCHAIN_PATH}/bgqos_0-xlc.cmake")
    elseif(SYS_TYPE MATCHES chaos_5_x86_64_ib)
      set(CMAKE_TOOLCHAIN_FILE "${CMAKE_TOOLCHAIN_PATH}/chaos_5_x86_64_ib-ic14.cmake")
      include("${CMAKE_TOOLCHAIN_PATH}/chaos_5_x86_64_ib-ic14.cmake")
    elseif(SYS_TYPE MATCHES ppc64el)
      set(CMAKE_TOOLCHAIN_FILE "${CMAKE_TOOLCHAIN_PATH}/ppc64el-gcc.cmake")
      include("${CMAKE_TOOLCHAIN_PATH}/ppc64el-gcc.cmake")
    endif()
    
  # If the user selected a toolchain (ex. ic14, gcc) then prepend the systype
  else()
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_TOOLCHAIN_PATH}/${SYS_TYPE}-${TOOLCHAIN}.cmake")
  endif()
  
  message(STATUS "Selected toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()


#
# Invoke the project, choose language, display some debug info
#


# Re-incluing the toolchain file, which allows us to modify the compiler flags
# since all of the flags are lost after calling project()
include(${CMAKE_TOOLCHAIN_FILE})


set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "tarball.py;kripke-tarball;.git;.cproject;.project;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

message(STATUS "Selected Build Type ${CMAKE_BUILD_TYPE}")
message(STATUS "Selected C++ Compiler ${CMAKE_CXX_COMPILER}")
message(STATUS "Selected C++ Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Selected NVCC Compiler ${CUDA_NVCC_EXECUTABLE}")
message(STATUS "Selected CUDA NVCC Flags: ${CMAKE_NVCC_FLAGS}")
message(STATUS "Selected BDIV Package Path ${PKG_PATH}")



#
# Libraries to include
#  

# Include macros to aid in BDIV packages
set(BDIV_LIBRARY_PATH "${KRIPKE_SOURCE_DIR}/cmake/Library")
include(cmake/bdiv_pkg.cmake)

# Typically used (turned on by default)
bdiv_opt_pkg(tcmalloc 0 OFF)

# Optional 
bdiv_opt_pkg(memP 0.0 OFF)
bdiv_opt_pkg(papi 4.2.0 OFF)
bdiv_opt_pkg(perftools 0 OFF)
bdiv_opt_pkg(SILO 4.9 OFF)
bdiv_opt_pkg(HDF5 gnu-serial-1.8.8 OFF)


# BG/Q Specific
bdiv_opt_pkg(BGPM 20030410 OFF)

# BLAS 
bdiv_opt_pkg(esslblas 0 OFF)

#HWLOC
bdiv_opt_pkg(hwloc 0 ON)

#HPC
bdiv_opt_pkg(hpclib 0 ON)



#
# Other Build Options
# 

# OpenMP build option
include(cmake/OpenMP.cmake)


#
# Process build files for Kripke
#
set(KRIPKE_LIBS kripkelib ${BDIV_LIBS})
set(KRIPKE_LIBS ${KRIPKE_LIBS} ${KRIPKE_LIBS})



add_subdirectory(src)

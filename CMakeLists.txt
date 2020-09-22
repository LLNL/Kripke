#
# Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0057 NEW)

project(KRIPKE LANGUAGES CXX)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})


set(KRIPKE_VERSION "1.2.4")


#
# Initialize the BLT build system
#

set(ENABLE_GTEST Off CACHE Bool "")
set(ENABLE_TESTS Off CACHE Bool "")
set(ENABLE_TBB Off CACHE Bool "")
set(ENABLE_EXAMPLES Off CACHE Bool "")
set(ENABLE_DOCUMENTATION Off CACHE Bool "")

# Use C++14 standard
set(BLT_CXX_STD "c++14" CACHE STRING "")

set(chai_DIR ${CURRENT_SOURCE_DIR}/tpl/chai)

include(blt/SetupBLT.cmake)



#
# Default Arch and Layout selection
# Sequential by default, but will be overriden if OpenMP or CUDA are enabled
#
set(KRIPKE_ARCH "Sequential")
set(KRIPKE_LAYOUT DGZ)


#
# Add CMake CUDA language support
#
if(ENABLE_CUDA)
  # Native CMake support  
  enable_language(CUDA)
  
  # Make sure we use nvcc with our selected host compiler to link executables
  set(CMAKE_CUDA_LINK_EXECUTABLE "${CMAKE_CUDA_COMPILER} -ccbin <CMAKE_CXX_COMPILER> <CMAKE_CUDA_LINK_FLAGS>  <FLAGS>  <LINK_FLAGS>  <OBJECTS> -o <TARGET>  <LINK_LIBRARIES>")
  # Use nvcc as the linker
  set(CMAKE_LINKER "${CMAKE_CUDA_COMPILER}" CACHE STRING "")
  
  # Pass nvcc the -cubin to point at our C++ compiler  
  set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

endif()



#
# Configure CHAI/Umpire for memory management
# 
option(ENABLE_CHAI "Enable CHAI/Umpire memory management" Off)
if(ENABLE_CHAI)

  # Add CHAI (which brings in Umpire as a submodule)
  add_subdirectory(tpl/chai)
  
  set(ENABLE_TESTS Off CACHE Bool "")
  
  set(KRIPKE_USE_CHAI 1)
  
  list(APPEND KRIPKE_DEPENDS chai)
  
endif()


#
# Configure RAJA (REQUIRED)
#

set(RAJA_ENABLE_TESTS Off CACHE Bool "")
add_subdirectory(tpl/raja)

list(APPEND KRIPKE_DEPENDS RAJA)


#
# Configure OpenMP (CPU threading, not target offload)
#

if(ENABLE_OPENMP)
  
  set(KRIPKE_USE_OPENMP 1)  
  
  list(APPEND KRIPKE_DEPENDS openmp)

  set(KRIPKE_ARCH "OpenMP")
  
endif()



#
# Configure CUDA
#

if(ENABLE_CUDA)

  set(KRIPKE_USE_CUDA 1)
    
  list(APPEND KRIPKE_DEPENDS cuda)
  
  set(KRIPKE_ARCH "CUDA")

  # Make sure that nvcc turns on the host compiler OpenMP flag	
  if(ENABLE_OPENMP)
    list(APPEND CUDA_NVCC_FLAGS -Xcompiler ${OpenMP_CXX_FLAGS})
  endif()
	
	
  # Pass on compiler info to KripkeConfig.h
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
  set(KRIPKE_NVCC_COMPILER "${CMAKE_CUDA_COMPILER}")
  set(KRIPKE_NVCC_FLAGS    "${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")

endif()





option(ENABLE_CALIPER "Enable Caliper support" FALSE)

if (ENABLE_CALIPER)
  find_package(caliper REQUIRED)

  set(KRIPKE_USE_CALIPER TRUE)
  
  blt_register_library(NAME caliper
    INCLUDES ${caliper_INCLUDE_DIR}
    LIBRARIES caliper)
  
  list(APPEND KRIPKE_DEPENDS caliper)
  
  if (ENABLE_MPI)
    list(APPEND KRIPKE_DEPENDS caliper-mpi)
  endif()
endif()

message(STATUS "Caliper support is ${KRIPKE_USE_CALIPER}")


#
# Configure MPI
# 
# Use ENABLE_MPI=On if you want CMake to automatically figure out MPI
# using FindMPI
#
# Use ENABLE_MPI_WRAPPER if you want to use mpi compiler wrappers for CC/CXX 
# 

if(ENABLE_MPI)
	set(KRIPKE_USE_MPI 1)	
	list(APPEND KRIPKE_DEPENDS mpi)
endif()

if(ENABLE_MPI_WRAPPER)
	set(KRIPKE_USE_MPI 1)
endif()









#
# Display/Configure our default ArchLayout
#
message(STATUS "Kripke selected default architecture: '${KRIPKE_ARCH}'")
message(STATUS "Kripke selected default layout:       '${KRIPKE_LAYOUT}'")

set(KRIPKE_ARCHV_DEFAULT ArchV_${KRIPKE_ARCH})
set(KRIPKE_LAYOUTV_DEFAULT LayoutV_${KRIPKE_LAYOUT})




#
# Pass compiler options to our KripkeConfig.h file so we can print them at
# runtime
#


set(KRIPKE_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
set(KRIPKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${KRIPKE_BUILD_TYPE}}")
set(KRIPKE_LINK_FLAGS   "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_LINK_FLAGS_${CMAKE_BUILD_TYPE}}")







#
# Create KripkeConfig.h file
#

configure_file(${PROJECT_SOURCE_DIR}/src/KripkeConfig.h.in
               ${PROJECT_BINARY_DIR}/include/KripkeConfig.h)






#
# Define build targets for Kripke
#
 


blt_add_library(
  NAME        kripke
  SOURCES     "src/Kripke/Core/BaseVar.cpp"
              "src/Kripke/Core/DataStore.cpp"
              "src/Kripke/Core/DomainVar.cpp"
              "src/Kripke/Generate.cpp"
              "src/Kripke/Generate/Data.cpp"
              "src/Kripke/Generate/Decomp.cpp"
              "src/Kripke/Generate/Energy.cpp"
              "src/Kripke/Generate/Quadrature.cpp"
              "src/Kripke/Generate/Space.cpp"
              "src/Kripke/InputVariables.cpp"
              "src/Kripke/Kernel/LPlusTimes.cpp"
              "src/Kripke/Kernel/LTimes.cpp"
              "src/Kripke/Kernel/Population.cpp"
              "src/Kripke/Kernel/Scattering.cpp"
              "src/Kripke/Kernel/Source.cpp"
              "src/Kripke/Kernel/SweepSubdomain.cpp"
              "src/Kripke/Kernel/ErrorNorms.cpp"
              "src/Kripke/ParallelComm/BlockJacobiComm.cpp"
              "src/Kripke/ParallelComm/SweepComm.cpp"
              "src/Kripke/ParallelComm.cpp"
              "src/Kripke/Core/PartitionSpace.cpp"
              "src/Kripke/Core/Set.cpp"
              "src/Kripke/SteadyStateSolver.cpp"
              "src/Kripke/SweepSolver.cpp"
              "src/Kripke/Timing.cpp"
  DEPENDS_ON  ${KRIPKE_DEPENDS}
)


target_include_directories(kripke PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>)

blt_add_executable( 
  NAME        kripke.exe
  SOURCES     "src/kripke.cpp"
  DEPENDS_ON  ${KRIPKE_DEPENDS} kripke 
)



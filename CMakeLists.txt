#
# Copyright (c) 2014-22, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0057 NEW)

project(KRIPKE LANGUAGES CXX)

set(KRIPKE_VERSION "1.2.5-dev")


#
# Initialize the BLT build system
#

set(ENABLE_GTEST Off CACHE BOOL "")
set(ENABLE_TESTS Off CACHE BOOL "")
set(ENABLE_TBB Off CACHE BOOL "")
set(ENABLE_EXAMPLES Off CACHE BOOL "")
set(ENABLE_DOCUMENTATION Off CACHE BOOL "")

if (ENABLE_CHAI)
  set(ENABLE_RAJA_PLUGIN On CACHE BOOL "")
  set(CHAI_ENABLE_RAJA_PLUGIN On CACHE BOOL "")
endif ()

# Use C++14 standard
set(BLT_CXX_STD "c++14" CACHE STRING "")

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
# Add googletest (needed by RAJA, etc)
# 
add_subdirectory(tpl/googletest)


#
# Configure CHAI/Umpire for memory management
# 
option(ENABLE_CHAI "Enable CHAI/Umpire memory management" Off)
if(ENABLE_CHAI)
  # Find camp
  if (NOT TARGET camp)
    find_package(camp QUIET NO_DEFAULT_PATH HINTS ${CAMP_DIR} ${camp_DIR})

    if (camp_FOUND)
      message(STATUS "Kripke: Using external CAMP")
      set_target_properties(camp PROPERTIES IMPORTED_GLOBAL TRUE)
    else ()
      message(STATUS "Kripke: Using CAMP submodule")

      if (NOT EXISTS ${PROJECT_SOURCE_DIR}/tpl/raja/tpl/camp/CMakeLists.txt)
        message(FATAL_ERROR "Kripke: CAMP submodule not initialized. Run 'git submodule update --init' in the git repository or set camp_DIR or CAMP_DIR to use an external build of CAMP.")
      else ()
        message(STATUS "Kripke: Using CAMP in Kripke TPL")
        add_subdirectory(${PROJECT_SOURCE_DIR}/tpl/raja/tpl/camp)
      endif ()
    endif()
  endif()

  # Set camp backends
  if (ENABLE_CUDA)
    message(STATUS "Kripke: Setting CAMP_HAVE_CUDA")
    blt_add_target_definitions(TO camp
                               SCOPE INTERFACE
                               TARGET_DEFINITIONS CAMP_HAVE_CUDA)
  endif ()

  if (ENABLE_HIP)
    message(STATUS "Kripke: Setting CAMP_HAVE_HIP")
    blt_add_target_definitions(TO camp
                               SCOPE INTERFACE
                               TARGET_DEFINITIONS CAMP_HAVE_HIP)
  endif ()

  set(UMPIRE_DEPENDS camp)

  if (DEFINED umpire_DIR)
    find_package(umpire REQUIRED)
    
    if (ENABLE_MPI)
      set(UMPIRE_DEPENDS "${UMPIRE_DEPENDS} mpi")
    endif()

    blt_register_library(
      NAME umpire
      INCLUDES ${UMPIRE_INCLUDE_DIRS}
      LIBRARIES umpire
      DEPENDS_ON ${UMPIRE_DEPENDS})
  else ()
    add_subdirectory(tpl/umpire)
  endif()

  if (DEFINED chai_DIR)
    find_package(chai REQUIRED)
  else ()
    add_subdirectory(tpl/chai)
  endif ()

  set(ENABLE_TESTS Off CACHE BOOL "")
  set(KRIPKE_USE_CHAI 1)
  list(APPEND KRIPKE_DEPENDS chai)
else ()
  #
  # Configure RAJA (REQUIRED)
  #
  set(ENABLE_RAJA_PLUGIN Off CACHE BOOL "")
  if (DEFINED RAJA_DIR)
    find_package(RAJA REQUIRED)
  else ()
    set(ENABLE_EXAMPLES Off CACHE BOOL "")
    set(ENABLE_EXERCISES Off CACHE BOOL "")
    add_subdirectory(tpl/raja)
  endif ()
  list(APPEND KRIPKE_DEPENDS RAJA)
endif ()

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


#
# Configure HIP
#

if(ENABLE_HIP)

  set(KRIPKE_USE_HIP 1)
    
  list(APPEND KRIPKE_DEPENDS hip)
  
  set(KRIPKE_ARCH "HIP")

  # Make sure that nvcc turns on the host compiler OpenMP flag	
  #if(ENABLE_OPENMP)
  #  list(APPEND CUDA_NVCC_FLAGS -Xcompiler ${OpenMP_CXX_FLAGS})
  #endif()
  #  
  #  
  ## Pass on compiler info to KripkeConfig.h
  string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
  #set(KRIPKE_NVCC_COMPILER "${CMAKE_CUDA_COMPILER}")
  #set(KRIPKE_NVCC_FLAGS    "${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")

endif()



option(ENABLE_CALIPER "Enable Caliper support" FALSE)

if (ENABLE_CALIPER)
	if(EXISTS "${PROJECT_SOURCE_DIR}/tpl/caliper/CMakeLists.txt")
		set(BUILD_SHARED_LIBS False CACHE BOOL "")
    add_subdirectory(tpl/caliper)	
		include_directories(tpl/caliper/include)
	else()
    find_package(caliper REQUIRED)
	
 	endif()

    blt_register_library(NAME caliper
      INCLUDES ${caliper_INCLUDE_DIR}
      LIBRARIES caliper)

	set(KRIPKE_USE_CALIPER TRUE)
  
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


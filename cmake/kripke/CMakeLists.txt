#
# Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

cmake_minimum_required(VERSION 3.23)

cmake_policy(SET CMP0057 NEW)

project(KRIPKE LANGUAGES CXX)

set(KRIPKE_VERSION "1.2.5-dev")

if(NOT DEFINED KRIPKE_SOURCE_ROOT)
  message(FATAL_ERROR "KRIPKE_SOURCE_ROOT is not set")
endif()

set(KRIPKE_SOURCE_DIR "${KRIPKE_SOURCE_ROOT}")
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

# Use C++17 standard
set(BLT_CXX_STD "c++17" CACHE STRING "")

if (NOT BLT_LOADED)
  if (DEFINED BLT_SOURCE_DIR)
    message(${BLT_SOURCE_DIR})
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
      message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
  else ()
    set (BLT_SOURCE_DIR "${KRIPKE_SOURCE_DIR}/blt" CACHE PATH "")

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
      message(FATAL_ERROR "The BLT submodule is not present.")
    endif ()
 endif ()

  set(BLT_CXX_STD "c++17" CACHE STRING "BLT requires C++17")
  set(ENABLE_ALL_WARNINGS OFF CACHE BOOL "")
  set(ENABLE_ASTYLE OFF CACHE BOOL "")
  set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")
  set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
  set(ENABLE_CPPCHECK OFF CACHE BOOL "")
  set(ENABLE_CLANGQUERY OFF CACHE BOOL "")
  set(ENABLE_GTEST ON CACHE BOOL "")
  set(ENABLE_TESTS ON CACHE BOOL "")
  set(ENABLE_FRUIT OFF CACHE BOOL "")
  set(ENABLE_FIND_MPI ON CACHE BOOL "")
  set(BLT_DOCS_TARGET_NAME "docsBLT" CACHE STRING "")

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR} CACHE PATH "")
  set(TEST_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR} CACHE PATH "")

  include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

endif()


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
  #enable_language(CUDA)
  
  # Make sure we use nvcc with our selected host compiler to link executables
  set(CMAKE_CUDA_LINK_EXECUTABLE "${CMAKE_CUDA_COMPILER} -ccbin <CMAKE_CXX_COMPILER> <CMAKE_CUDA_LINK_FLAGS>  <FLAGS>  <LINK_FLAGS>  <OBJECTS> -o <TARGET>  <LINK_LIBRARIES>")
  # Use nvcc as the linker
  set(CMAKE_LINKER "${CMAKE_CUDA_COMPILER}" CACHE STRING "")
  
  # Pass nvcc the -cubin to point at our C++ compiler  
  set(CUDA_HOST_COMPILER ${CMAKE_CUDA_HOST_COMPILER})

endif()

#
# Add googletest (needed by RAJA, etc)
# 
#add_subdirectory(${KRIPKE_SOURCE_DIR}/tpl/googletest)


#
# Configure CHAI/Umpire for memory management
# 
option(ENABLE_CHAI "Enable CHAI/Umpire memory management" Off)
option(ENABLE_GPU_AWARE_MPI "Pass GPU pointers directly to MPI communication buffers" Off)
if(ENABLE_CHAI)
  # Find camp
  if (NOT TARGET camp)
    find_package(camp QUIET NO_DEFAULT_PATH HINTS ${CAMP_DIR} ${camp_DIR})

    if (camp_FOUND)
      message(STATUS "Kripke: Using external CAMP")
      set_target_properties(camp PROPERTIES IMPORTED_GLOBAL TRUE)
      if (DEFINED RAJA_DIR)
        find_package(RAJA REQUIRED)
      endif()

      list(APPEND KRIPKE_DEPENDS RAJA)
    else ()
      message(STATUS "Kripke: Using CAMP submodule")

      if (NOT EXISTS ${KRIPKE_SOURCE_DIR}/tpl/raja/tpl/camp/CMakeLists.txt)
        message(FATAL_ERROR "Kripke: CAMP submodule not initialized. Run 'git submodule update --init' in the git repository or set camp_DIR or CAMP_DIR to use an external build of CAMP.")
      else ()
        message(STATUS "Kripke: Using CAMP in Kripke TPL")
        add_subdirectory(${KRIPKE_SOURCE_DIR}/tpl/raja/tpl/camp
                         ${PROJECT_BINARY_DIR}/camp)
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
      set(UMPIRE_DEPENDS "${UMPIRE_DEPENDS} blt::mpi")
    endif()

    blt_register_library(
      NAME umpire
      INCLUDES ${UMPIRE_INCLUDE_DIRS}
      LIBRARIES umpire
      DEPENDS_ON ${UMPIRE_DEPENDS})
  else ()
    add_subdirectory(${KRIPKE_SOURCE_DIR}/tpl/umpire
                     ${PROJECT_BINARY_DIR}/umpire)
  endif()

  if (DEFINED chai_DIR)
    find_package(chai REQUIRED)
  else ()
    add_subdirectory(${KRIPKE_SOURCE_DIR}/tpl/chai
                     ${PROJECT_BINARY_DIR}/chai)
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
    add_subdirectory(${KRIPKE_SOURCE_DIR}/tpl/raja
                     ${PROJECT_BINARY_DIR}/raja)
  endif ()
  list(APPEND KRIPKE_DEPENDS RAJA)
endif ()

#
# Configure OpenMP (CPU threading, not target offload)
#

if(ENABLE_OPENMP)
  
  set(KRIPKE_USE_OPENMP 1)  
  
  list(APPEND KRIPKE_DEPENDS blt::openmp)

  set(KRIPKE_ARCH "OpenMP")
  
endif()



#
# Configure CUDA
#

if(ENABLE_CUDA)

  set(KRIPKE_USE_CUDA 1)
    
  list(APPEND KRIPKE_DEPENDS blt::cuda)
  
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
    
  list(APPEND KRIPKE_DEPENDS blt::hip)
  
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
  find_package(adiak REQUIRED)
  find_package(caliper QUIET)

  if(NOT caliper_FOUND)
    set(BUILD_SHARED_LIBS False CACHE BOOL "")
    set(WITH_MPI TRUE CACHE BOOL "")
    set(WITH_ADIAK TRUE CACHE BOOL "")
    if(ENABLE_CUDA)
      set(WITH_NVTX TRUE CACHE BOOL "")
      set(WITH_CUPTI TRUE CACHE BOOL "")
    elseif(ENABLE_HIP)
      set(WITH_ROCM TRUE CACHE BOOL "")
    endif()
    add_subdirectory(${KRIPKE_SOURCE_DIR}/tpl/caliper
                     ${PROJECT_BINARY_DIR}/caliper)
    if(NOT TARGET Threads::Threads)
      find_package(Threads REQUIRED)
    endif()
    target_link_libraries(caliper PUBLIC Threads::Threads)
    include_directories(${KRIPKE_SOURCE_DIR}/tpl/caliper/include)
  endif()
  list(APPEND KRIPKE_DEPENDS caliper adiak::adiak)
  set(KRIPKE_USE_CALIPER TRUE)
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
	list(APPEND KRIPKE_DEPENDS blt::mpi)
endif()

if(ENABLE_MPI_WRAPPER)
	set(KRIPKE_USE_MPI 1)
endif()

if(ENABLE_GPU_AWARE_MPI)
  if(NOT (ENABLE_MPI OR ENABLE_MPI_WRAPPER))
    message(FATAL_ERROR "ENABLE_GPU_AWARE_MPI requires MPI. Enable either ENABLE_MPI=On or ENABLE_MPI_WRAPPER=On.")
  endif()
  if(NOT ENABLE_CHAI)
    message(FATAL_ERROR "ENABLE_GPU_AWARE_MPI requires ENABLE_CHAI=On.")
  endif()
  if(NOT ENABLE_CUDA AND NOT ENABLE_HIP)
    message(FATAL_ERROR "ENABLE_GPU_AWARE_MPI requires ENABLE_CUDA=On or ENABLE_HIP=On.")
  endif()
  set(KRIPKE_USE_GPU_AWARE_MPI 1)
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

configure_file(${KRIPKE_SOURCE_DIR}/src/KripkeConfig.h.in
               ${PROJECT_BINARY_DIR}/include/KripkeConfig.h)






#
# Define build targets for Kripke
#
 


blt_add_library(
  NAME        kripke
  SOURCES     "${KRIPKE_SOURCE_DIR}/src/Kripke/Core/BaseVar.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Core/MemoryManager.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Core/DataStore.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Core/DomainVar.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Generate.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Generate/Data.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Generate/Decomp.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Generate/Energy.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Generate/Quadrature.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Generate/Space.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/InputVariables.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Kernel/LPlusTimes.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Kernel/LTimes.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Kernel/Population.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Kernel/Scattering.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Kernel/Source.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Kernel/SweepSubdomain.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/ParallelComm/BlockJacobiComm.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/ParallelComm/SweepComm.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/ParallelComm.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Core/PartitionSpace.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Core/Set.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/SteadyStateSolver.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/SweepSolver.cpp"
              "${KRIPKE_SOURCE_DIR}/src/Kripke/Timing.cpp"
  DEPENDS_ON  ${KRIPKE_DEPENDS}
)


target_include_directories(kripke PUBLIC 
	$<BUILD_INTERFACE:${KRIPKE_SOURCE_DIR}/src>
	$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>)

blt_add_executable( 
  NAME        kripke.exe
  SOURCES     "${KRIPKE_SOURCE_DIR}/src/kripke.cpp"
  DEPENDS_ON  ${KRIPKE_DEPENDS} kripke 
)

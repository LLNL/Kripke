cmake_minimum_required(VERSION 3.1)

project(KRIPKE LANGUAGES CXX)

#
# Initialize the BLT build system
#

set(ENABLE_GTEST Off CACHE Bool "")
set(ENABLE_TESTS Off CACHE Bool "")
set(ENABLE_TBB Off CACHE Bool "")
set(ENABLE_EXAMPLES Off CACHE Bool "")
set(ENABLE_DOCUMENTATION Off CACHE Bool "")

include(blt/SetupBLT.cmake)



#
# Configure RAJA
#

set(RAJA_ENABLE_TESTS Off CACHE Bool "")
add_subdirectory(tpl/raja)# EXCLUDE_FROM_ALL)

get_property(RAJA_INCLUDE_DIRS TARGET RAJA PROPERTY INCLUDE_DIRECTORIES)
include_directories(${RAJA_INCLUDE_DIRS})




#
# Figure out what Arch we are selecting
#


set(${KRIPKE_ARCH_FOUND} 0)
if(NOT KRIPKE_ARCH)
  set(KRIPKE_ARCH "OpenMP")
endif()

if(${KRIPKE_ARCH} STREQUAL "Sequential")
  set(${KRIPKE_ARCH_FOUND} 1)
  set(KRIPKE_ARCH_SEQUENTIAL 1)
endif()


if(${KRIPKE_ARCH} STREQUAL "OpenMP")
  set(${KRIPKE_ARCH_FOUND} 1)
  set(KRIPKE_ARCH_OPENMP 1)
endif()


if(NOT ${KRIPKE_ARCH_FOUND})
  error("I don't know what architecture this is:'${KRIPKE_ARCH}'")
endif()


message(STATUS "Kripke selected architecture '${KRIPKE_ARCH}'")


#
# Setup variables to pass to KripkeConfig.h
# 



if(ENABLE_MPI)
	set(KRIPKE_USE_MPI 1)
	set(KRIPKE_DEPENDS ${KRIPKE_DEPENDS} mpi)
endif()

if(ENABLE_OPENMP)
	set(KRIPKE_USE_OPENMP 1)
	set(KRIPKE_DEPENDS ${KRIPKE_DEPENDS} openmp)
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} KRIPKE_BUILD_TYPE)
set(KRIPKE_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
set(KRIPKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${KRIPKE_BUILD_TYPE}}")
set(KRIPKE_LINK_FLAGS   "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_LINK_FLAGS_${CMAKE_BUILD_TYPE}}")

set(KRIPKE_VERSION "1.2.0")

configure_file(${PROJECT_SOURCE_DIR}/src/KripkeConfig.h.in
               ${PROJECT_BINARY_DIR}/include/KripkeConfig.h)




#
# Define build targets for Kripke
#
 

include_directories(./src ${PROJECT_BINARY_DIR}/include)

blt_add_library(
  NAME        kripke
  SOURCES     "src/Kripke/Core/DataStore.cpp"
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

add_dependencies(kripke RAJA)

blt_add_executable( 
  NAME        kripke.exe
  SOURCES     "src/kripke.cpp"
  DEPENDS_ON  ${KRIPKE_DEPENDS} kripke
)



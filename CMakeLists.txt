cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0007 NEW)
cmake_policy(SET CMP0057 NEW)

project(KRIPKE LANGUAGES CXX)



#
# Initialize the BLT build system
#

set(ENABLE_GTEST Off CACHE Bool "")
set(ENABLE_TESTS Off CACHE Bool "")
set(ENABLE_TBB Off CACHE Bool "")
set(ENABLE_EXAMPLES Off CACHE Bool "")
set(ENABLE_DOCUMENTATION Off CACHE Bool "")

set(BLT_CXX_STD c++14 CACHE STRING "")

set(RAJA_CXX_STANDARD_FLAG ${BLT_CXX_STD} CACHE STRING "")

include(blt/SetupBLT.cmake)


#
# Configure RAJA
#

set(RAJA_ENABLE_TESTS Off CACHE Bool "")
add_subdirectory(tpl/raja)

get_property(RAJA_INCLUDE_DIRS TARGET RAJA PROPERTY INCLUDE_DIRECTORIES)
include_directories(${RAJA_INCLUDE_DIRS})


#
# Configure CHAI/Umpire if CUDA was enabled
#

if(ENABLE_CUDA)
	add_subdirectory(tpl/chai)

	set(ENABLE_TESTS Off CACHE Bool "")
	set(KRIPKE_DEPENDS ${KRIPKE_DEPENDS} chai umpire)

	get_property(CHAI_INCLUDE_DIRS TARGET chai PROPERTY INCLUDE_DIRECTORIES)
	include_directories(${CHAI_INCLUDE_DIRS})

	get_property(UMPIRE_INCLUDE_DIRS TARGET umpire PROPERTY INCLUDE_DIRECTORIES)
	include_directories(${UMPIRE_INCLUDE_DIRS})
	
	if(ENABLE_OPENMP)
	 list(APPEND CUDA_NVCC_FLAGS -Xcompiler ${OpenMP_CXX_FLAGS})
	endif()

endif()








#
# Detect capabilities for Kripke to use
# 

set(KRIPKE_ARCH_AUTO "Sequential")

if(ENABLE_MPI)
	set(KRIPKE_USE_MPI 1)
	set(KRIPKE_DEPENDS ${KRIPKE_DEPENDS} mpi)	
endif()

if(ENABLE_MPI_WRAPPER)
	set(KRIPKE_USE_MPI 1)
endif()


if(ENABLE_OPENMP)
	set(KRIPKE_USE_OPENMP 1)
	set(KRIPKE_DEPENDS ${KRIPKE_DEPENDS} openmp)
	set(KRIPKE_ARCH_AUTO "OpenMP")
endif()

if(ENABLE_CUDA)
  set(KRIPKE_USE_CUDA 1)
  set(KRIPKE_DEPENDS ${KRIPKE_DEPENDS} cuda  chai)
  set(KRIPKE_ARCH_AUTO "CUDA")

endif()



#
# Setup default Arch (specified by user, or autodetected above)
# 



if(KRIPKE_ARCH)

	set(KRIPKE_ARCH_FOUND 0)

	if(${KRIPKE_ARCH} STREQUAL "Sequential")
		set(KRIPKE_ARCH_FOUND 1)
		set(KRIPKE_ARCH_SEQUENTIAL 1)
	endif()


	if(${KRIPKE_ARCH} STREQUAL "OpenMP")
		set(KRIPKE_ARCH_FOUND 1)
		set(KRIPKE_ARCH_OPENMP 1)
	endif()


	if(NOT ${KRIPKE_ARCH_FOUND})
		error("I don't know what architecture this is:'${KRIPKE_ARCH}'")
	endif()

else()
	set(KRIPKE_ARCH_FOUND 1)
	set(KRIPKE_ARCH ${KRIPKE_ARCH_AUTO})
endif()



if(KRIPKE_LAYOUT)
else()
	set(KRIPKE_LAYOUT DGZ)
endif()


message(STATUS "Kripke selected default architecture: '${KRIPKE_ARCH}'")
message(STATUS "Kripke selected default layout:       '${KRIPKE_LAYOUT}'")


set(KRIPKE_ARCHV_DEFAULT ArchV_${KRIPKE_ARCH})
set(KRIPKE_LAYOUTV_DEFAULT LayoutV_${KRIPKE_LAYOUT})



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

message(INFO, "KRIPKE_DEPENDS: ${KRIPKE_DEPENDS}")


blt_add_executable( 
  NAME        kripke.exe
  SOURCES     "src/kripke.cpp"
  DEPENDS_ON  ${KRIPKE_DEPENDS} kripke
)



cmake_minimum_required(VERSION 3.1)

project(KRIPKE LANGUAGES CXX)


#
# Configure RAJA
#

set(RAJA_ENABLE_TESTS Off)
set(RAJA_ENABLE_EXAMPLES Off)
set(RAJA_ENABLE_DOCUMENTATION Off)
add_subdirectory(tpl/raja EXCLUDE_FROM_ALL)

get_property(RAJA_INCLUDE_DIRS DIRECTORY tpl/raja PROPERTY INCLUDE_DIRECTORIES)
include_directories(${RAJA_INCLUDE_DIRS})



#
# Initialize the BLT build system
#


include(blt/SetupBLT.cmake)




#
# Setup variables to pass to KripkeConfig.h
# 

if(ENABLE_MPI)
	set(KRIPKE_USE_MPI 1)
endif()

if(ENABLE_OPENMP)
	set(KRIPKE_USE_OPENMP 1)
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
  SOURCES     "src/Kripke/Directions.cpp"
              "src/Kripke/Grid.cpp"
              "src/Kripke/Input_Variables.cpp"
              "src/Kripke/Kernel/Kernel_3d_DGZ.cpp"
              "src/Kripke/Kernel/Kernel_3d_DZG.cpp"
              "src/Kripke/Kernel/Kernel_3d_GDZ.cpp"
              "src/Kripke/Kernel/Kernel_3d_GZD.cpp"
              "src/Kripke/Kernel/Kernel_3d_ZDG.cpp"
              "src/Kripke/Kernel/Kernel_3d_ZGD.cpp"
              "src/Kripke/Kernel.cpp"
              "src/Kripke/Layout.cpp"
              "src/Kripke/ParallelComm/BlockJacobiComm.cpp"
              "src/Kripke/ParallelComm/SweepComm.cpp"
              "src/Kripke/ParallelComm.cpp"
              "src/Kripke/Subdomain.cpp"
              "src/Kripke/Sweep_Solver.cpp"
              "src/Kripke/Test/TestKernels.cpp"
              "src/Kripke/Timing.cpp"
  DEPENDS_ON  MPI 
)

add_dependencies(kripke RAJA)

blt_add_executable( 
  NAME        kripke.exe
  SOURCES     "src/kripke.cpp"
  DEPENDS_ON  kripke MPI 
)



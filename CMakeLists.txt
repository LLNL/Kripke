cmake_minimum_required(VERSION 3.0)

project(KRIPKE)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})
message(STATUS ${CMAKE_MODULE_PATH})


# MPI Support (Required)
# If the ENABLE_MPI is turned off, you must use the mpi* compiler wrappers, since
# Kripke will not build without MPI support
set(ENABLE_MPI On CACHE STRING "Turn this off if you are using mpi compiler wrappers")
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  set(KRIPKE_LIBS ${KRIPKE_LIBS} ${MPI_LIBRARIES})
endif()


# OpenMP Support (Optional)
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
    add_definitions (-DKRIPKE_USE_OPENMP)
    message(STATUS "OpenMP Enabled")
  else()
    message(WARNING "OpenMP NOT FOUND")
  endif()
endif()


# LLNL's Silo file support (Optional)
if(ENABLE_SILO)
  find_package(Silo)
  if(Silo_FOUND)
    add_definitions (-DKRIPKE_USE_SILO)
    include_directories(${Silo_INCLUDE_DIRS})
    set(KRIPKE_LIBS ${KRIPKE_LIBS} ${Silo_LIBRARIES})
    message(STATUS "Silo Enabled")
  else()
    message(WARNING "Silo NOT FOUND")
  endif()
endif()


# PAPI Support for Hardware Counters (Optional)
if(ENABLE_PAPI)
  find_package(PAPI)
  if(PAPI_FOUND)
    add_definitions(-DKRIPKE_USE_PAPI)
    include_directories(${PAPI_INCLUDE_DIRS})
    set(KRIPKE_LIBS ${KRIPKE_LIBS} ${PAPI_LIBRARIES})
    message(STATUS "PAPI Enabled")
  else()
    message(WARNING "PAPI NOT FOUND")
  endif()
endif()


# TCMalloc Support for better malloc/profiling (Optional)
if(ENABLE_TCMALLOC)
  find_package(TCMalloc)
  if(TCMalloc_FOUND)
    add_definitions(-DKRIPKE_USE_TCMALLOC)
    include_directories(${TCMalloc_INCLUDE_DIRS})
    set(KRIPKE_LIBS ${KRIPKE_LIBS} ${TCMalloc_LIBRARIES})
    message(STATUS "TCMalloc Enabled")
  else()
    message(WARNING "TCMalloc NOT FOUND")
  endif()
endif()


# BG/Q Performance tools (BGPM) Support (Optional, only on BG/Q)
if(ENABLED_BGPM)
  find_package(TCMalloc)
  if(TCMalloc_FOUND)
    add_definitions(-DKRIPKE_USE_TCMALLOC)
    include_directories(${TCMalloc_INCLUDE_DIRS})
    set(KRIPKE_LIBS ${KRIPKE_LIBS} ${TCMalloc_LIBRARIES})
    message(STATUS "TCMalloc Enabled")
  else()
    message(WARNING "TCMalloc NOT FOUND")
  endif()
endif()


# Tarball Creation Support
if(ENABLE_CPACK)
  set(CPACK_PACKAGE_VERSION_MAJOR "1")
  set(CPACK_PACKAGE_VERSION_MINOR "1")
  set(CPACK_PACKAGE_VERSION_PATCH "0")
  set(CPACK_SOURCE_GENERATOR "TGZ")
  set(CPACK_SOURCE_IGNORE_FILES "tarball.py;kripke-tarball;.git;.cproject;.project;${CPACK_SOURCE_IGNORE_FILES}")
  include(CPack)
endif()


# Optional 
#bdiv_opt_pkg(perftools 0 OFF)


#
# Other Build Options
# 

# Traverse subdirectories
include_directories(src)
add_subdirectory(src)

# Setup the kripke target
add_executable(kripke "src/kripke.cpp")
target_link_libraries(kripke ${KRIPKE_LIBS} ${KRIPKE_LIBS})


#
# Testing suite
#

enable_testing()

# Use testing harness to check all kernels, all nestings, with default problem
add_test(test_default_dgz kripke --nest dgz --test)
add_test(test_default_dzg kripke --nest dzg --test)
add_test(test_default_gdz kripke --nest gdz --test)
add_test(test_default_gzd kripke --nest gzd --test)
add_test(test_default_zdg kripke --nest zdg --test)
add_test(test_default_zgd kripke --nest zgd --test)

# Use testing harness to check use of a zone sets for all nestings
add_test(test_zset_dgz_block kripke --nest dgz --test --zset 2,3,5 --gset 3 --groups 6 --layout 0)
add_test(test_zset_dzg_block kripke --nest dzg --test --zset 2,3,5 --gset 3 --groups 6 --layout 0)
add_test(test_zset_gdz_block kripke --nest gdz --test --zset 2,3,5 --gset 3 --groups 6 --layout 0)
add_test(test_zset_gzd_block kripke --nest gzd --test --zset 2,3,5 --gset 3 --groups 6 --layout 0)
add_test(test_zset_zdg_block kripke --nest zdg --test --zset 2,3,5 --gset 3 --groups 6 --layout 0)
add_test(test_zset_zgd_block kripke --nest zgd --test --zset 2,3,5 --gset 3 --groups 6 --layout 0)

# Also check with scattered layout (just checking parallel algo, so just 1 nesting)
add_test(test_zset_dgz_scatter kripke --nest dgz --zset 2,3,5 --gset 3 --groups 6 --layout 1)

# Test the block-Jacobi parallel method
add_test(test_bj_dgz kripke --nest dgz --pmethod bj)




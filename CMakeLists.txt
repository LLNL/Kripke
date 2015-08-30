cmake_minimum_required(VERSION 3.0)

project(KRIPKE)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})
message(STATUS ${CMAKE_MODULE_PATH})


# MPI Support (Required)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(KRIPKE_LIBS ${KRIPKE_LIBS} ${MPI_LIBRARIES})


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

# BG/Q Specific
#bdiv_opt_pkg(BGPM 20030410 OFF)


#
# Other Build Options
# 

# Traverse subdirectories
include_directories(src)
add_subdirectory(src)

# Setup the kripke target
add_executable(kripke "src/kripke.cpp")
target_link_libraries(kripke ${KRIPKE_LIBS} ${KRIPKE_LIBS})




#
# Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)
#

if (ENABLE_HIP)
  cmake_minimum_required(VERSION 3.23)
else()
  cmake_minimum_required(VERSION 3.20)
endif()

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

# Use C++17 standard
set(BLT_CXX_STD "c++17" CACHE STRING "")

if (NOT BLT_LOADED)
  if (DEFINED BLT_SOURCE_DIR)
    message(${BLT_SOURCE_DIR})
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
      message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
  else ()
    set (BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/blt" CACHE PATH "")

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
#add_subdirectory(tpl/googletest)


#
# Configure CHAI/Umpire for memory management
# 
option(ENABLE_CHAI "Enable CHAI/Umpire memory management" Off)
option(ENABLE_CHAI_SINGLE_MEMORY "Enable single memory configuration in CHAI" Off)
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
    if(ENABLE_CHAI_SINGLE_MEMORY)
        # Set CHAI-related HIP variables for single memory address space.
        # Use the THIN_GPU_ALLOCATE mode of CHAI on MI250X.
        # This is to mirror the expected behavior of El Capitan.
        # Requires env variables HSA_XNACK set to 1 and MPICH_GPU_SUPPORT_ENABLED set to 1 to work on the MI250X.
        set (CHAI_ENABLE_UM OFF CACHE BOOL "")
        set (CHAI_ENABLE_PINNED ON CACHE BOOL "")
        set(CHAI_DISABLE_RM ON CACHE BOOL "")
        set(CHAI_THIN_GPU_ALLOCATE ON CACHE BOOL "")

        set(KRIPKE_USE_CHAI_SINGLE_MEMORY 1)
    endif()

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
  find_package(adiak QUIET)
  if (adiak_FOUND)
    # First, build Adiak with:
    # cd tpl/adiak/ ; mkdir build
    # cmake -DBUILD_SHARED_LIBS=Off -DENABLE_MPI=On -DMPI_C_COMPILER=mpiamdclang -DMPI_CXX_COMPILER=mpiamdclang++ -DCMAKE_C_COMPILER=mpiamdclang -DCMAKE_CXX_COMPILER=mpiamdclang++ -DCMAKE_INSTALL_PREFIX=./install ..
    # make -j
    # Then, set adiak_DIR=tpl/adiak/build/install/lib/cmake/adiak/ when building Kripke.
    list(APPEND KRIPKE_DEPENDS adiak::adiak)
  else ()
    # Stage a build-local Adiak copy and patch its MPI export handling without
    # modifying the vendored source tree.
    set(_kripke_adiak_staging_root "${PROJECT_BINARY_DIR}/tpl")
    set(_kripke_adiak_staging_dir "${_kripke_adiak_staging_root}/adiak")
    set(_kripke_adiak_build_dir "${PROJECT_BINARY_DIR}/tpl/adiak-build")
    set(_kripke_adiak_src_cmake "${_kripke_adiak_staging_dir}/src/CMakeLists.txt")
    set(_kripke_adiak_src_cmake_patched "${PROJECT_BINARY_DIR}/adiak-src-CMakeLists.txt")
    set(_kripke_adiak_patch_file "${PROJECT_BINARY_DIR}/adiak-mpi-export.patch")

    file(REMOVE_RECURSE "${_kripke_adiak_staging_dir}" "${_kripke_adiak_build_dir}")
    file(MAKE_DIRECTORY "${_kripke_adiak_staging_root}")
    file(COPY "${PROJECT_SOURCE_DIR}/tpl/adiak" DESTINATION "${_kripke_adiak_staging_root}")

    file(READ "${_kripke_adiak_src_cmake}" _kripke_adiak_src_contents)
    set(_kripke_adiak_src_original "${_kripke_adiak_src_contents}")

    string(REPLACE
      "set(adiak_public_headers"
      "if (MPI_FOUND AND TARGET mpi)\n  get_target_property(_adiak_mpi_imported mpi IMPORTED)\n  if (_adiak_mpi_imported AND NOT TARGET adiak_mpi_export)\n    add_library(adiak_mpi_export INTERFACE)\n    target_link_libraries(adiak_mpi_export INTERFACE mpi)\n  endif ()\nendif ()\n\nset(adiak_public_headers"
      _kripke_adiak_src_contents
      "${_kripke_adiak_src_contents}")

    string(REPLACE
      "  list(APPEND adiak_dependencies\n    mpi)\n  list(APPEND adiak_export_targets\n    mpi)"
      "  if (TARGET adiak_mpi_export)\n    list(APPEND adiak_dependencies\n      adiak_mpi_export)\n    list(APPEND adiak_export_targets\n      adiak_mpi_export)\n  else ()\n    list(APPEND adiak_dependencies\n      mpi)\n    list(APPEND adiak_export_targets\n      mpi)\n  endif()"
      _kripke_adiak_src_contents
      "${_kripke_adiak_src_contents}")

    if (_kripke_adiak_src_contents STREQUAL _kripke_adiak_src_original)
      message(FATAL_ERROR "Kripke: failed to prepare staged Adiak MPI export patch content")
    endif ()

    file(WRITE "${_kripke_adiak_src_cmake_patched}" "${_kripke_adiak_src_contents}")

    execute_process(
      COMMAND diff -u "src/CMakeLists.txt" "${_kripke_adiak_src_cmake_patched}"
      WORKING_DIRECTORY "${_kripke_adiak_staging_dir}"
      RESULT_VARIABLE _kripke_adiak_diff_result
      OUTPUT_FILE "${_kripke_adiak_patch_file}"
      ERROR_VARIABLE _kripke_adiak_diff_stderr
    )

    if (_kripke_adiak_diff_result LESS 0 OR _kripke_adiak_diff_result GREATER 1)
      message(FATAL_ERROR
        "Kripke: failed to generate staged Adiak MPI export patch\n"
        "stderr:\n${_kripke_adiak_diff_stderr}")
    endif ()

    execute_process(
      COMMAND patch -p0 -i "${_kripke_adiak_patch_file}"
      WORKING_DIRECTORY "${_kripke_adiak_staging_dir}"
      RESULT_VARIABLE _kripke_adiak_patch_result
      OUTPUT_VARIABLE _kripke_adiak_patch_stdout
      ERROR_VARIABLE _kripke_adiak_patch_stderr
    )

    if (NOT _kripke_adiak_patch_result EQUAL 0)
      message(FATAL_ERROR
        "Kripke: failed to patch staged Adiak MPI export handling\n"
        "stdout:\n${_kripke_adiak_patch_stdout}\n"
        "stderr:\n${_kripke_adiak_patch_stderr}")
    endif ()

		set(BUILD_SHARED_LIBS False CACHE BOOL "")
    add_subdirectory("${_kripke_adiak_staging_dir}" "${_kripke_adiak_build_dir}")
		include_directories("${_kripke_adiak_staging_dir}/include")
    set_target_properties(adiak PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${_kripke_adiak_build_dir}/lib"
      LIBRARY_OUTPUT_DIRECTORY "${_kripke_adiak_build_dir}/lib")
    set(_kripke_adiak_library
      "${_kripke_adiak_build_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}adiak${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(_kripke_adiak_pkg_dir "${_kripke_adiak_build_dir}/lib/cmake/adiak")
    file(MAKE_DIRECTORY "${_kripke_adiak_pkg_dir}")
    set(_kripke_adiak_config_template [=[
if (NOT TARGET adiak::adiak)
  add_library(adiak::adiak STATIC IMPORTED GLOBAL)
  set_target_properties(adiak::adiak PROPERTIES
    IMPORTED_LOCATION "@_kripke_adiak_library@"
    INTERFACE_INCLUDE_DIRECTORIES "@_kripke_adiak_staging_dir@/include")
  if (TARGET adiak)
    add_dependencies(adiak::adiak adiak)
  endif()
endif()
set(adiak_FOUND TRUE)
set(adiak_INCLUDE_DIR "@_kripke_adiak_staging_dir@/include")
set(adiak_INCLUDE_DIRS "@_kripke_adiak_staging_dir@/include")
set(adiak_LIBRARIES adiak::adiak)
]=])
    string(CONFIGURE "${_kripke_adiak_config_template}" _kripke_adiak_config @ONLY)
    file(WRITE "${_kripke_adiak_pkg_dir}/adiak-config.cmake" "${_kripke_adiak_config}")
    set(adiak_DIR "${_kripke_adiak_pkg_dir}" CACHE PATH
      "Build-tree Adiak package config for Caliper" FORCE)
    list(APPEND KRIPKE_DEPENDS adiak)

    unset(_kripke_adiak_staging_root)
    unset(_kripke_adiak_staging_dir)
    unset(_kripke_adiak_build_dir)
    unset(_kripke_adiak_config_template)
    unset(_kripke_adiak_config)
    unset(_kripke_adiak_library)
    unset(_kripke_adiak_pkg_dir)
    unset(_kripke_adiak_src_cmake)
    unset(_kripke_adiak_src_cmake_patched)
    unset(_kripke_adiak_src_contents)
    unset(_kripke_adiak_src_original)
    unset(_kripke_adiak_patch_file)
    unset(_kripke_adiak_diff_result)
    unset(_kripke_adiak_diff_stderr)
    unset(_kripke_adiak_patch_result)
    unset(_kripke_adiak_patch_stdout)
    unset(_kripke_adiak_patch_stderr)
  endif()

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
    add_subdirectory(tpl/caliper)	
    if (NOT TARGET Threads::Threads)
      find_package(Threads REQUIRED)
    endif()
    target_link_libraries(caliper PUBLIC Threads::Threads)
		include_directories(tpl/caliper/include)
 	endif()
  list(APPEND KRIPKE_DEPENDS caliper)

#    blt_register_library(NAME caliper
#      INCLUDES ${caliper_INCLUDE_DIR}
#      LIBRARIES caliper)

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
              "src/Kripke/Core/MemoryManager.cpp"
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

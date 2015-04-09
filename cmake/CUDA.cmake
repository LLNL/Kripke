#
# Cuda
# 
option(ENABLE_CUDA "Turn on compiler support for CUDA" ON)
message(STATUS "CUDA Support is ${ENABLE_CUDA}")

if(${ENABLE_CUDA})
  find_package(CUDA)
  if (CUDA_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUDA_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CUDA_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CUDA_EXE_LINKER_FLAGS}")
      add_definitions (-DKRIPKE_USE_CUDA)
  endif()
endif()


add_subdirectory(Kripke)

list(APPEND KRIPKE_LIBS ${KRIPKE_LIBS} PARENT_SCOPE)

#if(ENABLE_CUDA)
#	cuda_add_library(cukripkelib STATIC ${CUDA_KRIPKE_SOURCE})

#  list(APPEND KRIPKE_LIBS cukripkelib PARENT_SCOPE)
#  list(APPEND KRIPKE_LIBS ${CUDA_cublas_LIBRARY} ${CUDA_npps_LIBRARY} ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libnvToolsExt.so PARENT_SCOPE)

#endif()



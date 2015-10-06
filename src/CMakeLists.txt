include_directories(.)

add_subdirectory(Kripke)

add_library(kripkelib STATIC ${KRIPKE_SOURCE})

#add_subdirectory(Kripke)

if(ENABLE_CUDA)
	cuda_add_library(cukripkelib STATIC ${CUDA_KRIPKE_SOURCE})

	set(KRIPKE_LIBS ${KRIPKE_LIBS} cukripkelib)
        set(KRIPKE_LIBS ${KRIPKE_LIBS} /usr/local/cuda/lib64/libcublas.so)

endif()

add_subdirectory(tools)

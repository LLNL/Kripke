include_directories(.)

add_subdirectory(Kripke)

add_library(kripkelib STATIC ${KRIPKE_SOURCE})

#add_subdirectory(Kripke)


#cuda_add_library(cukripkelib STATIC Kripke/test.cu)
cuda_add_library(cukripkelib STATIC ${CUDA_KRIPKE_SOURCE})

set(KRIPKE_LIBS ${KRIPKE_LIBS} cukripkelib)

add_subdirectory(tools)

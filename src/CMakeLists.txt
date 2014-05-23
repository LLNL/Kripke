


#
# OpenMP
#





include_directories(.)

add_subdirectory(Kripke)
#add_subdirectory(test)

add_library(kripkelib STATIC ${KRIPKE_SOURCE})

add_subdirectory(tools)


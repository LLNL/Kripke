include_directories(.)

add_subdirectory(Kripke)

add_library(kripke STATIC ${KRIPKE_SOURCE})

add_subdirectory(tools)



#
# OpenMP
#





include_directories(.)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

add_subdirectory(Kripke)
#add_subdirectory(test)

add_library(kripkelib STATIC ${KRIPKE_SOURCE})

add_subdirectory(tools)


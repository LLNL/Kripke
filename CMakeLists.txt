cmake_minimum_required(VERSION 2.8)

project(KRIPKE)

enable_language(CXX)
#set(CMAKE_BUILD_TYPE Debug)

set(KRIPKE_BASE_LIBS kripkelib)
set(KRIPKE_LIBS ${KRIPKE_BASE_LIBS})
include_directories(${INCLUDE_PATHS})

add_subdirectory(src)


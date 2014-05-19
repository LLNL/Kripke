cmake_minimum_required(VERSION 2.8)

project(KRIPKE)

enable_language(CXX)

set(KRIPKE_BASE_LIBS kripke)
set(KRIPKE_LIBS ${KRIPKE_BASE_LIBS})
include_directories(${INCLUDE_PATHS})

add_subdirectory(src)


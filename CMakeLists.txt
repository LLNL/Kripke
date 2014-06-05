cmake_minimum_required(VERSION 2.8)

project(KRIPKE)

enable_language(CXX)

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "tarball.py;kripke-tarball;.git;.cproject;.project;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

#set(CMAKE_BUILD_TYPE Debug)

set(KRIPKE_BASE_LIBS kripkelib)
set(KRIPKE_LIBS ${KRIPKE_BASE_LIBS})
include_directories(${INCLUDE_PATHS})

add_subdirectory(src)


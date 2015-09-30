include_directories(.)


add_subdirectory(Kripke)
add_subdirectory(RAJA)

set(KRIPKE_LIBS ${KRIPKE_LIBS} PARENT_SCOPE)

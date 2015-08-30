include_directories(.)

add_subdirectory(Kripke)

set(KRIPKE_LIBS ${KRIPKE_LIBS} PARENT_SCOPE)

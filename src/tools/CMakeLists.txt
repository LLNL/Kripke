
include_directories(.)

add_executable(kripkeSweep "kripkeSweep.cpp")
target_link_libraries(kripkeSweep ${KRIPKE_LIBS} ${KRIPKE_LIBS})

add_executable(kripkeSearch "kripkeSearch.cpp")
target_link_libraries(kripkeSearch ${KRIPKE_LIBS} ${KRIPKE_LIBS})



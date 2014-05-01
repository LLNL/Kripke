
include_directories(.)

add_executable(sweep_kernel "sweep_kernel.cpp")
target_link_libraries(sweep_kernel ${KRIPKE_LIBS} ${KRIPKE_LIBS})



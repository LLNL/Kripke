include_directories(.)


add_definitions (-DRAJA_USE_DOUBLE -DRAJA_USE_RESTRICT_PTR -DRAJA_PLATFORM_X86_SSE -DRAJA_COMPILER_GNU)
include_directories(./RAJA/includes)

add_subdirectory(Kripke)
add_subdirectory(RAJA)

set(KRIPKE_LIBS ${KRIPKE_LIBS} PARENT_SCOPE)

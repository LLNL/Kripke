set(ICC_VER 14.0.174)

set(CMAKE_C_COMPILER mpiicc-${ICC_VER})
set(CMAKE_CXX_COMPILER mpiicpc-${ICC_VER})
set(CMAKE_LINKER mpiicpc-${ICC_VER})

set(MPI "-mpi=mvapich2-intel-1.9")

set(CMAKE_C_FLAGS "${MPI} -g -fno-omit-frame-pointer -unroll-aggressive -finline-functions -axAVX -msse4.2 -no-fma ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${MPI} -g -unroll-aggressive -finline-functions -axAVX -msse4.2  ${CMAKE_CXX_FLAGS}")

set(PKG_PATH "/usr/gapps/bdiv/${SYS_TYPE}/icc-14.0-opt-mvapich2-intel-1.9")



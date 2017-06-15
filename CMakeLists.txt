cmake_minimum_required(VERSION 3.1)

project(KRIPKE LANGUAGES CXX)

set(ENABLE_MPI "On")

include(blt/SetupBLT.cmake)

include_directories(./src)

blt_add_library(
  NAME        kripke
  SOURCES     "src/Kripke/Directions.cpp"
              "src/Kripke/Grid.cpp"
              "src/Kripke/Input_Variables.cpp"
              "src/Kripke/Kernel/Kernel_3d_DGZ.cpp"
              "src/Kripke/Kernel/Kernel_3d_DZG.cpp"
              "src/Kripke/Kernel/Kernel_3d_GDZ.cpp"
              "src/Kripke/Kernel/Kernel_3d_GZD.cpp"
              "src/Kripke/Kernel/Kernel_3d_ZDG.cpp"
              "src/Kripke/Kernel/Kernel_3d_ZGD.cpp"
              "src/Kripke/Kernel.cpp"
              "src/Kripke/Layout.cpp"
              "src/Kripke/ParallelComm/BlockJacobiComm.cpp"
              "src/Kripke/ParallelComm/SweepComm.cpp"
              "src/Kripke/ParallelComm.cpp"
              "src/Kripke/Subdomain.cpp"
              "src/Kripke/Sweep_Solver.cpp"
              "src/Kripke/Test/TestKernels.cpp"
              "src/Kripke/Timing.cpp"
  DEPENDS_ON  MPI
)

blt_add_executable( 
  NAME        kripke.exe
  SOURCES     "src/kripke.cpp"
  DEPENDS_ON  kripke MPI
)



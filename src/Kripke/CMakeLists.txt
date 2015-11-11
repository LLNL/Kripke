include_directories(..)

add_library(lib_kripke  
  Directions.cpp  
  Grid.cpp  
  Input_Variables.cpp
  Kernel.cpp
  Layout.cpp
  Subdomain.cpp  
  Sweep_Solver.cpp
  ParallelComm.cpp  
  Timing.cpp
  example.cpp  
    
  Kernel/Kernel_3d_GDZ.cpp
  Kernel/Kernel_3d_DGZ.cpp
  Kernel/Kernel_3d_ZDG.cpp
  Kernel/Kernel_3d_DZG.cpp
  Kernel/Kernel_3d_ZGD.cpp
  Kernel/Kernel_3d_GZD.cpp
  
  ParallelComm/BlockJacobiComm.cpp
  ParallelComm/SweepComm.cpp 
  
  Test/TestKernels.cpp
)

#target_compile_features(lib_kripke PUBLIC cxx_auto_type cxx_attributes cxx_lambdas)
                                                                                                          
set(KRIPKE_LIBS ${KRIPKE_LIBS} lib_kripke PARENT_SCOPE) 

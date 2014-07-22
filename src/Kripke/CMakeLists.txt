
set(KRIPKE_SOURCE ${KRIPKE_SOURCE}
  Kripke/Comm.cpp  
  Kripke/Directions.cpp  
  Kripke/Grid.cpp  
  Kripke/Kernel.cpp  
  Kripke/Sweep_Solver.cpp  
  Kripke/Timing.cpp  
  Kripke/User_Data.cpp
    
  Kripke/Kernel/Kernel_3d_GDZ.cpp
  Kripke/Kernel/Kernel_3d_DGZ.cpp
  Kripke/Kernel/Kernel_3d_ZDG.cpp
  Kripke/Kernel/Kernel_3d_DZG.cpp
  Kripke/Kernel/Kernel_3d_ZGD.cpp
  Kripke/Kernel/Kernel_3d_GZD.cpp
)
                                                                                                          
set(KRIPKE_SOURCE ${KRIPKE_SOURCE} PARENT_SCOPE) 


set(KRIPKE_SOURCE ${KRIPKE_SOURCE}
  Kripke/Directions.cpp  
  Kripke/Grid.cpp  
  Kripke/Kernel.cpp  
  Kripke/Layout.cpp
  Kripke/Subdomain.cpp
  Kripke/SweepComm.cpp
  Kripke/Sweep_Solver.cpp  
  Kripke/Timing.cpp  
 
  Kripke/Kernel/Kernel_3d_GDZ.cpp
  Kripke/Kernel/Kernel_3d_DGZ.cpp
  Kripke/Kernel/Kernel_3d_ZDG.cpp
  Kripke/Kernel/Kernel_3d_DZG.cpp
  Kripke/Kernel/Kernel_3d_ZGD.cpp
  Kripke/Kernel/Kernel_3d_GZD.cpp
)
                                                                                                          
set(KRIPKE_SOURCE ${KRIPKE_SOURCE} PARENT_SCOPE)


if(ENABLE_CUDA)
	set(CUDA_KRIPKE_SOURCE ${CUDA_KRIPKE_SOURCE}
	  Kripke/cu_utils.cu
	  Kripke/sweep_ZGD.cu
	  Kripke/sweep_ZDG.cu 
    Kripke/sweep_GZD.cu
	)

	set(CUDA_KRIPKE_SOURCE ${CUDA_KRIPKE_SOURCE} PARENT_SCOPE)
endif()



 

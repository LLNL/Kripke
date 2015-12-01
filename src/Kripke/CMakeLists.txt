include_directories(..)

add_library(lib_kripke
  Directions.cpp  
  Grid.cpp  
  Kernel.cpp  
  Layout.cpp
  Subdomain.cpp
  SweepComm.cpp
  Sweep_Solver.cpp  
  Timing.cpp  
 
  Kernel/Kernel_3d_GDZ.cpp
  Kernel/Kernel_3d_DGZ.cpp
  Kernel/Kernel_3d_ZDG.cpp
  Kernel/Kernel_3d_DZG.cpp
  Kernel/Kernel_3d_ZGD.cpp
  Kernel/Kernel_3d_GZD.cpp

	../tools/testKernels.cpp
)
                                                                                                          
#list(APPEND KRIPKE_LIBS lib_kripke PARENT_SCOPE) 
set(KRIPKE_LIBS ${KRIPKE_LIBS} lib_kripke PARENT_SCOPE) 


if(ENABLE_CUDA)
	cuda_add_library(cukripkelib
	  cu_utils.cu
	  sweep_ZGD.cu
	  sweep_ZDG.cu 
    sweep_GZD.cu
	)

	list(APPEND KRIPKE_LIBS cukripkelib PARENT_SCOPE) 
	list(APPEND KRIPKE_LIBS ${CUDA_cublas_LIBRARY} ${CUDA_npps_LIBRARY} ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libnvToolsExt.so PARENT_SCOPE)
endif()



 

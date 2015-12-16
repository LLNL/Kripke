include_directories(./includes)

if(ENABLE_CUDA)
	cuda_add_library(RAJA
    sources/AlignedRangeIndexSetBuilders.cxx
    sources/DepGraphNode.cxx
    sources/IndexSet.cxx
    sources/ListSegment.cxx
    sources/LockFreeIndexSetBuilders.cxx
    sources/MemUtils.cxx
    sources/RangeSegment.cxx
  )
else()

  add_library(RAJA  
    sources/AlignedRangeIndexSetBuilders.cxx
    sources/DepGraphNode.cxx
    sources/IndexSet.cxx
    sources/ListSegment.cxx
    sources/LockFreeIndexSetBuilders.cxx
    sources/MemUtils.cxx
    sources/RangeSegment.cxx
  )

endif()

#target_compile_features(RAJA PUBLIC cxx_auto_type cxx_attributes cxx_lambdas)
                                                                                                          
set(KRIPKE_LIBS ${KRIPKE_LIBS} RAJA PARENT_SCOPE) 

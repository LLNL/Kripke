include_directories(./includes)


add_library(RAJA  
  sources/AlignedRangeIndexSetBuilders.cxx
  sources/DepGraphNode.cxx
  sources/IndexSet.cxx
  sources/ListSegment.cxx
  sources/LockFreeIndexSetBuilders.cxx
  sources/MemUtils.cxx
  sources/RangeSegment.cxx
)

                                                                                                          
set(KRIPKE_LIBS ${KRIPKE_LIBS} RAJA PARENT_SCOPE) 

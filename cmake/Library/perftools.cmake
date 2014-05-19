
set(${PACKAGE}_PREFIX "${PKG_PATH}/google-perftools/${VERSION}")

set(perftools_DEFAULT_INCLUDE_DIRS "${perftools_PREFIX}/include")

set(perftools_DEFAULT_LIB_DIR "${perftools_PREFIX}/lib")
set(perftools_DEFAULT_LIBS "profiler")

set(perftools_DEFAULT_DEFINITIONS "-DKRIPKE_USE_PERFTOOLS")





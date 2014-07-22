
set(tcmalloc_PREFIX "${PKG_PATH}/google-perftools/2.1")
set(tcmalloc_DEFAULT_INCLUDE_DIRS "${tcmalloc_PREFIX}/include")
set(tcmalloc_DEFAULT_LIB_DIR "${tcmalloc_PREFIX}/lib")
set(tcmalloc_DEFAULT_LIBS "tcmalloc_minimal")

set(tcmalloc_DEFAULT_DEFINITIONS "-DKRIPKE_USE_TCMALLOC")




set(HDF5_PREFIX "/usr/gapps/silo/hdf5/1.8.10/bgqos_0_bgxlc")
set(HDF5_DEFAULT_INCLUDE_DIRS "${HDF5_PREFIX}/include")
set(HDF5_DEFAULT_LIB_DIR ${HDF5_PREFIX}/lib /usr/gapps/silo/zlib/1.2.3/bgqos_0_bgxlc/lib /usr/gapps/silo/szip/2.1/bgqos_0_bgxlc/lib)
set(HDF5_DEFAULT_LIBS "hdf5")
set(HDF5_EXTRALIBS sz z rt m)

set(HDF5_DEFAULT_DEFINITIONS "-DKRIPKE_USE_HDF5")





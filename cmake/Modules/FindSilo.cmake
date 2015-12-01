# (Slightly adapted from S. Johnson)
# - Find LLNL's Silo library
# This module defines
#   Silo_INCLUDE_DIR, where to find blitz/blitz.h, etc.
#   Silo_LIBRARIES, libraries to link against to use Silo.
#   Silo_FOUND, If false, do not try to use Silo.
# also defined, but not for general use are
#   Silo_LIBRARY, where to find the Silo library.
# The user should specify the head Silo director, Silo_DIR,
# or Silo_INC and Silo_LIB.  

find_path(Silo_INCLUDE_DIRS 
          NAMES silo.h
          PATHS ${Silo_DIR}/include 
                ${Silo_INC}
)

find_library(Silo_LIBRARY
             NAMES siloh5 silo siloxx
             PATHS ${Silo_DIR}/lib
                   ${Silo_LIB}
)

if (Silo_LIBRARY MATCHES "siloh5")
    FIND_PACKAGE(HDF5 REQUIRED)
endif()

# Requires ZLib
#find_package(ZLIB REQUIRED)

SET( Silo_FOUND "NO" )
IF(Silo_INCLUDE_DIRS)
  IF(Silo_LIBRARY)

    SET( Silo_LIBRARIES ${Silo_LIBRARY})
    SET( Silo_FOUND "YES" )

    #The following deprecated settings are for backwards compatibility with CMake1.4
    SET (Silo_INCLUDE_PATH ${Silo_INCLUDE_DIR})

  ENDIF(Silo_LIBRARY)
ENDIF(Silo_INCLUDE_DIRS)

MARK_AS_ADVANCED(
  Silo_INCLUDE_DIR
  Silo_LIBRARY
)





find_path(TCMalloc_INCLUDE_DIRS 
          NAMES gperftools/tcmalloc.h
          PATHS ${TCMalloc_DIR}/include 
                ${TCMalloc_INC}
)

find_library(TCMalloc_LIBRARY
             NAMES tcmalloc 
             PATHS ${TCMalloc_DIR}/lib
                   ${TCMalloc_LIB}
)

SET( TCMalloc_FOUND "NO" )
IF(TCMalloc_INCLUDE_DIRS)
  IF(TCMalloc_LIBRARY)

    SET( TCMalloc_LIBRARIES ${TCMalloc_LIBRARY})
    SET( TCMalloc_FOUND "YES" )

  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  TCMalloc_INCLUDE_DIR
  TCMalloc_LIBRARY
)




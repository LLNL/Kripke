# -*- cmake -*-

# - Find Google perftools
# Find the Google perftools includes and libraries
# This module defines
#  GOOGLE_PERFTOOLS_INCLUDE_DIR, where to find heap-profiler.h, etc.
#  GOOGLE_PERFTOOLS_FOUND, If false, do not try to use Google perftools.
# also defined for general use are
#  TCMALLOC_LIBRARIES, where to find the tcmalloc library.
#  STACKTRACE_LIBRARIES, where to find the stacktrace library.
#  PROFILER_LIBRARIES, where to find the profiler library.

FIND_PATH(GOOGLE_PERFTOOLS_INCLUDE_DIR google/heap-profiler.h
/usr/local/include
/usr/include
)

SET(TCMALLOC_NAMES ${TCMALLOC_NAMES} tcmalloc)
FIND_LIBRARY(TCMALLOC_LIBRARY
  NAMES ${TCMALLOC_NAMES}
  PATHS /usr/lib /usr/local/lib
  )

IF (TCMALLOC_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
    SET(TCMALLOC_LIBRARIES ${TCMALLOC_LIBRARY})
    SET(GOOGLE_PERFTOOLS_FOUND "YES")
ELSE (TCMALLOC_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
  SET(GOOGLE_PERFTOOLS_FOUND "NO")
ENDIF (TCMALLOC_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)

SET(STACKTRACE_NAMES ${STACKTRACE_NAMES} stacktrace)
FIND_LIBRARY(STACKTRACE_LIBRARY
  NAMES ${STACKTRACE_LIBRARY}
  PATHS /usr/lib /usr/local/lib
  )

IF (STACKTRACE_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
    SET(STACKTRACE_LIBRARIES ${STACKTRACE_LIBRARY})
ENDIF (STACKTRACE_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)

SET(PROFILER_NAMES ${PROFILER_NAMES} profiler)
FIND_LIBRARY(PROFILER_LIBRARY
  NAMES ${PROFILER_NAMES}
  PATHS /usr/lib /usr/local/lib
  )

IF (PROFILER_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)
    SET(PROFILER_LIBRARIES ${PROFILER_LIBRARY})
ENDIF (PROFILER_LIBRARY AND GOOGLE_PERFTOOLS_INCLUDE_DIR)

IF (GOOGLE_PERFTOOLS_FOUND)
   IF (NOT GOOGLE_PERFTOOLS_FIND_QUIETLY)
      MESSAGE(STATUS "Found Google perftools: ${GOOGLE_PERFTOOLS_LIBRARIES}")
   ENDIF (NOT GOOGLE_PERFTOOLS_FIND_QUIETLY)
ELSE (GOOGLE_PERFTOOLS_FOUND)
   IF (GOOGLE_PERFTOOLS_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Google perftools library")
   ENDIF (GOOGLE_PERFTOOLS_FIND_REQUIRED)
ENDIF (GOOGLE_PERFTOOLS_FOUND)

MARK_AS_ADVANCED(
  TCMALLOC_LIBRARY
  STACKTRACE_LIBRARY
  PROFILER_LIBRARY
  GOOGLE_PERFTOOLS_INCLUDE_DIR
  )

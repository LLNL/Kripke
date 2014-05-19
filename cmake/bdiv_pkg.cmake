# This provides an OPTIONAL package that lives in the BDIV package path
# But allows for redirecting to a different path

macro (bdiv_opt_pkg PACKAGE VERSION ENABLE_DEFAULT)

  # If enable isn't specified, assume it's enabled by default
  if(NOT ${ENABLE_DEFAULT})
    set(ENABLE_DEFAULT "yes")
  endif()
  
  # Use version from caller, unless specified on command line
  if(NOT DEFINED ${${PACKAGE}_VERSION})
    set(${PACKAGE}_VERSION ${VERSION})
  endif()
  
  # Create command line option, and display a status message
  option(ENABLE_${PACKAGE} "Turn on support for package ${PACKAGE}" ${ENABLE_DEFAULT})  
  
  # If support is requested, then try and find/support it
  if(${ENABLE_${PACKAGE}})
    # Setup the base package path
    set(${PACKAGE}_PREFIX "${PKG_PATH}/${PACKAGE}/${VERSION}")
 
    # Call the module to set the defaults, and define locations
    if(EXISTS "${BDIV_LIBRARY_PATH}/${PACKAGE}-${SYS_TYPE}-${VERSION}.cmake")      
      set(${PACKAGE}_PKG_FILE "${BDIV_LIBRARY_PATH}/${PACKAGE}-${SYS_TYPE}-${VERSION}.cmake")     
    elseif(EXISTS "${BDIV_LIBRARY_PATH}/${PACKAGE}-${SYS_TYPE}.cmake")
      set(${PACKAGE}_PKG_FILE "${BDIV_LIBRARY_PATH}/${PACKAGE}-${SYS_TYPE}.cmake")
    elseif(EXISTS "${BDIV_LIBRARY_PATH}/${PACKAGE}-${VERSION}.cmake")
      set(${PACKAGE}_PKG_FILE "${BDIV_LIBRARY_PATH}/${PACKAGE}-${VERSION}.cmake")
    elseif(EXISTS "${BDIV_LIBRARY_PATH}/${PACKAGE}.cmake")
      set(${PACKAGE}_PKG_FILE "${BDIV_LIBRARY_PATH}/${PACKAGE}.cmake")
    else()
      message(FATAL_ERROR "Couldn't find Library file for ${PACKAGE} in '${BDIV_LIBRARY_PATH}'")
    endif()    
    include(${${PACKAGE}_PKG_FILE})
      
    # Resolve include paths
    if(NOT DEFINED ${PACKAGE}_INCLUDE_DIRS)
      if(DEFINED ${PACKAGE}_DEFAULT_INCLUDE_DIRS)
        set(${PACKAGE}_INCLUDE_DIRS ${${PACKAGE}_DEFAULT_INCLUDE_DIRS})
      endif()
    endif()
    
       
    # Resolve library paths
    if(NOT DEFINED ${PACKAGE}_LIB_DIR)
      if(DEFINED ${PACKAGE}_DEFAULT_LIB_DIR)
        set(${PACKAGE}_LIB_DIR ${${PACKAGE}_DEFAULT_LIB_DIR})
      endif()
    endif()
    
    
    # Libraries
    if(NOT DEFINED ${PACKAGE}_LIBS)
      if(DEFINED ${PACKAGE}_DEFAULT_LIBS)
        set(${PACKAGE}_LIBS ${${PACKAGE}_DEFAULT_LIBS})
      endif()
    endif()

    find_library(${PACKAGE}_LIBS NAMES ${${PACKAGE}_LIBS} PATHS ${${PACKAGE}_LIB_DIR} NO_DEFAULT_PATH)
    if(NOT ${PACKAGE}_LIBS)
      message(WARNING "${PACKAGE} libs not found, disabling")
      set(ENABLE_${PACKAGE} OFF)
    endif()
    
    # Definitions
    if(NOT DEFINED ${PACKAGE}_DEFINITIONS)
      if(DEFINED ${PACKAGE}_DEFAULT_DEFINITIONS)
        set(${PACKAGE}_DEFINITIONS ${${PACKAGE}_DEFAULT_DEFINITIONS})
      endif()
    endif()
        
  endif()
  
  
  # If the packages is enabled now: it was enabled, and was found
  # Display some status information
  message(STATUS "Support for ${PACKAGE}-${${PACKAGE}_VERSION} is ${ENABLE_${PACKAGE}}")
  if(${ENABLE_${PACKAGE}})    
    message(STATUS "  ${PACKAGE}-${VERSION} pkgfile: ${${PACKAGE}_PKG_FILE}")    
    message(STATUS "  ${PACKAGE}-${VERSION} include: ${${PACKAGE}_INCLUDE_DIRS}")
    message(STATUS "  ${PACKAGE}-${VERSION} libdir:  ${${PACKAGE}_LIB_DIR}")
    message(STATUS "  ${PACKAGE}-${VERSION} libs:    ${${PACKAGE}_LIBS}")
    message(STATUS "  ${PACKAGE}-${VERSION} defines: ${${PACKAGE}_DEFINITIONS}")
    
    # Append Library Build Details
    include_directories(${${PACKAGE}_INCLUDE_DIRS})
    link_directories(${${PACKAGE}_LIB_DIR})
    add_definitions(${${PACKAGE}_DEFINITIONS})
    set(BDIV_LIBS ${BDIV_LIBS} ${${PACKAGE}_LIBS})
    
  endif()
  

endmacro (bdiv_opt_pkg)



# - Try to find Thread Building Blocks
# Once done this will define
#
#  tbb_FOUND - system has tbb
#  tbb_INCLUDE_DIR - the tbb include directory
#  tbb_LIBRARY - Link these to use tbb
#  tbb_DEFINITIONS - Compiler switches required for using tbb
#


FIND_PATH(tbb_INCLUDE_DIR NAMES tbb/tbb.h
  PATHS
  $ENV{tbb_INCLUDE_DIR}
  ENV CPATH
  /usr/include
  /usr/local/include
  NO_DEFAULT_PATH
)

FIND_LIBRARY(tbb_LIBRARY NAMES tbb
  PATHS
  $ENV{tbb_LIB_DIR}
  ENV LD_LIBRARY_PATH
  ENV LIBRARY_PATH
  /usr/lib64
  /usr/local/lib64
  /usr/lib
  /usr/local/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY(tbb_LIBRARY NAMES tbb)

IF(tbb_INCLUDE_DIR AND tbb_LIBRARY)
   SET(tbb_FOUND TRUE)
ENDIF(tbb_INCLUDE_DIR AND tbb_LIBRARY)

IF(tbb_FOUND)
  IF(NOT tbb_FIND_QUIETLY)
    MESSAGE(STATUS "Found QGLViewer: ${tbb_LIBRARY}")
  ENDIF(NOT tbb_FIND_QUIETLY)
ELSE(tbb_FOUND)
  IF(tbb_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find Thread Building blocks")
  ENDIF(tbb_FIND_REQUIRED)
ENDIF(tbb_FOUND)

# show the tbb_INCLUDE_DIR and tbb_LIBRARY variables only in the advanced view
MARK_AS_ADVANCED(tbb_INCLUDE_DIR tbb_LIBRARY )

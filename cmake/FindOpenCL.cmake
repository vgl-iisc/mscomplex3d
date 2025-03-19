# - Try to find OpenCL
# This module tries to find an OpenCL implementation on your system. It supports
# AMD / ATI, Apple and NVIDIA implementations, but shoudl work, too.
#
# Once done this will define
#  OPENCL_FOUND        - system has OpenCL
#  OPENCL_INCLUDE_DIRS  - the OpenCL include directory
#  OPENCL_LIBRARIES    - link these to use OpenCL
#
# WIN32 should work, but is untested
# 
#
# Nithin:: Copied froms http://www.gitorious.org/findopencl/findopencl/blobs/raw/master/FindOpenCL.cmake
# Nithin:: Not in a mood to make a submodule etc for a single file so just using as was.. 

FIND_PACKAGE( PackageHandleStandardArgs )

SET (OPENCL_VERSION_STRING "0.1.0")
SET (OPENCL_VERSION_MAJOR 0)
SET (OPENCL_VERSION_MINOR 1)
SET (OPENCL_VERSION_PATCH 0)

IF (APPLE)

  FIND_LIBRARY(OPENCL_LIBRARIES OpenCL DOC "OpenCL lib for OSX")
  FIND_PATH(OPENCL_INCLUDE_DIRS OpenCL/cl.h DOC "Include for OpenCL on OSX")
  FIND_PATH(_OPENCL_CPP_INCLUDE_DIRS OpenCL/cl.hpp DOC "Include for OpenCL CPP bindings on OSX")

ELSE (APPLE)
	IF (WIN32)
		if(NOT DEFINED CUDA_PATH)
    		message(STATUS "CUDA_PATH not defined, trying to detect automatically...")
    
    	# Check environment variable first
    		if(DEFINED ENV{CUDA_PATH})
        		set(CUDA_PATH $ENV{CUDA_PATH})
        		message(STATUS "Found CUDA_PATH from environment: ${CUDA_PATH}")
    		else()
				# Look in the default installation directory for the latest version
            	file(GLOB CUDA_PATH_CANDIDATES 
                	 "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v*")
            
	            if(CUDA_PATH_CANDIDATES)
    	            # Sort versions to get the latest one
        	        list(SORT CUDA_PATH_CANDIDATES)
            	    list(REVERSE CUDA_PATH_CANDIDATES)
                	list(GET CUDA_PATH_CANDIDATES 0 CUDA_PATH)
                	message(STATUS "Auto-detected CUDA installation: ${CUDA_PATH}")
            	else()
                	message(WARNING "Could not auto-detect CUDA installation. You may need to manually set CUDA_PATH")
                	message(STATUS "You can find it by going to cmd and typing 'echo %CUDA_PATH%'")
                	# Fallback to a common location, but warn the user
                	set(CUDA_PATH "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.5")
                	message(STATUS "Using fallback CUDA_PATH: ${CUDA_PATH} \nYou may need to set the correct version number or change the path entirely")
            	endif()
    		endif()
		endif()

		message(STATUS "Using CUDA path = ${CUDA_PATH}")

		set(CUDA_INCLUDE_DIRS "${CUDA_PATH}/include")
		set(OPENCL_LIBRARIES "${CUDA_PATH}/lib/x64/OpenCL.lib")
		set(_OPENCL_INC_CAND "${CUDA_PATH}/include")

		# On Win32 search relative to the library
		FIND_PATH(OPENCL_INCLUDE_DIRS CL/cl.h PATHS "${_OPENCL_INC_CAND}")
		FIND_PATH(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS "${_OPENCL_INC_CAND}")

	ELSE (WIN32)

            # Unix style platforms
            FIND_LIBRARY(OPENCL_LIBRARIES OpenCL
              ENV LD_LIBRARY_PATH
            )

            GET_FILENAME_COMPONENT(OPENCL_LIB_DIR ${OPENCL_LIBRARIES} PATH)
            GET_FILENAME_COMPONENT(_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include ABSOLUTE)

            # The AMD SDK currently does not place its headers
            # in /usr/include, therefore also search relative
            # to the library
            FIND_PATH(OPENCL_INCLUDE_DIRS CL/cl.h PATHS ${_OPENCL_INC_CAND})
            FIND_PATH(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS ${_OPENCL_INC_CAND})

	ENDIF (WIN32)

ENDIF (APPLE)

FIND_PACKAGE_HANDLE_STANDARD_ARGS( OpenCL DEFAULT_MSG OPENCL_LIBRARIES OPENCL_INCLUDE_DIRS )

IF( _OPENCL_CPP_INCLUDE_DIRS )
	SET( OPENCL_HAS_CPP_BINDINGS TRUE )
	LIST( APPEND OPENCL_INCLUDE_DIRS ${_OPENCL_CPP_INCLUDE_DIRS} )
	# This is often the same, so clean up
	LIST( REMOVE_DUPLICATES OPENCL_INCLUDE_DIRS )
ENDIF( _OPENCL_CPP_INCLUDE_DIRS )

MARK_AS_ADVANCED(
  OPENCL_INCLUDE_DIRS
)
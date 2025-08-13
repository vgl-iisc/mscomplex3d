#Define the project we are working on
project(mscomplex3d)

# Look for various software we need
cmake_minimum_required(VERSION 3.15)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Common Properties
set(MSCOMPLEX3D_LOG_LEVEL 2 CACHE STRING
      "log level value trace=0,debug=1,info=2,warning=3,error=4,fatal=5")

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fpic -DUTL_LOG_LEVEL=${MSCOMPLEX3D_LOG_LEVEL}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -DUTL_LOG_LEVEL=${MSCOMPLEX3D_LOG_LEVEL}")

# Handle installation paths
#set(MSCOMPLEX3D_INSTALL_DIR "" CACHE PATH "installation dir")

#if(MSCOMPLEX3D_INSTALL_DIR)
#  set(MSCOMPLEX3D_INSTALL_DIR_BIN   ${MSCOMPLEX3D_INSTALL_DIR})
#  set(MSCOMPLEX3D_INSTALL_DIR_LIB   ${MSCOMPLEX3D_INSTALL_DIR})
#  set(MSCOMPLEX3D_INSTALL_DIR_SHARE ${MSCOMPLEX3D_INSTALL_DIR})
#else(MSCOMPLEX3D_INSTALL_DIR)
#  set(MSCOMPLEX3D_INSTALL_DIR_BIN   bin)
#  set(MSCOMPLEX3D_INSTALL_DIR_LIB   lib)
#  set(MSCOMPLEX3D_INSTALL_DIR_SHARE share/mscomplex-tri)
#endif(MSCOMPLEX3D_INSTALL_DIR)

# setup OpenCL
set(OPENCL_SDK_BUILD_UTILITY_LIBRARIES FALSE CACHE BOOL "" FORCE)
set(OPENCL_SDK_BUILD_SAMPLES FALSE CACHE BOOL "" FORCE)
set(OPENCL_SDK_TEST_SAMPLES FALSE CACHE BOOL "" FORCE)

add_subdirectory(dep/OpenCL)

set(OPENCL_SDK_BASE ${CMAKE_CURRENT_SOURCE_DIR}/dep/OpenCL)
set(OPENCL_INCLUDE_DIRS ${OPENCL_SDK_BASE}/external/OpenCL-Headers ${OPENCL_SDK_BASE}/external/OpenCL-CLHPP/include CACHE STRING "" FORCE)
set(OPENCL_LIBRARIES OpenCL::OpenCL OpenCL::Headers OpenCL::HeadersCpp CACHE STRING "" FORCE)
set(OPENCL_HPP_SOURCE ${OPENCL_SDK_BASE}/external/OpenCL-CLHPP/include/CL/opencl.hpp CACHE STRING "" FORCE)

add_subdirectory(core)
# option(BUILD_VIEWER "build the mscomplex viewer" OFF)
# if(BUILD_VIEWER)
# add_subdirectory(viewer)
# endif(BUILD_VIEWER)
option(BUILD_PYMS3D "mstri python module" OFF)
if(BUILD_PYMS3D)
add_subdirectory(pyms3d)
endif(BUILD_PYMS3D)
#Define the project we are working on
project(mscomplex3d)

# Look for various software we need
cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 20)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Manually specify the OpenCL include directory
#set(CUDA_INCLUDE_DIRS "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.5\\include")


# Manually specify the OpenCL library path
#set(OPENCL_LIBRARIES "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.5\\lib\\x64\\OpenCL.lib")

# Automatically find OpenCL
find_package(OpenCL REQUIRED)


# Manually set the Boost include and library paths
#set(Boost_INCLUDE_DIR "C://Users//sachi//OneDrive//Documents//vcpckg//vcpkg//installed//x64-windows//include")
#set(Boost_LIBRARY_DIR "C://Users//sachi//OneDrive//Documents//vcpckg//vcpkg//installed//x64-windows//lib")

# Find Boost
#find_package(Boost REQUIRED COMPONENTS date_time system program_options)

#if(Boost_FOUND)
#    include_directories(${Boost_INCLUDE_DIRS})
#else()
#    message(FATAL_ERROR "Could not find Boost libraries!")
#endif()



# Common Properties
set(MSCOMPLEX3D_LOG_LEVEL 2 CACHE STRING
      "log level value trace=0,debug=1,info=2,warning=3,error=4,fatal=5")

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fpic -DUTL_LOG_LEVEL=${MSCOMPLEX3D_LOG_LEVEL}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -DUTL_LOG_LEVEL=${MSCOMPLEX3D_LOG_LEVEL}")

# Handle installation paths
set(MSCOMPLEX3D_INSTALL_DIR "" CACHE PATH "installation dir")

if(MSCOMPLEX3D_INSTALL_DIR)
  set(MSCOMPLEX3D_INSTALL_DIR_BIN   ${MSCOMPLEX3D_INSTALL_DIR})
  set(MSCOMPLEX3D_INSTALL_DIR_LIB   ${MSCOMPLEX3D_INSTALL_DIR})
  set(MSCOMPLEX3D_INSTALL_DIR_SHARE ${MSCOMPLEX3D_INSTALL_DIR})
else(MSCOMPLEX3D_INSTALL_DIR)
  set(MSCOMPLEX3D_INSTALL_DIR_BIN   bin)
  set(MSCOMPLEX3D_INSTALL_DIR_LIB   lib)
  set(MSCOMPLEX3D_INSTALL_DIR_SHARE share/mscomplex-tri)
endif(MSCOMPLEX3D_INSTALL_DIR)

add_subdirectory(core)

# option(BUILD_VIEWER "build the mscomplex viewer" OFF)
# if(BUILD_VIEWER)
# add_subdirectory(viewer)
# endif(BUILD_VIEWER)

option(BUILD_PYMS3D "mstri python module" OFF)
if(BUILD_PYMS3D)
add_subdirectory(pyms3d)
endif(BUILD_PYMS3D)




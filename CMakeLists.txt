#Define the project we are working on
project(mscomplex3d)

# Look for various software we need
cmake_minimum_required(VERSION 2.8.3)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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




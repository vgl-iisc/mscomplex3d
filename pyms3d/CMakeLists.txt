project(mscomplex3d-pyms3d)

#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Boost.Numpy/libs/numpy/cmake ${CMAKE_MODULE_PATH})


# Manually set the Boost include and library paths
#set(Boost_INCLUDE_DIR "C://Users//sachi//OneDrive//Documents//vcpckg//vcpkg//installed//x64-windows//include")
#set(Boost_LIBRARY_DIR "C://Users//sachi//OneDrive//Documents//vcpckg//vcpkg//installed//x64-windows//lib")
#set(Boost_PYTHON_LIBRARY "C://Users//sachi//OneDrive//Documents//vcpckg//vcpkg//installed//x64-windows//lib//boost_python311-vc143 mt-x64-1_85.lib")

#set(CMAKE_PREFIX_PATH "C://Python311//Lib//site-packages//pybind11//share//cmake//pybind11")

# Use FetchContent to download pybind11
include(FetchContent)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.13.6  # You can use a specific commit or branch if needed
)

# Make sure the pybind11 content is available
FetchContent_MakeAvailable(pybind11)

# Find pybind11
#find_package(pybind11 REQUIRED)

# Manually specify the OpenCL include directory
# set(CUDA_INCLUDE_DIRS "C://Program Files//NVIDIA GPU Computing Toolkit//CUDA//v12.5//include")


# Manually specify the OpenCL library path
# set(OPENCL_LIBRARIES "C://Program Files//NVIDIA GPU Computing Toolkit//CUDA//v12.5//lib//x64//OpenCL.lib")


#find_package(Boost COMPONENTS python system date_time REQUIRED)

find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)

#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
#find_package(NumPy REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../core/
#  ${Boost_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
#  ${CMAKE_CURRENT_SOURCE_DIR}/Boost.Numpy/
)

set(BOOST_NUMPY_SRCS

  # header files
#  Boost.Numpy/boost/numpy/dtype.hpp
#  Boost.Numpy/boost/numpy/internal.hpp
#  Boost.Numpy/boost/numpy/invoke_matching.hpp
#  Boost.Numpy/boost/numpy/matrix.hpp
#  Boost.Numpy/boost/numpy/ndarray.hpp
#  Boost.Numpy/boost/numpy/numpy_object_mgr_traits.hpp
#  Boost.Numpy/boost/numpy/scalars.hpp
#  Boost.Numpy/boost/numpy/ufunc.hpp

  # source files
#  Boost.Numpy/libs/numpy/src/dtype.cpp
#  Boost.Numpy/libs/numpy/src/scalars.cpp
#  Boost.Numpy/libs/numpy/src/ndarray.cpp
#  Boost.Numpy/libs/numpy/src/matrix.cpp
#  Boost.Numpy/libs/numpy/src/ufunc.cpp
#  Boost.Numpy/libs/numpy/src/numpy.cpp
)

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Combine everything into a shared library for Python (the .pyd file)
#add_library(pyms3d_core SHARED 
#    $<TARGET_OBJECTS:mscomplex3d-core>
#pyms3d.cpp 
#${BOOST_NUMPY_SRCS}
#)

# Link the shared library to dependencies
#target_link_libraries(pyms3d_core PRIVATE 
#    ${Boost_LIBRARIES}
#    ${OPENCL_LIBRARIES}
#    pybind11::module
#)

# Set properties for Python extension module (e.g., .pyd for Windows)
#set_target_properties(pyms3d_core PROPERTIES
#    PREFIX ""
#    SUFFIX ".pyd"
#)


# Add the Python module
pybind11_add_module(pyms3d_core 
    pyms3d.cpp
    $<TARGET_OBJECTS:mscomplex3d-core>
   
)
target_link_libraries(pyms3d_core PRIVATE ${OPENCL_LIBRARIES})

install(TARGETS pyms3d_core DESTINATION ${MSCOMPLEX3D_INSTALL_DIR_BIN})

# Optional executable for testing or standalone use
#add_executable(pyms3d-exe 
#    main.cpp
#    $<TARGET_OBJECTS:mscomplex3d-core>
#    $<TARGET_OBJECTS:pyms3d_core>
#)

# Link the executable to dependencies
#target_link_libraries(pyms3d-exe PRIVATE
#    ${Boost_LIBRARIES}
#    ${OPENCL_LIBRARIES}
#)


#install(TARGETS pyms3d-exe DESTINATION ${MSCOMPLEX3D_INSTALL_DIR_BIN})

# option(INSTALL_PYMSTRI_EXAMPLES "Install the sample scripts" OFF)
# 
# if(INSTALL_PYMSTRI_EXAMPLES)
# install(DIRECTORY examples/ DESTINATION 
# 	${MSCOMPLEX_TRI_INSTALL_DIR_SHARE})
# endif(INSTALL_PYMSTRI_EXAMPLES)

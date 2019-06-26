project(mscomplex3d-pyms3d)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Boost.Numpy/libs/numpy/cmake ${CMAKE_MODULE_PATH})

find_package(Boost 1.48 COMPONENTS python system date_time REQUIRED)

find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../core/
  ${Boost_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/Boost.Numpy/
)

set(BOOST_NUMPY_SRCS

  # header files
  Boost.Numpy/boost/numpy/dtype.hpp
  Boost.Numpy/boost/numpy/internal.hpp
  Boost.Numpy/boost/numpy/invoke_matching.hpp
  Boost.Numpy/boost/numpy/matrix.hpp
  Boost.Numpy/boost/numpy/ndarray.hpp
  Boost.Numpy/boost/numpy/numpy_object_mgr_traits.hpp
  Boost.Numpy/boost/numpy/scalars.hpp
  Boost.Numpy/boost/numpy/ufunc.hpp

  # source files
  Boost.Numpy/libs/numpy/src/dtype.cpp
  Boost.Numpy/libs/numpy/src/scalars.cpp
  Boost.Numpy/libs/numpy/src/ndarray.cpp
  Boost.Numpy/libs/numpy/src/matrix.cpp
  Boost.Numpy/libs/numpy/src/ufunc.cpp
  Boost.Numpy/libs/numpy/src/numpy.cpp
)

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(pyms3d-core OBJECT pyms3d.cpp ${BOOST_NUMPY_SRCS})


add_library(pyms3d SHARED $<TARGET_OBJECTS:mscomplex3d-core>
    $<TARGET_OBJECTS:pyms3d-core>)
target_link_libraries(pyms3d ${Boost_LIBRARIES} ${OPENCL_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(pyms3d PROPERTIES PREFIX "")
set(PYTHON_SITE_PACKAGES_INSTALL_DIR "" CACHE PATH
    "installation dir for the python interface module")


add_executable(pyms3d-exe $<TARGET_OBJECTS:mscomplex3d-core>
    $<TARGET_OBJECTS:pyms3d-core> main.cpp)
target_link_libraries(pyms3d-exe ${Boost_LIBRARIES} ${OPENCL_LIBRARIES}
    ${PYTHON_LIBRARIES})

  add_custom_command(TARGET pyms3d POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy #copy the 
    $<TARGET_FILE:pyms3d> ${CMAKE_CURRENT_SOURCE_DIR}/examples)
  


if(PYTHON_SITE_PACKAGES_INSTALL_DIR)
install(TARGETS pyms3d DESTINATION ${PYTHON_SITE_PACKAGES_INSTALL_DIR})
else(PYTHON_SITE_PACKAGES_INSTALL_DIR)
install(TARGETS pyms3d DESTINATION ${MSCOMPLEX3D_INSTALL_DIR_LIB})
endif(PYTHON_SITE_PACKAGES_INSTALL_DIR)

install(TARGETS pyms3d-exe DESTINATION ${MSCOMPLEX3D_INSTALL_DIR_BIN})

# option(INSTALL_PYMSTRI_EXAMPLES "Install the sample scripts" OFF)
# 
# if(INSTALL_PYMSTRI_EXAMPLES)
# install(DIRECTORY examples/ DESTINATION 
# 	${MSCOMPLEX_TRI_INSTALL_DIR_SHARE})
# endif(INSTALL_PYMSTRI_EXAMPLES)

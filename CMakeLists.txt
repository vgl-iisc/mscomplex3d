#Define the project we are working on
project(mscomplex3d)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenCL REQUIRED)

find_package(Boost 1.48 COMPONENTS program_options thread system REQUIRED)

set(NUM_PARALLEL_TASKS 4 CACHE STRING
      "number of tasks to do in parallel")

set(OPENCL_NUM_WORK_ITEMS_PER_GROUP 256 CACHE STRING
      "number of work items in an opencl workgroup" )

set(OPENCL_NUM_WORK_GROUPS 32 CACHE STRING
      "number of work items in an opencl workgroup")

option(BUILD_EXEC_TBB "build using thread building blocks" OFF)

if(BUILD_EXEC_TBB)
  find_package(tbb REQUIRED)
endif(BUILD_EXEC_TBB)

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${OPENCL_INCLUDE_DIRS}
  ${tbb_INCLUDE_DIR}

  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include
)

set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}

    ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/cpputils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/cpputils.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/n_vector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/n_vector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/aabb.h

    ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/timer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/timer.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset_bfs.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_cl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_cl.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset.clh
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_assigngradient.cl
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_markandcollect.cl
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_ownerextrema.cl

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/cl.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    )

set(BUILD_EXEC_OPENCL 1)

configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)

add_executable(${PROJECT_NAME}  ${${PROJECT_NAME}_SRCS})

target_link_libraries(${PROJECT_NAME}
  ${Boost_LIBRARIES}
  ${tbb_LIBRARY}
  ${OPENCL_LIBRARIES}
)

# option(BUILD_EXEC_CUDA "build with cuda support" OFF)
#
#   if(BUILD_EXEC_CUDA)
#
#     find_package(CUDA REQUIRED)
#
#     set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}
#         ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/CUDA/grid_dataset.cu
#         ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/CUDA/cutil.h
#     )
#
#     CUDA_ADD_EXECUTABLE(${PROJECT_NAME}  ${${PROJECT_NAME}_SRCS})
#
#   else(BUILD_EXEC_CUDA)
#
#     add_executable(${PROJECT_NAME}  ${${PROJECT_NAME}_SRCS})
#
#   endif(BUILD_EXEC_CUDA)

option(BUILD_VTK_VIEWER "build the mscomplex viewer" ON)

if(BUILD_VTK_VIEWER)
include(vtk_viewer/vtk_viewer.cmake)
endif(BUILD_VTK_VIEWER)

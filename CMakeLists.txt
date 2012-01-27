#Define the project we are working on
project(mscomplex3d)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(BUILD_EXEC "build standard exec" ON)

set(NUM_PARALLEL_TASKS 4 CACHE STRING "number of tasks to do in parallel")

if(BUILD_EXEC)

  find_package(Boost 1.48 COMPONENTS program_options thread REQUIRED)

  set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset_bfs.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_dataset.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_datamanager.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    )

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}
  )

  set(${PROJECT_NAME}_LIBS ${${PROJECT_NAME}_LIBS} ${Boost_LIBRARIES})

  option(BUILD_EXEC_CUDA "build with cuda support" OFF)

  option(BUILD_EXEC_OPENCL "build with cuda support" ON)

  option(BUILD_EXEC_GUI "build with qt gui" OFF)

  if(NOT BUILD_EXEC_GUI)
    set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}

        ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/cpputils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/cpputils.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/n_vector.h
        ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/n_vector.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/aabb.h

        ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/timer.h
        ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/timer.cpp
        )

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utls/include)

  else(NOT BUILD_EXEC_GUI)

    add_subdirectory(utls)

    set(${PROJECT_NAME}_LIBS ${${PROJECT_NAME}_LIBS} utls)

    include_directories(${utls_SOURCE_DIR}/include)

    set(QT_MIN_VERSION "4.5.0")
    set(QT_USE_QTGUI TRUE)
    set(QT_USE_QTOPENGL TRUE)
    set(QT_USE_QTXML TRUE)

    find_package(Qt4 REQUIRED)

    include(${QT_USE_FILE})

    find_package(QGLViewer REQUIRED)

    SET(QXT_USE_QXTCORE TRUE)
    SET(QXT_USE_QXTGUI TRUE)

    find_package(Qxt REQUIRED)

    find_package(PythonInterp REQUIRED)

    set(viewer_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer.h
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer.cpp

      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.h
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.cpp
      )

    QT4_WRAP_CPP(viewer_MOC_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.h
      )

    set(viewer_UIS
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/grid_viewer_mainwindow.ui
      )

    QT4_WRAP_UI(viewer_UI_SRCS ${viewer_UIS})

    file(GLOB viewer_shader_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/cell_shader_geom.glsl
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders/cell_shader_vert.glsl
      )

    set(viewer_shader_SRCS
      ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources.h
      ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources.cpp
      )

    add_custom_command(
      OUTPUT ${viewer_shader_SRCS}
      COMMAND ${PYTHON_EXECUTABLE} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
      ${CMAKE_CURRENT_BINARY_DIR}/viewer/shadersources ${viewer_shader_FILES}
      DEPENDS ${viewer_shader_FILES} ${utls_SOURCE_DIR}/copy_files_to_c_header.py
      )

    set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}
      ${viewer_SRCS}
      ${viewer_MOC_SRCS}
      ${viewer_UI_SRCS}
      ${viewer_shader_SRCS}
    )

    include_directories(
      ${CMAKE_CURRENT_SOURCE_DIR}/viewer/
      ${CMAKE_CURRENT_BINARY_DIR}/viewer/
      ${QT_ADDITIONAL_INCLUDE_PATHS}
      ${QGLVIEWER_INCLUDE_DIR}
      ${QXT_INCLUDE_DIRS}
      )

    set(${PROJECT_NAME}_LIBS ${${PROJECT_NAME}_LIBS}
      ${QT_LIBRARIES}
      ${QT_ADDITIONAL_LIBRARIES}
      ${QGLVIEWER_LIBRARY}
      ${QXT_LIBRARIES}
      )

  endif(NOT BUILD_EXEC_GUI)

  configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)

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

  if(BUILD_EXEC_OPENCL)

    find_package(OpenCL REQUIRED)

    set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}
        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_cl.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_cl.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset.clh
        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_assigngradient.cl
        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_markandcollect.cl
        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/grid_dataset_ownerextrema.cl

        ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL/cl.hpp)

    include_directories(${OPENCL_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/OpenCL)

    set(${PROJECT_NAME}_LIBS ${${PROJECT_NAME}_LIBS} ${OPENCL_LIBRARIES})


  endif(BUILD_EXEC_OPENCL)

  add_executable(${PROJECT_NAME}  ${${PROJECT_NAME}_SRCS})

  target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_LIBS})
endif(BUILD_EXEC)

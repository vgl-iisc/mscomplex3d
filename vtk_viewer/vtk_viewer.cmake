cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
  ${PROJECT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/vtk_viewer/
)

find_package(Boost 1.48 COMPONENTS program_options thread system REQUIRED)

# Set your files and resources here
SET(VTK_VIEWER_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_viewer/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_viewer/mainwindow.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/cpputils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/cpputils.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/n_vector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/n_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/aabb.h

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/timer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/timer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/grid_mscomplex.cpp
)
SET(VTK_VIEWER_UI
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_viewer/mainwindow.ui)
SET(VTK_VIEWER_MOC_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_viewer/mainwindow.h)
SET(VTK_VIEWER_RESOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_viewer/icons/icons.qrc)

QT4_WRAP_UI(VTK_VIEWER_UISrcs ${VTK_VIEWER_UI})

QT4_WRAP_CPP(VTK_VIEWER_MOCSrcs ${VTK_VIEWER_MOC_HEADERS}
  OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)

QT4_ADD_RESOURCES(VTK_VIEWER_RCSrcs ${VTK_VIEWER_RESOURCES})

SOURCE_GROUP("Resources" FILES
  ${VTK_VIEWER_UI}
)

SOURCE_GROUP("Generated" FILES
  ${VTK_VIEWER_UISrcs}
  ${VTK_VIEWER_MOCSrcs}
  ${VTK_VIEWER_RCSrcs}
)

ADD_EXECUTABLE(mscomplex3d-vtk-viewer
  ${VTK_VIEWER_SRCS}
  ${VTK_VIEWER_UISrcs}
  ${VTK_VIEWER_MOCSrcs}
  ${VTK_VIEWER_RCSrcs}
)
TARGET_LINK_LIBRARIES(mscomplex3d-vtk-viewer vtkHybrid QVTK ${Boost_LIBRARIES})

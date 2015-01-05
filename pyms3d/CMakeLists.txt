project(mscomplex3d-pyms3d)

find_package(Boost 1.48 COMPONENTS python serialization system REQUIRED)

find_package(PythonLibs REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../core/
  ${Boost_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
)

add_library(pyms3d SHARED $<TARGET_OBJECTS:mscomplex3d-core> pyms3d.cpp)

target_link_libraries(pyms3d ${Boost_LIBRARIES})

set_target_properties(pyms3d PROPERTIES PREFIX "")

set(PYTHON_SITE_PACKAGES_INSTALL_DIR "" CACHE PATH "installation dir for the python interface module")

if(PYTHON_SITE_PACKAGES_INSTALL_DIR)
install(TARGETS pyms3d DESTINATION ${PYTHON_SITE_PACKAGES_INSTALL_DIR})
else(PYTHON_SITE_PACKAGES_INSTALL_DIR)
install(TARGETS pyms3d DESTINATION ${MSCOMPLEX3D_INSTALL_DIR_LIB})
endif(PYTHON_SITE_PACKAGES_INSTALL_DIR)

# option(INSTALL_PYMSTRI_EXAMPLES "Install the sample scripts" OFF)
# 
# if(INSTALL_PYMSTRI_EXAMPLES)
# install(DIRECTORY examples/ DESTINATION 
# 	${MSCOMPLEX_TRI_INSTALL_DIR_SHARE})
# endif(INSTALL_PYMSTRI_EXAMPLES)
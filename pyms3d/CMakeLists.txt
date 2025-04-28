project(pyms3d_core)

#set(CMAKE_PREFIX_PATH "C://Python311//Lib//site-packages//pybind11//share//cmake//pybind11")

# Use FetchContent to download pybind11
include(FetchContent)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.13.0  # You can use a specific commit or branch if needed
)
# Make sure the pybind11 content is available
FetchContent_MakeAvailable(pybind11)


find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../core/
  ${PYTHON_INCLUDE_DIRS}
)

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

pybind11_add_module(pyms3d_core 
    pyms3d.cpp
    $<TARGET_OBJECTS:mscomplex3d-core>
   
)

target_include_directories(pyms3d_core PRIVATE ${OPENCL_INCLUDE_DIRS})
target_link_libraries(pyms3d_core PRIVATE ${OPENCL_LIBRARIES})

install(TARGETS pyms3d_core
        LIBRARY DESTINATION pyms3d_core
        ARCHIVE DESTINATION pyms3d_core
        RUNTIME DESTINATION pyms3d_core) 

# Optional executable for testing or standalone use
#add_executable(pyms3d-exe 
#    main.cpp
#    $<TARGET_OBJECTS:mscomplex3d-core>
#    $<TARGET_OBJECTS:pyms3d_core>
#)

# option(INSTALL_PYMSTRI_EXAMPLES "Install the sample scripts" OFF)
# 
# if(INSTALL_PYMSTRI_EXAMPLES)
# install(DIRECTORY examples/ DESTINATION 
# 	${MSCOMPLEX_TRI_INSTALL_DIR_SHARE})
# endif(INSTALL_PYMSTRI_EXAMPLES)

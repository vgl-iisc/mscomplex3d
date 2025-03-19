# **mscomplex3d** #

The project page is hosted [here](http://vgl.csa.iisc.ac.in/mscomplex/). More details on what the Morse-Smale complex and the algorithms etc are available on project page. 

The mscomplex3d computes the Morse-Smale complexes on 3d grids. Its available in two modules, either of which can be used. 

- A python loadable module named **pyms3d**, which is easy to build, install and run. See [pyms3d/examples](pyms3d/examples/).
- A command line tool named **mscomplex3d**. This is mainly meant for debugging purposes, and for testing the source code directly.

# Installation #

## Dependencies ##
- [Cmake](http://www.cmake.org/)
- [OpenCL](https://developer.nvidia.com/cuda-toolkit)
- [OpenMP](http://openmp.org/wp/)
- [Python 3.7](http://python.org)
- [Pybind11 2.13](https://github.com/pybind/pybind11/releases/tag/v2.13.0) 

## Fetch and Compile ##

### Linux ###

```bash
$ git clone https://bitbucket.org/vgl_iisc/mscomplex-3d.git
$ cd mscomplex3d
$ git submodule update --init --recursive
$ cd ..
$ mkdir build install
$ cd build
$ cmake ../mscomplex-3d/ -DMSCOMPLEX_3D_INSTALL_DIR=../install -DBUILD_PYMS3D=1  
$ make -j8
$ make -j8 install
```

### Windows ###

The following commands build the MSVC project as well the Python module immediately in the build release directory for PYMS3D

$ cmake -DBUILD_PYMS3D=ON ../MS_COMPLEX
$ cmake --build . --clean-first --config release

Run the above commands in bash shell. 

# OpenCL Installation #

OpenCL usually comes with the CUDA Computing Toolkit from Nvidia. It is recomended to install the entire toolkit to make use of the binary or python module.

# Pybind11 Installation #

The CMake file for the python module is configured to automatically download the Pybind11 dependency. The version to download might change depending on your version of Python. The currently set version works with any Python version from 3.7 onwards.

# Unit testing via Pytest #
Whenever we build the final module, we can test and verify the correctness of its functionality using Pytest

To use Pytest, the Pytest module is required. It is recommended to create a virtual environment in the root of the repository(where the pytest.ini file is available), with the proper Pytest and numpy version associated to that version of Python which the module is created and then running "pytest -v" to perform all the unit tests in test_hydrogen.py which processes the Hydrogen_128x128x128.raw file and compares it against correct pre-computed outputs.

$ py <python_version_number> -m venv venv
$ venv\Scripts\activate
$ pip install pytest
$ pip install numpy
$ pytest -v

The test should use the pyms3d_core module. While it is by default, assumed to be in ../build/pyms3d/Release, this can be edited if it being created elsewhere on your system.

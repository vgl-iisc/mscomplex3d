# mscomplex3d #

[![Build and Test with CMake](https://github.com/vgl-iisc/mscomplex3d/actions/workflows/cmake-build-test.yml/badge.svg)](https://github.com/vgl-iisc/mscomplex3d/actions/workflows/cmake-build-test.yml) [![Build and Test with Pip & Pytest](https://github.com/vgl-iisc/mscomplex3d/actions/workflows/pip-build-test.yml/badge.svg)](https://github.com/vgl-iisc/mscomplex3d/actions/workflows/pip-build-test.yml)

mscomplex3d computes the Morse-Smale complex of a scalar function defined on a 3D grid. It is available in two modules: 

- A C++ Library, **mscomplex3d-core**, for the construction, simplification, and general manipulation of Morse-Smale complexes.
- A python loadable module named **pyms3d**, which is easy to build, install and run. See [pyms3d/examples](pyms3d/examples/).
- A command line tool named **mscomplex3d**. This is primarily meant for debugging purposes, and for testing the source code directly.

The project page is hosted [here](http://vgl.csa.iisc.ac.in/mscomplex/). Please refer to the project page for an introduction to the Morse-Smale complex, description of algorithms and software design, and sample results. 

# References #

Please cite the following publications if you use this software in your work.

@article{shivashankar2011parallel,  
  title={Parallel computation of 2D Morse-Smale complexes},  
  author={Shivashankar, Nithin and Senthilnathan, M and Natarajan, Vijay},  
  journal={IEEE Transactions on Visualization and Computer Graphics},  
  volume={18},  
  number={10},  
  pages={1757--1770},  
  year={2011},  
  publisher={IEEE}  
}  
  
@article{shivashankar2012parallel,  
author = {Shivashankar, Nithin and Natarajan, Vijay},  
title = {Parallel Computation of 3D Morse-Smale Complexes},  
journal = {Computer Graphics Forum},  
volume = {31},  
number = {3pt1},  
pages = {965-974},  
keywords = {I.3.5 Computational Geometry and Object Modeling},  
doi = {https://doi.org/10.1111/j.1467-8659.2012.03089.x},  
url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-8659.2012.03089.x},  
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1467-8659.2012.03089.x},  
year = {2012}  
}  
  
@incollection{shivashankar2015efficient,  
  title={Efficient software for programmable visual analysis using Morse-Smale complexes},  
  author={Shivashankar, Nithin and Natarajan, Vijay},  
  booktitle={Topological Methods in Data Analysis and Visualization},  
  pages={317--331},  
  year={2015},  
  publisher={Springer}  
}  

# Documentation
The repository bundles Sphinx for documentation of the python module in the `docs/` directory. The documentation can be found [here](https://vgl-iisc.github.io/mscomplex3d/). 

# Building & Installing #

## Dependencies ##

Ensure the availability of the following packages on your system before proceeding:

- [Cmake 3.15+](http://www.cmake.org/)
- [OpenCL](https://developer.nvidia.com/cuda-toolkit)
- [OpenMP](http://openmp.org/wp/)
- [Python 3.8+](http://python.org), with development dependencies on linux (python3-dev)
- [Pybind11 2.13](https://github.com/pybind/pybind11/releases/tag/v2.13.0)
- C & C++20 Compilers (gcc-11+, MSBuild, etc.)

> ⚠️ If you install OpenCL from a method other than installing the entire Nvidia CUDA Toolkit, ensure that you install the associated [ICD](https://dooglz.github.io/gpu/2019/07/26/opencl.html). <br>For end-users, we recommend simply installing the Toolkit.

Besides CMake, the compilers, OpenMP, and OpenCL, these dependencies will be installed automatically. OpenMP will likely already be bundled in your compiler, with no further effort required from your side.

## Installing the Python Package ##

To install the python package directly, run the following command (ideally in a conda or virtual environment):

```bash
pip install git+https://github.com/vgl-iisc/mscomplex3d.git
```

Now, you should be able to import the module using `import pyms3d_core`. A variety of usage examples are available in the [`examples`](./examples/), [`testing`](./testing/), and [`py_utils`](./py_utils/) directories.

If you face any issues with this method of installation, please raise an issue.

## Building Manually ##
Run the following commands from the base repository directory.

```bash
$ git clone https://github.com/vgl-iisc/mscomplex3d.git --recursive
$ cd mscomplex3d
$ mkdir build
$ cd build
$ cmake -DBUILD_PYMS3D=ON -DCMAKE_CXX_STANDARD=20 ..
$ cmake --build . --config Release 
```

> Additional Build Options
> - If you'd like to build the C sample program (`core/main.cpp`), set -DBUILD_SAMPLE=ON.
> - If you'd only like to build the C library, don't set -DBUILD_PYMS3D to ON (default OFF).
> - If needed, set -DCMAKE_CXX_COMPILER="<path_to_your_compiler>", or leave it unspecified to let CMake handle compiler discovery.
> - If you'd like to build a development version with debugging symbols attached, replace `Release` with `Debug`.

After running these commands, you should be able to find the compiled binaries in `build/pyms3d/<config>`, `build/core/<config>`, and `build/core/mscomplex3d-core.dir/<config>`, depending on which components were built.

## Building Python Wheels ##

We use the modern `build` backend with `pyproject.toml` to build python wheels. Then, you can run the following commands to build wheels for your python version:
```bash
$ pip install build
$ python -m build . --outdir dist
```

This will invoke `scikit-build` to build the python library and roll it up into an installable wheel archive, which will be placed into the `dist/` directory.

Options for CMake, such as the compilers to use, or the build configuration to be used can be changed from `pyproject.toml`.

For advanced wheel building, support for `cibuildwheel` is also available.


# Development #
## 🐞 Debugging ##

Debugging can be performed normally using standard debuggers like `gdb` or Visual Studio debug, as long as the compile config used was `Debug`. 

### VSCode ###
The repository includes sample VSCode workspace config files in `vscode`, which can be used for setting up a development and debugging environment for VSCode users. To use these files as is, simply rename the `vscode` directory to `.vscode`.

Now, you can debug `core/main.cpp` by using the `(Windows) Launch` debug target.

If you want to follow the python binding's calls into the C++ library, do the following: 

- Build and install the wheel in debug mode (`pyproject.toml line 36`)
- Run an instance of the python interpreter which imports `pyms3d_core`
- Attach to the process using the `Debug Attach` debug target.

Your python calls should now be able to trigger breakpoints set in the C++ files in `core`.

## ✅ Unit Testing via Pytest

After installing the Python module, unit tests can be run using **Pytest** to verify correctness. Ensure that you run `git lfs install && git lfs pull` before the following:

1. Create and activate a virtual environment at the **root** of the repository (where `pytest.ini` is located):

  ```bash
  $ py <python_version_number> -m venv venv
  $ venv\Scripts\activate
  ```

2. Install test dependencies:

  ```bash
  $ pip install pytest
  ```

3. Run the tests:

  ```bash
  $ pytest -v
  ```
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
$ cd ..
$ mkdir build install
$ cd build
$cmake -DBUILD_PYMS3D=ON -DCMAKE_CXX_COMPILER=/usr/bin/g++-11 -DCMAKE_CXX_STANDARD=20 <project root path from build>
$ cmake --build . --config Release(or Debug)
```

### Windows ###

```bash
$ git clone https://bitbucket.org/vgl_iisc/mscomplex-3d.git
$ cd mscomplex3d
$ cd ..
$ mkdir build install
$ cd build
$cmake -DBUILD_PYMS3D=ON -DCMAKE_CXX_COMPILER=/usr/bin/g++-11 -DCMAKE_CXX_STANDARD=20 <project root path from build>
$ cmake --build . --config Release(or Debug)
```

> **Note:** Run the above commands in a **Bash shell** (e.g., Git Bash, MSYS2, or WSL) for proper execution of Unix-style command syntax.

This will generate the C++ library and the Python bindings directly into the `Release/` directory inside the build folder.

---

## ⚙️ OpenCL Installation

OpenCL is typically installed with the **CUDA Toolkit** from NVIDIA. You can download it from the official [NVIDIA Developer site](https://developer.nvidia.com/cuda-toolkit).

> **Note:** It is recommended to install the full CUDA toolkit to ensure compatibility with OpenCL binaries and drivers, especially if using an NVIDIA GPU.

---

## 📦 Pybind11 Dependency

The `CMakeLists.txt` file for the Python module is configured to automatically **download and integrate Pybind11** during the CMake configuration step.

- The Pybind11 version used is compatible with **Python 3.7 and above**.
- No manual installation of Pybind11 is required.

---

## ✅ Unit Testing via Pytest

After building the Python module, unit tests can be run using **Pytest** to verify correctness.

1. Create and activate a virtual environment at the **root** of the repository (where `pytest.ini` is located):

```bash
$ py <python_version_number> -m venv venv
$ venv\Scripts\activate
```

2. Install test dependencies:

```bash
$ pip install pytest numpy
```

3. Run the tests:

```bash
$ pytest -v
```

This will execute the tests in `test_hydrogen.py`, which uses the `pyms3d_core` module to process `Hydrogen_128x128x128.raw` and compare outputs against precomputed reference data.

> ⚠️ The default import path for `pyms3d_core` is assumed to be `../build/pyms3d/Debug`. Modify this in `pytest.ini` if your build directory structure differs.

---

## 📦 Building Python Wheels

### 🪟 Windows

We use **[cibuildwheel](https://github.com/pypa/cibuildwheel)** to build Python wheels on Windows for multiple Python versions.

- Configuration is handled through `pyproject.toml` and CI scripts.
- This automates the wheel generation process and ensures compliance with standard Python packaging practices.

From Project root

Build a single package for a single python version

```python3.xx -m build --outdir <wheel directory>```

Build multiple directories

```python3 -m cibuildwheel --output-dir <wheel directory>```

### 🐧 Linux

Building on Linux requires a workaround due to Docker limitations with OpenCL:

- **cibuildwheel** on Linux uses **Docker** to ensure compliance with the [manylinux](https://github.com/pypa/manylinux) standard.
- Due to issues with OpenCL inside Docker containers, **we currently do not use cibuildwheel on Linux**.
- Instead, we build manually for x86_64 Linux platforms.

Build a single package for a single python version. 
From project root: 

```python3.xx -m build --outdir <wheel directory>```

> **Note:** While this does not comply with the manylinux standard, it is an acceptable compromise for our current project requirements.

---

## 🐞 Debugging the C++ Source Code on Windows

To debug the `mscomplex` C++ source using **Visual Studio (not VSCode)**:

1. Download and install Visual Studio from:  
   [https://visualstudio.microsoft.com/](https://visualstudio.microsoft.com/)

2. Open the project’s `.sln` file.

3. In the **Solution Explorer**:
   - Right-click on the `mscomplex` project.
   - Select **"Set as Startup Project"**.

4. Build and run the project to execute the `main.cpp` file in `mscomplex`.

This setup allows for full debugging support, including breakpoints, watch variables, and step-through execution within the Visual Studio IDE.

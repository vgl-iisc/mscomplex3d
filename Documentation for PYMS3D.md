Documentation for PYMS3D 

This text document will hold the documentation for PYMS3D 

Project Structure

The project is made of 2 parts. 
One is the MSComplex software, which holds the C++ code that has the primary computations for processing a scalar field, as well as the OpenCL kernels. 
The second part is Pyms3d, which performs the python bindings for functions we want to expose in the Pyms3d module. These functions usually make calls to MSComplex, requiring MSComplex to be able to compile before compiling Pyms3d.

MSComplex exe
MSComplex also has a separate .exe standalone, which can be run without performing python bindings. This is useful for debugging purposes.

Input

File type
MSComplex requires a .raw file of type float-32 to process. Other .raw files of a different type can usually be converted into a float-32.

File Dimensions
A given .raw file usually represents a scalar field in a voxel grid of some given (X,Y,Z) dimensions along the 3 axes. These dimensions must be known beforehand and set before computation to ensure memory is properly allocated. 

File size
Files can often be too large to load on the computer, given limited RAM or GPU Video memory, computation may fail. Ensure that the workstation running computation has sufficient CPU and GPU memory.

Computation

OpenCL 
MSComplex utilises OpenCL C++ bindings in order to run OpenCL kernels on the GPU for parallel computation. A more thorough look at these ideas can be seen here:

Dependency
Given OpenCL utilisation, it is required to have OpenCL dependencies on the system. It is highly recommended (as of 2025), to have these dependencies installed via the Nvidia Computing Toolkit. Downloading the right one is very simple, and installation directories are usually uniform across systems. This also makes it easier for the CMake to identify where the installation directory is.

Kernel Workload Sizing
We can decide whether to use the CPU or GPU to run OpenCL kernels, given sufficient memory. We can usually decide how this workload is then allocated across threads. This can differ based on difference in parallel computation on CPU and GPU. 

Flag and flag values
When we load the data, each voxel is assigned a 8 bit flag. The

Unit testing via Pytest
Whenever we make changes to the final module, we can easily test and verify the correctness of its functionality using Pytest

To use Pytest, the Pytest module is required. It is recommended to create a virtual environment in the root of the repository(where the pytest.ini file is available), with the proper Pytest and numpy version associated to that version of Python which the module is created and then running "pytest -v" to perform all the unit tests in test_hydrogen.py which processes the Hydrogen_128x128x128.raw file and compares it against pre-computed outputs.






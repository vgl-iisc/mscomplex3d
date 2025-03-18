import sys
import os
import pytest
#sys.path.append("C:\\Users\\sachi\\AppData\\Local\\Packages\\PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0\\LocalCache\\local-packages\\Python310\\Scripts")
#sys.path.append("C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES")

# Get the parent directory
build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build', 'pyms3d', 'Release'))

# Add it to sys.path
sys.path.insert(0, build_dir)

# Now you can import your custom module
import pyms3d_core as pyms3d

def basic_hydrogen_dataset_test():
    try:
        os.chdir('testing')
        DataFile1 = "Hydrogen_128x128x128.raw"
        Dim1      = (128,128,128)
        pyms3d.get_hw_info()
        msc = pyms3d.MsComplexPyms3D()
        msc.compute_bin(DataFile1,Dim1)
    except Exception as e: 
         pytest.fail("Error on processing the hydrogen dataset..")

    try:
        msc.simplify_pers(thresh=0.05)
    except Exception as e:  
        pytest.fail("Error on simplifying the processed hydrogen data")


    #numerical test based on the known outcomes from the hydrogen dataset
    correctNumberOfMinima=1
    correctNumberOfSaddle1=1
    correctNumberOfSaddle2=5
    correctNumberOfMaxima=4

    #if any of these are triggered, it means the package is not calculating the number of critical points correctly
    assert(correctNumberOfMinima==len(msc.cps(0)))
    assert(correctNumberOfSaddle1==len(msc.cps(1)))
    assert(correctNumberOfSaddle2==len(msc.cps(2)))
    assert(correctNumberOfMaxima==len(msc.cps(3)))

def test_data():
    basic_hydrogen_dataset_test()

def test_has_required_hardware():
    try:
        pyms3d.get_hw_info()
    except Exception as e:  
         pytest.fail("Error on finding the required GPU and hardware for OpenCL initialisation..")

def test_msc_datastructure():
    try:
        msc = pyms3d.MsComplexPyms3D()
    except Exception as e:  
         pytest.fail("Error on creating the msc object")


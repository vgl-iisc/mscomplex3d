from skbuild import setup
from setuptools import find_packages

setup(
    name="pyms3d_core",
    version="0.1.0",
    description="Python bindings for the MS3D C++ library",
    author="Sachin",
    license="MIT",
    packages=find_packages(where="pyms3d_core"),
    package_dir={"": "pyms3d_core"},
    include_package_data=True,
    cmake_install_dir="pyms3d_core",
    install_requires=[
        "numpy",  
    ],
    zip_safe=False,
)



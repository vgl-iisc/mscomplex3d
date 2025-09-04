Example Usage
====

Here is a simple example of how to use PyMS3D to compute the Morse-Smale complex of a 3D scalar field:

.. literalinclude:: ../examples/basic_processing.py
    :language: python

Note that the package can take input either as a numpy array or as a raw binary file (with specified dimensions) in fortran order. 

In either case, the datatype should be float32. See the `dataset downloader script <https://github.com/vgl-iisc/mscomplex3d/blob/main/py_utils/test_dataset_downloader/downloader.py>`_ for an example of how datasets of other formats can be converted to f32.

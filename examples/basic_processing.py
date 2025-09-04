import pyms3d_core as pyms3d

# Identify the dataset and its dimensions
ds = "../testing/Hydrogen_128x128x128.raw" # f32 raw file in fortran order
dims = (128, 128, 128)

# Select the device to be used for processing
# NOTE: This must be called before any other pyms3d function
dev = pyms3d.GPU # Options: pyms3d.GPU (0) or pyms3d.CPU (1)
pyms3d.select_device(dev)

# Compute & save the MS complex
msc = pyms3d.compute_bin(ds, dims)
msc.save("Hydrogen_128x128x128_bin.msc") # extension .msc is recommended, but not mandatory

# Alternatively, compute from a numpy array
# import numpy as np
# data = np.fromfile(ds, dtype=np.float32).reshape(dims)
# msc = pyms3d.compute_arr(data)

# Simplify the MS complex using persistence simplification
msc.simplify_pers(0.1)
msc.save("Hydrogen_128x128x128_bin_simplified.msc")

# Load the MS complex from a file
msc2 = pyms3d.MSComplex()
msc2.load("Hydrogen_128x128x128_bin.msc")
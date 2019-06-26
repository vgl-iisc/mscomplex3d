import pyms3d
import numpy as np


X,Y,Z    = [2,15,12] # the X,Y,Z dim of the array

arr1 = np.arange(X*Y*Z).reshape([X,Y,Z],order="C")
arr2 = np.arange(X*Y*Z).reshape([X,Y,Z],order="F")


for arr in [arr1,arr2]:
    msc = pyms3d.mscomplex()
    msc.compute_arr(arr)

    for x in range(X):
        for y in range(Y):
            for z in range(Z):
                if (arr[x,y,z]-msc.vert_func(x,y,z)) > 0.000001:
                    raise "Something went wrong with your indexing"
                    
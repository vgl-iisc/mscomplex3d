"""
The following is a sample using VTK 8.1.1. The hydrogen dataset 
used is from http://lgdv.cs.fau.de/External/vollib/. 
The dataset was converted to a raw binary 
array of float32 using the pvm tools given in the same website.
"""

import pyms3d
import vtk
import vtk.util.numpy_support as nps
import numpy as np


#data info
DataFile = "Hydrogen_128x128x128.raw"
Dim      = (128,128,128)

print pyms3d.get_hw_info()
print "VTK_VERSION =", vtk.VTK_VERSION

# try to load from cache
msc = pyms3d.mscomplex()
msc.compute_bin(DataFile,Dim)
msc.simplify_pers(thresh=0.05)
msc.collect_geom(dim=2,dir=1)


# the ascending geometry of 2-saddles is defined
# on the dual grid. Get coorinates of dual points. 
# i.e. The centroids of cubes
dp = msc.dual_points()
pa = vtk.vtkPoints()
pa.SetData(nps.numpy_to_vtk(dp,"Pts"))

# Cache objects
setattr(pyms3d,"msc",msc)
setattr(pyms3d,"pa",pa)

#list of 2 saddle cps
cps_2sad = msc.cps(2)
#cps_2sad = [cps_2sad[i] for i in [0,1,3]]

#put the list in cache 
setattr(pyms3d,"cps_2sad",cps_2sad)

# create a vtk CellArray for the line segments
ca = vtk.vtkCellArray()
for s in cps_2sad:
    gm = msc.asc_geom(s)
    for a,b in gm:
        ca.InsertNextCell(2)
        ca.InsertCellPoint(a)
        ca.InsertCellPoint(b)

#Set the outputs
#pd = self.GetOutput()
pd = vtk.vtkPolyData()
pd.SetPoints(pa)
pd.SetLines(ca)
print(pd)

writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName("Hydrogen_2sad_asc.vtp");
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(pd)
else:
    writer.SetInputData(pd)
writer.Write()

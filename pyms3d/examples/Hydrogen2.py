"""
The following is a sample using VTK 8.1.1. The hydrogen dataset 
used is from http://lgdv.cs.fau.de/External/vollib/. 
The dataset was converted to a raw binary 
array of float32 using the pvm tools given in the same website.

"""

import vtk,pyms3d
from vtk.util.numpy_support import numpy_to_vtk

#data info
DataFile = "Hydrogen_128x128x128.raw"
Dim      = (128,128,128)

print pyms3d.get_hw_info()
print "VTK_VERSION =", vtk.VTK_VERSION


# compute the mscomplex simplify and collect the 
# ascending geometry of 2 saddles


msc = pyms3d.mscomplex()
msc.compute_bin(DataFile,Dim)
msc.simplify_pers(thresh=0.05)
msc.collect_geom(dim=2,dir=1)

# the ascending geometry of 2-saddles is defined
# on the dual grid. Get coorinates of dual points. 
# i.e. The centroids of cubes
dp   = msc.dual_points()

# pass on the dual points
pa = vtk.vtkPoints()
pa.SetData(numpy_to_vtk(dp))

# get filtered list of saddles
cps_2sad = msc.cps(2)
cps_2sad = [cps_2sad[i] for i in [0,1,3]]

# create a vtk CellArray for the line segments
ca = vtk.vtkCellArray()
for s in cps_2sad:
    gm = msc.asc_geom(s)
    for a,b in gm:
        ca.InsertNextCell(2)
        ca.InsertCellPoint(a)
        ca.InsertCellPoint(b)

# put points and line segments into a polydata
pd = vtk.vtkPolyData()
pd.SetPoints(pa)
pd.SetLines(ca)

# clean out unused points
cleaner = vtk.vtkCleanPolyData()
cleaner.SetInputDataObject(pd)

# create line strips
stripper = vtk.vtkStripper()
stripper.SetInputConnection(cleaner.GetOutputPort())

# make the line into tubes
tuber = vtk.vtkTubeFilter()
tuber.SetInputConnection(stripper.GetOutputPort())
tuber.SetRadius(1)

# map it 
mapper = vtk.vtkDataSetMapper()
mapper.SetInputConnection(tuber.GetOutputPort())

# actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.AddPosition(2, 0, 4)
actor.GetProperty().SetDiffuseColor(0.2, 0.4, 0.8)

# renderer
ren = vtk.vtkRenderer()
ren.SetBackground(1, 1, 1)
ren.AddActor(actor)
ren.ResetCamera()

# renderWindow
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(800, 600)

# interactor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
iren.Initialize()

# Good to go
iren.Start()
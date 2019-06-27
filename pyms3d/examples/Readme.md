# **pyms3d** examples : Hydrogen #


The following examples need the the Hydrogen dataset which is available [here](https://vgl.csa.iisc.ac.in/mscomplex/data/Hydrogen_128x128x128.raw) (raw float32 data - 8Mb).
Download the same and place it in this folder. 

After compilation of the mscomplex-3d source code (as per instructions [here](../../Readme.md)) a module named **```pyms3d.so```** (linux) should be automatically available in this folder. 


## Example 1: 2saddle-max connectivity ## 

This Example computes the mcomplex of the hydrogen dataset, simplifies it using persistence, and then prints the connectivity of the 2-saddle maxima subcomplex. Run the example using the following command 

```bash
mscomplex-3d/pyms3d/examples$ python Hydrogen1.py
```


## Example 2: 2saddle-max geometry  ##

This requires vtk (tested on version 8.1.1). This example extracts the geometry of the ascending manifolds of the 2saddles and saves into a file named ```Hydrogen_2sad_asc.vtp```. This file can be loaded in vtk/paraview. To execute, run the [Hydrogen2.py](Hydrogen2.py) file similar to previous example.

## Example 3: 2saddle-max geometry visualization ##

This requires vtk (tested on version 8.1.1). This example creates a visualization of the acending manifolds of 2 saddles in vtk.  To execute, run the [Hydrogen3.py](Hydrogen3.py) file similar to previous example.

## Example 4: Paraview Visualization ##

This requires vtk and paraview installed. This example creates a visualization (in paraview) of the acending manifolds of 2 saddles, 2-saddles, maxima, and a volume visualization. From command line the code can be executed as given below. 

```bash
mscomplex-3d/pyms3d/examples$ PYTHONPATH=./ paraview --state=Hydrogen4.pvsm
```

This should produce a visualization as shown below. Here a volume rendering of the dataset is shown . The 2saddles are in yellow spheres, maxima are in red spheres and the connecting manifolds are shown in green. 

![mscomplex3d-config.png](https://vgl.csa.iisc.ac.in/mscomplex/images/hydrogen_vren_2sadSel.png)

# References #

- The hydrogen dataset used is from http://lgdv.cs.fau.de/External/vollib/. The dataset was converted to a raw binary array of float32 using the pvm tools given in the same website.

- For technical details on the mscomplex-3d algorithms  see [here](https://vgl.csa.iisc.ac.in/mscomplex/)

- For a technical interpretation of the dataset and the visualizations, see [here](https://vgl.csa.iisc.ac.in/mscomplex/pyms3dEx.html)
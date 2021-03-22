# Time of Flight neutron Imaging
Status: [![Build Status](https://travis-ci.com/neutronimaging/ToFImaging.svg?branch=master)](https://travis-ci.com/neutronimaging/ToFImaging)

This package includes modules for data reduction and analysis of various Time Of Flight (TOF) neutron imaging data. Initially designed to reduce MCP/Timepix data, most of the modules are compatible with post-processed 3D matrixes with dimension (x,y,lambda/TOF).
Currently implemented packages:
  - Various tools for TOF imaging data processing / visualization
  - Advanced Bragg edge fitting of 1d-arrays, but also stacks of ToF images
  - Bragg edge fitting using derivative's Gaussian of 1d-arrays, but also stacks of ToF images
  - Estimation of phase fraction, from a linear combination of basis functions (requires a priori material cross sections)
  - Tools for Frame overlap Bragg edge imaging (FOBI) data reduction
 
For a detailed guide of how to use this package please visit https://neutronimaging.github.io/ToFImaging/.

For questions please contact matteo.busi@psi.ch or anders.kaestner@psi.ch

# How to Install
The package can be installed via pip using the command:
```python
pip install tofimaging
```
Then the modules can be imported using e.g. the commands:
```python
import tofimaging.EdgeFitting as efit
import tofimaging.PhaseFitting as pfit
import tofimaging.ReductionTools as rt
``` 

Alternatively, the package can be cloned via git or downloaded from the website to the local machine, then loaded by the commands:
```python
import sys  
sys.path.insert(0, "path-to-repository\\src")
import tofimaging.EdgeFitting as efit
import tofimaging.PhaseFitting as pfit
import tofimaging.ReductionTools as rt
```
Make sure to update the "path-to-repository" with the path to this downloaded package in the local machine and appending the "\\src" as shown.
This procedure, may require further installation of external modules, listed in requirements.txt

# How to Use
The functions can inspected in the documentation and called in the command prompt or jupyter notebook. E.g. if you installed using the above:
```python
import tofimaging.EdgeFitting as efit
efit.GaussianBraggEdgeFitting2D(ToFdata,spectrum)
```

For the software documentation please visit https://neutronimaging.github.io/ToFImaging/ 

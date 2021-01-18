# Time of Flight neutron Imaging
Status: [![Build Status](https://travis-ci.com/neutronimaging/ToFImaging.svg?branch=master)](https://travis-ci.com/neutronimaging/ToFImaging)

This package includes modules for data reduction and analysis of various Time Of Flight (TOF) neutron imaging data. Initially designed for MCP/Timepix data, most of the modules are compatible with 3D matrixes with dimension (x,y,lambda/TOF).
Currently implemented packages:
  - Advanced Bragg edge fitting of 1d-arrays, but also images
  - Bragg edge fitting using derivative's Gaussian of 1d-arrays, but also images
  - Estimation of phase fraction, from a linear combination of basis functions (a priori defined phases)
  - Tools for Frame overlap Bragg edge imaging (FOBI) data reduction
  - Various tools for TOF imaging data processing / visualization
 
For a detailed guide of how to use this package please visit https://neutronimaging.github.io/ToFImaging/.

For questions please contact matteo.busi@psi.ch or anders.kaestner@psi.ch

# How to Use
Clone or download the codes to your local machine, then insert the following two lines to your python codes:
```python
import sys  
sys.path.insert(0, "path-to-repository")
```
This may require further installation of external modules, listed in requirements.txt

For the software documentation please visit https://neutronimaging.github.io/ToFImaging/ 

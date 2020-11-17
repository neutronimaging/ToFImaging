# Time of Flight neutron Imaging
This package includes modules for data reduction and analysis of various Time Of Flight (TOF) neutron imaging data. Initially designed for MCP/Timepix data, most of the modules are compatible with 3D matrixes with dimension (x,y,lambda/TOF).
Currently implemented packages:
  - Advanced Bragg edge fitting of 1d-arrays, but also images
  - Bragg edge fitting using derivative's Gaussian of 1d-arrays, but also images
  - Estimation of phase fraction, from a linear combination of basis functions (a priori defined phases)
  - Tools for Frame overlap Bragg edge imaging (FOBI) data reduction
  - Various tools for TOF imaging data processing / visualization
 
For a detailed guid of how to use this package please open docs/index.md from the Github repository.
For questions please contact matteo.busi@psi.ch or anders.kaestner@psi.ch

# How to Use
Clone or download the codes to your local machine, move to the right folder then install the
dependencies using the following command

```
$ pip install -r requirements.txt
```

and then to install the package

```
$ pip install -e .
```

For the software documentation please visit https://neutronimaging.github.io/ToFImaging/ 

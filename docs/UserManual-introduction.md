[Return to table of contents](index.md)<br/>
# Introduction
This package includes modules for data reduction and analysis of various Time Of Flight (TOF) neutron imaging data. Initially designed for MCP/Timepix data, most of the modules are compatible with 3D matrixes with dimension (x,y,lambda/TOF).
Currently implemented packages:
  - Advanced Bragg edge fitting of 1d-arrays, but also images
  - Bragg edge fitting using derivative's Gaussian of 1d-arrays, but also images
  - Estimation of phase fraction, from a linear combination of basis functions (a priori defined phases)
  - Tools for Frame overlap Bragg edge imaging (FOBI) data reduction
  - Various tools for TOF imaging data processing / visualization
  
Please find the instruction for how to use the functions as comments in the modules.
All the relevant files, with python modules and example jupyter notebooks are in /scripts/.
For questions please contact matteo.busi@psi.ch or anders.kaestner@psi.ch
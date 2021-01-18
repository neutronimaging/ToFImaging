[Return to table of contents](index.md)<br/>
# Introduction
This package includes modules for data reduction and analysis of various Time Of Flight (TOF) neutron imaging data. Initially designed for MCP/Timepix data, most of the modules are compatible with 3Darrays with dimension (x,y,lambda/TOF/bin).
Currently implemented packages:
  - Advanced Bragg edge fitting of 1Darrays and 3Darrays.
  - Bragg edge fitting using derivative's Gaussian of 1Darrays and 3Darrays.
  - Estimation of phase fraction, from a linear combination of basis functions (a priori defined phases) of 1Darrays and 3Darrays.
  - Tools for Frame overlap Bragg edge imaging (FOBI) data reduction.
  - Various tools for TOF imaging data processing / visualization.

The python modules can be found in "python_modules/".
Example jupyter notebooks of how to process the data can be found in "jupyter_notebooks/".
For questions please contact [matteo.busi@psi.ch](matteo.busi@psi.ch) or [anders.kaestner@psi.ch](anders.kaestner@psi.ch).
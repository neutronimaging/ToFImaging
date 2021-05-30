[Return to table of contents](index.md)<br/>
# Installation
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
These are the three modules which are used for [edge fitting](UserManual-EdgeFitting.md), [crystallographic phase fitting](UserManual-PhaseFitting.md) and [tools for data reduction/normalization](UserManual-ReductionTools.md) respectively. To inspect the list with input and outputs of the fitting fuction please consult the respective pages, for examples of how to use them for data please inspect the [Jupiter Notebooks](UserManual-jupyter.md)


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
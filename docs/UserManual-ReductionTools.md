[Return to table of contents](index.md)<br/>
# Reduction Tools
The reduction_tools.py module contains various functions for the data processing/reduction/filtering.
Please inspect the file for the full list of the functions. Here we report on a few key functions.

## moving_average_2D
Carries out moving average with 2D kernel on an image or sequentially to a stack of ToF images if a 3D matrix is given. Either option selected for the kernel, this is normalized such as the sum of elements is one.

__INPUTS__:

|Parameter| Description|
|----------|------------|
| mysignal | 2D or 3D array with the image or stack of images [REQUIRED]|
| kernel_size | Size of a square kernel (integer) [Default = 3]|
| rect_kernel | Row and columns of a rectangular kernel ([rows,columns])  [Default = []]|
| custom_kernel | A user defined kernel kernel [Default = []]|

## savitzky_golay
Advance smoothing of a signal. For documentation see [here](https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html).

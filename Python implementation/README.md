# Python implementation
This folder contains the following files:
1) KinDrape.py (Python implementation of the original MATLAB KinDrape)
2) KinDrape_eff_NR.py (Python implementation of the MATLAB KinDrape_eff_NR)
3) CallKinDrape.py (script to define input and call KinDrape.py)

The motivation for the Python implementation was to enable users without a MATLAB license
to interact with the kinematic draping code. The code was developed using Python 3.8 and 
can e.g. be executed from the Spyder environment available with the Anaconda distribution 
(https://www.anaconda.com/). Here, interactive plots much like MATLAB are supported.

The code can also be accessed online through MyBinder, i.e. without the need to download and
install the Python environment. Please see the MyBinder badge at the bottom of the 
[primary readme file](../README.md) in the root directory. 

Please note that the focus in the implementation was to create a code similar in structure 
to the original MATLAB code and therefore not necessarily to create an elegant Python code.

## Brief overview of input and output
Input parameters to KinDrape.py:
- d: discretization distance (scalar)
- Grid: Dimensions of fabric grid (two-component vector with number of rows and columns)
- Org: x,y origin point on mold (two-component vector with x,y)
- Ang: Initial draping direction rel. to y-axis in degrees (scalar)
- OrgNode: Origin node of grid (two-component vector with row and column)

Output parameters from KinDrape:
- Node: 3D array with computed grid nodes (first two dimensions correspond to the location in 
the grid as row/column and the third dimension/page contains x,y,z-coordinates)
- Shear: 2D array with computed shear angles of the cells (first dimension is the cell number
and the second dimension is the four shear angles of the cell)

# Original code from paper
This folder contains the following files from the journal paper:
1) KinDrape.m (original draping code, Appendix A of the paper).
2) KinDrapeOptimization.m (optimization script, Appendix B of the paper)
3) KinDrape_with_ext.m (original draping code with extensions from Sec. 4.1-4.3 implemented).

## Brief overview of input and output
Input parameters to KinDrape:
- d: discretization distance (scalar)
- Grid: Dimensions of fabric grid (two-component vector with number of rows and columns)
- Org: x,y origin point on mold (two-component vector with x,y)
- Ang: Initial draping direction rel. to y-axis in degrees (scalar)
- OrgNode: Origin node of grid (two-component vector with row and column)

The following two inputs are implemented with the extensions of the code
- PreShear: Pre-shear angle in degrees (scalar)
- Plt: Variable to enable/disable plotting (true/false)

Output parameters from KinDrape:
- Node: 3D array with computed grid nodes (first two dimensions correspond to the location in 
the grid as row/column and the third dimension/page contains x,y,z-coordinates)
- P: 3D array with data for plotting the draped cells as colored patches (the first dimension is the cell
number, the second dimension is the vertices 1-4 of each cell and the third dimension has three coordinates
and a shear angle)

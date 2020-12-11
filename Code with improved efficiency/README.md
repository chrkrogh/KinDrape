# Code with improved effiency
This folder contains the following files: 
1) KinDrape_eff.m (an updated and more efficient version of the original draping code, KinDrape)

The updates concern the method for locating constrained nodes in Step 3. As discussed in the 
journal paper, it can be done more efficiently by setting up two spheres, centered in vertex 2
and vertex 4, respectively and with a radius equal to the discretization distance, d. The two 
spheres will intersect in a circle and the problem of locating vertex 3 thus reduces to finding 
the intersection between the intersection circle and the mold surface. This is achieved with a
simple bi-section algorithm. The bi-section method is also being used to locate the second node 
in Step 1. The new additions are implemented in a new auxiliary function, MoldCircIntersecFun.

Also, the PreShear angle has been included in the initial guess for the first cell in arm 1 and 3
in Step 2 (l. 20), which adds robustness for high values of pre-shear angles.

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

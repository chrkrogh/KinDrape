This folder contains the following files: 
1) KinDrape_eff.m: an updated and more efficient version of the original draping code, KinDrape.

The updates concern the method for locating constrained nodes in Step 3. As discussed in the 
journal paper, it can be done more efficiently by setting up two spheres, centered in vertex 2
and vertex 4, respectively and with a radius equal to the discretization distance, d. The two 
spheres will intersect in a circle and the problem of locating vertex 3 thus reduces to finding 
the intersection between the intersection circle and the mold surface. This is achieved with a
simple bi-section algorithm. The method as also been adapted to locate the second node in Step 1.
The new additions are implemented in a new auxiliary function, MoldCircIntersecFun.

Also, the PreShear angle has been included in the initial guess for the first cell in each arm
in Step 2 (l. 20), which adds robustness for high values of pre-shear angles.

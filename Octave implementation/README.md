# Octave implementation
This folder contains the following files:
1) KinDrape_eff_NR_octave.m (version of KinDrape_eff_NR.m that will run in Octave 6.2)

GNU Octave is a free high-level programming language that is mostly compatible with MATLAB.
The program can be downloaded from https://www.gnu.org/software/octave/download.html

While many toolboxes and built-in functions are available in Octave, there are some limitations
which have been addressed in this implementation of KinDrape. The main issue is the limited capability of the
optimization toolbox, in particular fmincon. This was the motivation for implementing the Newton-Raphson solver
such that the draping analysis can run independently of the optimization toolbox. A further benefit of the 
Newton-Raphson solver is that the code runs more efficiently. (please see the description 
of KinDrape_eff_NR.m in the " [Code with improved efficiency folder](../Code with improved effiency/README.md) "). 

Also a few minor things have been changed for this Octave KinDrape implementation:
- Limited 2D interpolation capabilities: Therefore an analytical expression for the hemisphere mold is used.
- Trigonometric functions in degrees are not supported: They have been rewritten with radians.
- A colorbar label is not supported: It has been ommitted from the plot.

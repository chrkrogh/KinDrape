# KinDrape
This project is about a simple implementation of a kinematic draping algorithm developed
for fiber-reinforced composites and described in the journal paper:

Krogh, C., Bak, B.L.V., Lindgaard, E. et al. A simple MATLAB draping code for fiber-reinforced 
composites with application to optimization of manufacturing process parameters. 
Struct Multidisc Optim 64, 457–471 (2021). https://doi.org/10.1007/s00158-021-02925-z

The paper can be accessed in a free read-only version at: https://rdcu.be/clPiQ

Given a mold definition (either in the form of a point cloud or an analytical expression),
dimensions of the fabric/ply to be draped and some draping parameters, the code will predict
the draped fabric/ply pattern on the mold under kinematic assumptions. Please refer to the 
paper for a more elaborate description.

The code was developed at Department of Materials and Production at Aalborg University,
Denmark by Christian Krogh, Brian L.V. Bak, Esben Lindgaard, Asbjørn M. Olesen, Sebastian
M. Hermansen, Peter H. Broberg, Jørgen A. Kepler, Erik Lund, and Johnny Jakobsen. 

Various versions of the code can be found in the different folders but all of the versions
are around 100 lines in length:
- Original code from paper (the original MATLAB version presented in the paper)
- Code with improved effiency (a more efficient version of the original code)
- Python implementation (a Python version of the original code)

Please feel free to use and adapt the codes but remember to give proper attribution,
i.e. cite the journal paper and the code's DOI number: 10.5281/zenodo.4316860

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4316860.svg)](https://doi.org/10.5281/zenodo.4316860)

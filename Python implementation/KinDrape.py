import numpy as np 
from scipy.interpolate import CloughTocher2DInterpolator as CT2DInt
from scipy.optimize import (fsolve, minimize)
from mpl_toolkits.mplot3d import (Axes3D, art3d)
import matplotlib as mpl
import matplotlib.pyplot as plt
def KinDrape(d, Grid, Org, Ang, OrgNode):
    ## Mold definition: Hemisphere
    The, Phi = np.meshgrid(np.linspace(0,2*np.pi,100),
                           np.linspace(1e-6,np.pi/2-np.pi/20,50))
    X = np.cos(The)*np.sin(Phi); Y = np.sin(The)*np.sin(Phi); Z = np.cos(Phi)
    F = CT2DInt((X.ravel(),Y.ravel()),Z.ravel(),fill_value = np.min(Z))
    ## Aux. variables, solver settings and initialization of Node, P and Shear
    Node = np.empty((*Grid, 3))*np.nan
    P = [list(((np.nan,np.nan,np.nan),)*4)]*(Grid[0]-1)*(Grid[1]-1)
    Shear = np.empty(((Grid[0]-1)*(Grid[1]-1),4))*np.nan 
    Dir1 = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])  
    Dir2 = np.array([[0, 1], [-1, 0], [0, -1], [1, 0]]) 
    ## STEP 1: Place org. node (1) and node (2) defined by ini. drape angle
    # Get indices for cell, place 1st node and solve for 2nd node
    Idx = CellIdx(OrgNode[0],OrgNode[1],Dir1,Dir2,Grid,0)[0]
    Node[tuple(Idx[:,0])] = [Org[0], Org[1], F(Org[0], Org[1])]
    a_sol = fsolve(DistFun, d, args=(Node[tuple(Idx[:,0:2])], F, d, Ang))
    Node[tuple(Idx[:,0:2])] = CellVertCoor(a_sol,Node[tuple(Idx[:,0:2])],F,Ang)
    ## STEP 2: Place generator cells (initial cells) while minimizing shear
    GenStart = OrgNode + np.array([[0, 0], [1, 1], [0, 1], [0,0]])
    nGen = np.hstack([Grid[0]-OrgNode[0]-1, Grid[1]-OrgNode[1]-2,  OrgNode])
    Opt = {'maxiter': 100, 'ftol': 5e-9, 'iprint': 0, 'disp': True}
    for i in range(4):  # arms
        a_0 = 3/4*d*np.tile(CosSin(Ang+i*90),(1,2))
        for k in GenStart[i,:] + np.arange(nGen[i]).reshape(-1,1)*Dir1[i,:]:
            # Define idx and solver input. Call optimizer, assign solution
            Idx, CellNo = CellIdx(k[0],k[1],Dir1,Dir2,Grid,i) 
            BndBox = 0.5*np.array([[-1,-1,-1,-1],[1,1,1,1]])
            Bnd = tuple(map(tuple,(a_0+np.linalg.norm(a_0[0:3])*BndBox).T))
            Arg = (Node[tuple(Idx)], F, d)
            Constr = {'type': 'eq', 'fun': DistFun, 'args': Arg}
            a_sol = minimize(ObjFun, a_0, args=Arg, bounds=Bnd,
                           method='SLSQP', constraints=Constr, options=Opt)
            # Put current cell coord. and shear in P and Shear and update a_0
            Node[tuple(Idx)] = CellVertCoor(a_sol.x, Node[tuple(Idx)], F)
            P[CellNo] = list(map(tuple,Node[tuple(Idx)]))
            Shear[CellNo,:] = ShearFun(Node[tuple(Idx)])
            a_0 = a_sol.x     
    ## STEP 3: Place remaining, constrained cells in 4 quadrants between arms
    ConStart = OrgNode + np.array([[1,1],[0,1],[0,0],[1,0]])
    nCon = (nGen[[0,1,2,1,2,3,0,3]]-[1,0,0,0,0,0,1,0]).reshape(4,2)
    for i in range(4): 
        for h in ConStart[i,0] + np.arange(nCon[i,0])*(Dir1[i,0]+Dir2[i,0]):
            for w in ConStart[i,1] + np.arange(nCon[i,1])*(Dir1[i,1]+Dir2[i,1]):  
                # Define idx and solver input. Call fsolve, assign solution
                Idx, CellNo = CellIdx(h,w,Dir1,Dir2,Grid,i) 
                a_0 = Node[tuple(Idx[:,3])][0:2] - Node[tuple(Idx[:,0])][0:2] 
                a_sol = fsolve(DistFun, a_0, args=(Node[tuple(Idx)], F, d),
                                factor=1,diag=(d,d))
                Node[tuple(Idx)] = CellVertCoor(a_sol, Node[tuple(Idx)], F)
                # Put curr. cell coord. and shear in P and Shear and upd. a_0
                P[CellNo] = list(map(tuple,Node[tuple(Idx)]))
                Shear[CellNo,:] = ShearFun(Node[tuple(Idx)])
    ## Plotting
    # Create 3D figure and plot surface
    fig = plt.figure(); ax = Axes3D(fig)
    ax.plot_surface(X,Y,Z,rstride=2,cstride=2,color=(0.64,0.71,0.8),shade=True)
    # Define colormap and map the mean shear of each cell to a list of colors
    cMin = np.nanmin(Shear); cMax = np.nanmax(Shear);  
    CMapName = 'jet'; CMap = mpl.cm.get_cmap(CMapName);
    C = list(map(list,CMap(mpl.colors.Normalize(cMin,cMax)(np.mean(Shear,1)))))
    # Plot cells as colored polygons, and create axis labels and colorbar
    pc = art3d.Poly3DCollection(P,cmap=CMapName)
    pc.set_facecolor(C); pc.set_edgecolor('k'); ax.add_collection3d(pc)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z'); 
    ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))
    fig.colorbar(pc,shrink=0.75,boundaries=np.linspace(cMin, cMax, 50),
                  label='Shear angle [deg]')
    plt.show() 
    return Node, Shear
## Auxiliary functions
def CosSin(Alpha):
    # Return a unit vector in direction of Alpha
    return np.array([np.cos(Alpha*np.pi/180), np.sin(Alpha*np.pi/180)])
def CellIdx(Row,Col,Dir1,Dir2,Grid,No):
    # Return all row and col indices for the cell given vertex 1
    Rows = Row + np.array([0, Dir2[No,0], Dir1[No,0]+Dir2[No,0], Dir1[No,0]])
    Cols = Col + np.array([0, Dir2[No,1], Dir1[No,1]+Dir2[No,1], Dir1[No,1]])
    CellNo = Rows[No] + Cols[No]*(Grid[0]-1);
    return np.vstack([Rows, Cols]), CellNo
def DistFun(a, Vert, F, d, Ang=0.0):  
    # Return the difference between the edge lengths and discretization d
    Vert = CellVertCoor(a, Vert, F, Ang)
    if a.size == 4:
        EdgeLength = np.linalg.norm(Vert[[2,3,0],:] - Vert[[1,2,3],:],2,1)
    elif a.size == 2:
        EdgeLength = np.linalg.norm(Vert[[2,3],:] - Vert[[1,2],:],2,1)
    elif a.size == 1:
        EdgeLength = np.linalg.norm(Vert[0,:] - Vert[1,:])
    return EdgeLength-d  
def ShearFun(Vert):  
    # Calculate shear angles in the cell using edge vector pairs u and v
    u = Vert[[1, 2, 3, 0], :] - Vert
    v = Vert[[3, 0, 1, 2], :] - Vert
    CellAng = np.arctan2(np.linalg.norm(np.cross(u,v),2,1),np.sum(u*v,1))
    return np.abs(CellAng*180/np.pi - 90)
def ObjFun(a, Vert, F, d):
    # The objectice function from STEP 2 as sum of shear angles
    Vert = CellVertCoor(a, Vert, F)
    return np.sum(ShearFun(Vert))
def CellVertCoor(a, Vert, F, Ang=0.0):  
    # Calculate vertices in a cell given known vertices and design vars.
    if a.size == 1: # STEP 1
        Vert[1,0:2] = Vert[0,0:2]+a*CosSin(Ang+90)
        Vert[1,2] = F(Vert[1, 0], Vert[1, 1]) 
    if a.size >= 2: # STEP 2 and 3
        Vert[2,0:2] = Vert[1,0:2] + a[0:2]
        Vert[2,2] = F(Vert[1,0] + a[0],Vert[1,1] + a[1])
    if a.size == 4: # STEP 2
        Vert[3,0:2] = Vert[0,0:2] + a[2:4]
        Vert[3,2] = F(Vert[0,0] + a[2],Vert[0,1] + a[3])
    assert ~np.any(np.isnan(Vert)), 'Interpolant F is evaluated outside domain'
    return Vert

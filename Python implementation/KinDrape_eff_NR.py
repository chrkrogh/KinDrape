import numpy as np 
from scipy.interpolate import CloughTocher2DInterpolator as CT2DInt
from mpl_toolkits.mplot3d import (Axes3D, art3d)
import matplotlib as mpl
import matplotlib.pyplot as plt
def KinDrape_eff_NR(d, Grid, Org, Ang, OrgNode, PreShear, Plt):
    ## Mold definition: Hemisphere
    The, Phi = np.meshgrid(np.linspace(0,2*np.pi,100),
                           np.linspace(1e-6,np.pi/2-np.pi/20,50))
    X = np.cos(The)*np.sin(Phi); Y = np.sin(The)*np.sin(Phi); Z = np.cos(Phi)
    F = CT2DInt((X.ravel(),Y.ravel()),Z.ravel(),fill_value = np.min(Z))
    ## Mold definition: Generic double-curved mold
    #X,Y = np.meshgrid(np.linspace(0,0.5,51),np.linspace(0,0.5,51));
    #F = lambda x,y: 1.004*x + 1.089*y - 3.667*x**2 -4.4*x*y - 3.75*y**2 + \
    #    3.086*x**3 + 8.889*x**2*y + 4.321*y**3; Z = F(X,Y);
    ## Aux. variables Dir1 and Dir1. Initialization of Node, P and CellShear
    Dir1 = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])  
    Dir2 = np.array([[0, 1], [-1, 0], [0, -1], [1, 0]]) 
    Node = np.empty((*Grid, 3))*np.nan
    P = [list(((np.nan,np.nan,np.nan),)*4)]*(Grid[0]-1)*(Grid[1]-1)
    CellShear = np.empty(((Grid[0]-1)*(Grid[1]-1)))*np.nan 
    ## STEP 1: Place org. node (1) and node (2) defined by ini. drape angle
    # Get indices for cell, place 1st node and solve for 2nd node
    Idx = CellIdx(OrgNode[0],OrgNode[1],Dir1,Dir2,Grid,0)[0]
    Node[Idx[0][0],Idx[1][0]] = [Org[0], Org[1], F(Org[0], Org[1])]
    Node[Idx] = MoldCircIntersec(Node[Idx],F,d,Ang+90,1)
    ## STEP 2: Place generator cells (initial cells) while minimizing shear
    GenStart = OrgNode + np.array([[0, 0], [1, 1], [0, 1], [0, 0]])
    nGen = np.hstack([Grid[0]-OrgNode[0]-1, Grid[1]-OrgNode[1]-2,  OrgNode])
    for i in range(4):  # Arms
        CellAng_0 = Ang+i*90 + PreShear*(1+(-1)**i)/2
        for j in GenStart[i,:] + np.arange(nGen[i]).reshape(-1,1)*Dir1[i,:]:
            # Get cell idx and no. Solve for Vert #2+4 using NR, upd. CellAng_0
            Idx, CellNo = CellIdx(j[0],j[1],Dir1,Dir2,Grid,i) 
            Node[Idx], CellAng_0 = NRSol(Node[Idx],CellAng_0,F,d,PreShear,i)
            # Put current cell vertex coord. and shear in P and CellShear
            P[CellNo] = list(map(tuple,Node[Idx]))
            CellShear[CellNo] = np.abs(np.mean(ShearFun(Node[Idx],i),0))
    ## STEP 3: Place remaining, constrained cells in 4 quadrants between arms
    ConStart = OrgNode + np.array([[1,1],[0,1],[0,0],[1,0]])
    nCon = (nGen[[0,1,2,1,2,3,0,3]]-[1,0,0,0,0,0,1,0]).reshape(4,2)
    for i in range(4): # Quadrants
        for j in ConStart[i,0] + np.arange(nCon[i,0])*(Dir1[i,0]+Dir2[i,0]):
            for k in ConStart[i,1] + np.arange(nCon[i,1])*(Dir1[i,1]+Dir2[i,1]):  
                # Get cell idx and no. Call MoldCircIntersec to get Vert #3
                Idx, CellNo = CellIdx(j,k,Dir1,Dir2,Grid,i)
                Node[Idx] = MoldCircIntersec(Node[Idx],F,d,[],2)
                # Put curr. cell coord. and shear in P and CellShear
                P[CellNo] = list(map(tuple,Node[Idx]))
                CellShear[CellNo] = np.abs(np.mean(ShearFun(Node[Idx],i),0))
    ## Plotting
    if Plt:
        # Create 3D figure and plot surface (offset in z by -1 mm)
        fig = plt.figure(); ax = Axes3D(fig,auto_add_to_figure=0); fig.add_axes(ax)
        ax.plot_surface(X,Y,Z-1e-3,color=(0.64,0.71,0.80),shade=True) 
        # Define colormap and map the shear of each cell to a list of colors
        cMin = np.nanmin(CellShear); cMax = np.nanmax(CellShear)
        CMapName = 'jet'; CMap = mpl.cm.get_cmap(CMapName)
        C = list(map(list,CMap(mpl.colors.Normalize(cMin,cMax)(CellShear))))
        # Plot cells as colored polygons, and create axis labels and colorbar
        pc = art3d.Poly3DCollection(P,cmap=CMapName)
        pc.set_facecolor(C); pc.set_edgecolor('k'); ax.add_collection3d(pc)
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
        ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))
        fig.colorbar(pc,shrink=0.75,boundaries=np.linspace(cMin, cMax, 50),
                      label='Shear angle [deg]'); plt.show() 
    return Node, CellShear, ax, fig
## Auxiliary functions
def CellIdx(Row,Col,Dir1,Dir2,Grid,No):
    # Return all row and col idx. of the cell + cell no. given vert. 1 idx.
    Rows = Row + np.array([0, Dir2[No,0], Dir1[No,0]+Dir2[No,0], Dir1[No,0]])
    Cols = Col + np.array([0, Dir2[No,1], Dir1[No,1]+Dir2[No,1], Dir1[No,1]])
    CellNo = Rows[No] + Cols[No]*(Grid[0]-1)
    return tuple(np.vstack([Rows, Cols])), CellNo
def NRSol(Vert,CellAng,F,d,PreShear,i):    
    # Newton-Raphson solver to find CellAng that min. shear in cell (Step 2)
    for j in range(100): # Max iter
        # Calculate the current objective and the vertex coord. Check converg.
        Obj_curr, Vert = Step2Obj(CellAng,Vert,F,d,i,PreShear)
        if np.abs(Obj_curr) < 1e-3 or j == 100:
            break
        # Calculate a perturbed objective and a forward finite diff. gradient
        Obj_pert = Step2Obj(CellAng+1e-8,Vert,F,d,i,PreShear)[0] 
        Grad = (Obj_pert - Obj_curr)/1e-8
        # Update the design variable using a scaled step
        CellAng = CellAng - 0.5*(Obj_curr / Grad)
    return Vert, CellAng
def Step2Obj(CellAng,Vert,F,d,i,PreShear):
    # Compute all cell vertices given the angle CellAng of cell edge 1-4
    # Place: node #4 based on CellAng and d, and node #3 based on d (kinematics)
    Vert = MoldCircIntersec(Vert,F,d,CellAng,3)
    Vert = MoldCircIntersec(Vert,F,d,[],2)
    # Evaluate shear and calculate the objective
    Shear = ShearFun(Vert,i)
    Obj = np.sum(Shear - PreShear)/4
    return Obj, Vert   
def MoldCircIntersec(Vert,F,d,Ang,UnknownVertIdx):
    # Location of unknown vertex based on intersec. of circle and mold surface
    # C, R: center and radius. Vec1, Vec2: perpend. unit vec. spanning circle
    if UnknownVertIdx in (1,3): # Step 1 (2nd node) and Step 2 iterative (Vert #4)
        # Circle is constructed along angle Ang, centered in Vert #1
        C = Vert[0,:]
        R = d
        Vec1 = np.array([-np.cos(Ang*np.pi/180), -np.sin(Ang*np.pi/180), 0])
        Vec2 = np.array([0, 0, 1])
    else: # Step 3 and Step 2 iterative (Vert # 3)
        # Circle is intersection of two spheres, centered in Vert 2 and Vert 4
        C = (Vert[1,:] + Vert[3,:])/2
        R = np.sqrt(d**2 - np.linalg.norm(C-Vert[1,:])**2)
        Vec1 = (Vert[0,:]-C)/np.linalg.norm(Vert[0,:]-C)
        CircAxis = (C-Vert[1,:])/np.linalg.norm(C-Vert[1,:])
        Vec2 = np.cross(Vec1,CircAxis)/np.linalg.norm(np.cross(Vec1,CircAxis))
    # Find the intersection between the circle and the surface using bisection 
    IntersecPt = np.array([np.NaN, np.NaN, np.NaN])
    Theta = np.array([60*np.pi/180, (360-60)*np.pi/180])
    for i in range(100): # Max iter
        # Compute middle point
        Theta_mid = np.sum(Theta)/2
        # Circle pt. in 3D based on center, radius, 2 perp. vectors and angle
        CircCoor = C + R*np.cos(Theta_mid)*Vec1 + R*np.sin(Theta_mid)*Vec2
        FunVal_mid = CircCoor[2] - F(CircCoor[0],CircCoor[1])
        # Stop or adjust interval based on function value at midpoint
        if np.abs(FunVal_mid) < 5e-6: # Stopping tolerance
            IntersecPt = np.hstack((CircCoor[0:2].T,F(CircCoor[0],CircCoor[1])))
            break
        elif FunVal_mid > 0:
            Theta[0] = Theta_mid
        else:
            Theta[1] = Theta_mid
    # Store the IntersecPt in Vert at the row specified by UnknownVertIdx        
    Vert[UnknownVertIdx,:] = IntersecPt
    return Vert
def ShearFun(Vert,i):  
    # Calc. 4 signed shear angles using edge vector pairs u and v + quad. # (i)
    u = Vert[[1, 2, 3, 0], :] - Vert
    v = Vert[[3, 0, 1, 2], :] - Vert
    Shear = np.arctan2(np.linalg.norm(np.cross(u,v),2,1),np.sum(u*v,1))-np.pi/2
    return 180/np.pi*Shear*np.array([-1, 1, -1, 1])*(-1)**i
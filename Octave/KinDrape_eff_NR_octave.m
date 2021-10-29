function [Node,P] = KinDrape_eff_NR_octave(d,Grid,Org,Ang,OrgNode,PreShear,Plt) 
%% Mold definition: Analytical Hemisphere
F =@(x,y) sqrt(1-x.^2 - y.^2);
[X,Y] = meshgrid(linspace(-1,1,100)); MaskIdx = X.^2 + Y.^2 > 0.99.^2; 
X(MaskIdx) = NaN; Y(MaskIdx) = NaN; Z = F(X,Y);
%% Mold definition: Single-curved parabolic cylinder
%[X,Y] = meshgrid(0:0.01:0.5);
%F =@(x,y) 3*(x-0.25).^2; Z = F(X,Y);
%% Mold definition: Double-curved mold
%[X,Y] = meshgrid(0:0.01:0.5);
%F =@(x,y) 1.004*x + 1.089*y - 3.667*x.^2 ...
%-4.4*x.*y - 3.75*y.^2 + 3.086*x.^3 + ...
%8.889*x.^2.*y + 4.321*y.^3; Z = F(X,Y);
%% Auxiliary variables, solver settings and initialization of Node and P
Dir1 = [1 0 ; 0 1 ; -1 0 ; 0 -1]; Dir2 = [-Dir1(:,2) Dir1(:,1)];
Node = NaN([Grid 3]); P = NaN(prod(Grid-1),4,4);
%% Step 1: Place org. node (1) and node (2) defined by ini. drape angle
% Define linear idx., place node 1 and locate node 2 by circle intersection
Idx = CellIdx(Grid,OrgNode(1),OrgNode(2),Dir1,Dir2,1);
Node(Idx(1,:)) = [Org(1), Org(2), F(Org(1), Org(2))]; 
Node(Idx(2,:)) = MoldCircIntersecFun(Node(Idx(1,:)),F,d,Ang+90);
%% Step 2: Place generator cells (initial cells) while minimizing shear
GenStart = OrgNode + [0 0 ; 1 1 ; 0 1 ; 0 0];
nGenCell = [Grid-OrgNode-[0 1]  OrgNode-1]; 
for i = 1:4
    CellAng_0 = Ang+(i-1)*90 + PreShear*(1-(-1)^i)/2;
    for j = GenStart(i,:)' + (0:nGenCell(i)-1).*Dir1(i,:)'
        % Get cell idx and no. Solve using NR, get new CellAng_0
        [Idx, CellNo] = CellIdx(Grid,j(1),j(2),Dir1,Dir2,i);
        [Node(Idx), CellAng_0] = NRSol(Node(Idx),CellAng_0,F,d,PreShear,i);
        Shear = ShearFun(Node(Idx),i);
        % Put current cell coord. and shear in P array
        P(CellNo,1:4,1:4) = [Node(Idx) Shear']; 
    end
end
%% Step 3: Place remaining, constrained cells
ConStart = OrgNode + [1 1 ; 0 1 ; 0 0 ; 1 0];
nConCell = nGenCell([1 2 ; 3 2 ; 3 4 ; 1 4]) - [1 0 ; 0 0 ; 0 0 ; 1 0];
for i = 1:4
    for j = ConStart(i,1) + (0:nConCell(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = ConStart(i,2) + (0:nConCell(i,2)-1)*(Dir1(i,2)+Dir2(i,2))
            % Get cell idx and call mold circle intersection function
            [Idx, CellNo] = CellIdx(Grid,j,k,Dir1,Dir2,i);
            Node(Idx(3,:)) = MoldCircIntersecFun(Node(Idx),F,d,[]);
            % Put current cell coord. and shear in P array for patch plot
            Shear = ShearFun(Node(Idx),i);
            P(CellNo,1:4,1:4) = [Node(Idx) Shear'];
        end
    end
end
%% Plot (surface offset in z by -0.5 mm)
if Plt
    figure; plot3(Org(1),Org(2),F(Org(1),Org(2)),'kx','LineWidth',5);
    hold on; xlabel('x'); ylabel('y'); zlabel('z');
    surf(X,Y,Z-5e-3,zeros(size(Z)),'EdgeColor','none','FaceColor',[.6,.7,.8] );
    patch(P(:,:,1)',P(:,:,2)',P(:,:,3)',P(:,:,4)');
    cb = colorbar; colormap('jet'); axis('equal','tight');
end
end
%% Aux. functions
function [Idx, CellNo] = CellIdx(Grid,Row,Col,Dir1,Dir2,No)
% For a cell, return linear ind. of vert. in Node (Idx) and # in P (CellNo)
Rows = Row + [0 Dir2(No,1) Dir1(No,1)+Dir2(No,1) Dir1(No,1)]';
Cols = Col + [0 Dir2(No,2) Dir1(No,2)+Dir2(No,2) Dir1(No,2)]';
Idx = Rows + (Cols-1)*Grid(1) + ((1:3)-1)*Grid(1)*Grid(2);
CellNo = Rows(No) + (Cols(No)-1)*(Grid(1)-1);
end
function [Vert,CellAng] = NRSol(Vert,CellAng,F,d,PreShear,i)
% Newton-Raphson solver to find CellAng that min. shear in cell (Step 2)
for j = 1:1e2 % Max iter
    % Calculate the current objective and the vertex coord. Check converg.
    [Obj_curr, Vert] = Step2Obj(CellAng,Vert,F,d,i,PreShear);
    if abs(Obj_curr) < 5e-3 || j == 1e2
        break
    end
    % Calculate a perturbed objective and a forward finite diff. gradient
    Obj_pert = Step2Obj(CellAng+sqrt(eps),Vert,F,d,i,PreShear);
    Grad = (Obj_pert - Obj_curr)/sqrt(eps);
    % Update the design variable using a scaled step.
    CellAng = CellAng - 0.5*(Obj_curr / Grad);
end
end
function [Obj,Vert] = Step2Obj(CellAng,Vert,F,d,i,PreShear)
% The location of all cell vert. given the angle CellAng of cell edge 1-4
% Place: node #4 based on CellAng and d + node #3 based on d (kinematics)
Vert(4,:) = MoldCircIntersecFun(Vert(1,:),F,d,CellAng);
Vert(3,:) = MoldCircIntersecFun(Vert,F,d,[]);
% Evaluate shear and calculate the objective
Shear = ShearFun(Vert,i);
Obj = sum(Shear - PreShear)/4;
end
function IntersecPt = MoldCircIntersecFun(Vert,F,d,Ang)
% Location of unkn. vertex based on intersection of circle and mold surface
% C, R: center and radius. Vec1, Vec2: perpend. unit vec. spanning circle
if size(Vert,1) == 1 % Step 1 (2nd node) and Step 2 iterative (Vert #4)
    % Circle is constructed along angle Ang, centered in Vert(1,:)
    C = Vert(1,:);
    R = d;
    Vec1 = -[cos(Ang*pi/180) sin(Ang*pi/180) 0];
    Vec2 = [0 0 1];
else % Step 3 and Step 2 iterative (Vert #3)
    % Circle is intersection of two spheres, centered in Vert 2 and Vert 4
    C = (Vert(2,:) + Vert(4,:))/2;
    R = sqrt(d^2 - norm(C-Vert(2,:))^2);
    Vec1 = (Vert(1,:)-C)/norm(Vert(1,:)-C);
    CircleAxis = (C-Vert(2,:))/norm(C-Vert(2,:));
    Vec2 = cross(Vec1,CircleAxis)/norm(cross(Vec1,CircleAxis));
end
% Find the intersection between the circle and the surface using bisection 
IntersecPt = [NaN NaN NaN];
Theta = [60 360-60]*pi/180;
for i = 1:1e2 % Max iter
    % Compute middle point
    Theta_mid = sum(Theta)/2;
    % Circle pt. in 3D based on center, radius, 2 perp. vectors and angle
    CircCoor = C + R*cos(Theta_mid).*Vec1 + R*sin(Theta_mid).*Vec2;
    FunVal_mid = CircCoor(3) - F(CircCoor(1),CircCoor(2));
    % Stop or adjust interval based on function value at midpoint
    if abs(FunVal_mid) < 5e-6 % Stopping tolerance
        IntersecPt = [CircCoor(1) CircCoor(2) F(CircCoor(1),CircCoor(2))];
        break
    elseif FunVal_mid > 0
        Theta(1) = Theta_mid;
    else
        Theta(2) = Theta_mid;
    end
end
end
function Shear = ShearFun(Vert,i)
% Calc. 4 signed shear angles using edge vector pairs u and v + quad. # (i)
u = Vert([2 3 4 1],:)' - Vert'; 
v = Vert([4 1 2 3],:)' - Vert';
Shear = (atan2d(vecnorm(cross(u,v),2,1),dot(u,v))-90).*[1 -1 1 -1]*(-1)^i;
end

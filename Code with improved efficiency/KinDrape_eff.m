function [Node,P] = KinDrape_eff(d,Grid,Org,Ang,OrgNode,PreShear,Plt) 
%% Mold definition: Hemisphere
[Theta,Phi] = meshgrid(linspace(0,2*pi,100),linspace(pi/20,pi/2-1e-10,50));
[X,Y,Z] = sph2cart(Theta,Phi,1); 
F = scatteredInterpolant(X(:),Y(:),Z(:),'linear','linear');
%% Auxiliary variables, solver settings and initialization of Node and P
Dir1 = [1 0 ; 0 1 ; -1 0 ; 0 -1]; Dir2 = [-Dir1(:,2) Dir1(:,1)];
Opt = optimoptions(@fmincon,'Algorithm','active-set','Display','notify'...
    ,'MaxFunctionEvaluations',5e3,'StepTolerance',1e-6);
Node = NaN([Grid 3]); P = NaN(prod(Grid-1),4,4);
%% Step 1: Place org. node (1) and node (2) defined by ini. drape angle
% Define linear idx., place node 1 and locate node 2 by circle intersection
Idx = CellIdx(Grid,OrgNode(1),OrgNode(2),Dir1,Dir2,1);
Node(Idx(1,:)) = [Org(1), Org(2), F(Org(1), Org(2))]; 
Node(Idx(2,:)) = MoldCircIntersecFun(Node(Idx),F,d,Ang);
%% Step 2: Place generator cells (initial cells) while minimizing shear
GenStart = OrgNode + [0 0 ; 1 1 ; 0 1 ; 0 0];
nGenCell = [Grid-OrgNode-[0 1]  OrgNode-1]; 
for i = 1:4
    ArmAng = Ang+(i-1)*90 + PreShear*(1-(-1)^i)/2;
    a_0 = repmat(3/4*d*[cosd(ArmAng) sind(ArmAng)],1,2);
    for j = GenStart(i,:)' + (0:nGenCell(i)-1).*Dir1(i,:)'
        % Get cell idx and def. solver input. Call fmincon, assign solution
        [Idx, CellNo] = CellIdx(Grid,j(1),j(2),Dir1,Dir2,i);
        Bnd = a_0 + 1/2*norm(a_0(1:2))*[-1 -1 -1 -1; 1 1 1 1];
        a_sol = fmincon(@(a)ShearFun(a,Node(Idx),F,PreShear,i),a_0,[],...
            [],[],[],Bnd(1,:),Bnd(2,:),@(a)DistFun(a,Node(Idx),F,d),Opt);
        [~,Node(Idx),Shear] = ShearFun(a_sol,Node(Idx),F,PreShear,i);
        % Put current cell coord. and shear in P array and update a_0
        P(CellNo,1:4,1:4) = [Node(Idx) Shear']; 
        a_0 = a_sol;
    end
end
%% Step 3: Place remaining, constrained cells
ConStart = OrgNode + [1 1 ; 0 1 ; 0 0 ; 1 0];
nConCell = nGenCell([1 2 ; 3 2 ; 3 4 ; 1 4]) - [1 0 ; 0 0 ; 0 0 ; 1 0];
for i = 1:4
    for j = ConStart(i,1) + (0:nConCell(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = ConStart(i,2) + (0:nConCell(i,2)-1)*(Dir1(i,2)+Dir2(i,2))
            % Get cell idx and call circle intersection function
            [Idx, CellNo] = CellIdx(Grid,j,k,Dir1,Dir2,i);
            Node(Idx(3,:)) = MoldCircIntersecFun(Node(Idx),F,d,[]);
            % Put current cell coord. and shear in P array for patch plot
            [~,~,Shear] = ShearFun([],Node(Idx),F,PreShear,i);
            P(CellNo,1:4,1:4) = [Node(Idx) Shear'];
        end
    end
end
%% Plot
if Plt
    figure; scatter3(Org(1),Org(2),F(Org(1),Org(2)),'kx','LineWidth',5);
    hold on; axis('equal','tight'); xlabel('x'); ylabel('y'); zlabel('z');
    surf(X,Y,Z,0,'EdgeColor','none','FaceColor',[.6,.7,.8],'FaceAlpha',.5);
    patch(P(:,:,1)',P(:,:,2)',P(:,:,3)',P(:,:,4)');
    cb = colorbar; cb.Label.String = 'Shear Angle [deg]'; colormap('jet');
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
function IntersecPt = MoldCircIntersecFun(Vert,F,d,Ang)
% Location of new vertex based on intersection of circle and mold surface
% C, R: center and radius. Vec1, Vec2: perpend. unit vec. spanning circle
if all(isnan(Vert(2,:))) % Step 1
    % Circle is constructed along initial angle, centered in Vert1
    C = Vert(1,:);
    R = d;
    Vec1 = -[cosd(90+Ang) sind(90+Ang) 0];
    Vec2 = [0 0 1];
else % Step 3
    % Circle is intersection of two spheres, centered in Vert 2 and Vert 4
    C = (Vert(2,:) + Vert(4,:))/2;
    R = sqrt(d^2 - norm(C-Vert(2,:))^2);
    Vec1 = (Vert(1,:)-C)/norm(Vert(1,:)-C);
    CircleAxis = (C-Vert(2,:))/norm(C-Vert(2,:));
    Vec2 = cross(Vec1,CircleAxis)/norm(cross(Vec1,CircleAxis));
end
% Find the intersection between the circle and the surface using bisection 
Theta = [45 360-45];
for i = 1:1e3 % Max iter
    % Compute middle point
    Theta_mid = mean(Theta);
    % Circle pt. in 3D based on center, radius, 2 perp. vectors and ang.
    CircCoor = C + R*cosd(Theta_mid).*Vec1 + R*sind(Theta_mid).*Vec2;
    FunVal_mid = CircCoor(3) - F(CircCoor(1),CircCoor(2));
    % Stop or adjust interval based on function value at midpoint
    IntersecPt = NaN;
    if abs(FunVal_mid) < 1e-5 % Stopping tolerance
        IntersecPt = [CircCoor(1) CircCoor(2) F(CircCoor(1),CircCoor(2))];
        break
    elseif FunVal_mid > 0
        Theta(1) = Theta_mid;
    else
        Theta(2) = Theta_mid;
    end
end
end
function [Out1, Out2] = DistFun(a,Vert,F,d)
% Constraint function for Step 2. Output: [] (inequality constr.) and 3 x
% dist. between unknown vertices
Vert(3,:) = [Vert(2,1:2)+a(1:2) F(Vert(2,1)+a(1),Vert(2,2)+a(2))];
Vert(4,:) = [Vert(1,1:2)+a(3:4) F(Vert(1,1)+a(3),Vert(1,2)+a(4))];
Out1 = [];
Out2 = vecnorm(Vert([3,4,4],:)-Vert([2,3,1],:),2,2)' - d;
end
function [Obj, Vert, Shear] = ShearFun(a,Vert,F,PreShear,i)
% Calculate unknown vertices and shear ang. Return shear ang. sum (obj. in 
% step 2), vertex coordinates and vector of shear ang. (for P array)
if length(a) == 4 % Step 2
    Vert(3,:) = [Vert(2,1:2)+a(1:2) F(Vert(2,1)+a(1),Vert(2,2)+a(2))];
    Vert(4,:) = [Vert(1,1:2)+a(3:4) F(Vert(1,1)+a(3),Vert(1,2)+a(4))];
end
% Calculate shear angles using edge vectors u and v. Calculate sum
u = Vert([2 3 4 1],:)' - Vert'; 
v = Vert([4 1 2 3],:)' - Vert';
Shear = (atan2d(vecnorm(cross(u,v),2,1),dot(u,v))-90).*[1 -1 1 -1]*(-1)^i;
Obj = sum(abs(Shear - PreShear)); 
end

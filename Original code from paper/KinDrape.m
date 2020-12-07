function [Node,P] = KinDrape(d,Grid,Org,Ang,OrgNode)
%% Mold definition: Hemisphere
[Theta,Phi] = meshgrid(linspace(0,2*pi,100),linspace(pi/20,pi/2-1e-10,50));
[X,Y,Z] = sph2cart(Theta,Phi,1); 
F = scatteredInterpolant(X(:),Y(:),Z(:),'linear','linear');
%% Auxiliary variables, solver settings and initialization of Node and P
Dir1 = [1 0 ; 0 1 ; -1 0 ; 0 -1]; Dir2 = [-Dir1(:,2) Dir1(:,1)];
Opt1 = optimoptions(@fsolve,'Display','off'); 
Opt2 = optimoptions(@fmincon,'Algorithm','active-set','Display','notify'...
    ,'MaxFunctionEvaluations',5e3);
Node = NaN([Grid 3]); P = NaN(prod(Grid-1),4,4);
%% Step 1: Place org. node (1) and node (2) defined by ini. drape angle
% Define linear indices for cell, place 1st node and solve for 2nd node
Idx = CellIdx(Grid,OrgNode(1),OrgNode(2),Dir1,Dir2,1);
Node(Idx(1,:)) = [Org(1), Org(2), F(Org(1), Org(2))]; 
a_sol = fsolve(@(a)DistFun(a,Node(Idx(1,:)),F,d,Ang),3/4*d,Opt1);
Node(Idx(1:2,:)) = CellVertCoor(a_sol,Node(Idx(1,:)),F,Ang);
%% Step 2: Place generator cells (initial cells) while minimizing shear
GenStart = OrgNode + [0 0 ; 1 1 ; 0 1 ; 0 0];
nGenCell = [Grid-OrgNode-[0 1]  OrgNode-1]; 
for i = 1:4
    a_0 = repmat(3/4*d*[cosd(Ang+(i-1)*90) sind(Ang+(i-1)*90)],1,2);
    for j = GenStart(i,:)' + (0:nGenCell(i)-1).*Dir1(i,:)'
        % Get cell idx and def. solver input. Call fmincon, assign solution
        [Idx, CellNo] = CellIdx(Grid,j(1),j(2),Dir1,Dir2,i);
        Bnd = a_0 + 1/2*norm(a_0(1:2))*[-1 -1 -1 -1; 1 1 1 1];
        a_sol = fmincon(@(a)ShearFun(a,Node(Idx),F),a_0,[],[],[],[],...
            Bnd(1,:),Bnd(2,:),@(a)DistFun(a,Node(Idx),F,d,[]),Opt2);
        [~,Node(Idx),Shear] = ShearFun(a_sol,Node(Idx),F);
        % Put current cell coord. and shear in P array and update a_0
        P(CellNo,1:4,1:4) = [Node(Idx) Shear']; 
        a_0 = a_sol;
    end
end
%% Step 3: Place remaining, constrained cells
ConStart = OrgNode + [1 1 ; 0 1 ; 0 0 ; 1 0];
nConCell = nGenCell([1 2 ; 3 2 ; 3 4 ; 1 4]) - [1 0 ; 0 0 ; 0 0 ; 1 0];
for i = 1:4
    a_0 = 3/4*d*[cosd(Ang+(i-1)*90) sind(Ang+(i-1)*90)];
    for j = ConStart(i,1) + (0:nConCell(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = ConStart(i,2) + (0:nConCell(i,2)-1)*(Dir1(i,2)+Dir2(i,2)) 
            % Get cell idx and def. solver input. Call fsolve, assign sol.
            [Idx, CellNo] = CellIdx(Grid,j,k,Dir1,Dir2,i);
            a_sol = fsolve(@(a)DistFun(a,Node(Idx),F,d,[]),a_0,Opt1);
            [~,Node(Idx),Shear] = ShearFun(a_sol,Node(Idx),F);
            % Put current cell coord. and shear in P array and update a_0
            P(CellNo,1:4,1:4) = [Node(Idx) Shear']; 
            a_0 = a_sol;
        end
    end
end
%% Plot
figure; scatter3(Org(1),Org(2),F(Org(1),Org(2)),'kx','LineWidth',5); 
hold on; axis('equal','tight'); xlabel('x'); ylabel('y'); zlabel('z');
surf(X,Y,Z,0,'EdgeColor','none','FaceColor',[0.6,0.7,0.8],'FaceAlpha',0.5);
patch(P(:,:,1)',P(:,:,2)',P(:,:,3)',P(:,:,4)');
cb = colorbar; cb.Label.String = 'Shear Angle [deg]'; colormap('jet');
end
%% Aux. functions
function [Idx, CellNo] = CellIdx(Grid,Row,Col,Dir1,Dir2,No)
% For a cell, return linear ind. of vert. in Node (Idx) and # in P (CellNo)
Rows = Row + [0 Dir2(No,1) Dir1(No,1)+Dir2(No,1) Dir1(No,1)]';
Cols = Col + [0 Dir2(No,2) Dir1(No,2)+Dir2(No,2) Dir1(No,2)]';
Idx = Rows + (Cols-1)*Grid(1) + ((1:3)-1)*Grid(1)*Grid(2);
CellNo = Rows(No) + (Cols(No)-1)*(Grid(1)-1);
end
function Vert = CellVertCoor(a,Vert,F,Ang)
% Calculte unknown vertices in cell depending on length of design var. a 
if length(a) == 1 % Step 1 (second node rel. to origin node)
    Vert(2,1:2) = Vert(1,1:2) + a*[cosd(Ang+90) sind(Ang+90)];
    Vert(2,3) = F(Vert(2,1),Vert(2,2));
elseif length(a) == 4 % Step 2 (Vert 3 rel. to 2 and Vert 4 rel. to 1) 
    Vert(3,:) = [Vert(2,1:2)+a(1:2) F(Vert(2,1)+a(1),Vert(2,2)+a(2))];
    Vert(4,:) = [Vert(1,1:2)+a(3:4) F(Vert(1,1)+a(3),Vert(1,2)+a(4))];
elseif length(a) == 2 % Step 3 (Vert 3 rel. to 2)
    Vert(3,:) = [Vert(2,1:2)+a(1:2) F(Vert(2,1)+a(1),Vert(2,2)+a(2))];
end
end
function [Out1, Out2] = DistFun(a,Vert,F,d,Ang)
% Get cell vertices and return distance between vertices depending on step
Vert = CellVertCoor(a,Vert,F,Ang);
if length(a) == 1 % Step 1 (1 x dist)
    Out1 = norm(Vert(2,:)-Vert(1,:)) - d;
elseif length(a) == 4 % Step 2 ([] (inequality constr.) and 3 x dist)
    Out1 = [];
    Out2 = vecnorm(Vert([3,4,1],:)-Vert([2,3,4],:),2,2)' - d;
elseif length(a) == 2 % Step 3 (2 x dist.)
    Out1 = vecnorm(Vert([3,4],:)-Vert([2,3],:),2,2)' - d;
end
end
function [Obj, Vert, Shear] = ShearFun(a,Vert,F)
% Get cell vertices and calc. shear angles. Return shear ang. sum (obj. in 
% step 2), vertex coordinates and vector of shear angles (for P array)
Vert = CellVertCoor(a,Vert,F);
% Calculate shear angles using cell edge vectors u and v
u = Vert([2 3 4 1],:)' - Vert'; 
v = Vert([4 1 2 3],:)' - Vert';
Shear = abs(atan2d(vecnorm(cross(u,v),2,1),dot(u,v))-90);
Obj = sum(Shear); 
end

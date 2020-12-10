clc; clear; close all; rng('default');
% Define variables, options and call ga function
d = 0.022; Grid = [21 21]; Org0 = [0.25 0.25]; OrgNode0 = [11 11]; p = 12;
Opt3 = optimoptions(@ga,'Display','iter','MaxGenerations',4);
nDesVar = 3; lb = [-10 -10 -10]; ub = [10 9 9]; IntegerCon = [2 3];
x_opt = ga(@(x)ObjFun(x,d,Grid,Org0,OrgNode0,p,false),...
    nDesVar,[],[],[],[],lb,ub,[],IntegerCon,Opt3);
% Evaluate function with x_opt. Plot and display result
[~,Node,Shear,AngDev] = ObjFun(x_opt,d,Grid,Org0,OrgNode0,p,true);
fprintf('\nMax shear: %g, Max angle dev.: %g \n',max(Shear),max(AngDev(:)))
function [Obj,Node,Shear,AngDev] = ObjFun(x,d,Grid,Org0,OrgNode0,p,Plt)
% Obj. fun. that evaluates KinDrape and calculates shear and angle dev.
% Create input variables to KinDrape based on design variables
Ang = x(1); 
OrgNode = OrgNode0 + [x(2) x(3)]; 
Org = Org0 + d*[x(2) x(3)];
[Node,P] = KinDrape(d,Grid,Org,Ang,OrgNode,Plt);
% Calculate shear angles (gamma) as mean of each cell
Shear = mean(P(:,:,4),2); 
% Calculate the warp fiber angle deviations (psi)
WarpVec = diff(Node(:,:,1:2),1,2); WarpVec(:,:,3) = 0.0;
NomVec = reshape([0 1 0],1,1,3) .* ones(size(WarpVec));
AngDev = atan2d(vecnorm(cross(WarpVec,NomVec),2,3),dot(WarpVec,NomVec,3));
% Calculate objective as sum of p-norms
Obj = sum(abs(Shear).^p)^(1/p) + sum(abs(AngDev(:)).^p)^(1/p);
end
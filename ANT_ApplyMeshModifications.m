function UserVar=ANT_ApplyMeshModifications(UserVar)

%% This function can be used to apply modifications to a basemesh before Ua is started.
%% This might be useful if, e.g., you have a basemesh that extends beyond the intended
%% domain boundary, and you would like to deactivate certain elements of the basemesh.

% Load base mesh and domain boundary
load(UserVar.BaseMesh);
load(UserVar.BCs);

xEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
yEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));

% take care of potential nans in boundary
variable = round(UserVar.Geometry);
BCx = MeshBoundaryCoordinates.(['yr',num2str(variable)])(:,1);
BCy = MeshBoundaryCoordinates.(['yr',num2str(variable)])(:,2);
Inan = find(isnan(BCx));
if ~isempty(Inan)
    BCx = BCx(1:Inan(1)-1);
    BCy = BCy(1:Inan(1)-1);
end

% remove elements that have 1 node or more outside the boundary
[cn1,on1] = inpoly([MUA.coordinates(MUA.connectivity,1),MUA.coordinates(MUA.connectivity,2)],[BCx(:) BCy(:)]);
NodesToBeDeactivated = ~cn1;
NodesToBeDeactivated = reshape(NodesToBeDeactivated,size(MUA.connectivity));
ElementsToBeDeactivated = find(sum(NodesToBeDeactivated,2)>0);

% then remove elements for which the center of gravity is outside the
% boundary
[cn2,on2] = inpoly([xEle(:) yEle(:)],[BCx(:) BCy(:)]);
Indtmp = find(~cn2);
ElementsToBeDeactivated = unique([ElementsToBeDeactivated(:); Indtmp(:)]);

MUA_tmp=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated);

% remove 'disconnected' elements: ensure that each element has at least 1 edge
% in common with 1 other element
TR = triangulation(MUA_tmp.connectivity,MUA_tmp.coordinates);
N = neighbors(TR); N(isnan(N))=0;
N = sum(N,2); ElementsToBeDeactivated_2 = find(N==0);

MUA_tmp = DeactivateMUAelements(CtrlVar,MUA_tmp,ElementsToBeDeactivated_2);

xB = MUA_tmp.Boundary.x;
yB = MUA_tmp.Boundary.y; 
Inan = find(isnan(xB));
if ~isempty(Inan)
	xB = xB(1:Inan(1)-1);
	yB = yB(1:Inan(1)-1);
end
I = find(~inpoly([xEle(:) yEle(:)],[xB(:) yB(:)]));
UserVar.BaseMesh.DeactivatedElements = I;

[In,On] = inpoly([x(:) y(:)],[xB(:) yB(:)]);
I = [1:MUA.Nnodes]; I(unique([find(In(:)==1);find(On(:)==1)]))=[];
UserVar.BaseMesh.DeactivatedNodes = I;

MUA = MUA_tmp;

% save updated boundary and mesh for Ua
MeshBoundaryCoordinates = [MUA.Boundary.x(:) MUA.Boundary.y(:)];

MeshBoundaryCoordinatesFile = "./"+UserVar.Experiment+"/"+UserVar.Experiment+"_MeshBoundaryCoordinates";
InitialMeshFileName = "./"+UserVar.Experiment+"/"+UserVar.Experiment+"_InitialMesh";

save(MeshBoundaryCoordinatesFile,"MeshBoundaryCoordinates");
save(InitialMeshFileName,"MUA","CtrlVar");

UserVar.MeshBoundaryCoordinatesFile = "./"+UserVar.Experiment+"_MeshBoundaryCoordinates";
UserVar.InitialMeshFileName = "./"+UserVar.Experiment+"_InitialMesh";

end 

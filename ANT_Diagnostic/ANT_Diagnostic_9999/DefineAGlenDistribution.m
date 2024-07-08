function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

AGlenFile = UserVar.NameOfFileForReadingAGlenEstimate;

tmp = load(AGlenFile,'MUA','AGlen','n'); % load results from inversion
n = tmp.n(1);

% check if we need to extrude AGlen by comparing the current mesh with the mesh from the inversion
if tmp.MUA.Nnodes ~= MUA.Nnodes % something is different between the meshes
    fprintf("Using AGlen from file %s and extruding unknown values along streamlines.\n",AGlenFile);
    % map AGlen onto basemesh, with zeros where there is no data
    basemesh = load(UserVar.BaseMesh.Mesh);
    % identify nodes of the basemesh that are within the tmp mesh (i.e.
    % inversion mesh) boundaries
    Ind = find(~inpoly(basemesh.MUA.coordinates,[tmp.MUA.Boundary.x tmp.MUA.Boundary.y]));
    % map from tmp mesh to basemesh
    CtrlVar.MapOldToNew.method = "ShapeAndScattered";
    AGlen_outside = 0;
    [~,AGlen] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,basemesh.MUA,AGlen_outside,tmp.AGlen);
    % now set any values outside the tmp mesh to zero
    AGlen(Ind) = 0;
    Izero = find(AGlen==0); % indices of nodes to be updated
    % now replace those zero values with streamline advected non-zero values.
    % This takes a bit of work because we need to regrid to a regular grid and
    % have access to a velocity field.
    
    
    
    % finally interpolate back onto required mesh
    FAGlenExtruded = scatteredInterpolant(basemesh.MUA.coordinates(:,1),basemesh.MUA.coordinates(:,2),AGlen_extruded,'nearest');
    AGlen = FAGlenExtruded(MUA.coordinates(:,1),MUA.coordinates(:,2));
else
    % meshes have the same number of nodes. Check that they are indeed the
    % same
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
    x_tmp = tmp.MUA.coordinates(:,1); y_tmp = tmp.MUA.coordinates(:,2);
    if sum(abs(x-x_tmp)+abs(y-y_tmp))<eps
        fprintf("Using AGlen from file %s.\n",AGlenFile);
        AGlen = tmp.AGlen;
    else
        error("ExpID "+UserVar.ExpID+": In DefineAGlenDistribution, new and original mesh have the same number of nodes "+
              "but they appear to have different coordinates: check!");
    end
end

save("AGlen_interpolated.mat","CtrlVar","MUA","AGlen");

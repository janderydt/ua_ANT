function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

AGlenFile = UserVar.NameOfFileForReadingAGlenEstimate;

tmp = load(AGlenFile,'MUA','AGlen','n'); % load results from inversion
n = tmp.n(1);

% check if we need to extrude AGlen by comparing the current mesh with the mesh from the inversion
if tmp.MUA.Nnodes ~= MUA.Nnodes % something is different between the meshes

    % which nodes are included in the original mesh and what are the corresponding node ids?
    Tarea=TriAreaFE(tmp.MUA.coordinates,tmp.MUA.connectivity);
    tol=1e-5*sqrt(2*min(Tarea)) ;  % I found that tol=1000*eps is not enough...    
    if isempty(tmp.MUA.TR)
        tmp.MUA.TR=CreateFEmeshTriRep(tmp.MUA.connectivity,tmp.MUA.coordinates);
    end
    [ID,d] = nearestNeighbor(tmp.MUA.TR,MUA.coordinates);      
    IdenticalNodes=d<tol ;

    if ~exist("Alen-Estimate_EXTRUDED.mat","file")
        fprintf("Using AGlen from file %s and extruding unknown values along streamlines.\n",AGlenFile);
        % load basemesh
        basemesh = load(UserVar.BaseMesh.Mesh);
        % interpolate AGlen onto regular grid for extrusion
        x_r=[min(MUA_base.coordinates(:,1)):5e3:max(MUA_base.coordinates(:,1))];
        y_r=[min(MUA_base.coordinates(:,2)):5e3:max(MUA_base.coordinates(:,2))];
        [X_r,Y_r]=ndgrid(x_r,y_r);
        FAGlen = scatteredInterpolant(tmp.MUA.coordinates(:,1),tmp.MUA.coordinates(:,2),tmp.AGlen);
        AGlen_r = FAGlen(X_r,Y_r); clearvars FAGlen;
        % set values outside original boundary to zero
        Ind = find(~inpoly2([X_r(:) Y_r(:)],[tmp.MUA.Boundary.x(:) tmp.MUA.Boundary.y(:)]));
        AGlen_r(Ind) = 0;
        % define gridded interpolant of original AGlen
        FAGlen_r = griddedInterpolant(X_r,Y_r,AGlen_r);
        Create_ExtrudedFields_GriddedInterpolants([],[],FAGlen_r,0,"-scalar-AGlen-"); clearvars FAGlen_r;
    else
        load("Alen-Estimate_EXTRUDED.mat");
    end

    AGlen = ScalarInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2));   

    % now map original values of identical nodes across
    AGlen(IdenticalNodes) = tmp.AGlen(ID(IdenticalNodes));

    clearvars ScalarInterpolant;

else
    % meshes have the same number of nodes. Check that meshes are indeed the
    % same
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
    x_tmp = tmp.MUA.coordinates(:,1); y_tmp = tmp.MUA.coordinates(:,2);
    if sum(abs(x-x_tmp)+abs(y-y_tmp))<eps
        fprintf("Using AGlen from file %s.\n",AGlenFile);
        AGlen = tmp.AGlen;
    else
        error("ExpID "+UserVar.ExpID+": In DefineAGlenDistribution, new and original mesh have the same number of nodes "+...
              "but they appear to have different coordinates: check!");
    end
end

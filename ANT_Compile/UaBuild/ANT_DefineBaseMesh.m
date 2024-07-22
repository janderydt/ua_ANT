function UserVar = ANT_DefineBaseMesh(UserVar,BaseMesh)

UserVar.BaseMesh.FixedBoundaryPoints = '';

switch char(BaseMesh)
    case {'2000_meshmin5000_meshmax100000','2000_meshmin1500_meshmax100000'}
        % adjust and copy mesh files  
        UserVar.BaseMesh.BCs = UserVar.datafolder+"ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+BaseMesh+"_extrudemesh0_variableboundaryres1.mat";
        UserVar.BaseMesh.Mesh = UserVar.datafolder+"ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+BaseMesh+"_extrudemesh0_variableboundaryres1.mat";
    case {'2000_2009_2014_2018_meshmin3000_meshmax100000'}
	% adjust and copy mesh files
	    UserVar.BaseMesh.BCs = UserVar.datafolder+"ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+BaseMesh+"_extrudemesh1_variableboundaryres1.mat";
        UserVar.BaseMesh.Mesh = UserVar.datafolder+"ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+BaseMesh+"_extrudemesh1_variableboundaryres1.mat";
    case {'2000_2009_2014_2018_meshmin3000_meshmax100000_refined'}
	% adjust and copy mesh files
	    UserVar.BaseMesh.BCs = UserVar.datafolder+"ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+BaseMesh+"_extrudemesh1_variableboundaryres1.mat";
        UserVar.BaseMesh.Mesh = UserVar.datafolder+"ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+BaseMesh+"_extrudemesh1_variableboundaryres1.mat";
    case {'AS_PROPHET_2000'}
    % adjust and copy mesh files
        froot = "/mnt/md0/Ua/cases/ASE_Inversions/";
	    UserVar.BaseMesh.BCs = froot + "BoundaryCoordinates_1996_v2.mat";
        UserVar.BaseMesh.Mesh = froot + "MeshFileAdapt_Ele159177_Nod3.mat";
        UserVar.BaseMesh.FixedBoundaryPoints = froot + "FixedBoundaryPoints_1996_v2.mat";
    otherwise
        error(['ExpID ',num2str(UserVar.ExpID),': Do not recognise Mesh flag in RunTable.']);
end
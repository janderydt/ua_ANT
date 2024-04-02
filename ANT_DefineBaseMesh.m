function UserVar = ANT_DefineBaseMesh(UserVar,BaseMesh)

switch char(BaseMesh)
    case {'2000_meshmin5000_meshmax100000','2000_meshmin1500_meshmax100000'}
        % adjust and copy mesh files  
        UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+BaseMesh+"_extrudemesh0_variableboundaryres1";
        UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+BaseMesh+"_extrudemesh0_variableboundaryres1";
    case {'2000_2009_2014_2018_meshmin3000_meshmax100000'}
	% adjust and copy mesh files
	    UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+BaseMesh+"_extrudemesh1_variableboundaryres1";
        UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+BaseMesh+"_extrudemesh1_variableboundaryres1";
    case {'2000_2009_2014_2018_meshmin3000_meshmax100000_refined'}
	% adjust and copy mesh files
	    UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+BaseMesh+"_extrudemesh1_variableboundaryres1";
        UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+BaseMesh+"_extrudemesh1_variableboundaryres1";
    otherwise
        error(['ExpID ',num2str(UserVar.ExpID),': Do not recognise Mesh flag in RunTable.']);
end
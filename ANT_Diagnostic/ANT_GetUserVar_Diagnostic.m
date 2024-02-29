function UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar)

cwd = pwd;

%% perturbation parameters
UserVar.Calv = RunTable{ind,"Calv"};

%% geometry interpolants
UserVar.ISGeometry = RunTable{ind,"ISthick"}; % this is the geometry BEFORE the perturbation
switch UserVar.ISGeometry
    case {2000,2009,2014,2018}
        UserVar.ISGeometryInterpolants = [cwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.ISGeometry),'_EXTRUDED.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise ISthick flag in RunTable.']);
end

UserVar.GIGeometry = RunTable{ind,"GIthick"}; % this is the geometry BEFORE the perturbation
switch UserVar.GIGeometry
    case {2000,2009,2014,2018}
        UserVar.GIGeometryInterpolants = [cwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.GIGeometry),'_EXTRUDED.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise GIthick flag in RunTable.']);
end

%% density interpolant
UserVar.DensityInterpolant = UserVar.GIGeometryInterpolants;

%% mesh
switch char(RunTable{ind,"Mesh"}{:})

    case {'2000_meshmin5000_meshmax100000','2000_meshmin1500_meshmax100000'}
        % adjust and copy mesh files  
        UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+RunTable{ind,"Mesh"}+"_extrudemesh0_variableboundaryres1";
        UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+RunTable{ind,"Mesh"}+"_extrudemesh0_variableboundaryres1";
        UserVar = ANT_ApplyMeshModifications(UserVar);
    case {'2000_2009_2014_2018_meshmin3000_meshmax100000'}
	    % adjust and copy mesh files
	    UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+ RunTable{ind,"Mesh"}+"_extrudemesh1_variableboundaryres1";
        UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+RunTable{ind,"Mesh"}+"_extrudemesh1_variableboundaryres1";
        UserVar = ANT_ApplyMeshModifications(UserVar);
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Mesh flag in RunTable.']);
end

%% sliding law
UserVar.InverseC = RunTable{ind,"InverseC"};

% copy relevant files   
UserVar.NameOfFileForReadingSlipperinessEstimate="ANT_Inverse_"+string(UserVar.InverseC)+"_C-Estimate.mat";
copyfile("../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseC)+...
        "/"+UserVar.NameOfFileForReadingSlipperinessEstimate,...
        "./ANT_Diagnostic_"+string(RunTable{ind,"ExpID"})+...
        "/"+UserVar.NameOfFileForReadingSlipperinessEstimate);
tmp = load("../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseC)+...
        "/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"CtrlVarInRestartFile");
UserVar.SlidingLaw = tmp.CtrlVarInRestartFile.SlidingLaw;

UserVar.InverseCFill = RunTable{ind,"InverseCFill"}; 
UserVar.NameOfFileForReadingSlipperinessEstimateFill="ANT_Inverse_"+string(UserVar.InverseCFill)+"_C-Estimate.mat";
copyfile("../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseCFill)+...
        "/"+UserVar.NameOfFileForReadingSlipperinessEstimateFill,...
        "./ANT_Diagnostic_"+string(RunTable{ind,"ExpID"})+...
        "/"+UserVar.NameOfFileForReadingSlipperinessEstimateFill);

%% Glen's exponent
UserVar.InverseA = RunTable{ind,"InverseA"};

% copy relevant files   
UserVar.NameOfFileForReadingAGlenEstimate="ANT_Inverse_"+string(UserVar.InverseA)+"_AGlen-Estimate.mat";
copyfile("../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseA)+...
        "/"+UserVar.NameOfFileForReadingAGlenEstimate,...
        "./ANT_Diagnostic_"+string(RunTable{ind,"ExpID"})+...
        "/"+UserVar.NameOfFileForReadingAGlenEstimate);

UserVar.InverseAFill = RunTable{ind,"InverseAFill"}; 
UserVar.NameOfFileForReadingAGlenEstimateFill="ANT_Inverse_"+string(UserVar.InverseAFill)+"_AGlen-Estimate.mat";
copyfile("../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseAFill)+...
        "/"+UserVar.NameOfFileForReadingAGlenEstimateFill,...
        "./ANT_Diagnostic_"+string(RunTable{ind,"ExpID"})+...
        "/"+UserVar.NameOfFileForReadingAGlenEstimateFill);

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

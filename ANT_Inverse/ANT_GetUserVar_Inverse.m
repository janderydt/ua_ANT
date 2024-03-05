function UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar)

%% This function deals with all the information in the RunTable, and feeds that information 
%% to Ua via the UserVar structure. The content of this function is entirely flexible, and
%% will need to be adapted based on which parameters you need to provide to Ua at runtime.

cwd = pwd;

%% Read inversion parameters from Runtable
UserVar.Inverse.Iterations = RunTable{ind,"Iterations"};
UserVar.Inverse.Measurements = RunTable{ind,"Measurements"}{:};
UserVar.Inverse.GradientCalculation = RunTable{ind,"GradientCalc"}{:};
UserVar.Inverse.InvertFor = RunTable{ind,"InvertFor"}{:};
UserVar.Inverse.logC.gs = RunTable{ind,"gsC"};
UserVar.Inverse.logAGlen.gs = RunTable{ind,"gsA"};
UserVar.Inverse.logC.ga = RunTable{ind,"gaC"};
UserVar.Inverse.logAGlen.ga = RunTable{ind,"gaA"};

%% Read geometry interpolants from Runtable
UserVar.Geometry =RunTable{ind,"Geometry"};
switch UserVar.Geometry
    case {2000,2009,2014,2018}
        UserVar.GeometryInterpolants = [cwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.Geometry),'_EXTRUDED.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Geometry flag in RunTable.']);
end

%% density interpolant
UserVar.DensityInterpolant = UserVar.GeometryInterpolants;

%% velocity interpolants
UserVar.Velocity = RunTable{ind,"Velocity"};
switch UserVar.Velocity
    case 2000
        UserVar.VelocityInterpolants = [cwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat'];
    case {2009,2014,2018} 
        UserVar.VelocityInterpolants = [cwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_',num2str(UserVar.Velocity),'-',num2str(UserVar.Velocity+1),'_MeaSUREs_ITSLIVE_Velocities.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Velocity flag in RunTable.']);
end

%% mesh
switch char(RunTable{ind,"Mesh"}{:})

    case {'2000_meshmin5000_meshmax100000','2000_meshmin1500_meshmax100000'}
        % adjust and copy mesh files  
        UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+RunTable{ind,"Mesh"}+"_extrudemesh0_variableboundaryres1";
        UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+RunTable{ind,"Mesh"}+"_extrudemesh0_variableboundaryres1";
        UserVar = ANT_ApplyMeshModifications(UserVar);
    case {'2000_2009_2014_2018_meshmin3000_meshmax100000'}
	% adjust and copy mesh files
	    UserVar.BaseMesh.Mesh = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_meshboundarycoordinates_"+ RunTable{ind,"Mesh"}+"_extrudemesh1_variableboundaryres1";
        UserVar.BaseMesh.BCs = "../ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_"+RunTable{ind,"Mesh"}+"_extrudemesh1_variableboundaryres1";
        UserVar = ANT_ApplyMeshModifications(UserVar);
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Mesh flag in RunTable.']);
end

%% sliding law
if RunTable{ind,"Cstart"} ~= 0
    % copy relevant files   
    UserVar.NameOfFileForReadingSlipperinessEstimate="ANT_Inverse_"+string(RunTable{ind,"Cstart"})+"_C-Estimate.mat";
    copyfile("./ANT_Inverse_"+string(RunTable{ind,"Cstart"})+...
            "/"+UserVar.NameOfFileForReadingSlipperinessEstimate,...
            "./ANT_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingSlipperinessEstimate);
else
    UserVar.NameOfFileForReadingSlipperinessEstimate="";
end
UserVar.NameOfFileForSavingSlipperinessEstimate="ANT_Inverse_"+string(RunTable{ind,"ExpID"})+"_C-Estimate.mat";
UserVar.SlidingLaw = RunTable{ind,"SlidingLaw"}{:};
UserVar.SlidingCoefficient =  RunTable{ind,"m"};

%% Glen's exponent
if RunTable{ind,"Aglenstart"} ~= 0
    % copy relevant files   
    UserVar.NameOfFileForReadingAGlenEstimate="ANT_Inverse_"+string(RunTable{ind,"Aglenstart"})+"_AGlen-Estimate.mat";
    copyfile("./ANT_Inverse_"+string(RunTable{ind,"Aglenstart"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate,...
            "./ANT_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate);
else
    UserVar.NameOfFileForReadingAGlenEstimate="";
end
UserVar.NameOfFileForSavingAGlenEstimate="ANT_Inverse_"+string(RunTable{ind,"ExpID"})+"_AGlen-Estimate.mat";
UserVar.n =  RunTable{ind,"n"};

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

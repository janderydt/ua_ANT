function UserVar = ANT_GetUserVar_Transient(RunTable,ind,UserVar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% © Qing Qin    October 2024 %%%
%%% qing.qin@northumbria.ac.uk %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UserVar.StartYear = RunTable{ind,"StartYear"};
UserVar.Inverse = RunTable{ind,"Inverse"};
UserVar.InverseCycle = RunTable{ind,"InverseCycle"};
UserVar.InverseRestartFile = UserVar.datafolder+"/ANT_InputsForTransientSimulations/ANT_nsmbl_Inverse_"+...
    string(UserVar.Inverse)+"/ANT_nsmbl_Inverse_"+string(UserVar.Inverse)+"-RestartFile_InverseCycle"+...
    string(UserVar.InverseCycle)+".mat";

%% mesh variables
UserVar.InitialMeshFileName = UserVar.InverseRestartFile;
UserVar.MeshBoundaryCoordinatesFile = "MeshBoundaryCoordinates.mat";

%% geometry interpolants: only needed for bedrock geometry
UserVar.GeometryInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.StartYear)+"_EXTRUDED.mat";

%% density interpolant: same file as geometry interpolants
UserVar.DensityInterpolant = UserVar.GeometryInterpolants;  

%% sliding law and rheology
% extract information from the relevant (inversion) restart file   
NameOfFiletoRead = UserVar.InverseRestartFile;

if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
    
    % copy C field
    C = F.C; m = F.m;
    UserVar.NameOfFileForReadingSlipperinessEstimate = "ANT_Inverse_"+string(UserVar.Inverse)+"_C-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"MUA","C","m");
    UserVar.SlidingLaw = CtrlVarInRestartFile.SlidingLaw;
    
    % copy AGlen field
    AGlen = F.AGlen; n = F.n;
    UserVar.NameOfFileForReadingAGlenEstimate = "ANT_Inverse_"+string(UserVar.Inverse)+"_AGlen-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimate,"MUA","AGlen","n");
    
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to copy C and AGlen fields: "+NameOfFiletoRead);
end

UserVar.AdaptMesh=RunTable{ind,"AdaptMesh"};

%% Basal Melt
UserVar.BasalMelt = string(RunTable{ind,"BasalMelt"});
UserVar.OceForcing = string(RunTable{ind,"OceForcing"});

switch UserVar.BasalMelt
    case "PICO"
         UserVar.PICOC1 = RunTable{ind,"PICO_C1"};
         UserVar.PICOgam = RunTable{ind,"PICO_gam"};   
    case "LQ"
         UserVar.LQgam= RunTable{ind,"LQ_gam"};
    otherwise
         error("ExpID "+UserVar.ExpID+": Could not find Basal Melting Parameterization "+UserVar.BaselMelt);
end

%% SMB
UserVar.SMB = string(RunTable{ind,"SMB"});

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

%% write log file
fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Transient: ExpID %s: Starting transient experiment on %s.\n",string(UserVar.ExpID),UserVar.hostname);



function UserVar = ANT_GetUserVar_Transient(RunTable,ind,UserVar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% © Qing Qin    October 2024 %%%
%%% qing.qin@northumbria.ac.uk %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time variables
UserVar.StartTime = datetime(RunTable{ind,"ExpStartDate"},format="yyyy-MM-dd"); % physical time (in yyyy-MM-dd) when the simulation should start
UserVar.EndTime = datetime(RunTable{ind,"ExpEndDate"},format="yyyy-MM-dd"); % physical time (in yyyy-MM-dd) when the simulation should end)
DeltaSeconds = seconds(UserVar.EndTime-UserVar.StartTime); % interval between end and start year in seconds
UserVar.TotalTime = double(DeltaSeconds/(365.25*24*60*60)); % to be consistent with Ua units, convert from seconds to years
UserVar.StartTime_DecimalYears = year(UserVar.StartTime)+years(UserVar.StartTime-dateshift(UserVar.StartTime,'start','year')); % start time in Ua units (decimal years)

%% information about inverse simulation to start from
UserVar.Inverse = RunTable{ind,"Inverse"};
UserVar.InverseCycle = RunTable{ind,"InverseCycle"};
UserVar.InverseRestartFile = UserVar.datafolder+"/ANT_InputsForTransientSimulations/ANT_nsmbl_Inverse_"+...
    string(UserVar.Inverse)+"/ANT_nsmbl_Inverse_"+string(UserVar.Inverse)+"-RestartFile_InverseCycle"+...
    string(UserVar.InverseCycle)+".mat";

%% is this a restart?
UserVar.YearsCompleted = RunTable{ind,"YearsCompleted"};
UserVar.Restart = RunTable{ind,"Restart"};

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
    xC = MUA.coordinates(:,1); yC = MUA.coordinates(:,2);
    m = F.m; muk = F.muk; q = F.q;
    UserVar.NameOfFileForReadingSlipperinessEstimate = "ANT_Inverse_"+string(UserVar.Inverse)+"_C-Estimate.mat";
    % C is not updated for floating ice, except for some
    % regularization-dependent smoothing. A more representative value for C
    % beneath the ice shelves could be obtained by extruding C from the
    % grounded ice along streamlines. We use streamlines based on the
    % velocity fields obtained through the inversion. Note that is very ad-hoc
    % and users might want to change this part of the code.
    C = ExtrudeC(CtrlVarInRestartFile,MUA,F);
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"xC","yC","C","m","muk","q");
    UserVar.SlidingLaw = CtrlVarInRestartFile.SlidingLaw;
    
    % copy AGlen field
    xA = MUA.coordinates(:,1); yA = MUA.coordinates(:,2);
    AGlen = F.AGlen; n = F.n;
    UserVar.NameOfFileForReadingAGlenEstimate = "ANT_Inverse_"+string(UserVar.Inverse)+"_AGlen-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimate,"xA","yA","AGlen","n");
    
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



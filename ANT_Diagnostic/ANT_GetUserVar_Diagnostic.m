function UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar)

%% target geometry interpolants
%% ICE SHELF
UserVar.ISGeometry = RunTable{ind,"ISthick"};
switch UserVar.ISGeometry
    case {2000,2009,2014,2018,2020}
        UserVar.ISGeometryInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.ISGeometry)+"_EXTRUDED.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise ISthick flag in RunTable.']);
end

%% GROUNDED ICE
UserVar.GIGeometry = RunTable{ind,"GIthick"};
switch UserVar.GIGeometry
    case {2000,2009,2014,2018,2020}
        UserVar.GIGeometryInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.GIGeometry)+"_EXTRUDED.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise GIthick flag in RunTable.']);
end

%% CALVING
UserVar.Calv = RunTable{ind,"Calv"};
UserVar = ANT_DefineBaseMesh(UserVar,RunTable{ind,"BaseMesh"}{:});
UserVar = ANT_ApplyMeshModifications(UserVar);
    
%% density interpolant: same as ice shelf run
UserVar.DensityInterpolant = UserVar.ISGeometryInterpolants;   

%% sliding law
UserVar.InverseC = RunTable{ind,"InverseC"};

% copy relevant files   
NameOfFiletoRead="../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(UserVar.InverseC)+"/"+UserVar.Domain+"_Inverse_"+string(UserVar.InverseC)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleC"})+".mat";
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
    C = F.C; m = F.m;
    UserVar.NameOfFileForReadingSlipperinessEstimate = UserVar.Domain+"_Inverse_"+string(UserVar.InverseC)+"_C-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"MUA","C","m");
    UserVar.SlidingLaw = CtrlVarInRestartFile.SlidingLaw;
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to copy C field: "+NameOfFiletoRead);
end

%% Glen's exponent
UserVar.InverseA = RunTable{ind,"InverseA"};

% copy relevant files   
NameOfFiletoRead="../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(UserVar.InverseA)+"/"+UserVar.Domain+"_Inverse_"+string(UserVar.InverseA)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleA"})+".mat";
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
    AGlen = F.AGlen; n = F.n;
    UserVar.NameOfFileForReadingAGlenEstimate = UserVar.Domain+"_Inverse_"+string(UserVar.InverseA)+"_AGlen-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimate,"MUA","AGlen","n");
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to copy AGlen field: "+NameOfFiletoRead);
end   

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

%% write log file
fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Diagnostic: ExpID %s: Starting diagnostic experiment on %s.\n",string(UserVar.ExpID),UserVar.hostname);



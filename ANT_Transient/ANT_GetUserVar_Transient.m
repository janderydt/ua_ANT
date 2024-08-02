function UserVar = ANT_GetUserVar_Transient(RunTable,ind,UserVar)

%% geometry interpolants: only needed for bedrock geometry
UserVar.BedGeometry = 2000;
UserVar.GeometryInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.BedGeometry)+"_EXTRUDED.mat";
    
%% density interpolant: same file as geometry interpolants
UserVar.DensityInterpolant = UserVar.GeometryInterpolants;   

%% sliding law
% extract information from the relevant inversion restart file   
NameOfFiletoRead=""; %% ADD CORRECT PATH
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
    C = F.C; m = F.m;
    UserVar.NameOfFileForReadingSlipperinessEstimate = "ANT_Inverse_"+string(UserVar.InverseC)+"_C-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"MUA","C","m");
    UserVar.SlidingLaw = CtrlVarInRestartFile.SlidingLaw;
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to copy C field: "+NameOfFiletoRead);
end

%% Glen's exponent
% extract information from the relevant inversion restart file   
NameOfFiletoRead=""; %% ADD CORRECT PATH
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
    AGlen = F.AGlen; n = F.n;
    UserVar.NameOfFileForReadingAGlenEstimate = "ANT_Inverse_"+string(UserVar.InverseA)+"_AGlen-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimate,"MUA","AGlen","n");
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to copy AGlen field: "+NameOfFiletoRead);
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD OTHER STUFF THAT IS RELEVANT FOR TRANSIENT RUNS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

%% write log file
fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Transient: ExpID %s: Starting transient experiment on %s.\n",string(UserVar.ExpID),UserVar.hostname);



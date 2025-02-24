function UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar)


%% target geometry interpolants
%% CALVING
UserVar.Calv = RunTable{ind,"Calv"};
UserVar = ANT_DefineBaseMesh(UserVar,RunTable{ind,"BaseMesh"}{:});
UserVar = ANT_ApplyMeshModifications(UserVar);

%% ICE SHELF and GROUNDED ICE
for geom=["IS","GI"]

    UserVar.(geom+"Geometry") = RunTable{ind,geom+"thick"};
    UserVar.("InverseCycle"+geom) = RunTable{ind,"InverseCycle"+geom};
    
    if UserVar.("InverseCycle"+geom) == 1
        switch UserVar.(geom+"Geometry")
            case {2000,2009,2014,2018,2020}
                UserVar.(geom+"BaseGeometry") = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+...
                    num2str(UserVar.(geom+"Geometry"))+"_EXTRUDED.mat";
                UserVar.("Adapt"+geom+"BaseGeometry") = 0;
            otherwise
                error("ExpID "+RunTable{ind,"ExpID"}+": Do not recognise "+geom+"thick flag in RunTable.");
        end
     
    elseif UserVar.InverseCycleIS > 1 
        UserVar.(geom+"BaseGeometry") = UserVar.casefolder+"../../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(RunTable{ind,"InverseA"})+"/"+UserVar.Domain+...
            "_Inverse_"+string(RunTable{ind,"InverseA"})+"-RestartFile_InverseCycle"+string(UserVar.("InverseCycle"+geom))+".mat";
    
        switch UserVar.(geom+"Geometry")
            case 2000
                UserVar.("Adapt"+geom+"BaseGeometry") = 0;
            case {2009,2014,2018,2020}
                UserVar.("Adapt"+geom+"BaseGeometry") = 1;
                UserVar.("InterpolantsToCalculateNew"+geom+"Geometry") = ...
                    [UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-2000_EXTRUDED.mat",...
                    UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+string(UserVar.(geom+"Geometry"))+"_EXTRUDED.mat"];
            otherwise
                error("ExpID "+RunTable{ind,"ExpID"}+": Do not recognise "+geom+"thick flag in RunTable.");
        end
    
    else
        error("ExpID "+RunTable{ind,"ExpID"}+": Do not recognise UserVar.InverseCycle"+geom+".");
    
    end
end
    
%% density interpolant: same as inverse run
UserVar.DensityInterpolant = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-2000_EXTRUDED.mat";   

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

%% velocity interpolants for fixed boundaries // we assume these do not change in the perturbation experiments
UserVar.Velocity = 2000;
UserVar.VelocityInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat";
    
%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

%% write log file
fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Diagnostic: ExpID %s: Starting diagnostic experiment on %s.\n",string(UserVar.ExpID),UserVar.hostname);



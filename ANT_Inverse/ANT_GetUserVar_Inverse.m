function UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar)

%% This function deals with all the information in the RunTable, and feeds that information 
%% to Ua via the UserVar structure. The content of this function is entirely flexible, and
%% will need to be adapted based on which parameters you need to provide to Ua at runtime.

% initialize variables: either start a new simulation with
% an inverse cycle, or continue an existing simulation with
% a spinup or inverse cycle

iter_tmp = RunTable{ind,"InverseIterations"}{:};
UserVar.Inverse.Iterations = str2double(split(iter_tmp,"+"));
UserVar.Inverse.IterationsDone = RunTable{ind,"InverseIterationsDone"};
UserVar.Inverse.Cycle = find([0; cumsum(UserVar.Inverse.Iterations)]==UserVar.Inverse.IterationsDone);
if isempty(UserVar.Inverse.Cycle) & UserVar.Restart == 0
    fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Inverse: Expecting InverseIterationsDone to be equal to [%s] but got %s instead.\n",...
            string(iter_tmp),string(UserVar.Inverse.IterationsDone));
    error("Unexpected number of InverseIterationsDone.");
elseif isempty(UserVar.Inverse.Cycle) & UserVar.Restart
    %% restarting inverse cycle 
    UserVar.InverseCycle = 1;
    UserVar.SpinupCycle = 0;
    Itmp = find(cumsum(UserVar.Inverse.Iterations)-UserVar.Inverse.IterationsDone>0);
    UserVar.Inverse.Cycle = Itmp(1);
end

spinup_tmp = RunTable{ind,"SpinupYears"}{:};
UserVar.Spinup.Years = str2double(split(spinup_tmp,"+"));
UserVar.Spinup.YearsDone = RunTable{ind,"SpinupYearsDone"};
UserVar.Spinup.Cycle = find([0; cumsum(UserVar.Spinup.Years)]==UserVar.Spinup.YearsDone);
UserVar.Spinup.Restart = 0;
if isempty(UserVar.Spinup.Cycle) & UserVar.Restart == 0
    fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Inverse: Expecting SpinupYearsDone to be equal to [%s] but got %s instead.\n",...
            string(spinup_tmp),string(UserVar.Spinup.YearsDone));
    error("Unexpected number of SpinupYearsDone.");
elseif isempty(UserVar.Spinup.Cycle) & UserVar.Restart == 1
    %% restarting spinup cycle 
    UserVar.Spinup.Restart = 1;
    Itmp = find(cumsum(UserVar.Spinup.Years)-UserVar.Spinup.YearsDone>0);
    UserVar.Spinup.Cycle = Itmp(1);
end

% Restart files
UserVar.NameOfRestartFiletoRead = UserVar.Experiment + "-RestartFile.mat";

% Density
UserVar.Geometry = RunTable{ind,"startGeometry"};
switch UserVar.Geometry
    case {2000,2009,2014,2018}
        UserVar.DensityInterpolant = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.Geometry)+"_EXTRUDED.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Geometry flag in RunTable.']);
end

% Sliding law
UserVar.SlidingLaw = RunTable{ind,"SlidingLaw"}{:};
UserVar.SlidingCoefficient = RunTable{ind,"m"};
UserVar.muk = RunTable{ind,"muk"};
UserVar.NameOfFileForSavingSlipperinessEstimate=string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+"_C-Estimate.mat";

% Rheology
UserVar.n =  RunTable{ind,"n"};
UserVar.NameOfFileForSavingAGlenEstimate=string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+"_AGlen-Estimate.mat";

% Mesh variables
%% !!!!!! This script assumes that the mesh is not adapted during the spinup cycles. !!!!!!!
%% Adapting in the mesh is not advised anyway, because it complicates the diagnostic perturbation experiments
UserVar = ANT_DefineBaseMesh(UserVar,RunTable{ind,"startMesh"}{:});
UserVar = ANT_ApplyMeshModifications(UserVar);

% Outputs
UserVar.UaOutputDirectory = './ResultsFiles';

% Initialize other variables
if UserVar.Inverse.Cycle == 1
    %% (re)start first inverse cycle
    UserVar.InverseCycle = 1;
    UserVar.SpinupCycle = 0;

    fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Inverse: ExpID %s: Start inverse cycle %s on %s.\n",string(UserVar.ExpID),string(UserVar.Inverse.Cycle),UserVar.hostname);

    UserVar = ANT_GetUserVar_FirstInverseRun(RunTable,ind,UserVar);

else
    if UserVar.Inverse.Cycle == UserVar.Spinup.Cycle
        %% start new inverse cycle
        UserVar.InverseCycle = 1;
        UserVar.SpinupCycle = 0;
      
        fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Inverse: ExpID %s: Start inverse cycle %s on %s.\n",string(UserVar.ExpID),string(UserVar.Inverse.Cycle),UserVar.hostname);

        UserVar = ANT_GetUserVar_InverseAfterSpinup(RunTable,ind,UserVar);

    elseif UserVar.Inverse.Cycle > UserVar.Spinup.Cycle
        %% start spinup cycle
        UserVar.SpinupCycle = 1;
        UserVar.InverseCycle = 0;
       
        fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Inverse: ExpID %s: Start spinup cycle %s on %s.\n",string(UserVar.ExpID),string(UserVar.Spinup.Cycle),UserVar.hostname);

        UserVar = ANT_GetUserVar_Spinup(RunTable,ind,UserVar);

    else
        error("Something odd happened: number of inverse cycles (%s) should be >= number of spinup cycles (%s).",...
            string(UserVar.Inverse.Cycle),string(UserVar.Spinup.Cycle));
    end
end


end


%% ======================== %%
%% === helper functions === %%
%% ======================== %%

function UserVar = ANT_GetUserVar_InverseAfterSpinup(RunTable,ind,UserVar)

UserVar.Restart = 0;

%% Read inversion parameters from Runtable
UserVar.Inverse.Measurements = RunTable{ind,"Measurements"}{:};
UserVar.Inverse.GradientCalculation = RunTable{ind,"GradientCalc"}{:};
UserVar.Inverse.InvertFor = RunTable{ind,"InvertFor"}{:};
UserVar.Inverse.logC.gs = RunTable{ind,"gsC"};
UserVar.Inverse.logAGlen.gs = RunTable{ind,"gsA"};
UserVar.Inverse.logC.ga = RunTable{ind,"gaC"};
UserVar.Inverse.logAGlen.ga = RunTable{ind,"gaA"};

%% geometry interpolants from spinup
NameOfRestartFiletoRead = ...
    strrep(UserVar.NameOfRestartFiletoRead,".mat","_SpinupCycle"+string(UserVar.Inverse.Cycle-1)+".mat");
load(UserVar.casefolder+"/"+UserVar.Experiment+"/"+NameOfRestartFiletoRead,"F","MUA");

UserVar.GeometryInterpolants = "ScatteredInterpolants_GeometryfromSpinupCycle"+string(UserVar.Inverse.Cycle-1)+".mat";
% Lines below are removed to save memory. The interpolants should not be
% needed because we don't change the mesh during the spinup
%FB = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.B,"linear");
%Fb = FB; Fb.Values = F.b;
%Fs = FB; Fs.Values = F.s;
%save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.GeometryInterpolants,"FB","Fb","Fs");
%clearvars FB Fb Fs; 

% We create a new file with the Ua geometry fields at the end of the spinup,
% so we can use these fields directly in the next inversion. This wastes a bit of 
% storage space, but for spinup simulation we do not adapt the mesh, so
% it is more memory efficient to read the Ua fields directly, rather than
% load interpolants
B=F.B; b=F.b; s=F.s; S=F.S; rho=F.rho;
save(UserVar.casefolder+"/"+UserVar.Experiment+"/GeometryfromSpinupCycle"+string(UserVar.Inverse.Cycle-1)+"_mesh_Nnodes"+string(MUA.Nnodes)+"_Nele"+string(MUA.Nele)+".mat","B","S","s","b","rho");

%% velocity interpolants
UserVar.Velocity = RunTable{ind,"Velocity"};
switch UserVar.Velocity
    case 2000
        UserVar.VelocityInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat";
    case {2009,2014,2018} 
        UserVar.VelocityInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_"+num2str(UserVar.Velocity)+"-"+num2str(UserVar.Velocity+1)+"_MeaSUREs_ITSLIVE_Velocities.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Velocity flag in RunTable.']);
end

%% dhdt interpolants
if contains(UserVar.Inverse.Measurements,"dhdt")
    UserVar.dhdtInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_dhdt_01-Jun-"+num2str(UserVar.Geometry)+".mat";
    UserVar.dhdt_err = RunTable{ind,"dhdt_err"};
end

%% sliding law
UserVar.NameOfFileForReadingSlipperinessEstimate=UserVar.NameOfFileForSavingSlipperinessEstimate;
UserVar.Inverse.priorC = RunTable{ind,"priorC"};

%% Glen's exponent
UserVar.NameOfFileForReadingAGlenEstimate=UserVar.NameOfFileForSavingAGlenEstimate;
UserVar.Inverse.priorAGlen = RunTable{ind,"priorAGlen"};

end

function UserVar = ANT_GetUserVar_Spinup(RunTable,ind,UserVar)

UserVar.Restart = 1;

% make sure to start from the correct restart file
if UserVar.Spinup.Restart == 0 
    % this means we start a new spinup, using the restart file from the previous inversion
    NameOfRestartInputFile = strrep(UserVar.NameOfRestartFiletoRead,".mat","_InverseCycle"+string(UserVar.Spinup.Cycle)+".mat");
    copyfile(UserVar.casefolder+"/"+string(UserVar.Experiment)+"/"+NameOfRestartInputFile,...
        UserVar.casefolder+"/"+string(UserVar.Experiment)+"/"+UserVar.NameOfRestartFiletoRead);
else
    % we restart the spinup in the current spinup cycle
    NameOfRestartInputFile = strrep(UserVar.NameOfRestartFiletoRead,".mat","_SpinupCycle"+string(UserVar.Spinup.Cycle)+".mat");
    copyfile(UserVar.casefolder+"/"+string(UserVar.Experiment)+"/"+NameOfRestartInputFile,...
        UserVar.casefolder+"/"+string(UserVar.Experiment)+"/"+UserVar.NameOfRestartFiletoRead);
end

%% Read geometry interpolants from Runtable
UserVar.Geometry = RunTable{ind,"startGeometry"};
switch UserVar.Geometry
    case {2000,2009,2014,2018}
        UserVar.GeometryInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.Geometry)+"_EXTRUDED.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Geometry flag in RunTable.']);
end

%% mesh boundary
UserVar.MeshBoundaryCoordinatesFile = UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.Experiment+"_MeshBoundaryCoordinates.mat";

%% sliding law
UserVar.NameOfFileForReadingSlipperinessEstimate=UserVar.NameOfFileForSavingSlipperinessEstimate;

%% Glen's exponent
UserVar.NameOfFileForReadingAGlenEstimate=UserVar.NameOfFileForSavingAGlenEstimate;
    
end


function UserVar = ANT_GetUserVar_FirstInverseRun(RunTable,ind,UserVar)

%% Read inversion parameters from Runtable
UserVar.Inverse.Measurements = "-uv-"; % don't use dhdt in first inverse cycle  
UserVar.Inverse.GradientCalculation = RunTable{ind,"GradientCalc"}{:};
UserVar.Inverse.InvertFor = RunTable{ind,"InvertFor"}{:};
UserVar.Inverse.logC.gs = RunTable{ind,"gsC"};
UserVar.Inverse.logAGlen.gs = RunTable{ind,"gsA"};
UserVar.Inverse.logC.ga = RunTable{ind,"gaC"};
UserVar.Inverse.logAGlen.ga = RunTable{ind,"gaA"};

%% Read geometry interpolants from Runtable
UserVar.Geometry = RunTable{ind,"startGeometry"};
switch UserVar.Geometry
    case {2000,2009,2014,2018}
        UserVar.GeometryInterpolants = UserVar.datafolder+"ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+num2str(UserVar.Geometry)+"_EXTRUDED.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Geometry flag in RunTable.']);
end

%% density interpolant
UserVar.DensityInterpolant = UserVar.GeometryInterpolants;

%% velocity interpolants
UserVar.Velocity = RunTable{ind,"Velocity"};
switch UserVar.Velocity
    case 2000
        UserVar.VelocityInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat";
    case {2009,2014,2018} 
        UserVar.VelocityInterpolants = UserVar.datafolder+"/ANT_Interpolants/GriddedInterpolants_"+num2str(UserVar.Velocity)+"-"+num2str(UserVar.Velocity+1)+"_MeaSUREs_ITSLIVE_Velocities.mat";
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Velocity flag in RunTable.']);
end

%% sliding law
UserVar.Inverse.startC = 0;
if RunTable{ind,"startC"} > 0
    % copy relevant files   
    RestartFile = UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startC"})+...
        "/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startC"})+"-RestartFile_InverseCycle1.mat";
    UserVar.NameOfFileForReadingSlipperinessEstimate=string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startC"})+"_C-Estimate.mat";
    if exist(RestartFile,"file")
        load(RestartFile,"F","MUA","CtrlVarInRestartFile");
        C = F.C; xC = MUA.coordinates(:,1); yC = MUA.coordinates(:,2); m = F.m;
        muk = F.muk; q = F.q;
        save(UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"MUA","CtrlVarInRestartFile","xC","yC","C","m","muk","q");
    else    
        copyfile(UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startC"})+...
                "/"+UserVar.NameOfFileForReadingSlipperinessEstimate,...
                UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+...
                "/"+UserVar.NameOfFileForReadingSlipperinessEstimate);
    end
elseif RunTable{ind,"startC"} == -9999
    UserVar.Inverse.startC = -9999;
    UserVar.NameOfFileForReadingSlipperinessEstimate="";
else
    UserVar.NameOfFileForReadingSlipperinessEstimate="";
end
UserVar.Inverse.priorC = RunTable{ind,"priorC"};


%% Glen's exponent
if RunTable{ind,"startAglen"} > 0
    % copy relevant files   
    RestartFile = UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startAglen"})+...
        "/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startAglen"})+"-RestartFile_InverseCycle1.mat";
    UserVar.NameOfFileForReadingAGlenEstimate=string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startAglen"})+"_AGlen-Estimate.mat";
    if exist(RestartFile,"file")
        load(RestartFile,"F","MUA","CtrlVarInRestartFile");
        AGlen = F.AGlen; xA = MUA.coordinates(:,1); yA = MUA.coordinates(:,2); n = F.n;
        save(UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate,"MUA","CtrlVarInRestartFile","xA","yA","AGlen","n");
    else    
        copyfile(UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startAglen"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate,...
            UserVar.casefolder+"/"+string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate);
    end 
elseif RunTable{ind,"startAGlen"} == -9999
    UserVar.Inverse.startAGlen = -9999;
    UserVar.NameOfFileForReadingAGlenEstimate="";  
else
    UserVar.NameOfFileForReadingAGlenEstimate="";
end
UserVar.Inverse.priorAGlen = RunTable{ind,"priorAGlen"};


end

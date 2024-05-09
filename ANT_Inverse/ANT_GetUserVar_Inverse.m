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
    fprintf(UserVar.fid,"Expecting InverseIterationsDone to be equal to [%s] but got %s instead.\n",...
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

% check stopping criteria for inverse run
%if UserVar.InverseCycle
%    UserVar.Inverse.MinGradNorm = RunTable{ind,"MinGradNorm"};
%end

% Restart files
UserVar.NameOfRestartFiletoRead = UserVar.Experiment + "-RestartFile.mat";

% Density
UserVar.Geometry = RunTable{ind,"startGeometry"};
switch UserVar.Geometry
    case {2000,2009,2014,2018}
        UserVar.DensityInterpolant = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.Geometry),'_EXTRUDED.mat'];
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
    %% start first inverse cycle
    UserVar.InverseCycle = 1;
    UserVar.SpinupCycle = 0;

    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,"> %s: Start inverse cycle %s.\n",UserVar.Experiment,string(UserVar.Inverse.Cycle));

    UserVar = ANT_GetUserVar_FirstInverseRun(RunTable,ind,UserVar);

else
    if UserVar.Inverse.Cycle == UserVar.Spinup.Cycle
        %% start new inverse cycle
        UserVar.InverseCycle = 1;
        UserVar.SpinupCycle = 0;

        fprintf(UserVar.fid,'============================\n');
        fprintf(UserVar.fid,string(datetime("now"))+"\n");
        fprintf(UserVar.fid,'============================\n');
        fprintf(UserVar.fid,"> %s: Start inverse cycle %s.\n",UserVar.Experiment,string(UserVar.Inverse.Cycle));

        UserVar = ANT_GetUserVar_InverseAfterSpinup(RunTable,ind,UserVar);

    elseif UserVar.Inverse.Cycle > UserVar.Spinup.Cycle
        %% start spinup cycle
        UserVar.SpinupCycle = 1;
        UserVar.InverseCycle = 0;

        fprintf(UserVar.fid,'============================\n');
        fprintf(UserVar.fid,string(datetime("now"))+"\n");
        fprintf(UserVar.fid,'============================\n');
        fprintf(UserVar.fid,"> %s: Start spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
        
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
load("./"+UserVar.Experiment+"/"+NameOfRestartFiletoRead,"F","MUA");
FB = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.B,"linear");
Fb = FB; Fb.Values = F.b;
Fs = FB; Fs.Values = F.s;
UserVar.GeometryInterpolants = "GeometryInterpolants_fromSpinupCycle"+string(UserVar.Inverse.Cycle-1)+".mat";
save("./"+UserVar.Experiment+"/"+UserVar.GeometryInterpolants,"FB","Fb","Fs");

%% velocity interpolants
UserVar.Velocity = RunTable{ind,"Velocity"};
switch UserVar.Velocity
    case 2000
        UserVar.VelocityInterpolants = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat'];
    case {2009,2014,2018} 
        UserVar.VelocityInterpolants = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_',num2str(UserVar.Velocity),'-',num2str(UserVar.Velocity+1),'_MeaSUREs_ITSLIVE_Velocities.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Velocity flag in RunTable.']);
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

% copy and rename inverse restart file
%UserVar.Spinup.NameOfRestartFiletoRead = ...
%    strrep(UserVar.Inverse.NameOfRestartInputFile,"-RestartFile","_SpinupCycle"+string(UserVar.Spinup.Cycle)+"-RestartFile");
%copyfile("./"+string(UserVar.Experiment)+"/"+UserVar.Inverse.NameOfRestartInputFile,...
%    "./"+string(UserVar.Experiment)+"/"+UserVar.Spinup.NameOfRestartFiletoRead);

%% Read geometry interpolants from Runtable
UserVar.Geometry = RunTable{ind,"startGeometry"};
switch UserVar.Geometry
    case {2000,2009,2014,2018}
        UserVar.GeometryInterpolants = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.Geometry),'_EXTRUDED.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Geometry flag in RunTable.']);
end

%% mesh boundary
UserVar.MeshBoundaryCoordinatesFile = "./"+UserVar.Experiment+"_MeshBoundaryCoordinates.mat";

%% sliding law
UserVar.NameOfFileForReadingSlipperinessEstimate=UserVar.NameOfFileForSavingSlipperinessEstimate;

%% Glen's exponent
UserVar.NameOfFileForReadingAGlenEstimate=UserVar.NameOfFileForSavingAGlenEstimate;
    
end


function UserVar = ANT_GetUserVar_FirstInverseRun(RunTable,ind,UserVar)

%% Read inversion parameters from Runtable
UserVar.Inverse.Measurements = RunTable{ind,"Measurements"}{:};
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
        UserVar.GeometryInterpolants = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.Geometry),'_EXTRUDED.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Geometry flag in RunTable.']);
end

%% density interpolant
UserVar.DensityInterpolant = UserVar.GeometryInterpolants;

%% velocity interpolants
UserVar.Velocity = RunTable{ind,"Velocity"};
switch UserVar.Velocity
    case 2000
        UserVar.VelocityInterpolants = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat'];
    case {2009,2014,2018} 
        UserVar.VelocityInterpolants = [pwd,'/../ANT_Data/ANT_Interpolants/GriddedInterpolants_',num2str(UserVar.Velocity),'-',num2str(UserVar.Velocity+1),'_MeaSUREs_ITSLIVE_Velocities.mat'];
    otherwise
        error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise Velocity flag in RunTable.']);
end

%% sliding law
if RunTable{ind,"startC"} ~= 0
    % copy relevant files   
    UserVar.NameOfFileForReadingSlipperinessEstimate=string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startC"})+"_C-Estimate.mat";
    copyfile(string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startC"})+...
            "/"+UserVar.NameOfFileForReadingSlipperinessEstimate,...
            string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingSlipperinessEstimate);
else
    UserVar.NameOfFileForReadingSlipperinessEstimate="";
end
UserVar.Inverse.priorC = RunTable{ind,"priorC"};


%% Glen's exponent
if RunTable{ind,"startAglen"} ~= 0
    % copy relevant files   
    UserVar.NameOfFileForReadingAGlenEstimate=string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startAglen"})+"_AGlen-Estimate.mat";
    copyfile(string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"startAglen"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate,...
            string(UserVar.Domain)+"_Inverse_"+string(RunTable{ind,"ExpID"})+...
            "/"+UserVar.NameOfFileForReadingAGlenEstimate);
else
    UserVar.NameOfFileForReadingAGlenEstimate="";
end
UserVar.Inverse.priorAGlen = RunTable{ind,"priorAGlen"};


end

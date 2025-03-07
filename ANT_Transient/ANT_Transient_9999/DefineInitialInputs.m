function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)

%% DEBUG OPTIONS ARE AT THE END OF THIS SCRIPT %%

%% Keep track of walltime
% These lines define a UI to keep track of the remaining walltime.
% If the walltime expires, the UI is set to false, which is picked
% up by the fmincon minimization algorithm or DefineRunStopCriterion.m and 
% used as a stopping % criteria to break out of the inversion or spinup.
% We apply a generous 60min buffer to allow fmincon and the runstep to 
% cleanly finish the current iteration. The remaining walltime should never 
% be negative, so we set it to a minimum value of 600s.
UserVar.walltime_remaining = max(UserVar.walltime_remaining-3600,600);
setappdata(0,'UAstopFlag',false); %stopping flag is false
T = timer('startdelay',UserVar.walltime_remaining,'timerfcn',@(src,evt)setappdata(0,'UAstopFlag',true)); %initialize timer to change value of uastopflag after wallclocktime
t0 = tic(); 
start(T); %start the timer
remainingTime = round(UserVar.walltime_remaining-toc(t0));
fprintf(UserVar.fid_experimentlog,"> At %s: remaining time on wallclock timer is %ss. Ua will be stopped when this time has been exceeded.\n",string(datetime("now")),num2str(remainingTime));

%%
CtrlVar.Experiment = UserVar.Experiment;

%% Type of run
%
CtrlVar.TimeDependentRun=1; 
CtrlVar.InverseRun=0;
CtrlVar.Restart=UserVar.Restart;
%
CtrlVar.StartTime=UserVar.StartTime_DecimalYears ; % start time (decimal year, e.g. 2000.5)
CtrlVar.EndTime=UserVar.EndTime_DecimalYears; % end time (decimal year, e.g. 2001.5)
CtrlVar.dt=0.001; % time step
CtrlVar.TotalNumberOfForwardRunSteps=1e10;  
%
UserVar.stoppedduetowalltime = 0;
CtrlVar.UseUserDefinedRunStopCriterion=1;

%% Grid options
CtrlVar.TriNodes=3;
CtrlVar.kH=10;
CtrlVar.nip=6;
CtrlVar.niph=6;
CtrlVar.ReadInitialMesh=1;
CtrlVar.ReadInitialMeshFileName=UserVar.InitialMeshFileName;
CtrlVar.RefineMeshOnStart=0;

%% Define domain boundary
if exist(UserVar.MeshBoundaryCoordinatesFile,"file")
    load(UserVar.MeshBoundaryCoordinatesFile,"MeshBoundaryCoordinates");
elseif exist(CtrlVar.ReadInitialMeshFileName,"file")
    tmp=load(CtrlVar.ReadInitialMeshFileName,"MUA");
    MeshBoundaryCoordinates = [tmp.MUA.Boundary.x(:) tmp.MUA.Boundary.y(:)];
    save(UserVar.MeshBoundaryCoordinatesFile,"MeshBoundaryCoordinates");
else
    error("Do not know what the domain boundary is. Could not find "+UserVar.MeshBoundaryCoordinatesFile+...
        " or "+CtrlVar.ReadInitialMeshFileName);
end

if UserVar.AdaptMesh == 1

    CtrlVar.InfoLevelAdaptiveMeshing=1;  
    CtrlVar.MeshSize=100e3;
    CtrlVar.MeshSizeMax=100e3;
    CtrlVar.MeshSizeMin=1.5e3;

    % CtrlVar.MeshAdapt.GLrange=[8000 4000 ; 3000  CtrlVar.MeshSizeMin];% (10000 range from GL) (5000 elm size within the range) ; (3000 set another range) (and another size)
    CtrlVar.MeshAdapt.GLrange=[10000 5000 ; 3000  CtrlVar.MeshSizeMin];
    CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed

    CtrlVar.AdaptMesh=1;
    CtrlVar.AdaptMeshMaxIterations=10;
    CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';

    CtrlVar.SaveAdaptMeshFileName='MeshFileAdapt';    %  file name for saving adapt mesh. If left empty, no file is written
    CtrlVar.AdaptMeshRunStepInterval=50 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0
    %CtrlVar.doAdaptMeshPlots=1; 
end

%% Physics
CtrlVar.SlidingLaw = UserVar.SlidingLaw;

%% Output options
CtrlVar.InfoLevelNonLinIt=1;
CtrlVar.InfoLevel=1;

if CtrlVar.Restart
    CtrlVar.CreateOutputsBeginningOfRun=0; 
    % capture outputs at CtrlVar.StartTime, but otherwise do not write
    % outputs at the start of a run
end
CtrlVar.CreateOutputsEndOfRun=0; % don't write outputs at the end of a run

CtrlVar.DefineOutputsDt=1/12; % write monthly output files

CtrlVar.NameOfRestartFiletoWrite=CtrlVar.Experiment+"-RestartFile.mat";

%% Plotting
CtrlVar.doplots=0;
CtrlVar.PlotWaitBar=0;
CtrlVar.PlotOceanLakeNodes=0;
CtrlVar.PlotMesh=0;  CtrlVar.PlotBCs=1;
CtrlVar.PlotXYscale=1;

%CtrlVar.AGlenmin=AGlenVersusTemp(-15)/1e4;
%CtrlVar.AGlenmax=AGlenVersusTemp(-15)*1e4;
CtrlVar.Cmin=1e-150;  CtrlVar.Cmax=1e150;                                                  

%% Other
CtrlVar.ThicknessConstraints=0;
CtrlVar.ResetThicknessToMinThickness=1;  % change this later on
CtrlVar.ThickMin=1;

%% Debugging options
if UserVar.debug == 1
    %CtrlVar.ExplicitEstimationMethod="-dhdt-" ;
    %CtrlVar.GuardAgainstWildExtrapolationInExplicit_uvh_Step=1;
    %CtrlVar.ATSdtMin=1e-10;
    %CtrlVar.NRitmax=25;
    CtrlVar.InfoLevel=100;
    CtrlVar.InfoLevelBackTrack=10;
    CtrlVar.InfoLevelNonLinIt=100;
    CtrlVar.doplots=1;   
    CtrlVar.PlotGLs=1;
    %CtrlVar.ThickMin=1;
end
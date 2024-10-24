function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)

CtrlVar.Experiment = UserVar.Experiment;

%% Type of run
%
CtrlVar.TimeDependentRun=1; 
CtrlVar.InverseRun=0;

CtrlVar.TotalNumberOfForwardRunSteps=inf;

CtrlVar.Restart=0;%UserVar.Restart;

CtrlVar.time=0 ; % start time
CtrlVar.dt=0.001; % time step
CtrlVar.TotalNumberOfForwardRunSteps=1e10;  
CtrlVar.TotalTime=1; % in years

%% Grid options
CtrlVar.TriNodes=3;
CtrlVar.kH=10;
CtrlVar.nip=6;
CtrlVar.niph=6;
CtrlVar.ReadInitialMesh=1;
CtrlVar.ReadInitialMeshFileName=UserVar.InitialMeshFileName;
CtrlVar.RefineMeshOnStart=0;

if UserVar.AdaptMesh ==1

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
    CtrlVar.AdaptMeshRunStepInterval=200 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0
    %CtrlVar.doAdaptMeshPlots=1; 
end


%% Boundary
load(UserVar.MeshBoundaryCoordinatesFile,"MeshBoundaryCoordinates");

%% Physics
CtrlVar.SlidingLaw = UserVar.SlidingLaw;

%% Output options
CtrlVar.InfoLevelNonLinIt=1;
CtrlVar.InfoLevel=1;

CtrlVar.CreateOutputsEndOfRun=1;

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
CtrlVar.ThicknessConstraints=1;
CtrlVar.ResetThicknessToMinThickness=0;  % change this later on
CtrlVar.ThickMin=1;

end

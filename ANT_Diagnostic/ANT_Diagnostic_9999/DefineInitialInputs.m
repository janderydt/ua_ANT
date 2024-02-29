function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)

CtrlVar.Experiment = UserVar.Experiment;

%% Type of run
%
CtrlVar.TimeDependentRun=0; 
CtrlVar.InverseRun=0;

CtrlVar.TotalNumberOfForwardRunSteps=1;

CtrlVar.Restart=0;%UserVar.Restart;

%% Grid options
CtrlVar.TriNodes=3;
CtrlVar.kH=1;
CtrlVar.nip=6;
CtrlVar.niph=6;
CtrlVar.ReadInitialMesh=1;
CtrlVar.ReadInitialMeshFileName=UserVar.InitialMeshFileName;
CtrlVar.RefineMeshOnStart=0;
CtrlVar.AdaptMesh=0;

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

CtrlVar.AGlenmin=AGlenVersusTemp(-15)/1e4;
CtrlVar.AGlenmax=AGlenVersusTemp(-15)*1e4;
CtrlVar.Cmin=1e-50;  CtrlVar.Cmax=1e10;                                                  

%% Other
CtrlVar.ThicknessConstraints=1;
CtrlVar.ResetThicknessToMinThickness=0;  % change this later on
CtrlVar.ThickMin=1;

end

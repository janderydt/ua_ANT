function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)

CtrlVar.Experiment = UserVar.Experiment;

%% Type of run
%
CtrlVar.TimeDependentRun=0; 
CtrlVar.doInverseStep=1;
CtrlVar.TotalNumberOfForwardRunSteps=inf; % an arbitrary large number
CtrlVar.TotalTime=1;
CtrlVar.time=0;
CtrlVar.Restart=UserVar.Restart;

CtrlVar.dt = 1e-3;
CtrlVar.RestartTime=0; 
CtrlVar.ResetTime=0;
CtrlVar.ResetTimeStep=0;    % perhaps this has to be reconsidered if model has issues converging

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
CtrlVar.InfoLevelInverse=1; % Overall level of information (inverse runs). 
                            % Note: generally good to combine with CtrlVar.InfoLevelNonLinIt=0;
                            % CtrlVar.InfoLevel=0; to suppress information related to the forward step. 
CtrlVar.Inverse.InfoLevelBackTrack=1000;  % info on backtracking within inverse step
CtrlVar.InfoLevelNonLinIt=0;
CtrlVar.InfoLevel=0;

CtrlVar.Inverse.WriteRestartFile=1;  % always a good idea to write a restart file. 
CtrlVar.Inverse.NameOfRestartOutputFile= [CtrlVar.Experiment,'-RestartFile.mat'];
CtrlVar.Inverse.NameOfRestartInputFile = CtrlVar.Inverse.NameOfRestartOutputFile;

CtrlVar.Inverse.SaveSlipperinessEstimateInSeperateFile=true ;
CtrlVar.Inverse.SaveAGlenEstimateInSeperateFile=true ;
CtrlVar.NameOfFileForSavingSlipperinessEstimate=UserVar.NameOfFileForSavingSlipperinessEstimate;
CtrlVar.NameOfFileForSavingAGlenEstimate=UserVar.NameOfFileForSavingAGlenEstimate;

CtrlVar.CreateOutputsEndOfRun=1;

%% Plotting
CtrlVar.doplots=0;
CtrlVar.PlotWaitBar=0;
CtrlVar.PlotOceanLakeNodes=0;
CtrlVar.PlotMesh=0;  CtrlVar.PlotBCs=1;
CtrlVar.PlotXYscale=1;

%% Inverse options
CtrlVar.Inverse.Iterations=UserVar.Iterations; % Number of inverse iterations

CtrlVar.Inverse.DataMisfit.GradientCalculation=UserVar.Inverse.GradientCalculation; % {'Adjoint','FixPoint'}

CtrlVar.Inverse.Measurements=UserVar.Inverse.Measurements;

CtrlVar.Inverse.InvertFor=UserVar.Inverse.InvertFor; % {'-C-','-logC-','-A-','-logA-'}
%CtrlVar.Inverse.InvertFor='-logC-' ;
CtrlVar.Inverse.Regularize.Field=CtrlVar.Inverse.InvertFor; 

CtrlVar.Inverse.DataMisfit.Multiplier=1;
CtrlVar.Inverse.Regularize.Multiplier=1;

% regularisation parameters
CtrlVar.Inverse.Regularize.logAGlen.ga=UserVar.Inverse.logAGlen.ga;
CtrlVar.Inverse.Regularize.logAGlen.gs=UserVar.Inverse.logAGlen.gs ;
CtrlVar.Inverse.Regularize.logC.ga=UserVar.Inverse.logC.ga;
CtrlVar.Inverse.Regularize.logC.gs=UserVar.Inverse.logC.gs; 

CtrlVar.AGlenmin=AGlenVersusTemp(-15)/1e4;
CtrlVar.AGlenmax=AGlenVersusTemp(-15)*1e4;
CtrlVar.Cmin=1e-50;  CtrlVar.Cmax=1e10;                                                  

%% Other
CtrlVar.ThicknessConstraints=1;
CtrlVar.ResetThicknessToMinThickness=0;  % change this later on
CtrlVar.ThickMin=1;

end

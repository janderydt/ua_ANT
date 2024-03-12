function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)

CtrlVar.Experiment = UserVar.Experiment;
CtrlVar.Restart=UserVar.Restart;

%% Type of run
if UserVar.InverseCycle
    CtrlVar.TimeDependentRun=0; 
    CtrlVar.doInverseStep=1;
elseif UserVar.SpinupCycle
    CtrlVar.TimeDependentRun=1; 
    CtrlVar.doInverseStep=0;
    CtrlVar.TotalNumberOfForwardRunSteps=inf; % an arbitrary large number
    CtrlVar.TotalTime=UserVar.Spinup.Years(UserVar.Spinup.Cycle);
    CtrlVar.time=0;
    CtrlVar.dt = 1e-3;
    CtrlVar.RestartTime=0; 
    CtrlVar.ResetTime=1;
    CtrlVar.ResetTimeStep=1;    % perhaps this has to be reconsidered if model has issues converging
    CtrlVar.InitialDiagnosticStep=1;    
else
    error("Unknown run case");
end

%% Grid options
CtrlVar.TriNodes=3;
CtrlVar.kH=1;
CtrlVar.nip=6;
CtrlVar.niph=6;    
CtrlVar.RefineMeshOnStart=0;

% always read intial mesh
CtrlVar.ReadInitialMesh=1;
CtrlVar.ReadInitialMeshFileName=UserVar.InitialMeshFileName;

% refine mesh around GL at first spinup cycle 
if UserVar.SpinupCycle && UserVar.Spinup.Cycle==1    
    CtrlVar.AdaptMesh=1;
else
    CtrlVar.AdaptMesh=0;
end

%% Boundary
load(UserVar.MeshBoundaryCoordinatesFile,"MeshBoundaryCoordinates");

%% Physics
CtrlVar.SlidingLaw = UserVar.SlidingLaw;

%% Adapt Mesh options
if UserVar.SpinupCycle
    CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';
    
    CtrlVar.MeshSizeMax=100e3;
    CtrlVar.MeshSize=100e3;
    CtrlVar.MeshSizeMin=1e3;
    
    CtrlVar.MeshAdapt.GLrange=[5e3 CtrlVar.MeshSizeMin*2; 2e3  CtrlVar.MeshSizeMin];
    
    CtrlVar.RefineMeshOnStart=0;
    CtrlVar.InfoLevelAdaptiveMeshing=1;                                            
    CtrlVar.AdaptMeshInitial=0  ; % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshInterval)~=0.
    CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                       % useful, for example, when trying out different remeshing options (then use CtrlVar.doAdaptMeshPlots=1 to get plots)
    CtrlVar.AdaptMeshMaxIterations=10;
    CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan = 20;
    CtrlVar.SaveAdaptMeshFileName='MeshFileAdapt';    %  file name for saving adapt mesh. If left empty, no file is written
    CtrlVar.AdaptMeshRunStepInterval=1e20 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0
    CtrlVar.doAdaptMeshPlots=1; 
end

%% Output options
if UserVar.InverseCycle
    CtrlVar.InfoLevelInverse=1; % Overall level of information (inverse runs). 
                                % Note: generally good to combine with CtrlVar.InfoLevelNonLinIt=0;
                                % CtrlVar.InfoLevel=0; to suppress information related to the forward step. 
    CtrlVar.Inverse.InfoLevelBackTrack=1000;  % info on backtracking within inverse step
    CtrlVar.InfoLevelNonLinIt=0;
    CtrlVar.InfoLevel=0;
    CtrlVar.Inverse.WriteRestartFile=1;  % always a good idea to write a restart file. 
    CtrlVar.Inverse.NameOfRestartInputFile = UserVar.Inverse.NameOfRestartInputFile;
    CtrlVar.Inverse.NameOfRestartOutputFile = CtrlVar.Inverse.NameOfRestartInputFile;
    CtrlVar.Inverse.SaveSlipperinessEstimateInSeperateFile=true;
    CtrlVar.Inverse.SaveAGlenEstimateInSeperateFile=true;
    CtrlVar.NameOfFileForSavingSlipperinessEstimate=UserVar.NameOfFileForSavingSlipperinessEstimate;
    CtrlVar.NameOfFileForSavingAGlenEstimate=UserVar.NameOfFileForSavingAGlenEstimate;
else
    CtrlVar.InfoLevelNonLinIt=1;
    CtrlVar.InfoLevel=1;
    CtrlVar.DefineOutputsDt=1;
    CtrlVar.WriteRestartFile = 1;
    CtrlVar.WriteRestartFileInterval = 100;
    CtrlVar.NameOfRestartFiletoWrite = UserVar.Spinup.NameOfRestartFiletoRead;
end
CtrlVar.CreateOutputsEndOfRun=1;

%% Plotting
CtrlVar.doplots=0;
CtrlVar.PlotWaitBar=0;
CtrlVar.PlotOceanLakeNodes=0;
CtrlVar.PlotMesh=0;  CtrlVar.PlotBCs=1;
CtrlVar.PlotXYscale=1;

%% Inverse options
if UserVar.InverseCycle
    CtrlVar.Inverse.Iterations=UserVar.TargetIterations; % Number of inverse iterations
    
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
end

CtrlVar.AGlenmin=AGlenVersusTemp(-15)/1e4;
CtrlVar.AGlenmax=AGlenVersusTemp(-15)*1e4;
CtrlVar.Cmin=1e-50;  CtrlVar.Cmax=1e10;                                                  

%% Minimum ice thickness
CtrlVar.ThicknessConstraints=0;
CtrlVar.ResetThicknessToMinThickness=1;  % change this later on
CtrlVar.ThickMin=1;

end

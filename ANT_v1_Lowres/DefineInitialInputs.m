function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)

UserVar.GeometryInterpolant='./Interpolants/Bedmap2GriddedInterpolantModifiedBathymetry.mat'; % this assumes you have downloaded the OneDrive folder `Interpolants'.
UserVar.DensityInterpolant='./Interpolants/DepthAveragedDensityGriddedInterpolant.mat';
UserVar.SurfaceVelocityInterpolant='./Interpolants/SurfVelMeasures990mInterpolants.mat';
UserVar.UaOutputDirectory = './ResultsFiles';

CtrlVar.Experiment = 'Antarctica_40km_RACMO_BalancedMelt';

%% Type of run
%
CtrlVar.TimeDependentRun=1; 
CtrlVar.doInverseStep=0;
CtrlVar.TotalNumberOfForwardRunSteps=inf; % an arbitrary large number
CtrlVar.TotalTime=10010;
CtrlVar.time=0;
CtrlVar.Restart=1;

CtrlVar.dt = 1e-3;
CtrlVar.RestartTime=0; 
CtrlVar.ResetTime=0;
CtrlVar.ResetTimeStep=0;    % perhaps this has to be reconsidered if model has issues converging

% Parallel options
%myCluster = parcluster('local') ;  
%myCluster.NumWorkers = 6;
%saveProfile(myCluster)

CtrlVar.Parallel.uvhAssembly.parfor.isOn=0;     % assembly over integration points done in parallel using parfor
CtrlVar.Parallel.uvhAssembly.spmd.isOn=0;       % assembly in parallel using spmd over sub-domain (domain decomposition)  
CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=[];

% Grid options
CtrlVar.TriNodes=3;
CtrlVar.kH=1;
CtrlVar.nip=6;
CtrlVar.niph=6;
CtrlVar.AdaptMesh=0;

% timestepping
%CtrlVar.ATStimeStepTarget = UserVar.UaMITgcm.ATStimeStepTarget; 
CtrlVar.ATSdtMax = 10;
CtrlVar.dtmin = 1e-10;
CtrlVar.ATStimeStepFactorUp=2 ;
CtrlVar.ATStimeStepFactorDown=1e2 ;
CtrlVar.ATSTargetIterations=3;

CtrlVar.InitialDiagnosticStep=1;
CtrlVar.InitialDiagnosticStepAfterRemeshing=0;
CtrlVar.Implicituvh=1;
CtrlVar.TG3=0 ; CtrlVar.Gamma=1;
CtrlVar.uvhTimeSteppingMethod='supg';

load BoundaryCoordinates MeshBoundaryCoordinates

CtrlVar.DefineOutputsDt = 10;

CtrlVar.WriteRestartFile=1;
CtrlVar.WriteRestartFileInterval=50;
CtrlVar.NameOfRestartFiletoWrite = [CtrlVar.Experiment,'-RestartFile.mat'];
CtrlVar.NameOfRestartFiletoRead = CtrlVar.NameOfRestartFiletoWrite;

CtrlVar.NameOfFileForReadingSlipperinessEstimate = 'Antarctica_40km_AGlen-Estimate.mat';
CtrlVar.NameOfFileForReadingAGlenEstimate = 'Antarctica_40km_C-Estimate.mat';
CtrlVar.NameOfFileForSavingSlipperinessEstimate='';
CtrlVar.NameOfFileForSavingAGlenEstimate='';

CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry='SB';

CtrlVar.doplots=0;
CtrlVar.PlotWaitBar=0;
CtrlVar.PlotOceanLakeNodes=0;
CtrlVar.PlotMesh=0;  CtrlVar.PlotBCs=1;

CtrlVar.MeltNodesDefinition='edge-wise';
CtrlVar.MassBalanceGeometryFeedback = 0;
CtrlVar.MeltRateFactor=1;
CtrlVar.MeltReductionTime=Inf;

CtrlVar.MeshSizeMax=100e3;
CtrlVar.MeshSize=20e3;
CtrlVar.MeshSizeMin=CtrlVar.MeshSize/5;
CtrlVar.MeshSizeBoundary=CtrlVar.MeshSize;

CtrlVar.ReadInitialMesh=1;
CtrlVar.ReadInitialMeshFileName='./Antarctica_ref_nEle168146.mat';
%CtrlVar.ReadInitialMeshFileName='./Antarctica_40km_nEle20177.mat';

CtrlVar.OnlyMeshDomainAndThenStop=0; % if true then only meshing is done and no further calculations. Usefull for checking if mesh is reasonable
CtrlVar.MaxNumberOfElements=150e4;

%% plotting
CtrlVar.PlotXYscale=1;

%%
CtrlVar.InfoLevelNonLinIt=1;

%% adapt mesh
%CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';
CtrlVar.MeshRefinementMethod='explicit:global';    % can have any of these values:
                                                   % 'explicit:global'
                                                   % 'explicit:local'
                                                   % 'implicit:global'  (broken at the moment, do not use)
                                                   % 'implicit:local'   (broken at the moment, do not use)

      
I=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='effective strain rates';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=1e-3;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1000;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='thickness gradient';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=1e-2;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1000;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='effective strain rates gradient';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=5e-7;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[0.7];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1000;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=false;
                                                   
CtrlVar.MeshAdapt.GLrange=[2000 CtrlVar.MeshSizeMin];

CtrlVar.RefineMeshOnStart=0;
CtrlVar.InfoLevelAdaptiveMeshing=1;                                            
CtrlVar.AdaptMeshInitial=1  ; % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshInterval)~=0.
CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                   % useful, for example, when trying out different remeshing options (then use CtrlVar.doAdaptMeshPlots=1 to get plots)
CtrlVar.AdaptMeshMaxIterations=5;
CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan = 20;
CtrlVar.SaveAdaptMeshFileName='MeshFileAdapt';    %  file name for saving adapt mesh. If left empty, no file is written
CtrlVar.AdaptMeshRunStepInterval=20 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0
CtrlVar.doAdaptMeshPlots=0; 

CtrlVar.Inverse.MinimisationMethod='MatlabOptimization'; % {,'UaOptimization'}
%CtrlVar.Inverse.MinimisationMethod='UaOptimization';
CtrlVar.CisElementBased=0;   
CtrlVar.AGlenisElementBased=0;
%CtrlVar.Inverse.CalcGradI=true;
%CtrlVar.AGlenmin=AGlenVersusTemp(-15)/1e4;
%CtrlVar.AGlenmax=AGlenVersusTemp(-15)*1e4;
CtrlVar.Cmin=1e-50;  CtrlVar.Cmax=1e20;      

CtrlVar.InitTrustRegionRadius = 1;
CtrlVar.InitBarrierParam = 1e-7;
CtrlVar.Inverse.Iterations = 8900;

% CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
%     'Algorithm','quasi-newton',...
%     'MaxIterations',CtrlVar.Inverse.Iterations,...
%     'MaxFunctionEvaluations',1000,...
%     'Display','iter-detailed',...
%     'OutputFcn',@fminuncOutfun,...
%     'Diagnostics','on',...
%     'OptimalityTolerance',1e-20,...
%     'StepTolerance',1e-20,...
%     'PlotFcn',{@optimplotfval,@optimplotstepsize},...
%     'SpecifyObjectiveGradient',CtrlVar.Inverse.CalcGradI);

CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        'CheckGradients',false,...
        'ConstraintTolerance',1e-10,...
        'HonorBounds',true,...
        'Diagnostics','on',...
        'DiffMaxChange',Inf,...
        'DiffMinChange',0,...
        'Display','iter-detailed',...
        'FunValCheck','off',...
        'MaxFunctionEvaluations',1e6,...
        'MaxIterations',CtrlVar.Inverse.Iterations,...,...
        'OptimalityTolerance',1e-10,...
        'OutputFcn',@fminuncOutfun,...
        'PlotFcn',{@optimplotlogfval,@optimplotstepsize},...
        'StepTolerance',1e-10,...
        'FunctionTolerance',1,...
        'UseParallel',true,...
        'HessianApproximation',{'lbfgs',30},...
        'HessianFcn',[],...
        'HessianMultiplyFcn',[],...
        'InitBarrierParam',CtrlVar.InitBarrierParam,...           % On a restart this might have to be reduced if objective function starts to increase
        'ScaleProblem','none',...
        'InitTrustRegionRadius',CtrlVar.InitTrustRegionRadius,...         % default:1, set to smaller value if the forward problem is not converging
        'SpecifyConstraintGradient',false,...
        'SpecifyObjectiveGradient',true,...
        'SubproblemAlgorithm','cg'); % 'cg' or 'factorization'

CtrlVar.Inverse.WriteRestartFile=1;  % always a good idea to write a restart file. 
CtrlVar.Inverse.NameOfRestartOutputFile=[CtrlVar.Experiment,'_InverseRestartFile.mat'];
CtrlVar.Inverse.NameOfRestartInputFile=[CtrlVar.Experiment,'_InverseRestartFile.mat'];
CtrlVar.NameOfFileForReadingSlipperinessEstimate=[CtrlVar.Experiment,'_C-Estimate.mat'];
CtrlVar.NameOfFileForReadingAGlenEstimate=[CtrlVar.Experiment,'_AGlen-Estimate.mat'];
CtrlVar.NameOfFileForSavingSlipperinessEstimate=[CtrlVar.Experiment,'_C-Estimate.mat'];
CtrlVar.NameOfFileForSavingAGlenEstimate=[CtrlVar.Experiment,'_AGlen-Estimate.mat'];

CtrlVar.Inverse.Measurements='-uv-';

% It is usually better to invert for log(A) and log(C) rather than A and C.
% The default is to invert for log(A) and log(C) simultaneously.
%CtrlVar.Inverse.InvertFor= 'logAGlenlogC'; % {'C','logC','AGlen','logAGlen','logAGlenlogC'}
CtrlVar.Inverse.InvertFor='logAGlenlogC';

% The gradient of the objective function is calculated using the adjoint method.
% When inverting for C only, one can also use a gradient based on a `FixPoint'
% iteration, which is often a very good initial approach. 
CtrlVar.Inverse.DataMisfit.GradientCalculation='Adjoint'; % {'Adjoint','FixPoint'}
%CtrlVar.Inverse.DataMisfit.GradientCalculation='Adjoint';

% The gradient of the objective function can be premultiplied with the inverse
% of the mass matrix. This creates a `mesh independent' gradient. This has both
% advantages and disadvantages. 
CtrlVar.Inverse.AdjointGradientPreMultiplier='I'; % {'I','M'}

% 
% CtrlVar.Inverse.TestAdjoint.isTrue=0; % If true then perform a brute force calculation 
%                                       % of the directinal derivative of the objective function.  
% %CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType='central' ; % {'central','forward'}
% CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize=1e-8 ;
% CtrlVar.Inverse.TestAdjoint.iRange=[] ;

% Regularisation can be applied on A and C or log(A) and log(C). Also possible
% to use a covariance matrix for A and C. 
%
% Select Bayesian motivated regularisation by setting 
% CtrlVar.Inverse.Regularize.Field='cov' and Tikhonov regularisation
% by setting CtrlVar.Inverse.Regularize.Field to either 'C','logC','AGlen','logAGlen','logAGlenlogC'
%
% Default is Tikhonov regularisation on log(A) and log(C)
CtrlVar.Inverse.Regularize.Field= 'logAGlenlogC'; % {'cov','C','logC','AGlen','logAGlen','logAGlenlogC'}

% [ -- Parameters specific to Tikhonov regularisation
% See the above definition of R in the case of Tikhonov regularisation.
% The values of these parameters can be expected to be highly problem dependent.
% By default regularisation is switched on, but can the switched off by setting
% the gs and the ga parameters to zero.
CtrlVar.Inverse.Regularize.C.gs=0; 
CtrlVar.Inverse.Regularize.C.ga=0;
CtrlVar.Inverse.Regularize.logC.gs=1e4;
CtrlVar.Inverse.Regularize.logC.ga=1;

CtrlVar.Inverse.Regularize.AGlen.gs=0;
CtrlVar.Inverse.Regularize.AGlen.ga=0;
CtrlVar.Inverse.Regularize.logAGlen.gs=1e4;
CtrlVar.Inverse.Regularize.logAGlen.ga=1;

% Some less often used parameters related to inversion 
CtrlVar.Inverse.InfoLevel=1;  % Set to 1 to get some basic information, >=2 for additional info on backtrackgin,
                              % >=100 for further info and plots
% In an inversion it it generally better to set other infolevels to a low value. So
% consider setting:
%CtrlVar.InfoLevelNonLinIt=1000; CtrlVar.InfoLevel=1000;                                           
                                                        

%%
CtrlVar.ThicknessConstraints=0;
CtrlVar.ResetThicknessToMinThickness=1;  % change this later on
CtrlVar.ThickMin=1;

end

function UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar)

%% target geometry interpolants
%% ICE SHELF
UserVar.ISGeometry = RunTable{ind,"ISthick"}; % this is the geometry BEFORE the perturbation
if ismember(UserVar.ISGeometry,[2000,2009,2014,2018])
    UserVar.ISGeometryInterpolants = [pwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.ISGeometry),'_EXTRUDED.mat'];
elseif floor(UserVar.ISGeometry/1e3)==1
    NameOfFiletoRead = "../ANT_Inverse/ANT_Inverse_"+string(UserVar.ISGeometry)+...
        "/ANT_Inverse_"+string(UserVar.ISGeometry)+"-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleIS"})+".mat";
    load(NameOfFiletoRead,"F","MUA");
    FB = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.B,"linear");
    Fb = FB; Fb.Values = F.b;
    Fs = FB; Fs.Values = F.s;
    Frho = FB; Frho.Values = F.rho;
    UserVar.ISGeometryInterpolants = "GeometryInterpolantsIS.mat";
    save("./"+UserVar.Experiment+"/"+UserVar.ISGeometryInterpolants,"FB","Fb","Fs","Frho");
else
    error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise ISthick flag in RunTable.']);
end

%% GROUNDED ICE
UserVar.GIGeometry = RunTable{ind,"GIthick"}; % this is the geometry BEFORE the perturbation
if ismember(UserVar.GIGeometry,[2000,2009,2014,2018])
    UserVar.GIGeometryInterpolants = [pwd,'../../ANT_Data/ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-',num2str(UserVar.GIGeometry),'_EXTRUDED.mat'];
elseif floor(UserVar.GIGeometry/1e3)==1
    NameOfFiletoRead = "../ANT_Inverse/ANT_Inverse_"+string(UserVar.GIGeometry)+...
        "/ANT_Inverse_"+string(UserVar.GIGeometry)+"-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleGI"})+".mat";
    load(NameOfFiletoRead,"F","MUA");
    FB = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.B,"linear");
    Fb = FB; Fb.Values = F.b;
    Fs = FB; Fs.Values = F.s;
    UserVar.GIGeometryInterpolants = "GeometryInterpolantsGI.mat";
    save("./"+UserVar.Experiment+"/"+UserVar.GIGeometryInterpolants,"FB","Fb","Fs");
else
    error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise GIthick flag in RunTable.']);
end

%% CALVING
UserVar.Calv = RunTable{ind,"Calv"};


%% density interpolant: same as ice shelf run
UserVar.DensityInterpolant = UserVar.ISGeometryInterpolants;

%% target mesh: take mesh from ice shelf run as base mesh 
if floor(UserVar.ISGeometry/1e3)==1
    UserVar.BaseMesh.BCs = "../ANT_Inverse/ANT_Inverse_"+string(UserVar.ISGeometry)+...
        "/ANT_Inverse_"+string(UserVar.ISGeometry)+"_MeshBoundaryCoordinates.mat";
    NameOfFiletoRead = "../ANT_Inverse/ANT_Inverse_"+string(UserVar.ISGeometry)+...
        "/ANT_Inverse_"+string(UserVar.ISGeometry)+"-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleIS"})+".mat";
    load(NameOfFiletoRead,"CtrlVarInRestartFile","MUA"); CtrlVar = CtrlVarInRestartFile;
    UserVar.BaseMesh.Mesh = "./"+UserVar.Experiment+"/"+"BaseMesh.mat";
    save(UserVar.BaseMesh.Mesh,"CtrlVar","MUA");
else
    error(['ExpID ',RunTable{ind,"ExpID"},': Do not recognise ISthick flag in RunTable.']);
end

UserVar = ANT_ApplyMeshModifications(UserVar);
   
%% sliding law
UserVar.InverseC = RunTable{ind,"InverseC"};
UserVar.InverseCFill = RunTable{ind,"InverseCFill"};

% copy relevant files   
NameOfFiletoRead="../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseC)+"/"+"ANT_Inverse_"+string(UserVar.InverseC)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleC"})+".mat";
load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
C = F.C; m = F.m;
UserVar.NameOfFileForReadingSlipperinessEstimate = "ANT_Inverse_"+string(UserVar.InverseC)+"_C-Estimate.mat";
save(UserVar.Experiment+"/"+UserVar.NameOfFileForReadingSlipperinessEstimate,"MUA","C","m");
UserVar.SlidingLaw = CtrlVarInRestartFile.SlidingLaw;

NameOfFiletoRead="../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseCFill)+"/"+"ANT_Inverse_"+string(UserVar.InverseCFill)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleC"})+".mat";
load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
C = F.C; m = F.m;
UserVar.NameOfFileForReadingSlipperinessEstimateFill = "ANT_Inverse_"+string(UserVar.InverseCFill)+"_C-Estimate.mat";
save(UserVar.Experiment+"/"+UserVar.NameOfFileForReadingSlipperinessEstimateFill,"MUA","C","m");

%% Glen's exponent
UserVar.InverseA = RunTable{ind,"InverseA"};
UserVar.InverseAFill = RunTable{ind,"InverseAFill"}; 

% copy relevant files   
NameOfFiletoRead="../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseA)+"/"+"ANT_Inverse_"+string(UserVar.InverseA)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleA"})+".mat";
load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
AGlen = F.AGlen; n = F.n;
UserVar.NameOfFileForReadingAGlenEstimate = "ANT_Inverse_"+string(UserVar.InverseA)+"_AGlen-Estimate.mat";
save(UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimate,"MUA","AGlen","n");

NameOfFiletoRead="../ANT_Inverse/ANT_Inverse_"+string(UserVar.InverseAFill)+"/"+"ANT_Inverse_"+string(UserVar.InverseAFill)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleA"})+".mat";
load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
AGlen = F.AGlen; n = F.n;
UserVar.NameOfFileForReadingAGlenEstimateFill = "ANT_Inverse_"+string(UserVar.InverseAFill)+"_AGlen-Estimate.mat";
save(UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimateFill,"MUA","AGlen","n");

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

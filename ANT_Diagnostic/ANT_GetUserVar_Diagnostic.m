function UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar)

%% target geometry interpolants
%% ICE SHELF
UserVar.ISGeometry = RunTable{ind,"ISthick"}; % this is the geometry BEFORE the perturbation
NameOfFiletoRead = "../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(UserVar.ISGeometry)+...
        "/"+UserVar.Domain+"_Inverse_"+string(UserVar.ISGeometry)+"-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleIS"})+".mat";
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA");
    FB = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.B,"linear");
    Fb = FB; Fb.Values = F.b;
    Fs = FB; Fs.Values = F.s;
    Frho = FB; Frho.Values = F.rho;
    UserVar.ISGeometryInterpolants = "GeometryInterpolantsIS.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.ISGeometryInterpolants,"FB","Fb","Fs","Frho");
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to read ISthick geometry: "+NameOfFiletoRead);
end

%% GROUNDED ICE
UserVar.GIGeometry = RunTable{ind,"GIthick"}; % this is the geometry BEFORE the perturbation
NameOfFiletoRead = "../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(UserVar.GIGeometry)+...
        "/"+UserVar.Domain+"_Inverse_"+string(UserVar.GIGeometry)+"-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleGI"})+".mat";
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA");
    FB = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.B,"linear");
    Fb = FB; Fb.Values = F.b;
    Fs = FB; Fs.Values = F.s;
    Frho = FB; Frho.Values = F.rho;
    UserVar.GIGeometryInterpolants = "GeometryInterpolantsGI.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.GIGeometryInterpolants,"FB","Fb","Fs","Frho");
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to read GIthick geometry: "+NameOfFiletoRead);
end

%% CALVING
UserVar.Calv = RunTable{ind,"Calv"};
NameOfFiletoRead = "../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(UserVar.Calv)+...
        "/"+UserVar.Domain+"_Inverse_"+string(UserVar.Calv)+"-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleGI"})+".mat";
if exist(NameOfFiletoRead,"file")
    % new boundary
    load(NameOfFiletoRead,"MUA");
    MeshBoundaryCoordinates = [MUA.Boundary.x(:) MUA.Boundary.y(:)];
    MeshBoundaryCoordinatesFile = UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.Experiment+"_MeshBoundaryCoordinates";
    save(MeshBoundaryCoordinatesFile,"MeshBoundaryCoordinates");
    UserVar.MeshBoundaryCoordinatesFile = UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.Experiment+"_MeshBoundaryCoordinates";
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to read calving front: "+NameOfFiletoRead);
end

% new mesh
InitialMeshFileName = UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.Experiment+"_InitialMesh";
save(InitialMeshFileName,"MUA");
UserVar.InitialMeshFileName = "./"+UserVar.Experiment+"_InitialMesh";

%%basemesh
UserVar = ANT_DefineBaseMesh(UserVar,RunTable{ind,"BaseMesh"}{:});
    
%% density interpolant: same as ice shelf run
UserVar.DensityInterpolant = UserVar.ISGeometryInterpolants;   

%% sliding law
UserVar.InverseC = RunTable{ind,"InverseC"};

% copy relevant files   
NameOfFiletoRead="../ANT_Inverse/cases/"+UserVar.Domain+"_Inverse_"+string(UserVar.InverseC)+"/"+UserVar.Domain+"_Inverse_"+string(UserVar.InverseC)+...
    "-RestartFile_InverseCycle"+string(RunTable{ind,"InverseCycleC"})+".mat";
if exist(NameOfFiletoRead,"file")
    load(NameOfFiletoRead,"F","MUA","CtrlVarInRestartFile");
    C = F.C; m = F.m;
    UserVar.NameOfFileForReadingSlipperinessEstimate = "ANT_Inverse_"+string(UserVar.InverseC)+"_C-Estimate.mat";
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
    UserVar.NameOfFileForReadingAGlenEstimate = "ANT_Inverse_"+string(UserVar.InverseA)+"_AGlen-Estimate.mat";
    save(UserVar.casefolder+"/"+UserVar.Experiment+"/"+UserVar.NameOfFileForReadingAGlenEstimate,"MUA","AGlen","n");
else
    error("ExpID "+RunTable{ind,"ExpID"}+": Could not find file to copy AGlen field: "+NameOfFiletoRead);
end   

%% outputs
UserVar.UaOutputDirectory = './ResultsFiles';

%% write log file
fprintf(UserVar.fid_experimentlog,"> ANT_GetUserVar_Diagnostic: ExpID %s: Starting diagnostic experiment on %s.\n",string(UserVar.ExpID),UserVar.hostname);



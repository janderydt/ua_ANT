function [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

CFile = UserVar.NameOfFileForReadingSlipperinessEstimate;

q=1 ;      % only needed for Budd sliding law
muk=0.5 ; 

tmp = load(CFile,'MUA','C','m');
m = tmp.m(1);

% map C onto basemesh, with zeros where there is no data
base = load("../"+UserVar.BaseMesh.Mesh);
CtrlVar.MapOldToNew.method = "ShapeAndScattered";
C_outside = 0;
[~,C_base] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,base.MUA,C_outside,tmp.C);
C_base(UserVar.BaseMesh.DeactivatedNodes) = 0;

[~,C] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],base.MUA,MUA,C_outside,C_base);
Izero = find(C==0);

% where no value for C is available, take it from UserVar.NameOfFileForReadingSlipperinessEstimateFill
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
if ~isempty(Izero)
    FillC = load(UserVar.NameOfFileForReadingSlipperinessEstimateFill);
    FCFill = scatteredInterpolant(FillC.MUA.coordinates(:,1),FillC.MUA.coordinates(:,2),FillC.C,'nearest');
    C(Izero) = FCFill(x(Izero),y(Izero));
    fprintf("Using C from file %s and filling holes with %s.\n",CFile,UserVar.NameOfFileForReadingSlipperinessEstimateFill);
else
    fprintf("Using C from file %s.\n",CFile);
end

Izero = find(C==0);

if ~isempty(Izero)
    error('Still zeros in C');
end


save("C_interpolated.mat","CtrlVar","MUA","C");

function [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

CFile = UserVar.NameOfFileForReadingSlipperinessEstimate;

q=1 ;      % only needed for Budd sliding law
muk=0.5 ; 

tmp = load(CFile,'MUA','C','m');
m = tmp.m(1);

% check if we need to interpolate C
if tmp.MUA.Nnodes ~= MUA.Nnodes % something is different between the meshes

    fprintf("Interpolating C from old to new mesh using MapNodalVariablesFromMesh1ToMesh2. " + ...
    "Nnodes old mesh = %s, Nnodes new mesh = %s\n",string(tmp.MUA.Nnodes),string(MUA.Nnodes));
    CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
    tau = 80; %kPa
    ub = 100; %m/yr
    C_outside = ub/tau^m;

    [~,AGlen] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,C_outside,tmp.C);

else
    % meshes have the same number of nodes. 
    C = tmp.C;

end

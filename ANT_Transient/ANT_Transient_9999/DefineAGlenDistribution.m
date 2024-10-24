function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

AGlenFile = UserVar.NameOfFileForReadingAGlenEstimate;

fprintf("Using AGlen from file %s\n",AGlenFile);

tmp = load(AGlenFile,'MUA','AGlen','n'); % load results from inversion
n = tmp.n(1);

% check if we need to interpolate AGlen 
if tmp.MUA.Nnodes ~= MUA.Nnodes % something is different between the meshes

    fprintf("Interpolating AGlen from old to new mesh using MapNodalVariablesFromMesh1ToMesh2. " + ...
    "Nnodes old mesh = %s, Nnodes new mesh = %s\n",string(tmp.MUA.Nnodes),string(MUA.Nnodes));
    CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
    tau = 80; 
    eps = 0.0026;
    AGlen_outside = eps/tau^n;

    [~,AGlen] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,AGlen_outside,tmp.AGlen);

else
    % meshes have the same number of nodes. 
    AGlen = tmp.AGlen;

end

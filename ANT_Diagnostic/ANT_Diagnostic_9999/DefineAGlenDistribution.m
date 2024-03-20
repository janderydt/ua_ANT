function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

AGlenFile = UserVar.NameOfFileForReadingAGlenEstimate;

tmp = load(AGlenFile,'MUA','AGlen','n');
n = tmp.n(1);

% map AGlen onto basemesh, with zeros where there is no data
CtrlVar.MapOldToNew.method = "ShapeAndScattered";
AGlen_outside = 0;
[~,AGlen] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,AGlen_outside,tmp.AGlen);

Izero = find(AGlen==0);

% where no value for AGlen is available, take it from UserVar.NameOfFileForReadingAGlenEstimateFill
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
if ~isempty(Izero)
    FillA = load(UserVar.NameOfFileForReadingAGlenEstimateFill);
    FAGlenFill = scatteredInterpolant(FillA.MUA.coordinates(:,1),FillA.MUA.coordinates(:,2),FillA.AGlen,'nearest');
    AGlen(Izero) = FAGlenFill(x(Izero),y(Izero));
    fprintf("Using AGlen from file %s and filling holes with %s.\n",AGlenFile,UserVar.NameOfFileForReadingAGlenEstimateFill);
else
    fprintf("Using AGlen from file %s.\n",AGlenFile);
end

Izero = find(AGlen==0);

if ~isempty(Izero)
    error('Still zeros in AGlen');
end

save("AGlen_interpolated.mat","CtrlVar","MUA","AGlen");

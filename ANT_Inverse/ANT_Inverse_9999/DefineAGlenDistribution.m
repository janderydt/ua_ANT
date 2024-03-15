function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FA;

n = UserVar.n;

AGlenFile = UserVar.NameOfFileForReadingAGlenEstimate;

if isempty(FA) & exist(AGlenFile,"file")

    tmp = load(AGlenFile,'MUA','AGlen');
    CtrlVar.MapOldToNew.method = "ShapeAndScattered";
    AGlen_outside = AGlenVersusTemp(-15);
    [~,AGlen] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,AGlen_outside,tmp.AGlen);

    save("AGlen_interpolated.mat","CtrlVar","MUA","AGlen");

    fprintf("Using AGlen from file %s.\n",AGlenFile);

    FA = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),AGlen,'linear');

elseif ~exist(AGlenFile,"file")

    AGlen=s*0+AGlenVersusTemp(-15);
    fprintf("Used %s as constant intital value for AGlen.\n",string(AGlen(1)));

elseif ~isempty(FA)

    AGlen = FA(MUA.coordinates(:,1),MUA.coordinates(:,2));

end

end

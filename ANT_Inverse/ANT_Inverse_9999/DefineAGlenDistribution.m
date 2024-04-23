function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FA;

n = UserVar.n;

AGlenFile = UserVar.NameOfFileForReadingAGlenEstimate;

if isempty(FA) & exist(AGlenFile,"file")

    tmp = load(AGlenFile,'MUA','AGlen','n');
    CtrlVar.MapOldToNew.method = "ShapeAndScattered";
    AGlen_outside = AGlenVersusTemp(-15);
    [~,AGlen] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,AGlen_outside,tmp.AGlen);

    save("AGlen_interpolated.mat","CtrlVar","MUA","AGlen");

    fprintf("Using AGlen from file %s.\n",AGlenFile);

    % check if AGlen needs to be rescaled to reflect new value of n
    n1 = tmp.n(1);
    n2 = n; % new

    if n1 ~= n2
        velfile = "../"+erase(AGlenFile,"_AGlen-Estimate.mat")+"/"+strrep(AGlenFile,"_AGlen-Estimate.mat","-RestartFile.mat");
        tmp = load(velfile,'MUA','F');
        [~,~,~,exxInt,eyyInt,exyInt] = calcStrainRatesEtaInt(CtrlVar,tmp.MUA,tmp.F.ub,tmp.F.vb,tmp.F.AGlen,tmp.F.n);
        [exx,eyy,exy]=ProjectFintOntoNodes(MUA,exxInt,eyyInt,exyInt);
        e_old = real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));
        [~,e_new] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,1e-4,e_old);
        AGlen = max((max(AGlen+CtrlVar.AGlenAdjointZero,CtrlVar.AGlenmin)).^(n2/n1).*(e_new+eps).^((n1-n2)/(n1))-CtrlVar.AGlenAdjointZero,CtrlVar.AGlenAdjointZero);
        fprintf("Rescaling AGlen values from n=%s to n=%s.\n",string(n1),string(n2));
    end

    FA = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),AGlen,'linear');

elseif ~exist(AGlenFile,"file")

    AGlen=s*0+UserVar.Inverse.priorAGlen;
    fprintf("Used %s as constant intital value for AGlen.\n",string(AGlen(1)));

elseif ~isempty(FA)

    AGlen = FA(MUA.coordinates(:,1),MUA.coordinates(:,2));

end

end

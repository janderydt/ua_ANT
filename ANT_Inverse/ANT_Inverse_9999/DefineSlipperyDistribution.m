function [UserVar,C,m,q,muk,V0]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FC;

CFile = UserVar.NameOfFileForReadingSlipperinessEstimate;

m = UserVar.SlidingCoefficient;
q = 1 ;      % only needed for Budd sliding law
V0 = [];

if ~isfield(UserVar,'muk')
    muk=0.5;
else
    muk = UserVar.muk ;   % only needed for mixed Weertman/Coulomb sliding law
end

if isempty(FC) & exist(CFile,"file")

    tmp = load(CFile,'MUA','C','m');
    CtrlVar.MapOldToNew.method = "ShapeAndScattered";
    ub=100; tau=80 ; % units meters, year , kPa
    C_outside = ub/tau^m;
    [~,C] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,C_outside,tmp.C);

    fprintf("Using C from file %s.\n",CFile);

    % check if C needs to be rescaled to reflect new value of m
    m1 = tmp.m(1);
    m2 = m; % new

    if m1 ~= m2
        velfile = "../"+erase(CFile,"_C-Estimate.mat")+"/"+strrep(CFile,"_C-Estimate.mat","-RestartFile.mat");
        tmp = load(velfile,'MUA','F');
        [~,v1] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,0,hypot(tmp.F.ub,tmp.F.vb));
        C = max((max(C+CtrlVar.Czero,CtrlVar.Cmin)).^(m2/m1).*(v1.^2+CtrlVar.SpeedZero^2).^((m1-m2)/(2*m1))-CtrlVar.Czero,CtrlVar.Cmin);
        fprintf("Rescaling C values from m=%s to m=%s.\n",string(m1),string(m2));
    end

    save("C_interpolated.mat","CtrlVar","MUA","C");

    FC = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),C,'linear');

elseif ~exist(CFile,"file")

    C=s*0+UserVar.Inverse.priorC;
    fprintf("Used %s as constant intital value for C.\n",string(C(1)));

elseif ~isempty(FC)

    C = FC(MUA.coordinates(:,1),MUA.coordinates(:,2));

end
    
   
end
    

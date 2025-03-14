function [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FC;

CFile = UserVar.NameOfFileForReadingSlipperinessEstimate;

m = UserVar.SlidingCoefficient;
q = 1 ;      % only needed for Budd sliding law
%V0 = [];

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

    velfile = UserVar.casefolder+"/"+erase(CFile,"_C-Estimate.mat")+"/"+strrep(CFile,"_C-Estimate.mat","-RestartFile_InverseCycle1.mat");
    if m1 ~= m2       
        fprintf("Rescaling C values from m=%s to m=%s using measured velocities from %s.\n",string(m1),string(m2),velfile);
        tmp = load(velfile,'MUA','Meas');
        [~,v1] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,0,hypot(tmp.Meas.us,tmp.Meas.vs));
        v1(isnan(v1))=1; v1(v1==0)=1;
        C = max((max(C+CtrlVar.Czero,CtrlVar.Cmin)).^(m2/m1).*(v1.^2+CtrlVar.SpeedZero^2).^((m1-m2)/(2*m1))-CtrlVar.Czero,CtrlVar.Cmin);
    end

    % this script is typically used after a fixpoint inversion, and the
    % resulting C field is used as an input to the adjoint inversion.
    % If so, the C field will have short wavelength - high amplitude 
    % variations, which we might not want to feed into the algorithm as a
    % starting point for the adjoint. Here we apply a Helmholtz smoothing
    % to remove the short wavelength variations
    tmp = load(velfile,"UserVarInRestartFile");
    if isfield(tmp.UserVarInRestartFile.Inverse,"GradientCalculation")
        if tmp.UserVarInRestartFile.Inverse.GradientCalculation == "FixPoint"
            L=10e3 ;  % Smoothing length scale 
            [~,log10C_smooth]=HelmholtzEquation([],CtrlVar,MUA,1,L^2,log10(C),0); 
            C = 10.^log10C_smooth;
            C(C < CtrlVar.Cmin) = CtrlVar.Cmin; 
            C(C > CtrlVar.Cmax) = CtrlVar.Cmax;
        end   
    end
    
    save("C_interpolated.mat","CtrlVar","MUA","C");

    FC = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),C,'linear');

elseif ~exist(CFile,"file")

    if UserVar.Inverse.startC == -9999
        fprintf("Using measured velocities to calculate initial guess for C.\n");
        tau = 80;
        tmp = load(UserVar.VelocityInterpolants);
        x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
        uMeas = tmp.Fus(x,y);
        vMeas = tmp.Fvs(x,y);
        clear tmp;
        ub = hypot(uMeas,vMeas);
        ub(isnan(ub))=1;
        C = ub/tau^m; 
        C = max(C,CtrlVar.Cmin);
    else
        C=s*0+UserVar.Inverse.priorC;
        fprintf("Used %s as constant intital value for C.\n",string(C(1)));
    end

elseif ~isempty(FC)

    C = FC(MUA.coordinates(:,1),MUA.coordinates(:,2));

end
    
   
end
    

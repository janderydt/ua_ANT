function [UserVar,C,m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FC

if exist(CtrlVar.NameOfFileForReadingSlipperinessEstimate)~=2
    
    m = 3;
    ub=750; tau=80 ; % units meters, year , kPa
    C=s*0+ub/tau^m;
    
else
    if isempty(FC)
        load(CtrlVar.NameOfFileForReadingSlipperinessEstimate,'xC','yC','C');
        FC = scatteredInterpolant(xC,yC,C,'linear');
        fprintf('\n Read slipperiness from file \n');
    end

    load(CtrlVar.NameOfFileForReadingSlipperinessEstimate,'m');
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    C = FC(x,y);
    m = m(1);

end

end

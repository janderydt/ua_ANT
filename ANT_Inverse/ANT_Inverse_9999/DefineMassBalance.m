function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent Fsmb_RACMO_climatology;

x=MUA.coordinates(:,1);
y=MUA.coordinates(:,2);
    
if UserVar.SpinupCycle

   
    %% Surface mass balance: RACMO 2000-2018 climatology
    
    if isempty(Fsmb_RACMO_climatology)
    
        load("../../ANT_Data/ANT_Interpolants/ScatteredInterpolants_SMB.mat","Fsmb_RACMO");
        
        % RACMO climatology between 2000 and 2018
        Istart = find(contains(string(Fsmb_RACMO.years),"2000"));
        Iend = find(contains(string(Fsmb_RACMO.years),"2018"));
        dn = Iend-Istart+1;
        RACMO_smb = 0;
        for ii=1:dn
            yrstr = "yr"+string(1999+ii);
            RACMO_smb = RACMO_smb + Fsmb_RACMO.(yrstr).Values;
        end
        Fsmb_RACMO_climatology = Fsmb_RACMO.climatology; Fsmb_RACMO_climatology.Values = RACMO_smb/dn;
    end
    
    as = Fsmb_RACMO_climatology(x,y);
    
    %% Basal mass balance: "balanced melt"
    if F.time>0
        ab=CalcIceShelfMeltRates(CtrlVar,MUA,F.ub,F.vb,F.s,F.b,F.S,F.B,F.rho,F.rhow,0*x,as,0*x);
    else
        ab = 0*x;
    end
    
    % apply melt on nodes downstream of GL
    GF =IceSheetIceShelves(CtrlVar,MUA,F.GF);
    [LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF);
    
    I = find(LakeNodes & GF.node==1);
    ab(I) = 0;

else

    as = 0*x;
    ab = 0*x;

end


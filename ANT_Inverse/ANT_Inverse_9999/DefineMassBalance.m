function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,F)

persistent Fsmb_RACMO_climatology Fdhdt;

x=MUA.coordinates(:,1);
y=MUA.coordinates(:,2);

if UserVar.SpinupCycle
    %% Surface mass balance: RACMO 2000-2020 climatology
    
    if isempty(Fsmb_RACMO_climatology)
    
        load(UserVar.datafolder+"ANT_Interpolants/ScatteredInterpolants_SMB.mat","Fsmb_RACMO");
        
        % RACMO climatology between 2000 and 2020
        Istart = find(contains(string(Fsmb_RACMO.years),"2000"));
        Iend = find(contains(string(Fsmb_RACMO.years),"2020"));
        dn = Iend-Istart+1;
        RACMO_smb = 0;
        for ii=1:dn
            yrstr = "yr"+string(1999+ii);
            RACMO_smb = RACMO_smb + Fsmb_RACMO.(yrstr).Values;
        end
        Fsmb_RACMO_climatology = Fsmb_RACMO.climatology; 
        Fsmb_RACMO_climatology.Values = RACMO_smb/dn; % in m w.e. / yr
    end
    
    as = Fsmb_RACMO_climatology(x,y)*1000./max(F.rho,500); % transfer to m ice equivalent / yr, putting some random cutoff on the density to avoid unrealistically high accumulation rates
    
    %% Basal mass balance: "balanced melt"
    ab = 0*x;
    if F.time>0
        F.as = as;
        F.ab = ab;
        BCs=BoundaryConditions;
        [~,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F);
        [~,dhdt]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F,BCs);
        ab = -dhdt;
        %ab=CalcIceShelfMeltRates(CtrlVar,MUA,F.ub,F.vb,F.s,F.b,F.S,F.B,F.rho,F.rhow,0*x,as,0*x);
    end
    
    % apply melt on nodes downstream of GL
    GF =IceSheetIceShelves(CtrlVar,MUA,F.GF);   
    %I = find(LakeNodes | GF.node>0.5); ab(I)=0;
    NodesDownstreamOfGroundingLines="Relaxed"; % strict
    [LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF,[],NodesDownstreamOfGroundingLines);
    ab(LakeNodes) = 0;
    ab(~OceanNodes) = 0;

    %% Make adjustments for second spinup cycle: adjust surface mass
    %% balance as -> as - dhdt and run to steady state with fixed ice shelf thickness
    if UserVar.SpinupCycle && UserVar.Spinup.Cycle > 1

        ab = 0*x; % ab can be zero because we keep ice shelf thickness fixed as boundary condition

        dhdt_filename = UserVar.datafolder+"/ANT_Interpolants/dhdt_"+string(UserVar.Geometry)+"_"+string(UserVar.Geometry+1)+".mat";

        if exist(dhdt_filename,"file")==2

            if isempty(Fdhdt)
                load(dhdt_filename);
                % only retain floating nodes from Paolo and grounded nodes from
                % Otosaka
                FGFmask = scatteredInterpolant(x,y,GF.node);
                GFo = FGFmask(Xo,Yo);
                GFp = FGFmask(Xp,Yp);
                Io = find(GFo > 0.5 & ~isnan(dhdto));
                Ip = find(GFp < 0.5 & ~isnan(dhdtp));
                xnew = [Xo(Io); Xp(Ip)]; ynew = [Yo(Io); Yp(Ip)];
                dhdt_obs = [dhdto(Io); dhdtp(Ip)];
                Fdhdt = scatteredInterpolant(xnew(:),ynew(:),dhdt_obs(:),'natural');
            end

        else

            error(dhdt_filename+" does not exist.");

        end

        as = as - Fdhdt(x,y);
        fprintf("Removed observed dhdt from as.\n");
            
    end


elseif UserVar.InverseCycle

    ab = 0*x;
    as = 0*x;

end

end


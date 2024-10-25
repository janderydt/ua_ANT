function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% © Qing Qin    October 2024 %%%
%%% qing.qin@northumbria.ac.uk %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch UserVar.BasalMelt
 
    case "PICO"
    
         %addpath(genpath('/home/qingqin/Documents/PICO_Ua-master/'))
         
         load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/BasinsInterpolant.mat");
    
         % Check if the file loaded successfully by checking the structure fields
         if isempty(Fbasins)
            error('Failed to load BasinsInterpolant.mat. The file is empty.');
         else
            disp('File loaded successfully.');
         end
         MUA=UpdateMUA(CtrlVar,MUA);
         PICO_opts = struct;
         PICO_opts.algorithm = 'watershed';%'polygon','oneshelf';
         PICO_opts.nmax = 10;% max number of boxes 5/10
         PICO_opts.minArea = 2e9;
         PICO_opts.minNumShelf = 20;
         PICO_opts.SmallShelfMelt = 0;
         PICO_opts.PICOres = 10000; % resolution in km (for watershed algorithm only)
         PICO_opts.BasinsInterpolant = Fbasins;
         PICO_opts.FloatingCriteria = 'GLthreshold'; %'GLthreshold' or 'StrictDownstream'
         PICO_opts.persistentBC = 0;
         PICO_opts.InfoLevel = 1; % 0,1,10,100
         % ---------------read PICO parameters from Runtable-------------------
         PICO_opts.C1 = UserVar.PICOC1;
         PICO_opts.gamTstar = UserVar.PICOgam;
         % ----------------oceanic forcing data for PICO-----------------------
         switch UserVar.OceForcing
             case "ismip6"
               ISMIP6_obs = load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/obs_oceanic_forcing_ISMIP6.mat");
               PICO_opts.Sbasins = ISMIP6_obs.salt_ID_30;
               PICO_opts.Tbasins = ISMIP6_obs.temp_ID_30;
         end
         % --------------------------------------------------------------------
         PICO_opts.MeshBoundaryCoordinates = CtrlVar.MeshBoundaryCoordinates;
         [Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(1,CtrlVar,MUA,GF,h,median(rho),rhow,PICO_opts);
         ab=Mk;
         NodesDownstreamOfGroundingLines="Relaxed"; % strict
         [LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF,[],NodesDownstreamOfGroundingLines);
         ab(LakeNodes)=0;
         ab(~OceanNodes) = 0;

    case "LQ"
         % ------------read Local Quadratic parameter from Runtable----------------
         gamma0 = UserVar.LQgam;
         % ------------------------------------------------------------------------
         ncfile = UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/coeff_gamma0_DeltaT_quadratic_local_median.nc";
         deltaT_basin = ncread(ncfile,'deltaT_basin');

         x = double(ncread(ncfile,'x')); % units kms, projection: 
         y = double(ncread(ncfile,'y'));

         [X,Y] = ndgrid(x,y);
         deltaT_basin1 = scatteredInterpolant(X(:),Y(:),deltaT_basin(:),'linear');
         deltaT_basin2 = deltaT_basin1(MUA.coordinates(:,1),MUA.coordinates(:,2));% ?

         rhoi_SI=918.0; % ice density (kg/m^3)
         rhosw_SI=1028.0; % sea water density
         rhofw_SI=1000.0; % freshwater density
         Lf_SI=3.34e5; % fusion latent heat of Ice (J/kg)
         cpw_SI=3974.0; % specific heat of sea water (J/kg/K)
         % ------------- oceanic forcing data for LQ-----------------------
         switch UserVar.OceForcing
             case "ismip6"
                  ISMIP6_obs_tf= load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/TF_climatology_ISMIP6.mat");
                  thermal_forcing = ISMIP6_obs_tf.TF([MUA.coordinates b]);
         end
         % ------------------------------------------------------------------------
         ab=-gamma0.*(rhosw_SI.*cpw_SI./rhoi_SI./Lf_SI).^2.*(max(thermal_forcing+deltaT_basin2,0.0)).^2; % here is mass balance, which needs the negative sign

         NodesDownstreamOfGroundingLines="Relaxed"; % strict
         [LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF,[],NodesDownstreamOfGroundingLines);
         ab(LakeNodes)=0;
         ab(~OceanNodes) = 0;
         clear tf_anomaly; 
         
     otherwise
     
         error("Unknown UserVar.BasalMelt: "+UserVar.BasalMelt);
end


switch UserVar.smb

    case "MAR"
    
	    ncfile = UserVar.datafolder+"/ANT_InputsForTransientSimulations/SMB/year-MAR_NorESM-1980-2100_zen.nc2";
	    fprintf(" Reading MAR SMB data from file %s \n",ncfile);

	    x = double(ncread(ncfile,'X')); % units kms, projection: 
	    y = double(ncread(ncfile,'Y'));
	    % t = ncread(ncfile,'TIME'); % days since 1901-01-15 00:00:00
	    smb = ncread(ncfile,'SMB'); %annual averages, units kg m-2 yr-1, missing value -1.e+34

	    % epochnum = datenum('1901-01-15 00:00:00');
	    % t = double(epochnum+t);
	    smb = double(smb)./1000; % units from kg m^-2 year^-1 to m/year
	    [X,Y] = ndgrid(x,y);
	    year_anom=floor(time+1);
	    yr=20+year_anom;
	    smb_yrly = smb(:,:,1,yr);
	    I = find(~isnan(smb_yrly(:)));
	    Fsmb_MAR = scatteredInterpolant(X(I).*1000,Y(I).*1000,smb_yrly(I),'linear');
	    as = Fsmb_MAR(MUA.coordinates(:,1),MUA.coordinates(:,2));

	    fprintf("Done. \n");
	  
    otherwise
    
        error("Unknown UserVar.smb: "+UserVar.smb);
   
end



end


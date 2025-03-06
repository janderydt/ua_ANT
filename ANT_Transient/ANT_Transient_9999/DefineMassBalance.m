function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% © Qing Qin    October 2024 %%%
%%% qing.qin@northumbria.ac.uk %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent Fbasins FdeltaT_basin ISMIP6_obs ISMIP6_obs_tf Fsmb_MAR

switch UserVar.BasalMelt
 
    case "PICO"
    
         if UserVar.hostname ~= "ARCHER2"
             addpath(genpath('/mnt/md0/Ua/UaPICO_master/'));
         end

         load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/BasinsInterpolant.mat");
    
         % Check if the file loaded successfully by checking the structure fields
         if isempty(Fbasins)
            load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/BasinsInterpolant.mat");
         else
            disp('File loaded successfully.');
         end
         MUA=UpdateMUA(CtrlVar,MUA); % make sure MUA.TR is non-empty
         PICO_opts = struct;
         PICO_opts.algorithm = 'watershed';%'polygon','oneshelf';
         PICO_opts.nmax = 10;% max number of boxes 5/10
         PICO_opts.ContinentArea = 5e10;
         PICO_opts.minArea = 1e9;
         PICO_opts.minNumShelf = 20;
         PICO_opts.SmallShelfMelt = 0;
         %PICO_opts.PICOres = 10000; % resolution in km (for watershed algorithm only)
         PICO_opts.BasinsInterpolant = Fbasins;
         PICO_opts.FloatingCriteria = 'GLthreshold'; %'GLthreshold' or 'StrictDownstream'
         PICO_opts.persistentBC = 0;
         PICO_opts.InfoLevel = 1; % 0,1,10,100
         % ----------------use PICO parameters from Runtable-------------------
         PICO_opts.C1 = UserVar.PICO_C1;
         PICO_opts.gamTstar = UserVar.PICO_gam;
         % ----------------oceanic forcing data for PICO-----------------------
         switch UserVar.OceForcing
             case "ismip6"
                 if isempty(ISMIP6_obs)
                      ISMIP6_obs = load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/obs_oceanic_forcing_ISMIP6.mat");
                 end
                 PICO_opts.Sbasins = ISMIP6_obs.salt_ID_30;
                 PICO_opts.Tbasins = ISMIP6_obs.temp_ID_30;
         end
         % --------------------------------------------------------------------
         PICO_opts.MeshBoundaryCoordinates = CtrlVar.MeshBoundaryCoordinates;
         [Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(1,CtrlVar,MUA,GF,h,median(rho),rhow,PICO_opts);
         ab=Mk;

    case "LQ"
         % -------------use Local Quadratic parameter from Runtable----------------
         gamma0 = UserVar.LQ_gam;
         % ------------------------------------------------------------------------
         if isempty(FdeltaT_basin)
            ncfile = UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/coeff_gamma0_DeltaT_quadratic_local_median.nc";
            deltaT_basin = ncread(ncfile,'deltaT_basin');

            x = double(ncread(ncfile,'x')); % units kms, projection: 
            y = double(ncread(ncfile,'y'));

            [X,Y] = ndgrid(x,y);
            FdeltaT_basin = scatteredInterpolant(X(:),Y(:),deltaT_basin(:),'linear');
         end

         deltaT_basin = FdeltaT_basin(MUA.coordinates(:,1),MUA.coordinates(:,2));

         rhoi_SI=918.0; % ice density (kg/m^3)
         rhosw_SI=1028.0; % sea water density
         rhofw_SI=1000.0; % freshwater density
         Lf_SI=3.34e5; % fusion latent heat of Ice (J/kg)
         cpw_SI=3974.0; % specific heat of sea water (J/kg/K)
         % ------------- oceanic forcing data for LQ-----------------------
         switch UserVar.OceForcing
             case "ismip6"
                 if isempty(ISMIP6_obs_tf)
                    ISMIP6_obs_tf= load(UserVar.datafolder+"/ANT_InputsForTransientSimulations/Basal_Melt/TF_climatology_ISMIP6.mat");
                 end
                 thermal_forcing = ISMIP6_obs_tf.TF([MUA.coordinates b]);
         end
         % ------------------------------------------------------------------------
         ab=-gamma0.*(rhosw_SI.*cpw_SI./rhoi_SI./Lf_SI).^2.*(max(thermal_forcing+deltaT_basin,0.0)).^2; % here is mass balance, which needs the negative sign
         
    case "PLUME"

        error("TO DO: add plume model to DefineMassBalance.");

     otherwise
     
         error("Unknown UserVar.BasalMelt: "+UserVar.BasalMelt);
end

NodesDownstreamOfGroundingLines="Relaxed"; % strict
[LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF,[],NodesDownstreamOfGroundingLines);
ab(LakeNodes) = 0;
ab(~OceanNodes) = 0;

switch UserVar.SMB

    case "MAR"
    
        if isempty(Fsmb_MAR)
	        ncfile = UserVar.datafolder+"/ANT_InputsForTransientSimulations/SMB/year-MAR_NorESM-1980-2100_zen.nc2";
	        fprintf(" Reading MAR SMB data from file %s \n",ncfile);
            % grid
	        x = double(ncread(ncfile,'X'))*1000; % convert units kms to m, projection: 
	        y = double(ncread(ncfile,'Y'))*1000;
            % time
	        t = ncread(ncfile,'TIME'); % days since 1901-01-15 00:00:00
            epochnum = datenum('1901-01-15 00:00:00');
	        t_years = year(double(epochnum+t));
            % remove years between start and end time of Ua simulation,
            % with buffer of 1 year on either side
            Ind_toremove = find(t_years<(UserVar.StartTime_DecimalYears-1) | t_years>(UserVar.EndTime_DecimalYears+1));
            t_years(Ind_toremove)=[];
            % smb
	        smb = ncread(ncfile,'SMB'); %annual averages, units kg m-2 yr-1, missing value -1.e+34
	        smb = squeeze(smb);
            smb(:,:,Ind_toremove)=[];
            % interpolant
            %% CLIMATOLOGY
            [X,Y] = ndgrid(x,y);
            smb_clim = mean(smb,3);
            Inan = find(~isnan(smb_clim));
            Fsmb_MAR = scatteredInterpolant(X(Inan),Y(Inan),smb_clim(Inan),'linear');
            %% ANNUALLY VARYING SMB
            %[X,Y,T] = ndgrid(x,y,t_years);
            %Inan = find(~isnan(smb));
            %Fsmb_MAR = scatteredInterpolant(X(Inan),Y(Inan),T(Inan),smb(Inan),'linear');
        end

	    as = Fsmb_MAR(MUA.coordinates(:,1),MUA.coordinates(:,2)); % in kg m-2 yr-1
        as = as./rho; % in m yr-1 ice equivalent
        % remove nans
        as = inpaint_nans(as,4);

	    fprintf("Done. \n");
	  
    otherwise
    
        error("Unknown UserVar.SMB: "+UserVar.SMB);
   
end



end


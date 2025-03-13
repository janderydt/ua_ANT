function out = loadvelocitydata(dataformat,years,only_grounded_ice)

% load measured velocity data
%tmp=load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
%Fus_2000=tmp.Fus; Fvs_2000=tmp.Fvs; Fxerr_2000=tmp.Fxerr; Fyerr_2000=tmp.Fyerr;
addpath(getenv("froot_data")+"Measures/Measures_annual");

% assemble velocity fields
years_tmp = strjoin(years,"-");
years_unique = unique(split(years_tmp,"-"));

for ii=1:numel(years_unique)

    % years_measures = string(years(ii))+"_"+string(years(ii)+1);
    % [vx,x_meas,y_meas]=measures_annual("vx",years_measures); 
    % vy=measures_annual("vy",years_measures);
    % stdx=measures_annual("vxstd",years_measures);
    % stdy=measures_annual("vystd",years_measures);
    % y_meas = flip(y_meas,1); [Xm,Ym]=ndgrid(x_meas,y_meas);
    % v(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(vx,vy),2));
    % std(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(stdx,stdy),2));
    if years_unique(ii)=="2000"
        fname = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    elseif years_unique(ii)=="2020"
        fname = "GriddedInterpolants_"+string(double(years_unique(ii))-1)+"-"+string(double(years_unique(ii)))+...
            "_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    else
        fname = "GriddedInterpolants_"+string(double(years_unique(ii)))+"-"+string(double(years_unique(ii))+1)+...
            "_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    end
    load("../ANT_Data/ANT_Interpolants/"+fname);
    v_tmp = hypot(Fus.Values,Fvs.Values);
    Fu = Fus; Fu.Values = v_tmp;
    std_tmp = hypot(Fxerr.Values,Fyerr.Values);
    Fstd = Fus; Fstd.Values = std_tmp;
    v.("yr"+years_unique(ii)).F = Fu;
    std.("yr"+years_unique(ii)).F = Fstd;

end

for dd=1:numel(dataformat)

    if dataformat(dd) == "du"
        
        years_tmp = split(years(dd),"-");
        yr1 = years_tmp(1);
        yr2 = years_tmp(2);

        load("Delta_u_AMUND_Weertman_"+yr1+"-"+yr2+".mat","MUA_yr2","GF_yr2");
        MUA = MUA_yr2; GF = GF_yr2;
    
        % interpolate initial and final speed onto Ua mesh
        u_init = v.("yr"+yr1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_init = std.("yr"+yr1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        u_target = v.("yr"+yr2).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_target = std.("yr"+yr2).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
    
        deltau = u_target-u_init; % measured change in speed (m/yr)
        deltau_err = hypot(std_init,std_target); % error estimate (m/yr)
        
        % remove nans
        Ind = find(isnan(deltau) | isnan(deltau_err));
        deltau(Ind) = eps(0);
        deltau_err(Ind) = 1e3;        
    
        %% what if we remove the floating ice?
        if only_grounded_ice
            deltau(GF.node<0.5) = eps(0);
            deltau_err(GF.node<0.5) = 1e-3;
        end
    
        out.("yr"+yr1+"_yr"+yr2).du = deltau(:);
        out.("yr"+yr1+"_yr"+yr2).stddu = deltau_err(:);

    elseif ismember(dataformat(dd),["u","LOGu"])
    
        load("u_AS_Calv_dh_cycle2_Weertman_"+years(dd)+".mat","MUA","GF");
        MUA = MUA.("yr"+years(dd)); GF = GF.("yr"+years(dd));

        % interpolate initial and final speed onto Ua mesh
        u_tmp = v.("yr"+years(dd)).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_tmp = std.("yr"+years(dd)).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        
        % remove nans
        %Ind = find(isnan(u_tmp) | isnan(std_tmp));
        %std_tmp(Ind) = 5e3;
        %u_tmp(Ind) = eps(0);

        % apply log if needed
        if dataformat(dd)=="LOGu"
            std_tmp = min(std_tmp./(u_tmp*log(10)),10);
            u_tmp = log10(u_tmp);    
        end

        out.("yr"+years(dd)).u = u_tmp(:);
        out.("yr"+years(dd)).stdu = std_tmp(:);
    
    end
end

end
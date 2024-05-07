function [dhdt,dhdt_err,dhdt_source] = Calc_ITSLIVE_CPOM_Paolo2023_dhdt(X,Y,targettime)

% X and Y are 2d matrices (ndgrid structure) for which dhdt will be returned
% 
% outputs will be returned for the targettime, which is a single number
% (datetime format 'dd-MMM-yyyy')

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

%% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %%
%% Read original data and linearly interpolate onto grid %%
%% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %%

%% ICE THICKNESS ELEVATION CHANGE
% 
% Original data grounded ice:
%----------------------------
% Otosaka et al. (2023) 
% Surface elevation change of the Amundsen Sea Embayment 1992-2019, CPOM
% https://zenodo.org/records/8117577 [Accessed on 14 November 2023]
%
% Original data floating ice:
%----------------------------
% Paolo et al. (2023)
% MEaSUREs ITS_LIVE Antarctic Ice Shelf Height Change and Basal Melt Rates, Version 1, National Snow and Ice Data Center 
% https://doi.org/10.5067/SE3XH9RXQWAM [Accessed on 14 November 2023]
% 
%
%% Otosaka et al. (2023) 
%
fprintf("Reading data from Otosaka et al. (2023).\n");
%
ncfileo = froot_data+"CPOM_dhdt/CPOM_altimetry_AIS_annual_dh_grids.nc";
%
xo = ncread(ncfileo,'x');
yo = ncread(ncfileo,'y');
timeo = ncread(ncfileo,'dhdt_time');
%
dhdto = ncread(ncfileo,'dhdt'); 
dhdto_std = ncread(ncfileo,'dhdt_total_err');
% In which 3-year wide window does the targettime fall?
years = floor(timeo);
Iwindow = find(years==year(targettime));

[Xo,Yo]=ndgrid(double(xo(:,1)),double(yo(1,:)));
dhdto = squeeze(dhdto(Iwindow,:,:));
Fdhdto = griddedInterpolant(Xo,Yo,dhdto,'linear','none');
Fdhdto_std = Fdhdto; Fdhdto_std.Values = squeeze(dhdto_std(Iwindow,:,:));
%
R = maprefcells([X(1,1) X(end,1)],[Y(1,1) Y(1,end)],[size(X,1) size(X,2)]);
geotiffwrite("./GeoTiffFiles/dhdto.tif",Fdhdto(X,Y)',R,'CoordRefSysCode','EPSG:3031');
%
%
%
%
%% Paolo et al. (2023)
%
fprintf("Reading data from Paolo et al. (2023).\n");
%
ncfilep = froot_data+"/Antarctica_dhdt/Paolo_2023/ANT_G1920V01_IceShelfMelt.nc";
%
xp=ncread(ncfilep,'x');
yp=flipdim(ncread(ncfilep,'y'),1);
[Xp,Yp] = ndgrid(unique(xp),unique(yp));
timep = ncread(ncfilep,'time'); % days since 1950-01-01
timep = datenum('1950-01-01 00:00')+timep;
% obtain fractional date
d=datetime(timep,'convertfrom','datenum');
startOfYear = dateshift(d, 'start', 'year');
endOfYear =  dateshift(d, 'end', 'year');
frac = (d-startOfYear)./(endOfYear-startOfYear);
timep_dec = year(d)+frac;
%
hp = ncread(ncfilep,'thickness');
hp = flipdim(permute(hp,[3 1 2]),3);
hp_err = ncread(ncfilep,'thickness_err');
hp_err = flipdim(permute(hp_err,[3 1 2]),3);
% calculate change in thickness and corresponding timestamp
dhp = hp(2:end,:,:)-hp(1:end-1,:,:);
dtp = timep_dec(2:end)-timep_dec(1:end-1);
dhdtp = dhp./repmat(dtp,1,size(dhp,2),size(dhp,3));
timep = (timep_dec(1:end-1)+timep_dec(2:end))/2;
% 3 year moving windows to be consistent with Otosaka
dhdtp_movmean = movmean(dhdtp,3*12,1);
%timep_movmean = movmean(timep,12*3,1);
[~,Itime] = min(abs(timep-(year(targettime)+month(targettime)/12)));
% interpolants
[Xp,Yp]=ndgrid(double(xp),double(yp));
dhdtp = squeeze(dhdtp_movmean(Itime,:,:));
Fdhdtp = griddedInterpolant(Xp,Yp,dhdtp,'linear','none');
%Fdhdtp_std = Fdhdtp; Fdhdtp_std.Values = squeeze(dhdt_movmean_std(Itime,:,:));
%
geotiffwrite("./GeoTiffFiles/dhdtp.tif",Fdhdtp(X,Y)',R,'CoordRefSysCode','EPSG:3031');
%
filename = "dhdt_"+string(year(targettime))+"_"+string(year(targettime)+1)+".mat";
save(filename,"Xo","Yo","Xp","Yp","dhdto","dhdtp");
%
clear hp hp_err dhp dhdtp dhdtp_movmean Xp Yp dhdto dhdto_std Xo Yo



return
%
%
%% Combine data for floating and grounded ice 
%
dhdt = zeros(size(X));
dhdt_err = zeros(size(X));
dhdt_source = zeros(size(X),'uint8');
%
dhdt_source_tmp = zeros(size(X),'uint8');
dhdt_source_string = '0: no data';
% indices of missing data: 
indo = ~isfinite(dso_ref) | ~isfinite(wo);
indp = ~isfinite(dsp_ref) | ~isfinite(wp);
% remove nans for calculation of averages
wo(isnan(wo)) = 0; dso_tmp = dso_ref; dso_tmp(isnan(dso_tmp))=0;
wn(isnan(wn)) = 0; dsn_tmp = dsn_ref; dsn_tmp(isnan(dsn_tmp))=0;
wp(isnan(wp)) = 0; dsp_tmp = dsp_ref; dsp_tmp(isnan(dsp_tmp))=0;
dso_err_tmp = dso_err_ref; dso_err_tmp(isnan(dso_err_tmp))=0;
dsn_err_tmp = dsn_err_ref; dsn_err_tmp(isnan(dsn_err_tmp))=0;
dsp_err_tmp = dsp_err_ref; dsp_err_tmp(isnan(dsp_err_tmp))=0;

% Replace the completely unknown values with NaNs: 
ds_tmp(indo & indn & indp) = nan;
ds_err_tmp(indo & indn & indp) = nan;
% source matrix
ds_source_tmp(~indo & indn & indp) = 1; 
ds_source_string = [ds_source_string,', 1: data from Otosaka et al. (2023) (zenodo.org/records/8117577)'];
ds_source_tmp(indo & indn & ~indp) = 2;
ds_source_string = [ds_source_string,', 2: data from Paolo et al. (2023) (doi.org/10.5067/SE3XH9RXQWAM )'];
ds_source_tmp(~indo & ~indn & indp) = 3;
ds_source_string = [ds_source_string,', 3: Blend of data from Otosaka et al. (2023) and Paolo et al. (2023)'];

ds = reshape(ds_tmp,size(X));
ds_err = reshape(ds_err_tmp,size(X));
ds_source = reshape(ds_source_tmp,size(X));

clear ds_err_tmp ds_source_tmp ds_tmp dsn_err_tmp dso_err_tmp dsp_err_tmp

end
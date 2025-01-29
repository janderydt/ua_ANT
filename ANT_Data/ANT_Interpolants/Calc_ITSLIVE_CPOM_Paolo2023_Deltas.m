function [ds,ds_err,ds_source] = Calc_ITSLIVE_CPOM_Paolo2023_Deltas(X,Y,referencetime,targettime)

% X and Y are 2d matrices (ndgrid structure) for which ds will be returned
% 
% referencetime is a single number or 2d matrix (same size as X) with 
% timestamps of the original datapoints (datetime format 'dd-MMM-yyyy')
%
% cumulative ds will be returned between referencetime and targettime. The
% latter is a single number (datetime format 'dd-MMM-yyyy')

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

useOtosaka = 0;
useNilsson = 1;
usePaolo = 1;

%% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %%
%% Read original data and linearly interpolate onto grid %%
%% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| %%

%% SURFACE ELEVATION CHANGE
% 
% Original data grounded ice:
%----------------------------
% Otosaka et al. (2023) 
% Surface elevation change of the Amundsen Sea Embayment 1992-2019, CPOM
% https://zenodo.org/records/8117577 [Accessed on 14 November 2023]
%
% Nilsson et al. (2023)
% MEaSUREs ITS_LIVE Antarctic Grounded Ice Sheet Elevation Change, Version 1, National Snow and Ice Data Center
% https://doi.org/10.5067/L3LSVDZS15ZV [Accessed on 14 November 2023]
%
% Original data floating ice:
%----------------------------
% Paolo et al. (2023)
% MEaSUREs ITS_LIVE Antarctic Ice Shelf Height Change and Basal Melt Rates, Version 1, National Snow and Ice Data Center 
% https://doi.org/10.5067/SE3XH9RXQWAM [Accessed on 14 November 2023]
% 
%
doy_target = day(targettime,'dayofyear')./day(datetime(year(targettime),12,31),'dayofyear');
doy_ref = day(referencetime,'dayofyear')./day(datetime(year(referencetime),12,31),'dayofyear');
referenceyear = year(referencetime)+doy_ref;
referenceyear = 0*X + referenceyear; % ensure correct dimensions
%
%
if useOtosaka
    %% Otosaka et al. (2023) 
    %
    fprintf("Reading data from Otosaka et al. (2023).\n");
    %
    ncfileo = froot_data+"CPOM_dhdt/CPOM_altimetry_AIS_annual_dh_grids.nc";
    %
    xo = ncread(ncfileo,'x');
    yo = ncread(ncfileo,'y');
    timeo = ncread(ncfileo,'dh_time');
    %
    dso = ncread(ncfileo,'dh_mean'); % this is the cumulative elevation change
    dso_std = ncread(ncfileo,'dh_stddev');
    % calculate elevation change ds with respect to targettime
    int = double(interp1(timeo,[1:numel(timeo)],year(targettime)+doy_target));
    fprintf("  > The target year is %s (%s-%s), which has index %s.\n",string(year(targettime)+doy_target),...
        string(timeo(floor(int))),string(timeo(ceil(int))),string(int));
    % calculate ds for target year (=dso_datum)
    fprintf("  > Calculate ds for target year and subtract from full timeseries.\n");
    a = year(targettime)+doy_target-timeo(floor(int));
    b = timeo(ceil(int))-timeo(floor(int));
    dso_datum = (1-a./b).*dso(floor(int),:,:) + a./b.*dso(ceil(int),:,:);
    dso_datum_std = sqrt((1-a./b).^2.*dso_std(floor(int),:,:).^2+...
        (a./b).^2.*dso_std(ceil(int),:,:).^2);
    % subtract ds_datum from timeseries (=dso_shifted)
    dso_shifted = dso - repmat(dso_datum,size(dso,1),1);
    dso_shifted_std = sqrt(dso_std.^2+repmat(dso_datum_std,size(dso,1),1).^2);
    % now interpolate ds to reference times    
    fprintf("  > Interpolate to reference times between %s and %s.\n",...
        string(min(referenceyear(:))),string(max(referenceyear(:))));
    [To,Xo,Yo]=ndgrid(double(timeo),double(xo(:,1)),double(yo(1,:)));
    Fdso = griddedInterpolant(To,Xo,Yo,dso_shifted,'linear','none');
    Fdso_std = Fdso; Fdso_std.Values = dso_shifted_std;
    dso_ref = Fdso(referenceyear,X,Y);
    dso_err_ref = Fdso_std(referenceyear,X,Y);
    %
    %R = maprefcells([X(1,1) X(end,1)],[Y(1,1) Y(1,end)],[size(X,1) size(X,2)]);
    %geotiffwrite("./GeoTiffFiles/dso_err.tif",flipdim(dso_err',1),R,'CoordRefSysCode','EPSG:3031');
    %
    clear Fdso Fdso_std dso dso_std To Xo Yo dso_shifted dso_shifted_std;
else
    dso_ref = nan*X;
    dso_err_ref = nan*X;
end
%
%
if useNilsson
    %% Nilsson et al. (2023)
    %
    fprintf("Reading data from Nilsson et al. (2023).\n");
    %
    ncfilen = froot_data+"ITS_LIVE/dhdt/ANT_G1920_GroundedIceHeight_v01.nc";
    % 
    xn = ncread(ncfilen,'x');
    yn = flipdim(ncread(ncfilen,'y'),1);
    [Xn,Yn] = ndgrid(xn,yn);
    timen = ncread(ncfilen,'time'); % days since 1950-01-01
    timen = datenum('1950-01-01 00:00')+timen;
    % obtain fractional date
    d=datetime(timen,'convertfrom','datenum');
    startOfYear = dateshift(d, 'start', 'year');
    endOfYear =  dateshift(d, 'end', 'year');
    frac = (d-startOfYear)./(endOfYear-startOfYear);
    timen_dec = year(d)+frac;
    %
    %dsn_scale = 0.01;
    dsn = ncread(ncfilen,'height_change'); % this is the cumulative elevation change relative to 16-12-2013
    dsn(dsn==-32767) = nan;
    dsn = flip(permute(dsn,[3 1 2]),3);
    %dsn_err_scale = 0.001;
    dsn_err = ncread(ncfilen,'height_change_rmse');
    dsn_err(dsn_err==-32767) = nan;
    dsn_err = flip(permute(dsn_err,[3 1 2]),3);
    % calculate elevation change ds with respect to targettime
    int = double(interp1(timen_dec,1:numel(timen_dec),year(targettime)+doy_target));
    fprintf("  > The target year is %s (%s-%s), which has index %s.\n",string(year(targettime)+doy_target),...
        string(timen_dec(floor(int))),string(timen_dec(ceil(int))),string(int));
    % calculate ds for target year (=dso_datum)
    fprintf("  > Calculate ds for target year and subtract from full timeseries.\n");
    a = year(targettime)+doy_target-timen_dec(floor(int));
    b = timen_dec(ceil(int))-timen_dec(floor(int));
    dsn_datum = (1-a./b).*dsn(floor(int),:,:) + a./b.*dsn(ceil(int),:,:);
    dsn_datum_std = sqrt((1-a./b).^2.*dsn_err(floor(int),:,:).^2+...
        (a./b).^2.*dsn_err(ceil(int),:,:).^2);
    dsn_shifted = dsn - repmat(dsn_datum,size(dsn,1),1);
    dsn_shifted_std = sqrt(dsn_err.^2+repmat(dsn_datum_std,size(dsn,1),1).^2);
    % now interpolate ds to reference times
    fprintf("  > Interpolate to reference times between %s and %s.\n",...
        string(min(referenceyear(:))),string(max(referenceyear(:))));
    [Tn,Xn,Yn]=ndgrid(double(timen_dec),double(xn),double(yn));
    Fdsn = griddedInterpolant(Tn,Xn,Yn,dsn_shifted,'linear','none');    
    Fdsn_std = Fdsn; Fdsn_std.Values = dsn_shifted_std;
    dsn_ref = Fdsn(referenceyear,X,Y);
    dsn_err_ref = Fdsn_std(referenceyear,X,Y);
    %
    %R = maprefcells([X(1,1) X(end,1)],[Y(1,1) Y(1,end)],[size(X,1) size(X,2)]);
    %geotiffwrite("./GeoTiffFiles/dsn_err.tif",flipdim(dsn_err',1),R,'CoordRefSysCode','EPSG:3031');
    %
    clear dsn dsn_err Fdsn Fdsn_std Tn Xn Yn dsn_shifted dsn_shifted_std
else
    dsn_ref = nan*X;
    dsn_err_ref = nan*X;
end
%
%
if usePaolo
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
    dsp = ncread(ncfilep,'height_change'); % cumulative elevation change with respect to 16-12-2013
    dsp = flip(permute(dsp,[3 1 2]),3);
    dsp_err = ncread(ncfilep,'height_change_err');
    dsp_err = flip(permute(repmat(dsp_err,1,1,size(dsp,1)),[3 1 2]),3);
    % calculate elevation change ds with respect to targettime
    int = double(interp1(timep_dec,1:numel(timep_dec),year(targettime)+doy_target));
    if ~isnan(int)
        fprintf("  > The target year is %s (%s-%s), which has index %s.\n",string(year(targettime)+doy_target),...
        string(timep_dec(floor(int))),string(timep_dec(ceil(int))),string(int));
        % calculate ds for target year (=dso_datum)
        
    elseif isnan(int)
        int = double(interp1(timep_dec,1:numel(timep_dec),year(targettime)+doy_target,'nearest','extrap'));
        fprintf("  > Requested target year %s falls outside the range of Paolo data (%s-%s),"+...
        " and values for %s are used instead.\n",string(year(targettime)+doy_target),string(timep_dec(1)),string(timep_dec(end)),string(timep_dec(int)));
    end
    fprintf("  > Calculate ds for target year and subtract from full timeseries.\n");
    a = year(targettime)+doy_target-timep_dec(floor(int));
    b = timep_dec(ceil(int))-timep_dec(floor(int));
    if b~=0
        dsp_datum = (1-a./b).*dsp(floor(int),:,:) + a./b.*dsp(ceil(int),:,:);
        dsp_datum_std = sqrt((1-a./b).^2.*dsp_err(floor(int),:,:).^2+...
            (a./b).^2.*dsp_err(ceil(int),:,:).^2);
    else
        dsp_datum = dsp(floor(int),:,:);
        dsp_datum_std = sqrt(dsp_err(floor(int),:,:).^2);
    end
    dsp_shifted = dsp - repmat(dsp_datum,size(dsp,1),1);
    dsp_shifted_std = sqrt(dsp_err.^2+repmat(dsp_datum_std,size(dsp,1),1).^2);
    % now interpolate ds to reference times
    fprintf("  > Interpolate to reference times between %s and %s.\n",...
        string(min(referenceyear(:))),string(max(referenceyear(:))));
    [Tp,Xp,Yp]=ndgrid(double(timep_dec),double(xp),double(yp));
    Fdsp = griddedInterpolant(Tp,Xp,Yp,dsp_shifted,'linear','none');
    Fdsp_std = Fdsp; Fdsp_std.Values = dsp_shifted_std;
    dsp_ref = Fdsp(referenceyear,X,Y);
    dsp_err_ref = Fdsp_std(referenceyear,X,Y);
    %
    clear dsp dsp_err Fdsp Fdsp_std Tp Xp Yp dsp_shifted dsp_shifted_std
    %
else
    dsp_ref = nan*X;
    dsp_err_ref = nan*X;
end
%
%
%% Combine data for floating and grounded ice 
%
ds = zeros(size(dsp_ref));
ds_err = zeros(size(dsp_ref));
ds_source = zeros(size(dsp_ref),'uint8');
%
ds_source_tmp = zeros(size(X),'uint8');
ds_source_string = '0: no data';
% define weights
wo = 1./dso_err_ref.^2;
wn = 1./dsn_err_ref.^2; 
wp = 1./dsp_err_ref.^2;
% indices of missing data: 
indo = ~isfinite(dso_ref) | ~isfinite(wo);
indn = ~isfinite(dsn_ref) | ~isfinite(wn);
indp = ~isfinite(dsp_ref) | ~isfinite(wp);
% remove nans for calculation of averages
wo(isnan(wo)) = 0; dso_tmp = dso_ref; dso_tmp(isnan(dso_tmp))=0;
wn(isnan(wn)) = 0; dsn_tmp = dsn_ref; dsn_tmp(isnan(dsn_tmp))=0;
wp(isnan(wp)) = 0; dsp_tmp = dsp_ref; dsp_tmp(isnan(dsp_tmp))=0;
dso_err_tmp = dso_err_ref; dso_err_tmp(isnan(dso_err_tmp))=0;
dsn_err_tmp = dsn_err_ref; dsn_err_tmp(isnan(dsn_err_tmp))=0;
dsp_err_tmp = dsp_err_ref; dsp_err_tmp(isnan(dsp_err_tmp))=0;
% weighted averages
ds_tmp = (wo.*dso_tmp + wn.*dsn_tmp + wp.*dsp_tmp)./(wo + wn + wp + eps);
ds_err_tmp = sqrt((dso_err_tmp.*wo).^2+...
    (dsn_err_tmp.*wn).^2+(dsp_err_tmp.*wp).^2)./(wo + wn + wp + eps);
% Replace the completely unknown values with NaNs: 
ds_tmp(indo & indn & indp) = nan;
ds_err_tmp(indo & indn & indp) = nan;
% source matrix
ds_source_tmp(~indo & indn & indp) = 1; 
ds_source_string = [ds_source_string,', 1: data from Otosaka et al. (2023) (zenodo.org/records/8117577)'];
ds_source_tmp(indo & ~indn & indp) = 2; 
ds_source_string = [ds_source_string,', 2: data from Nilsson et al. (2023) (doi.org/10.5067/L3LSVDZS15ZV)'];
ds_source_tmp(indo & indn & ~indp) = 3;
ds_source_string = [ds_source_string,', 3: data from Paolo et al. (2023) (doi.org/10.5067/SE3XH9RXQWAM )'];
ds_source_tmp(~indo & indn & ~indp) = 4;
ds_source_string = [ds_source_string,', 4: Blend of data from Otosaka et al. (2023) and Paolo et al. (2023)'];
ds_source_tmp(~indo & ~indn & indp) = 5;
ds_source_string = [ds_source_string,', 5: Blend of data from Otosaka et al. (2023) and Nilsson et al. (2023)'];
ds_source_tmp(indo & ~indn & ~indp) = 6;
ds_source_string = [ds_source_string,', 6: Blend of data from Nilsson et al. (2023) and Paolo et al. (2023)'];
ds_source_tmp(~indo & ~indn & ~indp) = 7; % 
ds_source_string = [ds_source_string,', 7: Blend of data from Otosaka et al. (2023), Nilsson et al. (2023) and Paolo et al. (2023)'];

ds = reshape(ds_tmp,size(X));
ds_err = reshape(ds_err_tmp,size(X));
ds_source = reshape(ds_source_tmp,size(X));

clear ds_err_tmp ds_source_tmp ds_tmp dsn_err_tmp dso_err_tmp dsp_err_tmp

end
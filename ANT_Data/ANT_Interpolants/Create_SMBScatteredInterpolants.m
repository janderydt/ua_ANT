function Create_SMBScatteredInterpolants(CreateGeotiff)

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

if nargin==0
    CreateGeotiff=1;
end

years_to_average = [2000:2020];

%% MAR
ncfile = froot_data+"/SMB_RACMO_MAR/MAR_Kittel2021/year-MAR_ERA5-1979-2019_zen.nc2";
fprintf(" Reading MAR SMB data from file %s \n",ncfile);

x = double(ncread(ncfile,'X'))*1e3; % units kms to m
y = double(ncread(ncfile,'Y'))*1e3;
t = ncread(ncfile,'TIME'); %days since 1976-03-01 00:00:00
smb = ncread(ncfile,'SMB'); %annual averages, units kg m-2 yr-1, missing value -1.e+34
smb(smb==-1e34)=nan;

fprintf("Done. \n");

fprintf("Creating scatteredInterpolants %s \n",ncfile);

epochnum = datenum('1976-03-01 00:00:00');
t = double(epochnum+t);
smb = squeeze(double(smb)/1000); % Ua wants m/yr in freshwater equivalent

[X,Y] = ndgrid(x,y);
Fsmb_MAR.years = [];
for ii=1:numel(t)
    field = ['yr',datestr(t(ii),'yyyy')];
    smb_tmp = squeeze(smb(:,:,ii));
    Fsmb_MAR.years = [Fsmb_MAR.years; datestr(t(ii),'yyyy')];
    Fsmb_MAR.(field) = griddedInterpolant(X,Y,smb_tmp,'cubic','nearest'); %scatteredInterpolant(X(:),Y(:),T(:),smb_year(:),'linear');
end

I = find(double(string(Fsmb_MAR.years))>=years_to_average(1) & double(string(Fsmb_MAR.years))<=years_to_average(end));
smb_climatology = mean(smb(:,:,I),3,"omitmissing");
Fsmb_MAR.climatology = griddedInterpolant(X,Y,smb_climatology,'cubic','nearest');

fprintf("Done. \n");

if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');
    
    R = maprefcells([x(1) x(end)],[y(1) y(end)],[numel(x),numel(y)]);
    smb_mean = mean(smb(:,:,I),3);
    geotiffwrite("./GeoTiffFiles/smb_MAR_mean"+string(years_to_average(1))+"-"+string(years_to_average(end))+".tif",...
        smb_mean',R,'CoordRefSysCode','EPSG:3031');

    fprintf('Done.\n');

end

%% RACMO
% yearly data
% ncfile = froot_data+"/SMB_RACMO_MAR/RACMO2.3/SMB_RACMO2.3p2_yearly_ANT27_1979_2016.nc";
% monthly data
ncfile = froot_data+"/SMB_RACMO_MAR/RACMO2.3/smb_monthlyS_ANT27_ERA5-3H_RACMO2.3p2_197901_202212.nc";

fprintf(" Reading RACMO SMB data from file %s \n",ncfile);

lon = ncread(ncfile,'lon'); 
lat = ncread(ncfile,'lat');
t = ncread(ncfile,'time'); %days since 1950-01-01 00:00:00.0
height =  ncread(ncfile,'height'); %height above the surface
%smb = ncread(ncfile,'smb'); % yearly average rates, units kg m-2 s-1, missing value -9999
smb = ncread(ncfile,'smb'); % monthly integrated values, units kg m-2, missing value -9999
smb(smb==-9999)=nan;

fprintf("Done. \n");

fprintf("Creating scatteredInterpolants %s \n",ncfile);

[x,y] = ll2psxy(lat,lon,-71,0);
epochnum = datenum('1950-01-01 00:00:00');
t = double(epochnum+t);
%smb_year = smb/1000*(365.25*24*60*60); % Ua want m/yr in freshwater equivalent
smb_month = squeeze(smb/1000); % in m/mo
% add monthly values in any given year to obtain yearly rates:
smb_year = [];
for ii=0:numel(t)/12-1
    smb_year(:,:,ii+1) = sum(smb_month(:,:,12*ii+1:12*(ii+1)),3);
    years(ii+1) = round(mean(t(12*ii+1:12*(ii+1))));
end
%t = years;
% X = repmat(x,[1,1,numel(t)]);
% Y = repmat(y,[1,1,numel(t)]);
% T = repmat(t,[1,size(x,1),size(x,2)]); T = permute(T,[2 3 1]);

Fsmb_RACMO.years = [];
for ii=1:numel(years)
    field = ['yr',datestr(years(ii),'yyyy')];
    smb_tmp = smb_year(:,:,ii);
    Fsmb_RACMO.years = [Fsmb_RACMO.years; datestr(years(ii),'yyyy')];
    Fsmb_RACMO.(field) = scatteredInterpolant(x(:),y(:),smb_tmp(:),'natural','nearest'); %scatteredInterpolant(X(:),Y(:),T(:),smb_year(:),'linear');
end

I = find(year(years)>=years_to_average(1) & year(years)<=years_to_average(end));
smb_climatology = mean(smb_year(:,:,I),3,"omitmissing");
Fsmb_RACMO.climatology = scatteredInterpolant(x(:),y(:),smb_climatology(:),'natural','nearest');

fprintf("Done. \n");

fprintf(' Saving SMB interpolants with Fsmb_MAR and Fsmb_RACMO \n')

save('ScatteredInterpolants_SMB','Fsmb_MAR','Fsmb_RACMO','-v7.3');

if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');

    smb_mean = mean(smb(:,:,I),3);

    R = georasterref('RasterSize',size(smb_mean), ...
        'LatitudeLimits',[min(lat(1,:)) max(lat(1,:))],'LongitudeLimits',[min(lon(:,1)) max(lon(:,1))]);
       
    geotiffwrite("./GeoTiffFiles/smb_RACMO_mean"+string(years_to_average(1))+"-"+string(years_to_average(end))+".tif",...
        smb_mean',R);

    fprintf('Done.\n');

end


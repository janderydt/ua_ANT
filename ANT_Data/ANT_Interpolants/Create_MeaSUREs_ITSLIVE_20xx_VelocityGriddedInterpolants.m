function Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(year,CreateGeotiff)

%%  
% Creates gridded interpolants velocities based on a combination of 
% MEaSUREs and ITSLive data between 20xx and 20xx+1 with some gaps filled by
% MEaSUREs and ITSLive mosaics that span multiple years.

% step1.1 MEaSUREs annual 20xx-20xx+1 - remove data with errors > error_cutoff m/yr
% step1.2 ITSLIVE annual 20xx - remove data with errors > error_cutoff
% step1.3 weighted average of MEaSUREs and ITSLIVE data 

% step2 fill gaps with MEaSUREs 450m v2 and ITSLIVE 240m mosaic
% and remove data with errors > error_cutoff m/yr

% step3 fill remaining gaps with MEaSUREs Central Antarctica (dated 1997)
% and remove data with errors > error_cutoff m/yr
%

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));
addpath(getenv("froot_data")+"/Measures/Measures_annual");
addpath(getenv("froot_data")+"/ITS_LIVE/velocities");

relative_error_cutoff = 1; % ERR_V./V

if nargin==0
    year = 2019;
end

if nargin<2
    CreateGeotiff = 1;
end

%% STEP1
% Load MeaSUREs annual data
fprintf("Reading MeaSUREs velocity data for %s-%s",string(year),string(year+1));
[vxm,x,y] = measures_annual('vx',[num2str(year),'_',num2str(year+1)]);
vym = measures_annual('vy',[num2str(year),'_',num2str(year+1)]);
xerrm = measures_annual('vxerr',[num2str(year),'_',num2str(year+1)]); xerrm(isnan(xerrm))=0; % change nan to zero for calculation of  total error.
yerrm = measures_annual('vyerr',[num2str(year),'_',num2str(year+1)]); yerrm(isnan(yerrm))=0; % change nan to zero for calculation of  total error.
vm = hypot(vxm,vym);
errm = hypot(vxm./vm.*xerrm,vym./vm.*yerrm);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wm = 1./errm.^2; % weights

[X,Y] = ndgrid(x,y);
    
fprintf("...done.\n")

badm = errm./hypot(vxm,vym)>relative_error_cutoff;
vxm(badm) = nan; vym(badm) = nan; 

% Load ITSLIVE annual data
fprintf("Reading ITSLIVE velocity data for %s",string(year));
vxi = itslive_annual('vx',num2str(year),X,Y);
vyi = itslive_annual('vy',num2str(year),X,Y);
xerri = itslive_annual('vxerr',num2str(year),X,Y); xerri(isnan(xerri)) = 0;
yerri = itslive_annual('vyerr',num2str(year),X,Y); yerri(isnan(yerri)) = 0;
vi = hypot(vxi,vyi);
erri = hypot(vxi./vi.*xerri,vyi./vi.*yerri);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wi = 1./erri.^2; % weights
    
fprintf("...done.\n")

badi = erri./hypot(vxi,vyi)>relative_error_cutoff;
vxi(badi) = nan; vyi(badi) = nan; 
% Indices of missing data: 
indm = ~isfinite(vxm) | ~isfinite(vym) | ~isfinite(wm); 
vxm(indm) = 0; 
vym(indm) = 0; 
wm(indm) = 0; 
indi = ~isfinite(vxi) | ~isfinite(vyi) | ~isfinite(wi); 
vxi(indi) = 0; 
vyi(indi) = 0; 
wi(indi) = 0;

% Weighted mean velocity components: 
vx = (vxi.*wi + vxm.*wm)./(wi + wm + eps); 
vy = (vyi.*wi + vym.*wm)./(wi + wm + eps); 
xerr = hypot(xerri.*wi,xerrm.*wm)./(wi + wm + eps);
yerr = hypot(yerri.*wi,yerrm.*wm)./(wi + wm + eps);

% Replace the completely unknown values with NaNs: 
vx(indi & indm) = nan; 
vy(indi & indm) = nan; 
xerr(indi & indm) = nan;
yerr(indi & indm) = nan;

vx = double(vx); 
vy = double(vy); 
xerr = double(xerr);
yerr = double(yerr);

source_count = 1;
v_source = zeros(size(vx),'uint8');
v_source(~indi & indm) = source_count; % ITS_LIVE
v_source(indi & ~indm) = source_count+1; % Measures v2
v_source(~indi & ~indm) = source_count+2; % Blend of ITS_LIVE and MEasures 
v_source_string(source_count) =  string(source_count)+": ITS_LIVE Antarctica annual "+string(year);
v_source_string(source_count+1) =  string(source_count+1)+": MeaSUREs Antarctica annual "+string(year)+"-"+string(year+1);
v_source_string(source_count+2) =  string(source_count+2)+": Blend of ITS_LIVE and MeaSUREs Antarctica annual "+string(year)+"-"+string(year+1);
source_count = source_count + 3;

nanind = isnan(vx);

%% STEP2
% MeaSUREs data
ncfile = froot_data+"/Measures/Measures_v2/antarctica_ice_velocity_450m_v2.nc";
fprintf("Reading Velocity data from MeaSURES 450m v2");

xtmp = double(ncread(ncfile,'x'));
ytmp = double(ncread(ncfile,'y'));
ytmp = flipdim(ytmp,1);
[Xtmp,Ytmp] = ndgrid(xtmp,ytmp);
xerrtmp = ncread(ncfile,'ERRX'); xerrtmp = flipdim(xerrtmp,2);
yerrtmp = ncread(ncfile,'ERRY'); yerrtmp = flipdim(yerrtmp,2);
vxtmp = ncread(ncfile,'VX'); vxtmp = flipdim(vxtmp,2);
vytmp = ncread(ncfile,'VY'); vytmp = flipdim(vytmp,2);
%vtmp = hypot(vxtmp,vytmp);
%errtmp = hypot(vxtmp./vtmp.*xerrtmp,vytmp./vtmp.*yerrtmp);
Fvx = griddedInterpolant(Xtmp,Ytmp,vxtmp,'cubic','none');
Fvy = griddedInterpolant(Xtmp,Ytmp,vytmp,'cubic','none');
Fxerr = griddedInterpolant(Xtmp,Ytmp,xerrtmp,'cubic','none');
Fyerr = griddedInterpolant(Xtmp,Ytmp,yerrtmp,'cubic','none');
            
vxm = Fvx(X,Y);
vym = Fvy(X,Y);
xerrm = Fxerr(X,Y); xerrm(isnan(xerrm))=0; % replace nans by 0 for calculation of total error 
yerrm = Fyerr(X,Y); yerrm(isnan(yerrm))=0; % replace nans by 0 for calculation of total error 

vm = hypot(vxm,vym);
errm = hypot(vxm./vm.*xerrm,vym./vm.*yerrm);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wm = 1./errm.^2;

vxm(errm./vm>relative_error_cutoff) = nan;
vym(errm./vm>relative_error_cutoff) = nan;

fprintf("...done.\n");

% Load ITSLIVE data
fprintf("Processing ITS_LIVE 240m data");

vxi = 0*X + itslive_annual('vx',"0000",X,Y);
vyi = 0*X + itslive_annual('vy',"0000",X,Y);
xerri = 0*X + itslive_annual('vxerr',"0000",X,Y); xerri(isnan(xerri)) = 0;
yerri = 0*X + itslive_annual('vyerr',"0000",X,Y); yerri(isnan(yerri)) = 0;
vi = hypot(vxi,vyi);
erri = hypot(vxi./vi.*xerri,vyi./vi.*yerri);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wi = 1./erri.^2; % weights

vxi(erri./vi>relative_error_cutoff) = nan;
vyi(erri./vi>relative_error_cutoff) = nan;

fprintf("...done.\n")

% Indices of missing data: 
indm = ~isfinite(vxm) | ~isfinite(vym) | ~isfinite(wm); 
vxm(indm) = 0; 
vym(indm) = 0; 
wm(indm) = 0; 
indi = ~isfinite(vxi) | ~isfinite(vyi) | ~isfinite(wi); 
vxi(indi) = 0; 
vyi(indi) = 0; 
wi(indi) = 0;

% Weighted mean velocity components: 
vxtmp = (vxi.*wi + vxm.*wm)./(wi + wm + eps); 
vytmp = (vyi.*wi + vym.*wm)./(wi + wm + eps); 
xerrtmp = hypot(xerri.*wi,xerrm.*wm)./(wi + wm + eps);
yerrtmp = hypot(yerri.*wi,yerrm.*wm)./(wi + wm + eps);

% Replace the completely unknown values with NaNs: 
vxtmp(indi & indm) = nan; 
vytmp(indi & indm) = nan; 
xerrtmp(indi & indm) = nan;
yerrtmp(indi & indm) = nan;

vxtmp = double(vxtmp); 
vytmp = double(vytmp); 
xerrtmp = double(xerrtmp);
yerrtmp = double(yerrtmp);

% Now add to final datasets
vx(nanind) = vxtmp(nanind);
vy(nanind) = vytmp(nanind);
xerr(nanind) = xerrtmp(nanind);
yerr(nanind) = yerrtmp(nanind);

v_source(~indi & indm & nanind) = source_count; % ITS_LIVE
fprintf("   > added %s data points from ITS_LIVE 240m mosaic.\n",string(numel(find(~indi & indm & nanind))));
v_source(indi & ~indm & nanind) = source_count+1; % Measures v2
fprintf("   > added %s data points from Measures 450m v2 mosaic.\n",string(numel(find(indi & ~indm & nanind))));
v_source(~indi & ~indm & nanind) = source_count+2; % Blend of ITS_LIVE and MEasures 
fprintf("   > added %s data points from ITS_LIVE + Measures mosaic.\n",string(numel(find(~indi & ~indm & nanind))));
v_source_string(source_count) =  string(source_count)+": ITS_LIVE Antarctica 240m mosaic ";
v_source_string(source_count+1) =  string(source_count+1)+": MeaSUREs Antarctica 450m v2 mosaic";
v_source_string(source_count+2) =  string(source_count+2)+": Blend of ITS_LIVE 240m and MeaSUREs 450m v2 mosaics";
source_count = source_count+3;

nanind = isnan(vx);

%% STEP3
% Fill gaps with MeaSUREs data Central Antarctica
ncfile = froot_data+"/Measures/Measures_CentralAntarctica/1997/Central_Antarctica_ice_velocity_1997.nc";
fprintf("Reading Velocity data from MeaSURES Central Antarctia");

nx = double(ncreadatt(ncfile,'/','nx')); ny = double(ncreadatt(ncfile,'/','ny'));
xmin = str2double(erase(ncreadatt(ncfile,'/','xmin'),[" ","m"])); ymax = str2double(erase(ncreadatt(ncfile,'/','ymax'),[" ","m"]));
dl = str2double(erase(ncreadatt(ncfile,'/','spacing'),[" ","m"]));
xtmp = double([xmin:dl:xmin+(nx-1)*dl]);
ytmp = double([ymax-(ny-1)*dl:dl:ymax]);
[Xtmp,Ytmp] = ndgrid(xtmp,ytmp);
errtmp = double(ncread(ncfile,'err')); errtmp = flipdim(errtmp,2);
errtmp(errtmp==0) = nan;
vxtmp = double(ncread(ncfile,'vx')); vxtmp = flipdim(vxtmp,2);
vxtmp(vxtmp==0) = nan;
vytmp = double(ncread(ncfile,'vy')); vytmp = flipdim(vytmp,2);
vytmp(vytmp==0) = nan;
Fvx = griddedInterpolant(Xtmp,Ytmp,vxtmp,'cubic','none');
Fvy = griddedInterpolant(Xtmp,Ytmp,vytmp,'cubic','none');
Ferr = griddedInterpolant(Xtmp,Ytmp,errtmp,'cubic','none');
            
vxm = Fvx(X,Y);
vym = Fvy(X,Y);
xerrm = Ferr(X,Y); xerrm(isnan(xerrm))=0; % replace nans by 0 for calculation of total error 
yerrm = xerrm; 

vm = hypot(vxm,vym);
errm = hypot(vxm./vm.*xerrm,vym./vm.*yerrm);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wm = 1./errm.^2;

vxm(errm./vm>relative_error_cutoff) = nan;
vym(errm./vm>relative_error_cutoff) = nan;

fprintf("...done.\n");

% Indices of missing data: 
indm = ~isfinite(vxm) | ~isfinite(vym) | ~isfinite(wm); 
vxm(indm) = nan; 
vym(indm) = nan; 
xerrm(indm) = nan;
yerrm(indm) = nan;

vxm = double(vxm);
vym = double(vym); 
xerrm = double(xerrm);
yerrm = double(yerrm);

% Now add to final datasets
vx(nanind) = vxm(nanind);
vy(nanind) = vym(nanind);
xerr(nanind) = xerrm(nanind);
yerr(nanind) = yerrm(nanind);

v_source(~indm & nanind) = source_count; % ITS_LIVE
fprintf("   > added %s data points from MeaSURES Central Antarctica 1997.\n",string(numel(find(~indm & nanind))));
v_source_string(source_count) =  string(source_count)+": MeaSURES Central Antarctica 1997";
source_count = source_count+1;

nanind = isnan(vx);

if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');
    
    R = maprefcells([x(1) x(end)],[y(end) y(1)],[numel(x),numel(y)]);
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_vx.tif",flipdim(vx',1),R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_vy.tif",flipdim(vy',1),R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_xerr.tif",flipdim(xerr',1),R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_yerr.tif",flipdim(yerr',1),R,'CoordRefSysCode','EPSG:3031');
    v = sqrt(vx.^2 + vy.^2); err = sqrt((vx./v.*xerr).^2+(vy./v.*yerr).^2);
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_err.tif",flipdim(err',1),R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_source.tif",flipdim(v_source',1),R,'CoordRefSysCode','EPSG:3031');

    fprintf('Done.\n');

end

fprintf(' Creating gridded interpolants...');

y = flipdim(y,1);
[X,Y] = ndgrid(x,y);

Fus = griddedInterpolant(X,Y,flipdim(vx,2),'linear');
Fvs = griddedInterpolant(X,Y,flipdim(vy,2),'linear');
Fxerr = griddedInterpolant(X,Y,flipdim(xerr,2),'linear');
Fyerr = griddedInterpolant(X,Y,flipdim(yerr,2),'linear');
Fsource = griddedInterpolant(X,Y,flipdim(double(v_source),2),'nearest');

fprintf('done.\n')

fprintf(' Saving VelocityInterpolants with Fus, Fvs, Fxerr, Fyerr. \n')
save("GriddedInterpolants_"+string(year)+"-"+string(year+1)+"_MeaSUREs_ITSLIVE_Velocities",'Fus','Fvs','Fxerr','Fsource','Fyerr','-v7.3');

return

%% Test gridded interpolants
% load in Sainan's ISMIP6 mesh of Antarctica and map on this mesh. 
% Just done for testing purposes and to see if the gridded data looks sensible.
fprintf(' Testing interpolants by mapping on a FE grid...')

load("./MeshAntISMIP6_2_SainanSun.mat"); 
xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 

usFE=Fus(xFE,yFE);
vsFE=Fvs(xFE,yFE);
uerrFE=Fuerr(xFE,yFE);

uFE=sqrt(usFE.^2 + vsFE.^2);

CtrlVar.PlotXYscale=1000;

ufig=FindOrCreateFigure('u') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,uFE);
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('u') ; title(cbar,'m/yr')

uerrfig=FindOrCreateFigure('uerr') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,uerrFE);
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('uerr') ; title(cbar,'m/yr')





%%  
% Reads nc files and creates Úa gridded interpolants with s, b, B, velocity and rho.
%
%

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

%%% VELOCITIES
ncfile = froot_data+"/BedMachine_Antarctica/vel_nsidc.CF16_2.nc";
fprintf(" Reading Velocity data from file %s \n",ncfile);

xdata = double(ncread(ncfile,'x'));
ydata = double(ncread(ncfile,'y'));

ERR_x  = double(ncread(ncfile,'ERRX'))';
ERR_y  = double(ncread(ncfile,'ERRY'))';
ERR = sqrt(ERR_x.^2 + ERR_y.^2);
vx  = double(ncread(ncfile,'VX'))';
vy =  double(ncread(ncfile,'VY'))';
V = sqrt(vx.^2+vy.^2);
    
fprintf("...done.\n")

xdata=unique(xdata);
ydata=unique(ydata);
[Xdata,Ydata]=ndgrid(xdata,ydata);

fprintf(' Creating gridded interpolants...')

Fus = griddedInterpolant(Xdata,Ydata,flipdim(vx',2),'cubic');
Fvs = griddedInterpolant(Xdata,Ydata,flipdim(vy',2),'cubic');
Fuerr = griddedInterpolant(Xdata,Ydata,flipdim(ERR',2),'cubic');
fprintf('done.\n')

fprintf(' Saving VelocityInterpolants with Fus, Fvs, Fuerr. \n')
tic
    save('GriddedInterpolants_MeaSUREs_Fus_2016','Fus','-v7.3','-nocompression');
    save('GriddedInterpolants_MeaSUREs_Fvs_2016','Fvs','-v7.3','-nocompression');
    save('GriddedInterpolants_MeaSUREs_Fuerr_2016','Fuerr','-v7.3','-nocompression');
toc

if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');
    
    R = makerefmat(xdata(1),ydata(1),xdata(2)-xdata(1),ydata(2)-ydata(1));
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2016.tif",V',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocityError_2016.tif",ERR',R,'CoordRefSysCode','EPSG:3031');

    fprintf('Done.\n');

end

return 

%% Test gridded interpolants
% load in Sainan's ISMIP6 mesh of Antarctica and map on this mesh. 
% Just done for testing purposes and to see if the gridded data looks sensible.
fprintf(' Testing interpolants by mapping on a FE grid...')

load("./MeshAntISMIP6_2_SainanSun.mat"); 
xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 

usFE=Fus(xFE,yFE);
vsFE=Fvs(xFE,yFE);
uerrFE=Fuerr(xFE,yFE);

uFE=sqrt(usFE.^2 + vsFE.^2);

CtrlVar.PlotXYscale=1000;

ufig=FindOrCreateFigure('u') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,uFE);
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('u') ; title(cbar,'m/yr')

uerrfig=FindOrCreateFigure('uerr') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,uerrFE);
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('uerr') ; title(cbar,'m/yr')

fprintf('done.\n')
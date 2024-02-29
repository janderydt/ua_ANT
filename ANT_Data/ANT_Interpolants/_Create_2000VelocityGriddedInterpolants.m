function Create2000VelocityGriddedInterpolants(use_Paolo_for_AmundsenSea,CreateGeotiff)

%%  
% Creates gridded interpolants for mosaic of early 2000s velocities.
% This is based on a combination of MEaSUREs and ITSLive data between 1996
% and 2003.

% step1.1 MEaSUREs annual 2000-2001 - remove data with errors > error_cutoff m/yr
% step1.2 ITSLIVE annual 2000-2001 - remove data with errors > error_cutoff
% step1.3 weighted average of MEaSUREs and ITSLIVE data 

% step2 fill gaps with with weighted average of MEaSUREs Amundsen Sea and 
% ITSLIVE data from 2000:1:2003 and 1996, removing data with errors > 
% error_cutoff m/yr at each step

% step3. fill remaining gaps with MEaSUREs 450m v2 and ITSLIVE 240m mosaic
% and remove data with errors > error_cutoff m/yr

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));
addpath(getenv("froot_data")+"/Measures/Measures_annual");
addpath(getenv("froot_data")+"/ITS_LIVE/velocities");
addpath(froot_data+"/Measures/Measures_AmundsenSea/ALOS, ERS-1, RADARSAT-1, RADARSAT-2, TDX");

relative_error_cutoff = 1; % ERR_V./V

if nargin==0
    use_Paolo_for_AmundsenSea = 0;
    CreateGeotiff=1;
elseif nargin ==1
    CreateGeotiff=1;
end

%% STEP1
% Load MeaSUREs annual data
fprintf("Reading MeaSUREs velocity data for 2000");
[vxm,x,y] = measures_annual('vx','2000_2001');
vym = measures_annual('vy','2000_2001');
xerrm = measures_annual('vxerr','2000_2001'); xerrm(isnan(xerrm))=0; % change nan to zero for calculation of  total error.
yerrm = measures_annual('vyerr','2000_2001'); yerrm(isnan(yerrm))=0; % change nan to zero for calculation of  total error.
vm = sqrt(vxm.^2+vym.^2);
errm = sqrt((vxm./vm.*xerrm).^2+(vym./vm.*yerrm).^2);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wm = 1./errm.^2; % weights

[X,Y] = ndgrid(x,y);
    
fprintf("...done.\n")

badm = errm./sqrt(vxm.^2+vym.^2)>relative_error_cutoff;
vxm(badm) = nan; vym(badm) = nan; 

% Load ITSLIVE annual data
fprintf("Reading ITSLIVE velocity data for 2000");
vxi = itslive_annual('vx','2000',X,Y);
vyi = itslive_annual('vy','2000',X,Y);
xerri = itslive_annual('vxerr','2000',X,Y); xerri(isnan(xerri)) = 0;
yerri = itslive_annual('vyerr','2000',X,Y); yerri(isnan(yerri)) = 0;
vi = sqrt(vxi.^2+vyi.^2);
erri = sqrt((vxi./vi.*xerri).^2+(vyi./vi.*yerri).^2);% propagation of error for velocity v = sqrt(vx^2+vy^2)
wi = 1./erri.^2; % weights
    
fprintf("...done.\n")

badi = erri./sqrt(vxi.^2+vyi.^2)>relative_error_cutoff;
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
xerr = sqrt((xerri.*wi).^2+(xerrm.*wm).^2)./(wi + wm + eps);
yerr = sqrt((yerri.*wi).^2+(yerrm.*wm).^2)./(wi + wm + eps);

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
v_source_string(source_count) =  string(source_count)+": ITS_LIVE Antarctica annual 2000";
v_source_string(source_count+1) =  string(source_count+1)+": MeaSUREs Antarctica annual 2000-2001";
v_source_string(source_count+2) =  string(source_count+2)+": Blend of ITS_LIVE and MeaSUREs Antarctica annual 2000-2001";
source_count = source_count + 3;

nanind = isnan(vx);

%% STEP2_DEFAULT
if  use_Paolo_for_AmundsenSea ==0

    years = ["2000","2001","2002","2003","1996"]; % provide in order you want the infill to happen
    
    for ii=1:numel(years)

        % load MeaSUREs data
        fprintf("Processing MeaSUREs Amundsen Sea data for %s",years(ii));

        vxm = 0*X+measures_amundsen(['vx',char(years(ii))],X,Y);
        vym = 0*X+measures_amundsen(['vy',char(years(ii))],X,Y);
        xerrm = 0*X+measures_amundsen(['err',char(years(ii))],X,Y); xerrm(isnan(xerrm))=0; % replace nans by 0 for calculation of total error 
        yerrm = yerrm; 

        vm = sqrt(vxm.^2+vym.^2);
        errm = xerrm;
        wm = 1./errm.^2;

        vxm(errm./vm>relative_error_cutoff) = nan;
        vym(errm./vm>relative_error_cutoff) = nan;

        fprintf("...done.\n");
        
        % Load ITSLIVE data
        fprintf("Processing ITS_LIVE data for %s",years(ii));

        vxi = 0*X + itslive_annual('vx',years(ii),X,Y);
        vyi = 0*X + itslive_annual('vy',years(ii),X,Y);
        xerri = 0*X + itslive_annual('vxerr',years(ii),X,Y); xerri(isnan(xerri)) = 0;
        yerri = 0*X + itslive_annual('vyerr',years(ii),X,Y); yerri(isnan(yerri)) = 0;
        vi = sqrt(vxi.^2+vyi.^2);
        erri = sqrt((vxi./vi.*xerri).^2+(vyi./vi.*yerri).^2);% propagation of error for velocity v = sqrt(vx^2+vy^2)
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
        xerrtmp = sqrt((xerri.*wi).^2+(xerrm.*wm).^2)./(wi + wm + eps);
        yerrtmp = sqrt((yerri.*wi).^2+(yerrm.*wm).^2)./(wi + wm + eps);
        
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
        fprintf("   > added %s data points from ITS_LIVE %s.\n",string(numel(find(~indi & indm & nanind))),years(ii));
        v_source(indi & ~indm & nanind) = source_count+1; % Measures v2
        fprintf("   > added %s data points from Measures Amundsen %s.\n",string(numel(find(indi & ~indm & nanind))),years(ii));
        v_source(~indi & ~indm & nanind) = source_count+2; % Blend of ITS_LIVE and MEasures 
        fprintf("   > added %s data points from ITS_LIVE + Measures Amundsen %s.\n",string(numel(find(~indi & ~indm & nanind))),years(ii));
        v_source_string(source_count) =  string(source_count)+": ITS_LIVE Antarctica annual "+years(ii);
        v_source_string(source_count+1) =  string(source_count+1)+": MeaSUREs Amundsen annual 2000-2001"+years(ii);
        v_source_string(source_count+2) =  string(source_count+2)+": Blend of ITS_LIVE and MeaSUREs Amundsen annual "+years(ii);
        source_count = source_count+3;

        nanind = isnan(vx);

    end     

else

%% STEP2 ALTERNATIVE
    
    error('TO DO');

end

%% STEP3
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
vtmp = sqrt(vxtmp.^2 + vytmp.^2);
errtmp = sqrt((vxtmp./vtmp.*xerrtmp).^2+(vytmp./vtmp.*yerrtmp).^2);
Fvx = griddedInterpolant(Xtmp,Ytmp,vxtmp,'cubic','none');
Fvy = griddedInterpolant(Xtmp,Ytmp,vytmp,'cubic','none');
Fxerr = griddedInterpolant(Xtmp,Ytmp,xerrtmp,'cubic','none');
Fyerr = griddedInterpolant(Xtmp,Ytmp,yerrtmp,'cubic','none');
            
vxm = Fvx(X,Y);
vym = Fvy(X,Y);
xerrm = Fxerr(X,Y); xerrm(isnan(xerrm))=0; % replace nans by 0 for calculation of total error 
yerrm = Fyerr(X,Y); yerrm(isnan(yerrm))=0; % replace nans by 0 for calculation of total error 

vm = sqrt(vxm.^2+vym.^2);
errm = sqrt((vxm./vm.*xerrm).^2+(vym./vm.*yerrm).^2);% propagation of error for velocity v = sqrt(vx^2+vy^2)
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
vi = sqrt(vxi.^2+vyi.^2);
erri = sqrt((vxi./vi.*xerri).^2+(vyi./vi.*yerri).^2);% propagation of error for velocity v = sqrt(vx^2+vy^2)
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
xerrtmp = sqrt((xerri.*wi).^2+(xerrm.*wm).^2)./(wi + wm + eps);
yerrtmp = sqrt((yerri.*wi).^2+(yerrm.*wm).^2)./(wi + wm + eps);

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
vxtmp = double(ncread(ncfile,'vx')); vxtmp = flipdim(vxtmp,2);
vytmp = double(ncread(ncfile,'vy')); vytmp = flipdim(vytmp,2);
vtmp = sqrt(vxtmp.^2 + vytmp.^2);
Fvx = griddedInterpolant(Xtmp,Ytmp,vxtmp,'cubic','none');
Fvy = griddedInterpolant(Xtmp,Ytmp,vytmp,'cubic','none');
Ferr = griddedInterpolant(Xtmp,Ytmp,errtmp,'cubic','none');
            
vxm = Fvx(X,Y);
vym = Fvy(X,Y);
xerrm = Ferr(X,Y); xerrm(isnan(xerrm))=0; % replace nans by 0 for calculation of total error 
yerrm = xerrm; 

vm = sqrt(vxm.^2+vym.^2);
errm = sqrt((vxm./vm.*xerrm).^2+(vym./vm.*yerrm).^2);% propagation of error for velocity v = sqrt(vx^2+vy^2)
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

v_source(indm & nanind) = source_count; % ITS_LIVE
fprintf("   > added %s data points from MeaSURES Central Antarctica 1997.\n",string(numel(find(~indi & indm & nanind))));
v_source_string(source_count) =  string(source_count)+": MeaSURES Central Antarctica 1997";
source_count = source_count+1;

nanind = isnan(vx);

if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');
    
    R = makerefmat(x(1),y(1),x(2)-x(1),y(2)-y(1));
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_vx.tif",vx',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_vy.tif",vy',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_xerr.tif",xerr',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_yerr.tif",yerr',R,'CoordRefSysCode','EPSG:3031');
    v = sqrt(vx.^2 + vy.^2); err = sqrt((vx./v.*xerr).^2+(vy./v.*yerr).^2);
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_err.tif",err',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/AntarcticVelocity_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_source.tif",v_source',R,'CoordRefSysCode','EPSG:3031');

    fprintf('Done.\n');

end

fprintf(' Creating gridded interpolants...');

y = flipdim(y,1);
[X,Y] = ndgrid(x,y);

Fus = griddedInterpolant(X,Y,flipdim(vx,2),'linear');
Fvs = griddedInterpolant(X,Y,flipdim(vy,2),'linear');
Fxerr = griddedInterpolant(X,Y,flipdim(xerr,2),'linear');
Fyerr = griddedInterpolant(X,Y,flipdim(yerr,2),'linear');
Fsource = griddedInterpolant(X,Y,flipdim(single(v_source),2),'nearest');

fprintf('done.\n')

fprintf(' Saving VelocityInterpolants with Fus, Fvs, Fxerr, Fyerr. \n')
save('GriddedInterpolants_2000-2001-2002-2003-1996_MeaSUREs_ITSLIVE_velocities','Fus','Fvs','Fxerr','Fsource','Fyerr','-v7.3');

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
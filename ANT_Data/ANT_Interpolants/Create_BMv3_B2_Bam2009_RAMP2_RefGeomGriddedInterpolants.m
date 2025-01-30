function Create_BMv3_B2_Bam2009_RAMP2_RefGeomGriddedInterpolants(CreateGeotiff)

%%  
% Create a reference geometry for 2015, based on Bedmachine v3. As we want
% to use this product as a reference for earlier ice-sheet geometries, we
% will extend the ice front beyond its current extent. By default, the 
% Bedmachine surface and density will be extruded alogn flow lines, beyond 
% the current ice front. Alternatively, one can include data from Bedmap2, 
% Bamber 2009 and RAMP2 in areas that were previously ice covered, e.g. the
% Larsen B ice Shelf, but the quality of those products near the edges of 
% these ice shelves is generally rather poor.

ExtrudeBedMachine = 0;
UseBedmap2 = 1;
UseBamber2004 = 1;
UseRAMP2 = 1;

if nargin==0 
    CreateGeotiff = 1;
end

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

%%% GEOMETRY
%% Read/process Bedmachine data
ncfile = froot_data+"/BedMachine_Antarctica/BedMachineAntarctica-v3.nc";
fprintf(" Reading BedMachine data from file %s...",ncfile);
ncfile_pigdig = froot_data+"/BedMachine_Antarctica/BedMachineAntarctica-v3-PIGdig_PHolland.nc";

x_bm = double(ncread(ncfile,'x'));
y_bm = double(ncread(ncfile,'y'));
x_bm = unique(x_bm);
y_bm = unique(y_bm);
[X_bm,Y_bm] = ndgrid(x_bm,y_bm);

%bed  = double(flipdim(ncread(ncfile,'bed'),2));
% use the modified bed from P. Holland. This is identical to BMv3, except
% for some digging downstream of the PIG grounding line
bed = double(flipdim(ncread(ncfile_pigdig,'bed_pigdig'),2));
surface = double(flipdim(ncread(ncfile,'surface'),2));
thickness = double(flipdim(ncread(ncfile,'thickness'),2));
firn  = double(flipdim(ncread(ncfile,'firn'),2));
mask  = double(flipdim(ncread(ncfile,'mask'),2));
% geoid = double(flipdim(ncread(ncfile,'geoid'),2));
% err =  double(flipdim(ncread(ncfile,'errbed'),2));
% source_bm =  double(flipdim(ncread(ncfile,'source'),2));
% R = makerefmat(x_bm(1),y_bm(1),x_bm(2)-x_bm(1),y_bm(2)-y_bm(1));
% geotiffwrite("./GeoTiffFiles/Bedmachine_v3_error",err',R,'CoordRefSysCode','EPSG:3031');
% geotiffwrite("./GeoTiffFiles/Bedmachine_v3_source",source_bm',R,'CoordRefSysCode','EPSG:3031');
fprintf("Done.\n")

% All heights are referenced to mean sea level using the geod EIGEN-6D4
% To obtain height with respect to WGS84 add geoid:
% Rema = Surface + Firn + Geoid
% Surface=bed+thickness over the grounded areas
% However this surface is not the REMA surface
% So the actual surface (with respect to sea level) is: s = surface+firn
%                                                       b=bed=surface-thickness=s-firn-thickness
%                                           therefore   h=s-b=surface+firn-(surface-thickness)=firn+thickness
%
%  ocean/ice/land mask
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice 
%  4 = Lake Vostok
%
% I'm assuming here that the user wants the variable s to be surface+firn, in which case the vertically averaged ice density must be
% modfied accordingly.
% 
s_bm=surface+firn ;
b_bm=surface-thickness ;
h_bm=firn+thickness ; % or just h=s-b 
bed_bm=bed ; 
rhoi= 917 ; 
rhow = 1027;
rhoMin = 100; 
%
% The 'firn' field in the Bedmachine data is not the firn thickness, but
% `firn air content'
%
% $$ \mathrm{firn} =  \frac{1}{\rho_i}  \int_{firn layer} (\rho_i - \rho_f) \, dz $$
%
% The firn air content has units of distance.
%
% \rho = (1-firn/h) \rho_i   where h is the total ice thickness (Ua definition).
%
%
rho=(1-firn./(h_bm+eps)).*rhoi ;  % vertically averaged ice density over the deph h=s-b=firn+thicknes
% Note: one would expect that the minimum value of rho would be equal to the firn density,
%       but there are a number of values where rho=0, that is when firn=h. Not sure how this can happen
%       when the ice thickness is larger than zero... Possibly some inconsistencies between the data sets providing air content (firn) and
%       ice thickness. Note: firn is 'air content'
I=h_bm>0 ; %figure(190) ; histogram(rho(I),'Normalization','probability') ; title('histogram of densities where h>0')
fprintf('Over glaciated areas %f%% of densities are smaller than %f kg/m^3. These densities are set to %f kg/m^3. \n',...
    100*numel(find(rho(I)<rhoMin))/numel(find(I)),rhoMin,rhoMin);
I=h_bm>0 & rho<rhoMin ; rho(I)=rhoMin ; 
% run some consistency checks
b_bm(mask==0)=0; s_bm(mask==0)=0; % ocean has zero ice draft and surface height
b_bm(mask==1)=bed_bm(mask==1); s_bm(mask==1)=bed_bm(mask==1); % ice-free land has ice draft and surface height equal to the bed height
b_bm(mask==2)=bed_bm(mask==2); % grounded ice has ice draft equal to bed height
bed_bm(mask==4)=b_bm(mask==4); % get rid of Lake Vostok by setting lower ice surface (b) equal to bedrock elevation (B) over the lake as defined my mask.
mask_bm = mask;

grounded = ismember(mask_bm,[1 2 4]); 

% Start with BedMachine as the basis: 
source = zeros(size(h_bm),'uint8'); 
source(h_bm>0 | grounded) = 1; 

if CreateGeotiff
    
    fprintf('Writing origial Bedmachine GeoTiff files...');
    
    R = maprefcells([x_bm(1) x_bm(end)],[y_bm(1) y_bm(end)],[numel(x_bm),numel(y_bm)]);
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_bed",bed_bm',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_mask",mask_bm',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_icethickness",h_bm',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_watercolumnthickness",bed_bm'-b_bm',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_surface",s_bm',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_draft",b_bm',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_density",rho',R,'CoordRefSysCode','EPSG:3031');
     
    fprintf('Done.\n');

end

% Distance in kilometers to BedMachine ocean data: 
D_ocean = bwdist(mask_bm==0); % euclidean distance to the nearest non-zero pixel
% narrow strips of floating ice with thin water column (<15m) exist along
% the ice front. Remove these by looking in a 1.5km strip along the open ocean
% for cells with mask_bm==3. Then set these cells to open ocean
D_grounded = bwdist(mask_bm==2);
I = D_ocean<1.5 & D_grounded<1.5 & mask_bm==3;% & b_bm-bed_bm<15;
b_bm(I) = 0; % ocean
mask_bm(I) = 0;
rho(I) = rhoi;
s_bm(I) = 0;
h_bm(I) = 0;

source(I) = 2; 

s_b2 = 0*X_bm+nan; s_bam = s_b2; s_ramp = s_b2;
if UseBedmap2
    %% Load Bedmap2 surface
    fprintf('Adding Bedmap2 data...');
    addpath(froot_data+"/Bedmap2/bedmap2_tif");
    s_b2 = bedmap2_interp(X_bm,Y_bm,'surface');
    isn = ~(s_b2>0);  
    s_b2(isn) = bedmap2_interp(X_bm(isn),Y_bm(isn),'surface','nearest'); 
    s_b2(bwdist(~isfinite(s_b2))*.5<1) = NaN; % Eliminates 1 km perimeter to avoid narrow strip of bedmap2 data along edges
    fprintf('Done. \n');
else
    s_b2 = nan*X_bm;
end

if UseBamber2004
    %% Load Bamber2004 surface
    fprintf('Adding Bamber2004 data...');
    addpath(froot_data+"/Bamber_SurfDEM_2004");
    [~,~,s_bam] = read_Bamber_2004(X_bm,Y_bm,1);
    s_bam(bwdist(~isfinite(s_bam))*.5<5) = NaN; % Eliminates 10 km perimeter b/c drooping ice sheet edges in Bamber and RAMP dems
    fprintf('Done. \n');
else
    s_bam = nan*X_bm;
end

if UseRAMP2
    %% Load RAMP2 surface
    fprintf('Adding RAMP2 data...');
    [Iramp,xramp,yramp] = geoimread(froot_data+"/RAMP2/osu91a200m.tif"); % The OSU version of the file is referenced to the geoid whereas the RAM2_DEM.tif is referenced to wgs84. Use OSU.  
    Iramp = double(Iramp); 
    Iramp(Iramp==0) = NaN; % Convert zeros to NaN so interpolation won't produce thin ice shelf edges.  
    s_ramp = interp2(xramp,yramp,Iramp,X_bm,Y_bm); 
    s_ramp(bwdist(~isfinite(s_ramp))*.5<5) = NaN;  % Eliminates 10 km perimeter b/c drooping ice sheet edges in Bamber and RAMP dems
    fprintf('Done. \n');
else
    s_ramp = nan*X_bm;
end

if UseBedmap2 || UseBamber2004 || UseRAMP2
    %% convert surface elevation to thickness, assuming floatation
    % first do some adjustments to rho: remove all rho=917 for ocean cells, and
    % interpolate linearly from grounded/floating ice.
    Frho = scatteredInterpolant(X_bm(mask_bm~=0),Y_bm(mask_bm~=0),rho(mask_bm~=0),'natural');
    rho = Frho(X_bm,Y_bm);
    hydro = rhow./(rhow-rho);
    h_b2 = s_b2.*hydro;
    h_bam = s_bam.*hydro;
    h_ramp = s_ramp.*hydro;
    
    %% Prevent new grounding:
    h_max = -rhow./rho.*bed_bm; % maximum thickness that could exist at neutral buoyancy (no new grounding).
    h_b2(h_b2>h_max) = h_max(h_b2>h_max); 
    h_bam(h_bam>h_max) = h_max(h_bam>h_max); 
    h_ramp(h_ramp>h_max) = h_max(h_ramp>h_max); 
    
    % Eliminate negative thickness: 
    h_b2(h_b2<0) = nan; 
    h_bam(h_bam<0) = nan; 
    h_ramp(h_ramp<0) = nan;
    
    %% Combine surface heights 
    % Distance in kilometers to BedMachine ocean data: 
    D = bwdist(mask_bm==0)*0.5; % euclidean distance to the nearest non-zero pixel
    D_thresh_thinice = 10; % replace thin ice within 10 km of ocean (for example, where Mertz Glacier tongue broke off and is now thin on its end in the BedMachine data, overwrite that thin ice with the full thickness from Bedmap2.) 
    
    h_bm_adj = h_bm;
    % Start overwriting BedMachine: 
    source(h_bm_adj<(h_b2/2) & D<D_thresh_thinice & ~grounded) = 3; 
    h_bm_adj(source==3) = h_b2(source==3); 
    
    source(h_bm_adj<(h_bam/2) & D<D_thresh_thinice & ~grounded) = 4; 
    h_bm_adj(source==4) = h_bam(source==4); 
    
    I = find(h_bm_adj<(h_ramp/2) & D<D_thresh_thinice & ~grounded);
    % remove ramp2 data from FRISS region due to odd regrounding
    J = find(Y_bm(I)>840441 & Y_bm(I)<934255 & X_bm(I)>-1182106 & X_bm(I)<-970374); I(J)=[];
    source(I) = 5;
    h_bm_adj(source==5) = h_ramp(source==5); 
    
    h_bm_adj(h_bm_adj==0 & ~grounded) = nan;
    
    L = bwlabel(isnan(h_bm_adj)); 
    h_bm_adj = regionfill(h_bm_adj,L>1);
    source(L>1) = 6; % interpolated
end

if ExtrudeBedMachine

    % We will extrude floating ice along flowlines. To know the ice
    % thickness, we need to extrude surface elevation and ice density.

    % Temporarily create and save interpolants
    Fs = griddedInterpolant(X_bm,Y_bm,s_bm,'cubic');
    Fb = Fs; Fb.Values = b_bm; 
    FB = Fs; FB.Values = bed_bm;
    Frho = Fs; Frho.Values = rho;
    Fmask = Fs; Fmask.Values = mask_bm;

    Velinterpolantfile = "./GriddedInterpolants_2015-2016_MeaSUREs_ITSLIVE_Velocities.mat";
    Geominterpolantfile =  "BMv3_tmp.mat"; save(Geominterpolantfile,"Fs","Fb","FB","Frho","Fmask","-v7.3");
    ScalarInterpolant = [];
    fields_to_extrude = "-geom-";
    Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,ScalarInterpolant,0,fields_to_extrude);

    load(Geominterpolantfile+"_EXTRUDED.mat","Fs","Fb","Frho");
    s_bm = Fs.Values;
    rho_bm = Frho.Values;

    bfloat = s_bm.*(rho_bm./(rho_bm-rhow));
    b_bm = max(bfloat,bed_bm);
    h_bm_adj = s_bm - b_bm;

end

% Calculate s and b from h and B
S = 0*X_bm;
hf=rhow.*(0*S-bed_bm)./rho;
GF = HeavisideApprox(100,h_bm-hf,1); % 1 if grounded, 0 if afloat;
bfloat=S-rho.*h_bm_adj./rhow;
b_bm_adj=GF.*bed_bm + (1-GF) .* bfloat ;

% Enforce b Above B
% because the grounding line is `smeared out' a bit for a finite CtrlVar.kH one can
% have situations where b<B. For CtrlVar.kH>0.1 this is not really much of an issue
% and anyhow this is a direct consequence of using a smooth step function (explained in detail in UaCompendium)
I=b_bm_adj<bed_bm ; b_bm_adj(I)=bed_bm(I);
s_bm_adj=b_bm_adj+h_bm_adj;

% To avoid grounding along ice front margins, impose that any ice seaward of 
% BedmapMachine iceshelf mask is floating
I = find(b_bm_adj<=bed_bm & mask_bm==0);
b_bm_adj(I) = bed_bm(I) + 10;
% recalculate s and h
s_bm_adj(I) = b_bm_adj(I).*(rho(I)-rhow)./rho(I);
h_bm_adj(I) = s_bm_adj(I) - b_bm_adj(I);

mask_bm_adj = 0*X_bm;
%  ocean/ice/land mask
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice 
%  4 = Lake Vostok
mask_bm_adj(mask_bm==1)=1;
mask_bm_adj(mask_bm==4)=4;
mask_bm_adj(b_bm_adj==bed_bm)=2;
mask_bm_adj(b_bm_adj>bed_bm)=3;

if CreateGeotiff
    
    fprintf('Writing GeoTiff files for new geometry...');
    
    R = maprefcells([x_bm(1) x_bm(end)],[y_bm(1) y_bm(end)],[numel(x_bm),numel(y_bm)]);
    geotiffwrite("./GeoTiffFiles/Bedmachinev3_Bedmap2_Bamber2009_RAMP2_mask",mask_bm_adj',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachinev3_Bedmap2_Bamber2009_RAMP2_icethickness",h_bm_adj',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachinev3_Bedmap2_Bamber2009_RAMP2_watercolumnthickness",bed_bm'-b_bm_adj',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachinev3_Bedmap2_Bamber2009_RAMP2_surface",s_bm_adj',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachinev3_Bedmap2_Bamber2009_RAMP2_draft",b_bm_adj',R,'CoordRefSysCode','EPSG:3031');      
    geotiffwrite("./GeoTiffFiles/Bedmap2_surface",s_b2',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bamber2004_surface",s_bam',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/RAMP2_surface",s_ramp',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachinev3_Bedmap2_Bamber2009_RAMP2_source",source',R,'CoordRefSysCode','EPSG:3031');
   
    %geotiffwrite("./GeoTiffFiles/Bedmachine_v3_s",flipdim(s',2)',R,'CoordRefSysCode','EPSG:3031');
    %geotiffwrite("./GeoTiffFiles/Bedmachine_v3_b",flipdim(b',2)',R,'CoordRefSysCode','EPSG:3031');
    %geotiffwrite("./GeoTiffFiles/Bedmachine_v3_B",flipdim(B',2)',R,'CoordRefSysCode','EPSG:3031');
    %geotiffwrite("./GeoTiffFiles/Bedmachine_v3_rho_adjusted",rho',R,'CoordRefSysCode','EPSG:3031');
    
    fprintf('Done.\n');

end

%% create gridded interpolants
fprintf(' Creating gridded interpolants for unmodified Bedmachine...')

FB = griddedInterpolant(X_bm,Y_bm,bed_bm,'cubic');
Fs = griddedInterpolant(X_bm,Y_bm,s_bm_adj,'cubic');
Fb = griddedInterpolant(X_bm,Y_bm,b_bm_adj,'cubic');
Frho = griddedInterpolant(X_bm,Y_bm,rho,'cubic'); 
Fmask = griddedInterpolant(X_bm,Y_bm,mask_bm_adj,'nearest');
Fsource = griddedInterpolant(X_bm,Y_bm,single(source),'nearest');

fprintf('done.\n')

fprintf(' Saving BedMachineGriddedInterpolants with Fs, Fb, FB, Frho, Fsource and Fmask. \n');

save('GriddedInterpolants_Bedmachinev3_Bedmap2_Bamber2009_RAMP2_RefGeometry','Fs','Fb','FB','Frho','Fmask','Fsource','-v7.3');

fprintf('done.\n')

return


%% Be carefull here! Bedmachine-v3 has floating ice in thin strips along
%% most of the Antarctic coastline. This creates many tiny ice shelves when
%% interpolated to Ua. I modify the Bedmachine-v3 mask and geometry fields
%% to remove these thin strips. The criteria are a bit ad-hoc but seem to do
%% a reasonable job: 
% 1. look for cells (i,j) with mask = 3 (floating ice)
% 2. check 8 neigbouring cells (i-1:i+1,j-1:j+1). If at least 1 is open
% ocean and 2 are grounded, then modify cell (i,j) to be open ocean
% 3. fill cell (i,j) with s=0, b=0

kk=1; I_flag = []; mask_mod = mask;
while kk<5
    %% STEP1.
    I = find(mask_mod == 3);
    %% STEP2.
    Mim1j = [zeros(1,size(mask_mod,2)); mask_mod(1:end-1,:)];
    Mim1jp1 =  [zeros(1,size(mask_mod,2)); [mask_mod(1:end-1,2:end) zeros(size(mask_mod,1)-1,1)]];
    Mim1jm1 = [zeros(1,size(mask_mod,2)); [zeros(size(mask_mod,1)-1,1) mask_mod(1:end-1,1:end-1)]];
    Mijp1 = [mask_mod(1:end,2:end) zeros(size(mask_mod,1),1)];
    Mijm1 = [zeros(size(mask_mod,1),1) mask_mod(1:end,1:end-1)];
    Mip1j = [mask_mod(2:end,1:end); zeros(1,size(mask_mod,2))];
    Mip1jp1 = [[mask_mod(2:end,2:end); zeros(1,size(mask_mod,2)-1)] zeros(size(mask_mod,1),1)];
    Mip1jm1 = [[mask_mod(2:end,1:end-1) zeros(size(mask_mod,1)-1,1)]; zeros(1,size(mask_mod,2))];
    Neighbours = [Mim1j(I) Mim1jp1(I) Mim1jm1(I) Mijp1(I) Mijm1(I) ...
        Mip1j(I) Mip1jp1(I) Mip1jm1(I)];
    Neighbours_OpenOcean = 0*Neighbours; Neighbours_OpenOcean(Neighbours==0)=1;
    Neighbours_Grounded = 0*Neighbours; Neighbours_Grounded(Neighbours==2)=1;
    Neighbours_Floating = 0*Neighbours; Neighbours_Floating(Neighbours==3)=1;
    if kk<3     
        J_OpenOcean = sum(Neighbours_OpenOcean,2)>1;
        J_Grounded = sum(Neighbours_Grounded,2)>1;
        Itmp = I(J_OpenOcean(:) & J_Grounded(:));
    else
        J_OpenOcean = sum(Neighbours_OpenOcean,2)>0;
        %J_Grounded = sum(Neighbours_Grounded,2)>0;
        J_Floating = sum(Neighbours_Floating,2)<2;
        Itmp = I(J_OpenOcean(:) & J_Floating(:));
    end   
    
    I_flag = [I_flag; Itmp(:)];  
    mask_mod(Itmp) = 0;

    kk=kk+1;
end

%% STEP3.
I = unique(I_flag);
b_mod = b; s_mod = s; rho_mod = rho;
b_mod(I) = 0;
s_mod(I) = 0;
rho_mod(I) = rhoi; % rhoi over open ocean

%% create gridded interpolants
fprintf(' Creating gridded interpolants for modified Bedmachine...')

FB = griddedInterpolant(Xdata,Ydata,flipdim(B',2),'cubic');
Fs = griddedInterpolant(Xdata,Ydata,flipdim(s_mod',2),'cubic');
Fb = griddedInterpolant(Xdata,Ydata,flipdim(b_mod',2),'cubic');
Frho = griddedInterpolant(Xdata,Ydata,flipdim(rho_mod',2),'cubic'); 
Fmask = griddedInterpolant(Xdata,Ydata,flipdim(mask_mod',2),'nearest');

fprintf('done.\n')

fprintf(' Saving BedMachineGriddedInterpolants with Fs, Fb, FB, Frho and Fmask. \n');
save('GriddedInterpolants_Bedmachine_v3_FastIceRemoved','Fs','Fb','FB','Frho','Fmask','-v7.3');


%% Write Geotiff files
if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');
    
    R = makerefmat(xdata(1),ydata(1),xdata(2)-xdata(1),ydata(2)-ydata(1));
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_mask_FastIceRemoved",flipdim(mask_mod',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_s_FastIceRemoved",flipdim(s_mod',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_b_FastIceRemoved",flipdim(b_mod',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_B_FastIceRemoved",flipdim(B',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_rho_FastIceRemoved",flipdim(rho_mod',2)',R,'CoordRefSysCode','EPSG:3031');
    
    fprintf('Done.\n');

end

%% I generate scatteredInterpolants for the modified geometry that do not include 
%% open ocean gridpoints. This allows extrapolation of nearest data to areas 
%% where Ua might think there is ice, but Bedmachine indicates there is ocean
% fprintf(' Creating scattered interpolants for s and b...')
% 
% I = find(mask_mod~=0);
% s_mod = flipdim(s_mod',2);
% b_mod = flipdim(b_mod',2);
% B = flipdim(B',2);
% Fs = scatteredInterpolant(Xdata(I),Ydata(I),s_mod(I),'natural','nearest');
% Fb = Fs; Fb.Values = b_mod(I);
% FB = Fs; FB.Values = B(I);
% 
% fprintf('done.\n')
% 
% fprintf(' Saving BedMachineGriddedInterpolants with Fs, Fb\n');
% save('ScatteredInterpolants_Bedmachine_v3_FastIceRemoved_NoOpenOcean','Fs','Fb','FB','-v7.3');
% 
% fprintf('done.\n')
%

%% Rather than use scatteredInterpolants, which is very slow, I reduce the 
%% risk of potential extrapolation errors by filling s, b, h and rho 
%% with nearest neighbouring values in a 5km-wide buffer over the open ocean.
%% These extruded fields can then be stored as gridded interpolants
b_buffer = b_mod;
s_buffer = s_mod;
h_buffer = s_mod-b_mod;
rho_buffer = rho;
mask_buffer = mask_mod;

% 1. find nearest neighbour land or shelf ice cells to all ocean cells
% within a 10km buffer
[Ix_Ocean,Iy_Ocean] = find(mask_mod==0);
[Ix_notOcean,Iy_notOcean] = find(mask_mod~=0);
for ii=1:numel(Ix_Ocean)
    Ix_min = max(Ix_Ocean(ii)-10,1); % buffer is 10 cells wide (10 x 500m = 5km)
    Ix_max = min(Ix_Ocean(ii)+10,nx);
    Iy_min = max(Iy_Ocean(ii)-10,1); % buffer is 10 cells wide (10 x 500m = 5km)
    Iy_max = min(Iy_Ocean(ii)+10,ny);
    Ix = [Ix_min:Ix_max]; Iy = [Iy_min:Iy_max];
 %   if any(ismember(Ix,Ix_notOcean)) && any(ismember(Iy,Iy_notOcean))
        dist = sqrt((Xdata(Ix_Ocean(ii),Iy_Ocean(ii))-Xdata(Ix,Iy)).^2+...
            (Ydata(Ix_Ocean(ii),Iy_Ocean(ii))-Ydata(Ix,Iy)).^2);
        dist = dist.*(mask(Ix,Iy)~=0); dist(dist==0)=NaN;
        [M,Imin] = min(dist,[],'all','omitnan');
        if ~isnan(M)
            [Iminx,Iminy] = ind2sub(size(dist),Imin);
            b_buffer(Ix_Ocean(ii),Iy_Ocean(ii)) = b_buffer(Ix(Iminx),Iy(Iminy));
            s_buffer(Ix_Ocean(ii),Iy_Ocean(ii)) = s_buffer(Ix(Iminx),Iy(Iminy));
            h_buffer(Ix_Ocean(ii),Iy_Ocean(ii)) = h_buffer(Ix(Iminx),Iy(Iminy));
            rho_buffer(Ix_Ocean(ii),Iy_Ocean(ii)) = rho(Ix(Iminx),Iy(Iminy));
            mask_buffer(Ix_Ocean(ii),Iy_Ocean(ii)) = mask_buffer(Ix(Iminx),Iy(Iminy));
        end

 %   end
    if mod(ii,1e5)==0
        fprintf('Done %i out of %i \n',ii,numel(Ix_Ocean));
    end
end 

%% create gridded interpolants
fprintf(' Creating gridded interpolants with buffer...')

FB = griddedInterpolant(Xdata,Ydata,flipdim(B',2),'cubic');
Fs = griddedInterpolant(Xdata,Ydata,flipdim(s_buffer',2),'cubic');
Fb = griddedInterpolant(Xdata,Ydata,flipdim(b_buffer',2),'cubic');
Frho = griddedInterpolant(Xdata,Ydata,flipdim(rho_buffer',2),'cubic'); 
Fmask = griddedInterpolant(Xdata,Ydata,flipdim(mask_buffer',2),'cubic'); 

%Fmask = griddedInterpolant(Xdata,Ydata,flipdim(mask',2),'nearest');
fprintf('done.\n')

%% adjust rho to avoid noise near mountain ranges
% fprintf(' Cleaning up rho ...')
% rhotemp = Frho.Values;
% rhotemp(rhotemp==917)=NaN;
% I = find(isnan(rhotemp));
% Xtemp=Xdata; Ytemp=Ydata;
% Xtemp(I)=[]; Ytemp(I)=[]; rhotemp(I)=[];
% Frhotemp=scatteredInterpolant(Xtemp(:),Ytemp(:),rhotemp(:),'natural');
% Frho.Values=Frhotemp(Xdata,Ydata);
% Frho.Values(Fs.Values==0)=917;
% Frho.Values(Frho.Values<100)=100;
% fprintf('done.\n')

fprintf(' Saving BedMachineGriddedInterpolants with Fs, Fb, FB and Frho. \n');

save('GriddedInterpolants_Bedmachine_v3_FastIceRemoved_WithBuffer','Fs','Fb','FB','Frho','Fmask','-v7.3');

fprintf('Done.\n');

%% Write Geotiff files
if CreateGeotiff
    
    fprintf('Writing GeoTiff files \n');
    
    R = makerefmat(xdata(1),ydata(1),xdata(2)-xdata(1),ydata(2)-ydata(1));
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_FastIceRemoved_s_WithBuffer",flipdim(s_buffer',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_FastIceRemoved_b_WithBuffer",flipdim(b_buffer',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_FastIceRemoved_B_WithBuffer",flipdim(B',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_FastIceRemoved_rho_WithBuffer",flipdim(rho_buffer',2)',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/Bedmachine_v3_FastIceRemoved_mask_WithBuffer",flipdim(mask_buffer',2)',R,'CoordRefSysCode','EPSG:3031');
    
    fprintf('Done.\n');

end

return


%% Test gridded interpolants
% load in Sainan's ISMIP6 mesh of Antarctica and map on this mesh. 
% Just done for testing purposes and to see if the gridded data looks sensible.
fprintf(' Testing interpolants by mapping on a FE grid...')

load("./MeshAntISMIP6_2.mat"); 
xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 

sFE=Fs(xFE,yFE);
bFE=Fb(xFE,yFE);
BFE=FB(xFE,yFE);
rhoFE=Frho(xFE,yFE);

CtrlVar.PlotXYscale=1000;
 
bfig=FindOrCreateFigure('b') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,bFE);
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('b') ; title(cbar,'m a.s.l')

sfig=FindOrCreateFigure('s') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,sFE);
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('s') ; title(cbar,'m a.s.l')

Bfig=FindOrCreateFigure('B') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,BFE);
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('B') ; title(cbar,'m a.s.l')

rhofig=FindOrCreateFigure('rho') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,rhoFE); 
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('rho') ; title(cbar,'kg/m^3')

wctfig=FindOrCreateFigure('bFE-BFE') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,bFE-BFE); 
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('b-B') ; title(cbar,'m')

fprintf('done.\n')
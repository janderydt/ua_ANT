function Create_Geometry_from_reference_and_ds(targettime,CreateGeotiff)

% This script takes the reference geometry, which is a mosaic of Bedmachine
% v3, Bedmap2, Bamber2009 and RAMP2 data, and applies a surface correction
% from ITS-LIVE, CPOM and Paolo2023 to obtain a geometry that is
% representative for the targettime below.

if nargin==0
    targettime = datetime(2020,06,01); % format (YYYY,MM,DD)
    CreateGeotiff=1;
elseif nargin==1
    CreateGeotiff=1;
end

fprintf("Loading gridded interpolants for reference geometry...")

load('GriddedInterpolants_Bedmachinev3_Bedmap2_Bamber2009_RAMP2_RefGeometry');
FB_ref = FB; Fb_ref = Fb; Fmask_ref = Fmask; Frho_ref = Frho; Fs_ref = Fs; Fsource_ref = Fsource;

clear FB Fb Fmask Frho Fs Fsource;

fprintf("done.\n");

x = Fb_ref.GridVectors{1};
y = Fb_ref.GridVectors{2};
s = Fs_ref.Values;
B = FB_ref.Values;
rho = Frho_ref.Values;
mask = Fmask_ref.Values;

[X,Y] = ndgrid(x,y);

ds_file = "ds_Nilsson_Paolo_"+string(targettime)+".mat";

if ~exist(ds_file,"file")

    fprintf("Calculating surface correction.\n");

    source = Fsource_ref.Values;
    referencetime = datetime(2015,01,01)+0*X;
    %referencetime(source<3) = datetime(2015,01,01); % Bedmachine v3
    referencetime(source==3) = datetime(1995,01,01); % Bedmap2: Griggs and Bamber for the ice shelves
    referencetime(source==4) = datetime(2004,01,01); % Bamber2009
    referencetime(source==5) = datetime(1990,01,01); % RAMP2 best guess
    
    [ds,ds_err,ds_source] = Calc_ITSLIVE_CPOM_Paolo2023_Deltas(X,Y,referencetime,targettime);    

    if CreateGeotiff
    
        fprintf('Writing GeoTiff files \n');
        
        R = maprefcells([x(1) x(end)],[y(1) y(end)],[numel(x),numel(y)]);
        geotiffwrite("./GeoTiffFiles/ds_"+string(targettime)+".tif",ds',R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/ds_err_"+string(targettime)+".tif",ds_err',R,'CoordRefSysCode','EPSG:3031');
       
        fprintf('done.\n');

    end

    save(ds_file,"x","y","ds","ds_source","ds_err");

else
    
    fprintf("Loading existing surface correction for %s...",string(targettime));

    load(ds_file);

    fprintf('done.\n');

end

fprintf("Applying surface correction and calculating new draft and mask...");
% apply mask: THIS DOES NOT HELP TO REDUCE NOISE
%dsmask = load("ds_TREND_mask","mask");
%dsmask.mask(dsmask.mask==0)=nan;
%ds_tmp = dsmask.mask.*ds;
%clear dsmask;
% remove ds with errors > 2.5m
%ds_tmp(ds_err./ds>0.1)=nan;: THIS REMOVES TOO MANY DATA POINTS
ds(ds_err>2.5)=nan;

% replace nans with interpolated values if nearest non-nan is not more than 10km away
ind = ~isnan(ds);
Fds = scatteredInterpolant(X(ind),Y(ind),double(ds(ind)),'natural','none');
dx = X(2,1)-X(1,1);
ind = find((bwdist(~isnan(ds))*dx)<10e3);
ds(ind) = Fds(X(ind),Y(ind));
ds(isnan(ds) & mask~=0) = 0;
ds(mask==0)=nan;

s_new = s - ds; % subtract ds here because it is negative for thinning
rhow = 1027;

[b_new,h_new] = Calc_dh_from_ds(X,Y,s_new,B,rho,rhow);

% To avoid grounding along ice front margins, impose that any ice seaward of 
% BedmapMachine iceshelf mask is floating
ncfile = getenv("froot_data")+"/BedMachine_Antarctica/BedMachineAntarctica-v3.nc";
fprintf("Reading mask from BedMachine data...");
mask_bm  = double(flipdim(ncread(ncfile,'mask'),2));
I = find(b_new<=B & mask_bm==0);
b_new(I) = B(I) + 10;
% recalculate s and h
s_new(I) = b_new(I).*(rho(I)-rhow)./rho(I);
h_new(I) = s_new(I) - b_new(I);

%  ocean/ice/land mask
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice 
mask_new = 0*X;
mask_new(h_new==0 & B>0) = 1;
mask_new(b_new<=B) = 2;
mask_new(b_new>B) = 3;

fprintf('Done.\n')

%% Save interpolants
fprintf("Creating and saving gridded interpolants for modified geometry...");

% make sure all fields are doubles.
FB = griddedInterpolant(double(X),double(Y),double(B),'cubic');
Fb = FB; Fb.Values = double(b_new);
Fs = FB; Fs.Values = double(s_new);
Frho = FB; Frho.Values = double(rho);
Fmask = griddedInterpolant(double(X),double(Y),mask_new,'nearest');

save("GriddedInterpolants_Geometry_"+string(targettime),'Fs','Fb','FB','Frho','Fmask','-v7.3');

fprintf('done.\n')


%% Write Geotiff files
if CreateGeotiff
    
    fprintf("Writing GeoTiff files...");
    
    R = maprefcells([x(1) x(end)],[y(1) y(end)],[numel(x),numel(y)]);
    geotiffwrite("./GeoTiffFiles/ds_"+string(targettime)+".tif",ds',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/ds_err_"+string(targettime)+".tif",ds_err',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/surface_"+string(targettime)+".tif",s_new',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/draft_"+string(targettime)+".tif",b_new',R,'CoordRefSysCode','EPSG:3031');
    geotiffwrite("./GeoTiffFiles/mask_"+string(targettime)+".tif",mask_new',R,'CoordRefSysCode','EPSG:3031');

    fprintf("Done.\n");

end

clear Fds FB Fb Fs Frho Fmask Fsource ds B mask referencetime rho s source ds_err s_new b_new mask_new 

return




%% SOME PLOTTING
% Change in ice thickness
figure; hold on; pcolor(X,Y,deltah); shading flat; axis equal; caxis([-30 30]);
colormap(othercolor('RdYlBu7'));

% Water column thickness and Grounding lines
%figure; hold on; pcolor(X,Y,b_new-B); shading flat; axis equal; caxis([-1 1]);
S_insar = shaperead(froot_data+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
S_insar_mililloTW = shaperead(froot_data+"/AmundsenSea/Milillo_GL_Thwaites/Thwaites_2016-2017.shp");
S_insar_mililloDC = shaperead(froot_data+"/AmundsenSea/Milillo_GL_CrossonDotson/Grounding_Lines_CSK_Pope_Smith_Kohler.shp");
for ii=1:length(S_insar)
    date=S_insar(ii).DATE1;
    year_insar = date(1:4);
    [x_insar,y_insar] = ll2psxy(S_insar(ii).Y,S_insar(ii).X,-71,0);
    if str2double(year_insar)>=2011
        %ginsar=plot(x_insar,y_insar,"color",[255, 165, 0]/255,"linewidth",1); 
    else
        plot(x_insar,y_insar,"color",[0, 0, 0]/255,"linewidth",1); 
    end
end

% STD for Change in ice thickness
figure; hold on; pcolor(X,Y,deltah_err); shading flat; axis equal; caxis([0 10]);
colormap('parula');

% Water column thickness and Grounding lines
%figure; hold on; pcolor(X,Y,b_new-B); shading flat; axis equal; caxis([-1 1]);
S_insar = shaperead(froot_data+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
S_insar_mililloTW = shaperead(froot_data+"/AmundsenSea/Milillo_GL_Thwaites/Thwaites_2016-2017.shp");
S_insar_mililloDC = shaperead(froot_data+"/AmundsenSea/Milillo_GL_CrossonDotson/Grounding_Lines_CSK_Pope_Smith_Kohler.shp");
for ii=1:length(S_insar)
    date=S_insar(ii).DATE1;
    year_insar = date(1:4);
    [x_insar,y_insar] = ll2psxy(S_insar(ii).Y,S_insar(ii).X,-71,0);
    if str2double(year_insar)>=2011
        %ginsar=plot(x_insar,y_insar,"color",[255, 165, 0]/255,"linewidth",1); 
    else
        plot(x_insar,y_insar,"color",[0, 0, 0]/255,"linewidth",1); 
    end
end

% for ii=1:length(S_insar_mililloTW)
%     plot(S_insar_mililloTW(ii).X,S_insar_mililloTW(ii).Y,"color",[255, 165, 0]/255,"linewidth",2.5);  
% end
% 
% for ii=1:length(S_insar_mililloDC)
%     plot(S_insar_mililloDC(ii).X,S_insar_mililloDC(ii).Y,"color",[255, 165, 0]/255,"linewidth",2.5);                
% end
% 
end

function [b,h] = sB_to_bh(X,Y,s,B,rho,rhow)

    S = 0*X;
    
    % get a rough and a reasonable initial estimate for b if none is provided
    % The lower surface b is 
    %
    %
    %   b=max( B , (rhow S - rho s)/(rhow-rho) ) 
    %   where
    %
    %  h_f = rhow (S-B) / rho
    %
    %  b=s-h_f =

    hf=rhow*(S-B)./rho ;
    
    b0 =  max(B,(rho.*s-rhow.*S)./(rho-rhow)) ; 
    
    b=b0;
    h=s-b;
    
    % iteration
    ItMax=30 ; tol=100*eps ;  J=Inf ; I=0 ;
    JVector=zeros(ItMax,1)+NaN ;        
    
    CtrlVar.kH = 1;
    CtrlVar.Hh0 = 0;
    
    while I< ItMax && J > tol
        I=I+1;

        G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1
        dGdb=-DiracDelta(CtrlVar.kH,h-hf,CtrlVar.Hh0) ;
        
        F0=    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
        dFdb = 1 - dGdb.* (B -  (rho.*s-rhow.*S)./(rho-rhow)) ;
        
        db= -F0./dFdb ;
        
        b=b+db ;
        h=s-b ;
        
        F1 =    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
        
        JLast=J ;
        J=sum(F1.^2,'all','omitnan')/2 ;
        
        JVector(I)=J ;
        
    end

    fprintf("%s iterations needed for bh calculation.\n",string(I));

end

% return
% load('GriddedInterpolants_Bedmachine_v3_withbuffer.mat');
% load('MeshAntISMIP6_2_SainanSun.mat');
% B = FB(MUA.coordinates(:,1),MUA.coordinates(:,2));
% s1 = Fs(MUA.coordinates(:,1),MUA.coordinates(:,2));
% s2 = s1 + FDeltas(MUA.coordinates(:,1),MUA.coordinates(:,2));
% h1 = s1-Fb(MUA.coordinates(:,1),MUA.coordinates(:,2));
% rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
% rhow = 1027;
% [b2,h2,GF2]=Calc_bh_From_sBS([],MUA,s2,B,0*B,rho,rhow);
%  
% figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,h2-h1);
% return
% %% deleted some CPOM data on floating ice
% S_GL = shaperead('InSAR_GL_92-96_Combined.shp');
% I = find(inpoly([X_CPOM(:) Y_CPOM(:)],[S_GL(1).X(1:end-1)' S_GL(1).Y(1:end-1)']));
% X_CPOM(I)=[]; Y_CPOM(I)=[]; dh_CPOM(I)=[]; dh_err(I)=[];
% 
% % I = find(Y_CPOM(:)>-345e3 & Y_CPOM(:)<-310e3 & X_CPOM(:)>-1670e3 & X_CPOM(:)<-1660e3);
% % dhdt_CPOM(I) = NaN;
% % 
% % I = find(Y_CPOM(:)>-290e3 & Y_CPOM(:)<-280e3 & X_CPOM(:)>-1670e3 & X_CPOM(:)<-1655e3);
% % dhdt_CPOM(I) = NaN;
% % 
% % I = find(Y_CPOM(:)>-345e3 & Y_CPOM(:)<-315e3 & X_CPOM(:)>-1595e3 & X_CPOM(:)<-1565e3);
% % dhdt_CPOM(I) = NaN;
% 
% I = find(dh_err(:)>3 | isnan(dh_CPOM(:))); %dh_err(:)<-10 | dh_err(:)>8 | 
% X_CPOM(I)=[]; Y_CPOM(I)=[];
% dh_CPOM(I)=[]; dh_err(I)=[];
% 
% %% F. Paolo data for ice shelf
% if strfind(Paolo_version,'PIG_Paolo_trimmedv2')
%     dhdt_file = 'dhdt_PIG_Paolo_trimmedv2.shp';
% else
%     dhdt_file = ['dhdt_PIG_1996-2016_Paolo_',Paolo_version,'.shp'];
% end
% if ~exist(dhdt_file)
%     error([dhdt_file,' not found']);
% else
%     S = shaperead(dhdt_file);     
% end
% 
% X_Paolo = [S(:).X];
% Y_Paolo = [S(:).Y];
% dhdt_Paolo = [S(:).dhdt];
% 
% %% remove any data in excess of 15m/yr.
% %dhdt_Paolo(dhdt_Paolo<-15)=NaN;
% 
% dh_Paolo = dhdt_Paolo*20; % multiply average dhdt by number of years (1996-2016) to obtain dh
% I = find(isnan(dh_Paolo));
% X_Paolo(I)=[]; Y_Paolo(I)=[]; dh_Paolo(I)=[];
% 
% x = MUA.coordinates(:,1);
% y = MUA.coordinates(:,2);
% Fdh = scatteredInterpolant(double([X_CPOM(:);X_Paolo(:)]),double([Y_CPOM(:);Y_Paolo(:)]),double([dh_CPOM(:);dh_Paolo(:)]),'natural');
% dh = Fdh(x,y);
% 
% N=1;
% 
% for ii=1:numel(fields)
%     
%     datafield = fields(ii)  ;
%     
%     geotiff_file = froot_data+"/ICESat1_ICESat2_mass_change_updated_2_2021/"+extractBetween(datafield,1,4)+"/ais_"+datafield+".tif";
%     
%     fprintf(" Reading %s data from file %s \n",datafield,geotiff_file);
%     
%     [geotiff_data,~,refmap,bbox]=geotiffread(geotiff_file);
%     dvardt=squeeze(geotiff_data(1:N:end,1:N:end,1));
%     dvardt=double(dvardt'); dvardt=flipdim(dvardt,2);
%     x=[bbox(1,1)+refmap(2,1)/2:refmap(2,1):bbox(2,1)-refmap(2,1)/2];
%     y=[bbox(2,2)+refmap(1,2)/2:refmap(1,2):bbox(1,2)-refmap(1,2)/2];
%     y=flipdim(y,2);
%     [X,Y] = ndgrid(x,y);
% 
%     fprintf("Done. \n");
% 
%     geotiff_file = froot_data+"/ICESat1_ICESat2_mass_change_updated_2_2021/"+extractBetween(datafield,1,4)+"/ais_"+datafield+"_rmse.tif";
%     
%     fprintf(" Reading %s data from file %s \n",datafield,geotiff_file);
% 
%     [geotiff_data,~,refmap,bbox]=geotiffread(geotiff_file);
%     dvardt_err=squeeze(geotiff_data(1:N:end,1:N:end,1));
%     dvardt_err=double(dvardt_err'); dvardt_err=flipdim(dvardt_err,2);
%     x_err=[bbox(1,1)+refmap(2,1)/2:refmap(2,1):bbox(2,1)-refmap(2,1)/2];
%     y_err=[bbox(2,2)+refmap(1,2)/2:refmap(1,2):bbox(1,2)-refmap(1,2)/2];
%     y_err=flipdim(y_err,2);
%     [X_err,Y_err] = ndgrid(x_err,y_err);
% 
%     fprintf("Done. \n");
% 
%     data.(datafield).X = X;
%     data.(datafield).Y = Y;
%     data.(datafield).data = dvardt;
%     data.(datafield).X_err = X_err;
%     data.(datafield).Y_err = Y_err;
%     data.(datafield).rmse = dvardt_err;
% 
% end
% 
% fprintf("Creating gridded Interpolants \n");
% %%
% for ii=[1 3]
% 
%    X = [data.(fields(ii)).X(:); data.(fields(ii+1)).X(:)];
%    Y = [data.(fields(ii)).Y(:); data.(fields(ii+1)).Y(:)];
%    data_full = [data.(fields(ii)).data(:); data.(fields(ii+1)).data(:)];
%    I = find(~isnan(data_full));
% 
%    X_err = [data.(fields(ii)).X_err(:); data.(fields(ii+1)).X_err(:)];
%    Y_err = [data.(fields(ii)).Y_err(:); data.(fields(ii+1)).Y_err(:)];
%    rmse_full = [data.(fields(ii)).rmse(:); data.(fields(ii+1)).rmse(:)];
%    J = find(~isnan(rmse_full));
% 
%    if ii==1
%         Fdmdt = scatteredInterpolant(X(I),Y(I),data_full(I),'linear');
%         Fdmdt_err = scatteredInterpolant(X_err(J),Y_err(J),rmse_full(J),'linear');
%    elseif ii==3
%         Fdhdt = scatteredInterpolant(X(I),Y(I),data_full(I),'linear');
%         Fdhdt_err = scatteredInterpolant(X_err(J),Y_err(J),rmse_full(J),'linear');
%    end
% 
% end
% 
% fprintf("Done \n");
% 
% fprintf(' Saving dm/dt and dh/dt Interpolants with Fdmdt, Fdmdt_err, Fdhdt, Fdhdt_err. \n');
% 
% save('scatteredInterpolants_dhdt_2003-2019_Smith-etal-2000','Fdmdt','Fdmdt_err','Fdhdt','Fdhdt_err','-v7.3');
% 
% fprintf("Done \n");
% 
% return
% %% Test gridded interpolants
% % load in Sainan's ISMIP6 mesh of Antarctica and map on this mesh. 
% % Just done for testing purposes and to see if the gridded data looks sensible.
% fprintf(' Testing interpolants by mapping on a FE grid...')
% 
% load("./MeshAntISMIP6_2.mat"); 
% xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 
% 
% dmdt_FE=Fdmdt(xFE,yFE);
% dmdt_err_FE=Fdmdt_err(xFE,yFE);
% dhdt_FE=Fdhdt(xFE,yFE);
% dhdt_err_FE=Fdhdt_err(xFE,yFE);
% 
% CtrlVar.PlotXYscale=1000;
%  
% % thinning rates in ice equivalent, 917 kg/m3
% dmdtfig=FindOrCreateFigure('dmdt') ; 
% [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dmdt_FE*16); 
% clim([-60 60]); colormap(othercolor('RdYlBu4'));
% xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('\Deltam between 2003 and 2019 from IceSat (Smith et al 2020)') ; title(cbar,'m');
% 
% % convert dm/dt to dh/dt by multiplying with a factor 917/rho_ua
% load('GriddedInterpolants_Bedmachine_v3.mat','Frho');
% rho_FE = Frho(xFE,yFE);
% dhdt_ua_fig=FindOrCreateFigure('dhdt_ua') ; 
% [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dmdt_FE*917./rho_FE*16);
% clim([-60 60]); colormap(othercolor('RdYlBu4'));
% xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('\Deltah between 2003 and 2019 for Ua') ; title(cbar,'m');
% 
% % thinning rates as measured by IceSat
% dhdtfig=FindOrCreateFigure('dhdt') ; 
% [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dhdt_FE*16);
% clim([-60 60]); colormap(othercolor('RdYlBu4'));
% xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('\Deltah from IceSat (Smith et al 2020)') ; title(cbar,'m');
% 
% % difference between both methods
% dhdtdifffig=FindOrCreateFigure('dhdt difference') ; 
% [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,(dhdt_FE-dmdt_FE*917./rho_FE)*16);
% clim([-60 60]); colormap(othercolor('RdYlBu4'));
% xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title({'difference between \Deltah from IceSat (not corrected for FAC) and';' \Deltam converted to \Delath using Ua ice densities'}) ; title(cbar,'m');
% 
% 
% fprintf('done.\n')





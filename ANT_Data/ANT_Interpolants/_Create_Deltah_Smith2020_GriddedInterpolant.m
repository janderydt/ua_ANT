function Create_Deltah_Smith2020_GriddedInterpolant(numberofyears)

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

fields = ["dmdt_grounded", "dmdt_floating",...
 "dhdt_grounded", "dhdt_floating"];

N=1;

for ii=1:numel(fields)
    
    datafield = fields(ii)  ;
    
    geotiff_file = froot_data+"/ICESat1_ICESat2_mass_change_updated_2_2021/"+extractBetween(datafield,1,4)+"/ais_"+datafield+".tif";
    
    fprintf(" Reading %s data from file %s \n",datafield,geotiff_file);
    
    [geotiff_data,~,refmap,bbox]=geotiffread(geotiff_file);
    dvardt=squeeze(geotiff_data(1:N:end,1:N:end,1));
    dvardt=double(dvardt'); dvardt=flipdim(dvardt,2);
    x=[bbox(1,1)+refmap(2,1)/2:refmap(2,1):bbox(2,1)-refmap(2,1)/2];
    y=[bbox(2,2)+refmap(1,2)/2:refmap(1,2):bbox(1,2)-refmap(1,2)/2];
    y=flipdim(y,2);
    [X,Y] = ndgrid(x,y);

    fprintf("Done. \n");

    geotiff_file = froot_data+"/ICESat1_ICESat2_mass_change_updated_2_2021/"+extractBetween(datafield,1,4)+"/ais_"+datafield+"_rmse.tif";
    
    fprintf(" Reading %s data from file %s \n",datafield,geotiff_file);

    [geotiff_data,~,refmap,bbox]=geotiffread(geotiff_file);
    dvardt_err=squeeze(geotiff_data(1:N:end,1:N:end,1));
    dvardt_err=double(dvardt_err'); dvardt_err=flipdim(dvardt_err,2);
    x_err=[bbox(1,1)+refmap(2,1)/2:refmap(2,1):bbox(2,1)-refmap(2,1)/2];
    y_err=[bbox(2,2)+refmap(1,2)/2:refmap(1,2):bbox(1,2)-refmap(1,2)/2];
    y_err=flipdim(y_err,2);
    [X_err,Y_err] = ndgrid(x_err,y_err);

    fprintf("Done. \n");

    data.(datafield).X = X;
    data.(datafield).Y = Y;
    data.(datafield).data = dvardt;
    data.(datafield).X_err = X_err;
    data.(datafield).Y_err = Y_err;
    data.(datafield).rmse = dvardt_err;

end

fprintf("Creating gridded Interpolants \n");
%%
% load rho form Bedmachine
load("GriddedInterpolants_Bedmachine_v3_withbuffer.mat","Frho");

for ii=[1 3]

   X = [data.(fields(ii)).X(:); data.(fields(ii+1)).X(:)];
   Y = [data.(fields(ii)).Y(:); data.(fields(ii+1)).Y(:)];
   data_full = [data.(fields(ii)).data(:); data.(fields(ii+1)).data(:)];
   I = find(~isnan(data_full));

   X_err = [data.(fields(ii)).X_err(:); data.(fields(ii+1)).X_err(:)];
   Y_err = [data.(fields(ii)).Y_err(:); data.(fields(ii+1)).Y_err(:)];
   rmse_full = [data.(fields(ii)).rmse(:); data.(fields(ii+1)).rmse(:)];
   J = find(~isnan(rmse_full));

   if ii==1
        load("GriddedInterpolants_Bedmachine_v3_withbuffer.mat","Frho");
        % mask out any areas with rho<800 to avoid noise in Deltah
        Frho.Values(Frho.Values<800 | isnan(Frho.Values))=800;
        rho = Frho(X,Y); rho_err = Frho(X_err,Y_err);
        Deltah_from_dmdt = numberofyears*data_full*917./rho;
        Deltah_from_dmdt_err =  numberofyears*rmse_full*917./rho_err;
        FDeltah_from_dmdt = scatteredInterpolant(X(I),Y(I),Deltah_from_dmdt(I),'linear');
        FDeltah_from_dmdt_err = scatteredInterpolant(X_err(J),Y_err(J),Deltah_from_dmdt_err(J),'linear');

   elseif ii==3
        Deltah_from_dhdt = numberofyears*data_full;
        Deltah_from_dhdt_err = numberofyears*rmse_full;
        FDeltah_from_dhdt = scatteredInterpolant(X(I),Y(I),Deltah_from_dhdt(I),'linear');
        FDeltah_from_dhdt_err = scatteredInterpolant(X_err(J),Y_err(J),Deltah_from_dhdt_err(J),'linear');
   end

end

fprintf("Done \n");

fprintf(' Saving dm/dt and dh/dt Interpolants with FDeltah_from_dmdt, FDeltah_from_dmdt_err, FDeltah_from_dhdt, FDeltah_from_dhdt_err. \n');

save('scatteredInterpolants_Deltah_2003-2019_Smith-etal-2000','FDeltah_from_dmdt', 'FDeltah_from_dmdt_err', 'FDeltah_from_dhdt', 'FDeltah_from_dhdt_err','-v7.3');

fprintf("Done \n");

return
%% Test gridded interpolants
% load in Sainan's ISMIP6 mesh of Antarctica and map on this mesh. 
% Just done for testing purposes and to see if the gridded data looks sensible.
fprintf(' Testing interpolants by mapping on a FE grid...')

load("./MeshAntISMIP6_2.mat"); 
xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 

dmdt_FE=Fdmdt(xFE,yFE);
dmdt_err_FE=Fdmdt_err(xFE,yFE);
dhdt_FE=Fdhdt(xFE,yFE);
dhdt_err_FE=Fdhdt_err(xFE,yFE);

CtrlVar.PlotXYscale=1000;
 
% thinning rates in ice equivalent, 917 kg/m3
dmdtfig=FindOrCreateFigure('dmdt') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dmdt_FE*16); 
clim([-60 60]); colormap(othercolor('RdYlBu4'));
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('\Deltam between 2003 and 2019 from IceSat (Smith et al 2020)') ; title(cbar,'m');

% convert dm/dt to dh/dt by multiplying with a factor 917/rho_ua
load('GriddedInterpolants_Bedmachine_v3.mat','Frho');
rho_FE = Frho(xFE,yFE);
dhdt_ua_fig=FindOrCreateFigure('dhdt_ua') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dmdt_FE*917./rho_FE*16);
clim([-60 60]); colormap(othercolor('RdYlBu4'));
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('\Deltah between 2003 and 2019 for Ua') ; title(cbar,'m');

% thinning rates as measured by IceSat
dhdtfig=FindOrCreateFigure('dhdt') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dhdt_FE*16);
clim([-60 60]); colormap(othercolor('RdYlBu4'));
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title('\Deltah from IceSat (Smith et al 2020)') ; title(cbar,'m');

% difference between both methods
dhdtdifffig=FindOrCreateFigure('dhdt difference') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,(dhdt_FE-dmdt_FE*917./rho_FE)*16);
clim([-60 60]); colormap(othercolor('RdYlBu4'));
xlabel('xps (km)' ) ; ylabel('yps (km)' ) ; title({'difference between \Deltah from IceSat (not corrected for FAC) and';' \Deltam converted to \Delath using Ua ice densities'}) ; title(cbar,'m');


fprintf('done.\n')



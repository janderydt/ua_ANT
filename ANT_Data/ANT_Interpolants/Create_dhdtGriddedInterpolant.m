function Create_dhdtGriddedInterpolant(targettime,CreateGeotiff)

% This script creates dhdt fields from CPOM and ITSLive data on the
% Bedmachine grid. The output is an average over 9 years. 
% It also produces an interpolant for dhdt of grounded ice, with nans 
% removed using inpaint_nans.

if nargin==0
    targettime = datetime(2020,06,01); % format (YYYY,MM,DD)
    CreateGeotiff=1;
elseif nargin==1
    CreateGeotiff=1;
end

fprintf("Loading Bedmachine grid...")

ncfile = getenv("froot_data")+"/BedMachine_Antarctica/BedMachineAntarctica-v3.nc";

x = unique(double(ncread(ncfile,'x')));
y = unique(double(ncread(ncfile,'y')));

fprintf("done.\n");

[X,Y] = ndgrid(x,y);

if ~exist("dhdt_"+string(targettime)+".mat","file")

    fprintf("Calculating surface elevation change.\n");

    [dhdt,dhdt_CPOM,dhdt_err,dhdt_err_CPOM,dhdt_source] = Calc_ITSLIVE_CPOM_Paolo2023_dhdt(X,Y,targettime);

    if CreateGeotiff
    
        fprintf('Writing GeoTiff files \n');
        
        R = maprefcells([x(1) x(end)],[y(1) y(end)],[numel(x),numel(y)]);
        geotiffwrite("./GeoTiffFiles/dhdt_"+string(targettime)+".tif",dhdt',R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/dhdt_err_"+string(targettime)+".tif",dhdt_err',R,'CoordRefSysCode','EPSG:3031');
       
        fprintf('done.\n');

    end

    save("dhdt_"+string(targettime)+".mat","x","y","dhdt","dhdt_CPOM","dhdt_err","dhdt_err_CPOM","dhdt_source");

else

    fprintf("Loading existing dhdt for %s...",string(targettime));
    
    load("dhdt_"+string(targettime)+".mat","dhdt_CPOM","dhdt_err_CPOM");
    
    fprintf('done.\n');

end

fprintf("Replacing nans...");
% inpaint nans
dhdt_CPOM_inpaint=inpaint_nans(dhdt_CPOM,4);
dhdt_err_CPOM_inpaint=inpaint_nans(dhdt_err_CPOM,4);
fprintf('done.\n');

fprintf("Creating and saving gridded interpolants...");
% create gridded interpolants
Fdhdt_CPOM = griddedInterpolant(X,Y,dhdt_CPOM_inpaint);
Fdhdt_err_CPOM = Fdhdt_CPOM; 
Fdhdt_err_CPOM.Values = dhdt_err_CPOM_inpaint;

% save
save("GriddedInterpolants_dhdt_"+string(targettime)+".mat","Fdhdt_CPOM","Fdhdt_err_CPOM");
fprintf('done.\n');
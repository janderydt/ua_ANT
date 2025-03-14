function Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,ScalarInterpolant,CreateGeotiff,fields_to_extrude)

persistent Fus Fvs Fxerr Fyerr Fsource Frho Fb Fs

%addpath("../ANT_HelperFunctions/");

% This function extrudes surface velocities, ice geometry (surface and draft), 
% densities or any other scalar quantity from the present-day ice edge of Antarctica along flowlines
% of ice flux (H*q). The processing chain is largely based on a script by
% C. Greene
% (https://github.com/chadagreene/ice-shelf-geometry/blob/main/code/flow_dem_extend.m).

if nargin==0
    Velinterpolantfile = "GriddedInterpolants_2009-2010_MeaSUREs_ITSLIVE_Velocities";
    Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2009";
    ScalarInterpolant = [];   
    CreateGeotiff = 1;
    fields_to_extrude = '-v-geom-';
end

if nargin==1
    error('Specify input Geominterpolantfile');
end

if nargin==2
    ScalarInterpolant = [];
    CreateGeotiff = 0;
    fields_to_extrude = '-v-geom-';
end

if nargin==3
    fields_to_extrude = '-scalar-';
    CreateGeotiff = 0;
end

if nargin==4
    fields_to_extrude = '-scalar-';
end

fprintf("Interpolate data onto correct grids...");

if contains(fields_to_extrude,'-v-')

    fprintf("Load velocity and geometry interpolants...");

    if isempty(Fus)
        load(Velinterpolantfile,"Fus","Fvs","Fxerr","Fyerr","Fsource");
        load(Geominterpolantfile,"Fb","Fs");
    end
    
    fprintf("done.\n");

    x_v = Fus.GridVectors{1};
    y_v = Fus.GridVectors{2};
    [X_v,Y_v] = meshgrid(x_v,y_v);
    
    vx_v = Fus.Values';
    vy_v = Fvs.Values';
    v_source_v = Fsource.Values';
    err_v = hypot(Fxerr.Values',Fyerr.Values');
    H_v = Fs(X_v,Y_v) - Fb(X_v,Y_v);

    Fus=[]; Fvs=[]; Fxerr=[]; Fyerr=[]; Fsource=[]; Fb=[]; Fs=[];

end

if contains(fields_to_extrude,'-geom-')

    fprintf("Load velocity and geometry interpolants...");

    if isempty(Fus)
        load(Velinterpolantfile,"Fus","Fvs","Fxerr","Fyerr");
        load(Geominterpolantfile,"Frho","Fb","Fs","Fmask","FB");
    end
    
    fprintf("done.\n");

    x_g = Fb.GridVectors{1};
    y_g = Fs.GridVectors{2};
    [X_g,Y_g] = meshgrid(x_g,y_g);

    s = Fs.Values'; b = Fb.Values';
    rho = Frho.Values';

    rho(s==0) = nan;
    s(s==0) = nan; b(b==0) = nan;
    
    H_g = s-b;

    % remove velocities with large errors || they will not be used to
    % extrude the geometry
    err_tmp = hypot(Fxerr.Values',Fyerr.Values');
    Fus.Values(err_tmp>15)=nan;
    Fvs.Values(err_tmp>15)=nan;
    vx_g = Fus(X_g,Y_g);
    vy_g = Fvs(X_g,Y_g);
    %v_source_g =  Fsource(X_g,Y_g);

end

if contains(fields_to_extrude,'-scalar-')

    [nx,ny] = size(ScalarInterpolant.Values);
    filetoread = "Fields_to_extrude_scalar_nx"+string(nx)+"_ny"+string(ny)+".mat";

    if ~exist(filetoread,"file")

        x_scal = ScalarInterpolant.GridVectors{1};
        y_scal = ScalarInterpolant.GridVectors{2};
        [X_scal,Y_scal] = meshgrid(x_scal,y_scal);    

        fprintf("Load geometry interpolants...");
    
        if isempty(Fs)
            load(Geominterpolantfile,"Fs");
        end

        s = Fs(X_scal,Y_scal);
        clear Fs

        if isempty(Fb)
            load(Geominterpolantfile,"Fb");
        end

        b = Fb(X_scal,Y_scal);
        clear Fb

        fprintf("done.\n");

        s(s==0) = nan; b(b==0) = nan;    
        H_scal = s-b;

        fprintf("Load velocity interpolants...");

        if isempty(Fxerr)
            load(Velinterpolantfile,"Fxerr","Fyerr");
        end

        % remove velocities with large errors || they will not be used to
        % extrude the scalar field
        err_tmp = hypot(Fxerr.Values',Fyerr.Values');

        clear Fxerr Fyerr

        if isempty(Fus)
            load(Velinterpolantfile,"Fus","Fvs");
        end

        Fus.Values(err_tmp>15)=nan;
        Fvs.Values(err_tmp>15)=nan;
        vx_scal = Fus(X_scal,Y_scal);
        vy_scal =  Fvs(X_scal,Y_scal);
        %v_source_g =  Fsource(X_g,Y_g);

        clear Fus Fvs
    
        save(filetoread,"x_scal","y_scal","H_scal","vx_scal","vy_scal","-v7.3");

    else

        load(filetoread);

    end

    [X_scal,Y_scal] = meshgrid(x_scal,y_scal);
    scal = ScalarInterpolant.Values';
    scal(scal==0) = nan;

end

clear X_g X_v Y_g Y_v

fprintf("done.\n");

%% -v-
if contains(fields_to_extrude,'-v-')
    
    fprintf("Processing velocity fields.\n");

    fprintf("  > Inpaint nan regions...");

    % Remove any velocities with high errors as we don't want these to be
    % used to interpolate/extrude the velocity field. The values will be 
    % added back at the end
    I_higherr = find(err_v>15);
    vx_v_orig = vx_v; vx_v(I_higherr) = nan;
    vy_v_orig = vy_v; vy_v(I_higherr) = nan;

    % Label the nan regions, and the ocean will be L=1. 
    L = bwlabel(isnan(vx_v)); 

    % Fill any nan regions that are not ocean, and assign label
    vx_v = regionfill(vx_v,L>1); 
    vy_v = regionfill(vy_v,L>1);
    v_source_v(L>1) = max(v_source_v(:))+1; % interpolated
    %v = hypot(vx_v,vy_v);

    clear L;
    
    % subsample velocities for flow *directions*
    sc = 1/4;%/2; % scale for resizing velocity
    [vx_r,x_r,y_r] = demresize(vx_v,x_v,y_v,sc); 
    vy_r = imresize(vy_v,sc);
    H_r = imresize(H_v,sc);
    
    figure; hold on; 
    subplot(1,3,1); imagesc(H_r); caxis([-1 1]); title("H_r");
    subplot(1,3,2); imagesc(vx_r); caxis([-1 1]); title("vx_r");
    subplot(1,3,3); imagesc(vy_r); caxis([-1 1]); title("vy_r");

    Hvx_r = inpaint_nans(H_r.*vx_r,4);
    Hvy_r = inpaint_nans(H_r.*vy_r,4);

    clear H_r vx_r vy_r;

    fprintf("done.\n");

    fprintf("  > Extrude velocity field...");
    % Extrude velocity
    % This coarse resolution round is only to constrain the faraway bits later
    % when the entire grid is inpainted. There's no real practical scientific
    % purpose to this coarse resolution, but it does ultimately help us create
    % a beautiful velocity map, and so it's well worth the extra ten minutes or
    % so that it takes for this section to run. 
    % 
    % Displacement vectors will be used for the flow directions when extrapolating terminus speeds. 
    %v_r = hypot(vx_r,vy_r);
    %dx = interp2(x_r,y_r,vx_r./v_r,X_v,Y_v); 
    %dy = interp2(x_r,y_r,vy_r./v_r,X_v,Y_v); 
    
    vx_ext = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_v,y_v,vx_v,0.03,35000);
    fprintf("done vx...");
    vy_ext = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_v,y_v,vy_v,0.03,35000);
    fprintf("done vy.\n");

    clear x_r y_r Hvx_r Hvy_r

    fprintf("  > Combine with original data and fill remaining nan regions..."); 
    % Wherever original data is nan and gridbin gives us data, 
    % overwrite vx_v, vy_v:  
    isf = isfinite(vx_ext) & ~isfinite(vx_v); 
    vx_v(isf) = vx_ext(isf); 
    isf = isfinite(vy_ext) & ~isfinite(vy_v);
    vy_v(isf) = vy_ext(isf);

    clear vx_ext vy_ext isf
    
    % Fill holes in the velocity data: 
    vx_v = regionfill(vx_v,isnan(vx_v)); 
    vy_v = regionfill(vy_v,isnan(vy_v)); 

    % Restore original values with large errors
    vx_v(I_higherr) = vx_v_orig(I_higherr);
    vy_v(I_higherr) = vy_v_orig(I_higherr);
     
    v_source_v(v_source_v==0 & isfinite(vx_v)) = max(v_source_v(:))+1; 

    fprintf("done.\n");

    clear vx_v_orig vy_v_orig I_higherr

    if CreateGeotiff
    
        fprintf('Writing GeoTiff files...');
        
        R = maprefcells([x_v(1) x_v(end)],[y_v(1) y_v(end)],[numel(x_v),numel(y_v)]);
        velfile_prim = "AntarcticVelocity"+erase(Velinterpolantfile,["GriddedInterpolants",...
            "Velocities"]);
        geotiffwrite("./GeoTiffFiles/"+velfile_prim+"_vx_EXTRUDED.tif",vx_v,R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/"+velfile_prim+"_vy_EXTRUDED.tif",vy_v,R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/"+velfile_prim+"_source_EXTRUDED.tif",v_source_v,R,'CoordRefSysCode','EPSG:3031');
    
        fprintf('Done.\n');
    
    end

    fprintf(' Creating gridded interpolants...');

    if isempty(Fus)
        load(Velinterpolantfile,"Fus","Fvs","Fsource","Fxerr","Fyerr");
    end

    Fus.Values = vx_v';
    Fvs.Values = vy_v';
    Fsource.Values = v_source_v';

    save(Velinterpolantfile+"_EXTRUDED.mat","Fus","Fvs","Fsource","Fxerr","Fyerr","-v7.3");

    clear vx_v vy_v v_source_v

    fprintf("done.\n");

end

%% -geom-
if contains(fields_to_extrude,'-geom-')

    fprintf("Processing geometry fields.\n");

    fprintf("  > Inpaint nan regions...");

    L = bwlabel(isnan(vx_g)); % Label the nan regions, and the ocean will be L=1.  
    vx_g = regionfill(vx_g,L>1); 
    vy_g = regionfill(vy_g,L>1);
    %v_source(L>1) = 4; % interpolated
    %v = hypot(vx_g,vy_g);

    clear L;
    
    sc = 1;%1/4; % scale for resizing velocity (for flow *directions* only) 
    [vx_r,x_r,y_r] = demresize(vx_g,x_g,y_g,sc);
    vy_r = imresize(vy_g,sc); 
    H_r = imresize(H_g,sc);

    figure; hold on; 
    subplot(1,3,1); imagesc(H_r); caxis([-1 1]); title("H_r");
    subplot(1,3,2); imagesc(vx_r); caxis([-1 1]); title("vx_r");
    subplot(1,3,3); imagesc(vy_r); caxis([-1 1]); title("vy_r");
    
    Hvx_r = inpaint_nans(H_r.*vx_r,4);
    Hvy_r = inpaint_nans(H_r.*vy_r,4);

    clear H_r vx_r vy_r

    fprintf("done.\n");

    %% Extrude surface elevation
    fprintf("  > Extrude geometry fields...");

    sf = filt2(s,x_g(2)-x_g(1),5e3,'lp'); 
    sg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,sf,0.1,10000);
    tmps = s; isf = isfinite(sg); tmps(isf) = sg(isf); 
    sg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,sf,0.01,10000);
    isf = isfinite(sg); tmps(isf) = sg(isf); 
    fprintf("done s...");

    clear sf sg isf

    bf = filt2(b,x_g(2)-x_g(1),5e3,'lp'); 
    bg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,bf,0.1,10000);
    tmpb = b; isf = isfinite(bg); tmpb(isf) = bg(isf); 
    bg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,bf,0.01,10000);
    isf = isfinite(bg); tmpb(isf) = bg(isf); 
    fprintf("done b. \n");

    clear bf bg isf

    %% commented out because it is better to interpolate/extrapolate rho 
    %% from bedmachine
    % rhof = filt2(rho,x_g(2)-x_g(1),5e3,'lp'); 
    % rhof(Fmask.Values'==0)=nan;
    % rhog = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,rhof,0.1,10000);
    % tmprho = rho; isf = isfinite(rhog); tmprho(isf) = rhog(isf); 
    % rhog = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,rhof,0.01,10000);
    % isf = isfinite(rhog); tmprho(isf) = rhog(isf); 
    % fprintf("done rho...");
    % 
    % clear rhof isf

    fprintf("  > Combine with original data and fill remaining nan regions...");

    s = regionfill(tmps,isnan(tmps)); 
    b = regionfill(tmpb,isnan(tmpb)); 
    %rho = regionfill(tmprho,isnan(tmprho));

    clear tmp*

    fprintf("done.\n");
    
    if CreateGeotiff
    
        fprintf('Writing GeoTiff files...');
        
        geomfile_prim = erase(Geominterpolantfile,"GriddedInterpolants_Geometry");

        R = maprefcells([x_g(1) x_g(end)],[y_g(1) y_g(end)],[numel(x_g),numel(y_g)]);
        geotiffwrite("./GeoTiffFiles/surface"+geomfile_prim+"_EXTRUDED.tif",s,R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/draft"+geomfile_prim+"_EXTRUDED.tif",b,R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/density"+geomfile_prim+"_EXTRUDED.tif",rho,R,'CoordRefSysCode','EPSG:3031');
    
        fprintf('Done.\n');
    
    end

    fprintf(' Creating gridded interpolants...');

    Fs.Values = s';
    Fb.Values = b';
    Frho.Values = rho';

    save(Geominterpolantfile+"_EXTRUDED.mat","Fs","Fb","Frho","FB","Fmask","-v7.3");

    clear tmp* s b rho

    fprintf("Done.\n");


end

%% -scalar-
if contains(fields_to_extrude,'-scalar-')

    fprintf("Processing scalar field.\n");

    fprintf("  > Inpaint nan regions...");

    L = bwlabel(isnan(vx_scal)); % Label the nan regions, and the ocean will be L=1.  
    vx_scal = regionfill(vx_scal,L>1); 
    vy_scal = regionfill(vy_scal,L>1);
    %v_source(L>1) = 4; % interpolated
    %v = hypot(vx_g,vy_g);

    clearvars L;
    
    sc = 1/4; % scale for resizing velocity (for flow *directions* only) 
    [vx_r,x_r,y_r] = demresize(vx_scal,x_scal,y_scal,sc);
    vy_r = imresize(vy_scal,sc); 
    H_r = imresize(H_scal,sc);

    % figure; hold on; 
    % subplot(1,3,1); imagesc(H_r); caxis([-1 1]); title("H_r");
    % subplot(1,3,2); imagesc(vx_r); caxis([-1 1]); title("vx_r");
    % subplot(1,3,3); imagesc(vy_r); caxis([-1 1]); title("vy_r");
    
    Hvx_r = inpaint_nans(H_r.*vx_r,4);
    Hvy_r = inpaint_nans(H_r.*vy_r,4);

    clearvars vx_r vy_r H_r

    fprintf("done.\n");

    %% Extrude scalar field
    fprintf("  > Extrude scalar field...");

    sf = filt2(scal,x_scal(2)-x_scal(1),5e3,'lp');
    sg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_scal,y_scal,sf,0.1,10000);
    tmps = scal; isf = isfinite(sg); tmps(isf) = sg(isf); 
    sg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_scal,y_scal,sf,0.01,10000);
    isf = isfinite(sg); tmps(isf) = sg(isf); 
    fprintf("done...");

    clearvars sf sg isf

    fprintf("  > Combine with original data and fill remaining nan regions...");

    scal = regionfill(tmps,isnan(tmps)); 

    clear tmp*

    fprintf("done.\n");
    
    if CreateGeotiff

        fprintf('Writing GeoTiff files \n');

        R = maprefcells([x_scal(1) x_scal(end)],[y_scal(1) y_scal(end)],[numel(x_scal),numel(y_scal)]);
        geotiffwrite("./GeoTiffFiles/scalarfield_EXTRUDED.tif",scal,R,'CoordRefSysCode','EPSG:3031');

        fprintf('Done.\n');

    end

    fprintf(' Creating gridded interpolants...');

    ScalarInterpolant.Values = scal';
   
    save(erase(fields_to_extrude,"-scalar-")+"_EXTRUDED.mat","ScalarInterpolant","-v7.3");

    fprintf("done.\n");

end







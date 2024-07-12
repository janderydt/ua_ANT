function Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,ScalarInterpolant,CreateGeotiff,fields_to_extrude)

% This function extrudes surface velocities, ice geometry (surface and draft), 
% densities or any other scalar quantity from the present-day ice edge of Antarctica along flowlines
% of ice flux (H*q). The processing chain is largely based on a script by
% C. Greene
% (https://github.com/chadagreene/ice-shelf-geometry/blob/main/code/flow_dem_extend.m).

%Velinterpolantfile = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities";
%Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jan-2000";
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
    CreateGeotiff = 1;
    fields_to_extrude = '-v-geom-';
end

if nargin==3
    fields_to_extrude = '-scalar-';
    CreateGeotiff = 1;
end

if nargin==4
    fields_to_extrude = '-scalar-';
end

fprintf("Interpolate data onto correct grids...");

if contains(fields_to_extrude,'-v-')

    fprintf("Load velocity and geometry interpolants...");

    load(Velinterpolantfile);
    load(Geominterpolantfile);
    
    fprintf("done.\n");

    x_v = Fus.GridVectors{1};
    y_v = Fus.GridVectors{2};
    [X_v,Y_v] = meshgrid(x_v,y_v);
    
    vx_v = Fus.Values';
    vy_v = Fvs.Values';
    v_source_v = Fsource.Values';
    H_v = Fs(X_v,Y_v)' - Fb(X_v,Y_v)';

end

if contains(fields_to_extrude,'-geom-')

    fprintf("Load velocity and geometry interpolants...");

    load(Velinterpolantfile);
    load(Geominterpolantfile);
    
    fprintf("done.\n");

    x_g = Fb.GridVectors{1};
    y_g = Fs.GridVectors{2};
    [X_g,Y_g] = meshgrid(x_g,y_g);

    s = Fs.Values'; b = Fb.Values';
    rho = Frho.Values';

    rho(s==0) = nan;
    s(s==0) = nan; b(b==0) = nan;
    
    H_g = s-b;

    vx_g = Fus(X_g,Y_g);
    vy_g =  Fvs(X_g,Y_g);
    %v_source_g =  Fsource(X_g,Y_g);

end

if contains(fields_to_extrude,'-scalar-')

    [nx,ny] = size(ScalarInterpolant.Values);
    filetoread = "Fields_to_extrude_scalar_nx"+string(nx)+"_ny"+string(ny)+".mat";

    if ~exist(filetoread,"file")
        
        fprintf("Load velocity and geometry interpolants...");
    
        load(Velinterpolantfile);
        load(Geominterpolantfile);
        
        fprintf("done.\n");
    
        x_scal = ScalarInterpolant.GridVectors{1};
        y_scal = ScalarInterpolant.GridVectors{2};
        [X_scal,Y_scal] = meshgrid(x_scal,y_scal);
     
        s = Fs(X_scal,Y_scal); b = Fb(X_scal,Y_scal);
        s(s==0) = nan; b(b==0) = nan;    
        H_scal = s-b;
        
        vx_scal = Fus(X_scal,Y_scal);
        vy_scal =  Fvs(X_scal,Y_scal);
        %v_source_g =  Fsource(X_g,Y_g);
    
        save(filetoread,"x_scal","y_scal","H_scal","vx_scal","vy_scal");
    else

        load(filetoread);

    end

    [X_scal,Y_scal] = meshgrid(x_scal,y_scal);
    scal = ScalarInterpolant.Values';
    scal(scal==0) = nan;

end

fprintf("done.\n");

%% -v-
if contains(fields_to_extrude,'-v-')
    fprintf("Processing velocity fields.\n");

    fprintf("  > Inpaint nan regions...");
    L = bwlabel(isnan(vx_v)); % Label the nan regions, and the ocean will be L=1.  
    vx_v = regionfill(vx_v,L>1); 
    vy_v = regionfill(vy_v,L>1);
    %v_source(L>1) = 4; % interpolated
    v = hypot(vx_v,vy_v);
    
    sc = 1/4; % scale for resizing velocity (for flow *directions* only) 
    [vx_r,x_r,y_r] = demresize(vx_v,x_v,y_v,sc); 
    vy_r = imresize(vy_v,sc); 
    H_r = imresize(H_v,sc);
    
    Hvx_r = inpaint_nans(H_r.*vx_r,4);
    Hvy_r = inpaint_nans(H_r.*vy_r,4);
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
    
    vx_ext = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_v,y_v,vx_v,0.03,25000);
    fprintf("done vx...");
    vy_ext = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_v,y_v,vy_v,0.03,25000);
    fprintf("done vy.\n");

    fprintf("  > Combine with original data and fill remaining nan regions...");
    % Create temporary grids so we don't accidentally do anything dumb and
    % overwrite any good data: 
    tmpvx = vx_v; 
    tmpvy = vy_v; 
    % Wherever gridbin gives us data, overwrite tmpvx, tmpvy:  
    isf = isfinite(vx_ext);
    %tmpvx(isf) = dx(isf).*vx_g(isf); 
    tmpvx(isf) = vx_ext(isf);
    %tmpvy(isf) = dy(isf).*vg(isf); 
    isf = isfinite(vy_ext);
    tmpvy(isf) = vy_ext(isf);
    
    % Fill holes in the velocity data: 
    tmpvx2 = regionfill(tmpvx,isnan(tmpvx)); 
    tmpvy2 = regionfill(tmpvy,isnan(tmpvy)); 
    
    vx_v = tmpvx2; 
    vy_v = tmpvy2; 
    v_source_v(v_source_v==0 & isfinite(vx_v)) = max(v_source_v(:))+1; 

    fprintf("done.\n");

    clear tmp*

    if CreateGeotiff
    
        fprintf('Writing GeoTiff files \n');
        
        R = maprefcells([x_v(1) x_v(end)],[y_v(1) y_v(end)],[numel(x_v),numel(y_v)]);
        velfile_prim = "AntarcticVelocity"+erase(Velinterpolantfile,["GriddedInterpolants",...
            "Velocities"]);
        geotiffwrite("./GeoTiffFiles/"+velfile_prim+"_vx_EXTRUDED.tif",vx_v,R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/"+velfile_prim+"_vy_EXTRUDED.tif",vy_v,R,'CoordRefSysCode','EPSG:3031');
        geotiffwrite("./GeoTiffFiles/"+velfile_prim+"_source_EXTRUDED.tif",v_source_v,R,'CoordRefSysCode','EPSG:3031');
    
        fprintf('Done.\n');
    
    end

    fprintf(' Creating gridded interpolants...');

    Fus.Values = vx_v';
    Fvs.Values = vy_v';
    Fsource.Values = v_source_v';

    save(Velinterpolantfile+"_EXTRUDED.mat","Fus","Fvs","Fsource","Fxerr","Fyerr","-v7.3");

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
    
    sc = 1/4; % scale for resizing velocity (for flow *directions* only) 
    [vx_r,x_r,y_r] = demresize(vx_g,x_g,y_g,sc);
    vy_r = imresize(vy_g,sc); 
    H_r = imresize(H_g,sc);
    
    Hvx_r = inpaint_nans(H_r.*vx_r,4);
    Hvy_r = inpaint_nans(H_r.*vy_r,4);

    fprintf("done.\n");

    %% Extrude surface elevation
    fprintf("  > Extrude geometry fields...");

    sf = filt2(s,x_g(2)-x_g(1),5e3,'lp'); 
    sg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,sf,0.1,10000);
    tmps = s; isf = isfinite(sg); tmps(isf) = sg(isf); 
    sg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,sf,0.01,10000);
    isf = isfinite(sg); tmps(isf) = sg(isf); 
    fprintf("done s...");

    clear sf isf

    bf = filt2(b,x_g(2)-x_g(1),5e3,'lp'); 
    bg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,bf,0.1,10000);
    tmpb = b; isf = isfinite(bg); tmpb(isf) = bg(isf); 
    bg = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,bf,0.01,10000);
    isf = isfinite(bg); tmpb(isf) = bg(isf); 
    fprintf("done b...");

    clear bf isf

    rhof = filt2(rho,x_g(2)-x_g(1),5e3,'lp'); 
    rhof(Fmask.Values'==0)=nan;
    rhog = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,rhof,0.1,10000);
    tmprho = rho; isf = isfinite(rhog); tmprho(isf) = rhog(isf); 
    rhog = ExtrudeField(x_r,y_r,Hvx_r,Hvy_r,x_g,y_g,rhof,0.01,10000);
    isf = isfinite(rhog); tmprho(isf) = rhog(isf); 
    fprintf("done rho...");

    clear rhof isf

    fprintf("  > Combine with original data and fill remaining nan regions...");

    s = regionfill(tmps,isnan(tmps)); 
    b = regionfill(tmpb,isnan(tmpb)); 
    rho = regionfill(tmprho,isnan(tmprho));

    clear tmp*

    fprintf("done.\n");
    
    if CreateGeotiff
    
        fprintf('Writing GeoTiff files \n');
        
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

    fprintf("done.\n");


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







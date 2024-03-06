function plot_MassBalanceComponents

%% This function plots mass balance components for grounded and/or floating ice
%% Grounded ice: surface mass balance, ice-shelf flux and Ocean flux for individual basins
%% Floating ice: surface mass balance, grounding line flux and calving flux for indididual ice shelves

grounded = 1;
floating = 0;

addpath(getenv("froot_tools"));

froot_ua = getenv("froot_ua")+"/cases/ANT/";

years = [1795 1317 1034 1381]; % 
year_datestr = ["2000","2009","2014","2018"];
year_datenum = datenum("0106" + year_datestr,"ddmmyyyy");

load("../ANT_Data/ANT_Interpolants/ScatteredInterpolants_SMB.mat","Fsmb_RACMO","Fsmb_MAR");

% RACMO climatology between year_datestr(1) and year_datestr(end)
Istart = find(contains(string(Fsmb_RACMO.years),year_datestr(1)));
Iend = find(contains(string(Fsmb_RACMO.years),year_datestr(end)));
dn = Iend-Istart+1;
RACMO_smb = 0;
for ii=1:dn
    yrstr = "yr"+string(double(year_datestr(1))-1+ii);
    RACMO_smb = RACMO_smb + Fsmb_RACMO.(yrstr).Values;
end
Fsmb_RACMO_climatology = Fsmb_RACMO.climatology; Fsmb_RACMO_climatology.Values = RACMO_smb/dn;

% MAR climatology between year_datestr(1) and year_datestr(end)
Istart = find(contains(string(Fsmb_MAR.years),year_datestr(1)));
Iend = find(contains(string(Fsmb_MAR.years),year_datestr(end)));
dn = Iend-Istart+1;
MAR_smb = 0;
for ii=1:dn
    yrstr = "yr"+string(double(year_datestr(1))-1+ii);
    MAR_smb = MAR_smb + Fsmb_MAR.(yrstr).Values;
end
Fsmb_MAR_climatology = Fsmb_MAR.climatology; Fsmb_MAR_climatology.Values = MAR_smb/dn;


%% GROUNDED ICE SHEET
if grounded

%% drainage basins from IMBIE - uses C Greene routines below
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);

% Only keep ice rises and islands larger than 250km2. We don't treat
% smaller ice rises and/or islands seperately
Bx = B.x{1}; Bxtmp = [];
By = B.y{1}; Bytmp = [];
I = [0; find(isnan(Bx)); numel(Bx)+1];
for ii=1:numel(I)-1
    xtmp = Bx(I(ii)+1:I(ii+1)-1);
    ytmp = By(I(ii)+1:I(ii+1)-1);
    P = polyshape(xtmp(:),ytmp(:));
    A = area(P)/1e6; % square km
    if A>250
        Bxtmp = [Bxtmp(:);xtmp(:);NaN];
        Bytmp = [Bytmp(:);ytmp(:);NaN];
    end
end
B.x{1} = Bxtmp;
B.y{1} = Bytmp;

%% To compare modelled fluxes with those from the literature, we need the flux gates that were used
%% to calculate drainage. 
%% flux gates from Gardner et al NSIDC
froot_gates = getenv("froot_data")+"/Measures/FluxGates/line_shp/";
GL0 = shaperead(froot_gates+"GL0.shp");
FG1 = shaperead(froot_gates+"FG1.shp");
FG2 = shaperead(froot_gates+"FG2.shp");

for yy=1:numel(years)

    load(froot_ua+"/ANT_Inverse/ANT_Inverse_"+string(years(yy))+"/ANT_Inverse_"+string(years(yy))+"-RestartFile.mat");

    % Surface mass balance grounded (RACMO climatology)
    smb = Fsmb_RACMO_climatology(MUA.coordinates(:,1),MUA.coordinates(:,2));
    %smb = Fsmb_MAR_climatology(MUA.coordinates(:,1),MUA.coordinates(:,2));
    smb_grounded = FEintegrate2D(CtrlVarInRestartFile,MUA,smb.*F.GF.node);
    %smb_floating = FEintegrate2D(CtrlVarInRestartFile,MUA,smb.*(1-F.GF.node));

    % Fluxes through gates used in Gardner et al. 2018
    Fub=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
    Fvb=Fub; Frho=Fub; Fh=Fub; 
    Fvb.Values=F.vb;
    Frho.Values=F.rho;
    Fh.Values=F.h;
    [GL0.x,GL0.y,GL0.flux,~,~]=FluxAcrossBoundary([GL0.X(:) GL0.X([2:end 1])'],[GL0.Y(:) GL0.Y([2:end 1])'],Fub,Fvb,Fh,Frho);
    [FG1.x,FG1.y,FG1.flux,~,~]=FluxAcrossBoundary([FG1.X(:) FG1.X([2:end 1])'],[FG1.Y(:) FG1.Y([2:end 1])'],Fub,Fvb,Fh,Frho);
    [FG2.x,FG2.y,FG2.flux,~,~]=FluxAcrossBoundary([FG2.X(:) FG2.X([2:end 1])'],[FG2.Y(:) FG2.Y([2:end 1])'],Fub,Fvb,Fh,Frho);

    % Order Fluxes per basin
    [~,GL0.ind] = Define_Quantity_PerBasin(GL0.x,GL0.y,B);
    [~,FG1.ind] = Define_Quantity_PerBasin(FG1.x,FG1.y,B);
    [~,FG2.ind] = Define_Quantity_PerBasin(FG2.x,FG2.y,B);

    % Order Ua elements per basin
    basin = Define_Quantity_PerBasin(MUA.xEle,MUA.yEle,B);
    for ii=1:numel(B.x)
        B.UaEle{ii} = basin(ii).ind;
    end

    % Obtain Ua fluxes across the grounding line (qGL) into floating areas
    [B,GL] = Calc_UaGLFlux_PerBasin(MUA,F,GF,B,CtrlVarInRestartFile);

    % Obtain Ua fluxes across the open boundary (qOB) into the ocean
    B = Calc_UaOBFlux_PerBasin(MUA,F,GF,B,CtrlVarInRestartFile); 
    
    % Sum values of SMB, qGL and qOB for each basin
    for ii=1:numel(B.x)
        B.SMB{ii} = sum(smb_grounded(B.UaEle{ii}),'omitmissing')/1e9; 
        B.qGL_tot{ii} = sum(B.qGL{ii},'omitmissing')/1e12; 
        B.qOB_tot{ii} = sum(B.qOB{ii},'omitmissing')/1e12; 
    end

    % Initalize some variables for plotting
    if yy==1
        CM1 = crameri('lajolla',64);
        cmin1 = 0;
        cmax1 = 350;
        CMtmp = flipdim(crameri('vik',18),1);
        CM2 = [CMtmp(1:9,:);0.95*ones(2,3);CMtmp(10:18,:)];
        cmin2 = -150;
        cmax2 = 150;
        CMtmp = crameri('vik',14);
        CM3 = [CMtmp(1:7,:);0.95*ones(2,3);CMtmp(8:14,:)];
        cmin3 = -100;
        cmax3 = 100;
        % save data for first year in seperate variable
        B0 = B;
        Fh0 = Fh;
        Fub0 = Fub;
        Fvb0 = Fvb;
    end

    %% -------- %%
    %% PLOTTING %%
    %% -------- %%
    figure; hold on;
    PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,F.h-Fh0(MUA.coordinates(:,1),MUA.coordinates(:,2))); 
    title("\Deltah "+year_datestr(yy)+"-"+year_datestr(1)); cb = colorbar; caxis([-100 100]);
    PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF);

    speed0 = hypot(Fub0(MUA.coordinates(:,1),MUA.coordinates(:,2)),Fvb0(MUA.coordinates(:,1),MUA.coordinates(:,2)));
    speed = hypot(F.ub,F.vb);
    figure; hold on; PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,speed-speed0); 
    title("\Deltau "+year_datestr(yy)+"-"+year_datestr(1)); cb = colorbar; caxis([-1000 1000]);
    PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF);

    H=fig('units','inches','width',130*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');

    tlo_fig = tiledlayout(1,4,"TileSpacing","compact");
    for i = 1:4
        ax_fig(i) = nexttile(tlo_fig,i); hold on;
    end

    % To exclude values for the ice rises and islands, set startB=2
    startB=2;

    %% SMB
    for ii=startB:numel(B.x)
        pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
        [xc,yc] = centroid(pgon); 
        h1(ii) = plot(ax_fig(1),pgon,'FaceColor',CM1(max(min(round(1+63/(cmax1-cmin1)*(B.SMB{ii}-cmin1)),64),1),:),'FaceAlpha',1,'EdgeColor',[0.75 0.75 0.75]); 
        h1(ii).LineStyle = '-';
        text(ax_fig(1),xc,yc,{B.name{ii};string(round(B.SMB{ii}))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h1(ii).FaceColor*.25)
    end
    colormap(ax_fig(1),CM1); cb1=colorbar(ax_fig(1)); 
    cb1.Label.String = "Flux [Gt/yr]"; 
    caxis(ax_fig(1),[cmin1 cmax1]);
    title(ax_fig(1),"SMB 2000-2018: "+string(round(sum([B.SMB{startB:end}])))+" Gt/yr");


    %% GL Flux
    for ii=startB:numel(B.x)
        pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
        [xc,yc] = centroid(pgon); 
        h2(ii) = plot(ax_fig(2),pgon,'FaceColor',CM1(max(min(round(1+63/(cmax1-cmin1)*(B.qGL_tot{ii}-cmin1)),64),1),:),'FaceAlpha',1,'EdgeColor',[0.75 0.75 0.75]); 
        h2(ii).LineStyle = '-';
        text(ax_fig(2),xc,yc,{B.name{ii};string(round(B.qGL_tot{ii}))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h2(ii).FaceColor*.25)
    end
    colormap(ax_fig(2),CM1); 
    caxis(ax_fig(2),[cmin1 cmax1]);  
    title(ax_fig(2),"GL flux "+year_datestr(yy)+": "+string(round(sum([B.qGL_tot{startB:end}])))+" Gt/yr");

    %% Ocean Flux
    for ii=startB:numel(B.x)
        pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
        [xc,yc] = centroid(pgon); 
        h3(ii) = plot(ax_fig(3),pgon,'FaceColor',CM1(max(min(round(1+63/(cmax1-cmin1)*(B.qOB_tot{ii}-cmin1)),64),1),:),'FaceAlpha',1,'EdgeColor',[0.75 0.75 0.75]); 
        h3(ii).LineStyle = '-';
        text(ax_fig(3),xc,yc,{B.name{ii};string(round(B.qOB_tot{ii}))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h3(ii).FaceColor*.25)
    end
    colormap(ax_fig(3),CM1); 
    caxis(ax_fig(3),[cmin1 cmax1]);  
    title(ax_fig(3),"Ocean flux "+year_datestr(yy)+": "+string(round(sum([B.qOB_tot{startB:end}])))+" Gt/yr");


    %% SMB - GL Flux - Ocean Flux
    for ii=startB:numel(B.x)
        pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
        [xc,yc] = centroid(pgon); 
        h4(ii) = plot(ax_fig(4),pgon,'FaceColor',CM2(max(min(round(1+19/(cmax2-cmin2)*(B.SMB{ii}-B.qGL_tot{ii}-B.qOB_tot{ii}-cmin2)),20),1),:),'FaceAlpha',1,'EdgeColor',[0.75 0.75 0.75]); 
        h4(ii).LineStyle = '-';
        text(ax_fig(4),xc,yc,{B.name{ii};string(round(B.SMB{ii}-B.qGL_tot{ii}-B.qOB_tot{ii}))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h4(ii).FaceColor*.25)
    end
    colormap(ax_fig(4),CM2); cb4=colorbar(ax_fig(4)); 
    cb4.Label.String = "\DeltaFlux [Gt/yr]";
    caxis(ax_fig(4),[cmin2 cmax2]);
    title(ax_fig(4),"Net mass balance "+year_datestr(yy)+": "+string(round(sum([B.SMB{startB:end}]-[B.qGL_tot{startB:end}]-[B.qOB_tot{startB:end}])))+" Gt/yr");

    %% DEAL WITH AXIS ETC
    for i=1:4
        CtrlVarInRestartFile.PlotGLs=0;
        [xGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF);%,[],[],[],'color','k');
        plot(ax_fig(i),xGL/1e3,yGL/1e3,'-k');
        if i>1
            yticklabels(ax_fig(i),{});
        end
        plot(ax_fig(i),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
        axis(ax_fig(i),"equal"); axis(ax_fig(i),"tight");
        grid(ax_fig(i),"on"); box(ax_fig(i),"on");  
    end
    xlabel(tlo_fig,"psx [km]"); ylabel(tlo_fig,"psy [km]"); 

    figure; hold on;

    %% Change in GL Flux
    if yy>1
        figure; hold on;
        for ii=startB:numel(B.x)
            pgon = polyshape(B.x{ii},B.y{ii},'Simplify',false); 
            [xc,yc] = centroid(pgon); 
            h = plot(pgon,'FaceColor',CM3(max(min(round(1+15/(cmax3-cmin3)*(B.qGL_tot{ii}-B0.qGL_tot{ii}-cmin3)),16),1),:),'FaceAlpha',1,'EdgeColor',[0.75 0.75 0.75]); 
            h.LineStyle = '-';
            text(xc,yc,{B.name{ii};string(round(B.qGL_tot{ii}-B0.qGL_tot{ii}))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h.FaceColor*.25)
        end
        PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF,[],[],[],'color','k');
        plot(MUA.Boundary.x,MUA.Boundary.y,'-k');
        axis equal; grid on; box on;
        colormap(CM3); cb=colorbar; 
        cb.Label.String = "Change in GL Flux [Gt/yr]"; 
        caxis([cmin3 cmax3]);
        title("Change in GL Flux: "+year_datestr(yy)+"-"+year_datestr(1));
    end

    %% For each year, compare Ua flux through Ua GL to Ua flux through GL0, FG1, FG2 from Gardner et al 2018
    % Assemble data in n*4 matrix, where n is the number of basins:
    % [qGLUa,qGL0,qFG1,qFG2]
    M = zeros(4,numel(B.x)); xlabels=[];
    for ii=1:numel(B.x)
        M(1,ii) = B.qGL_tot{ii}+B.qOB_tot{ii}; % Ua GL flux
        M(2,ii) = sum(GL0.flux(GL0.ind==ii)/1e12);  % GL0 flux
        M(3,ii) = sum(FG1.flux(FG1.ind==ii)/1e12); % FG1 flux
        M(4,ii) = sum(FG2.flux(FG2.ind==ii)/1e12); % FG2 flux
        xlabels = [xlabels string(B.name{ii})];
    end
    % plot data
    
    H2=fig('units','inches','width',130*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');

    bar([1:numel(B.x)],M');
    xticks(1:numel(B.x)); xticklabels(xlabels);
    title(year_datestr(yy));
    
end


end

if floating

    addpath(getenv("froot_ua")+"/UaPICO_master/");

    for yy=1%:numel(years)

        load(froot_ua+"/ANT_Inverse/ANT_Inverse_"+string(years(yy))+"/ANT_Inverse_"+string(years(yy))+"-RestartFile.mat");
        PICO_opts.FloatingCriteria='GLthreshold';
        %PICO_opts.algorithm='graph';
        %PICO_opts.ContinentArea=5e10;
        %PICO_opts.PICOres=1000;

        PICO_opts = PICO_DefaultParameters(MUA,PICO_opts);
        [MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
        [ShelfID,PBOX,Ak,floating] = PICO_IdentifyIceShelvesWatershedOption(UserVarInRestartFile,CtrlVarInRestartFile,MUA,GF,PICO_opts);

        x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
        figure(yy);
        ShelfIDtemp = ShelfID; ShelfIDtemp(ShelfID==0) = nan;
        PlotMeshScalarVariable(CtrlVarInRestartFile, MUA, ShelfIDtemp);
        hold on;
        PlotGroundingLines(CtrlVarInRestartFile,MUA,GF,[],[],[],'k');
        colormap(lines(max(ShelfID)))
        colorbar off
        
        for shelf_i=1:max(ShelfID)
            average_x_loc_per_shelf(shelf_i) = mean(x(ShelfID==shelf_i));
            average_y_loc_per_shelf(shelf_i) = mean(y(ShelfID==shelf_i));
            text(mean(x(ShelfID==shelf_i))/CtrlVarInRestartFile.PlotXYscale,mean(y(ShelfID==shelf_i))/CtrlVarInRestartFile.PlotXYscale,num2str(shelf_i),'Color','k','FontWeight','bold','FontSize',12,'BackgroundColor','w','Margin',0.1,'EdgeColor','r');
        end
        title('Shelf IDs');
    end

end

 return
%% ICE SHELVES
% Surface mass balance floating
SMBFloating(ii) = sum(FEintegrate2D(CtrlVarInRestartFile,MUA,smb.*(1-F.GF.node)),'omitmissing')/1e9; %m3/yr water equivalent to Gt/yr
    
% Calving flux
[~,~,qIF,~]= FluxAcrossModelBoundary(CtrlVarInRestartFile,MUA,F.ub,F.vb,F.h,F.b,F.B,F.rho);
IFFlux(ii) = -abs(sum(qIF,'all')/1e12); %kg/yr to Gt/yr

% Basal melt (assuming balanced melt)
BasalMelt(ii) = -(SMBFloating(ii)+GLFlux(ii)+IFFlux(ii));

figure; hold on;
b1 = bar([1:4],[SMBGrounded(:) -GLFlux(:) NetGrounded(:)]);   
xticks(gca,[1:4]);
xticklabels(gca,["2000","2009","2014","2018"]);
ylabel('Flux [Gt/yr]');
legend(["smb grounded MAR","grounding line flux","net"]);
grid on;
box on;


figure; hold on;
b2 = bar([1:4],[GLFlux(:) SMBFloating(:) IFFlux(:) BasalMelt(:)]);
xticks(gca,[1:4]);
xticklabels(gca,["2000","2009","2014","2018"]);
ylabel('Flux [Gt/yr]');
legend(["grounding line flux","MAR smb shelves","calving flux","balanced melt flux"]);
grid on;
box on;

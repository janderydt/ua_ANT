function plot_Perturbations

addpath(getenv("froot_tools"));

froot_ua = getenv("froot_ua")+"/cases/ANT/";

baseline = 1814; % n=3, m=3, 2000
%baseline = 1421; % n=3, m=3, 2009
%baseline = 1915; % n=3, m=3, 2014
%baseline = 1792; % n=3, m=3, 2018

%target =  1421; % n=3, m=3, 2009
%target =  1915; % n=3, m=3, 2014
target =  1792; % n=3, m=3, 2018

% n=3, m=3, 2018
perturbations = [1905 1126 1913 1632];%,...
    %1097 1278 1141 nan,...
    %1800 nan nan nan];

Table = readtable("../ANT_Diagnostic/RunTable.csv");

%% basins from C Greene
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);

%% load data baseline simulation
load(froot_ua+"/ANT_Diagnostic/ANT_Diagnostic_"+string(baseline)+"/ANT_Diagnostic_"+string(baseline)+"-RestartFile.mat");

Fbaseline.speed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb));
Fbaseline.h = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
Fbaseline.AGlen = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.AGlen);
Fbaseline.C = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.C);

baselineMUA = MUA;
baselineF = F;

% Obtain Ua fluxes across the grounding line (qGL) into floating areas
[Bbaseline,~] = Calc_UaGLFlux_PerBasin(MUA,F,GF,B,CtrlVarInRestartFile);

% Obtain Ua fluxes across the open boundary (qOB) into the ocean
Bbaseline = Calc_UaOBFlux_PerBasin(MUA,F,GF,Bbaseline,CtrlVarInRestartFile); 
    
% Sum values of SMB, qGL and qOB for each basin
for ii=1:numel(Bbaseline.x) 
    Bbaseline.qGL_tot{ii} = sum(Bbaseline.qGL{ii},'omitmissing')/1e12; 
    Bbaseline.qOB_tot{ii} = sum(Bbaseline.qOB{ii},'omitmissing')/1e12; 
    Bbaseline.qtot{ii} = Bbaseline.qGL_tot{ii}+Bbaseline.qOB_tot{ii}; 
end

%% load data target simulation
load(froot_ua+"/ANT_Diagnostic/ANT_Diagnostic_"+string(target)+"/ANT_Diagnostic_"+string(target)+"-RestartFile.mat");

Ftarget.speed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb));
Ftarget.h = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
Ftarget.AGlen = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.AGlen);
Ftarget.C = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.C);

targetMUA = MUA;
targetF = F;

% Obtain Ua fluxes across the grounding line (qGL) into floating areas
[Btarget,~] = Calc_UaGLFlux_PerBasin(MUA,F,GF,B,CtrlVarInRestartFile);

% Obtain Ua fluxes across the open boundary (qOB) into the ocean
Btarget = Calc_UaOBFlux_PerBasin(MUA,F,GF,Btarget,CtrlVarInRestartFile); 
    
% Sum values of SMB, qGL and qOB for each basin
for ii=1:numel(Btarget.x) 
    Btarget.qGL_tot{ii} = sum(Btarget.qGL{ii},'omitmissing')/1e12; 
    Btarget.qOB_tot{ii} = sum(Btarget.qOB{ii},'omitmissing')/1e12; 
    Btarget.qtot{ii} = Btarget.qGL_tot{ii}+Btarget.qOB_tot{ii}; 
end

np = numel(perturbations);

%% load data diagnostic geometry perturbation experiments
for pp = 1:np
    load(froot_ua+"/ANT_Diagnostic/ANT_Diagnostic_"+string(perturbations(pp))+"/ANT_Diagnostic_"+string(perturbations(pp))+"-RestartFile.mat");
    
    Fperturb(pp).speed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb));
    Fperturb(pp).h = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
    Fperturb(pp).AGlen = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.AGlen);
    Fperturb(pp).C = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.C);
    
    % Obtain Ua fluxes across the grounding line (qGL) into floating areas
    [Bperturb(pp).B,~] = Calc_UaGLFlux_PerBasin(MUA,F,GF,B,CtrlVarInRestartFile);
    
    % Obtain Ua fluxes across the open boundary (qOB) into the ocean
    Bperturb(pp).B = Calc_UaOBFlux_PerBasin(MUA,F,GF,Bperturb(pp).B,CtrlVarInRestartFile); 
        
    % Sum values of SMB, qGL and qOB for each basin
    for ii=1:numel(Bperturb(pp).B.x) 
        Bperturb(pp).B.qGL_tot{ii} = sum(Bperturb(pp).B.qGL{ii},'omitmissing')/1e12; 
        Bperturb(pp).B.qOB_tot{ii} = sum(Bperturb(pp).B.qOB{ii},'omitmissing')/1e12; 
        Bperturb(pp).B.qtot{ii} = Bperturb(pp).B.qGL_tot{ii}+Bperturb(pp).B.qOB_tot{ii}; 
    end

    I = find(Table.ExpID == perturbations(pp));
    perturbationdata(pp).description = string(Table.Comments(I)); 

end

%% -------- %%
%% PLOTTING %%
%% -------- %%

%% 1- Difference in flux across the GL and into the open ocean between perturbation and baseline
% Initalize some variables for plotting
CMtmp = crameri('vik',18);
CM = [CMtmp(1:9,:);0.9*ones(2,3);CMtmp(10:18,:)];
cmin = -125;
cmax = 125;

H=fig("units","inches","width",130*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(1,np,"TileSpacing","compact");
for i = 1:np
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

for pp = 1:np

    startB=2;

    for ii=startB:numel(B.x)

        % plotting 100*(qGL(perturb)-qGL(baseline))./(qGL(target)-qGL(baseline))
        data(ii) = (Bperturb(pp).B.qtot{ii}-Bbaseline.qtot{ii});%./(Btarget.qtot{ii}-Bbaseline.qtot{ii});

        pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
        [xc,yc] = centroid(pgon); 
        h(ii) = plot(ax_fig(pp),pgon,'FaceColor',CM(max(min(round(1+17/(cmax-cmin)*(data(ii)-cmin)),18),1),:),'FaceAlpha',1,'EdgeColor',[0.75 0.75 0.75]); 
        h(ii).LineStyle = '-';
        text(ax_fig(pp),xc,yc,{B.name{ii};string(round(data(ii)))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h(ii).FaceColor*.25);

    end

    if pp==1 || pp==4
        plot(ax_fig(pp),targetMUA.Boundary.x/1e3,targetMUA.Boundary.y/1e3,'-m');
    end

    if pp==3 || pp==4
        CtrlVarInRestartFile.PlotGLs=0;
        [xGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,targetMUA,targetF.GF);%,[],[],[],'color','k');
        plot(ax_fig(pp),xGL/1e3,yGL/1e3,'-m');
    end

    colormap(ax_fig(pp),CM);    
    caxis(ax_fig(pp),[cmin cmax]);
    title(ax_fig(pp),perturbationdata(pp).description);

end

for i=1:np
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

cb=colorbar(ax_fig(np)); 
cb.Layout.Tile = 'east'; 
cb.Label.String = "qGL_{pert}+qOB_{pert}-qGL_{base}-qOB_{base} [Gt/yr]"; 

%% 2- Difference in surface speed between perturbation and baseline
% Initalize some variables for plotting
CMtmp = crameri('vik',14);
CM = [CMtmp(1:7,:);0.9*ones(2,3);CMtmp(8:14,:)];
cmin = -100;
cmax = 100;

H=fig("units","inches","width",130*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(1,np,"TileSpacing","compact");
for i = 1:np
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

for pp = 1:np

    CtrlVarInRestartFile.PlotXYscale = 1e3;
    PlotNodalBasedQuantities_JDR(ax_fig(pp),targetMUA.connectivity,targetMUA.coordinates,Fperturb(pp).speed(targetMUA.coordinates(:,1),targetMUA.coordinates(:,2)) -...
        Fbaseline.speed(targetMUA.coordinates(:,1),targetMUA.coordinates(:,2)),CtrlVarInRestartFile); hold on;

    for ii=startB:numel(B.x)

        % plotting 100*(qGL(perturb)-qGL(baseline))./(qGL(target)-qGL(baseline))
        data(ii) = (Bperturb(pp).B.qtot{ii}-Bbaseline.qtot{ii});%./(Btarget.qtot{ii}-Bbaseline.qtot{ii});

        plot(ax_fig(pp),B.x{ii}/1e3,B.y{ii}/1e3,'color',[0.75 0.75 0.75]); 
        
        pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
        [xc,yc] = centroid(pgon); 
        text(ax_fig(pp),xc,yc,{B.name{ii};string(round(data(ii)))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h(ii).FaceColor*.25);

    end

    if pp==1 || pp==4
        plot(ax_fig(pp),targetMUA.Boundary.x/1e3,targetMUA.Boundary.y/1e3,'-m');
    end

    if pp==3 || pp==4
        CtrlVarInRestartFile.PlotGLs=0;
        [xGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,targetMUA,targetF.GF);%,[],[],[],'color','k');
        plot(ax_fig(pp),xGL/1e3,yGL/1e3,'-m');
    end

    colormap(ax_fig(pp),CM);    
    caxis(ax_fig(pp),[cmin cmax]);

    title(ax_fig(pp),perturbationdata(pp).description);

end

% deal with axes etc
for i=1:np
    CtrlVarInRestartFile.PlotGLs=0;
    [xGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,baselineMUA,baselineF.GF);%,[],[],[],'color','k');
    plot(ax_fig(i),xGL/1e3,yGL/1e3,'-k');
    if i>1
        yticklabels(ax_fig(i),{});
    end
    plot(ax_fig(i),baselineMUA.Boundary.x/1e3,baselineMUA.Boundary.y/1e3,'-k');
    axis(ax_fig(i),"equal"); axis(ax_fig(i),"tight");
    grid(ax_fig(i),"on"); box(ax_fig(i),"on");
end
xlabel(tlo_fig,"psx [km]"); ylabel(tlo_fig,"psy [km]"); 

cb=colorbar(ax_fig(np)); 
cb.Layout.Tile = 'east'; 
cb.Label.String = "Change in speed"; 

%% 3- Difference in surface speed between target and baseline
% Initalize some variables for plotting
CMtmp = crameri('vik',14);
CM = [CMtmp(1:7,:);0.9*ones(2,3);CMtmp(8:14,:)];
cmin = -100;
cmax = 100;

H=fig("units","inches","width",55*12/72.27,"height",55*12/72.27,"fontsize",14,"font","Helvetica");

subplot("position",[0.1 0.1 0.85 0.85]); hold on;

CtrlVarInRestartFile.PlotXYscale = 1e3;
PlotNodalBasedQuantities_JDR(gca,targetMUA.connectivity,targetMUA.coordinates,Ftarget.speed(targetMUA.coordinates(:,1),targetMUA.coordinates(:,2)) -...
        Fbaseline.speed(targetMUA.coordinates(:,1),targetMUA.coordinates(:,2)),CtrlVarInRestartFile); hold on;

CtrlVarInRestartFile.PlotGLs=0;
[xGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,targetMUA,targetF.GF);%,[],[],[],'color','k');
g(1)=plot(gca,xGL/1e3,yGL/1e3,'-k');

[xGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,baselineMUA,baselineF.GF);%,[],[],[],'color','k');
g(2)=plot(gca,xGL/1e3,yGL/1e3,'-m');

for ii=startB:numel(B.x)

    plot(gca,B.x{ii}/1e3,B.y{ii}/1e3,'color',[0.75 0.75 0.75]); 

    pgon = polyshape(B.x{ii}/1e3,B.y{ii}/1e3,'Simplify',false); 
    [xc,yc] = centroid(pgon); 
    text(gca,xc,yc,{B.name{ii};string(round(data(ii)))+" Gt/yr"},'vert','middle','horiz','center','fontsize',8,'color',h(ii).FaceColor*.25);

end

plot(gca,targetMUA.Boundary.x/1e3,targetMUA.Boundary.y/1e3,'-k');
plot(gca,baselineMUA.Boundary.x/1e3,baselineMUA.Boundary.y/1e3,'-m');
   
colormap(gca,CM);    
caxis(gca,[cmin cmax]);

title(gca,"Target - baseline");

axis(gca,"equal"); axis(gca,"tight");
grid(gca,"on"); box(gca,"on");

xlabel(gca,"psx [km]"); ylabel(gca,"psy [km]"); 

cb=colorbar; 
cb.Label.String = "Change in speed"; 

function plotPerturbations

addpath(getenv("froot_tools"));

froot_ua = getenv("froot_ua")+"/cases/ANT/";


baseline = 1814; % n=3, m=3, 2000
%baseline = 1421; % n=3, m=3, 2009
%baseline = 1915; % n=3, m=3, 2014
%baseline = 1792; % n=3, m=3, 2018

target =  1421; % n=3, m=3, 2009
target =  1915; % n=3, m=3, 2014
%target =  1792; % n=3, m=3, 2018


% n=3, m=3, 2018
perturbations = [1905 1126 1913 1632,...
    1097 1278 1141 nan,...
    1800 nan nan nan];

Table = readtable("../ANT_Diagnostic/RunTable.csv");

%CM = flipdim(othercolor('RdYlBu5'),1);
CM=crameri('vik',11); CM(6,:)=[1 1 1];

%% baseline
load(froot_ua+"/ANT_Diagnostic/ANT_Diagnostic_"+string(baseline)+"/ANT_Diagnostic_"+string(baseline)+"-RestartFile.mat");

Fspeed_baseline = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb));
Fh_baseline = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
FAGlen_baseline = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.AGlen);
FC_baseline = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.C);

%% basins from C Greene
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);

%% baseline vs target
H=fig("units","inches","width",90*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");

% tlo_fig = tiledlayout(2,3,"TileSpacing","compact");
% for i = 1:6
%     ax_fig(i) = nexttile(tlo_fig,i); hold on;
% end
ax_fig=gca; hold on;

load(froot_ua+"/ANT_Diagnostic/ANT_Diagnostic_"+string(target)+"/ANT_Diagnostic_"+string(target)+"-RestartFile.mat");

CtrlVarInRestartFile.PlotXYscale=1e3;

% velocity
ii=1;
PlotNodalBasedQuantities_JDR(ax_fig(ii),MUA.connectivity,MUA.coordinates,hypot(F.ub,F.vb) - Fspeed_baseline(MUA.coordinates(:,1),MUA.coordinates(:,2)),CtrlVarInRestartFile); hold on;
colormap(ax_fig(ii),CM); caxis(ax_fig(ii),[-300 300]);
PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVarInRestartFile.PlotXYscale,MUA.Boundary.y/CtrlVarInRestartFile.PlotXYscale,'-k')
%set(ax_fig(ii),'ColorScale','log');
grid(ax_fig(ii),"off"); 
box(ax_fig(ii),"off");
axis(ax_fig(ii),"equal","off");
title(ax_fig(ii),string(target));

% icefront

% ice shelf thickness

% grounded ice thickness

% AGlen

% C
return


%% perturbations
H=fig("units","inches","width",120*12/72.27,"height",70*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(3,4,"TileSpacing","tight");
for i = 1:numel(perturbations)
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

for ii=1:numel(perturbations)

    if ~isnan(perturbations(ii))
        load(froot_ua+"/ANT_Diagnostic/ANT_Diagnostic_"+string(perturbations(ii))+"/ANT_Diagnostic_"+string(perturbations(ii))+"-RestartFile.mat");
        
        perturbationdata(ii).dspeed = hypot(F.ub,F.vb) - Fspeed_baseline(MUA.coordinates(:,1),MUA.coordinates(:,2));
        perturbationdata(ii).dh = F.h - Fh_baseline(MUA.coordinates(:,1),MUA.coordinates(:,2));
        perturbationdata(ii).dAGlen = F.AGlen - FAGlen_baseline(MUA.coordinates(:,1),MUA.coordinates(:,2));
        perturbationdata(ii).dC = F.C - FC_baseline(MUA.coordinates(:,1),MUA.coordinates(:,2));
    
        I = find(Table.ExpID == perturbations(ii));
        perturbationdata(ii).description = string(Table.Comments(I)); 
    
        CtrlVarInRestartFile.PlotsXaxisLabel=' ';
        CtrlVarInRestartFile.PlotsYaxisLabel=' ';
        PlotNodalBasedQuantities_JDR(ax_fig(ii),MUA.connectivity,MUA.coordinates,perturbationdata(ii).dspeed,CtrlVarInRestartFile); hold on;
        
        %shading(ax_fig(ii),"flat");
        colormap(ax_fig(ii),CM); caxis(ax_fig(ii),[-300 300]);
        %set(ax_fig(ii),'ColorScale','log');
        
        %axfig(ii).XAxis.Visible = "off";
        %axfig(ii).YAxis.Visible = "off";
        %grid(ax_fig(ii),"off"); 
        %box(ax_fig(ii),"off");
    
        title(ax_fig(ii),perturbationdata(ii).description);
        %if ii~=1 && ii~=5 && ii~=9
        %    yticklabels(ax_fig(ii),'');
        %end
        %if ii<10
        %    xticklabels(ax_fig(ii),'');
        %end
    end

    axis(ax_fig(ii),"equal","off");

end

%xlabel(tlo_fig,"psx [km]");
%ylabel(tlo_fig,"psy [km]");
cb=colorbar(ax_fig(numel(perturbationdata)),"Location","eastoutside","Ticks",[-300:100:300]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.025;
cb.FontSize=16;
cb.Label.String = "";


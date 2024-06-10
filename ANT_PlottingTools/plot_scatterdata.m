function plot_scatterdata

variable_to_plot = 'qGL'; %options: qGL, niter

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_1.csv";
UserVar.idrange = [2000,2999];

addpath("/mnt/md0/Ua/cases/ANT/");

%% read run table
RunTable = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

%% load basins
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);

%% ExpIDs
ExpID = RunTable{:,"ExpID"};
I = find(ExpID>=UserVar.idrange(1) & ExpID<=UserVar.idrange(2));

%% Gather data
for ii=1:numel(I)
    folder = UserVar.home+"ANT_"+UserVar.type+"/ANT_nsmbl_Inverse_"+ExpID(I(ii));
    restartfile = folder+"/ANT_nsmbl_Inverse_"+ExpID(I(ii))+"-RestartFile.mat";
    if exist(restartfile,"file")
        load(restartfile,"UserVarInRestartFile","CtrlVarInRestartFile","F","MUA");
        m(ii) = F.m(1);
        n(ii) = F.n(1);
        gaA(ii) = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.ga;
        gaC(ii) = CtrlVarInRestartFile.Inverse.Regularize.logC.ga;
        gsA(ii) = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.gs;
        gsC(ii) = CtrlVarInRestartFile.Inverse.Regularize.logC.gs;
        % number of iterations done
        niter(ii) = UserVarInRestartFile.Inverse.IterationsDone;
        % Obtain Ua fluxes across the grounding line (qGL) into floating areas
        %[B,GL] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
        % qGL(ii) = 0;
        % for jj=1:numel(GL)
        %     qGL(ii) = qGL(ii)+sum(GL(jj).qGL);
        % end            
        GL=FluxAcrossGroundingLine(CtrlVarInRestartFile,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho);
        qGL(ii) = sum(GL);
     else
        table_ind = I(ii);
        m(ii) = RunTable{table_ind,"m"};
        n(ii) = RunTable{table_ind,"n"};
        gaA(ii) = RunTable{table_ind,"gaA"};
        gaC(ii) = RunTable{table_ind,"gaC"};
        gsA(ii) = RunTable{table_ind,"gsA"};
        gsC(ii) = RunTable{table_ind,"gsC"};
        niter(ii) = 0;
        qGL(ii) = 0;
    end
    
    fprintf("Done %s out of %s.\n",string(ii),string(numel(I)));
end

save("scatterdata.mat","m","n","gaA","gaC","gsA","gsC","niter","qGL");

qGL(niter<100)=nan;
qGL = qGL/1e12; % convert to Gt/yr


%% Plotting
H=fig('units','inches','width',120*12/72.27,'height',80*12/72.27,'fontsize',14,'font','Helvetica');

switch variable_to_plot
    case 'qGL'
        plotdata = qGL;        
        cmin = 1000;
        cmax = 3000;
    case 'niter' 
        plotdata = niter;        
        cmin = 0;
        cmax = max(niter(:));
end

dummydata = nan*plotdata;

plotdata(plotdata<cmin)=cmin; 
plotdata(plotdata>cmax)=cmax; 

colormap('jet');
marker_size = 25;


tlo_fig = tiledlayout(6,6,"TileSpacing","compact");
for i = 1:36
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

ind = 1;
% row 1: m
scatter(ax_fig(ind),m,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"m"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),n,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaA,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaC,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsA,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsC,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;

% row 2: n
scatter(ax_fig(ind),m,n,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"n"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),n,n,marker_size,dummydata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaA,n,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaC,n,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsA,n,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsC,n,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;

% row 3: gaA
scatter(ax_fig(ind),m,gaA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"gaA"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),n,gaA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaA,gaA,marker_size,dummydata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaC,gaA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsA,gaA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsC,gaA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;


% row 4: gaC
scatter(ax_fig(ind),m,gaC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"gaC"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),n,gaC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaA,gaC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaC,gaC,marker_size,dummydata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsA,gaC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsC,gaC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;


% row 5: gsA
scatter(ax_fig(ind),m,gsA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"gsA"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),n,gsA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaA,gsA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaC,gsA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsA,gsA,marker_size,dummydata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsC,gsA,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;

% row 6: gsC
scatter(ax_fig(ind),m,gsC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"gsC"); xlabel(ax_fig(ind),"m"); 
ind = ind+1;
scatter(ax_fig(ind),n,gsC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"n"); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaA,gsC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gaA"); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gaC,gsC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gaC"); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsA,gsC,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gsA"); yticklabels(ax_fig(ind),"");
ind = ind+1;
scatter(ax_fig(ind),gsC,gsC,marker_size,dummydata,'filled'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gsC"); yticklabels(ax_fig(ind),"");
ind = ind+1;

for i = 1:36
    grid(ax_fig(i),"on");
    box(ax_fig(i),"on");
end

cb=colorbar;
cb.Layout.Tile='east';
clim([cmin cmax]);
%cb.Label.String = "Iterations";
cb.Label.String = "Grounding line flux [Gt/yr]";







function plot_InverseResults_Ensemble

variable_to_plot = 'misfit'; %options: qGL, niter, misfit

load("inversiondata.mat");

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
UserVar.cycle = 2;
UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_"+string([2 5 8])+".csv";
UserVar.idrange = [3000,3999;6000,6999;9000,9999];

%% Plotting
H=fig('units','inches','width',120*12/72.27,'height',80*12/72.27,'fontsize',14,'font','Helvetica');

switch variable_to_plot
    case 'qGL'
        for ii=1:numel(data)
            plotdata(ii,:) = data(ii).qGL(:)'/1e12;  
        end
        cmin = 1000;
        cmax = 3000;
        cbLabel = "Grounding line flux [Gt/yr]";
    case 'niter' 
        for ii=1:numel(data)
            plotdata(ii,:) = data(ii).niter(:)';  
        end       
        cmin = 0;
        cmax = max(niter(:));
        cbLabel = "Number of iterations";
    case 'misfit'
        for ii=1:numel(data)
            plotdata(ii,:) = data(ii).misfit(:)';  
        end 
        cmin = 0;
        cmax = 1000;
        cbLabel = "Misfit";
end
m = [data(:).m];
n = [data(:).n];
gaA = [data(:).gaA];
gaC = [data(:).gaC];
gsA = [data(:).gsA];
gsC = [data(:).gsC];


plotdata = plotdata(:,UserVar.cycle);
dummydata = nan*plotdata;

%plotdata(plotdata<cmin)=cmin; 
%plotdata(plotdata>cmax)=cmax; 

colormap('jet');
marker_size = 30;%;1+100*niter/(max(niter)-min(niter));

tlo_fig = tiledlayout(6,6,"TileSpacing","compact");
for i = 1:36
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

ind = 1;
% row 1: m
%scatter(ax_fig(ind),m,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"m"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
% edges = linspace(min(m),max(m),10);
% [~,~,bin] = histcounts(m,edges);
% meany = accumarray(bin(:),plotdata(:))./accumarray(bin(:),1);
% xmid = 0.5*(edges(1:end-1)+edges(2:end));
% bar(ax_fig(ind),xmid,meany);
% ylabel(ax_fig(ind),"qGL [Gt/yr]"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
% ind = ind+1;
scatter(ax_fig(ind),n,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"m"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
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
cb.Label.String = cbLabel;
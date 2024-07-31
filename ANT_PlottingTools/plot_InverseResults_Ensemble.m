function plot_InverseResults_Ensemble

variable_to_plot = 'misfit'; %options: qGL, niter, misfit

load("inversiondata.mat");

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
UserVar.cycle = 1;
UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_"+string([2 5 8])+".csv";
UserVar.idrange = [3000,3999;6000,6999;9000,9999];

% define size for markers: resize with misfit
for ii=1:numel(data)
    misfit(ii) = data(ii).misfit(UserVar.cycle);
    ind_finished(ii) = ismember(data(ii).niter(UserVar.cycle),[5000,15000]);
end 
misfit = misfit(:);
N = misfit-min(misfit);
N = N./max(N);
alphavalue = 1-N;
marker_size = 50*alphavalue+eps;

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

tlo_fig = tiledlayout(6,6,"TileSpacing","compact");
for i = 1:36
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

ind = 1;
% row 1: m
s(ind)=scatter(ax_fig(ind),m,m,marker_size,dummydata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),"m"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
% edges = linspace(min(m),max(m),10);
% [~,~,bin] = histcounts(m,edges);
% meany = accumarray(bin(:),plotdata(:))./accumarray(bin(:),1);
% xmid = 0.5*(edges(1:end-1)+edges(2:end));
% bar(ax_fig(ind),xmid,meany);
% ylabel(ax_fig(ind),"qGL [Gt/yr]"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),n,m,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaA,m,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaC,m,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsA,m,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsC,m,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;

% row 2: n
s(ind)=scatter(ax_fig(ind),m,n,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),"n"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),n,n,marker_size,dummydata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaA,n,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaC,n,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsA,n,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsC,n,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;

% row 3: gaA
s(ind)=scatter(ax_fig(ind),m,gaA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),"gaA"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),n,gaA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaA,gaA,marker_size,dummydata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaC,gaA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsA,gaA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsC,gaA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;


% row 4: gaC
s(ind)=scatter(ax_fig(ind),m,gaC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),"gaC"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),n,gaC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaA,gaC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaC,gaC,marker_size,dummydata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsA,gaC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsC,gaC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;


% row 5: gsA
s(ind)=scatter(ax_fig(ind),m,gsA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),"gsA"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),n,gsA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaA,gsA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaC,gsA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsA,gsA,marker_size,dummydata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsC,gsA,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),""); yticklabels(ax_fig(ind),"");
ind = ind+1;

% row 6: gsC
s(ind)=scatter(ax_fig(ind),m,gsC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),"gsC"); xlabel(ax_fig(ind),"m"); 
ind = ind+1;
s(ind)=scatter(ax_fig(ind),n,gsC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"n"); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaA,gsC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gaA"); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gaC,gsC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gaC"); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsA,gsC,marker_size,plotdata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gsA"); yticklabels(ax_fig(ind),"");
ind = ind+1;
s(ind)=scatter(ax_fig(ind),gsC,gsC,marker_size,dummydata,'filled','MarkerEdgeColor','none'); ylabel(ax_fig(ind),""); xlabel(ax_fig(ind),"gsC"); yticklabels(ax_fig(ind),"");
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

for ii=1:numel(s)
    s(ii).AlphaData = alphavalue;
    s(ii).MarkerFaceAlpha='flat';
end

%% histogram
figure; hold on;
bins = linspace(min(plotdata),max(plotdata),20);
h1=histogram(plotdata,bins); 
h2=histogram(plotdata(ind_finished==1),bins);
xlabel(cbLabel);
grid on; box on;

return

%% weighted histogram
nbins = 20;
minV  = min(plotdata);
maxV  = max(plotdata);
delta = (maxV-minV)/nbins;
vinterval = linspace(minV, maxV, nbins)-delta/2.0;
bincounter = zeros(nbins,1);
bincounter_finished = zeros(nbins,1);
histw = zeros(nbins, 1);
histw_finished = zeros(nbins, 1);
for i=1:length(plotdata)
    ind = find(vinterval < plotdata(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + 1./misfit(i);
      bincounter(ind) = bincounter(ind) + 1;

      histw_finished(ind) = histw_finished(ind) + ind_finished(i)*(1./misfit(i));
      bincounter_finished(ind) = bincounter_finished(ind) + 1*ind_finished(i);
    end
    
end
figure; bar(vinterval,histw./bincounter); hold on;
bar(vinterval,histw_finished./bincounter_finished);
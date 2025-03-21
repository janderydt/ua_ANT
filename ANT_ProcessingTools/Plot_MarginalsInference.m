function Plot_MarginalsInference

Klear

%% initialize UQLAB
addpath(genpath(getenv("froot_matlabfunctions")+"/../UQLab_Rel2.0.0"));

rng(1,'twister'); % set the random number generator for reproducible results
uqlab; % initialize uqlab

%% load sample from existing Bayesian Analysis
[fname,location] = uigetfile("/mnt/md0/Ua/cases/ANT/ANT_ProcessingTools/BayesianAnalysis/*.mat","MultiSelect","on");
fname = string(fname);

for ff=1:numel(fname)

    load(location+"/"+fname(ff));
    
    if ~exist("myPosteriorDist","var")
        error(fname(FF)+":  myPoseriorDist does not exist.");
    end

    data(ff).myPosteriorDist = myPosteriorDist;
    data(ff).l = UserVar.l;
    data(ff).fac_su = UserVar.fac_su;
    data(ff).mpsrf = myBayesianAnalysis.Results.PostProc.MPSRF;
    tmp = split(UserVar.fname,"_");
    data(ff).nametag = erase(tmp(end),".mat");

end

% reorder data
[~,I] = sort([data(:).fac_su],'ascend');
data = data(I);
fac_su = unique([data(:).fac_su]);
for ii=1:numel(fac_su)
    I=find([data(:).fac_su]==fac_su(ii));
    % reorder data
    [~,I2] = sort([data(I).l],'ascend');
    data(I) = data(I(I2));
end

l = unique([data(:).l]);
ndiscrete = numel(data(1).myPosteriorDist);
nmarginals = numel(data(1).myPosteriorDist(1).dist.Marginals);

linewidth = linspace(2,0.75,numel(fac_su));
linestyle = ["-","--",":","-.","-","--",":","-."];
CM = flipud(cmocean('speed',numel(l)+1)); CM=CM(1:end-1,:);
if ndiscrete==2
    CM = [CM; flipud(cmocean('matter',numel(l)+1))]; CM=CM(1:end-1,:);
end

for ff=1:numel(data)
    data(ff).linewidthind = find(fac_su==data(ff).fac_su);
    data(ff).linestyleind = find(fac_su==data(ff).fac_su);
    data(ff).colorind = find(l==data(ff).l);
end


%% PLOT
H=fig('units','inches','width',120*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');

tlo_fig = tiledlayout(1,nmarginals,"TileSpacing","compact");
for i = 1:nmarginals
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

nn=1;
for jj=1:nmarginals
    for ii=1:numel(data)
        for kk=1:ndiscrete
            if ~isempty(data(ii).myPosteriorDist(kk).dist)
                X = linspace(data(ii).myPosteriorDist(kk).dist.Marginals(jj).Bounds(1),data(ii).myPosteriorDist(kk).dist.Marginals(jj).Bounds(2),50);            
                f = uq_all_pdf(X(:),data(ii).myPosteriorDist(kk).dist.Marginals(jj))*data(ii).myPosteriorDist(kk).frac;
            else
                X = [nan nan];
                f = [nan nan];
            end
            if jj==1 && kk==1
                g(nn) = plot(ax_fig(jj),X,f,linestyle(data(ii).linestyleind),'color',CM((kk-1)*numel(fac_su)+data(ii).colorind,:),...
                    'LineWidth',linewidth(data(ii).linewidthind),'markersize',4);
                label(nn) = "("+data(ii).nametag+") fac="+string(data(ii).fac_su)+", l="+string(data(ii).l)+", MPSRF="+sprintf("%.1f",data(ii).mpsrf);
                nn=nn+1;
            else
                plot(ax_fig(jj),X,f,linestyle(data(ii).linestyleind),'color',CM((kk-1)*numel(fac_su)+data(ii).colorind,:),...
                    'LineWidth',linewidth(data(ii).linewidthind),'markersize',4);
            end
        end
    end
end

for ii=1:nmarginals
    grid(ax_fig(ii),'on'); 
    box(ax_fig(ii),'on');
    xlabel(ax_fig(ii),myBayesianAnalysis.PriorDist.Marginals(ii).Name);
    if ii==1
        leg=legend(g,label,'Orientation','Horizontal','NumColumns',5);
        leg.Layout.Tile = 'north';
    end
end
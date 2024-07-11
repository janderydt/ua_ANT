function plot_scatterdata

variable_to_plot = 'misfit'; %options: qGL, niter, misfit

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
UserVar.cycle = 1;
UserVar.Table = UserVar.home+"ANT_"+UserVar.type+["/RunTable_ARCHER2_2.csv"];%,"/RunTable_ARCHER2_5.csv"];
UserVar.idrange = [3000,3999];

addpath("/mnt/md0/Ua/cases/ANT/");

%% load basins
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);

kk=0;

for tt=1:numel(UserVar.Table)

    %% read run table
    RunTable = ANT_ReadWritetable(UserVar,UserVar.Table(tt),[],'read');
    
    %% ExpIDs
    ExpID = RunTable{:,"ExpID"};
    Ind = find(ExpID>=UserVar.idrange(tt,1) & ExpID<=UserVar.idrange(tt,2));
    
    %% Gather data
    for ii=1:numel(Ind)
        folder = UserVar.home+"/ANT_"+UserVar.type+"/cases/ANT_nsmbl_Inverse_"+ExpID(Ind(ii));
        restartfile = folder+"/ANT_nsmbl_Inverse_"+ExpID(Ind(ii))+"-RestartFile_InverseCycle"+...
            string(UserVar.cycle)+".mat";
        if exist(restartfile,"file")
            load(restartfile,"UserVarInRestartFile","CtrlVarInRestartFile","F","MUA","InvFinalValues");
            if UserVarInRestartFile.Inverse.IterationsDone == 15000
                m(kk+ii) = F.m(1);
                n(kk+ii) = F.n(1);
                gaA(kk+ii) = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.ga;
                gaC(kk+ii) = CtrlVarInRestartFile.Inverse.Regularize.logC.ga;
                gsA(kk+ii) = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.gs;
                gsC(kk+ii) = CtrlVarInRestartFile.Inverse.Regularize.logC.gs;
                % number of iterations done
                niter(kk+ii) = UserVarInRestartFile.Inverse.IterationsDone;
                fprintf("(%s/%s) ExpID %s: done %s iterations.\n",string(ii),string(numel(Ind)),string(UserVarInRestartFile.ExpID),string(niter(ii)));
                % Obtain Ua fluxes across the grounding line (qGL) into floating areas
                %[B,GL] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
                % qGL(ii) = 0;
                % for jj=1:numel(GL)
                %     qGL(ii) = qGL(ii)+sum(GL(jj).qGL);
                % end            
                GL=FluxAcrossGroundingLine(CtrlVarInRestartFile,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho);
                qGL(kk+ii) = sum(GL);
                I(kk+ii) = InvFinalValues.I; % calculated as 
            else
                table_ind = Ind(ii);
                m(kk+ii) = RunTable{table_ind,"m"};
                n(kk+ii) = RunTable{table_ind,"n"};
                gaA(kk+ii) = RunTable{table_ind,"gaA"};
                gaC(kk+ii) = RunTable{table_ind,"gaC"};
                gsA(kk+ii) = RunTable{table_ind,"gsA"};
                gsC(kk+ii) = RunTable{table_ind,"gsC"};
                niter(kk+ii) = 0;
                qGL(kk+ii) = nan;
                I(kk+ii) = nan;
            end
         else
            table_ind = Ind(ii);
            m(kk+ii) = RunTable{table_ind,"m"};
            n(kk+ii) = RunTable{table_ind,"n"};
            gaA(kk+ii) = RunTable{table_ind,"gaA"};
            gaC(kk+ii) = RunTable{table_ind,"gaC"};
            gsA(kk+ii) = RunTable{table_ind,"gsA"};
            gsC(kk+ii) = RunTable{table_ind,"gsC"};
            niter(kk+ii) = 0;
            qGL(kk+ii) = nan;
            I(kk+ii) = nan;
        end
              
        %fprintf("Done %s out of %s.\n",string(ii),string(numel(I)));
    end
    kk = numel(I);
end

save("scatterdata.mat","m","n","gaA","gaC","gsA","gsC","niter","qGL","I");

%qGL(niter<100)=nan;
qGL = qGL/1e12; % convert to Gt/yr


%% Plotting
H=fig('units','inches','width',120*12/72.27,'height',80*12/72.27,'fontsize',14,'font','Helvetica');

switch variable_to_plot
    case 'qGL'
        plotdata = qGL;        
        cmin = 1000;
        cmax = 3000;
        cbLabel = "Grounding line flux [Gt/yr]";
    case 'niter' 
        plotdata = niter;        
        cmin = 0;
        cmax = max(niter(:));
        cbLabel = "Number of iterations";
    case 'misfit'
        plotdata = I;
        cmin = 0;
        cmax = 1000;
        cbLabel = "Misfit";
end

dummydata = nan*plotdata;

%plotdata(plotdata<cmin)=cmin; 
%plotdata(plotdata>cmax)=cmax; 

colormap('jet');
marker_size = 1+100*niter/(max(niter)-min(niter));


tlo_fig = tiledlayout(6,6,"TileSpacing","compact");
for i = 1:36
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

ind = 1;
% row 1: m
%scatter(ax_fig(ind),m,m,marker_size,plotdata,'filled'); ylabel(ax_fig(ind),"m"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
edges = linspace(min(m),max(m),10);
[~,~,bin] = histcounts(m,edges);
meany = accumarray(bin(:),plotdata(:))./accumarray(bin(:),1);
xmid = 0.5*(edges(1:end-1)+edges(2:end));
bar(ax_fig(ind),xmid,meany);
ylabel(ax_fig(ind),"qGL [Gt/yr]"); xlabel(ax_fig(ind),""); xticklabels(ax_fig(ind),"");
ind = ind+1;
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







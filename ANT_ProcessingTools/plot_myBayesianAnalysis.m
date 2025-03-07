function plot_myBayesianAnalysis(myBayesianAnalysis)

if nargin==0
    [file,location] = uigetfile("/mnt/md0/Ua/cases/ANT/ANT_ProcessingTools/BayesianAnalysis/*.mat");
    load(location+"/"+file);
end

mySolver = myBayesianAnalysis.Options.Solver;
myPriorDist = myBayesianAnalysis.PriorDist;
myResults = myBayesianAnalysis.Results;

%% gather some data
MAP=myBayesianAnalysis.Results.PostProc.PointEstimate.X{:};

nDiscrete = 0;
for ii=1:numel(myPriorDist.Marginals)

    % Sample equilibrium distribution from Markov Chains; we use the final 
    % 20% of the Markov Chains, i.e. discard the first 80% as burn-in
    Post(ii).Sample = [];
    for jj = 1:mySolver.MCMC.NChains
        Post(ii).Sample = [Post(ii).Sample; ...
            myResults.Sample(mySolver.MCMC.Steps-round(0.20*mySolver.MCMC.Steps):end,ii,jj)];
    end

    % Calculate spread of sample, discarding outliers that are outside the
    % 1% and 99% percentiles. This is used for plotting purposes ONLY
    Outliers = isoutlier(Post(ii).Sample,"percentiles",[1 99]);
    Post(ii).Sample_NoOutliers = Post(ii).Sample;
    Post(ii).Sample_NoOutliers(find(Outliers))=[];
    Post(ii).min = min(Post(ii).Sample_NoOutliers);
    Post(ii).max = max(Post(ii).Sample_NoOutliers);
    dsample = Post(ii).max-Post(ii).min;
    Post(ii).min = Post(ii).min-0.4*dsample;
    Post(ii).max = Post(ii).max+0.4*dsample;
    Post(ii).Name = myPriorDist.Marginals(ii).Name;

    % Identify discrete variables
    if Post(ii).Name=="Cycle"
        D1 = "no spinup";
        D2 = "3yr spinup";
        Post(ii).isDiscrete = 1;
        nDiscrete = nDiscrete+1;
    elseif Post(ii).Name=="SlidingLaw"
        D1 = "Weertman";
        D2 = "Weertman+Coulomb";
        Post(ii).isDiscrete = 1;
        nDiscrete = nDiscrete+1;
    else
        Post(ii).isDiscrete = 0;
    end
end

%% plot results
if nDiscrete <= 1
    CMtmp = cmocean('curl',64);
    CM(1).map = CMtmp(5:30,:);
    CM(2).map = CMtmp(35:60,:);
else
    error("Implement additional colormaps");
end

H=fig('units','inches','width',120*12/72.27,'height',80*12/72.27,'fontsize',14,'font','Helvetica');

tlo_fig = tiledlayout(numel(Post),numel(Post),"TileSpacing","compact");
for i = 1:numel(Post)^2
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

for iy=1:numel(Post)

    for ix=1:numel(Post)
        
        Ind = (iy-1)*numel(Post)+ix;
        edges=[]; xmid=[]; ymid=[];

        if ix==iy
            % histograms
            if Post(ix).isDiscrete
                edges = linspace(0,1,3);
                xmid = 0.5*(edges(1:end-1)+edges(2:end));
                [n,~,~] = histcounts(Post(ix).Sample,edges); 
                b=bar(ax_fig(Ind),xmid,n/sum(n));
                b.FaceColor='flat';
                b.CData(1,:) = CM(1).map(16,:);
                b.CData(2,:) = CM(2).map(16,:);
                MAPx = round(MAP(ix)/0.5)*0.5-0.25;
            else
                edges = linspace(Post(ix).min,Post(ix).max,30);
                xmid = 0.5*(edges(1:end-1)+edges(2:end));
                if nDiscrete == 1
                    % split data into two parts
                    Ind1 = find(Post(end).Sample<=0.5);
                    Ind2 = find(Post(end).Sample>0.5);
                    [n1,~,~] = histcounts(Post(ix).Sample(Ind1),edges); 
                    [n2,~,~] = histcounts(Post(ix).Sample(Ind2),edges); 
                    b = bar(ax_fig(Ind),xmid,[n1(:)'; n2(:)']/sum([n1,n2]),'stacked');
                    b(1).FaceColor = CM(1).map(16,:);
                    b(2).FaceColor = CM(2).map(16,:);                   
                elseif nDiscrete > 1
                    error("plotting so far only implemented for 1 discrete variable.");
                else
                    bar(ax_fig(Ind),xmid,n/sum(n),'FaceColor',CM(1).map(end,:));     
                end
                MAPx = MAP(ix);
            end           
            plot(ax_fig(Ind),[MAPx MAPx],[0 1],'-','Color',[0.929 0.694 0.125],'linewidth',1.5);
        else
            % scatter of data
            if Post(ix).isDiscrete
                edges{2} = linspace(0,1,3);
                MAPx = round(MAP(ix)/0.5)*0.5-0.25;
            else
                edges{2} = linspace(Post(ix).min,Post(ix).max,50); % edges in second (x) dimension
                MAPx = MAP(ix);
            end
            if Post(iy).isDiscrete
                edges{1} = linspace(0,1,3);
                MAPy = round(MAP(iy)/0.5)*0.5-0.25;
            else
                edges{1} = linspace(Post(iy).min,Post(iy).max,50); % edges in second (x) dimension
                MAPy = MAP(iy);
            end            
            xmid = 0.5*(edges{2}(1:end-1)+edges{2}(2:end));
            ymid = 0.5*(edges{1}(1:end-1)+edges{1}(2:end));
            [Xedge,Yedge] = ndgrid(edges{2},edges{1});
            if nDiscrete == 1
                % split data into two parts
                Ind1 = find(Post(end).Sample<=0.5);
                Ind2 = find(Post(end).Sample>0.5);
                values1 = hist3([Post(iy).Sample(Ind1) Post(ix).Sample(Ind1)],'Edges',edges);
                values2 = hist3([Post(iy).Sample(Ind2) Post(ix).Sample(Ind2)],'Edges',edges);
                sum_values = sum([values1,values2],"all");
                values1 = -values1/sum_values; 
                values2 = values2/sum_values; 
                values1(values1==0)=nan; values2(values2==0)=nan;
                pcolor(ax_fig(Ind),Xedge,Yedge,values1');
                pcolor(ax_fig(Ind),Xedge,Yedge,values2'); 
                shading(ax_fig(Ind),'flat'); 
                colormap(ax_fig(Ind),[CM(1).map; CM(2).map]);
            else
                %scatter(ax_fig(Ind),Post(jj).Sample,Post(ii).Sample,15,...
                    %'MarkerFaceColor',[1 0.5 0.5],'MarkerFaceAlpha',.4,'MarkerEdgeColor','none');
                values = values./sum(values,"all");
                values(values==0)=nan;      
                pcolor(ax_fig(Ind),Xedge,Yedge,values'); shading(ax_fig(Ind),'flat'); 
                colormap(ax_fig(Ind),CM(1).map);
            end
            caxis(ax_fig(Ind),[-0.01 0.01]);
            plot(ax_fig(Ind),MAPx,MAPy,'x','Color',[0.929 0.694 0.125],'MarkerFaceColor',[0.929 0.694 0.125],'MarkerSize',10,'linewidth',3);
        end
        if ix==1
            if Post(iy).isDiscrete
                yticks(ax_fig(Ind),ymid);
                yticklabels(ax_fig(Ind),[D1,D2]);
                ylabel(ax_fig(Ind),[]);
            else
                ylabel(ax_fig(Ind),Post(iy).Name);
            end
        else
            yticklabels(ax_fig(Ind),[]);
        end
        if iy==numel(Post)
            if Post(ix).isDiscrete
                xticks(ax_fig(Ind),xmid);
                xticklabels(ax_fig(Ind),[D1,D2]);
                xlabel(ax_fig(Ind),[]);
            else
                xlabel(ax_fig(Ind),Post(ix).Name);
            end
        else
            xticklabels(ax_fig(Ind),[]); 
        end
        xlim(ax_fig(Ind),[Post(ix).min Post(ix).max]);
        if ix~=iy
            ylim(ax_fig(Ind),[Post(iy).min Post(iy).max]);
        else
            if Post(ix).isDiscrete
                ylim(ax_fig(Ind),[0 1]);
            else
                ylim(ax_fig(Ind),[0 0.3]);
            end
        end
        grid(ax_fig(Ind),'on'); box(ax_fig(Ind),'on');
    end
end

return

for ii=1:numel(years)
    years_tmp = split(years(ii),"-");
    yr1 = years_tmp(1);
    yr2 = years_tmp(2);
    load("Delta_u_AMUND_Weertman_"+yr1+"-"+yr2+".mat");
    MUA = MUA_yr2; GF = GF_yr2;

    % nearest predictor and target
    X = [log10(gaA(:)) log10(gaC(:)) log10(gsA(:)) log10(gsC(:)) m(:) n(:)];
    if cycle==2 % add dhdt_err to predictors
        X = [X log10(dhdt_err(:))];
    end
    X = double(X);
    Ind_MAP=knnsearch(X,MAP);
    deltau_ua = Delta_u.Calv_dh.map(Ind_MAP,:,1);

    % emulated values at MAP
    % -- load emulator
    if dataformat == "LOGu"
        yearstr = "LOGu"+years(ii); 
    else
        yearstr = years(ii);
    end
    switch NN(ii)
        case "RNN"
            net_tmp=importNetworkFromTensorFlow("./RNN/TF_files/tuned_model_SVD-nodata_Calv_dh_"+...
                yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii)));
            Net = initialize(net_tmp);
            load("./RNN/mat_files/data_SVD-nodata_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat",...
                "X_train_C","X_train_S","T_train_C","T_train_S");
            load("./RNN/mat_files/SVD-nodata_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat",...
                "B_trunc","T_mean","T_reproj");        
        case "FNN"
            load("./FNN/mat_files/FNN_trainscg_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat",...
                "Net_opt");
            load("./FNN/mat_files/SVD_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat",...
                "B_trunc","T_mean","T_reproj");
            Net = Net_opt.trained;
            X_train_C = Net_opt.X_train_C(:)';
            X_train_S = Net_opt.X_train_S(:)';
            T_train_C = Net_opt.T_train_C(:)';
            T_train_S = Net_opt.T_train_S(:)';
    end
    % -- apply normalization to parameters before feeding into emulator
    predictors = (MAP-X_train_C)./X_train_S;
    % -- evaluate emulator 
    switch NN(ii)
        case "RNN"
            modelRun = double(predict(Net,predictors)); 
        case "FNN"
            modelRun = double(Net(predictors')');
    end
    % -- Undo normalization of the output and project to nodal basis
    modelRun = modelRun.*T_train_S+T_train_C;
    deltau_emul = modelRun*B_trunc+T_mean;
    
    % measurements
    out = loadvelocitydata(dataformat(ii),years(ii),0);
    deltau_meas = out.("yr"+yr1+"_yr"+yr2).du;
    
    % measurements in SVD basis
    deltau_meas_SVD = deltau_meas(:)'*T_reproj'*B_trunc;
    
    if only_grounded_ice(ii)
        deltau_ua(GF.node<0.5) = 0;
        deltau_emul(GF.node<0.5) = 0;
        deltau_meas(GF.node<0.5) = 0;
        deltau_meas_SVD(GF.node<0.5) = 0;
    end
    
    figure; tlo = tiledlayout(1,4,"TileSpacing","tight");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_ua(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([0 750]); colorbar("off");
    axis tight; xlim([-19 -9]*1e5); axis off;
    title("Ua closest to MAP");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_emul(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([0 750]); colorbar("off")
    axis tight;   xlim([-19 -9]*1e5); axis off;
    title("Emulator at MAP");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_meas(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([0 750]); colorbar("off")
    axis tight;  xlim([-19 -9]*1e5); axis off;
    title("Measurements");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_meas_SVD(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([0 750]); colorbar("off"); colormap(cmocean('deep'));
    axis tight;  xlim([-19 -9]*1e5); axis off;
    title(["Measurements ";"in SVD basis"]);
    
    cb = colorbar;
    cb.Layout.Tile = "east";
    cb.Label.String = "\Deltau [m/yr]";

    title(tlo,["Comparison of MAP with measurements "+yearstr;...
        "MAP: "+...
        sprintf('log_{10}(gaA)=%.1f, log_{10}(gaC)=%.1f, log_{10}(gsA)=%.1f, log_{10}(gsC)=%.1f, m=%.1f, n=%.1f',MAP(1:6))]);
    
    figure; tlo = tiledlayout(1,4,"TileSpacing","tight");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_ua(:)-deltau_meas(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([-250 250]); colorbar("off")
    axis tight;   xlim([-19 -9]*1e5); axis off;
    title("Ua - Meas");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_ua(:)-deltau_meas_SVD(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([-250 250]); colorbar("off")
    axis tight;   xlim([-19 -9]*1e5); axis off;
    title("Ua - Meas SVD");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_emul(:)-deltau_meas(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([-250 250]); colorbar("off")
    axis tight;   xlim([-19 -9]*1e5); axis off;
    title("Emul - Meas");
    
    nexttile; PlotMeshScalarVariable([],MUA_yr2,deltau_emul(:)-deltau_meas_SVD(:)); hold on;
    PlotGroundingLines([],MUA,GF,[],[],[],'-k','linewidth',1);
    caxis([-250 250]); colorbar("off"); colormap(cmocean('balance'));
    axis tight;   xlim([-19 -9]*1e5); axis off;
    title("Emul - Meas SVD");
    
    cb = colorbar;
    cb.Layout.Tile = "east";
    cb.Label.String = "\Deltau [m/yr]";
    %subplot(1,2,1); hold on; PlotMeshScalarVariable([],MUA_yr2,deltau_ua(:)-deltau_meas(:)); title("Ua nearest to MAP - Measured");
    %subplot(1,2,2); hold on; PlotMeshScalarVariable([],MUA_yr2,deltau_emul(:)-deltau_meas(:)); title("Emulator at MAP - Measured");

    title(tlo,["Comparison of MAP with measurements "+yearstr;...
        "MAP: "+...
        sprintf('log_{10}(gaA)=%.1f, log_{10}(gaC)=%.1f, log_{10}(gsA)=%.1f, log_{10}(gsC)=%.1f, m=%.1f, n=%.1f',MAP(1:6))]);
end
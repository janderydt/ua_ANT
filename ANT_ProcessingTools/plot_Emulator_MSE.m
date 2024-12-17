function plot_Emulator_MSE

plotRNN=0;
plotFNN=1;

perturbation='Calv_dh';
startyear="2000";
targetyear=""; % this can be an empty string to plot results for the u emulator, rather than the du emulator
slidinglaw="Weertman";
cycle=2;
n_test=1; % index of test data to plot (integer between 1 and size of test dataset)

if targetyear ~= ""
    years=startyear+"-"+targetyear;
    datafile = "Delta_u_AS_"+slidinglaw+"_"+years+".mat";
    load(datafile);
    
    MUA = MUA_yr2; GF = GF_yr2;
    T = Delta_u.(perturbation).map(:,:,cycle);
else
    years = startyear;
    datafile = "u_AS_"+perturbation+"_cycle"+string(cycle)+"_"+slidinglaw+"_2000_2009_2014_2018.mat";
    load(datafile);
    T = speed.("yr"+string(years));
    MUA = MUA.("yr"+string(years)); GF = GF.("yr"+string(years));
end

Ind_toremove=find(gaC>50);
% check for nans and inf
if any(isnan(T))
    warning("Removing nan from T");
    [rows,~]=find(isnan(T));
    rows = unique(rows);
    T(rows,:)=[];
end
T(Ind_toremove,:)=[];
T_mean=mean(T,1);
T = T-repmat(T_mean(:)',size(T,1),1);

pct=string([950 960 970 980 990 995]);
%pct = string([995]);

%% RNN Tensorflow
if plotRNN

    for ii=1:numel(pct)
    
        net=importNetworkFromTensorFlow("./RNN/TF_files/tuned_model_"+perturbation+"_"+years+"_"+slidinglaw+...
            "_cycle"+string(cycle)+"_N0k"+pct(ii));
        netUpdated=initialize(net);
        load("./RNN/mat_files/data_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+string(cycle)+"_N0k"+pct(ii)+".mat");
        load("./RNN/mat_files/SVD_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+string(cycle)+"_N0k"+pct(ii)+".mat");
        num_train = size(X_train,1);
        num_val = size(X_val,1);
        num_test = size(X_test,1);
        
        % emulate predictors
        Y_train=double(predict(netUpdated,X_train));
        Y_val=double(predict(netUpdated,X_val));
        Y_test=double(predict(netUpdated,X_test));
    
        % for m=1:size(T_full,2)
        %     figure; hold on;
        %     plot(T_train(:,m),Y_train(:,m),'ok');
        %     plot(T_val(:,m),Y_val(:,m),'xk');
        %     plot(T_test(:,m),Y_test(:,m),'dk');
        %     plot(T_test(10,m),Y_test(10,m),'dr');
        %     %plotregression(T_full,Y(:,m));
        %     title("mode "+num2str(m));
        % end
    
        % reconstruct spatial maps of emulated targets
        Y_test = Y_test.*repmat(T_train_S,num_test,1)+...
            repmat(T_train_C,num_test,1); % undo normalization
        Y_reproj = Y_test*B_trunc;% reproject onto nodal basis
        
        % ua outputs
        Ind_test_orig = seq(num_train+num_val+1:end);
        Ua_orig = T(Ind_test_orig,:);
        Ua_proj = (Ua_orig-repmat(T_mean,num_test,1))*T_reproj';
    
        % reconstruct original target data
        T_test = T_test.*repmat(T_train_S,num_test,1)+...
            repmat(T_train_C,num_test,1); % undo normalization
        T_test = T_test*B_trunc;% reproject onto nodal basis
        
        % calc mse of test data
        % in nodal basis
        mse(ii,:)=1/num_test*sum((Y_reproj-Ua_orig).^2,1);
        % in truncated svd basis
        mse_RNN=1/num_test*sum((Y_test-Ua_proj).^2,1);
        save("./RNN/mat_files/MSE_RNN_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+"_N0k"+pct(ii),"mse_RNN");

        % plot some maps for particular test case
        predictorvalues = X_test(n_test,:).*X_train_S + X_train_C;
        Ua_orig = Ua_orig(n_test,:);
        Y_reproj = Y_reproj(n_test,:);
        T_test = T_test(n_test,:);
        plot_comparison_emulator_ua(T_mean,Ua_orig,T_test,Y_reproj,MUA,GF,pct(ii),predictorvalues);
    
    end

    % plot mse for different svd trunctions
    plot_mse(mse,pct,MUA,GF);

end

%% FNN Matlab
if plotFNN
    for ii=1:numel(pct)
    
        tmp=load("./FNN/mat_files/FNN_trainscg_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+"_N0k"+pct(ii)+".mat");
        load("./FNN/mat_files/SVD_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+"_N0k"+pct(ii)+".mat","B_trunc","T_reproj","seq");
        net=tmp.Net_opt;
        
        X_train = net.X_train';
        X_val = net.X_val';
        X_test = net.X_test';
        num_train = size(X_train,1);
        num_val = size(X_val,1);
        num_test = size(X_test,1);
        
        % emulate predictors
        Y_train=double(net.trained(X_train')');
        Y_val=double(net.trained(X_val')');
        Y_test=double(net.trained(X_test')');
    
        % for m=1:size(T_full,2)
        %     figure; hold on;
        %     plot(T_train(:,m),Y_train(:,m),'ok');
        %     plot(T_val(:,m),Y_val(:,m),'xk');
        %     plot(T_test(:,m),Y_test(:,m),'dk');
        %     plot(T_test(10,m),Y_test(10,m),'dr');
        %     %plotregression(T_full,Y(:,m));
        %     title("mode "+num2str(m));
        % end
    
        % reconstruct spatial maps of emulated targets
        Y_test = Y_test.*repmat(net.T_train_S',num_test,1)+...
            repmat(net.T_train_C',num_test,1); % undo normalization
        Y_reproj = Y_test*B_trunc;% reproject onto nodal basis
        
        % ua output
        Ind_test_orig = seq(num_train+num_val+1:end);
        Ua_orig = T(Ind_test_orig,:);
        Ua_proj = (Ua_orig-repmat(T_mean,num_test,1))*T_reproj';
    
        % reconstruct original target data
        T_test = net.T_test'.*repmat(net.T_train_S',num_test,1)+...
            repmat(net.T_train_C',num_test,1); % undo normalization
        T_test = T_test*B_trunc;% reproject onto nodal basis
        
        % calc mse of test data
        % in nodal basis
        mse(ii,:)=1/num_test*sum((Y_reproj-Ua_orig).^2,1);
        % in truncated svd basis
        mse_FNN=1/num_test*sum((Y_test-Ua_proj).^2,1);
        save("./FNN/mat_files/MSE_FNN_trainscg_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+"_N0k"+pct(ii),"mse_FNN");

        % plot some maps for particular test case
        predictorvalues = X_test(n_test,:).*net.X_train_S' + net.X_train_C';
        Ua_orig = Ua_orig(n_test,:);
        Y_reproj = Y_reproj(n_test,:);
        T_test = T_test(n_test,:);
        %plot_comparison_emulator_ua(T_mean,Ua_orig,T_test,Y_reproj,MUA,GF,pct(ii),predictorvalues);
    
    end

    % plot mse for different svd trunctions
    plot_mse(mse,pct,MUA,GF);

end

end

function plot_comparison_emulator_ua(T_mean,T_orig,Y_reproj_orig,Y_reproj,MUA,GF,pct,predictorvalues)

CM = othercolor('Reds8',15);
CM2 = flipdim(othercolor('RdBu11',15),1);
CtrlVar=Ua2D_DefaultParameters;
CtrlVar.PlotXYscale=1e3;
CtrlVar.PlotsXaxisLabel='' ; 
CtrlVar.PlotsYaxisLabel='' ;

H=fig('units','inches','width',120*12/72.27,'height',60*12/72.27,'fontsize',14,'font','Helvetica');

tlo=tiledlayout(1,3,"TileSpacing","tight");    

ax(1)=nexttile; hold on;
PlotMeshScalarVariable(CtrlVar,MUA,T_orig'+T_mean');
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);
colormap(ax(1),CM); 
cb=colorbar; cb.Label.String='\Delta u [m/yr]';
caxis(ax(1),[0 2000]); 
axis equal; axis tight; 
grid(ax(1),"off"); box(ax(1),"off"); axis(ax(1),"off");
yticklabels([]); xticklabels([]); 
title('Ua');

ax(2)=nexttile; hold on;
PlotMeshScalarVariable(CtrlVar,MUA,Y_reproj_orig'-T_orig');
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);
colormap(ax(2),CM2); colorbar('off')
caxis(ax(2),[-300 300]); axis equal; axis tight;
grid(ax(2),"off"); box(ax(2),"off"); axis(ax(2),"off");
yticklabels([]); xticklabels([]);
title("Truncated SVD ("+pct+"%) - Ua");

ax(3)=nexttile; hold on;
PlotMeshScalarVariable(CtrlVar,MUA,Y_reproj'-T_orig');
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);
colormap(ax(3),CM2); cb2=colorbar; cb2.Label.String='\Delta u [m/yr]';
caxis(ax(3),[-300 300]); axis equal; axis tight;
grid(ax(3),"off"); box(ax(3),"off"); axis(ax(3),"off");
yticklabels([]); xticklabels([]);
title("Emulator - Ua");

%xlabel(tlo,'psx [km]');
%ylabel(tlo,'psy [km]');

if numel(predictorvalues)==6
    title(tlo,sprintf('log10(gaA)=%3.1f, log10(gaC)=%3.1f, log10(gsA)=%3.1f, log10(gsC)=%3.1f, m=%2.2f, n=%2.2f',...
        predictorvalues));
else
    title(tlo,sprintf(['log10(gaA)=%3.1f, log10(gaC)=%3.1f, log10(gsA)=%3.1f, log10(gsC)=%3.1f, m=%2.2f, n=%2.2f, log10(dhdterr)=%3.1f'],...
        predictorvalues));
end

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Deltau_example";
print(H,fname,"-dpng","-r400");

end

function plot_mse(mse,pct,MUA,GF)

H=fig('units','inches','width',120*12/72.27,'height',30*12/72.27,'fontsize',14,'font','Helvetica');

tlo=tiledlayout(1,numel(pct),"TileSpacing","tight");    

CM = flipdim(othercolor('RdYlBu7',15),1);
CM2 = flipdim(othercolor('RdBu11',15),1);
CtrlVar=Ua2D_DefaultParameters;
CtrlVar.PlotXYscale=1e3;
CtrlVar.PlotsXaxisLabel='' ; 
CtrlVar.PlotsYaxisLabel='' ;

rmse = sqrt(mse);

for ii=1:numel(pct)
    ax(ii)=nexttile; hold on;
    
    if ii==1
        PlotMeshScalarVariable(CtrlVar,MUA,log10(rmse(ii,:)'));
        caxis(ax(ii),[0 4]); colormap(ax(ii),CM); 
        cb=colorbar; cb.Label.String='log_{10}(rmse) between emulator and Ua [m^2/yr^2]';
    else
        PlotMeshScalarVariable(CtrlVar,MUA,log10(rmse(ii,:)')-log10(rmse(1,:)'));
        caxis(ax(ii),[-1 1]); colormap(ax(ii),CM2); 
        if ii==numel(pct)
            cb=colorbar; cb.Label.String='log_{10}(rmse_{pct}) - log_{10}(rmse_{95})';
        else
            colorbar("off");
        end
    end
    PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
    plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);        
            
    axis equal; axis tight; 
    grid(ax(ii),"off"); box(ax(ii),"off"); axis(ax(ii),"off");
    yticklabels([]); xticklabels([]); 
    title(pct(ii)+"%");
end

end
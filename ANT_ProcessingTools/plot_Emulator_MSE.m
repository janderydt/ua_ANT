function plot_Emulator_MSE

plotRNN=1;
plotFNN=0;

domain="AMUND";
perturbation="Calv_dh";
trainingdataformat="u"; % u or LOGu
startyear="2000";
targetyear="2020"; % this can be an empty string to plot results for the u emulator, rather than the du emulator
slidinglaw="Umbi";
cycle=2;
only_grounded_ice=1;
include_measurements = 0;
n_test=[];%[1:20]; % indices of test data to plot (integers between 1 and size of test dataset)

if targetyear ~= ""
    years=startyear+"-"+targetyear;
    datafile = "Delta_u_"+domain+"_"+slidinglaw+"_"+years+".mat";
    load(datafile);    
    MUA = MUA_yr2; GF = GF_yr2;
    T = Delta_u.(perturbation).map(:,:,cycle);
else
    years = trainingdataformat+startyear;
    datafile = "u_"+domain+"_"+perturbation+"_cycle"+string(cycle)+"_"+slidinglaw+"_"+startyear+".mat";
    load(datafile);
    T = speed.("yr"+string(startyear));
    MUA = MUA.("yr"+string(startyear)); GF = GF.("yr"+string(startyear));
    if trainingdataformat=="LOGu"
        T = log10(T);
    end 
end

% check for nans and inf
if any(isnan(T))
    warning("Removing nan from T");
    [rows,~]=find(isnan(T));
    rows = unique(rows);
    T(rows,:)=[];
end
T_mean = mean(T,1);
T = T-repmat(T_mean(:)',size(T,1),1);

% remove floating ice?
if only_grounded_ice
    % remove floating nodes
    load("Delta_u_"+domain+"_"+slidinglaw+"_"+startyear+"-"+targetyear+".mat","MUA_yr2","GF_yr2");
    MUA = MUA_yr2; GF = GF_yr2;
    Nodes_floating = find(GF.node<0.5);
    T(:,Nodes_floating) = 0;
    T_mean(Nodes_floating) = 0;  
end

%% RNN Tensorflow
if plotRNN
        
    TF_dir = dir("./RNN/TF_files/tuned_model_"+domain+"_"+perturbation+"_"+years+...
        "_"+slidinglaw+"_cycle"+string(cycle)+"_floatingice"+string(1-only_grounded_ice)+...
        "_includemeasurements"+string(include_measurements)+"*");
    net=importNetworkFromTensorFlow(TF_dir.folder+"/"+TF_dir.name);
    netUpdated=initialize(net);

    data_file = dir("./RNN/mat_files/data_"+domain+"_"+perturbation+"_"+years+...
        "_"+slidinglaw+"_cycle"+string(cycle)+"_floatingice"+string(1-only_grounded_ice)+...
        "_includemeasurements"+string(include_measurements)+"*.mat");
    load(data_file.folder+"/"+data_file.name);

    SVD_file = dir("./RNN/mat_files/SVD_"+domain+"_"+perturbation+"_"+years+...
        "_"+slidinglaw+"_cycle"+string(cycle)+"_floatingice"+string(1-only_grounded_ice)+...
        "_includemeasurements"+string(include_measurements)+"*.mat");
    load(SVD_file.folder+"/"+SVD_file.name);

    num_train = size(X_train,1);
    num_val = size(X_val,1);
    num_test = size(X_test,1);
    
    % emulate predictors
    Y_train=double(predict(netUpdated,X_train));
    Y_val=double(predict(netUpdated,X_val));
    Y_test=double(predict(netUpdated,X_test));

    figure;
    ntiles = size(T_train,2); 
    tlo = tiledlayout(ceil(sqrt(ntiles)),ceil(sqrt(ntiles)),"TileSpacing","tight");
    for m=1:size(T_train,2)
        tile(m)=nexttile; hold on;
        plot(tile(m),T_train(:,m),Y_train(:,m),'ok');
        plot(tile(m),T_val(:,m),Y_val(:,m),'xg');
        plot(tile(m),T_test(:,m),Y_test(:,m),'dm');
        %plot(T_test(10,m),Y_test(10,m),'dr');
        %plotregression(T_full,Y(:,m));
        title(tile(m),"mode "+num2str(m));
        xvec = xlim(tile(m));           
        plot(tile(m),[-10 10],[-10 10],'-r');
        xlim(tile(m),xvec); ylim(tile(m),xvec);
    end

    % reconstruct spatial maps of emulated targets
    Y_test = Y_test.*repmat(T_train_S,num_test,1)+...
        repmat(T_train_C,num_test,1); % undo normalization
    Y_reproj = Y_test*B_trunc;% reproject onto nodal basis
    
    % ua outputs
    Ind_test_orig = seq(num_train+num_val+1:end);
    Ua_orig = T(Ind_test_orig,:);
    Ua_proj = Ua_orig*T_reproj';

    % reconstruct original target data
    T_test = T_test.*repmat(T_train_S,num_test,1)+...
        repmat(T_train_C,num_test,1); % undo normalization
    T_test = T_test*B_trunc;% reproject onto nodal basis
    
    % calc mse of test data
    % in nodal basis
    mse = 1/(num_test-1)*sum((Y_reproj-Ua_orig).^2,1);
    %mse_tmp(mse_tmp>1e6) = 1e6;
    %nNodes = size(Ua_orig,2);
    %mse_tmp=spdiags(mse_tmp(:),0,nNodes,nNodes);

    % in truncated svd basis
    %mse_RNN = T_reproj*mse_tmp*T_reproj';   
    mse_RNN=1/(num_test-1)*sum((Y_test-Ua_proj).^2,1);

    MSE_file.name = strrep(SVD_file.name,"SVD","MSE_RNN");
    MSE_file.folder = SVD_file.folder;
    save(MSE_file.folder+"/"+MSE_file.name,"mse_RNN");

    % plot some maps for particular test case
    for nn=n_test
        predictorvalues = X_test(nn,:).*X_train_S + X_train_C;
        Ua_orig_tmp = Ua_orig(nn,:);
        Y_reproj_tmp = Y_reproj(nn,:);
        T_test_tmp = T_test(nn,:);
        plot_comparison_emulator_ua(T_mean,Ua_orig_tmp,T_test_tmp,Y_reproj_tmp,MUA,GF,"",predictorvalues,trainingdataformat);
    end

    % plot mse for different svd trunctions
    plot_mse(mse,"",MUA,GF);

end

%% FNN Matlab
if plotFNN
    for ii=1:numel(pct)
    
        tmp=load("./FNN/mat_files/FNN_trainscg_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+...
            "_floatingice"+string(1-only_grounded_ice)+"_N0k"+pct(ii)+".mat");
        load("./FNN/mat_files/SVD_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+...
            "_floatingice"+string(1-only_grounded_ice)+"_N0k"+pct(ii)+".mat","B_trunc","T_reproj","seq");
        net=tmp.Net_opt;
        
        X_train = net.X_train';
        X_val = net.X_val';
        X_test = net.X_test';
        num_train = size(X_train,1);
        num_val = size(X_val,1);
        num_test = size(X_test,1);
        T_train = net.T_train';
        T_val = net.T_val';
        T_test = net.T_test';

        % emulate predictors
        Y_train=double(net.trained(X_train')');
        Y_val=double(net.trained(X_val')');
        Y_test=double(net.trained(X_test')');
    
        figure;
        ntiles = size(T_train,2); 
        tlo = tiledlayout(ceil(sqrt(ntiles)),ceil(sqrt(ntiles)),"TileSpacing","tight");
        for m=1:size(T_train,2)
            tile(m)=nexttile; hold on;
            plot(tile(m),T_train(:,m),Y_train(:,m),'ok');
            plot(tile(m),T_val(:,m),Y_val(:,m),'xg');
            plot(tile(m),T_test(:,m),Y_test(:,m),'dm');
            %plot(T_test(10,m),Y_test(10,m),'dr');
            %plotregression(T_full,Y(:,m));
            title(tile(m),"mode "+num2str(m));
            xvec = xlim(tile(m));           
            plot(tile(m),[-10 10],[-10 10],'-r');
            xlim(tile(m),xvec); ylim(tile(m),xvec);
        end
    
        % reconstruct spatial maps of emulated targets
        Y_test = Y_test.*repmat(net.T_train_S',num_test,1)+...
            repmat(net.T_train_C',num_test,1); % undo normalization
        Y_reproj = Y_test*B_trunc;% reproject onto nodal basis
        
        % ua output
        Ind_test_orig = seq(num_train+num_val+1:end);
        Ua_orig = T(Ind_test_orig,:);        
        Ua_proj = Ua_orig*T_reproj';
    
        % reconstruct original target data
        T_test = net.T_test'.*repmat(net.T_train_S',num_test,1)+...
            repmat(net.T_train_C',num_test,1); % undo normalization
        T_test = T_test*B_trunc;% reproject onto nodal basis
        
        % calc mse of test data
        % in nodal basis
        mse_tmp = 1/(num_test-1)*sum((Y_reproj-Ua_orig).^2,1);        
        %mse_tmp(mse_tmp>1e6)=1e6;
        mse(ii,:)= mse_tmp;

        std_orig = std(Ua_orig,1);
        figure; PlotMeshScalarVariable([],MUA,std_orig(:)./mse_tmp(:));
        caxis([0 50]);
        title("std\_orig/mse\_tmp "+pct(ii));

        % in truncated svd basis
        nNodes = size(Ua_orig,2);
        mse_tmp=spdiags(mse_tmp(:),0,nNodes,nNodes);

        % in truncated svd basis
        mse_FNN = T_reproj*mse_tmp*T_reproj';   

        %mse_FNN=1/(num_test-1)*sum((Y_test-Ua_proj).^2,1);
        
        save("./FNN/mat_files/MSE_FNN_trainscg_"+perturbation+"_"+years+"_"+slidinglaw+"_cycle"+cycle+...
            "_floatingice"+string(1-only_grounded_ice)+"_N0k"+pct(ii),"mse_FNN");

        % plot some maps for particular test case
        for nn=n_test
            predictorvalues = X_test(nn,:).*net.X_train_S(:)' + net.X_train_C(:)';
            Ua_orig_tmp = Ua_orig(nn,:);
            Y_reproj_tmp = Y_reproj(nn,:);
            T_test_tmp = T_test(nn,:);
            plot_comparison_emulator_ua(T_mean,Ua_orig_tmp,T_test_tmp,Y_reproj_tmp,MUA,GF,pct(ii),predictorvalues,trainingdataformat);
        end
    end

    % plot mse for different svd trunctions
    plot_mse(mse,pct,MUA,GF);

end

end

function plot_comparison_emulator_ua(T_mean,T_orig,Y_reproj_orig,Y_reproj,MUA,GF,pct,predictorvalues,trainingdataformat)

%CM = othercolor('Reds8',15);
CM = crameri('batlow',15);
CM2 = flipdim(othercolor('RdYlBu11',15),1);

CtrlVar=Ua2D_DefaultParameters;
CtrlVar.PlotXYscale=1e3;
CtrlVar.PlotsXaxisLabel='' ; 
CtrlVar.PlotsYaxisLabel='' ;

H=fig('units','inches','width',120*12/72.27,'height',60*12/72.27,'fontsize',14,'font','Helvetica');

tlo=tiledlayout(H,1,3,"TileSpacing","tight");    

ax(1)=nexttile; hold on;
data_to_plot = T_orig(:);
PlotMeshScalarVariable(CtrlVar,MUA,data_to_plot);
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);
colormap(ax(1),CM2); 
cb=colorbar; 
if trainingdataformat == "LOGu"
    caxis(ax(1),[-1 1]); 
    cb.Label.String='log_{10}u [m/yr]';
else
    caxis(ax(1),[-300 300]); 
    cb.Label.String='u [m/yr]';
end
axis equal; axis tight; 
grid(ax(1),"off"); box(ax(1),"off"); axis(ax(1),"off");
yticklabels([]); xticklabels([]); 
title('Ua');

ax(2)=nexttile; hold on;
data_to_plot = Y_reproj_orig(:);
PlotMeshScalarVariable(CtrlVar,MUA,data_to_plot);
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);
colormap(ax(2),CM2); colorbar('off');
if trainingdataformat == "LOGu" 
    caxis(ax(2),[-1 1]); 
else
    caxis(ax(2),[-300 300]); 
end
axis equal; axis tight;
grid(ax(2),"off"); box(ax(2),"off"); axis(ax(2),"off");
yticklabels([]); xticklabels([]);
title("Truncated SVD ("+pct+"%)");

ax(3)=nexttile; hold on;
data_toplot = Y_reproj(:);
PlotMeshScalarVariable(CtrlVar,MUA,data_toplot);
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','k');
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','LineWidth',0.5);
colormap(ax(3),CM2); cb2=colorbar; 
if trainingdataformat == "LOGu" 
    cb2.Label.String='u [m/yr]'; 
    caxis(ax(3),[-1 1]);     
else
    cb2.Label.String='\Delta u [m/yr]';
    caxis(ax(3),[-300 300]); 
end
title("Emulator");

axis equal; axis tight;
grid(ax(3),"off"); box(ax(3),"off"); axis(ax(3),"off");
yticklabels([]); xticklabels([]);

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
function prepare_data_for_Delta_u_emulators

%% Script reads maps of speed change from perturbation experiments
%% ('Calv', 'dhIS', 'dh' or 'Calv_dh') between startyear and targetyear,
%% uses SVD to reduce the dimensionality for a range of truncated 
%% decompositions (describing 95% - 99.5% of the variance), then carries out
%% the FNN training and writes data for the Tensorflow RNN training.

addpath(getenv("froot_tools"));

%% -- perturbation data to use
domain = "AMUND";
perturbation = 'Calv_dh'; % can be 'Calv', 'dhIS', 'dh' or 'Calv_dh'
startyear = "2000";
targetyear = "2009";
slidinglaw = "Umbi"; % can be "Weertman" or "Umbi"
cycle = 2; % cycle 1: inversion without spinup, cycle 2: inversion after spinup
only_grounded_ice = 1;

%% -- SVD settings
add_measurements_to_SVD = 0;
use_existing_SVD_basis = 1;
SVD_basistoread = "Wc2_Uc2";
SVD_filetoread = "SVD_"+domain+"_"+perturbation+"_"+startyear+"-"+targetyear+...
    "_"+SVD_basistoread+"_floatingice"+string(1-only_grounded_ice)+"_includemeasurements0.mat";
max_SVD_components = 25; % for a large number of SVD components the ML 
% algorithms might struggle to learn the variability in the highest 
% compoments. From basic test 25 seems like a reasonable number.
SVD_targetpercentage = 0.999; % this needs to be high enough so we can
% reliably represent the measurements in the truncated SVD basis

%% -- FNN settings
FNN_train = 0; % train FNN emulator?
FNN_type = "feedforwardnet"; % feedforwardnet or cascadeforwardnet, a simple feedforwardnet seems to perform just fine
FNN_trainFcn = "trainscg"; % trainlm is fast on CPU and seems to perform just fine, 
% alternatives that seem to work well are trainscg
FNN_UseGPU = 0; % UseGPU=1 only works with trainFcn="trainscg"

%% -- outputs
do_plots = 1;
write_outputs_for_TF = 1;

%% ----------------------------------------- %%
%% No need to make changes beyond this point %%
%% ----------------------------------------- %%
% load X and T for perturbation data
data.domain = domain;
data.perturbation = perturbation;
data.startyear = startyear;
data.targetyear = targetyear;
data.slidinglaw = slidinglaw;
data.cycle = cycle;
data.only_grounded_ice = only_grounded_ice;
SVD_settings.add_measurements_to_SVD = add_measurements_to_SVD;

strtmp = char(slidinglaw);
strtmp = strtmp(1)+"c"+string(cycle(1))+"_";
SVD_filetowrite = "SVD_"+domain+"_"+perturbation+"_"+startyear+"-"+targetyear+"_"+...
    strtmp+"floatingice"+string(1-only_grounded_ice)+"_includemeasurements"+string(add_measurements_to_SVD)+".mat";

if ~exist(SVD_filetowrite)
    Calv_SVD_basis(data,SVD_settings,SVD_filetowrite);
end

SVD_new = load(SVD_filetowrite);

% load T for existing basis
if use_existing_SVD_basis
    SVD_existing = load(SVD_filetoread);
    S = SVD_existing.S;
    T = SVD_existing.T;
    T_mean = SVD_existing.T_mean;
else
    S = SVD_new.S;
    T = SVD_new.T;
    T_mean = SVD_new.T_mean;
end

%% DIMENSIONALITY REDUCTION: use singular value decomposition of targets to
%% express each map of \Delta u in terms of its k dominant principal components
% do the full svd to calculate how many components are required for each
% increment in pct
tmp = (cumsum(diag(S).^2)./sum(diag(S).^2));
T_nComp = find(tmp>SVD_targetpercentage,1,'first');

if T_nComp > max_SVD_components
    SVD_targetpercentage = tmp(max_SVD_components);
    T_nComp = max_SVD_components;
end

fprintf("Retaining %s components in the SVD with %s percentage of the variance explained.\n",...
    string(T_nComp),string(SVD_targetpercentage));

T_pct1 = diag(S).^2./sum(diag(S).^2);
T_pct = T_pct1(1:T_nComp);

% return the left singular vectors U, diagonal matrix S of singular 
% values, and right singular vectors V for T_nComp largest singular values
[U_trunc,S_trunc,V_trunc] = svds(T-repmat(T_mean(:)',size(T,1),1),T_nComp);

% note that T is transposed compared to the canonical treatment of the
% svd. 
% [T]=n-by-m with n number of samples and m size of each sample
% [U]=n-by-n 
% [S]=n-by-m
% [V']=m-by-m
% T = U*S*V' = U*B with B=S*V' with [B]=n-by-m
% truncate decomposation at rank r:
% [T_trunc]=n-by-m
% [U_trunc]=n-by-r
% [S_trunc]=r-by-r: the eigen values
% [V_trunc']=r-by-m: the eigen basis
% T_trunc = U_trunc*S_trunc*V_trunc' = U_trunc*B_trunc with 
% B_trunc=S_trunc*V_trunc' with [B_trunc]=r-by-m
% So T_trunc = U_trunc*B_trunc
% => T_trunc*B_trunc' = U_trunc*B_trunc*B_trunc'
% ==> T_trunc*B_trunc'*inv(B_trunc*B_trunc') = U_trunc
% and define T_reproj = inv(B_trunc*B_trunc')*B_trunc
% then T_reproj' = B_trunc'*inv(B_trunc*B_trunc')
% with [T_reproj']=m-by-r
% and U_trunc = T_trunc*T_reproj';
% Instead of training the emulator with U_trunc, we train the emulator
% with T_hat = T*T_reproj'
% To reproject T_hat to the original basis, simply multiply with 
% B_trunc: T = T_hat*B_trunc
B_trunc = S_trunc*V_trunc';
T_reproj = (B_trunc*B_trunc')\B_trunc; % equivalent to inv(B*B')*B

%% TARGETS (training/input data). In this case, input data consists of 
%% simulated instantaneous changes in surface speed in response to changes
%% in ice-sheet geometry (ice thickness, calving front location). The input
%% data consists of num_nodes nodal values for nun_exp experiments.
if add_measurements_to_SVD
    T = SVD_new.T(2:end,:)
    T_hat = (T-repmat(T_mean(:)',size(T,1),1))*T_reproj';
    % we remove the first line here, which contains the measurements
    % because we don't want to include them in the training data        
else
    T = SVD_new.T;
    T_hat = (T-repmat(T_mean(:)',size(T,1),1))*T_reproj';
end

num_exp = size(SVD_new.T,1); num_nodes = size(SVD_new.T,2);

%% PREDICTORS
X = SVD_new.X;

if do_plots

    % load Ua mesh
    load("Delta_u_"+domain+"_"+slidinglaw+"_"+startyear+"-"+targetyear+".mat","MUA_yr2","GF_yr2");
    MUA = MUA_yr2; GF = GF_yr2;

    %% Can measurements be represented adequately in truncated basis?
    T_meas = SVD_new.T_meas(:)'-T_mean(:)';
    T_meas_trunc = T_meas*T_reproj';
    T_meas_reproj = T_meas_trunc*B_trunc; % reproject truncated field onto nodal basis

    figure; tlo=tiledlayout(1,3,"TileSpacing","compact");
    
    t(1)=nexttile; hold on; PlotMeshScalarVariable([],MUA,T_meas(:)+T_mean(:)); 
    title(t(1),'\Delta U_{meas}');
    axis(t(1),"equal","tight"); caxis(t(1),[0 1000]); colormap(t(1),flipud(cmocean('haline')));
    
    t(2)=nexttile; hold on; PlotMeshScalarVariable([],MUA,T_meas_reproj(:)+T_mean(:)); 
    title(t(2),'Reprojected \Delta U_{meas}');
    axis(t(2),"equal","tight"); caxis(t(2),[0 1000]); colormap(t(2),flipud(cmocean('haline')));
    
    t(3)=nexttile; hold on; PlotMeshScalarVariable([],MUA,T_meas(:)-T_meas_reproj(:)); 
    title(t(3),'Original - Reprojected');
    axis(t(3),"equal","tight"); caxis(t(3),[-250 250]); colormap(t(3),flipud(cmocean('balance')));
    
    title(tlo,strrep(SVD_filetoread+", SVD truncation: "+string(T_nComp)+" components, "+string(SVD_targetpercentage*100)+"%.","_","\_"));
    
    %% Are model data adequately represented in truncated basis?
    if use_existing_SVD_basis
    
        T_reproj_tmp = T_hat*B_trunc;
    
        figure; tlo=tiledlayout(1,2,"TileSpacing","compact");
    
        nexttile; hold on; PlotMeshScalarVariable([],MUA,mean(T-repmat(T_mean(:)',size(T,1),1)-T_reproj_tmp,1)'); title('mean(Original - Reprojected)');
        axis equal; axis tight;  caxis([-25 25]); colormap(cmocean('balance'));
    
        nexttile; hold on; PlotMeshScalarVariable([],MUA,std(T-repmat(T_mean(:)',size(T,1),1)-T_reproj_tmp,1)'); title('std(Original - Reprojected)');
        axis equal; axis tight; caxis([-150 150]); colormap(cmocean('balance'));
    
        title(tlo,[strrep(domain+", "+perturbation+", "+startyear+"-"+targetyear+", "+slidinglaw+", cycle: "+string(cycle)+...
            ", floating ice: "+string(1-only_grounded_ice),"_","\_");...
            "SVD: "+strrep(SVD_filetoread,"_","\_");...
            "SVD truncation: "+string(T_nComp)+" components, "+string(SVD_targetpercentage*100)+"%."]);
    
    end
end

%% scramble data
seq = randperm(num_exp);
data = T_hat(seq,:);
predictors = X(seq,:);

%% Now simulate FeedForward NN
if FNN_train
    year = startyear + "-" + targetyear;
    fname1 = sprintf('./FNN/mat_files/data_%s_%s_%s-%s_%s_cycle%s_floatingice%s_includemeasurements%s_Comp%s_Var%s',...
        domain,perturbation,startyear,targetyear,...
        slidinglaw,string(cycle),string(1-only_grounded_ice),string(add_measurements_to_SVD),string(T_nComp),...
        strrep(string(round(SVD_targetpercentage*1e4)/1e2),".","k"));
    save(fname1, 'T_mean','V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');
    fname2 = sprintf('./FNN/mat_files/SVD_%s_%s_%s-%s_%s_cycle%s_floatingice%s_includemeasurements%s_Comp%s_Var%s',...
        domain,perturbation,startyear,targetyear,...
        slidinglaw,string(cycle),string(1-only_grounded_ice),string(add_measurements_to_SVD),string(T_nComp),...
        strrep(string(round(SVD_targetpercentage*1e4)/1e2),".","k"));
    addpath("./FNN");
    TrainFNN(predictors',data',FNN_type,FNN_trainFcn,FNN_UseGPU,fname2,do_plots);
end

% Emulate full dataset
% X_full = [Net.X_train Net.X_val Net.X_test];
% Y = Net.trained(X_full);
% T_full = [Net.T_train Net.T_val Net.T_test];
% 
% % Plot emulator vs targets
% for ii=1:size(T_full,1)
%     figure;
%     plotregression(T_full(ii,:),Y(ii,:));
%     title("mode "+num2str(ii));
% end

%% WRITE OUTPUTS for tensorflow
if write_outputs_for_TF
    % split targets and predictors into training, cross-validation and test
    % datasets
    val_idx = floor(num_exp*0.8);
    test_idx = floor(num_exp*0.1) + val_idx;
    
    T_train = data(1:val_idx,:); num_train = size(T_train,1);
    T_val = data(val_idx+1:test_idx,:); num_val = size(T_val,1);
    T_test = data(test_idx+1:end,:); num_test = size(T_test,1);
    
    X_train = predictors(1:val_idx,:);
    X_val = predictors(val_idx+1:test_idx,:);
    X_test = predictors(test_idx+1:end,:);

    % need to feed normalized data to tensorflow
    [X_train,X_train_C,X_train_S]=normalize(X_train);
    X_val = (X_val-repmat(X_train_C,num_val,1))./repmat(X_train_S,num_val,1);
    X_test = (X_test-repmat(X_train_C,num_test,1))./repmat(X_train_S,num_test,1);

    [T_train,T_train_C,T_train_S]=normalize(T_train);
    T_val = (T_val-repmat(T_train_C,num_val,1))./repmat(T_train_S,num_val,1);
    T_test = (T_test-repmat(T_train_C,num_test,1))./repmat(T_train_S,num_test,1);

    fname1 = sprintf('./RNN/mat_files/data_%s_%s_%s-%s_%s_cycle%s_floatingice%s_includemeasurements%s_Basis_%s-Comp%s-Var%s',...
        domain,perturbation,startyear,targetyear,...
        slidinglaw,string(cycle),string(1-only_grounded_ice),string(add_measurements_to_SVD),...
        strrep(SVD_basistoread,"_","-"),...
        string(T_nComp),strrep(string(round(SVD_targetpercentage*1e4)/1e2),".","k"));
    fname2 = sprintf('./RNN/mat_files/SVD_%s_%s_%s-%s_%s_cycle%s_floatingice%s_includemeasurements%s_Basis_%s-Comp%s-Var%s',...
        domain,perturbation,startyear,targetyear,...
        slidinglaw,string(cycle),string(1-only_grounded_ice),string(add_measurements_to_SVD), ...
        strrep(SVD_basistoread,"_","-"),...
        string(T_nComp),strrep(string(round(SVD_targetpercentage*1e4)/1e2),".","k"));
    save(fname1,'X_train','X_val','X_test','T_train','T_val','T_test',...
        'X_train_C','X_train_S','T_train_C','T_train_S');
    save(fname2, 'T_mean','V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');
end

%% PLOTTING

if do_plots

    CtrlVar=Ua2D_DefaultParameters;
    CtrlVar.PlotXYscale = 1e3;

    cumulative_variance = cumsum(diag(S).^2)/sum(diag(S).^2);
    
    % figure(111), tlo1=tiledlayout(2,2,'TileSpacing','tight'); title(tlo1,'A');
    % nexttile;
    % imagesc(A), axis off; colormap(slanCM('YlGnBu')); cb1=colorbar;
    % title('Original');
    
    % figure(222), tlo2=tiledlayout(2,2,'TileSpacing','tight'); title(tlo2,'t = 1 year');
    % nexttile; 
    % PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,A(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu')); cb2=colorbar(gca);
    % title('Original');
    % 
    % figure(333), tlo3=tiledlayout(2,2,'TileSpacing','tight'); title(tlo3,'t = 850 years');
    % nexttile; 
    % PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,data(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu')); cb3=colorbar(gca);
    % title('Original');
    
    figure; tlo4=tiledlayout(2,5,'TileSpacing','tight');
    for ii=1:10
        nexttile;
        PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,V_trunc(:,ii),CtrlVar);
        colormap(othercolor('RdYlBu8'));
        variance =100*[cumulative_variance(1); cumulative_variance(2:end)-cumulative_variance(1:end-1)];
        title(['mode',num2str(ii),' (',num2str(variance(ii),'%2.2f'),'%  /  ',num2str(100*cumulative_variance(ii),'%2.2f'),'%)']); caxis([-0.03 0.03]);
        axis tight; axis off;
    
        %cb=colormap; cb.visible='off';
    end
    cb=colorbar(gca);
    plotind = 2;
    for r = [10 50 100]
        %Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
        %Xapprox(Xapprox<0)=0; 
        %Xapprox(Xapprox>1)=1;
        %figure(111); nexttile, imagesc(Xapprox), axis off, colormap(slanCM('YlGnBu'));
        %figure(222); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Xapprox(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
        %title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
        %figure(333); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Xapprox(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
        %title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
        %plotind = plotind + 1;
    end

    cb.Layout.Tile = 'east';
    
    %% singular values
    figure, subplot(1,2,1)
    semilogy(diag(S),'k','linewidth',2); grid on;
    xlabel('r');
    ylabel('singlar value, \sigma_r');
    set(gca,'fontsize',14);
    subplot(1,2,2)
    plot(cumulative_variance,'k','linewidth',2); grid on;
    xlabel('r');
    ylabel('cumulative energy');
    set(gca,'fontsize',14);
end
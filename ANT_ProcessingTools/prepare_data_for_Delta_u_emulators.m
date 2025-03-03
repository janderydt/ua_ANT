function prepare_data_for_Delta_u_emulators

%% Script reads maps of speed change from perturbation experiments
%% ('Calv', 'dhIS', 'dh' or 'Calv_dh') between startyear and targetyear,
%% uses SVD to reduce the dimensionality for a range of truncated 
%% decompositions (describing 95% - 99.5% of the variance), then carries out
%% the FNN training and writes data for the Tensorflow RNN training.

addpath(getenv("froot_tools"));

perturbation = 'Calv_dh'; % can be 'Calv', 'dhIS', 'dh' or 'Calv_dh'
startyear = "2000";
targetyear = "2020";
slidinglaw = "Umbi";
cycle=1; % cycle 1: inversion without spinup, cycle 2: inversion after spinup
FNNtype = "feedforwardnet"; % feedforwardnet or cascadeforwardnet, a simple feedforwardnet seems to perform just fine
trainFcn = "trainscg"; % trainlm is fast on CPU and seems to perform just fine, 
% alternatives that seem to work well are trainscg
UseGPU = 0; % UseGPU=1 only works with trainFcn="trainscg"
doplots = 0;
write_outputs_for_TF = 1;
add_measurements_to_SVD = 0;
only_grounded_ice = 1;
if only_grounded_ice
    pct = 0.99; %% chose high enough value so we can adequately represent 
% the measurements in the truncated basis
else
    pct = 0.999; %% chose high enough value so we can adequately represent 
% the measurements in the truncated basis
end

% load data file || this file is produced by the
% plot_PerturbationResults_Ensemble.m function
load("Delta_u_AMUND_"+slidinglaw+"_"+startyear+"-"+targetyear+".mat");

%% PREDICTORS
Ind_toremove = find(gaC>50);
X = [log10(gaA(:)) log10(gaC(:)) log10(gsA(:)) log10(gsC(:)) m(:) n(:)];
if cycle==2 % add dhdt_err to predictors
    X = [X log10(dhdt_err(:))];
end
X = double(X);
X(Ind_toremove,:)=[];

%% TARGETS (training/input data). In this case, input data consists of 
%% simulated instantaneous changes in surface speed in response to changes
%% in ice-sheet geometry (ice thickness, calving front location). The input
%% data consists of num_nodes nodal values for nun_exp experiments.
T = Delta_u.(perturbation).map(:,:,cycle);
% remove gaA>50
T(Ind_toremove,:)=[];
% check for nans and inf
if any(isnan(T) | isinf(T))
    warning("Removing nan from T");
    [rows,~]=find(isnan(T));
    rows = unique(rows);
    X(rows,:)=[];
    T(rows,:)=[];
end
% check dimensions of input data; rows: experiments, columns: nodes
if numel(size(T))==2
    num_exp = size(T,1); num_nodes = size(T,2);
    if num_nodes < num_exp % number of experiments is highly likely to be smaller than number of nodes
        T = T'; % transpose data
        num_exp = size(T,1); num_nodes = size(T,2);
    end
else
    error("check dimensions of input data");
end
% remove mean
T_mean = mean(T,1);
T = T-repmat(T_mean(:)',size(T,1),1);

%% Load observations
for yy=[startyear targetyear]
    if yy=="2000"
        fname = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    else
        fname = "GriddedInterpolants_"+string(double(yy)-1)+"-"+string(double(yy))+...
        "_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    end
    load("../ANT_Data/ANT_Interpolants/"+fname);
    v_tmp = hypot(Fus.Values,Fvs.Values);
    Fu = Fus; Fu.Values = v_tmp;
    std_tmp = hypot(Fxerr.Values,Fyerr.Values);
    Fstd = Fus; Fstd.Values = std_tmp;
    v.("yr"+yy).F = Fu;
    std.("yr"+yy).F = Fstd;
end
load("Delta_u_AMUND_Weertman_"+startyear+"-"+targetyear+".mat","MUA_yr2","GF_yr2");
MUA = MUA_yr2; GF = GF_yr2;
% interpolate initial and final speed onto Ua mesh
u_init = v.("yr"+startyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
std_init = std.("yr"+startyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
u_target = v.("yr"+targetyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
std_target = std.("yr"+targetyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltau = u_target-u_init; % measured change in speed (m/yr)
% subtract model ensemble mean
T_meas = deltau(:) - T_mean(:);
if any(isnan(T_meas))
    error("dT_meas contains nans. Consider using extruded field or inpaint_nans.");
end

if only_grounded_ice
    %% remove floating nodes
    Nodes_floating = find(GF.node<0.5);
    T(:,Nodes_floating) = 0;
    T_meas(Nodes_floating) = 0;
    T_mean(Nodes_floating) = 0;
end

if add_measurements_to_SVD
    T = [T_meas(:)' ; T];
    datastr = "withdata";
else
    datastr = "nodata";
end    

%% DIMENSIONALITY REDUCTION: use singular value decomposition of targets to
%% express each map of \Delta u in terms of its k dominant principal components
% do the full svd to calculate how many components are required for each
% increment in pct
[~,S,~] = svd(T,'econ');

seq = randperm(num_exp);

for ii=1:numel(pct)
    tmp = (cumsum(diag(S).^2)./sum(diag(S).^2));
    T_nComp = find(tmp>pct(ii),1,'first');

    if T_nComp>25
        pct(ii) = tmp(25);
        T_nComp = 25;
    end

    fprintf("Retaining %s components in the SVD with %s percentage of the variance explained.\n",...
        string(T_nComp),string(pct(ii)));

    T_pct1 = diag(S).^2./sum(diag(S).^2);
    T_pct = T_pct1(1:T_nComp);

    % return the left singular vectors U, diagonal matrix S of singular 
    % values, and right singular vectors V for T_nComp largest singular values
    [U_trunc,S_trunc,V_trunc] = svds(T,T_nComp);

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

    %% Plot eigen basis
    %for jj=1:T_nComp
     %   figure; PlotMeshScalarVariable([],MUA,V_trunc(:,jj));...
     %   axis equal; axis tight; axis off; 
     %   title(string(pct)+" SVD eigenfunction "+string(jj)+" ("+string(T_pct(jj)*100)+"%)");
    %end
   
    if add_measurements_to_SVD
        T_hat = T(2:end,:)*T_reproj'; 
        % we remove the first line here, which contains the measurements
        % because we don't want to include them in the training data        
    else
        T_hat = T*T_reproj';
    end

    %% Can measurements be represented adequately in truncated basis?
    T_meas_trunc = T_meas(:)'*T_reproj';
    T_meas_reproj = T_meas_trunc*B_trunc; % reproject truncated field onto nodal basis

    figure; subplot(1,3,1); hold on; PlotMeshScalarVariable([],MUA,T_meas(:)+T_mean(:)); title('\Delta U_{meas}');
    axis equal; axis tight; caxis([0 1000]);
    
    subplot(1,3,2); hold on; PlotMeshScalarVariable([],MUA,T_meas_reproj(:)+T_mean(:)); title('Reprojected \Delta U_{meas}');
    axis equal; axis tight; caxis([0 1000]);

    subplot(1,3,3); hold on; PlotMeshScalarVariable([],MUA,T_meas(:)-T_meas_reproj(:)); title('Original - Reprojected');
    axis equal; axis tight; caxis([-250 250]);

    data = T_hat(seq,:);
    predictors = X(seq,:);

    %% Now simulate FeedForward NN
    year = startyear + "-" + targetyear;
    fname1 = sprintf("./FNN/mat_files/SVD-%s_%s_%s_%s_cycle%s_floatingice%s_N0k%.4g",datastr,perturbation,year,...
        slidinglaw,string(cycle),string(1-only_grounded_ice),10000*pct(ii));
    save(fname1, 'T_mean','V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');
    fname2 = sprintf("./FNN/mat_files/FNN_SVD-%s_%s_%s_%s_%s_cycle%s_floatingice%s_N0k%.4g",datastr,trainFcn,perturbation,year,...
        slidinglaw,string(cycle),string(1-only_grounded_ice),10000*pct(ii));
    addpath("./FNN");
    TrainFNN(predictors',data',FNNtype,trainFcn,UseGPU,fname2,doplots);

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

        fname1 = sprintf('./RNN/mat_files/data_SVD-%s_%s_%s_%s_cycle%s_floatingice%s_N0k%.4g',datastr,perturbation,year,...
            slidinglaw,string(cycle),string(1-only_grounded_ice),pct(ii)*10000);
        fname2 = sprintf('./RNN/mat_files/SVD-%s_%s_%s_%s_cycle%s_floatingice%s_N0k%.4g',datastr,perturbation,year,...
            slidinglaw,string(cycle),string(1-only_grounded_ice),pct(ii)*10000);
        save(fname1,'X_train','X_val','X_test','T_train','T_val','T_test',...
            'X_train_C','X_train_S','T_train_C','T_train_S');
        save(fname2, 'T_mean','V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');
    end
end

%% PLOTTING

if doplots

    load("perturbation_grids.mat");
    if targetyear == "2009"
        MUA=MUA_2009;
    elseif targetyear == "2014"
        MUA=MUA_2014;
    elseif targetyear == "2018"
        MUA=MUA_2018;
    elseif targetyear == "2020"
        MUA=MUA_2018;
    end

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
        PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,U(:,ii),CtrlVar);
        colormap(othercolor('RdYlBu8'));
        variance =100*[cumulative_variance(1); cumulative_variance(2:end)-cumulative_variance(1:end-1)];
        title(['mode',num2str(ii),' (',num2str(variance(ii),'%2.2f'),'%  /  ',num2str(100*cumulative_variance(ii),'%2.2f'),'%)']); caxis([-0.03 0.03]);
        axis tight; axis off;
    
        %cb=colormap; cb.visible='off';
    end
    cb4=colorbar(gca);
    plotind = 2;
    for r = [10 50 100]
        Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
        %Xapprox(Xapprox<0)=0; 
        %Xapprox(Xapprox>1)=1;
        %figure(111); nexttile, imagesc(Xapprox), axis off, colormap(slanCM('YlGnBu'));
        %figure(222); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Xapprox(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
        %title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
        %figure(333); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Xapprox(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
        %title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
        %plotind = plotind + 1;
    end
    
    cb1.Layout.Tile = 'east';
    cb2.Layout.Tile = 'east';
    cb3.Layout.Tile = 'east';
    cb4.Layout.Tile = 'east';
    
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
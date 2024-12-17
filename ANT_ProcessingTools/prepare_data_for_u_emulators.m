function prepare_data_for_u_emulators

%% Script reads maps of speed from inversion and perturbation experiments
%% ('Calv', 'dhIS', 'dh' or 'Calv_dh') for the years specified below.
%% For each year, a SVD is used to reduce the dimensionality for a range of 
%% truncated decompositions (describing 95% - 99.5% of the variance), then 
%% a simple FNN model is trained, and data is written for Tensorflow RNN training.

addpath(getenv("froot_tools"));
addpath("..");

years_to_analyze = [2000 2009 2014 2018];
perturbation = 'Calv_dh'; % can be 'Calv', 'dhIS', 'dh' or 'Calv_dh'
basins_to_analyze = {'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbott
slidinglaw = "Weertman";
cycle=2; % cycle 1: inversion without spinup, cycle 2: inversion after spinup
FNNtype = "feedforwardnet"; % feedforwardnet or cascadeforwardnet, a simple feedforwardnet seems to perform just fine
trainFcn = "trainscg"; % trainlm is fast on CPU and seems to perform just fine, 
% alternatives that seem to work well are trainlm
UseGPU = 1; % UseGPU=1 only works with trainFcn="trainscg"
doplots = 0;
writeoutputsforTF = 1;

file_with_perturbation_data_to_read = "perturbationdata_"+slidinglaw+".mat";
file_with_speedmaps_to_read = "u_AS_"+perturbation+"_cycle"+string(cycle)+"_"+slidinglaw+"_"+strjoin(string(years_to_analyze),"_")+".mat";

%% Gather speed maps
if ~exist(file_with_speedmaps_to_read,"file")
    % load perturbation data
    load(file_with_perturbation_data_to_read);
    load("perturbation_grids.mat");
    % load basins
    filename = 'basins_IMBIE_v2.mat'; 
    B = load(filename);
    B = RemoveIceRisesAndIslands(B,1e10); % cutoff=1e10: remove all ice rises and islands
    % prepare meshes
    CtrlVar = Ua2D_DefaultParameters;
    UserVar.casefolder=pwd;
    UserVar.datafolder=pwd+"/../ANT_Data/";
    UserVar.Experiment="";
    BaseMesh='2000_2009_2014_2018_meshmin3000_meshmax100000_refined';
    UserVar = ANT_DefineBaseMesh(UserVar,BaseMesh);
    for yy=1:numel(years_to_analyze)
        yrstring = "yr"+string(years_to_analyze(yy));
        UserVar.Geometry=double(years_to_analyze(yy));
        UserVar=ANT_ApplyMeshModifications(UserVar);
        tmp = load(UserVar.InitialMeshFileName);
        MUA.(yrstring) = tmp.MUA; 
        % identify basin id of each MUA node
        [MUA.(yrstring).basins,~] = Define_Quantity_PerBasin(MUA.(yrstring).coordinates(:,1),MUA.(yrstring).coordinates(:,2),B,0);
        MUA_basinnames = erase({MUA.(yrstring).basins(:).name},'-');
        basinnodes_all = [];
        for bb=1:numel(basins_to_analyze) 
            basin = char(erase(basins_to_analyze{bb},'-'));
            [~,BasinInd] = ismember(basin,MUA_basinnames);
            basinnodes = MUA.(yrstring).basins(BasinInd).ind;
            basinnodes_all = [basinnodes_all; basinnodes];
        end
        ElementsToBeDeactivated=any(~ismember(MUA.(yrstring).connectivity,basinnodes_all),2);
        [MUA.(yrstring),MUA.(yrstring).k,MUA.(yrstring).l]=DeactivateMUAelements(CtrlVar,MUA.(yrstring),ElementsToBeDeactivated);
        Ind_nan = find(isnan(MUA.(yrstring).Boundary.x));
        if ~isempty(Ind_nan)
            MUA.(yrstring).Boundary.x = MUA.(yrstring).Boundary.x(1:Ind_nan(1)-1);
            MUA.(yrstring).Boundary.y = MUA.(yrstring).Boundary.y(1:Ind_nan(1)-1);
        end
        original_node_numbers = MUA.(yrstring).k(find(~isnan(MUA.(yrstring).k)));
        switch yrstring
            case "yr2000"
                GF.(yrstring).node = GF_2000.node(original_node_numbers);
            case "yr2009"
                GF.(yrstring).node = GF_2009.node(original_node_numbers);
            case "yr2014"
                GF.(yrstring).node = GF_2014.node(original_node_numbers);
            case "yr2018"
                GF.(yrstring).node = GF_2018.node(original_node_numbers);
        end
        for ii=1:numel(data)
            if yrstring=="yr2000"                
                speed.(yrstring)(ii,:) = data(ii).Original.speed(original_node_numbers,cycle)';
            else
                Ind_tmp = find(data(ii).(perturbation).(yrstring).cycle==cycle);
                if ~isempty(Ind_tmp)
                    speed.(yrstring)(ii,:) = data(ii).(perturbation).(yrstring).speed(original_node_numbers,Ind_tmp)';
                else
                    speed.(yrstring)(ii,:) = nan*original_node_numbers';
                end
            end   
        end
    end

    m=[data(:).m];
    n=[data(:).n];
    for ii=1:numel(data)

        dhdt_err(ii)=data(ii).Inverse.dhdt_err;
        gaA(ii)=data(ii).Inverse.gaA;
        gaC(ii)=data(ii).Inverse.gaC;
        gsA(ii)=data(ii).Inverse.gsA;
        gsC(ii)=data(ii).Inverse.gsC;

    end

    save(file_with_speedmaps_to_read,"speed","GF","MUA","dhdt_err","m","n","gaA","gaC","gsA","gsC","-v7.3");

else

    load(file_with_speedmaps_to_read);

end

%% PREDICTORS
Ind_toremove = find(gaC>50);
X = [log10(gaA(:)) log10(gaC(:)) log10(gsA(:)) log10(gsC(:)) m(:) n(:)];
if cycle==2 % add dhdt_err to predictors
    X = [X log10(dhdt_err(:))];
end
X = double(X);
X(Ind_toremove,:)=[];

X_orig = X;

%% TARGETS (training/input data). In this case, input data consists of 
%% simulated surface speed at the end of the inversion (yr2000) or in 
%% response to changes in ice-sheet geometry (ice thickness, calving front location). 
%% The input data consists of num_nodes nodal values for nun_exp experiments.
for yy=1:numel(years_to_analyze)
    X = X_orig;
    yrstring = "yr"+string(years_to_analyze(yy));
    T = speed.(yrstring);
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
    T_mean=mean(T,1);
    T = T-repmat(T_mean(:)',size(T,1),1);

    %% DIMENSIONALITY REDUCTION: use singular value decomposition of targets to
    %% express each map of \Delta u in terms of its k dominant principal components
    %% do the full svd to calculate how many components are required for each
    %% increment in pct
    [~,S,~] = svd(T,'econ');
    
    seq = randperm(num_exp);
    pct = [0.95 0.96 0.97 0.98 0.99 0.995];
    
    for ii=1:numel(pct)
        T_nComp = find((cumsum(diag(S).^2)./sum(diag(S).^2))>pct(ii),1,'first');
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
        % [S_trunc]=r-by-r
        % [V_trunc']=r-by-m
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
        T_hat = T*T_reproj';
        
        data = T_hat(seq,:);
        predictors = X(seq,:);
    
        %% Now simulate FeedForward NN
        year = string(years_to_analyze(yy));
        fname1 = sprintf("./FNN/mat_files/SVD_%s_%s_%s_cycle%s_N0k%.3g",perturbation,year,slidinglaw,string(cycle),1000*pct(ii));
        save(fname1, 'T_mean','V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');
        fname2 = sprintf("./FNN/mat_files/FNN_%s_%s_%s_%s_cycle%s_N0k%.3g",trainFcn,perturbation,year,slidinglaw,string(cycle),1000*pct(ii));
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
        if writeoutputsforTF
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
    
            fname1 = sprintf('./RNN/mat_files/data_%s_%s_%s_cycle%s_N0k%.3g',perturbation,year,slidinglaw,string(cycle),pct(ii)*1000);
            fname2 = sprintf('./RNN/mat_files/SVD_%s_%s_%s_cycle%s_N0k%.3g',perturbation,year,slidinglaw,string(cycle),pct(ii)*1000);
            save(fname1,'X_train','X_val','X_test','T_train','T_val','T_test',...
                'X_train_C','X_train_S','T_train_C','T_train_S');
            save(fname2, 'T_mean','V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');
        end
    end
end

return
%% PLOTTING

if doplots

    load("perturbation_grids.mat");
    if targetyear == "2009"
        MUA=MUA_2009;
    elseif targetyear == "2014"
        MUA=MUA_2014;
    elseif targetyear == "2018"
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
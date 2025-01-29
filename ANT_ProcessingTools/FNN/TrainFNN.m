function Net_opt = TrainFNN(X,T,FNNtype,trainFcn,UseGPU,filename,doplots)

%% ----------------------------------------------------------------------- %%
% This script takes data to train a simple feedforward neural 
% network (FNN) with 2 hidden layers. The input data is split into a 
% training set (80%) which is fed to the neural network, and a 
% cross-validation set (20%), which is used to find the optimal size of 
% each hidden layer. We use a k-fold cross-validation approach, where we 
% split the data into k subsets, and train on k-1 different subsets, k 
% times. By default we use k=10. The size of layer 1 (nL1) ranges from 
% 2 to 15, and the size of layer 2 (nL2) ranging from 1 to nL1. The 
% performance is measured using a simple cost function (misfit)
% that calculates the mean square error between FNN outputs and targets: 
% J=1/m*sum((trained-targets)^2)
%% ----------------------------------------------------------------------- %%
% 
%% INPUTS: 
% X: n x m vector with n the number of predictor parameters and m the 
% number of samples
% T: k x m vector with k the number of target parameters and m the number
% of samples
% FNNtype: feedforwardnet or cascadeforwardnet
% UseGPU: 0 or 1
% trainFcn: training function (examples: 'trainlm', 'trainbr', 'trainscg')
% filename: where to save the network output
% doplots: 0 or 1
%
%% OUTPUTS:
% Net_opt: optimal network

kfold = 10;
layersizemax = 10;

addpath(getenv("froot_tools"));

% step 1: remove test data (10%)
num_exp = size(X,2);
test_idx = floor(num_exp*0.9);
X_test = X(:,test_idx+1:end); num_test = size(X_test,2);
T_test = T(:,test_idx+1:end); 
X_res = X(:,1:test_idx);
T_res = T(:,1:test_idx);

% step 2: k-fold partitioning into training and cross-validation
if ~exist(filename,"file")

    % create k-fold partition 
    c = cvpartition(test_idx,"KFold",kfold);

    % train model for k partitions
    for it=1:kfold
    
        idx = training(c,it);
        train_idx = find(idx==1);
        val_idx = find(idx==0);

        % split data into training and cross-validation segments
        X_train = X_res(:,train_idx); num_train = size(X_train,2);
        X_val = X_res(:,val_idx); num_val = size(X_val,2);
        T_train = T_res(:,train_idx);
        T_val = T_res(:,val_idx);

        % normalize predictors
        [X_train,X_train_C,X_train_S] = normalize(X_train,2); % X_norm = (X-X_C)/X_S
        X_val = (X_val-repmat(X_train_C,1,num_val))./repmat(X_train_S,1,num_val);
        % normalize targets
        [T_train,T_train_C,T_train_S] = normalize(T_train,2);
        T_val = (T_val-repmat(T_train_C,1,num_val))./repmat(T_train_S,1,num_val);

        kk=1;
        
        for ii=2:layersizemax
        
            for jj=1:ii
        
                if FNNtype == "feedforwardnet"
                    net = feedforwardnet([ii jj],trainFcn);
                elseif FNNtype == "cascadeforwardnet"
                    net = cascadeforwardnet([ii jj],trainFcn);
                else
                    error("Unknown FNN type "+FNNtype);
                end

                % set early stopping parameters
                net.divideFcn= 'dividerand';
                net.divideParam.trainRatio = 0.9; % training set [%]
                net.divideParam.valRatio = 0.1; % validation set [%]
                net.divideParam.testRatio = 0; % test set [%]   
                %net.inputs{1}.processFcns = {'mapstd'}; % Normalize inputs/targets to have zero mean and unity variance
                net.trainParam.showWindow = 0;
                net.performFcn = 'mse'; % for other option use >>help nnperformance

                net = configure(net,X_train,T_train);

                if UseGPU==1
                    X_train_gpu = gpuArray(X_train);
                    T_train_gpu = gpuArray(T_train);                 
                    Net(kk).it(it).trained = train(net,X_train_gpu,T_train_gpu,'showResources','no','useGPU','only');
                else
                    Net(kk).it(it).trained = train(net,X_train,T_train,'showResources','no');
                end

                Net(kk).nL1 = ii;
                Net(kk).nL2 = jj;
                fprintf("("+num2str(ii)+","+num2str(jj)+")\n");

                kk=kk+1;
        
            end
            
        end
        
        % Now calculate cost functions
        for kk=1:numel(Net)

            % performance on training data using mse
            Y = Net(kk).it(it).trained(X_train);
            Net(kk).J_train(it) = 0.5/size(X_train,2)*sum((T_train(:)-Y(:)).^2);
                
            Ptmp = polyfit(T_train(:),Y(:),1);
            Net(kk).slope(it)=Ptmp(1);
            Net(kk).intercept(it)=Ptmp(2);
        
            [Rtmp,~]=corrcoef(T_train(:),Y(:));
            Net(kk).R(it)=Rtmp(2,1);
        
            % performance on cross-validation data
            Y = Net(kk).it(it).trained(X_val);
            Net(kk).J_val(it) = 0.5/size(X_val,2)*sum((T_val(:)-Y(:)).^2);
        end

        fold(it).X_train = X_train;
        fold(it).X_train_C = X_train_C;
        fold(it).X_train_S = X_train_S;
        fold(it).T_train = T_train;
        fold(it).T_train_C = T_train_C;
        fold(it).T_train_S = T_train_S;
        fold(it).X_val = X_val;
        fold(it).T_val = T_val;

        fprintf("Done "+num2str(it)+"/"+num2str(kfold)+" iterations in k-fold cross-validation.\n")
    
    end
    

    for kk=1:numel(Net)
    
        J_train = Net(kk).J_train;
        J_val = Net(kk).J_val;
    
        Net(kk).J_train_min = min(J_train);
        Net(kk).J_train_max = max(J_train);
        Net(kk).J_train_mean = mean(J_train);
        Net(kk).J_train_std = std(J_train);
    
        Net(kk).J_val_min = min(J_val);
        Net(kk).J_val_max = max(J_val);
        Net(kk).J_val_mean = mean(J_val);
        Net(kk).J_val_std = std(J_val);
    
    end

    % Choose 'optimal' architecture based on mean of cost function
    [~,ind]=min([Net(:).J_val_mean]);
    [~,ind2] = min(Net(ind).J_val);
    Net_opt.trained = Net(ind).it(ind2).trained;
    Net_opt.X_train = fold(ind2).X_train; % normalized
    Net_opt.X_train_C = fold(ind2).X_train_C;
    Net_opt.X_train_S = fold(ind2).X_train_S;
    Net_opt.T_train = fold(ind2).T_train; % normalized
    Net_opt.T_train_C = fold(ind2).T_train_C;
    Net_opt.T_train_S = fold(ind2).T_train_S;
    Net_opt.X_val = fold(ind2).X_val; % normalized
    Net_opt.T_val = fold(ind2).T_val; % normalized
    Net_opt.X_test = (X_test-repmat(Net_opt.X_train_C,1,num_test))./repmat(Net_opt.X_train_S,1,num_test);
    Net_opt.T_test = (T_test-repmat(Net_opt.T_train_C,1,num_test))./repmat(Net_opt.T_train_S,1,num_test);

    save(filename,"X","T","Net_opt","kfold","layersizemax","FNNtype","trainFcn");

else

    load(filename);
    
end

%% Plotting
if doplots
    
    % Plot dependency of performance on size of hidden layers
    figure; tlo=tiledlayout(2,3,"TileSpacing","tight");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).J_train_mean],'filled');
    clim([0.05 0.5]); grid on; axis on;
    title("J training mean");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).J_train_min],'filled');
    clim([0.05 0.5]); grid on; axis on;
    title("J training min");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).J_train_max],'filled');
    clim([0.05 0.5]); grid on; axis on;
    title("J training max");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).J_val_mean],'filled');
    [Jmin,I]=min([Net(:).J_val_mean]);
    plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
    clim([0.05 0.5]); grid on; axis on;
    yticklabels("");
    title("J validation mean");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).J_val_min],'filled');
    [Jmin,I]=min([Net(:).J_val_mean]);
    plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
    clim([0.05 0.5]); grid on; axis on;
    yticklabels("");
    title("J validation min");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).J_val_max],'filled');
    [Jmin,I]=min([Net(:).J_val_mean]);
    plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
    clim([0.05 0.5]); grid on; axis on;
    yticklabels("");
    title("J validation max");
    
    xlabel(tlo,"nL1");
    ylabel(tlo,"nL2");
    
    cb=colorbar;
    cb.Layout.Tile='east';
end


function Net_opt = TrainANN(X,V,trainFcn,UseGPU,filename,doplots)

%% ----------------------------------------------------------------------- %%
% This script takes normalized data to train a simple feedforward neural 
% network (ANN) with 2 hidden layers. The input data is split into a 
% training set (80%) which is fed to the neural network, and a 
% cross-validation set (20%), which is used to find the optimal size of 
% each hidden layer, with the size of layer 1 (nL1) ranging from 2 to 10, 
% and the size of layer 2 (nL2) ranging from 2 to nL1. The performance is
% measured using a simple cost function that calculates the mean square
% error between ANN outputs and targets: J=1/m*sum((trained-targets)^2)
%% ----------------------------------------------------------------------- %%
% 
%% INPUTS: 
% X: n x m vector with n the number of training parameters and m the number
% of samples
% V: k x m vector with k the number of output parameters and m the number
% of samples
% UseGPU: 0 or 1
% trainFcn: training function (examples: 'trainlm', 'trainbr', 'trainscg')
% filename: where to save the network output
% doplots: 0 or 1
%
%% OUTPUTS:
% Net_opt: optimal network

warning("Make sure to normalize each row of the X and V vectors.");

addpath(getenv("froot_tools"));

if ~exist(filename,"file")

    for it=1:10 % split data into training and cross-validation segments in 10 different ways
    
        ntot = size(X,2);
        ind = randi(ntot,[1 round(ntot*0.2)]);
        XTrain = X; XTrain(:,ind(:))=[];
        XVal = X(:,ind(1,:));
        XTest = [];
        VTrain = V; VTrain(:,ind(:))=[];
        VVal = V(:,ind(1,:));
        VTest = [];
        
        nmax = 10;
        
        kk=1;
        
        for ii=2:nmax
        
            for jj=2:ii
        
                %net = feedforwardnet([ii jj],trainFcn);
                net = cascadeforwardnet([ii jj],trainFcn);
        
                % set early stopping parameters
                net.divideFcn= 'dividerand';
                net.divideParam.trainRatio = 0.8; % training set [%]
                net.divideParam.valRatio = 0; % validation set [%]
                net.divideParam.testRatio = 0.2; % test set [%]   
                %net.inputs{1}.processFcns = {'mapstd'}; % Normalize inputs/targets to have zero mean and unity variance
                net.trainParam.showWindow = 0;
                
                net = configure(net,XTrain,VTrain);
                
                if UseGPU==1
                    Net(kk).it(it).trained = train(net,XTrain,VTrain,'showResources','no','useGPU','only');
                else
                    Net(kk).it(it).trained = train(net,XTrain,VTrain,'showResources','no');
                end

                Net(kk).nL1 = ii;
                Net(kk).nL2 = jj;
        
                kk=kk+1;
        
            end
        end
        
        % Now calculate cost functions
        for kk=1:numel(Net)

            % performance on training data using simple quadratic misfit
            Y = Net(kk).it(it).trained(XTrain);
            Net(kk).JTrain(it) = 0.5/size(XTrain,2)*sum((VTrain(:)-Y(:)).^2);
        
            %Net(kk).Perf = perform(Net(kk).trained,Y,VTrain);
        
            Ptmp = polyfit(VTrain(:),Y(:),1);
            Net(kk).slope(it)=Ptmp(1);
            Net(kk).intercept(it)=Ptmp(2);
        
            [Rtmp,~]=corrcoef(VTrain(:),Y(:));
            Net(kk).R(it)=Rtmp(2,1);
        
            % performance on cross-validation data
            Y = Net(kk).it(it).trained(XVal);
            Net(kk).JVal(it) = 0.5/size(XVal,2)*sum((VVal(:)-Y(:)).^2);
        
            % performance on test data
            %Y = Net(kk).it(it).trained(XTest);
            %Net(kk).JTest(it) = 0.5/size(XTest,2)*sum((VTest(:)-Y(:)).^2);
        end

        fprintf("Done "+num2str(it)+"/10 iterations of the train/validate split of the data.\n")
    
    end
    

    for kk=1:numel(Net)
    
        JTrain = Net(kk).JTrain;
        JVal = Net(kk).JVal;
        %JTest = Net(kk).JTest;
    
        Net(kk).JTrain_min = min(JTrain);
        Net(kk).JTrain_max = max(JTrain);
        Net(kk).JTrain_mean = mean(JTrain);
        Net(kk).JTrain_std = std(JTrain);
    
        Net(kk).JVal_min = min(JVal);
        Net(kk).JVal_max = max(JVal);
        Net(kk).JVal_mean = mean(JVal);
        Net(kk).JVal_std = std(JVal);
    
        %Net(kk).JTest_min = min(JTest);
        %Net(kk).JTest_max = max(JTest);
        %Net(kk).JTest_mean = mean(JTest);
        %Net(kk).JTest_std = std(JTest);

    end

    save(filename,"Net");

else

    load(filename);
    
end

%% Choose 'optimal' architecture based on mean of cost function
[~,ind]=min([Net(:).JVal_mean]);
[~,ind2] = min(Net(ind).JVal);
Net_opt = Net(ind).it(ind2).trained;

%% Plotting
if doplots
    % Plot emulator vs data
    tlo=tiledlayout(1,2,"TileSpacing","tight");
    
    figure;
    plotregression(log10(I),log10(Y(1,:)));
    title('Misfit');
    
    figure;
    plotregression(log10(R),log10(Y(2,:)));
    title('Regularization');
    
    % Plot dependency of performance on size of hidden layers
    figure; tlo=tiledlayout(2,3,"TileSpacing","tight");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTrain_mean],'filled');
    clim([0.05 0.5]); grid on; axis on;
    title("J training mean");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTrain_min],'filled');
    clim([0.05 0.5]); grid on; axis on;
    title("J training min");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTrain_max],'filled');
    clim([0.05 0.5]); grid on; axis on;
    title("J training max");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JVal_mean],'filled');
    [Jmin,I]=min([Net(:).JVal_mean]);
    plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
    clim([0.05 0.5]); grid on; axis on;
    yticklabels("");
    title("J validation mean");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JVal_min],'filled');
    [Jmin,I]=min([Net(:).JVal_mean]);
    plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
    clim([0.05 0.5]); grid on; axis on;
    yticklabels("");
    title("J validation min");
    
    nexttile; hold on;
    scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JVal_max],'filled');
    [Jmin,I]=min([Net(:).JVal_mean]);
    plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
    clim([0.05 0.5]); grid on; axis on;
    yticklabels("");
    title("J validation max");
    
    xlabel(tlo,"nL1");
    ylabel(tlo,"nL2");
    
    cb=colorbar;
    cb.Layout.Tile='east';
end


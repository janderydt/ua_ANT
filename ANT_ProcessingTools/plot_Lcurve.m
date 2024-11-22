function [X,V]=plot_Lcurve

addpath("/mnt/md0/Ua/cases/ANT/");
addpath(getenv("froot_tools"));

UserVar.cycle = 1;

load("inversiondata_Weertman.mat");

ExpID = [data(:).InverseExpID];
gaA = [data(:).gaA];
gaC = [data(:).gaC];
gsA = [data(:).gsA];
gsC = [data(:).gsC];
m = [data(:).m];
n = [data(:).n];
for ii=1:numel(ExpID)
    I_tmp = data(ii).misfit;
    I(ii) = I_tmp(UserVar.cycle);
    R_tmp = data(ii).regularization;
    R(ii) = R_tmp(UserVar.cycle);
end

% 
% 
% %RunTable = "RunTable_ARCHER2_08-10-2024_"+string([14 15 16 17])+".csv";
% %UserVar.idrange = [14000,14999; 15000, 15999; 16000, 16999; 17000, 17999];
% 
% %RunTable = "RunTable_ARCHER2_"+string([3 6 9])+".csv";
% 
% gaA=[]; gaC=[]; gsA=[]; gsC=[]; m=[]; n=[]; finished=[]; error=[];
% Ind_finished=[]; Ind_error=[]; ExpID=[]; nIt=[]; nIt_tmp=[];
% 
% for tt=1:numel(RunTable)
% 
%     UserVar.home = pwd; 
%     UserVar.type = "Inverse";
%     UserVar.Table = UserVar.home+"/../ANT_"+UserVar.type+"/"+RunTable(tt);
% 
%     RunTable_tmp = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');
% 
%     gaA = [gaA; RunTable_tmp{:,'gaA'}];
%     gaC = [gaC; RunTable_tmp{:,'gaC'}];
%     gsA = [gsA; RunTable_tmp{:,'gsA'}];
%     gsC = [gsC; RunTable_tmp{:,'gsC'}];
%     nIt_tmp = RunTable_tmp{:,'InverseIterationsDone'};
%     nIt = [nIt; nIt_tmp];
%     m = [m; RunTable_tmp{:,'m'}];
%     n = [n; RunTable_tmp{:,'n'}];
%     finished_tmp = RunTable_tmp{:,'Finished'};
%     finished = [finished; finished_tmp];
%     error_tmp = RunTable_tmp{:,'Error'};
%     error = [error; error_tmp];
%     ExpID_tmp = RunTable_tmp{:,'ExpID'};
%     ExpID = [ExpID; ExpID_tmp];
% 
%     %Ind_finished_tmp = find(finished_tmp==1);
%     %Ind_finished = [Ind_finished; Ind_finished_tmp];
%     Ind_finished_tmp = find(nIt_tmp>=1000);
%     Ind_finished = [Ind_finished; Ind_finished_tmp];
%     Ind_error_tmp = find(error_tmp==1);
%     Ind_error = [Ind_error; Ind_error_tmp];
% 
% end
% 
% gaA = gaA(Ind_finished);
% gaC = gaC(Ind_finished);
% gsA = gsA(Ind_finished);
% gsC = gsC(Ind_finished);
% m = m(Ind_finished);
% n = n(Ind_finished);
% ExpID = ExpID(Ind_finished);

cm = crameri('roma',numel(ExpID));

% for i = 1:numel(ExpID)
% 
%     load("../ANT_Inverse/cases/ANT_nsmbl_Inverse_"+string(ExpID(i))+"/ANT_nsmbl_Inverse_"+string(ExpID(i))+...
%         "-RestartFile_InverseCycle"+string(UserVar.cycle)+".mat",...
%         "InvFinalValues");
% 
%     I(i) = InvFinalValues.I;
%     R(i) = InvFinalValues.R;
%     R_gsC(i) = InvFinalValues.R/(gsC(i)^2);
%     R_gsA(i) = InvFinalValues.R/(gsA(i)^2);
%     R_gaC(i) = InvFinalValues.R/(gaC(i)^2);
%     R_gaA(i) = InvFinalValues.R/(gaA(i)^2);
% 
%     %figure(111); hold on; 
%     %yyaxis left
%     %g(i)=plot(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-x','LineWidth',2,'color',cm(i,:));
%     %ylabel('J','interpreter','latex');
% 
% %     if ~isempty(RunInfo.Inverse.GradNorm)  && ~all(isnan(RunInfo.Inverse.GradNorm)) ...
% %             &&  numel(RunInfo.Inverse.Iterations) == numel(RunInfo.Inverse.GradNorm)
% % 
% %         hold off
% %         yyaxis right
% %         semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.GradNorm,'-r+')
% %         ylabel('Norm of gradient ($l_2$)','interpreter','latex')
% %         legend('Objective function','$\| \hbox{Gradient} \|$','Location','northeast','interpreter','latex')
% %         
% %     end
% 
% %     yyaxis left
%     %xlabel('Inverse iteration','interpreter','latex');
%     %hold off
% 
%     disp("Done "+string(i)+" out of "+string(numel(ExpID)));
% 
% end

%figure(111); legend(g(:),string(ysC))

% R_gsC = R./gsC;
% R_gsA = R./gsA;
% R_gaC = R./gaC;
% R_gaA = R./gaA;

% figure; scatter(log10(R_gsC),I,20,log10(gsC),'filled'); title("\gamma_sC"); grid on; box on; xlabel("log_{10}(R/\gamma_{sC}^2)"); ylabel("I");
% figure; scatter(log10(R_gsA),I,20,log10(gsA),'filled'); title("\gamma_sA"); grid on; box on; xlabel("log_{10}(R/\gamma_{sA}^2)"); ylabel("I");
% figure; scatter(log10(R_gaC),I,20,gaC,'filled'); title("\gamma_aC"); grid on; box on; xlabel("log_{10}(R/\gamma_{aC}^2)"); ylabel("I");
% figure; scatter(log10(R_gaA),I,20,gaA,'filled'); title("\gamma_aA"); grid on; box on; xlabel("log_{10}(R/\gamma_{aA}^2)"); ylabel("I");
% 

%% Prepare data
Ind = find(I>400);
gaA(Ind)=[]; gaC(Ind)=[]; gsA(Ind)=[]; gsC(Ind)=[];
R(Ind)=[]; I(Ind)=[];

[gaA_norm,gaA_C,gaA_S] = normalize(gaA); % gaA_norm = (gaA-gaA_C)/gaA_S
[gaC_norm,gaC_C,gaC_S] = normalize(gaC); 
[gsA_norm,gsA_C,gsA_S] = normalize(log10(gsA));
[gsC_norm,gsC_C,gsC_S] = normalize(log10(gsC));

[I_norm,I_C,I_S] = normalize(log10(I)); % log improves the performance of the emulator
[R_norm,R_C,R_S] = normalize(log10(R)); % log improves the performance of the emulator

Xtmp=[];
Xtmp(1,:)=gaA_norm(:)'; 
Xtmp(2,:)=gaC_norm(:)';
Xtmp(3,:)=gsA_norm(:)';
Xtmp(4,:)=gsC_norm(:)';
X=double(Xtmp);
%X = gpuArray(double(Xtmp));

Vtmp=[];
Vtmp(1,:) = I_norm(:)';
Vtmp(2,:) = R_norm(:)';
V=double(Vtmp);
%V = gpuArray(double(Vtmp));

%% Split data into training, test, validation sets to find optimal NN architecture
if ~exist("Lcurve_UNN_it10.mat")

    for it=1:10
    
        ntot = size(X,2);
        ind = randi(ntot,[2 round(ntot*0.15)]);
        XTrain = X; XTrain(:,ind(:))=[];
        XVal = X(:,ind(1,:));
        XTest = X(:,ind(2,:));
        VTrain = V; VTrain(:,ind(:))=[];
        VVal = V(:,ind(1,:));
        VTest = V(:,ind(2,:));
        
        %% UNN - find optimal values of 2-layer network
        nmax = 10;
        
        kk=1;
        
        for ii=2:nmax
        
            for jj=2:ii
        
                net = feedforwardnet([ii jj],'trainbr');
        
                % set early stopping parameters
                net.divideFcn= 'dividerand';
                net.divideParam.trainRatio = 0.7; % training set [%]
                net.divideParam.valRatio = 0.15; % validation set [%]
                net.divideParam.testRatio = 0.15; % test set [%]   
                %net.inputs{1}.processFcns = {'mapstd'}; % Normalize inputs/targets to have zero mean and unity variance
                net.trainParam.showWindow = 0;
                
                net = configure(net,XTrain,VTrain);
                
                Net(kk).it(it).trained = train(net,XTrain,VTrain,'showResources','no');%'useGPU','only');
                Net(kk).nL1 = ii;
                Net(kk).nL2 = jj;
        
                kk=kk+1;
        
                fprintf("ii="+num2str(ii)+", jj="+num2str(jj)+"\n");
                
                %figure; plotregression(V,Y,'Regression'); hold on;
                %plot([-1.5 3],[-1.5 3]);
        
            end
        end
        
        %% Now add cost functions
        for kk=1:numel(Net)
            % performance on training data
            Y = Net(kk).it(it).trained(XTrain);
            Net(kk).JTrain(it) = 0.5/size(XTrain,2)*sum((VTrain(:)-Y(:)).^2);
        
            %Net(kk).Perf = perform(Net(kk).trained,Y,VTrain);
        
            Ptmp = polyfit(VTrain(:),Y(:),1);
            Net(kk).slope(it)=Ptmp(1);
            Net(kk).intercept(it)=Ptmp(2);
        
            [Rtmp,~]=corrcoef(VTrain(:),Y(:));
            Net(kk).R(it)=Rtmp(2,1);
        
            % performance on validation data
            Y = Net(kk).it(it).trained(XVal);
            Net(kk).JVal(it) = 0.5/size(XVal,2)*sum((VVal(:)-Y(:)).^2);
        
            % performance on test data
            Y = Net(kk).it(it).trained(XTest);
            Net(kk).JTest(it) = 0.5/size(XTest,2)*sum((VTest(:)-Y(:)).^2);
        end

        fprintf("Done "+num2str(it)+" iterations of the train/validate/test split of the data.\n")
    
    end
    
    for kk=1:numel(Net)
    
        JTrain = Net(kk).JTrain;
        JVal = Net(kk).JVal;
        JTest = Net(kk).JTest;
    
        Net(kk).JTrain_min = min(JTrain);
        Net(kk).JTrain_max = max(JTrain);
        Net(kk).JTrain_mean = mean(JTrain);
        Net(kk).JTrain_std = std(JTrain);
    
        Net(kk).JVal_min = min(JVal);
        Net(kk).JVal_max = max(JVal);
        Net(kk).JVal_mean = mean(JVal);
        Net(kk).JVal_std = std(JVal);
    
        Net(kk).JTest_min = min(JTest);
        Net(kk).JTest_max = max(JTest);
        Net(kk).JTest_mean = mean(JTest);
        Net(kk).JTest_std = std(JTest);

    end

    save("Lcurve_UNN_it10.mat","Net");

else

    load("Lcurve_UNN_it10.mat");
    
end

%% Choose 'optimal' architecture
[Jmin,ind]=min([Net(:).JVal_mean]);
[JVal_opt,ind2] = min(Net(ind).JVal);
Net_opt = Net(ind).it(ind2).trained;

%% Emulate full dataset
Y = Net_opt(X); % 
Y(1,:) = 10.^(Y(1,:)*I_S + I_C); % undo normalization
Y(2,:) = 10.^(Y(2,:)*R_S + R_C); % undo normalization and log

%% Plot emulator vs data
tlo=tiledlayout(1,2,"TileSpacing","tight");

figure;
plotregression(I,Y(1,:));
title('Misfit');

figure;
plotregression(log10(R),log10(Y(2,:)));
title('Regularization');

%% Produce L-curve with new data
CM = parula(10);
for ga=1:10
    X_new(1,:) = ones(1,15); %gaA
    X_new(2,:) = ones(1,15); %gaC
    X_new(3,:) = 10.^linspace(3,6,15); %gsA
    X_new(4,:) = 1e4*ones(1,15); %gsC
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-gaA_C)/gaA_S;
    X_new_norm(2,:) = (X_new(2,:)-gaC_C)/gaC_S;
    X_new_norm(3,:) = (log10(X_new(3,:))-gsA_C)/gsA_S;
    X_new_norm(4,:) = (log10(X_new(4,:))-gsC_C)/gsC_S;
    
    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    figure(999); tlo=tiledlayout(1,4,'TileSpacing','tight');
    
    nexttile; hold on;
    plot(log10(Y_new(2,:)./X_new(3,:).^2),Y_new(1,:),'o-');
    title("gsA");

    %%
    X_new(1,:) = ones(1,15); %gaA
    X_new(2,:) = ones(1,15); %gaC
    X_new(3,:) = 1e4*ones(1,15); %gsA
    X_new(4,:) = 10.^linspace(3,6,15); %gsC
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-gaA_C)/gaA_S;
    X_new_norm(2,:) = (X_new(2,:)-gaC_C)/gaC_S;
    X_new_norm(3,:) = (log10(X_new(3,:))-gsA_C)/gsA_S;
    X_new_norm(4,:) = (log10(X_new(4,:))-gsC_C)/gsC_S;
    
    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    nexttile; hold on;
    plot(log10(Y_new(2,:)./X_new(4,:).^2),Y_new(1,:),'o-');
    title("gsC");

    %%
    X_new(1,:) = 10.^linspace(-1,2,15); %gaA
    X_new(2,:) = ones(1,15); %gaC
    X_new(3,:) = 1e4*ones(1,15); %gsA
    X_new(4,:) = 1e4*ones(1,15); %gsC
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-gaA_C)/gaA_S;
    X_new_norm(2,:) = (X_new(2,:)-gaC_C)/gaC_S;
    X_new_norm(3,:) = (log10(X_new(3,:))-gsA_C)/gsA_S;
    X_new_norm(4,:) = (log10(X_new(4,:))-gsC_C)/gsC_S;
    
    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    nexttile; hold on;
    plot(log10(Y_new(2,:)./X_new(1,:).^2),Y_new(1,:),'o-');
    title("gaA");

    %%
    X_new(1,:) = ones(1,15); %gaA
    X_new(2,:) = 10.^linspace(-1,2,15); %gaC
    X_new(3,:) = 1e4*ones(1,15); %gsA
    X_new(4,:) = 1e4*ones(1,15); %gsC
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-gaA_C)/gaA_S;
    X_new_norm(2,:) = (X_new(2,:)-gaC_C)/gaC_S;
    X_new_norm(3,:) = (log10(X_new(3,:))-gsA_C)/gsA_S;
    X_new_norm(4,:) = (log10(X_new(4,:))-gsC_C)/gsC_S;
    
    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    nexttile; hold on;
    plot(log10(Y_new(2,:)./X_new(2,:).^2),Y_new(1,:),'o-');
    title("gaC");

end


%% PLOTTING
figure; tlo=tiledlayout(3,3,"TileSpacing","tight");

nexttile; hold on;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTrain_mean],'filled');
clim([0.05 0.5]); grid on; axis on;
title("J training");

nexttile; hold on;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JVal_mean],'filled');
[Jmin,I]=min([Net(:).JVal_mean]);
plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
clim([0.05 0.5]); grid on; axis on;
yticklabels("");
title("J validation");

nexttile;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTest_mean],'filled');
clim([0.05 0.5]); grid on; axis on;
yticklabels("");
title("J test");

nexttile; hold on;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTrain_min],'filled');
clim([0.05 0.5]); grid on; axis on;

nexttile; hold on;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JVal_min],'filled');
[Jmin,I]=min([Net(:).JVal_mean]);
plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
clim([0.05 0.5]); grid on; axis on;
yticklabels("");

nexttile;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTest_min],'filled');
clim([0.05 0.5]); grid on; axis on;
yticklabels("");

nexttile; hold on;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTrain_max],'filled');
clim([0.05 0.5]); grid on; axis on;

nexttile; hold on;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JVal_max],'filled');
[Jmin,I]=min([Net(:).JVal_mean]);
plot(Net(I).nL1,Net(I).nL2,'ok','markersize',10);
clim([0.05 0.5]); grid on; axis on;
yticklabels("");

nexttile;
scatter([Net(:).nL1],[Net(:).nL2],50,[Net(:).JTest_max],'filled');
clim([0.05 0.5]); grid on; axis on;
yticklabels("");

xlabel(tlo,"nL1");
ylabel(tlo,"nL2");

cb=colorbar;
cb.Layout.Tile='east';
    

    %[Jmin,I]=min([Net(:).JVal])
    %hold on
    %plot(Net(I).nL1,Net(I).nL2],'ok','markersize',10);

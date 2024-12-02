function [U,S,V] = SVD_deltau

addpath(getenv("froot_tools"));

cycle=1;
ANNtype = "feedforwardnet"; % feedforwardnet or cascadeforwardnet, a simple feedforwardnet seems to perform just fine
trainFcn = "trainlm"; % trainlm is fast on CPU and seems to perform just fine
UseGPU=0;
doplots=1;

% load data file
load("Delta_u_AS_Weertman.mat");

%% PREDICTORS
X = [m(:) n(:) gaA(:) gaC(:) log10(gsA(:)) log10(gsC(:))];
if cycle==2 % add dhdt_err to predictors
    X = [X log10(dhdt_err(:))];
end
X = double(X);

%% TARGETS (training/input data). In this case, input data consists of 
%% simulated instanteneous changes in surface speed in response to changes
%% in ice-sheet geometry (ice thickness, calving front location). The input
%% data consists of num_nodes nodal values for nun_exp experiments.
T = Delta_u.Calv_dh.map(:,:,cycle);
% check for nans and inf
if any(isnan(T) | isinf(T))
    error("Remove nan and inf from T");
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
MUA=MUA_2018;

%% DIMENSIONALITY REDUCTION: use singular value decomposition of targets to
%% express each map of \Delta u in terms of its k dominant principal components
%% do the full svd to calculate how many components are required for each
%% increment in pct
[~,S,~] = svd(T,'econ');

seq = randperm(num_exp);
pct = 0.99;
T_nComp = find((cumsum(diag(S).^2)./sum(diag(S).^2))>pct,1,'first');
T_pct1 = diag(S).^2./sum(diag(S).^2);
T_pct = T_pct1(1:T_nComp);

% return the left singular vectors U, diagonal matrix S of singular 
% values, and right singular vectors V for T_nComp largest singular values
[U_trunc,S_trunc,V_trunc] = svds(T,T_nComp);

B_trunc = S_trunc*V_trunc';
T_reproj = (B_trunc*B_trunc')\B_trunc; % equivalent to inv(B*B')*B
T_hat = T*T_reproj';

data = T_hat(seq,:);
predictors = X(seq,:);

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

%% WRITE OUTPUTS for tensorflow
fname1 = sprintf('./mat_files/data_N0k%.2g',pct*100);
fname2 = sprintf('./mat_files/SVD_N0k%.2g',pct*100);
save(fname1,'X_test','X_val','X_train','T_train','T_val','T_test');
save(fname2, 'V_trunc', 'S_trunc', 'B_trunc', 'T_reproj', 'T_pct','seq');

%% Now simulate FeedForward NN
filename = "Perturbation_Calv_dh_UNN_cycle"+num2str(cycle)+"_"+ANNtype+"_"+trainFcn+".mat";
Net = TrainANN(predictors',data',ANNtype,trainFcn,UseGPU,filename,doplots);

%% Emulate full dataset
X_full = [Net.X_train Net.X_val Net.X_test];
Y = Net.trained(X_full);
T_full = [Net.T_train Net.T_val Net.T_test];

%% Plot emulator vs targets
for ii=1:size(X_full,1)
    figure;
    plotregression(T_full(ii,:),Y(ii,:));
    title("mode "+num2str(ii));
end

return

%% PLOTTING

if doplots

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
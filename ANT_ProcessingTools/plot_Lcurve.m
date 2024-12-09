function plot_Lcurve

%% emulate L-surface for predictors (m,n,gaA,gaC,gsA,gsC,dhdt) using a 
%% feedforward neural network and plot some example L-curves

addpath("/mnt/md0/Ua/cases/ANT/");
addpath(getenv("froot_tools"));

%% Inverse cycle: either 1 (no spin-up) or 2 (with spin-up)
cycle = 1;
FNNtype = "feedforwardnet"; % feedforwardnet or cascadeforwardnet, a simple feedforwardnet seems to perform just fine
trainFcn = "trainscg"; % trainlm is fast on CPU and seems to perform just fine
UseGPU = 1;

%% Load data - need to run Calc_InverseResults_Ensemble first
load("inversiondata_Weertman.mat");

ExpID = [data(:).InverseExpID];
gaA = [data(:).gaA];
gaC = [data(:).gaC];
gsA = [data(:).gsA];
gsC = [data(:).gsC];
m = [data(:).m];
n = [data(:).n];
dhdt_err = [data(:).dhdt_err];

for ii=1:numel(ExpID)
    I_tmp = data(ii).misfit;
    I(ii) = I_tmp(cycle);
    R_tmp = data(ii).regularization;
    R(ii) = R_tmp(cycle);
end

%% Prepare data for training 
% remove outliers for misfit
if cycle == 1
    Ind = find(I>400 | I==0);
else
    Ind = find(I>400 | I==0 | isnan(dhdt_err));
end

% predictors
m(Ind)=[]; n(Ind)=[];
gaA(Ind)=[]; gaC(Ind)=[]; gsA(Ind)=[]; gsC(Ind)=[];
R(Ind)=[]; I(Ind)=[]; dhdt_err(Ind)=[];

X = [m(:) n(:) gaA(:) gaC(:) log10(gsA(:)) log10(gsC(:))];
if cycle==2 % add dhdt_err to predictors
    X = [X log10(dhdt_err(:))];
end
X = double(X');

% targets
T = [log10(I(:)) log10(R(:))];
T = double(T');

%% Start training
filename = "./FNN/mat_files/Lcurve_FNN_cycle"+num2str(cycle)+"_"+FNNtype+"_"+trainFcn+".mat";
doplots = 0;
addpath("FNN/");
Net = TrainFNN(X,T,FNNtype,trainFcn,UseGPU,filename,doplots);

%% Emulate full dataset
X_full = [Net.X_train Net.X_val Net.X_test];
Y = Net.trained(X_full);
T_full = [Net.T_train Net.T_val Net.T_test];

%% Plot emulator vs targets
% figure;
% plotregression(T_full(1,:),Y(1,:));
% title('Misfit');
% 
% figure;
% plotregression(T_full(2,:),Y(2,:));
% title('Regularization');

%% Produce L-curve for new data
CM = parula(10);

if cycle==1
    tiles=4;     
    axestag=["gaA","gaC","gsA","gsC"];
else
    tiles=5;
    axestag=["gaA","gaC","gsA","gsC","dhdt"];
end

fh=findobj('Type','Figure','Name','L-curves'); 

if isempty(fh)
    fh=fig('units','inches','width',120*12/72.27,'height',45*12/72.27,...
        'Name','L-curves','fontsize',14,'font','Helvetica');
    tlo=tiledlayout(fh,1,tiles,'TileSpacing','tight');
    for tt=1:tiles
        ax(tt)=nexttile(tlo); hold on;
        ax(tt).Tag=axestag(tt);
    end
else
    for tt=1:numel(axestag)
        ax(tt)=findobj(fh,'type','axes','tag',axestag(tt));
    end
end

%% gaA
X_new(1,:) = 3*ones(1,15); %m
X_new(2,:) = 3*ones(1,15); %n
X_new(3,:) = 10.^linspace(-1,2,15); %gaA
X_new(4,:) = ones(1,15); %gaC
X_new(5,:) = log10(1e4*ones(1,15)); %gsA
X_new(6,:) = log10(1e4*ones(1,15)); %gsC
if cycle>1
    X_new(7,:)=log10(0.5*ones(1,15)); %dhdt_err
end

% find nearest in original dataset
Ind=knnsearch(X',X_new','K',1); % X contains non-normalized values
X_orig = X(:,Ind);
T_orig = T(:,Ind);

% apply normalization and predict using trained network
X_new_norm = (X_new-repmat(Net.X_train_C,1,15))./repmat(Net.X_train_S,1,15);
Y_new = Net.trained(X_new_norm);

% undo normalization and log
for ii=1:size(Y_new,1)
    Y_new(ii,:) = 10.^(Y_new(ii,:)*Net.T_train_S(ii) + Net.T_train_C(ii)); 
end

% undo log of original predictor and target data
X_orig(5,:)=10.^X_orig(5,:);
X_orig(6,:)=10.^X_orig(6,:);
if cycle>1
    X_orig(7,:)=10.^X_orig(7,:);
end
T_orig(1,:) = 10.^T_orig(1,:);
T_orig(2,:) = 10.^T_orig(2,:);

% plot
x=log10(Y_new(2,:)./X_new(3,:).^2);
y=Y_new(1,:);
plot(ax(1),x,y,'+-','LineWidth',2,'MarkerSize',8);
for kk=1:2:size(Y_new,2)
    text(ax(1),x(kk),y(kk)+2,erase(num2str(X_new(3,kk),'%3.2e'),'+0'),'FontSize',8);
end
plot(ax(1),log10(T_orig(2,:)./X_orig(3,:).^2),T_orig(1,:),'xr');
title(ax(1),"gaA"); grid(ax(1),"on"); box(ax(1),"on");

%% gaC
X_new(1,:) = 3*ones(1,15); %m
X_new(2,:) = 3*ones(1,15); %n
X_new(3,:) = 10*ones(1,15); %gaA
X_new(4,:) = 10.^linspace(-1,2,15); %gaC
X_new(5,:) = log10(1e4*ones(1,15)); %gsA
X_new(6,:) = log10(1e4*ones(1,15)); %gsC
if cycle>1
    X_new(7,:)=log10(0.5*ones(1,15)); %dhdt_err
end

% find nearest in original dataset
Ind=knnsearch(X',X_new','K',1); % X contains non-normalized values
X_orig = X(:,Ind);
T_orig = T(:,Ind);

% apply normalization and predict using trained network
X_new_norm = (X_new-repmat(Net.X_train_C,1,15))./repmat(Net.X_train_S,1,15);
Y_new = Net.trained(X_new_norm);

% undo normalization and log
for ii=1:size(Y_new,1)
    Y_new(ii,:) = 10.^(Y_new(ii,:)*Net.T_train_S(ii) + Net.T_train_C(ii)); 
end

% undo log of original predictor and target data
X_orig(5,:)=10.^X_orig(5,:);
X_orig(6,:)=10.^X_orig(6,:);
if cycle>1
    X_orig(7,:)=10.^X_orig(7,:);
end
T_orig(1,:) = 10.^T_orig(1,:);
T_orig(2,:) = 10.^T_orig(2,:);
        
x=log10(Y_new(2,:)./X_new(4,:).^2);
y=Y_new(1,:);
plot(ax(2),x,y,'+-','LineWidth',2,'MarkerSize',8);
for kk=1:2:size(Y_new,2)
    text(ax(2),x(kk),y(kk)+2,erase(num2str(X_new(4,kk),'%3.2e'),'+0'),'FontSize',8);
end
plot(ax(2),log10(T_orig(2,:)./X_orig(4,:).^2),T_orig(1,:),'xr');
title(ax(2),"gaC"); grid(ax(2),"on"); box(ax(2),"on");

%% gsA
X_new(1,:) = 3*ones(1,15); %m
X_new(2,:) = 3*ones(1,15); %n
X_new(3,:) = 10*ones(1,15); %gaA
X_new(4,:) = 10*ones(1,15); %gaC
X_new(5,:) = log10(10.^linspace(3,6,15)); %gsA
X_new(6,:) = log10(1e4*ones(1,15)); %gsC
if cycle>1
    X_new(7,:)=log10(0.5*ones(1,15)); %dhdt_err
end

% find nearest in original dataset
Ind=knnsearch(X',X_new','K',1); % X contains non-normalized values
X_orig = X(:,Ind);
T_orig = T(:,Ind);

% apply normalization and predict using trained network
X_new_norm = (X_new-repmat(Net.X_train_C,1,15))./repmat(Net.X_train_S,1,15);
Y_new = Net.trained(X_new_norm);

% undo normalization and log
for ii=1:size(Y_new,1)
    Y_new(ii,:) = 10.^(Y_new(ii,:)*Net.T_train_S(ii) + Net.T_train_C(ii)); 
end

% undo log of original predictor and target data
X_orig(5,:)=10.^X_orig(5,:); X_new(5,:)=10.^X_new(5,:);
X_orig(6,:)=10.^X_orig(6,:); X_new(6,:)=10.^X_new(6,:);
if cycle>1
    X_orig(7,:)=10.^X_orig(7,:); X_new(7,:)=10.^X_new(7,:);
end
T_orig(1,:) = 10.^T_orig(1,:);
T_orig(2,:) = 10.^T_orig(2,:);

x=log10(Y_new(2,:)./X_new(5,:).^2);
y=Y_new(1,:);
plot(ax(3),x,y,'+-','LineWidth',2,'MarkerSize',8);
for kk=1:2:size(Y_new,2)
    text(ax(3),x(kk),y(kk)+2,erase(num2str(X_new(5,kk),'%3.2e'),'+0'),'FontSize',8);
end
plot(ax(3),log10(T_orig(2,:)./X_orig(5,:).^2),T_orig(1,:),'xr');
title(ax(3),"gsA"); grid(ax(3),"on"); box(ax(3),"on");

%% gsC
X_new(1,:) = 3*ones(1,15); %m
X_new(2,:) = 3*ones(1,15); %n
X_new(3,:) = 10*ones(1,15); %gaA
X_new(4,:) = 10*ones(1,15); %gaC
X_new(5,:) = log10(1e4*ones(1,15)); %gsA
X_new(6,:) = log10(10.^linspace(3,6,15)); %gsC
if cycle>1
    X_new(7,:)=log10(0.5*ones(1,15)); %dhdt_err
end

% find nearest in original dataset
Ind=knnsearch(X',X_new','K',1); % X contains non-normalized values
X_orig = X(:,Ind);
T_orig = T(:,Ind);

% apply normalization and predict using trained network
X_new_norm = (X_new-repmat(Net.X_train_C,1,15))./repmat(Net.X_train_S,1,15);
Y_new = Net.trained(X_new_norm);

% undo normalization and log
for ii=1:size(Y_new,1)
    Y_new(ii,:) = 10.^(Y_new(ii,:)*Net.T_train_S(ii) + Net.T_train_C(ii)); 
end

% undo log of original predictor and target data
X_orig(5,:)=10.^X_orig(5,:); X_new(5,:)=10.^X_new(5,:);
X_orig(6,:)=10.^X_orig(6,:); X_new(6,:)=10.^X_new(6,:);
if cycle>1
    X_orig(7,:)=10.^X_orig(7,:); X_new(7,:)=10.^X_new(7,:);
end
T_orig(1,:) = 10.^T_orig(1,:);
T_orig(2,:) = 10.^T_orig(2,:);

x=log10(Y_new(2,:)./X_new(6,:).^2);
y=Y_new(1,:);
plot(ax(4),x,y,'+-','LineWidth',2,'MarkerSize',8);
for kk=1:2:size(Y_new,2)
    text(ax(4),x(kk),y(kk)+2,erase(num2str(X_new(6,kk),'%3.2e'),'+0'),'FontSize',8);
end
plot(ax(4),log10(T_orig(2,:)./X_orig(6,:).^2),T_orig(1,:),'xr');
title(ax(4),"gsC"); grid(ax(4),"on"); box(ax(4),"on");

%% dhdt
if cycle>1
    X_new(1,:) = 3*ones(1,15); %m
    X_new(2,:) = 3*ones(1,15); %n
    X_new(3,:) = 10*ones(1,15); %gaA
    X_new(4,:) = 10*ones(1,15); %gaC
    X_new(5,:) = log10(1e4*ones(1,15)); %gsA
    X_new(6,:) = log10(1e4*ones(1,15)); %gsC
    X_new(7,:) = log10(linspace(0.1,0.5,15)); %dhdt_err
    
    % find nearest in original dataset
    Ind=knnsearch(X',X_new','K',1); % X contains non-normalized values
    X_orig = X(:,Ind);
    T_orig = T(:,Ind);
    
    % apply normalization and predict using trained network
    X_new_norm = (X_new-repmat(Net.X_train_C,1,15))./repmat(Net.X_train_S,1,15);
    Y_new = Net.trained(X_new_norm);
    
    % undo normalization and log
    for ii=1:size(Y_new,1)
        Y_new(ii,:) = 10.^(Y_new(ii,:)*Net.T_train_S(ii) + Net.T_train_C(ii)); 
    end
    
    % undo log of original predictor and target data
    X_orig(5,:)=10.^X_orig(5,:); X_new(5,:)=10.^X_new(5,:);
    X_orig(6,:)=10.^X_orig(6,:); X_new(6,:)=10.^X_new(6,:);
    if cycle>1
        X_orig(7,:)=10.^X_orig(7,:); X_new(7,:)=10.^X_new(7,:);
    end
    T_orig(1,:) = 10.^T_orig(1,:);
    T_orig(2,:) = 10.^T_orig(2,:);

    x=log10(Y_new(2,:));
    y=Y_new(1,:);
    plot(ax(5),x,y,'+-');
    for kk=1:2:size(Y_new,2)
        text(ax(5),x(kk),y(kk)+4,erase(num2str(X_new(7,kk),'%3.2e'),'+0'));
    end
    plot(ax(5),log10(Y_orig(2,:)),Y_orig(1,:),'xr');
    %plot(ax(5),log10(Y_orig(2,:)),Y_orig(1,:),'xr');
    title(ax(5),"dhdt_err"); grid(ax(5),"on"); box(ax(5),"on");
end

xlabel(tlo,'Misfit');
ylabel(tlo,'Regularization');

pos = get(fh,"Position");
set(fh,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Lcurve";
print(fh,fname,"-dpng","-r400");
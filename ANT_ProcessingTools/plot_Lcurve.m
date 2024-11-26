function [X,V]=plot_Lcurve

addpath("/mnt/md0/Ua/cases/ANT/");
addpath(getenv("froot_tools"));

%% Inverse cycle: either 1 (no spin-up) or 2 (with spin-up)
UserVar.cycle = 1;
trainFcn = "trainlm";
UseGPU = 0;

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
    I(ii) = I_tmp(UserVar.cycle);
    R_tmp = data(ii).regularization;
    R(ii) = R_tmp(UserVar.cycle);
end

%% Prepare data for training (remove outliers for misfit and normalize)
if UserVar.cycle == 1
    Ind = find(I>400 | I==0);
else
    Ind = find(I>400 | I==0 | isnan(dhdt_err));
end
m(Ind)=[]; n(Ind)=[];
gaA(Ind)=[]; gaC(Ind)=[]; gsA(Ind)=[]; gsC(Ind)=[];
R(Ind)=[]; I(Ind)=[]; dhdt_err(Ind)=[];

% input training parameters
[m_norm,m_C,m_S] = normalize(m);
[n_norm,n_C,n_S] = normalize(n);
[gaA_norm,gaA_C,gaA_S] = normalize(gaA); % gaA_norm = (gaA-gaA_C)/gaA_S
[gaC_norm,gaC_C,gaC_S] = normalize(gaC); 
[gsA_norm,gsA_C,gsA_S] = normalize(log10(gsA));
[gsC_norm,gsC_C,gsC_S] = normalize(log10(gsC));

Xtmp=[];
Xtmp(1,:)=m_norm(:)';
Xtmp(2,:)=n_norm(:)';
Xtmp(3,:)=gaA_norm(:)'; 
Xtmp(4,:)=gaC_norm(:)';
Xtmp(5,:)=gsA_norm(:)';
Xtmp(6,:)=gsC_norm(:)';

if UserVar.cycle>1
    [dhdt_err_norm,dhdt_err_C,dhdt_err_S] = normalize(log10(dhdt_err));
    Xtmp(7,:)=dhdt_err_norm(:)';
end

X=double(Xtmp);

% targets
[I_norm,I_C,I_S] = normalize(log10(I)); % log improves the performance of the emulator
[R_norm,R_C,R_S] = normalize(log10(R)); % log improves the performance of the emulator

Vtmp=[];
Vtmp(1,:) = I_norm(:)';
Vtmp(2,:) = R_norm(:)';
V=double(Vtmp);

if UseGPU
    X = gpuArray(X);
    V = gpuArray(V);
end

%% Start training
filename = "Lcurve_UNN_cycle"+num2str(UserVar.cycle)+"_cascadeforwardnet_"+trainFcn+".mat";
doplots = 1;
Net = TrainANN(X,V,trainFcn,UseGPU,filename,doplots);

%% Emulate full dataset
Y = Net(X); %
if UseGPU
    X=gather(X);
    Y=gather(Y);
end
Y(1,:) = 10.^(Y(1,:)*I_S + I_C); % undo normalization
Y(2,:) = 10.^(Y(2,:)*R_S + R_C); % undo normalization and log

%% Produce L-curve for new data
CM = parula(10);

if UserVar.cycle==1
    tiles=4;       
else
    tiles=5;
end
figure(999); tlo=tiledlayout(1,tiles,'TileSpacing','tight');
for tt=1:tiles
    ax(tt)=nexttile(tlo); hold on;
end

for mexp=3

    %% gaA
    X_new(1,:) = mexp*ones(1,15); %m
    X_new(2,:) = 3*ones(1,15); %n
    X_new(3,:) = 10.^linspace(-1,2,15); %gaA
    X_new(4,:) = ones(1,15); %gaC
    X_new(5,:) = 1e4*ones(1,15); %gsA
    X_new(6,:) = 1e4*ones(1,15); %gsC
    if UserVar.cycle>1
        X_new(7,:)=0.5*ones(1,15); %dhdt_err
    end
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-m_C)/m_S;
    X_new_norm(2,:) = (X_new(2,:)-n_C)/n_S;
    X_new_norm(3,:) = (X_new(3,:)-gaA_C)/gaA_S;
    X_new_norm(4,:) = (X_new(4,:)-gaC_C)/gaC_S;
    X_new_norm(5,:) = (log10(X_new(5,:))-gsA_C)/gsA_S;
    X_new_norm(6,:) = (log10(X_new(6,:))-gsC_C)/gsC_S;
    if UserVar.cycle>1
        X_new_norm(7,:)=(log10(X_new(7,:))-dhdt_err_C)/dhdt_err_S;
    end

    % find nearest to original dataset
    Ind=knnsearch(X',X_new_norm','K',1);
    Y_orig = V(:,Ind);
    Y_orig(1,:) = 10.^(Y_orig(1,:)*I_S + I_C); % undo normalization
    Y_orig(2,:) = 10.^(Y_orig(2,:)*R_S + R_C); % undo normalization and log
    X_orig = X(:,Ind);
    X_orig(1,:) = X_orig(1,:)*m_S + m_C; % undo normalization
    X_orig(2,:) = X_orig(2,:)*n_S + n_C; % undo normalization and log
    X_orig(3,:) = X_orig(3,:)*gaA_S + gaA_C; % undo normalization
    X_orig(4,:) = X_orig(4,:)*gaC_S + gaC_C; % undo normalization and log
    X_orig(5,:) = 10.^(X_orig(5,:)*gsA_S + gsA_C); % undo normalization
    X_orig(6,:) = 10.^(X_orig(6,:)*gsC_S + gsC_C); % undo normalization and log
    if UserVar.cycle>1
        X_orig(7,:)=10.^(X_orig(7,:)*dhdt_err_S+dhdt_err_C); % undo normalization and log
    end

    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log

    x=log10(Y_new(2,:)./X_new(3,:).^2);
    y=Y_new(1,:);
    plot(ax(1),x,y,'+-');
    if mexp==3
        for kk=1:2:size(Y_new,2)
            text(ax(1),x(kk),y(kk)+4,erase(num2str(X_new(3,kk),'%3.2e'),'+0'));
        end
    end
    plot(ax(1),log10(Y_orig(2,:)./X_orig(3,:).^2),Y_orig(1,:),'xr');
    title(ax(1),"gaA"); grid(ax(1),"on"); box(ax(1),"on");

    %% gaC
    X_new(1,:) = mexp*ones(1,15); %m
    X_new(2,:) = 3*ones(1,15); %n
    X_new(3,:) = ones(1,15); %gaA
    X_new(4,:) = 10.^linspace(-1,2,15); %gaC
    X_new(5,:) = 1e4*ones(1,15); %gsA
    X_new(6,:) = 1e4*ones(1,15); %gsC
    if UserVar.cycle>1
        X_new(7,:)=0.5*ones(1,15); %dhdt_err
    end
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-m_C)/m_S;
    X_new_norm(2,:) = (X_new(2,:)-n_C)/n_S;
    X_new_norm(3,:) = (X_new(3,:)-gaA_C)/gaA_S;
    X_new_norm(4,:) = (X_new(4,:)-gaC_C)/gaC_S;
    X_new_norm(5,:) = (log10(X_new(5,:))-gsA_C)/gsA_S;
    X_new_norm(6,:) = (log10(X_new(6,:))-gsC_C)/gsC_S;
    if UserVar.cycle>1
        X_new_norm(7,:)=(log10(X_new(7,:))-dhdt_err_C)/dhdt_err_S;
    end

    % find nearest to original dataset
    Ind=knnsearch(X',X_new_norm','K',1);
    Y_orig = V(:,Ind);
    Y_orig(1,:) = 10.^(Y_orig(1,:)*I_S + I_C); % undo normalization
    Y_orig(2,:) = 10.^(Y_orig(2,:)*R_S + R_C); % undo normalization and log
    X_orig = X(:,Ind);
    X_orig(1,:) = X_orig(1,:)*m_S + m_C; % undo normalization
    X_orig(2,:) = X_orig(2,:)*n_S + n_C; % undo normalization and log
    X_orig(3,:) = X_orig(3,:)*gaA_S + gaA_C; % undo normalization
    X_orig(4,:) = X_orig(4,:)*gaC_S + gaC_C; % undo normalization and log
    X_orig(5,:) = 10.^(X_orig(5,:)*gsA_S + gsA_C); % undo normalization
    X_orig(6,:) = 10.^(X_orig(6,:)*gsC_S + gsC_C); % undo normalization and log
    if UserVar.cycle>1
        X_orig(7,:) = 10.^(X_orig(7,:)*dhdt_err_S+dhdt_err_C); % undo normalization and log
    end

    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    x=log10(Y_new(2,:)./X_new(4,:).^2);
    y=Y_new(1,:);
    plot(ax(2),x,y,'+-');
    if mexp==3
        for kk=1:2:size(Y_new,2)
            text(ax(2),x(kk),y(kk)+4,erase(num2str(X_new(4,kk),'%3.2e'),'+0'));
        end
    end
    plot(ax(2),log10(Y_orig(2,:)./X_orig(4,:).^2),Y_orig(1,:),'xr');
    title(ax(2),"gaC"); grid(ax(2),"on"); box(ax(2),"on");

    %% gsA
    X_new(1,:) = mexp*ones(1,15); %m
    X_new(2,:) = 3*ones(1,15); %n
    X_new(3,:) = ones(1,15); %gaA
    X_new(4,:) = ones(1,15); %gaC
    X_new(5,:) = 10.^linspace(3,6,15); %gsA
    X_new(6,:) = 1e4*ones(1,15); %gsC
    if UserVar.cycle>1
        X_new(7,:)=0.5*ones(1,15); %dhdt_err
    end
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-m_C)/m_S;
    X_new_norm(2,:) = (X_new(2,:)-n_C)/n_S;
    X_new_norm(3,:) = (X_new(3,:)-gaA_C)/gaA_S;
    X_new_norm(4,:) = (X_new(4,:)-gaC_C)/gaC_S;
    X_new_norm(5,:) = (log10(X_new(5,:))-gsA_C)/gsA_S;
    X_new_norm(6,:) = (log10(X_new(6,:))-gsC_C)/gsC_S;
    if UserVar.cycle>1
        X_new_norm(7,:)=(log10(X_new(7,:))-dhdt_err_C)/dhdt_err_S;
    end

    % find nearest to original dataset
    Ind=knnsearch(X',X_new_norm','K',1);
    Y_orig = V(:,Ind);
    Y_orig(1,:) = 10.^(Y_orig(1,:)*I_S + I_C); % undo normalization
    Y_orig(2,:) = 10.^(Y_orig(2,:)*R_S + R_C); % undo normalization and log
    X_orig = X(:,Ind);
    X_orig(1,:) = X_orig(1,:)*m_S + m_C; % undo normalization
    X_orig(2,:) = X_orig(2,:)*n_S + n_C; % undo normalization and log
    X_orig(3,:) = X_orig(3,:)*gaA_S + gaA_C; % undo normalization
    X_orig(4,:) = X_orig(4,:)*gaC_S + gaC_C; % undo normalization and log
    X_orig(5,:) = 10.^(X_orig(5,:)*gsA_S + gsA_C); % undo normalization
    X_orig(6,:) = 10.^(X_orig(6,:)*gsC_S + gsC_C); % undo normalization and log
    if UserVar.cycle>1
        X_orig(7,:)=10.^(X_orig(7,:)*dhdt_err_S+dhdt_err_C); % undo normalization and log
    end

    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    x=log10(Y_new(2,:)./X_new(5,:).^2);
    y=Y_new(1,:);
    plot(ax(3),x,y,'+-');
    if mexp==3
        for kk=1:2:size(Y_new,2)
            text(ax(3),x(kk),y(kk)+4,erase(num2str(X_new(5,kk),'%3.2e'),'+0'));
        end
    end
    plot(ax(3),log10(Y_orig(2,:)./X_orig(5,:).^2),Y_orig(1,:),'xr');
    title(ax(3),"gsA"); grid(ax(3),"on"); box(ax(3),"on");

    %% gsC
    X_new(1,:) = mexp*ones(1,15); %m
    X_new(2,:) = 3*ones(1,15); %n
    X_new(3,:) = ones(1,15); %gaA
    X_new(4,:) = ones(1,15); %gaC
    X_new(5,:) = 1e4*ones(1,15); %gsA
    X_new(6,:) = 10.^linspace(3,6,15); %gsC
    if UserVar.cycle>1
        X_new(7,:)=0.5*ones(1,15); %dhdt_err
    end
    
    % apply normalization
    X_new_norm(1,:) = (X_new(1,:)-m_C)/m_S;
    X_new_norm(2,:) = (X_new(2,:)-n_C)/n_S;
    X_new_norm(3,:) = (X_new(3,:)-gaA_C)/gaA_S;
    X_new_norm(4,:) = (X_new(4,:)-gaC_C)/gaC_S;
    X_new_norm(5,:) = (log10(X_new(5,:))-gsA_C)/gsA_S;
    X_new_norm(6,:) = (log10(X_new(6,:))-gsC_C)/gsC_S;
    if UserVar.cycle>1
        X_new_norm(7,:)=(log10(X_new(7,:))-dhdt_err_C)/dhdt_err_S;
    end

    % find nearest to original dataset
    Ind=knnsearch(X',X_new_norm','K',1);
    Y_orig = V(:,Ind);
    Y_orig(1,:) = 10.^(Y_orig(1,:)*I_S + I_C); % undo normalization
    Y_orig(2,:) = 10.^(Y_orig(2,:)*R_S + R_C); % undo normalization and log
    X_orig = X(:,Ind);
    X_orig(1,:) = X_orig(1,:)*m_S + m_C; % undo normalization
    X_orig(2,:) = X_orig(2,:)*n_S + n_C; % undo normalization and log
    X_orig(3,:) = X_orig(3,:)*gaA_S + gaA_C; % undo normalization
    X_orig(4,:) = X_orig(4,:)*gaC_S + gaC_C; % undo normalization and log
    X_orig(5,:) = 10.^(X_orig(5,:)*gsA_S + gsA_C); % undo normalization
    X_orig(6,:) = 10.^(X_orig(6,:)*gsC_S + gsC_C); % undo normalization and log
    if UserVar.cycle>1
        X_orig(7,:)=10.^(X_orig(7,:)*dhdt_err_S+dhdt_err_C); % undo normalization and log
    end

    Y_new = Net_opt(X_new_norm);
    Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
    Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
    
    x=log10(Y_new(2,:)./X_new(6,:).^2);
    y=Y_new(1,:);
    plot(ax(4),x,y,'+-');
    if mexp==3
        for kk=1:2:size(Y_new,2)
            text(ax(4),x(kk),y(kk)+4,erase(num2str(X_new(6,kk),'%3.2e'),'+0'));
        end
    end
    plot(ax(4),log10(Y_orig(2,:)./X_orig(6,:).^2),Y_orig(1,:),'xr');
    title(ax(4),"gsC"); grid(ax(4),"on"); box(ax(4),"on");

    %% dhdt
    if UserVar.cycle>1
        X_new(1,:) = mexp*ones(1,15); %m
        X_new(2,:) = 3*ones(1,15); %n
        X_new(3,:) = ones(1,15); %gaA
        X_new(4,:) = ones(1,15); %gaC
        X_new(5,:) = 1e4*ones(1,15); %gsA
        X_new(6,:) = 1e4*ones(1,15); %gsC
        X_new(7,:) = linspace(0.1,0.5,15); %dhdt_err
        
        % apply normalization
        X_new_norm(1,:) = (X_new(1,:)-m_C)/m_S;
        X_new_norm(2,:) = (X_new(2,:)-n_C)/n_S;
        X_new_norm(3,:) = (X_new(3,:)-gaA_C)/gaA_S;
        X_new_norm(4,:) = (X_new(4,:)-gaC_C)/gaC_S;
        X_new_norm(5,:) = (log10(X_new(5,:))-gsA_C)/gsA_S;
        X_new_norm(6,:) = (log10(X_new(6,:))-gsC_C)/gsC_S;
        X_new_norm(7,:) = (log10(X_new(7,:))-dhdt_err_C)/dhdt_err_S;
    
        % find nearest to original dataset
        Ind=knnsearch(X',X_new_norm','K',1);
        Y_orig = V(:,Ind);
        Y_orig(1,:) = 10.^(Y_orig(1,:)*I_S + I_C); % undo normalization
        Y_orig(2,:) = 10.^(Y_orig(2,:)*R_S + R_C); % undo normalization and log
        X_orig = X(:,Ind);
        X_orig(1,:) = X_orig(1,:)*m_S + m_C; % undo normalization
        X_orig(2,:) = X_orig(2,:)*n_S + n_C; % undo normalization and log
        X_orig(3,:) = X_orig(3,:)*gaA_S + gaA_C; % undo normalization
        X_orig(4,:) = X_orig(4,:)*gaC_S + gaC_C; % undo normalization and log
        X_orig(5,:) = 10.^(X_orig(5,:)*gsA_S + gsA_C); % undo normalization
        X_orig(6,:) = 10.^(X_orig(6,:)*gsC_S + gsC_C); % undo normalization and log
        X_orig(7,:) = 10.^(X_orig(7,:)*dhdt_err_S + dhdt_err_C); % undo normalization and log
    
        Y_new = Net_opt(X_new_norm);
        Y_new(1,:) = 10.^(Y_new(1,:)*I_S + I_C); % undo normalization
        Y_new(2,:) = 10.^(Y_new(2,:)*R_S + R_C); % undo normalization and log
        
        x=log10(Y_new(2,:));
        y=Y_new(1,:);
        plot(ax(5),x,y,'+-');
        if mexp==3
            for kk=1:2:size(Y_new,2)
                text(ax(5),x(kk),y(kk)+4,erase(num2str(X_new(7,kk),'%3.2e'),'+0'));
            end
        end
        plot(ax(5),log10(Y_orig(2,:)),Y_orig(1,:),'xr');
        %plot(ax(5),log10(Y_orig(2,:)),Y_orig(1,:),'xr');
        title(ax(5),"dhdt_err"); grid(ax(5),"on"); box(ax(5),"on");
    end
end

function TrainANN

addpath(getenv("froot_tools"));

geometry = 2018;
cycle = 1;
% truncate svd at r basis functions
r=1;

net = feedforwardnet([10 2],'trainlm');

% set early stopping parameters
net.divideFcn= 'dividerand';
net.divideParam.trainRatio = 0.8; % training set [%]
net.divideParam.valRatio = 0.1; % validation set [%]
net.divideParam.testRatio = 0.1; % test set [%]
 
plot_PerturbationResults_Ensemble('Delta_u','m',[1 2]);
load("Delta_u.mat");

% parameter space
% [X(1,:),XC(1),XS(1)] = normalize(gaA);
% [X(2,:),XC(2),XS(2)] = normalize(gaC); 
% [X(3,:),XC(3),XS(3)] = normalize(gsA); 
% [X(4,:),XC(4),XS(4)] = normalize(gsC); 
% [X(5,:),XC(5),XS(5)] = normalize(m); 
% [X(6,:),XC(6),XS(6)] = normalize(n); 

X(:,1)=gaA(:)';
X(:,2)=gaC(:)';
X(:,3)=gsA(:)';
X(:,4)=gsC(:)';
X(:,5)=m(:)';
X(:,6)=n(:)';
% 
% %X = gpuArray(double(X));            
% X=double(X);
% if numel(size(X))==2
%     nx = size(X,1); ny = size(X,2);
%     if nx < ny % number of experiments is highly likely to be larger than number of model parameters
%         X = X'; % transpose data
%     end
% else
%     error("check dimensions of input data");
% end


% training data
data=Delta_u.Calv_dh.map(:,:,cycle);

% this script takes as input 2 matrices: Yv, and X, where 
% - Yv (num_experiments x num_nodes) is the change in ice speed between the
% start and end of each simulation, at each node in the global mesh
% - X (num_experiments x num_parameters) is the set of model parameters
% (randomly sampled) that were used to run each ice sheet model experiment

% rows: nodes, columns: experiments
if numel(size(data))==2
    nx = size(data,1); ny = size(data,2);
    if nx < ny % number of experiments is highly likely to be smaller than number of nodes
        data = data'; % transpose data
    end
else
    error("check dimensions of input data");
end
data(isnan(data))=0;

% Yv=data;
% 
% save ("/mnt/SSD1/Documents/Projects/2024_ANT/UQ_ASE-1.0_Rosier/RNN_surrogate/model_results.mat","Yv","X");


[U,S,V] = SVD_deltau(data,1);


return
U_trunc = U(:,1:r);
S_trunc = S(1:r,1:r);
V_trunc = V(:,1:r);
Vtmp = V_trunc';

[Vtmp2,VC,VS] = normalize(Vtmp(:)); % V = (A-VC)/VS
%V = gpuArray(double(reshape(Vtmp2,size(Vtmp))));
V=double(reshape(Vtmp2,size(Vtmp)));
net = configure(net,X,V);
%view(net);

trainedNet = train(net,X,V,'useParallel','yes','showResources','yes');%,'useGPU','only');
Y = trainedNet(X);
perf = perform(trainedNet,Y,V);

figure; plot(Y,V,'*k'); hold on;
plot([-1.5 3],[-1.5 3]);
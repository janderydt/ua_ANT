function Calc_Emulator_Perturbations

load("Delta_u.mat");


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

return

[U,S,V] = SVD_deltau(data,1);

U_trunc = U(:,1:r);
S_trunc = S(1:r,1:r);
V_trunc = V(:,1:r);
Vtmp = V_trunc';

[Vtmp2,VC,VS] = normalize(Vtmp(:)); % V = (A-VC)/VS
%V = gpuArray(double(reshape(Vtmp2,size(Vtmp))));
V=double(reshape(Vtmp2,size(Vtmp)));
net = configure(net,X,V);
%view(net);
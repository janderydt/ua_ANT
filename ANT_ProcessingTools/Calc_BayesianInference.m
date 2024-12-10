function Calc_BayesianInference

Klear;

%% user defined parameters
cycle = 2; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
pct = 995; % trunction of SVD in emulator (98% seems optimal)
years = [2000 2018]; % a vector with years for which velocity data is used
NN = "RNN"; % which emulator? FNN or RNN

%% initialize UQLAB
addpath(genpath(getenv("froot_matlabfunctions")+"/../UQLab_Rel2.0.0"));

rng(1,'twister'); % set the random number generator for reproducible results
uqlab; % initialize uqlab

%% define prior distribution
% prior distributions for gaA, gaC, gsA, gsC are based on an L-curve 
% analysis (plot_Lcurve.m). Other parameters are given a uniform
% distribution over their support. The bounds might be set automatically,
% but just to be sure we add them here. They can be optained from 
% plot_Inverse_ParameterDistribution.m 
% !!! Make sure that the order of priors is consistent with how the
% customLogLikelihood function expects the parameters as input. Check
% that these are also the parameters in emulator prepare_perturbationresults_for_emulators.m 
% for which the emulators have been trained. By default, we use
% X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n];
PriorOpts.Name = 'Model parameters prior';

ind = 1;
PriorOpts.Marginals(ind).Name = 'log10(gaA)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [1 1]; % mean and std
PriorOpts.Marginals(ind).Bounds = [-1 2];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gaC)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [1 1];
PriorOpts.Marginals(ind).Bounds = [-1 2];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gsA)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [4 0.5];
PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gsC)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [4 0.5];
PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'm';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 9];
PriorOpts.Marginals(ind).Bounds = [2 9];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'n';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 4];
PriorOpts.Marginals(ind).Bounds = [2 4];

if cycle > 1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'log(dhdt_err)';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [0.05 0.5];
    PriorOpts.Marginals(ind).Bounds = [0.05 0.5];
end

myPriorDist = uq_createInput(PriorOpts);

%% define data
myData.Name = 'velocity observations and errors';
myData.y = loadvelocitydata(years); % y is a structure with size n, where n
% is the number of years, and for each year a N-by-2 dimensional matrix, 
% data, with N the number of nodes in the Ua mesh for that year, and 
% 2 columns corresponding to \delta_u and err_\delta_u.

%% define custom likelihood function
% function takes parameter values (m, n, gaA, gaC, log10(gsA), 
% log10(gsC), dhdt_err) and data (y), and returns the log-likelihood
% function at these points. The forward model and covariance matrix are 
% defined within the function.
myLogLikelihood = @(params,y) customLogLikelihood(params, y, [years(2:end),cycle,pct,NN]);

%% define solver options
mySolver.Type = 'MCMC';
mySolver.MCMC.Sampler = 'AIES'; % AM, HMS or AIES (default)
mySolver.MCMC.Steps = 5000; % T=300 default
mySolver.MCMC.NChains = 25; % C=100 default
mySolver.MCMC.Visualize.Parameters = 1:numel(PriorOpts.Marginals);
mySolver.MCMC.Visualize.Interval = 20;

if mySolver.MCMC.Sampler == "AIES"
    mySolver.MCMC.a = 2; % default a=2
else
    error("Define parameters for sampler.")
end

%% define the Bayesian inference model
BayesOpts.Name = 'Bayesian inference for Ua parameters';
BayesOpts.Type = 'Inversion';
BayesOpts.Prior = myPriorDist;
BayesOpts.Data = myData;
BayesOpts.LogLikelihood = myLogLikelihood;
BayesOpts.Solver = mySolver;
BayesOpts.Display = 'verbose'; % quiet, standard or verbose

%% perform analysis
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%% postprocessing and save
uq_postProcessInversionMCMC(myBayesianAnalysis,'pointEstimate','MAP','gelmanRubin','true');

fname = "./BayesianAnalysis/BayesianAnalysis_Steps"+string(mySolver.MCMC.Steps)+"_NChains"+...
    string(mySolver.MCMC.NChains)+"_Calv_dh_"+...
            strjoin(string(years),"_")+"_Weertman_cycle"+string(cycle)+"_"+NN+"_N0k"+string(pct);
save(fname,"myBayesianAnalysis");

%% print some results
uq_print(myBayesianAnalysis);
uq_display(myBayesianAnalysis);

end

function y = loadvelocitydata(years)

% load measured velocity data
%tmp=load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
%Fus_2000=tmp.Fus; Fvs_2000=tmp.Fvs; Fxerr_2000=tmp.Fxerr; Fyerr_2000=tmp.Fyerr;
addpath(getenv("froot_data")+"Measures/Measures_annual");

% assemble velocity fields
for ii=1:numel(years)
    
    years_measures = string(years(ii))+"_"+string(years(ii)+1);
    [vx,x_meas,y_meas]=measures_annual("vx",years_measures); 
    vy=measures_annual("vy",years_measures);
    stdx=measures_annual("vxstd",years_measures);
    stdy=measures_annual("vystd",years_measures);
    y_meas = flip(y_meas,1); [Xm,Ym]=ndgrid(x_meas,y_meas);
    v(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(vx,vy),2));
    std(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(stdx,stdy),2));

end

% assemble velocity differences
for ii=1:numel(years)-1
    yr = string(years(ii+1));
    % Ua mesh
    fname = "Delta_u_"+yr+"_AS_Weertman.mat";
    if yr=="2009"
        load(fname,"MUA_2009","GF_2009");
        MUA = MUA_2009; GF = GF_2009;
    elseif yr=="2014"
        load(fname,"MUA_2014","GF_2014");
        MUA = MUA_2014; GF = GF_2014;
    elseif yr=="2018"
        load(fname,"MUA_2018","GF_2018");
        MUA = MUA_2018; GF = GF_2018;
    else
        error("unknown year")
    end
    % interpolate initial and final speed onto Ua mesh
    u_init = v(ii).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
    std_init = std(ii).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
    u = v(ii+1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
    std = std(ii+1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));

    deltau = u-u_init; % measured change in speed (m/yr)
    deltau_err = hypot(std_init,std); % error estimate (m/yr)
    
    % remove nans
    Ind = find(isnan(deltau) | isnan(deltau_err));
    deltau_err(Ind) = 1e3;
    deltau(Ind) = 0;

    % what if we remove the floating ice?
    deltau(GF.node<0.5) = 0;
    deltau_err(GF.node<0.5) = 1e3;

    y(ii).data = [deltau(:) deltau_err(:)];
end

end

function logL = customLogLikelihood(params, measurements, mymodel)

persistent M

% evaluates the log likelihood for parameter values specified in params
% INPUTS:
% ** params is C-by-M matrix, where C is the number of MCMC chains and M is
% the number of parameters
% ** measurements is a data structure, as specified by the loadvelocitydata 
% function
% ** years is a vector with years for which data/fwd model is required
% OUTPUT:
% ** logL is a C-by-1 vector

% Initialization
nReal = size(params,1); % number of queried realizations
nNodes = numel(measurements(1).data(:,1));
years = mymodel(1:end-3);
cycle = mymodel(end-2);
pct = mymodel(end-1);
NN = mymodel(end);

% Load forward model(s)
if isempty(M)
    for ii=1:numel(years)
        % Check prepare_perturbationresults_for_emulators.m to see
        % what parameters the emulator is trained for, and in what order
        % The default is X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n dhdt_err];
        if NN=="RNN"
            net_tmp=importNetworkFromTensorFlow("./RNN/TF_files/tuned_model_forAll_Calv_dh_"+...
            string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct));
            M(ii).net=initialize(net_tmp);
            load("./RNN/mat_files/data_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "X_train_C","X_train_S","T_train_C","T_train_S");
            load("./RNN/mat_files/SVD_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "T_reproj","T_mean");
            load("./RNN/mat_files/MSE_RNN_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat");
            MSE = mse_RNN;
        elseif NN=="FNN"
            load("./FNN/mat_files/FNN_trainscg_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "Net_opt");
            load("./FNN/mat_files/SVD_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "T_reproj","T_mean");
            M(ii).net=Net_opt.trained;
            X_train_C = Net_opt.X_train_C(:)';
            X_train_S = Net_opt.X_train_S(:)';
            T_train_C = Net_opt.T_train_C(:)';
            T_train_S = Net_opt.T_train_S(:)';
            load("./FNN/mat_files/MSE_FNN_trainscg_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat");
            MSE = mse_FNN;
        else
            error("Unknown emulator type");
        end
        Nmodes = size(T_reproj,1);
        M(ii).X_train_C = X_train_C;
        M(ii).X_train_S = X_train_S;
        M(ii).T_train_C = T_train_C;
        M(ii).T_train_S = T_train_S;
        % Subtract mean and project measurements onto truncated basis.
        % These are the unnormalized values.
        M(ii).data = (measurements(1).data(:,1)-T_mean(:))'*T_reproj';
        % Assemble covariance matrices
        % 1. measurement errors
        S_meas = T_reproj*spdiags(measurements(1).data(:,2).^2,0,nNodes,nNodes)*T_reproj';
        % 2. Ua errors: obtained from  FIXME
        S_ua = 0*S_meas;
        % 3. Emulator errors: obtained from plot_Emulator_MSE
        S_emulator = spdiags(MSE(:),0,Nmodes,Nmodes);
        M(ii).S = S_meas + S_ua + S_emulator;
        M(ii).detS = det(M(ii).S);
        M(ii).Sinv = inv(M(ii).S);
    end    
end

% Loop through years and realizations
logL = zeros(nReal,1);

for ii=1:numel(years)
    % Apply normalization to parameters before feeding into emulator
    predictors = (params-repmat(M(ii).X_train_C,nReal,1))./repmat(M(ii).X_train_S,nReal,1);

    % Evaluate forward model. Undo normalization of the output but keep in the
    % projected basis to make it compatible with data
    if NN=="RNN"
        modelRun = double(predict(M(ii).net,predictors)); 
    elseif NN=="FNN"
        modelRun = double(M(ii).net(predictors')');
    end
    modelRun = modelRun.*repmat(M(ii).T_train_S,nReal,1)+repmat(M(ii).T_train_C,nReal,1);

    for jj = 1:nReal
      % Evaluate log-likelihood
      logLikeli = - 1/2*log(2*pi*M(ii).detS) - 1/2*(M(ii).data...
        -modelRun(jj,:))*M(ii).Sinv*(M(ii).data-modelRun(jj,:))';
      % Assign to logL vector
      logL(jj) = logL(jj)+logLikeli;
    end
end

end




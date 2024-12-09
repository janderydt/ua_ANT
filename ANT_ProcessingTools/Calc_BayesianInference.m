function Calc_BayesianInference

Klear;

%% user defined parameters
cycle = 1; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
pct = 98; % trunction of SVD in emulator (98% seems optimal)
years = 2018; % a vector with years for which velocity data is used

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
    PriorOpts.Marginals(ind).Name = 'dhdt_err';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = 'Uniform';
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
myLogLikelihood = @(params,y) customLogLikelihood(params, y, [years,cycle,pct]);

%% define solver options
mySolver.Type = 'MCMC';
mySolver.MCMC.Sampler = 'AIES'; % AM, HMS or AIES (default)
mySolver.MCMC.Steps = 1000; % T=300 default
mySolver.MCMC.NChains = 50; % C=100 default
mySolver.MCMC.Visualize.Parameters = 1:6;
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


%% perforn analysis
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%% print some results
uq_print(myBayesianAnalysis);
uq_display(myBayesianAnalysis);

end


function y = loadvelocitydata(years)

% load measured velocity data
tmp=load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
Fus_2000=tmp.Fus; Fvs_2000=tmp.Fvs; Fxerr_2000=tmp.Fxerr; Fyerr_2000=tmp.Fyerr;
       
for ii=1:numel(years)
    
    yr = string(years(ii));

    % Ua mesh
    fname = "Delta_u_"+yr+"_AS_Weertman.mat";
    if yr=="2009"
        load(fname,"MUA_2009");
        MUA = MUA_2009;
        load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_2009-2010_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
    elseif yr=="2014"
        load(fname,"MUA_2014");
        MUA = MUA_2014;
        load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_2014-2015_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
    elseif yr=="2018"
        load(fname,"MUA_2018");
        MUA = MUA_2018;
        load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_2018-2019_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
    else
        error("unknown year")
    end
    
    % interpolate initial and final speed onto Ua mesh
    ux2000 = Fus_2000(MUA.coordinates(:,1),MUA.coordinates(:,2));
    uxerr2000 = Fxerr_2000(MUA.coordinates(:,1),MUA.coordinates(:,2));
    uy2000 = Fvs_2000(MUA.coordinates(:,1),MUA.coordinates(:,2));
    uyerr2000 = Fyerr_2000(MUA.coordinates(:,1),MUA.coordinates(:,2));
    u2000 = hypot(ux2000,uy2000);
    uerr2000 = hypot(uxerr2000,uyerr2000); % as defined in https://nsidc.org/sites/default/files/nsidc-0720-v001-userguide_0.pdf
        
    ux = Fus(MUA.coordinates(:,1),MUA.coordinates(:,2));
    uxerr = Fxerr(MUA.coordinates(:,1),MUA.coordinates(:,2));
    uy = Fvs(MUA.coordinates(:,1),MUA.coordinates(:,2));
    uyerr = Fyerr(MUA.coordinates(:,1),MUA.coordinates(:,2));
    u = hypot(ux,uy);
    uerr = hypot(uxerr,uyerr);

    deltau = u-u2000; % measured change in speed (m/yr)
    deltau_err = hypot(uerr2000,uerr); % error estimate (m/yr)
    
    deltau_err(isnan(deltau)) = 1e5;
    deltau(isnan(deltau)) = 0;

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
years = mymodel(1:end-2);
cycle = mymodel(end-1);
pct = mymodel(end);

% Load forward model(s)
if isempty(M)
    for ii=1:numel(years)
        net_tmp=importNetworkFromTensorFlow("./RNN/TF_files/tuned_model_forAll_Calv_dh_"+...
            string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N"+string(pct));
        % Check prepare_perturbationresults_for_emulators.m to see
        % what parameters the emulator is trained for, and in what order
        % The default is X = [log10(gaA(:)) log10(gaC(:)) log10(gsA(:)) log10(gsC(:)) m(:) n(:)];
        M(ii).net=initialize(net_tmp);
        load("./RNN/mat_files/data_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
            "X_train_C","X_train_S","T_train_C","T_train_S");
        load("./RNN/mat_files/SVD_Calv_dh_"+string(years(ii))+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
            "T_reproj","T_mean");
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
        % 3. Emulator errors: obtained from plot_Emulator_MSE FIXME
        S_emulator = 0*S_meas;
        M(ii).S = S_meas + S_ua + S_emulator;
        M(ii).detS = det(M(ii).S);
        M(ii).Sinv = inv(M(ii).S);
    end    
end

% Apply normalization to parameters before feeding into emulator
predictors = (params-repmat(M(1).X_train_C,nReal,1))./repmat(M(1).X_train_S,nReal,1);

% Evaluate forward model. Undo normalization of the output but keep in the
% projected basis to make it compatible with data
modelRun = double(predict(M(1).net,predictors)); 
modelRun = modelRun.*repmat(M(1).T_train_S,nReal,1)+repmat(M(1).T_train_C,nReal,1);

% Loop through realizations
logL = zeros(nReal,1);
for ii = 1:nReal
  % Evaluate log-likelihood
  logLikeli = - 1/2*log(2*pi*M(1).detS) - 1/2*(M(1).data...
    -modelRun(ii,:))*M(1).Sinv*(M(1).data-modelRun(ii,:))';
  % Assign to logL vector
  logL(ii) = logLikeli;
end

end




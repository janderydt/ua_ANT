function Calc_BayesianInference(UserVar)

if nargin==0
    Klear
    %% user defined parameters
    UserVar.a = 2;
    UserVar.l = 50e3; % range in the ua semivariogram
    UserVar.fac_su = 3; % multiplication factor for model errors compared to observation errors
    UserVar.domain = "AMUND";
    UserVar.slidinglaw = ["Weertman" "Weertman" "Weertman" "Umbi" "Umbi" "Umbi"];
    UserVar.cycle = [1 1 1 1 1 1]; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
    UserVar.dataformat = ["du" "du" "du" "du" "du" "du"]; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
    UserVar.only_grounded_ice = [1 1 1 1 1 1];
    UserVar.years = ["2000-2009" "2000-2014" "2000-2020" "2000-2009" "2000-2014" "2000-2020"]; % a vector with years for which velocity data is used
    UserVar.NN = ["RNN" "RNN" "RNN" "RNN" "RNN" "RNN"]; % which emulator? FNN or RNN
    UserVar.do_plots = 1;
end

% dataformat = ["du"]; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
% cycle = [1]; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
% pct = [9897]; % trunction of SVD in emulator
% only_grounded_ice = [1];
% years = ["2000-2009"]; % a vector with years for which velocity data is used
% NN = ["RNN"]; % which emulator? FNN or RNN

%% initialize UQLAB
addpath(genpath(getenv("froot_matlabfunctions")+"/../UQLab_Rel2.0.0"));

rng(1,'twister'); % set the random number generator for reproducible results
clear uqlab;
uqlab; % initialize uqlab

%% define prior distribution
% prior distributions for gaA, gaC, gsA, gsC are based on an L-curve 
% analysis (plot_Lcurve.m). Other parameters are given a uniform
% distribution over their support. The bounds might be set automatically,
% but just to be sure we add them here. They can be optained from 
% plot_Inverse_ParameterDistribution.m 
% !!! Make sure that the order of priors is consistent with how the
% customLogLikelihood function expects the parameters as input. Check
% that these are also the parameters in emulator prepare_data_for_***_emulators.m 
% for which the emulators have been trained. By default, we use
% X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n log10(dhdt_err)];
PriorOpts.Name = 'Model parameters prior';

ind = 1;
PriorOpts.Marginals(ind).Name = 'log10(gaA)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [0.7 0.51]; % mean and std
% mean is based on L curve
% std is chosen such that integral under distribution between mean-1 
% and mean+1 (in log space) captures 95% of the variability:
% -0.5*(1+erf((-0.7+(0.7-1))/(sqrt(2)*0.51)))+0.5*(1+erf(((0.7+1)-0.7)/(sqrt(2)*0.51)))=0.95
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [-1 log10(200)];
PriorOpts.Marginals(ind).Bounds = [-1 log10(200)];
ind = ind+1;

PriorOpts.Marginals(ind).Name = 'log10(gaC)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [0.7 0.51];
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [-1 log10(50)];
PriorOpts.Marginals(ind).Bounds = [-1 log10(50)];
ind = ind+1;

PriorOpts.Marginals(ind).Name = 'log10(gsA)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [4.3 0.51];
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [3 6];
PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;

PriorOpts.Marginals(ind).Name = 'log10(gsC)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [4.3 0.51];
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [3 6];
PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;

PriorOpts.Marginals(ind).Name = 'm';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 9];
PriorOpts.Marginals(ind).Bounds = [2 9];
ind = ind+1;

PriorOpts.Marginals(ind).Name = 'n';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 5];
PriorOpts.Marginals(ind).Bounds = [2 5];

% add dhdt if any cycle>1
if max(UserVar.cycle) > 1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'log(dhdt_err)';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [log10(0.05) log10(0.5)];
    PriorOpts.Marginals(ind).Bounds = [log10(0.05) log10(0.5)];
end

% add discrete variables
nslidinglaw = numel(unique(UserVar.slidinglaw));
ncycle = numel(unique(UserVar.cycle));
if ncycle>1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'Cycle';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [0 1];
    PriorOpts.Marginals(ind).Bounds = [0 1];
end
if nslidinglaw>1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'SlidingLaw';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [0 1];
    PriorOpts.Marginals(ind).Bounds = [0 1];
end

myPriorDist = uq_createInput(PriorOpts);

%% define data
myData.y = loadvelocitydata(UserVar.dataformat,UserVar.years,UserVar.only_grounded_ice); % myData.y.(year).u and myData.y.(year).stdu 
    % are N-by-1 matrices with N the number of nodes in the Ua mesh for that year
myData.Name = 'velocity observations and errors';

%% define custom likelihood function
% function takes parameter values (log10(gaA), log19(gaC), log10(gsA), 
% log10(gsC), m, n, log10(dhdt_err)) and data (y), and returns the log-likelihood
% function at these points. The forward model and covariance matrix are 
% defined within the function.

myLogLikelihood = @(params,y) customLogLikelihood(params, y, UserVar);

%% define solver options
mySolver.Type = 'MCMC';
mySolver.MCMC.Sampler = 'AIES'; % AM, HMS or AIES (default)
mySolver.MCMC.Steps = 20000; % T=300 default
mySolver.MCMC.NChains = 25; % C=100 default
mySolver.MCMC.Visualize.Parameters = 1:numel(PriorOpts.Marginals);
mySolver.MCMC.Visualize.Interval = 20;

if mySolver.MCMC.Sampler == "AIES"
    mySolver.MCMC.a = UserVar.a; % default a=2, Rosier et al.: a=1.5
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

UserVar.location = pwd+"/BayesianAnalysis/";
fname = "BayesianAnalysis_Steps"+string(mySolver.MCMC.Steps)+"_NChains"+...
    string(mySolver.MCMC.NChains)+"_"+UserVar.dataformat(1)+"_Calv_dh_"+...
            strjoin(string(UserVar.years),"_")+"_"+strjoin(string(UserVar.slidinglaw),"_")+"_cycle"+string(UserVar.cycle(1))+...
            "_floatingice"+string(1-UserVar.only_grounded_ice(1));
n_exist = numel(dir(UserVar.location+fname(1)+"*.mat"));
UserVar.fname = fname(1)+"_v"+string(n_exist+1)+".mat";

save(UserVar.location+UserVar.fname,"myBayesianAnalysis","UserVar");

Calc_MarginalsInference(UserVar.location,UserVar.fname);

%% print some results
uq_print(myBayesianAnalysis);
%uq_display(myBayesianAnalysis);

if UserVar.do_plots
    plot_myBayesianAnalysis(myBayesianAnalysis,UserVar);
end

end


function logL = customLogLikelihood(params, measurements, UserVar)

persistent M N

% evaluates the log likelihood for parameter values specified in params

% INPUTS:
% ** params is C-by-M matrix, where C is the number of MCMC chains and M is
% the number of parameters
% ** measurements is a data structure, as specified by the loadvelocitydata 
% function
% ** mymodel is a vector with years for which data/fwd model is required,
% as well as the cycle, pct and NN to be used

% OUTPUT:
% ** logL is a C-by-1 vector

% Initialization
nReal = size(params,1); % number of queried realizations
% mymodel has format [domain,slidinglaw,cycle,dataformat,only_grounded_ice,years,NN],
% each of which is an 1xn vector apart from domain. 
domain = UserVar.domain;
slidinglaw = UserVar.slidinglaw;
cycle = UserVar.cycle;
dataformat = UserVar.dataformat;
only_grounded_ice = double(UserVar.only_grounded_ice);
years = UserVar.years;
NN = UserVar.NN;
l = UserVar.l; % range in the semivariogram
fac_su = UserVar.fac_su; % multiplication factor for model errors compared to observation errors
N = numel(slidinglaw);

% Load forward model(s)
if isempty(M)
    for ii=1:N
        if dataformat(ii) == "LOGu"
            yearstr = "LOGu"+years(ii); 
        else
            yearstr = years(ii);
        end

        % Check prepare_perturbationresults_for_emulators.m to see
        % what parameters the emulator is trained for, and in what order
        % The default is X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n log10(dhdt_err)];
        if NN(ii)=="RNN"

            TF_dir = dir("./RNN/TF_files/tuned_model_"+domain+"_Calv_dh_"+years(ii)+...
                "_"+slidinglaw(ii)+"_cycle"+string(cycle(ii))+"_floatingice"+string(1-only_grounded_ice(ii))+...
                "_includemeasurements0*");
            net_tmp=importNetworkFromTensorFlow(TF_dir.folder+"/"+TF_dir.name);
            M(ii).net=initialize(net_tmp);

            data_file = dir("./RNN/mat_files/data_"+domain+"_Calv_dh_"+years(ii)+...
                "_"+slidinglaw(ii)+"_cycle"+string(cycle(ii))+"_floatingice"+string(1-only_grounded_ice(ii))+...
                "_includemeasurements0*.mat");
            load(data_file.folder+"/"+data_file.name,"X_train_C","X_train_S","T_train_C","T_train_S");
        
            SVD_file = dir("./RNN/mat_files/SVD_"+domain+"_Calv_dh_"+years(ii)+...
                "_"+slidinglaw(ii)+"_cycle"+string(cycle(ii))+"_floatingice"+string(1-only_grounded_ice(ii))+...
                "_includemeasurements0*.mat");
            load(SVD_file.folder+"/"+SVD_file.name,"T_reproj","T_mean");

            MSE_file = dir("./RNN/mat_files/MSE_RNN_"+domain+"_Calv_dh_"+years(ii)+...
                "_"+slidinglaw(ii)+"_cycle"+string(cycle(ii))+"_floatingice"+string(1-only_grounded_ice(ii))+...
                "_includemeasurements0*.mat");
            load(MSE_file.folder+"/"+MSE_file.name);
            MSE = mse_RNN;

        elseif NN(ii)=="FNN"
            load("./FNN/mat_files/FNN_trainscg_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat",...
                "Net_opt");
            load("./FNN/mat_files/SVD_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat",...
                "T_reproj","T_mean");
            M(ii).net=Net_opt.trained;
            X_train_C = Net_opt.X_train_C(:)';
            X_train_S = Net_opt.X_train_S(:)';
            T_train_C = Net_opt.T_train_C(:)';
            T_train_S = Net_opt.T_train_S(:)';
            load("./FNN/mat_files/MSE_FNN_trainscg_Calv_dh_"+yearstr+"_Weertman_cycle"+string(cycle(ii))+...
                "_floatingice"+string(1-only_grounded_ice(ii))+"_N0k"+string(pct(ii))+".mat");
            MSE = mse_FNN;
        else
            error("Unknown emulator type");
        end
        [nModes,nNodes] = size(T_reproj);
        M(ii).X_train_C = X_train_C;
        M(ii).X_train_S = X_train_S;
        M(ii).T_train_C = T_train_C;
        M(ii).T_train_S = T_train_S;
        
        if dataformat(ii) == "du"
            stdstr = "stddu";
            ustr = "du";
            yearstr2 = "yr"+replace(years(ii),"-","_yr");
        else
            stdstr = "stdu";
            ustr = "u";
            yearstr2 = "yr"+years(ii);
        end

        % Subtract mean and deal with nans in data
        data_tmp = measurements(1).(yearstr2).(ustr)(:)-T_mean(:);
        std_tmp = measurements(1).(yearstr2).(stdstr)(:);
        Indnan = find(isnan(data_tmp) | isnan(std_tmp));
        if dataformat(ii) == "LOGu"
            data_tmp(Indnan) = 0; 
            std_tmp(Indnan) = 3;
        elseif dataformat(ii) == "u"
            data_tmp(Indnan) = 0; 
            std_tmp(Indnan) = 1e3;
        elseif dataformat(ii) == "du"
            data_tmp(Indnan) = 0; 
            std_tmp(Indnan) = 1e3;
        end

        % Project data onto truncated space. Make sure that this does not
        % distort the data too much, i.e. the data can be represented well in the
        % truncated base. To check this, you need to reproject M(ii).data back onto
        % the nodal base by multiplying with B_trunc, which is saved in SVD_Calv_dh_****
        % If the data cannot be faithfully represented using the truncated
        % basis functions then you might want to try adding the data to the SVD
        % decomposition, or increase the number of basis functions (both
        % are options in prepare_data_for_***_emulators.m)
        M(ii).data = data_tmp'*T_reproj';

        % Assemble covariance matrices
        % 1. measurement errors
        S_meas = T_reproj*diag(std_tmp.^2)*T_reproj';

        % 2. Ua errors: obtained from  FIXME
        su = fac_su*mean(std_tmp); %m/yr
        % if dataformat == "LOGu"
        %     s_u = log10(s_u);
        % end
        load("Delta_u_"+domain+"_"+slidinglaw(1)+"_"+yearstr+".mat","MUA_yr2","GF_yr2");
        MUA = MUA_yr2;
        GF= GF_yr2;
        D_tmp = pdist2(MUA.coordinates,MUA.coordinates,"squaredeuclidean");
        su_nodal = su*ones(MUA.Nnodes,1);
        if only_grounded_ice(ii)==0
            su_nodal(GF.node<0.5)=2*su; % give floating nodes a larger error
        end
        D_tmp = tensorprod(su_nodal,su_nodal,2).*exp(-D_tmp/(2*l^2));

        S_ua = T_reproj*D_tmp*T_reproj';
        clear D_tmp

        % 3. Emulator errors: obtained from plot_Emulator_MSE
        S_emulator = diag(MSE);

        M(ii).S = S_meas + S_emulator + S_ua;
        M(ii).detS = det(M(ii).S);
        M(ii).Sinv = inv(M(ii).S);
    end    
end

% Loop through years and realizations
logL = zeros(nReal,1);

% Deal with discrete variables and which emulator to use
ncycles = numel(unique(cycle));
nslidinglaws = numel(unique(slidinglaw));
% make a copy of the MCMC parameters
params_tmp = params;
% set up some dummy variables
ind=1;
% which emulator to use for the cycles?
if ncycles > 1    % cycles is a discrete variable that takes multiple values
    if max(cycle)>1
        param_discr(:,ind)=round(params(:,8))+1; 
        params_tmp(:,8)=[];
    else
        param_discr(:,ind)=round(params(:,7))+1;
        params_tmp(:,7)=[];
    end
    ind=ind+1;
end
% which emulator to use for the sliding law?
if nslidinglaws > 1  % sliding law is a discrete variable that takes multiple values
    param_discr(:,ind)=round(params(:,end));
    param_discr = string(param_discr);
    params_tmp(:,end)=[];
    param_discr=strrep(param_discr,"0","Weertman");
    param_discr=strrep(param_discr,"1","Umbi");
end
params=params_tmp; % with discrete parameters removed

% deal with multiple years
EmulatorIndex = zeros(nReal,N);
for ii=1:N
    if ncycles > 1
        Itmp = find(param_discr == UserVar.cycle(ii));
    elseif nslidinglaws > 1
        Itmp = find(param_discr == UserVar.slidinglaw(ii));
    else
        Itmp = 1:nReal;
    end
    EmulatorIndex(Itmp,ii) = ii;
end

% unique combinations of emulators
Emulators_Unique = unique(EmulatorIndex,'rows');

% now assemble likelihood
for ii=1:size(Emulators_Unique,1)

    Itmp = find(ismember(EmulatorIndex,Emulators_Unique(ii,:),'rows'));
    nItmp = numel(Itmp);
    Ind_Emulators = find(Emulators_Unique(ii,:)~=0);

    for ee=Ind_Emulators

        % Apply normalization to parameters before feeding into emulator
        if cycle(ii)==1
            predictors = (params(Itmp,1:6)-repmat(M(ee).X_train_C,nItmp,1))./repmat(M(ee).X_train_S,nItmp,1);
        else
            predictors = (params(Itmp,:)-repmat(M(ee).X_train_C,nItmp,1))./repmat(M(ee).X_train_S,nItmp,1);
        end

        % Evaluate forward model. 
        if NN(ee)=="RNN"
            modelRun = double(predict(M(ee).net,predictors)); 
        elseif NN(ee)=="FNN"
            modelRun = double(M(ee).net(predictors')');
        end

        % Undo normalization of the output but keep in the projected basis to 
        % make it compatible with data
        modelRun = modelRun.*repmat(M(ee).T_train_S,nItmp,1)+repmat(M(ee).T_train_C,nItmp,1);
        Nout = size(modelRun,2);

        % Assemble log likelihood
        for jj = 1:nItmp
          % Evaluate log-likelihood
          prefac = -0.5*log((2*pi)^Nout*M(ee).detS);
          Ndistrib = - 0.5*(M(ee).data...
            -modelRun(jj,:))*M(ee).Sinv*(M(ee).data-modelRun(jj,:))';
          logLikeli = prefac + Ndistrib;
          % Assign to logL vector
          logL(Itmp(jj)) = logL(Itmp(jj))+logLikeli;
        end

        %figure(115); hold on; plot(params(1,5),prefac(1),'.',Color=CM(ii,:)); title('m');
        %plot(params(1,5),Ndistrib(1),'.',Color=CM(ii+2,:)); title('m');

    end

end

%figure(111); hold on; plot(params(1,1),logL(1),'.k'); title('gaA');
%figure(112); hold on; plot(params(1,2),logL(1),'.k'); title('gaC');
%figure(113); hold on; plot(params(1,3),logL(1),'.k'); title('gsA');
%figure(114); hold on; plot(params(1,4),logL(1),'.k'); title('gsC');
%figure(116); hold on; plot(params(1,6),logL(1),'.k'); title('n');

end




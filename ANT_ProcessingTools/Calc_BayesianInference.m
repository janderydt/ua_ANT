function Calc_BayesianInference

Klear;

%% user defined parameters
dataformat = "LOGu"; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
cycle = 2; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
pct = 960; % trunction of SVD in emulator
years = [2000]; % a vector with years for which velocity data is used
NN = "FNN"; % which emulator? FNN or RNN

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
% X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n log10(dhdt_err)];
PriorOpts.Name = 'Model parameters prior';

ind = 1;
PriorOpts.Marginals(ind).Name = 'log10(gaA)';
% PriorOpts.Marginals(ind).Type = 'Gaussian';
% PriorOpts.Marginals(ind).Parameters = [1 0.5]; % mean and std
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [-1 log10(200)];
% PriorOpts.Marginals(ind).Bounds = [-1 log10(200)];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gaC)';
% PriorOpts.Marginals(ind).Type = 'Gaussian';
% PriorOpts.Marginals(ind).Parameters = [1 0.5];
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [-1 log10(50)];
%PriorOpts.Marginals(ind).Bounds = [-1 log10(50)];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gsA)';
% PriorOpts.Marginals(ind).Type = 'Gaussian';
% PriorOpts.Marginals(ind).Parameters = [4 0.5];
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [3 6];
%PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gsC)';
%PriorOpts.Marginals(ind).Type = 'Gaussian';
%PriorOpts.Marginals(ind).Parameters = [4 0.5];
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [3 6];
%PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'm';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 9];
%PriorOpts.Marginals(ind).Bounds = [2 9];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'n';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 4];
%PriorOpts.Marginals(ind).Bounds = [2 4];

if cycle > 1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'log(dhdt_err)';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [log10(0.05) log10(0.5)];
    %PriorOpts.Marginals(ind).Bounds = [log10(0.05) log10(0.5)];
end

myPriorDist = uq_createInput(PriorOpts);

%% define data
myData.y = loadvelocitydata(dataformat,years); % myData.y.(year).u and myData.y.(year).stdu 
% are N-by-1 matrices with N the number of nodes in the Ua mesh for that year
myData.Name = 'velocity observations and errors';

%% define custom likelihood function
% function takes parameter values (log10(gaA), log19(gaC), log10(gsA), 
% log10(gsC), m, n, log10(dhdt_err)) and data (y), and returns the log-likelihood
% function at these points. The forward model and covariance matrix are 
% defined within the function.
myLogLikelihood = @(params,y) customLogLikelihood(params, y, [dataformat,years,cycle,pct,NN]);

%% define solver options
mySolver.Type = 'MCMC';
mySolver.MCMC.Sampler = 'AIES'; % AM, HMS or AIES (default)
mySolver.MCMC.Steps = 20000; % T=300 default
mySolver.MCMC.NChains = 20; % C=100 default
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
    string(mySolver.MCMC.NChains)+"_"+dataformat+"_Calv_dh_"+...
            strjoin(string(years),"_")+"_Weertman_cycle"+string(cycle)+"_"+NN+"_N0k"+string(pct);
save(fname,"myBayesianAnalysis");

%% print some results
uq_print(myBayesianAnalysis);
uq_display(myBayesianAnalysis);

% MAP=myBayesianAnalysis.Results.PostProc.PointEstimate.X{:};
% % nearest predictor and target
% load("u_AS_Calv_dh_cycle2_Weertman_2000_2009_2014_2018.mat");
% Ind_toremove = find(gaC>50);
% X = [log10(gaA(:)) log10(gaC(:)) log10(gsA(:)) log10(gsC(:)) m(:) n(:)];
% if cycle==2 % add dhdt_err to predictors
%     X = [X log10(dhdt_err(:))];
% end
% X = double(X);
% X(Ind_toremove,:)=[];
% Ind_MAP=knnsearch(X,MAP);
% MAP
% X(Ind_MAP,:)
% deltau = speed.yr2000(Ind_MAP,:);
% deltau(deltau<0)=nan;
% out = loadvelocitydata("u",[2000]);
% figure; PlotMeshScalarVariable([],MUA.yr2000,log10(deltau(:))); 
% deltau_meas = out.yr2000.u;
% deltau_meas(deltau_meas<0)=nan;
% figure; PlotMeshScalarVariable([],MUA.yr2000,log10(deltau_meas(:))); 
% figure; PlotMeshScalarVariable([],MUA.yr2000,deltau_meas(:)-deltau(:));



end

function out = loadvelocitydata(dataformat,years)

% load measured velocity data
%tmp=load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
%Fus_2000=tmp.Fus; Fvs_2000=tmp.Fvs; Fxerr_2000=tmp.Fxerr; Fyerr_2000=tmp.Fyerr;
addpath(getenv("froot_data")+"Measures/Measures_annual");

% assemble velocity fields
for ii=1:numel(years)

    % years_measures = string(years(ii))+"_"+string(years(ii)+1);
    % [vx,x_meas,y_meas]=measures_annual("vx",years_measures); 
    % vy=measures_annual("vy",years_measures);
    % stdx=measures_annual("vxstd",years_measures);
    % stdy=measures_annual("vystd",years_measures);
    % y_meas = flip(y_meas,1); [Xm,Ym]=ndgrid(x_meas,y_meas);
    % v(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(vx,vy),2));
    % std(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(stdx,stdy),2));
    if years(ii)==2000
        fname = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    else
        fname = "GriddedInterpolants_"+string(years(ii))+"-"+string(years(ii)+1)+"_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    end
    load("../ANT_Data/ANT_Interpolants/"+fname);
    v_tmp = hypot(Fus.Values,Fvs.Values);
    Fu = Fus; Fu.Values = v_tmp;
    std_tmp = hypot(Fxerr.Values,Fyerr.Values);
    Fstd = Fus; Fstd.Values = std_tmp;
    v(ii).F = Fu;
    std(ii).F = Fstd;

end

if dataformat == "du"
    % assemble velocity differences
    for ii=1:numel(years)-1
        % Ua mesh
        yr1 = string(years(ii));
        yr2 = string(years(ii+1));
        load("Delta_u_AS_Weertman_"+yr1+"-"+yr2+".mat","MUA_yr2","GF_yr2");
        MUA = MUA_yr2; GF = GF_yr2;
    
        % interpolate initial and final speed onto Ua mesh
        u_init = v(ii).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_init = std(ii).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        u_target = v(ii+1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_target = std(ii+1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
    
        deltau = u_target-u_init; % measured change in speed (m/yr)
        deltau_err = hypot(std_init,std_target); % error estimate (m/yr)
        
        % remove nans
        Ind = find(isnan(deltau) | isnan(deltau_err));
        deltau_err(Ind) = 1e3;
        deltau(Ind) = 0;
    
        % what if we remove the floating ice?
        %deltau(GF.node<0.5) = 0;
        %deltau_err(GF.node<0.5) = 1e3;
    
        out.("yr"+years(ii)+"_yr"+years(ii+1)).du = deltau(:);
        out.("yr"+years(ii)+"_yr"+years(ii+1)).stddu = deltau_err(:);
    end
elseif ismember(dataformat,["u","LOGu"])
    for ii=1:numel(years)
        load("u_AS_Calv_dh_cycle2_Weertman_2000_2009_2014_2018.mat","MUA","GF");
        MUA = MUA.("yr"+years(ii)); GF = GF.("yr"+years(ii));

        % interpolate initial and final speed onto Ua mesh
        u_tmp = v(ii).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_tmp = std(ii).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        
        % remove nans
        Ind = find(isnan(u_tmp) | isnan(std_tmp));
        std_tmp(Ind) = 5e3;
        u_tmp(Ind) = eps(0);

        % apply log if needed
        if dataformat=="LOGu"
            u_tmp = log10(u_tmp);
            std_tmp = std_tmp./(u_tmp*log(10));
        end
        out.("yr"+years(ii)).u = u_tmp(:);
        out.("yr"+years(ii)).stdu = std_tmp(:);

    end
end

end

function logL = customLogLikelihood(params, measurements, mymodel)

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
dataformat = mymodel(1);
years = mymodel(2:end-3);
cycle = mymodel(end-2);
pct = mymodel(end-1);
NN = mymodel(end);

% Load forward model(s)
if isempty(M)
    if dataformat == "du"
        N = numel(years)-1;
        for ii=1:N
            yearstr(ii) = years(ii)+"-"+years(ii+1);
        end
    elseif dataformat == "u"
        N = numel(years);
        yearstr = "u"+years;
    elseif dataformat == "LOGu"
        N = numel(years);
        yearstr = "LOGu"+years;
    end
    for ii=1:N
        % Check prepare_perturbationresults_for_emulators.m to see
        % what parameters the emulator is trained for, and in what order
        % The default is X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n log10(dhdt_err)];
        if NN=="RNN"
            net_tmp=importNetworkFromTensorFlow("./RNN/TF_files/tuned_model_Calv_dh_"+...
            yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct));
            M(ii).net=initialize(net_tmp);
            load("./RNN/mat_files/data_Calv_dh_"+yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "X_train_C","X_train_S","T_train_C","T_train_S");
            load("./RNN/mat_files/SVD_Calv_dh_"+yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "T_reproj","T_mean");
            load("./RNN/mat_files/MSE_RNN_Calv_dh_"+yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat");
            MSE = mse_RNN;
        elseif NN=="FNN"
            load("./FNN/mat_files/FNN_trainscg_Calv_dh_"+yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "Net_opt");
            load("./FNN/mat_files/SVD_Calv_dh_"+yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat",...
                "T_reproj","T_mean");
            M(ii).net=Net_opt.trained;
            X_train_C = Net_opt.X_train_C(:)';
            X_train_S = Net_opt.X_train_S(:)';
            T_train_C = Net_opt.T_train_C(:)';
            T_train_S = Net_opt.T_train_S(:)';
            load("./FNN/mat_files/MSE_FNN_trainscg_Calv_dh_"+yearstr(ii)+"_Weertman_cycle"+string(cycle)+"_N0k"+string(pct)+".mat");
            MSE = mse_FNN;
        else
            error("Unknown emulator type");
        end
        [nModes,nNodes] = size(T_reproj);
        M(ii).X_train_C = X_train_C;
        M(ii).X_train_S = X_train_S;
        M(ii).T_train_C = T_train_C;
        M(ii).T_train_S = T_train_S;
        % Subtract mean and project measurements onto truncated basis.
        % These are the unnormalized values.
        if dataformat == "du"
            yearstr2 = "yr"+years(ii)+"_yr"+years(ii+1);
            stdstr = "stddu";
            ustr = "du";
        else
            yearstr2 = "yr"+years(ii);
            stdstr = "stdu";
            ustr = "u";
        end
        M(ii).data = (measurements(1).(yearstr2).(ustr)(:)-T_mean(:))'*T_reproj';
        % Assemble covariance matrices
        % 1. measurement errors
        S_meas = T_reproj*spdiags(measurements(1).(yearstr2).(stdstr)(:).^2,0,nNodes,nNodes)*T_reproj';
        % 2. Ua errors: obtained from  FIXME
        s_u = 200; %m/yr
        if dataformat == "LOGu"
            s_u = log10(s_u);
        end
        l = 200e3; % range in the semivariogram
        load("u_AS_Calv_dh_cycle2_Weertman_2000_2009_2014_2018.mat","MUA");
        MUA = MUA.(yearstr2);
        D_tmp = pdist2(MUA.coordinates,MUA.coordinates,"squaredeuclidean");
        D_tmp = s_u^2*exp(-D_tmp/(2*l^2));
        D_tmp = D_tmp*T_reproj';
        S_ua = T_reproj*D_tmp;
        clear D_tmp
        % 3. Emulator errors: obtained from plot_Emulator_MSE
        S_emulator = spdiags(MSE(:),0,nModes,nModes);
        M(ii).S = S_meas + S_ua + S_emulator;
        M(ii).detS = det(M(ii).S);
        M(ii).Sinv = inv(M(ii).S);
    end    
end

% Loop through years and realizations
logL = zeros(nReal,1);

for ii=1:N
    % Apply normalization to parameters before feeding into emulator
    predictors = (params-repmat(M(ii).X_train_C,nReal,1))./repmat(M(ii).X_train_S,nReal,1);

    % Evaluate forward model. 
    if NN=="RNN"
        modelRun = double(predict(M(ii).net,predictors)); 
    elseif NN=="FNN"
        modelRun = double(M(ii).net(predictors')');
    end
    % Undo normalization of the output but keep in the projected basis to 
    % make it compatible with data
    modelRun = modelRun.*repmat(M(ii).T_train_S,nReal,1)+repmat(M(ii).T_train_C,nReal,1);

    % Assemble log likelihood
    for jj = 1:nReal
      % Evaluate log-likelihood
      logLikeli = - 1/2*log(2*pi*M(ii).detS) - 1/2*(M(ii).data...
        -modelRun(jj,:))*M(ii).Sinv*(M(ii).data-modelRun(jj,:))';
      % Assign to logL vector
      logL(jj) = logL(jj)+logLikeli;
    end
end

end




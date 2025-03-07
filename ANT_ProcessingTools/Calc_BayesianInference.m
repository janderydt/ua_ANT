function Calc_BayesianInference

Klear;

%% user defined parameters
domain = "AMUND";
slidinglaw = ["Weertman" "Umbi"];
cycle = [2 2]; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
dataformat = ["du" "du"]; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
only_grounded_ice = [1 1];
years = ["2000-2020" "2000-2020"]; % a vector with years for which velocity data is used
NN = ["RNN" "RNN"]; % which emulator? FNN or RNN

do_plots = 1;

% dataformat = ["du"]; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
% cycle = [1]; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
% pct = [9897]; % trunction of SVD in emulator
% only_grounded_ice = [1];
% years = ["2000-2009"]; % a vector with years for which velocity data is used
% NN = ["RNN"]; % which emulator? FNN or RNN

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
% that these are also the parameters in emulator prepare_data_for_***_emulators.m 
% for which the emulators have been trained. By default, we use
% X = [log10(gaA) log10(gaC) log10(gsA) log10(gsC) m n log10(dhdt_err)];
PriorOpts.Name = 'Model parameters prior';

ind = 1;
PriorOpts.Marginals(ind).Name = 'log10(gaA)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [0 0.5]; % mean and std
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [-1 log10(200)];
PriorOpts.Marginals(ind).Bounds = [-1 log10(200)];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gaC)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [0 0.5];
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [-1 log10(50)];
PriorOpts.Marginals(ind).Bounds = [-1 log10(50)];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gsA)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [4 1];
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [3 6];
PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'log10(gsC)';
PriorOpts.Marginals(ind).Type = 'Gaussian';
PriorOpts.Marginals(ind).Parameters = [4 1];
%PriorOpts.Marginals(ind).Type = 'Uniform';
%PriorOpts.Marginals(ind).Parameters = [3 6];
PriorOpts.Marginals(ind).Bounds = [3 6];
ind = ind+1;
PriorOpts.Marginals(ind).Name = 'm';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 9];

ind = ind+1;
PriorOpts.Marginals(ind).Name = 'n';
PriorOpts.Marginals(ind).Type = 'Uniform';
PriorOpts.Marginals(ind).Parameters = [2 5];

% add dhdt if any cycle>1
if max(cycle) > 1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'log(dhdt_err)';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [log10(0.05) log10(0.5)];
    %PriorOpts.Marginals(ind).Bounds = [log10(0.05) log10(0.5)];
end

% add discrete variables
nslidinglaw = numel(unique(slidinglaw));
ncycle = numel(unique(cycle));
if ncycle>1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'Cycle';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [0 1];
end
if nslidinglaw>1
    ind = ind+1;
    PriorOpts.Marginals(ind).Name = 'SlidingLaw';
    PriorOpts.Marginals(ind).Type = 'Uniform';
    PriorOpts.Marginals(ind).Parameters = [0 1];
end

myPriorDist = uq_createInput(PriorOpts);

%% define data
myData.y = loadvelocitydata(dataformat,years,only_grounded_ice); % myData.y.(year).u and myData.y.(year).stdu 
    % are N-by-1 matrices with N the number of nodes in the Ua mesh for that year
myData.Name = 'velocity observations and errors';

%% define custom likelihood function
% function takes parameter values (log10(gaA), log19(gaC), log10(gsA), 
% log10(gsC), m, n, log10(dhdt_err)) and data (y), and returns the log-likelihood
% function at these points. The forward model and covariance matrix are 
% defined within the function.

myLogLikelihood = @(params,y) customLogLikelihood(params, y, [domain,slidinglaw,cycle,dataformat,only_grounded_ice,years,NN]);

%% define solver options
mySolver.Type = 'MCMC';
mySolver.MCMC.Sampler = 'AIES'; % AM, HMS or AIES (default)
mySolver.MCMC.Steps = 5000; % T=300 default
mySolver.MCMC.NChains = 250; % C=100 default
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
    string(mySolver.MCMC.NChains)+"_"+dataformat(1)+"_Calv_dh_"+...
            strjoin(string(years),"_")+"_"+strjoin(string(slidinglaw),"_")+"_cycle"+string(cycle(1))+...
            "_floatingice"+string(1-only_grounded_ice(1));
save(fname(1),"myBayesianAnalysis");

%% print some results
uq_print(myBayesianAnalysis);
%uq_display(myBayesianAnalysis);

if do_plots
    plot_myBayesianAnalysis(myBayesianAnalysis);
end

end

function out = loadvelocitydata(dataformat,years,only_grounded_ice)

% load measured velocity data
%tmp=load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat","Fus","Fvs","Fxerr","Fyerr");
%Fus_2000=tmp.Fus; Fvs_2000=tmp.Fvs; Fxerr_2000=tmp.Fxerr; Fyerr_2000=tmp.Fyerr;
addpath(getenv("froot_data")+"Measures/Measures_annual");

% assemble velocity fields
years_tmp = strjoin(years,"-");
years_unique = unique(split(years_tmp,"-"));

for ii=1:numel(years_unique)

    % years_measures = string(years(ii))+"_"+string(years(ii)+1);
    % [vx,x_meas,y_meas]=measures_annual("vx",years_measures); 
    % vy=measures_annual("vy",years_measures);
    % stdx=measures_annual("vxstd",years_measures);
    % stdy=measures_annual("vystd",years_measures);
    % y_meas = flip(y_meas,1); [Xm,Ym]=ndgrid(x_meas,y_meas);
    % v(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(vx,vy),2));
    % std(ii).F = griddedInterpolant(Xm,Ym,flip(hypot(stdx,stdy),2));
    if years_unique(ii)=="2000"
        fname = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    elseif years_unique(ii)=="2020"
        fname = "GriddedInterpolants_"+string(double(years_unique(ii))-1)+"-"+string(double(years_unique(ii)))+...
            "_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    else
        fname = "GriddedInterpolants_"+string(double(years_unique(ii)))+"-"+string(double(years_unique(ii))+1)+...
            "_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    end
    load("../ANT_Data/ANT_Interpolants/"+fname);
    v_tmp = hypot(Fus.Values,Fvs.Values);
    Fu = Fus; Fu.Values = v_tmp;
    std_tmp = hypot(Fxerr.Values,Fyerr.Values);
    Fstd = Fus; Fstd.Values = std_tmp;
    v.("yr"+years_unique(ii)).F = Fu;
    std.("yr"+years_unique(ii)).F = Fstd;

end

for dd=1:numel(dataformat)

    if dataformat(dd) == "du"
        
        years_tmp = split(years(dd),"-");
        yr1 = years_tmp(1);
        yr2 = years_tmp(2);

        load("Delta_u_AMUND_Weertman_"+yr1+"-"+yr2+".mat","MUA_yr2","GF_yr2");
        MUA = MUA_yr2; GF = GF_yr2;
    
        % interpolate initial and final speed onto Ua mesh
        u_init = v.("yr"+yr1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_init = std.("yr"+yr1).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        u_target = v.("yr"+yr2).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_target = std.("yr"+yr2).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
    
        deltau = u_target-u_init; % measured change in speed (m/yr)
        deltau_err = hypot(std_init,std_target); % error estimate (m/yr)
        
        % remove nans
        Ind = find(isnan(deltau) | isnan(deltau_err));
        deltau(Ind) = eps(0);
        deltau_err(Ind) = 1e3;        
    
        %% what if we remove the floating ice?
        if only_grounded_ice
            deltau(GF.node<0.5) = eps(0);
            deltau_err(GF.node<0.5) = 1e-3;
        end
    
        out.("yr"+yr1+"_yr"+yr2).du = deltau(:);
        out.("yr"+yr1+"_yr"+yr2).stddu = deltau_err(:);

    elseif ismember(dataformat(dd),["u","LOGu"])
    
        load("u_AS_Calv_dh_cycle2_Weertman_"+years(dd)+".mat","MUA","GF");
        MUA = MUA.("yr"+years(dd)); GF = GF.("yr"+years(dd));

        % interpolate initial and final speed onto Ua mesh
        u_tmp = v.("yr"+years(dd)).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        std_tmp = std.("yr"+years(dd)).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
        
        % remove nans
        %Ind = find(isnan(u_tmp) | isnan(std_tmp));
        %std_tmp(Ind) = 5e3;
        %u_tmp(Ind) = eps(0);

        % apply log if needed
        if dataformat(dd)=="LOGu"
            std_tmp = min(std_tmp./(u_tmp*log(10)),10);
            u_tmp = log10(u_tmp);    
        end

        out.("yr"+years(dd)).u = u_tmp(:);
        out.("yr"+years(dd)).stdu = std_tmp(:);
    
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
% mymodel has format [domain,slidinglaw,cycle,dataformat,only_grounded_ice,years,NN],
% each of which is an 1xn vector apart from domain. 
domain = mymodel(1);
N = (numel(mymodel)-1)/6;
slidinglaw = mymodel(2:1+N);
cycle = double(mymodel(2+N:1+2*N));
dataformat = mymodel(2+2*N:1+3*N);
only_grounded_ice = double(mymodel(2+3*N:1+4*N));
years = mymodel(2+4*N:1+5*N);
NN = mymodel(2+5*N:1+6*N);

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
        s_u = 200; %m/yr
        % if dataformat == "LOGu"
        %     s_u = log10(s_u);
        % end
        l = 200e3; % range in the semivariogram
        load("Delta_u_AMUND_Weertman_"+yearstr+".mat","MUA_yr2");
        MUA = MUA_yr2;
        D_tmp = pdist2(MUA.coordinates,MUA.coordinates,"squaredeuclidean");
        D_tmp = s_u^2*exp(-D_tmp/(2*l^2));
        D_tmp = D_tmp*T_reproj';
        S_ua = 0*T_reproj*D_tmp;
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
params_tmp = params;
C = zeros(nReal,1); ind=1;

if ncycles > 1    
    if max(cycle)>1
        C(:,ind)=round(params(:,8));
        params_tmp(:,8)=[];
    else
        C(:,ind)=round(params(:,7));
        params_tmp(:,7)=[];
    end
    ind=ind+1;
end

if nslidinglaws > 1
    C(:,ind)=round(params(:,end));
    params_tmp(:,end)=[];
end

EmulatorIndex = 1;
for ii=1:size(C,2)
    EmulatorIndex = EmulatorIndex+2^(size(C,2)-1)*C(:,ii);
end
params=params_tmp;

for ii=unique(EmulatorIndex(:)')

    Itmp = find(EmulatorIndex==ii);
    nItmp = numel(Itmp);

    % Apply normalization to parameters before feeding into emulator
    if cycle(ii)==1
        predictors = (params(Itmp,1:6)-repmat(M(ii).X_train_C,nItmp,1))./repmat(M(ii).X_train_S,nItmp,1);
    else
        predictors = (params(Itmp,:)-repmat(M(ii).X_train_C,nItmp,1))./repmat(M(ii).X_train_S,nItmp,1);
    end

    % Evaluate forward model. 
    if NN(ii)=="RNN"
        modelRun = double(predict(M(ii).net,predictors)); 
    elseif NN(ii)=="FNN"
        modelRun = double(M(ii).net(predictors')');
    end

    % Undo normalization of the output but keep in the projected basis to 
    % make it compatible with data
    modelRun = modelRun.*repmat(M(ii).T_train_S,nItmp,1)+repmat(M(ii).T_train_C,nItmp,1);
    Nout = size(modelRun,2);

    % Assemble log likelihood
    for jj = 1:nItmp
      % Evaluate log-likelihood
      prefac = -0.5*log((2*pi)^Nout*M(ii).detS);
      Ndistrib = - 0.5*(M(ii).data...
        -modelRun(jj,:))*M(ii).Sinv*(M(ii).data-modelRun(jj,:))';
      logLikeli = prefac + Ndistrib;
      % Assign to logL vector
      logL(Itmp(jj)) = logL(Itmp(jj))+logLikeli;
    end

    %figure(115); hold on; plot(params(1,5),prefac(1),'.',Color=CM(ii,:)); title('m');
    %plot(params(1,5),Ndistrib(1),'.',Color=CM(ii+2,:)); title('m');

end

%figure(111); hold on; plot(params(1,1),logL(1),'.k'); title('gaA');
%figure(112); hold on; plot(params(1,2),logL(1),'.k'); title('gaC');
%figure(113); hold on; plot(params(1,3),logL(1),'.k'); title('gsA');
%figure(114); hold on; plot(params(1,4),logL(1),'.k'); title('gsC');
%figure(116); hold on; plot(params(1,6),logL(1),'.k'); title('n');

end




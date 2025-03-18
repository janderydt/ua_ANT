function Driver_BayesianInference

l = [1 50e3*ones(1,4) 100e3*ones(1,4) 150e3*ones(1,4) 200e3*ones(1,4)];
fac_su = [0 1:4 1:4 1:4 1:4];

% for ii=1:numel(l)
%     UserVar.l = l(ii); % range in the ua semivariogram
%     UserVar.fac_su = fac_su(ii); % multiplication factor for model errors compared to observation errors
%     UserVar.domain = "AMUND";
%     UserVar.slidinglaw = ["Weertman" "Weertman" "Weertman" "Umbi" "Umbi" "Umbi"];
%     UserVar.cycle = [1 1 1 1 1 1]; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
%     UserVar.dataformat = ["du" "du" "du" "du" "du" "du"]; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
%     UserVar.only_grounded_ice = [1 1 1 1 1 1];
%     UserVar.years = ["2000-2009" "2000-2014" "2000-2020" "2000-2009" "2000-2014" "2000-2020"]; % a vector with years for which velocity data is used
%     UserVar.NN = ["RNN" "RNN" "RNN" "RNN" "RNN" "RNN"]; % which emulator? FNN or RNN
%     UserVar.do_plots = 1;
%     Calc_BayesianInference(UserVar);
%     clear Calc_BayesianInference;
%     close all;
% end

for ii=1:numel(l)
    UserVar.l = l(ii); % range in the ua semivariogram
    UserVar.fac_su = fac_su(ii); % multiplication factor for model errors compared to observation errors
    UserVar.domain = "AMUND";
    UserVar.slidinglaw = ["Weertman" "Weertman" "Weertman" "Umbi" "Umbi" "Umbi"];
    UserVar.cycle = [2 2 2 2 2 2]; % without spinup (cycle=1) or with spinup and dhdt (cycle=2)
    UserVar.dataformat = ["du" "du" "du" "du" "du" "du"]; % use speed ("u"), log of speed ("LOGu") or change in speed ("du")
    UserVar.only_grounded_ice = [1 1 1 1 1 1];
    UserVar.years = ["2000-2009" "2000-2014" "2000-2020" "2000-2009" "2000-2014" "2000-2020"]; % a vector with years for which velocity data is used
    UserVar.NN = ["RNN" "RNN" "RNN" "RNN" "RNN" "RNN"]; % which emulator? FNN or RNN
    UserVar.do_plots = 1;
    Calc_BayesianInference(UserVar);
    clear Calc_BayesianInference;
    close all;
end


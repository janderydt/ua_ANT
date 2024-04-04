function ANT_GenerateSampleSpace

addpath(genpath('/mnt/md0/Matlab/UQLab_Rel2.0.0'));

clearvars;
rng(1,'twister'); % set the random number generator for reproducible results
uqlab; % initialize uqlab

Input.Name = 'Parameter array for inverse simulations';
Input.Marginals = uq_Marginals(8,'Uniform',[0]);
ind = 1;

%% -------------------- %%
%% CONTINUOUS VARIABLES %%
%% -------------------- %%

%% Regularization
Input.Marginals(ind).Name = 'gsC';
Input.Marginals(ind).Parameters = [10e3 500e3];
Input.Marginals(ind).Bounds = [10e3 500e3];
ind = ind + 1;

Input.Marginals(ind).Name = 'gsA';
Input.Marginals(ind).Parameters = [10e3 500e3]; % parameters are bounds
Input.Marginals(ind).Bounds = [10e3 500e3];
ind = ind + 1;

Input.Marginals(ind).Name = 'gaC';
Input.Marginals(ind).Parameters = [1 200];
Input.Marginals(ind).Bounds = [1 200];
ind = ind + 1;

Input.Marginals(ind).Name = 'gaA';
Input.Marginals(ind).Parameters = [1 200]; % parameters are bounds
Input.Marginals(ind).Bounds = [1 200];
ind = ind + 1;

%% Sliding law
Input.Marginals(ind).Name = 'm';
Input.Marginals(ind).Parameters = [2 9];
Input.Marginals(ind).Bounds = [2 9];
ind = ind + 1;

Input.Marginals(ind).Name = 'ubprior';
Input.Marginals(ind).Parameters = [10 200];
Input.Marginals(ind).Bounds = [10 200];
ind = ind + 1;

%% Ice Rheology
Input.Marginals(ind).Name = 'n';
Input.Marginals(ind).Parameters = [2 4];
Input.Marginals(ind).Bounds = [2 4];
ind = ind + 1;

Input.Marginals(ind).Name = 'epsprior';
Input.Marginals(ind).Parameters = [1 6]*1e-3;
Input.Marginals(ind).Bounds = [1 6]*1e-3;

myInput = uq_createInput(Input);

uq_print(myInput);

uq_selectInput(myInput);
X = uq_getSample(numel(Input.Marginals)*10,'LHS');

%uq_display(myInput);

ANT_CreateNewRunsTable(X);

return
%% ------------------ %%
%% DISCRETE VARIABLES %%
%% ------------------ %%

%% Spinup
Name = 'spinup';
Values = [0 1];

%% Cost function
Name = 'dhdt';
Values = [0 1];

%% Mesh
Name = 'meshrefinement';
Values = [0 1];

%% Coulomb
Name = 'coulomb';
Values = [0 1];

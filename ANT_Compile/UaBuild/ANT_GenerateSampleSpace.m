function ANT_GenerateSampleSpace

addpath(genpath(getenv("froot_matlabfunctions")+"/../UQLab_Rel2.0.0"));

clearvars;
rng(1,'twister'); % set the random number generator for reproducible results
uqlab; % initialize uqlab

GradientCalc = "Adjoint"; % options: FixPoint or Adjoint
SlidingLaw = "Weertman"; 
Enrich = 0; % enriches the existing Latin Hypercube
EnrichSampleSize = 180;

Input.Name = 'Parameter array for inverse simulations';
ind = 1;

%% -------------------- %%
%% CONTINUOUS VARIABLES %%
%% -------------------- %%
switch GradientCalc
    case "Adjoint"

        Input.Marginals = uq_Marginals(6,'Uniform',[0]);

        %% Regularization
        Input.Marginals(ind).Name = 'gsC';
        %Input.Marginals(ind).Parameters = log10([25e3 50e3]); %v1: [50e3 1000e3]
        %Input.Marginals(ind).Bounds = log10([25e3 50e3]); %v1: [50e3 1000e3]
        Input.Marginals(ind).Parameters = [0 (log(50e3)-log(1e3))/log(2)];
        Input.Marginals(ind).Bounds = [0 (log(50e3)-log(1e3))/log(2)];
        ind = ind + 1;
        
        Input.Marginals(ind).Name = 'gsA';
        %Input.Marginals(ind).Parameters = log10([25e3 50e3]);  %v1: [50e3 1000e3]
        %Input.Marginals(ind).Bounds = log10([25e3 50e3]); %v1: [50e3 1000e3]
        Input.Marginals(ind).Parameters = [0 (log(50e3)-log(1e3))/log(2)];
        Input.Marginals(ind).Bounds =[0 (log(50e3)-log(1e3))/log(2)];
        ind = ind + 1;
        
        Input.Marginals(ind).Name = 'gaC';
        %Input.Marginals(ind).Parameters = log10([1 100]);
        %Input.Marginals(ind).Bounds = log10([1 100]);
        Input.Marginals(ind).Parameters = [0 (log(100)-log(0.1))/log(2)];
        Input.Marginals(ind).Bounds = [0 (log(100)-log(0.1))/log(2)];
        ind = ind + 1;
        
        Input.Marginals(ind).Name = 'gaA';
        %Input.Marginals(ind).Parameters = log10([1 250]); % parameters are bounds
        %Input.Marginals(ind).Bounds = log10([1 250]);
        Input.Marginals(ind).Parameters = [0 (log(250)-log(0.1))/log(2)];
        Input.Marginals(ind).Bounds = [0 (log(250)-log(0.1))/log(2)];
        ind = ind + 1;
        
        %% Sliding law
        Input.Marginals(ind).Name = 'm';
        Input.Marginals(ind).Parameters = [2 9];
        Input.Marginals(ind).Bounds = [2 9];
        %Input.Marginals(ind).Moments = [0];
        ind = ind + 1;
        
        %Input.Marginals(ind).Name = 'ubprior';
        %Input.Marginals(ind).Parameters = [100];
        %Input.Marginals(ind).Bounds = [100];
        %ind = ind + 1;
        
        %% Ice Rheology
        Input.Marginals(ind).Name = 'n';
        Input.Marginals(ind).Parameters = [2 4];
        Input.Marginals(ind).Bounds = [2 4];
        %Input.Marginals(ind).Moments = [0];
        %ind = ind + 1;
        
        %Input.Marginals(ind).Name = 'epsprior';
        %Input.Marginals(ind).Parameters = [2.6]*1e-3;
        %Input.Marginals(ind).Bounds = [2.6]*1e-3;
        % 
        

        %% Latin hypercube
        myInput = uq_createInput(Input);
        uq_print(myInput);
        
        uq_selectInput(myInput);
        X = uq_getSample(numel(Input.Marginals)*30,'LHS');

        X(:,1) = 1e3.*2.^X(:,1);
        X(:,2) = 1e3.*2.^X(:,2); 
        X(:,3) = 0.1.*2.^X(:,3);
        X(:,4) = 0.1.*2.^X(:,4);
        figure; hold on;
        for ii=1:6
            for jj=1:6
                subplot(6,6,(ii-1)*6+jj);
                if ii==jj
                    histogram(X(:,ii),15)
                else
                    plot(X(:,ii),X(:,jj),'ok');
                end
            end
        end

    case "FixPoint"

        %% Regularization: does not matter here, but need to provide some numbers
        Input.Marginals(ind).Name = 'gsC';
        gsC = 1;
        ind = ind + 1;
        
        Input.Marginals(ind).Name = 'gsA';
        gsA = 1;
        ind = ind + 1;
        
        Input.Marginals(ind).Name = 'gaC';
        gaC = 1;
        ind = ind + 1;
        
        Input.Marginals(ind).Name = 'gaA';
        gaA = 1;
        ind = ind + 1;


        %% Sliding law      
        Input.Marginals(ind).Name = 'm';
        m = [2:7];
        ind = ind + 1;

        %% Ice Rheology
        Input.Marginals(ind).Name = 'n';
        n = [2.5:0.5:3.5];

        %% ND grid
        [m,n]=ndgrid(m,n); ny=numel(m(:));
        X = [gsC*ones(ny,1) gsA*ones(ny,1) gaC*ones(ny,1) gaA*ones(ny,1) m(:) n(:)];

end

if Enrich

    Xnew = uq_enrichLHS(X,EnrichSampleSize,myInput);
    X = Xnew;

end

T=array2table(X,'VariableNames',{Input.Marginals(:).Name});

ANT_GenerateRunTable(T,GradientCalc,SlidingLaw,Enrich);

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

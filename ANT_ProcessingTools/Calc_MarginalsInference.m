function Calc_MarginalsInference

Klear

%% initialize UQLAB
addpath(genpath(getenv("froot_matlabfunctions")+"/../UQLab_Rel2.0.0"));

rng(1,'twister'); % set the random number generator for reproducible results
uqlab; % initialize uqlab

%% load sample from existing Bayesian Analysis
[fname,location] = uigetfile("/mnt/md0/Ua/cases/ANT/ANT_ProcessingTools/BayesianAnalysis/*.mat","MultiSelect","on");
fname = string(fname);

for ff=1:numel(fname)

    if numel(fname)>1
        load(location+"/"+fname(ff));
    else
        load(location+"/"+fname);
    end
    
    if ~exist("myPosteriorDist","var")

        mySolver = myBayesianAnalysis.Options.Solver;
        myPriorDist = myBayesianAnalysis.PriorDist;
        myResults = myBayesianAnalysis.Results;
        isdiscrete = 0;
        
        for ii=1:numel(myPriorDist.Marginals)
            Sample_tmp = [];
            % Sample equilibrium distribution from Markov Chains; we use the final 
            % 1% of the Markov Chains, i.e. discard the first 99% as burn-in
            for jj = 1:mySolver.MCMC.NChains
                Sample_tmp = [Sample_tmp; ...
                    myResults.Sample(mySolver.MCMC.Steps-round(0.01*mySolver.MCMC.Steps):end,ii,jj)];
            end
            X(:,ii) = Sample_tmp(:);        
            if ismember(myPriorDist.Marginals(ii).Name,["Cycle","SlidingLaw"])
                I_discr1 = find(X(:,ii)<0.5);
                I_discr2 = find(X(:,ii)>=0.5);
                X(:,ii)=[];
                isdiscrete = 1;
            end
        end
        
        %% For discrete variables, split distributions
        X1=[]; X2=[];
        if isdiscrete
            for ii=1:size(X,2)
                X1(:,ii) = X(I_discr1,ii);
                X2(:,ii) = X(I_discr2,ii);
            end
        end

        %% Calculate marginals/copulas fully automated
        iOpts1=[]; iOpts2=[];
        if isdiscrete
            % discrete set 1
            iOpts1.Inference.Data = X1;
            iOpts1.Copula.Type = 'auto';
            for ii=1:numel(myPriorDist.Marginals)-1
                iOpts1.Marginals(ii).Bounds = myPriorDist.Marginals(ii).Bounds;
                iOpts1.Marginals(ii).Type = {'Gaussian','Laplace'};
            end            
            if ~isempty(X1)
                myPosteriorDist(1).dist = uq_createInput(iOpts1);
                myPosteriorDist(1).frac = numel(I_discr1)/size(X,1);
            else
                myPosteriorDist(1).dist = [];
                myPosteriorDist(1).frac = [];
            end
            %  discrete set 2
            iOpts2.Inference.Data = X2;
            iOpts2.Copula.Type = 'auto';
            for ii=1:numel(myPriorDist.Marginals)-1
                iOpts2.Marginals(ii).Bounds = myPriorDist.Marginals(ii).Bounds;
                iOpts2.Marginals(ii).Type = {'Gaussian','Laplace'};
            end           
            if ~isempty(X2)
                myPosteriorDist(2).dist = uq_createInput(iOpts2);   
                myPosteriorDist(2).frac = numel(I_discr2)/size(X,1);
            else
                myPosteriorDist(2).dist = [];
                myPosteriorDist(2).frac = [];
            end
        else
            iOpts.Inference.Data = X;
            iOpts.Copula.Type = 'auto';
            for ii=1:numel(myPriorDist.Marginals)
                iOpts.Marginals(ii).Bounds = myPriorDist.Marginals(ii).Bounds;
                iOpts.Marginals(ii).Type = {'Gaussian','Laplace'};
            end            
            myPosteriorDist.dist = uq_createInput(iOpts);     
            myPosteriorDist.frac = 1;
        end

        save(location+"/"+fname(ff),"myPosteriorDist","-append");  

    end

    clear myBayesianAnalysis UserVar myPosteriorDist X f

end
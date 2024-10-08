function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

fprintf('Reading inputs for inverse run \n');

load(UserVar.VelocityInterpolants);

x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);

if strfind(CtrlVar.Inverse.Measurements,'dhdt')
    
    % to speed up the dhdt assembly, we check if the fields have
    % previously been constructed for the same mesh. If not, we read the
    % interpolants and construct a new file for this particular mesh.
    filename_dhdtfields = UserVar.dhdtInterpolants;
    filename_dhdtfields = erase(filename_dhdtfields,["GriddedInterpolants_","ScatteredInterpolants_",".mat"]);
    filename_dhdtfields = filename_dhdtfields + "_mesh_Nnodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

    fprintf('Loading dhdt fields for grounded ice.\n');

    fprintf('Trying to read fields from %s...',filename_dhdtfields);

    if exist(filename_dhdtfields,"file")

        load(filename_dhdtfields);

    else

        load(UserVar.dhdtInterpolants);
        dhdtMeas=Fdhdt_CPOM(x,y);
        dhdtError=Fdhdt_err_CPOM(x,y);

        save(filename_dhdtfields,"dhdtMeas","dhdtError");

    end
     
%     if strfind(CtrlVar.Experiment,'GL')
%         % small errors for dhdt at GL and on floating parts
%         fprintf('Setting small error for dhdt close to GL and for floating parts \n');
%         ds=10000;
%         CtrlVar.PlotGLs=0;
%         [xGL,yGL,~] = PlotGroundingLines(CtrlVar,MUA,GF,[],[],[]);
%         ID=FindAllNodesWithinGivenRangeFromGroundingLine([],MUA,xGL,yGL,ds);
%         I = find(GF.node<1);
%         dhdtError([I(:);ID(:)])=1e-1;
%     end

    % set zero values and high errors for floating areas
    dhdtMeas(GF.node<0.5) = 0;
    dhdtError(GF.node<0.5) = 1e10;
     
    Meas.dhdt=dhdtMeas;

     % Apply scaling of the dhdt errors and set minimum value of the error 
     % to 0.1. This is a somewhat arbitrary number, but previous 
     % experiments have shown that for smaller errors, the inversion tends 
     % to overly constrain dhdt, which leads to 'noisy' AGlen and C fields,
     % and a larger misfit to the surface velocities. It is questionable 
     % anyhow if dhdt_err<0.1m can be achieved by satellite instruments.
     dhdtError = max(UserVar.Inverse.dhdt_err*dhdtError,0.1);
    
     Meas.dhdtCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,dhdtError.^2,MUA.Nnodes,MUA.Nnodes);
     if any(isnan(dhdtError))
         error('NaN in dhdtError'); 
     end
 end

uMeas = Fus(x,y);
vMeas = Fvs(x,y);
xerrMeas = Fxerr(x,y);
yerrMeas = Fyerr(x,y);

clearvars Fus Fvs Fxerr Fyerr

% we don't want zero errors:
xerrMeas(xerrMeas==0) = 0.1;
yerrMeas(yerrMeas==0) = 0.1;

if contains(CtrlVar.Inverse.DataMisfit.GradientCalculation,'FixPoint','IgnoreCase',true) || CtrlVar.SlidingLaw ~= "Weertman"
    %% for fixpoint inversion
    % Siple Coast
    I = find(x>-600942 & x<8245 & y<-286360 & y>-800000 & ~isnan(uMeas));
    xerrMeas(I) = 0.01;
    yerrMeas(I) = 0.01;
    I = find(x>-800000 & x<-450000 & y<-700000 & y>-950000 & ~isnan(uMeas));
    xerrMeas(I) = 0.01;
    yerrMeas(I) = 0.01;
    % Shackleton Range
    I = find(x>-6e5 & x<-1.5e5 & y>7.5e5 & y<10.5e5);
    xerrMeas(I) = 0.1;
    yerrMeas(I) = 0.1;
    % Pensacola Mountains
    I = find(x>-10e5 & x<-8e5 & y>1.5e5 & y<3.5e5);
    xerrMeas(I) = 0.1;
    yerrMeas(I) = 0.1;
    % Elsworth Mountains and Rutford
    I = find(x>-14.5e5 & x<-12e5 & y>0.5e5 & y<2.5e5);
    xerrMeas(I) = 0.1;
    yerrMeas(I) = 0.1;
    % Pine Island
    I = find(x>-1.7e6 & x<-1.5e6 & y>-2.3e5 & y<0);
    xerrMeas(I) = 0.1;
    yerrMeas(I) = 0.1;
    % Foundation Ice Stream
    I = find(x>-7e5 & x<-5e5 & y>2e5 & y<4e5);
    xerrMeas(I) = 0.1;
    yerrMeas(I) = 0.1;
    % Foundation Ice Stream
    I = find(x>-1.6e6 & x<-1.4e6 & y>2.5e5 & y<4.5e5);
    xerrMeas(I) = 0.1;
    yerrMeas(I) = 0.1;
else
    if UserVar.Geometry == 2000
        % Thwaites
        I = find(x>-1.65e6 & x<-1.45e6 & y>-5.5e5 & y<-4e5);
        xerrMeas(I) = xerrMeas(I)+1;
        yerrMeas(I) = yerrMeas(I)+1;
    end
end

if contains(UserVar.Domain,'AS_PROPHET')
    xerrMeas = 0*x+1;
    yerrMeas = 0*x+1;
end

% Because I put NaN in where there is no data, this will give NaN at location not surounded by four data points
% I can then afterwards find the NaN and put in some data with very high errors.
% I need values everywhere, so I now set vel to zero where no measurments are available, 
% and the error to a very large value

Mask=isnan(uMeas);
uMeas(Mask)=0;
xerrMeas(Mask)=1e10;

Mask=isnan(vMeas);
vMeas(Mask)=0;
yerrMeas(Mask)=1e10;

Meas.us = uMeas;
Meas.vs = vMeas;

% figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(sqrt(Meas.us.^2+Meas.vs.^2))); hold on;
% drawnow;
% figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,errMeas); hold on;
% drawnow;

usError = xerrMeas;
vsError = yerrMeas;

if any(isnan(Meas.us))
	error('NaN in Meas.us');
elseif any(isnan(Meas.vs))
    error('NaN in Meas.vs');
elseif any(isnan(Meas.dhdt))
    error('NaN in Meas.dhdt');
elseif any(isnan(usError))
    error('NaN in usError');
 elseif any(isnan(vsError))
    error('NaN in vsError');   
end

Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);

clearvars xerrMeas yerrMeas uMeas vMeas

fprintf('Reading start values for AGlen and C...\n');
[UserVar,InvStartValues.C,InvStartValues.m,InvStartValues.q,InvStartValues.muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
[UserVar,InvStartValues.AGlen,InvStartValues.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
fprintf('done \n');

% original=[CtrlVar.Experiment,'_InverseResults_Intermediate.mat'] ;
% new=[CtrlVar.Experiment,'_InverseResults_Intermediate_',date,'.mat'];
% if exist(original)==2
%     copyfile(original,new);
%     Frestart=load(original,'F');
%     InvStartValues.C = Frestart.F.C;
%     InvStartValues.AGlen = Frestart.F.AGlen;
% end

listingCC=dir('CC.mat') ; listingCA=dir('CAGlen.mat') ;
CC = []; CAGlen = []; 

%%  Covariance matrices of priors
% 
% if CtrlVar.AGlenisElementBased
%     CAGlen=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
% else
%     CAGlen=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
% end
% 
% if strcmpi(CtrlVar.Inverse.Regularize.Field,'cov')
% 
%     Err=1e-2 ; Sigma=1e3 ; DistanceCutoff=10*Sigma;
% 
%     if CtrlVar.CisElementBased
%         [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
%     else
%         [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
%     end
% 
% else
%     if CtrlVar.CisElementBased
%         CC=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
%     else
%         CC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
%     end
% end

%% Define Priors
Priors.B=F.B;

if ~isfield(UserVar.Inverse,'priorC')
    ub=100; tau=80 ; % units meters, year , kPa
    Priors.C = x*0 + ub/tau^InvStartValues.m(1);
else
    Priors.C = x*0 + UserVar.Inverse.priorC;
end
Priors.m = InvStartValues.m(1);
Priors.muk = InvStartValues.muk(1);

if ~isfield(UserVar.Inverse,'priorAGlen')
    Priors.AGlen = x*0 + AGlenVersusTemp(-15); % this is for n=3 only
else
    Priors.AGlen = x*0 + UserVar.Inverse.priorAGlen;
end
Priors.n=InvStartValues.n(1);

disp(['Constant prior value for C: ',num2str(Priors.C(1))]);
disp(['Constant prior value for AGlen: ',num2str(Priors.AGlen(1))]);

Priors.rho = F.rho;
Priors.rhow = F.rhow;
% Priors.CovAGlen = CAGlen;
% Priors.CovC = CC;

end

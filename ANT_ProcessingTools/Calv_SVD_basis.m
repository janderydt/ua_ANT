function Calv_SVD_basis(data,SVD_settings,SVD_filetowrite)

if nargin==0
    %% -- perturbation data to use
    domain = "AMUND";
    perturbation = "Calv_dh"; % valid choices are "Calv", "dhIS", "dh" or "Calv_dh"
    startyear = "2000";
    targetyear = "2014";
    slidinglaw = ["Weertman" "Umbi"]; % valid choices are "Weertman" or "Umbi"
    cycle = [2 2]; % cycle 1: inversion without spinup, cycle 2: inversion after spinup
    only_grounded_ice = 0;
    
    %% -- SVD settings
    add_measurements_to_SVD = 0;

    %% filename
    strtmp = "";
    for ii=1:numel(slidinglaw)
        strtmp2 = char(slidinglaw(ii));
        strtmp = strtmp+strtmp2(1)+"c"+string(cycle(ii))+"_";
    end
    SVD_filetowrite = "SVD_"+domain+"_"+perturbation+"_"+startyear+"-"+targetyear+"_"+...
        strtmp+"floatingice"+string(1-only_grounded_ice)+"_includemeasurements"+string(add_measurements_to_SVD)+".mat";

else
    if nargin<2
        error("Need 2 inputs.");
    else
        domain = data.domain;
        perturbation = data.perturbation;
        startyear = data.startyear;
        targetyear = data.targetyear;
        slidinglaw = data.slidinglaw;
        cycle = data.cycle;
        only_grounded_ice = data.only_grounded_ice;
        add_measurements_to_SVD = SVD_settings.add_measurements_to_SVD;
    end
end

X = []; T = [];

for ii=1:numel(slidinglaw)

    load("Delta_u_"+domain+"_"+slidinglaw(ii)+"_"+startyear+"-"+targetyear+".mat");

    %% PREDICTORS
    Xtmp = [log10(gaA(:)) log10(gaC(:)) log10(gsA(:)) log10(gsC(:)) m(:) n(:)];
    if cycle(ii)==2 % add dhdt_err to predictors
        Xtmp = [Xtmp log10(dhdt_err(:))];
    elseif cycle(ii)==1 && max(cycle)>1
        Xtmp = [Xtmp 0*m(:)];
    end
    Xtmp = double(Xtmp);

    %% TARGETS (training/input data). In this case, input data consists of 
    %% simulated instantaneous changes in surface speed in response to changes
    %% in ice-sheet geometry (ice thickness, calving front location). The input
    %% data consists of num_nodes nodal values for nun_exp experiments.
    Ttmp = Delta_u.(perturbation).map(:,:,cycle(ii));

    % check for nans and inf
    if any(isnan(Ttmp) | isinf(Ttmp))
        warning("Removing nan from T");
        [rows_to_delete,~]=find(isnan(Ttmp));
        rows_to_delete = unique(rows_to_delete);
        Xtmp(rows_to_delete,:)=[];
        Ttmp(rows_to_delete,:)=[];
    end
    % check dimensions of input data; rows: experiments, columns: nodes
    if numel(size(Ttmp))==2
        num_exp = size(Ttmp,1); num_nodes = size(Ttmp,2);
        if num_nodes < num_exp % number of experiments is highly likely to be smaller than number of nodes
            Ttmp = Ttmp'; % transpose data
            num_exp = size(Ttmp,1); num_nodes = size(Ttmp,2);
        end
    else
        error("check dimensions of input data");
    end

    %% ASSEMBLE FULL MATRICES
    X = [X; Xtmp];
    T = [T; Ttmp];
end

% remove mean
T_mean = mean(T,1);
T_meanremoved = T-repmat(T_mean(:)',size(T,1),1);

%% Load Ua target grid
load("Delta_u_"+domain+"_Weertman_"+startyear+"-"+targetyear+".mat","MUA_yr2","GF_yr2");
MUA = MUA_yr2; GF = GF_yr2;

%% Load observations 
for yy=[startyear targetyear]
    if yy=="2000"
        fname = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    else
        fname = "GriddedInterpolants_"+string(double(yy))+"-"+string(double(yy)+1)+...
        "_MeaSUREs_ITSLIVE_Velocities_EXTRUDED.mat";
    end
    load("../ANT_Data/ANT_Interpolants/"+fname);
    v_tmp = hypot(Fus.Values,Fvs.Values);
    Fu = Fus; Fu.Values = v_tmp;
    std_tmp = hypot(Fxerr.Values,Fyerr.Values);
    Fstd = Fus; Fstd.Values = std_tmp;
    v.("yr"+yy).F = Fu;
    std.("yr"+yy).F = Fstd;
end

% interpolate initial and final speed onto Ua mesh
u_init = v.("yr"+startyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
std_init = std.("yr"+startyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
u_target = v.("yr"+targetyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
std_target = std.("yr"+targetyear).F(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltau = u_target-u_init; % measured change in speed (m/yr)

% subtract model ensemble mean
T_meas = deltau(:);

if any(isnan(T_meas))
    error("dT_meas contains nans. Consider using extruded field or inpaint_nans.");
end

if add_measurements_to_SVD    
    T = [T_meas(:)' T];
    T_meanremoved = T-repmat(T_mean(:)',size(T,1),1);
end

%% remove floating nodes   
if only_grounded_ice   
    Nodes_floating = find(GF.node<0.5);
    T(:,Nodes_floating) = 0;
    T_meas(Nodes_floating) = 0;
    T_mean(Nodes_floating) = 0;
end

[~,S,~] = svd(T_meanremoved,'econ');

%% save data
save(SVD_filetowrite,"S","X","T","T_mean","T_meanremoved","T_meas");

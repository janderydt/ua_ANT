function ab = CalcMelt_UPP(UserVar,MUA,F,CtrlVar,doplots)

if nargin<4
    folder = "/mnt/md0/Ua/cases/ANT/ANT_Inverse/cases/ANT_nsmbl_Inverse_14423";
    file = "ANT_nsmbl_Inverse_14423-RestartFile.mat";
    load(folder+"/"+file);
    CtrlVar = CtrlVarInRestartFile;
    UserVar = UserVarInRestartFile;
    doplots=1;
end

% definition of physics parameters
betaT=3.87e-5; % thermal expansion coefficient
betaS=7.86e-4; % haline contraction coefficient
c=3.974e3; % specific heat capacity of water
Cd=2.5e-3; % drag coefficient
g=9.81; % gravity
l1=-5.73e-2; % freezing point salinity coefficient
l2=8.32e-2; % freezing point offset
l3=7.61e-4; % freezing point-depth coefficient
L=3.35e5; % latent heat of fusion for ice

% fitting parameters
if isfield(UserVar,'E0')
    E0=UserVar.E0;
else
    E0=0.036; % entrainment factor default value from Lazeroms et al 2019
end
if isfield(UserVar,'GammaTS')
    GammTS=UserVar.GammaTS;
else
    GammTS=5.9e-4; % exchange coefficient default value from Lazeroms et al 2019
end

%% obtain plume geometry
[N_m,~,ind_meltnodes,grad_mean,grad_local,zGL]=CalcPlumeGeometry(MUA,F,CtrlVar,1);
a_local = atan(grad_local); % local slope
a_mean = atan(grad_mean); % mean slope
z = N_m(:,3); % depth of melt node

% find melt nodes that are at the same depth as GL nodes
ind_samedepth = find(isapprox(z,zGL));

%% define a simple piecewise-linear temperature profile and constant salinity that is 
%% representative for AS ice shelves and apply everywhere
% THIS NEEDS TO BE REPLACED BY ICE SHELF-DEPENDENT PROFILES

% typical ocean properties in AS, values taken from Fig3 in De Rydt et al 2013
zT = -700; % thermocline depth
Ts = -1; % surface temperature
Td = 1.2; % typical mCDW temperature in AS
Sd = 34.7; % typical mCDW salinity in AS

T = Ts+z*(Td-Ts)/zT; T(z<zT)=Td; % piecewise linear temperature profile similar to Rosier et al (2024)
Tgl = Ts+zGL*(Td-Ts)/zT; Tgl(zGL<zT)=Td; % temperature at GL
% piecewise linear averaging to obtain mean temperature:
T_belowTz = max(zT-zGL,0)*Td; % this is equal to zT-zGL for zT>zGL, or 0 for zGL>zT
T_aboveTz = (z-max(zGL,zT)).*(T+Td*(zGL<zT)+Tgl.*(zGL>=zT))/2;
Tm = (T_belowTz+T_aboveTz)./(z-zGL);
Tm(ind_samedepth) = T(ind_samedepth); % remove infinite values

% freezing temperatures
Tf = l1*Sd + l2 + l3*z; % local freezing temperature at ice base
Tglf = l1*Sd + l2 + l3*zGL; % freezing temperature at GL

%% now assemble melt rates - following Rosier et al (2024)
% scaled distance
x = l3*(z-zGL)./(Tm-Tglf).*(1+0.6*(E0*sin(a_mean)./(GammTS*(1-l1*(Sd*c/L))+E0*sin(a_mean))).^(3/4)).^(-1);

% dimensionless velocity and thermal forcing functions
fU = 1/sqrt(2)*(1-x).^(1/3).*(1-(1-x).^(4/3)).^(1/2);
fdT = 1/2*(3*(1-x)-1./(1-x).^(1/3));

% velocity scale 
U = (betaT*g*E0*sin(a_local)./(l3*(Cd+E0*sin(a_local)))).^(1/2).*...
    ((GammTS*(betaS/betaT*(Sd*c/L)-1))./(GammTS*(1-l1*(Sd*c/L))+E0*sin(a_local))).^(1/2).*(T-Tf); % we use local values here

% thermal forcing scale
dT = (E0.*sin(a_mean)./(GammTS*(1-l1*(Sd*c/L))+E0*sin(a_mean))).*(Tm-Tglf); % we use mean values here

% melt rate
ab_tmp = (GammTS/(L/c)).*U.*fU.*dT.*fdT;

% project back onto Ua mesh
ab = zeros(MUA.Nnodes,1)+nan;
ab(ind_meltnodes) = -ab_tmp*365.25*24*60*60;

%% plotting
if doplots
    figure; hold on;
    PlotMeshScalarVariable(CtrlVar,MUA,ab);
    CtrlVar.PlotGLs=1;
    PlotGroundingLines(CtrlVar,MUA,F.GF); 
    title("Rosier et al 2024");
end

end
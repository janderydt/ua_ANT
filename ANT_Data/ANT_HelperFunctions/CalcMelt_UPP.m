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
    GammaTS=UserVar.GammaTS;
else
    GammaTS=0.0118; % exchange coefficient default value from Lazeroms et al 2019
end
Stanton = sqrt(Cd)*GammaTS; % effective thermal stanton number

% plume geometry
[N_m,~,ind_meltnodes,a_mean,a_local,zGL]=CalcPlumeGeometry(MUA,F,CtrlVar,0);
z = N_m(:,3);

% melt nodes that are at the same depth as GL nodes
ind_samedepth = find(isapprox(z,zGL));

% define simple piecewise-linear profile that is representative for AS and apply everywhere
% THIS NEEDS TO BE REPLACED BY ICE SHELF-DEPENDENT PROFILES

% typical ocean properties in AS, values taken from Fig3 in De Rydt et al 2013
zT = -800; % thermocline depth
Ts = -1; % surface temperature
Td = 1.2; % typical mCDW temperature in AS
Ss = 33.9; % surface temperature
Sd = 34.7; % typical mCDW salinity in AS

% GL properties
Tgl = Ts+zGL*(Td-Ts)/zT; Tgl(zGL<zT)=Td; % temperature at GL
Sgl = Ss+zGL*(Sd-Ss)/zT; Sgl(zGL<zT)=Sd; % salinity at GL
Tfgl = l1*Sgl + l2 + l3*zGL; % freezing temperature at GL

Ta = Ts+z*(Td-Ts)/zT; Ta(z<zT)=Td; % piecewise linear temperature profile similar to Rosier et al (2024)
% piecewise linear integration to obtain mean temperature:
T_belowTz = max(zT-zGL,0)*Td;
T_aboveTz = (z-max(zGL,zT)).*(Ta+Td*(zGL<zT)+Tgl.*(zGL>=zT))/2;
Tmean = (T_belowTz+T_aboveTz)./(z-zGL);
Tmean(ind_samedepth) = Ta(ind_samedepth); % remove inf

Sa =  Ss+z*(Sd-Ss)/zT; Sa(z<zT)=Sd; % piecewise linear salinity profile similar to Rosier et al (2024)
% piecewise linear integration to obtain mean salinity:
S_belowTz = max(zT-zGL,0)*Sd;
S_aboveTz = (z-max(zGL,zT)).*(Sa+Sd*(zGL<zT)+Sgl.*(zGL>=zT))/2;
Smean = (S_belowTz+S_aboveTz)./(z-zGL);
Smean(ind_samedepth) = Sa(ind_samedepth); % remove inf

% now assemble melt rates - following Rosier et al (2024)

% scaled distance
x = l3*(z-zGL)./(Tmean-Tfgl).*(1+0.6*(E0*sin(a_mean)./(Stanton*(1-l1*(Smean*c/L))+E0*sin(a_mean))).^(3/4)).^(-1);

% dimensionless velocity and thermal forcing functions
fU = 1/sqrt(2)*(1-x).^(1/3).*(1-(1-x).^(4/3)).^(1/2);
fdT = 1/2*(3*(1-x)-1./(1-x).^(1/3));

% velocity scale 
U = (betaT*g*E0*sin(a_local)./(l3*(Cd+E0*sin(a_local)))).^(1/2).*...
    ((Stanton*(betaS/betaT*(Sa*c/L)-1))./(Stanton*(1-l1*(Sa*c/L))+E0*sin(a_local))).^(1/2).*(Ta-Tfgl);

% thermal forcing scale
dT = (E0.*sin(a_local)./(Stanton*(1-l1*(Smean*c/L))+E0*sin(a_local))).*(Tmean-Tfgl);

% melt rate
ab_tmp = (Stanton/(L/c)).*U.*fU.*dT.*fdT;
ab = zeros(MUA.Nnodes,1)+nan;
ab(ind_meltnodes) =  -ab_tmp*365.25*24*60*60;

% plotting
if doplots
    figure; hold on;
    PlotMeshScalarVariable(CtrlVar,MUA,ab);
    CtrlVar.PlotGLs=1;
    PlotGroundingLines(CtrlVar,MUA,F.GF);
end
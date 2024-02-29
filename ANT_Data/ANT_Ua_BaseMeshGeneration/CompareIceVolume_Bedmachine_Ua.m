function CompareIceVolume_Bedmachine_Ua

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

%% calculate total ice volume native Bedmachine dataset
ncfile = froot_data+"/BedMachine_Antarctica/BedMachineAntarctica-v3.nc";
xdata = double(ncread(ncfile,'x'));
resolution = xdata(2)-xdata(1);
thickness = double(ncread(ncfile,'thickness')');

IceVolume_Bedmachine = sum(thickness*resolution^2,'all');

%% calculate total ice volume of Bedmachine interpolated onto Ua grid
load("../ANT_Interpolants/GriddedInterpolants_Bedmachine_v3_unmodified.mat","Fs","Fb");
load("ANT_basemesh_2003");
%load("MeshAntISMIP6_2_SainanSun");
[p,t] = refine(MUA.coordinates,MUA.connectivity);
[p,t] = refine(p,t);
[p,t] = refine(p,t);
%[p,t] = smooth2(p,[],t);
MUA = CreateMUA(CtrlVar,t,p);

%load("Antarctica_ref_nEle168146.mat");
%load("MeshAntISMIP6_2_SainanSun");
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
h_Ua = Fs(x,y) - Fb(x,y);

Int = FEintegrate2D([],MUA,h_Ua);
IceVolume_Ua = sum(Int,'all');

RelDiff = (IceVolume_Ua-IceVolume_Bedmachine)/IceVolume_Bedmachine*100;

fprintf("Difference between Ua and Bedmachine ice volume = %s%%\n",string(RelDiff));
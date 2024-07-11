function test

load('/mnt/md0/Ua/cases/ANT/ANT_Inverse/cases/ANT_nsmbl_Inverse_3145/ANT_nsmbl_Inverse_3145_AGlen-Estimate.mat');
FAGlen = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),AGlen);
MUA_orig = MUA;

load('/mnt/md0/Ua/cases/ANT/ANT_Inverse/cases/ANT_nsmbl_Inverse_3145/ANT_nsmbl_Inverse_3145_C-Estimate.mat');
load('/mnt/md0/Ua/cases/ANT/ANT_Inverse/cases/ANT_nsmbl_Inverse_3145/ANT_nsmbl_Inverse_3145-RestartFile_InverseCycle1.mat','GF');
GF_orig = GF.node; C(GF_orig<0.99)=0;
FC =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),C);

load('/mnt/md0/Ua/cases/ANT/ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_2000_2009_2014_2018_meshmin3000_meshmax100000_refined_extrudemesh1_variableboundaryres1.mat')
MUA_base = MUA;
x_r=[min(MUA_base.coordinates(:,1)):5e3:max(MUA_base.coordinates(:,1))];
y_r=[min(MUA_base.coordinates(:,2)):5e3:max(MUA_base.coordinates(:,2))];
[X_r,Y_r]=ndgrid(x_r,y_r);

AGlen_r = FAGlen(X_r,Y_r);
C_r = FC(X_r,Y_r);

Ind = find(~inpoly2([X_r(:) Y_r(:)],[MUA_orig.Boundary.x(:) MUA_orig.Boundary.y(:)]));
AGlen_r(Ind) = 0;
C_r(Ind) = 0;

FAGlen_r = griddedInterpolant(X_r,Y_r,AGlen_r);
FC_r = griddedInterpolant(X_r,Y_r,C_r);

Velinterpolantfile = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities";
Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2000";

Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,FAGlen_r,0,"-scalar-AGlen-");
Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,FC_r,0,"-scalar-C-");




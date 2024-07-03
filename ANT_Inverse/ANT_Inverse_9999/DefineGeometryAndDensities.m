function [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

if nargin<5
    FieldsToBeDefined='sbSBrho';
end

s=[]; b=[]; S=[]; B=[]; rho = [];
alpha=0 ;

% to speed up the geometry assembly, we check if the geometry fields have
% previously been constructed for the same mesh. If not, we read the
% interpolants and construct a new file for this particular mesh.
filename_geometryfields = UserVar.GeometryInterpolants;
filename_geometryfields = erase(filename_geometryfields,["GriddedInterpolants_","ScatteredInterpolants_",".mat"]);
filename_geometryfields = filename_geometryfields + "_mesh_Nodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

fprintf('Loading geometry and density fields %s...',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

if exist(filename_geometryfields,"file")
    load(filename_geometryfields,"B","b","S","s","rho");
else
    load(UserVar.GeometryInterpolants,'FB');
    B = FB(x,y);
    clearvars FB;
    B = inpaint_nans(B,4);
    load(UserVar.GeometryInterpolants,'Fs');
    s = Fs(x,y);
    clearvars Fs; 
    s = inpaint_nans(s,4);
    load(UserVar.GeometryInterpolants,'Fb'); 
    b = Fb(x,y);
    clearvars Fb;
    b = inpaint_nans(b,4);
    S = 0*x;
    load(UserVar.DensityInterpolant,'Frho');
    rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
    clearvars Frho; 
    rho(rho<100)=100;
    rho(rho>917)=917;
    save(filename_geometryfields,"B","b","S","s","rho");
end

rhow=1027; 
            
g=9.81/1000;

fprintf('done.\n');
   


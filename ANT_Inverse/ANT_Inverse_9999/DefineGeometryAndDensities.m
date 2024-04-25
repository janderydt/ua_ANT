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
filename_geometryfields = erase(filename_geometryfields,"GriddedInterpolants_");
filename_geometryfields = erase(filename_geometryfields,".mat");
filename_geometryfields = filename_geometryfields + "_mesh_Nodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

fprintf('Loading geometry and density fields %s...',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

if exist(filename_geometryfields,"file")
    load(filename_geometryfields,"B","b","S","s","rho");
else
    load(UserVar.GeometryInterpolants,'FB','Fs','Fb');
    B = FB(x,y);
    B = inpaint_nans(B,4);
    s = Fs(x,y);
    s = inpaint_nans(s,4);
    b = Fb(x,y);
    b = inpaint_nans(b,4);
    S = 0*x;
    load(UserVar.DesnityInterpolant,'Frho');
    rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
    rho(rho<100)=100;
    rho(rho>917)=917;
    save(filename_geometryfields,"B","b","S","s","rho");
end

rhow=1027; 
            
g=9.81/1000;

fprintf('done.\n');
   


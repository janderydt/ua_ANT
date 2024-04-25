function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

if nargin<5
    FieldsToBeDefined='sbSB';
end

s=[]; b=[]; S=[]; B=[];
alpha=0 ;

% to speed up the geometry assembly, we check if the geometry fields have
% previously been constructed for the same mesh. If not, we read the
% interpolants and construct a new file for this particular mesh.
filename_geometryfields = UserVar.GeometryInterpolants;
filename_geometryfields = erase(filename_geometryfields,"GriddedInterpolants_");
filename_geometryfields = erase(filename_geometryfields,".mat");
filename_geometryfields = filename_geometryfields + "_mesh_Nodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

fprintf('Loading geometry fields %s...',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

if contains(FieldsToBeDefined,'S')
    S = 0*x;
end

if contains(FieldsToBeDefined,'B')
    if exist(filename_geometryfields,"file")
        load(filename_geometryfields,"B")
    else
        load(UserVar.GeometryInterpolants,'FB');
        B = FB(x,y);
        B = inpaint_nans(B,4);
        save(filename_geometryfields,"B","-append");
    end
end

if contains(FieldsToBeDefined,'s')
    if exist(filename_geometryfields,"file")
        load(filename_geometryfields,"s")
    else
        load(UserVar.GeometryInterpolants,'Fs');
        s = Fs(x,y);
        s = inpaint_nans(s,4);
        save(filename_geometryfields,"s","-append");
    end
end

if contains(FieldsToBeDefined,'b')
    if exist(filename_geometryfields,"file")
        load(filename_geometryfields,"b")
    else
        load(UserVar.GeometryInterpolants,'Fb');
        b = Fb(x,y);
        b = inpaint_nans(b,4);
        save(filename_geometryfields,"b","-append");
    end
end

fprintf('done.\n');
   


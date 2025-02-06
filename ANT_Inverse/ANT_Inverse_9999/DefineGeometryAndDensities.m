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
filename_geometryfields = filename_geometryfields + "_mesh_Nnodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

fprintf('Loading geometry and density fields %s.\n',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
S = 0*x;

fprintf('Trying to read fields from %s...',filename_geometryfields);

if exist(filename_geometryfields,"file")
    if contains(FieldsToBeDefined,"B")
        load(filename_geometryfields,"B");
    end
    if contains(FieldsToBeDefined,"b")
        load(filename_geometryfields,"b");
    end
    if contains(FieldsToBeDefined,"s")
        load(filename_geometryfields,"s");
    end
    if contains(FieldsToBeDefined,"rho")
        load(filename_geometryfields,"rho");
    end
    fprintf('done.\n');
else
    fprintf('file does not exist, try interpolants instead...');
    if exist(UserVar.GeometryInterpolants)
        if contains(FieldsToBeDefined,"B")
            load(UserVar.GeometryInterpolants,'FB');
            B = FB(x,y);
            clearvars FB;
            B = inpaint_nans(B,4);
        end
        if contains(FieldsToBeDefined,"s")
            load(UserVar.GeometryInterpolants,'Fs');
            s = Fs(x,y);
            clearvars Fs; 
            s = inpaint_nans(s,4);
        end
        if contains(FieldsToBeDefined,"b")
            load(UserVar.GeometryInterpolants,'Fb'); 
            b = Fb(x,y);
            clearvars Fb;
            b = inpaint_nans(b,4);
        end
        if contains(FieldsToBeDefined,"rho")
            load(UserVar.DensityInterpolant,'Frho');
            rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
            clearvars Frho;
        end    
        save(filename_geometryfields,"B","b","S","s","rho");
        fprintf('done.\n');
        fprintf('Used geometry interpolants from %s.\n',UserVar.GeometryInterpolants);
    else
        error("File with interpolants from previous spinup ("+UserVar.GeometryInterpolants+") does not exist. "+...
            "Check if it is created by ANT_GetUserVar_Inverse. Breaking out.");
    end
end

rho(rho<100)=100;
rho(rho>917)=917;
rhow=1027; 
            
g=9.81/1000;

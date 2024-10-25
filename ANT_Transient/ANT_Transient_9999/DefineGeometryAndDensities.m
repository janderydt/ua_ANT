function [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

persistent Ftmp MUAtmp

if nargin<5
    FieldsToBeDefined='sbSBrho';
end

s=[]; b=[]; B=[]; rho=[];
g=9.81/1000;
rhow = 1027;

% If remeshing is required, we read the bed geometry from the restart file
% This is to speed up the geometry assembly, and minimize the memory
% requirements on ARCHER2
%% THIS METHOD IGNORES AVAILABLE BED DATA AND COULD BE IMPROVED

fprintf('Loading geometry and density fields %s.\n',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
S = 0*x;

if exist(UserVar.InverseRestartFile,"file")

    if contains(FieldsToBeDefined,["B","s","b","rho"]) & isempty(Ftmp)
        tmp = load(UserVar.InverseRestartFile,"F","MUA");
        Ftmp = tmp.F;
        MUAtmp = tmp.MUA;
    end

    if contains(FieldsToBeDefined,"B")
        B = Ftmp.B;
        if numel(B) ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,B] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],MUAtmp,MUA,0,B);
        end
    end

    if contains(FieldsToBeDefined,"b")   
        b = Ftmp.b;
        if numel(b) ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,b] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],MUAtmp,MUA,0,b);
        end
    end

    if contains(FieldsToBeDefined,"s")
        s = Ftmp.s;
        if numel(s) ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,s] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],MUAtmp,MUA,0,s);
        end
    end

    rho = Ftmp.rho;
    if numel(rho) ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,rho] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,917,rho);
    end
    
    fprintf('Using fields from %s. Done.\n',UserVar.InverseRestartFile);

else

    fprintf('%s does not exist, try interpolants instead...',UserVar.InverseRestartFile);

    if exist(UserVar.GeometryInterpolants,"file")

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

    else

        error("Neither "+UserVar.InverseRestartFile+" nor "+UserVar.GeometryInterpolants+" exist"+...
            " so I don't know what geometry fields to use. Breaking out.");

    end

    if exist(UserVar.DensityInterpolants,"file")
	    load(UserVar.DensityInterpolant,'Frho');
        rho = Frho(x,y);
        clearvars Frho;
    else

        error("Neither "+UserVar.InverseRestartFile+" nor "+UserVar.DensityInterpolant+" exist"+...
            " so I don't know what ice density to use. Breaking out.");
    end

    %save(filename_geometryfields,"B","b","S","s","rho");
    fprintf('done.\n');
    fprintf('Used geometry interpolants from %s for grounded ice.\n',UserVar.GeometryInterpolants);

end

rho(rho<100)=100;
rho(rho>917)=917;
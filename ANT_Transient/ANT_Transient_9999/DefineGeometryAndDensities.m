function [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

if nargin<5
    FieldsToBeDefined='sbSBrho';
end

s=[]; b=[]; S=[]; B=[]; rho=[];
g=9.81/1000;
rhow = 1027;

% If remeshing is required, we read the bed geometry from the restart file
% This is to speed up the geometry assembly, and minimize the memory requirements
%% THIS NEEDS TO BE IMPROVED

fprintf('Loading geometry and density fields %s.\n',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
S = 0*x;

if exist(UserVar.RestartFile,"file")
    tmp=load(FieldsToBeDefined,"MUA");
    if contains(FieldsToBeDefined,"B")
        load(UserVar.RestartFile,"B");
        if tmp.MUA.Nnodes ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,B] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,0,B);
        end
    end
    if contains(FieldsToBeDefined,"b")   
        load(UserVar.RestartFile,"b");
        if tmp.MUA.Nnodes ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,b] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,0,b);
        end
    end
    if contains(FieldsToBeDefined,"s")
        load(UserVar.RestartFile,"s");
        if tmp.MUA.Nnodes ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,s] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,0,s);
        end
    end
    
    
    load(UserVar.RestartFile,"rho");
    if tmp.MUA.Nnodes ~= MUA.Nnodes % Nnodes in restartfile does not equal number of nodes in new mesh. Interpolate.
            CtrlVar.MapOldToNew.method = "ShapeAndScattered"; 
            [~,rho] = MapNodalVariablesFromMesh1ToMesh2(CtrlVar,[],tmp.MUA,MUA,917,rho);
    end
    
    fprintf('Using fields from %s. Done.\n',UserVar.RestartFile);

else

    fprintf('%s does not exist, try interpolants instead...',UserVar.RestartFile);
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
        
	load(UserVar.DensityInterpolant,'Frho');
        rho = Frho(x,y);
        clearvars Frho;
        
        %save(filename_geometryfields,"B","b","S","s","rho");
        fprintf('done.\n');
        fprintf('Used geometry interpolants from %s for grounded ice.\n',UserVar.GIGeometryInterpolants);
    else
        error("File with interpolants ("+UserVar.GeometryInterpolants+") does not exist. Breaking out.");
    end
end

rho(rho<100)=100;
rho(rho>917)=917;

%% where is the grounding line?
h = s-b;
[b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow);

%% find floating nodes that are not lakes
[GF,~,~,~]=IceSheetIceShelves(CtrlVar,MUA,GF);
[LakeNodes,~,~,~] = LakeOrOcean3(CtrlVar,MUA,GF);
ISnodes = find(GF.node<0.5 & GF.NodesDownstreamOfGroundingLines & ~LakeNodes);

function [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

%if nargin<5
    FieldsToBeDefined='sbSBrho';
%end

s=[]; b=[]; S=[]; B=[]; rho = [];
alpha=0 ;

% to speed up the geometry assembly, we check if the geometry fields have
% previously been constructed for the same mesh. If not, we read the
% interpolants and construct a new file for this particular mesh.
filename_ISgeometryfields = UserVar.ISGeometryInterpolants;
filename_ISgeometryfields = erase(filename_ISgeometryfields,["GriddedInterpolants_","ScatteredInterpolants_",".mat"]);
filename_ISgeometryfields = filename_ISgeometryfields + "_mesh_Nnodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

filename_GIgeometryfields = UserVar.GIGeometryInterpolants;
filename_GIgeometryfields = erase(filename_GIgeometryfields,["GriddedInterpolants_","ScatteredInterpolants_",".mat"]);
filename_GIgeometryfields = filename_GIgeometryfields + "_mesh_Nnodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

fprintf('Loading geometry and density fields %s.\n',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
S = 0*x;
rhow = 1027;

%% first deal with the grounded ice and densities
fprintf('Grounded ice: trying to read fields from %s...',filename_GIgeometryfields);

if exist(filename_GIgeometryfields,"file")
    if contains(FieldsToBeDefined,"B")
        load(filename_GIgeometryfields,"B");
    end
    if contains(FieldsToBeDefined,"b")
        load(filename_GIgeometryfields,"b");
    end
    if contains(FieldsToBeDefined,"s")
        load(filename_GIgeometryfields,"s");
    end
    
    load(filename_GIgeometryfields,"rho");
    fprintf('done.\n');

else

    fprintf('file does not exist, try interpolants instead...');
    if exist(UserVar.GIGeometryInterpolants)
        if contains(FieldsToBeDefined,"B")
            load(UserVar.GIGeometryInterpolants,'FB');
            B = FB(x,y);
            clearvars FB;
            B = inpaint_nans(B,4);
        end
        if contains(FieldsToBeDefined,"s")
            load(UserVar.GIGeometryInterpolants,'Fs');
            s = Fs(x,y);
            clearvars Fs; 
            s = inpaint_nans(s,4);
        end
        if contains(FieldsToBeDefined,"b")
            load(UserVar.GIGeometryInterpolants,'Fb'); 
            b = Fb(x,y);
            clearvars Fb;
            b = inpaint_nans(b,4);
        end
        
	    load(UserVar.DensityInterpolant,'Frho');
        rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
        clearvars Frho;
        
        %save(filename_GIgeometryfields,"B","b","S","s","rho");
        fprintf('done.\n');
        fprintf('Used geometry interpolants from %s for grounded ice.\n',UserVar.GIGeometryInterpolants);
    else
        error("File with interpolants ("+UserVar.GIGeometryInterpolants+") does not exist. Breaking out.");
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

%% now deal with floating ice
fprintf('Floating ice: trying to read fields from %s...',filename_ISgeometryfields);

if exist(filename_ISgeometryfields,"file")
    if contains(FieldsToBeDefined,"b")
        tmp=load(filename_ISgeometryfields,"b");
        b(ISnodes) = tmp.b(ISnodes);
    end
    if contains(FieldsToBeDefined,"s")
        tmp=load(filename_ISgeometryfields,"s");
        s(ISnodes) = tmp.s(ISnodes);
    end
    fprintf('done.\n');
else
    fprintf('file does not exist, try interpolants instead...');
    if exist(UserVar.ISGeometryInterpolants)
        if contains(FieldsToBeDefined,"s")
            load(UserVar.ISGeometryInterpolants,'Fs');
            s_IS = Fs(x(ISnodes),y(ISnodes));
            clearvars Fs; 
            s_IS = inpaint_nans(s_IS,4);
            s(ISnodes) = s_IS;
        end
        if contains(FieldsToBeDefined,"b")
            load(UserVar.ISGeometryInterpolants,'Fb'); 
            b_IS = Fb(x(ISnodes),y(ISnodes));
            clearvars Fb;
            b_IS = inpaint_nans(b_IS,4);
            b(ISnodes) = b_IS;
        end
        fprintf('done.\n');
        fprintf('Used geometry interpolants from %s for floating ice.\n',UserVar.ISGeometryInterpolants);
    else
        error("File with interpolants ("+UserVar.ISGeometryInterpolants+") does not exist. Breaking out.");
    end
end

g=9.81/1000;


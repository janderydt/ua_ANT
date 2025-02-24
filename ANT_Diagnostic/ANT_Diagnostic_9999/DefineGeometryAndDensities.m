function [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

%if nargin<5
    FieldsToBeDefined='sbSBrho';
%end

s=[]; b=[]; S=[]; B=[]; rho=[];
alpha=0 ;

fprintf('Loading geometry and density fields %s.\n',FieldsToBeDefined);

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
S = 0*x;
rhow = 1027;
g=9.81/1000;

%% first deal with the grounded ice and densities
if contains(UserVar.GIBaseGeometry,["GriddedInterpolants_","ScatteredInterpolants_"])
    % this will always be the case when the geometry from the first inverse
    % cycle is used
    
    % to speed up the geometry assembly, we check if the geometry fields have
    % previously been constructed for the same mesh. If not, we read the
    % interpolants and construct a new file for this particular mesh.
    filename_GIgeometryfields = UserVar.GIBaseGeometry;
    filename_GIgeometryfields = erase(filename_GIgeometryfields,["GriddedInterpolants_","ScatteredInterpolants_",".mat"]);
    filename_GIgeometryfields = filename_GIgeometryfields + "_mesh_Nnodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

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
        if exist(UserVar.GIBaseGeometry)
            if contains(FieldsToBeDefined,"B")
                load(UserVar.GIBaseGeometry,'FB');
                B = FB(x,y);
                clearvars FB;
                B = inpaint_nans(B,4);
            end
            if contains(FieldsToBeDefined,"s")
                load(UserVar.GIBaseGeometry,'Fs');
                s = Fs(x,y);
                clearvars Fs; 
                s = inpaint_nans(s,4);
            end
            if contains(FieldsToBeDefined,"b")
                load(UserVar.GIBaseGeometry,'Fb'); 
                b = Fb(x,y);
                clearvars Fb;
                b = inpaint_nans(b,4);
            end
            
            load(UserVar.DensityInterpolant,'Frho');
            rho = Frho(x,y);
            clearvars Frho;
            
            save(filename_GIgeometryfields,"B","b","S","s","rho");
            fprintf('done.\n');
            fprintf('Used geometry interpolants from %s for grounded ice.\n',UserVar.GIBaseGeometry);
        else
            error("File with interpolants ("+UserVar.GIBaseGeometry+") does not exist. Breaking out.");
        end
    end

    h = s-b;

else
    % this will be the case when the geometry from a relaxation run is used; 
    % in this case we want to construct the updated geometry by extrapolting 
    % the relaxed geometry to potentially new cells.
    
    fprintf('Grounded ice: reading fields from %s...',UserVar.GIBaseGeometry);

    tmp = load(UserVar.GIBaseGeometry);

    if contains(FieldsToBeDefined,"s")
        Fs = scatteredInterpolant(tmp.MUA.coordinates(:,1),tmp.MUA.coordinates(:,2),tmp.F.s,'natural');
        s = Fs(x,y);
        clearvars Fs;
    end

    if contains(FieldsToBeDefined,"b")
        Fb = scatteredInterpolant(tmp.MUA.coordinates(:,1),tmp.MUA.coordinates(:,2),tmp.F.b,'natural');
        b = Fb(x,y);
        clearvars Fb;
    end

    clearvars tmp

    if contains(FieldsToBeDefined,"rho")
        load(UserVar.DensityInterpolant,'Frho');
        rho = Frho(x,y);
        clearvars Frho;
        rho = inpaint_nans(rho,4);
    end

    if contains(FieldsToBeDefined,"B")
        load(UserVar.DensityInterpolant,'FB');
        B = FB(x,y);
        clearvars FB;
        B = inpaint_nans(B,4);
    end

    fprintf('done.\n');
   
    % now compute correction to ice surface and base
    if UserVar.AdaptGIBaseGeometry

        fprintf('Apply thickness corrections for grounded and floating ice from %s and %s...',...
           UserVar.InterpolantsToCalculateNewGIGeometry(1),...
           UserVar.InterpolantsToCalculateNewGIGeometry(2));

        tmp1=load(UserVar.InterpolantsToCalculateNewGIGeometry(1),"Fs","Fb");
        h1 = tmp1.Fs(x,y)-tmp1.Fb(x,y);
        clearvars tmp1
        tmp2=load(UserVar.InterpolantsToCalculateNewGIGeometry(2),"Fs","Fb");
        h2 = tmp2.Fs(x,y)-tmp2.Fb(x,y);
        clearvars tmp2
        dh = h2-h1;

        fprintf('done. \n');

    else
        dh = 0*x;
    end

    h = s-b; 
    h = h+dh;
    
end

rho(rho<100)=100;
rho(rho>917 | ~isfinite(rho))=917;

%% where is the grounding line?
[b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow);

%% find floating nodes that are not lakes
[GF,~,~,~]=IceSheetIceShelves(CtrlVar,MUA,GF);
[LakeNodes,~,~,~] = LakeOrOcean3(CtrlVar,MUA,GF);
ISnodes = find(GF.node<0.5 & GF.NodesDownstreamOfGroundingLines & ~LakeNodes); 

%% now deal with floating ice
if UserVar.AdaptISBaseGeometry && ~UserVar.AdaptGIBaseGeometry
    
    if contains(UserVar.ISBaseGeometry,["GriddedInterpolants_","ScatteredInterpolants_"])
        % this will always be the case when the geometry from the first inverse
        % cycle is used
        
        % to speed up the geometry assembly, we check if the geometry fields have
        % previously been constructed for the same mesh. If not, we read the
        % interpolants and construct a new file for this particular mesh.
        filename_ISgeometryfields = UserVar.ISBaseGeometry;
        filename_ISgeometryfields = erase(filename_ISgeometryfields,["GriddedInterpolants_","ScatteredInterpolants_",".mat"]);
        filename_ISgeometryfields = filename_ISgeometryfields + "_mesh_Nnodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";
    
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
            if exist(UserVar.ISBaseGeometry)
                if contains(FieldsToBeDefined,"s")
                    load(UserVar.ISBaseGeometry,'Fs');
                    s_IS = Fs(x(ISnodes),y(ISnodes));
                    clearvars Fs; 
                    s_IS = inpaint_nans(s_IS,4);
                    s(ISnodes) = s_IS;
                end
                if contains(FieldsToBeDefined,"b")
                    load(UserVar.ISBaseGeometry,'Fb');
                    b_IS = Fb(x(ISnodes),y(ISnodes));
                    clearvars Fb; 
                    b_IS = inpaint_nans(b_IS,4);
                    b(ISnodes) = b_IS;
                end
                fprintf('done.\n');
                fprintf('Used geometry interpolants from %s for grounded ice.\n',UserVar.GIBaseGeometry);
            else
                error("File with interpolants ("+UserVar.ISBaseGeometry+") does not exist. Breaking out.");
            end
        end
    
    else
    
        % this will be the case when the geometry from a relaxation run is used; 
        % in this case we want to construct the updated geometry by extrapolting 
        % the relaxed geometry to potentially new cells.
    
        fprintf('Floating ice: reading fields from %s...',UserVar.ISBaseGeometry);
    
        tmp = load(UserVar.ISBaseGeometry);
    
        if contains(FieldsToBeDefined,["s","b"])
            Fs = scatteredInterpolant(tmp.MUA.coordinates(:,1),tmp.MUA.coordinates(:,2),tmp.F.s,'natural');
            s_IS = Fs(x(ISnodes),y(ISnodes));
            clearvars Fs;
            s_IS = inpaint_nans(s_IS,4);
        
            Fb = scatteredInterpolant(tmp.MUA.coordinates(:,1),tmp.MUA.coordinates(:,2),tmp.F.b,'natural');
            b_IS = Fb(x(ISnodes),y(ISnodes));
            clearvars Fb;
            b_IS = inpaint_nans(b_IS,4);
    
            clearvars tmp
        
            fprintf('done.\n');
        
            % now compute correction to ice surface
            fprintf('Apply thickness corrections for floating ice from %s and %s...',...
               UserVar.InterpolantsToCalculateNewISGeometry(1),...
               UserVar.InterpolantsToCalculateNewISGeometry(2));
        
            tmp1=load(UserVar.InterpolantsToCalculateNewISGeometry(1),"Fs","Fb");
            h1_IS = tmp1.Fs(x(ISnodes),y(ISnodes))-tmp1.Fb(x(ISnodes),y(ISnodes));
            clearvars tmp1
            tmp2=load(UserVar.InterpolantsToCalculateNewISGeometry(2),"Fs","Fb");
            h2_IS = tmp2.Fs(x(ISnodes),y(ISnodes))-tmp2.Fb(x(ISnodes),y(ISnodes));
            clearvars tmp2
            dh_IS = h2_IS-h1_IS;
        
            h_IS = s_IS-b_IS; 
            h_IS = h_IS+dh_IS;
            h(ISnodes) = h;
            [b,s,~,~]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow);
    
            fprintf('done.\n');

        end
        
    end

end
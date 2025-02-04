function ANT_RefineMesh

addpath("/mnt/md0/Matlab/Matlab_Functions/poly_stuff/");

basemesh = "ANT_basemesh_2000_2009_2014_2020_meshmin1500_meshmax100000_extrudemesh1_variableboundaryres1";
years = ["2000" "2009" "2014" "2020"];

% element size of base mesh
load(basemesh);
baseMUA = MUA;
EleArea=TriAreaFE(baseMUA.coordinates,baseMUA.connectivity);
M=Ele2Nodes(baseMUA.connectivity,baseMUA.Nnodes);
EleSizeCurrent=sqrt(M*EleArea);  % Elesizes around nodes

CtrlVar.PlotGLs=0; CtrlVar.GLsubdivide=1; CtrlVar.LineUpGLs=0;
CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';
CtrlVar.MeshSizeMax=100e3;
CtrlVar.MeshSize=100e3;
CtrlVar.MeshSizeMin=1.5e3;
CtrlVar.MeshAdapt.GLrange=[10e3 1 CtrlVar.MeshSizeMin*3; 5e3 10 CtrlVar.MeshSizeMin*2;...
    2e3 25 CtrlVar.MeshSizeMin];
CtrlVar.AdaptMeshMaxIterations=10;

% identify elements to be refined
EleSizeDesired = RefineGL(CtrlVar,EleSizeCurrent,baseMUA,years);

eRatio=EleSizeDesired./EleSizeCurrent;
eRatio=Nodes2EleMean(baseMUA.connectivity,eRatio);

% do not refine a greater number of elements than CtrlVar.MeshRefinementRatio*CurrentNumberOfElements
% at any given refinement step
test=sort(eRatio);
Ratio=0.9;
ElementsToBeRefined=eRatio<=test(ceil(numel(eRatio)*CtrlVar.LocalAdaptMeshRatio)) & eRatio<Ratio;

% have to make sure that if an element has just been refined that it will not
% then afterwards be a candidate for coarsening. If an element was refined, the
% size decreased by about a factor of 2 so if the ratio was R+eps it is now
% R/2+eps and I must set eRatio>2*R at the very least, for coarsening
ElementsToBeCoarsened=eRatio>=test(floor(numel(eRatio)*CtrlVar.LocalAdaptMeshRatio)) & eRatio>(2.1*Ratio);

kk=1; MUAold = baseMUA;
% iteratively continue to refine elements
while any(ElementsToBeRefined)
    MUAnew = LocalMeshRefinement(CtrlVar,[],MUAold,ElementsToBeRefined,ElementsToBeCoarsened);
    EleArea=TriAreaFE(MUAnew.coordinates,MUAnew.connectivity);
    M=Ele2Nodes(MUAnew.connectivity,MUAnew.Nnodes);
    EleSizeCurrent=sqrt(M*EleArea);  % Elesizes around nodes
    EleSizeDesired = RefineGL(CtrlVar,EleSizeCurrent,MUAnew,years);
    eRatio=EleSizeDesired./EleSizeCurrent;
    eRatio=Nodes2EleMean(MUAnew.connectivity,eRatio);
    test=sort(eRatio);
    Ratio=0.9;
    ElementsToBeRefined=eRatio<=test(ceil(numel(eRatio)*CtrlVar.LocalAdaptMeshRatio)) & eRatio<Ratio;
    ElementsToBeCoarsened=eRatio>=test(floor(numel(eRatio)*CtrlVar.LocalAdaptMeshRatio)) & eRatio>(2.1*Ratio);
    MUAold = MUAnew;
    fprintf("Remeshing iteration %s\n",string(kk));
    kk=kk+1;
end

MUA = MUAnew;

save(strrep(basemesh,"_extrude","_refined_extrude"),"MUA","CtrlVar");

BCfile = strrep(basemesh,"basemesh","meshboundarycoordinates");
copyfile(BCfile+".mat",strrep(BCfile,"_extrude","_refined_extrude")+".mat");


end



%% Helper function

function EleSizeDesired = RefineGL(CtrlVar,EleSizeDesired,MUAold,years)

persistent Interp

% part of this script is taken from NewDesiredEleSizesAndElementsToRefineOrCoarsen2.m
x = MUAold.coordinates(:,1); y = MUAold.coordinates(:,2);

if isempty(Interp)
    for ii=1:numel(years)
        geometryinterpolant = "../ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-"+years(ii)+"_EXTRUDED.mat";
        load(geometryinterpolant,"Fs","FB","Frho","Fb");
        Interp.("yr"+years(ii)).Fs=Fs;
        Interp.("yr"+years(ii)).FB=FB;
        Interp.("yr"+years(ii)).Frho=Frho;
        Interp.("yr"+years(ii)).Fb=Fb;
    end
    velocityinterpolant = "../ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED";
    load(velocityinterpolant,"Fus","Fvs");
    Interp.Fus = Fus;
    Interp.Fvs = Fvs;
      
end

for ii=1:numel(years)
    % initialize geometry
    F.s = Interp.("yr"+years(ii)).Fs(x,y); F.b = Interp.("yr"+years(ii)).Fb(x,y);
    F.S = 0*F.s; F.B = Interp.("yr"+years(ii)).FB(x,y); 
    F.rho = Interp.("yr"+years(ii)).Frho(x,y); F.rhow=1027;
    F.h=F.s-F.b;
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUAold,F.h,F.S,F.B,F.rho,F.rhow); 
    U = hypot(Interp.Fus(x,y),Interp.Fvs(x,y));
    
    cooA=[x y];
    KdTree=KDTreeSearcher(cooA) ;
    EleSizeIndicator =zeros(MUAold.Nnodes,1)+CtrlVar.MeshSizeMax ;
    
    [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAold,F.GF);  % no need to align GL.
    for I=1:size(CtrlVar.MeshAdapt.GLrange,1)   
        ds=CtrlVar.MeshAdapt.GLrange(I,1);
        du=CtrlVar.MeshAdapt.GLrange(I,2);
        dh=CtrlVar.MeshAdapt.GLrange(I,3);
        if dh<CtrlVar.MeshSizeMin
            if CtrlVar.InfoLevelAdaptiveMeshing>=1
                fprintf('---> Warning: CtrlVar.MeshAdapt.GLrange(%i,2)=%g<CtrlVar.MeshSizeMin=%g \n',I,dh,CtrlVar.MeshSizeMin)
                fprintf('              Setting CtrlVar.MeshAdapt.GLrange(%i,2)=%g \n',I,CtrlVar.MeshSizeMin)
            end
            dh=CtrlVar.MeshSizeMin;
        end
        if CtrlVar.InfoLevelAdaptiveMeshing>=10
            fprintf('Nodes within the distance of %g from the grounding line are given the target element size %g \n',ds,dh)
        end
        [IDds,~,~,KdTree]=FindAllNodesWithinGivenRangeFromGroundingLine(CtrlVar,MUAold,xGL,yGL,ds,KdTree);
        % only refine elements where velocity is above certain threshold.
        % We chose 10m/yr here, which is arbitrary, but it limits the number
        % of elements to be refined.
        IDdu = find(U>du);
        ID = intersect(IDds,IDdu);
        EleSizeIndicator(ID)=dh;

        EleSizeDesired=min(EleSizeDesired,EleSizeIndicator);
    end
end

end



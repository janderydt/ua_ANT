function ANT_RefineMesh

addpath("/mnt/md0/Matlab/Matlab_Functions/poly_stuff/");

basemesh = "ANT_basemesh_2000_2009_2014_2018_meshmin3000_meshmax100000_extrudemesh1_variableboundaryres1";

refinedmeshes = [basemesh+"_Refined2000" basemesh+"_Refined2009" basemesh+"_Refined2014" basemesh+"_Refined2018"];

% element size of base mesh
load(basemesh);
baseMUA = MUA;
baseEleSize = FEintegrate2D(CtrlVar,baseMUA,ones(baseMUA.Nele,1));

% identify elements to be refined
NewEleSize = calcNewEleSize(CtrlVar,baseMUA,refinedmeshes);
ElementsToBeRefined = find(NewEleSize<baseEleSize-1);
ElementsToBeCoarsened = find(NewEleSize>baseEleSize+1);

CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';
CtrlVar.InfoLevelAdaptiveMeshing = 10;
CtrlVar.MeshSizeMax=100e3;
CtrlVar.MeshSize=100e3;
CtrlVar.MeshSizeMin=1e3;

kk=1; MUAold = baseMUA;
% iteratively continue to refine elements
while ~isempty(ElementsToBeRefined)
    MUAnew = LocalMeshRefinement(CtrlVar,[],MUAold,ElementsToBeRefined,ElementsToBeCoarsened);
    EleSize = FEintegrate2D(CtrlVar,MUAnew,ones(MUAnew.Nele,1));
    NewEleSize = calcNewEleSize(CtrlVar,MUAnew,refinedmeshes);
    ElementsToBeRefined = find(NewEleSize<EleSize-1);
    ElementsToBeCoarsened = find(NewEleSize>EleSize+1);
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

function NewEleSize = calcNewEleSize(CtrlVar,MUAold,refinedmeshes)

NewEleSize = [];

for ii=1:numel(refinedmeshes)
    EleSizeold = FEintegrate2D(CtrlVar,MUAold,ones(MUAold.Nele,1));

    load(refinedmeshes(ii));
    EleSizerefined = FEintegrate2D(CtrlVar,MUA,ones(MUA.Nele,1));
    xEle = MUAold.xEle(:); yEle = MUAold.yEle(:);
    Inan = find(isnan(MUA.Boundary.x(:)));
    if ~isempty(Inan)
        xB = MUA.Boundary.x(1:Inan(1)-1);
        yB = MUA.Boundary.y(1:Inan(1)-1);
    else
        xB = MUA.Boundary.x;
        yB = MUA.Boundary.y;
    end
    Iele = find(inpoly([xEle yEle],[xB yB]));
    modifiedbaseEleSize = EleSizeold(:);
    FEleSize = scatteredInterpolant(MUA.xEle(:),MUA.yEle(:),EleSizerefined,'nearest');
    modifiedbaseEleSize(Iele) = FEleSize(xEle(Iele),yEle(Iele));
    NewEleSize = [NewEleSize modifiedbaseEleSize(:)];

end

NewEleSize = min(NewEleSize,[],2);

end










function AMUND_GenerateBaseMesh

addpath(getenv("froot_tools"));
addpath(getenv("froot_ua")+"cases/ANT");
addpath(getenv("froot_ua")+"cases/ANT/ANT_ProcessingTools");

basins_to_analyze = {'','F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbott
CtrlVar=Ua2D_DefaultParameters; 
ANT_BaseMesh = "2000_2009_2014_2020_meshmin1500_meshmax100000_refined";
years = ["2000" "2009" "2014" "2020"];
AMUND_BaseMesh = years+"_meshmin1500_meshmax100000_refined";

%% load basins
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
%B = RemoveIceRisesAndIslands(B);

UserVar.casefolder = pwd;
UserVar.datafolder = pwd+"/../";
UserVar.Experiment = "";

UserVar = ANT_DefineBaseMesh(UserVar,ANT_BaseMesh);

for yy=1:numel(years)
    UserVar.Geometry = double(years(yy));
    UserVar = ANT_ApplyMeshModifications(UserVar);
    tmp = load(UserVar.InitialMeshFileName);
    MUA = tmp.MUA; 
    % identify basin id of each MUA node
    [MUA.basins,~] = Define_Quantity_PerBasin(MUA.coordinates(:,1),MUA.coordinates(:,2),B,0);
    MUA_basinnames = erase({MUA.basins(:).name},'-');
    basinnodes_all = [];
    for bb=1:numel(basins_to_analyze) 
        basin = char(erase(basins_to_analyze{bb},'-'));
        [~,BasinInd] = ismember(basin,MUA_basinnames);
        basinnodes = MUA.basins(BasinInd).ind;
        basinnodes_all = [basinnodes_all; basinnodes];
    end
    ElementsToBeDeactivated=any(~ismember(MUA.connectivity,basinnodes_all),2);
    [MUA,MUA.k,MUA.l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated);
    Ind_nan = find(isnan(MUA.Boundary.x));

    % only keep largest fully connected domain: this throws away ice rises and
    % islands that are not connected to the basins of interest
    MUA.Boundary.x = MUA.Boundary.x(1:Ind_nan(1)-1);
    MUA.Boundary.y = MUA.Boundary.y(1:Ind_nan(1)-1);
    % deactivate any elements outside the new boundary
    basinnodes = find(~inpoly(MUA.coordinates,[MUA.Boundary.x MUA.Boundary.y]));
    ElementsToBeDeactivated=any(ismember(MUA.connectivity,basinnodes),2);
    [MUA,MUA.k,MUA.l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated);

    MeshBoundaryCoordinates = [MUA.Boundary.x(:) MUA.Boundary.y(:)];

    switch years(yy)
        case "2000"
            xOuter = MUA.Boundary.x(858:1064);
            yOuter = MUA.Boundary.y(858:1064);
        case "2009"
            xOuter = MUA.Boundary.x(766:971);
            yOuter = MUA.Boundary.y(766:971);
        case "2014"
            xOuter = MUA.Boundary.x(830:1036);
            yOuter = MUA.Boundary.y(830:1036);
        case "2020"
            xOuter = MUA.Boundary.x(627:833);
            yOuter = MUA.Boundary.y(627:833);
    end

    figure; PlotMuaMesh(CtrlVar,MUA,[],'k'); hold on;
    plot(MUA.Boundary.x,MUA.Boundary.y,'-m');
    plot(xOuter,yOuter,'-og');
    title(years(yy));

    save("AMUND_basemesh_"+AMUND_BaseMesh(yy)+"_extrudemesh1_variableboundaryres1.mat","MUA");
    save("AMUND_meshboundarycoordinates_"+AMUND_BaseMesh(yy)+"_extrudemesh1_variableboundaryres1.mat","MeshBoundaryCoordinates");
    save("AMUND_fixedboundarypoints_"+AMUND_BaseMesh(yy)+"_extrudemesh1_variableboundaryres1.mat","xOuter","yOuter");

end
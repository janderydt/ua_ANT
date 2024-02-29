function plotGroundingLines_mod_obs

runID = "1381";

Table="../ANT_Inverse/RunTable.csv";

addpath(getenv("froot_tools"));
froot_data = getenv("froot_data");

%% read table
if exist(Table,'file')
    RunTable=readtable(Table); 
    ind = find(RunTable{:,'ExpID'}==double(runID));
else    
    error("Runtable does not exist"); 
end

year = RunTable{ind,'Geometry'};
year = 2018

%% model GL
load("../ANT_Inverse/ANT_Inverse_"+runID+"/ANT_Inverse_"+runID+"-RestartFile.mat","CtrlVarInRestartFile","MUA","GF");

figure; hold on;
CtrlVar = CtrlVarInRestartFile;
CtrlVar.PlotXYscale = 1e3;
PlotGroundingLines(CtrlVar,MUA,GF);
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k');


%% observed GL
S_insar = shaperead(froot_data+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
%S_insar_mililloTW = shaperead(froot_data+"/PIG_Thwaites_DotsonCrosson/Milillo_GL_Thwaites/Thwaites_2016-2017.shp");
%S_insar_mililloDC = shaperead(froot_data+"/PIG_Thwaites_DotsonCrosson/Milillo_GL_CrossonDotson/Grounding_Lines_CSK_Pope_Smith_Kohler.shp");

for ii=1:length(S_insar)
    date=S_insar(ii).DATE1;
    year_insar = date(1:4);
    [x_insar,y_insar] = ll2psxy(S_insar(ii).Y,S_insar(ii).X,-71,0);
    %if str2double(year_insar)>year-10 && str2double(year_insar)<year+10
    ginsar=plot(x_insar/CtrlVar.PlotXYscale,y_insar/CtrlVar.PlotXYscale,"color",[255, 165, 0]/255,"linewidth",2.5); 
    %end
end

% for ii=1:length(S_insar_mililloTW)
%     plot(ax_fig(bb),S_insar_mililloTW(ii).X/1e3,S_insar_mililloTW(ii).Y/1e3,"color",[255, 165, 0]/255,"linewidth",2.5);  
% end
% 
% for ii=1:length(S_insar_mililloDC)
%     plot(ax_fig(bb),S_insar_mililloDC(ii).X/1e3,S_insar_mililloDC(ii).Y/1e3,"color",[255, 165, 0]/255,"linewidth",2.5);                
% end


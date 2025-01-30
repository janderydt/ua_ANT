function plot_GroundingLines_FluxGates

runID = "1034"; %(2000: 1795, 2009: 1317, 2014: 1034, 2018:1381)

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
year = 2014;

%% model GL
load("../ANT_Inverse/ANT_Inverse_"+runID+"/ANT_Inverse_"+runID+"-RestartFile.mat","CtrlVarInRestartFile","MUA","GF");

figure; hold on;
CtrlVar = CtrlVarInRestartFile;
CtrlVar.PlotXYscale = 1e3;
CtrlVar.PlotGLs=0;
[xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,GF);
glua=plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'-k','linewidth',2);
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-','color',[0.5 0.5 0.5]);

%% observed GL
S_insar = shaperead(froot_data+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
S_insar_mililloTW = shaperead(froot_data+"/AmundsenSea/Milillo_GL_Thwaites/Thwaites_2016-2017.shp");
S_insar_mililloDC = shaperead(froot_data+"/AmundsenSea/Milillo_GL_CrossonDotson/Grounding_Lines_CSK_Pope_Smith_Kohler.shp");

for ii=1:length(S_insar)
    date=S_insar(ii).DATE1;
    year_insar = date(1:4);
    [x_insar,y_insar] = ll2psxy(S_insar(ii).Y,S_insar(ii).X,-71,0);
    %if str2double(year_insar)>year-10 && str2double(year_insar)<year+10
    %ginsar=plot(x_insar/CtrlVar.PlotXYscale,y_insar/CtrlVar.PlotXYscale,"color",[255, 165, 0]/255,"linewidth",2.5); 
    %end
end
% for ii=1:length(S_insar_mililloTW)
%     plot(ax_fig(bb),S_insar_mililloTW(ii).X/1e3,S_insar_mililloTW(ii).Y/1e3,"color",[255, 165, 0]/255,"linewidth",2.5);  
% end
% 
% for ii=1:length(S_insar_mililloDC)
%     plot(ax_fig(bb),S_insar_mililloDC(ii).X/1e3,S_insar_mililloDC(ii).Y/1e3,"color",[255, 165, 0]/255,"linewidth",2.5);                
% end

%% Gardner et al 2018 Flux Gates GL0, FG1 and FG2
froot_gates = getenv("froot_data")+"/Measures/FluxGates/line_shp/";
GL0 = shaperead(froot_gates+"GL0.shp");
FG1 = shaperead(froot_gates+"FG1.shp");
FG2 = shaperead(froot_gates+"FG2.shp");

gl0=plot(GL0.X/CtrlVar.PlotXYscale,GL0.Y/CtrlVar.PlotXYscale,'--','linewidth',1.5);
fg1=plot(FG1.X/CtrlVar.PlotXYscale,FG1.Y/CtrlVar.PlotXYscale,':','linewidth',1.5);
fg2=plot(FG2.X/CtrlVar.PlotXYscale,FG2.Y/CtrlVar.PlotXYscale,':','linewidth',1.5);

grid on; box on; xlabel('psx [km]'); ylabel('psy [km]'); axis equal;

legend([glua gl0 fg1 fg2],["Ua","GL0 (Gardner et al., 2018)","FG1","FG2"]);


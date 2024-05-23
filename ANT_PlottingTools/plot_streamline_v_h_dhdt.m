function plot_streamline_v_h_dhdt

glacier = 'PIG';

RunID = "1223"; % Weertman
RunID = "1585"; % Cornford

froot = "/mnt/md0/Ua/cases/ANT/ANT_Inverse/";
files = ["AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_InverseCycle1.mat",...
    "AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_SpinupCycle1.mat",...
    "AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_InverseCycle2.mat",...
    "AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_SpinupCycle2.mat",...
    "AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_InverseCycle3.mat",...
    "AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_SpinupCycle3.mat",...
    "AS_PROPHET_Inverse_"+RunID+"/AS_PROPHET_Inverse_"+RunID+"-RestartFile_InverseCycle4.mat"];

switch glacier
    case 'PIG'
        xstart = -1579830;
        ystart = -193393;
end

CM = jet(4);

%% set up figure
H=fig('units','inches','width',120*12/72.27,'height',65*12/72.27,'fontsize',14,'font','Helvetica');

tlo_fig = tiledlayout(1,3,"TileSpacing","compact");
for i = 1:3
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

ylabels=["Speed [m/yr]","Ice Thickness [m]","Change in ice thickness [m/yr]"];


%% set up streamline data
load(froot+files(1));
x=MUA.coordinates(:,1); y=MUA.coordinates(:,2); xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
[X,Y]=meshgrid([xmin:1e3:xmax],[ymin:1e3:ymax]);
Fu=scatteredInterpolant(x,y,F.ub);
Fv=scatteredInterpolant(x,y,F.vb);
U=Fu(X,Y);
V=Fv(X,Y);
U(find(~inpoly2([X(:),Y(:)],[MUA.Boundary.x,MUA.Boundary.y])))=nan;
V(find(~inpoly2([X(:),Y(:)],[MUA.Boundary.x,MUA.Boundary.y])))=nan;
SL=streamline(X,Y,U,V,xstart,ystart);
SLx=SL.XData; SLy=SL.YData;
SLd=cumsum(hypot(SLx(2:end)-SLx(1:end-1),SLy(2:end)-SLy(1:end-1)));

%% load and plot data
for ii=1:numel(files)

    load(froot+files(ii));
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    Fv = scatteredInterpolant(x,y,hypot(F.ub,F.vb));
    Fh = scatteredInterpolant(x,y,F.h);
    MUA=UpdateMUA(CtrlVarInRestartFile,MUA);
    [~,dhdt]=dhdtExplicitSUPG(UserVarInRestartFile,CtrlVarInRestartFile,MUA,F,BCs);
    Fdhdt = scatteredInterpolant(x,y,dhdt);
    FGF = scatteredInterpolant(x,y,F.GF.node);

    V = Fv(SLx,SLy);
    H = Fh(SLx,SLy);
    dHdT = Fdhdt(SLx,SLy);
    GF = FGF(SLx,SLy);
    I = find(GF<0.9); SLd_GL = SLd(I(1));

    if contains(files(ii),'InverseCycle')
        style = '--';
        width = 2;
    else
        style = '-';
        width = 1;
    end

    g(ii) = plot(ax_fig(1),SLd/1e3,V(2:end),LineStyle=style,LineWidth=width,Color=CM(ceil(ii/2),:));
    plot(ax_fig(1),[SLd_GL,SLd_GL]/1e3,[1000 4000],'-k');
    if ii==1
        FVmeas = scatteredInterpolant(x,y,hypot(Meas.us,Meas.vs));
        VMeas=FVmeas(SLx,SLy);
        plot(ax_fig(1),SLd/1e3,VMeas(2:end),LineStyle=style,LineWidth=width,Color='k');
    end

    plot(ax_fig(2),SLd/1e3,H(2:end),LineStyle=style,LineWidth=width,Color=CM(ceil(ii/2),:));
    plot(ax_fig(2),[SLd_GL,SLd_GL]/1e3,[0 2000],'-k');

    plot(ax_fig(3),SLd/1e3,dHdT(2:end),LineStyle=style,LineWidth=width,Color=CM(ceil(ii/2),:));
    plot(ax_fig(3),[SLd_GL,SLd_GL]/1e3,[-1000 1000],'-k');

end

for ii=1:3
    xlim(ax_fig(ii),[0 150]);
    ylabel(ax_fig(ii),ylabels(ii));
    grid(ax_fig(ii),"on");
    box(ax_fig(ii),"on");
end
xlabel(tlo_fig,'Distance [km]');






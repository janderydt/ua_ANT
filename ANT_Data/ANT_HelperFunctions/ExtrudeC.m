function Cext_Ua = ExtrudeC(CtrlVar,MUA,F)

addpath("/mnt/md0/Ua/cases/ANT/ANT_Data/ANT_Interpolants");

xUa = MUA.coordinates(:,1);
yUa = MUA.coordinates(:,2);
C = F.C;
vx = F.ub;
vy = F.vb;

% replace C for floating nodes with nan
C(F.GF.node<0.9)=nan;

% create interpolants for gridding later
FC = scatteredInterpolant(xUa ,yUa ,C);
Fvx = FC; Fvx.Values = vx;
Fvy = FC; Fvy.Values = vy;

% ExtrudeField requires gridded coordinates
xmin = min(xUa); ymin = min(yUa);
xmax = max(xUa); ymax = max(yUa);
dx = 500; dy = 500;
x = [xmin:dx:xmax]; y = [ymin:dy:ymax];
[X,Y]=ndgrid(x,y);
Iout = find(~inpoly2([X(:) Y(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));


Cgrid = FC(X,Y);
Cgrid(Iout) = nan;
vxgrid = Fvx(X,Y);
%vxgrid(Iout) = nan;
vygrid = Fvy(X,Y);
%vygrid(Iout) = nan;

% Extrude C and fill nan regions
if any(contains({'ASE','AS_PROPHET'},CtrlVar.Experiment)
    step = 0.01;
    n = 25000;
else
    error('Unknown experiment name for extruding C');
end

Cext = ExtrudeField(x,y,vxgrid',vygrid',x,y,Cgrid',step,n);
Cext = Cext';
I = find(~isnan(Cext));

% Reproject onto Ua mesh
FCext = scatteredInterpolant(X(I),Y(I),Cext(I),'natural');

Cext_Ua = FCext(xUa,yUa);
Cext_Ua(F.GF.node>=0.9) = F.C(F.GF.node>=0.9);


% figure; PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,log10(Cext_Ua));
% hold on;
% PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF);
% 
% 
% figure; PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,Cext_Ua-F.C);
% hold on;
% PlotGroundingLines(CtrlVarInRestartFile,MUA,F.GF);



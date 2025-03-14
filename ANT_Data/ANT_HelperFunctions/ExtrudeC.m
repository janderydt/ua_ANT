function Cext_Ua = ExtrudeC(CtrlVar,MUA,F)

%addpath("/mnt/md0/Ua/cases/ANT/ANT_Data/ANT_Interpolants");

xUa = MUA.coordinates(:,1);
yUa = MUA.coordinates(:,2);
C = F.C;
vx = F.ub;
vy = F.vb;

% replace C for floating nodes with nan
C(F.GF.node<0.5)=nan; 

% ExtrudeField requires gridded coordinates
xmin = min(xUa); ymin = min(yUa);
xmax = max(xUa); ymax = max(yUa);
if contains(CtrlVar.Experiment,["ASE","AS_PROPHET","AMUND"])
    dx = 2e3; dy = 2e3;
elseif contains(CtrlVar.Experiment,"ANT")
    dx = 5e3; dy = 5e3;
else
    error('Unknown experiment name for extruding C');
end
x = [xmin:dx:xmax]; y = [ymin:dy:ymax];
[X,Y]=ndgrid(x,y);
Iout = find(~inpoly2([X(:) Y(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));

% create interpolants for gridding later
FC = scatteredInterpolant(xUa ,yUa ,C);
Cgrid = FC(X,Y);
Cgrid(Iout) = nan;

Fvx = FC; Fvx.Values = vx;
vxgrid = Fvx(X,Y);

clear Fvx;

Fvy = FC; Fvy.Values = vy;
vygrid = Fvy(X,Y);

clear Fvy FC;

% Extrude C and fill nan regions
step = 0.01;
n = 10000;

Cext = ExtrudeField(x,y,vxgrid',vygrid',x,y,Cgrid',step,n);
Cext = Cext';
I = find(~isnan(Cext));

% Reproject onto Ua mesh
FCext = scatteredInterpolant(X(I),Y(I),Cext(I),'natural');

Cext_Ua = FCext(xUa,yUa);
Cext_Ua(F.GF.node>=0.5) = F.C(F.GF.node>=0.5);
Cext_Ua(Cext_Ua<CtrlVar.Cmin)=CtrlVar.Cmin;
Cext_Ua(Cext_Ua>CtrlVar.Cmax)=CtrlVar.Cmax;

% figure; PlotMeshScalarVariable(CtrlVar,MUA,imag(Cext_Ua));%log10(Cext_Ua));
% hold on;
% CtrlVar.PlotGLs=1;
% PlotGroundingLines(CtrlVar,MUA,F.GF);

%figure; %PlotMeshScalarVariable(CtrlVar,MUA,Cext_Ua-F.C);
%hold on;
%scatter(MUA.coordinates(:,1),MUA.coordinates(:,2),10,Cext_Ua-F.C,'filled');
%PlotGroundingLines(CtrlVar,MUA,F.GF);



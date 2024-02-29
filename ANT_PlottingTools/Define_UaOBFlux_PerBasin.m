function B = Define_UaOBFlux_PerBasin(MUA,F,GF,B,CtrlVar)

% We define OB nodes as Boundary nodes that are not afloat. 
x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);
% grounded fraction of edges
OB_GF = GF.node(MUA.Boundary.Edges);
% find edges that are not fully floating:
IOB = find(sum(OB_GF,2)>0.5);
% calculate midpoint of those edges
xOB = [x(MUA.Boundary.Edges(IOB,1)),x(MUA.Boundary.Edges(IOB,2))];
yOB = [y(MUA.Boundary.Edges(IOB,1)),y(MUA.Boundary.Edges(IOB,2))];

Fub=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
Fvb=Fub; Frho=Fub; Fh=Fub;
Fvb.Values=F.vb;
Frho.Values=F.rho;
Fh.Values=F.h;

[xOB,yOB,qOB]=FluxAcrossBoundary(xOB,yOB,Fub,Fvb,Fh,Frho);

%% check proximity to basins
PQ = [xOB(:) yOB(:)]; dist = [];
for ii=1:numel(B.x)
    P = [B.x{ii} B.y{ii}];
    [~,dist(ii,:)] = dsearchn(P,PQ);
end
[~,I] = min(dist,[],1); 

for ii=1:numel(B.x)
    B.OBpoints{ii} = find(I==ii);
    B.xOB{ii} = xOB(B.OBpoints{ii});
    B.yOB{ii} = yOB(B.OBpoints{ii});
    B.qOB{ii} = qOB(B.OBpoints{ii});
end

function [N_m,N_GL,theta_global,theta_local]=CalcPlumeGeometry(MUA,F,CtrlVar)

%% OUTPUTS:
% N_m: n x 3 matrix, n rows are melt nodes, 3 columns are x,y,z coordinats of each melt node
% N_GL: m x 3 matrix, m rows are GL nodes, 3 columns are x,y,z coordinats of each GL node
% idx: n x 10 matrix, n rows are melt nodes, 10 columns are indices of corresponding GL nodes that form the plume origins
% theta_global: n x 10 matrix, n rows are melt nodes, 10 columns are global slopes of the plumes
% theta_local: n x 10 matrix, n rows are melt nodes, 10 columns are local slopes of the plumes

if nargin<3
    folder = "/mnt/md0/Ua/cases/ANT/ANT_Inverse/cases/ANT_nsmbl_Inverse_14423";
    file = "ANT_nsmbl_Inverse_14423-RestartFile.mat";
    load(folder+"/"+file);
    CtrlVar = CtrlVarInRestartFile;
end

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);
b = F.b;
Fb = scatteredInterpolant(x,y,b);

%% find GL nodes
[GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
% this slows the code down a lot, but is needed to remove the smallest ice rises
[x_GL,y_GL] = ArrangeGroundingLinePos(CtrlVar,GLgeo);

%% define plume source nodes // remove small ice rises
x_GL = [nan; x_GL(:); nan]; y_GL = [nan; y_GL(:); nan]; 
Inan = find(isnan(x_GL));
x_GL_main = []; y_GL_main = []; 
for ii=1:numel(Inan)-1
    x_seg = [x_GL(Inan(ii)+1:Inan(ii+1)-1)];
    y_seg = y_GL(Inan(ii)+1:Inan(ii+1)-1);
    L_seg = sum(hypot(x_seg(2:end)-x_seg(1:end-1),y_seg(2:end)-y_seg(1:end-1)));
    if L_seg>100e3 % arbitrary cut off: any continuous GL segment with length<100km is removed.
        x_GL_main = [x_GL_main; x_seg(:); nan];
        y_GL_main = [y_GL_main; y_seg(:); nan];
    end
end

N_GL = [x_GL_main(:) y_GL_main(:) Fb(x_GL_main(:),y_GL_main(:))];

%% find melt nodes
NodesDownstreamOfGroundingLines="Relaxed";%"Strickt"; % strict
[LakeNodes,OceanNodes] = LakeOrOcean3_JDR(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele,[],NodesDownstreamOfGroundingLines);
x_m = x; y_m = y; b_m = b;
I_toremove = find(LakeNodes | ~OceanNodes);
x(I_toremove) = [];
y(I_toremove) = [];
b_m(I_toremove) = [];
N_m = [x(:) y(:) b_m(:)];

%% offset z-coordinate of N_m with large negative number (following Rosier et al. 2024)
N_m_adj = N_m;
N_m_adj(:,3) = N_m_adj(:,3)-1e6;

%% multiply z-coordinate of N_GL by two (following Rosier et al. 2024)
N_GL_adj = N_GL;
N_GL_adj(:,3) = N_GL_adj(:,3)*2;

%% k nearest neighbour search with default Euclidean metric
k=10;
[idx, dist] = knnsearch(N_GL_adj,N_m_adj,'K',k);

%% calculate global & local slope
[dbdx,dbdy,xint,yint]=calcFEderivativesMUA(b,MUA,CtrlVar); % local slope at integration points
[dbdx,dbdy]=ProjectFintOntoNodes(MUA,dbdx,dbdy); % local slope at nodes
dbdx_m = dbdx; dbdy_m = dbdy;
dbdx_m(I_toremove)=[]; % retain only melt nodes
dbdy_m(I_toremove)=[];

theta_global = zeros(size(idx));
theta_local = zeros(size(idx));
for ii=1:size(N_m,1) 
    dx = N_m(ii,1)-N_GL(idx(ii,:),1);
    dy = N_m(ii,2)-N_GL(idx(ii,:),2);
    dz = N_m(ii,3)-N_GL(idx(ii,:),3);
    dl = hypot(dx,dy);
    theta_global(ii,:) = atan(dz./dl);
    theta_local(ii,:) = atan(hypot(dbdx_m(ii),dbdy_m(ii)));
end

%% PLOTTING

figure; hold on;

PlotMuaMesh(CtrlVarInRestartFile,MUA,[],'color',[0.8 0.8 0.8]); hold on;
g(1)=scatter(N_m(:,1),N_m(:,2),10,mean(theta_global,2),'filled'); % global slope
g(2)=plot(x_GL,y_GL,'-m'); % original GL
g(3)=plot(x_GL_main,y_GL_main,'xb','markersize',2); % GL pinning points removed - these are
%the origins of the plume

% XLim=[-1.695e+06 -1.524e+06];
% YLim=[-3.844e+05 -2.498e+05];
% xlim(XLim); ylim(YLim);
% 
% % plot some pathways
% for ii=1:50:size(idx,1)
%     for jj=1:k
%         plot([N_m(ii,1) N_GL(idx(ii,jj),1)],[N_m(ii,2) N_GL(idx(ii,jj),2)],'--k');
%     end
% end

%scatter(N_m(:,1),N_m(:,2),10,dz_mean,'filled');
%figure; scatter(N_m(:,1),N_m(:,2),10,dl_mean,'filled');

caxis([-0.025 0.025])
CM=othercolor('RdYlBu7');
colormap(CM);

title('Global slope');
legend(g(:),["Global slope","GL","Possible plume origins",]);


figure; hold on;

PlotMuaMesh(CtrlVarInRestartFile,MUA,[],'color',[0.8 0.8 0.8]); hold on;
g(1)=scatter(N_m(:,1),N_m(:,2),10,mean(theta_local,2),'filled'); % global slope
g(2)=plot(x_GL,y_GL,'-m'); % original GL
g(3)=plot(x_GL_main,y_GL_main,'xb','markersize',2); % GL pinning points removed - these are
%the origins of the plume

% XLim=[-1.695e+06 -1.524e+06];
% YLim=[-3.844e+05 -2.498e+05];
% xlim(XLim); ylim(YLim);
% 
% % plot some pathways
% for ii=1:50:size(idx,1)
%     for jj=1:k
%         plot([N_m(ii,1) N_GL(idx(ii,jj),1)],[N_m(ii,2) N_GL(idx(ii,jj),2)],'--k');
%     end
% end

%scatter(N_m(:,1),N_m(:,2),10,dz_mean,'filled');
%figure; scatter(N_m(:,1),N_m(:,2),10,dl_mean,'filled');

caxis([-0.025 0.025])
CM=othercolor('RdYlBu7');
colormap(CM);

title('Local slope');
legend(g(:),["Local slope","GL","Possible plume origins",]);
%figure; scatter(N_m(:,1),N_m(:,2),10,theta_local,'filled'); % local slope
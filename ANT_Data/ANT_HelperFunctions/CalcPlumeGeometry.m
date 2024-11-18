function [N_m,N_GL,ind_meltnodes,theta_mean,theta_local,z_GL]=CalcPlumeGeometry(MUA,F,CtrlVar,doplots)

%% OUTPUTS:
% N_m: n x 3 matrix, n melt nodes, 3 columns are x,y,z coordinates of each melt node
% N_GL: m x 3 matrix, m GL nodes, 3 columns are x,y,z coordinates of each GL node
% ind_meltnode: n x 3 matrix, nodal indices of the melt nodes
% theta_mean: n x 1 matrix, mean slope of the plume path for each melt node
% theta_local: n x 1 matrix, local slope of the ice draft for each melt node
% z_GL: n x 1 matrix, grounding line depth of plume source

if nargin<4
    folder = "/mnt/md0/Ua/cases/ANT/ANT_Inverse/cases/ANT_nsmbl_Inverse_14423";
    file = "ANT_nsmbl_Inverse_14423-RestartFile.mat";
    load(folder+"/"+file);
    CtrlVar = CtrlVarInRestartFile;
    doplots=0;
end

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);
b = F.b;
Fb = scatteredInterpolant(x,y,b);

%% find GL nodes
[GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
% this slows the code down a lot, but is needed to remove ice rises and
% pinning points
[x_GL,y_GL] = ArrangeGroundingLinePos(CtrlVar,GLgeo);

%% define plume source nodes // remove small ice rises
x_GL = [nan; x_GL(:); nan]; y_GL = [nan; y_GL(:); nan]; 
Inan = find(isnan(x_GL));
x_GL_main = []; y_GL_main = []; 
for ii=1:numel(Inan)-1
    x_seg = [x_GL(Inan(ii)+1:Inan(ii+1)-1)];
    y_seg = [y_GL(Inan(ii)+1:Inan(ii+1)-1)];
    % calc length of GL segment
    L_seg = sum(hypot(x_seg(2:end)-x_seg(1:end-1),y_seg(2:end)-y_seg(1:end-1)));
    % calc distance between start and end point of GL segment
    d_start_end = hypot(x_seg(end)-x_seg(1),y_seg(end)-y_seg(1));

    % now impose 2 criteria to identify segments that will be removed from
    % the possible plume source nodes:
    % 1) cut-off for length of continuous GL segment: segments with
    % length<100km are removed. This deals with the smallest pinning points.    
    % 2) remove GL segments that form an approximately closed loop. This 
    % removes most ice rises and larger pinning points, but retains some 
    % larger ones near the domain bounadry, sucha as Berkner Island or Bear
    % Peninsula
    if L_seg>100e3 && d_start_end>5e3
        x_GL_main = [x_GL_main; x_seg(:); nan];
        y_GL_main = [y_GL_main; y_seg(:); nan];
    end
end

N_GL = [x_GL_main(:) y_GL_main(:) Fb(x_GL_main(:),y_GL_main(:))];

%% find melt nodes
NodesDownstreamOfGroundingLines="Relaxed";%"Strickt" or "Relaxed";
[LakeNodes,OceanNodes] = LakeOrOcean3_JDR(CtrlVar,MUA,F.GF,GLgeo,GLnodes,GLele,[],NodesDownstreamOfGroundingLines);
x_m = x; y_m = y; b_m = b;
ind_toremove = find(LakeNodes | ~OceanNodes);
ind_meltnodes = find(~LakeNodes & OceanNodes);
x_m(ind_toremove) = [];
y_m(ind_toremove) = [];
b_m(ind_toremove) = [];
N_m = [x_m(:) y_m(:) b_m(:)];

%% offset z-coordinate of N_m with large negative number (following Rosier et al. 2024)
N_m_adj = N_m;
N_m_adj(:,3) = N_m_adj(:,3)-1e6; % CAUTION: the plume paths are sensitive to the offset, 
% which was chosen to give reasonable-looking result for Antarctica. This
% number can be changed.

%% multiply z-coordinate of N_GL to vertically stretch GL geometry (following Rosier et al. 2024)
N_GL_adj = N_GL;
N_GL_adj(:,3) = N_GL_adj(:,3)*2; % CAUTION: the plume paths are sensitive to the multiplier, 
% which was chosen to give reasonable-looking result for Antarctica. This
% number can be changed.

%% k nearest neighbour search with default Euclidean metric. 
% Choose k sufficiently large such that a positive mean slope can (almost) always
% be found
k=100;
[idx, dist] = knnsearch(N_GL_adj,N_m_adj,'K',k);

%% calculate local & global slopes
[dbdx,dbdy,xint,yint]=calcFEderivativesMUA(b,MUA,CtrlVar); % local slope at integration points
[dbdx,dbdy]=ProjectFintOntoNodes(MUA,dbdx,dbdy); % local slope at nodes
dbdx_m = dbdx; dbdy_m = dbdy;
dbdx_m(ind_toremove)=[]; % retain only melt nodes
dbdy_m(ind_toremove)=[];
theta_local = atan(hypot(dbdx_m,dbdy_m));

theta_mean = 0*theta_local;
k_GL = 0*theta_local;
z_GL = 0*theta_local;

for ii=1:size(N_m,1) 
    dx = N_m(ii,1)-N_GL(idx(ii,:),1);
    dy = N_m(ii,2)-N_GL(idx(ii,:),2);
    dz = N_m(ii,3)-N_GL(idx(ii,:),3);
    dl = hypot(dx,dy);
    theta_mean_tmp = atan(dz./dl);
    % find GL nodes for which mean slope is positive.
    Ipos = find(theta_mean_tmp>0);
    if isempty(Ipos)
        theta_mean(ii) = eps; % no positive mean slope could be found
        k_GL(ii) = nan;
        z_GL(ii) = N_m(ii,3); % plume source has same depth as melt node
    else
        theta_mean(ii) = theta_mean_tmp(Ipos(1));
        k_GL(ii) = idx(ii,Ipos(1));
        z_GL(ii) = N_GL(k_GL(ii),3); % depth of plume source
    end
end

%% PLOTTING
if doplots
    
    figure; hold on;

    slope_tmp = 0*x+nan;
    slope_tmp(ind_meltnodes) = theta_mean;
    PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,slope_tmp); hold on;
    %PlotMuaMesh(CtrlVar,MUA,[],'color',[0.8 0.8 0.8]); hold on;
    g(1)=plot(x_GL,y_GL,'-m'); % original GL
    g(2)=plot(x_GL_main,y_GL_main,'xb','markersize',2); % GL pinning points removed - these are
    %the origins of the plume
    axis off;
    
    % XLim=[-1.695e+06 -1.524e+06];
    % YLim=[-3.844e+05 -2.498e+05];
    % xlim(XLim); ylim(YLim);
    % 
    % % plot some pathways
    % I_toplot = find(~isnan(k_GL));
    % for ii=1:10:numel(I_toplot)
    %     plot([N_m(I_toplot(ii),1) N_GL(k_GL(I_toplot(ii)),1)],[N_m(I_toplot(ii),2) N_GL(k_GL(I_toplot(ii)),2)],':k');
    % end
    
    caxis([0 0.025])
    %CM=othercolor('RdBu11');
    %colormap(CM);
    cb.Label.String='slope';
    
    title('Global slope');
    legend(g(:),["GL","Possible plume origins"],"Location","southwest");
    
    figure; hold on;
    
    slope_tmp = 0*x+nan;
    slope_tmp(ind_meltnodes) = theta_local;
    PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,slope_tmp); hold on;
    %PlotMuaMesh(CtrlVar,MUA,[],'color',[0.8 0.8 0.8]); hold on;
    g(1)=plot(x_GL,y_GL,'-m'); % original GL
    g(2)=plot(x_GL_main,y_GL_main,'xb','markersize',2); % GL pinning points removed - these are
    %the origins of the plume
    axis off;
    
    % XLim=[-1.695e+06 -1.524e+06];
    % YLim=[-3.844e+05 -2.498e+05];
    % xlim(XLim); ylim(YLim);
    % 
    
    caxis([0 0.04])
    %CM=othercolor('RdBu11');
    %colormap(CM);
    cb=colorbar;
    cb.Label.String='slope';
    
    title('Local slope');
    legend(g(:),["GL","Possible plume origins"],"Location","southwest");
end
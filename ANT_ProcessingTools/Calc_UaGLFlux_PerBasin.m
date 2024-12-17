function [B,GL] = Calc_UaGLFlux_PerBasin(MUA,F,GF,B,CtrlVar)

%% INPUTS: MUA, F, GF, B, CtrlVar
%% B is a structure with basin information. One way to construct B is to use 
%% basin outlines from IMBIE (routines from C Greene):

%% 
%% filename = 'basins_IMBIE_v2.mat'; 
%% B = load(filename);
%% B = RemoveIceRisesAndIslands(B);

% GL flux
GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ; ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
[xGL,yGL]=LineUpEdges2(CtrlVar,xa,xb,ya,yb);

% add breaks in middle of Ross, FilchnerRonne and Peninsula:
breakx = [-2.11e6 -6.13e5 -2.6e5];
breaky = [1.065e6 3.37e5 -4.35e5];
for ii=1:numel(breakx)
    [~,I]=min(hypot(xGL-breakx(ii),yGL-breaky(ii)),[],'all');
    xGL = [xGL(1:I); nan; xGL(I+1:end)];
    yGL = [yGL(1:I); nan; yGL(I+1:end)];
end

% GL fluxes will be defined at midpoints
[xGLmid,yGLmid,qGL,qGLx,qGLy,Fub,Fvb,Fr,Fh,LakeNodes,GLgeo,ubGL,vbGL] = FluxAcrossGroundingLine_JDR(CtrlVar,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho,[],[],[],[],[],xGL,yGL,GLgeo);

%figure(111); hold on; quiver(xGLmid,yGLmid,qGLx/3e7,qGLy/3e7,'off');

%% keep n longest GL segments
I = find(isnan(xGL)); I = [0;I(:);numel(xGL)+1];
% make structure for individual GL segements
for ii=1:numel(I)-1
    GL(ii).x = xGL(I(ii)+1:I(ii+1)-1);
    GL(ii).y = yGL(I(ii)+1:I(ii+1)-1);
    if numel(GL(ii).x)>1
        GL(ii).xmid = xGLmid(I(ii)+1:I(ii+1)-2);
        GL(ii).ymid = yGLmid(I(ii)+1:I(ii+1)-2);
        GL(ii).qGL = qGL(I(ii)+1:I(ii+1)-2);
        GL(ii).l = hypot(GL(ii).x(2:end)-GL(ii).x(1:end-1),GL(ii).y(2:end)-GL(ii).y(1:end-1)); 
        L = cumsum(GL(ii).l);
        GL(ii).L = L(end);
        GL(ii).ind = 0*GL(ii).x+ii;        
    else
        GL(ii).xmid = GL(ii).x;
        GL(ii).ymid = GL(ii).y;
        GL(ii).qGL = 0;
        GL(ii).l = 0;
        GL(ii).L = 0;
        GL(ii).ind = ii;
    end
    GL(ii).basinindex = [];
end
n_tokeep = numel(GL); % keep all segments
[~,I] = maxk([GL(:).L],n_tokeep);

%for ii=1:numel(I)
%    GLtmp = GL(I(ii));
%    plot(GLtmp.x,GLtmp.y,'-r','linewidth',2); 
%end

GL = GL(I);

%% check proximity to basins
xGL = cell2mat({GL(:).xmid}'); yGL = cell2mat({GL(:).ymid}'); qGL = cell2mat({GL(:).qGL}');
PQ = [xGL yGL]; dist = [];
for ii=1:numel(B.x)
    P = [B.x{ii} B.y{ii}];
    [~,dist(ii,:)] = dsearchn(P,PQ);
end
[~,I] = min(dist,[],1); 

%CM = jet(numel(B.x));
%[~,CM] = cmpermute([],CM);
for ii=1:numel(B.x)
    B.GLpoints{ii} = find(I==ii);
    B.xGL{ii} = xGL(B.GLpoints{ii});
    B.yGL{ii} = yGL(B.GLpoints{ii});
    B.qGL{ii} = qGL(B.GLpoints{ii});
    %figure(111); hold on;
    %plot(B.xGL{ii},B.yGL{ii},'.','color',CM(ii,:))
end

for ii=1:numel(GL)
    basinindex = I(1:numel(GL(ii).xmid));        
    GL(ii).basinindex = mode(basinindex);
    I(1:numel(GL(ii).xmid))=[];
end

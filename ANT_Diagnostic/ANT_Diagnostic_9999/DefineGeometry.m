function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

persistent FB Fs Fb;

if nargin<5
    FieldsToBeDefined='sbSB';
end

alpha=0 ;
fprintf('Loading s, b, S and B \n');

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

%% first deal with the grounded ice
S = 0*x;
load(UserVar.GIGeometryInterpolants,'FB','Fs','Fb');
B = FB(x,y);
B = inpaint_nans(B,4);
s = Fs(x,y);
s = inpaint_nans(s,4);
b = Fb(x,y);
b = inpaint_nans(b,4);
%% where is the grounding line?
h = s-b;
[~,rho,rhow,~]=DefineDensities(UserVar,CtrlVar,MUA,[],[],[],[],[],[]);
[b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow);
%% find floating nodes that are not lakes
[GF,~,~,~]=IceSheetIceShelves(CtrlVar,MUA,GF);
[LakeNodes,~,~,~] = LakeOrOcean3(CtrlVar,MUA,GF);
ISnodes = find(GF.node==0 & GF.NodesDownstreamOfGroundingLines & ~LakeNodes); 

%% now deal with floating ice
load(UserVar.ISGeometryInterpolants,'Fs','Fb');
xIS = x(ISnodes); yIS = y(ISnodes);
s(ISnodes) = Fs(xIS,yIS);
s(ISnodes) = inpaint_nans(s(ISnodes),4);
b(ISnodes) = Fb(xIS,yIS);
b(ISnodes) = inpaint_nans(b(ISnodes),4);


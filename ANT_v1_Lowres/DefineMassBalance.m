function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent DTacc SMB_av Fab

%as = 0*x;
%ab = 0*x;

x_new = MUA.coordinates(:,1) ;
y_new = MUA.coordinates(:,2);

if isempty(DTacc)
    %RACMO
    load('./Interpolants/SMB_RACMO_1979_2013.mat');
    DTacc = DelaunayTri(xRACMO,yRACMO);
end
as=Grid1toGrid2(DTacc,SMB_av,x_new,y_new,CtrlVar);

%figure; scatter(x_new,y_new,[],as);
% define your geographical area: 
%I = find(x_new>1000km);
% change as for those coordinates:
%as(I) = 2*as(I);

if isempty(Fab)
    load('./Interpolants/BalancedMeltRates.mat');
    Fab = scatteredInterpolant(x,y,ab);
end
ab = Fab(x_new,y_new);
ab = ab;
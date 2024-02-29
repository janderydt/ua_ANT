function [Ind,NearestPoints] = ShiftFront(p, coo, FdesiredEleSize)

%INPUTS
% p: icefront 1, to be adjusted
% coo: icefront 2, stay fixed

% OUTPUTS
% Ind: nx1 index vector of points in p to be adjusted
% NearestPoint: nx2 matrix of (x,y) coordinates for nearest points in coo

% DEFINE TOLERANCES
nTol = FdesiredEleSize(p(:,1),p(:,2));
Ind = [];
NearestPoints = [];

% FIND NEAREST POINTS ON COO AND CHECK EUCLIDIAN DISTANCE AGAINST TOLERANCES
for ii=1:size(p,1)
    L = hypot(p(ii,1)-coo(:,1),p(ii,2)-coo(:,2));
    ind = find(L<nTol(ii));
    if ~isempty(ind)
        dist=[];
        for jj=1:numel(ind)
            dist(jj) = L(ind(jj));
        end
        [~,mindistI] = min(dist);
        Ind = [Ind; ii];
        NearestPoints = [NearestPoints; coo(ind(mindistI),:)];
    end
end
function [Ind,NearestPoints] = Node_to_face_distance(p, coo, FdesiredEleSize,meshmin)

%INPUTS
% p: icefront 1, to be adjusted
% coo: icefront 2, stay fixed

% OUTPUTS
% Ind: nx1 index matrix of segments in p that need to be adjusted. 
% NearestPoint: nx2 matrix of (x,y) coordinates in coo that need to be
% added to p
% E.g. NearestPoint(i,:) needs to be added to segment Ind(i)

% DEFINE SEGMENTS
seg.start = p;
seg.end = [p(2:end,:); p(1,:)];
seg.mid = 0.5*(seg.start+seg.end);
seg.d = [seg.end - seg.start];
seg.n = [-seg.d(:,2) seg.d(:,1)];
seg.n = seg.n./repmat(vecnorm(seg.n,2,2),1,2);

% DEFINE TOLERANCES
seg.nTol = max(FdesiredEleSize(seg.mid(:,1),seg.mid(:,2))/3,meshmin);

% DEFINE TRANSVERSE STRIP TO SEGMENTS with WIDTH = seg.nTol+1
strp.p1 = seg.start - (repmat(seg.nTol,1,2)+1).*seg.n;
strp.p2 = seg.start + (repmat(seg.nTol,1,2)+1).*seg.n;
strp.p3 = seg.end + (repmat(seg.nTol,1,2)+1).*seg.n;
strp.p4 = seg.end - (repmat(seg.nTol,1,2)+1).*seg.n;

% FIND POINTS ON COO in the transverse strip, calculate normal distance,
% and add to output vector if within tolerance
Ind = []; NearestPoints = [];

for ii=1:size(strp.p1,1) % loop through all strips

    stripx = [strp.p1(ii,1) strp.p2(ii,1) strp.p3(ii,1) strp.p4(ii,1) strp.p1(ii,1)];
    stripy = [strp.p1(ii,2) strp.p2(ii,2) strp.p3(ii,2) strp.p4(ii,2) strp.p1(ii,2)];
    Icoo = find(inpoly2([coo(:,1) coo(:,2)],[stripx(:) stripy(:)])); % find elements of coo within strip

    a = [seg.start(ii,:) - seg.end(ii,:) 0]; % direction of line segment
    for pp=1:numel(Icoo)    
        b = [coo(Icoo(pp),:) - seg.end(ii,:) 0];
        d = norm(cross(a,b))/norm(a); % normal distance of coo(Icoo(pp),:) to line segment
        if d>eps(0) && d< seg.nTol(ii) % if normal distance to line segment is smaller than tolerance...
            Ind = [Ind; ii];% add to output vector
            NearestPoints = [NearestPoints; coo(Icoo(pp),:)];
        end
    end
end
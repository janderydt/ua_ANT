
%%
%
% Determines if a point p is close to a line segment defined by the end points 
% defined in coo
%
% Also calculates normal and along distances from a point to a line segment.
%
%
% The input format is as follows:
%
%   [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, coo ,
%   [],FdesireEleSize)
%
% where coo is a vector defining the end points of individual connected line segments.
%
% A point is considered to be 'close to' the line segment if:
%
% # if the distance along the line segment is between 0 and 1 (normalized
% units)
% # the distance to one or both of the end points of the line segment is less
% than nTol, where nTol is calculated based on the desired element size around 
% point p, using FdesiredEleSize
%
% Inputs:
%
% p  : Nx2 matrix of with the (x,y) cooridnates of points
%
% coo   : Mx2 matrix defining the (x,y) coordinates of connected line segments.
%
% FdesiredEleSize is an interpolant that specifies the spatially variable tolerances. 
%
% Outputs:
%
% Ind : points in p that are within tolerance to a line segment.
% NearestSegment : for each point that is within tolerance of a line segment, this output
%   provides the corresponding segment
% NearestPointOnSegment: (x,y) coordinates of the nearest point on the B-A
%   segment. This is calculated as follows: 
% 

%%% OLD APPROACH, DO NOT USE
%%%   1. If the AlongDirection is less than NormalTolerance, shift the point to the nearest point on the yr2 boundary
%%%   2. If the AlongDistance is more than NormalTolerance, shift the point to that line segment in the normal direction
% 

B=[coo([2:end 1],1) coo([2:end 1],2) ];
A=[coo(1:end,1) coo(1:end,2)] ;

M=size(A,1);
N=size(p,1);

Ind=[]; NearestSegment=[]; NearestPointOnSegment=[];

for J=1:N % cycle through all the points in p

    nTol = FdesiredEleSize(p(J,1),p(J,2));
    
    d=B-A; % segments in coo (boundary to compare to)
    l=vecnorm(d,2,2); % length of the segments
    
    f = 0*A; e = 0*A;
    f(:,1)=d(:,2) ; f(:,2)=-d(:,1);
    e(:,1)=p(J,1)-A(:,1);
    e(:,2)=p(J,2)-A(:,2);
    
    dalong=dot(e,d,2)./l.^2;    % distance along the segments, if point p within the `transverse' strip to b-a, then 0<dalong<1
    dnorm=abs(dot(e,f,2)./l);   % absolute distance normal to each segment in B-A

    % find segments for which p falls within the transverse strip
    Itemp = find(dalong>=0 & dalong<=1);

    % find the nearest segment to point p in the normal direction
    [MINdnorm,I2temp] = min(dnorm(Itemp));
    if isempty(I2temp)
        MINdnorm = nTol+1e10;
    end

    % define nTol, depending on the region of the domain
%     [latp,lonp] = psxy2ll(p(J,1),p(J,2),-71,0);
%     nTol = IF.defaultds ; % default
%     nTol_adjusted = 0; rr=1;
%     while nTol_adjusted == 0 && rr<=numel(IF.segments)
%         lonlimits = [IF.segments(rr).lonlim(1), IF.segments(rr).lonlim(1),...
%                 IF.segments(rr).lonlim(2) IF.segments(rr).lonlim(2)];
%         latlimits = [IF.segments(rr).latlim(1), IF.segments(rr).latlim(2),...
%                 IF.segments(rr).latlim(2) IF.segments(rr).latlim(1)];
%         if inpoly([lonp latp],[lonlimits(:) latlimits(:)])
%             nTol = IF.segments(rr).ds;
%             nTol_adjusted = 1;
%         end
%         rr = rr+1;
%     end
        
    % check if dnorm is less than the tolerance. If so, calculate new
    % position for that point, otherwise leave unmodified
    %%if MINdnorm <= nTol

    if ~isempty(I2temp)

        SInd = Itemp(I2temp);
        lA = hypot(A(SInd,1)-p(J,1),A(SInd,2)-p(J,2)); 
        lB = hypot(B(SInd,1)-p(J,1),B(SInd,2)-p(J,2));

    else

        lA=0; lB=0;
    end

    if lA<nTol/2 || lB<nTol/2

        Ind = [Ind; J];
        NearestSegment = [NearestSegment; SInd];

        if lA<lB
            % move point to A
            d

        else
            % move point to B
            NearestPointOnSegment = [NearestPointOnSegment; B(SInd,:)];
        end
        %%% OLD APPROACH, DO NOT USE
        %AlongDist = dalong(SInd)*l(SInd);
        %if AlongDist <= nTol 
            % move point to A
        %    NearestPointOnSegment = [NearestPointOnSegment; A(SInd,:)];
        %elseif l(SInd)-AlongDist <= nTol
            % move point to B
        %    NearestPointOnSegment = [NearestPointOnSegment; B(SInd,:)];

       
       %%% else
            % move to intersection of segment and normal
       %%%     NearestPointOnSegment = [NearestPointOnSegment; p(J,:)+dnorm(SInd)*f(SInd,:)/norm(f(SInd,:))];
        end

    end
    
end

function vg = ExtrudeField(x,y,vx,vy,x_orig,y_orig,v_orig,step,nsteps)

% Adapted from Chad Greene
[X,Y] = meshgrid(x_orig,y_orig);

perim = bwperim(isfinite(v_orig));
xseed = X(perim); 
yseed = Y(perim); 
[row,col] = find(perim); 
XY = stream2(x,y,vx,vy,xseed,yseed,[step nsteps]);

clear x y vx vy X Y 

% Remove empty entries
IndEmpty = find(cellfun(@isempty,XY));
XY(IndEmpty)=[];
% xseed(IndEmpty)=[];
% yseed(IndEmpty)=[];
row(IndEmpty)=[];
col(IndEmpty)=[];

clear IndEmpty

% ORIGINAL CODE, which is a bit wasteful on memory
% Loop through each seed location to extrapolate seed value along flowlines:  
% V = XY;
% for k = 1:length(V) 
%     V{k}(:,1) = v_orig(row(k),col(k));
% end
% 
% M = cell2nancat(XY);
% V = cell2nancat(V); 
% 
% % Are the seed locations captured in M? I'm not sure, so for good measure we'll create these arrays:  
% xx = [X(perim);M(:,1)] ;
% yy = [Y(perim);M(:,2)] ;
% vv = [v_orig(perim);V(:,1)] ;

V = XY;
for k = 1:length(V) 
     V{k}(:,1) = v_orig(row(k),col(k));
end
V = cell2nancat(V);
vv = [v_orig(perim);V(:,1)] ;

clear V v_orig;

XY = cell2nancat(XY);
xx = [xseed;XY(:,1)] ;
yy = [yseed;XY(:,2)] ;

clear XY perim

isf = isfinite(vv);
xx = xx(isf);
yy = yy(isf);
vv = vv(isf);

clear isf

% Grid up the streamline data:  
vg = gridbin(xx,yy,vv,x_orig,y_orig);

clear xx yy vv x_orig y_orig
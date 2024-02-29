function vg = ExtrudeField(x,y,vx,vy,x_orig,y_orig,v_orig,step,nsteps)

v = hypot(vx,vy); 
[X,Y] = meshgrid(x_orig,y_orig);

perim = bwperim(isfinite(v_orig));
xseed = X(perim); 
yseed = Y(perim); 
[row,col] = find(perim); 
XY = stream2(x,y,vx,vy,xseed,yseed,[step nsteps]);

% Loop through each seed location to extrapolate terminus velocity and ice shelf name along flowlines:  
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = v_orig(row(k),col(k));
end

M = cell2nancat(XY); % cell2nancat is a function at the bottom of this script. 
V = cell2nancat(V); 

% Are the seed locations captured in M? I'm not sure, so for good measure we'll create these arrays:  
xx = [X(perim);M(:,1)] ;
yy = [Y(perim);M(:,2)] ;
vv = [v_orig(perim);V(:,1)] ;

% Remove NaNs: 
isf = isfinite(vv);
xx = xx(isf);
yy = yy(isf);
vv = vv(isf);

% Grid up the streamline data:  
vg = gridbin(xx,yy,vv,x_orig,y_orig); % gridbin is on my GitHub. It bins the average velocity observation in each grid cell. It's not as beautiful as gridfit, but it's the only way to deal with nearly a billion scattered data points.   

clear v X Y perim xseed yseed XY V M xx yy vv isf;
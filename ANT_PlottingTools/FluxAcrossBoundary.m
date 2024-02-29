function [xB,yB,qB]=FluxAcrossBoundary(X,Y,Fub,Fvb,Fh,Fr)

% INPUT: X and Y are (n,2) dimensional arrays with x/y coordinates of the endpoints 
% of n edges through which the flux should be calculated. 
% OUTPUT: xB and yB are the midpoints of the edges, qB the flux in kg/yr

Inan = find(isnan(X(:,1))|isnan(X(:,2))|isnan(Y(:,1))|isnan(Y(:,2))); 
fprintf("Removing %s nan values.\n",string(numel(Inan)));

if ~isempty(Inan)
    X(Inan,:)=[];
    Y(Inan,:)=[];
end

ds=hypot(X(:,2)-X(:,1),Y(:,2)-Y(:,1));

Izerods = find(ds<eps);
fprintf("Removing %s segments with zero length.\n",string(numel(Izerods)));

ds(Izerods)=[];
X(Izerods,:)=[];
Y(Izerods,:)=[];
    
nx = (Y(:,2)-Y(:,1))./ds;
ny = -(X(:,2)-X(:,1))./ds;

xB = 0.5*(X(:,1)+X(:,2));
yB = 0.5*(Y(:,1)+Y(:,2));

rhoB=Fr(xB,yB);
ubB=Fub(xB,yB); 
vbB=Fvb(xB,yB); 
hB=Fh(xB,yB); 

qB = ds.*hB.*rhoB.*(ubB.*nx+vbB.*ny);  % units: m kg/m^2 (m/yr) = kg/yr

qBx=qB.*nx; qBy=qB.*ny;


end

addpath(genpath('/media/wchm8/data1-JDeRydt/Antarctic_datasets'));

[x,y,~] = bedmap2_data('surface','xy');

mask = bedmap2_data('icemask');

[C,h] = contour(x,y,mask,[127 127]);

ii=1; kk=1;
while ii<size(C,2)
   step = C(2,ii);
   Cont(kk).step = step;
   Cont(kk).x = C(1,ii+1:ii+step);
   Cont(kk).y = C(2,ii+1:ii+step);
   ii = ii + step + 1;
   kk = kk+1;   
end

[Maxstepsize,I] = max([Cont(:).step]);

CtrlVar.GLtension=0.5;
CtrlVar.GLds=40e3;
[X,Y,nx,ny,~] = Smooth2dPos(Cont(I).x(:),Cont(I).y(:),CtrlVar);

figure; hold on;
plot(Cont(I).x(:),Cont(I).y(:),'ok');
plot(X,Y,'+r');

MeshBoundaryCoordinates = [X(:) Y(:)];

save('./BoundaryCoordinates.mat','MeshBoundaryCoordinates');


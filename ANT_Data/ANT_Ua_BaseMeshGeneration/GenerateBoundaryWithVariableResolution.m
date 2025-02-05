function [xB,yB,FdesiredEleSize] = GenerateBoundaryWithVariableResolution(velocityfile,geometryfile,xB_init,yB_init,FdesiredEleSize,elesizemetric,meshmin,meshmax)

% Create boundary coordinates with variable spacing, starting from xB_init
% and yB_init as an initial boundary, and using the gridded interpolant
% FdesiredEleSize(X,Y,EleSize) to identify what resolution is required.
% FdesiredEleSize is an optional input. If it is not provided, a new
% gridded interpolant will be constructed, based on a range of criteria
% such as speed, floatation, thickness_gradient etc, defined in elesizemetric.

persistent Fus Fvs Fmask;

if nargin < 2
    FdesiredEleSize = [];
    geometry = 'antarctica'; % square, antarctica
else
    if isempty(xB_init)
        geometry = 'antarctica';
    else
        geometry = 'user';
    end
end

if nargin < 4
    elesizemetric = ["speed","floatation","thickness_gradient","effective_strain_rate_gradient"]; 
    % any combination of {speed, floatation, thickness_gradient, effective_strain_rate_graident}
end

if nargin < 5
    meshmin = 5e3;
    meshmax = 100e3;
end

%% define a base mesh and matrix with desired element sizes
switch geometry
    case 'square'
        % something very simple
        fprintf('A simple square.\n');
        dx=1e-2; dy=1e-2;
        x = [0:dx:1];
        y = [0:dy:1];
    case 'antarctica'
        % something a bit more complicated
        fprintf('  > Coarse resolution Antarctic boundary.\n');
        load('/mnt/md0/Ua/cases/ANT/ANT_Lowres/BoundaryCoordinates.mat');
        dx=1e3; dy=1e3;
        x = [min(MeshBoundaryCoordinates(:,1))-10*dx:dx:max(MeshBoundaryCoordinates(:,1))+10*dx];
        y = [min(MeshBoundaryCoordinates(:,2))-10*dy:dy:max(MeshBoundaryCoordinates(:,2))+10*dy];
    case 'user'
        fprintf('  > User has specified initial boundary.\n');
        dx=1e3; dy=1e3;
        x = [min(xB_init)-100*dx:dx:max(xB_init)+100*dx];
        y = [min(yB_init)-100*dy:dy:max(yB_init)+100*dy];
    otherwise
        fprintf('  > Geometry unknown.\n');
        error('Stopping.');
end

%% define boundary with some nominal resolution
switch geometry
    case 'square'
        nx=2; ny=2;
        xB = linspace(x(1),x(end),nx);
        yB = linspace(y(1),y(end),ny);
        boundary = [xB(1)+0*yB(2:end) xB(1:end-1) xB(end)+0*yB(2:end) flipdim(xB(1:end),2);...
            yB(1:end-1) yB(end)+0*xB(2:end) flipdim(yB(2:end),2) yB(1)+0*xB(1:end)]';
        elesizemetric = "checkerboard";
    case 'antarctica'
        boundary = MeshBoundaryCoordinates;
    case 'user'
        boundary = [xB_init(:) yB_init(:)];
    otherwise 
        fprintf('  > Geometry unknown.\n');
        error('Stopping.');
end

[X,Y] = meshgrid(x,y);

if isempty(FdesiredEleSize)
    %% Define desired element sizes
    desiredEleSize = 0*X+meshmax;
    
    for ii=1:numel(elesizemetric)

        metric = elesizemetric(ii);
    
        switch metric
            case 'checkerboard'
                desiredEleSize(1:numel(x)/2,1:numel(y)/2)=dx;
                desiredEleSize(numel(x)/2:end,1:numel(y)/2)=4*dx;
                desiredEleSize(1:numel(x)/2,numel(y)/2:end)=4*dx;
                desiredEleSize(numel(x)/2:end,numel(y)/2:end)=dx;
    
            case 'speed'
                fprintf('  > Implementing speed for desiredEleSize matrix...');
                if isempty(Fus)
                    load(velocityfile,"Fus","Fvs");
                    if isempty(Fus)
                        Fu = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb));
                        U = Fu(X,Y);
                    else
                        U = hypot(Fus(X,Y),Fvs(X,Y));
                    end
                end
%                 if ~exist('../ANT_Interpolants/ScatteredInterpolant_Nearest_FU_2016.mat','file')
%                     fprintf('Generating ScatteredInterpolant_Nearest_MeaSUREs_FU_2016.mat\n');
%                     load("../ANT_Interpolants/GriddedInterpolants_2015-2016_MeaSUREs_ITSLIVE_velocities.mat","Fus","Fvs");
%                     Fx = Fus.GridVectors{1}; Fy = Fus.GridVectors{2};   
%                     %Fus.Method = 'nearest'; Fus.ExtrapolationMethod = 'nearest';
%                     %Fvs.Method = 'nearest'; Fvs.ExtrapolationMethod = 'nearest';
%                     U = sqrt(Fus.Values.^2+Fvs.Values.^2);
%                     I = find(~isnan(U)); [FX,FY]=ndgrid(Fx,Fy);
%                     FU = scatteredInterpolant(FX(I),FY(I),U(I),'nearest','nearest');
%                     save('../ANT_Interpolants/ScatteredInterpolant_Nearest_MeaSUREs_FU_2016.mat','FU','-v7.3','-nocompression');
%                 else
%                     load('../ANT_Interpolants/ScatteredInterpolant_Nearest_MeaSUREs_FU_2016.mat');
%                 end
                 
                % evolve mesh resolution from meshmax km where U<1m/yr to meshmin km for
                % U>1km/yr
                desiredEleSize(U<1) = meshmax;
                % log dependency of mesh size on speed
                I = find(U>=1 & U<15);
                desiredEleSize(I) = meshmax + (10e3-meshmax)/log10(15)*log10(U(I));
                I = find(U>=15 & U<250);
                desiredEleSize(I) = 10e3 + (2*meshmin-10e3)/(log10(250)-log10(15))*(log10(U(I))-log10(15));               
                desiredEleSize(U>=250) = 2*meshmin;
                desiredEleSize(desiredEleSize==0) = meshmax;
                fprintf('Done.\n');
    
            case 'floatation' % floating areas have meshmin km resolution or more
                fprintf('  > Implementing floatation for desiredEleSize matrix...');
                if isempty(Fmask)
                    load(geometryfile,"Fmask");
                end
                Fmask.Method = 'nearest'; Fmask.ExtrapolationMethod = 'nearest';
                Mask = Fmask(X,Y); % 3 is floating ice, 0 is open ocean
                desiredEleSize((Mask==3 | Mask==0) & desiredEleSize>2*meshmin) = 2*meshmin; % meshmin km resolution for ice shelves
                fprintf('Done.\n');
    
            case 'thickness_gradient'  % unrefine some areas that are very flat, such as the large ice shelves   
                fprintf('  > Implementing thickness gradients for desiredEleSize matrix...');
                if isempty(Fmask)
                    load(geometryfile,"Fmask");
                end
                load(geometryfile,"Fs","Fb");
                Fs.Method = 'nearest'; Fs.ExtrapolationMethod = 'nearest';
                Fb.Method = 'nearest'; Fb.ExtrapolationMethod = 'nearest';
                
                Mask = Fmask(X,Y);
                h = Fs(X,Y)-Fb(X,Y); 
                [Dhx,Dhy] = gradient(h,x,y);
                Dh = 0.5*(Dhx+Dhy);
                % do some smoothing
                Dh = imgaussfilt(abs(Dh),'FilterSize',201);
                I = find(Dh<0.005 & (Mask==3 | Mask==0));
                desiredEleSize(I) = max(desiredEleSize(I),meshmin*4); % lower resolution for floating regions with low thickness gradient
                fprintf('Done.\n');

            case 'effective_strain_rate_gradient'

                fprintf('  > Implementing effective strain rate gradients for desiredEleSize matrix...');
                
                % calc effective strain rate
                [dudx,dudy] = gradient(Fus(X,Y),x,y);
                [dvdx,dvdy] = gradient(Fvs(X,Y),x,y);
                e=real(sqrt(dudx.^2+dvdy.^2+(dudy+dvdx).^2/2));
                % do some smoothing
                e = imgaussfilt(abs(e),'FilterSize',501);
                I = find(e>0.03);
                desiredEleSize(I) = min(desiredEleSize(I),meshmin);
                         
                fprintf('Done.\n');
                
            otherwise
                error('unkown case \n');
        end
    
    
    end
    
    desiredEleSize = movmean(desiredEleSize,[10,10]);
    [Xn,Yn] = ndgrid(x,y);
    FdesiredEleSize = griddedInterpolant(Xn,Yn,desiredEleSize','linear','nearest');

end

desiredEleSize = FdesiredEleSize(X,Y);

%% construct resampled boundary
fprintf('  > Constructing resampled boundary...');
dL = [0; cumsum(hypot(diff(boundary(:,1),1,1),diff(boundary(:,2),1,1)))];
% now use boundary points as interpolant, and regrid with desired spacing
% define start point as boundary(1,2)
dx0x1=1e10;
x0 = boundary(1,1); y0 = boundary(1,2);
cutoff = 0.99*interp2(X,Y,desiredEleSize,x0,y0);
ds = 0; xnew = x0; ynew = y0;

while dx0x1 > cutoff
    step = interp2(X,Y,desiredEleSize,x0,y0);
    ds = ds+step;
    % walk distance ds along boundary and add increments of 100m to ds as long as
    % hypot(xnew-x0,ynew-y0)<step
    Ltmp = 0;
    while Ltmp < step
        finalPathXY = interp1(dL, boundary, ds);
        xtmp = finalPathXY(1);
        ytmp = finalPathXY(2);   
        Ltmp = hypot(xtmp-x0,ytmp-y0);
        ds = ds + 100;
    end
    ds = ds-100;
    xnew = [xnew xtmp];
    ynew = [ynew ytmp];
    dx0x1 = hypot(xtmp-xnew(1),ytmp-ynew(1));
    x0 = xtmp; y0 = ytmp;
end
fprintf('Done.\n');

xB = xnew;
yB = ynew;

%plot(boundary(:,1),boundary(:,2),'.-k');
%plot(xB,yB,'x-r');
%axis equal;

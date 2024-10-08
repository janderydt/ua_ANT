function ANT_GenerateBaseMesh

warning('off','all');

%velocityfile = "../ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED";
%velocityfile = "/mnt/md0/Ua/cases/ANT/ANT_Inverse/ANT_Inverse_1035/ANT_Inverse_1035-RestartFile.mat";
velocityfile = "../ANT_Interpolants/GriddedInterpolants_2018-2019_MeaSUREs_ITSLIVE_Velocities_EXTRUDED";
%geometryfile = "../ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-2000_EXTRUDED";
geometryfile = "../ANT_Interpolants/GriddedInterpolants_Geometry_01-Jun-2018_EXTRUDED";

% This function generates boundary and internal coordinates for Antarctica,
% based on calving front outlines from C. Greene et al. (2022) for years 
% specified below. You can choose multiple years, in which case you need to
% set ExtrudeMesh=1. The icefronts for different years will be imposed as internal
% boundaries in the extruded mesh. This is probably only useful if you want to 
% study calving processes. Mesh refinement is based on a number of criteria
% that are specified in GenerateBoundaryWithVariableResolution.m. If no
% refinement criteria are defined then a mesh with uniform resolution 
% meshav is produced.

years_to_include = 2018;%[2000 2009 2014 2018]; %tested with 2000
ExtrudeMesh = 0;
VariableBoundaryResolution = 1;
meshmin = 5e3;
meshav = 10e3;
meshmax = 100e3;

if numel(years_to_include)>1 & ExtrudeMesh == 0

        error("This script only works with multiple years and internal boundaries if you extrude the mesh");
end

% first we read data from Green et al. (2022)
froot_data = getenv("froot_data");
load(froot_data+"Antarctica_IceFronts/icemask_composite_CGreene_2022.mat");

year = floor(year);
year = year(:)'; % need row format

figure; hold on;

% assemble ice front locations for different years 
FdesiredEleSize = [];

for ii=1:numel(years_to_include)

    fprintf("======================================\n");
    fprintf("Collecting initial ice front data for %s...",string(years_to_include(ii)));

    I = find(year == years_to_include(ii));
    xtmp = cx{I(1)}; ytmp = cy{I(1)};
    J = find(isnan(xtmp)); J = [0; J; numel(xtmp)+1];
    % extract longest contour line
    dJ = J(2:end)-J(1:end-1); [~,K]=max(dJ);
    icefront(ii).x = xtmp(J(K)+1:J(K+1)-1);
    icefront(ii).y = ytmp(J(K)+1:J(K+1)-1);
    icefront(ii).year = year(I(1));

    fprintf("Done.\n");
    fprintf("======================================\n");

    if VariableBoundaryResolution

        % resample according to DesiredEleSize
        fprintf("Sampling Boundary at variable resolution...\n")
    
        [xtmp,ytmp,FdesiredEleSize] = GenerateBoundaryWithVariableResolution(velocityfile,geometryfile,icefront(ii).x,icefront(ii).y,...
            FdesiredEleSize,{'speed','floatation','thickness_gradient'},meshmin,meshmax);

        fprintf("Done.\n");
        fprintf("======================================\n");

    else
    
        % resample according to DesiredEleSize
        fprintf("Sampling boundary at fixed resolution %sm...",string(meshav));

        xtmp = icefront(ii).x;
        ytmp = icefront(ii).y;
        [xtmp,ytmp] = Smooth2dPos(xtmp,ytmp,1,meshav);
        xtmp = xtmp(:)'; ytmp = ytmp(:)';
        dx = 1e3; dy = 1e3;
        x = [min(xtmp)-100*dx:dx:max(xtmp)+100*dx];
        y = [min(ytmp)-100*dy:dy:max(ytmp)+100*dy];
        [Xn,Yn] = ndgrid(x,y);
        FdesiredEleSize = griddedInterpolant(Xn,Yn,meshav+zeros(size(Xn)),'linear','nearest');
        
        fprintf("Done.\n");
        fprintf("======================================\n");

    end

    % do some smoothing for Thwaites ice front   
    % [~,Imin1] = min(sqrt((xtmp+1557580).^2+(ytmp+553192).^2));
    % [~,Imin2] = min(sqrt((xtmp+1574000).^2+(ytmp+568115).^2));
    % Imin = sort([Imin1,Imin2]);
    % dx = xtmp(Imin(2))-xtmp(Imin(1)); dy = ytmp(Imin(2))-ytmp(Imin(1)); 
    % N = round(sqrt(dx.^2+dy.^2)/5e3); dL = dx/N;
    % xnew = xtmp(Imin(1))+[1:N-1]*dL;
    % ynew = ytmp(Imin(1))+dy/dx*(xnew-xtmp(Imin(1)));
    % xtmp = [xtmp(1:Imin(1)) xnew xtmp(Imin(2):end)];
    % ytmp = [ytmp(1:Imin(1)) ynew ytmp(Imin(2):end)];

    %plot(xtmp,ytmp,'-ok','linewidth',0.5,'MarkerSize',1);

    fprintf("Looking for self-interesecting boundary segments...");

    [~,~,segments]=selfintersect(xtmp,ytmp);
    while ~isempty(segments)
        xtmp([segments(1,1):segments(1,2)])=[];
        ytmp([segments(1,1):segments(1,2)])=[];  
        [~,~,segments]=selfintersect(xtmp,ytmp);
    end

    %fprintf("Looking for and removing local clusters...");

    % D = pdist2([xtmp(:) ytmp(:)],[xtmp(:) ytmp(:)]);
    % D = tril(D,-1); 
    % D(D==0)=nan;
    % 
    % nTol = repmat(FdesiredEleSize(xtmp(:),ytmp(:)),1,numel(xtmp));
    % 
    % mask = D<nTol/2;
    % 
    % [C,h] = contour([1:numel(xtmp)],[1:numel(xtmp)],mask,[0.5 0.5]);
    % kk=1;
    % while size(C,2)>0
    %     l = C(2,1);
    %     cont(kk).x = C(1,2:l+1);
    %     cont(kk).y = C(2,2:l+1);
    %     C(:,1:l+1) = [];
    %     kk=kk+1;
    % end
    % 
    % addpath('/mnt/md0/Matlab/Matlab_Functions/poly_stuff');
    % 
    % figure(111); hold on; plot(xtmp,ytmp,'-xk',LineWidth=0.5);
    % nTol=FdesiredEleSize(xtmp,ytmp);
    % L = hypot(xtmp(2:end)-xtmp(1:end-1),ytmp(2:end)-ytmp(1:end-1));
    % for ii=2:numel(xtmp)-1
    %     if L(ii)<(nTol(ii-1)+nTol(ii)+nTol(ii+1))/6
    %         text(xtmp(ii)/2+xtmp(ii+1)/2,ytmp(ii)/2+ytmp(ii+1)/2,num2str(L(ii),2));
    %         plot([xtmp(ii-1) xtmp(ii+1)],[ytmp(ii-1) ytmp(ii+1)],'--r');
    %     end
    % end
    % for cc=1:numel(cont)   
    %         xctmp = cont(cc).x(:); yctmp = cont(cc).y(:);
    %         [X,Y] = ndgrid([min(floor(xctmp)):max(ceil(xctmp))],...
    %             [min(floor(yctmp)):max(ceil(yctmp))]);
    %         I = find(inpoly([X(:) Y(:)],[xctmp(:) yctmp(:)]));
    %     if numel(I)>1
    %         IndInPoly = unique([X(I),Y(I)]);
    %         figure(111);
    %         plot(xtmp(IndInPoly),ytmp(IndInPoly),'o-r',LineWidth=2);
    %         nTol = mean(FdesiredEleSize(xtmp(IndInPoly),ytmp(IndInPoly)))/2;
    %         dn = floor(hypot(xtmp(IndInPoly(end))-xtmp(IndInPoly(1)),...
    %             ytmp(IndInPoly(end))-ytmp(IndInPoly(1)))/nTol);
    %         if dn>1
    %             xnew = linspace(xtmp(IndInPoly(1)),xtmp(IndInPoly(end)),dn);
    %             ynew = linspace(ytmp(IndInPoly(1)),ytmp(IndInPoly(end)),dn);
    %         else
    %             xnew = [xtmp(IndInPoly(1)),xtmp(IndInPoly(end))];
    %             ynew = [ytmp(IndInPoly(1)),ytmp(IndInPoly(end))];
    %         end
    %         figure(111); plot(xnew,ynew,'--m');
    %     end
    % end

    icefront(ii).x_vards = xtmp;
    icefront(ii).y_vards = ytmp;

    %plot(icefront(ii).x_vards,icefront(ii).y_vards,'-r','linewidth',0.5,'MarkerSize',1);
    icefront(ii).mask = 1+0*icefront(ii).x_vards;

    fprintf("Done.\n");
    fprintf("======================================\n");
end

% now remove any icefront locations that are within tolerance of any other
% icefront location

if numel(years_to_include)>1

    years_to_compare = flipdim(nchoosek([1:numel(years_to_include)],2),2);
    
    for yr=1:size(years_to_compare,1)
        
        yr1 = years_to_compare(yr,1); % year to adjust
        yr2 = years_to_compare(yr,2); % year to compare

        fprintf("Removing points from succesive ice fronts that are within tolerance of each other: %s - %s...\n",string(years_to_include(yr1)),string(years_to_include(yr2)));

        p = [icefront(yr1).x_vards(:) icefront(yr1).y_vards(:)];
        coo = [icefront(yr2).x_vards(:) icefront(yr2).y_vards(:)];
    
        % Now find line segments from yr2 that form the transverse strip to 
        % each
        % boundary point along yr1. Then find the nearest normal distance to those
        % line segments. If that distance is less than IF.segments.ds in that region, then 
        % 1. If the AlongDirection is less than NormalTolerance, shift the point to the nearest point on the yr2 boundary
        % 2. If the AlongDistance is more than NormalTolerance, shift the point to that line segment in the normal direction
        [Ind,NearestPoints] = ShiftFront(p, coo, FdesiredEleSize);

        fprintf("Year %s: found %s points to within tolerance of the %s ice front.\n",...
            string(years_to_include(yr1)),string(numel(Ind)),string(years_to_include(yr2)));
        % shift Ind to point on nearest segment of the other icefront, year(1) ice front
        % remains intact, other ice fronts are adjusted. A mask keeps track of
        % which segments remain intact
        %figure; hold on;
        %plot(icefront(yr2).x_vards(:),icefront(yr2).y_vards(:),'-xb');
        %plot(icefront(yr1).x_vards(:),icefront(yr1).y_vards(:),'--r');

        for ii=1:numel(Ind)
            %[~,J] = min(sqrt((icefront(yr2).x_vards(:)-icefront(yr1).x_vards(Ind(ii))).^2+...
            %    (icefront(yr2).y_vards(:)-icefront(yr1).y_vards(Ind(ii))).^2));
            %icefront(yr1).x_vards(Ind(ii)) = icefront(yr2).x_vards(J(1));
            %icefront(yr1).y_vards(Ind(ii)) = icefront(yr2).y_vards(J(1));  
            %plot(icefront(yr1).x_vards(Ind(ii)),icefront(yr1).y_vards(Ind(ii)),'og');
            %plot(NearestPointOnSegment(ii,1),NearestPointOnSegment(ii,2),'ob');
            %plot([icefront(yr1).x_vards(Ind(ii)),NearestPointOnSegment(ii,1)],...
            %    [icefront(yr1).y_vards(Ind(ii)),NearestPointOnSegment(ii,2)],'-g');
            icefront(yr1).x_vards(Ind(ii)) = NearestPoints(ii,1);
            icefront(yr1).y_vards(Ind(ii)) = NearestPoints(ii,2);
            icefront(yr1).mask(Ind(ii)) = NaN;
        end
        %plot(icefront(yr1).x_vards(:),icefront(yr1).y_vards(:),'-or');

        fprintf("Done.\n");
        fprintf("======================================\n");
        
    end

end

if ExtrudeMesh

    fprintf("Define new ice front that allows for future advance, and/or as an efficient" + ...
        "way to implement internal boundaries along all desired ice fronts...")
    
    % define new mesh boundary of Antarctica, which allows advancing ice
    % fronts. Start from the present-day ice front location, and use
    % present-day velocities to estimate a hypothetical future ice front,
    % assuming that the calving rate is 0. The algorithm takes each point p(x,y)
    % along the present-day ice front, constructs a polyshape circle with 
    % radius v(p)*dt where v is the present-day speed and dt = 100 years. 
    % The algorithm then calculates the union of all the polyshapes, and the
    % outer contour of this union defines the newly expanded outline of
    % Antarctia.
    x = icefront(1).x_vards;
    y = icefront(1).y_vards;
    % resample at regular 5km intervals
    CtrlVar.GLtension=1; % tension of spline, 1: natural smoothing; 0: straight line
    CtrlVar.GLds=5e3 ; 
    [x,y,~,~] = Smooth2dPos(x,y,CtrlVar);
    
    vfile = "../ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED";
    load(vfile,"Fus","Fvs");
    N=5;
    xu = Fus.GridVectors{1}; yu = Fus.GridVectors{2}; 
    xu=xu(1:N:end); yu=yu(1:N:end);
    [Xu,Yu] = ndgrid(xu,yu);
    us = Fus.Values(1:N:end,1:N:end); vs = Fvs.Values(1:N:end,1:N:end);
    I = find(isnan(us));
    Xu(I)=[]; Yu(I)=[]; us(I)=[]; vs(I)=[];
    Fus = scatteredInterpolant(Xu(:),Yu(:),us(:),'nearest');
    Fvs = Fus; Fvs.Values = vs(:);
    vx = Fus(x,y); vy = Fvs(x,y);
    v = hypot(vx,vy);
    
    n=20; dt = 1.5*(years_to_include(end)-years_to_include(1)); % advance ice front at present-day speed for 25 years 

    for ii=1:numel(x)
    
        r = max(v(ii)*dt,15e3); % impose minimum 15km halo
        % xout(ii) = x(ii) + r*vx(ii)/v(ii);
        % yout(ii) = y(ii) + r*vy(ii)/v(ii);
        theta = (0:n-1)*(2*pi/n);
        xcircle = x(ii) + r*cos(theta);
        ycircle = y(ii) + r*sin(theta);
        P(ii) = polyshape(xcircle,ycircle);
    
    end
    
    U = union(P);
    I=find(isnan(U.Vertices));
    xout = U.Vertices(1:I(1)-1,1);
    yout = U.Vertices(1:I(1)-1,2);

    CtrlVar.GLtension=0.5; % tension of spline, 1: natural smoothing; 0: straight line
    CtrlVar.GLds=25e3 ; 
    [xout,yout,~,~] = Smooth2dPos(xout,yout,CtrlVar);
    CtrlVar.GLtension=1; % tension of spline, 1: natural smoothing; 0: straight line
    CtrlVar.GLds=25e3 ; 
    [xouter,youter,~,~] = Smooth2dPos(xout,yout,CtrlVar);
    plot(xouter,youter,'-m');


    fprintf("Done.\n");
    fprintf("======================================\n");


else

    if numel(years_to_include)==1

        xouter = icefront(1).x_vards;
        youter = icefront(1).y_vards;
        I = find(isnan(xouter));
        xouter(I) = []; youter(I) = [];

        %% do some smoothing
        %x=xouter(:) ; y=youter(:);
        %xy=[x';y']; df=diff(xy,1,2);
        %t = cumsum([0, sqrt([1 1]*(df.*df))]);
        %cv = csaps(t,xy,1e-12);
        %X=fnval(cv,t);
        %xouter = X(1,:); youter = X(2,:);

    end

end

% now assemble all the necessary data for meshing:
% 1. outer boundary: xouter, youter
% 2. internal boundaries: icefront(:).x_vards,icefront(:).y_vards
% and generate a mesh with internal constraints
nBoundary = numel(xouter);

node = [xouter(:) youter(:)];
edge = [[1:nBoundary]' [2:nBoundary 1]'];

save("dump.mat");

figure(111); plot(xouter,youter,'-k'); hold on;

if ExtrudeMesh

    fprintf("Assembling segments for internal constraints...");

    for ii=1:numel(years_to_include)
        icefront(ii).internalnodes=[];
        kk=1; ll=1; mm=1; firstnan=0;
        while kk<numel(icefront(ii).x_vards(:))
            if ~isnan(icefront(ii).mask(kk))
                firstnan=1;
                icefront(ii).segment(ll).x(mm) = icefront(ii).x_vards(kk);
                icefront(ii).segment(ll).y(mm) = icefront(ii).y_vards(kk);
                icefront(ii).segment(ll).ind(mm) = kk;
                mm=mm+1;
            else
                if firstnan
                    ll=ll+1;
                    mm=1;
                    firstnan=0;
                end
            end
            kk=kk+1;
        end

    end

    fprintf("Done.\n");
    fprintf("======================================\n");

    fprintf("Checking for intersecting segments...");

    pairs = nchoosek([1:numel(years_to_include)],2);

    for ii=1:size(pairs,1)  
        yr1 = pairs(ii,1); yr2 = pairs(ii,2);
        for l1=1:numel(icefront(yr1).segment)
            x1 = icefront(yr1).segment(l1).x(:)';
            y1 = icefront(yr1).segment(l1).y(:)';
            ind1 = icefront(yr1).segment(l1).ind(:)';
            %plot(x1,y1,'-b');
            for l2=1:numel(icefront(yr2).segment)
                x2 = icefront(yr2).segment(l2).x(:)';
                y2 = icefront(yr2).segment(l2).y(:)';
                ind2 = icefront(yr2).segment(l2).ind(:)';
                %plot(x2,y2,'-b');
                [xi,yi,segind] = polyxpoly(x1',y1',x2',y2');
                indi = [ind1(segind(:,1))' ind2(segind(:,2))'];

                plot(xi,yi,'or');
                
                % reassemble line segments with intersection points at
                % correct locations

                if ~isempty(xi) & numel(x1)>1 & numel(x2)>1
                    seg1=[]; seg2=[]; orig1=[]; orig2=[];
                    for i=1:numel(xi)
                        % check that xi is not too close to any of the
                        % already existing nodes
                        L1 = hypot(xi(i) - [x1(segind(i,1):segind(i,1)+1)],...
                            yi(i)-[y1(segind(i,1):segind(i,1)+1)]);
                        L2 = hypot(xi(i)-[x2(segind(i,2):segind(i,2)+1)],...
                            yi(i)-[y2(segind(i,2):segind(i,2)+1)]);
                        L = [L1(:); L2(:)];
                        [minL,indL] = min(L);
                        if minL<FdesiredEleSize(xi(i),yi(i))
                            if indL==1
                                xi(i) = x1(segind(i,1));
                                yi(i) = y1(segind(i,1));
                            elseif indL==2
                                xi(i) = x1(min(segind(i,1),numel(x1)));
                                yi(i) = y1(min(segind(i,1),numel(y1)));
                            elseif indL==3
                                xi(i) = x2(segind(i,2));
                                yi(i) = y2(segind(i,2));
                            elseif indL==4
                                xi(i) = x2(min(segind(i,2),numel(x2)));
                                yi(i) = y2(min(segind(i,2),numel(y2)));
                            end
                        end
                        if i==1
                            seg1(i).x = [x1(1:segind(i,1)),xi(i)];
                            seg1(i).y = [y1(1:segind(i,1)),yi(i)];
                            seg2(i).x = [x2(1:segind(i,2)),xi(i)];
                            seg2(i).y = [y2(1:segind(i,2)),yi(i)];  
                            orig1(i).x = [icefront(yr1).x_vards(1:indi(i,1)),xi(i)];
                            orig1(i).y = [icefront(yr1).y_vards(1:indi(i,1)),yi(i)];
                            orig2(i).x = [icefront(yr2).x_vards(1:indi(i,2)),xi(i)];
                            orig2(i).y = [icefront(yr2).y_vards(1:indi(i,2)),yi(i)];
                        else
                            seg1(i).x = [x1(segind(i-1,1)+1:segind(i,1)),xi(i)];
                            seg1(i).y = [y1(segind(i-1,1)+1:segind(i,1)),yi(i)];
                            seg2(i).x = [x2(segind(i-1,2)+1:segind(i,2)),xi(i)];
                            seg2(i).y = [y2(segind(i-1,2)+1:segind(i,2)),yi(i)];
                            orig1(i).x = [icefront(yr1).x_vards(indi(i-1,1)+1:indi(i,1)),xi(i)];
                            orig1(i).y = [icefront(yr1).y_vards(indi(i-1,1)+1:indi(i,1)),yi(i)];
                            orig2(i).x = [icefront(yr2).x_vards(indi(i-1,2)+1:indi(i,2)),xi(i)];
                            orig2(i).y = [icefront(yr2).y_vards(indi(i-1,2)+1:indi(i,2)),yi(i)];
                        end
                        plot([x1(segind(i,1)),xi(i),x1(segind(i,1)+1)],...
                                [y1(segind(i,1)),yi(i),y1(segind(i,1)+1)],'x-b');
                        plot([x2(segind(i,2)),xi(i),x2(segind(i,2)+1)],...
                            [y2(segind(i,2)),yi(i),y2(segind(i,2)+1)],'+-g');
                        % insert xi at correct place in original icefront 

                    end  
                    icefront(yr1).segment(l1).x = [seg1(:).x, x1(segind(end,1)+1:end)];
                    icefront(yr1).segment(l1).y = [seg1(:).y, y1(segind(end,1)+1:end)];
                    icefront(yr2).segment(l2).x = [seg2(:).x, x2(segind(end,2)+1:end)];
                    icefront(yr2).segment(l2).y = [seg2(:).y, y2(segind(end,2)+1:end)];

                    icefront(yr1).x_vards = [orig1(:).x, icefront(yr1).x_vards(indi(end,1)+1:end)];
                    icefront(yr1).y_vards = [orig1(:).y, icefront(yr1).y_vards(indi(end,1)+1:end)];
                    icefront(yr2).x_vards = [orig2(:).x, icefront(yr2).x_vards(indi(end,2)+1:end)];
                    icefront(yr2).y_vards = [orig2(:).y, icefront(yr2).y_vards(indi(end,2)+1:end)];

                    %plot(icefront(yr1).segment(l1).x,icefront(yr1).segment(l1).y,'-.m');
                    %plot(icefront(yr2).segment(l2).x,icefront(yr2).segment(l2).y,'-xg');
                end
            end
        end
    end

    fprintf("Done.\n");
    fprintf("======================================\n");

    fprintf("Restructuring data to appropriate mesh2d format...");

    to_remove = [-1765 -1760 996 997.2; ...
        -1725 -1724 984.8 985.4;...
        -1452 -1450 788 790;...
        -470.4 -470 1913.6 1914;...
        800 801 2064 2064.5]*1e3; %xmin xmax ymin ymax

    for ii=1:numel(years_to_include)

        for kk=1:numel(icefront(ii).segment)

            newsegment = [icefront(ii).segment(kk).x(:), icefront(ii).segment(kk).y(:)];
            
            % points to remove
            I = find(~inpoly(newsegment,[xouter(:) youter(:)]));
            for rr=1:size(to_remove,1)
                Itmp = find(icefront(ii).segment(kk).x(:)>to_remove(rr,1) & ...
                    icefront(ii).segment(kk).x(:)<to_remove(rr,2) & ...
                    icefront(ii).segment(kk).y(:)>to_remove(rr,3) & ...
                    icefront(ii).segment(kk).y(:)<to_remove(rr,4));
                I = [I(:); Itmp(:)];
            end
            newsegment(I,:)=[];

            n = size(newsegment,1);

            if n>1
                
                %plot(newsegment,'.--k');

                node = [node; newsegment];
                edge = [edge; nBoundary + [[1:n-1]' [2:n]']];
                nBoundary = nBoundary + n;     
                
            end
  
        end 

        % remove self-interesections and duplicate boundary coordinates
        tmp = [icefront(ii).x_vards(:) icefront(ii).y_vards(:)];
        [~,~,segtmp] = selfintersect(tmp(:,1),tmp(:,2));

        for kk=1:size(segtmp,1)
            if years_to_include(ii)~=2018 && kk~=1
                ind = [segtmp(kk,1)+1:segtmp(kk,2)-1];
                tmp(ind,:)=nan;
            end
        end
        for rr=1:size(to_remove,1)
            Itmp = find(tmp(:,1)>to_remove(rr,1) & ...
                tmp(:,1)<to_remove(rr,2) & ...
                tmp(:,2)>to_remove(rr,3) & ...
                tmp(:,2)<to_remove(rr,4));
            tmp(Itmp,:) = nan;
        end
        tmp(isnan(tmp(:,1)),:)=[];
        tmp = unique(tmp,"rows","stable");

        %D=pdist2(tmp,tmp);
        %D(D==0)=nan;
        %min(D(:))

        MeshBoundaryCoordinates.(['yr',num2str(years_to_include(ii))]) = tmp;

    end
    
    fprintf("Done.\n");
    fprintf("======================================\n");

end

save("dump2.mat");

% for ii=1:size(node,1)
%     if ii==1
%         [k,dist] = dsearchn(node(2:end,:),node(ii,:));
%     elseif ii==size(node,1)
%         [k,dist] = dsearchn([node(1:end-1,:)],node(ii,:));        
%     else
%         [k,dist] = dsearchn([node(1:ii-1,:);node(ii+1:end,:)],node(ii,:));
%     end
%     Lmin(ii) = dist;
% end
% Lmin(Lmin==0)=1; 
% figure; hold on; scatter(node(:,1),node(:,2),[],abs(Lmin(:)-FdesiredEleSize(node(:,1),node(:,2)))/500,'filled');
% ind_x = find(Lmin(:) < 0.9*FdesiredEleSize(node(:,1),node(:,2)));
% scatter(node(ind_x,1),node(ind_x,2),'LineWidth',1.5,'MarkerEdgeColor','b');
% caxis([-5 5]);
% hold on
% plot(node(:,1),node(:,2),'xk');
% colormap(jet)
% axis equal;
% axis tight;

part{1} = [1:numel(xouter)];

opts.kind = 'delfront'; %'delaunay'
opts.rho2 = 1;%1.025;
opts.ref1 = 'preserve';%'refine';%'preserve';
opts.siz1 = 1.333;%1.333
opts.disp = 1;
%hfun = 6e3;
hfun = @(X) FdesiredEleSize(X(:,1),X(:,2));

[vert,etri,...
    tria,tnum] = refine2(node,edge,part,opts,hfun);

[vert,etri,...
    tria,tnum] = smooth2(vert,etri,tria,tnum);

CtrlVar = Ua2D_DefaultParameters;
CtrlVar.MeshBoundaryCoordinates = [xouter(:) youter(:)];
CtrlVar.TriNodes=3;

[coordinates,connectivity] = ChangeElementType(vert,tria,CtrlVar.TriNodes);

if CtrlVar.sweep
    [coordinates,connectivity] = ElementSweep(coordinates,connectivity,CtrlVar.SweepAngle);
    [coordinates,connectivity] = NodalSweep(coordinates,connectivity,CtrlVar.SweepAngle);
end

if CtrlVar.CuthillMcKee
    M=connectivity2adjacency(connectivity);
    [coordinates,connectivity] = CuthillMcKeeFE(coordinates,connectivity,M);
end

connectivity = TestAndCorrectForInsideOutElements(CtrlVar,coordinates,connectivity);

MUA = CreateMUA(CtrlVar,connectivity,coordinates);

Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity);
Tlength=sqrt(2*Tarea) ;

figure; hold on;
%PlotMeshScalarVariable(CtrlVar,MUA,FdesiredEleSize(MUA.coordinates(:,1),MUA.coordinates(:,2))); hold on;
%caxis([0 5e3]);
%PlotMuaMesh(CtrlVar,MUA);
for ii=1:size(edge,1)
    plot(node(edge(ii,:),1)/1e3,node(edge(ii,:),2)/1e3,'-k','LineWidth',2);
end
xEle = Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1)/1e3);
yEle = Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2)/1e3);
I = find(Tlength<100);
hold on
plot(node(:,1)/1e3,node(:,2)/1e3,'xg');
plot(xEle(I),yEle(I),'or');

fname = "";
for ii=1:numel(years_to_include)
    fname = fname + "_" +string(years_to_include(ii));
end
fname = fname + "_meshmin" + string(meshmin) + "_meshmax" + string(meshmax) + ...
    "_extrudemesh" + string(ExtrudeMesh) + "_variableboundaryres" + ...
    string(VariableBoundaryResolution) ;

fname_MUA = "ANT_basemesh" + fname;
save(fname_MUA,"MUA","CtrlVar");

fname_Boundary = "ANT_meshboundarycoordinates" + fname;

MeshBoundaryCoordinates.outer = [MUA.Boundary.x(:) MUA.Boundary.y(:)];

save(fname_Boundary,"MeshBoundaryCoordinates");

figure; hold on;
patch('faces',tria(:,1:3),'vertices',vert, ...
    'facecolor','w', ...
    'edgecolor',[0.2 0.2 0.2]);
patch('faces',tria(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[0.1 0.1 0.1], ...
    'linewidth',1.5);

% for ii=1:length(out)
%    plot(out(ii).points(:,1),out(ii).points(:,2),'-b'); 
% end

for ii=1:numel(years)
    plot(icefront(ii).x_vards,icefront(ii).y_vards,'-*r');
end
axis equal;

fprintf("Done.\n");

warning('on','all');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ANYTHING BELOW HERE ARE OLD SCRIPTS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% now do spatially variable resampling
    % Peninsula, East Antarctica: 15km
    % Ross and FRIS: 5km
    % Amundsen and Bellingshausen, Amery, Totten: 2km
    
%     IF.tension=1;
%     IF.defaultds = 15e3;
%     %
%     IF.segments(1).name = 'Weddell Sea';
%     IF.segments(1).lonlim = [-70 -30];
%     IF.segments(1).latlim = [-80 -75];
%     IF.segments(1).ds = 5e3; 
%     %
%     IF.segments(2).name = 'Ross Sea';
%     IF.segments(2).lonlim = [150 210];
%     IF.segments(2).latlim = [-80 -75];
%     IF.segments(2).ds = 5e3; 
%     %
%     IF.segments(3).name = 'Amery';
%     IF.segments(3).lonlim = [65 80];
%     IF.segments(3).latlim = [-70 -65];
%     IF.segments(3).ds = 2e3; 
%     %
%     IF.segments(4).name = 'Totten';
%     IF.segments(4).lonlim = [110 120];
%     IF.segments(4).latlim = [-70 -65];
%     IF.segments(4).ds = 2e3; 
%     %
%     IF.segments(5).name = 'Amundsen';
%     IF.segments(5).lonlim = [-120 -80];
%     IF.segments(5).latlim = [-77.5 -70];
%     IF.segments(5).ds = 2e3;
% 
%     [xtmp,ytmp] = IFresample(xtmp,ytmp,IF);


%%
% step 1. create outlines from Bedmap2 (2003) and Bedmachine (2016)
% step 2. take the union of both shapes, which will be the maximum grid extent
% step 3. ice fronts from 2003 and 2016 will be imposed as internal
% boundaries in the mesh


froot_data = getenv("froot_data");
addpath(froot_data+"/Bedmap2/bedmap2_tif");

BoundaryResolution = 2e3;
SaveOutputs = 1;

%% Meshboundary coordinates Bedmap2 (2003)

[x,y,~] = bedmap2_data('surface','xy');
mask = bedmap2_data('icemask');

fprintf('Creating MeshBoundary coordinates Bedmap2...')
IceMask=(mask>10) ;  IceMask=double(IceMask);
DataResolution=x(2,1)-x(1,1); % the resolution of the data set
NN=min([round(BoundaryResolution/DataResolution) 1]) ;

fc=FindOrCreateFigure('contour') ;  
[M]=contour(x(1,1:NN:end),y(1:NN:end,1),IceMask(1:NN:end,1:NN:end),1) ; axis equal
hold on; plot(M(1,:), M(2, :), 'r.');

Boundary = BoundaryFromContour(M);

plot(Boundary(1,:),Boundary(2,:),'o-k')  ; axis equal

Boundary=Boundary';
fprintf('done.\n')

if SaveOutputs
    
    fprintf('Saving MeshBoundaryCoordinates Bedmap2. \n ')
    save('MeshBoundaryCoordinatesForAntarcticaBasedOnBedmap2_2003','Boundary');

    S=mapshape();
    S.Geometry='line';
    %S.BoundingBox=[min(Boundary(1,:)) min(Boundary(2,:));max(Boundary(1,:)) max(Boundary(2,:))];
    S.X=Boundary(:,1);
    S.Y=Boundary(:,2);

    S.Name='MeshBoundaryCoordinatesForAntarcticaBasedOnBedmap2_2003';
    shapewrite(S,'MeshBoundaryCoordinatesForAntarcticaBasedOnBedmap2_2003.shp');

end

%% Meshboundary coordinates Bedmachine 2016
% Find the boundary by extracting the 0.5 contour line of the IceMask Then extract the largest single contourline. Depending on the
% situation, this may or may not be what the user wants. But this appears a reasonable guess as to what most users might want most of the
% time. 

ncfile = froot_data+"/BedMachine_Antarctica/BedMachineAntarctica-v3.nc";

xdata = double(ncread(ncfile,'x'));
ydata = double(ncread(ncfile,'y'));
%xdata = unique(xdata);
%ydata = unique(ydata);

mask  = double(ncread(ncfile,'mask')');

fprintf('Creating MeshBoundary coordinates Bedmachine...')
IceMask=(mask~=0) ;  IceMask=double(IceMask);
DataResolution=xdata(2)-xdata(1); % the resolution of the data set is 500 m
NN=min([round(BoundaryResolution/DataResolution) 1]) ;

fc=FindOrCreateFigure('contour') ;  
[M]=contour(xdata(1:NN:end),ydata(1:NN:end),IceMask(1:NN:end,1:NN:end),1) ; axis equal
hold on; plot(M(1,:), M(2, :), 'm.');

Boundary = BoundaryFromContour(M);

plot(Boundary(1,:),Boundary(2,:),'o-b')  ; axis equal

Boundary=Boundary';
fprintf('done.\n')

if SaveOutputs
    
    fprintf('Saving MeshBoundaryCoordinates Bedmachinev3. \n ');
    save('MeshBoundaryCoordinatesForAntarcticaBasedOnBedmachinev3_2016','Boundary');

    S=mapshape();
    S.Geometry='line';
    %S.BoundingBox=[min(Boundary(1,:)) min(Boundary(2,:));max(Boundary(1,:)) max(Boundary(2,:))];
    S.X=Boundary(:,1);
    S.Y=Boundary(:,2);

    S.Name='MeshBoundaryCoordinatesForAntarcticaBasedOnBedmachinev3_2016';
    shapewrite(S,'MeshBoundaryCoordinatesForAntarcticaBasedOnBedmachinev3_2016.shp');

end




return
%%  Testing mesh boundary coordinates and creating  a new one with different spacing between points and some level of smoothing

CtrlVar.GLtension=1e-12; % tension of spline, 1: no smoothing; 0: straight line
CtrlVar.GLds=5e3 ; 

[xB,yB,nx,ny] = Smooth2dPos(Boundary(:,1),Boundary(:,2),CtrlVar);
MeshBoundaryCoordinates=[xB(:) yB(:)] ;
fc=FindOrCreateFigure('MeshBoundaryCoordinates') ;  
plot(MeshBoundaryCoordinates(:,1),MeshBoundaryCoordinates(:,2),'.-') ; axis equal
title("Example of a smoothed and resampled boundary")
end

function Boundary = BoundaryFromContour(M)
    % now find longest contourline
    level=0.5 ;  % this contour level must be in M
    I=find(M(1,:)==level) ; [~,J]=max(M(2,I)) ; 
    fprintf(' %i points in the longest contour line segment.\n',M(2,I(J)) );
    Boundary=M(:,I(J)+1:I(J)+M(2,I(J))) ; 
end

function [x,y] = IFresample(xinit,yinit,IF)
    [latinit,loninit] = psxy2ll(xinit,yinit,-71,0);
    % default
    CtrlVar.GLtension=IF.tension; % tension of spline, 1: no smoothing; 0: straight line
    CtrlVar.GLds=IF.defaultds ; 
    [xtmp,ytmp,~,~] = Smooth2dPos(xinit(:),yinit(:),CtrlVar);
    [lattmp,lontmp] = psxy2ll(xtmp,ytmp,-71,0);
    % refine regions
    for ii=1:numel(IF.segments)
        CtrlVar.GLds=IF.segments(ii).ds ; 
        lonlimits = [IF.segments(ii).lonlim(1), IF.segments(ii).lonlim(1),...
            IF.segments(ii).lonlim(2) IF.segments(ii).lonlim(2)];
        latlimits = [IF.segments(ii).latlim(1), IF.segments(ii).latlim(2),...
            IF.segments(ii).latlim(2) IF.segments(ii).latlim(1)];
        I = find(inpoly([loninit(:) latinit(:)],[lonlimits(:) latlimits(:)]));
        J = find(inpoly([lontmp(:) lattmp(:)],[lonlimits(:) latlimits(:)]));
        xtmp(J)=[]; ytmp(J)=[];
        [IF.segments(ii).xtmp,IF.segments(ii).ytmp,~,~] = Smooth2dPos(xinit(I(:)),yinit(I(:)),CtrlVar);
        xtmp = [xtmp(1:J(1)-1); IF.segments(ii).xtmp; xtmp(J(1):end)];
        ytmp = [ytmp(1:J(1)-1); IF.segments(ii).ytmp; ytmp(J(1):end)];
        [lattmp,lontmp] = psxy2ll(xtmp,ytmp,-71,0);
    end
    x = xtmp;
    y = ytmp;
end


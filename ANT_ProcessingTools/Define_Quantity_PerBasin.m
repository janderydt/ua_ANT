function [basin,ind] = Define_Quantity_PerBasin(x,y,B,doplots)

%% detect basin for each coordinate pair
if nargin<4
    doplots=0;
end

% 1. inpoly
xtmp = x;
ytmp = y;
for ii=1:numel(B.x)
    Bx = B.x{ii};
    By = B.y{ii};
    basin(ii).ind = [];
    Inan = [0; find(isnan(Bx)); numel(Bx)];
    for jj=1:numel(Inan)-2
        basin(ii).ind = [basin(ii).ind; find(inpoly2([xtmp(:) ytmp(:)],...
            [Bx(Inan(jj)+1:Inan(jj+1)-1) By(Inan(jj)+1:Inan(jj+1)-1)]))];
    end
    % avoid counting double
    xtmp(basin(ii).ind) = 1e12;
    ytmp(basin(ii).ind) = 1e12;
end

% 2. any remaining nodes are allocated to basin with nearest boundary
remaining = find(xtmp<1e12-eps);
PQ = [xtmp(remaining) ytmp(remaining)]; dist = [];
for ii=1:numel(B.x)
    P = [B.x{ii} B.y{ii}];
    [~,dist(ii,:)] = dsearchn(P,PQ);
end
[~,I] = min(dist,[],1);
   
% group indices
ind = zeros(numel(x),1)+nan;
for ii=1:numel(B.x)
    basin(ii).ind = [basin(ii).ind; remaining(find(I==ii))];
    ind(basin(ii).ind(:)) = ii;
    basin(ii).name = B.name{ii};
end

if doplots

    CM = parula(numel(B.x));
    I = randperm(numel(B.x));
    CM = CM(I,:);

    figure; hold on;
    for ii=1:numel(B.x)
        plot(B.x{ii},B.y{ii},'-k');
        plot(x(basin(ii).ind),y(basin(ii).ind),'o','markersize',1,'color',CM(ii,:));
    end

end
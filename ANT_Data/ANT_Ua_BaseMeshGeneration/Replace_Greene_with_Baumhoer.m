function [xnew,ynew] = Replace_Greene_with_Baumhoer(xold,yold,lims,year)

xmin = lims(1);
xmax = lims(2);
ymin = lims(3);
ymax = lims(4);

%% load Icelines data
datafolder = "/mnt/md0/Antarctic_datasets/IceLines/Data/";

if year == 2020
    files = ["2020noQ1_mean-Thwaites1.shp","2021noQ1_mean-Thwaites2.shp",...
       "2020noQ1_mean-Crosson.shp","2020noQ1_mean-PineIsland.shp"];
end

for ii=1:numel(files)
    tmp = shaperead(datafolder+files(ii));
    Icelines(ii).x = [];
    Icelines(ii).y = [];
    for jj=1:numel(tmp)
        Ind = find(~isnan(tmp(jj).X));
        Icelines(ii).x = [Icelines(ii).x(:); tmp(jj).X(Ind)'];
        Icelines(ii).y = [Icelines(ii).y(:); tmp(jj).Y(Ind)'];
    end
end

% combine Thwaites ice lines
Icelines(1).x = [Icelines(1).x ; Icelines(2).x];
Icelines(1).y = [Icelines(1).y ; Icelines(2).y];
Icelines(2).x = Icelines(3).x; Icelines(2).y = Icelines(3).y;
Icelines(3).x = Icelines(4).x; Icelines(3).y = Icelines(4).y;
Icelines = Icelines(1:3);
xtmp=xold;
ytmp=yold;
for ii=1:numel(Icelines)
    tmp = unique([Icelines(ii).x(:) Icelines(ii).y(:)],"rows","stable");
    Icelines(ii).x = tmp(:,1); Icelines(ii).y = tmp(:,2);
    [Icelines(ii).x,Icelines(ii).y,~] = Arrange2dPos(Icelines(ii).x,Icelines(ii).y);
    Ind = find(xmin<Icelines(ii).x & xmax>Icelines(ii).x &...
        ymin<Icelines(ii).y & ymax>Icelines(ii).y);
    Icelines(ii).x = Icelines(ii).x(Ind);
    Icelines(ii).y = Icelines(ii).y(Ind);
    %plot(Icelines(ii).x,Icelines(ii).y);
    % find points in [xold,yold] that are nearest to the start and end
    % point of [Icelines(ii).x Icelines(ii).y]
    [~,Ind_start]=min(hypot(xold-Icelines(ii).x(1),yold-Icelines(ii).y(1)));
    [~,Ind_end]=min(hypot(xold-Icelines(ii).x(end),yold-Icelines(ii).y(end)));

    if ii==2 || ii==3
        xold = [xold(1:Ind_end); flipdim(Icelines(ii).x,1); xold(Ind_start:end)];
        yold = [yold(1:Ind_end); flipdim(Icelines(ii).y,1); yold(Ind_start:end)];
    elseif ii==1
        xold = [xold(1:Ind_start); Icelines(ii).x; xold(Ind_end:end)];
        yold = [yold(1:Ind_start); Icelines(ii).y; yold(Ind_end:end)];
    end
    
end

xnew = xold;
ynew = yold;



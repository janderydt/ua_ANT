function [xtmp,ytmp] = Replace_Greene_with_Baumhoer(xold,yold,lims,year)

xmin = lims(1);
xmax = lims(2);
ymin = lims(1);
ymax = lims(2);

%% load Icelines data
datafolder = "/mnt/md0/Antarctic_datasets/IceLines/Data/";

if year == 2020
    files = ["2020noQ1_mean-Thwaites1.shp","2021noQ1_mean-Thwaites2.shp",...
        "2020noQ1_mean-Dotson.shp","2020noQ1_mean-Crosson.shp","2020noQ1_mean-PineIsland.shp"];
end

for ii=1:numel(files)
    tmp = shaperead(datafolder+files(ii));
    Icelines(ii).x = [];
    Icelines(ii).y = [];
    for jj=1:numel(tmp)
        Icelines(ii).x = [Icelines(ii).x; tmp(jj).X(:)];
        Icelines(ii).y = [Icelines(ii).y; tmp(jj).Y(:)];
    end
end



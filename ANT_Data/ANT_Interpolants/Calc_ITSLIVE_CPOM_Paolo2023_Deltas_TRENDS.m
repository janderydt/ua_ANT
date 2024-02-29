function [ds_trend,ds_significance] = Calc_ITSLIVE_CPOM_Paolo2023_Deltas_TRENDS(X,Y)

% X and Y are optinal 2d matrices (ndgrid structure) for which ds_trend and 
% ds_significance can be returned
% if X and Y are not provide, coordinates of the original ds files are
% assumed.
% linear trends and associated sifnificance of ds over the available 
% observational period will be returned and saved

froot_data = getenv("froot_data");
addpath(getenv("froot_tools"));

%% |||||||||||||||||||||| %%
%% Read original ds files %%
%% |||||||||||||||||||||| %%

years = [1995:2018];

if nargin==0
    load("ds_01-Jun-1995.mat","x","y");
    [X,Y] = ndgrid(x,y);
end

ds = zeros(13333,13333,numel(years));

for ii=1:numel(years)
    ds_tmp = load("ds_01-Jun-"+string(years(ii))+".mat","ds");
    ds(:,:,ii) = ds_tmp.ds;
    fprintf("Done %s\n",string(years(ii)));
end

[m,b] = trend(ds,years,3);  

time = permute(repmat(years',[1,size(m)]),[2 3 1]);
m_tmp = repmat(m,[1,1,numel(years)]);
b_tmp = repmat(b,[1,1,numel(years)]);
ds_pred = m_tmp.*time + b_tmp;
clear time m_tmp b_tmp;
SSE = sum((ds-ds_pred).^2,3);
clear ds_pred;
se = sqrt(SSE./(numel(years)-2));
clear SSE;
avtime = mean(years);
sb = se./sqrt(sum((years-avtime).^2));
t = m./sb;
df = numel(years)-2;
p = tpdf(t,df);

m(p>1e-2)=nan;

mask = 0*m+1;
mask(isnan(mask))=0;

save("ds_TREND_mask.mat","mask","x","y");

end
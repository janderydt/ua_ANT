function dqGL = Calc_Delta_qGL(cycle,basins_to_analyze)

if nargin==0
    UserVar.cycle = 2;
    basins_to_analyze = {'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbott
else
    UserVar.cycle = cycle;
end

% gather basin data
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveIceRisesAndIslands(B);
[~,BasinInd] = ismember(basins_to_analyze,B.name);

% load inversion data
load("inversiondata.mat");

m = [data(:).m];
n = [data(:).n];
gaA = [data(:).gaA];
gaC = [data(:).gaC];
gsA = [data(:).gsA];
gsC = [data(:).gsC];
expID = [data(:).InverseExpID];

for ii=1:numel(data)
    if UserVar.cycle==1
        ind_finished(ii) = ismember(data(ii).niter(UserVar.cycle),[5000,15000]);
        qGL_tmp = data(ii).qGL(:,UserVar.cycle);
    elseif UserVar.cycle==2
        if numel(data(ii).niter)>1
            ind_finished(ii) = ismember(data(ii).niter(UserVar.cycle),[6000,16000]);
            qGL_tmp = data(ii).qGL(:,UserVar.cycle);
        else
            ind_finished(ii) = 0;
            qGL_tmp = nan*ones(numel(B.name),1);
        end
    end   
    qGL(ii) = sum(qGL_tmp(BasinInd));
end

% find 2018 inversions
inv_2018=find(expID>5000 & expID<6000 & ind_finished==1);

% find 2000 inversions
ind_2000=find(expID>3000 & expID<4000 & ind_finished==1);

% for each 2018 inversion, find corresponding 2000 inversion
for ii=1:numel(inv_2018)
    ind_2018 = inv_2018(ii);
    metric = sqrt((m(ind_2018)-m(ind_2000)).^2+...
        (n(ind_2018)-n(ind_2000)).^2+...
        (gsA(ind_2018)-gsA(ind_2000)).^2+...
        (gsC(ind_2018)-gsC(ind_2000)).^2+...
        (gaA(ind_2018)-gaA(ind_2000)).^2+...
        (gaC(ind_2018)-gaC(ind_2000)).^2);
    [output,tmp] = min(metric);
    if output < eps
        ind_2000_matching(ii) = ind_2000(tmp);
    else
        ind_2000_matching(ii) = nan;
    end
end

inv_2018(isnan(ind_2000_matching))=[];
ind_2000_matching(isnan(ind_2000_matching))=[];
% calculate difference in grounding line flux between both
dqGL = qGL(inv_2018)-qGL(ind_2000_matching);

% plot distribution
edges = linspace(min(qGL),max(qGL),25);
[N_2018,~] = histcounts(qGL(inv_2018),edges,'Normalization','percentage');
[N_2000,~] =  histcounts(qGL(ind_2000_matching),edges,'Normalization','percentage');
figure; bar(0.5*(edges(1:end-1)+edges(2:end)),[N_2000(:)'; N_2018(:)']);

figure; histogram(dqGL,25,'Normalization','percentage');
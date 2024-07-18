function plot_ensemble_perturbation

variable_to_plot = 'm'; % n, m, gsA, gsC, gaA, gaC

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Diagnostic";
UserVar.Table = UserVar.home+"ANT_Diagnostic/"+["RunTable_ARCHER2_Diagnostic_2.csv"];
UserVar.idrange = [20000 29999];

addpath("/mnt/md0/Ua/cases/ANT/");

% load basins
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);

kk=0;

if exist("perturbationdata.mat","file")
    load("perturbationdata.mat");
else
    data=[];    
    perturbation_experiments_analyzed = [];
end
tmp = load("inversiondata.mat");
data_inverse = tmp.data;

for tt=1:numel(UserVar.Table)

    % read run table

    RunTable = ANT_ReadWritetable(UserVar,UserVar.Table(tt),[],'read');
    
    % ExpIDs
    ExpID = RunTable{:,"ExpID"};
    Ind = find(ExpID>=UserVar.idrange(tt,1) & ExpID<=UserVar.idrange(tt,2));
    % only keep experiments that have not been analyzed yet
    if ~isempty(data)
        Ind_ignore = ismember(Ind,perturbation_experiments_analyzed);
    else
        Ind_ignore = 0*Ind;
    end
    % only keep experiments that have finished
    Ind_finished = RunTable{Ind,"Finished"}==1;
    Ind = Ind(Ind_ignore==0 & Ind_finished==1);
    Comments = RunTable{Ind,"Comments"};   
    
    %% Gather data
    for ii=1:numel(Ind)

        perturbation_experiments_analyzed(end+1) = Ind(ii);
        InverseExpID = RunTable{Ind(ii),"InverseA"};
        InverseCycle = RunTable{Ind(ii),"InverseCycleA"};        
        Ind_inverse = find([data_inverse(:).InverseExpID]==InverseExpID);
        folder = UserVar.home+"/ANT_Diagnostic/cases/ANT_nsmbl_Diagnostic_"+ExpID(Ind(ii));
        outputfiles = dir(folder+"/ResultsFiles/*.mat");

        if ~isempty(outputfiles)
            outputfile = outputfiles(1).folder + "/" + outputfiles(1).name;
            expinfo = Comments{Ind(ii)};
            load(outputfile,"CtrlVar","F","MUA");             
            GL=FluxAcrossGroundingLine(CtrlVar,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho);
            qGL = sum(GL);    
            % store in data array
            if isempty(data)
                data_ind = 1;            
            else
                ExpID_list=[];
                for nn=1:numel(data)
                    ExpID_list(nn)=data(nn).Inverse.ExpID;
                end
                [~,data_ind] = ismember(InverseExpID,ExpID_list);
                if data_ind==0 % add new element to data structure
                    ndata = numel(data);
                    data_ind = ndata+1;
                end
            end

            data(data_ind).Inverse.ExpID = InverseExpID;
            data(data_ind).m = F.m(1);
            data(data_ind).n = F.n(1);
            data(data_ind).SlidingLaw = CtrlVar.SlidingLaw;
            data(data_ind).Inverse.gaC = data_inverse(Ind_inverse).gaC;
            data(data_ind).Inverse.gaA = data_inverse(Ind_inverse).gaA;
            data(data_ind).Inverse.gsC = data_inverse(Ind_inverse).gsC;
            data(data_ind).Inverse.gsA = data_inverse(Ind_inverse).gsA;
            data(data_ind).Inverse.misfit(InverseCycle) = data_inverse(Ind_inverse).misfit(InverseCycle);
            data(data_ind).startgeometry = data_inverse(Ind_inverse).startgeometry;

            %initialize structure
            initialize = 0;
            if ~isfield(data(data_ind),'Original')
                initialize = 1;
            else
                if ~isfield(data(data_ind).Original,'geometry')
                    initialize = 1;
                end
            end
            if initialize
                geomfields = {'Original','Calv','dhIS','dh','Calv_dh'};
                for ff=1:numel(geomfields)
                    data(data_ind).(geomfields{ff}).geometry=[];
                    data(data_ind).(geomfields{ff}).qGL=[];
                    data(data_ind).(geomfields{ff}).cycle=[];
                end
            end

            if contains(expinfo,"Original")
                year = RunTable{Ind(ii),"Calv"};              
                data(data_ind).Original.geometry(end+1) = year;
                data(data_ind).Original.qGL(end+1) = qGL;
                data(data_ind).Original.cycle(end+1) = InverseCycle;

            elseif contains(expinfo,"Ice front geometry")
                year = RunTable{Ind(ii),"Calv"};
                data(data_ind).Calv.geometry(end+1) = year;
                data(data_ind).Calv.qGL(end+1) = qGL;
                data(data_ind).Calv.cycle(end+1) = InverseCycle;
                
            elseif contains(expinfo,"Ice shelf thickness")
                year = RunTable{Ind(ii),"ISthick"};
                data(data_ind).dhIS.geometry(end+1) = year;
                data(data_ind).dhIS.qGL(end+1) = qGL;
                data(data_ind).dhIS.cycle(end+1) = InverseCycle;

            elseif contains(expinfo,"Ice thickness")
                year = RunTable{Ind(ii),"ISthick"};
                data(data_ind).dh.geometry(end+1) = year;
                data(data_ind).dh.qGL(end+1) = qGL;
                data(data_ind).dh.cycle(end+1) = InverseCycle;

            elseif contains(expinfo,"Ice front and thickness")
                year = RunTable{Ind(ii),"Calv"};
                data(data_ind).Calv_dh.geometry(end+1) = year;
                data(data_ind).Calv_dh.qGL(end+1) = qGL;
                data(data_ind).Calv_dh.cycle(end+1) = InverseCycle;

            else
                error("Unknown experiment info "+expinfo);
            end
        end             
        fprintf("Done %s out of %s.\n",string(ii),string(numel(Ind)));
    end
end

save("perturbationdata.mat","data","perturbation_experiments_analyzed");

%% gather data in userfriendly format 
for ii=1:numel(data)

    orig(ii,:)=data(ii).Original.qGL(:)'/1e12; % convert from kg/yr to Gt/yr
    calv(ii,:)=data(ii).Calv.qGL(:)'/1e12;
    dhIS(ii,:)=data(ii).dhIS.qGL(:)'/1e12;
    dh(ii,:)=data(ii).dh.qGL(:)'/1e12;
    calvdh(ii,:)=data(ii).Calv_dh.qGL(:)'/1e12;
    misfit(ii,:)=data(ii).Inverse.misfit(:)';
    gsA(ii) = [data(ii).Inverse.gsA];
    gsC(ii) = [data(ii).Inverse.gsC];
    gaA(ii) = [data(ii).Inverse.gaA];
    gaC(ii) = [data(ii).Inverse.gaC];

end

m = [data(:).m];
n = [data(:).n];

N = misfit-min(misfit);
N = N./max(N);
alphavalue = 1-N;
markersize = 50*alphavalue+eps;
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];

switch variable_to_plot
    case 'n'
        xdata = n;
    case 'm'
        xdata = m;
    case 'gsA'
        xdata = gsA;
    case 'gsC'
        xdata = gsC;
    case 'gaA'
        xdata = gaA;
    case 'gaC'
        xdata = gaC;
end

%% Plotting
figure; hold on;

s(1)=scatter(xdata,calv(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');
s(2)=scatter(xdata,dhIS(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');
s(3)=scatter(xdata,dh(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');
s(4)=scatter(xdata,calvdh(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');

s(5)=scatter(xdata,calv(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');
s(6)=scatter(xdata,dhIS(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');
s(7)=scatter(xdata,dh(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');
s(8)=scatter(xdata,calvdh(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');

colororder([newcolors;newcolors]);

xlim([floor(min(xdata)) ceil(max(xdata))]);

for ii=0:7
    s(ii+1).AlphaData = alphavalue(:,floor(ii/4)+1);
    s(ii+1).MarkerFaceAlpha='flat';
end

for ii=1:4
    s(ii)=plot(0,0,'o','color',newcolors(ii,:),'MarkerSize',10,'MarkerFaceColor',newcolors(ii,:));
end
s(5)=plot(0,0,'ok','markersize',10);
s(6)=plot(0,0,'sk','markersize',10);

legend(s(1:6),{'Calving','Ice Shelf thickness','Ice thickness','Calving + Ice thickness','no spinup','with spinup'},...
    'NumColumns',2,'Location','northwest');

xlabel(variable_to_plot); ylabel('\Delta q_{GL} [Gt/yr]');
grid on;
box on;



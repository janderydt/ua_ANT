function Calc_PerturbationResults_Ensemble

addpath(getenv("froot_ua")+"cases/ANT");

basins_to_analyze = {'A-Ap',...  % Queen Maud Land
    'Ap-B',... % Enderby Land
    'B-C',...  % Amery
    'C-Cp',... % 
    'Cp-D',... % Totton/Wilkes Land
    'D-Dp',... % George V Coast
    'Dp-E',... % Victoria Land
    'E-Ep',... % Ross west
    'Ep-F',... % Ross east
    'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp',... % Abbot
    'Hp-I',... % English Coast
    'I-Ipp',... % Northern Peninsula
    'Ipp-J',... % Eastern Peninsula
    'J-Jpp',... % Ronne 
    'Jpp-K',... % Filchner
    'K-A'}; % Caird Coast

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Diagnostic";
UserVar.Table = UserVar.home+"ANT_Diagnostic/"+["RunTable_ARCHER2_Diagnostic_2.csv"];
UserVar.idrange = [20000 29999];

% load basins
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);
Ind = contains(B.name,basins_to_analyze);
Bfields= fields(B);
for ii=1:numel(Bfields)
    if ~contains(Bfields{ii},'note')
        tmp = B.(Bfields{ii});
        B.(Bfields{ii}) = tmp(Ind);
    end
end

% load mesh for interpolation of speed
% tmp = load(UserVar.home+"ANT_Data/ANT_Ua_BaseMeshGeneration/ANT_basemesh_2018_meshmin5000_meshmax100000_extrudemesh0_variableboundaryres1.mat");
% MUA_coarse = tmp.MUA;
% 
% % identify basin id of each MUA node
% [MUA_coarse.basins,~] = Define_Quantity_PerBasin(MUA_coarse.coordinates(:,1),MUA_coarse.coordinates(:,2),B,0);

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
        Ind_ignore = ismember(ExpID(Ind),perturbation_experiments_analyzed);
    else
        Ind_ignore = 0*Ind;
    end
    % only keep experiments that have finished
    Ind_finished = RunTable{Ind,"Finished"}==1;
    Ind = Ind(Ind_ignore==0 & Ind_finished==1);
    Comments = RunTable{Ind,"Comments"};
    
    %% Gather data
    for ii=1:numel(Ind)

        perturbation_experiments_analyzed(end+1) = ExpID(Ind(ii));
        InverseExpID = RunTable{Ind(ii),"InverseA"};
        InverseCycle = RunTable{Ind(ii),"InverseCycleA"};        
        Ind_inverse = find([data_inverse(:).InverseExpID]==InverseExpID);
        folder = UserVar.home+"/ANT_Diagnostic/cases/ANT_nsmbl_Diagnostic_"+ExpID(Ind(ii));
        outputfiles = dir(folder+"/ResultsFiles/*.mat");

        if ~isempty(outputfiles)

            outputfile = outputfiles(1).folder + "/" + outputfiles(1).name;
            expinfo = Comments{ii};
            load(outputfile,"CtrlVar","F","MUA");      

            %GL=FluxAcrossGroundingLine(CtrlVar,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho);
            %qGL = sum(GL);
            % Obtain Ua fluxes across the grounding line (qGL) into floating areas
            [B,~] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVar);
            
            % Obtain Ua fluxes across the open boundary (qOB) into the ocean
            B = Calc_UaOBFlux_PerBasin(MUA,F,F.GF,B,CtrlVar); 
    
            % Sum values of SMB, qGL and qOB for each basin
            for bb=1:numel(B.x) 
                B.qGL_tot{bb} = sum(B.qGL{bb},'omitmissing')/1e12; 
                B.qOB_tot{bb} = sum(B.qOB{bb},'omitmissing')/1e12; 
                B.qtot{bb} = B.qGL_tot{bb}+B.qOB_tot{bb}; 
            end
          
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
                    for bb=1:numel(basins_to_analyze)
                        basin = char(erase(basins_to_analyze(bb),'-'));
                        data(data_ind).(geomfields{ff}).qGL.(basin)=[];
                        data(data_ind).(geomfields{ff}).qOB.(basin)=[];
                    end
                    data(data_ind).(geomfields{ff}).cycle=[];
                    data(data_ind).(geomfields{ff}).ExpID=[];
                    data(data_ind).(geomfields{ff}).speed=[];
                end
            end

            if contains(expinfo,"Original")
                year = RunTable{Ind(ii),"Calv"};   
                fieldname = 'Original';
            elseif contains(expinfo,"Ice front geometry")
                year = RunTable{Ind(ii),"Calv"};
                fieldname = 'Calv';
            elseif contains(expinfo,"Ice shelf thickness")
                year = RunTable{Ind(ii),"ISthick"};
                fieldname = 'dhIS';
            elseif contains(expinfo,"Ice thickness")
                year = RunTable{Ind(ii),"ISthick"};
                fieldname = 'dh';
            elseif contains(expinfo,"Ice front and thickness")
                year = RunTable{Ind(ii),"Calv"};
                fieldname = 'Calv_dh';
            else
                error("Unknown experiment info "+expinfo);
            end

            % save data
            data(data_ind).(fieldname).geometry(end+1) = year;
            for bb=1:numel(basins_to_analyze)
                basin = char(erase(basins_to_analyze(bb),'-'));
                data(data_ind).(fieldname).qGL.(basin)(end+1) = cell2mat(B.qGL_tot(bb));
                data(data_ind).(fieldname).qOB.(basin)(end+1) = cell2mat(B.qOB_tot(bb));               
            end           
            data(data_ind).(fieldname).cycle(end+1) = InverseCycle;
            data(data_ind).(fieldname).ExpID(end+1) = ExpID(Ind(ii));
            
            % Interpolate speed to relevant grid
            % Ind_nan = find(isnan(MUA.Boundary.x));
            % if ~isempty(Ind_nan)
            %     MUA.Boundary.x = MUA.Boundary.x(1:Ind_nan(1)-1);
            %     MUA.Boundary.y = MUA.Boundary.y(1:Ind_nan(1)-1);
            % end
            % Ind_out = find(~inpoly2(MUA_coarse.coordinates,[MUA.Boundary.x MUA.Boundary.y]));
            % 
            % if ismember(fieldname,["Original","dhIS","dh"])
            % 
            % 
            % 
            %     if isempty(Fspeed_2000)
            %         Fspeed_2000 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb),"natural");
            %     else
            %         Fspeed_2000.Values = hypot(F.ub,F.vb);
            %     end 
            %     speed = Fspeed_2000(MUA_coarse.coordinates(:,1),MUA_coarse.coordinates(:,2));
            % 
            % else
            %     if isempty(Fspeed_2018)
            %         Fspeed_2018 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),hypot(F.ub,F.vb),"natural");
            %     else
            %         Fspeed_2018.Values = hypot(F.ub,F.vb);
            %     end 
            %     speed = Fspeed_2018(MUA_coarse.coordinates(:,1),MUA_coarse.coordinates(:,2));
            % end
            % 
            % speed(Ind_out) = NaN;
            data(data_ind).(fieldname).speed(:,end+1) = hypot(F.ub(:),F.vb(:));

        end             
        fprintf("Done %s out of %s.\n",string(ii),string(numel(Ind)));
    end
end

save("perturbationdata.mat","data","perturbation_experiments_analyzed","-v7.3");
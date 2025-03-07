function Calc_InverseResults_Ensemble

addpath("/mnt/md0/Ua/cases/ANT/");
addpath(getenv("froot_tools"));

only_finished=1;

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
%UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_"+string([2 4 5 8])+".csv";
%UserVar.idrange = [3000,3999;5000,5999;6000,6999;9000,9999];
%UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_10-10-2024_"+string([10 11 12 13])+".csv";
%UserVar.idrange = [10000,10999; 11000, 11999; 12000, 12999; 13000, 13999];

%UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_"+string([3 6 9])+".csv";
%UserVar.idrange = [3000,3999; 6000, 6999; 9000, 9999];

UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_17-02-2025_"+string([24 25 26 27])+".csv";
UserVar.idrange = [24000 24999; 25000 25999; 26000 26999; 27000 27999];

inversiondata_filename = "inversiondata_AMUND_Umbi.mat";

if exist(inversiondata_filename,"file")
    load(inversiondata_filename);
else
    data=[];    
    inverse_experiments_analyzed = [];
end

%% ALL
% basins_to_analyze = {'A-Ap',...  % Queen Maud Land
%     'Ap-B',... % Enderby Land
%     'B-C',...  % Amery
%     'C-Cp',... % 
%     'Cp-D',... % Totton/Wilkes Land
%     'D-Dp',... % George V Coast
%     'Dp-E',... % Victoria Land
%     'E-Ep',... % Ross west
%     'Ep-F',... % Ross east
%     'F-G',...  % Getz
%     'G-H',...  % PIG, Thwaites
%     'H-Hp',... % Abbot
%     'Hp-I',... % English Coast
%     'I-Ipp',... % Northern Peninsula
%     'Ipp-J',... % Eastern Peninsula
%     'J-Jpp',... % Ronne 
%     'Jpp-K',... % Filchner
%     'K-A'}; % Caird Coast
%% AMUND
basins_to_analyze = {'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbott

filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveIceRisesAndIslands(B);
Ind = contains(B.name,basins_to_analyze);
Bfields= fields(B);
for ii=1:numel(Bfields)
    if ~contains(Bfields{ii},'note')
        tmp = B.(Bfields{ii});
        B.(Bfields{ii}) = tmp(Ind);
    end
end

for tt=1:numel(UserVar.Table)

    %% read run table
    RunTable = ANT_ReadWritetable(UserVar,UserVar.Table(tt),[],'read');
    
    %% ExpIDs
    ExpID = RunTable{:,"ExpID"};
    Ind = find(ExpID>=UserVar.idrange(tt,1) & ExpID<=UserVar.idrange(tt,2));
    % only keep experiments that have not been analyzed yet
    if ~isempty(data)
        Ind_ignore = ismember(ExpID(Ind),inverse_experiments_analyzed);
    else
        Ind_ignore = 0*Ind;
    end
    % only keep experiments that have finished
    Ind_finished = RunTable{Ind,"Finished"}==1;
    if only_finished
        Ind = Ind(Ind_ignore==0 & Ind_finished==1);
    else
        Ind = Ind(Ind_ignore==0);
    end

    %% Gather data
    for ii=1:numel(Ind)

        UserVar.domain = RunTable{Ind(ii),"Domain"};
        folder = UserVar.home+"/ANT_"+UserVar.type+"/cases/"+UserVar.domain+"_Inverse_"+ExpID(Ind(ii));
        
        % store in data array
        if isempty(data)
            data_ind = 1;            
        else
            [~,data_ind] = ismember(Ind(ii),[data(:).InverseExpID]);
            if data_ind==0 % add new element to data structure
                ndata = numel(data);
                data_ind = ndata+1;
            end
        end

        for cc=1:2

            restartfile = folder+"/"+UserVar.domain+"_"+UserVar.type+"_"+ExpID(Ind(ii))+"-RestartFile_InverseCycle"+...
                string(cc)+".mat";

            if exist(restartfile,"file")

                inverse_experiments_analyzed(end+1) = ExpID(Ind(ii));

                try   

                    load(restartfile,"UserVarInRestartFile","CtrlVarInRestartFile","F","MUA","InvFinalValues");
                               
                    data(data_ind).InverseExpID = ExpID(Ind(ii));
                    data(data_ind).cycle(cc) = cc;
                    data(data_ind).m = F.m(1);
                    data(data_ind).n = F.n(1);
                    data(data_ind).SlidingLaw = CtrlVarInRestartFile.SlidingLaw;
                    data(data_ind).gaA = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.ga;
                    data(data_ind).gaC = CtrlVarInRestartFile.Inverse.Regularize.logC.ga;
                    data(data_ind).gsA = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.gs;
                    data(data_ind).gsC = CtrlVarInRestartFile.Inverse.Regularize.logC.gs;
                    if isfield(UserVarInRestartFile.Inverse,"dhdt_err")
                        data(data_ind).dhdt_err = UserVarInRestartFile.Inverse.dhdt_err;
                    else
                        data(data_ind).dhdt_err = nan;
                    end
                    data(data_ind).startgeometry = RunTable{Ind(ii),"startGeometry"};
    
                    if ExpID(Ind(ii))<10000 && cc>1
                        
                        nvec_basins = numel(B);
                        nvec_nodes = MUA.Nnodes;
    
                        data(data_ind).niter(cc) = nan;
    
                        data(data_ind).qGL(:,cc) = nan*ones(nvec_basins,1);
                        data(data_ind).qOB(:,cc) = nan*ones(nvec_basins,1);
                        data(data_ind).TotalBalancedMelt(:,cc) = nan*ones(nvec_basins,1);
    
                        data(data_ind).BalancedMeltMap(:,cc) = nan*ones(nvec_nodes,1);
        
                        data(data_ind).misfit(cc) = nan;
                        data(data_ind).regularization(cc) = nan;
    
                    else
    

                        [B,~] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
                        B = Calc_UaOBFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
                        % Sum values of SMB, qGL and qOB for each basin
                        for bb=1:numel(B.x) 
                            B.qGL_tot{bb} = sum(B.qGL{bb},'omitmissing')/1e12;   
                            B.qOB_tot{bb} = sum(B.qOB{bb},'omitmissing')/1e12;   
                        end
                        qGL = cell2mat(B.qGL_tot);
                        qOB = cell2mat(B.qOB_tot);
                        % 
                        % 
                        ab = CalcIceShelfMeltRates(CtrlVarInRestartFile,MUA,F.ub,F.vb,F.s,F.b,F.S,F.B,F.rho,F.rhow,0*F.ub,F.as,0*F.ub);
                        x = MUA.coordinates(:,1);
                        y = MUA.coordinates(:,2);
                        for bb=1:numel(B.x)
                            xB = B.x{bb};
                            yB = B.y{bb};
                            ab_tot = 0;
                            Indnan = [0; find(isnan(xB)); numel(xB)+1];
                            for nn=1:numel(Indnan)-1
                                if ~isempty([Indnan(nn)+1:Indnan(nn+1)-1])
                                    Indpoly = find(inpoly([x y],[xB(Indnan(nn)+1:Indnan(nn+1)-1) yB(Indnan(nn)+1:Indnan(nn+1)-1)]));
                                    ab_tot = ab_tot + sum(ab(Indpoly));
                                end
                            end
                            BalancedMelt(bb) = ab_tot;
                        end
                        %
                        %
                        data(data_ind).niter(cc) = UserVarInRestartFile.Inverse.IterationsDone;
        
                        data(data_ind).qGL(:,cc) = qGL(:);
                        data(data_ind).qOB(:,cc) = qOB(:);
                        data(data_ind).TotalBalancedMelt(:,cc) = BalancedMelt(:);
                        data(data_ind).BalancedMeltMap(:,cc) = ab(:);
        
                        data(data_ind).misfit(cc) = InvFinalValues.I;
                        data(data_ind).regularization(cc) = InvFinalValues.R;
    
                    end
    
                    % Obtain Ua fluxes across the grounding line (qGL) into floating areas
                    %[B,GL] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
                    % qGL(ii) = 0;
                    % for jj=1:numel(GL)
                    %     qGL(ii) = qGL(ii)+sum(GL(jj).qGL);
                    % end            
                     % calculated as 
    
                    % store in data array

                catch

                    warning("Unable to load file "+restartfile+".");

                end

            end

        end

        fprintf("done %s our of %s.\n",string(ii),string(numel(Ind)));

    end

    if ~isempty(Ind)
        GF = F.GF;
        inverse_experiments_analyzed = unique(inverse_experiments_analyzed,"stable");
        save(inversiondata_filename,"data","inverse_experiments_analyzed","MUA","GF","-v7.3");
    end
        
end







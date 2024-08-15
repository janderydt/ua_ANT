function Calc_InverseResults_Ensemble

only_finished=0;

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_"+string([2 4 5 8])+".csv";
UserVar.idrange = [3000,3999;5000,5999;6000,6999;9000,9999];

addpath("/mnt/md0/Ua/cases/ANT/");

if exist("inversiondata.mat","file")
    load("inversiondata.mat");
else
    data=[];    
    inverse_experiments_analyzed = [];
end

filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);

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
    % only keep experiments that have finished?
    Ind_finished = RunTable{Ind,"Finished"}==1;
    if only_finished
        Ind = Ind(Ind_ignore==0 & Ind_finished==1);
    else
        Ind = Ind(Ind_ignore==0);
    end

    %% Gather data
    for ii=1:numel(Ind)

        inverse_experiments_analyzed(end+1) = ExpID(Ind(ii));
        folder = UserVar.home+"/ANT_"+UserVar.type+"/cases/ANT_nsmbl_Inverse_"+ExpID(Ind(ii));
        
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
            restartfile = folder+"/ANT_nsmbl_Inverse_"+ExpID(Ind(ii))+"-RestartFile_InverseCycle"+...
                string(cc)+".mat";

            if exist(restartfile,"file")

                load(restartfile,"UserVarInRestartFile","CtrlVarInRestartFile","F","MUA","InvFinalValues");   
                
                [B,~] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
                % Sum values of SMB, qGL and qOB for each basin
                for bb=1:numel(B.x) 
                    B.qGL_tot{bb} = sum(B.qGL{bb},'omitmissing')/1e12;   
                end
                qGL = cell2mat(B.qGL_tot);

                data(data_ind).InverseExpID = ExpID(Ind(ii));
                data(data_ind).cycle(cc) = cc;
                data(data_ind).m = F.m(1);
                data(data_ind).n = F.n(1);
                data(data_ind).SlidingLaw = CtrlVarInRestartFile.SlidingLaw;
                data(data_ind).gaA = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.ga;
                data(data_ind).gaC = CtrlVarInRestartFile.Inverse.Regularize.logC.ga;
                data(data_ind).gsA = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.gs;
                data(data_ind).gsC = CtrlVarInRestartFile.Inverse.Regularize.logC.gs;
                data(data_ind).startgeometry = RunTable{Ind(ii),"startGeometry"};
                data(data_ind).niter(cc) = UserVarInRestartFile.Inverse.IterationsDone;
                data(data_ind).misfit(cc) = InvFinalValues.I;
                data(data_ind).qGL(:,cc) = qGL(:);

                
                % Obtain Ua fluxes across the grounding line (qGL) into floating areas
                %[B,GL] = Calc_UaGLFlux_PerBasin(MUA,F,F.GF,B,CtrlVarInRestartFile);
                % qGL(ii) = 0;
                % for jj=1:numel(GL)
                %     qGL(ii) = qGL(ii)+sum(GL(jj).qGL);
                % end            
                 % calculated as 

                % store in data array
            end
        end

        fprintf("done %s our of %s.\n",string(ii),string(numel(Ind)));
        
    end
end

save("inversiondata.mat","data","inverse_experiments_analyzed");







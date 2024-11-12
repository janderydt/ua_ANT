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

UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_"+string([3 6 9])+".csv";
UserVar.idrange = [3000,3999; 6000, 6999; 9000, 9999];

UserVar.Table = [UserVar.Table UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_08-10-2024_"+string([14 15 16 17])+".csv"];
UserVar.idrange = [UserVar.idrange; 14000,14999; 15000, 15999; 16000, 16999; 17000, 17999];

inversiondata_filename = "inversiondata_Weertman.mat";

if exist(inversiondata_filename,"file")
    load(inversiondata_filename);
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

                inverse_experiments_analyzed(end+1) = ExpID(Ind(ii));

                load(restartfile,"UserVarInRestartFile","CtrlVarInRestartFile","F","MUA","InvFinalValues");
                
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
                data(data_ind).niter(cc) = UserVarInRestartFile.Inverse.IterationsDone;

                data(data_ind).qGL(:,cc) = qGL(:);
                data(data_ind).qOB(:,cc) = qOB(:);
                data(data_ind).TotalBalancedMelt(:,cc) = BalancedMelt(:);
                data(data_ind).BalancedMeltMap(:,cc) = ab(:);

                data(data_ind).misfit(cc) = InvFinalValues.I;
                data(data_ind).regularization(cc) = InvFinalValues.R;

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

GF = F.GF;
save(inversiondata_filename,"data","inverse_experiments_analyzed","MUA","GF","-v7.3");







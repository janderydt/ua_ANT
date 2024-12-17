function Calc_dhdt_Inverse

persistent MUA BCs

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

UserVar.Table = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_08-10-2024_"+string([14 15 16 17])+".csv";
UserVar.idrange = [14000,14999; 15000, 15999; 16000, 16999; 17000, 17999];

dhdt_filename = "dhdt_Weertman.mat";

if exist(dhdt_filename,"file")
    load(dhdt_filename);
end

filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveIceRisesAndIslands(B);

kk=1;

for tt=1:numel(UserVar.Table)

    %% read run table
    RunTable = ANT_ReadWritetable(UserVar,UserVar.Table(tt),[],'read');
    
    %% ExpIDs
    ExpID = RunTable{:,"ExpID"};
    dhdt_err_table = RunTable{:,"dhdt_err"};

    Ind = find(ExpID>=UserVar.idrange(tt,1) & ExpID<=UserVar.idrange(tt,2));
    
    % only keep experiments that have finished
    Ind = find(RunTable{Ind,"Finished"}==1);

    %% Gather data
    for ii=1:numel(Ind)

        folder = UserVar.home+"/ANT_"+UserVar.type+"/cases/ANT_nsmbl_Inverse_"+ExpID(Ind(ii));

        for cc=1
            
            restartfile = folder+"/ANT_nsmbl_Inverse_"+ExpID(Ind(ii))+"-RestartFile_InverseCycle"+...
                string(cc)+".mat";

            if exist(restartfile,"file")

                load(restartfile,"F","UserVarInRestartFile","CtrlVarInRestartFile");

                if isempty(MUA)
                    load(restartfile,"MUA","BCs");
                    MUA=UpdateMUA(CtrlVarInRestartFile,MUA);
                    x = MUA.coordinates(:,1);
                    y = MUA.coordinates(:,2);
                    PIG_interestregion = find(x>-1.7e6 & x<-1.5e6 & y>-3e5 & y<-1e5);
                end

                [~,dhdt]=dhdtExplicitSUPG(UserVarInRestartFile,CtrlVarInRestartFile,MUA,F,BCs);
                dhdt = dhdt.*F.GF.node; dhdt(dhdt==0)=nan;

                %% PIG
                dhdt_PIG(kk) = mean(dhdt(PIG_interestregion),"omitmissing");
                dhdt_err(kk) = dhdt_err_table(Ind(ii));
                kk=kk+1;

                fprintf("done "+num2str(ii)+" out of "+num2str(numel(Ind))+".\n");

            end
        end

    end

end

save("dhdt_PIG","dhdt_err","dhdt_PIG");

[Y,~]=discretize(dhdt_err,15);
boxplot(dhdt_PIG,Y);
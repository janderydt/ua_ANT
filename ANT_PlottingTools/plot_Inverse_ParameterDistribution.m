function plot_Inverse_ParameterDistribution

addpath("..");

RunTable = "RunTable_ARCHER2_"+string([2 5 8])+".csv";

gaA=[]; gaC=[]; gsA=[]; gsC=[]; m=[]; n=[]; finished=[]; error=[];
Ind_finished=[]; Ind_error=[];

for tt=1:numel(RunTable)

    UserVar.home = pwd; 
    UserVar.type = "Inverse";
    UserVar.Table = UserVar.home+"/../ANT_"+UserVar.type+"/"+RunTable(tt);

    RunTable_tmp = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

    gaA = [gaA; RunTable_tmp{:,'gaA'}];
    gaC = [gaC; RunTable_tmp{:,'gaC'}];
    gsA = [gsA; RunTable_tmp{:,'gsA'}];
    gsC = [gsC; RunTable_tmp{:,'gsC'}];
    m = [m; RunTable_tmp{:,'m'}];
    n = [n; RunTable_tmp{:,'n'}];
    finished_tmp = RunTable_tmp{:,'Finished'};
    finished = [finished; finished_tmp];
    error_tmp = RunTable_tmp{:,'Error'};
    error = [error; error_tmp];

    Ind_finished_tmp = find(finished_tmp==1);
    Ind_finished = [Ind_finished; Ind_finished_tmp];
    Ind_error_tmp = find(error_tmp==1);
    Ind_error = [Ind_error; Ind_error_tmp];

    fprintf("%s: Finished %s out of %s. %s ended with an error message. %s still running.\n",...
        RunTable(tt),...
        string(numel(Ind_finished_tmp)),string(numel(finished_tmp)),...
        string(numel(Ind_error_tmp)),...
        string(numel(finished_tmp)-numel(Ind_finished_tmp)-numel(Ind_error_tmp)));


end

Ind_finished = find(finished==1);
Ind_error = find(error==1);

fprintf("TOTAL: Finished %s out of %s. %s ended with an error message. %s still running.\n",...
    string(numel(Ind_finished)),string(numel(finished)),string(numel(Ind_error)),...
    string(numel(finished)-numel(Ind_error)-numel(Ind_finished)));

figure;  hold on;
g(1)=histogram(log10(gaA),linspace(min(log10(gaA)),max(log10(gaA)),50));
g(2)=histogram(log10(gaA(Ind_finished)),linspace(min(log10(gaA)),max(log10(gaA)),50));
g(3)=histogram(log10(gaA(Ind_error)),linspace(min(log10(gaA)),max(log10(gaA)),50));
title("gaA (min: "+string(min(gaA))+", max: "+string(max(gaA))+")");
xlabel('log10(gaA)'); ylabel('count');
grid on; box on;
legend(g(:),["all","finished","error"],'Location','northwest');

figure;  hold on;
g(1)=histogram(log10(gaC),linspace(min(log10(gaC)),max(log10(gaC)),50));
g(2)=histogram(log10(gaC(Ind_finished)),linspace(log10(min(gaC)),max(log10(gaC)),50));
g(3)=histogram(log10(gaC(Ind_error)),linspace(min(log10(gaC)),max(log10(gaC)),50));
title("gaC (min: "+string(min(gaC))+", max: "+string(max(gaC))+")");
xlabel('log10(gaC)'); ylabel('count');
grid on; box on;
legend(g(:),["all","finished","error"],'Location','northwest');

figure;  hold on;
g(1)=histogram(log10(gsA),linspace(min(log10(gsA)),max(log10(gsA)),50));
g(2)=histogram(log10(gsA(Ind_finished)),linspace(min(log10(gsA)),max(log10(gsA)),50));
g(3)=histogram(log10(gsA(Ind_error)),linspace(min(log10(gsA)),max(log10(gsA)),50));
title("gsA (min: "+string(min(gsA))+", max: "+string(max(gsA))+")");
xlabel('log10(gsA)'); ylabel('count');
grid on; box on;
legend(g(:),["all","finished","error"],'Location','northwest');

figure;  hold on;
g(1)=histogram(log10(gsC),linspace(min(log10(gsC)),max(log10(gsC)),50));
g(2)=histogram(log10(gsC(Ind_finished)),linspace(min(log10(gsC)),max(log10(gsC)),50));
g(3)=histogram(log10(gsC(Ind_error)),linspace(min(log10(gsC)),max(log10(gsC)),50));
title("gsC (min: "+string(min(gsC))+", max: "+string(max(gsC))+")");
xlabel('log10(gsC)'); ylabel('count');
grid on; box on;
legend(g(:),["all","finished","error"],'Location','northwest');

figure;  hold on;
g(1)=histogram(m,linspace(min(m),max(m),50));
g(2)=histogram(m(Ind_finished),linspace(min(m),max(m),50));
g(3)=histogram(m(Ind_error),linspace(min(m),max(m),50));
title("m (min: "+string(min(m))+", max: "+string(max(m))+")");
xlabel('m'); ylabel('count');
grid on; box on;
legend(g(:),["all","finished","error"],'Location','northwest');

figure;  hold on;
g(1)=histogram(n,linspace(min(n),max(n),50));
g(2)=histogram(n(Ind_finished),linspace(min(n),max(n),50));
g(3)=histogram(n(Ind_error),linspace(min(n),max(n),50));
title("n (min: "+string(min(n))+", max: "+string(max(n))+")");
xlabel('n'); ylabel('count');
grid on; box on;
legend(g(:),["all","finished","error"],'Location','northwest');

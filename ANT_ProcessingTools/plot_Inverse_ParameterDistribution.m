function plot_Inverse_ParameterDistribution

addpath("..");

%RunTable = "RunTable_ARCHER2_"+string([2 5 8])+".csv";
%RunTable = "RunTable_ARCHER2_"+string([6])+".csv";
%RunTable = "RunTable_ARCHER2_10-10-2024_"+string([10 11 12 13])+".csv";
RunTable = "RunTable_ARCHER2_08-10-2024_"+string([14 15 16 17])+".csv";

gaA=[]; gaC=[]; gsA=[]; gsC=[]; m=[]; n=[]; dhdt=[];
finished=[]; error=[]; running=[]; ExpID=[];
Ind_finished=[]; Ind_error=[]; Ind_running=[];

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
    dhdt = [dhdt; RunTable_tmp{:,'dhdt_err'}];    
    ExpID = [ExpID; RunTable_tmp{:,'ExpID'}];

    finished_tmp = RunTable_tmp{:,'Finished'};
    finished = [finished; finished_tmp];
    error_tmp = RunTable_tmp{:,'Error'};
    error = [error; error_tmp];
    running_tmp = ~finished_tmp & ~error_tmp;
    running = [running; running_tmp];

    Ind_finished_tmp = find(finished_tmp==1);
    Ind_finished = [Ind_finished; Ind_finished_tmp];

    Ind_error_tmp = find(error_tmp==1);
    Ind_error = [Ind_error; Ind_error_tmp];

    Ind_running_tmp = find(finished_tmp==0 & error_tmp==0);
    Ind_running = [Ind_running; Ind_running_tmp];

    fprintf("%s: Finished %s out of %s. %s ended with an error message. %s still running.\n",...
        RunTable(tt),...
        string(numel(Ind_finished_tmp)),string(numel(finished_tmp)),...
        string(numel(Ind_error_tmp)),...
        string(numel(Ind_running_tmp)));

end

Ind_finished = find(finished==1);
Ind_error = find(error==1);
Ind_running = find(running==1);

fprintf("TOTAL: Finished %s out of %s. %s ended with an error message. %s still running.\n",...
    string(numel(Ind_finished)),string(numel(finished)),string(numel(Ind_error)),...
    string(numel(finished)-numel(Ind_error)-numel(Ind_finished)));

H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
%counts(1,:)=histc(log10(gaA),linspace(min(log10(gaA)),max(log10(gaA)),50));
binrng = linspace(min((gaA)),max(log10(gaA)),50);
counts(1,:)=histc(log10(gaA(Ind_error)),binrng);
counts(2,:)=histc(log10(gaA(Ind_finished)),binrng);
counts(3,:)=histc(log10(gaA(Ind_running)),binrng);
b = bar(binrng,counts,'stacked','EdgeColor','none');
title("gaA (min: "+string(min(gaA))+", max: "+string(max(gaA))+")");
xlabel('log10(gaA)'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northeast');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_gaA_distribution";
print(H,fname,"-dpng","-r400");


H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
binrng = linspace(min(log10(gaC)),max(log10(gaC)),50);
counts(1,:)=histc(log10(gaC(Ind_error)),binrng);
counts(2,:)=histc(log10(gaC(Ind_finished)),binrng);
counts(3,:)=histc(log10(gaC(Ind_running)),binrng);
b = bar(binrng,counts,'stacked','EdgeColor','none');
title("gaC (min: "+string(min(gaC))+", max: "+string(max(gaC))+")");
xlabel('log10(gaC)'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northeast');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_gaC_distribution";
print(H,fname,"-dpng","-r400");

H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
binrng = linspace(min(log10(gsA)),max(log10(gsA)),50);
counts(1,:)=histc(log10(gsA(Ind_error)),binrng);
counts(2,:)=histc(log10(gsA(Ind_finished)),binrng);
counts(3,:)=histc(log10(gsA(Ind_running)),binrng);
b = bar(binrng,counts,'stacked','EdgeColor','none');
title("gsA (min: "+string(min(gsA))+", max: "+string(max(gsA))+")");
xlabel('log10(gsA)'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northwest');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_gsA_distribution";
print(H,fname,"-dpng","-r400");

H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
binrng = linspace(min(log10(gsC)),max(log10(gsC)),50);
counts(1,:)=histc(log10(gsC(Ind_error)),binrng);
counts(2,:)=histc(log10(gsC(Ind_finished)),binrng);
counts(3,:)=histc(log10(gsC(Ind_running)),binrng);
b = bar(binrng,counts,'stacked','EdgeColor','none');
title("gsC (min: "+string(min(gsC))+", max: "+string(max(gsC))+")");
xlabel('log10(gsC)'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northwest');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_gsC_distribution";
print(H,fname,"-dpng","-r400");

H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
binrng = linspace(min(m),max(m),50);
counts(1,:)=histc(m(Ind_error),binrng);
counts(2,:)=histc(m(Ind_finished),binrng);
counts(3,:)=histc(m(Ind_running),binrng);
b = bar(binrng,counts,'stacked','EdgeColor','none');
title("m (min: "+string(min(m))+", max: "+string(max(m))+")");
xlabel('m'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northwest');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_m_distribution";
print(H,fname,"-dpng","-r400");

H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
binrng = linspace(min(n),max(n),50);
counts(1,:)=histc(n(Ind_error),binrng);
counts(2,:)=histc(n(Ind_finished),binrng);
counts(3,:)=histc(n(Ind_running),binrng);
b = bar(binrng,counts,"stacked",'EdgeColor','none');
title("n (min: "+string(min(n))+", max: "+string(max(n))+")");
xlabel('n'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northwest');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_n_distribution";
print(H,fname,"-dpng","-r400");

H=fig('units','inches','width',50*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');
subplot("position",[0.1 0.1 0.85 0.85]); hold on;
binrng = linspace(min(dhdt),max(dhdt),50);
counts(1,:)=histc(dhdt(Ind_error),binrng);
counts(2,:)=histc(dhdt(Ind_finished),binrng);
counts(3,:)=histc(dhdt(Ind_running),binrng);
b = bar(binrng,counts,"stacked",'EdgeColor','none');
title("dhdt_{err} (min: "+string(min(dhdt))+", max: "+string(max(dhdt))+")");
xlabel('dhdt_{err}'); ylabel('count');
grid on; box on;
legend(b,["error","finished","running"],'Location','northwest');
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "./Figures/Ensemble_dhdt_err_distribution";
print(H,fname,"-dpng","-r400");